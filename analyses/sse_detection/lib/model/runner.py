"""CLI orchestration for Bayesian socio-geodemographic regression models.

The preparation and formula-building logic lives in ``model/prep.py``.
This module adds the repeatable fitting layer used by the domain-specific
``model/mixing.py`` and ``model/composition.py`` scripts.
"""

from __future__ import annotations

import argparse
import contextlib
import os
import sys
import time
import traceback
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import TextIO

import pandas as pd

from ..concurrent_io import (
    LockAlreadyHeldError,
    atomic_write_csv,
    exclusive_create_lock,
    exclusive_file_lock,
)
from ..sse.config import BAYESIAN_OUTPUT_DIR, PROJECT_ROOT, SSE_OUTPUT_DIR
from ..sse.io import load_sse_outputs
from .bayesian import (
    BayesianFitConfig,
    fit_prepared_model,
    print_diagnostic_report,
    save_prepared_model_result,
)
from .prep import (
    SCORE_OUTCOMES,
    Domain,
    PreparedModelFrame,
    PreparedRegressionRun,
    RegressionDataBundle,
    RegressionFamily,
    SampleSpec,
    default_model_specs,
    model_output_files,
    prepare_regression_data,
    prepare_regression_run,
)

ALL_FAMILIES: tuple[RegressionFamily, ...] = ("logistic", "linear")
ALL_OUTCOMES = ("candidate", *SCORE_OUTCOMES)


@dataclass(frozen=True)
class RegressionCliConfig:
    """Resolved command-line configuration for one domain script."""

    project_root: Path
    sse_output_dir: Path
    result_dir: Path
    domain: Domain
    families: tuple[RegressionFamily, ...]
    model_sets: tuple[str, ...] | None
    outcomes: tuple[str, ...] | None
    sample_rows: int | None
    sample_fraction: float | None
    positive_fraction: float | None
    no_sample: bool
    fit_config: BayesianFitConfig
    save_idata: bool = False
    dry_run: bool = False
    write_tables: bool = True
    skip_existing: bool = False
    continue_on_error: bool = False
    display_tables: bool = False
    print_diagnostics: bool = True
    log_file_name: str = "fit.log"
    live_progress: bool = False


def find_project_root(start: Path | str | None = None) -> Path:
    """Find the repository root from a file, directory, or cwd."""
    path = Path.cwd() if start is None else Path(start).expanduser()
    path = path.resolve()
    if path.is_file():
        path = path.parent

    for candidate in (path, *path.parents):
        if (candidate / "utils").exists() and (
            candidate / "analyses/sse_detection"
        ).exists():
            return candidate
    raise RuntimeError("Could not find the scotland repository root.")


def default_regression_result_dir(project_root: Path | str) -> Path:
    """Return the combined regression result directory used by the scripts."""
    relative = BAYESIAN_OUTPUT_DIR.relative_to(PROJECT_ROOT)
    return Path(project_root) / relative


def available_model_sets(domain: Domain) -> tuple[str, ...]:
    """Return valid model-set names for a regression domain."""
    return tuple(spec.model_set for spec in default_model_specs(domain))


def run_domain_models(config: RegressionCliConfig) -> pd.DataFrame:
    """Prepare, optionally fit, save, and manifest selected domain models."""
    _ensure_project_on_path(config.project_root)
    config.result_dir.mkdir(parents=True, exist_ok=True)

    data = load_regression_data(config.sse_output_dir)
    atomic_write_csv(
        data.eligibility_summary,
        config.result_dir / "eligibility_summary.csv",
        index=False,
    )

    print(f"Minimum candidate cluster size: {data.min_candidate_size}")
    print(_format_counts(data))

    selected_model_sets = resolve_model_sets(config.domain, config.model_sets)
    manifest_rows: list[dict[str, object]] = []
    selected_grid_rows: list[dict[str, object]] = []

    for family in config.families:
        outcomes = outcomes_for_family(family, config.outcomes)
        if not outcomes:
            print(f"Skipping {family}: no selected outcomes apply.")
            continue

        prepared = prepare_domain_regression_run(
            data,
            family=family,
            domain=config.domain,
            result_dir=config.result_dir / family,
            outcomes=outcomes,
            sample=sample_for_domain(data, config),
            write_tables=config.write_tables,
        )
        frames = [
            frame
            for frame in prepared.frames.values()
            if frame.model_set in selected_model_sets
        ]
        selected_grid_rows.extend(
            {
                **frame.grid_row(),
                "log_file": str(model_log_path(frame, config.log_file_name)),
            }
            for frame in frames
        )
        print(
            f"Prepared {len(frames)} {config.domain} {family} model(s): "
            + ", ".join(frame.result_key for frame in frames)
        )

        if config.dry_run:
            continue

        for frame in frames:
            try:
                manifest_rows.append(fit_save_and_log_frame(frame, config=config))
            except Exception:
                if manifest_rows:
                    write_saved_model_manifest(manifest_rows, config.result_dir)
                raise

    selected_grid = pd.DataFrame(selected_grid_rows)
    if not selected_grid.empty:
        atomic_write_csv(
            selected_grid,
            config.result_dir / f"{config.domain}_selected_model_grid.csv",
            index=False,
        )

    manifest = write_saved_model_manifest(manifest_rows, config.result_dir)
    if config.dry_run:
        print(
            "Dry run only. Selected model grid written to "
            f"{config.result_dir / f'{config.domain}_selected_model_grid.csv'}"
        )
    elif not manifest.empty:
        print(
            f"Saved model manifest to {config.result_dir / 'saved_model_manifest.csv'}"
        )
    return manifest


def load_regression_data(sse_output_dir: Path | str) -> RegressionDataBundle:
    """Load SSE outputs and align node/sequence regression frames."""
    return prepare_regression_data(load_sse_outputs(sse_output_dir))


def prepare_domain_regression_run(
    data: RegressionDataBundle,
    *,
    family: RegressionFamily,
    domain: Domain,
    result_dir: Path | str,
    outcomes: Sequence[str],
    sample: SampleSpec | None,
    write_tables: bool = True,
) -> PreparedRegressionRun:
    """Prepare model frames for one family/domain combination."""
    return prepare_regression_run(
        data,
        family=family,
        result_dir=result_dir,
        outcomes=outcomes,
        domains=(domain,),
        mixing_sample=sample if domain == "mixing" else None,
        composition_sample=sample if domain == "composition" else None,
        write_tables=write_tables,
    )


def sample_for_domain(
    data: RegressionDataBundle,
    config: RegressionCliConfig,
) -> SampleSpec | None:
    """Build the sample spec for the selected domain."""
    if config.no_sample:
        return None
    source = (
        data.eligible_nodes
        if config.domain == "mixing"
        else data.eligible_sequence_data
    )
    positive_fraction = config.positive_fraction
    if positive_fraction is None:
        positive_fraction = float(source["candidate"].mean())
    return SampleSpec(
        rows=config.sample_rows,
        fraction=config.sample_fraction,
        positive_fraction=positive_fraction,
        random_state=config.fit_config.random_seed,
    )


def runner_fit_config(config: RegressionCliConfig) -> BayesianFitConfig:
    """Return fit settings adjusted for the runner's backend-output policy."""
    return replace(
        config.fit_config,
        progressbar=config.live_progress,
        quiet=not config.live_progress,
    )


def resolve_model_sets(
    domain: Domain,
    requested: Sequence[str] | None,
) -> tuple[str, ...]:
    """Validate requested model-set names for a domain."""
    available = available_model_sets(domain)
    if requested is None:
        return available
    missing = [name for name in requested if name not in available]
    if missing:
        raise ValueError(
            f"Unknown {domain} model set(s): {', '.join(missing)}. "
            f"Available: {', '.join(available)}"
        )
    return tuple(requested)


def outcomes_for_family(
    family: RegressionFamily,
    requested: Sequence[str] | None,
) -> tuple[str, ...]:
    """Return selected outcomes that apply to one regression family."""
    defaults = ("candidate",) if family == "logistic" else SCORE_OUTCOMES
    if requested is None:
        return tuple(defaults)
    unknown = [outcome for outcome in requested if outcome not in ALL_OUTCOMES]
    if unknown:
        raise ValueError(
            f"Unknown outcome(s): {', '.join(unknown)}. "
            f"Available: {', '.join(ALL_OUTCOMES)}"
        )
    return tuple(outcome for outcome in requested if outcome in defaults)


def fit_save_and_log_frame(
    frame: PreparedModelFrame,
    *,
    config: RegressionCliConfig,
) -> dict[str, object]:
    """Fit one prepared frame, capture prints, save outputs, and return manifest."""
    log_path = model_log_path(frame, config.log_file_name)
    started_at = _now_iso()
    start = time.monotonic()
    try:
        with exclusive_create_lock(
            model_lock_path(frame),
            details=f"model={frame_key(frame)}\nlog={log_path}",
        ):
            if config.skip_existing and model_outputs_exist(
                frame,
                save_idata=config.save_idata,
            ):
                print(f"Skipping existing model: {frame_key(frame)}")
                return {
                    **base_manifest_row(frame, log_path=log_path),
                    "status": "skipped_existing",
                    "started_at": started_at,
                    "finished_at": _now_iso(),
                    "elapsed_seconds": time.monotonic() - start,
                }

            print(f"Fitting {frame_key(frame)}")
            print(f"  log: {log_path}")
            with log_path.open("w", encoding="utf-8") as log_handle:
                with (
                    contextlib.redirect_stdout(log_handle),
                    contextlib.redirect_stderr(log_handle),
                ):
                    write_model_log_header(frame, project_root=config.project_root)
                    if config.live_progress:
                        print("Backend stdout/stderr: terminal only (--live-progress).")
                    else:
                        print("Backend stdout/stderr: suppressed from fit.log.")
                    print()
                log_handle.flush()

            with suppress_backend_output(live=config.live_progress):
                result = fit_prepared_model(
                    frame,
                    config=runner_fit_config(config),
                    display_tables=False,
                    print_diagnostics=False,
                )

            with log_path.open("a", encoding="utf-8") as log_handle:
                with (
                    contextlib.redirect_stdout(log_handle),
                    contextlib.redirect_stderr(log_handle),
                ):
                    if config.print_diagnostics:
                        print_diagnostic_report(
                            result.diagnostics,
                            result.summary,
                            display_tables=config.display_tables,
                        )
                    manifest_row = save_prepared_model_result(
                        result,
                        frame,
                        save_idata=config.save_idata,
                    )
                    print()
                    print("Saved outputs")
                    for label, path in model_output_files(frame.output_dir).items():
                        if path.exists():
                            print(f"  {label}: {path}")
                log_handle.flush()
    except LockAlreadyHeldError as exc:
        elapsed = time.monotonic() - start
        row = {
            **base_manifest_row(frame, log_path=log_path),
            "status": "locked",
            "error": str(exc),
            "started_at": started_at,
            "finished_at": _now_iso(),
            "elapsed_seconds": elapsed,
        }
        if config.continue_on_error:
            print(f"Locked {frame_key(frame)}; {exc}")
            return row
        raise
    except Exception as exc:
        elapsed = time.monotonic() - start
        with log_path.open("a", encoding="utf-8") as log_handle:
            print(file=log_handle)
            print("Model failed", file=log_handle)
            traceback.print_exc(file=log_handle)
        row = {
            **base_manifest_row(frame, log_path=log_path),
            "status": "failed",
            "error": repr(exc),
            "started_at": started_at,
            "finished_at": _now_iso(),
            "elapsed_seconds": elapsed,
        }
        if config.continue_on_error:
            print(f"Failed {frame_key(frame)}; see {log_path}")
            return row
        raise

    elapsed = time.monotonic() - start
    print(f"Finished {frame_key(frame)} in {elapsed / 60:.1f} min")
    return {
        **manifest_row,
        "status": "fit",
        "log_file": str(log_path),
        "started_at": started_at,
        "finished_at": _now_iso(),
        "elapsed_seconds": elapsed,
    }


class _TeeTextIO:
    """Write text to multiple streams while preserving terminal-like behavior."""

    def __init__(self, *streams: TextIO) -> None:
        self._streams = streams
        self.encoding = getattr(streams[0], "encoding", None) if streams else None
        self.errors = getattr(streams[0], "errors", None) if streams else None

    def write(self, text: str) -> int:
        for stream in self._streams:
            stream.write(text)
        return len(text)

    def flush(self) -> None:
        for stream in self._streams:
            stream.flush()

    def isatty(self) -> bool:
        return any(
            getattr(stream, "isatty", lambda: False)() for stream in self._streams
        )


@contextlib.contextmanager
def redirect_model_output(log_path: Path, *, echo: bool = False) -> Iterator[TextIO]:
    """Redirect stdout and stderr to one per-model log file, optionally echoing live."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log_handle:
        stdout: TextIO = log_handle
        stderr: TextIO = log_handle
        if echo:
            stdout = _TeeTextIO(log_handle, sys.stdout)  # type: ignore[assignment]
            stderr = _TeeTextIO(log_handle, sys.stderr)  # type: ignore[assignment]
        with (
            contextlib.redirect_stdout(stdout),
            contextlib.redirect_stderr(stderr),
        ):
            yield log_handle


@contextlib.contextmanager
def suppress_backend_output(*, live: bool) -> Iterator[None]:
    """Show backend output live or suppress it so fit logs remain readable."""
    if live:
        yield
        return

    with (
        Path(os.devnull).open("w", encoding="utf-8") as sink,
        contextlib.redirect_stdout(sink),
        contextlib.redirect_stderr(sink),
    ):
        yield


def write_model_log_header(
    frame: PreparedModelFrame,
    *,
    project_root: Path | None = None,
) -> None:
    """Write model metadata before fitting so failures still have context."""
    print("=" * 100)
    print(frame_key(frame))
    print(frame.formula)
    print(f"Started: {_now_iso()}")
    print(f"Fit rows: {len(frame.fit_df):,}")
    print(f"Complete-case rows: {len(frame.full_df):,}")
    output_dir = frame.output_dir
    if project_root is not None:
        try:
            output_dir = output_dir.relative_to(project_root)
        except ValueError:
            pass
    print(f"Output dir: {output_dir}")
    print("=" * 100)
    print()


def model_outputs_exist(frame: PreparedModelFrame, *, save_idata: bool) -> bool:
    """Return True when the standard saved outputs already exist."""
    files = model_output_files(frame.output_dir)
    required = ("summary", "diagnostics", "metadata")
    if save_idata:
        required = (*required, "idata")
    return all(files[name].exists() for name in required)


def model_log_path(frame: PreparedModelFrame, log_file_name: str = "fit.log") -> Path:
    """Return the per-model log path."""
    return frame.output_dir / log_file_name


def model_lock_path(frame: PreparedModelFrame) -> Path:
    """Return the fail-fast lock path for one model output directory."""
    return frame.output_dir / ".fit.lock"


def frame_key(frame: PreparedModelFrame) -> str:
    """Return a stable family/domain/outcome/model identifier."""
    return f"{frame.family}:{frame.domain}:{frame.outcome}:{frame.model_set}"


def base_manifest_row(
    frame: PreparedModelFrame,
    *,
    log_path: Path,
) -> dict[str, object]:
    """Return model metadata shared by fit and skipped/failed manifest rows."""
    return {
        "family": frame.family,
        "domain": frame.domain,
        "outcome": frame.outcome,
        "model_set": frame.model_set,
        "model_dir": str(frame.output_dir),
        "n_rows": len(frame.fit_df),
        "candidate_rate": (
            float(frame.fit_df[frame.outcome].mean())
            if frame.family == "logistic"
            else None
        ),
        "outcome_mean": float(frame.fit_df[frame.outcome].mean()),
        "outcome_sd": (
            None
            if frame.family == "logistic"
            else float(frame.fit_df[frame.outcome].std())
        ),
        "use_sample": len(frame.fit_df) != len(frame.full_df),
        "log_file": str(log_path),
    }


def write_saved_model_manifest(
    manifest_rows: Sequence[dict[str, object]],
    result_dir: Path | str,
) -> pd.DataFrame:
    """Write current-run and accumulated manifests.

    Existing accumulated manifests are updated by
    family/domain/outcome/model_set, so running the mixing and composition
    scripts separately keeps both domains represented.
    """
    manifest = pd.DataFrame(manifest_rows)
    if manifest.empty:
        return manifest
    result_dir = Path(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)
    with exclusive_file_lock(result_dir / ".saved_model_manifest.lock"):
        atomic_write_csv(
            manifest,
            result_dir / "last_saved_model_manifest.csv",
            index=False,
        )

        combined = _merge_existing_manifest(
            result_dir / "saved_model_manifest.csv",
            manifest,
        )
        atomic_write_csv(
            combined,
            result_dir / "saved_model_manifest.csv",
            index=False,
        )
        for family, family_table in combined.groupby("family", sort=False):
            atomic_write_csv(
                family_table,
                result_dir / f"{family}_saved_model_manifest.csv",
                index=False,
            )
            family_dir = result_dir / str(family)
            if family_dir.exists():
                atomic_write_csv(
                    family_table,
                    family_dir / "saved_model_manifest.csv",
                    index=False,
                )
    return combined


def _merge_existing_manifest(path: Path, current: pd.DataFrame) -> pd.DataFrame:
    if not path.exists():
        return current
    existing = pd.read_csv(path)
    key_cols = ["family", "domain", "outcome", "model_set"]
    if not all(col in existing.columns and col in current.columns for col in key_cols):
        return pd.concat([existing, current], ignore_index=True)

    current_keys = set(map(tuple, current[key_cols].astype(str).to_numpy()))
    existing_keys = list(map(tuple, existing[key_cols].astype(str).to_numpy()))
    keep_existing = [key not in current_keys for key in existing_keys]
    return pd.concat([existing.loc[keep_existing], current], ignore_index=True)


def build_domain_arg_parser(domain: Domain) -> argparse.ArgumentParser:
    """Create the command-line parser for one domain script."""
    parser = argparse.ArgumentParser(
        description=f"Fit Bayesian SSE {domain} regression models.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--project-root", type=Path, default=None)
    parser.add_argument("--sse-output-dir", type=Path, default=None)
    parser.add_argument("--result-dir", type=Path, default=None)
    parser.add_argument(
        "--family",
        nargs="+",
        default=["all"],
        help="Regression family: logistic, linear, or all.",
    )
    parser.add_argument(
        "--model-set",
        "--model",
        dest="model_sets",
        nargs="+",
        default=["all"],
        help="Model-set name(s), comma-separated values, or all.",
    )
    parser.add_argument(
        "--outcome",
        nargs="+",
        default=["all"],
        help="Outcome name(s), comma-separated values, or all.",
    )
    parser.add_argument(
        "--sample-fraction",
        type=float,
        default=1.0,
        help="Fraction of complete-case rows to fit; ignored with --sample-rows.",
    )
    parser.add_argument("--sample-rows", type=int, default=None)
    parser.add_argument(
        "--positive-fraction",
        type=float,
        default=None,
        help="Target positive fraction for logistic samples; defaults to observed.",
    )
    parser.add_argument("--no-sample", action="store_true")
    parser.add_argument("--draws", type=int, default=2_000)
    parser.add_argument("--tune", type=int, default=2_000)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--target-accept", type=float, default=0.99)
    parser.add_argument("--random-seed", type=int, default=123)
    parser.add_argument("--inference-method", default="pymc")
    parser.add_argument("--fixed-prior-sigma", type=float, default=None)
    parser.add_argument("--intercept-prior-sigma", type=float, default=None)
    parser.add_argument("--random-effect-sigma", type=float, default=None)
    parser.add_argument("--residual-sigma", type=float, default=None)
    parser.add_argument("--no-log-likelihood", action="store_true")
    parser.add_argument("--centered", action="store_true")
    parser.add_argument("--save-idata", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument("--display-tables", action="store_true")
    parser.add_argument("--quiet-diagnostics", action="store_true")
    parser.add_argument("--no-write-tables", action="store_true")
    parser.add_argument("--log-file-name", default="fit.log")
    parser.add_argument(
        "--live-progress",
        action="store_true",
        help="Show backend model fit output in the terminal; keep fit.log clean.",
    )
    parser.add_argument(
        "--jax-platforms",
        default="cpu",
        help="Default JAX_PLATFORMS value set before fitting.",
    )
    parser.add_argument("--list-models", action="store_true")
    return parser


def main_for_domain(domain: Domain, argv: Sequence[str] | None = None) -> int:
    """Script entrypoint shared by mixing and composition modules."""
    parser = build_domain_arg_parser(domain)
    args = parser.parse_args(argv)

    if args.list_models:
        print_available_models(domain)
        return 0

    project_root = find_project_root(args.project_root)
    _ensure_project_on_path(project_root)
    if args.jax_platforms:
        os.environ.setdefault("JAX_PLATFORMS", args.jax_platforms)

    result_dir = _resolve_path(
        args.result_dir,
        base=project_root,
        default=default_regression_result_dir(project_root),
    )
    sse_output_dir = _resolve_path(
        args.sse_output_dir,
        base=project_root,
        default=project_root / SSE_OUTPUT_DIR.relative_to(PROJECT_ROOT),
    )
    model_sets = _parse_multi_arg(args.model_sets)
    outcomes = _parse_multi_arg(args.outcome)
    families = _resolve_families(_parse_multi_arg(args.family))
    _validate_sample_args(args)

    fit_config = BayesianFitConfig(
        draws=args.draws,
        tune=args.tune,
        chains=args.chains,
        cores=args.cores,
        target_accept=args.target_accept,
        inference_method=args.inference_method,
        random_seed=args.random_seed,
        fixed_prior_sigma=args.fixed_prior_sigma,
        intercept_prior_sigma=args.intercept_prior_sigma,
        random_effect_sigma=args.random_effect_sigma,
        residual_sigma=args.residual_sigma,
        log_likelihood=not args.no_log_likelihood,
        noncentered=not args.centered,
        progressbar=args.live_progress,
        quiet=not args.live_progress,
    )
    config = RegressionCliConfig(
        project_root=project_root,
        sse_output_dir=sse_output_dir,
        result_dir=result_dir,
        domain=domain,
        families=families,
        model_sets=model_sets,
        outcomes=outcomes,
        sample_rows=args.sample_rows,
        sample_fraction=None if args.sample_rows is not None else args.sample_fraction,
        positive_fraction=args.positive_fraction,
        no_sample=args.no_sample,
        fit_config=fit_config,
        save_idata=args.save_idata,
        dry_run=args.dry_run,
        write_tables=not args.no_write_tables,
        skip_existing=args.skip_existing,
        continue_on_error=args.continue_on_error,
        display_tables=args.display_tables,
        print_diagnostics=not args.quiet_diagnostics,
        log_file_name=args.log_file_name,
        live_progress=args.live_progress,
    )
    run_domain_models(config)
    return 0


def print_available_models(domain: Domain) -> None:
    """Print valid model selectors for one domain."""
    print(f"Available {domain} model sets:")
    for model_set in available_model_sets(domain):
        print(f"  - {model_set}")
    print()
    print("Families and outcomes:")
    print("  - logistic: candidate")
    print(f"  - linear: {', '.join(SCORE_OUTCOMES)}")


def _resolve_families(requested: Sequence[str] | None) -> tuple[RegressionFamily, ...]:
    if requested is None:
        return ALL_FAMILIES
    normalised = tuple(item.lower() for item in requested)
    if "all" in normalised:
        return ALL_FAMILIES
    allowed = set(ALL_FAMILIES)
    unknown = [item for item in normalised if item not in allowed]
    if unknown:
        raise ValueError(
            f"Unknown family value(s): {', '.join(unknown)}. "
            "Use logistic, linear, or all."
        )
    return normalised  # type: ignore[return-value]


def _parse_multi_arg(values: Sequence[str] | None) -> tuple[str, ...] | None:
    if not values:
        return None
    parsed: list[str] = []
    saw_all = False
    for value in values:
        for part in value.split(","):
            token = part.strip()
            if not token:
                continue
            if token.lower() == "all":
                saw_all = True
            else:
                parsed.append(token)
    if saw_all and not parsed:
        return None
    return tuple(dict.fromkeys(parsed)) if parsed else None


def _validate_sample_args(args: argparse.Namespace) -> None:
    if args.no_sample and (args.sample_rows is not None or args.sample_fraction != 1.0):
        raise ValueError("--no-sample cannot be combined with sample size settings.")
    if args.sample_rows is not None and args.sample_rows < 1:
        raise ValueError("--sample-rows must be at least 1.")
    if args.sample_fraction is not None and not 0 < args.sample_fraction <= 1:
        raise ValueError("--sample-fraction must be in (0, 1].")
    if args.positive_fraction is not None and not 0 < args.positive_fraction < 1:
        raise ValueError("--positive-fraction must be in (0, 1).")


def _resolve_path(
    value: Path | None,
    *,
    base: Path,
    default: Path,
) -> Path:
    path = default if value is None else value.expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _ensure_project_on_path(project_root: Path) -> None:
    project = str(project_root)
    if project not in sys.path:
        sys.path.insert(0, project)


def _format_counts(data: RegressionDataBundle) -> str:
    rows = []
    for row in data.eligibility_summary.to_dict("records"):
        rows.append(
            "{dataset}: {rows:,} rows, {candidates:,} candidates, "
            "{candidate_rate:.3%} candidate rate".format_map(row)
        )
    return "\n".join(rows)


def _now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


__all__ = [
    "RegressionCliConfig",
    "available_model_sets",
    "base_manifest_row",
    "build_domain_arg_parser",
    "default_regression_result_dir",
    "find_project_root",
    "fit_save_and_log_frame",
    "frame_key",
    "load_regression_data",
    "main_for_domain",
    "model_lock_path",
    "model_log_path",
    "model_outputs_exist",
    "outcomes_for_family",
    "prepare_domain_regression_run",
    "print_available_models",
    "redirect_model_output",
    "resolve_model_sets",
    "run_domain_models",
    "runner_fit_config",
    "sample_for_domain",
    "write_model_log_header",
    "write_saved_model_manifest",
]
