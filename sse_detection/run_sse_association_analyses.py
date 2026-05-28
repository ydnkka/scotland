"""Command-line runner for the SSE association analyses.

The default invocation runs the primary socio-geodemographic association
analysis plus the clade, window, observed-entropy, policy, and vaccination
sensitivity/context analyses.
Each run writes its CSV outputs to a dedicated subdirectory under
``sse_detection/results``.

Model rationale, inputs, output table definitions, and interpretation notes
live in ``sse_sociodemographic_association_notes.md``. Pass one or more
analysis names on the command line to run a subset.

Examples
--------
From the repository root:
    ``python -m sse_detection.run_sse_association_analyses``

From ``sse_detection/``:
    ``python -m run_sse_association_analyses main observed-entropy``
"""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
import sys
import time
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


DEFAULT_MODEL_METHOD = "firth_glm"
DEFAULT_MIXING_REFERENCE = "per 1 null-model SD increase in entropy"
_SSELIB: Any | None = None


def load_association_lib() -> Any:
    """Import the analysis library only when an analysis is about to run."""
    global _SSELIB
    if _SSELIB is None:
        from sse_detection import lib as sselib 

        _SSELIB = sselib
    return _SSELIB


@dataclass(frozen=True)
class AnalysisSpec:
    """Configuration for one analysis run."""

    key: str
    label: str
    result_subdir: str
    run: Callable[..., dict[str, Any]]


def result_dir_for(
    *,
    results_root: Path | None,
    result_subdir: str,
) -> Path:
    """Return the concrete output directory for one analysis."""
    root = (
        results_root
        if results_root is not None
        else PROJECT_ROOT / "sse_detection" / "results"
    )
    return root / result_subdir


def run_main_analysis(
    *,
    output_dir: Path | None,
    result_dir: Path,
    model_method: str,
    window_stride: int,
) -> dict[str, Any]:
    """Run the primary overall composition and mixing association analysis."""
    sselib = load_association_lib()
    return sselib.run_main_association_analysis(
        output_dir=output_dir,
        result_dir=result_dir,
        model_method=model_method,
        variant_adjuster="clade",
        window_adjustment="fixed_effects",
        mixing_reference=DEFAULT_MIXING_REFERENCE,
        window_stride=window_stride,
    )


def run_clade_sensitivity(
    *,
    output_dir: Path | None,
    result_dir: Path,
    model_method: str,
    window_stride: int,
) -> dict[str, Any]:
    """Run the clade-stratified fixed-window sensitivity analysis."""
    sselib = load_association_lib()
    model_sets = sselib.default_model_sets(
        variant_adjuster=None,
        window_adjustment="fixed_effects",
    )
    return sselib.run_association_pipeline(
        output_dir=output_dir,
        result_dir=result_dir,
        model_method=model_method,
        variant_adjuster=None,
        window_adjustment="fixed_effects",
        composition_model_sets=model_sets,
        mixing_model_sets=model_sets,
        group_by_clade=True,
        window_stride=window_stride,
    )


def run_window_sensitivity(
    *,
    output_dir: Path | None,
    result_dir: Path,
    model_method: str,
    window_stride: int,
) -> dict[str, Any]:
    """Run the clade-stratified surveillance-adjusted window sensitivity."""
    sselib = load_association_lib()
    model_sets = sselib.default_model_sets(
        variant_adjuster=None,
        window_adjustment="surveillance",
    )
    return sselib.run_association_pipeline(
        output_dir=output_dir,
        result_dir=result_dir,
        model_method=model_method,
        variant_adjuster=None,
        window_adjustment="surveillance",
        composition_model_sets=model_sets,
        mixing_model_sets=model_sets,
        group_by_clade=True,
        window_stride=window_stride,
    )


def run_observed_entropy_sensitivity(
    *,
    output_dir: Path | None,
    result_dir: Path,
    model_method: str,
    window_stride: int,
) -> dict[str, Any]:
    """Run the observed-normalised entropy mixing sensitivity analysis."""
    sselib = load_association_lib()
    mixing_model_sets = sselib.default_model_sets(
        variant_adjuster="clade",
        window_adjustment="fixed_effects",
    )
    return sselib.run_association_pipeline(
        output_dir=output_dir,
        result_dir=result_dir,
        model_method=model_method,
        variant_adjuster="clade",
        window_adjustment="fixed_effects",
        mixing_model_sets=mixing_model_sets,
        mixing_features=sselib.OBSERVED_MIXING_FEATURES_X10,
        mixing_reference=sselib.OBSERVED_MIXING_REFERENCE_X10,
        run_composition=False,
        run_mixing=True,
        window_stride=window_stride,
    )


def run_policy_analysis(
    *,
    output_dir: Path | None,
    result_dir: Path,
    model_method: str,
    window_stride: int,
) -> dict[str, Any]:
    """Run the policy-era association analysis."""
    sselib = load_association_lib()
    return sselib.run_policy_analysis(
        output_dir=output_dir,
        result_dir=result_dir,
        model_method=model_method,
        window_stride=window_stride,
    )


def run_vaccination_analysis(
    *,
    output_dir: Path | None,
    result_dir: Path,
    model_method: str,
    window_stride: int,
) -> dict[str, Any]:
    """Run the vaccination-context association analysis."""
    sselib = load_association_lib()
    return sselib.run_vaccination_analysis(
        output_dir=output_dir,
        result_dir=result_dir,
        model_method=model_method,
        window_stride=window_stride,
    )


ANALYSES: dict[str, AnalysisSpec] = {
    "main": AnalysisSpec(
        key="main",
        label="Main association analysis",
        result_subdir="association_outputs",
        run=run_main_analysis,
    ),
    "clade": AnalysisSpec(
        key="clade",
        label="Clade sensitivity analysis",
        result_subdir="sensitivity_clade",
        run=run_clade_sensitivity,
    ),
    "window": AnalysisSpec(
        key="window",
        label="Window sensitivity analysis",
        result_subdir="sensitivity_window",
        run=run_window_sensitivity,
    ),
    "observed-entropy": AnalysisSpec(
        key="observed-entropy",
        label="Observed entropy sensitivity analysis",
        result_subdir="sensitivity_observed_entropy",
        run=run_observed_entropy_sensitivity,
    ),
    "policy": AnalysisSpec(
        key="policy",
        label="Policy-era context analysis",
        result_subdir="policy_outputs",
        run=run_policy_analysis,
    ),
    "vaccination": AnalysisSpec(
        key="vaccination",
        label="Vaccination context analysis",
        result_subdir="vaccination_outputs",
        run=run_vaccination_analysis,
    ),
}
DEFAULT_ANALYSIS_ORDER = tuple(ANALYSES)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the main SSE socio-geodemographic association analysis and "
            "its sensitivity analyses."
        )
    )
    parser.add_argument(
        "analyses",
        nargs="*",
        default=None,
        metavar="{main,clade,window,observed-entropy,policy,vaccination,all}",
        help=(
            "Analyses to run. Defaults to all. Available analyses: "
            + ", ".join(ANALYSES.keys())
            + "."
        ),
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Repository root containing config.yaml.",
    )
    parser.add_argument(
        "--sse-output-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing SSE parquet outputs. Defaults to "
            "<project-root>/sse_detection/results/sse_outputs."
        ),
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=None,
        help=(
            "Parent directory for result subdirectories. Defaults to "
            "<project-root>/sse_detection/results."
        ),
    )
    parser.add_argument(
        "--model-method",
        default=DEFAULT_MODEL_METHOD,
        choices=["firth_glm", "glm_clustered", "conditional_logit_by_window"],
        help="Regression fitter to use for all selected analyses.",
    )
    parser.add_argument(
        "--window-stride",
        type=int,
        default=2,
        help="Rolling-window stride passed to sequence-level frame loading.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue with later analyses if one selected analysis raises.",
    )
    parser.add_argument(
        "--fail-on-model-failure",
        action="store_true",
        help="Exit non-zero when any model-level failures are recorded.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available analyses and exit.",
    )
    return parser.parse_args(argv)


def selected_analyses(values: Sequence[str] | None) -> tuple[str, ...]:
    allowed = {*ANALYSES.keys(), "all"}
    unknown = [value for value in values or [] if value not in allowed]
    if unknown:
        expected = ", ".join([*ANALYSES.keys(), "all"])
        raise ValueError(
            f"unknown analysis {unknown[0]!r}; expected one of: {expected}"
        )
    if not values:
        return DEFAULT_ANALYSIS_ORDER
    if "all" in values:
        if len(values) > 1:
            raise ValueError("'all' cannot be combined with named analyses.")
        return DEFAULT_ANALYSIS_ORDER
    return tuple(dict.fromkeys(values))


def print_available_analyses() -> None:
    for key in DEFAULT_ANALYSIS_ORDER:
        spec = ANALYSES[key]
        print(f"{key}: {spec.label} -> {spec.result_subdir}")


def print_result_summary(
    spec: AnalysisSpec,
    result: dict[str, Any],
    elapsed: float,
) -> None:
    result_dir = result["result_dir"]
    print(f"\nCompleted {spec.label} in {elapsed:.1f}s")
    print(f"Results saved to: {result_dir}")

    for filename, table in result["summary_tables"].items():
        print(f"  {filename}: {len(table):,} rows")

    diagnostics = result["cluster_diagnostics"]
    print(f"  cluster_diagnostics.csv: {len(diagnostics):,} rows")

    failures = result["failures"]
    if failures.empty:
        print("  model failures: none")
    else:
        print(f"  model_failures.csv: {len(failures):,} rows")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    if args.list:
        print_available_analyses()
        return 0
    
    output_dir = args.sse_output_dir.resolve() if args.sse_output_dir else None
    results_root = args.results_root.resolve() if args.results_root else None

    try:
        analysis_keys = selected_analyses(args.analyses)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    recorded_failures = False
    raised_failures: list[tuple[str, Exception]] = []

    for key in analysis_keys:
        spec = ANALYSES[key]
        result_dir = result_dir_for(
            results_root=results_root,
            result_subdir=spec.result_subdir,
        )

        print(f"\nRunning {spec.label}")
        print(f"Result directory: {result_dir}")

        started = time.perf_counter()
        try:
            result = spec.run(
                output_dir=output_dir,
                result_dir=result_dir,
                model_method=args.model_method,
                window_stride=args.window_stride,
            )
        except Exception as exc:
            print(f"FAILED {spec.label}: {exc}", file=sys.stderr)
            raised_failures.append((key, exc))
            if not args.keep_going:
                return 1
            continue

        elapsed = time.perf_counter() - started
        print_result_summary(spec, result, elapsed)
        recorded_failures = recorded_failures or not result["failures"].empty

    if raised_failures:
        print("\nAnalyses with uncaught failures:", file=sys.stderr)
        for key, exc in raised_failures:
            print(f"  {key}: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    if args.fail_on_model_failure and recorded_failures:
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
