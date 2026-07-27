"""Build compatibility-network mixing matrices and assortativity summaries.

Pairwise compatibility edge scans can be large. For development, pass a small
window set, for example:

    python -m chapter_analyses.genomic_networks.build_mixing --windows W080 W081

For the full Chapter 4 run:

    python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4 \
        --include-giants --giant-workers 1

Use --dry-run first to inspect the small-file and giant-file schedule.
Giant files are skipped by default.

This refactored version uses the NumPy-based assortativity implementation in
``lib.mixing`` and estimates uncertainty with multiplier bootstrap rather than
node jackknife. Analysis rows can be skipped/marked with NaN estimates when the
retained graph has fewer than ``--min-edges`` edges.
"""

from __future__ import annotations

import argparse
import logging
import os
import tempfile
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import pandas as pd

from utils import load_pairwise_edges

from .lib.config import (
    DEFAULT_MIXING_ATTRIBUTES,
    INTERMEDIATE_DIR,
    PROJECT_ROOT,
    SPARSIFICATION_THRESHOLD,
)
from .lib.io import (
    ensure_results_dirs,
    load_sequence_data,
    write_table,
)
from .lib.mixing import (
    build_degree_assortativity_for_edge_table,
    build_mixing_for_edge_table,
    specs_by_name,
)

LOGGER = logging.getLogger(__name__)

PAIRWISE_DATASET_DIR = PROJECT_ROOT / "data/processed/pairwise_distances_dataset"
EDGE_MANIFEST_PATH = (
    PROJECT_ROOT / "data/processed/sparsified_edge_counts_by_window_lineage.parquet"
)

# New analysis-file directories for the multiplier-bootstrap implementation.
# If you want to overwrite the original analysis names, change these back to:
#   mixing_matrix, comp_assortativity, deg_assortativity
INTERMEDIATE_TABLE_DIRS = {
    "matrix": INTERMEDIATE_DIR / "mixing_matrix_bootstrap",
    "summary": INTERMEDIATE_DIR / "comp_assortativity_bootstrap",
    "topology": INTERMEDIATE_DIR / "deg_assortativity_bootstrap",
}

FINAL_TABLES = {
    "matrix": "compatibility_mixing_matrix_bootstrap",
    "summary": "compatibility_assortativity_bootstrap",
    "topology": "compatibility_degree_assortativity_bootstrap",
}

PAIRWISE_COLUMNS = [
    "id1",
    "id2",
    "epilink_compatibility",
]


@dataclass(frozen=True)
class PairwiseMetadata:
    window_id: str
    pango_lineage: str
    sequence_count: int


@dataclass(frozen=True)
class ScheduledTask:
    args: tuple
    edge_cost: int | None
    edge_cost_source: str

    @property
    def pairwise_path(self) -> Path:
        return self.args[0]

    @property
    def stem(self) -> str:
        return self.pairwise_path.stem


def _format_task_progress(processed: int, total: int) -> str:
    if total <= 0:
        return "0/0 tasks (0.0%)"
    return f"{processed:,}/{total:,} tasks ({processed / total:.1%})"


def _format_status_counts(statuses: dict[str, int]) -> str:
    return (
        f"{statuses.get('ok', 0):,} ok, "
        f"{statuses.get('skipped', 0):,} skipped, "
        f"{statuses.get('no_edges', 0):,} no_edges"
    )


def _normalise_window(value: str) -> str:
    value = str(value).strip()
    upper = value.upper()
    if upper.startswith("W") and upper[1:].isdigit():
        return f"W{int(upper[1:]):03d}"
    if upper.isdigit():
        return f"W{int(upper):03d}"
    return value


def _pairwise_metadata_from_path(path: Path) -> PairwiseMetadata:
    # Stems are {window}_{lineage}_{count}; Pango lineages use dots, not
    # underscores, but join the middle fields so the parser is robust.
    parts = path.stem.split("_")
    if len(parts) < 3 or not parts[-1].isdigit():
        raise ValueError(
            f"Pairwise filename must look like {{window}}_{{lineage}}_{{count}}: "
            f"{path.name}"
        )

    lineage = "_".join(parts[1:-1])
    if not lineage:
        raise ValueError(
            f"Could not parse Pango lineage from pairwise filename: {path.name}"
        )

    return PairwiseMetadata(
        window_id=_normalise_window(parts[0]),
        pango_lineage=lineage,
        sequence_count=int(parts[-1]),
    )


def _window_from_pairwise_path(path: Path) -> str:
    return _pairwise_metadata_from_path(path).window_id


def _select_pairwise_files(
    pairwise_dir: Path,
    *,
    windows: Sequence[str],
    max_windows: int | None,
) -> list[Path]:
    files = sorted(pairwise_dir.glob("*.parquet"), key=lambda path: path.stem)

    if windows:
        wanted = set(windows)
        files = [path for path in files if _window_from_pairwise_path(path) in wanted]

    if max_windows is not None:
        retained_windows = sorted({_window_from_pairwise_path(path) for path in files})
        retained_windows = retained_windows[:max_windows]
        wanted = set(retained_windows)
        files = [path for path in files if _window_from_pairwise_path(path) in wanted]

    return files


def _output_paths_for_stem(stem: str) -> dict[str, Path]:
    return {
        kind: table_dir / f"{stem}.parquet"
        for kind, table_dir in INTERMEDIATE_TABLE_DIRS.items()
    }


def _all_outputs_exist(paths: dict[str, Path]) -> bool:
    return all(path.exists() for path in paths.values())


def _ensure_intermediate_dirs() -> None:
    for table_dir in INTERMEDIATE_TABLE_DIRS.values():
        table_dir.mkdir(parents=True, exist_ok=True)


def _read_pairwise_edges(
    pairwise_path: Path,
    *,
    compatibility_threshold: float | None,
) -> pd.DataFrame:
    return load_pairwise_edges(
        PAIRWISE_COLUMNS,
        compatibility_threshold=compatibility_threshold,
        pairwise_dataset=pairwise_path,
    )


def _with_pairwise_metadata(
    df: pd.DataFrame,
    metadata: PairwiseMetadata,
    stem: str,
    *,
    kind: str,
) -> pd.DataFrame:
    out = df.copy()
    out["window_id"] = metadata.window_id
    out["pango_lineage"] = metadata.pango_lineage
    out["pairwise_stem"] = stem

    metadata_cols = ["window_id", "pango_lineage", "pairwise_stem"]
    non_metadata = [col for col in out.columns if col not in metadata_cols]

    if kind == "topology":
        ordered = ["window_id", "pango_lineage", *non_metadata, "pairwise_stem"]
    elif kind == "summary" and "attribute_label" in non_metadata:
        insert_at = non_metadata.index("attribute_label") + 1
        ordered = [
            *non_metadata[:insert_at],
            "window_id",
            "pango_lineage",
            *non_metadata[insert_at:],
            "pairwise_stem",
        ]
    else:
        ordered = [*non_metadata, "window_id", "pango_lineage", "pairwise_stem"]

    return out[ordered]


def _empty_outputs_with_metadata(
    metadata: PairwiseMetadata,
    stem: str,
) -> dict[str, pd.DataFrame]:
    return {
        "matrix": _with_pairwise_metadata(
            pd.DataFrame(),
            metadata,
            stem,
            kind="matrix",
        ),
        "summary": _with_pairwise_metadata(
            pd.DataFrame(),
            metadata,
            stem,
            kind="summary",
        ),
        "topology": _with_pairwise_metadata(
            pd.DataFrame(),
            metadata,
            stem,
            kind="topology",
        ),
    }


def _write_parquet_atomic(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent,
            prefix=f".{path.stem}.",
            suffix=f".tmp{path.suffix}",
            delete=False,
        ) as tmp_file:
            tmp_path = Path(tmp_file.name)
            df.to_parquet(tmp_path, index=False)
        os.replace(tmp_path, path)
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink()


def _write_missing_outputs(
    outputs: dict[str, pd.DataFrame],
    paths: dict[str, Path],
    *,
    force: bool,
) -> None:
    for kind, df in outputs.items():
        path = paths[kind]
        if path.exists() and not force:
            continue
        _write_parquet_atomic(df, path)


def _process_pairwise_file(
    pairwise_path: Path,
    nodes: pd.DataFrame,
    attributes,
    compatibility_threshold: float,
    missing_label: str | None,
    bootstrap_replicates: int,
    bootstrap_alpha: float,
    bootstrap_seed: int,
    min_edges: int,
    force: bool,
) -> tuple[str, str]:
    """Build and write the mixing outputs for a single pairwise parquet file.

    Runs inside a worker process, so it must only use picklable arguments and
    module-level imports. Returns the pairwise file stem and a status string so
    the parent process can do all the logging with the configured handlers.
    """
    stem = pairwise_path.stem
    metadata = _pairwise_metadata_from_path(pairwise_path)
    output_paths = _output_paths_for_stem(stem)

    if not force and _all_outputs_exist(output_paths):
        return stem, "skipped"

    edges = _read_pairwise_edges(
        pairwise_path,
        compatibility_threshold=compatibility_threshold,
    )

    if edges.empty:
        _write_missing_outputs(
            _empty_outputs_with_metadata(metadata, stem),
            output_paths,
            force=force,
        )
        return stem, "no_edges"

    topology = build_degree_assortativity_for_edge_table(
        edges,
        nodes,
        node_id_col="sequence_id",
        source_col="id1",
        target_col="id2",
        weight_col="epilink_compatibility",
        group_cols=None,
        bootstrap_replicates=bootstrap_replicates,
        bootstrap_alpha=bootstrap_alpha,
        bootstrap_seed=bootstrap_seed,
        min_edges=min_edges,
    )

    matrix, summary = build_mixing_for_edge_table(
        edges,
        nodes,
        attributes=attributes,
        node_id_col="sequence_id",
        source_col="id1",
        target_col="id2",
        weight_col="epilink_compatibility",
        group_cols=None,
        symmetric=True,
        missing_label=missing_label,
        bootstrap_replicates=bootstrap_replicates,
        bootstrap_alpha=bootstrap_alpha,
        bootstrap_seed=bootstrap_seed,
        min_edges=min_edges,
    )

    _write_missing_outputs(
        {
            "matrix": _with_pairwise_metadata(matrix, metadata, stem, kind="matrix"),
            "summary": _with_pairwise_metadata(summary, metadata, stem, kind="summary"),
            "topology": _with_pairwise_metadata(
                topology,
                metadata,
                stem,
                kind="topology",
            ),
        },
        output_paths,
        force=force,
    )

    return stem, "ok"


def _sort_output_table(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    preferred = [
        "window_id",
        "pango_lineage",
        "pairwise_stem",
        "attribute",
        "source_category",
        "target_category",
    ]
    sort_cols = [col for col in preferred if col in df.columns]

    if not sort_cols:
        return df.reset_index(drop=True)

    return df.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)


def _concat_intermediate_table(kind: str, stems: Sequence[str]) -> pd.DataFrame:
    table_dir = INTERMEDIATE_TABLE_DIRS[kind]
    frames: list[pd.DataFrame] = []

    for stem in stems:
        path = table_dir / f"{stem}.parquet"
        if not path.exists():
            continue

        df = pd.read_parquet(path)
        if not df.empty:
            frames.append(df)

    if not frames:
        return pd.DataFrame()

    return _sort_output_table(pd.concat(frames, ignore_index=True, sort=False))


def _load_edge_costs(edge_manifest: Path) -> dict[str, int]:
    if not edge_manifest.exists():
        LOGGER.warning(
            "Edge manifest not found at %s; falling back to parquet file sizes",
            edge_manifest,
        )
        return {}

    try:
        manifest = pd.read_parquet(
            edge_manifest,
            columns=["pairwise_stem", "sparse_edges"],
        )
    except ValueError as exc:
        raise SystemExit(
            "Edge manifest must contain pairwise_stem and sparse_edges: "
            f"{edge_manifest}"
        ) from exc

    costs: dict[str, int] = {}

    for row in manifest.itertuples(index=False):
        sparse_edges = row.sparse_edges
        if pd.isna(sparse_edges):
            continue
        costs[str(row.pairwise_stem)] = int(cast(int, sparse_edges))

    return costs


def _edge_cost_for_pairwise_path(
    pairwise_path: Path,
    manifest_costs: dict[str, int],
) -> tuple[int | None, str]:
    manifest_cost = manifest_costs.get(pairwise_path.stem)
    if manifest_cost is not None:
        return manifest_cost, "manifest"

    # Missing manifest rows still need a deterministic cost; file size is a
    # conservative local proxy, and stat failures are treated as giant later.
    try:
        return pairwise_path.stat().st_size, "file_size"
    except OSError:
        return None, "unknown"


def _is_giant_task(task: ScheduledTask, threshold: int) -> bool:
    return task.edge_cost is None or task.edge_cost >= threshold


def _descending_edge_cost(task: ScheduledTask) -> float:
    return float("inf") if task.edge_cost is None else float(task.edge_cost)


def _format_edge_cost(edge_cost: int | None) -> str:
    return "unknown" if edge_cost is None else f"{edge_cost:,}"


def _phase_edge_summary(tasks: Sequence[ScheduledTask]) -> str:
    known_total = sum(task.edge_cost for task in tasks if task.edge_cost is not None)
    unknown = sum(1 for task in tasks if task.edge_cost is None)
    suffix = f", {unknown:,} unknown" if unknown else ""
    return f"{known_total:,} summed edges/cost{suffix}"


def _print_dry_run_schedule(
    small_tasks: Sequence[ScheduledTask],
    giant_tasks: Sequence[ScheduledTask],
    *,
    workers: int,
    giant_workers: int,
    giant_threshold: int,
    include_giants: bool,
) -> None:
    total = len(small_tasks) + len(giant_tasks)
    processing_total = len(small_tasks) + (len(giant_tasks) if include_giants else 0)

    print("Dry run: no pairwise files will be processed.")
    print(f"Total selected files: {total:,}")
    print(f"Files that would be processed: {processing_total:,}")
    print(
        "Small phase: "
        f"{len(small_tasks):,} files, {_phase_edge_summary(small_tasks)}, "
        f"workers={workers:,}"
    )

    if include_giants:
        print(
            "Giant phase: "
            f"{len(giant_tasks):,} files, {_phase_edge_summary(giant_tasks)}, "
            f"workers={giant_workers:,}, threshold={giant_threshold:,}"
        )
    else:
        print(
            "Giant phase: "
            f"{len(giant_tasks):,} files skipped by default, "
            f"{_phase_edge_summary(giant_tasks)}, threshold={giant_threshold:,}; "
            "pass --include-giants to process them"
        )

    print("Giants in processing order:")
    if not giant_tasks:
        print("  (none)")
        return

    for task in giant_tasks:
        print(
            "  "
            f"{task.stem}\t{_format_edge_cost(task.edge_cost)}\t"
            f"{task.edge_cost_source}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--windows",
        nargs="*",
        help="Specific windows to process, e.g. W080 W081 or 80 81.",
    )
    parser.add_argument(
        "--all-windows",
        action="store_true",
        help="Process every available retained analysis window.",
    )
    parser.add_argument(
        "--max-windows",
        type=int,
        default=None,
        help="Optional cap for development runs after window selection.",
    )
    parser.add_argument(
        "--attributes",
        nargs="*",
        default=None,
        help=(
            "Attribute names to process. Defaults to all: "
            + ", ".join(spec.name for spec in DEFAULT_MIXING_ATTRIBUTES)
        ),
    )
    parser.add_argument(
        "--compatibility-threshold",
        type=float,
        default=SPARSIFICATION_THRESHOLD,
        help="Retain pairwise edges with epilink_compatibility above this threshold.",
    )
    parser.add_argument(
        "--pairwise-dir",
        type=Path,
        default=PAIRWISE_DATASET_DIR,
        help=(
            "Directory of per-window-lineage pairwise parquet files. "
            f"Default: {PAIRWISE_DATASET_DIR}"
        ),
    )
    parser.add_argument(
        "--missing-label",
        default=None,
        help="Optional label for missing node attributes. Default drops missing pairs.",
    )
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=500,
        help=(
            "Number of multiplier-bootstrap replicates for assortativity uncertainty. "
            "Use 0 to skip uncertainty. Default: 500."
        ),
    )
    parser.add_argument(
        "--bootstrap-alpha",
        type=float,
        default=0.05,
        help="Bootstrap interval alpha. Default 0.05 gives a 95% percentile interval.",
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=123,
        help="Base random seed for multiplier bootstrap. Default: 123.",
    )
    parser.add_argument(
        "--min-edges",
        type=int,
        default=0,
        help=(
            "Minimum number of retained compatibility edges required before reporting "
            "assortativity estimates. Analyses below this threshold are retained with "
            "NaN estimates and skipped_reason. Default: 0."
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help=(
            "Worker processes for the small-file phase. Default use 1. "
            "Use 1 to run serially, which is easier for debugging/profiling."
        ),
    )
    parser.add_argument(
        "--giant-threshold",
        type=int,
        default=50_000_000,
        help=(
            "Files with cost >= this many sparse edges are processed in the "
            "giant phase. Manifest misses fall back to file size; unknown costs "
            "are treated as giant. Default: 50,000,000."
        ),
    )
    parser.add_argument(
        "--giant-workers",
        type=int,
        default=1,
        help=(
            "Worker processes for the giant-file phase when --include-giants "
            "is passed. Default: 1."
        ),
    )
    parser.add_argument(
        "--include-giants",
        action="store_true",
        help=(
            "Process files at or above --giant-threshold. By default these "
            "memory-heavy files are skipped."
        ),
    )
    parser.add_argument(
        "--edge-manifest",
        type=Path,
        default=EDGE_MANIFEST_PATH,
        help=(
            "Parquet with pairwise_stem and sparse_edges columns used for scheduling. "
            f"Default: {EDGE_MANIFEST_PATH}"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the edge-count-aware schedule and exit without processing files.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute per-file intermediate outputs even when they already exist.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help=(
            "Log aggregate INFO progress every N completed pairwise files. "
            "Per-file messages are DEBUG. Default: 100."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=args.log_level, format="%(levelname)s: %(message)s")

    if not args.all_windows and not args.windows:
        raise SystemExit("Specify --windows or --all-windows.")

    if args.workers < 1:
        raise SystemExit("--workers must be at least 1.")

    if args.giant_workers < 1:
        raise SystemExit("--giant-workers must be at least 1.")

    if args.giant_threshold < 0:
        raise SystemExit("--giant-threshold must be non-negative.")

    if args.progress_every < 1:
        raise SystemExit("--progress-every must be at least 1.")

    if args.bootstrap_replicates < 0:
        raise SystemExit("--bootstrap-replicates must be non-negative.")

    if not 0 < args.bootstrap_alpha < 1:
        raise SystemExit("--bootstrap-alpha must be between 0 and 1.")

    if args.min_edges < 0:
        raise SystemExit("--min-edges must be non-negative.")

    pairwise_dir = args.pairwise_dir
    if not pairwise_dir.exists():
        raise SystemExit(f"Pairwise dataset directory not found: {pairwise_dir}")

    if args.bootstrap_replicates > 0:
        LOGGER.info(
            "Computing multiplier-bootstrap uncertainty with %s replicates, alpha=%s",
            f"{args.bootstrap_replicates:,}",
            args.bootstrap_alpha,
        )
    else:
        LOGGER.info("Skipping multiplier-bootstrap uncertainty.")

    if args.min_edges > 0:
        LOGGER.info(
            "Applying minimum-edge filter: at least %s retained edges per analysis",
            f"{args.min_edges:,}",
        )

    attributes = specs_by_name(args.attributes)
    attr_cols = [spec.column for spec in attributes]

    columns = ["window_id", "window_idx", "sequence_id", *attr_cols]
    sequence_df = load_sequence_data(columns=columns)

    if args.all_windows:
        windows = sorted(sequence_df["window_id"].dropna().unique())
    else:
        windows = [_normalise_window(window) for window in args.windows]

    pairwise_files = _select_pairwise_files(
        pairwise_dir,
        windows=windows,
        max_windows=args.max_windows,
    )

    if not pairwise_files:
        raise SystemExit("No pairwise parquet files matched the requested windows.")

    node_by_window = {
        window_id: group.copy()
        for window_id, group in sequence_df.groupby("window_id", sort=False)
    }

    edge_costs = _load_edge_costs(args.edge_manifest)

    # Build the per-file tasks up front. Each task corresponds to one
    # window-lineage pairwise parquet file.
    tasks: list[ScheduledTask] = []

    for pairwise_path in pairwise_files:
        window_id = _window_from_pairwise_path(pairwise_path)
        nodes = node_by_window.get(window_id)

        if nodes is None or nodes.empty:
            LOGGER.warning("Skipping %s: no node rows found", pairwise_path.stem)
            continue

        edge_cost, edge_cost_source = _edge_cost_for_pairwise_path(
            pairwise_path,
            edge_costs,
        )

        tasks.append(
            ScheduledTask(
                args=(
                    pairwise_path,
                    nodes,
                    attributes,
                    args.compatibility_threshold,
                    args.missing_label,
                    args.bootstrap_replicates,
                    args.bootstrap_alpha,
                    args.bootstrap_seed,
                    args.min_edges,
                    args.force,
                ),
                edge_cost=edge_cost,
                edge_cost_source=edge_cost_source,
            )
        )

    total_tasks = len(tasks)
    if total_tasks == 0:
        raise SystemExit("No pairwise tasks could be built from the requested inputs.")

    small_tasks = [
        task for task in tasks if not _is_giant_task(task, args.giant_threshold)
    ]
    giant_tasks = sorted(
        (task for task in tasks if _is_giant_task(task, args.giant_threshold)),
        key=_descending_edge_cost,
        reverse=True,
    )

    if args.dry_run:
        _print_dry_run_schedule(
            small_tasks,
            giant_tasks,
            workers=args.workers,
            giant_workers=args.giant_workers,
            giant_threshold=args.giant_threshold,
            include_giants=args.include_giants,
        )
        return 0

    if giant_tasks and not args.include_giants:
        LOGGER.info(
            "Skipping %s giant pairwise files at threshold %s; pass "
            "--include-giants to process them",
            f"{len(giant_tasks):,}",
            f"{args.giant_threshold:,}",
        )

    process_tasks = [*small_tasks, *giant_tasks] if args.include_giants else small_tasks
    total_process_tasks = len(process_tasks)

    if total_process_tasks == 0:
        raise SystemExit(
            "No pairwise tasks left after skipping giant files. "
            "Pass --include-giants to process the selected giant files."
        )

    ensure_results_dirs()
    _ensure_intermediate_dirs()

    LOGGER.info(
        "Processing compatibility mixing for %s pairwise files (%s)",
        f"{total_process_tasks:,}",
        _format_task_progress(0, total_process_tasks),
    )

    statuses: dict[str, int] = {"ok": 0, "skipped": 0, "no_edges": 0}

    def _record(stem: str, status: str) -> None:
        statuses[status] = statuses.get(status, 0) + 1
        if status == "no_edges":
            LOGGER.warning("Skipping %s: no compatibility edges found", stem)

    def _log_progress(processed: int) -> None:
        is_progress_interval = processed % args.progress_every == 0
        if not is_progress_interval and processed != total_process_tasks:
            return

        LOGGER.info(
            "Processed %s pairwise files (%s)",
            _format_task_progress(processed, total_process_tasks),
            _format_status_counts(statuses),
        )

    processed_total = 0

    def _run_phase(
        phase_name: str,
        phase_tasks: Sequence[ScheduledTask],
        *,
        workers: int,
    ) -> None:
        nonlocal processed_total

        if not phase_tasks:
            LOGGER.info("Skipping %s phase: no pairwise files", phase_name)
            return

        LOGGER.info(
            "Starting %s phase: %s pairwise files, workers=%s",
            phase_name,
            f"{len(phase_tasks):,}",
            f"{workers:,}",
        )

        if workers == 1:
            for task in phase_tasks:
                pairwise_path = task.pairwise_path

                LOGGER.debug(
                    "Processing compatibility mixing for %s (%s)",
                    pairwise_path.stem,
                    _format_task_progress(processed_total, total_process_tasks),
                )

                _record(*_process_pairwise_file(*task.args))
                processed_total += 1

                LOGGER.debug(
                    "Finished compatibility mixing for %s (%s)",
                    pairwise_path.stem,
                    _format_task_progress(processed_total, total_process_tasks),
                )

                _log_progress(processed_total)
            return

        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_stem = {
                executor.submit(_process_pairwise_file, *task.args): task.stem
                for task in phase_tasks
            }

            for future in as_completed(future_to_stem):
                stem = future_to_stem[future]

                _record(*future.result())
                processed_total += 1

                LOGGER.debug(
                    "Finished compatibility mixing for %s (%s)",
                    stem,
                    _format_task_progress(processed_total, total_process_tasks),
                )

                _log_progress(processed_total)

    # Run small files first, then the memory-heavy files. This avoids clumping
    # the largest BA.2-era parquets on 16 GB machines using macOS spawn.
    _run_phase("small", small_tasks, workers=args.workers)

    if args.include_giants:
        _run_phase("giant", giant_tasks, workers=args.giant_workers)

    LOGGER.info(
        "Pairwise mixing chunks complete: %s",
        _format_status_counts(statuses),
    )

    LOGGER.info("Writing compatibility mixing outputs")
    stems = [task.stem for task in process_tasks]

    matrix_table = _concat_intermediate_table("matrix", stems)
    summary_table = _concat_intermediate_table("summary", stems)
    topology_table = _concat_intermediate_table("topology", stems)

    write_table(matrix_table, FINAL_TABLES["matrix"], formats=("parquet",))
    write_table(
        summary_table,
        FINAL_TABLES["summary"],
        formats=("csv", "parquet"),
    )
    write_table(
        topology_table,
        FINAL_TABLES["topology"],
        formats=("csv", "parquet"),
    )

    LOGGER.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
