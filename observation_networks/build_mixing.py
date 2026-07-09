"""Build compatibility-network mixing matrices and assortativity summaries.

Pairwise compatibility edge scans can be large. For development, pass a small
window set, for example:

    python -m observation_networks.build_mixing --windows W080 W081

For the full Chapter 4 run:

    python -m observation_networks.build_mixing --all-windows --workers 5
"""

from __future__ import annotations

import argparse
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Sequence

import pandas as pd

from .lib.config import (
    DEFAULT_MIXING_ATTRIBUTES,
    INTERMEDIATE_DIR,
    PROJECT_ROOT,
)
from .lib.io import (
    ensure_results_dirs,
    load_chapter4_sequence_data,
    write_table,
)
from .lib.mixing import (
    build_degree_assortativity_for_edge_table,
    build_mixing_for_edge_table,
    specs_by_name,
)
from utils import load_pairwise_edges


LOGGER = logging.getLogger(__name__)

PAIRWISE_DATASET_DIR = PROJECT_ROOT / "data/processed/pairwise_distances_dataset"
INTERMEDIATE_TABLE_DIRS = {
    "matrix": INTERMEDIATE_DIR / "mixing_matrix",
    "summary": INTERMEDIATE_DIR / "comp_assortativity",
    "topology": INTERMEDIATE_DIR / "deg_assortativity",
}
FINAL_TABLES = {
    "matrix": "compatibility_mixing_matrix",
    "summary": "compatibility_assortativity",
    "topology": "compatibility_degree_assortativity",
}
PAIRWISE_COLUMNS = [
    "window_id",
    "pango_lineage",
    "id1",
    "id2",
    "epilink_compatibility",
]


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


def _window_from_pairwise_path(path: Path) -> str:
    return _normalise_window(path.stem.split("_", maxsplit=1)[0])


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


def _add_pairwise_stem(df: pd.DataFrame, stem: str) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["pairwise_stem"] = stem
    return out


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
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(path, index=False)


def _process_pairwise_file(
    pairwise_path: Path,
    nodes: pd.DataFrame,
    attributes,
    compatibility_threshold: float,
    missing_label: str | None,
    n_permutations: int,
    permutation_seed: int | None,
    force: bool,
) -> tuple[str, str]:
    """Build and write the mixing outputs for a single pairwise parquet file.

    Runs inside a worker process, so it must only use picklable arguments and
    module-level imports. Returns the pairwise file stem and a status string so
    the parent process can do all the logging with the configured handlers.
    """
    stem = pairwise_path.stem
    output_paths = _output_paths_for_stem(stem)
    if not force and _all_outputs_exist(output_paths):
        return stem, "skipped"

    edges = _read_pairwise_edges(
        pairwise_path,
        compatibility_threshold=compatibility_threshold,
    )
    if edges.empty:
        _write_missing_outputs(
            {
                "matrix": pd.DataFrame(),
                "summary": pd.DataFrame(),
                "topology": pd.DataFrame(),
            },
            output_paths,
            force=force,
        )
        return stem, "no_edges"

    topology = build_degree_assortativity_for_edge_table(
        edges,
        source_col="id1",
        target_col="id2",
        weight_col="epilink_compatibility",
        group_cols=("window_id", "pango_lineage"),
    )

    matrix, summary = build_mixing_for_edge_table(
        edges,
        nodes,
        attributes=attributes,
        node_id_col="sequence_id",
        source_col="id1",
        target_col="id2",
        weight_col="epilink_compatibility",
        group_cols=("window_id", "pango_lineage"),
        symmetric=True,
        missing_label=missing_label,
        n_permutations=n_permutations,
        seed=permutation_seed,
    )
    _write_missing_outputs(
        {
            "matrix": _add_pairwise_stem(matrix, stem),
            "summary": _add_pairwise_stem(summary, stem),
            "topology": _add_pairwise_stem(topology, stem),
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
        default=0.0001,
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
        "--n-permutations",
        "--permutations",
        dest="n_permutations",
        type=int,
        default=0,
        help=(
            "Number of vertex-label permutations for empirical p-values. "
            "Default 0 computes only observed assortativity."
        ),
    )
    parser.add_argument(
        "--permutation-seed",
        type=int,
        default=42,
        help="Base random seed for deterministic permutation p-values.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help=(
            "Number of worker processes. Default use 1. "
            "Use 1 to run serially (easier debugging/profiling)."
        ),
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
    ensure_results_dirs()
    _ensure_intermediate_dirs()

    if not args.all_windows and not args.windows:
        raise SystemExit("Specify --windows or --all-windows.")
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1.")
    if args.progress_every < 1:
        raise SystemExit("--progress-every must be at least 1.")
    if args.n_permutations < 0:
        raise SystemExit("--n-permutations must be non-negative.")
    pairwise_dir = args.pairwise_dir
    if not pairwise_dir.exists():
        raise SystemExit(f"Pairwise dataset directory not found: {pairwise_dir}")
    if args.n_permutations > 0:
        LOGGER.info(
            "Computing permutation p-values with %s permutations",
            f"{args.n_permutations:,}",
        )

    attributes = specs_by_name(args.attributes)
    attr_cols = [spec.column for spec in attributes]
    columns = ["window_id", "window_idx", "sequence_id", *attr_cols]
    sequence_df = load_chapter4_sequence_data(columns=columns, add_policy=False)

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

    # Build the per-file tasks up front. Each task corresponds to one
    # window-lineage pairwise parquet file.
    tasks: list[tuple] = []
    for pairwise_path in pairwise_files:
        window_id = _window_from_pairwise_path(pairwise_path)
        nodes = node_by_window.get(window_id)
        if nodes is None or nodes.empty:
            LOGGER.warning("Skipping %s: no node rows found", pairwise_path.stem)
            continue
        tasks.append(
            (
                pairwise_path,
                nodes,
                attributes,
                args.compatibility_threshold,
                args.missing_label,
                args.n_permutations,
                args.permutation_seed,
                args.force,
            )
        )

    total_tasks = len(tasks)
    if total_tasks == 0:
        raise SystemExit("No pairwise tasks could be built from the requested inputs.")

    LOGGER.info(
        "Processing compatibility mixing for %s pairwise files (%s)",
        f"{total_tasks:,}",
        _format_task_progress(0, total_tasks),
    )

    statuses: dict[str, int] = {"ok": 0, "skipped": 0, "no_edges": 0}

    def _record(stem: str, status: str) -> None:
        statuses[status] = statuses.get(status, 0) + 1
        if status == "no_edges":
            LOGGER.warning("Skipping %s: no compatibility edges found", stem)

    def _log_progress(processed: int) -> None:
        if processed % args.progress_every != 0 and processed != total_tasks:
            return
        LOGGER.info(
            "Processed %s pairwise files (%s)",
            _format_task_progress(processed, total_tasks),
            _format_status_counts(statuses),
        )

    if args.workers == 1:
        # Serial path.
        for processed, task in enumerate(tasks, start=1):
            pairwise_path = task[0]
            LOGGER.debug(
                "Processing compatibility mixing for %s (%s)",
                pairwise_path.stem,
                _format_task_progress(processed - 1, total_tasks),
            )
            _record(*_process_pairwise_file(*task))
            LOGGER.debug(
                "Finished compatibility mixing for %s (%s)",
                pairwise_path.stem,
                _format_task_progress(processed, total_tasks),
            )
            _log_progress(processed)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            future_to_stem = {
                executor.submit(_process_pairwise_file, *task): task[0].stem
                for task in tasks
            }
            for processed, future in enumerate(as_completed(future_to_stem), start=1):
                stem = future_to_stem[future]
                _record(*future.result())
                LOGGER.debug(
                    "Finished compatibility mixing for %s (%s)",
                    stem,
                    _format_task_progress(processed, total_tasks),
                )
                _log_progress(processed)

    LOGGER.info(
        "Pairwise mixing chunks complete: %s",
        _format_status_counts(statuses),
    )

    LOGGER.info("Writing compatibility mixing outputs")
    stems = [task[0].stem for task in tasks]
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
