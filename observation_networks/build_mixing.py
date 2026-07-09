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

import pandas as pd

from .lib.config import DEFAULT_MIXING_ATTRIBUTES
from .lib.io import (
    ensure_results_dirs,
    load_chapter4_sequence_data,
    load_pairwise_compatibility_edges,
    write_table,
)
from .lib.mixing import (
    build_degree_assortativity_for_edge_table,
    build_mixing_for_edge_table,
    specs_by_name,
)


LOGGER = logging.getLogger(__name__)


def _format_window_progress(processed: int, total: int) -> str:
    if total <= 0:
        return "0/0 windows (0.0%)"
    return f"{processed:,}/{total:,} windows ({processed / total:.1%})"


def _normalise_window(value: str) -> str:
    value = str(value).strip()
    upper = value.upper()
    if upper.startswith("W") and upper[1:].isdigit():
        return f"W{int(upper[1:]):03d}"
    if upper.isdigit():
        return f"W{int(upper):03d}"
    return value


def _process_window(
    window_id: str,
    nodes: pd.DataFrame,
    attributes,
    compatibility_threshold: float,
    missing_label: str | None,
    n_permutations: int,
    permutation_seed: int | None,
) -> tuple[str, pd.DataFrame | None, pd.DataFrame | None, pd.DataFrame | None]:
    """Build the mixing matrix/summary for a single window.

    Runs inside a worker process, so it must only use picklable arguments and
    module-level imports. Returns a status string plus the result frames so the
    parent process can do all the logging with the configured handlers.
    """
    edges = load_pairwise_compatibility_edges(
        windows=window_id,
        compatibility_threshold=compatibility_threshold,
    )
    if edges.empty:
        return "no_edges", None, None, None

    topology = build_degree_assortativity_for_edge_table(
        edges,
        source_col="id1",
        target_col="id2",
        weight_col="epilink_compatibility",
        group_cols=("window_id",),
    )

    matrix, summary = build_mixing_for_edge_table(
        edges,
        nodes,
        attributes=attributes,
        node_id_col="sequence_id",
        source_col="id1",
        target_col="id2",
        weight_col="epilink_compatibility",
        group_cols=("window_id",),
        symmetric=True,
        missing_label=missing_label,
        n_permutations=n_permutations,
        seed=permutation_seed,
    )
    return "ok", matrix, summary, topology


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
        help="EpiLink compatibility threshold passed to utils.load_pairwise_edges.",
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
        default=None,
        help=(
            "Number of worker processes. Default uses all CPUs. "
            "Use 1 to run serially (easier debugging/profiling)."
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

    if not args.all_windows and not args.windows:
        raise SystemExit("Specify --windows or --all-windows.")
    if args.n_permutations < 0:
        raise SystemExit("--n-permutations must be non-negative.")
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

    if args.max_windows is not None:
        windows = windows[: args.max_windows]

    # Build the per-window tasks up front (slice node rows once each).
    tasks: list[tuple] = []
    for window_id in windows:
        nodes = sequence_df.loc[sequence_df["window_id"].eq(window_id)].copy()
        if nodes.empty:
            LOGGER.warning("Skipping %s: no node rows found", window_id)
            continue
        tasks.append(
            (
                window_id,
                nodes,
                attributes,
                args.compatibility_threshold,
                args.missing_label,
                args.n_permutations,
                args.permutation_seed,
            )
        )

    results: dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]] = {}
    total_tasks = len(tasks)
    LOGGER.info(
        "Processing compatibility mixing for %s",
        _format_window_progress(0, total_tasks),
    )

    def _record(window_id: str, status: str, matrix, summary, topology) -> None:
        if status == "ok":
            results[window_id] = (matrix, summary, topology)
        elif status == "no_edges":
            LOGGER.warning("Skipping %s: no compatibility edges found", window_id)

    if args.workers == 1:
        # Serial path.
        for processed, task in enumerate(tasks, start=1):
            window_id = task[0]
            LOGGER.info(
                "Processing compatibility mixing for %s (%s)",
                window_id,
                _format_window_progress(processed - 1, total_tasks),
            )
            _record(window_id, *_process_window(*task))
            LOGGER.info(
                "Finished compatibility mixing for %s (%s)",
                window_id,
                _format_window_progress(processed, total_tasks),
            )
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            future_to_window = {
                executor.submit(_process_window, *task): task[0] for task in tasks
            }
            for processed, future in enumerate(as_completed(future_to_window), start=1):
                window_id = future_to_window[future]
                _record(window_id, *future.result())
                LOGGER.info(
                    "Finished compatibility mixing for %s (%s)",
                    window_id,
                    _format_window_progress(processed, total_tasks),
                )

    # Reassemble in the original window order for deterministic output.
    matrix_parts = [results[w][0] for w in windows if w in results]
    summary_parts = [results[w][1] for w in windows if w in results]
    topology_parts = [results[w][2] for w in windows if w in results]

    matrix_table = (
        pd.concat(matrix_parts, ignore_index=True, sort=False)
        if matrix_parts
        else pd.DataFrame()
    )
    summary_table = (
        pd.concat(summary_parts, ignore_index=True, sort=False)
        if summary_parts
        else pd.DataFrame()
    )
    topology_table = (
        pd.concat(topology_parts, ignore_index=True, sort=False)
        if topology_parts
        else pd.DataFrame()
    )

    LOGGER.info("Writing compatibility mixing outputs")
    write_table(matrix_table, "compatibility_mixing_matrix", formats=("parquet",))
    write_table(
        summary_table,
        "compatibility_assortativity",
        formats=("csv", "parquet"),
    )
    write_table(
        topology_table,
        "compatibility_degree_assortativity",
        formats=("csv", "parquet"),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
