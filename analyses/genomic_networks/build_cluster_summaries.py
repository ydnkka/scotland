"""Build all genomic-network cohort, context, and cluster summary tables.

Run from the Scotland repository root:

    python -m analyses.genomic_networks.build_cluster_summaries
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .lib.cluster_distance_summary import (
    MIN_PAIRWISE_ROWS,
    build_cluster_pairwise_distance_overall_summary,
)
from .lib.cluster_pairwise_distances import (
    build_cluster_pairwise_distance_summary,
)
from .lib.clusters import (
    build_cluster_period_summary,
    build_cluster_table,
    build_cluster_window_summary,
)
from .lib.cohort import (
    build_cohort_summary,
    build_denominator_contrasts,
    build_sequence_composition,
    build_test_reason_by_policy_era,
    build_vaccination_window_context,
    build_window_coverage,
)
from .lib.config import TABLES_DIR
from .lib.io import (
    TABLE_OUTPUT_FORMATS,
    ensure_results_dirs,
    load_sequence_data,
    read_table,
    write_table,
)
from .lib.windows import normalise_windows

LOGGER = logging.getLogger(__name__)

PAIRWISE_SUMMARY_NAME = "cluster_pairwise_genetic_temporal_distance"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--max-windows",
        type=int,
        default=None,
        help=("Development cap on the number of main-analysis windows after loading."),
    )

    parser.add_argument(
        "--windows",
        nargs="*",
        help=(
            "Optional windows for pairwise-distance summaries, "
            "for example W080 W081. Default processes all windows."
        ),
    )

    parser.add_argument(
        "--table-dir",
        type=Path,
        default=TABLES_DIR,
        help=f"Input/output table directory. Default: {TABLES_DIR}.",
    )

    parser.add_argument(
        "--reuse-input-tables",
        action="store_true",
        help=(
            "Reuse existing cluster_table, window_coverage, and "
            "cluster_pairwise_distance_summary tables instead of "
            "rebuilding those tables."
        ),
    )

    parser.add_argument(
        "--skip-pairwise-distances",
        action="store_true",
        help=(
            "Build cohort and cluster summaries, but skip all "
            "pairwise-derived summaries."
        ),
    )

    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    )

    args = parser.parse_args()

    if args.reuse_input_tables and args.max_windows is not None:
        parser.error(
            "--max-windows cannot be combined with "
            "--reuse-input-tables because the reused tables may contain "
            "a different set of windows."
        )

    return args


def load_main_data(max_windows: int | None = None):
    """Load sequence-window data and optionally apply a window cap."""
    LOGGER.info("Loading sequence-window data")
    df = load_sequence_data()

    if max_windows is not None:
        keep = sorted(df["window_idx"].dropna().unique())[:max_windows]
        df = df.loc[df["window_idx"].isin(keep)].copy()

    LOGGER.info("Loaded %s sequence-window rows", f"{len(df):,}")
    return df


def write_table_jobs(table_jobs, table_dir: Path) -> None:
    """Write a collection of named tables."""
    for name, (table, formats) in table_jobs.items():
        LOGGER.info(
            "Writing %s (%s rows)",
            name,
            f"{len(table):,}",
        )
        write_table(
            table,
            name,
            table_dir=table_dir,
            formats=formats,
        )


def main() -> int:
    args = parse_args()

    logging.basicConfig(
        level=args.log_level,
        format="%(levelname)s: %(message)s",
    )

    ensure_results_dirs()
    args.table_dir.mkdir(parents=True, exist_ok=True)

    # The cohort tables always need sequence-window data.
    df = load_main_data(args.max_windows)

    if args.reuse_input_tables:
        LOGGER.info("Reading existing cluster_table and window_coverage")

        cluster_table = read_table(
            "cluster_table",
            table_dir=args.table_dir,
        )

        window_coverage = read_table(
            "window_coverage",
            table_dir=args.table_dir,
        )
    else:
        LOGGER.info("Building cluster_table and window_coverage")
        window_coverage = build_window_coverage(df)
        cluster_table = build_cluster_table(df)

    # ------------------------------------------------------------------
    # Cohort and context tables
    # ------------------------------------------------------------------
    cohort_table_jobs = {
        "cohort_summary": (
            build_cohort_summary(df),
            TABLE_OUTPUT_FORMATS,
        ),
        "window_coverage": (
            window_coverage,
            TABLE_OUTPUT_FORMATS,
        ),
        "window_denominator_contrasts": (
            build_denominator_contrasts(window_coverage),
            TABLE_OUTPUT_FORMATS,
        ),
        "sequence_composition_by_policy": (
            build_sequence_composition(df),
            TABLE_OUTPUT_FORMATS,
        ),
        "test_reason_by_policy_era": (
            build_test_reason_by_policy_era(df),
            TABLE_OUTPUT_FORMATS,
        ),
        "vaccination_window_context": (
            build_vaccination_window_context(df),
            TABLE_OUTPUT_FORMATS,
        ),
    }

    write_table_jobs(
        cohort_table_jobs,
        table_dir=args.table_dir,
    )

    # ------------------------------------------------------------------
    # Cluster tables
    # ------------------------------------------------------------------
    cluster_period = build_cluster_period_summary(cluster_table)

    cluster_table_jobs = {
        "cluster_table": (
            cluster_table,
            ("parquet",),
        ),
        "cluster_window_summary": (
            build_cluster_window_summary(cluster_table),
            TABLE_OUTPUT_FORMATS,
        ),
        "cluster_period_summary": (
            cluster_period,
            TABLE_OUTPUT_FORMATS,
        ),
    }

    write_table_jobs(
        cluster_table_jobs,
        table_dir=args.table_dir,
    )

    # ------------------------------------------------------------------
    # Pairwise-derived cluster tables
    # ------------------------------------------------------------------
    if args.skip_pairwise_distances:
        LOGGER.info("Skipping pairwise-derived cluster summaries")
        return 0

    if args.reuse_input_tables:
        LOGGER.info(
            "Reading existing %s",
            PAIRWISE_SUMMARY_NAME,
        )
        pairwise_summary = read_table(
            PAIRWISE_SUMMARY_NAME,
            table_dir=args.table_dir,
        )
    else:
        LOGGER.info("Building pairwise distance summary")
        pairwise_summary = build_cluster_pairwise_distance_summary(
            windows=normalise_windows(args.windows),
        )

    pairwise_overall = build_cluster_pairwise_distance_overall_summary(
        pairwise_summary,
        min_pairwise_rows=MIN_PAIRWISE_ROWS,
    )

    pairwise_table_jobs = {
        PAIRWISE_SUMMARY_NAME: (
            pairwise_summary,
            TABLE_OUTPUT_FORMATS,
        ),
        "cluster_pairwise_distance_summary": (
            pairwise_overall,
            TABLE_OUTPUT_FORMATS,
        ),
    }

    write_table_jobs(
        pairwise_table_jobs,
        table_dir=args.table_dir,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
