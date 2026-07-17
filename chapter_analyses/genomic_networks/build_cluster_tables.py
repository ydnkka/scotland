"""Build core Chapter 4 observation/network analysis tables.

Run from the Scotland repository root:

    python -m chapter_analyses.genomic_networks.build_cluster_tables
"""

from __future__ import annotations

import argparse
import logging

from .lib.clusters import (
    build_cluster_attribute_composition,
    build_cluster_period_summary,
    build_cluster_table,
    build_cluster_window_summary,
)
from .lib.cohort import (
    build_clade_window_counts,
    build_cohort_summary,
    build_denominator_contrasts,
    build_sequence_composition,
    build_vaccination_context_by_policy,
    build_vaccination_window_context,
    build_window_coverage,
)
from .lib.io import ensure_results_dirs, load_sequence_data, write_table


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-transition",
        action="store_true",
        help=(
            "Deprecated no-op. Transition graph tables are now built by "
            "chapter_analyses.sse_detection.lib.sse.detection."
        ),
    )
    parser.add_argument(
        "--max-windows",
        type=int,
        default=None,
        help="Development cap on the number of main-analysis windows after loading.",
    )
    parser.add_argument(
        "--max-transition-windows",
        type=int,
        default=None,
        help=(
            "Deprecated no-op. Transition graph tables are now built by "
            "chapter_analyses.sse_detection.lib.sse.detection."
        ),
    )
    parser.add_argument(
        "--transition-window-stride",
        type=int,
        default=2,
        help=(
            "Deprecated no-op. Transition graph tables are now built by "
            "chapter_analyses.sse_detection.lib.sse.detection."
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

    LOGGER.info("Loading Chapter 4 sequence-window data")
    df = load_sequence_data()
    if args.max_windows is not None:
        keep = sorted(df["window_idx"].dropna().unique())[: args.max_windows]
        df = df.loc[df["window_idx"].isin(keep)].copy()
    LOGGER.info("Loaded %s sequence-window rows", f"{len(df):,}")

    window_coverage = build_window_coverage(df)
    cluster_table = build_cluster_table(df)

    table_jobs = {
        "cohort_summary": (build_cohort_summary(df), ("csv", "parquet")),
        "window_coverage": (window_coverage, ("csv", "parquet")),
        "window_denominator_contrasts": (
            build_denominator_contrasts(window_coverage),
            ("csv", "parquet"),
        ),
        "clade_window_counts": (build_clade_window_counts(df), ("parquet",)),
        "sequence_composition_by_policy": (
            build_sequence_composition(df, group_cols=("policy_period",)),
            ("parquet", "csv"),
        ),
        "vaccination_context_by_policy": (
            build_vaccination_context_by_policy(df),
            ("csv", "parquet"),
        ),
        "vaccination_window_context": (
            build_vaccination_window_context(df),
            ("csv", "parquet"),
        ),
        "cluster_table": (cluster_table, ("parquet",)),
        "cluster_window_summary": (
            build_cluster_window_summary(cluster_table),
            ("csv", "parquet"),
        ),
        "cluster_period_summary": (
            build_cluster_period_summary(cluster_table),
            ("csv", "parquet"),
        ),
        "cluster_attribute_composition": (
            build_cluster_attribute_composition(cluster_table),
            ("parquet", "csv"),
        ),
    }

    for name, (table, formats) in table_jobs.items():
        LOGGER.info("Writing %s (%s rows)", name, f"{len(table):,}")
        write_table(table, name, formats=formats)

    LOGGER.info(
        "Transition graph tables are built by "
        "python -m chapter_analyses.sse_detection.lib.sse.detection"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
