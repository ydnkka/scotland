"""Build core Chapter 4 observation/network analysis tables.

Run from the Scotland repository root:

    python -m observation_networks.build_tables
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
    build_window_coverage,
)
from .lib.config import TRANSITION_WINDOW_STRIDE
from .lib.io import ensure_results_dirs, load_chapter4_sequence_data, write_table
from .lib.mixing import build_transition_mixing
from .lib.transition_graph import build_transition_outputs


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-transition",
        action="store_true",
        help="Skip alternate-window transition-graph table construction.",
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
        help="Development cap on alternate-window transition input after loading.",
    )
    parser.add_argument(
        "--transition-window-stride",
        type=int,
        default=TRANSITION_WINDOW_STRIDE,
        help="Window stride for the transition graph input. Default: 2.",
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
    df = load_chapter4_sequence_data()
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

    if args.skip_transition:
        return 0

    LOGGER.info(
        "Loading alternate-window data for transition graph (stride=%s)",
        args.transition_window_stride,
    )
    transition_df = load_chapter4_sequence_data(
        window_stride=args.transition_window_stride,
    )
    if args.max_transition_windows is not None:
        keep = sorted(transition_df["window_idx"].dropna().unique())[
            : args.max_transition_windows
        ]
        transition_df = transition_df.loc[transition_df["window_idx"].isin(keep)].copy()
    LOGGER.info("Building transition graph baseline tables")
    transition = build_transition_outputs(transition_df)
    transition_matrix, transition_assortativity = build_transition_mixing(
        transition.edge_table,
        transition.node_table,
    )

    transition_jobs = {
        "transition_edge_table": (transition.edge_table, ("parquet",)),
        "transition_node_table": (transition.node_table, ("parquet",)),
        "transition_graph_summary": (transition.graph_summary, ("csv", "parquet")),
        "transition_window_summary": (transition.window_summary, ("csv", "parquet")),
        "transition_component_summary": (
            transition.component_summary,
            ("csv", "parquet"),
        ),
        "transition_mixing_matrix": (transition_matrix, ("parquet",)),
        "transition_assortativity": (
            transition_assortativity,
            ("csv", "parquet"),
        ),
    }
    for name, (table, formats) in transition_jobs.items():
        LOGGER.info("Writing %s (%s rows)", name, f"{len(table):,}")
        write_table(table, name, formats=formats)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
