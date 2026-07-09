"""Regenerate Chapter 4 figures from saved observation_networks tables."""

from __future__ import annotations

import argparse
import logging

from .lib.figures import (
    plot_assortativity_over_time,
    plot_clade_frequencies,
    plot_cluster_size_summary,
    plot_degree_assortativity_over_time,
    plot_transition_window_summary,
    plot_window_coverage,
)
from .lib.io import ensure_results_dirs, read_table


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="Skip figures whose input tables have not been built yet.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    )
    return parser.parse_args()


def _run(name: str, func, table_name: str, *, skip_missing: bool) -> None:
    try:
        table = read_table(table_name)
    except FileNotFoundError:
        if skip_missing:
            LOGGER.warning("Skipping %s: missing %s", name, table_name)
            return
        raise
    LOGGER.info("Writing figure %s", name)
    func(table)


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=args.log_level, format="%(levelname)s: %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    ensure_results_dirs()

    jobs = [
        ("window_coverage", plot_window_coverage, "window_coverage"),
        ("clade_window_frequencies", plot_clade_frequencies, "clade_window_counts"),
        ("cluster_size_summary", plot_cluster_size_summary, "cluster_window_summary"),
        (
            "transition_graph_window_summary",
            plot_transition_window_summary,
            "transition_window_summary",
        ),
        (
            "compatibility_assortativity",
            plot_assortativity_over_time,
            "compatibility_assortativity",
        ),
        (
            "compatibility_degree_assortativity",
            plot_degree_assortativity_over_time,
            "compatibility_degree_assortativity",
        ),
    ]
    for name, func, table_name in jobs:
        _run(name, func, table_name, skip_missing=args.skip_missing)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
