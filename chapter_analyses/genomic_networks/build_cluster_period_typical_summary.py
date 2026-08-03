"""Build a combined period-level non-singleton cluster summary table.

Run from the Scotland repository root after building ``cluster_period_summary``,
``window_coverage``, and ``cluster_pairwise_distance_summary``:

    python -m chapter_analyses.genomic_networks.build_cluster_period_typical_summary
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .lib.config import TABLES_DIR
from .lib.io import ensure_results_dirs, read_table, write_table
from .lib.period_summary import (
    MIN_PAIRWISE_ROWS,
    build_cluster_period_typical_summary,
)

LOGGER = logging.getLogger(__name__)
DEFAULT_TABLE_NAME = "cluster_period_typical_summary"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=TABLES_DIR,
        help=f"Input/output table directory. Default: {TABLES_DIR}.",
    )
    parser.add_argument(
        "--output-name",
        default=DEFAULT_TABLE_NAME,
        help=f"Output table stem under --table-dir. Default: {DEFAULT_TABLE_NAME}.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=("parquet", "csv"),
        help="Output formats passed to write_table. Default: parquet csv.",
    )
    parser.add_argument(
        "--min-pairwise-rows",
        type=int,
        default=MIN_PAIRWISE_ROWS,
        help=(
            "Minimum observed within-cluster pairwise rows per window-lineage "
            f"summary. Default: {MIN_PAIRWISE_ROWS}."
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

    LOGGER.info("Reading input tables from %s", args.table_dir)
    cluster_period = read_table("cluster_period_summary", table_dir=args.table_dir)
    pairwise_distance = read_table(
        "cluster_pairwise_distance_summary",
        table_dir=args.table_dir,
    )
    window_coverage = read_table("window_coverage", table_dir=args.table_dir)

    summary = build_cluster_period_typical_summary(
        cluster_period,
        pairwise_distance,
        window_coverage,
        min_pairwise_rows=args.min_pairwise_rows,
    )
    written = write_table(
        summary,
        args.output_name,
        table_dir=args.table_dir,
        formats=args.formats,
    )
    for fmt, path in written.items():
        LOGGER.info("Wrote %s output to %s", fmt, path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
