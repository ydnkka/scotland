"""Build compact SIMD population-weighting validation tables for the appendix.

Run from the Scotland repository root:

    python -m chapter_analyses.genomic_networks.build_simd_validation
"""

from __future__ import annotations

import argparse
import logging

from .lib.io import ensure_results_dirs, write_table
from .lib.simd import build_simd_validation_tables


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--n-groups",
        type=int,
        default=5,
        choices=(5, 10, 20),
        help="SIMD grouping granularity. Default: 5 quintiles.",
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

    LOGGER.info("Building SIMD population-weighting validation tables")
    tables = build_simd_validation_tables(n_groups=args.n_groups)

    outputs = {
        "simd_population_weighting_datazone_assignments": (
            tables.datazone_assignments,
            ("parquet",),
        ),
        "simd_population_weighting_group_summary": (
            tables.group_summary,
            ("csv", "parquet"),
        ),
        "simd_population_weighting_movement": (
            tables.movement_table,
            ("csv", "parquet"),
        ),
        "simd_population_weighting_change_summary": (
            tables.change_summary,
            ("csv", "parquet"),
        ),
        "simd_population_weighting_diagnostics": (
            tables.diagnostics,
            ("csv", "parquet"),
        ),
    }

    for name, (table, formats) in outputs.items():
        LOGGER.info("Writing %s (%s rows)", name, f"{len(table):,}")
        write_table(table, name, formats=formats)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
