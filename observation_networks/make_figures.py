"""Regenerate Chapter 4 figures and LaTeX tables from saved result tables."""

from __future__ import annotations

import argparse
import logging

from .lib.figures import build_all_figures, build_all_tables
from .lib.io import ensure_results_dirs


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="Skip artifacts whose input tables have not been built yet.",
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
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    ensure_results_dirs()
    build_all_figures(skip_missing=args.skip_missing, logger=LOGGER)
    build_all_tables(skip_missing=args.skip_missing, logger=LOGGER)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
