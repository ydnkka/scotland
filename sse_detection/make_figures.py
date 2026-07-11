"""Regenerate SSE detection figures from saved result tables."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .lib.figs.common import DEFAULT_TABLE_DIR, FIGURE_DIR
from .lib.figures import build_all_figures


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=DEFAULT_TABLE_DIR,
        help="Directory containing SSE output tables.",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=FIGURE_DIR,
        help="Directory for generated SSE figures.",
    )
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
    build_all_figures(
        table_dir=args.table_dir,
        figure_dir=args.figure_dir,
        skip_missing=args.skip_missing,
        logger=LOGGER,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
