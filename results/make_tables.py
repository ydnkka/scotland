"""Build publication table artifacts."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .lib.config import FIGURES_DIR, TABLES_DIR
from .lib.registry import DOMAINS, build_tables, list_builders, table_builders

LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "builders",
        nargs="*",
        help=(
            "Optional builder names. Use either 'name' if unique or "
            "'domain:name' for an exact builder."
        ),
    )
    parser.add_argument(
        "--domain",
        action="append",
        choices=DOMAINS,
        help="Restrict to one analysis domain. May be repeated.",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=FIGURES_DIR,
        help="Unused figure directory kept for a consistent builder context.",
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=TABLES_DIR,
        help="Directory for generated top-level LaTeX tables.",
    )
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="Skip builders whose input tables or model outputs are missing.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available table builders and exit.",
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

    if args.list:
        print(list_builders(table_builders(), domains=args.domain))
        return 0

    build_tables(
        names=args.builders,
        domains=args.domain,
        figure_dir=args.figure_dir,
        table_dir=args.table_dir,
        skip_missing=args.skip_missing,
        logger=LOGGER,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
