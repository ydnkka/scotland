"""Build policy index correlation table."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..config import TABLES_DIR
from ..io import write_table
from .policy_data import build_correlation_summary, load_policy_indices

TABLE_NAME = "tab_policy_index_correlation"


def build(*, table_dir: Path = TABLES_DIR, start_date=None, end_date=None):
    table = build_correlation_summary(
        load_policy_indices(start_date=start_date, end_date=end_date)
    )
    return write_table(table, TABLE_NAME, table_dir=table_dir)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-dir", type=Path, default=TABLES_DIR)
    args = parser.parse_args()
    build(table_dir=args.table_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
