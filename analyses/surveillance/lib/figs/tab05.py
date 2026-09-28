"""Build clade overtake events table."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..config import SEQUENCE_WINDOW_STRIDE, TABLES_DIR
from ..io import write_table
from .sequence_data import compute_lineage_frequency_tables, load_sequences

TABLE_NAME = "tab_clade_overtake_events"


def build(*, table_dir: Path = TABLES_DIR, window_stride: int = SEQUENCE_WINDOW_STRIDE):
    sequences = load_sequences(window_stride=window_stride)
    table = compute_lineage_frequency_tables(sequences, clade_col="variant")[4]
    return write_table(table, TABLE_NAME, table_dir=table_dir)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-dir", type=Path, default=TABLES_DIR)
    args = parser.parse_args()
    build(table_dir=args.table_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
