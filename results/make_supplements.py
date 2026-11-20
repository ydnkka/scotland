"""Build supplementary CSV exports from saved analysis results."""

from __future__ import annotations

import argparse
from pathlib import Path

from analyses.sse_detection.lib.concurrent_io import atomic_write_csv
from analyses.sse_detection.lib.figs.tab05 import build_fixed_effects_full_table
from analyses.sse_detection.lib.sse.config import BAYESIAN_OUTPUT_DIR

from .lib.config import SUPPLEMENTS_DIR


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bayesian-result-dir",
        type=Path,
        default=BAYESIAN_OUTPUT_DIR,
        help="Directory containing consolidated Bayesian model results.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SUPPLEMENTS_DIR,
        help="Directory for generated supplementary CSV files.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    table = build_fixed_effects_full_table(args.bayesian_result_dir)
    output_path = args.output_dir / "full_bayesian_estimates.csv"
    atomic_write_csv(table, output_path, index=False)
    print(f"Wrote {len(table):,} rows -> {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
