"""Orchestrate every figure for Paper 1.

Run from the repository root:

    python -m manuscripts.paper1_socioeconomic.make_figures [--output DIR] [--only N]
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from pathlib import Path

FIGURES = [
    "manuscripts.paper1_socioeconomic.figures.fig1_sequences_by_simd_over_time",
    "manuscripts.paper1_socioeconomic.figures.fig2_cluster_size_by_simd",
    "manuscripts.paper1_socioeconomic.figures.fig3_simd_domain_forest",
    "manuscripts.paper1_socioeconomic.figures.fig4_singleton_odds_by_epoch",
    "manuscripts.paper1_socioeconomic.figures.fig5_deprivation_lineage_heatmap",
    "manuscripts.paper1_socioeconomic.figures.fig6_domain_decomposition",
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=None,
                    help="Output directory (default: manuscripts/paper1_socioeconomic/output)")
    ap.add_argument("--only", type=int, nargs="+", default=None,
                    help="Run only specific figures by 1-based index")
    args = ap.parse_args()

    select = set(args.only) if args.only else set(range(1, len(FIGURES) + 1))
    for idx, module_name in enumerate(FIGURES, start=1):
        if idx not in select:
            continue
        t0 = time.time()
        print(f"[{idx}/{len(FIGURES)}] {module_name} …", flush=True)
        try:
            mod = importlib.import_module(module_name)
            path = mod.main(args.output)
            print(f"    -> {path}   ({time.time() - t0:.1f}s)")
        except Exception as e:  # noqa: BLE001
            print(f"    !! FAILED: {e}", file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
