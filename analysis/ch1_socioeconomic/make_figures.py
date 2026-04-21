"""Orchestrate every figure for Chapter 1.

Run from the repository root:

    python -m analysis.ch1_socioeconomic.make_figures [--figures DIR] [--only N]
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from pathlib import Path

SCRIPTS = [
    "analysis.ch1_socioeconomic.scripts.fig1_sequences_by_simd_over_time",
    "analysis.ch1_socioeconomic.scripts.fig2_cluster_size_by_simd",
    "analysis.ch1_socioeconomic.scripts.fig3_simd_domain_forest",
    "analysis.ch1_socioeconomic.scripts.fig4_singleton_odds_by_epoch",
    "analysis.ch1_socioeconomic.scripts.fig5_deprivation_lineage_heatmap",
    "analysis.ch1_socioeconomic.scripts.fig6_domain_decomposition",
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures", type=Path, default=None,
                    help="Output directory (default: analysis/ch1_socioeconomic/figures)")
    ap.add_argument("--only", type=int, nargs="+", default=None,
                    help="Run only specific scripts by 1-based index")
    args = ap.parse_args()

    select = set(args.only) if args.only else set(range(1, len(SCRIPTS) + 1))
    for idx, module_name in enumerate(SCRIPTS, start=1):
        if idx not in select:
            continue
        t0 = time.time()
        print(f"[{idx}/{len(SCRIPTS)}] {module_name} …", flush=True)
        try:
            mod = importlib.import_module(module_name)
            path = mod.main(args.figures)
            print(f"    -> {path}   ({time.time() - t0:.1f}s)")
        except Exception as e:  # noqa: BLE001
            print(f"    !! FAILED: {e}", file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
