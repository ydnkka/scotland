"""Orchestrate every figure for Paper 3."""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from pathlib import Path

FIGURES = [
    "manuscripts.paper3_spatial.scripts.fig1_cluster_centroid_map",
    "manuscripts.paper3_spatial.scripts.fig2_within_cluster_spread",
    "manuscripts.paper3_spatial.scripts.fig3_distance_decay",
    "manuscripts.paper3_spatial.scripts.fig4_urban_rural_footprint",
    "manuscripts.paper3_spatial.scripts.fig5_spread_by_deprivation",
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures", type=Path, default=None)
    ap.add_argument("--only", type=int, nargs="+", default=None)
    args = ap.parse_args()
    select = set(args.only) if args.only else set(range(1, len(FIGURES) + 1))
    for idx, name in enumerate(FIGURES, start=1):
        if idx not in select:
            continue
        t0 = time.time()
        print(f"[{idx}/{len(FIGURES)}] {name} …", flush=True)
        try:
            mod = importlib.import_module(name)
            p = mod.main(args.output)
            print(f"    -> {p}   ({time.time() - t0:.1f}s)")
        except Exception as e:  # noqa: BLE001
            print(f"    !! FAILED: {e}", file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
