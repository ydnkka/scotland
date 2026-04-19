#!/usr/bin/env python3
"""
Characterise socioeconomic (SIMD) factors associated with clustering.

For each (window, resolution), computes per-cluster SIMD summaries then
aggregates to a cluster-level dataset suitable for regression.

Cluster-level features computed:
  - simd_quintile_mode: modal SIMD quintile of members' datazones
  - simd_decile_mode: modal SIMD decile
  - simd_rank_mean/median: continuous deprivation measure
  - frac_deprived_q1: fraction of members in SIMD quintile 1 (most deprived)
  - n_datazones: number of distinct datazones represented
  - n_sequences: cluster size
  - is_singleton: 1 if cluster size == 1

Writes: data/processed/cluster_simd_features.parquet
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def _mode(s: pd.Series):
    m = s.mode()
    return m.iloc[0] if len(m) > 0 else float("nan")


def compute_cluster_simd(ds: pd.DataFrame) -> pd.DataFrame:
    grp = ds.groupby(["window_id", "resolution", "cluster_id"])
    return (
        grp.agg(
            n_sequences=("sequence_id", "nunique"),
            n_datazones=("datazone", "nunique"),
            simd_rank_mean=("dz_simd_rank", "mean"),
            simd_rank_median=("dz_simd_rank", "median"),
            simd_quintile_mode=("dz_simd_quintile", _mode),
            simd_decile_mode=("dz_simd_decile", _mode),
            frac_deprived_q1=("dz_simd_quintile", lambda x: (x == 1).mean()),
            wn_mid_date=("wn_mid_date", "first"),
            window_idx=("window_idx", "first"),
            pango_lineage=("pango_lineage", lambda x: x.mode().iloc[0] if len(x) > 0 else None),
        )
        .reset_index()
        .assign(is_singleton=lambda df: (df["n_sequences"] == 1).astype(int))
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument("--root", type=Path, default=Path("."))
    args = ap.parse_args()

    cfg = load_config(args.root / args.config)
    proc = {k: args.root / v for k, v in cfg["data"]["processed"].items()}

    ds = pd.read_parquet(proc["analysis_dataset"])
    features = compute_cluster_simd(ds)

    out = args.root / "data/processed/cluster_simd_features.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    features.to_parquet(out, index=False, compression="zstd")
    print(f"SIMD cluster features: {len(features)} rows → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
