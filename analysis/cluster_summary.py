#!/usr/bin/env python3
"""
Summarise cluster characteristics per (window, resolution).

We only consider good-quality sequences (i.e. those with Nextclade QC = good).

For each (window_id, resolution) group computes:
  - n_clusters: number of distinct clusters
  - n_sequences: total sequences
  - cluster_size_mean/median/max: distribution of cluster sizes
  - n_singletons: clusters of size 1
  - singleton_frac: fraction of sequences in singleton clusters
  - largest_cluster_frac: fraction of sequences in the largest cluster

Writes: data/processed/cluster_summary.parquet
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


def summarise_clusters(ds: pd.DataFrame) -> pd.DataFrame:
    def _agg(g: pd.DataFrame) -> pd.Series:
        sizes = g.groupby("cluster_id")["sequence_id"].nunique()
        n_seq = sizes.sum()
        return pd.Series({
            "n_clusters": len(sizes),
            "n_sequences": int(n_seq),
            "cluster_size_mean": sizes.mean(),
            "cluster_size_median": sizes.median(),
            "cluster_size_max": int(sizes.max()),
            "n_singletons": int((sizes == 1).sum()),
            "singleton_frac": float((sizes == 1).sum() / len(sizes)),
            "largest_cluster_frac": float(sizes.max() / n_seq) if n_seq > 0 else float("nan"),
        })

    return (
        ds.groupby(["window_idx", "window_id", "wn_mid_date", "resolution"])
        .apply(_agg, include_groups=False)
        .reset_index()
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument("--root", type=Path, default=Path("."))
    args = ap.parse_args()

    cfg = load_config(args.root / args.config)
    proc = {k: args.root / v for k, v in cfg["data"]["processed"].items()}

    ds = pd.read_parquet(proc["analysis_dataset"])
    ds = ds[ds["nextclade_qc"] == "good"]  # restrict to good-quality sequences for this summary
    summary = summarise_clusters(ds)

    out = args.root / "data/processed/cluster_summary.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_parquet(out, index=False, compression="zstd")
    print(f"Cluster summary: {len(summary)} rows → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
