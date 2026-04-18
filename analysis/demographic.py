#!/usr/bin/env python3
"""
Characterise demographic factors (age, sex, vaccination) associated with clustering.

For each (window, resolution, cluster_id) computes:
  - n_sequences: cluster size
  - median_age: median age midpoint of members
  - frac_female: fraction female
  - frac_vaccinated: fraction with any prior vaccination
  - mean_vacc_dose: mean vacc_dose_number among members
  - age_diversity: std of age midpoints (within-cluster age spread)

Writes: data/processed/cluster_demographic_features.parquet
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


def compute_cluster_demographics(ds: pd.DataFrame) -> pd.DataFrame:
    grp = ds.groupby(["window_id", "resolution", "cluster_id"])
    return (
        grp.agg(
            n_sequences=("sequence_id", "nunique"),
            median_age=("age_midpoint", "median"),
            age_diversity=("age_midpoint", "std"),
            frac_female=("is_female", "mean"),
            frac_vaccinated=("is_vaccinated", "mean"),
            mean_vacc_dose=("vacc_dose_number", "mean"),
            wn_mid_date=("wn_mid_date", "first"),
            window_idx=("window_idx", "first"),
            pango_lineage=("pango_lineage", "first"),
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
    features = compute_cluster_demographics(ds)

    out = args.root / "data/processed/cluster_demographic_features.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    features.to_parquet(out, index=False, compression="zstd")
    print(f"Demographic cluster features: {len(features)} rows → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
