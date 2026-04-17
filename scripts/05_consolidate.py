#!/usr/bin/env python3
"""
Consolidate per-group cluster parquets and merge with all metadata.

Steps:
1. Concatenate all per-group long parquets from cluster_long_dir.
2. Add singleton assignments (lineage groups of size 1, not processed by tn93).
3. Merge with sequence metadata, SIMD, testing, vaccination.
4. Compute window-level aggregates and date bounds.
5. Anonymise patient_id.
6. Write final analysis dataset parquet.

Output columns (46 total):
  window_idx, window_id, wn_start/mid/end_date, wn_no_sequences, wn_positive_tests,
  wn_prop_sequenced, sequence_id, patient_id, resolution, cluster_id,
  collection_date, datazone, dz_x, dz_y, sex, is_female, age_band, age_midpoint,
  is_vaccinated, vacc_dose_number, pango_lineage, clade, who_voc, nextclade_qc,
  dz_population, dz_working_age_population, dz_simd_rank/quintile/decile/vigintile,
  dz_simd_*_rank (×7), dz_total/positive/negative_tests, dz_total_vaccinated, dz_prop_vaccinated

Usage:
  python3 scripts/05_consolidate.py
  python3 scripts/05_consolidate.py --config config.yaml --root /path/to/repo
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from datetime import timedelta
from pathlib import Path

import pandas as pd
import yaml


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def build_windows(
    min_date: pd.Timestamp, max_date: pd.Timestamp,
    window_size: timedelta, step: timedelta,
) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    if max_date < min_date + window_size:
        return []
    starts = pd.date_range(start=min_date, end=max_date - window_size, freq=step)
    return [(s, s + window_size) for s in starts]


def load_cluster_parquets(cluster_long_dir: Path) -> pd.DataFrame:
    files = sorted(cluster_long_dir.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"No parquet files in {cluster_long_dir}")
    logging.info("Loading %d cluster parquets from %s", len(files), cluster_long_dir)
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


def build_singletons(
    metadata: pd.DataFrame,
    windows: list[tuple[pd.Timestamp, pd.Timestamp]],
    resolutions: list[float],
) -> pd.DataFrame:
    rows = []
    for i, (start, end) in enumerate(windows, start=1):
        window_id = f"W{str(i).zfill(3)}"
        wdf = metadata[(metadata["collection_date"] >= start) & (metadata["collection_date"] < end)]
        for lin, gdf in wdf.groupby("pango_lineage", sort=False):
            if gdf["sequence_id"].nunique() == 1:
                sid = gdf["sequence_id"].iloc[0]
                for res in resolutions:
                    rows.append({
                        "sequence_id": sid,
                        "window_id": window_id,
                        "resolution": res,
                        "cluster_id": f"{window_id}|{lin}|R{res}|C0",
                    })
    return pd.DataFrame(rows)


def build_window_info(
    metadata: pd.DataFrame,
    testing: pd.DataFrame,
    windows: list[tuple[pd.Timestamp, pd.Timestamp]],
) -> pd.DataFrame:
    rows = []
    for i, (start, end) in enumerate(windows, start=1):
        window_id = f"W{str(i).zfill(3)}"
        wdf = metadata[(metadata["collection_date"] >= start) & (metadata["collection_date"] < end)]
        tdf = testing[(testing["collection_date"] >= start) & (testing["collection_date"] < end)]
        pos_tests = tdf["dz_positive_tests"].sum()
        rows.append({
            "window_id": window_id,
            "wn_no_sequences": wdf["sequence_id"].nunique(),
            "wn_positive_tests": pos_tests,
            "wn_prop_sequenced": wdf["sequence_id"].nunique() / pos_tests if pos_tests > 0 else float("nan"),
        })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description="Consolidate cluster parquets into analysis dataset.")
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()

    setup_logging(args.log_level)
    cfg = load_config(args.root / args.config)
    pipe = cfg["pipeline"]
    proc = {k: args.root / v for k, v in cfg["data"]["processed"].items()}

    resolutions = [float(x) for x in pipe["leiden_resolutions"]]
    window_size = timedelta(weeks=pipe["window_size_weeks"])
    step = timedelta(weeks=pipe["step_weeks"])

    metadata = pd.read_parquet(proc["metadata"])
    metadata["collection_date"] = pd.to_datetime(metadata["collection_date"])
    testing = pd.read_parquet(proc["testing"])
    testing["collection_date"] = pd.to_datetime(testing["collection_date"])
    vaccination = pd.read_parquet(proc["vaccination"])
    vaccination["vaccination_date"] = pd.to_datetime(vaccination["vaccination_date"])
    simd = pd.read_parquet(proc["simd"])

    windows = build_windows(
        metadata["collection_date"].min(), metadata["collection_date"].max(), window_size, step
    )
    logging.info("%d time windows", len(windows))

    clustering = load_cluster_parquets(proc["cluster_long_dir"])
    singletons = build_singletons(metadata, windows, resolutions)
    logging.info(
        "Cluster rows: %d + %d singletons = %d total",
        len(clustering), len(singletons), len(clustering) + len(singletons),
    )
    clustering = pd.concat([clustering, singletons], ignore_index=True)

    window_info = build_window_info(metadata, testing, windows)

    full_meta = metadata.merge(simd, on="datazone", how="inner")
    full_meta = full_meta.drop(columns=["vaccination_date"], errors="ignore")

    ds = (
        clustering
        .merge(full_meta, on="sequence_id", how="inner")
        .merge(window_info, on="window_id", how="inner")
        .merge(testing, on=["collection_date", "datazone"], how="left")
        .merge(
            vaccination.rename(columns={"vaccination_date": "collection_date"}),
            on=["collection_date", "datazone"], how="left",
        )
    )

    ds["window_idx"] = ds["window_id"].str.extract(r"(\d+)", expand=False).astype(int)

    window_dates = (
        ds.groupby("window_idx")
        .agg(wn_start_date=("collection_date", "min"), wn_end_date=("collection_date", "max"))
        .reset_index()
    )
    window_dates["wn_mid_date"] = (
        window_dates["wn_start_date"]
        + (window_dates["wn_end_date"] - window_dates["wn_start_date"]) / 2
    )
    ds = ds.merge(window_dates, on="window_idx", how="left")

    ds["dz_prop_vaccinated"] = ds["dz_total_vaccinated"] / ds["dz_population"]

    # Anonymise patient_id
    uid = ds["patient_id"].unique()
    id_map = {old: f"P{i:06d}" for i, old in enumerate(uid, start=1)}
    ds["patient_id"] = ds["patient_id"].map(id_map)

    column_order = [
        "window_idx", "window_id", "wn_start_date", "wn_mid_date", "wn_end_date",
        "wn_no_sequences", "wn_positive_tests", "wn_prop_sequenced",
        "sequence_id", "patient_id", "resolution", "cluster_id",
        "collection_date", "datazone", "dz_x", "dz_y",
        "sex", "is_female", "age_band", "age_midpoint",
        "is_vaccinated", "vacc_dose_number",
        "pango_lineage", "clade", "who_voc", "nextclade_qc",
        "dz_population", "dz_working_age_population",
        "dz_simd_rank", "dz_simd_quintile", "dz_simd_decile", "dz_simd_vigintile",
        "dz_simd_income_rank", "dz_simd_employment_rank", "dz_simd_education_rank",
        "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank", "dz_simd_housing_rank",
        "dz_total_tests", "dz_positive_tests", "dz_negative_tests",
        "dz_total_vaccinated", "dz_prop_vaccinated",
    ]
    ds = ds[[c for c in column_order if c in ds.columns]].reset_index(drop=True)

    out_path: Path = proc["analysis_dataset"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_parquet(out_path, index=False, compression="zstd")
    logging.info("Analysis dataset: %s rows × %s cols → %s", len(ds), len(ds.columns), out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
