#!/usr/bin/env python3
"""
Prepare all processed metadata files from raw PHS and COG-UK data.

Outputs (parquet, relative to repo root unless --root is given):
  data/processed/scotland_sequence_metadata.parquet
  data/processed/scotland_testing.parquet
  data/processed/scotland_datazone_vaccinations.parquet
  data/processed/scotland_datazone_simd_data.parquet
  data/processed/scotland_geography.parquet

Usage:
  python3 method/01_prep_metadata.py
  python3 method/01_prep_metadata.py --config config.yaml --root /path/to/repo
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import yaml

# Fix PROJ_LIB for conda environments
_conda = os.environ.get("CONDA_PREFIX") or sys.prefix
os.environ.setdefault("PROJ_LIB", os.path.join(_conda, "share", "proj"))


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def band_to_midpoints(bands: pd.Series) -> pd.Series:
    s = bands.astype("string")
    lower = s.str.extract(r"(\d+)")[0].astype(float)
    upper = s.str.extract(r"-(\d+)")[0].astype(float)
    open_ended = s.str.endswith("+").fillna(False)
    upper = np.where(open_ended, lower + 5.0, upper)
    with np.errstate(invalid="ignore"):
        mid = (lower + upper) / 2.0
    return pd.to_numeric(pd.Series(mid, index=bands.index), errors="coerce")


def prep_testing(csv_path: Path, out_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path, parse_dates=["date_ecoss_specimen"], date_format="%Y%m%d")
    df.rename(columns={
        "date_ecoss_specimen": "collection_date",
        "PatientID": "patient_id",
        "datazone2011": "datazone",
    }, inplace=True)
    df.sort_values("collection_date", inplace=True)
    df.drop_duplicates(subset=["patient_id", "specimen_id", "datazone"], keep="first", inplace=True)

    total = df.groupby(["collection_date", "datazone"]).size().reset_index(name="dz_total_tests")
    pos = (df[df["test_result"] == "POSITIVE"]
           .groupby(["collection_date", "datazone"]).size()
           .reset_index(name="dz_positive_tests"))
    neg = (df[df["test_result"] == "NEGATIVE"]
           .groupby(["collection_date", "datazone"]).size()
           .reset_index(name="dz_negative_tests"))

    out = total.merge(pos, on=["collection_date", "datazone"], how="left")
    out = out.merge(neg, on=["collection_date", "datazone"], how="left")
    out.fillna(0, inplace=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_path, index=False, compression="zstd")
    logging.info("Testing data: %d rows → %s", len(out), out_path)
    return out


def prep_simd(csv_path: Path, out_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    col_map = {
        "DZ": "datazone",
        "Population": "dz_population",
        "Working_Age_Population": "dz_working_age_population",
        "SIMD2020v2_Rank": "dz_simd_rank",
        "SIMD2020v2_Quintile": "dz_simd_quintile",
        "SIMD2020v2_Decile": "dz_simd_decile",
        "SIMD2020v2_Vigintile": "dz_simd_vigintile",
        "SIMD2020v2_Income_Domain_Rank": "dz_simd_income_rank",
        "SIMD2020_Employment_Domain_Rank": "dz_simd_employment_rank",
        "SIMD2020_Education_Domain_Rank": "dz_simd_education_rank",
        "SIMD2020_Health_Domain_Rank": "dz_simd_health_rank",
        "SIMD2020_Access_Domain_Rank": "dz_simd_access_rank",
        "SIMD2020_Crime_Domain_Rank": "dz_simd_crime_rank",
        "SIMD2020_Housing_Domain_Rank": "dz_simd_housing_rank",
    }
    df.rename(columns=col_map, inplace=True)
    df = df[list(col_map.values())]
    assert df.notna().all().all(), "SIMD data contains unexpected nulls"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False, compression="zstd")
    logging.info("SIMD data: %d datazones → %s", len(df), out_path)
    return df


def prep_geography(shp_path: Path, simd: pd.DataFrame, out_path: Path) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(shp_path).set_index("DataZone")
    gdf.index.name = "datazone"
    gdf["dz_centroid"] = gdf.geometry.centroid
    gdf["dz_xcoord"] = gdf["dz_centroid"].x
    gdf["dz_ycoord"] = gdf["dz_centroid"].y

    simd_idx = simd.set_index("datazone")
    gdf = gdf[["dz_xcoord", "dz_ycoord"]].merge(
        simd_idx, how="left", left_index=True, right_index=True
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_parquet(out_path, index=True, compression="zstd")
    logging.info("Geography: %d datazones → %s", len(gdf), out_path)
    return gdf


def prep_vaccination(csv_path: Path, out_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path, low_memory=False)
    df["vacc_occurence_time"] = pd.to_datetime(df["vacc_occurence_time"], format="%Y%m%d", errors="coerce")
    df.dropna(subset=["vacc_occurence_time", "age_band"], inplace=True)
    df.rename(columns={
        "PatientID": "patient_id",
        "vacc_occurence_time": "vaccination_date",
        "datazone2011": "datazone",
    }, inplace=True)

    df["age_midpoint"] = band_to_midpoints(df["age_band"])

    out = (
        df.groupby(["vaccination_date", "datazone"])
        .agg(
            dz_total_vaccinated=("patient_id", "nunique"),
            dz_mean_vacc_age=("age_midpoint", "mean"),
            dz_median_vacc_age=("age_midpoint", "median"),
            dz_mean_vdose_number=("vacc_dose_number", "mean"),
            dz_median_vdose_number=("vacc_dose_number", "median"),
        )
        .reset_index()
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_path, index=False, compression="zstd")
    logging.info("Vaccination data: %d rows → %s", len(out), out_path)
    return out


def prep_sequence_metadata(
    metadata_csv: Path,
    nextclade_tsv: Path,
    vaccination_csv: Path,
    geography: gpd.GeoDataFrame,
    out_path: Path,
) -> pd.DataFrame:
    nextclade = pd.read_table(nextclade_tsv, low_memory=False)
    assert not nextclade["seqName"].duplicated(keep=False).any(), "Duplicate seqName in nextclade TSV"
    nextclade["seq_id"] = nextclade["seqName"].str.split("/").apply(lambda x: x[1])
    nextclade = nextclade.set_index("seq_id")

    meta = pd.read_csv(metadata_csv, parse_dates=["Collection_Date"])
    meta.rename(columns={
        "Collection_Date": "collection_date",
        "subject_sex": "sex",
        "SequenceID": "seq_id",
        "PatientID": "patient_id",
        "datazone2011": "datazone",
    }, inplace=True)

    assert set(meta["seq_id"]).issubset(set(nextclade.index)), "seq_ids not in nextclade"

    meta["sequence_id"] = nextclade.loc[meta["seq_id"].values, "seqName"].values
    meta["clade"] = nextclade.loc[meta["seq_id"].values, "clade"].values
    meta["who_voc"] = nextclade.loc[meta["seq_id"].values, "clade_who"].values
    meta["pango_lineage"] = nextclade.loc[meta["seq_id"].values, "Nextclade_pango"].values
    meta["nextclade_qc"] = nextclade.loc[meta["seq_id"].values, "qc.overallStatus"].values
    meta["age_midpoint"] = band_to_midpoints(meta["age_band"])

    meta.sort_values("collection_date", inplace=True)
    meta.drop_duplicates("specimen_id", keep="first", inplace=True)

    geo_cols = geography.reset_index()[["datazone", "dz_xcoord", "dz_ycoord"]]
    meta = meta.merge(geo_cols, on="datazone", how="left")

    # Attach latest vaccination before collection date
    vacc = pd.read_csv(vaccination_csv, low_memory=False)
    vacc["vacc_occurence_time"] = pd.to_datetime(vacc["vacc_occurence_time"], format="%Y%m%d", errors="coerce")
    vacc.dropna(subset=["vacc_occurence_time", "age_band"], inplace=True)
    vacc.rename(columns={"PatientID": "patient_id", "vacc_occurence_time": "vaccination_date"}, inplace=True)

    seq_vacc = meta.merge(vacc[["patient_id", "vaccination_date", "vacc_dose_number"]], on="patient_id", how="left")
    seq_vacc = seq_vacc[seq_vacc["vaccination_date"] <= seq_vacc["collection_date"]]
    latest_vacc = (
        seq_vacc.sort_values(["patient_id", "vaccination_date"])
        .groupby("patient_id")
        .tail(1)[["patient_id", "vaccination_date", "vacc_dose_number"]]
    )
    meta = meta.merge(latest_vacc, on="patient_id", how="left")
    meta["is_female"] = (meta["sex"] == "Female").astype(float)
    meta["is_vaccinated"] = (meta["vacc_dose_number"] > 0).astype(float)

    required = ["datazone", "collection_date", "patient_id", "sex", "age_band",
                "sequence_id", "clade", "pango_lineage", "nextclade_qc"]
    meta.dropna(subset=required, inplace=True)
    assert meta[required].notna().all().all()

    meta = meta[[
        "datazone", "dz_xcoord", "dz_ycoord",
        "collection_date", "patient_id", "sex", "is_female",
        "age_band", "age_midpoint", "sequence_id",
        "clade", "who_voc", "pango_lineage", "nextclade_qc",
        "vaccination_date", "vacc_dose_number", "is_vaccinated",
    ]].copy()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    meta.to_parquet(out_path, index=False, compression="zstd")
    logging.info("Sequence metadata: %d rows → %s", len(meta), out_path)
    return meta


def main() -> int:
    ap = argparse.ArgumentParser(description="Prepare all processed metadata from raw inputs.")
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument("--root", type=Path, default=Path("."),
                    help="Repo root (paths in config.yaml are relative to this)")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()

    setup_logging(args.log_level)
    cfg = load_config(args.root / args.config)
    raw = {k: args.root / v for k, v in cfg["data"]["raw"].items()}
    proc = {k: args.root / v for k, v in cfg["data"]["processed"].items()}

    _ = prep_testing(raw["testing_csv"], proc["testing"])
    simd = prep_simd(raw["simd_csv"], proc["simd"])
    geography = prep_geography(raw["geography_shp"], simd, proc["geography"])
    prep_vaccination(raw["vaccination_csv"], proc["vaccination"])
    prep_sequence_metadata(
        raw["metadata_csv"], raw["nextclade_tsv"], raw["vaccination_csv"],
        geography, proc["metadata"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
