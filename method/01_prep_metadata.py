"""
Prepare all processed metadata files from raw PHS and COG-UK data.

Outputs (parquet, relative to repo root unless --root is given):
    data/processed/scotland_sequence_metadata.parquet
    data/processed/scotland_testing.parquet
    data/processed/scotland_datazone_vaccinations.parquet
    data/processed/scotland_datazone_simd_data.parquet
    data/processed/scotland_geography.parquet
    data/processed/scotland_hb_daily_trends.parquet
    data/processed/scotland_policy.parquet

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


AGE_GROUP_MAP = {
    "00-04": "00-04",
    "05-09": "05-14",
    "10-14": "05-14",
    "15-19": "15-24",
    "20-24": "15-24",
    "25-29": "25-64",
    "30-34": "25-64",
    "35-39": "25-64",
    "40-44": "25-64",
    "45-49": "25-64",
    "50-54": "25-64",
    "55-59": "25-64",
    "60-64": "25-64",
    "65-69": "65-74",
    "70-74": "65-74",
    "75+": "75+",
}

TEST_REASON_MAP = {
    "symptomatic-citizen": "symptomatic_citizen",
    "I have coronavirus symptoms": "symptomatic_citizen",
    "I live~ work or study in a lockdown area with a coronavirus outbreak": "symptomatic_citizen",
    "symptomatic-essential-worker": "symptomatic_essential_worker",
    "Im an essential worker": "symptomatic_essential_worker",
    "scotland-wales-keyworker": "symptomatic_essential_worker",
    "wales-keyworker": "symptomatic_essential_worker",
    "test-for-contact-tracing": "contact_tracing",
    "test-for-contact-tracing-app": "contact_tracing",
    "test-for-contact-self-referral": "contact_tracing",
    "for-symptomatic-household-member": "contact_tracing",
    "Ive been in contact with a person who has tested positive for coronavirus and Ive been asked to take a test by a contact tracer (Northern Ireland and Scotland)": "contact_tracing",
    "Ive been in contact with a person who has tested positive for coronavirus and have since developed symptoms": "contact_tracing",
    "confirmatory-positive-test": "confirmatory",
    "confirmatory-other-reason": "confirmatory",
    "confirmatory-test-unclear": "confirmatory",
    "confirmatory-test-borders": "confirmatory",
    "told-to-order-repeat-test": "confirmatory",
    "self-isolation-support-grant": "isolation_scheme",
    "isolation-testing-home": "isolation_scheme",
    "isolation-testing-facility": "isolation_scheme",
    "gp-healthcare-request": "clinical",
    "antiviral-order": "clinical",
    "dental-patient-testing": "clinical",
    "I have been told to have a test before I go into hospital~ for example~ for surgery": "clinical",
    "zoe-symptom-study": "surveillance_research",
    "contact-testing-study": "surveillance_research",
    "events-research-programme": "surveillance_research",
    "serial-testing": "surveillance_research",
    "ntrg-member": "surveillance_research",
    "local-council-request": "local_outbreak",
    "attended-outbreak-venue": "local_outbreak",
    "community-testing": "local_outbreak",
    "scotland-university": "local_outbreak",
    "wales-university": "local_outbreak",
    "green-traveller": "travel",
    "other": "other",
    "Other": "other",
    "none": "other",
    "do-not-know": "other",
    "general-cta-referral": "other",
    "personal-assistant": "other",
    "Im a visiting professional": "other",
    "asymptomatic-home-order": "other",
    "school-trained-staff-test": "other",
    "school-self-supervised-test": "other",
}


def setup_logging(level: str = "INFO") -> None:
    """Configure timestamped console logging for the metadata prep script."""
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def load_config(path: Path) -> dict:
    """Load the YAML pipeline configuration file."""
    with open(path) as f:
        return yaml.safe_load(f)


def get_age_midpoints(bands: pd.Series) -> np.ndarray:
    """Convert age-band strings (e.g. '30-34', '85+') to numeric midpoints.

    Open-ended bands (e.g. '85+') are given a half-width of 5 years, placing
    their midpoint at lower + 2.5 — consistent with closed 5-year bands.
    """
    s = bands.astype("string")
    lower = s.str.extract(r"(\d+)")[0].astype(float)
    upper = s.str.extract(r"-(\d+)")[0].astype(float)
    open_ended = s.str.endswith("+").fillna(False)
    # Give open-ended bands (e.g. '85+') a synthetic upper bound of lower + 5
    upper = np.where(open_ended, lower + 5.0, upper)
    with np.errstate(invalid="ignore"):
        mid = (lower + upper) / 2.0
    return pd.to_numeric(pd.Series(mid, index=bands.index), errors="coerce").to_numpy()


def _load_oxcgrt_index(
    csv_path: Path,
    *,
    region_name: str,
    value_col: str,
) -> pd.DataFrame:
    """Select one region and reshape a wide OxCGRT index to daily rows."""
    table = pd.read_csv(csv_path)
    if "RegionName" not in table.columns:
        raise KeyError(f"OxCGRT table needs 'RegionName': {csv_path}")

    selected = table.loc[table["RegionName"].eq(region_name)]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one OxCGRT row for {region_name!r} in {csv_path}; "
            f"found {len(selected)}."
        )

    parsed_dates = pd.to_datetime(table.columns, format="%d%b%Y", errors="coerce")
    date_mask = parsed_dates.notna()
    date_columns = table.columns[date_mask]
    return pd.DataFrame(
        {
            "date": parsed_dates[date_mask],
            value_col: pd.to_numeric(
                selected.iloc[0][date_columns], errors="coerce"
            ).to_numpy(dtype=float),
        }
    ).sort_values("date", ignore_index=True)


def prep_policy(
    stringency_csv: Path,
    containment_csv: Path,
    period_spec: list[dict],
    out_path: Path,
    *,
    region_name: str = "Scotland",
) -> pd.DataFrame:
    """Build the daily Scotland policy lookup from OxCGRT and period metadata."""
    periods = pd.DataFrame(period_spec).rename(
        columns={
            "code": "period_code",
            "label": "period_label",
            "start_date": "period_start_date",
            "end_date": "period_end_date",
            "era": "policy_era",
        }
    )
    required = {
        "period_code",
        "period_label",
        "period_start_date",
        "period_end_date",
        "policy_era",
    }
    missing = sorted(required - set(periods.columns))
    if missing:
        raise KeyError(f"Policy period specification is missing columns: {missing}")
    if periods.empty:
        raise ValueError("At least one policy period must be configured.")

    periods["period_start_date"] = pd.to_datetime(
        periods["period_start_date"], errors="raise"
    ).dt.normalize()
    periods["period_end_date"] = pd.to_datetime(
        periods["period_end_date"], errors="raise"
    ).dt.normalize()
    periods = periods.sort_values("period_start_date", ignore_index=True)
    periods["period_order"] = np.arange(len(periods), dtype=int)

    if periods["period_code"].duplicated().any():
        duplicates = periods.loc[
            periods["period_code"].duplicated(keep=False), "period_code"
        ].tolist()
        raise ValueError(f"Duplicate policy period code(s): {duplicates}")
    if periods["period_start_date"].duplicated().any():
        raise ValueError("Policy periods contain duplicate start dates.")
    if periods["period_end_date"].lt(periods["period_start_date"]).any():
        raise ValueError("A policy period ends before it starts.")

    expected_starts = periods["period_end_date"].shift(1) + pd.Timedelta(days=1)
    discontinuous = periods.index[1:][
        periods.loc[1:, "period_start_date"].ne(expected_starts.loc[1:])
    ]
    if len(discontinuous):
        raise ValueError(
            "Policy periods must be contiguous and non-overlapping; check rows "
            f"{discontinuous.tolist()}."
        )

    stringency = _load_oxcgrt_index(
        stringency_csv,
        region_name=region_name,
        value_col="stringency_index",
    )
    containment = _load_oxcgrt_index(
        containment_csv,
        region_name=region_name,
        value_col="containment_index",
    )
    indices = stringency.merge(
        containment,
        on="date",
        how="outer",
        validate="one_to_one",
    )

    daily = pd.DataFrame(
        {
            "date": pd.date_range(
                periods["period_start_date"].min(),
                periods["period_end_date"].max(),
                freq="D",
            )
        }
    ).merge(indices, on="date", how="left", validate="one_to_one")

    for row in periods.itertuples(index=False):
        mask = daily["date"].between(
            row.period_start_date,
            row.period_end_date,
            inclusive="both",
        )
        daily.loc[mask, "period_code"] = str(row.period_code)
        daily.loc[mask, "period_label"] = str(row.period_label)
        daily.loc[mask, "period_start_date"] = row.period_start_date
        daily.loc[mask, "period_end_date"] = row.period_end_date
        daily.loc[mask, "period_order"] = int(row.period_order) # type: ignore
        daily.loc[mask, "policy_era"] = str(row.policy_era)

    if daily["period_code"].isna().any():
        raise AssertionError("The daily policy calendar contains unassigned dates.")
    daily["period_order"] = daily["period_order"].astype(int)
    daily = daily[
        [
            "date",
            "stringency_index",
            "containment_index",
            "period_code",
            "period_label",
            "period_start_date",
            "period_end_date",
            "period_order",
            "policy_era",
        ]
    ]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    daily.to_parquet(out_path, index=False, compression="zstd")
    logging.info(
        "Policy data: %d daily Scotland rows (%s to %s) \u2192 %s",
        len(daily),
        daily["date"].min().date(),
        daily["date"].max().date(),
        out_path,
    )
    missing_indices = daily[["stringency_index", "containment_index"]].isna().all(axis=1)
    if missing_indices.any():
        logging.warning(
            "Policy data has %d calendar dates after/between available OxCGRT values.",
            int(missing_indices.sum()),
        )
    return daily


def prep_testing(csv_path: Path, out_path: Path) -> pd.DataFrame:
    """Aggregate raw test records to daily datazone-level counts by result and test type.

    Output columns per (collection_date, datazone): 
        dz_total_tests          — all tests
        dz_positive_tests       — POSITIVE results (PCR + LFD)
        dz_negative_tests       — NEGATIVE results
        dz_pcr_positive_tests   — PCR-confirmed positives
        dz_lfd_positive_tests   — LFD-only positives (test_type == ANTIGEN)
        dz_care_home_tests      — tests linked to a care home facility (care_home_id not null)
    """
    df = pd.read_csv(
        csv_path, parse_dates=["date_ecoss_specimen"], date_format="%Y%m%d"
    )
    df.rename(
        columns={
            "date_ecoss_specimen": "collection_date",
            "PatientID": "patient_id",
            "datazone2011": "datazone",
        },
        inplace=True,
    )
    df.sort_values("collection_date", inplace=True)

    # Drop duplicate records for the same patient/specimen/datazone on the same day,
    # keeping the first occurrence after sorting by date.
    df.drop_duplicates(
        subset=["patient_id", "specimen_id", "datazone"], keep="first", inplace=True
    )

    # Aggregate separately then merge so that zeros are preserved for
    # datazones that had tests but no positives/negatives on a given day.
    grp = ["collection_date", "datazone"]
    total = df.groupby(grp).size().reset_index(name="dz_total_tests")
    pos = (
        df[df["test_result"] == "POSITIVE"]
        .groupby(grp)
        .size()
        .reset_index(name="dz_positive_tests")
    )
    neg = (
        df[df["test_result"] == "NEGATIVE"]
        .groupby(grp)
        .size()
        .reset_index(name="dz_negative_tests")
    )
    pcr_pos = (
        df[(df["test_result"] == "POSITIVE") & (df["test_type"] == "PCR")]
        .groupby(grp)
        .size()
        .reset_index(name="dz_pcr_positive_tests")
    )
    lfd_pos = (
        df[(df["test_result"] == "POSITIVE") & (df["test_type"] == "ANTIGEN")]
        .groupby(grp)
        .size()
        .reset_index(name="dz_lfd_positive_tests")
    )
    care_hm = (
        df[df["care_home_id"].notna()]
        .groupby(grp)
        .size()
        .reset_index(name="dz_care_home_tests")
    )

    out = total
    for part in (pos, neg, pcr_pos, lfd_pos, care_hm):
        out = out.merge(part, on=grp, how="left")
    out.fillna(0, inplace=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_path, index=False, compression="zstd")
    logging.info("Testing data: %d rows \u2192 %s", len(out), out_path)
    return out


def prep_simd(csv_path: Path, out_path: Path) -> pd.DataFrame:
    """Select and rename SIMD 2020v2 domain ranks/deciles for each datazone.

    Includes Health Board and Local Authority codes (needed to join daily HB trends).
    """
    df = pd.read_csv(csv_path)
    col_map = {
        "DZ": "datazone",
        "Population": "dz_population",
        "Working_Age_Population": "dz_working_age_population",
        "URname": "dz_urban_rural_class",
        "LAcode": "dz_local_authority_code",
        "LAname": "dz_local_authority",
        "HBcode": "dz_health_board_code",
        "HBname": "dz_health_board",
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

    # Keep only columns that are present (LAcode/HBcode may be absent in some SIMD releases)
    available = {k: v for k, v in col_map.items() if k in df.columns}
    missing = set(col_map) - set(available)
    if missing:
        logging.warning("SIMD CSV is missing expected columns: %s", sorted(missing))

    df.rename(columns=available, inplace=True)
    df = df[list(available.values())]

    null_counts = df.isna().sum()
    null_cols = null_counts[null_counts > 0]
    if not null_cols.empty:
        n_before = len(df)
        df = df.dropna()
        logging.warning(
            "SIMD data has null values in %d column(s); dropped %d/%d datazones. "
            "Affected columns: %s",
            len(null_cols),
            n_before - len(df),
            n_before,
            null_cols.index.tolist(),
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    logging.info("SIMD data: %d datazones \u2192 %s", len(df), out_path)
    return df


def prep_geography(
    shp_path: Path, simd: pd.DataFrame, out_path: Path
) -> gpd.GeoDataFrame:
    """Compute datazone centroids and join SIMD attributes; retain geometry for geoparquet output."""
    gdf = gpd.read_file(shp_path).set_index("DataZone")
    gdf.index.name = "datazone"

    # Compute centroids in the native CRS (OSGB36 / EPSG:27700)
    gdf["dz_centroid"] = gdf.geometry.centroid
    gdf["dz_xcoord"] = gdf["dz_centroid"].x
    gdf["dz_ycoord"] = gdf["dz_centroid"].y

    simd_idx = simd.set_index("datazone")

    # Area in km² from the shapefile; used downstream to compute population density.
    gdf["dz_area_km2"] = gdf["StdAreaKm2"]
    gdf = gdf[["geometry", "dz_xcoord", "dz_ycoord", "dz_area_km2"]].merge(
        simd_idx, how="left", left_index=True, right_index=True
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_parquet(out_path, index=True)
    logging.info("Geography: %d datazones \u2192 %s", len(gdf), out_path)
    return gdf


def prep_vaccination(csv_path: Path, out_path: Path) -> pd.DataFrame:
    """Aggregate raw vaccination records to daily datazone-level summary statistics."""
    df = pd.read_csv(csv_path, low_memory=False)
    df["vacc_occurence_time"] = pd.to_datetime(
        df["vacc_occurence_time"], format="%Y%m%d", errors="coerce"
    )
    df.dropna(subset=["vacc_occurence_time", "age_band"], inplace=True)
    df.rename(
        columns={
            "PatientID": "patient_id",
            "vacc_occurence_time": "vaccination_date",
            "datazone2011": "datazone",
        },
        inplace=True,
    )

    df["age_midpoint"] = get_age_midpoints(df["age_band"])

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
    logging.info("Vaccination data: %d rows \u2192 %s", len(out), out_path)
    return out


def prep_hb_trends(csv_path: Path, out_path: Path) -> pd.DataFrame:
    """Clean and select daily health board COVID trends.

    Keeps one row per (Health Board, date) with hospital, ICU, reinfection,
    and case-count metrics.  The Scotland-wide aggregate row (HB = S92000003)
    is dropped; per-HB rows are retained.

    Output columns:
        date, hb_code,
        hb_daily_positive, hb_cumulative_positive,
        hb_hospital_admissions, hb_hospital_occupancy,
        hb_icu_admissions, hb_icu_occupancy_lt28d, hb_icu_occupancy_ge28d,
        hb_daily_reinfections, hb_reinfection_rate,
        hb_total_tests, hb_positive_tests
    """
    df = pd.read_csv(csv_path, low_memory=False)
    df["Date"] = pd.to_datetime(
        df["Date"].astype(str), format="%Y%m%d", errors="coerce"
    )
    df.dropna(subset=["Date", "HB"], inplace=True)

    # Drop Scotland-wide aggregate
    df = df[df["HB"] != "S92000003"].copy()

    rename = {
        "Date": "date",
        "HB": "hb_code",
        "DailyPositive": "hb_daily_positive",
        "CumulativePositive": "hb_cumulative_positive",
        "HospitalAdmissions": "hb_hospital_admissions",
        "HospitalOccupancy": "hb_hospital_occupancy",
        "ICUAdmissions": "hb_icu_admissions",
        "ICUOccupancy28daysless": "hb_icu_occupancy_lt28d",
        "ICUOccupancy28daysmore": "hb_icu_occupancy_ge28d",
        "Reinfections": "hb_daily_reinfections",
        "PercentReinfections": "hb_reinfection_rate",
        "TotalTests": "hb_total_tests",
        "PositiveTests": "hb_positive_tests",
    }
    df.rename(
        columns={k: v for k, v in rename.items() if k in df.columns}, inplace=True
    )

    keep = [v for v in rename.values() if v in df.columns]
    df = df[keep].copy()

    # Coerce all metric columns to numeric
    for col in keep:
        if col not in ("date", "hb_code"):
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df.sort_values(["hb_code", "date"], inplace=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False, compression="zstd")
    logging.info("HB trends: %d rows \u2192 %s", len(df), out_path)
    return df


def prep_sequence_metadata(
    metadata_csv: Path,
    nextclade_tsv: Path,
    vaccination_csv: Path,
    testing_csv: Path,
    geography: gpd.GeoDataFrame,
    policy: pd.DataFrame,
    out_path: Path,
    required_nextclade_qc: str,
) -> pd.DataFrame:
    """Join sequence metadata with Nextclade QC/lineage calls, geography, vaccination
    history, and test-level attributes (test type, reason, S-gene status).

    The Nextclade QC filter is applied immediately after annotations are attached, so
    the written metadata and every downstream pipeline stage share one QC cohort.
    """

    # --- Nextclade annotations ---
    nextclade = pd.read_table(nextclade_tsv, low_memory=False)
    assert not nextclade["seqName"].duplicated(keep=False).any(), (
        "Duplicate seqName in nextclade TSV"
    )
    # Extract the numeric PHS sequence ID from the COGUK seqName (e.g. 'Scotland/LIVE-XXXXX/2021')
    nextclade["seq_id"] = nextclade["seqName"].str.split("/").apply(lambda x: x[1])
    nextclade = nextclade.set_index("seq_id")

    # --- Core sequence metadata ---
    meta = pd.read_csv(metadata_csv, parse_dates=["Collection_Date"])
    meta.rename(
        columns={
            "Collection_Date": "collection_date",
            "subject_sex": "sex",
            "SequenceID": "seq_id",
            "PatientID": "patient_id",
            "datazone2011": "datazone",
        },
        inplace=True,
    )
    meta["specimen_id"] = meta["specimen_id"].astype("string").str.strip()

    assert set(meta["seq_id"]).issubset(set(nextclade.index)), (
        "seq_ids not in nextclade"
    )

    # Attach Nextclade fields using seq_id as the alignment key
    nextclade_aligned = nextclade.reindex(meta["seq_id"])
    meta["sequence_id"] = nextclade_aligned["seqName"].to_numpy()
    meta["clade"] = nextclade_aligned["clade"].to_numpy()
    meta["who_voc"] = nextclade_aligned["clade_who"].to_numpy()
    meta["pango_lineage"] = nextclade_aligned["Nextclade_pango"].to_numpy()
    meta["nextclade_qc"] = nextclade_aligned["qc.overallStatus"].to_numpy()

    qc_status = str(required_nextclade_qc).strip()
    if not qc_status:
        raise ValueError("tn93.nextclade_qc must be a non-empty status.")
    qc_matches = (
        meta["nextclade_qc"]
        .astype("string")
        .str.strip()
        .str.casefold()
        .eq(qc_status.casefold())
    )
    n_before_qc = len(meta)
    meta = meta.loc[qc_matches].copy()
    logging.info(
        "Nextclade QC filter (%s): retained %d/%d sequences; dropped %d.",
        qc_status,
        len(meta),
        n_before_qc,
        n_before_qc - len(meta),
    )
    if meta.empty:
        raise ValueError(
            f"No sequence metadata rows have Nextclade QC status {qc_status!r}."
        )

    meta["age_midpoint"] = get_age_midpoints(meta["age_band"])
    meta["age_group"] = meta["age_band"].map(AGE_GROUP_MAP)
    unmapped_age_bands = sorted(
        meta.loc[meta["age_band"].notna() & meta["age_group"].isna(), "age_band"]
        .astype(str)
        .unique()
    )
    if unmapped_age_bands:
        raise ValueError(f"Unmapped sequence age band(s): {unmapped_age_bands}")

    meta.sort_values("collection_date", inplace=True)
    # A small number of specimens have duplicate records; keep the earliest entry
    meta.drop_duplicates("specimen_id", keep="first", inplace=True)

    # --- Geography: attach datazone centroid coordinates and area ---
    _geo = geography.reset_index()
    _geo_want = ["datazone", "dz_xcoord", "dz_ycoord", "dz_area_km2"]
    geo_cols = _geo[[c for c in _geo_want if c in _geo.columns]]
    meta = meta.merge(geo_cols, on="datazone", how="left")

    # --- Test-level attributes and reinfection flag ---
    # Load test records with both the per-specimen fields (type/reason/S-gene) and the
    # per-patient/date fields needed to derive an individual-level reinfection indicator.
    tests_raw = pd.read_csv(
        testing_csv,
        usecols=[
            "specimen_id",
            "PatientID",
            "date_ecoss_specimen",
            "test_result",
            "test_type",
            "test_reason",
            "test_result_s_gene_status",
        ],
        low_memory=False,
    )
    tests_raw["date_ecoss_specimen"] = pd.to_datetime(
        tests_raw["date_ecoss_specimen"], format="%Y%m%d", errors="coerce"
    )
    tests_raw.rename(
        columns={
            "PatientID": "patient_id",
            "date_ecoss_specimen": "collection_date",
            "test_result_s_gene_status": "s_gene_status",
        },
        inplace=True,
    )
    tests_raw["specimen_id"] = tests_raw["specimen_id"].astype("string").str.strip()

    # Specimen-level attributes: one row per specimen (keep earliest if duplicates).
    test_attr_cols = ["test_type", "test_reason", "s_gene_status"]
    specimen_attrs = tests_raw[["specimen_id"] + test_attr_cols].drop_duplicates(
        "specimen_id", keep="first"
    )
    # The raw sequence metadata can already carry test_type/test_reason, which would
    # otherwise force pandas to suffix the testing-derived columns and silently drop
    # them from the final output schema.
    meta = meta.merge(
        specimen_attrs, on="specimen_id", how="left", suffixes=("", "_testing")
    )
    for col in test_attr_cols:
        testing_col = f"{col}_testing"
        if testing_col in meta.columns:
            if col in meta.columns:
                meta[col] = meta[testing_col].combine_first(meta[col])
            else:
                meta[col] = meta[testing_col]
            meta.drop(columns=[testing_col], inplace=True)

    # Preserve the source value and expose one stable analytical grouping. Missing
    # reasons remain distinguishable from observed reasons reported as "other".
    meta.rename(columns={"test_reason": "test_reason_raw"}, inplace=True)
    meta["test_reason"] = meta["test_reason_raw"].map(TEST_REASON_MAP)
    missing_reason = meta["test_reason_raw"].isna()
    unmapped_reason = meta["test_reason_raw"].notna() & meta["test_reason"].isna()
    if unmapped_reason.any():
        values = sorted(meta.loc[unmapped_reason, "test_reason_raw"].astype(str).unique())
        logging.warning(
            "Mapped %d records with unknown test reasons to 'other': %s",
            int(unmapped_reason.sum()),
            values,
        )
    meta.loc[unmapped_reason, "test_reason"] = "other"
    meta.loc[missing_reason, "test_reason"] = "missing"

    # Attach the sample-level policy classification from the processed daily lookup.
    policy_lookup = policy[
        ["date", "period_code", "period_label", "policy_era"]
    ].rename(
        columns={
            "date": "collection_date",
            "period_code": "policy_period",
            "period_label": "policy_period_label",
        }
    )
    meta = meta.merge(
        policy_lookup,
        on="collection_date",
        how="left",
        validate="many_to_one",
    )
    missing_policy = meta["policy_period"].isna()
    if missing_policy.any():
        missing_dates = meta.loc[missing_policy, "collection_date"]
        raise ValueError(
            f"No policy period for {int(missing_policy.sum())} sequence record(s), "
            f"covering {missing_dates.min()} to {missing_dates.max()}."
        )

    # Reinfection flag: is_reinfection = 1 if this positive test occurred ≥ 90 days
    # after the same patient's most-recent prior positive test.  First positives and
    # negative/void tests receive 0.  The 90-day threshold follows standard PHS/ECDC
    # definitions for SARS-CoV-2 reinfection.
    positives = (
        tests_raw[tests_raw["test_result"] == "POSITIVE"][
            ["patient_id", "collection_date"]
        ]
        .dropna()
        .drop_duplicates()
        .sort_values(["patient_id", "collection_date"])
    )
    positives["_prev_pos"] = positives.groupby("patient_id")["collection_date"].shift(1)
    positives["is_reinfection"] = (
        (positives["collection_date"] - positives["_prev_pos"]).dt.days >= 90
    ).astype(float)
    positives = positives[["patient_id", "collection_date", "is_reinfection"]]
    meta = meta.merge(positives, on=["patient_id", "collection_date"], how="left")
    meta["is_reinfection"] = meta["is_reinfection"].fillna(0)

    # --- Vaccination history: find the most recent dose before each collection date ---
    vacc = pd.read_csv(vaccination_csv, low_memory=False)
    vacc["vacc_occurence_time"] = pd.to_datetime(
        vacc["vacc_occurence_time"], format="%Y%m%d", errors="coerce"
    )
    vacc.dropna(subset=["vacc_occurence_time", "age_band"], inplace=True)
    vacc.rename(
        columns={"PatientID": "patient_id", "vacc_occurence_time": "vaccination_date"},
        inplace=True,
    )
    vacc["vacc_booster"] = np.where(vacc["vacc_booster"] == "TRUE", 1, 0)

    # Cross-join sequences with vaccination records for the same patient, then
    # keep only vaccinations that occurred on or before the sequence collection date.
    # This produces one row per (sequence, prior vaccination) combination.
    vacc_cols = [
        "patient_id",
        "vaccination_date",
        "vacc_dose_number",
        "vacc_product_name",
        "vacc_booster",
    ]
    seq_vacc = meta.merge(vacc[vacc_cols], on="patient_id", how="left")
    seq_vacc = seq_vacc[seq_vacc["vaccination_date"] <= seq_vacc["collection_date"]]

    # Group by both patient_id AND collection_date so that patients with
    # multiple sequenced specimens each get the correct latest-prior-dose for
    # their own collection date.
    latest_vacc_cols = [
        "patient_id",
        "collection_date",
        "vaccination_date",
        "vacc_dose_number",
    ]
    for extra in ("vacc_product_name", "vacc_booster"):
        if extra in seq_vacc.columns:
            latest_vacc_cols.append(extra)

    latest_vacc = (
        seq_vacc.sort_values(["patient_id", "collection_date", "vaccination_date"])
        .groupby(["patient_id", "collection_date"])
        .tail(1)[latest_vacc_cols]
    )
    meta = meta.merge(latest_vacc, on=["patient_id", "collection_date"], how="left")

    meta["is_female"] = (meta["sex"] == "Female").astype(float)
    meta["vacc_dose_number"] = meta["vacc_dose_number"].fillna(0)
    meta["is_vaccinated"] = (meta["vacc_dose_number"] > 0).astype(float)

    # Drop rows missing any field required for downstream modelling
    required = [
        "datazone",
        "collection_date",
        "patient_id",
        "sex",
        "age_band",
        "sequence_id",
        "clade",
        "pango_lineage",
        "nextclade_qc",
    ]
    meta.dropna(subset=required, inplace=True)
    assert meta[required].notna().all().all()

    output_cols = [
        "datazone",
        "dz_xcoord",
        "dz_ycoord",
        "dz_area_km2",
        "collection_date",
        "patient_id",
        "sex",
        "is_female",
        "age_band",
        "age_group",
        "age_midpoint",
        "specimen_id",
        "sequence_id",
        "clade",
        "who_voc",
        "pango_lineage",
        "nextclade_qc",
        "test_type",
        "test_reason_raw",
        "test_reason",
        "s_gene_status",
        "policy_period",
        "policy_period_label",
        "policy_era",
        "is_reinfection",
        "vaccination_date",
        "vacc_dose_number",
        "is_vaccinated",
        "vacc_product_name",
        "vacc_booster",
    ]
    meta = meta[[c for c in output_cols if c in meta.columns]].copy()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    meta.to_parquet(out_path, index=False, compression="zstd")
    logging.info("Sequence metadata: %d rows \u2192 %s", len(meta), out_path)
    return meta


def main() -> int:
    """Run all raw-to-processed metadata preparation steps."""
    ap = argparse.ArgumentParser(
        description="Prepare all processed metadata from raw inputs."
    )
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="Repo root (paths in config.yaml are relative to this)",
    )
    ap.add_argument("--log-level", default="INFO")
    
    args = ap.parse_args()

    setup_logging(args.log_level)
    cfg = load_config(args.root / args.config)
    raw = {k: args.root / v for k, v in cfg["data"]["raw"].items()}
    proc = {k: args.root / v for k, v in cfg["data"]["processed"].items()}

    policy_cfg = cfg["policy"]
    policy = prep_policy(
        raw["oxcgrt_stringency_csv"],
        raw["oxcgrt_containment_csv"],
        policy_cfg["periods"],
        proc["policy"],
        region_name=policy_cfg.get("region_name", "Scotland"),
    )

    _ = prep_testing(raw["testing_csv"], proc["testing"])
    _ = prep_vaccination(raw["vaccination_csv"], proc["vaccination"])
    simd = prep_simd(raw["simd_csv"], proc["simd"])
    geography = prep_geography(raw["geography_shp"], simd, proc["geography"])
    _ = prep_hb_trends(raw["daily_hb_trends_csv"], proc["hb_trends"])
    _ = prep_sequence_metadata(
        metadata_csv=raw["metadata_csv"],
        nextclade_tsv=raw["nextclade_tsv"],
        vaccination_csv=raw["vaccination_csv"],
        testing_csv=raw["testing_csv"],
        geography=geography,
        policy=policy,
        out_path=proc["metadata"],
        required_nextclade_qc=cfg["tn93"]["nextclade_qc"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
