"""Part 2 descriptive cluster-characterisation tables.

This script rebuilds primary-resolution cluster aggregates from the sequence-
level analysis dataset, adding vaccination, demographic, SIMD-domain, and
cluster-category summaries that are not present in the cached Part 1 cluster
table.

Run from the repository root with:

    conda run -n PhD python part2/cluster_characterisation.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


PRIMARY_RESOLUTION = 0.3
QC_DEFAULT = "good"
DEFAULT_SIMD_RANK_MAX = 6976.0

AGE_BAND_ORDER = [
    "00-04",
    "05-09",
    "10-14",
    "15-19",
    "20-24",
    "25-29",
    "30-34",
    "35-39",
    "40-44",
    "45-49",
    "50-54",
    "55-59",
    "60-64",
    "65-69",
    "70-74",
    "75+",
]

VACCINATION_AGE_GROUP_ORDER = [
    "00-14",
    "15-19",
    "20-29",
    "30-39",
    "40-49",
    "50-54",
    "55-59",
    "60-64",
    "65-69",
    "70-74",
    "75+",
]
AGE_GROUP_ORDER = VACCINATION_AGE_GROUP_ORDER
SEX_ORDER = ["Female", "Male"]
SIMD_LABELS = {
    1: "1 most deprived",
    2: "2",
    3: "3",
    4: "4",
    5: "5 least deprived",
}
SIMD_LABEL_ORDER = [SIMD_LABELS[i] for i in range(1, 6)]
SIMD_DOMAIN_SPECS = {
    "overall": {
        "label": "Overall SIMD",
        "rank_col": "dz_simd_rank",
        "mean_col": "mean_simd_rank",
        "quintile_col": "simd_quintile",
        "label_col": "simd_quintile_label",
    },
    "income": {
        "label": "Income",
        "rank_col": "dz_simd_income_rank",
        "mean_col": "mean_simd_income_rank",
        "quintile_col": "simd_income_quintile",
        "label_col": "simd_income_quintile_label",
    },
    "employment": {
        "label": "Employment",
        "rank_col": "dz_simd_employment_rank",
        "mean_col": "mean_simd_employment_rank",
        "quintile_col": "simd_employment_quintile",
        "label_col": "simd_employment_quintile_label",
    },
    "education": {
        "label": "Education",
        "rank_col": "dz_simd_education_rank",
        "mean_col": "mean_simd_education_rank",
        "quintile_col": "simd_education_quintile",
        "label_col": "simd_education_quintile_label",
    },
    "health": {
        "label": "Health",
        "rank_col": "dz_simd_health_rank",
        "mean_col": "mean_simd_health_rank",
        "quintile_col": "simd_health_quintile",
        "label_col": "simd_health_quintile_label",
    },
    "access": {
        "label": "Geographic access",
        "rank_col": "dz_simd_access_rank",
        "mean_col": "mean_simd_access_rank",
        "quintile_col": "simd_access_quintile",
        "label_col": "simd_access_quintile_label",
    },
    "crime": {
        "label": "Crime",
        "rank_col": "dz_simd_crime_rank",
        "mean_col": "mean_simd_crime_rank",
        "quintile_col": "simd_crime_quintile",
        "label_col": "simd_crime_quintile_label",
    },
    "housing": {
        "label": "Housing",
        "rank_col": "dz_simd_housing_rank",
        "mean_col": "mean_simd_housing_rank",
        "quintile_col": "simd_housing_quintile",
        "label_col": "simd_housing_quintile_label",
    },
}
SIMD_DOMAIN_ORDER = list(SIMD_DOMAIN_SPECS)

WAVE_ORDER = [
    "B.1.177",
    "Alpha",
    "Delta",
    "BA.1",
    "BA.2",
    "BA.4",
    "BA.5",
    "BQ.1",
    "XBB",
    "Other",
]

VACCINATION_MIXING_LABELS = ("homogeneous", "baseline", "mixed", "not available")
SEQUENCE_COLUMNS = [
    "cluster_id",
    "sequence_id",
    "resolution",
    "window_id",
    "window_idx",
    "wn_mid_date",
    "collection_date",
    "datazone",
    "pango_lineage",
    "nextclade_qc",
    "age_band",
    "age_midpoint",
    "sex",
    "is_female",
    "dz_simd_rank",
    "dz_simd_quintile",
    "dz_simd_income_rank",
    "dz_simd_employment_rank",
    "dz_simd_education_rank",
    "dz_simd_health_rank",
    "dz_simd_access_rank",
    "dz_simd_crime_rank",
    "dz_simd_housing_rank",
    "is_vaccinated",
    "vacc_dose_number",
    "vacc_date_prior",
    "vacc_product_name",
    "vacc_booster",
    "days_since_vaccination",
]

CLUSTER_COMPACT_COLUMNS = [
    "cluster_id",
    "resolution",
    "window_id",
    "window_idx",
    "wn_mid_date",
    "cluster_week",
    "cluster_start_date",
    "cluster_end_date",
    "pango_lineage",
    "wave_group",
    "cluster_size",
    "cluster_size_category",
    "cluster_n_datazones",
    "geographic_dispersion_category",
    "mean_simd_rank",
    "simd_quintile",
    "simd_quintile_label",
    "mean_simd_income_rank",
    "simd_income_quintile",
    "simd_income_quintile_label",
    "mean_simd_employment_rank",
    "simd_employment_quintile",
    "simd_employment_quintile_label",
    "mean_simd_education_rank",
    "simd_education_quintile",
    "simd_education_quintile_label",
    "mean_simd_health_rank",
    "simd_health_quintile",
    "simd_health_quintile_label",
    "mean_simd_access_rank",
    "simd_access_quintile",
    "simd_access_quintile_label",
    "mean_simd_crime_rank",
    "simd_crime_quintile",
    "simd_crime_quintile_label",
    "mean_simd_housing_rank",
    "simd_housing_quintile",
    "simd_housing_quintile_label",
    "mean_age_midpoint",
    "predominant_age_band",
    "predominant_sex",
    "cluster_prop_female",
    "n_vaccination_known",
    "n_vaccinated",
    "cluster_prop_vaccinated",
    "cluster_vaccination_profile",
    "vaccination_n_valid",
    "vaccination_discordance",
    "vaccination_expected_discordance",
    "vaccination_excess_discordance",
    "vaccination_mixing_category",
    "mean_dose_all_members",
    "mean_dose_vaccinated_members",
    "median_dose_vaccinated_members",
    "n_boosted",
    "cluster_prop_boosted_all_members",
    "cluster_prop_boosted_vaccinated_members",
    "mean_days_since_vaccination",
    "median_days_since_vaccination",
    "index_is_vaccinated",
    "index_vacc_dose_number",
    "index_days_since_vaccination",
]


def _bootstrap_repo_root_for_utils() -> Path:
    """Ensure the repository root is importable for direct script execution."""
    here = Path(__file__).resolve()
    for candidate in [here, *here.parents]:
        if (candidate / "config.yaml").exists():
            root_str = str(candidate)
            if root_str not in sys.path:
                sys.path.insert(0, root_str)
            return candidate
    raise FileNotFoundError("Could not locate config.yaml.")


_bootstrap_repo_root_for_utils()

from utils.data import load_analysis_columns_pandas, repo_root as data_repo_root  # noqa: E402
from categorise_clusters import (  # noqa: E402
    GEOGRAPHY_LABELS,
    SIZE_LABELS,
    categorise_count,
    resolve_thresholds,
    simd_quintile_from_rank,
)


ROOT = data_repo_root()


def resolve_repo_path(root: Path, path: Path | None, default: Path) -> Path:
    target = path if path is not None else default
    return target if target.is_absolute() else root / target


def assign_wave(lineage: str) -> str:
    """Assign a broad variant wave from Pango lineage, matching Part 1."""
    if not isinstance(lineage, str):
        return "Other"
    if lineage.startswith("B.1.177"):
        return "B.1.177"
    if lineage == "B.1.1.7" or lineage.startswith("B.1.1.7."):
        return "Alpha"
    if lineage.startswith("AY.") or lineage == "B.1.617.2":
        return "Delta"
    if lineage.startswith("BA.1"):
        return "BA.1"
    if lineage.startswith("BA.2"):
        return "BA.2"
    if lineage.startswith("BA.4"):
        return "BA.4"
    if lineage.startswith("BA.5") or lineage.startswith("BE."):
        return "BA.5"
    if lineage.startswith("BQ."):
        return "BQ.1"
    if lineage.startswith("XBB"):
        return "XBB"
    return "Other"


def vaccination_age_group(age_band: object) -> str | pd.NA:
    """Approximate JCVI rollout age groups from the available age-band field."""
    if not isinstance(age_band, str):
        return pd.NA
    if age_band in {"00-04", "05-09", "10-14"}:
        return "00-14"
    if age_band == "15-19":
        return "15-19"
    if age_band in {"20-24", "25-29"}:
        return "20-29"
    if age_band in {"30-34", "35-39"}:
        return "30-39"
    if age_band in {"40-44", "45-49"}:
        return "40-49"
    if age_band in {
        "50-54",
        "55-59",
        "60-64",
        "65-69",
        "70-74",
        "75+",
    }:
        return age_band
    return pd.NA


def age_group(age_band: object) -> str | pd.NA:
    """Backward-compatible alias for the vaccination rollout grouping."""
    return vaccination_age_group(age_band)


def mode_or_na(values: pd.Series) -> object:
    clean = values.dropna()
    if clean.empty:
        return pd.NA
    counts = clean.astype(str).value_counts(sort=True)
    return counts.index[0]


def week_start(values: pd.Series) -> pd.Series:
    return pd.to_datetime(values).dt.to_period("W-SUN").dt.start_time


def simd_quintile_labels(quintile: pd.Series) -> pd.Categorical:
    return pd.Categorical(
        quintile.map(SIMD_LABELS),
        categories=SIMD_LABEL_ORDER,
        ordered=True,
    )


def add_sequence_simd_domain_quintiles(
    seq: pd.DataFrame,
    *,
    simd_rank_max: float,
) -> pd.DataFrame:
    for spec in SIMD_DOMAIN_SPECS.values():
        quintile = simd_quintile_from_rank(seq[spec["rank_col"]], simd_rank_max)
        seq[spec["quintile_col"]] = quintile
        seq[spec["label_col"]] = simd_quintile_labels(quintile)
    return seq


def add_cluster_simd_domain_quintiles(
    clusters: pd.DataFrame,
    *,
    simd_rank_max: float,
) -> pd.DataFrame:
    for spec in SIMD_DOMAIN_SPECS.values():
        quintile = simd_quintile_from_rank(clusters[spec["mean_col"]], simd_rank_max)
        clusters[spec["quintile_col"]] = quintile
        clusters[spec["label_col"]] = simd_quintile_labels(quintile)
    return clusters


def read_sequence_rows(
    qc: str | None,
    primary_resolution: float,
    simd_rank_max: float,
) -> pd.DataFrame:
    seq = load_analysis_columns_pandas(
        columns=SEQUENCE_COLUMNS,
        resolution=primary_resolution,
        qc=qc,
    )
    seq["collection_date"] = pd.to_datetime(seq["collection_date"])
    seq["wn_mid_date"] = pd.to_datetime(seq["wn_mid_date"])
    seq["vacc_date_prior"] = pd.to_datetime(seq["vacc_date_prior"])
    seq["is_vaccinated"] = pd.to_numeric(seq["is_vaccinated"], errors="coerce")
    seq["vacc_dose_number"] = pd.to_numeric(seq["vacc_dose_number"], errors="coerce")
    seq["vacc_booster"] = pd.to_numeric(seq["vacc_booster"], errors="coerce")
    seq["days_since_vaccination"] = pd.to_numeric(
        seq["days_since_vaccination"],
        errors="coerce",
    )
    seq["wave_group"] = pd.Categorical(
        seq["pango_lineage"].astype(str).map(assign_wave),
        categories=WAVE_ORDER,
        ordered=True,
    )
    seq["case_week"] = week_start(seq["collection_date"])
    seq["vaccination_age_group"] = pd.Categorical(
        seq["age_band"].map(vaccination_age_group),
        categories=VACCINATION_AGE_GROUP_ORDER,
        ordered=True,
    )
    seq["age_group"] = seq["vaccination_age_group"]
    seq = add_sequence_simd_domain_quintiles(seq, simd_rank_max=simd_rank_max)
    return seq


def deduplicate_cases(seq: pd.DataFrame) -> pd.DataFrame:
    """Return one row per sequence/case for case-level vaccination summaries."""
    return (
        seq.sort_values(["sequence_id", "window_idx", "collection_date"])
        .drop_duplicates("sequence_id")
        .copy()
    )


def add_cluster_categories(
    clusters: pd.DataFrame,
    *,
    large_size_min: int,
    very_large_size_min: int,
    large_geography_min: int,
    very_large_geography_min: int,
    simd_rank_max: float,
) -> pd.DataFrame:
    out = clusters.copy()
    out["cluster_size_category"] = categorise_count(
        out["cluster_size"],
        large_min=large_size_min,
        very_large_min=very_large_size_min,
        labels=SIZE_LABELS,
    )
    out["geographic_dispersion_category"] = categorise_count(
        out["cluster_n_datazones"],
        large_min=large_geography_min,
        very_large_min=very_large_geography_min,
        labels=GEOGRAPHY_LABELS,
    )
    singleton_mask = out["cluster_size"] <= 1
    out.loc[singleton_mask, "cluster_size_category"] = pd.NA
    out.loc[singleton_mask, "geographic_dispersion_category"] = pd.NA
    out = add_cluster_simd_domain_quintiles(out, simd_rank_max=simd_rank_max)
    return out


def pairwise_discordance_from_counts(
    counts: pd.DataFrame,
    group_cols: list[str],
) -> pd.DataFrame:
    totals = counts.groupby(group_cols, observed=True)["n"].sum().rename("n_valid")
    same_pairs = (
        counts.assign(same_pairs=counts["n"] * (counts["n"] - 1))
        .groupby(group_cols, observed=True)["same_pairs"]
        .sum()
    )
    out = pd.concat([totals, same_pairs], axis=1).reset_index()
    denom = out["n_valid"] * (out["n_valid"] - 1)
    out["discordance"] = np.where(denom > 0, 1 - out["same_pairs"] / denom, np.nan)
    return out.drop(columns=["same_pairs"])


def observed_vaccination_discordance(seq: pd.DataFrame) -> pd.DataFrame:
    counts = (
        seq.dropna(subset=["is_vaccinated"])
        .assign(is_vaccinated=lambda df: df["is_vaccinated"].astype(int))
        .groupby(["cluster_id", "is_vaccinated"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = pairwise_discordance_from_counts(counts, ["cluster_id"])
    return out.rename(
        columns={
            "n_valid": "vaccination_n_valid",
            "discordance": "vaccination_discordance",
        }
    )


def expected_vaccination_discordance(seq: pd.DataFrame) -> pd.DataFrame:
    strata = ["window_id", "pango_lineage"]
    counts = (
        seq.dropna(subset=["is_vaccinated"])
        .assign(is_vaccinated=lambda df: df["is_vaccinated"].astype(int))
        .groupby([*strata, "is_vaccinated"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = pairwise_discordance_from_counts(counts, strata)
    return out.rename(
        columns={
            "n_valid": "vaccination_stratum_n_valid",
            "discordance": "vaccination_expected_discordance",
        }
    )


def categorise_vaccination_mixing(
    excess_discordance: pd.Series,
    baseline_band: float,
) -> pd.Categorical:
    """Classify vaccination-status mixing relative to lineage-window expectation."""
    if baseline_band < 0:
        raise ValueError("vaccination_mixing_baseline_band must be non-negative.")
    values = pd.to_numeric(excess_discordance, errors="coerce")
    out = pd.Series("not available", index=values.index, dtype="object")
    out[values < -baseline_band] = "homogeneous"
    out[values.abs() <= baseline_band] = "baseline"
    out[values > baseline_band] = "mixed"
    return pd.Categorical(out, categories=list(VACCINATION_MIXING_LABELS), ordered=True)


def build_cluster_table(
    seq: pd.DataFrame,
    *,
    large_size_min: int,
    very_large_size_min: int,
    large_geography_min: int,
    very_large_geography_min: int,
    simd_rank_max: float,
    vaccination_mixing_baseline_band: float,
) -> pd.DataFrame:
    grouped = seq.groupby("cluster_id", observed=True, sort=False)
    cluster_aggs = {
        "cluster_size": ("sequence_id", "nunique"),
        "cluster_n_datazones": ("datazone", "nunique"),
        "cluster_start_date": ("collection_date", "min"),
        "cluster_end_date": ("collection_date", "max"),
        "resolution": ("resolution", "first"),
        "window_id": ("window_id", "first"),
        "window_idx": ("window_idx", "first"),
        "wn_mid_date": ("wn_mid_date", "first"),
        "pango_lineage": ("pango_lineage", "first"),
        "wave_group": ("wave_group", "first"),
        "mean_age_midpoint": ("age_midpoint", "mean"),
        "predominant_age_band": ("age_band", mode_or_na),
        "predominant_sex": ("sex", mode_or_na),
        "cluster_prop_female": ("is_female", "mean"),
        "predominant_simd_quintile": ("dz_simd_quintile", mode_or_na),
        "n_vaccination_known": ("is_vaccinated", "count"),
        "n_vaccinated": ("is_vaccinated", "sum"),
        "cluster_prop_vaccinated": ("is_vaccinated", "mean"),
        "mean_dose_all_members": ("vacc_dose_number", "mean"),
        "n_boosted": ("vacc_booster", "sum"),
    }
    for spec in SIMD_DOMAIN_SPECS.values():
        cluster_aggs[spec["mean_col"]] = (spec["rank_col"], "mean")
    clusters = grouped.agg(**cluster_aggs).reset_index()

    vaccinated = seq.loc[seq["is_vaccinated"] == 1].copy()
    vacc_agg = (
        vaccinated.groupby("cluster_id", observed=True)
        .agg(
            mean_dose_vaccinated_members=("vacc_dose_number", "mean"),
            median_dose_vaccinated_members=("vacc_dose_number", "median"),
            mean_days_since_vaccination=("days_since_vaccination", "mean"),
            median_days_since_vaccination=("days_since_vaccination", "median"),
            cluster_prop_boosted_vaccinated_members=("vacc_booster", "mean"),
        )
        .reset_index()
    )
    clusters = clusters.merge(vacc_agg, on="cluster_id", how="left")

    obs_mix = observed_vaccination_discordance(seq)
    exp_mix = expected_vaccination_discordance(seq)
    clusters = clusters.merge(obs_mix, on="cluster_id", how="left")
    clusters = clusters.merge(exp_mix, on=["window_id", "pango_lineage"], how="left")
    clusters["vaccination_excess_discordance"] = (
        clusters["vaccination_discordance"]
        - clusters["vaccination_expected_discordance"]
    )
    clusters["vaccination_mixing_category"] = categorise_vaccination_mixing(
        clusters["vaccination_excess_discordance"],
        baseline_band=vaccination_mixing_baseline_band,
    )

    index_rows = (
        seq.sort_values(["cluster_id", "collection_date", "sequence_id"])
        .groupby("cluster_id", observed=True)
        .first()
        .reset_index()
    )
    index_rows = index_rows[
        [
            "cluster_id",
            "is_vaccinated",
            "vacc_dose_number",
            "vacc_booster",
            "days_since_vaccination",
        ]
    ].rename(
        columns={
            "is_vaccinated": "index_is_vaccinated",
            "vacc_dose_number": "index_vacc_dose_number",
            "vacc_booster": "index_vacc_booster",
            "days_since_vaccination": "index_days_since_vaccination",
        }
    )
    clusters = clusters.merge(index_rows, on="cluster_id", how="left")

    clusters["cluster_week"] = week_start(clusters["wn_mid_date"])
    clusters["n_vaccinated"] = clusters["n_vaccinated"].astype(int)
    clusters["n_boosted"] = clusters["n_boosted"].fillna(0).astype(int)
    clusters["cluster_prop_boosted_all_members"] = (
        clusters["n_boosted"] / clusters["cluster_size"]
    )

    clusters["cluster_vaccination_profile"] = pd.Categorical(
        np.select(
            [
                clusters["cluster_prop_vaccinated"].eq(0),
                clusters["cluster_prop_vaccinated"].eq(1),
                clusters["cluster_prop_vaccinated"].between(0, 1, inclusive="neither"),
            ],
            ["none vaccinated", "all vaccinated", "mixed vaccination"],
            default="vaccination unknown",
        ),
        categories=[
            "none vaccinated",
            "mixed vaccination",
            "all vaccinated",
            "vaccination unknown",
        ],
        ordered=True,
    )

    return add_cluster_categories(
        clusters,
        large_size_min=large_size_min,
        very_large_size_min=very_large_size_min,
        large_geography_min=large_geography_min,
        very_large_geography_min=very_large_geography_min,
        simd_rank_max=simd_rank_max,
    )


def summarise_vaccination(grouped, id_col: str) -> pd.DataFrame:
    base = grouped.agg(
        n_cases=(id_col, "nunique"),
        n_vaccination_known=("is_vaccinated", "count"),
        n_vaccinated=("is_vaccinated", "sum"),
        mean_dose_all_cases=("vacc_dose_number", "mean"),
        n_boosted=("vacc_booster", "sum"),
    ).reset_index()
    base["n_vaccinated"] = base["n_vaccinated"].astype(int)
    base["n_boosted"] = base["n_boosted"].fillna(0).astype(int)
    base["prop_vaccinated"] = base["n_vaccinated"] / base["n_vaccination_known"]
    base["prop_boosted_all_cases"] = base["n_boosted"] / base["n_cases"]
    base["prop_boosted_vaccinated_cases"] = np.divide(
        base["n_boosted"],
        base["n_vaccinated"],
        out=np.full(len(base), np.nan, dtype=float),
        where=base["n_vaccinated"].to_numpy() > 0,
    )
    return base


def weekly_case_summary(cases: pd.DataFrame) -> pd.DataFrame:
    frames = []
    specs = [
        ("overall", None, None),
        ("age_band", "age_band", AGE_BAND_ORDER),
        (
            "vaccination_age_group",
            "vaccination_age_group",
            VACCINATION_AGE_GROUP_ORDER,
        ),
        ("sex", "sex", SEX_ORDER),
        ("simd_quintile", "simd_quintile_label", SIMD_LABEL_ORDER),
    ]
    for stratum_type, col, order in specs:
        if col is None:
            work = cases.assign(stratum_type=stratum_type, stratum="Overall")
            grouped = work.groupby(["case_week", "stratum_type", "stratum"], observed=True)
        else:
            work = cases.dropna(subset=[col]).copy()
            work["stratum_type"] = stratum_type
            work["stratum"] = work[col].astype(str)
            grouped = work.groupby(["case_week", "stratum_type", "stratum"], observed=True)
        summary = summarise_vaccination(grouped, "sequence_id")
        if order is not None:
            summary["stratum"] = pd.Categorical(summary["stratum"], categories=order, ordered=True)
        frames.append(summary)
    return pd.concat(frames, ignore_index=True).sort_values(
        ["stratum_type", "stratum", "case_week"]
    )


def case_weekly_simd_domain_summary(cases: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for domain in SIMD_DOMAIN_ORDER:
        spec = SIMD_DOMAIN_SPECS[domain]
        work = cases.dropna(subset=[spec["label_col"]]).copy()
        work["simd_domain"] = domain
        work["simd_domain_label"] = spec["label"]
        work["simd_quintile"] = work[spec["quintile_col"]]
        work["simd_quintile_label"] = work[spec["label_col"]]
        grouped = work.groupby(
            [
                "case_week",
                "simd_domain",
                "simd_domain_label",
                "simd_quintile",
                "simd_quintile_label",
            ],
            observed=True,
        )
        frames.append(summarise_vaccination(grouped, "sequence_id"))
    return pd.concat(frames, ignore_index=True).sort_values(
        ["simd_domain", "simd_quintile", "case_week"]
    )


def cluster_wave_category_summary(clusters: pd.DataFrame) -> pd.DataFrame:
    frames = []
    specs = [
        ("cluster_size_category", "cluster_size_category"),
        ("geographic_dispersion_category", "geographic_dispersion_category"),
        ("cluster_vaccination_profile", "cluster_vaccination_profile"),
        ("vaccination_mixing_category", "vaccination_mixing_category"),
        ("simd_quintile", "simd_quintile_label"),
    ]
    for category_variable, col in specs:
        work = clusters.dropna(subset=[col]).copy()
        work["category_variable"] = category_variable
        work["category"] = work[col].astype(str)
        grouped = work.groupby(["wave_group", "category_variable", "category"], observed=True)
        summary = grouped.agg(
            n_clusters=("cluster_id", "size"),
            n_cases=("cluster_size", "sum"),
            mean_cluster_size=("cluster_size", "mean"),
            median_cluster_size=("cluster_size", "median"),
            mean_datazones=("cluster_n_datazones", "mean"),
            median_datazones=("cluster_n_datazones", "median"),
            mean_prop_vaccinated=("cluster_prop_vaccinated", "mean"),
            median_prop_vaccinated=("cluster_prop_vaccinated", "median"),
            mean_dose_vaccinated_members=("mean_dose_vaccinated_members", "mean"),
            mean_days_since_vaccination=("mean_days_since_vaccination", "mean"),
            mean_prop_boosted_vaccinated_members=(
                "cluster_prop_boosted_vaccinated_members",
                "mean",
            ),
        ).reset_index()
        frames.append(summary)
    return pd.concat(frames, ignore_index=True).sort_values(
        ["category_variable", "wave_group", "category"]
    )


def cluster_simd_domain_wave_summary(clusters: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for domain in SIMD_DOMAIN_ORDER:
        spec = SIMD_DOMAIN_SPECS[domain]
        work = clusters.dropna(subset=[spec["label_col"]]).copy()
        work["simd_domain"] = domain
        work["simd_domain_label"] = spec["label"]
        work["simd_quintile"] = work[spec["quintile_col"]]
        work["simd_quintile_label"] = work[spec["label_col"]]
        grouped = work.groupby(
            [
                "wave_group",
                "simd_domain",
                "simd_domain_label",
                "simd_quintile",
                "simd_quintile_label",
            ],
            observed=True,
        )
        summary = grouped.agg(
            n_clusters=("cluster_id", "size"),
            n_cases=("cluster_size", "sum"),
            mean_cluster_size=("cluster_size", "mean"),
            median_cluster_size=("cluster_size", "median"),
            mean_datazones=("cluster_n_datazones", "mean"),
            median_datazones=("cluster_n_datazones", "median"),
            mean_prop_vaccinated=("cluster_prop_vaccinated", "mean"),
            median_prop_vaccinated=("cluster_prop_vaccinated", "median"),
            mean_dose_vaccinated_members=("mean_dose_vaccinated_members", "mean"),
            mean_days_since_vaccination=("mean_days_since_vaccination", "mean"),
            mean_prop_boosted_vaccinated_members=(
                "cluster_prop_boosted_vaccinated_members",
                "mean",
            ),
        ).reset_index()
        frames.append(summary)
    return pd.concat(frames, ignore_index=True).sort_values(
        ["simd_domain", "wave_group", "simd_quintile"]
    )


def cluster_weekly_category_summary(clusters: pd.DataFrame) -> pd.DataFrame:
    frames = []
    specs = [
        ("cluster_size_category", "cluster_size_category"),
        ("geographic_dispersion_category", "geographic_dispersion_category"),
        ("cluster_vaccination_profile", "cluster_vaccination_profile"),
        ("vaccination_mixing_category", "vaccination_mixing_category"),
    ]
    for category_variable, col in specs:
        work = clusters.dropna(subset=[col]).copy()
        work["category_variable"] = category_variable
        work["category"] = work[col].astype(str)
        grouped = work.groupby(["cluster_week", "category_variable", "category"], observed=True)
        summary = grouped.agg(
            n_clusters=("cluster_id", "size"),
            n_cases=("cluster_size", "sum"),
            mean_cluster_size=("cluster_size", "mean"),
            mean_prop_vaccinated=("cluster_prop_vaccinated", "mean"),
            median_prop_vaccinated=("cluster_prop_vaccinated", "median"),
        ).reset_index()
        summary["fraction_clusters_within_week"] = summary["n_clusters"] / summary.groupby(
            ["cluster_week", "category_variable"],
            observed=True,
        )["n_clusters"].transform("sum")
        frames.append(summary)
    return pd.concat(frames, ignore_index=True).sort_values(
        ["category_variable", "category", "cluster_week"]
    )


def cluster_descriptives(clusters: pd.DataFrame, cases: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {"measure": "primary_resolution_case_rows", "value": int(len(cases))},
        {"measure": "primary_resolution_windowed_sequence_rows", "value": int(clusters["cluster_size"].sum())},
        {"measure": "clusters", "value": int(len(clusters))},
        {"measure": "singleton_clusters", "value": int((clusters["cluster_size"] == 1).sum())},
        {"measure": "mean_case_prop_vaccinated", "value": float(cases["is_vaccinated"].mean())},
        {
            "measure": "mean_cluster_prop_vaccinated",
            "value": float(clusters["cluster_prop_vaccinated"].mean()),
        },
        {
            "measure": "mean_cluster_prop_boosted_all_members",
            "value": float(clusters["cluster_prop_boosted_all_members"].mean()),
        },
        {
            "measure": "mean_cluster_dose_vaccinated_members",
            "value": float(clusters["mean_dose_vaccinated_members"].mean()),
        },
        {
            "measure": "median_cluster_days_since_vaccination",
            "value": float(clusters["median_days_since_vaccination"].median()),
        },
    ]
    for category, count in clusters["vaccination_mixing_category"].value_counts(sort=False).items():
        rows.append({"measure": f"vaccination_mixing_clusters_{category}", "value": int(count)})
    for wave, count in clusters["wave_group"].value_counts(sort=False).items():
        rows.append({"measure": f"clusters_wave_{wave}", "value": int(count)})
    return pd.DataFrame(rows)


def key_questions_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "priority": 1,
                "question": (
                    "Among sequenced SARS-CoV-2 cases, how did breakthrough-case "
                    "frequency change over calendar time, and did this differ by "
                    "vaccination-rollout age group, sex, overall SIMD quintile, "
                    "or SIMD domain quintile?"
                ),
                "why_it_matters": (
                    "Separates vaccine rollout/booster timing from demographic and "
                    "deprivation differences in the infected sequenced population."
                ),
                "suggested_outputs": (
                    "Weekly vaccinated-case proportions by JCVI rollout-informed age "
                    "group, sex, overall SIMD quintile, and SIMD domain quintile."
                ),
                "main_caveat": (
                    "This is among sequenced cases, not vaccine effectiveness in the population."
                ),
            },
            {
                "priority": 2,
                "question": (
                    "Were clusters more homogeneous or mixed by vaccination status than expected "
                    "for cases from the same lineage and calendar window?"
                ),
                "why_it_matters": (
                    "Treats vaccination as a cluster-characterisation variable rather than "
                    "assuming binary vaccination proportion should explain size."
                ),
                "suggested_outputs": (
                    "Observed-minus-expected vaccination-status discordance and homogeneous/"
                    "baseline/mixed categories by wave, size, geography, overall SIMD, "
                    "SIMD domains, and age."
                ),
                "main_caveat": (
                    "Expected mixing should be interpreted as a descriptive benchmark, not a "
                    "causal counterfactual."
                ),
            },
            {
                "priority": 3,
                "question": (
                    "Did vaccine dose profile and recency among infected cases vary across "
                    "overall SIMD quintiles, SIMD domain quintiles, or epidemic waves?"
                ),
                "why_it_matters": (
                    "Can reveal inequality in vaccine exposure among cases and distinguish "
                    "early-dose from booster-era breakthrough clusters."
                ),
                "suggested_outputs": (
                    "Mean dose among vaccinated cluster members, booster proportion, and days "
                    "since latest prior dose by wave/SIMD domain/cluster category."
                ),
                "main_caveat": (
                    "Dose timing is both a biological exposure and a proxy for rollout phase."
                ),
            },
            {
                "priority": 4,
                "question": (
                    "Do all-vaccinated, mixed-vaccination, and unvaccinated clusters have "
                    "different social or geographic mixing profiles?"
                ),
                "why_it_matters": (
                    "Connects vaccination with the Part 1 transmission-structure results."
                ),
                "suggested_outputs": (
                    "Join this vaccination cluster table to mixing metrics or recompute mixing in "
                    "Part 2, then compare categories by wave."
                ),
                "main_caveat": (
                    "Requires careful handling of singleton clusters where mixing is undefined."
                ),
            },
        ]
    )


def compact_cluster_table(clusters: pd.DataFrame) -> pd.DataFrame:
    cols = [col for col in CLUSTER_COMPACT_COLUMNS if col in clusters.columns]
    return clusters.loc[:, cols].copy()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build Part 2 cluster-characterisation tables from "
            "primary-resolution sequence data."
        ),
    )
    parser.add_argument("--qc", default=QC_DEFAULT, help="Nextclade QC filter. Use 'none' to skip.")
    parser.add_argument("--primary-resolution", type=float, default=PRIMARY_RESOLUTION)
    parser.add_argument("--out-dir", type=Path, default=None, help="Output directory. Default: part2.")
    parser.add_argument("--large-size-min", type=int, default=None)
    parser.add_argument("--very-large-size-min", type=int, default=None)
    parser.add_argument("--large-size-quantile", type=float, default=0.90)
    parser.add_argument("--very-large-size-quantile", type=float, default=0.99)
    parser.add_argument("--large-geography-min", type=int, default=None)
    parser.add_argument("--very-large-geography-min", type=int, default=None)
    parser.add_argument("--large-geography-quantile", type=float, default=0.90)
    parser.add_argument("--very-large-geography-quantile", type=float, default=0.99)
    parser.add_argument("--simd-rank-max", type=float, default=DEFAULT_SIMD_RANK_MAX)
    parser.add_argument(
        "--vaccination-mixing-baseline-band",
        type=float,
        default=0.01,
        help=(
            "Absolute observed-minus-expected vaccination discordance classified "
            "as baseline. Default: 0.01, i.e. +/-1 percentage point."
        ),
    )
    return parser.parse_args()


def run(args: argparse.Namespace) -> None:
    out_dir = resolve_repo_path(ROOT, args.out_dir, Path("part2"))
    cache_dir = out_dir / "cache"
    tables_dir = out_dir / "tables"
    cache_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    qc = None if str(args.qc).lower() == "none" else args.qc
    print("Reading primary-resolution sequence rows with vaccination fields", flush=True)
    seq = read_sequence_rows(
        qc=qc,
        primary_resolution=args.primary_resolution,
        simd_rank_max=args.simd_rank_max,
    )
    cases = deduplicate_cases(seq)
    print(
        f"Loaded {len(seq):,} windowed sequence rows and {len(cases):,} unique cases.",
        flush=True,
    )

    cluster_sizes = seq.groupby("cluster_id", observed=True)["sequence_id"].nunique()
    cluster_datazones = seq.groupby("cluster_id", observed=True)["datazone"].nunique()
    non_singleton_cluster_ids = cluster_sizes.index[cluster_sizes > 1]
    if len(non_singleton_cluster_ids) == 0:
        raise ValueError("No non-singleton clusters available for categorisation.")
    non_singleton_sizes = cluster_sizes.loc[non_singleton_cluster_ids]
    non_singleton_datazones = cluster_datazones.loc[non_singleton_cluster_ids]
    large_size_min, very_large_size_min = resolve_thresholds(
        non_singleton_sizes,
        large_min=args.large_size_min,
        very_large_min=args.very_large_size_min,
        large_quantile=args.large_size_quantile,
        very_large_quantile=args.very_large_size_quantile,
    )
    large_geography_min, very_large_geography_min = resolve_thresholds(
        non_singleton_datazones,
        large_min=args.large_geography_min,
        very_large_min=args.very_large_geography_min,
        large_quantile=args.large_geography_quantile,
        very_large_quantile=args.very_large_geography_quantile,
    )
    print(
        "Building cluster-level vaccination aggregates "
        f"(size thresholds {large_size_min}/{very_large_size_min}; "
        f"geography thresholds {large_geography_min}/{very_large_geography_min}).",
        flush=True,
    )
    clusters = build_cluster_table(
        seq,
        large_size_min=large_size_min,
        very_large_size_min=very_large_size_min,
        large_geography_min=large_geography_min,
        very_large_geography_min=very_large_geography_min,
        simd_rank_max=args.simd_rank_max,
        vaccination_mixing_baseline_band=args.vaccination_mixing_baseline_band,
    )

    categorised_clusters = clusters.loc[clusters["cluster_size"] > 1].copy()
    weekly = weekly_case_summary(cases)
    weekly_simd_domains = case_weekly_simd_domain_summary(cases)
    wave_summary = cluster_wave_category_summary(categorised_clusters)
    wave_simd_domains = cluster_simd_domain_wave_summary(categorised_clusters)
    weekly_cluster = cluster_weekly_category_summary(categorised_clusters)
    descriptives = cluster_descriptives(clusters, cases)
    questions = key_questions_table()

    clusters.to_parquet(cache_dir / "vaccination_cluster_table.parquet", index=False)
    compact_cluster_table(clusters).to_csv(
        tables_dir / "vaccination_cluster_table.csv",
        index=False,
    )
    weekly.to_csv(tables_dir / "vaccination_case_weekly_summary.csv", index=False)
    weekly_simd_domains.to_csv(
        tables_dir / "vaccination_case_weekly_simd_domain_summary.csv",
        index=False,
    )
    wave_summary.to_csv(
        tables_dir / "vaccination_cluster_wave_category_summary.csv",
        index=False,
    )
    wave_simd_domains.to_csv(
        tables_dir / "vaccination_cluster_wave_simd_domain_summary.csv",
        index=False,
    )
    weekly_cluster.to_csv(
        tables_dir / "vaccination_cluster_weekly_category_summary.csv",
        index=False,
    )
    descriptives.to_csv(tables_dir / "vaccination_descriptives.csv", index=False)
    questions.to_csv(tables_dir / "vaccination_key_questions.csv", index=False)

    print(
        f"Wrote vaccination Part 2 tables/cache under {out_dir} "
        f"({len(clusters):,} clusters).",
        flush=True,
    )


if __name__ == "__main__":
    run(parse_args())
