from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .config import (
    AnalysisConfig,
    DATA_PATH,
    DATE_COLUMNS,
    PRIMARY_CATEGORICAL_FEATURES,
    READ_COLUMNS,
    SENSITIVITY_NUMERIC_FEATURES,
    SIMD_DOMAIN_FEATURES_BY_MODE,
    SIMD_DOMAIN_RANK_FEATURES,
    SIMD_OVERALL_FEATURE_BY_MODE,
    TEST_REASON_MAP,
)


@dataclass(slots=True)
class FeatureSpec:
    simd_overall_feature: str
    simd_domain_mode: str
    simd_domain_features: list[str]
    numeric_features: list[str]
    categorical_features: list[str]
    all_features: list[str]
    decomposition_numeric_features: list[str]
    decomposition_features: list[str]


def read_analysis_data(
    data_path: str | Path = DATA_PATH,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """Read the analysis parquet with only the columns needed for modelling."""
    return pd.read_parquet(data_path, columns=READ_COLUMNS if columns is None else columns)


def build_sequence_table(
    source: pd.DataFrame,
    config: AnalysisConfig,
    resolution: float | None = None,
    large_min: int | None = None,
) -> pd.DataFrame:
    """Filter to one resolution/QC pass and de-duplicate overlapping windows."""
    resolution = config.primary_resolution if resolution is None else resolution
    large_min = config.primary_large_min if large_min is None else large_min

    df = source.loc[
        source["resolution"].eq(resolution)
        & source["nextclade_qc"].isin(["good", "mediocre"])
    ].copy()

    for col in DATE_COLUMNS:
        df[col] = pd.to_datetime(df[col])

    df["dist_to_mid"] = (df["collection_date"] - df["wn_mid_date"]).abs()
    df = (
        df.sort_values(["sequence_id", "dist_to_mid", "window_idx"], kind="mergesort")
        .drop_duplicates("sequence_id", keep="first")
        .reset_index(drop=True)
    )

    df["cluster_type"] = np.select(
        [df["cluster_size"].eq(1), df["cluster_size"].between(2, large_min - 1)],
        ["singleton", "small_cluster"],
        default="large_cluster",
    )
    df["cluster_type"] = pd.Categorical(
        df["cluster_type"],
        categories=["singleton", "small_cluster", "large_cluster"],
        ordered=True,
    )
    df["is_large_cluster"] = df["cluster_type"].eq("large_cluster").astype("int8")
    df["is_singleton"] = df["cluster_type"].eq("singleton").astype("int8")
    df["large_min_threshold"] = large_min
    return df


def engineer_features(df: pd.DataFrame) -> pd.DataFrame:
    """Create derived modelling variables and grouped categorical labels."""
    out = df.copy()

    out["days_since_epoch"] = (
        out["collection_date"] - pd.Timestamp("2020-01-01")
    ).dt.days.astype("float64")
    out["log_dz_population_density"] = np.log1p(out["dz_population_density"].clip(lower=0))
    out["log_dz_cum_incidence_per_capita"] = np.log1p(
        out["dz_cum_incidence_per_capita"].clip(lower=0)
    )
    out["log_wn_no_sequences"] = np.log1p(out["wn_no_sequences"].clip(lower=0))
    out["log_hb_hospital_occupancy"] = np.log1p(out["hb_hospital_occupancy"].clip(lower=0))

    out["is_unvaccinated"] = out["is_vaccinated"].fillna(0).eq(0).astype("int8")
    out["days_since_vaccination_model"] = out["days_since_vaccination"].where(
        out["is_unvaccinated"].eq(0), -1
    )
    out["vacc_booster_model"] = out["vacc_booster"].where(out["is_unvaccinated"].eq(0), 0)

    out["vacc_product_name_model"] = out["vacc_product_name"].astype("object")
    out.loc[out["is_unvaccinated"].eq(1), "vacc_product_name_model"] = "None"
    out["vacc_product_name_model"] = out["vacc_product_name_model"].fillna("Unknown")

    out["who_voc_model"] = out["who_voc"].fillna("Other/Non-VOC")
    out["test_reason_category"] = out["test_reason"].map(TEST_REASON_MAP)
    out["test_reason_unmapped"] = (
        out["test_reason"].notna() & out["test_reason_category"].isna()
    )
    out["test_reason_model"] = out["test_reason_category"].fillna("unknown")
    out["test_type_model"] = out["test_type"].fillna("Missing")
    out["dz_health_board_model"] = out["dz_health_board"].fillna("Missing")
    out["dz_urban_rural_class_model"] = out["dz_urban_rural_class"].fillna("Missing")

    # Domain ranks use 1 = most deprived and 6,976 = least deprived.
    # These derived bands preserve that direction: quintile/decile 1 = most deprived.
    for rank_col in SIMD_DOMAIN_RANK_FEATURES:
        stem = rank_col.removesuffix("_rank")
        out[f"{stem}_quintile"] = np.ceil(out[rank_col] / (6976 / 5)).clip(1, 5).astype(int)
        out[f"{stem}_decile"] = np.ceil(out[rank_col] / (6976 / 10)).clip(1, 10).astype(int)
    return out


def build_feature_spec(config: AnalysisConfig) -> FeatureSpec:
    """Return main and SIMD-domain feature lists for a config."""
    if config.simd_overall_mode not in SIMD_OVERALL_FEATURE_BY_MODE:
        raise ValueError(
            "simd_overall_mode must be one of "
            f"{sorted(SIMD_OVERALL_FEATURE_BY_MODE)}; got {config.simd_overall_mode!r}"
        )
    if config.simd_domain_mode not in SIMD_DOMAIN_FEATURES_BY_MODE:
        raise ValueError(
            "simd_domain_mode must be one of "
            f"{sorted(SIMD_DOMAIN_FEATURES_BY_MODE)}; got {config.simd_domain_mode!r}"
        )

    simd_overall_feature = SIMD_OVERALL_FEATURE_BY_MODE[config.simd_overall_mode]
    simd_domain_features = SIMD_DOMAIN_FEATURES_BY_MODE[config.simd_domain_mode]
    numeric_features = [
        "age_midpoint",
        "is_female",
        simd_overall_feature,
        "log_dz_population_density",
        "vacc_dose_number",
        "is_unvaccinated",
        "days_since_vaccination_model",
        "vacc_booster_model",
        "is_reinfection",
        "dz_cum_prop_sequenced",
        "log_dz_cum_incidence_per_capita",
        "dz_7d_test_positivity",
        "wn_prop_sequenced",
        "log_wn_no_sequences",
        "days_since_epoch",
    ]
    if config.use_sensitivity_context_controls:
        numeric_features = numeric_features + SENSITIVITY_NUMERIC_FEATURES

    categorical_features = list(PRIMARY_CATEGORICAL_FEATURES)
    all_features = numeric_features + categorical_features
    decomposition_numeric_features = [
        f for f in numeric_features if f != simd_overall_feature
    ] + list(simd_domain_features)
    decomposition_features = decomposition_numeric_features + categorical_features

    return FeatureSpec(
        simd_overall_feature=simd_overall_feature,
        simd_domain_mode=config.simd_domain_mode,
        simd_domain_features=list(simd_domain_features),
        numeric_features=numeric_features,
        categorical_features=categorical_features,
        all_features=all_features,
        decomposition_numeric_features=decomposition_numeric_features,
        decomposition_features=decomposition_features,
    )


def make_outcome_tables(seq: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Build reusable outcome diagnostic tables."""
    outcome_counts = (
        seq["cluster_type"]
        .value_counts(sort=False)
        .rename_axis("cluster_type")
        .reset_index(name="n_sequences")
    )
    outcome_counts["prop_sequences"] = (
        outcome_counts["n_sequences"] / outcome_counts["n_sequences"].sum()
    )

    tmp = seq.copy()
    tmp["collection_quarter"] = tmp["collection_date"].dt.to_period("Q").astype(str)
    tmp["who_voc_display"] = tmp["who_voc"].fillna("Other/Non-VOC")

    by_voc = pd.crosstab(tmp["who_voc_display"], tmp["cluster_type"])
    by_voc["total"] = by_voc.sum(axis=1)
    by_voc["large_pct"] = 100 * by_voc["large_cluster"] / by_voc["total"]

    by_quarter = pd.crosstab(tmp["collection_quarter"], tmp["cluster_type"])
    by_quarter["total"] = by_quarter.sum(axis=1)
    by_quarter["large_pct"] = 100 * by_quarter["large_cluster"] / by_quarter["total"]

    test_reason_counts = (
        tmp["test_reason_model"]
        .value_counts(dropna=False)
        .rename_axis("test_reason_model")
        .reset_index(name="n_sequences")
    )
    unmapped_test_reasons = (
        tmp.loc[tmp["test_reason_unmapped"], "test_reason"].dropna().drop_duplicates().sort_values()
    )

    return {
        "outcome_counts": outcome_counts,
        "outcome_by_who_voc": by_voc,
        "outcome_by_quarter": by_quarter,
        "test_reason_counts": test_reason_counts,
        "unmapped_test_reasons": unmapped_test_reasons.to_frame("unmapped_test_reason"),
    }


def prepare_model_inputs(
    raw: pd.DataFrame,
    config: AnalysisConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, FeatureSpec, dict[str, pd.Series | pd.DataFrame]]:
    """Full data preparation from raw rows to sequence table and model inputs."""
    seq = build_sequence_table(raw, config)
    model_df = engineer_features(seq)
    feature_spec = build_feature_spec(config)
    targets = {
        "large": model_df["is_large_cluster"].copy(),
        "singleton": model_df["is_singleton"].copy(),
        "multiclass": model_df["cluster_type"].astype(str).copy(),
        "groups": model_df["cluster_id"].copy(),
    }
    X = model_df[feature_spec.all_features].copy()
    return model_df, X, feature_spec, targets
