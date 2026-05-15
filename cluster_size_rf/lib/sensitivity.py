from __future__ import annotations

import pandas as pd

from .config import AnalysisConfig
from .data import build_feature_spec, build_sequence_table, engineer_features
from .models import fit_binary_model, make_split_data


def run_threshold_sensitivity(
    raw: pd.DataFrame,
    config: AnalysisConfig,
    large_mins: tuple[int, ...] | None = None,
) -> pd.DataFrame:
    """Refit the primary binary model at alternative large-cluster thresholds."""
    rows = []
    large_mins = config.sensitivity_large_mins if large_mins is None else large_mins
    feature_spec = build_feature_spec(config)
    for large_min in large_mins:
        seq = build_sequence_table(raw, config, large_min=large_min)
        model_df = engineer_features(seq)
        X = model_df[feature_spec.all_features].copy()
        y = model_df["is_large_cluster"].copy()
        groups = model_df["cluster_id"].copy()
        split = make_split_data(X, y, groups, config)
        result = fit_binary_model(
            name=f"large_min_{large_min}",
            X_train=split.X_train,
            y_train=split.y_train,
            X_test=split.X_test,
            y_test=split.y_test,
            config=config,
            numeric_features=feature_spec.numeric_features,
            categorical_features=feature_spec.categorical_features,
            negative_name="not_large",
            positive_name="large",
        )
        rows.append({**result.metrics, "large_min": large_min, "splitter": split.splitter_name})
    return pd.DataFrame(rows)


def run_resolution_sensitivity(
    raw: pd.DataFrame,
    config: AnalysisConfig,
    resolutions: tuple[float, ...] | None = None,
) -> pd.DataFrame:
    """Refit the primary binary model at alternative Leiden resolutions."""
    rows = []
    resolutions = config.sensitivity_resolutions if resolutions is None else resolutions
    feature_spec = build_feature_spec(config)
    for resolution in resolutions:
        seq = build_sequence_table(raw, config, resolution=resolution)
        model_df = engineer_features(seq)
        X = model_df[feature_spec.all_features].copy()
        y = model_df["is_large_cluster"].copy()
        groups = model_df["cluster_id"].copy()
        split = make_split_data(X, y, groups, config)
        result = fit_binary_model(
            name=f"resolution_{resolution}",
            X_train=split.X_train,
            y_train=split.y_train,
            X_test=split.X_test,
            y_test=split.y_test,
            config=config,
            numeric_features=feature_spec.numeric_features,
            categorical_features=feature_spec.categorical_features,
            negative_name="not_large",
            positive_name="large",
        )
        rows.append({**result.metrics, "resolution": resolution, "splitter": split.splitter_name})
    return pd.DataFrame(rows)
