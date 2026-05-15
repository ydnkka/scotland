from __future__ import annotations

import json
from pathlib import Path

import joblib
import pandas as pd

from .config import AnalysisConfig, OUT_DIR
from .data import FeatureSpec
from .models import BinaryResult, MulticlassResult


def ensure_output_dir(out_dir: str | Path = OUT_DIR) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_outcome_tables(tables: dict[str, pd.DataFrame], out_dir: str | Path = OUT_DIR) -> None:
    out = ensure_output_dir(out_dir)
    mapping = {
        "outcome_counts": "primary_outcome_counts.csv",
        "outcome_by_who_voc": "primary_outcome_by_who_voc.csv",
        "outcome_by_quarter": "primary_outcome_by_quarter.csv",
        "test_reason_counts": "test_reason_category_counts.csv",
        "unmapped_test_reasons": "test_reason_unmapped_values.csv",
    }
    for key, filename in mapping.items():
        if key in tables:
            keep_index = key in {"outcome_by_who_voc", "outcome_by_quarter"}
            tables[key].to_csv(out / filename, index=keep_index)


def save_binary_result(
    result: BinaryResult,
    out_dir: str | Path = OUT_DIR,
    save_model: bool = True,
) -> None:
    out = ensure_output_dir(out_dir)
    prefix = result.name
    pd.DataFrame([result.metrics]).to_csv(out / f"{prefix}_metrics.csv", index=False)
    result.confusion_matrix.to_csv(out / f"{prefix}_confusion_matrix.csv")
    result.classification_report.to_csv(out / f"{prefix}_classification_report.csv")
    if result.importance is not None:
        result.importance.to_csv(out / f"{prefix}_permutation_importance.csv", index=False)
    if save_model:
        joblib.dump(result.model, out / f"{prefix}_random_forest.joblib")


def save_multiclass_result(
    result: MulticlassResult,
    out_dir: str | Path = OUT_DIR,
    save_model: bool = True,
) -> None:
    out = ensure_output_dir(out_dir)
    prefix = result.name
    pd.DataFrame([result.metrics]).to_csv(out / f"{prefix}_metrics.csv", index=False)
    result.confusion_matrix.to_csv(out / f"{prefix}_confusion_matrix.csv")
    result.classification_report.to_csv(out / f"{prefix}_classification_report.csv")
    if save_model:
        joblib.dump(result.model, out / f"{prefix}_random_forest.joblib")


def save_run_metadata(
    config: AnalysisConfig,
    feature_spec: FeatureSpec,
    out_dir: str | Path = OUT_DIR,
    extra: dict | None = None,
) -> None:
    out = ensure_output_dir(out_dir)
    metadata = {
        **config.to_dict(),
        "simd_overall_feature": feature_spec.simd_overall_feature,
        "simd_domain_mode": feature_spec.simd_domain_mode,
        "simd_band_weighting": feature_spec.simd_band_weighting,
        "simd_domain_features": feature_spec.simd_domain_features,
        "numeric_features": feature_spec.numeric_features,
        "categorical_features": feature_spec.categorical_features,
        "decomposition_numeric_features": feature_spec.decomposition_numeric_features,
        "decomposition_features": feature_spec.decomposition_features,
    }
    if extra:
        metadata.update(extra)
    with open(out / "run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)


def save_sensitivity_metrics(
    threshold_results: pd.DataFrame | None = None,
    resolution_results: pd.DataFrame | None = None,
    out_dir: str | Path = OUT_DIR,
) -> None:
    out = ensure_output_dir(out_dir)
    if threshold_results is not None and len(threshold_results):
        threshold_results.to_csv(out / "threshold_sensitivity_metrics.csv", index=False)
    if resolution_results is not None and len(resolution_results):
        resolution_results.to_csv(out / "resolution_sensitivity_metrics.csv", index=False)
