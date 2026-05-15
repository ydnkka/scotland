from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import OUT_DIR


def transformed_matrix_for_shap(model, X_raw: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Transform raw feature rows to dense numeric arrays accepted by SHAP."""
    X_transformed = model.named_steps["preprocessor"].transform(X_raw)
    if hasattr(X_transformed, "toarray"):
        X_shap = X_transformed.toarray()
    else:
        X_shap = np.asarray(X_transformed)
    X_shap = np.asarray(X_shap, dtype=np.float64)
    feature_names = model.named_steps["preprocessor"].get_feature_names_out()
    return X_shap, feature_names


def compute_tree_shap_for_positive_class(model, X_raw: pd.DataFrame):
    """Compute SHAP values for the positive class of a binary RF pipeline."""
    import shap

    X_shap, feature_names = transformed_matrix_for_shap(model, X_raw)
    explainer = shap.TreeExplainer(model.named_steps["rf"])
    shap_values = explainer.shap_values(X_shap, check_additivity=False)
    if isinstance(shap_values, list):
        positive_shap = shap_values[1]
    elif getattr(shap_values, "ndim", None) == 3:
        positive_shap = shap_values[:, :, 1]
    else:
        positive_shap = shap_values
    return np.asarray(positive_shap), X_shap, feature_names


def original_feature_name(encoded_name: str, all_features: list[str]) -> str:
    for feature in all_features:
        if encoded_name == feature or encoded_name.startswith(feature + "_"):
            return feature
    if encoded_name.startswith("missingindicator_"):
        return encoded_name.replace("missingindicator_", "")
    return encoded_name


def summarise_shap(
    shap_values: np.ndarray,
    feature_names: np.ndarray,
    all_features: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    encoded = (
        pd.DataFrame(
            {
                "encoded_feature": feature_names,
                "mean_abs_shap": np.abs(shap_values).mean(axis=0),
                "mean_shap": shap_values.mean(axis=0),
                "std_shap": shap_values.std(axis=0),
            }
        )
        .sort_values("mean_abs_shap", ascending=False)
        .reset_index(drop=True)
    )
    encoded["feature"] = encoded["encoded_feature"].map(
        lambda value: original_feature_name(value, all_features)
    )
    grouped = (
        encoded.groupby("feature", as_index=False)
        .agg(
            mean_abs_shap=("mean_abs_shap", "sum"),
            mean_shap=("mean_shap", "sum"),
            n_encoded_features=("encoded_feature", "size"),
        )
        .sort_values("mean_abs_shap", ascending=False)
        .reset_index(drop=True)
    )
    return encoded, grouped


def save_shap_artifacts(
    shap_values: np.ndarray,
    X_shap: np.ndarray,
    X_raw: pd.DataFrame,
    feature_names: np.ndarray,
    all_features: list[str],
    out_dir: str | Path = OUT_DIR,
    prefix: str = "binary_large_cluster",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Save expensive SHAP arrays and compact summaries."""
    shap_dir = Path(out_dir) / "shap"
    shap_dir.mkdir(parents=True, exist_ok=True)

    shap_array = np.asarray(shap_values, dtype=np.float32)
    X_array = np.asarray(X_shap, dtype=np.float32)
    feature_names_array = np.asarray(feature_names, dtype=str)

    np.savez_compressed(
        shap_dir / f"{prefix}_shap_values_encoded.npz",
        large_shap=shap_array,
        X_shap=X_array,
        feature_names=feature_names_array,
        sample_index=X_raw.index.to_numpy(),
    )
    X_raw.to_parquet(shap_dir / f"{prefix}_shap_sample_rows.parquet", index=True)

    encoded, grouped = summarise_shap(shap_array, feature_names_array, all_features)
    encoded.to_csv(shap_dir / f"{prefix}_shap_encoded_feature_summary.csv", index=False)
    grouped.to_csv(shap_dir / f"{prefix}_shap_grouped_feature_summary.csv", index=False)
    return encoded, grouped


def load_shap_artifacts(
    out_dir: str | Path = OUT_DIR,
    prefix: str = "binary_large_cluster",
) -> dict[str, np.ndarray]:
    shap_path = Path(out_dir) / "shap" / f"{prefix}_shap_values_encoded.npz"
    data = np.load(shap_path, allow_pickle=True)
    return {key: data[key] for key in data.files}
