"""Summarise pairwise distances within clusters."""
from __future__ import annotations

import numpy as np
import pandas as pd

MIN_PAIRWISE_ROWS = 10

POLICY_COLUMNS = ["policy_era", "policy_period"]


def numeric(values: pd.Series) -> pd.Series:
    """Return a numeric version of a pandas series."""
    return pd.to_numeric(values, errors="coerce")


def weighted_quantile(
    values: pd.Series,
    weights: pd.Series,
    quantile: float,
) -> float:
    """Return a weighted quantile, dropping missing and non-positive weights."""
    work = pd.DataFrame(
        {
            "value": numeric(values),
            "weight": numeric(weights),
        }
    ).replace([np.inf, -np.inf], np.nan)
    work = work.dropna()
    work = work.loc[work["weight"].gt(0)]
    if work.empty:
        return np.nan

    work = work.sort_values("value", kind="mergesort")
    values_array = work["value"].to_numpy(dtype=float)
    weights_array = work["weight"].to_numpy(dtype=float)
    total_weight = float(weights_array.sum())
    if total_weight <= 0:
        return np.nan

    positions = (np.cumsum(weights_array) - 0.5 * weights_array) / total_weight
    return float(np.interp(quantile, positions, values_array))


def _require_columns(
    frame: pd.DataFrame,
    columns: set[str],
    *,
    table_name: str,
) -> None:
    missing = columns - set(frame.columns)
    if missing:
        raise KeyError(f"Missing {table_name} columns: {sorted(missing)}")


def prepare_pairwise_distances(
    pairwise_distances: pd.DataFrame,
    *,
    min_pairwise_rows: int = MIN_PAIRWISE_ROWS,
    distance_columns: tuple[str, ...] = (
        "snp_distance_median",
        "temporal_distance_median",
    ),
    require_non_empty: bool = False,
) -> pd.DataFrame:
    """Filter pairwise distance summaries to supported within-cluster rows."""
    required = {
        "window_idx",
        "pango_lineage",
        "status",
        "n_pairwise_rows",
        *distance_columns,
    }
    _require_columns(
        pairwise_distances,
        required,
        table_name="cluster_pairwise_distance_summary",
    )
    if min_pairwise_rows < 1:
        raise ValueError("min_pairwise_rows must be at least 1.")

    work = pairwise_distances.copy()
    numeric_columns = ("window_idx", "n_pairwise_rows", *distance_columns)
    for col in numeric_columns:
        work[col] = numeric(work[col])
    work = work.replace([np.inf, -np.inf], np.nan)
    work = work.loc[
        work["status"].eq("ok")
        & work["n_pairwise_rows"].ge(min_pairwise_rows)
        & work["n_pairwise_rows"].gt(0)
    ].copy()
    work = work.dropna(subset=["window_idx", "n_pairwise_rows", *distance_columns])

    if require_non_empty and work.empty:
        raise ValueError(
            "No supported pairwise-distance rows remain after filtering "
            f"status == 'ok' and n_pairwise_rows >= {min_pairwise_rows}."
        )
    return work


def build_cluster_pairwise_distance_overall_summary(
    pairwise_distances: pd.DataFrame,
    *,
    min_pairwise_rows: int = MIN_PAIRWISE_ROWS,
) -> pd.DataFrame:
    """Summarise typical group-level distances with and without pair weights."""
    summary = prepare_pairwise_distances(
        pairwise_distances,
        min_pairwise_rows=min_pairwise_rows,
        require_non_empty=True,
    )
    weights = numeric(summary["n_pairwise_rows"])
    common = {
        "n_window_lineage_summaries": len(summary),
        "n_windows": summary["window_idx"].nunique(),
        "n_pairwise_rows": int(weights.sum()),
    }
    rows = []
    for distance_metric, value_col, unit in (
        ("snp_distance", "snp_distance_median", "SNPs"),
        ("temporal_distance", "temporal_distance_median", "days"),
    ):
        values = numeric(summary[value_col]).replace([np.inf, -np.inf], np.nan)
        valid_values = values.dropna()
        rows.append(
            {
                "distance_metric": distance_metric,
                "unit": unit,
                "weighting": "unweighted",
                "q25": float(valid_values.quantile(0.25)),
                "median": float(valid_values.median()),
                "q75": float(valid_values.quantile(0.75)),
                **common,
            }
        )
        rows.append(
            {
                "distance_metric": distance_metric,
                "unit": unit,
                "weighting": "pair_count_weighted",
                "q25": weighted_quantile(values, weights, 0.25),
                "median": weighted_quantile(values, weights, 0.50),
                "q75": weighted_quantile(values, weights, 0.75),
                **common,
            }
        )
    return pd.DataFrame(rows)
