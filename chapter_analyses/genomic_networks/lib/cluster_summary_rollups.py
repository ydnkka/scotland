"""Reusable cluster-summary helpers for Chapter 4."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .policy import sort_by_policy_period

MIN_PAIRWISE_ROWS = 10

POLICY_COLUMNS = ["policy_era", "policy_period"]
PERIOD_TYPICAL_OUTPUT_COLUMNS = [
    *POLICY_COLUMNS,
    "n_windows",
    "n_clusters",
    "n_singleton_clusters",
    "n_non_singleton_clusters",
    "median_non_singleton_cluster_size",
    "p90_non_singleton_cluster_size",
    "max_non_singleton_cluster_size",
    "q25_non_singleton_duration_days",
    "median_non_singleton_duration_days",
    "q75_non_singleton_duration_days",
    "median_non_singleton_datazones",
    "p90_non_singleton_datazones",
    "max_non_singleton_datazones",
    "q25_non_singleton_spatial_distance_km",
    "median_non_singleton_spatial_distance_km",
    "q75_non_singleton_spatial_distance_km",
    "n_pairwise_windows",
    "n_pairwise_window_lineage_summaries",
    "n_pairwise_rows",
    "within_cluster_distance_weighting",
    "min_pairwise_rows_per_window_lineage",
    "within_cluster_snp_distance_q25",
    "within_cluster_snp_distance_median",
    "within_cluster_snp_distance_q75",
    "within_cluster_temporal_distance_q25_days",
    "within_cluster_temporal_distance_median_days",
    "within_cluster_temporal_distance_q75_days",
]


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


def prepare_pairwise_distance_summary(
    pairwise_distance_summary: pd.DataFrame,
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
        pairwise_distance_summary,
        required,
        table_name="cluster_pairwise_distance_summary",
    )
    if min_pairwise_rows < 1:
        raise ValueError("min_pairwise_rows must be at least 1.")

    work = pairwise_distance_summary.copy()
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
    pairwise_distance_summary: pd.DataFrame,
    *,
    min_pairwise_rows: int = MIN_PAIRWISE_ROWS,
) -> pd.DataFrame:
    """Summarise typical group-level distances with and without pair weights."""
    summary = prepare_pairwise_distance_summary(
        pairwise_distance_summary,
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


def _pairwise_distance_by_period(
    pairwise_distance_summary: pd.DataFrame,
    window_coverage: pd.DataFrame,
    *,
    min_pairwise_rows: int,
) -> pd.DataFrame:
    pairwise_required = {
        "window_id",
        "window_idx",
        "pango_lineage",
        "status",
        "n_pairwise_rows",
        "snp_distance_median",
        "temporal_distance_median",
    }
    window_required = [
        "window_id",
        "window_idx",
        "policy_period",
        "policy_era",
    ]
    _require_columns(
        pairwise_distance_summary,
        pairwise_required,
        table_name="cluster_pairwise_distance_summary",
    )
    _require_columns(
        window_coverage,
        set(window_required),
        table_name="window_coverage",
    )

    window_policy = (
        window_coverage[window_required]
        .drop_duplicates(["window_id", "window_idx"])
        .copy()
    )
    work = prepare_pairwise_distance_summary(
        pairwise_distance_summary,
        min_pairwise_rows=min_pairwise_rows,
        distance_columns=("snp_distance_median", "temporal_distance_median"),
    )
    work = work.merge(
        window_policy,
        on=["window_id", "window_idx"],
        how="left",
        validate="many_to_one",
    ).dropna(subset=POLICY_COLUMNS)

    rows: list[dict[str, Any]] = []
    for keys, group in work.groupby(POLICY_COLUMNS, dropna=False, sort=False):
        policy_period, policy_era = keys
        rows.append(
            {
                "policy_era": policy_era,
                "policy_period": policy_period,
                "n_pairwise_windows": group["window_id"].nunique(),
                "n_pairwise_window_lineage_summaries": len(group),
                "n_pairwise_rows": int(group["n_pairwise_rows"].sum()),
                "within_cluster_distance_weighting": "pair_count_weighted",
                "min_pairwise_rows_per_window_lineage": min_pairwise_rows,
                "within_cluster_snp_distance_q25": weighted_quantile(
                    group["snp_distance_median"], group["n_pairwise_rows"], 0.25
                ),
                "within_cluster_snp_distance_median": weighted_quantile(
                    group["snp_distance_median"], group["n_pairwise_rows"], 0.50
                ),
                "within_cluster_snp_distance_q75": weighted_quantile(
                    group["snp_distance_median"], group["n_pairwise_rows"], 0.75
                ),
                "within_cluster_temporal_distance_q25_days": weighted_quantile(
                    group["temporal_distance_median"],
                    group["n_pairwise_rows"],
                    0.25,
                ),
                "within_cluster_temporal_distance_median_days": weighted_quantile(
                    group["temporal_distance_median"],
                    group["n_pairwise_rows"],
                    0.50,
                ),
                "within_cluster_temporal_distance_q75_days": weighted_quantile(
                    group["temporal_distance_median"],
                    group["n_pairwise_rows"],
                    0.75,
                ),
            }
        )

    distance_columns = [
        "policy_era",
        "policy_period",
        "n_pairwise_windows",
        "n_pairwise_window_lineage_summaries",
        "n_pairwise_rows",
        "within_cluster_distance_weighting",
        "min_pairwise_rows_per_window_lineage",
        "within_cluster_snp_distance_q25",
        "within_cluster_snp_distance_median",
        "within_cluster_snp_distance_q75",
        "within_cluster_temporal_distance_q25_days",
        "within_cluster_temporal_distance_median_days",
        "within_cluster_temporal_distance_q75_days",
    ]
    return pd.DataFrame(rows, columns=distance_columns)


def build_cluster_period_typical_summary(
    cluster_period_summary: pd.DataFrame,
    pairwise_distance_summary: pd.DataFrame,
    window_coverage: pd.DataFrame,
    *,
    min_pairwise_rows: int = MIN_PAIRWISE_ROWS,
) -> pd.DataFrame:
    """Return one period-level table of non-singleton cluster characteristics."""
    cluster_required = {
        "policy_era",
        "policy_period",
        "n_windows",
        "n_clusters",
        "n_singleton_clusters",
        "n_non_singleton_clusters",
        "median_non_singleton_cluster_size",
        "p90_non_singleton_cluster_size",
        "max_non_singleton_cluster_size",
        "median_non_singleton_duration_days",
        "median_non_singleton_datazones",
    }
    _require_columns(
        cluster_period_summary,
        cluster_required,
        table_name="cluster_period_summary",
    )

    cluster = cluster_period_summary.copy()
    for col in (
        "q25_non_singleton_spatial_distance_km",
        "median_non_singleton_spatial_distance_km",
        "q75_non_singleton_spatial_distance_km",
    ):
        if col not in cluster.columns:
            cluster[col] = np.nan

    distance = _pairwise_distance_by_period(
        pairwise_distance_summary,
        window_coverage,
        min_pairwise_rows=min_pairwise_rows,
    )
    out = cluster.merge(distance, on=POLICY_COLUMNS, how="left")
    out["within_cluster_distance_weighting"] = out[
        "within_cluster_distance_weighting"
    ].fillna("pair_count_weighted")
    out["min_pairwise_rows_per_window_lineage"] = out[
        "min_pairwise_rows_per_window_lineage"
    ].fillna(min_pairwise_rows)
    for col in (
        "n_pairwise_windows",
        "n_pairwise_window_lineage_summaries",
        "n_pairwise_rows",
    ):
        out[col] = numeric(out[col]).fillna(0).astype(int)

    for col in PERIOD_TYPICAL_OUTPUT_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    return sort_by_policy_period(out[PERIOD_TYPICAL_OUTPUT_COLUMNS])
