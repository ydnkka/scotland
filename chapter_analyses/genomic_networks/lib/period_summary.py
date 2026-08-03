"""Combined period-level summaries for Chapter 4 cluster characteristics."""

from __future__ import annotations

import sys
from typing import Any

import numpy as np
import pandas as pd

from .config import PROJECT_ROOT

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import policy_order

MIN_PAIRWISE_ROWS = 10

POLICY_COLUMNS = ["policy_period", "policy_era"]
OUTPUT_COLUMNS = [
    "policy_period",
    "policy_era",
    "n_windows",
    "n_clusters",
    "n_singleton_clusters",
    "n_non_singleton_clusters",
    "median_non_singleton_cluster_size",
    "p90_non_singleton_cluster_size",
    "max_non_singleton_cluster_size",
    "median_non_singleton_duration_days",
    "median_non_singleton_residential_datazones",
    "n_pairwise_windows",
    "n_pairwise_window_lineage_summaries",
    "n_pairwise_rows",
    "within_cluster_snp_distance_weighting",
    "min_pairwise_rows_per_window_lineage",
    "within_cluster_snp_distance_q25",
    "within_cluster_snp_distance_median",
    "within_cluster_snp_distance_q75",
]


def _numeric(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce")


def _weighted_quantile(
    values: pd.Series,
    weights: pd.Series,
    quantile: float,
) -> float:
    work = pd.DataFrame(
        {
            "value": _numeric(values),
            "weight": _numeric(weights),
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


def _sort_by_policy(df: pd.DataFrame) -> pd.DataFrame:
    period_order = {
        period: idx for idx, period in enumerate(policy_order("policy_period"))
    }
    era_order = {era: idx for idx, era in enumerate(policy_order("policy_era"))}
    out = df.copy()
    out["_policy_sort"] = out["policy_period"].astype(str).map(period_order)
    out["_era_sort"] = out["policy_era"].astype(str).map(era_order)
    return (
        out.sort_values(
            ["_policy_sort", "_era_sort", "policy_period", "policy_era"],
            kind="mergesort",
            na_position="last",
        )
        .drop(columns=["_policy_sort", "_era_sort"])
        .reset_index(drop=True)
    )


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
    }
    window_required = {
        "window_id",
        "window_idx",
        "policy_period",
        "policy_era",
    }
    _require_columns(
        pairwise_distance_summary,
        pairwise_required,
        table_name="cluster_pairwise_distance_summary",
    )
    _require_columns(window_coverage, window_required, table_name="window_coverage")
    if min_pairwise_rows < 1:
        raise ValueError("min_pairwise_rows must be at least 1.")

    window_policy = (
        window_coverage[list(window_required)]
        .drop_duplicates(["window_id", "window_idx"])
        .copy()
    )
    work = pairwise_distance_summary.merge(
        window_policy,
        on=["window_id", "window_idx"],
        how="left",
        validate="many_to_one",
    )
    work["n_pairwise_rows"] = _numeric(work["n_pairwise_rows"])
    work["snp_distance_median"] = _numeric(work["snp_distance_median"])
    work = work.replace([np.inf, -np.inf], np.nan)
    work = work.loc[
        work["status"].eq("ok")
        & work["n_pairwise_rows"].ge(min_pairwise_rows)
        & work["n_pairwise_rows"].gt(0)
        & work["snp_distance_median"].notna()
    ].copy()
    work = work.dropna(subset=POLICY_COLUMNS)

    rows: list[dict[str, Any]] = []
    for keys, group in work.groupby(POLICY_COLUMNS, dropna=False, sort=False):
        policy_period, policy_era = keys
        rows.append(
            {
                "policy_period": policy_period,
                "policy_era": policy_era,
                "n_pairwise_windows": group["window_id"].nunique(),
                "n_pairwise_window_lineage_summaries": len(group),
                "n_pairwise_rows": int(group["n_pairwise_rows"].sum()),
                "within_cluster_snp_distance_weighting": "pair_count_weighted",
                "min_pairwise_rows_per_window_lineage": min_pairwise_rows,
                "within_cluster_snp_distance_q25": _weighted_quantile(
                    group["snp_distance_median"], group["n_pairwise_rows"], 0.25
                ),
                "within_cluster_snp_distance_median": _weighted_quantile(
                    group["snp_distance_median"], group["n_pairwise_rows"], 0.50
                ),
                "within_cluster_snp_distance_q75": _weighted_quantile(
                    group["snp_distance_median"], group["n_pairwise_rows"], 0.75
                ),
            }
        )

    distance_columns = [
        "policy_period",
        "policy_era",
        "n_pairwise_windows",
        "n_pairwise_window_lineage_summaries",
        "n_pairwise_rows",
        "within_cluster_snp_distance_weighting",
        "min_pairwise_rows_per_window_lineage",
        "within_cluster_snp_distance_q25",
        "within_cluster_snp_distance_median",
        "within_cluster_snp_distance_q75",
    ]
    return pd.DataFrame(rows, columns=distance_columns)


def build_cluster_period_typical_summary(
    cluster_period_summary: pd.DataFrame,
    pairwise_distance_summary: pd.DataFrame,
    window_coverage: pd.DataFrame,
    *,
    min_pairwise_rows: int = MIN_PAIRWISE_ROWS,
) -> pd.DataFrame:
    """Return one period-level table of non-singleton cluster characteristics.

    Cluster size, duration, and residential spread come from
    ``cluster_period_summary`` and are already restricted to non-singleton
    clusters. Genetic distance is summarised from window-lineage median
    within-cluster SNP distances after joining pairwise rows to policy periods
    through ``window_coverage``.
    """
    cluster_required = {
        "policy_period",
        "policy_era",
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
    cluster = cluster.rename(
        columns={
            "median_non_singleton_datazones": (
                "median_non_singleton_residential_datazones"
            )
        }
    )
    distance = _pairwise_distance_by_period(
        pairwise_distance_summary,
        window_coverage,
        min_pairwise_rows=min_pairwise_rows,
    )
    out = cluster.merge(distance, on=POLICY_COLUMNS, how="left")
    out["within_cluster_snp_distance_weighting"] = out[
        "within_cluster_snp_distance_weighting"
    ].fillna("pair_count_weighted")
    out["min_pairwise_rows_per_window_lineage"] = out[
        "min_pairwise_rows_per_window_lineage"
    ].fillna(min_pairwise_rows)
    for col in (
        "n_pairwise_windows",
        "n_pairwise_window_lineage_summaries",
        "n_pairwise_rows",
    ):
        out[col] = _numeric(out[col]).fillna(0).astype(int)

    for col in OUTPUT_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    return _sort_by_policy(out[OUTPUT_COLUMNS])
