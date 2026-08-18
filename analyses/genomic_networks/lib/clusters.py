"""Window-level EpiLink cluster summaries"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from typing import Any

import numpy as np
import pandas as pd
from pandas.core.groupby.generic import DataFrameGroupBy

from .config import (
    ASSORTATIVITY_ATTRIBUTES,
    CLUSTER_COLUMNS,
    AttributeSpec,
)
from .policy import canonicalise_policy_fields, sort_by_policy_period

CLUSTER_KEYS = ["cluster_id", "window_id", "window_idx"]
SPATIAL_COORD_COLUMNS = ["dz_xcoord", "dz_ycoord"]
SPATIAL_SUMMARY_COLUMNS = [
    "n_spatial_sequences",
    "n_spatial_pairs",
    "q25_pairwise_residential_distance_km",
    "median_pairwise_residential_distance_km",
    "q75_pairwise_residential_distance_km",
    "max_pairwise_residential_distance_km",
]


def safe_mode(values: pd.Series) -> object:
    """Return the first modal non-missing value, or NaN if none exists."""
    values = values.dropna()
    if values.empty:
        return np.nan
    return values.mode().iloc[0]


def modal_fraction(values: pd.Series) -> float:
    """Return the fraction of non-missing values in the modal category."""
    values = values.dropna()
    if values.empty:
        return np.nan
    counts = values.value_counts(dropna=True)
    return float(counts.iloc[0] / counts.sum())


def _pairwise_residential_distances_km(coords: pd.DataFrame) -> np.ndarray:
    arr = coords[SPATIAL_COORD_COLUMNS].dropna().to_numpy(dtype=float)
    n_coords = len(arr)
    if n_coords < 2:
        return np.array([], dtype=float)

    chunks = []
    for idx in range(n_coords - 1):
        delta = arr[idx + 1 :] - arr[idx]
        chunks.append(np.sqrt((delta * delta).sum(axis=1)) / 1000)
    if not chunks:
        return np.array([], dtype=float)
    return np.concatenate(chunks)


def _spatial_stats_for_group(group: pd.DataFrame) -> dict[str, float | int]:
    coords = group[["sequence_id", *SPATIAL_COORD_COLUMNS]].dropna(
        subset=SPATIAL_COORD_COLUMNS
    )
    coords = coords.drop_duplicates("sequence_id")
    distances = _pairwise_residential_distances_km(coords)
    if distances.size == 0:
        return {
            "n_spatial_sequences": len(coords),
            "n_spatial_pairs": 0,
            "q25_pairwise_residential_distance_km": np.nan,
            "median_pairwise_residential_distance_km": np.nan,
            "q75_pairwise_residential_distance_km": np.nan,
            "max_pairwise_residential_distance_km": np.nan,
        }

    return {
        "n_spatial_sequences": len(coords),
        "n_spatial_pairs": int(distances.size),
        "q25_pairwise_residential_distance_km": float(np.quantile(distances, 0.25)),
        "median_pairwise_residential_distance_km": float(np.median(distances)),
        "q75_pairwise_residential_distance_km": float(np.quantile(distances, 0.75)),
        "max_pairwise_residential_distance_km": float(np.max(distances)),
    }


def build_cluster_spatial_summary(sequence_rows: pd.DataFrame) -> pd.DataFrame:
    """Return residential centroid distance summaries by cluster."""
    required = {*CLUSTER_KEYS, "sequence_id", *SPATIAL_COORD_COLUMNS}
    missing = required - set(sequence_rows.columns)
    if missing:
        raise KeyError(f"Missing cluster spatial columns: {sorted(missing)}")

    work = sequence_rows[[*CLUSTER_KEYS, "sequence_id", *SPATIAL_COORD_COLUMNS]].copy()
    work = work.dropna(subset=SPATIAL_COORD_COLUMNS)
    work = work.drop_duplicates([*CLUSTER_KEYS, "sequence_id"])
    eligible = (
        work.groupby(CLUSTER_KEYS, dropna=False, sort=False)
        .size()
        .rename("n_spatial_sequences")
        .reset_index()
        .loc[lambda x: x["n_spatial_sequences"].ge(2), CLUSTER_KEYS]
    )
    if eligible.empty:
        return pd.DataFrame(columns=[*CLUSTER_KEYS, *SPATIAL_SUMMARY_COLUMNS])
    work = work.merge(eligible, on=CLUSTER_KEYS, how="inner")

    rows = []
    for keys, group in work.groupby(CLUSTER_KEYS, dropna=False, sort=False):
        cluster_id, window_id, window_idx = keys
        rows.append(
            {
                "cluster_id": cluster_id,
                "window_id": window_id,
                "window_idx": window_idx,
                **_spatial_stats_for_group(group),
            }
        )
    return pd.DataFrame(rows, columns=[*CLUSTER_KEYS, *SPATIAL_SUMMARY_COLUMNS])


def attach_cluster_spatial_summary(
    cluster_table: pd.DataFrame,
    sequence_rows: pd.DataFrame,
) -> pd.DataFrame:
    """Attach residential spatial-distance summaries to a cluster table."""
    spatial = build_cluster_spatial_summary(sequence_rows)
    out = cluster_table.drop(columns=SPATIAL_SUMMARY_COLUMNS, errors="ignore")
    return out.merge(spatial, on=CLUSTER_KEYS, how="left", validate="one_to_one")


def _cluster_attribute_stats(
    grouped: DataFrameGroupBy,
    spec: AttributeSpec,
) -> pd.DataFrame:
    """Summarise one categorical attribute within each cluster."""
    stats = grouped[spec.column].agg(
        **{
            f"modal_{spec.name}": safe_mode,
            f"{spec.name}_modal_fraction": modal_fraction,
            f"n_{spec.name}_levels": lambda x: x.dropna().nunique(),
        }
    )
    return stats.reset_index()


def build_cluster_table(
    df: pd.DataFrame,
    *,
    attributes: Iterable[AttributeSpec] = ASSORTATIVITY_ATTRIBUTES,
) -> pd.DataFrame:
    """Collapse sequence-window rows to one row per EpiLink cluster."""
    base = df[list(CLUSTER_COLUMNS)].copy()
    base = base.drop_duplicates(keep="first")
    agg = {
        "clade": ("clade", safe_mode),
        "who_voc": ("who_voc", safe_mode),
        "proportion_vaccinated": ("is_vaccinated", "mean"),
        "proportion_reinfection": ("is_reinfection", "mean"),
        "policy_era": ("policy_era", safe_mode),
        "policy_period": ("policy_period", safe_mode),
        "policy_period_label": ("policy_period_label", safe_mode),
        "test_reason": ("test_reason", safe_mode),
    }

    grouped = df.groupby(CLUSTER_KEYS, dropna=False)
    out = grouped.agg(**agg).reset_index()

    for spec in attributes:
        if spec.column not in df.columns:
            continue
        out = out.merge(
            _cluster_attribute_stats(grouped, spec),
            on=CLUSTER_KEYS,
            how="left",
        )

    out = canonicalise_policy_fields(out)
    out = attach_cluster_spatial_summary(out, df)
    out = base.merge(out, on=CLUSTER_KEYS, how="left", validate="one_to_one")
    return out.sort_values(["window_idx", "cluster_size"], ascending=[True, False])


def _with_singleton_flags(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach reusable singleton/non-singleton flags based on cluster size."""
    work = frame.copy()
    cluster_size = pd.to_numeric(work["cluster_size"], errors="coerce")
    work["_is_singleton"] = cluster_size.eq(1)
    work["_is_non_singleton"] = cluster_size.gt(1)
    return work


def _with_non_singleton_summary_values(cluster_table: pd.DataFrame) -> pd.DataFrame:
    """Attach reusable non-singleton-only summary columns."""
    work = _with_singleton_flags(canonicalise_policy_fields(cluster_table))
    cluster_size = pd.to_numeric(work["cluster_size"], errors="coerce")
    work["_non_singleton_cluster_size"] = cluster_size.where(work["_is_non_singleton"])
    work["_non_singleton_duration_days"] = pd.to_numeric(
        work["cluster_duration_days"], errors="coerce"
    ).where(work["_is_non_singleton"])
    work["_non_singleton_datazones"] = pd.to_numeric(
        work["cluster_n_datazones"], errors="coerce"
    ).where(work["_is_non_singleton"])
    for source_col, work_col in (
        (
            "median_pairwise_residential_distance_km",
            "_non_singleton_pairwise_residential_distance_km",
        ),
        (
            "max_pairwise_residential_distance_km",
            "_non_singleton_spatial_span_km",
        ),
    ):
        if source_col in work.columns:
            values = pd.to_numeric(work[source_col], errors="coerce")
        else:
            values = pd.Series(np.nan, index=work.index)
        work[work_col] = values.where(work["_is_non_singleton"])
    return work


AggSpec = dict[str, tuple[str, str | Callable[..., Any]]]


_BASE_CLUSTER_SUMMARY_AGG: AggSpec = {
    "n_sequence_memberships": ("cluster_size", "sum"),
    "n_clusters": ("cluster_id", "nunique"),
    "n_singleton_clusters": ("_is_singleton", "sum"),
    "n_non_singleton_clusters": ("_is_non_singleton", "sum"),
    "median_non_singleton_cluster_size": (
        "_non_singleton_cluster_size",
        "median",
    ),
    "p90_non_singleton_cluster_size": (
        "_non_singleton_cluster_size",
        lambda x: x.quantile(0.90),
    ),
    "max_non_singleton_cluster_size": (
        "_non_singleton_cluster_size",
        "max",
    ),
    "q25_non_singleton_duration_days": (
        "_non_singleton_duration_days",
        lambda x: x.quantile(0.25),
    ),
    "median_non_singleton_duration_days": (
        "_non_singleton_duration_days",
        "median",
    ),
    "q75_non_singleton_duration_days": (
        "_non_singleton_duration_days",
        lambda x: x.quantile(0.75),
    ),
    "median_non_singleton_datazones": (
        "_non_singleton_datazones",
        "median",
    ),
    "p90_non_singleton_datazones": (
        "_non_singleton_datazones",
        lambda x: x.quantile(0.90),
    ),
    "max_non_singleton_datazones": (
        "_non_singleton_datazones",
        "max",
    ),
    "q25_non_singleton_spatial_distance_km": (
        "_non_singleton_pairwise_residential_distance_km",
        lambda x: x.quantile(0.25),
    ),
    "median_non_singleton_spatial_distance_km": (
        "_non_singleton_pairwise_residential_distance_km",
        "median",
    ),
    "q75_non_singleton_spatial_distance_km": (
        "_non_singleton_pairwise_residential_distance_km",
        lambda x: x.quantile(0.75),
    ),
}


def _aggregate_cluster_summary(
    work: pd.DataFrame,
    group_columns: Sequence[str],
    extra_aggregations: AggSpec | None = None,
) -> pd.DataFrame:
    """Apply the common cluster-summary aggregation backbone."""
    aggregations = {
        **_BASE_CLUSTER_SUMMARY_AGG,
        **(extra_aggregations or {}),
    }

    return (
        work.groupby(list(group_columns), dropna=False)
        .agg(**aggregations)
        .reset_index()
    )


def build_cluster_window_summary(cluster_table: pd.DataFrame) -> pd.DataFrame:
    """Summarise cluster-size and spread distributions by rolling window."""
    work = _with_non_singleton_summary_values(cluster_table)

    out = _aggregate_cluster_summary(
        work,
        group_columns=("window_id", "window_idx", "wn_no_sequences", "wn_positive_tests", "wn_prop_sequenced"),
    ).sort_values("window_idx")

    sequence_memberships = out["n_sequence_memberships"].replace(0, np.nan)

    out["clusters_per_1000_sequences"] = 1000 * out["n_clusters"] / sequence_memberships
    out["non_singleton_clusters_per_1000_sequences"] = (
        1000 * out["n_non_singleton_clusters"] / sequence_memberships
    )

    return out


def build_cluster_period_summary(cluster_table: pd.DataFrame) -> pd.DataFrame:
    """Summarise cluster distributions by policy period where available."""

    work = _with_non_singleton_summary_values(cluster_table)

    out = _aggregate_cluster_summary(
        work,
        group_columns=("policy_era", "policy_period"),
        extra_aggregations={
            "n_windows": ("window_id", "nunique"),
            "median_window_prop_sequenced": ("wn_prop_sequenced", "median")
        },
    )


    return sort_by_policy_period(out)
