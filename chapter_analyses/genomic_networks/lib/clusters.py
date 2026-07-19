"""Window-level EpiLink cluster summaries for Chapter 4."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from pandas.core.groupby.generic import DataFrameGroupBy

from .config import DEFAULT_MIXING_ATTRIBUTES, DISCLOSURE_MIN_CELL, AttributeSpec


CLUSTER_KEYS = ["cluster_id", "window_id", "window_idx"]


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
    attributes: Iterable[AttributeSpec] = DEFAULT_MIXING_ATTRIBUTES,
) -> pd.DataFrame:
    """Collapse sequence-window rows to one row per EpiLink cluster."""
    missing = [col for col in CLUSTER_KEYS if col not in df.columns]
    if missing:
        raise KeyError(f"Missing cluster key columns: {missing}")

    agg = {
        "wn_start_date": ("wn_start_date", "first"),
        "wn_mid_date": ("wn_mid_date", "first"),
        "wn_end_date": ("wn_end_date", "first"),
        "wn_no_sequences": ("wn_no_sequences", "first"),
        "wn_positive_tests": ("wn_positive_tests", "first"),
        "wn_prop_sequenced": ("wn_prop_sequenced", "first"),
        "cluster_size": ("cluster_size", "first"),
        "cluster_n_datazones": ("cluster_n_datazones", "first"),
        "cluster_duration_days": ("cluster_duration_days", "first"),
        "first_collection_date": ("collection_date", "min"),
        "last_collection_date": ("collection_date", "max"),
        "n_sequences": ("sequence_id", "nunique"),
        "n_datazones_observed": ("datazone", "nunique"),
        "pango_lineage": ("pango_lineage", safe_mode),
        "clade": ("clade", safe_mode),
        "who_voc": ("who_voc", safe_mode),
    }

    optional_mode_cols = [
        "policy_period",
        "policy_era",
        "policy_period_label",
        "test_reason",
        "is_vaccinated",
        "is_reinfection",
    ]
    for col in optional_mode_cols:
        if col in df.columns:
            agg[col] = (col, safe_mode)

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

    date_cols = [
        "wn_start_date",
        "wn_mid_date",
        "wn_end_date",
        "first_collection_date",
        "last_collection_date",
    ]
    for col in date_cols:
        if col in out.columns:
            out[col] = pd.to_datetime(out[col], errors="coerce")

    return out.sort_values(["window_idx", "cluster_size"], ascending=[True, False])


def build_cluster_window_summary(cluster_table: pd.DataFrame) -> pd.DataFrame:
    """Summarise cluster-size and spread distributions by rolling window."""
    required = {"window_id", "window_idx", "cluster_id", "cluster_size"}
    missing = required - set(cluster_table.columns)
    if missing:
        raise KeyError(f"Missing cluster-window summary columns: {sorted(missing)}")

    work = cluster_table.copy()
    cluster_size = pd.to_numeric(work["cluster_size"], errors="coerce")
    work["_is_singleton"] = cluster_size.eq(1)
    work["_is_non_singleton"] = cluster_size.gt(1)
    work["_non_singleton_cluster_size"] = cluster_size.where(
        work["_is_non_singleton"]
    )
    work["_non_singleton_duration_days"] = pd.to_numeric(
        work["cluster_duration_days"], errors="coerce"
    ).where(work["_is_non_singleton"])
    work["_non_singleton_datazones"] = pd.to_numeric(
        work["cluster_n_datazones"], errors="coerce"
    ).where(work["_is_non_singleton"])

    out = (
        work.groupby(["window_id", "window_idx"], dropna=False)
        .agg(
            n_clusters=("cluster_id", "nunique"),
            n_singleton_clusters=("_is_singleton", "sum"),
            n_non_singleton_clusters=("_is_non_singleton", "sum"),
            n_sequence_memberships=("cluster_size", "sum"),
            median_cluster_size=("cluster_size", "median"),
            p90_cluster_size=("cluster_size", lambda x: x.quantile(0.90)),
            max_cluster_size=("cluster_size", "max"),
            median_non_singleton_cluster_size=(
                "_non_singleton_cluster_size",
                "median",
            ),
            p90_non_singleton_cluster_size=(
                "_non_singleton_cluster_size",
                lambda x: x.quantile(0.90),
            ),
            max_non_singleton_cluster_size=(
                "_non_singleton_cluster_size",
                "max",
            ),
            median_non_singleton_duration_days=(
                "_non_singleton_duration_days",
                "median",
            ),
            median_datazones=("cluster_n_datazones", "median"),
            max_datazones=("cluster_n_datazones", "max"),
            median_non_singleton_datazones=(
                "_non_singleton_datazones",
                "median",
            ),
            max_non_singleton_datazones=("_non_singleton_datazones", "max"),
        )
        .reset_index()
        .sort_values("window_idx")
    )
    out["clusters_per_1000_sequences"] = (
        1000 * out["n_clusters"] / out["n_sequence_memberships"].replace(0, np.nan)
    )
    out["non_singleton_clusters_per_1000_sequences"] = (
        1000
        * out["n_non_singleton_clusters"]
        / out["n_sequence_memberships"].replace(0, np.nan)
    )
    return out


def build_cluster_period_summary(cluster_table: pd.DataFrame) -> pd.DataFrame:
    """Summarise cluster distributions by policy period where available."""
    if "policy_period" not in cluster_table.columns and "policy_era" not in cluster_table.columns:
        return pd.DataFrame()

    work = cluster_table.copy()
    cluster_size = pd.to_numeric(work["cluster_size"], errors="coerce")
    work["_is_singleton"] = cluster_size.eq(1)
    work["_is_non_singleton"] = cluster_size.gt(1)
    work["_non_singleton_cluster_size"] = cluster_size.where(
        work["_is_non_singleton"]
    )
    work["_non_singleton_duration_days"] = pd.to_numeric(
        work["cluster_duration_days"], errors="coerce"
    ).where(work["_is_non_singleton"])
    work["_non_singleton_datazones"] = pd.to_numeric(
        work["cluster_n_datazones"], errors="coerce"
    ).where(work["_is_non_singleton"])

    out = (
        work.groupby(["policy_period", "policy_era"], dropna=False)
        .agg(
            n_clusters=("cluster_id", "nunique"),
            n_singleton_clusters=("_is_singleton", "sum"),
            n_non_singleton_clusters=("_is_non_singleton", "sum"),
            median_non_singleton_cluster_size=(
                "_non_singleton_cluster_size",
                "median",
            ),
            p90_non_singleton_cluster_size=(
                "_non_singleton_cluster_size",
                lambda x: x.quantile(0.90),
            ),
            max_non_singleton_cluster_size=(
                "_non_singleton_cluster_size",
                "max",
            ),
            median_non_singleton_duration_days=(
                "_non_singleton_duration_days",
                "median",
            ),
            median_non_singleton_datazones=(
                "_non_singleton_datazones",
                "median",
            ),
            n_windows=("window_id", "nunique"),
        )
        .reset_index()
    )
    return out


def build_cluster_attribute_composition(
    cluster_table: pd.DataFrame,
    *,
    attributes: Iterable[AttributeSpec] = DEFAULT_MIXING_ATTRIBUTES,
    min_cell: int = DISCLOSURE_MIN_CELL,
) -> pd.DataFrame:
    """Count modal cluster attributes for appendix mixing tables."""
    rows: list[pd.DataFrame] = []
    group_cols = [col for col in ("policy_period", "policy_era") if col in cluster_table.columns]

    for spec in attributes:
        modal_col = f"modal_{spec.name}"
        if modal_col not in cluster_table.columns:
            continue

        cols = [*group_cols, modal_col, "cluster_id", "cluster_size"]
        work = cluster_table[cols].copy()
        work[modal_col] = work[modal_col].astype("string").fillna("Missing")
        cluster_size = pd.to_numeric(work["cluster_size"], errors="coerce")
        work["_is_singleton"] = cluster_size.eq(1)
        work["_is_non_singleton"] = cluster_size.gt(1)
        by = [*group_cols, modal_col] if group_cols else [modal_col]
        counts = (
            work.groupby(by, dropna=False)
            .agg(
                n_clusters=("cluster_id", "nunique"),
                n_singleton_clusters=("_is_singleton", "sum"),
                n_non_singleton_clusters=("_is_non_singleton", "sum"),
            )
            .reset_index()
            .rename(columns={modal_col: "category"})
        )
        if group_cols:
            totals = (
                counts.groupby(group_cols, dropna=False)
                .agg(
                    group_clusters=("n_clusters", "sum"),
                    group_singleton_clusters=("n_singleton_clusters", "sum"),
                    group_non_singleton_clusters=(
                        "n_non_singleton_clusters",
                        "sum",
                    ),
                )
                .reset_index()
            )
            counts = counts.merge(totals, on=group_cols, how="left")
        else:
            counts["group_clusters"] = counts["n_clusters"].sum()
            counts["group_singleton_clusters"] = counts[
                "n_singleton_clusters"
            ].sum()
            counts["group_non_singleton_clusters"] = counts[
                "n_non_singleton_clusters"
            ].sum()
        counts.insert(0, "attribute", spec.name)
        counts.insert(1, "attribute_label", spec.label)
        counts["proportion"] = counts["n_clusters"] / counts["group_clusters"].replace(
            0, np.nan
        )
        counts["singleton_proportion"] = counts[
            "n_singleton_clusters"
        ] / counts["group_singleton_clusters"].replace(0, np.nan)
        counts["non_singleton_proportion"] = counts[
            "n_non_singleton_clusters"
        ] / counts["group_non_singleton_clusters"].replace(0, np.nan)
        counts["small_cell"] = counts["n_clusters"].between(1, min_cell - 1)
        counts["singleton_small_cell"] = counts["n_singleton_clusters"].between(
            1, min_cell - 1
        )
        counts["non_singleton_small_cell"] = counts[
            "n_non_singleton_clusters"
        ].between(1, min_cell - 1)
        rows.append(counts)

    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True, sort=False)
