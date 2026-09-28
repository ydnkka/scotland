"""Shared sequence counts, sampling proportions, and lineage summaries."""

from __future__ import annotations

import pandas as pd

from analyses.surveillance.lib.config import SEQUENCE_WINDOW_STRIDE
from utils import CLADE_PALETTE, CLADES, attach_policy_calendar, load_analysis_columns

OTHER_CLADE_LABEL = "Other"


def _ordered_clade_groups(observed_groups: pd.Index) -> list:
    """Order observed clade groups using the shared raw-clade palette order."""
    observed = {group for group in observed_groups if pd.notna(group)}
    clade_order = []

    for raw_clade in CLADE_PALETTE:
        for group in dict.fromkeys([raw_clade, CLADES.get(raw_clade, raw_clade)]):
            if group in observed and group not in clade_order:
                clade_order.append(group)

    remaining = sorted(
        observed - set(clade_order) - {OTHER_CLADE_LABEL},
        key=str,
    )
    clade_order.extend(remaining)

    if OTHER_CLADE_LABEL in observed:
        clade_order.append(OTHER_CLADE_LABEL)

    return clade_order


def build_daily_sequence_counts(
    sequences: pd.DataFrame,
    *,
    smooth_window: int = 7,
) -> pd.DataFrame:
    """Return gap-filled daily counts with a centered rolling mean."""
    daily_counts = (
        sequences.groupby("collection_date")["sequence_id"]
        .nunique()
        .rename("count")
        .reset_index()
        .sort_values("collection_date")
    )

    if daily_counts.empty:
        return pd.DataFrame(columns=["collection_date", "count", "smoothed_count"])

    all_dates = pd.DataFrame(
        {
            "collection_date": pd.date_range(
                daily_counts["collection_date"].min(),
                daily_counts["collection_date"].max(),
                freq="D",
            )
        }
    )

    df_full = all_dates.merge(daily_counts, on="collection_date", how="left")
    df_full["count"] = df_full["count"].fillna(0).astype(int)
    df_full["smoothed_count"] = (
        df_full["count"]
        .astype(float)
        .rolling(window=smooth_window, min_periods=1, center=True)
        .mean()
    )
    return df_full


def attach_policy_timeline(df_full: pd.DataFrame) -> pd.DataFrame:
    """Join policy metadata directly from the processed daily calendar."""
    return attach_policy_calendar(df_full, "collection_date")


def compute_sequencing_proportion(
    df: pd.DataFrame,
    *,
    date_col: str = "collection_date",
    prop_sequenced_col: str = "wn_prop_sequenced",
    time_freq: str = "W",
    smooth_window: int | None = 3,
) -> pd.DataFrame:
    """Return the proportion of positive cases (PCR + antigen tests) sequenced over time."""
    dd = df.dropna(subset=[date_col]).copy()

    if prop_sequenced_col not in dd.columns:
        raise ValueError(f"Missing required column: {prop_sequenced_col}")

    if time_freq == "MS":
        dd["time_period"] = dd[date_col].dt.to_period("M").dt.start_time
    else:
        dd["time_period"] = dd[date_col].dt.to_period(time_freq).dt.start_time

    sampling_df = (
        dd.groupby("time_period")[prop_sequenced_col]
        .mean()
        .reset_index(name="prop_cases_sequenced")
        .sort_values("time_period")
        .set_index("time_period")
    )

    if smooth_window is not None:
        sampling_df["plot_prop_cases_sequenced"] = (
            sampling_df["prop_cases_sequenced"]
            .rolling(window=smooth_window, min_periods=1)
            .mean()
        )
    else:
        sampling_df["plot_prop_cases_sequenced"] = sampling_df["prop_cases_sequenced"]

    return sampling_df


def compute_lineage_frequency_tables(
    df: pd.DataFrame,
    *,
    date_col: str = "collection_date",
    clade_col: str = "variant",
    sequence_col: str = "sequence_id",
    time_freq: str = "W",
    smooth_window: int | None = 3,
    min_sequences_per_period: int = 1,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    list,
]:
    """Return clade frequency, count, and dominance tables for plotting/reporting."""
    dd = df.copy()

    required_cols = {date_col, clade_col, sequence_col}
    missing_cols = required_cols - set(dd.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    dd = dd.dropna(subset=[date_col, clade_col])

    if time_freq == "MS":
        dd["time_period"] = dd[date_col].dt.to_period("M").dt.start_time
    else:
        dd["time_period"] = dd[date_col].dt.to_period(time_freq).dt.start_time

    counts = (
        dd.groupby(["time_period", clade_col])[sequence_col]
        .nunique()
        .reset_index(name="n")
    )

    totals = counts.groupby("time_period")["n"].sum().reset_index(name="total")

    counts = counts.merge(totals, on="time_period")
    counts["frequency"] = counts["n"] / counts["total"]
    counts = counts[counts["total"] >= min_sequences_per_period]

    clade_freq = (
        counts.pivot(index="time_period", columns=clade_col, values="frequency")
        .fillna(0)
        .sort_index()
    )
    clade_order = _ordered_clade_groups(pd.Index(clade_freq.columns))

    clade_freq = clade_freq.reindex(columns=clade_order)

    if clade_freq.empty or not clade_order:
        raise ValueError("No sequences matched the selected clade groups.")

    clade_counts = (
        counts.pivot(index="time_period", columns=clade_col, values="n")
        .fillna(0)
        .astype(int)
        .sort_index()
        .reindex(columns=clade_order, fill_value=0)
    )
    clade_counts.insert(0, "total_sequences", clade_counts.sum(axis=1))

    if smooth_window is not None:
        plot_freq = clade_freq.rolling(window=smooth_window, min_periods=1).mean()
    else:
        plot_freq = clade_freq.copy()

    dominance_df = pd.DataFrame(
        {
            "time_period": plot_freq.index,
            "dominant_clade": plot_freq.idxmax(axis=1).values,
            "dominant_frequency": plot_freq.max(axis=1).values,
        }
    )
    dominance_df["previous_dominant_clade"] = dominance_df["dominant_clade"].shift()

    overtakes = (
        dominance_df[
            dominance_df["dominant_clade"] != dominance_df["previous_dominant_clade"]
        ]
        .dropna(subset=["previous_dominant_clade"])
        .copy()
    )

    return clade_freq, plot_freq, clade_counts, dominance_df, overtakes, clade_order


def load_sequences(*, window_stride: int = SEQUENCE_WINDOW_STRIDE) -> pd.DataFrame:
    if window_stride < 1:
        raise ValueError("window_stride must be at least 1.")
    sequences = load_analysis_columns(
        ["sequence_id", "collection_date", "clade", "wn_prop_sequenced"],
        window_stride=window_stride,
    )
    sequences["variant"] = sequences["clade"].map(CLADES).fillna(OTHER_CLADE_LABEL)
    return sequences
