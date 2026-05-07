"""Plot daily sequenced cases with policy periods and lineage-group frequency."""

from __future__ import annotations

from typing import Callable

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import pandas as pd
import polars as pl

from utils import data, policy, style


WAVE_GROUPS: dict[str, Callable[[str], bool]] = {
    "B.1.177": lambda lineage: lineage.startswith("B.1.177"),
    "Alpha": lambda lineage: lineage == "B.1.1.7" or lineage.startswith("B.1.1.7."),
    "Delta": lambda lineage: lineage.startswith("AY.") or lineage == "B.1.617.2",
    "BA.1": lambda lineage: lineage.startswith("BA.1"),
    "BA.2": lambda lineage: lineage.startswith("BA.2"),
    "BA.4": lambda lineage: lineage.startswith("BA.4"),
    "BA.5": lambda lineage: lineage.startswith("BA.5") or lineage.startswith("BE."),
    "BQ.1": lambda lineage: lineage.startswith("BQ."),
    "XBB": lambda lineage: lineage.startswith("XBB"),
}

WAVE_GROUP_PALETTE: dict[str, str] = {
    "B.1.177": "#4e79a7",
    "Alpha": "#f28e2b",
    "Delta": "#e15759",
    "BA.1": "#76b7b2",
    "BA.2": "#59a14f",
    "BA.4": "#edc948",
    "BA.5": "#b07aa1",
    "BQ.1": "#ff9da7",
    "XBB": "#9c755f",
}


def assign_wave_group(lineage: str | None) -> str | None:
    """Map a Pango lineage to one of the selected wave groups."""
    if lineage is None:
        return None

    lineage = str(lineage).strip()
    if not lineage:
        return None

    for group_name, matcher in WAVE_GROUPS.items():
        if matcher(lineage):
            return group_name

    return None


def build_daily_sequence_counts(
    sequences: pl.DataFrame,
    *,
    smooth_window: int = 7,
) -> pl.DataFrame:
    """Return gap-filled daily counts with a centered rolling mean."""
    sequences = sequences.with_columns(pl.col("collection_date").cast(pl.Date))

    daily_counts = (
        sequences
        .group_by("collection_date")
        .agg(pl.col("sequence_id").n_unique().alias("count"))
        .sort("collection_date")
    )

    all_dates = pl.DataFrame({
        "collection_date": pl.date_range(
            daily_counts["collection_date"].min(),
            daily_counts["collection_date"].max(),
            interval="1d",
            eager=True,
        )
    })

    df_full = (
        all_dates
        .join(daily_counts, on="collection_date", how="left")
        .with_columns(pl.col("count").fill_null(0))
    )

    half_window = smooth_window // 2
    return (
        df_full
        .with_columns(
            pl.col("count")
            .cast(pl.Float64)
            .rolling_mean(window_size=smooth_window, min_samples=1)
            .shift(-half_window)
            .alias("smoothed_count")
        )
        .with_columns(pl.col("smoothed_count").forward_fill().backward_fill())
    )


def attach_policy_timeline(df_full: pl.DataFrame) -> pl.DataFrame:
    """Attach policy metadata for each daily date."""
    period_lookup = (
        policy.POLICY_PERIODS
        .select(pl.col("period_code").alias("policy_period"))
        .with_row_index("policy_level")
    )

    return (
        policy.attach_period(df_full, "collection_date")
        .join(period_lookup, on="policy_period", how="left")
    )


def configure_date_axis(ax: plt.Axes, dates: object) -> None:
    """Use denser half-year ticks for shorter date ranges."""
    date_index = pd.to_datetime(pd.Index(dates)).dropna()
    if date_index.empty:
        return

    start = date_index.min()
    end = date_index.max()
    span_days = (end - start).days

    ax.set_xlim(
        (start - pd.Timedelta(days=7)).to_pydatetime(),
        (end + pd.Timedelta(days=7)).to_pydatetime(),
    )

    if span_days <= 365 * 4:
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
        ax.tick_params(axis="x")
        for label in ax.get_xticklabels():
            label.set_ha("center")
    else:
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))


def add_policy_background(
    ax: plt.Axes,
    dates: object,
    *,
    show_labels: bool = False,
) -> None:
    """Shade alternating policy periods behind the plotted data."""
    date_index = pd.to_datetime(pd.Index(dates)).dropna()
    if date_index.empty:
        return

    plot_start = date_index.min().normalize()
    plot_end = date_index.max().normalize()
    policy_periods = policy.POLICY_PERIODS.sort("start_date").to_pandas()

    for idx, row in policy_periods.iterrows():
        start = max(pd.Timestamp(row["start_date"]), plot_start)
        end = min(pd.Timestamp(row["end_date"]), plot_end)
        if start > end:
            continue

        if idx % 2 == 0:
            ax.axvspan(
                start,
                end + pd.Timedelta(days=1),
                color="#000000",
                alpha=0.035,
                lw=0,
                zorder=0,
            )

        ax.axvline(
            start,
            color="#b3b3b3",
            lw=0.45,
            alpha=0.35,
            zorder=1,
        )

        if show_labels:
            midpoint = start + (end - start) / 2
            ax.text(
                midpoint,
                1.01,
                str(row["period_code"]),
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=6.5,
                color="#666666",
                clip_on=False,
                zorder=4,
            )


def plot_sequences_with_policy(
    timeline: pl.DataFrame,
    *,
    ax: plt.Axes,
    show_xlabel: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot daily sequences with subtle policy-period shading."""
    fig = ax.figure

    dates = timeline["collection_date"].to_list()
    counts = timeline["count"].to_list()
    smoothed = timeline["smoothed_count"].to_list()
    add_policy_background(ax, dates, show_labels=True)

    ax.bar(
        dates,
        counts,
        width=1.0,
        color="#d9e2ec",
        edgecolor="none",
        alpha=0.95,
        label="Daily sequences",
        zorder=2,
    )
    ax.plot(
        dates,
        smoothed,
        color="#0b1f3b",
        lw=1.8,
        label="7-day smoothed count",
        zorder=3,
    )
    ax.set_ylabel("Number of sequences")
    if show_xlabel:
        ax.set_xlabel("Collection date")
    ax.set_facecolor("white")
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.5)
    ax.margins(x=0.01)

    configure_date_axis(ax, dates)
    ax.legend(loc="upper left", ncol=2, frameon=False)

    return fig, ax

def plot_lineage_frequency_and_overtakes(
    df: pl.DataFrame,
    *,
    date_col: str = "collection_date",
    clade_col: str = "wave_group",
    sequence_col: str = "sequence_id",
    prop_sequenced_col: str = "wn_prop_sequenced",
    time_freq: str = "W",
    smooth_window: int | None = 3,
    min_sequences_per_period: int = 1,
    ax: plt.Axes,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, plt.Axes]:
    """Plot selected lineage-group frequency and sequencing coverage over time.

    Returns
    -------
    lineage_freq
        Unsmoothed lineage-group frequency table.
    plot_freq
        Smoothed lineage-group frequency table used for plotting.
    dominance_df
        Dominant lineage group per period.
    overtakes
        Periods where the dominant lineage group changed.
    sampling_df
        Proportion of cases sequenced per period.
    ax2
        Twin y-axis for proportion of cases sequenced.
    """
    dd = df.to_pandas()

    required_cols = {date_col, clade_col, sequence_col, prop_sequenced_col}
    missing_cols = required_cols - set(dd.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    dd[date_col] = pd.to_datetime(dd[date_col], errors="coerce")
    dd[prop_sequenced_col] = pd.to_numeric(dd[prop_sequenced_col], errors="coerce")
    dd = dd.dropna(subset=[date_col, clade_col])

    if time_freq == "W":
        dd["time_period"] = dd[date_col].dt.to_period("W").dt.start_time
    elif time_freq == "MS":
        dd["time_period"] = dd[date_col].dt.to_period("M").dt.start_time
    else:
        dd["time_period"] = dd[date_col].dt.to_period(time_freq).dt.start_time

    counts = (
        dd.groupby(["time_period", clade_col])[sequence_col]
        .nunique()
        .reset_index(name="n")
    )

    totals = (
        counts.groupby("time_period")["n"]
        .sum()
        .reset_index(name="total")
    )

    counts = counts.merge(totals, on="time_period")
    counts["frequency"] = counts["n"] / counts["total"]
    counts = counts[counts["total"] >= min_sequences_per_period]

    lineage_freq = (
        counts.pivot(index="time_period", columns=clade_col, values="frequency")
        .fillna(0)
        .sort_index()
    )
    lineage_order = [group for group in WAVE_GROUPS if group in lineage_freq.columns]
    lineage_freq = lineage_freq.reindex(columns=lineage_order)

    if lineage_freq.empty or not lineage_order:
        raise ValueError("No sequences matched the selected lineage groups.")

    sampling_df = (
        dd.groupby("time_period")[prop_sequenced_col]
        .mean()
        .reset_index(name="prop_cases_sequenced")
        .sort_values("time_period")
        .set_index("time_period")
    )

    # Keep sampling data aligned to periods that pass the sequence-count threshold.
    sampling_df = sampling_df.reindex(lineage_freq.index)

    if smooth_window is not None:
        plot_freq = lineage_freq.rolling(window=smooth_window, min_periods=1).mean()
        plot_sampling = (
            sampling_df["prop_cases_sequenced"]
            .rolling(window=smooth_window, min_periods=1)
            .mean()
        )
    else:
        plot_freq = lineage_freq.copy()
        plot_sampling = sampling_df["prop_cases_sequenced"].copy()

    sampling_df = sampling_df.copy()
    sampling_df["plot_prop_cases_sequenced"] = plot_sampling

    dominance_df = pd.DataFrame({
        "time_period": plot_freq.index,
        "dominant_lineage_group": plot_freq.idxmax(axis=1).values,
        "dominant_frequency": plot_freq.max(axis=1).values,
    })
    dominance_df["previous_dominant_lineage_group"] = (
        dominance_df["dominant_lineage_group"].shift()
    )

    overtakes = dominance_df[
        dominance_df["dominant_lineage_group"]
        != dominance_df["previous_dominant_lineage_group"]
    ].dropna(subset=["previous_dominant_lineage_group"]).copy()

    ax.set_facecolor("white")
    add_policy_background(ax, plot_freq.index)
    stack_colors = [WAVE_GROUP_PALETTE.get(group, "#999999") for group in lineage_order]
    stack_handles = ax.stackplot(
        plot_freq.index,
        *[plot_freq[group].values for group in lineage_order],
        labels=lineage_order,
        colors=stack_colors,
        alpha=0.92,
        linewidth=0.5,
        edgecolor="white",
        zorder=2,
    )

    ax2 = ax.twinx()
    ax2.plot(
        plot_sampling.index,
        plot_sampling.values,
        linestyle=":",
        linewidth=1.6,
        alpha=0.85,
        color="#303030",
        label="Cases sequenced",
    )
    ax2.set_ylabel("Proportion of cases sequenced")
    ax2.set_ylim(0, 1)
    ax2.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax2.tick_params(axis="y", colors="#4d4d4d")
    ax2.grid(False)

    ax.set_xlabel("Collection date")
    ax.set_ylabel("Lineage-group frequency")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.5)
    ax.margins(x=0.01)

    configure_date_axis(ax, dd[date_col])

    coverage_handle = Line2D(
        [0], [0],
        color="#303030",
        linestyle=":",
        linewidth=1.6,
        label="Cases sequenced",
    )
    ax.legend(
        [*stack_handles, coverage_handle],
        [*lineage_order, "Cases sequenced"],
        ncol=5,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        frameon=False,
        columnspacing=1.1,
        handlelength=1.5,
    )

    return lineage_freq, plot_freq, dominance_df, overtakes, sampling_df, ax2


def main() -> None:
    style.set_theme(context="paper")
    paths = data.Paths.from_config()
    out_dir = paths.root / "figures"
    table_dir = paths.root / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    qc_statuses: list[data.QCStatus] = ["good", "mediocre", "bad"]

    sequences = data.load_analysis_columns(
        ["sequence_id", "collection_date", "pango_lineage", "wn_prop_sequenced"],
        resolution=data.PRIMARY_RESOLUTION,
        qc=qc_statuses,
    )

    sequences = (
        sequences.group_by("sequence_id")
        .agg([
            pl.first("collection_date").alias("collection_date"),
            pl.first("pango_lineage").alias("pango_lineage"),
            pl.mean("wn_prop_sequenced").alias("wn_prop_sequenced"),
        ])
        .with_columns(
            pl.col("pango_lineage")
            .map_elements(assign_wave_group, return_dtype=pl.Utf8)
            .alias("wave_group")
        )
    )

    timeline = attach_policy_timeline(
        build_daily_sequence_counts(sequences, smooth_window=7)
    )

    fig, axes = style.new_figure(
        width="double",
        height_in=6.4,
        nrows=2, ncols=1,
    )
    fig.subplots_adjust(hspace=0.08)

    axes = list(axes)

    ax_top =axes[0]
    ax_bottom = axes[1]

    fig, _ = plot_sequences_with_policy(timeline, ax=ax_top, show_xlabel=False)
    ax_top.tick_params(axis="x", labelbottom=False)

    lineage_freq, plot_freq, dominance_df, overtakes, sampling_df, _ = (
        plot_lineage_frequency_and_overtakes(sequences, ax=ax_bottom)
    )

    lineage_freq.to_csv(table_dir / "lineage_frequency_by_period.csv")
    plot_freq.to_csv(table_dir / "lineage_frequency_by_period_smoothed.csv")
    dominance_df.to_csv(table_dir / "lineage_dominance_by_period.csv", index=False)
    overtakes.to_csv(table_dir / "lineage_overtake_events.csv", index=False)
    sampling_df.to_csv(table_dir / "sequencing_proportion_by_period.csv")

    fig.align_ylabels([ax_top, ax_bottom])
    style.add_panel_labels([ax_top, ax_bottom], x=-0.1, y=1)

    _ = style.save_figure(
        fig,
        out_dir / "policy_sequences_over_time",
        width="double",
        save_png=True,
        save_pdf=True,
    )


if __name__ == "__main__":
    main()
