"""Plot daily sequenced cases with policy periods and lineage-group frequency."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import colors
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import PercentFormatter

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from chapter_analyses.surveillance.lib.config import (  # noqa: E402
    DAILY_SMOOTH_WINDOW,
    FIGURE_NAME,
    FIGURES_DIR,
    SEQUENCE_WINDOW_STRIDE,
    TABLES_DIR,
)
from chapter_analyses.surveillance.lib.io import write_table  # noqa: E402

from utils import (  # noqa: E402
    CLADES,
    CLADE_PALETTE,
    attach_policy_calendar,
    set_theme,
    load_analysis_columns,
    load_policy_calendar,
    add_panel_labels,
    new_figure,
    save_figure,
    policy_era_labels,
)

POLICY_STRINGENCY_CMAP = plt.get_cmap("RdYlGn_r")
POLICY_STRINGENCY_NORM = colors.Normalize(
    vmin=1,
    vmax=100,
)
LOGGER = logging.getLogger(__name__)


POLICY_ERA_LABELS = policy_era_labels()
POLICY_CALENDAR = load_policy_calendar()
POLICY_PERIODS = (
    POLICY_CALENDAR[
        [
            "policy_period",
            "policy_period_label",
            "policy_period_start_date",
            "policy_period_end_date",
            "policy_period_order",
            "policy_era",
        ]
    ]
    .drop_duplicates()
    .merge(
        POLICY_CALENDAR.groupby("policy_period", sort=False, observed=True)[
            ["stringency_index", "containment_index"]
        ]
        .mean()
        .rename(
            columns={
                "stringency_index": "policy_stringency",
                "containment_index": "policy_containment",
            }
        )
        .reset_index(),
        on="policy_period",
        how="left",
        validate="one_to_one",
    )
    .rename(
        columns={
            "policy_period": "period_code",
            "policy_period_label": "period_label",
            "policy_period_start_date": "start_date",
            "policy_period_end_date": "end_date",
            "policy_period_order": "period_order",
        }
    )
    .sort_values("period_order", ignore_index=True)
)

# Coarse epidemic eras derived from the individual policy periods.
POLICY_ERAS = (
    POLICY_PERIODS.groupby("policy_era", sort=False)
    .agg(
        start_date=("start_date", "min"),
        end_date=("end_date", "max"),
        era_order=("period_order", "min"),
    )
    .reset_index()
    .sort_values("era_order", ignore_index=True)
)


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


def configure_date_axis(ax: Axes, dates: pd.Series):
    """Use quarterly date ticks with month above year."""
    dates = dates.dropna()
    if dates.empty:
        return

    start = dates.min()
    end = dates.max()
    span_days = (end - start).days
    ax.set_xlim(
        (start - pd.Timedelta(days=7)),
        (end + pd.Timedelta(days=7)),
    )

    if span_days <= 365 * 4:
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
        ax.tick_params(axis="x")
        for label in ax.get_xticklabels():
            label.set_horizontalalignment("center")
    else:
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))


def add_epidemic_eras(
    ax: Axes,
    dates: pd.Series,
    *,
    label: bool = True,
    label_y: float = 0.4,
    zorder: int = 1,
) -> None:
    """Mark coarse epidemic-era boundaries with vertical lines (and labels)."""
    dates = dates.dropna()
    if dates.empty:
        return

    plot_start = dates.min().normalize()
    plot_end = dates.max().normalize()

    eras = POLICY_ERAS.sort_values("start_date")
    blended = ax.get_xaxis_transform()

    for _, row in eras.iterrows():
        start = max(row.start_date, plot_start)
        end = min(row.end_date, plot_end)
        if start > end:
            continue

        ax.axvline(
            row.start_date,
            color="#4d4d4d",
            lw=0.8,
            ls="--",
            alpha=0.55,
            zorder=zorder,
        )

        if label:
            era_name = POLICY_ERA_LABELS.get(
                str(row.policy_era),
                str(row.policy_era).replace("_", " ").upper(),
            )
            ax.text(
                start + (pd.Timedelta(weeks=1)),
                label_y,
                era_name,
                transform=blended,
                ha="left",
                va="top",
                rotation=90,
                rotation_mode="anchor",
                fontsize=6,
                fontweight="bold",
                color="#2304EF",
                clip_on=True,
                zorder=7,
                path_effects=[pe.withStroke(linewidth=1.6, foreground="white")],
            )


def add_policy_strip(ax: Axes, dates: pd.Series):
    """Draw a continuous daily-stringency strip with period-boundary lines."""
    dates = dates.dropna()
    if dates.empty:
        return

    plot_start = dates.min().normalize()
    plot_end = dates.max().normalize()

    daily = (
        POLICY_CALENDAR[["date", "stringency_index"]]
        .dropna(subset=["date"])
        .sort_values("date")
    )
    daily = daily[(daily["date"] >= plot_start) & (daily["date"] <= plot_end)]

    if not daily.empty:
        values = daily["stringency_index"].to_numpy(dtype=float).reshape(1, -1)
        x0 = float(mdates.date2num(daily["date"].min()))
        x1 = float(mdates.date2num(daily["date"].max() + pd.Timedelta(days=1)))
        ax.imshow(
            values,
            aspect="auto",
            cmap=POLICY_STRINGENCY_CMAP,
            norm=POLICY_STRINGENCY_NORM,
            extent=(x0, x1, 0.0, 1.0),
            origin="lower",
            interpolation="nearest",
            zorder=0,
        )

    # Vertical lines marking individual policy-period boundaries, plus code labels.
    for _, row in POLICY_PERIODS.sort_values("start_date").iterrows():
        start = max(row.start_date, plot_start)
        end = min(row.end_date, plot_end)
        if start > end:
            continue

        # Boundary line at the period start.
        if plot_start <= row.start_date <= plot_end:
            ax.axvline(
                row.start_date,
                color="black",
                lw=0.5,
                alpha=0.8,
                zorder=2,
            )

        width_days = (end - start).days + 1
        if width_days >= 18:
            midpoint = start + (end - start) / 2
            ax.text(
                midpoint,
                0.5,
                str(row.period_code),
                ha="center",
                va="center",
                color="black",
                fontsize=6,
                fontweight="bold",
                clip_on=True,
                zorder=3,
            )

    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xlim(
        float(mdates.date2num(plot_start - pd.Timedelta(days=7))),
        float(mdates.date2num(plot_end + pd.Timedelta(days=7))),
    )
    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)


def add_policy_stringency_colorbar(
    fig: Figure,
    ax_policy: Axes,
) -> Axes:
    """Add a slim horizontal stringency colour bar centred above the policy strip."""
    policy_box = ax_policy.get_position()

    bar_width = policy_box.width * 0.28
    bar_height = 0.010
    cax = fig.add_axes(
        (
            policy_box.x0 + (policy_box.width - bar_width) / 2,
            policy_box.y1 + 0.006,
            bar_width,
            bar_height,
        )
    )

    scalar = plt.cm.ScalarMappable(
        norm=POLICY_STRINGENCY_NORM,
        cmap=POLICY_STRINGENCY_CMAP,
    )
    scalar.set_array([])
    cbar = fig.colorbar(scalar, cax=cax, orientation="horizontal")
    cbar.set_label("Restriction stringency", fontsize=6.5, labelpad=3)
    cbar.set_ticks([10, 30, 55, 75, 95])
    cbar.ax.tick_params(labelsize=6.0, length=2.0, width=0.6, pad=1.5)
    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.xaxis.set_label_position("top")
    outline = getattr(cbar, "outline", None)
    if outline is not None:
        try:
            outline.set_linewidth(0.4)
        except TypeError:
            try:
                outline.set_lw(0.4)
            except Exception:
                pass
    return cax


def place_policy_strip_flush(
    ax_policy: Axes,
    ax_top: Axes,
) -> None:
    """Make the policy strip nearly flush with panel A."""
    policy_box = ax_policy.get_position()
    top_box = ax_top.get_position()

    policy_gap = 0.006

    ax_policy.set_position(
        (
            top_box.x0,
            top_box.y1 + policy_gap,
            top_box.width,
            policy_box.height,
        )
    )


def plot_sequences_with_policy(
    timeline: pd.DataFrame,
    sampling_df: pd.DataFrame,
    *,
    ax: Axes,
    show_xlabel: bool = True,
) -> Axes:
    """Plot daily sequences (main axis) and proportion sequenced (twin axis)."""

    dates = timeline["collection_date"]
    counts = timeline["count"]
    smoothed = timeline["smoothed_count"]

    add_epidemic_eras(ax, dates, label=True)

    bar_handle = ax.bar(
        dates,
        counts,
        width=1.0,
        color="#b8b9ba",
        edgecolor="none",
        alpha=0.95,
        label="Daily sequences",
        zorder=5,
    )
    (mean_handle,) = ax.plot(
        dates,
        smoothed,
        color="#0b1f3b",
        lw=1.8,
        label="7-day rolling mean",
        zorder=6,
    )
    ax.set_ylabel("Number of sequences")
    if show_xlabel:
        ax.set_xlabel("Collection date")
    ax.set_facecolor("white")
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.5)
    ax.margins(x=0.01)

    # Twin axis: proportion of positive cases (PCR + antigen tests) sequenced.
    ax2 = ax.twinx()
    (prop_handle,) = ax2.plot(
        sampling_df.index,
        sampling_df["plot_prop_cases_sequenced"].to_numpy(),
        linestyle=":",
        linewidth=1.6,
        alpha=0.9,
        color="#c1272d",
        label="Proportion sequenced",
        zorder=7,
    )
    ax2.set_ylabel("Positive tests sequenced")
    ax2.set_ylim(0, 1)
    ax2.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax2.tick_params(axis="y", colors="#c1272d")
    ax2.grid(False)

    configure_date_axis(ax, dates)

    ax.legend(
        [bar_handle, mean_handle, prop_handle],
        ["Daily sequences", "7-day rolling mean", "Proportion sequenced"],
        loc="upper right",
        ncol=1,
        frameon=False,
    )
    return ax2


def plot_lineage_frequency_and_overtakes(
    df: pd.DataFrame,
    *,
    date_col: str = "collection_date",
    clade_col: str = "variant",
    sequence_col: str = "sequence_id",
    time_freq: str = "W",
    smooth_window: int | None = 3,
    min_sequences_per_period: int = 1,
    ax: Axes,
    legend_ax: Axes | None = None,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    list,
    list,
]:
    """Plot selected lineage-group frequency over time.

    Returns
    -------
    clade_freq
        Unsmoothed clade frequency table.
    plot_freq
        Smoothed clade frequency table used for plotting.
    clade_counts
        Unsmoothed sequence counts by clade and period.
    dominance_df
        Dominant clade group per period.
    overtakes
        Periods where the dominant clade group changed.
    """
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
    clade_order = [clade for clade in CLADE_PALETTE if clade in clade_freq.columns]

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

    ax.set_facecolor("white")
    stack_colors = [
        "white" if clade == "Other" else CLADE_PALETTE.get(clade, "#999999")
        for clade in clade_order
    ]
    stack_handles = ax.stackplot(
        plot_freq.index,
        *[plot_freq[clade].to_numpy() for clade in clade_order],
        labels=clade_order,
        colors=stack_colors,
        alpha=0.92,
        linewidth=0.5,
        edgecolor="white",
        zorder=2,
    )

    # Epidemic-era boundaries (labels live on panel A).
    add_epidemic_eras(ax, dd[date_col], label=False, zorder=5)

    ax.set_xlabel("Collection date", labelpad=3)
    ax.set_ylabel("Clade frequency")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.5)
    ax.margins(x=0.01)

    configure_date_axis(ax, dd[date_col])

    legend_handles = [
        handle for handle, clade in zip(stack_handles, clade_order) if clade != "Other"
    ]
    legend_labels = [clade for clade in clade_order if clade != "Other"]
    if legend_ax is not None:
        legend_ax.axis("off")
        legend_ax.legend(
            legend_handles,
            legend_labels,
            ncol=4,
            loc="center",
            frameon=False,
            columnspacing=1.3,
            handlelength=1.5,
            borderaxespad=0.0,
        )

    return (
        clade_freq,
        plot_freq,
        clade_counts,
        dominance_df,
        overtakes,
        legend_handles,
        legend_labels,
    )


def main() -> int:
    """Load processed metadata and write policy/clade surveillance outputs."""
    # ---- Configuration (previously command-line arguments) ------------------
    smooth_window = DAILY_SMOOTH_WINDOW
    window_stride = SEQUENCE_WINDOW_STRIDE
    figure_dir = FIGURES_DIR
    table_dir = TABLES_DIR
    log_level = "INFO"
    # -------------------------------------------------------------------------

    logging.basicConfig(level=log_level, format="%(levelname)s: %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    if smooth_window < 1:
        raise SystemExit("smooth_window must be at least 1.")
    if window_stride < 1:
        raise SystemExit("window_stride must be at least 1.")

    set_theme(context="paper")
    out_dir = figure_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Loading surveillance sequence data")
    sequences = load_analysis_columns(
        ["sequence_id", "collection_date", "clade", "wn_prop_sequenced"],
        window_stride=window_stride,
    )
    sequences["variant"] = sequences["clade"].map(CLADES).fillna("Other")

    timeline = attach_policy_timeline(
        build_daily_sequence_counts(sequences, smooth_window=smooth_window)
    )

    sampling_df = compute_sequencing_proportion(sequences)

    fig, axes = new_figure(
        width="double",
        height_in=7,
        nrows=3,
        ncols=1,
        gridspec_kw={"height_ratios": [0.15, 2.5, 2.5], "hspace": 0.25},
    )

    axes = axes.ravel()

    ax_policy = axes[0]
    ax_top = axes[1]
    ax_bottom = axes[2]

    add_policy_strip(ax_policy, timeline["collection_date"])
    plot_sequences_with_policy(timeline, sampling_df, ax=ax_top, show_xlabel=False)
    ax_top.tick_params(axis="x", labelbottom=False)

    (
        clade_freq,
        plot_freq,
        clade_counts,
        dominance_df,
        overtakes,
        legend_handles,
        legend_labels,
    ) = plot_lineage_frequency_and_overtakes(
        sequences,
        ax=ax_bottom,
        clade_col="variant",
    )

    place_policy_strip_flush(ax_policy, ax_top)
    add_policy_stringency_colorbar(fig, ax_policy)

    fig.legend(
        legend_handles,
        legend_labels,
        ncol=4,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.05),
        bbox_transform=fig.transFigure,
        frameon=False,
        columnspacing=2,
        handlelength=1.5,
        borderaxespad=0.0,
    )

    tables = {
        "clade_frequency_by_period": clade_freq,
        "clade_frequency_by_period_smoothed": plot_freq,
        "clade_counts_by_period": clade_counts,
        "clade_dominance_by_period": dominance_df,
        "clade_overtake_events": overtakes,
        "sequencing_proportion_by_period": sampling_df,
    }
    for name, table in tables.items():
        LOGGER.info("Writing %s (%s rows)", name, f"{len(table):,}")
        write_table(table, name, table_dir=table_dir)

    fig.align_ylabels([ax_top, ax_bottom])
    add_panel_labels([ax_top, ax_bottom], x=-0.1, y=1.1)

    _ = save_figure(
        fig,
        out_dir / FIGURE_NAME,
        width="double",
        save_png=True,
        save_pdf=True,
    )
    LOGGER.info("Wrote surveillance outputs under %s", out_dir.parent)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
