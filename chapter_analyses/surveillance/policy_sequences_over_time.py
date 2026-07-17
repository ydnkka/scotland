"""Plot daily sequenced cases with policy periods and lineage-group frequency."""

from __future__ import annotations

import argparse
import logging
import matplotlib.dates as mdates
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import pandas as pd
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

from pathlib import Path
import sys

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
    set_theme,
    load_analysis_columns,
    load_policy_data,
    add_panel_labels,
    new_figure,
    save_figure,
)

POLICY_STRINGENCY_CMAP = plt.get_cmap("RdYlGn_r")
POLICY_STRINGENCY_NORM = colors.Normalize(
    vmin=1,
    vmax=100,
)
LOGGER = logging.getLogger(__name__)


POLICY_DAILY = load_policy_data()
POLICY_PERIODS = (
    POLICY_DAILY[
        [
            "period_code",
            "period_label",
            "period_start_date",
            "period_end_date",
            "period_order",
            "policy_era",
        ]
    ]
    .drop_duplicates()
    .merge(
        POLICY_DAILY.groupby("period_code", sort=False, observed=True)[
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
        on="period_code",
        how="left",
        validate="one_to_one",
    )
    .rename(
        columns={
            "period_start_date": "start_date",
            "period_end_date": "end_date",
        }
    )
    .sort_values("period_order", ignore_index=True)
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
    lookup = POLICY_DAILY.rename(
        columns={
            "date": "collection_date",
            "period_code": "policy_period",
            "period_label": "policy_period_label",
            "period_start_date": "policy_period_start_date",
            "period_end_date": "policy_period_end_date",
            "period_order": "policy_level",
        }
    )
    return df_full.merge(
        lookup,
        on="collection_date",
        how="left",
        validate="one_to_one",
    )


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


def policy_stringency_color(stringency) -> tuple[float, float, float, float]:
    """Map policy stringency to the shared policy strip colour scale."""
    return POLICY_STRINGENCY_CMAP(POLICY_STRINGENCY_NORM(float(stringency)))


def add_policy_background(
    ax: Axes,
    dates: pd.Series,
) -> None:
    """Shade policy periods behind plotted data."""
    dates = dates.dropna()
    if dates.empty:
        return

    plot_start = dates.min().normalize()
    plot_end = dates.max().normalize()
    policy_periods = POLICY_PERIODS.sort_values("start_date")

    for _, row in enumerate(policy_periods.itertuples(index=False)):
        start = max(row.start_date, plot_start)
        end = min(row.end_date, plot_end)
        if start > end:
            continue

        shade_color = policy_stringency_color(row.policy_stringency)
        shade_alpha = 0.25

        if shade_alpha > 0:
            ax.axvspan(
                start,
                end + pd.Timedelta(days=1),
                color=shade_color,
                alpha=shade_alpha,
                lw=0,
                zorder=-20,
            )

        ax.axvline(start, color="#b3b3b3", lw=0.45, alpha=0.3, zorder=-10)


def add_policy_strip(ax: Axes, dates: pd.Series):
    """Draw a horizontal policy-period strip with stringency-coloured segments."""
    dates = dates.dropna()
    if dates.empty:
        return

    plot_start = dates.min().normalize()
    plot_end = dates.max().normalize()
    policy_periods = POLICY_PERIODS.sort_values("start_date")

    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(False)

    for _, row in policy_periods.iterrows():
        start = max(row.start_date, plot_start)
        end = min(row.end_date, plot_end)
        if start > end:
            continue

        width_days = (end - start).days + 1
        ax.broken_barh(
            [(float(mdates.date2num(start)), float(width_days))],
            (0.08, 0.84),
            facecolors=[policy_stringency_color(row["policy_stringency"])],
            edgecolors="white",
            linewidth=0.45,
        )

        if width_days >= 18:
            midpoint = start + (end - start) / 2
            ax.text(
                midpoint,
                0.5,
                str(row["period_code"]),
                ha="center",
                va="center",
                color="white",
                fontsize=5.5,
                fontweight="bold",
                clip_on=True,
                path_effects=[pe.withStroke(linewidth=0.9, foreground="#333333")],
            )

    ax.set_xlim(
        (plot_start - pd.Timedelta(days=7)).to_pydatetime(),
        (plot_end + pd.Timedelta(days=7)).to_pydatetime(),
    )
    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)


def add_policy_stringency_colorbar(
    fig: Figure,
    ax_policy: Axes,
    ax_top: Axes,
) -> Axes:
    """Add a slim stringency colour bar spanning the policy strip and panel A."""
    policy_box = ax_policy.get_position()
    top_box = ax_top.get_position()
    full_height = policy_box.y1 - top_box.y0
    bar_height = full_height * 0.95
    bar_bottom = top_box.y0 + (full_height - bar_height) / 2
    cax = fig.add_axes(
        (
            policy_box.x1 + 0.006,
            bar_bottom,
            0.010,
            bar_height,
        )
    )

    scalar = plt.cm.ScalarMappable(
        norm=POLICY_STRINGENCY_NORM,
        cmap=POLICY_STRINGENCY_CMAP,
    )
    scalar.set_array([])
    cbar = fig.colorbar(scalar, cax=cax, orientation="vertical")
    cbar.set_label("Restriction stringency", fontsize=7, labelpad=5)
    cbar.set_ticks([10, 30, 55, 75, 95])
    cbar.ax.tick_params(labelsize=6.5, length=2.2, width=0.6, pad=1.5)
    outline = getattr(cbar, "outline", None)
    if outline is not None:
        try:
            outline.set_linewidth(0.4)
        except TypeError:
            try:
                outline.set_lw(0.4)
            except Exception:
                # best-effort: ignore if neither method exists
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
    *,
    ax: Axes,
    show_xlabel: bool = True,
):
    """Plot daily sequences with subtle policy-period shading."""

    dates = timeline["collection_date"]
    counts = timeline["count"]
    smoothed = timeline["smoothed_count"]
    add_policy_background(ax, dates)

    ax.bar(
        dates,
        counts,
        width=1.0,
        color="#b8b9ba",
        edgecolor="none",
        alpha=0.95,
        label="Daily sequences",
        zorder=5,
    )
    ax.plot(
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

    configure_date_axis(ax, dates)
    ax.legend(loc="upper left", ncol=2, frameon=False)


def plot_lineage_frequency_and_overtakes(
    df: pd.DataFrame,
    *,
    date_col: str = "collection_date",
    clade_col: str = "variant",
    sequence_col: str = "sequence_id",
    prop_sequenced_col: str = "wn_prop_sequenced",
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
    pd.DataFrame,
    Axes,
    list,
    list,
]:
    """Plot selected lineage-group frequency and sequencing coverage over time.

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
    sampling_df
        Proportion of cases sequenced per period.
    ax2
        Twin y-axis for proportion of cases sequenced.
    """
    dd = df.copy()

    required_cols = {date_col, clade_col, sequence_col, prop_sequenced_col}
    missing_cols = required_cols - set(dd.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

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

    sampling_df = (
        dd.groupby("time_period")[prop_sequenced_col]
        .mean()
        .reset_index(name="prop_cases_sequenced")
        .sort_values("time_period")
        .set_index("time_period")
    )

    # Keep sampling data aligned to periods that pass the sequence-count threshold.
    sampling_df = sampling_df.reindex(clade_freq.index)

    if smooth_window is not None:
        plot_freq = clade_freq.rolling(window=smooth_window, min_periods=1).mean()
        plot_sampling = (
            sampling_df["prop_cases_sequenced"]
            .rolling(window=smooth_window, min_periods=1)
            .mean()
        )
    else:
        plot_freq = clade_freq.copy()
        plot_sampling = sampling_df["prop_cases_sequenced"].copy()

    sampling_df = sampling_df.copy()
    sampling_df["plot_prop_cases_sequenced"] = plot_sampling

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
    stack_colors = [CLADE_PALETTE.get(clade, "#999999") for clade in clade_order]
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

    ax2 = ax.twinx()
    coverage_line = ax2.plot(
        plot_sampling.index,
        plot_sampling.to_numpy(),
        linestyle=":",
        linewidth=1.6,
        alpha=0.85,
        color="#303030",
        label="Sequenced cases",
    )[0]
    ax2.set_ylabel("Proportion of cases sequenced")
    ax2.set_ylim(0, 1)
    ax2.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax2.tick_params(axis="y", colors="#4d4d4d")
    ax2.grid(False)

    ax.set_xlabel("Collection date", labelpad=3)
    ax.set_ylabel("Clade frequency")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.5)
    ax.margins(x=0.01)

    configure_date_axis(ax, dd[date_col])

    legend_handles = [*stack_handles, coverage_line]
    legend_labels = [*clade_order, "Sequenced cases"]
    if legend_ax is not None:
        legend_ax.axis("off")
        legend_ax.legend(
            legend_handles,
            legend_labels,
            ncol=5,
            loc="center",
            frameon=False,
            columnspacing=1.1,
            handlelength=1.5,
            borderaxespad=0.0,
        )

    return (
        clade_freq,
        plot_freq,
        clade_counts,
        dominance_df,
        overtakes,
        sampling_df,
        ax2,
        legend_handles,
        legend_labels,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=DAILY_SMOOTH_WINDOW,
        help="Centered daily sequence-count smoothing window.",
    )
    parser.add_argument(
        "--window-stride",
        type=int,
        default=SEQUENCE_WINDOW_STRIDE,
        help="Processed analysis-window stride used while loading sequences.",
    )
    parser.add_argument("--figure-dir", type=Path, default=FIGURES_DIR)
    parser.add_argument("--table-dir", type=Path, default=TABLES_DIR)
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    )
    return parser.parse_args()


def main() -> int:
    """Load processed metadata and write policy/clade surveillance outputs."""
    args = parse_args()
    logging.basicConfig(level=args.log_level, format="%(levelname)s: %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    if args.smooth_window < 1:
        raise SystemExit("--smooth-window must be at least 1.")
    if args.window_stride < 1:
        raise SystemExit("--window-stride must be at least 1.")

    set_theme(context="paper")
    out_dir = args.figure_dir
    table_dir = args.table_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Loading surveillance sequence data")
    sequences = load_analysis_columns(
        ["sequence_id", "collection_date", "clade", "wn_prop_sequenced"],
        window_stride=args.window_stride,
        qc=None,
    )
    sequences["variant"] = sequences["clade"].map(CLADES).fillna("Other")

    timeline = attach_policy_timeline(
        build_daily_sequence_counts(sequences, smooth_window=args.smooth_window)
    )

    fig, axes = new_figure(
        width="double",
        height_in=6,
        nrows=3,
        ncols=1,
        gridspec_kw={"height_ratios": [0.16, 2.5, 2.5], "hspace": 0.20},
    )

    axes = axes.ravel()

    ax_policy = axes[0]
    ax_top = axes[1]
    ax_bottom = axes[2]

    add_policy_strip(ax_policy, timeline["collection_date"])
    plot_sequences_with_policy(timeline, ax=ax_top, show_xlabel=False)
    ax_top.tick_params(axis="x", labelbottom=False)
    place_policy_strip_flush(ax_policy, ax_top)
    add_policy_stringency_colorbar(fig, ax_policy, ax_top)

    (
        clade_freq,
        plot_freq,
        clade_counts,
        dominance_df,
        overtakes,
        sampling_df,
        _,
        legend_handles,
        legend_labels,
    ) = plot_lineage_frequency_and_overtakes(
        sequences,
        ax=ax_bottom,
        clade_col="variant",
    )

    fig.legend(
        legend_handles,
        legend_labels,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.075),
        bbox_transform=fig.transFigure,
        frameon=False,
        columnspacing=1.1,
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
    add_panel_labels([ax_top, ax_bottom], x=-0.05, y=1.1)

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
