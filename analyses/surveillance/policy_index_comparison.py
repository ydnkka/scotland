"""Compare Scotland's OxCGRT Stringency and Containment and Health indices."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib.dates as mdates
import numpy as np
import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from analyses.surveillance.lib.config import (
    FIGURES_DIR,
    POLICY_INDEX_FIGURE_NAME,
    TABLES_DIR,
)
from analyses.surveillance.lib.io import write_table
from utils import (
    add_panel_labels,
    load_daily_policy_data,
    new_figure,
    save_figure,
)

LOGGER = logging.getLogger(__name__)
STRINGENCY_COLOR = "#2166ac"
CONTAINMENT_COLOR = "#b2182b"


def load_policy_indices(
    start_date: str | pd.Timestamp | None = None,
    end_date: str | pd.Timestamp | None = None,
) -> pd.DataFrame:
    """Load and optionally date-filter Scotland's processed policy table."""
    daily = load_daily_policy_data()
    if start_date is not None:
        daily = daily.loc[daily["date"].ge(pd.Timestamp(start_date))]
    if end_date is not None:
        daily = daily.loc[daily["date"].le(pd.Timestamp(end_date))]
    return daily.sort_values("date", ignore_index=True)


def build_correlation_summary(
    daily: pd.DataFrame,
) -> pd.DataFrame:
    """Summarise complete-day agreement without independence-based p-values."""
    complete = daily.dropna(subset=["stringency_index", "containment_index"]).copy()
    if len(complete) < 2:
        raise ValueError("At least two complete daily index pairs are required.")

    x = complete["stringency_index"].to_numpy(dtype=float)
    y = complete["containment_index"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, deg=1)
    pearson = float(complete["stringency_index"].corr(complete["containment_index"]))
    spearman = float(
        complete["stringency_index"].corr(
            complete["containment_index"],
            method="spearman",
        )
    )
    return pd.DataFrame(
        [
            {
                "start_date": complete["date"].min(),
                "end_date": complete["date"].max(),
                "n_complete_days": len(complete),
                "pearson_r": pearson,
                "spearman_rho": spearman,
                "linear_slope": float(slope),
                "linear_intercept": float(intercept),
                "pearson_r_squared": pearson**2,
            }
        ]
    )


def configure_date_axis(ax, dates: pd.Series) -> None:
    """Use the quarterly month/year ticks from the surveillance timeline."""
    dates = pd.to_datetime(dates, errors="coerce").dropna()
    if dates.empty:
        return
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("center")


def build_figure(daily: pd.DataFrame, summary: pd.DataFrame):
    """Build vertically stacked time-series and correlation panels."""
    complete = daily.dropna(subset=["stringency_index", "containment_index"])
    row = summary.iloc[0]

    fig, axes = new_figure(
        width="double",
        height_in=6.6,
        nrows=2,
        ncols=1,
        gridspec_kw={"hspace": 0.34},
    )
    ax_time, ax_corr = np.asarray(axes).ravel()
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.08, top=0.94)

    # --- Panel A: time series -------------------------------------------------
    ax_time.plot(
        daily["date"],
        daily["stringency_index"],
        color=STRINGENCY_COLOR,
        linewidth=1.25,
        label="Stringency Index",
        zorder=3,
    )
    ax_time.plot(
        daily["date"],
        daily["containment_index"],
        color=CONTAINMENT_COLOR,
        linewidth=1.25,
        label="Containment and Health Index",
        zorder=3,
    )
    ax_time.set_xlabel("Date")
    ax_time.set_ylabel("OxCGRT index")
    ax_time.set_ylim(-2, 102)
    configure_date_axis(ax_time, daily["date"])
    margin = pd.Timedelta(days=7)
    ax_time.set_xlim(
        daily["date"].min().normalize() - margin,
        daily["date"].max().normalize() + margin,
    )
    ax_time.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.6)
    ax_time.legend(loc="upper right", frameon=False)

    # --- Panel B: correlation -------------------------------------------------
    ax_corr.scatter(
        complete["stringency_index"],
        complete["containment_index"],
        s=9,
        color="#525252",
        alpha=0.28,
        linewidths=0,
        label="Daily values",
    )
    x_line = np.linspace(0, 100, 200)
    y_line = row["linear_intercept"] + row["linear_slope"] * x_line
    ax_corr.plot(
        x_line,
        y_line,
        color="#000000",
        linewidth=1.25,
        label="Linear fit",
    )
    ax_corr.plot(
        x_line,
        x_line,
        color="#969696",
        linewidth=0.8,
        linestyle="--",
        label="Identity",
    )
    ax_corr.set_xlim(-2, 102)
    ax_corr.set_ylim(-2, 102)
    ax_corr.set_xlabel("Stringency Index")
    ax_corr.set_ylabel("Containment and Health Index")
    ax_corr.grid(color="#e5e5e5", linewidth=0.5, alpha=0.55)
    ax_corr.text(
        0.04,
        0.96,
        f"Pearson $r$ = {row['pearson_r']:.3f}\n"
        f"Spearman $\\rho$ = {row['spearman_rho']:.3f}\n"
        f"$ n $ = {int(row['n_complete_days']):,} days",
        transform=ax_corr.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2.5},
    )
    ax_corr.legend(loc="lower right", frameon=False)

    add_panel_labels([ax_time, ax_corr], x=-0.1, y=1.1)
    return fig


def main() -> int:
    # ---- Configuration (previously command-line arguments) ------------------
    policy_dates = load_daily_policy_data(["date"])["date"]
    start_date = str(policy_dates.min().date())
    end_date = str(policy_dates.max().date())
    figure_dir = FIGURES_DIR
    table_dir = TABLES_DIR
    log_level = "INFO"
    # -------------------------------------------------------------------------

    logging.basicConfig(level=log_level, format="%(levelname)s: %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)

    daily = load_policy_indices(
        start_date=start_date,
        end_date=end_date,
    )
    summary = build_correlation_summary(daily)
    LOGGER.info(
        "Index correlation over %s complete days: Pearson r=%.3f; Spearman rho=%.3f",
        summary.loc[0, "n_complete_days"],
        summary.loc[0, "pearson_r"],
        summary.loc[0, "spearman_rho"],
    )

    write_table(daily, "policy_indices_daily", table_dir=table_dir)
    write_table(summary, "policy_index_correlation", table_dir=table_dir)
    fig = build_figure(daily, summary)
    saved = save_figure(
        fig,
        figure_dir / POLICY_INDEX_FIGURE_NAME,
        width="double",
        save_png=True,
        save_pdf=True,
    )
    LOGGER.info("Wrote policy-index figure: %s", ", ".join(map(str, saved.values())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())










# """Compare Scotland's OxCGRT Stringency and Containment and Health indices."""

# from __future__ import annotations

# import argparse
# import logging
# from pathlib import Path
# import sys

# import matplotlib.dates as mdates
# import matplotlib.patheffects as pe
# import matplotlib.pyplot as plt
# from matplotlib import colors
# import numpy as np
# import pandas as pd

# _PROJECT_ROOT = Path(__file__).resolve().parents[2]
# if str(_PROJECT_ROOT) not in sys.path:
#     sys.path.insert(0, str(_PROJECT_ROOT))

# from analyses.surveillance.lib.config import (  # noqa: E402
#     FIGURES_DIR,
#     POLICY_INDEX_FIGURE_NAME,
#     TABLES_DIR,
# )
# from analyses.surveillance.lib.io import write_table  # noqa: E402

# from utils import (  # noqa: E402
#     add_panel_labels,
#     load_daily_policy_data,
#     new_figure,
#     save_figure,
# )


# LOGGER = logging.getLogger(__name__)
# STRINGENCY_COLOR = "#2166ac"
# CONTAINMENT_COLOR = "#b2182b"
# POLICY_STRINGENCY_CMAP = plt.get_cmap("RdYlGn_r")
# POLICY_STRINGENCY_NORM = colors.Normalize(vmin=0, vmax=100)


# def load_policy_indices(
#     start_date: str | pd.Timestamp | None = None,
#     end_date: str | pd.Timestamp | None = None,
# ) -> pd.DataFrame:
#     """Load and optionally date-filter Scotland's processed policy table."""
#     daily = load_daily_policy_data()
#     if start_date is not None:
#         daily = daily.loc[daily["date"].ge(pd.Timestamp(start_date))]
#     if end_date is not None:
#         daily = daily.loc[daily["date"].le(pd.Timestamp(end_date))]
#     return daily.sort_values("date", ignore_index=True)


# def build_correlation_summary(
#     daily: pd.DataFrame,
# ) -> pd.DataFrame:
#     """Summarise complete-day agreement without independence-based p-values."""
#     complete = daily.dropna(subset=["stringency_index", "containment_index"]).copy()
#     if len(complete) < 2:
#         raise ValueError("At least two complete daily index pairs are required.")

#     x = complete["stringency_index"].to_numpy(dtype=float)
#     y = complete["containment_index"].to_numpy(dtype=float)
#     slope, intercept = np.polyfit(x, y, deg=1)
#     pearson = float(complete["stringency_index"].corr(complete["containment_index"]))
#     spearman = float(
#         complete["stringency_index"].corr(
#             complete["containment_index"],
#             method="spearman",
#         )
#     )
#     return pd.DataFrame(
#         [
#             {
#                 "start_date": complete["date"].min(),
#                 "end_date": complete["date"].max(),
#                 "n_complete_days": len(complete),
#                 "pearson_r": pearson,
#                 "spearman_rho": spearman,
#                 "linear_slope": float(slope),
#                 "linear_intercept": float(intercept),
#                 "pearson_r_squared": pearson**2,
#             }
#         ]
#     )


# def build_policy_period_summary() -> pd.DataFrame:
#     """Return one ordered row per period with its mean daily stringency."""
#     policy = load_daily_policy_data()
#     descriptors = policy[
#         [
#             "period_code",
#             "period_label",
#             "period_start_date",
#             "period_end_date",
#             "period_order",
#         ]
#     ].drop_duplicates()
#     means = (
#         policy.groupby("period_code", sort=False, observed=True)["stringency_index"]
#         .mean()
#         .rename("policy_stringency")
#         .reset_index()
#     )
#     return (
#         descriptors.merge(means, on="period_code", how="left", validate="one_to_one")
#         .sort_values("period_order")
#         .reset_index(drop=True)
#     )


# def policy_stringency_color(value: float):
#     """Map mean period stringency to the shared policy colour scale."""
#     return POLICY_STRINGENCY_CMAP(POLICY_STRINGENCY_NORM(float(value)))


# def add_policy_background(
#     ax,
#     dates: pd.Series,
#     periods: pd.DataFrame,
# ) -> None:
#     """Shade the time-series background by policy period."""
#     plot_start = pd.to_datetime(dates, errors="coerce").min().normalize()
#     plot_end = pd.to_datetime(dates, errors="coerce").max().normalize()
#     for row in periods.itertuples(index=False):
#         start = max(row.period_start_date, plot_start)
#         end = min(row.period_end_date, plot_end)
#         if start > end:
#             continue
#         ax.axvspan(
#             start,
#             end + pd.Timedelta(days=1),
#             color=policy_stringency_color(row.policy_stringency),
#             alpha=0.20,
#             linewidth=0,
#             zorder=-20,
#         )
#         ax.axvline(start, color="#b3b3b3", linewidth=0.45, alpha=0.35, zorder=-10)


# def add_policy_strip(
#     ax,
#     dates: pd.Series,
#     periods: pd.DataFrame,
# ) -> None:
#     """Draw the policy-period strip above the time-series panel."""
#     plot_start = pd.to_datetime(dates, errors="coerce").min().normalize()
#     plot_end = pd.to_datetime(dates, errors="coerce").max().normalize()
#     ax.set_ylim(0, 1)
#     ax.set_yticks([])
#     ax.set_facecolor("white")
#     for spine in ax.spines.values():
#         spine.set_visible(False)

#     for row in periods.itertuples(index=False):
#         start = max(row.period_start_date, plot_start)
#         end = min(row.period_end_date, plot_end)
#         if start > end:
#             continue
#         width_days = (end - start).days + 1
#         ax.broken_barh(
#             [(float(mdates.date2num(start)), float(width_days))],
#             (0.08, 0.84),
#             facecolors=[policy_stringency_color(row.policy_stringency)],
#             edgecolors="white",
#             linewidth=0.45,
#         )
#         if width_days >= 18:
#             midpoint = start + (end - start) / 2
#             ax.text(
#                 midpoint,
#                 0.5,
#                 str(row.period_code),
#                 ha="center",
#                 va="center",
#                 color="white",
#                 fontsize=6,
#                 fontweight="bold",
#                 clip_on=True,
#                 path_effects=[pe.withStroke(linewidth=0.9, foreground="#333333")],
#             )

#     margin = pd.Timedelta(days=7)
#     ax.set_xlim(plot_start - margin, plot_end + margin)
#     ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)


# def configure_date_axis(ax, dates: pd.Series) -> None:
#     """Use the quarterly month/year ticks from the surveillance timeline."""
#     dates = pd.to_datetime(dates, errors="coerce").dropna()
#     if dates.empty:
#         return
#     ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
#     ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
#     ax.xaxis.set_minor_locator(mdates.MonthLocator())
#     for label in ax.get_xticklabels():
#         label.set_horizontalalignment("center")


# def build_figure(daily: pd.DataFrame, summary: pd.DataFrame):
#     """Build vertically stacked time-series and correlation panels."""
#     complete = daily.dropna(subset=["stringency_index", "containment_index"])
#     row = summary.iloc[0]
#     periods = build_policy_period_summary()
#     fig, axes = new_figure(
#         width="double",
#         height_in=6.6,
#         nrows=2,
#         ncols=1,
#         gridspec_kw={"hspace": 0.34},
#     )
#     ax_time, ax_corr = np.asarray(axes).ravel()
#     fig.subplots_adjust(left=0.10, right=0.98, bottom=0.08, top=0.90)
#     time_box = ax_time.get_position()
#     ax_policy = fig.add_axes((time_box.x0, time_box.y1 + 0.006, time_box.width, 0.022))

#     add_policy_strip(ax_policy, daily["date"], periods)
#     add_policy_background(ax_time, daily["date"], periods)

#     ax_time.plot(
#         daily["date"],
#         daily["stringency_index"],
#         color=STRINGENCY_COLOR,
#         linewidth=1.25,
#         label="Stringency Index",
#         zorder=3,
#     )
#     ax_time.plot(
#         daily["date"],
#         daily["containment_index"],
#         color=CONTAINMENT_COLOR,
#         linewidth=1.25,
#         label="Containment and Health Index",
#         zorder=3,
#     )
#     ax_time.set_xlabel("Date")
#     ax_time.set_ylabel("OxCGRT index")
#     ax_time.set_ylim(-2, 102)
#     configure_date_axis(ax_time, daily["date"])
#     ax_time.set_xlim(ax_policy.get_xlim())
#     ax_time.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.6)
#     ax_time.legend(loc="upper right", frameon=False)

#     ax_corr.scatter(
#         complete["stringency_index"],
#         complete["containment_index"],
#         s=9,
#         color="#525252",
#         alpha=0.28,
#         linewidths=0,
#         label="Daily values",
#     )
#     x_line = np.linspace(0, 100, 200)
#     y_line = row["linear_intercept"] + row["linear_slope"] * x_line
#     ax_corr.plot(
#         x_line,
#         y_line,
#         color="#000000",
#         linewidth=1.25,
#         label="Linear fit",
#     )
#     ax_corr.plot(
#         x_line,
#         x_line,
#         color="#969696",
#         linewidth=0.8,
#         linestyle="--",
#         label="Identity",
#     )
#     ax_corr.set_xlim(-2, 102)
#     ax_corr.set_ylim(-2, 102)
#     ax_corr.set_xlabel("Stringency Index")
#     ax_corr.set_ylabel("Containment and Health Index")
#     ax_corr.grid(color="#e5e5e5", linewidth=0.5, alpha=0.55)
#     ax_corr.text(
#         0.04,
#         0.96,
#         f"Pearson $r$ = {row['pearson_r']:.3f}\nSpearman $\\rho$ = {row['spearman_rho']:.3f}\n$n$ = {int(row['n_complete_days']):,} days",
#         transform=ax_corr.transAxes,
#         ha="left",
#         va="top",
#         fontsize=8,
#         bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2.5},
#     )
#     ax_corr.legend(loc="lower right", frameon=False)

#     add_panel_labels([ax_time, ax_corr], x=-0.07, y=1.02)
#     return fig


# def parse_args() -> argparse.Namespace:
#     parser = argparse.ArgumentParser(description=__doc__)
#     policy_dates = load_daily_policy_data(["date"])["date"]
#     parser.add_argument(
#         "--start-date",
#         default=str(policy_dates.min().date()),
#         help="First included date; defaults to the first defined policy period.",
#     )
#     parser.add_argument(
#         "--end-date",
#         default=str(policy_dates.max().date()),
#         help="Last included date; defaults to the last defined policy period.",
#     )
#     parser.add_argument("--figure-dir", type=Path, default=FIGURES_DIR)
#     parser.add_argument("--table-dir", type=Path, default=TABLES_DIR)
#     parser.add_argument(
#         "--log-level",
#         default="INFO",
#         choices=("DEBUG", "INFO", "WARNING", "ERROR"),
#     )
#     return parser.parse_args()


# def main() -> int:
#     args = parse_args()
#     logging.basicConfig(level=args.log_level, format="%(levelname)s: %(message)s")
#     logging.getLogger("fontTools").setLevel(logging.WARNING)

#     daily = load_policy_indices(
#         start_date=args.start_date,
#         end_date=args.end_date,
#     )
#     summary = build_correlation_summary(daily)
#     LOGGER.info(
#         "Index correlation over %s complete days: Pearson r=%.3f; Spearman rho=%.3f",
#         summary.loc[0, "n_complete_days"],
#         summary.loc[0, "pearson_r"],
#         summary.loc[0, "spearman_rho"],
#     )

#     write_table(daily, "policy_indices_daily", table_dir=args.table_dir)
#     write_table(summary, "policy_index_correlation", table_dir=args.table_dir)
#     fig = build_figure(daily, summary)
#     saved = save_figure(
#         fig,
#         args.figure_dir / POLICY_INDEX_FIGURE_NAME,
#         width="double",
#         save_png=True,
#         save_pdf=True,
#     )
#     LOGGER.info("Wrote policy-index figure: %s", ", ".join(map(str, saved.values())))
#     return 0


# if __name__ == "__main__":
#     raise SystemExit(main())
