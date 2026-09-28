"""Compare Scotland's OxCGRT Stringency and Containment and Health indices."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib.dates as mdates
import numpy as np
import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from analyses.surveillance.lib.config import (
    FIGURES_DIR,
    POLICY_INDEX_FIGURE_NAME,
)
from utils import (
    add_panel_labels,
    load_daily_policy_data,
    new_figure,
    save_figure,
)

from .policy_data import build_correlation_summary, load_policy_indices

LOGGER = logging.getLogger(__name__)
STRINGENCY_COLOR = "#2166ac"
CONTAINMENT_COLOR = "#b2182b"
FIGURE_NAME = POLICY_INDEX_FIGURE_NAME


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

    add_panel_labels([ax_time, ax_corr])
    return fig


def build(
    *,
    figure_dir: Path = FIGURES_DIR,
    start_date: str | pd.Timestamp | None = None,
    end_date: str | pd.Timestamp | None = None,
) -> dict[str, object]:
    """Build the policy-index comparison figure."""
    policy_dates = load_daily_policy_data(["date"])["date"]
    if start_date is None:
        start_date = str(policy_dates.min().date())
    if end_date is None:
        end_date = str(policy_dates.max().date())

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

    fig = build_figure(daily, summary)
    saved = save_figure(
        fig,
        figure_dir / POLICY_INDEX_FIGURE_NAME,
        width="double",
        save_png=True,
        save_pdf=True,
    )
    outputs = saved
    LOGGER.info(
        "Wrote policy-index figure: %s",
        ", ".join(map(str, saved.values())),
    )
    return {"figures": outputs}


def main() -> int:
    logging.basicConfig(level="INFO", format="%(levelname)s: %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    build()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
