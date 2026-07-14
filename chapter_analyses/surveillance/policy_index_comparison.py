"""Compare Scotland's OxCGRT Stringency and Containment and Health indices."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib.dates as mdates
import numpy as np
import pandas as pd

from utils import (
    POLICY_PERIODS,
    add_panel_labels,
    load_oxcgrt_containment,
    load_oxcgrt_stringency,
    new_figure,
    save_figure,
    set_theme,
)

from .lib.config import FIGURES_DIR, POLICY_INDEX_FIGURE_NAME, TABLES_DIR
from .lib.io import write_table


LOGGER = logging.getLogger(__name__)
STRINGENCY_COLOR = "#2166ac"
CONTAINMENT_COLOR = "#b2182b"


def load_policy_indices(
    *,
    region_name: str = "Scotland",
    start_date: str | pd.Timestamp | None = None,
    end_date: str | pd.Timestamp | None = None,
) -> pd.DataFrame:
    """Load, align, and optionally date-filter the two daily OxCGRT indices."""
    stringency = load_oxcgrt_stringency(region_name=region_name)
    containment = load_oxcgrt_containment(region_name=region_name)
    daily = stringency.merge(containment, on="date", how="inner", validate="one_to_one")
    daily["date"] = pd.to_datetime(daily["date"], errors="coerce")
    if start_date is not None:
        daily = daily.loc[daily["date"].ge(pd.Timestamp(start_date))]
    if end_date is not None:
        daily = daily.loc[daily["date"].le(pd.Timestamp(end_date))]
    return daily.sort_values("date", ignore_index=True)


def build_correlation_summary(
    daily: pd.DataFrame,
    *,
    region_name: str = "Scotland",
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
                "region_name": region_name,
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


def build_figure(daily: pd.DataFrame, summary: pd.DataFrame):
    """Build a two-column time-series and correlation figure."""
    complete = daily.dropna(subset=["stringency_index", "containment_index"])
    row = summary.iloc[0]
    set_theme(context="paper")
    fig, axes = new_figure(width="double", height_in=3.1, nrows=1, ncols=2)
    ax_time, ax_corr = np.asarray(axes).ravel()

    ax_time.plot(
        daily["date"],
        daily["stringency_index"],
        color=STRINGENCY_COLOR,
        linewidth=1.25,
        label="Stringency Index",
    )
    ax_time.plot(
        daily["date"],
        daily["containment_index"],
        color=CONTAINMENT_COLOR,
        linewidth=1.25,
        label="Containment and Health Index",
    )
    ax_time.set_xlabel("Date")
    ax_time.set_ylabel("OxCGRT index")
    ax_time.set_ylim(-2, 102)
    ax_time.xaxis.set_major_locator(mdates.YearLocator())
    ax_time.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax_time.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.6)
    ax_time.legend(loc="upper right", frameon=False)

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
        f"Pearson $r$ = {row['pearson_r']:.3f}\nSpearman $\\rho$ = {row['spearman_rho']:.3f}\n$n$ = {int(row['n_complete_days']):,} days",
        transform=ax_corr.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2.5},
    )
    ax_corr.legend(loc="lower right", frameon=False)

    add_panel_labels([ax_time, ax_corr], x=-0.12, y=1.06)
    fig.tight_layout(w_pad=2.0)
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-name", default="Scotland")
    parser.add_argument(
        "--start-date",
        default=str(POLICY_PERIODS["start_date"].min().date()),
        help="First included date; defaults to the first defined policy period.",
    )
    parser.add_argument(
        "--end-date",
        default=str(POLICY_PERIODS["end_date"].max().date()),
        help="Last included date; defaults to the last defined policy period.",
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
    args = parse_args()
    logging.basicConfig(level=args.log_level, format="%(levelname)s: %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)

    daily = load_policy_indices(
        region_name=args.region_name,
        start_date=args.start_date,
        end_date=args.end_date,
    )
    summary = build_correlation_summary(daily, region_name=args.region_name)
    LOGGER.info(
        "Index correlation over %s complete days: Pearson r=%.3f; Spearman rho=%.3f",
        summary.loc[0, "n_complete_days"],
        summary.loc[0, "pearson_r"],
        summary.loc[0, "spearman_rho"],
    )

    write_table(daily, "policy_indices_daily", table_dir=args.table_dir)
    write_table(summary, "policy_index_correlation", table_dir=args.table_dir)
    fig = build_figure(daily, summary)
    saved = save_figure(
        fig,
        args.figure_dir / POLICY_INDEX_FIGURE_NAME,
        width="double",
        save_png=True,
        save_pdf=True,
    )
    LOGGER.info("Wrote policy-index figure: %s", ", ".join(map(str, saved.values())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
