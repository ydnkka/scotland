"""Publication-ready figures for the Part 3 policy period analysis.

Outputs are written to ``part3/manuscript/figures/`` as PDF, PNG, and TIFF.
The script uses the shared project style module at ``utils/style.py`` and reads
pre-computed tables from ``part3/tables/``.

Main figures
------------
fig1  Weekly cluster outcomes and policy context — two-panel time series
      Panel A: median log cluster size (non-singletons), with policy-period
               background shading colour-coded by intensity.
      Panel B: policy intensity as a stepped line with period code labels.
fig2  Interrupted-time-series plots at three policy transitions — 3×2 panel
      showing pre/post trends in log cluster size (left) and log datazones
      (right) with fitted ITS regression lines.
fig3  Policy-period cluster outcome comparison — dot plot of median log cluster
      size and median log datazones per policy period, annotated with intensity.

Supplementary figures
---------------------
supp_fig1  Weekly mixing metric evolution with policy overlay — 2×2 panel for
           SIMD excess discordance, age excess discordance, mean log datazones,
           and policy intensity.

Run from the repository root:

    conda run -n PhD python part3/manuscript/make_figures.py
"""

from __future__ import annotations

import sys
import os
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Bootstrap repo root
# ---------------------------------------------------------------------------

def _bootstrap_root() -> Path:
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "config.yaml").exists():
            root = str(cand)
            if root not in sys.path:
                sys.path.insert(0, root)
            return cand
    raise FileNotFoundError("Cannot locate config.yaml.")


ROOT = _bootstrap_root()
TABLE_DIR = ROOT / "part3" / "tables"
OUT_DIR   = ROOT / "part3" / "manuscript" / "figures"

import utils.style as style  # noqa: E402

from utils.policy import (  # noqa: E402
    POLICY_PERIODS_PD,
    PERIOD_ORDER,
    PERIOD_LABELS,
    PERIOD_INTENSITY,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# ITS transitions (must match part3_analysis.py)
ITS_TRANSITIONS = [
    ("T1_onset",  pd.Timestamp("2020-10-02"), "P3", "T1",
     "Phase 3 → Pre-tier\n(2020-10-02)"),
    ("L2_to_SL",  pd.Timestamp("2021-04-02"), "L2", "SL",
     "2nd Lockdown → Stay-local\n(2021-04-02)"),
    ("NN_onset",  pd.Timestamp("2021-08-09"), "L0", "NN",
     "Level 0 → Near-normal\n(2021-08-09)"),
]

ITS_OUTCOMES_LEFT  = "log_cluster_size"
ITS_OUTCOMES_RIGHT = "log_datazones"

# Colour map for policy intensity (0–100): dark blue (low) to dark red (high).
_INTENSITY_CMAP = plt.cm.RdYlBu_r

def intensity_colour(v: float) -> tuple:
    """Map a 0–100 intensity value to a colour."""
    return _INTENSITY_CMAP(v / 100.0)

# Background shading alpha for period bands.
PERIOD_ALPHA = 0.18

# Periods present in the study data (excludes pre-July 2020 periods).
STUDY_PERIODS = ["P3", "T1", "F5", "L2", "SL", "L3", "L21", "L0",
                 "NN", "OM", "FE", "PR"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def save_all(fig: plt.Figure, stem: Path, width: str = "double",
             height_in: float = 4.0) -> None:
    """Save as PDF + PNG + TIFF using the shared save_figure helper."""
    style.save_figure(
        fig, stem, width,
        height_in=height_in,
        save_pdf=True, save_png=True, save_tiff=True,
    )


def _add_period_bands(
    ax: plt.Axes,
    ymin: float,
    ymax: float,
    periods: list[str] | None = None,
) -> None:
    """Add colour-coded vertical shading bands for each policy period."""
    if periods is None:
        periods = STUDY_PERIODS
    for code in periods:
        row = POLICY_PERIODS_PD[POLICY_PERIODS_PD["period_code"] == code]
        if row.empty:
            continue
        row = row.iloc[0]
        col = intensity_colour(row["intensity"])
        ax.axvspan(
            row["start_date"], row["end_date"],
            ymin=0, ymax=1,
            color=col, alpha=PERIOD_ALPHA, linewidth=0,
        )


def _period_legend_handles() -> list[mpatches.Patch]:
    """Return legend patches for the intensity scale reference."""
    levels = [(0, "Low restriction (0–25)"),
              (40, "Moderate (25–60)"),
              (80, "High restriction (60–100)")]
    return [
        mpatches.Patch(
            color=intensity_colour(v), alpha=0.5, label=lbl
        )
        for v, lbl in levels
    ]


# ---------------------------------------------------------------------------
# Figure 1 — Weekly time series with policy context
# ---------------------------------------------------------------------------

def plot_weekly_time_series(weekly: pd.DataFrame, out_dir: Path) -> None:
    """Two-panel weekly time series: cluster size and policy intensity."""
    fig, axes = style.new_figure(
        width="double", height_in=5.0, nrows=2, ncols=1,
        font_scale=0.9, gridspec_kw={"height_ratios": [2, 1], "hspace": 0.12},
    )
    ax_top, ax_bot = axes.ravel()

    dates = pd.to_datetime(weekly["week_start"])

    # Panel A: median log cluster size with period shading
    _add_period_bands(ax_top, 0, 1)
    ax_top.plot(
        dates, weekly["median_log_cluster_size"],
        color="#333333", linewidth=1.2, zorder=3, label="Median log cluster size",
    )
    ax_top.set_ylabel("Median log cluster size\n(non-singletons)", fontsize=8)
    ax_top.set_xlim(dates.min(), dates.max())
    ax_top.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax_top.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax_top.tick_params(labelbottom=False)

    # Add ITS transition vlines
    for _, tdate, _, _, desc in ITS_TRANSITIONS:
        ax_top.axvline(tdate, color="#c44e52", linestyle="--", linewidth=0.9,
                       alpha=0.8, zorder=4)
        # short label at top
        ax_top.text(
            tdate, 1.05,
            " " + tdate.strftime("%d %b %Y"),
            fontsize=5.5, color="#c44e52", va="top", rotation=0, ha="center",
            transform=ax_top.get_xaxis_transform(),
        )

    # Legend for period bands
    handles = _period_legend_handles()
    ax_top.legend(handles=handles, loc="upper right", fontsize=6,
                  title="Policy intensity", title_fontsize=6)

    # Panel B: stepped policy intensity line
    _add_period_bands(ax_bot, 0, 1)
    ax_bot.step(
        dates, weekly["dominant_intensity"],
        where="mid", color="#2b2b2b", linewidth=1.3, zorder=3,
    )
    ax_bot.set_ylabel("Policy\nintensity", fontsize=8)
    ax_bot.set_ylim(0, 105)
    ax_bot.set_yticks([0, 25, 50, 75, 100])
    ax_bot.set_xlim(dates.min(), dates.max())
    ax_bot.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax_bot.xaxis.set_major_locator(mdates.MonthLocator(interval=3))

    # Period code labels on intensity panel
    for code in STUDY_PERIODS:
        row = POLICY_PERIODS_PD[POLICY_PERIODS_PD["period_code"] == code]
        if row.empty:
            continue
        row = row.iloc[0]
        mid = row["start_date"] + (row["end_date"] - row["start_date"]) / 2
        ax_bot.text(
            mid, row["intensity"] + 3, code,
            ha="center", va="bottom", fontsize=5, color="#333333",
        )

    style.add_panel_labels([ax_top, ax_bot], x=-0.07, y=1.05)

    fig.subplots_adjust(left=0.10, right=0.97, top=0.93, bottom=0.09)

    save_all(fig, out_dir / "fig1_weekly_time_series", "double", 5.0)
    print("  fig1 saved.")


# ---------------------------------------------------------------------------
# Figure 2 — ITS transition plots
# ---------------------------------------------------------------------------

def plot_its_transitions(
    its_data: dict[str, pd.DataFrame],
    coef_table: pd.DataFrame,
    out_dir: Path,
) -> None:
    """3×2 ITS panel: log cluster size (left) and log datazones (right)."""

    outcomes = [
        (ITS_OUTCOMES_LEFT,  "Median log\ncluster size"),
        (ITS_OUTCOMES_RIGHT, "Median log\ndatazones"),
    ]

    fig, axes = style.new_figure(
        width="double", height_in=7.0, nrows=3, ncols=2,
        font_scale=0.85,
        # gridspec_kw={"hspace": 0.46, "wspace": 0.30},
        layout="constrained",
        sharex=True,
    )

    for row_idx, (label, tdate, pre_code, post_code, desc) in enumerate(ITS_TRANSITIONS):
        its_df = its_data[label]

        for col_idx, (outcome_col, outcome_label) in enumerate(outcomes):
            ax = axes[row_idx, col_idx]
            valid = its_df[["t", "post", "t_post", outcome_col]].dropna()

            # Observed weekly points
            pre  = valid[valid["post"] == 0]
            post = valid[valid["post"] == 1]
            ax.scatter(pre[outcome_col].index.map(lambda i: valid.loc[i, "t"])
                       if False else pre["t"],
                       pre[outcome_col],
                       s=14, color="#4e79a7", zorder=3, label="Pre-transition")
            ax.scatter(post["t"], post[outcome_col],
                       s=14, color="#e15759", zorder=3, label="Post-transition")

            # Fitted ITS regression lines
            row_coef = coef_table[
                (coef_table["transition"] == label) &
                (coef_table["outcome"] == outcome_col)
            ]
            if not row_coef.empty and "coef_const" in row_coef.columns:
                rc = row_coef.iloc[0]
                t_pre  = np.linspace(valid["t"].min(), -0.5, 30)
                t_post = np.linspace(0, valid["t"].max(), 30)
                y_pre  = rc["coef_const"] + rc["coef_t"] * t_pre
                y_post = (rc["coef_const"]
                          + rc["coef_t"] * t_post
                          + rc["coef_post"]
                          + rc["coef_t_post"] * t_post)
                ax.plot(t_pre,  y_pre,  color="#4e79a7", linewidth=1.4, zorder=4)
                ax.plot(t_post, y_post, color="#e15759", linewidth=1.4, zorder=4)

                # Annotate with level change β_post
                p  = rc.get("pval_post", np.nan)
                β  = rc.get("coef_post", np.nan)
                ci_lo = rc.get("ci_lo_post", np.nan)
                ci_hi = rc.get("ci_hi_post", np.nan)
                if not np.isnan(β):
                    sig = ("*" if p < 0.05 else "")
                    ann = f"Δ={β:+.2f} [{ci_lo:+.2f},{ci_hi:+.2f}]{sig}"
                    ax.text(0.03, 0.97, ann, transform=ax.transAxes,
                            fontsize=6, va="top", ha="left",
                            bbox=dict(boxstyle="round,pad=0.2",
                                      fc="white", ec="none", alpha=0.8))

            # Transition line
            ax.axvline(x=-0.5, color="#999999", linestyle="--",
                       linewidth=0.8, zorder=2)

            ax.set_ylabel(outcome_label, fontsize=7)

            # Title only on top row
            if row_idx == 0:
                ax.set_title(outcome_label.replace("\n", " "), fontsize=8, pad=4)

            # Transition descriptor on left column
            if col_idx == 0:
                short = {
                    "T1_onset":  "T1-onset\n(Oct 2020)",
                    "L2_to_SL":  "L2→SL\n(Apr 2021)",
                    "NN_onset":  "NN-onset\n(Aug 2021)",
                }[label]
                ax.set_ylabel(f"{short}\n\n{outcome_label}", fontsize=6.5)

    style.add_panel_labels(axes.ravel(), x=-0.18, y=1.07)

    # Figure-level legend at the top of the figure
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, fontsize=6,
               framealpha=0.7, bbox_to_anchor=(0.5, 1.05),
               bbox_transform=fig.transFigure)

    axes[2, 0].set_xlabel("Week from transition", fontsize=7)
    axes[2, 1].set_xlabel("Week from transition", fontsize=7)


    # fig.subplots_adjust(left=0.15, right=0.97, top=0.93, bottom=0.07)
    save_all(fig, out_dir / "fig2_its_transitions", "double", 7.0)
    print("  fig2 saved.")


# ---------------------------------------------------------------------------
# Figure 3 — Period-level cluster outcome dot plot
# ---------------------------------------------------------------------------

def plot_period_dot_chart(period_desc: pd.DataFrame, out_dir: Path) -> None:
    """Dot chart of median log cluster size and log datazones by policy period."""
    # Compute log-scale medians for plotting
    pd_ = period_desc.copy()
    pd_["log_med_size"]     = np.log(pd_["median_cluster_size"].clip(lower=1))
    pd_["log_med_datazones"] = np.log(pd_["median_datazones"].clip(lower=1))

    # Order by appearance (chronological = PERIOD_ORDER intersection)
    order = [c for c in PERIOD_ORDER if c in pd_["period_code"].values]
    pd_ = pd_.set_index("period_code").loc[order].reset_index()

    fig, axes = style.new_figure(
        width="double", height_in=4.8, nrows=1, ncols=2,
        font_scale=0.9, layout="constrained", sharey=True,
    )
    ax_size, ax_dz = axes.ravel()

    yticks = range(len(pd_))
    ylabels = [f"{row['period_code']}\n{PERIOD_LABELS[row['period_code']]}"
               for _, row in pd_.iterrows()]

    # Colour points by intensity
    colours = [intensity_colour(row["policy_intensity"]) for _, row in pd_.iterrows()]

    for ax, col, xlabel in [
        (ax_size, "log_med_size",     "Log median cluster size (non-singletons)"),
        (ax_dz,   "log_med_datazones", "Log median datazones (non-singletons)"),
    ]:
        for i, (y, val, col_, row) in enumerate(
            zip(yticks, pd_[col], colours, pd_.itertuples())
        ):
            ax.scatter(val, y, color=col_, s=60, zorder=3)
            ax.annotate(
                f"n={row.n_clusters_nonsingleton:,}",
                (val, y), textcoords="offset points",
                xytext=(6, 0), fontsize=5, va="center", color="#666666",
            )

        ax.set_yticks(list(yticks))
        ax.set_yticklabels(ylabels, fontsize=6)
        ax.set_xlabel(xlabel, fontsize=7.5)
        ax.invert_yaxis()
        ax.axvline(0, color="#cccccc", linewidth=0.7, zorder=1)

    # Intensity colourbar
    sm = plt.cm.ScalarMappable(
        cmap=_INTENSITY_CMAP,
        norm=plt.Normalize(vmin=0, vmax=100),
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.6, pad=0.03)
    cbar.set_label("Policy intensity", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    style.add_panel_labels(axes.ravel(), x=-0.15, y=1.04)
    # fig.subplots_adjust(left=0.22, right=0.88, top=0.95, bottom=0.10)

    save_all(fig, out_dir / "fig3_period_outcomes", "double", 4.8)
    print("  fig3 saved.")


# ---------------------------------------------------------------------------
# Supplementary figure 1 — Weekly mixing metrics with policy overlay
# ---------------------------------------------------------------------------

def plot_supp_weekly_mixing(weekly: pd.DataFrame, out_dir: Path) -> None:
    """2×2 weekly series: SIMD mixing, age mixing, log datazones, intensity."""
    dates = pd.to_datetime(weekly["week_start"])

    panels = [
        ("mean_simd_excess",         "Mean SIMD excess\ndiscordance",   "#4e79a7"),
        ("mean_age_excess",          "Mean age excess\ndiscordance",    "#59a14f"),
        ("median_log_datazones",     "Median log\ndatazones",           "#e15759"),
        ("dominant_intensity",       "Policy intensity",                "#2b2b2b"),
    ]

    fig, axes = style.new_figure(
        width="double", height_in=5.6, nrows=2, ncols=2,
        font_scale=0.85, layout="constrained",
        sharex=True,
        # gridspec_kw={"hspace": 0.40, "wspace": 0.35},
    )

    for idx, (col, ylabel, colour) in enumerate(panels):
        ax = axes.ravel()[idx]
        _add_period_bands(ax, 0, 1)

        if col == "dominant_intensity":
            ax.step(dates, weekly[col], where="mid",
                    color=colour, linewidth=1.2, zorder=3)
        else:
            ax.plot(dates, weekly[col], color=colour,
                    linewidth=1.1, zorder=3)

        for _, tdate, _, _, _ in ITS_TRANSITIONS:
            ax.axvline(tdate, color="#c44e52", linestyle="--",
                       linewidth=0.8, alpha=0.7, zorder=4)

        ax.set_ylabel(ylabel, fontsize=7.5)
        ax.set_xlim(dates.min(), dates.max())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=4))

    style.add_panel_labels(axes.ravel(), x=-0.18, y=1.05)
    # fig.subplots_adjust(left=0.12, right=0.97, top=0.93, bottom=0.08)

    save_all(fig, out_dir / "supp_fig1_weekly_mixing", "double", 5.6)
    print("  supp_fig1 saved.")


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def remove_stale_figures(out_dir: Path) -> None:
    """Remove old figure files so stale outputs don't persist."""
    expected = {
        "fig1_weekly_time_series",
        "fig2_its_transitions",
        "fig3_period_outcomes",
        "supp_fig1_weekly_mixing",
    }
    for f in out_dir.glob("*"):
        if f.stem in expected:
            try:
                f.unlink()
            except OSError:
                pass


def run() -> None:
    print("Part 3 figures")
    print("=" * 40)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    remove_stale_figures(OUT_DIR)

    # Load tables
    weekly       = pd.read_csv(TABLE_DIR / "weekly_summaries.csv",
                                parse_dates=["week_start"])
    period_desc  = pd.read_csv(TABLE_DIR / "period_descriptives.csv")
    coef_table   = pd.read_csv(TABLE_DIR / "its_coefficients.csv")

    its_data = {}
    for label, *_ in ITS_TRANSITIONS:
        path = TABLE_DIR / f"its_weekly_{label}.csv"
        its_data[label] = pd.read_csv(path, parse_dates=["week_start"])

    # Generate figures
    plot_weekly_time_series(weekly, OUT_DIR)
    plot_its_transitions(its_data, coef_table, OUT_DIR)
    plot_period_dot_chart(period_desc, OUT_DIR)
    plot_supp_weekly_mixing(weekly, OUT_DIR)

    print(f"\nAll figures written to {OUT_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    run()
