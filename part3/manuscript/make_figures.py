"""Create publication-ready figures for Part 3 policy-period analyses.

Outputs are written to ``part3/manuscript/figures`` as PDF, PNG, and TIFF.

Run from the repository root:

    conda run -n PhD python part3/manuscript/make_figures.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib import colors as mpl_colors
from matplotlib import patheffects
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils import policy, style


TABLE_DIR = ROOT / "part3" / "tables"
FIGURE_DIR = ROOT / "part3" / "manuscript" / "figures"

FOCUS_PERIODS = {"P3", "T1", "F5", "L2", "SL", "L0", "NN"}
TRANSITIONS = {
    "t1_onset": ("T1", pd.Timestamp("2020-10-02")),
    "l2_to_sl": ("SL", pd.Timestamp("2021-04-02")),
    "nn_onset": ("NN", pd.Timestamp("2021-08-09")),
}
TRANSITION_ROW_LABELS = {
    "t1_onset": "P3 -> T1",
    "l2_to_sl": "L2 -> SL",
    "nn_onset": "L0 -> NN",
}

POLICY_INTENSITY_CMAP = plt.get_cmap("RdYlGn_r")
POLICY_INTENSITY_NORM = mpl_colors.Normalize(
    vmin=policy.POLICY_PERIODS_PD["intensity"].min(),
    vmax=policy.POLICY_PERIODS_PD["intensity"].max(),
)

OUTCOME_LABELS = {
    "median_log_cluster_size": "Median log cluster size",
    "median_log_datazones": "Median log datazones",
    "mean_simd_excess_discordance": "Mean SIMD excess discordance",
    "mean_age_excess_discordance": "Mean age excess discordance",
}


def setup_environment() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def save_all(fig: plt.Figure, out_base: Path, *, width: str = "double", height_in: float = 3.8) -> None:
    style.save_figure(
        fig,
        out_base,
        width=width,
        height_in=height_in,
        dpi=600,
        save_pdf=True,
        save_png=True,
        save_tiff=True,
    )


def read_csv(name: str, date_cols: list[str] | None = None) -> pd.DataFrame:
    path = TABLE_DIR / name
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run conda run -n PhD python part3/part3_analysis.py first."
        )
    df = pd.read_csv(path)
    for col in date_cols or []:
        if col in df:
            df[col] = pd.to_datetime(df[col])
    return df


def policy_intensity_color(intensity: float) -> tuple[float, float, float, float]:
    return POLICY_INTENSITY_CMAP(POLICY_INTENSITY_NORM(float(intensity)))


def add_policy_spans(
    ax: plt.Axes,
    start: pd.Timestamp,
    end: pd.Timestamp,
    *,
    label_codes: bool = False,
    y_text: float = 0.98,
    alpha: float = 0.075,
) -> None:
    periods = policy.POLICY_PERIODS_PD.copy()
    for _, row in periods.iterrows():
        left = max(pd.Timestamp(row["start_date"]), start)
        right = min(pd.Timestamp(row["end_date"]), end)
        if right < left:
            continue
        code = row["period_code"]
        ax.axvspan(
            left,
            right + pd.Timedelta(days=1),
            color=policy_intensity_color(row["intensity"]),
            alpha=alpha if code in FOCUS_PERIODS else alpha * 0.65,
            lw=0,
            zorder=-20,
        )
        ax.axvline(left, color="#b8b8b8", lw=0.45, alpha=0.30, zorder=-10)
        if label_codes and code in FOCUS_PERIODS:
            ax.text(
                left + (right - left) / 2,
                y_text,
                code,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=6,
                fontweight="bold",
                color="#222222",
                path_effects=[patheffects.withStroke(linewidth=1.5, foreground="white")],
            )
    ax.set_axisbelow(True)


def add_policy_strip(
    ax: plt.Axes,
    start: pd.Timestamp,
    end: pd.Timestamp,
    *,
    label_codes: bool = True,
) -> None:
    """Draw an intensity-coloured policy strip, matching surveillance figures."""
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(False)

    periods = policy.POLICY_PERIODS_PD.copy()
    for _, row in periods.iterrows():
        left = max(pd.Timestamp(row["start_date"]), start.normalize())
        right = min(pd.Timestamp(row["end_date"]), end.normalize())
        if right < left:
            continue

        width_days = (right - left).days + 1
        ax.broken_barh(
            [(mdates.date2num(left), width_days)],
            (0.10, 0.80),
            facecolors=[policy_intensity_color(row["intensity"])],
            edgecolors="white",
            linewidth=0.45,
        )

        if label_codes and width_days >= 18:
            midpoint = left + (right - left) / 2
            ax.text(
                midpoint,
                0.50,
                str(row["period_code"]),
                ha="center",
                va="center",
                fontsize=6.5,
                fontweight="bold",
                color="white",
                clip_on=True,
                path_effects=[patheffects.withStroke(linewidth=0.9, foreground="#333333")],
            )

    ax.set_xlim(start, end)
    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)


def add_policy_intensity_colorbar(
    fig: plt.Figure,
    ax_top: plt.Axes,
    ax_bottom: plt.Axes,
    *,
    label: str = "Restriction intensity",
) -> plt.Axes:
    """Add a slim shared policy intensity colour bar beside a set of axes."""
    top_box = ax_top.get_position()
    bottom_box = ax_bottom.get_position()

    full_height = top_box.y1 - bottom_box.y0
    bar_height = full_height * 0.85
    bar_bottom = bottom_box.y0 + (full_height - bar_height) / 2
    cax = fig.add_axes([
        top_box.x1 + 0.006,
        bar_bottom,
        0.010,
        bar_height,
    ])

    scalar = plt.cm.ScalarMappable(
        norm=POLICY_INTENSITY_NORM,
        cmap=POLICY_INTENSITY_CMAP,
    )
    scalar.set_array([])
    cbar = fig.colorbar(scalar, cax=cax, orientation="vertical")
    cbar.set_label(label, fontsize=7, labelpad=5)
    cbar.set_ticks([10, 30, 55, 75, 95])
    cbar.ax.tick_params(labelsize=6.5, length=2.2, width=0.6, pad=1.5)
    cbar.outline.set_linewidth(0.4)
    return cax


def place_policy_strip_flush(ax_policy: plt.Axes, ax_top: plt.Axes, *, gap: float = 0.004) -> None:
    """Position a policy strip directly above its primary axis."""
    policy_box = ax_policy.get_position()
    top_box = ax_top.get_position()
    ax_policy.set_position([
        top_box.x0,
        top_box.y1 + gap,
        top_box.width,
        policy_box.height,
    ])


def format_date_axis(ax: plt.Axes, *, monthly: bool = False) -> None:
    locator = mdates.MonthLocator(interval=1 if monthly else 3)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax.tick_params(axis="x", labelrotation=0)


def plot_policy_timeline() -> None:
    weekly = read_csv("weekly_summaries.csv", ["wn_mid_date"])
    start = weekly["wn_mid_date"].min() - pd.Timedelta(days=7)
    end = weekly["wn_mid_date"].max() + pd.Timedelta(days=7)

    fig, axes = style.new_figure(
        width="double",
        height_in=3.55,
        nrows=2,
        ncols=1,
        font_scale=0.9,
        sharex=False,
        gridspec_kw={"height_ratios": [0.12, 1.0]},
    )
    fig.subplots_adjust(hspace=0.02, right=0.90)
    ax_policy, ax = axes
    add_policy_strip(ax_policy, start, end)
    add_policy_spans(ax, start, end, label_codes=False, alpha=0.055)

    ax.plot(
        weekly["wn_mid_date"],
        weekly["median_log_cluster_size"],
        color="#1f4e79",
        lw=1.8,
        marker="o",
        ms=2.4,
        zorder=3,
    )
    ax.set_ylabel("Median log cluster size")
    ax.set_xlabel("")
    ax.set_xlim(start, end)
    format_date_axis(ax)

    for label, date in [
        ("T1", pd.Timestamp("2020-10-02")),
        ("L2", pd.Timestamp("2021-01-05")),
        ("SL", pd.Timestamp("2021-04-02")),
        ("NN", pd.Timestamp("2021-08-09")),
    ]:
        ax.axvline(date, color="#333333", lw=0.7, ls="--", alpha=0.70, zorder=6)
        ax.text(
            date,
            0.98,
            label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=6.7,
            color="#222222",
            path_effects=[patheffects.withStroke(linewidth=1.4, foreground="white")],
        )

    place_policy_strip_flush(ax_policy, ax)
    add_policy_intensity_colorbar(fig, ax_policy, ax)
    save_all(fig, FIGURE_DIR / "fig1_policy_timeline_cluster_structure", height_in=3.55)


def plot_selected_its() -> None:
    fig, axes = style.new_figure(
        width="double",
        height_in=6.65,
        nrows=3,
        ncols=2,
        font_scale=0.82,
        sharex=False,
    )
    fig.subplots_adjust(top=0.91, hspace=0.44, wspace=0.24)
    outcomes = ["median_log_cluster_size", "median_log_datazones"]

    for r, (slug, (short_label, date)) in enumerate(TRANSITIONS.items()):
        its = read_csv(f"its_weekly_{slug}.csv", ["wn_mid_date"])
        start = its["wn_mid_date"].min() - pd.Timedelta(days=3)
        end = its["wn_mid_date"].max() + pd.Timedelta(days=3)
        for c, outcome in enumerate(outcomes):
            ax = axes[r, c]
            add_policy_spans(ax, start, end, label_codes=False, alpha=0.060)
            ax.scatter(
                its["wn_mid_date"],
                its[outcome],
                s=16,
                color="#1f4e79",
                alpha=0.85,
                zorder=5,
            )
            fitted = f"fitted_{outcome}"
            ax.plot(
                its["wn_mid_date"],
                its[fitted],
                color="#b23a2e",
                lw=1.7,
                zorder=6,
            )
            ax.axvline(date, color="#333333", lw=0.8, ls="--", zorder=7)
            if c == 0:
                ax.set_ylabel(OUTCOME_LABELS[outcome])
                ax.text(
                    0.015,
                    0.96,
                    TRANSITION_ROW_LABELS.get(slug, short_label),
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=6.8,
                    fontweight="bold",
                    color="#222222",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.6},
                    zorder=8,
                )
            else:
                ax.set_ylabel(OUTCOME_LABELS[outcome])
            ax.set_xlim(start, end)
            format_date_axis(ax, monthly=True)
            if r == len(TRANSITIONS) - 1:
                ax.set_xlabel("Week midpoint")
            if r == 0:
                ax.set_title(OUTCOME_LABELS[outcome])

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#1f4e79",
               markeredgecolor="#1f4e79", markersize=4.2, label="Weekly outcome"),
        Line2D([0], [0], color="#b23a2e", lw=1.7, label="Segmented fit"),
        Line2D([0], [0], color="#333333", lw=0.8, ls="--", label="Transition date"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=3,
        frameon=False,
    )
    style.add_panel_labels(axes.ravel(), x=-0.14, y=1.08, size=9)
    save_all(fig, FIGURE_DIR / "fig2_selected_policy_transitions", height_in=6.65)


def plot_alpha_emergence() -> None:
    traj = read_csv("alpha_mutation_trajectories.csv", ["wn_mid_date"])
    hb = read_csv("alpha_health_board_weekly.csv", ["wn_mid_date"])

    start = pd.Timestamp("2020-10-01")
    end = pd.Timestamp("2021-04-20")
    traj = traj[traj["wn_mid_date"].between(start, end)].copy()
    hb = hb[hb["wn_mid_date"].between(start, end)].copy()

    fig, axes = style.new_figure(
        width="double",
        height_in=5.15,
        nrows=3,
        ncols=1,
        font_scale=0.85,
        sharex=False,
        gridspec_kw={"height_ratios": [0.15, 1.0, 1.15]},
    )
    fig.subplots_adjust(hspace=0.13, right=0.91)
    ax_policy = axes[0]
    ax = axes[1]
    add_policy_strip(ax_policy, start, end)
    add_policy_spans(ax, start, end, label_codes=False, alpha=0.060)
    ax.plot(
        traj["wn_mid_date"],
        traj["freq_s_n501y"],
        color="#4e79a7",
        marker="o",
        ms=3,
        lw=1.8,
        label="S:N501Y",
    )
    ax.plot(
        traj["wn_mid_date"],
        traj["freq_s_a222v"],
        color="#d55e00",
        marker="o",
        ms=3,
        lw=1.8,
        label="S:A222V",
    )
    l2_start = pd.Timestamp("2021-01-05")
    ax.axvline(l2_start, color="#333333", lw=0.8, ls="--", zorder=7)
    ax.text(
        l2_start,
        0.98,
        "L2",
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=6.7,
        color="#222222",
        path_effects=[patheffects.withStroke(linewidth=1.4, foreground="white")],
    )
    ax.set_ylabel("Mutation frequency")
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlim(start, end)
    format_date_axis(ax, monthly=True)
    ax.tick_params(axis="x", labelbottom=False)
    ax.text(
        0.985,
        0.88,
        "S:N501Y",
        transform=ax.transAxes,
        ha="right",
        va="center",
        color="#4e79a7",
        fontsize=6.8,
        fontweight="bold",
        path_effects=[patheffects.withStroke(linewidth=1.4, foreground="white")],
    )
    ax.text(
        0.985,
        0.12,
        "S:A222V",
        transform=ax.transAxes,
        ha="right",
        va="center",
        color="#d55e00",
        fontsize=6.8,
        fontweight="bold",
        path_effects=[patheffects.withStroke(linewidth=1.4, foreground="white")],
    )

    ax_occ = ax.twinx()
    if "hb_hospital_occupancy_total" in traj:
        ax_occ.plot(
            traj["wn_mid_date"],
            traj["hb_hospital_occupancy_total"],
            color="#555555",
            lw=1.0,
            alpha=0.35,
            label="Hospital occupancy",
        )
        ax_occ.set_ylabel("Hospital occupancy")
        ax_occ.spines["right"].set_visible(True)

    ax = axes[2]
    add_policy_spans(ax, start, end, label_codes=False, alpha=0.060)
    if not hb.empty:
        totals = hb.groupby("health_board")["n_alpha_sequences"].sum().sort_values(ascending=False)
        top = list(totals.head(7).index)
        hb["health_board_plot"] = np.where(hb["health_board"].isin(top), hb["health_board"], "Other")
        pivot = (
            hb.groupby(["wn_mid_date", "health_board_plot"])["n_alpha_sequences"]
            .sum()
            .unstack(fill_value=0)
            .sort_index()
        )
        order = [name for name in top if name in pivot.columns]
        if "Other" in pivot.columns:
            order.append("Other")
        pivot = pivot.reindex(columns=order, fill_value=0)
        palette = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#b07aa1", "#edc948", "#9c755f", "#bab0ac"]
        ax.stackplot(
            pivot.index,
            [pivot[col].to_numpy() for col in pivot.columns],
            labels=pivot.columns,
            colors=palette[: len(pivot.columns)],
            alpha=0.86,
        )
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(0.008, 0.985),
            ncol=4,
            fontsize=5.8,
            frameon=True,
            facecolor="white",
            edgecolor="none",
            framealpha=0.78,
            borderpad=0.25,
            handlelength=1.4,
            columnspacing=0.8,
        )
    ax.axvline(l2_start, color="#333333", lw=0.8, ls="--", zorder=7)
    ax.set_ylabel("Alpha sequences")
    ax.set_xlim(start, end)
    format_date_axis(ax, monthly=True)
    ax.set_xlabel("Week midpoint")

    place_policy_strip_flush(ax_policy, axes[1])
    style.add_panel_labels([axes[1], axes[2]], x=-0.08, y=1.04, size=9)
    save_all(fig, FIGURE_DIR / "fig3_alpha_emergence_f5_l2", height_in=5.15)


def plot_counterfactuals() -> None:
    traj = read_csv("alpha_mutation_trajectories.csv", ["wn_mid_date"])
    cf = read_csv("alpha_counterfactual_trajectories.csv", ["date", "requested_switch_date"])
    params = read_csv("alpha_growth_params.csv")

    start = pd.Timestamp("2020-11-01")
    end = pd.Timestamp("2021-03-10")
    traj = traj[traj["wn_mid_date"].between(start, end)].copy()
    cf = cf[cf["date"].between(start, end)].copy()

    fig, axes = style.new_figure(
        width="double",
        height_in=3.85,
        nrows=1,
        ncols=2,
        font_scale=0.85,
        gridspec_kw={"width_ratios": [1.35, 1.0]},
    )
    fig.subplots_adjust(wspace=0.28, top=0.90, bottom=0.18)

    ax = axes[0]
    add_policy_spans(ax, start, end, label_codes=False, alpha=0.060)
    ax.scatter(
        traj["wn_mid_date"],
        traj["freq_s_n501y"],
        color="#222222",
        s=18,
        zorder=5,
        label="Observed S:N501Y",
    )
    scenario_order = [
        "actual_l2_start",
        "expansion_date_2020_12_08",
        "nearest_w021_2020_12_02",
        "f5_start_2020_11_02",
    ]
    colors = {
        "actual_l2_start": "#4e79a7",
        "expansion_date_2020_12_08": "#f28e2b",
        "nearest_w021_2020_12_02": "#59a14f",
        "f5_start_2020_11_02": "#e15759",
    }
    for scenario in scenario_order:
        dat = cf[cf["scenario"] == scenario]
        if dat.empty:
            continue
        ax.plot(
            dat["date"],
            dat["projected_n501y_frequency"],
            color=colors[scenario],
            lw=1.5,
            label=dat["scenario_label"].iloc[0],
        )
    l2_start = pd.Timestamp("2021-01-05")
    ax.axhline(0.5, color="#333333", lw=0.8, ls=":", zorder=2)
    ax.axvline(l2_start, color="#333333", lw=0.8, ls="--", alpha=0.75, zorder=3)
    ax.text(
        l2_start,
        0.98,
        "L2",
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=6.7,
        color="#222222",
        path_effects=[patheffects.withStroke(linewidth=1.4, foreground="white")],
    )
    ax.text(
        0.985,
        0.505,
        "50%",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=6.4,
        color="#333333",
        path_effects=[patheffects.withStroke(linewidth=1.4, foreground="white")],
    )
    ax.set_title("Counterfactual S:N501Y timing")
    ax.set_ylabel("S:N501Y frequency")
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlim(start, end)
    format_date_axis(ax, monthly=True)
    ax.legend(
        loc="lower right",
        fontsize=5.8,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.80,
        borderpad=0.25,
        handlelength=1.5,
    )

    ax = axes[1]
    order = ["alpha_f5_n501y", "alpha_l2_n501y", "b1177_l2_a222v"]
    labels = ["Alpha\nF5", "Alpha\nL2", "S:A222V\nL2"]
    sub = params[params["analysis"].isin(order)].copy()
    sub["order"] = sub["analysis"].map({k: i for i, k in enumerate(order)})
    sub = sub.sort_values("order")
    x = np.arange(len(sub))
    y = sub["slope_per_week"].to_numpy()
    yerr = np.vstack(
        [
            y - sub["slope_ci_low_per_week"].to_numpy(),
            sub["slope_ci_high_per_week"].to_numpy() - y,
        ]
    )
    ax.axhline(0, color="#333333", lw=0.8)
    ax.errorbar(
        x,
        y,
        yerr=yerr,
        fmt="o",
        color="#1f4e79",
        ecolor="#1f4e79",
        capsize=3,
        lw=1.4,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels[: len(sub)])
    ax.set_title("Growth-rate comparison")
    ax.set_ylabel("Log-odds slope per week")
    ax.set_xlim(-0.5, max(2.5, len(sub) - 0.5))

    style.add_panel_labels(axes.ravel(), x=-0.12, y=1.06, size=9)
    save_all(fig, FIGURE_DIR / "fig4_alpha_counterfactual_timing", height_in=3.85)


def plot_supplementary_its_mixing() -> None:
    fig, axes = style.new_figure(
        width="double",
        height_in=6.65,
        nrows=3,
        ncols=2,
        font_scale=0.82,
        sharex=False,
    )
    fig.subplots_adjust(top=0.91, hspace=0.44, wspace=0.24)
    outcomes = ["mean_simd_excess_discordance", "mean_age_excess_discordance"]
    for r, (slug, (short_label, date)) in enumerate(TRANSITIONS.items()):
        its = read_csv(f"its_weekly_{slug}.csv", ["wn_mid_date"])
        start = its["wn_mid_date"].min() - pd.Timedelta(days=3)
        end = its["wn_mid_date"].max() + pd.Timedelta(days=3)
        for c, outcome in enumerate(outcomes):
            ax = axes[r, c]
            add_policy_spans(ax, start, end, label_codes=False, alpha=0.060)
            ax.scatter(its["wn_mid_date"], its[outcome], s=16, color="#1f4e79", alpha=0.85, zorder=5)
            ax.plot(its["wn_mid_date"], its[f"fitted_{outcome}"], color="#b23a2e", lw=1.7, zorder=6)
            ax.axvline(date, color="#333333", lw=0.8, ls="--", zorder=7)
            if c == 0:
                ax.set_ylabel(OUTCOME_LABELS[outcome])
                ax.text(
                    0.015,
                    0.96,
                    TRANSITION_ROW_LABELS.get(slug, short_label),
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=6.8,
                    fontweight="bold",
                    color="#222222",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.6},
                    zorder=8,
                )
            else:
                ax.set_ylabel(OUTCOME_LABELS[outcome])
            ax.set_xlim(start, end)
            format_date_axis(ax, monthly=True)
            if r == len(TRANSITIONS) - 1:
                ax.set_xlabel("Week midpoint")
            if r == 0:
                ax.set_title(OUTCOME_LABELS[outcome])

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#1f4e79",
               markeredgecolor="#1f4e79", markersize=4.2, label="Weekly outcome"),
        Line2D([0], [0], color="#b23a2e", lw=1.7, label="Segmented fit"),
        Line2D([0], [0], color="#333333", lw=0.8, ls="--", label="Transition date"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=3,
        frameon=False,
    )
    style.add_panel_labels(axes.ravel(), x=-0.14, y=1.08, size=9)
    save_all(fig, FIGURE_DIR / "supp_fig1_its_mixing_outcomes", height_in=6.65)


def main() -> None:
    setup_environment()
    plot_policy_timeline()
    plot_selected_its()
    plot_alpha_emergence()
    plot_counterfactuals()
    plot_supplementary_its_mixing()
    print(f"Wrote Part 3 manuscript figures to {FIGURE_DIR}")


if __name__ == "__main__":
    main()
