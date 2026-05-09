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

POLICY_COLORS = {
    "P3": "#8fbf88",
    "T1": "#d6a03a",
    "F5": "#c97a32",
    "L2": "#9d3c32",
    "SL": "#d6a03a",
    "L0": "#71a9c9",
    "NN": "#8fbf88",
    "OM": "#b7a4d8",
    "FE": "#a7c6a1",
    "PR": "#d0d0d0",
}

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


def add_policy_spans(
    ax: plt.Axes,
    start: pd.Timestamp,
    end: pd.Timestamp,
    *,
    label_codes: bool = False,
    y_text: float = 0.98,
) -> None:
    periods = policy.POLICY_PERIODS_PD.copy()
    for _, row in periods.iterrows():
        left = max(pd.Timestamp(row["start_date"]), start)
        right = min(pd.Timestamp(row["end_date"]), end)
        if right < left:
            continue
        code = row["period_code"]
        color = POLICY_COLORS.get(code, "#d5d5d5")
        alpha = 0.12 if code in FOCUS_PERIODS else 0.055
        ax.axvspan(left, right, color=color, alpha=alpha, lw=0, zorder=0)
        if label_codes and code in FOCUS_PERIODS:
            ax.text(
                left + (right - left) / 2,
                y_text,
                code,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=6.5,
                color="#333333",
            )


def format_date_axis(ax: plt.Axes, *, monthly: bool = False) -> None:
    locator = mdates.MonthLocator(interval=1 if monthly else 2)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax.tick_params(axis="x", labelrotation=0)


def plot_policy_timeline() -> None:
    weekly = read_csv("weekly_summaries.csv", ["wn_mid_date"])
    start = weekly["wn_mid_date"].min() - pd.Timedelta(days=7)
    end = weekly["wn_mid_date"].max() + pd.Timedelta(days=7)

    fig, ax = style.new_figure(width="double", height_in=3.25, font_scale=0.9)
    add_policy_spans(ax, start, end, label_codes=True)

    ax.plot(
        weekly["wn_mid_date"],
        weekly["median_log_cluster_size"],
        color="#1f4e79",
        lw=1.8,
        marker="o",
        ms=2.4,
        label="Median log cluster size",
        zorder=3,
    )
    ax.set_ylabel("Median log cluster size")
    ax.set_xlabel("")
    ax.set_xlim(start, end)
    format_date_axis(ax)

    ax2 = ax.twinx()
    ax2.step(
        weekly["wn_mid_date"],
        weekly["policy_intensity"],
        where="mid",
        color="#333333",
        alpha=0.45,
        lw=1.0,
        label="Policy intensity",
    )
    ax2.set_ylim(0, 105)
    ax2.set_ylabel("Policy intensity")
    ax2.spines["right"].set_visible(True)

    for label, date in [
        ("T1", pd.Timestamp("2020-10-02")),
        ("L2", pd.Timestamp("2021-01-05")),
        ("SL", pd.Timestamp("2021-04-02")),
        ("NN", pd.Timestamp("2021-08-09")),
    ]:
        ax.axvline(date, color="#333333", lw=0.7, ls="--", alpha=0.65)
        ax.text(
            date,
            1.04,
            label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=7,
        )

    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc="upper left", ncol=2)
    save_all(fig, FIGURE_DIR / "fig1_policy_timeline_cluster_structure", height_in=3.25)


def plot_selected_its() -> None:
    fig, axes = style.new_figure(
        width="double",
        height_in=6.4,
        nrows=3,
        ncols=2,
        font_scale=0.82,
        sharex=False,
    )
    outcomes = ["median_log_cluster_size", "median_log_datazones"]

    for r, (slug, (short_label, date)) in enumerate(TRANSITIONS.items()):
        its = read_csv(f"its_weekly_{slug}.csv", ["wn_mid_date"])
        start = its["wn_mid_date"].min() - pd.Timedelta(days=3)
        end = its["wn_mid_date"].max() + pd.Timedelta(days=3)
        for c, outcome in enumerate(outcomes):
            ax = axes[r, c]
            add_policy_spans(ax, start, end, label_codes=False)
            ax.scatter(
                its["wn_mid_date"],
                its[outcome],
                s=16,
                color="#1f4e79",
                alpha=0.85,
                zorder=3,
            )
            fitted = f"fitted_{outcome}"
            ax.plot(
                its["wn_mid_date"],
                its[fitted],
                color="#b23a2e",
                lw=1.7,
                zorder=4,
            )
            ax.axvline(date, color="#333333", lw=0.8, ls="--")
            if c == 0:
                ax.set_ylabel(f"{short_label}\n{OUTCOME_LABELS[outcome]}")
            else:
                ax.set_ylabel(OUTCOME_LABELS[outcome])
            ax.set_xlim(start, end)
            format_date_axis(ax, monthly=True)
            if r == 0:
                ax.set_title(OUTCOME_LABELS[outcome])

    style.add_panel_labels(axes.ravel(), x=-0.08, y=1.08, size=9)
    save_all(fig, FIGURE_DIR / "fig2_selected_policy_transitions", height_in=6.4)


def plot_alpha_emergence() -> None:
    traj = read_csv("alpha_mutation_trajectories.csv", ["wn_mid_date"])
    hb = read_csv("alpha_health_board_weekly.csv", ["wn_mid_date"])

    start = pd.Timestamp("2020-10-01")
    end = pd.Timestamp("2021-04-20")
    traj = traj[traj["wn_mid_date"].between(start, end)].copy()
    hb = hb[hb["wn_mid_date"].between(start, end)].copy()

    fig, axes = style.new_figure(
        width="double",
        height_in=4.9,
        nrows=2,
        ncols=1,
        font_scale=0.85,
        sharex=False,
        gridspec_kw={"height_ratios": [1.0, 1.15]},
    )
    ax = axes[0]
    add_policy_spans(ax, start, end, label_codes=True)
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
    ax.axvline(pd.Timestamp("2021-01-05"), color="#333333", lw=0.8, ls="--")
    ax.set_ylabel("Mutation frequency")
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlim(start, end)
    format_date_axis(ax, monthly=True)
    ax.legend(loc="upper left", ncol=2)

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

    ax = axes[1]
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
        ax.legend(loc="upper left", ncol=3, fontsize=6.2)
    ax.axvline(pd.Timestamp("2021-01-05"), color="#333333", lw=0.8, ls="--")
    add_policy_spans(ax, start, end, label_codes=False)
    ax.set_ylabel("Alpha sequences")
    ax.set_xlim(start, end)
    format_date_axis(ax, monthly=True)

    style.add_panel_labels(axes.ravel(), x=-0.08, y=1.08, size=9)
    save_all(fig, FIGURE_DIR / "fig3_alpha_emergence_f5_l2", height_in=4.9)


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
        height_in=3.7,
        nrows=1,
        ncols=2,
        font_scale=0.85,
        gridspec_kw={"width_ratios": [1.35, 1.0]},
    )

    ax = axes[0]
    add_policy_spans(ax, start, end, label_codes=True)
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
    ax.axhline(0.5, color="#333333", lw=0.8, ls=":")
    ax.set_ylabel("S:N501Y frequency")
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlim(start, end)
    format_date_axis(ax, monthly=True)
    ax.legend(loc="upper left", fontsize=6.1)

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
    ax.set_ylabel("Log-odds slope per week")
    ax.set_xlim(-0.5, max(2.5, len(sub) - 0.5))

    style.add_panel_labels(axes.ravel(), x=-0.12, y=1.08, size=9)
    save_all(fig, FIGURE_DIR / "fig4_alpha_counterfactual_timing", height_in=3.7)


def plot_supplementary_its_mixing() -> None:
    fig, axes = style.new_figure(
        width="double",
        height_in=6.4,
        nrows=3,
        ncols=2,
        font_scale=0.82,
        sharex=False,
    )
    outcomes = ["mean_simd_excess_discordance", "mean_age_excess_discordance"]
    for r, (slug, (short_label, date)) in enumerate(TRANSITIONS.items()):
        its = read_csv(f"its_weekly_{slug}.csv", ["wn_mid_date"])
        start = its["wn_mid_date"].min() - pd.Timedelta(days=3)
        end = its["wn_mid_date"].max() + pd.Timedelta(days=3)
        for c, outcome in enumerate(outcomes):
            ax = axes[r, c]
            add_policy_spans(ax, start, end, label_codes=False)
            ax.scatter(its["wn_mid_date"], its[outcome], s=16, color="#1f4e79", alpha=0.85)
            ax.plot(its["wn_mid_date"], its[f"fitted_{outcome}"], color="#b23a2e", lw=1.7)
            ax.axvline(date, color="#333333", lw=0.8, ls="--")
            if c == 0:
                ax.set_ylabel(f"{short_label}\n{OUTCOME_LABELS[outcome]}")
            else:
                ax.set_ylabel(OUTCOME_LABELS[outcome])
            ax.set_xlim(start, end)
            format_date_axis(ax, monthly=True)
            if r == 0:
                ax.set_title(OUTCOME_LABELS[outcome])
    style.add_panel_labels(axes.ravel(), x=-0.08, y=1.08, size=9)
    save_all(fig, FIGURE_DIR / "supp_fig1_its_mixing_outcomes", height_in=6.4)


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
