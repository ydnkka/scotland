"""Publication-ready manuscript figures for Part 4.

Outputs are written to ``part4/manuscript/figures/`` as PDF, PNG, and TIFF.

Main figures
------------
fig1  Alpha regional emergence across three phases using a health-board
      expansion summary and heatmap.
fig2  Alpha/B.1.177 temporal transition: mutation frequencies, lineage
      displacement, policy periods, and hospital occupancy.
fig3  Counterfactual restriction timing and growth-rate comparison.

Run from the repository root:

    conda run -n PhD python part4/manuscript/make_figures.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Keep Matplotlib/fontconfig cache writes inside the writable sandbox.
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-codex")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd
import pyarrow.parquet as pq


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
TABLE_DIR = ROOT / "part4" / "tables"
OUT_DIR = ROOT / "part4" / "manuscript" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATA_PARQUET = ROOT / "data" / "processed" / "scotland_clustering_analysis_dataset.parquet"

import utils.style as style  # noqa: E402
from utils.policy import POLICY_PERIODS_PD  # noqa: E402


PHASE_SPECS = [
    {
        "phase": "Cryptic GGC chain",
        "windows": [f"W{i:03d}" for i in range(16, 22)],
        "label": "W016-W021\n4 Nov-11 Dec 2020",
        "short": "Cryptic phase",
    },
    {
        "phase": "Multi-region expansion",
        "windows": [f"W{i:03d}" for i in range(22, 25)],
        "label": "W022-W024\n30 Nov-1 Jan",
        "short": "Expansion",
    },
    {
        "phase": "F5/L2 bridge",
        "windows": ["W025"],
        "label": "W025\n19 Dec-8 Jan",
        "short": "F5/L2 bridge",
    },
]

PERIOD_COLOURS = {
    "P3": "#9ecae1",
    "T1": "#fdae6b",
    "F5": "#fb6a4a",
    "L2": "#cb181d",
    "SL": "#74c476",
}

SCENARIO_STYLES = [
    ("Actual", "#c0392b", "-", "Actual"),
    ("L2 from 2 Nov", "#2c7fb8", "--", "L2 from 2 Nov"),
    ("L2 from 2 Dec", "#8856a7", "--", "L2 from 2 Dec"),
    ("L2 from 8 Dec", "#31a354", "--", "L2 from 8 Dec"),
]


def save_all(fig: plt.Figure, stem: str, height_in: float) -> None:
    style.save_figure(
        fig,
        OUT_DIR / stem,
        "double",
        height_in=height_in,
        save_pdf=True,
        save_png=True,
        save_tiff=True,
    )


def add_period_bands(ax: plt.Axes, start: str, end: str, alpha: float = 0.12) -> None:
    start_ts = pd.Timestamp(start)
    end_ts = pd.Timestamp(end)
    periods = POLICY_PERIODS_PD[
        (POLICY_PERIODS_PD["end_date"] >= start_ts) &
        (POLICY_PERIODS_PD["start_date"] <= end_ts)
    ].copy()
    for _, row in periods.iterrows():
        code = row["period_code"]
        ax.axvspan(
            max(row["start_date"], start_ts),
            min(row["end_date"], end_ts),
            color=PERIOD_COLOURS.get(code, "#dddddd"),
            alpha=alpha,
            linewidth=0,
            zorder=0,
        )


def add_policy_lines(ax: plt.Axes, *, show_labels: bool = True) -> None:
    for date, code in [
        ("2020-10-02", "T1"),
        ("2020-11-02", "F5"),
        ("2021-01-05", "L2"),
        ("2021-04-02", "SL"),
    ]:
        ax.axvline(
            pd.Timestamp(date),
            color=PERIOD_COLOURS.get(code, "#444444"),
            linestyle="--",
            linewidth=0.9,
            alpha=0.9,
            zorder=5,
        )
        if show_labels:
            ax.text(
                pd.Timestamp(date),
                0.98,
                code,
                transform=ax.get_xaxis_transform(),
                ha="left",
                va="top",
                fontsize=6.2,
                color=PERIOD_COLOURS.get(code, "#444444"),
                fontweight="bold",
            )


def fmt_date_axis(ax: plt.Axes) -> None:
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax.tick_params(axis="x", labelsize=7)


def load_tables() -> dict[str, pd.DataFrame]:
    return {
        "wpm": pd.read_csv(TABLE_DIR / "part4_window_period_map.csv", parse_dates=["wn_mid_date"]),
        "traj": pd.read_csv(TABLE_DIR / "part4_mutation_trajectories.csv", parse_dates=["wn_mid_date"]),
        "growth": pd.read_csv(TABLE_DIR / "part4_growth_params.csv", parse_dates=["t0"]),
        "proj": pd.read_csv(TABLE_DIR / "part4_counterfactual_projections.csv", parse_dates=["date"]),
        "lineage": pd.read_csv(TABLE_DIR / "part4_lineage_composition.csv", parse_dates=["wn_mid_date"]),
        "hospital": pd.read_csv(TABLE_DIR / "part4_scotland_hospital.csv", parse_dates=["date"]),
        "demo": pd.read_csv(TABLE_DIR / "part4_alpha_phase_demographic_summary.csv"),
    }


def load_alpha_sequences_for_phase_plots() -> pd.DataFrame:
    cols = [
        "sequence_id",
        "window_id",
        "pango_lineage",
        "resolution",
        "nextclade_qc",
        "datazone",
        "dz_health_board",
        "dz_local_authority",
    ]
    df = pq.read_table(str(DATA_PARQUET), columns=cols).to_pandas()
    df = df[
        (df["resolution"] == 0.3) &
        (df["nextclade_qc"] == "good") &
        (df["pango_lineage"].str.startswith("B.1.1.7", na=False))
    ].copy()
    return df


def plot_regional_emergence(tables: dict[str, pd.DataFrame]) -> None:
    alpha = load_alpha_sequences_for_phase_plots()
    demo = tables["demo"].set_index("phase")
    count_rows = []
    summary_rows = []
    max_count = 1

    for spec in PHASE_SPECS:
        sub = alpha[alpha["window_id"].isin(spec["windows"])].drop_duplicates("sequence_id")
        counts = (
            sub.groupby("dz_health_board")["sequence_id"]
            .nunique()
            .rename("n_alpha_sequences")
            .reset_index()
        )
        for _, row in counts.iterrows():
            count_rows.append({
                "health_board": row["dz_health_board"],
                "phase": spec["short"],
                "n_alpha_sequences": int(row["n_alpha_sequences"]),
            })
        summary_rows.append({
            "phase": spec["short"],
            "phase_full": spec["phase"],
            "window": spec["label"].split("\n")[0],
            "dates": spec["label"].split("\n")[1],
            "unique_sequences": int(demo.loc[spec["phase"], "unique_sequences"]),
            "health_boards": int(counts.shape[0]),
            "local_authorities": int(sub["dz_local_authority"].nunique()),
        })
        if not counts.empty:
            max_count = max(max_count, int(counts["n_alpha_sequences"].max()))

    phase_order = [spec["short"] for spec in PHASE_SPECS]
    heat = (
        pd.DataFrame(count_rows)
        .pivot(index="health_board", columns="phase", values="n_alpha_sequences")
        .reindex(columns=phase_order)
        .fillna(0)
        .astype(int)
    )
    first_seen = heat.gt(0).replace(False, np.nan).apply(lambda s: s.first_valid_index(), axis=1)
    first_seen_rank = first_seen.map({phase: i for i, phase in enumerate(phase_order)}).fillna(99)
    heat = heat.assign(
        _first_seen=first_seen_rank,
        _total=heat[phase_order].sum(axis=1),
    ).sort_values(["_first_seen", "_total"], ascending=[True, False])
    heat = heat.drop(columns=["_first_seen", "_total"])
    summary = pd.DataFrame(summary_rows).set_index("phase").loc[phase_order].reset_index()

    bounds = [-0.5, 0.5, 10, 30, 60, 120, max(max_count + 1, 200)]
    tick_locs = [(lo + hi) / 2 for lo, hi in zip(bounds[:-1], bounds[1:])]
    tick_labels = ["0", "1-9", "10-29", "30-59", "60-119", "120+"]
    cmap = ListedColormap(["#f0f0f0", "#fff4b8", "#fed976", "#fd8d3c", "#e31a1c", "#800026"])
    norm = BoundaryNorm(bounds, cmap.N)

    fig, axes = style.new_figure(
        "double",
        height_in=4.8,
        nrows=2,
        ncols=1,
        font_scale=0.82,
        gridspec_kw={"height_ratios": [0.9, 2.6], "hspace": 0.36},
    )
    ax_top, ax_heat = np.ravel(axes)

    x = np.arange(len(summary))
    bar_colours = ["#fdbf6f", "#fb6a4a", "#cb181d"]
    bars = ax_top.bar(
        x,
        summary["unique_sequences"],
        color=bar_colours,
        width=0.58,
        linewidth=0,
    )
    for bar, n_seq in zip(bars, summary["unique_sequences"]):
        if n_seq < 100:
            y_text = bar.get_height() + 14
            va = "bottom"
            text_colour = "#222222"
        else:
            y_text = bar.get_height() - 18
            va = "top"
            text_colour = "white"
        ax_top.text(
            bar.get_x() + bar.get_width() / 2,
            y_text,
            f"{n_seq:,}",
            ha="center",
            va=va,
            fontsize=6.5,
            fontweight="bold",
            color=text_colour,
        )
    ax_top.set_ylabel("Unique Alpha\nsequences")
    ax_top.set_ylim(0, max(summary["unique_sequences"]) * 1.25)
    ax_top.set_xticks(x)
    ax_top.set_xticklabels(
        [f"{row.phase}\n{row.window}" for row in summary.itertuples()],
        fontsize=6.6,
    )
    ax_top.set_title("A  Expansion in sequence count and regional footprint", loc="left")
    ax_top.grid(axis="y", linewidth=0.4, alpha=0.25)

    ax_hb = ax_top.twinx()
    ax_hb.plot(
        x,
        summary["health_boards"],
        color="#2c7fb8",
        marker="o",
        linewidth=1.2,
        markersize=4,
    )
    for xi, hb in zip(x, summary["health_boards"]):
        ax_hb.text(
            xi,
            hb + 0.25,
            f"{hb} HB",
            color="#2c7fb8",
            ha="center",
            va="bottom",
            fontsize=6.2,
        )
    ax_hb.set_ylabel("Health boards", color="#2c7fb8")
    ax_hb.set_ylim(0, max(summary["health_boards"]) + 2)
    ax_hb.tick_params(axis="y", colors="#2c7fb8", labelsize=7)
    ax_hb.spines["right"].set_visible(True)

    matrix = heat[phase_order].values
    im = ax_heat.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    ax_heat.set_title("B  Unique Alpha sequences by health board and phase", loc="left")
    ax_heat.set_yticks(np.arange(len(heat.index)))
    ax_heat.set_yticklabels(heat.index, fontsize=6.4)
    ax_heat.set_xticks(np.arange(len(phase_order)))
    ax_heat.set_xticklabels(
        [f"{row.phase}\n{row.dates}" for row in summary.itertuples()],
        fontsize=6.4,
    )
    ax_heat.set_xticks(np.arange(-0.5, len(phase_order), 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, len(heat.index), 1), minor=True)
    ax_heat.grid(which="minor", color="#ffffff", linewidth=0.8)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    for yi in range(matrix.shape[0]):
        for xi in range(matrix.shape[1]):
            val = int(matrix[yi, xi])
            if val == 0:
                continue
            colour = "white" if val >= 60 else "#222222"
            ax_heat.text(
                xi,
                yi,
                str(val),
                ha="center",
                va="center",
                fontsize=6.1,
                color=colour,
                fontweight="bold" if val >= 60 else "normal",
            )
    ax_heat.set_xlabel("Alpha phase")
    ax_heat.set_ylabel("Health board")

    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.035, pad=0.015)
    cbar.ax.set_title("Unique\nsequences", fontsize=6.1, pad=4)
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(tick_labels)
    cbar.ax.tick_params(labelsize=6)
    fig.suptitle(
        "Figure 1. Regional expansion of Alpha across the seeding, expansion, and F5/L2 bridge phases",
        fontsize=7.7,
        fontweight="bold",
        y=0.98,
    )
    fig.subplots_adjust(left=0.21, right=0.92, bottom=0.10, top=0.89)
    save_all(fig, "fig1_alpha_geographic_emergence", height_in=4.8)


def plot_temporal_transition(tables: dict[str, pd.DataFrame]) -> None:
    traj = tables["traj"]
    lineage = tables["lineage"].copy()
    hospital = tables["hospital"].copy()

    n501 = traj[(traj["mutation"] == "S:N501Y") & (traj["mut_type"] == "amino_acid")].copy()
    s222 = traj[(traj["mutation"] == "S:A222V") & (traj["mut_type"] == "amino_acid")].copy()
    n501 = n501[(n501["wn_mid_date"] >= "2020-10-01") & (n501["wn_mid_date"] <= "2021-04-15")]
    s222 = s222[(s222["wn_mid_date"] >= "2020-10-01") & (s222["wn_mid_date"] <= "2021-04-15")]

    lineage = lineage[
        (lineage["wn_mid_date"] >= "2020-10-01") &
        (lineage["wn_mid_date"] <= "2021-04-15")
    ]
    lineage["lineage_group"] = lineage["pango_lineage"].astype(str).map(
        lambda x: "Alpha" if x.startswith("B.1.1.7")
        else ("B.1.177" if x.startswith("B.1.177") else "Other")
    )
    comp = (
        lineage.groupby(["wn_mid_date", "lineage_group"])["n"]
        .sum()
        .reset_index()
        .pivot(index="wn_mid_date", columns="lineage_group", values="n")
        .fillna(0)
    )
    comp = comp.div(comp.sum(axis=1), axis=0)
    for col in ["Other", "B.1.177", "Alpha"]:
        if col not in comp.columns:
            comp[col] = 0.0
    comp = comp[["Other", "B.1.177", "Alpha"]]

    fig, axes = style.new_figure(
        "double",
        height_in=5.2,
        nrows=2,
        ncols=1,
        font_scale=0.84,
        sharex=True,
        gridspec_kw={"height_ratios": [1.05, 1.0], "hspace": 0.23},
    )
    ax_top, ax_bottom = np.ravel(axes)

    add_period_bands(ax_top, "2020-10-01", "2021-04-15")
    add_policy_lines(ax_top)
    ax_top.plot(
        n501["wn_mid_date"],
        n501["frequency"],
        color="#4e79a7",
        marker="o",
        markersize=3.2,
        linewidth=1.2,
        label="S:N501Y (Alpha marker)",
    )
    ax_top.plot(
        s222["wn_mid_date"],
        s222["frequency"],
        color="#e15759",
        marker="s",
        markersize=3.0,
        linewidth=1.2,
        label="S:A222V (B.1.177 marker)",
    )
    ax_top.annotate(
        "8 Dec\n3.2% to 17.7%",
        xy=(pd.Timestamp("2020-12-08"), 0.177),
        xytext=(pd.Timestamp("2020-11-22"), 0.31),
        arrowprops={"arrowstyle": "->", "linewidth": 0.8, "color": "#333333"},
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.4},
        fontsize=6.5,
        ha="center",
        va="center",
    )
    ax_top.text(
        pd.Timestamp("2021-03-02"),
        0.965,
        "S:N501Y",
        color="#4e79a7",
        fontsize=6.8,
        fontweight="bold",
        ha="left",
        va="center",
    )
    ax_top.text(
        pd.Timestamp("2021-02-03"),
        0.085,
        "S:A222V",
        color="#e15759",
        fontsize=6.8,
        fontweight="bold",
        ha="left",
        va="center",
    )
    ax_top.axhline(0.5, color="#666666", linestyle=":", linewidth=0.8)
    ax_top.set_ylim(-0.03, 1.03)
    ax_top.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
    ax_top.set_ylabel("Mutation frequency")
    ax_top.set_title("A  Marker mutation frequencies show Alpha replacing B.1.177", loc="left")
    ax_top.grid(axis="y", linewidth=0.4, alpha=0.25)

    add_period_bands(ax_bottom, "2020-10-01", "2021-04-15")
    add_policy_lines(ax_bottom, show_labels=False)
    ax_bottom.stackplot(
        comp.index,
        comp["Other"].values,
        comp["B.1.177"].values,
        comp["Alpha"].values,
        colors=["#d9d9d9", "#f28e2b", "#4e79a7"],
        alpha=0.88,
        linewidth=0,
    )
    ax_bottom.set_ylim(0, 1)
    ax_bottom.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
    ax_bottom.set_ylabel("Lineage frequency")
    ax_bottom.set_title("B  Sequenced lineage composition and hospital pressure", loc="left")
    ax_bottom.text(pd.Timestamp("2020-10-10"), 0.14, "Other", color="#666666", fontsize=6.8)
    ax_bottom.text(pd.Timestamp("2020-11-10"), 0.72, "B.1.177", color="#7a4300", fontsize=6.8, fontweight="bold")
    ax_bottom.text(pd.Timestamp("2021-02-03"), 0.80, "Alpha", color="#24476a", fontsize=6.8, fontweight="bold")

    hosp = hospital[(hospital["date"] >= "2020-10-01") & (hospital["date"] <= "2021-04-15")].copy()
    hosp["occupancy_7d"] = hosp["hb_hospital_occupancy"].rolling(7, center=True, min_periods=1).mean()
    ax_hosp = ax_bottom.twinx()
    ax_hosp.plot(
        hosp["date"],
        hosp["occupancy_7d"],
        color="#7f1d1d",
        linewidth=1.1,
        label="Hospital occupancy",
    )
    ax_hosp.set_ylabel("Hospital occupancy", color="#7f1d1d")
    ax_hosp.tick_params(axis="y", colors="#7f1d1d", labelsize=7)
    ax_hosp.spines["right"].set_visible(True)

    fmt_date_axis(ax_bottom)
    ax_bottom.set_xlabel("Window mid-date")
    fig.suptitle(
        "Figure 2. Alpha replaced B.1.177 as hospital pressure rose",
        fontsize=8.6,
        fontweight="bold",
        y=0.98,
    )
    fig.subplots_adjust(top=0.90, bottom=0.11)
    save_all(fig, "fig2_alpha_b1177_transition", height_in=5.2)


def plot_counterfactual_and_growth(tables: dict[str, pd.DataFrame]) -> None:
    traj = tables["traj"]
    growth = tables["growth"].set_index("label")
    proj = tables["proj"]

    n501 = traj[(traj["mutation"] == "S:N501Y") & (traj["mut_type"] == "amino_acid")].copy()
    n501 = n501[(n501["wn_mid_date"] >= "2020-11-01") & (n501["wn_mid_date"] <= "2021-03-15")]

    fig, axes = style.new_figure(
        "double",
        height_in=4.2,
        nrows=1,
        ncols=2,
        font_scale=0.82,
        gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.35},
    )
    ax_left, ax_right = np.ravel(axes)

    add_period_bands(ax_left, "2020-11-01", "2021-03-15")
    add_policy_lines(ax_left)
    ax_left.scatter(
        n501["wn_mid_date"],
        n501["frequency"],
        color="#222222",
        s=14,
        zorder=8,
        label="Observed S:N501Y",
    )
    for scenario_prefix, colour, linestyle, label in SCENARIO_STYLES:
        sub = proj[proj["scenario"].str.startswith(scenario_prefix)]
        ax_left.plot(
            sub["date"],
            sub["frequency"],
            color=colour,
            linestyle=linestyle,
            linewidth=1.4,
            label=label,
        )
    ax_left.axhline(0.5, color="#666666", linestyle=":", linewidth=0.8)
    ax_left.text(pd.Timestamp("2020-11-05"), 0.515, "50% dominance", fontsize=6.5, color="#666666")
    ax_left.set_ylim(-0.03, 1.03)
    ax_left.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
    ax_left.set_ylabel("Alpha frequency")
    ax_left.set_title("A  Counterfactual L2 timing", loc="left")
    ax_left.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=3,
        fontsize=6.3,
        columnspacing=1.1,
        handlelength=2.2,
    )
    ax_left.grid(axis="y", linewidth=0.4, alpha=0.25)
    fmt_date_axis(ax_left)

    labels = ["Alpha\nF5", "Alpha\nL2", "B.1.177\nL2 decline"]
    slopes = [
        float(growth.loc["Alpha_F5", "slope"]),
        float(growth.loc["Alpha_L2", "slope"]),
        abs(float(growth.loc["B1177_L2_decline", "slope"])),
    ]
    colours = ["#fb6a4a", "#cb181d", "#2c7fb8"]
    yerr = np.array([
        [
            slopes[0] - float(growth.loc["Alpha_F5", "ci_lo"]),
            slopes[1] - float(growth.loc["Alpha_L2", "ci_lo"]),
            0.0,
        ],
        [
            float(growth.loc["Alpha_F5", "ci_hi"]) - slopes[0],
            float(growth.loc["Alpha_L2", "ci_hi"]) - slopes[1],
            0.0,
        ],
    ])
    bars = ax_right.bar(labels, slopes, color=colours, yerr=yerr, capsize=3, linewidth=0)
    doubling = [
        float(growth.loc["Alpha_F5", "doubling_days"]),
        float(growth.loc["Alpha_L2", "doubling_days"]),
        float(growth.loc["B1177_L2_decline", "doubling_days"]),
    ]
    for bar, days, text in zip(bars, doubling, ["doubling", "doubling", "halving"]):
        ax_right.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.006,
            f"{days:.1f}d\n{text}",
            ha="center",
            va="bottom",
            fontsize=6.8,
        )
    reduction = (1 - slopes[1] / slopes[0]) * 100
    ax_right.text(
        0.02,
        0.95,
        f"L2 growth {reduction:.0f}% slower than F5",
        transform=ax_right.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
    )
    ax_right.set_ylabel("Absolute growth/decline rate per day")
    ax_right.set_ylim(0, max(slopes) * 1.55)
    ax_right.set_title("B  Estimated growth rates", loc="left")
    ax_right.grid(axis="y", linewidth=0.4, alpha=0.25)

    fig.suptitle(
        "Figure 3. Earlier L2 timing delays fitted Alpha dominance",
        fontsize=8.6,
        fontweight="bold",
        y=0.98,
    )
    fig.subplots_adjust(top=0.86, bottom=0.25)
    save_all(fig, "fig3_counterfactual_growth", height_in=4.2)


def main() -> None:
    tables = load_tables()
    print("Generating Part 4 manuscript figures...")
    plot_regional_emergence(tables)
    print("  saved fig1_alpha_geographic_emergence")
    plot_temporal_transition(tables)
    print("  saved fig2_alpha_b1177_transition")
    plot_counterfactual_and_growth(tables)
    print("  saved fig3_counterfactual_growth")
    print(f"Done. Figures written to {OUT_DIR}")


if __name__ == "__main__":
    main()
