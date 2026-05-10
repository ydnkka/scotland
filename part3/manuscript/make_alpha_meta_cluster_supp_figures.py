"""Create supplementary manuscript figures for the pre-L2 Alpha meta-cluster analysis.

These figures extend the Part 3 Alpha/F5-L2 case study with the exploratory
connected-component analysis from ``part3/notebooks``.  They are descriptive:
meta-clusters represent continuity of rolling-window genomic clusters through
shared sequence membership, not proven single transmission chains.

Run from the repository root:

    conda run -n PhD python part3/manuscript/make_alpha_meta_cluster_supp_figures.py
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
import seaborn as sns
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils import style


FIGURE_DIR = ROOT / "part3" / "manuscript" / "figures"
META_TABLE_DIR_CANDIDATES = [
    ROOT / "part3" / "notebooks" / "tables",
    ROOT / "part3" / "tables" / "alpha_meta_clusters",
]

L2_START = pd.Timestamp("2021-01-05")
TOP_N = 6


def setup_environment() -> Path:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for candidate in META_TABLE_DIR_CANDIDATES:
        if (candidate / "alpha_pre_l2_meta_cluster_summary.csv").exists():
            return candidate
    raise FileNotFoundError(
        "Could not find Alpha meta-cluster tables. Run "
        "part3/notebooks/alpha_pre_l2_meta_cluster_network.ipynb and "
        "part3/notebooks/alpha_top6_meta_cluster_demographics_over_time.ipynb first."
    )


def save_all(fig: plt.Figure, out_base: Path, *, height_in: float) -> None:
    style.save_figure(
        fig,
        out_base,
        width="double",
        height_in=height_in,
        dpi=600,
        save_pdf=True,
        save_png=True,
        save_tiff=True,
    )


def read_meta_csv(table_dir: Path, name: str, date_cols: list[str] | None = None) -> pd.DataFrame:
    path = table_dir / name
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path)
    for col in date_cols or []:
        if col in df:
            df[col] = pd.to_datetime(df[col])
    return df


def week_start(dates: pd.Series) -> pd.Series:
    dates = pd.to_datetime(dates)
    return dates - pd.to_timedelta(dates.dt.weekday, unit="D")


def format_week_axis(ax: plt.Axes, *, interval: int = 2) -> None:
    ax.xaxis.set_major_locator(mdates.WeekdayLocator(interval=interval))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%Y"))
    ax.tick_params(axis="x", labelrotation=0)


def format_month_axis(ax: plt.Axes) -> None:
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax.tick_params(axis="x", labelrotation=0)


def meta_label(meta: str, sizes: pd.Series) -> str:
    return f"{meta} (n={int(sizes.loc[meta])})"


def plot_alpha_meta_cluster_amplification(table_dir: Path) -> None:
    summary = read_meta_csv(
        table_dir,
        "alpha_pre_l2_meta_cluster_summary.csv",
        ["first_collection_date", "last_collection_date"],
    ).sort_values("n_alpha_sequences", ascending=False)
    membership = read_meta_csv(
        table_dir,
        "alpha_pre_l2_meta_cluster_sequence_membership.csv",
        ["collection_date"],
    )
    impact = read_meta_csv(table_dir, "alpha_pre_l2_meta_cluster_signature_impact_summary.csv")
    trajectories = read_meta_csv(
        table_dir,
        "alpha_pre_l2_signature_mutation_trajectories.csv",
        ["collection_week"],
    )

    top_meta = summary.head(TOP_N)["meta_cluster_id"].tolist()
    sizes = summary.set_index("meta_cluster_id")["n_alpha_sequences"]
    top_colors = dict(zip(top_meta, sns.color_palette("tab10", TOP_N)))
    other_color = "#c9c9c9"

    membership = membership.copy()
    membership["collection_week"] = week_start(membership["collection_date"])
    membership["plot_group"] = membership["meta_cluster_id"].where(
        membership["meta_cluster_id"].isin(top_meta),
        f"Other {len(summary) - TOP_N} meta-clusters",
    )

    fig, axes = style.new_figure(
        width="double",
        height_in=6.2,
        nrows=2,
        ncols=2,
        font_scale=0.78,
        gridspec_kw={"height_ratios": [1.0, 1.05], "width_ratios": [1.0, 1.18]},
    )
    fig.subplots_adjust(hspace=0.42, wspace=0.30)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    ranked = summary.reset_index(drop=True).copy()
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    for _, row in ranked.iterrows():
        meta = row["meta_cluster_id"]
        color = top_colors.get(meta, other_color)
        alpha = 0.95 if meta in top_meta else 0.42
        lw = 1.0 if meta in top_meta else 0.45
        ax_a.vlines(row["rank"], 0.9, row["n_alpha_sequences"], color=color, alpha=alpha, lw=lw)
        ax_a.scatter(
            row["rank"],
            row["n_alpha_sequences"],
            color=color,
            edgecolor="white" if meta in top_meta else "none",
            linewidth=0.35,
            s=24 if meta in top_meta else 8,
            alpha=alpha,
            zorder=3,
        )
    am001_row = ranked[ranked["meta_cluster_id"] == "AM001"].iloc[0]
    ax_a.text(
        am001_row["rank"] + 1.0,
        am001_row["n_alpha_sequences"],
        "AM001",
        va="center",
        fontsize=6.8,
        fontweight="bold",
        color="#222222",
    )
    ax_a.set_yscale("log")
    ax_a.set_yticks([1, 2, 5, 10, 20, 50, 100, 250])
    ax_a.set_yticklabels(["1", "2", "5", "10", "20", "50", "100", "250"])
    ax_a.set_xlabel("Meta-cluster rank")
    ax_a.set_ylabel("Unique pre-L2 Alpha sequences")
    ax_a.set_title("Six largest of 78 inferred Alpha meta-clusters")
    ax_a.grid(axis="y", color="#dddddd", lw=0.45, alpha=0.7)

    weekly = (
        membership.groupby(["collection_week", "plot_group"], as_index=False)["sequence_id"]
        .nunique()
        .rename(columns={"sequence_id": "n_sequences"})
    )
    groups = [f"Other {len(summary) - TOP_N} meta-clusters"] + list(reversed(top_meta))
    pivot = weekly.pivot(index="collection_week", columns="plot_group", values="n_sequences").fillna(0)
    for group in groups:
        if group not in pivot:
            pivot[group] = 0
    pivot = pivot[groups].sort_index()
    stack_colors = [other_color] + [top_colors[meta] for meta in reversed(top_meta)]
    ax_b.stackplot(pivot.index, [pivot[col].to_numpy() for col in pivot.columns], colors=stack_colors, alpha=0.92)
    ax_b.axvline(L2_START, color="#333333", lw=0.8, ls="--")
    ax_b.text(L2_START, ax_b.get_ylim()[1] * 0.95, "L2", ha="center", va="top", fontsize=6.5)
    ax_b.set_ylabel("Unique sequences")
    ax_b.set_title("Weekly pre-L2 Alpha burden by meta-cluster")
    format_month_axis(ax_b)
    ax_b.grid(axis="y", color="#dddddd", lw=0.45, alpha=0.7)

    cumulative = (
        membership.groupby(["collection_week", "meta_cluster_id"], as_index=False)["sequence_id"]
        .nunique()
        .rename(columns={"sequence_id": "n_sequences"})
    )
    total_weekly = cumulative.groupby("collection_week")["n_sequences"].sum().sort_index().cumsum()
    top6_weekly = (
        cumulative[cumulative["meta_cluster_id"].isin(top_meta)]
        .groupby("collection_week")["n_sequences"]
        .sum()
        .reindex(total_weekly.index, fill_value=0)
        .cumsum()
    )
    am001_weekly = (
        cumulative[cumulative["meta_cluster_id"] == "AM001"]
        .groupby("collection_week")["n_sequences"]
        .sum()
        .reindex(total_weekly.index, fill_value=0)
        .cumsum()
    )
    ax_c.plot(total_weekly.index, total_weekly, color="#222222", lw=1.5, label="All 78 meta-clusters")
    ax_c.plot(top6_weekly.index, top6_weekly, color="#4e79a7", lw=1.8, label="Top six")
    ax_c.plot(am001_weekly.index, am001_weekly, color="#b23a2e", lw=2.0, label="AM001")
    ax_c.axvline(L2_START, color="#333333", lw=0.8, ls="--")
    ax_c.set_ylabel("Cumulative unique sequences")
    ax_c.set_xlabel("Collection week")
    ax_c.set_title("Cumulative pre-L2 expansion")
    ax_c.legend(loc="upper left", frameon=False, fontsize=6.5)
    format_month_axis(ax_c)
    ax_c.grid(axis="y", color="#dddddd", lw=0.45, alpha=0.7)

    if trajectories.empty:
        ax_d.text(0.5, 0.5, "No signature trajectories available", transform=ax_d.transAxes, ha="center")
    else:
        signature_order = (
            impact[impact["meta_cluster_id"].isin(top_meta)]
            .sort_values(["meta_cluster_id", "q_value", "prevalence_difference"], ascending=[True, True, False])
            .drop_duplicates("meta_cluster_id")
        )
        signature_order["meta_cluster_id"] = pd.Categorical(
            signature_order["meta_cluster_id"],
            categories=top_meta,
            ordered=True,
        )
        signature_order = signature_order.sort_values("meta_cluster_id")
        plotted_labels: set[str] = set()
        for _, row in signature_order.iterrows():
            meta = str(row["meta_cluster_id"])
            mutation = row["mutation"]
            dat = trajectories[
                (trajectories["mutation"] == mutation)
                & (trajectories["collection_week"].between(pd.Timestamp("2020-11-01"), pd.Timestamp("2021-04-30")))
            ].sort_values("collection_week")
            if dat.empty:
                continue
            label = f"{meta}: {mutation}"
            if label in plotted_labels:
                continue
            plotted_labels.add(label)
            is_am001 = meta == "AM001"
            ax_d.plot(
                dat["collection_week"],
                dat["freq_alpha_sequences"],
                color=top_colors.get(meta, "#777777"),
                lw=2.1 if is_am001 else 1.35,
                alpha=0.98 if is_am001 else 0.78,
                label=label,
            )
        ax_d.axvline(L2_START, color="#333333", lw=0.8, ls="--")
        ax_d.set_ylim(-0.03, 1.03)
        ax_d.set_ylabel("Frequency among Alpha")
        ax_d.set_xlabel("Collection week")
        ax_d.set_title("Candidate signature trajectories after L2")
        ax_d.legend(loc="upper right", frameon=False, fontsize=5.7)
        format_month_axis(ax_d)
        ax_d.grid(axis="y", color="#dddddd", lw=0.45, alpha=0.7)

    legend_handles = [
        Line2D([0], [0], color=top_colors[meta], lw=2, label=meta_label(meta, sizes))
        for meta in top_meta
    ] + [
        Line2D([0], [0], color=other_color, lw=2, label=f"Other {len(summary) - TOP_N} meta-clusters")
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.52, 1.015),
        ncol=4,
        frameon=False,
        fontsize=6.2,
    )

    style.add_panel_labels([ax_a, ax_b, ax_c, ax_d], x=-0.14, y=1.08, size=9)
    save_all(fig, FIGURE_DIR / "supp_fig2_alpha_meta_cluster_amplification", height_in=6.2)


def proportional_table(
    df: pd.DataFrame,
    category_col: str,
    categories: list[str],
    *,
    meta_order: list[str],
) -> pd.DataFrame:
    tmp = df[["meta_cluster_id", "sequence_id", category_col]].copy()
    tmp[category_col] = tmp[category_col].fillna("Missing").astype(str)
    tmp[category_col] = tmp[category_col].where(tmp[category_col].isin(categories), "Other")
    categories = [category for category in categories if category in set(tmp[category_col])] + (
        ["Other"] if "Other" in set(tmp[category_col]) and "Other" not in categories else []
    )
    counts = (
        tmp.groupby(["meta_cluster_id", category_col], as_index=False)["sequence_id"]
        .nunique()
        .rename(columns={"sequence_id": "n_sequences", category_col: "category"})
    )
    counts = counts.pivot(index="meta_cluster_id", columns="category", values="n_sequences").fillna(0)
    for category in categories:
        if category not in counts:
            counts[category] = 0
    counts = counts.reindex(meta_order)[categories]
    return counts.div(counts.sum(axis=1), axis=0).fillna(0)


def plot_prop_stackedh(
    ax: plt.Axes,
    prop: pd.DataFrame,
    colors: dict[str, str | tuple[float, float, float]],
    *,
    title: str,
    sizes: pd.Series,
) -> None:
    y = np.arange(len(prop.index))
    left = np.zeros(len(prop.index))
    for category in prop.columns:
        vals = prop[category].to_numpy()
        ax.barh(y, vals, left=left, color=colors.get(category, "#bdbdbd"), label=category, height=0.72)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels([meta_label(meta, sizes) for meta in prop.index])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Proportion of unique sequences")
    ax.set_title(title)
    ax.xaxis.set_major_formatter(lambda x, _: f"{x:.0%}")
    ax.grid(axis="x", color="#dddddd", lw=0.45, alpha=0.75)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.18), ncol=2, frameon=False, fontsize=5.8)


def plot_alpha_top6_context(table_dir: Path) -> None:
    sequence_meta = read_meta_csv(
        table_dir,
        "alpha_top6_meta_cluster_sequence_metadata.csv",
        ["collection_date", "collection_week"],
    )
    summary = read_meta_csv(
        table_dir,
        "alpha_pre_l2_meta_cluster_summary.csv",
        ["first_collection_date", "last_collection_date"],
    ).sort_values("n_alpha_sequences", ascending=False)
    top_meta = summary.head(TOP_N)["meta_cluster_id"].tolist()
    sizes = summary.set_index("meta_cluster_id")["n_alpha_sequences"]

    sequence_meta = sequence_meta[sequence_meta["meta_cluster_id"].isin(top_meta)].copy()
    sequence_meta["health_board"] = sequence_meta["health_board"].fillna("Missing").astype(str)
    sequence_meta["age_broad"] = sequence_meta["age_broad"].fillna("Missing").astype(str)
    sequence_meta["simd_quintile"] = sequence_meta["simd_quintile"].fillna("Missing").astype(str)
    sequence_meta["test_reason_group"] = sequence_meta["test_reason_group"].fillna("Missing").astype(str)

    top_health_boards = (
        sequence_meta["health_board"]
        .value_counts()
        .loc[lambda s: s >= 8]
        .head(7)
        .index.tolist()
    )
    health_prop = proportional_table(sequence_meta, "health_board", top_health_boards, meta_order=top_meta)
    age_categories = ["0-17", "18-39", "40-64", "65+", "Missing"]
    age_prop = proportional_table(sequence_meta, "age_broad", age_categories, meta_order=top_meta)
    simd_categories = ["Q1 most deprived", "Q2", "Q3", "Q4", "Q5 least deprived", "Missing"]
    simd_prop = proportional_table(sequence_meta, "simd_quintile", simd_categories, meta_order=top_meta)
    reason_categories = [
        "Symptomatic",
        "Contact/isolation",
        "Confirmatory/repeat",
        "Outbreak/local request",
        "Other/unclear",
        "Missing",
    ]
    reason_prop = proportional_table(sequence_meta, "test_reason_group", reason_categories, meta_order=top_meta)

    fig, axes = style.new_figure(
        width="double",
        height_in=7.4,
        nrows=2,
        ncols=2,
        font_scale=0.76,
        gridspec_kw={"hspace": 1.02, "wspace": 0.34},
    )
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    hb_palette = dict(zip(health_prop.columns, sns.color_palette("tab10", len(health_prop.columns))))
    age_palette = {
        "0-17": "#7b3294",
        "18-39": "#008837",
        "40-64": "#80cdc1",
        "65+": "#c2a5cf",
        "Missing": "#d9d9d9",
    }
    simd_palette = {
        "Q1 most deprived": "#d7191c",
        "Q2": "#fdae61",
        "Q3": "#ffffbf",
        "Q4": "#abdda4",
        "Q5 least deprived": "#2b83ba",
        "Missing": "#d9d9d9",
    }
    reason_palette = {
        "Symptomatic": "#66c2a5",
        "Contact/isolation": "#fc8d62",
        "Confirmatory/repeat": "#8da0cb",
        "Outbreak/local request": "#e78ac3",
        "Other/unclear": "#a6d854",
        "Missing": "#bdbdbd",
    }

    plot_prop_stackedh(ax_a, health_prop, hb_palette, title="Health-board composition", sizes=sizes)
    plot_prop_stackedh(ax_b, age_prop, age_palette, title="Age composition", sizes=sizes)
    plot_prop_stackedh(ax_c, simd_prop, simd_palette, title="SIMD quintile composition", sizes=sizes)
    plot_prop_stackedh(ax_d, reason_prop, reason_palette, title="Testing-reason composition", sizes=sizes)

    style.add_panel_labels([ax_a, ax_b, ax_c, ax_d], x=-0.14, y=1.08, size=9)

    save_all(fig, FIGURE_DIR / "supp_fig3_alpha_top6_meta_cluster_context", height_in=7.4)


def main() -> None:
    table_dir = setup_environment()
    print(f"Reading Alpha meta-cluster tables from {table_dir}")
    plot_alpha_meta_cluster_amplification(table_dir)
    plot_alpha_top6_context(table_dir)
    print(f"Wrote supplementary Alpha meta-cluster figures to {FIGURE_DIR}")


if __name__ == "__main__":
    main()
