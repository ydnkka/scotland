"""Build Chapter 4 Figure 2: window-level EpiLink cluster landscape."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import PercentFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    add_panel_labels,
    add_policy_bands,
    date_axis,
    new_figure,
    paths_from_args,
    read_table,
    styled_save_figure,
)

FIGURE_NAME = "fig_ch4_cluster_landscape"

WINDOW_CORRELATION_METRICS = (
    "wn_no_sequences",
    "wn_positive_tests",
    "wn_prop_sequenced",
)

CLUSTER_CORRELATION_METRICS = (
    "n_clusters",
    "n_non_singleton_clusters",
    "non_singleton_clusters_per_1000_sequences",
    "median_non_singleton_cluster_size",
    "p90_non_singleton_cluster_size",
    "max_non_singleton_cluster_size",
)

CLUSTER_SIZE_METRICS = (
    "median_non_singleton_cluster_size",
    "p90_non_singleton_cluster_size",
    "max_non_singleton_cluster_size",
)

CLUSTER_SPREAD_METRICS = (
    "q25_non_singleton_duration_days",
    "median_non_singleton_duration_days",
    "q75_non_singleton_duration_days",
    "median_non_singleton_datazones",
    "p90_non_singleton_datazones",
    "max_non_singleton_datazones",
    "q25_non_singleton_spatial_distance_km",
    "median_non_singleton_spatial_distance_km",
    "q75_non_singleton_spatial_distance_km",
)

PLOT_LABELS = {
    "wn_no_sequences": "Sequences",
    "wn_positive_tests": "Positive tests",
    "wn_prop_sequenced": "Seq. coverage",
    "n_clusters": "All clusters",
    "n_non_singleton_clusters": "Non-singleton\nclusters",
    "non_singleton_clusters_per_1000_sequences": "Non-singleton\nclusters /\n1,000 seq.",
    "median_non_singleton_cluster_size": "Median\ncluster size",
    "p90_non_singleton_cluster_size": "90th pct.\ncluster size",
    "max_non_singleton_cluster_size": "Maximum\ncluster size",
    "q25_non_singleton_duration_days": "25th pct.\nduration",
    "median_non_singleton_duration_days": "Median\nduration",
    "q75_non_singleton_duration_days": "75th pct.\nduration",
    "median_non_singleton_datazones": "Median\ndatazones",
    "p90_non_singleton_datazones": "90th pct.\ndatazones",
    "max_non_singleton_datazones": "Maximum\ndatazones",
    "q25_non_singleton_spatial_distance_km": "25th pct.\ndistance",
    "median_non_singleton_spatial_distance_km": "Median\ndistance",
    "q75_non_singleton_spatial_distance_km": "75th pct.\ndistance",
}


def _correlation_matrix(
    data: pd.DataFrame,
    row_vars: tuple[str, ...],
    col_vars: tuple[str, ...],
) -> pd.DataFrame:
    missing = [col for col in [*row_vars, *col_vars] if col not in data.columns]
    if missing:
        raise KeyError(f"Missing columns needed for correlation heatmap: {missing}")

    columns = list(dict.fromkeys([*row_vars, *col_vars]))
    numeric = data[columns].apply(pd.to_numeric, errors="coerce")
    numeric = numeric.replace([np.inf, -np.inf], np.nan)
    return numeric.corr(method="pearson").loc[list(row_vars), list(col_vars)]


def _draw_correlation_heatmap(
    ax,
    correlations: pd.DataFrame,
    *,
    xlabel: str,
    ylabel: str,
):
    values = correlations.to_numpy(dtype=float)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#f2f2f2")
    image = ax.imshow(
        values,
        aspect="auto",
        cmap=cmap,
        norm=TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0),
    )

    ax.set_xticks(np.arange(correlations.shape[1]))
    ax.set_xticklabels([PLOT_LABELS[col] for col in correlations.columns])
    ax.set_yticks(np.arange(correlations.shape[0]))
    ax.set_yticklabels([PLOT_LABELS[row] for row in correlations.index])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # x_rotation = 35 if correlations.shape[1] <= 6 else 55
    # x_fontsize = 6.5 if correlations.shape[1] <= 6 else 6.0
    plt.setp(
        ax.get_xticklabels(),
        rotation=90,
        ha="right",
        va="center",
        rotation_mode="anchor",
        # fontsize=x_fontsize,
    )
    ax.tick_params(axis="y")
    ax.tick_params(axis="both", length=0)
    ax.set_xticks(np.arange(-0.5, correlations.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, correlations.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    for row_idx in range(values.shape[0]):
        for col_idx in range(values.shape[1]):
            value = values[row_idx, col_idx]
            if not np.isfinite(value):
                continue
            ax.text(
                col_idx,
                row_idx,
                f"{value:.2f}",
                ha="center",
                va="center",
                color="white" if abs(value) >= 0.65 else "#1f1f1f",
                fontsize=6.5,
            )

    return image


def build(paths: Paths) -> None:
    window_coverage = read_table(paths, "window_coverage")
    cluster_window = read_table(paths, "cluster_window_summary")
    cluster_table = read_table(paths, "cluster_table")
    window_coverage["wn_mid_date"] = pd.to_datetime(
        window_coverage["wn_mid_date"], errors="coerce"
    )
    merge_cols = ["window_idx", "wn_mid_date"]
    merge_cols.extend(
        col
        for col in WINDOW_CORRELATION_METRICS
        if col not in cluster_window.columns and col in window_coverage.columns
    )
    cluster_window = cluster_window.merge(
        window_coverage[merge_cols],
        on="window_idx",
        how="left",
    )

    sizes = pd.to_numeric(cluster_table["cluster_size"], errors="coerce").dropna()
    sizes = sizes[sizes > 0].sort_values().to_numpy()
    ccdf_y = 1.0 - np.arange(len(sizes)) / len(sizes)

    fig, axes = new_figure(
        width="double",
        height_in=8,
        nrows=3,
        ncols=2,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [0.35, 0.35, 0.3]},
    )
    ax = axes[0, 0]
    ax.plot(sizes, ccdf_y, color="#1f4e79", lw=1.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cluster size (all clusters)")
    ax.set_ylabel("Pr(cluster size >= x)")

    ax = axes[0, 1]
    add_policy_bands(ax, window_coverage)
    ax.plot(
        cluster_window["wn_mid_date"],
        cluster_window["median_non_singleton_cluster_size"],
        label="Median",
        color="#1f4e79",
        lw=1.2,
    )
    ax.plot(
        cluster_window["wn_mid_date"],
        cluster_window["p90_non_singleton_cluster_size"],
        label="90th percentile",
        color="#d95f02",
        lw=1.2,
    )
    ax.plot(
        cluster_window["wn_mid_date"],
        cluster_window["max_non_singleton_cluster_size"],
        label="Maximum",
        color="#1b9e77",
        lw=1.1,
    )
    ax.set_yscale("log")
    ax.set_ylabel("Non-singleton cluster size")
    ax.legend(loc="upper left")
    date_axis(ax)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
    ax.tick_params(axis="x")

    ax = axes[1, 0]
    ax.scatter(
        cluster_window["wn_prop_sequenced"],
        cluster_window["max_non_singleton_cluster_size"],
        s=22,
        color="#1f4e79",
        alpha=0.75,
        edgecolor="none",
    )
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_xlabel("Window sequencing coverage")
    ax.set_ylabel("Maximum cluster size")

    ax = axes[1, 1]
    spread = cluster_table[["cluster_size", "median_pairwise_residential_distance_km"]].copy()
    spread = spread.replace([np.inf, -np.inf], np.nan).dropna()
    spread = spread.loc[
        spread["cluster_size"].gt(1) & spread["median_pairwise_residential_distance_km"].gt(0)
    ]
    x = np.log10(spread["cluster_size"].to_numpy())
    y = np.log10(spread["median_pairwise_residential_distance_km"].to_numpy())
    hb = ax.hexbin(x, y, gridsize=42, mincnt=1, cmap="viridis", bins="log")
    fig.colorbar(hb, ax=ax, label="Clusters")
    ticks = [0, 1, 2, 3]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{10**tick:g}" for tick in ticks])
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{10**tick:g}" for tick in ticks])
    ax.set_xlabel("Non-singleton cluster size")
    ax.set_ylabel("Median residential distance (km)")

    ax = axes[2, 0]
    window_cluster_corr = _correlation_matrix(
        cluster_window,
        WINDOW_CORRELATION_METRICS,
        CLUSTER_CORRELATION_METRICS,
    )
    heatmap = _draw_correlation_heatmap(
        ax,
        window_cluster_corr,
        xlabel="Cluster count and size summary",
        ylabel="Window denominator",
    )

    ax = axes[2, 1]
    cluster_spread_corr = _correlation_matrix(
        cluster_window,
        CLUSTER_SIZE_METRICS,
        CLUSTER_SPREAD_METRICS,
    )
    _draw_correlation_heatmap(
        ax,
        cluster_spread_corr,
        xlabel="Spatial and temporal summary",
        ylabel="Cluster size summary",
    )
    cbar = fig.colorbar(heatmap, ax=axes[2, :], shrink=0.88, pad=0.02)
    cbar.set_label("Pearson correlation coefficient")

    add_panel_labels(axes.ravel())
    styled_save_figure(fig, paths, FIGURE_NAME)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    print(f"Wrote {FIGURE_NAME} to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
