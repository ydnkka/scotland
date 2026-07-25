"""Build Chapter 4 Figure 2: window-level EpiLink cluster landscape."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    add_policy_bands,
    date_axis,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)

FIGURE_NAME = "fig_ch4_cluster_landscape"


def build(paths: Paths) -> None:
    window_coverage = read_table(paths, "window_coverage")
    cluster_window = read_table(paths, "cluster_window_summary")
    cluster_table = read_table(paths, "cluster_table")
    window_coverage["wn_mid_date"] = pd.to_datetime(
        window_coverage["wn_mid_date"], errors="coerce"
    )
    cluster_window = cluster_window.merge(
        window_coverage[["window_idx", "wn_mid_date", "wn_prop_sequenced"]],
        on="window_idx",
        how="left",
    )

    sizes = pd.to_numeric(cluster_table["cluster_size"], errors="coerce").dropna()
    sizes = sizes[sizes > 0].sort_values().to_numpy()
    ccdf_y = 1.0 - np.arange(len(sizes)) / len(sizes)

    fig, axes = styled_new_figure(
        width="double",
        height_in=6.7,
        nrows=2,
        ncols=2,
        constrained_layout=True,
    )
    ax = axes[0, 0]
    ax.plot(sizes, ccdf_y, color="#1f4e79", lw=1.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cluster size (all clusters)")
    ax.set_ylabel("Pr(cluster size >= x)")
    panel_label(ax, "A")

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
    panel_label(ax, "B")

    ax = axes[1, 0]
    ax.scatter(
        cluster_window["wn_prop_sequenced"],
        cluster_window["max_cluster_size"],
        s=22,
        color="#1f4e79",
        alpha=0.75,
        edgecolor="none",
    )
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_xlabel("Window sequencing coverage")
    ax.set_ylabel("Maximum cluster size")
    panel_label(ax, "C")

    ax = axes[1, 1]
    spread = cluster_table[["cluster_size", "cluster_n_datazones"]].copy()
    spread = spread.replace([np.inf, -np.inf], np.nan).dropna()
    spread = spread.loc[
        spread["cluster_size"].gt(1) & spread["cluster_n_datazones"].gt(0)
    ]
    x = np.log10(spread["cluster_size"].to_numpy())
    y = np.log10(spread["cluster_n_datazones"].to_numpy())
    hb = ax.hexbin(x, y, gridsize=42, mincnt=1, cmap="viridis", bins="log")
    fig.colorbar(hb, ax=ax, label="Clusters")
    ticks = [0, 1, 2, 3]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{10**tick:g}" for tick in ticks])
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{10**tick:g}" for tick in ticks])
    ax.set_xlabel("Non-singleton cluster size")
    ax.set_ylabel("Observed Data Zones")
    panel_label(ax, "D")
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
