"""Build Chapter 4 Supplementary Figure 2: compatibility topology diagnostics."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.dates as mdates
import numpy as np
import pandas as pd

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
    window_idx_from_id,
)

FIGURE_NAME = "fig_ch4_compatibility_topology"


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return np.nan
    return float(np.average(values.loc[mask], weights=weights.loc[mask]))


def build(paths: Paths) -> None:
    degree = read_table(paths, "compatibility_degree_assortativity_bootstrap")
    degree["window_idx"] = window_idx_from_id(degree["window_id"])
    degree = degree.loc[degree["n_edges_used"].gt(0)].copy()
    metrics = [
        ("n_nodes", "Nodes"),
        ("n_edges_used", "Edges"),
        ("edge_weight_total", "Total edge weight"),
        ("degree_assortativity", "Degree assortativity (equal edge weights)"),
        ("weighted_degree_assortativity", "Degree assortativity (EpiLink-weighted edges)"),
        ("strength_assortativity", "Strength assortativity (EpiLink-weighted edges)"),
    ]
    rows = []
    for window_idx, group in degree.groupby("window_idx"):
        weights = group["edge_weight_total"].clip(lower=0)
        row = {"window_idx": window_idx}
        for metric, _ in metrics:
            row[metric] = weighted_mean(group[metric], weights)
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values("window_idx")
    window_coverage = read_table(paths, "window_coverage")
    window_coverage["wn_mid_date"] = pd.to_datetime(
        window_coverage["wn_mid_date"], errors="coerce"
    )
    summary = summary.merge(
        window_coverage[["window_idx", "wn_mid_date", "policy_era"]],
        on="window_idx",
        how="left",
    )
    x = summary["wn_mid_date"]

    fig, axes = new_figure(
        width="double",
        height_in=7.0,
        nrows=3,
        ncols=1,
        sharex=True,
        constrained_layout=True,
    )
    for ax in axes:
        add_policy_bands(ax, window_coverage)

    axes[0].plot(x, summary["n_nodes"], label="Nodes")
    axes[0].plot(x, summary["n_edges_used"], label="Edges")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Weighted mean")
    axes[0].legend(loc="upper left")

    axes[1].plot(x, summary["edge_weight_total"], color="#1f4e79")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Total edge weight")

    colors = ["#1f4e79", "#d95f02", "#1b9e77"]
    for (metric, label), color in zip(metrics[3:], colors):
        axes[2].plot(x, summary[metric], label=label, color=color)
    axes[2].axhline(0, color="#777777", lw=0.8, ls=":")
    axes[2].set_xlabel("Window midpoint date")
    axes[2].set_ylabel("Assortativity")
    axes[2].legend(loc="upper right", bbox_to_anchor=(0.95, 1.0))
    date_axis(axes[2])
    axes[2].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    add_panel_labels(axes)
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
