"""Build Chapter 4 Supplementary Figure 2: compatibility topology diagnostics."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
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
    degree = read_table(paths, "compatibility_degree_assortativity")
    degree["window_idx"] = window_idx_from_id(degree["window_id"])
    degree = degree.loc[degree["n_edges_used"].gt(0)].copy()
    metrics = [
        ("n_nodes", "Nodes"),
        ("n_edges_used", "Edges"),
        ("edge_weight_total", "Total edge weight"),
        ("degree_assortativity", "Degree assortativity"),
        ("weighted_degree_assortativity", "Weighted degree assortativity"),
        ("strength_assortativity", "Strength assortativity"),
    ]
    rows = []
    for window_idx, group in degree.groupby("window_idx"):
        weights = group["edge_weight_total"].clip(lower=0)
        row = {"window_idx": window_idx}
        for metric, _ in metrics:
            row[metric] = weighted_mean(group[metric], weights)
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values("window_idx")

    fig, axes = styled_new_figure(
        width="double",
        height_in=7.0,
        nrows=3,
        ncols=1,
        sharex=True,
    )
    axes[0].plot(summary["window_idx"], summary["n_nodes"], label="Nodes")
    axes[0].plot(summary["window_idx"], summary["n_edges_used"], label="Edges")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Weighted mean")
    axes[0].legend(loc="upper left")
    panel_label(axes[0], "A")

    axes[1].plot(summary["window_idx"], summary["edge_weight_total"], color="#1f4e79")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Total edge weight")
    panel_label(axes[1], "B")

    colors = ["#1f4e79", "#d95f02", "#1b9e77"]
    for (metric, label), color in zip(metrics[3:], colors):
        axes[2].plot(summary["window_idx"], summary[metric], label=label, color=color)
    axes[2].axhline(0, color="#777777", lw=0.8, ls=":")
    axes[2].set_xlabel("Window")
    axes[2].set_ylabel("Assortativity")
    axes[2].legend(loc="upper right")
    panel_label(axes[2], "C")
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
