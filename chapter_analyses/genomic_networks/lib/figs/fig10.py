"""Build the topology-denominator correlation heatmap."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    new_figure,
    paths_from_args,
    read_table,
    styled_save_figure,
    window_idx_from_id,
)

FIGURE_NAME = "fig_compatibility_topology_correlations"

WINDOW_DENOMINATOR_METRICS = (
    "wn_no_sequences",
    "wn_positive_tests",
    "wn_prop_sequenced",
)

TOPOLOGY_METRICS = (
    "n_nodes",
    "n_edges_used",
    "edge_weight_total",
    "mean_degree",
    "max_degree",
    "mean_strength",
    "max_strength",
    "degree_assortativity",
    "weighted_degree_assortativity",
    "strength_assortativity",
)

PLOT_LABELS = {
    "wn_no_sequences": "Sequences",
    "wn_positive_tests": "Positive tests",
    "wn_prop_sequenced": "Seq. coverage",
    "n_nodes": "Nodes",
    "n_edges_used": "Edges",
    "edge_weight_total": "Total edge weight",
    "mean_degree": "Mean degree",
    "max_degree": "Max. degree",
    "mean_strength": "Mean strength",
    "max_strength": "Max. strength",
    "degree_assortativity": "Degree assort.\n(equal edge weights)",
    "weighted_degree_assortativity": "Degree assort.\n(EpiLink-weighted edges)",
    "strength_assortativity": "Strength assort.\n(EpiLink-weighted edges)",
}


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return np.nan
    return float(np.average(values.loc[mask], weights=weights.loc[mask]))


def _build_topology_window_summary(paths: Paths) -> pd.DataFrame:
    degree = read_table(paths, "compatibility_degree_assortativity_bootstrap")
    degree["window_idx"] = window_idx_from_id(degree["window_id"])
    degree = degree.loc[degree["n_edges_used"].gt(0)].copy()

    rows = []
    for window_idx, group in degree.groupby("window_idx"):
        weights = group["edge_weight_total"].clip(lower=0)
        row = {"window_idx": window_idx}
        for metric in TOPOLOGY_METRICS:
            row[metric] = weighted_mean(group[metric], weights)
        rows.append(row)

    summary = pd.DataFrame(rows).sort_values("window_idx")
    window_coverage = read_table(paths, "window_coverage")
    return summary.merge(
        window_coverage[
            [
                "window_idx",
                *WINDOW_DENOMINATOR_METRICS,
            ]
        ],
        on="window_idx",
        how="left",
    )


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


def _draw_correlation_heatmap(ax, correlations: pd.DataFrame):
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
    ax.set_xlabel("Compatibility-topology summary")
    ax.set_ylabel("Window denominator")
    plt.setp(
        ax.get_xticklabels(),
        rotation=90,
        ha="right",
        va="center",
        rotation_mode="anchor",
    )
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
    summary = _build_topology_window_summary(paths)
    correlations = _correlation_matrix(
        summary,
        WINDOW_DENOMINATOR_METRICS,
        TOPOLOGY_METRICS,
    )

    fig, ax = new_figure(
        width="double",
        height_in=3.8,
        constrained_layout=True,
    )
    heatmap = _draw_correlation_heatmap(ax, correlations)
    cbar = fig.colorbar(heatmap, ax=ax, shrink=0.9, pad=0.02)
    cbar.set_label("Pearson correlation coefficient")
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
