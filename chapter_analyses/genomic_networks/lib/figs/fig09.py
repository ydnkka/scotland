"""Build the compatibility topology diagnostics figure."""

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

FIGURE_NAME = "fig_compatibility_topology"


def _date_values(values: pd.Series) -> np.ndarray:
    dates = pd.to_datetime(values, errors="coerce")
    return mdates.date2num(dates.dt.to_pydatetime())


def _numeric_values(values: pd.Series) -> np.ndarray:
    return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return np.nan
    return float(np.average(values.loc[mask], weights=weights.loc[mask]))


def weighted_mean_ci_from_se(
    values: pd.Series,
    weights: pd.Series,
    standard_errors: pd.Series,
) -> dict[str, float]:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return {
            "weighted_mean": np.nan,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }

    values = values.loc[mask].astype(float)
    weights = weights.loc[mask].astype(float)
    weighted_mean_value = float(np.average(values, weights=weights))

    se = standard_errors.loc[mask].astype(float)
    se_mask = se.notna() & np.isfinite(se)
    if not se_mask.any():
        return {
            "weighted_mean": weighted_mean_value,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }

    normalized = weights.loc[se_mask] / weights.sum()
    combined_se = float(np.sqrt(np.sum((normalized * se.loc[se_mask]) ** 2)))
    return {
        "weighted_mean": weighted_mean_value,
        "combined_se": combined_se,
        "ci_low": weighted_mean_value - 1.96 * combined_se,
        "ci_high": weighted_mean_value + 1.96 * combined_se,
    }


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
        strength_ci = weighted_mean_ci_from_se(
            group["strength_assortativity"],
            weights,
            group["strength_assortativity_se"],
        )
        row["strength_assortativity"] = strength_ci["weighted_mean"]
        row["strength_assortativity_se"] = strength_ci["combined_se"]
        row["strength_assortativity_ci_low"] = strength_ci["ci_low"]
        row["strength_assortativity_ci_high"] = strength_ci["ci_high"]
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
    x = _date_values(summary["wn_mid_date"])

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

    axes[0].plot(x, _numeric_values(summary["n_nodes"]), label="Nodes")
    axes[0].plot(x, _numeric_values(summary["n_edges_used"]), label="Edges")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Weighted mean")
    axes[0].legend(loc="upper left")

    axes[1].plot(x, _numeric_values(summary["edge_weight_total"]), color="#1f4e79")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Total edge weight")

    colors = ["#1f4e79", "#d95f02", "#1b9e77"]
    strength_ci_low = _numeric_values(summary["strength_assortativity_ci_low"])
    strength_ci_high = _numeric_values(summary["strength_assortativity_ci_high"])
    ribbon_mask = (
        np.isfinite(x) & np.isfinite(strength_ci_low) & np.isfinite(strength_ci_high)
    )
    axes[2].fill_between(
        x,
        strength_ci_low,
        strength_ci_high,
        where=ribbon_mask,
        color="#1b9e77",
        alpha=0.18,
        linewidth=0,
        zorder=1,
        label="Strength assortativity 95% CI",
    )
    for (metric, label), color in zip(metrics[3:], colors):
        axes[2].plot(
            x,
            _numeric_values(summary[metric]),
            label=label,
            color=color,
            zorder=2,
        )
    axes[2].axhline(0, color="#777777", lw=0.8, ls=":")
    axes[2].set_xlabel("Window midpoint date")
    axes[2].set_ylabel("Assortativity")
    axes[2].legend(
        loc="upper right",
        bbox_to_anchor=(0.95, 1.0),
        frameon=True,
        facecolor="white",
        framealpha=0.82,
        edgecolor="#dddddd",
    )
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
