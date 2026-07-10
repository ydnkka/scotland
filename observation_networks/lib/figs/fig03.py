"""Build Chapter 4 Figure 3: weighted assortativity baseline."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

from matplotlib.image import AxesImage
from matplotlib.axes import Axes
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    ATTRIBUTE_ORDER,
    Paths,
    add_common_args,
    new_blank_figure,
    panel_label,
    paths_from_args,
    read_table,
    save_figure,
    window_idx_from_id,
)


def weighted_mean_ci_from_se(
    values: pd.Series,
    weights: pd.Series,
    standard_errors: pd.Series,
) -> dict[str, float]:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return {"weighted_mean": np.nan, "combined_se": np.nan}

    values = values.loc[mask].astype(float)
    weights = weights.loc[mask].astype(float)
    weighted_mean_value = float(np.average(values, weights=weights))

    se_mask = standard_errors.loc[mask].notna()
    if not se_mask.any():
        return {"weighted_mean": weighted_mean_value, "combined_se": np.nan}

    ci_weights = weights.loc[se_mask]
    ci_standard_errors = standard_errors.loc[mask].loc[se_mask].astype(float)
    normalized = ci_weights / weights.sum()
    combined_se = float(np.sqrt(np.sum((normalized * ci_standard_errors) ** 2)))
    return {"weighted_mean": weighted_mean_value, "combined_se": combined_se}


def compatibility_window_assortativity(paths: Paths) -> pd.DataFrame:
    assort = read_table(paths, "compatibility_assortativity")
    assort["window_idx"] = window_idx_from_id(assort["window_id"])
    work = assort.loc[
        assort["assortativity"].notna()
        & assort["edge_weight_total"].gt(0)
        & assort["n_categories"].gt(1)
        & assort["n_edge_contributions_used"].ge(20)
    ].copy()

    rows = []
    for (window_idx, attribute, label), group in work.groupby(
        ["window_idx", "attribute", "attribute_label"], dropna=False
    ):
        ci = weighted_mean_ci_from_se(
            group["assortativity"],
            group["edge_weight_total"],
            group["assortativity_se"],
        )
        rows.append(
            {
                "window_idx": window_idx,
                "attribute": attribute,
                "attribute_label": label,
                "assortativity": ci["weighted_mean"],
                "assortativity_se": ci["combined_se"],
                "edge_weight_total": group["edge_weight_total"].sum(),
            }
        )
    return pd.DataFrame(rows)


def transition_window_assortativity(paths: Paths) -> pd.DataFrame:
    assort = read_table(paths, "transition_assortativity")
    assort["window_idx"] = window_idx_from_id(assort["source_window_id"])
    return assort.loc[
        assort["assortativity"].notna()
        & assort["edge_weight_total"].gt(0)
        & assort["n_categories"].gt(1)
    ].copy()


def ordered_attribute_pivot(
    df: pd.DataFrame, value_col: str = "assortativity"
) -> pd.DataFrame:
    pivot = df.pivot_table(
        index="attribute_label",
        columns="window_idx",
        values=value_col,
        aggfunc="mean",
    )
    rows = [label for label in ATTRIBUTE_ORDER if label in pivot.index]
    rows.extend([label for label in pivot.index if label not in rows])
    return pivot.loc[rows].sort_index(axis=1)


def heatmap_with_window_ticks(
    ax: Axes,
    pivot: pd.DataFrame,
    *,
    title: str,
    vmin: float = -0.6,
    vmax: float = 0.6,
) -> AxesImage:
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    image = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="RdBu_r", norm=norm)
    ax.set_title(title)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    columns = pivot.columns.to_numpy()
    if len(columns) > 0:
        tick_positions = np.linspace(0, len(columns) - 1, min(8, len(columns))).round()
        tick_positions = np.unique(tick_positions.astype(int))
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([f"W{int(columns[pos]):03d}" for pos in tick_positions])
    ax.set_xlabel("Window")
    return image


def build(paths: Paths) -> None:
    comp = compatibility_window_assortativity(paths)
    trans = transition_window_assortativity(paths)
    comp_pivot = ordered_attribute_pivot(comp)
    trans_pivot = ordered_attribute_pivot(trans)

    fig = new_blank_figure(
        "double",
        width_in=8.7,
        height_in=5.8,
        constrained_layout=True,
    )
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.0, 0.035],
        height_ratios=[1.15, 1.0],
    )
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[1, 0])]
    colorbar_axis = fig.add_subplot(grid[:, 1])
    image = heatmap_with_window_ticks(
        axes[0],
        comp_pivot,
        title="Compatibility graph, edge-weighted over Pango lineage networks",
        vmax=1.0,
    )
    panel_label(axes[0], "A")
    heatmap_with_window_ticks(
        axes[1],
        trans_pivot,
        title="Temporal transition graph",
        vmax=1.0,
    )
    panel_label(axes[1], "B")
    cbar = fig.colorbar(image, cax=colorbar_axis)
    cbar.set_label("Assortativity")
    save_figure(fig, paths, "fig03_assortativity_baseline", tight=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    print(f"Wrote fig03_assortativity_baseline to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
