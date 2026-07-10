"""Build Chapter 4 Figure 4: selected compatibility mixing matrices."""

from __future__ import annotations

from pathlib import Path
import argparse
import re
import sys
from typing import Any

from matplotlib.image import AxesImage
from matplotlib.axes import Axes
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    new_blank_figure,
    panel_label,
    paths_from_args,
    read_table,
    save_figure,
)


def display_matrix_category(attribute: str, value: Any) -> str:
    if attribute == "simd_quintile":
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return str(value)
        if numeric.is_integer():
            return str(int(numeric))
    return str(value)


def matrix_order(
    categories: list[str],
    *,
    attribute: str,
    grouped: pd.DataFrame | None = None,
) -> list[str]:
    if attribute == "age_group":
        preferred = ["00-04", "05-14", "15-24", "25-64", "65-74", "75+"]
        return [value for value in preferred if value in categories]
    if attribute == "urban_rural":
        preferred = [
            "Large Urban Areas",
            "Other Urban Areas",
            "Accessible Small Towns",
            "Remote Small Towns",
            "Accessible Rural",
            "Remote Rural",
        ]
        return [value for value in preferred if value in categories]
    if attribute == "health_board" and grouped is not None:
        row_totals = grouped.groupby("source_category", dropna=False)[
            "edge_weight"
        ].sum()
        col_totals = grouped.groupby("target_category", dropna=False)[
            "edge_weight"
        ].sum()
        totals = row_totals.add(col_totals, fill_value=0.0)
        return totals.sort_values(ascending=False).index.astype(str).tolist()

    def key(value: str) -> tuple[int, object]:
        text = str(value)
        match = re.search(r"\d+", text)
        if match:
            return (0, int(match.group()))
        return (1, text)

    return sorted(categories, key=key)


def aggregate_row_matrix(
    matrix_table: pd.DataFrame,
    attribute: str,
) -> tuple[pd.DataFrame, str]:
    work = matrix_table.loc[matrix_table["attribute"].eq(attribute)].copy()
    label = work["attribute_label"].dropna().iloc[0] if not work.empty else attribute
    work["source_category"] = work["source_category"].map(
        lambda value: display_matrix_category(attribute, value)
    )
    work["target_category"] = work["target_category"].map(
        lambda value: display_matrix_category(attribute, value)
    )
    grouped = (
        work.groupby(["source_category", "target_category"], dropna=False)[
            "edge_weight"
        ]
        .sum()
        .reset_index()
    )
    cats = matrix_order(
        sorted(
            set(grouped["source_category"].astype(str))
            | set(grouped["target_category"].astype(str))
        ),
        attribute=attribute,
        grouped=grouped,
    )
    matrix = grouped.pivot_table(
        index="source_category",
        columns="target_category",
        values="edge_weight",
        aggfunc="sum",
        fill_value=0.0,
    ).reindex(index=cats, columns=cats, fill_value=0.0)
    row_totals = matrix.sum(axis=1).replace(0, np.nan)
    row_share = matrix.div(row_totals, axis=0)
    return row_share, str(label)


def draw_matrix_heatmap(
    ax: Axes,
    matrix: pd.DataFrame,
    *,
    title: str,
    vmax: float | None = None,
    label_size: float = 7.0,
) -> AxesImage:
    values = matrix.to_numpy(dtype=float)
    vmax = vmax or np.nanquantile(values, 0.98)
    image = ax.imshow(values, cmap="Blues", vmin=0, vmax=vmax, aspect="auto")
    ax.set_title(title)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=label_size)
    ax.set_yticklabels(matrix.index, fontsize=label_size)
    return image


def build(paths: Paths) -> None:
    compatibility = read_table(paths, "compatibility_mixing_matrix")

    panels = [
        ("age_group", "Age group", "A", (0, 0), 7.0),
        ("simd_quintile", "SIMD quintile", "B", (1, 0), 7.0),
        ("urban_rural", "Urban/rural class", "C", (2, 0), 6.7),
        ("health_board", "Health board", "D", (slice(None), 1), 6.0),
    ]
    matrices = [aggregate_row_matrix(compatibility, attr)[0] for attr, *_ in panels]
    vmax = max(np.nanquantile(matrix.to_numpy(), 0.98) for matrix in matrices)

    fig = new_blank_figure(
        "double",
        width_in=9.4,
        height_in=8.2,
        constrained_layout=True,
    )
    grid = fig.add_gridspec(3, 3, width_ratios=[1.0, 1.65, 0.045])
    colorbar_axis = fig.add_subplot(grid[:, 2])
    for attr, title, label, grid_position, label_size in panels:
        ax = fig.add_subplot(grid[grid_position])
        matrix, _ = aggregate_row_matrix(compatibility, attr)
        image = draw_matrix_heatmap(
            ax,
            matrix,
            title=title,
            vmax=vmax,
            label_size=label_size,
        )
        panel_label(ax, label)
    cbar = fig.colorbar(image, cax=colorbar_axis)
    cbar.set_label("Row share of weighted edge mass")
    fig.supxlabel("Linked endpoint category")
    fig.supylabel("Reference endpoint category")
    save_figure(fig, paths, "fig04_mixing_matrices", tight=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    print(f"Wrote fig04_mixing_matrices to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
