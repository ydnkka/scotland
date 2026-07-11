"""Build Chapter 4 Supplementary Figure 3: SIMD population-weighting validation."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)


FIGURE_NAME = "fig_ch4_simd_population_weighting"


def build(paths: Paths) -> None:
    group_summary = read_table(paths, "simd_population_weighting_group_summary")
    movement = read_table(paths, "simd_population_weighting_movement")
    movement = movement.loc[movement["comparison_method"].eq("equal_datazone")]

    fig, axes = styled_new_figure(
        width="double", height_in=3.5, nrows=1, ncols=2
    )
    ax = axes[0]
    for method, group in group_summary.groupby("grouping_method_label"):
        group = group.sort_values("simd_group")
        ax.plot(
            group["simd_group"],
            group["pct_population"],
            marker="o",
            lw=1.2,
            label=method,
        )
    ax.axhline(20, color="#777777", lw=0.8, ls=":")
    ax.set_xlabel("SIMD quintile")
    ax.set_ylabel("Population share (%)")
    ax.legend(loc="best")
    panel_label(ax, "A")

    matrix = (
        movement.pivot_table(
            index="comparison_group",
            columns="population_weighted_group",
            values="pct_population",
            fill_value=0.0,
        )
        .sort_index()
        .sort_index(axis=1)
    )
    image = axes[1].imshow(matrix.to_numpy(), cmap="Blues", vmin=0)
    axes[1].set_xticks(np.arange(matrix.shape[1]))
    axes[1].set_xticklabels(matrix.columns)
    axes[1].set_yticks(np.arange(matrix.shape[0]))
    axes[1].set_yticklabels(matrix.index)
    axes[1].set_xlabel("Population-weighted group")
    axes[1].set_ylabel("Equal-Data-Zone group")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i, j]
            if value > 0: # type: ignore
                axes[1].text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=7)
    fig.colorbar(image, ax=axes[1], label="Population share (%)")
    panel_label(axes[1], "B")
    styled_save_figure(fig, paths, FIGURE_NAME, tight=False)


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
