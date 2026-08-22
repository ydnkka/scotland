"""Build the detector-threshold robustness figure."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from ..sse.io import write_table
from ..sse.scoring import add_sse_node_metrics
from .common import (
    HIGH_PRIORITY,
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
    read_table,
    styled_save_figure,
)

FILE_NAME = "threshold_robustness"
MIN_SIZE_VALUES = tuple(range(2, 21, 2))
ALPHA_VALUES = tuple(i / 100 for i in range(1, 11))
PRIMARY_MIN_SIZE = 6
PRIMARY_ALPHA = 0.05
CELL_FONT_SIZE = 5.2
COLORBAR_LABEL_SIZE = 6.5


def _contrast_text_color(image, value: float) -> str:
    """Choose black or white text against a heatmap cell."""
    rgba = image.cmap(image.norm(value))
    r, g, b = rgba[:3]
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "black" if luminance > 0.58 else "white"


def build_threshold_robustness(nodes: pd.DataFrame) -> pd.DataFrame:
    baseline = set(
        nodes.loc[nodes["candidate_tier"].isin(HIGH_PRIORITY), "cluster_id"]
    )
    rows = []
    for min_size in MIN_SIZE_VALUES:
        # Re-score in memory so threshold sensitivity does not alter primary outputs.
        scored = add_sse_node_metrics(nodes, min_cluster_size=min_size)
        tested = scored.loc[scored["sse_tested"].fillna(False)].copy()
        size_ok = tested["cluster_size"].ge(min_size)
        for alpha in ALPHA_VALUES:
            selected = tested.loc[
                size_ok
                & (
                    tested["burst_score_upper_p"].le(alpha)
                    | (
                        tested["burden_eligible"].fillna(False)
                        & tested["burden_score_upper_p"].le(alpha)
                    )
                ),
                "cluster_id",
            ]
            chosen = set(selected)
            union = chosen | baseline
            rows.append(
                {
                    "min_cluster_size": min_size,
                    "alpha": alpha,
                    "eligible_n": int(size_ok.sum()),
                    "candidate_n": len(chosen),
                    "baseline_overlap_n": len(chosen & baseline),
                    "baseline_jaccard": len(chosen & baseline) / len(union)
                    if union
                    else 1.0,
                }
            )
    return pd.DataFrame(rows)


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_threshold_robustness(nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")
    sizes = sorted(table["min_cluster_size"].unique())
    alphas = sorted(table["alpha"].unique())
    fig, axes = new_figure(
        width="double",
        height_in=5.2,
        nrows=2,
        ncols=2,
        sharex="col",
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.55, 1.0]},
    )
    heatmap_specs = (
        (
            axes[0, 0],
            "candidate_n",
            "Candidates retained",
            "Blues",
            None,
            "Count",
        ),
        (
            axes[1, 0],
            "baseline_jaccard",
            "Agreement with primary",
            "Greens",
            1,
            "Jaccard",
        ),
    )
    alpha_ticks = [0.01, 0.05, 0.10]
    alpha_tick_pos = [alphas.index(alpha) for alpha in alpha_ticks if alpha in alphas]
    alpha_tick_labels = [f"{alpha:.2f}" for alpha in alpha_ticks if alpha in alphas]
    size_tick_pos = list(range(len(sizes)))
    size_tick_labels = [str(size) for size in sizes]
    primary_x = alphas.index(PRIMARY_ALPHA)
    primary_y = sizes.index(PRIMARY_MIN_SIZE)

    for row_idx, (ax, value, title, cmap, vmax, cbar_label) in enumerate(
        heatmap_specs
    ):
        matrix = table.pivot(
            index="min_cluster_size", columns="alpha", values=value
        ).reindex(index=sizes, columns=alphas)
        image = ax.imshow(
            matrix,
            aspect="auto",
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            interpolation="nearest",
            origin="lower",
        )
        ax.add_patch(
            Rectangle(
                (primary_x - 0.5, primary_y - 0.5),
                1,
                1,
                fill=False,
                edgecolor="red",
                linewidth=1.2,
            )
        )
        ax.set_xticks(alpha_tick_pos, alpha_tick_labels)
        ax.set_yticks(size_tick_pos, size_tick_labels)
        ax.set_xticks(np.arange(-0.5, len(alphas), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(sizes), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.6)
        ax.tick_params(which="both", length=0)
        ax.tick_params(which="minor", bottom=False, left=False)
        if row_idx == 0:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("Upper-tail p threshold")
        ax.set_ylabel("Minimum cluster size")
        ax.set_title(title)
        for y, row in enumerate(matrix.to_numpy()):
            for x, cell in enumerate(row):
                if pd.isna(cell):
                    continue
                label = f"{cell:.2f}" if value == "baseline_jaccard" else f"{cell:.0f}"
                ax.text(
                    x,
                    y,
                    label,
                    ha="center",
                    va="center",
                    fontsize=CELL_FONT_SIZE,
                    color=_contrast_text_color(image, float(cell)),
                )
        cbar = fig.colorbar(
            image,
            ax=ax,
            fraction=0.024,
            pad=0.018,
            shrink=0.9,
            aspect=16,
        )
        cbar.set_label(cbar_label, fontsize=COLORBAR_LABEL_SIZE, labelpad=2)
        cbar.ax.tick_params(length=2, pad=2, labelsize=COLORBAR_LABEL_SIZE)

    primary = table.loc[np.isclose(table["alpha"], PRIMARY_ALPHA)].sort_values(
        "min_cluster_size"
    )
    line_specs = (
        (
            axes[0, 1],
            "candidate_n",
            f"Candidates at p <= {PRIMARY_ALPHA:.2f}",
            "Candidates retained",
        ),
        (
            axes[1, 1],
            "baseline_jaccard",
            f"Jaccard at p <= {PRIMARY_ALPHA:.2f}",
            "Jaccard similarity",
        ),
    )
    for row_idx, (ax, value, title, ylabel) in enumerate(line_specs):
        ax.plot(
            primary["min_cluster_size"],
            primary[value],
            marker="o",
            markersize=3,
            linewidth=1.2,
            color="black",
        )
        primary_row = primary.loc[primary["min_cluster_size"].eq(PRIMARY_MIN_SIZE)]
        if not primary_row.empty:
            ax.scatter(
                primary_row["min_cluster_size"],
                primary_row[value],
                s=22,
                color="red",
                zorder=3,
            )
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlim(min(sizes) - 0.5, max(sizes) + 0.5)
        ax.set_xticks(sizes)
        ax.grid(axis="y", color="0.9", linewidth=0.6)
        if value == "baseline_jaccard":
            ax.set_ylim(0, 1.05)
        else:
            ax.set_ylim(bottom=0)
        if row_idx == 0:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("Minimum cluster size")

    add_panel_labels(axes.ravel())
    outputs = styled_save_figure(fig, paths, f"fig_{FILE_NAME}")
    return {"figure": fig, "outputs": outputs}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    paths = paths_from_args(parser.parse_args())
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
