"""Shared size-profile aggregates, category labels, and plotting helpers."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from utils import load_analysis_columns
from utils.style import FIG_WIDTHS_IN, set_theme

from ..characterisation import (
    ATTRIBUTES,
    ENTROPY_ATTRIBUTES,
    build_summaries,
    validate_summaries,
)
from ..sse.config import TRANSITION_WINDOW_STRIDE
from ..sse.io import write_table
from .common import Paths, add_panel_labels, new_figure, read_table, styled_save_figure

COLORS = {"CSC": "#0072B2", "Background": "#D55E00"}
SCENARIO_ORDER = ("p_0.010", "p_0.025", "primary", "p_0.100", "size_10", "size_20")
SCENARIO_LABELS = ("p≤0.01", "p≤0.025", "Primary", "p≤0.10", "Size≥10", "Size≥20")
CATEGORY_ORDERS = {
    "sex": ["Female", "Male"],
    "age": ["00-04", "05-14", "15-24", "25-64", "65-74", "75+"],
    "simd": ["1", "2", "3", "4", "5"],
    "urban_rural": [
        "Large Urban Areas",
        "Other Urban Areas",
        "Accessible Small Towns",
        "Remote Small Towns",
        "Accessible Rural",
        "Remote Rural",
    ],
}
SHORT_LABELS = {
    "Large Urban Areas": "Large urban",
    "Other Urban Areas": "Other urban",
    "Accessible Small Towns": "Accessible towns",
    "Remote Small Towns": "Remote towns",
    "Accessible Rural": "Accessible rural",
    "Remote Rural": "Remote rural",
    "Greater Glasgow and Clyde": "Glasgow and Clyde",
    "Dumfries and Galloway": "Dumfries and Galloway",
}


def prepare(paths: Paths) -> dict[str, pd.DataFrame]:
    nodes = read_table(paths, "cluster_table")
    columns = [
        "cluster_id",
        "window_id",
        "sequence_id",
        *[column for _, column, _ in ATTRIBUTES],
    ]
    records = load_analysis_columns(columns, window_stride=TRANSITION_WINDOW_STRIDE)
    tables = build_summaries(nodes, records)
    validate_summaries(tables, scotland_baseline=True)
    for key, table in tables.items():
        write_table(table, paths.result_table_dir, f"tab_size_profile_{key}")
    return tables


def load_saved_summaries(paths: Paths) -> dict[str, pd.DataFrame]:
    """Load existing aggregates for presentation-only rendering."""
    names = (
        "size_summary",
        "size_cumulative",
        "composition",
        "composition_contrasts",
        "entropy",
    )
    tables = {}
    for name in names:
        try:
            tables[name] = read_table(paths, f"tab_size_profile_{name}")
        except FileNotFoundError:
            # Existing aggregates remain usable after the output-name migration.
            tables[name] = read_table(paths, f"tab_ch6_{name}")
    validate_summaries(tables, scotland_baseline=True)
    return tables


def _finish(fig, paths, name):
    outputs = styled_save_figure(fig, paths, name)
    plt.close(fig)
    return outputs


def _heatmap(fig, ax, matrix, limit, title, *, annotate=False):
    cmap = plt.get_cmap("RdBu").copy()
    cmap.set_bad("#eeeeee")
    mesh = ax.imshow(
        matrix.to_numpy(float),
        cmap=cmap,
        norm=TwoSlopeNorm(0, -limit, limit),
        aspect="auto",
    )
    ax.set_yticks(
        range(len(matrix)),
        [SHORT_LABELS.get(str(x), str(x)) for x in matrix.index],
        fontsize=7,
    )
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, fontsize=8)
    ax.set_title(title, loc="left")
    ax.tick_params(length=0)
    if annotate:
        for row in range(len(matrix)):
            for col in range(len(matrix.columns)):
                value = matrix.iloc[row, col]
                ax.text(
                    col,
                    row,
                    f"{value:+.1f}" if pd.notna(value) else "–",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if abs(value) > 0.65 * limit else "#222222",
                )
    return mesh


def draw_composition(data, *, xcol, order, labels, rotate_labels=False):
    set_theme()
    fig = plt.figure(figsize=(FIG_WIDTHS_IN["double"], 7.0), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 1.8])
    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1:]),
    ]
    limit = max(1, np.ceil(data.difference_pp.abs().max()))
    for (key, _, label), ax in zip(ATTRIBUTES, axes):
        block = data.loc[data.attribute.eq(key)]
        matrix = block.pivot(
            index="category", columns=xcol, values="difference_pp"
        ).reindex(columns=order)
        matrix = matrix.reindex(CATEGORY_ORDERS.get(key, sorted(matrix.index)))
        matrix.columns = labels
        if key == "simd":
            matrix.index = ["Q" + x for x in matrix.index]
        mesh = _heatmap(fig, ax, matrix, limit, label, annotate=True)
        if rotate_labels:
            ax.tick_params(axis="x", labelrotation=45)
    add_panel_labels(axes)
    fig.colorbar(
        mesh,
        ax=axes,
        location="bottom",
        shrink=0.7,
        pad=0.025,
        aspect=40,
        label="Candidate minus background (percentage points)",
    )
    return fig


def draw_entropy(data, *, xcol, order, labels, sensitivity=False):
    fig, axes = new_figure(width="double", nrows=6, ncols=2, height_in=9.2)
    for row, (key, _, label) in enumerate(ENTROPY_ATTRIBUTES):
        for col, scale in enumerate(("obs", "z")):
            ax = axes[row, col]
            for offset, group in [(-0.12, "Background"), (0.12, "CSC")]:
                block = (
                    data.loc[
                        data.attribute.eq(key)
                        & data.scale.eq(scale)
                        & data.group.eq(group)
                    ]
                    .set_index(xcol)
                    .reindex(order)
                )
                x = np.arange(len(order)) + offset
                ax.errorbar(
                    x,
                    block["median"],
                    yerr=np.vstack(
                        [block["median"] - block.q25, block.q75 - block["median"]]
                    ),
                    fmt="o",
                    markersize=3,
                    capsize=2,
                    lw=1,
                    color=COLORS[group],
                    label=group,
                )
            ax.set_title(
                label + (" · observed" if scale == "obs" else " · null-standardised"),
                fontsize=8,
                loc="left",
            )
            ax.set_xticks(
                range(len(order)), labels if row == 5 else [""] * len(order), fontsize=7
            )
            if sensitivity and row == 5:
                ax.tick_params(axis="x", labelrotation=40)
            if scale == "obs":
                ax.set_ylim(-0.03, 1.03)
            else:
                ax.axhline(0, color="#777777", ls=":", lw=0.7)
            ax.grid(axis="y", alpha=0.18)
            ax.tick_params(axis="y", labelsize=9)
    add_panel_labels(axes)
    if not sensitivity:
        handles, legend_labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(
            handles,
            legend_labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.005),
            fontsize=9,
            ncol=2,
            frameon=False,
        )
        fig.tight_layout(h_pad=0.8, rect=(0, 0.04, 1, 1))
    else:
        axes[0, 0].legend(loc="lower left", fontsize=8, ncol=2)
        fig.tight_layout(h_pad=0.8)
    return fig
