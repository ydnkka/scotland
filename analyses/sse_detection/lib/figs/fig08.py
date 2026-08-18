"""Build the transition-graph diagnostics and role counts figure."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle

from ..sse.io import write_table
from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    add_policy_bands,
    date_axis,
    new_figure,
    paths_from_args,
    read_table,
    styled_save_figure,
)
from .fig01 import CLUSTER_ROLE_GROUPS, ROLE_COLORS

FILE_NAME = "transition_graph_characteristics"

PREFERRED_ROLES = [
    "Isolated",
    "Source",
    "Sink",
    "Linear continuation",
    "Branching source",
    "Merging sink",
    "Internal branching",
    "Internal merging",
    "Merge and branch",
    "Other",
]

INCOMING_LABELS = ["No incoming", "One incoming", "Multiple incoming"]
OUTGOING_LABELS = ["No outgoing", "One outgoing", "Multiple outgoing"]
ROLE_GRID = [
    ("Isolated", 0, 0),
    ("Source", 0, 1),
    ("Branching source", 0, 2),
    ("Sink", 1, 0),
    ("Linear continuation", 1, 1),
    ("Internal branching", 1, 2),
    ("Merging sink", 2, 0),
    ("Internal merging", 2, 1),
    ("Merge and branch", 2, 2),
]


def _degree_mask(degree: pd.Series, level: int) -> pd.Series:
    if level == 0:
        return degree.eq(0)
    if level == 1:
        return degree.eq(1)
    return degree.gt(1)


def build_role_count_grid(nodes: pd.DataFrame) -> pd.DataFrame:
    in_degree = pd.to_numeric(nodes["in_degree"], errors="coerce").fillna(0)
    out_degree = pd.to_numeric(nodes["out_degree"], errors="coerce").fillna(0)
    rows = []
    for role, row, col in ROLE_GRID:
        mask = _degree_mask(in_degree, row) & _degree_mask(out_degree, col)
        rows.append(
            {
                "incoming": INCOMING_LABELS[row],
                "outgoing": OUTGOING_LABELS[col],
                "role_group": role,
                "row": row,
                "col": col,
                "n": int(mask.sum()),
            }
        )
    return pd.DataFrame(rows)


def build_transition_context(
    window: pd.DataFrame, nodes: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    nodes = nodes.copy()
    nodes["wn_mid_date"] = pd.to_datetime(nodes["wn_mid_date"], errors="coerce")
    window_dates = (
        nodes[["window_id", "window_idx", "wn_mid_date"]]
        .drop_duplicates(["window_id", "window_idx"])
        .sort_values("window_idx")
    )
    window = window.merge(
        window_dates[["window_id", "window_idx", "wn_mid_date"]],
        on=["window_id", "window_idx"],
        how="left",
    )
    policy_context_cols = ["window_id", "window_idx", "wn_mid_date", "policy_era"]
    policy_context = (
        nodes[[col for col in policy_context_cols if col in nodes.columns]]
        .drop_duplicates(["window_id", "window_idx"])
        .sort_values("window_idx")
    )

    nodes["role_group"] = (
        nodes["primary_graph_role"].map(CLUSTER_ROLE_GROUPS).fillna("Other")
    )
    role_counts = (
        nodes.groupby(["window_idx", "wn_mid_date", "role_group"], dropna=False)
        .size()
        .rename("n")
        .reset_index()
    )
    role_pivot = role_counts.pivot_table(
        index="wn_mid_date",
        columns="role_group",
        values="n",
        fill_value=0,
    ).sort_index()
    role_pivot = role_pivot[[col for col in PREFERRED_ROLES if col in role_pivot]]
    return window, policy_context, role_pivot


def draw_transition_counts(ax, window: pd.DataFrame, policy_context: pd.DataFrame) -> None:
    add_policy_bands(ax, policy_context)
    nodes_line = ax.plot(
        window["wn_mid_date"], window["n_nodes"], color="#1f4e79", lw=1.3
    )[0]
    edges_line = ax.plot(
        window["wn_mid_date"], window["n_out_edges"], color="#d95f02", lw=1.3
    )[0]
    ax.set_ylabel("Count")
    ax.legend([nodes_line, edges_line], ["Nodes", "Outgoing edges"], loc="upper left")
    date_axis(ax)


def draw_role_distribution(
    ax, role_pivot: pd.DataFrame, policy_context: pd.DataFrame
) -> None:
    add_policy_bands(ax, policy_context)
    colors = [ROLE_COLORS.get(role, "#bab0ac") for role in role_pivot.columns]
    ax.stackplot(
        role_pivot.index,
        role_pivot.T.to_numpy(),
        labels=role_pivot.columns,
        colors=colors,
        linewidth=0,
    )
    ax.set_ylabel("Nodes")
    ax.legend(loc="upper left")
    date_axis(ax)


def draw_component_size_distribution(ax, components: pd.DataFrame) -> None:
    comp_sizes = components["n_nodes"].dropna().sort_values().to_numpy()
    if comp_sizes.size:
        ccdf = 1.0 - np.arange(len(comp_sizes)) / len(comp_sizes)
        ax.plot(comp_sizes, ccdf, color="#1f4e79", lw=1.4)
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.set_xlabel("Weak-component size (nodes)")
    ax.set_ylabel("Pr(component size >= x)")


def _contrast_text_color(color: str) -> str:
    red, green, blue = to_rgb(color)
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "#ffffff" if luminance < 0.48 else "#222222"


def draw_role_count_grid(ax, table: pd.DataFrame) -> None:
    for row in table.itertuples(index=False):
        y0 = 2 - row.row
        color = ROLE_COLORS.get(row.role_group, "#bab0ac")
        ax.add_patch(
            Rectangle(
                (row.col + 0.03, y0 + 0.04),
                0.94,
                0.86,
                facecolor=color,
                edgecolor="white",
                lw=1.0,
                zorder=0,
            )
        )
        ax.text(
            row.col + 0.50,
            y0 + 0.47,
            f"{int(row.n):,}",
            ha="center",
            va="center",
            color=_contrast_text_color(color),
            fontweight="bold",
            fontsize=11,
        )

    for col, label in enumerate(OUTGOING_LABELS):
        ax.text(
            col + 0.5,
            3.06,
            label.replace(" ", "\n"),
            ha="center",
            va="bottom",
            fontsize=9,
            linespacing=0.9,
        )
    for row, label in enumerate(INCOMING_LABELS):
        ax.text(-0.08, 2.47 - row, label, ha="right", va="center", fontsize=9)

    ax.set_xlim(-0.48, 3.02)
    ax.set_ylim(-0.03, 3.34)
    ax.set_axis_off()


def build(paths: Paths) -> dict[str, object]:
    transition_window = read_table(paths, "transition_window_summary")
    transition_nodes = read_table(paths, "transition_node_table")
    components = read_table(paths, "transition_component_summary")

    table = build_role_count_grid(transition_nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")

    window, policy_context, role_pivot = build_transition_context(
        transition_window, transition_nodes
    )

    fig, axes = new_figure(
        nrows=2,
        ncols=2,
        width="double",
        height_in=6,
        constrained_layout=True,
    )

    draw_transition_counts(axes[0, 0], window, policy_context)
    draw_component_size_distribution(axes[0, 1], components)
    draw_role_distribution(axes[1, 0], role_pivot, policy_context)
    draw_role_count_grid(axes[1, 1], table)

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
