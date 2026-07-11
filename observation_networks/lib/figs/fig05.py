"""Build Chapter 4 Figure 5: temporal cluster-transition graph baseline."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

from matplotlib.patches import Rectangle
from matplotlib.axes import Axes
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    add_policy_bands,
    date_axis,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)


FIGURE_NAME = "fig_ch4_transition_graph_baseline"


CLUSTER_ROLE_GROUPS = {
    "isolated": "Isolated",
    "single_outgoing_source": "Source",
    "source_branching": "Branching source",
    "single_incoming_sink": "Sink",
    "merging_sink": "Merging sink",
    "linear_continuation": "Linear continuation",
    "internal_branching": "Internal branching",
    "internal_merging": "Internal merging",
    "merge_and_branch": "Merge and branch",
    "other": "Other",
}

ROLE_COLORS = {
    "Isolated": "#4c78a8",
    "Source": "#f58518",
    "Branching source": "#e45756",
    "Sink": "#54a24b",
    "Merging sink": "#e15ebc",
    "Linear continuation": "#8e6c8a",
    "Internal branching": "#ff9da6",
    "Internal merging": "#9d755d",
    "Merge and branch": "#72b7b2",
    "Other": "#bab0ac",
}


def draw_role_edge(
    ax: Axes,
    start: tuple[float, float],
    end: tuple[float, float],
) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "-|>",
            "color": "#666666",
            "lw": 0.9,
            "mutation_scale": 8,
            "shrinkA": 2,
            "shrinkB": 2,
        },
        zorder=2,
    )


def draw_role_cell(
    ax: Axes,
    *,
    col: int,
    row: int,
    label: str,
    incoming: int,
    outgoing: int,
) -> None:
    y0 = 2 - row
    x0 = col
    color = ROLE_COLORS[label]
    ax.add_patch(
        Rectangle(
            (x0 + 0.03, y0 + 0.04),
            0.94,
            0.86,
            facecolor="#f7f7f7",
            edgecolor="#d0d0d0",
            lw=0.55,
            zorder=0,
        )
    )

    cx = x0 + 0.50
    cy = y0 + 0.55
    if incoming:
        incoming_shifts = [0.0] if incoming == 1 else [-0.12, 0.12]
        for shift in incoming_shifts:
            source = (x0 + 0.18, cy + shift)
            ax.scatter(*source, s=20, color="#9a9a9a", zorder=3)
            draw_role_edge(ax, (source[0] + 0.04, source[1]), (cx - 0.09, cy))
    if outgoing:
        outgoing_shifts = [0.0] if outgoing == 1 else [-0.12, 0.12]
        for shift in outgoing_shifts:
            target = (x0 + 0.82, cy + shift)
            ax.scatter(*target, s=20, color="#9a9a9a", zorder=3)
            draw_role_edge(ax, (cx + 0.09, cy), (target[0] - 0.04, target[1]))

    ax.scatter(
        cx,
        cy,
        s=96,
        color=color,
        edgecolor="white",
        linewidth=0.9,
        zorder=4,
    )
    ax.text(
        cx,
        y0 + 0.15,
        label,
        ha="center",
        va="center",
        fontsize=7.1,
        color="#222222",
    )


def draw_graph_role_schematic(ax: Axes) -> None:
    roles = [
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
    for idx, (label, incoming, outgoing) in enumerate(roles):
        row, col = divmod(idx, 3)
        draw_role_cell(
            ax,
            col=col,
            row=row,
            label=label,
            incoming=incoming,
            outgoing=outgoing,
        )

    for col, text in enumerate(["No outgoing", "One outgoing", "Multiple outgoing"]):
        ax.text(col + 0.5, 3.10, text, ha="center", va="bottom", fontsize=8.0)
    for row, text in enumerate(["No incoming", "One incoming", "Multiple incoming"]):
        ax.text(-0.08, 2.45 - row, text, ha="right", va="center", fontsize=8.0)

    ax.set_xlim(-0.42, 3.02)
    ax.set_ylim(-0.03, 3.25)
    ax.set_title("Graph role defined by continuity-edge degree", pad=8)
    ax.set_axis_off()
    panel_label(ax, "A")


def build(paths: Paths) -> None:
    window = read_table(paths, "transition_window_summary")
    nodes = read_table(paths, "transition_node_table")
    components = read_table(paths, "transition_component_summary")
    window_coverage = read_table(paths, "window_coverage")
    nodes["wn_mid_date"] = pd.to_datetime(nodes["wn_mid_date"], errors="coerce")
    window_coverage["wn_mid_date"] = pd.to_datetime(
        window_coverage["wn_mid_date"], errors="coerce"
    )
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
        index="wn_mid_date", columns="role_group", values="n", fill_value=0
    ).sort_index()
    preferred_roles = [
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
    role_pivot = role_pivot[[col for col in preferred_roles if col in role_pivot]]

    fig, placeholder_ax = styled_new_figure(
        width="double",
        height_in=8.4,
        constrained_layout=True,
    )
    placeholder_ax.remove()
    grid = fig.add_gridspec(3, 2, height_ratios=[1.35, 1.0, 1.0])

    ax = fig.add_subplot(grid[0, :])
    draw_graph_role_schematic(ax)

    ax = fig.add_subplot(grid[1, 0])
    add_policy_bands(ax, window_coverage)
    nodes_line = ax.plot(
        window["wn_mid_date"], window["n_nodes"], color="#1f4e79", lw=1.3
    )[0]
    edges_line = ax.plot(
        window["wn_mid_date"], window["n_out_edges"], color="#d95f02", lw=1.3
    )[0]
    ax.set_title("Transition graph size")
    ax.set_ylabel("Count")
    ax.legend([nodes_line, edges_line], ["Nodes", "Outgoing edges"], loc="upper left")
    date_axis(ax)
    panel_label(ax, "B")

    ax = fig.add_subplot(grid[1, 1])
    add_policy_bands(ax, window_coverage)
    colors = [ROLE_COLORS.get(role, "#bab0ac") for role in role_pivot.columns]
    ax.stackplot(
        role_pivot.index,
        role_pivot.T.to_numpy(),
        labels=role_pivot.columns,
        colors=colors,
        linewidth=0,
    )
    ax.set_title("Graph-role composition")
    ax.set_ylabel("Nodes")
    date_axis(ax)
    panel_label(ax, "C")

    ax = fig.add_subplot(grid[2, 0])
    comp_sizes = components["n_nodes"].dropna().sort_values().to_numpy()
    ccdf = 1.0 - np.arange(len(comp_sizes)) / len(comp_sizes)
    ax.plot(comp_sizes, ccdf, color="#1f4e79", lw=1.4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Weak-component size distribution")
    ax.set_xlabel("Weak-component size (nodes)")
    ax.set_ylabel("Pr(component size >= x)")
    panel_label(ax, "D")

    ax = fig.add_subplot(grid[2, 1])
    add_policy_bands(ax, window_coverage)
    isolates_line = ax.plot(
        window["wn_mid_date"],
        window["n_isolates"],
        color=ROLE_COLORS["Isolated"],
        lw=1.1,
    )[0]
    branching_line = ax.plot(
        window["wn_mid_date"],
        window["n_branching"],
        color=ROLE_COLORS["Branching source"],
        lw=1.1,
    )[0]
    merging_line = ax.plot(
        window["wn_mid_date"],
        window["n_merging"],
        color=ROLE_COLORS["Merging sink"],
        lw=1.1,
    )[0]
    ax.set_title("Selected degree-class counts")
    ax.set_ylabel("Nodes")
    ax.legend(
        [isolates_line, branching_line, merging_line],
        ["Isolated", "Multiple outgoing", "Multiple incoming"],
        loc="upper left",
    )
    date_axis(ax)
    panel_label(ax, "E")
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
