"""Build the transition-graph role schematic."""

from __future__ import annotations

import argparse

from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

from .common import (
    Paths,
    add_common_args,
    new_figure,
    paths_from_args,
    styled_save_figure,
)

FIGURE_NAME = "fig_transition_graph_roles"


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
        ax.text(col + 0.5, 3.10, text, ha="center", va="bottom")
    for row, text in enumerate(["No incoming", "One incoming", "Multiple incoming"]):
        ax.text(-0.08, 2.45 - row, text, ha="right", va="center")

    ax.set_xlim(-0.42, 3.02)
    ax.set_ylim(-0.03, 3.25)
    ax.set_axis_off()

def build(paths: Paths) -> None:
    fig, ax = new_figure(width="double", height_in=4.2, constrained_layout=True)
    draw_graph_role_schematic(ax)
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
