"""Build Chapter 5 Figure 8: transition graph diagnostics and selection flow."""

from __future__ import annotations

import argparse
from textwrap import fill

import numpy as np
import pandas as pd

from ..sse.config import MIN_CLUSTER_SIZE
from ..sse.io import write_table
from .common import (
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
from .fig01 import CLUSTER_ROLE_GROUPS, ROLE_COLORS

FILE_NAME = "ch5_transition_graph_characteristics"

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


def build_selection_funnel(nodes: pd.DataFrame) -> pd.DataFrame:
    tier = nodes["candidate_tier"]
    rows = [
        ("All transition-network clusters", len(nodes)),
        (
            f"Size eligible (cluster size ≥ {MIN_CLUSTER_SIZE})",
            int(nodes["sse_tested"].fillna(False).sum()),
        ),
        (
            "Background or low information",
            int(tier.eq("background_or_low_information").sum()),
        ),
        ("Possible review", int(tier.eq("possible_review").sum())),
        ("High-priority local burst", int(tier.eq("high_priority_burst").sum())),
        ("High-priority onward burden", int(tier.eq("high_priority_burden").sum())),
        ("High priority on both axes", int(tier.eq("high_priority_both_axes").sum())),
    ]
    out = pd.DataFrame(rows, columns=["stage", "n"])
    out["order"] = np.arange(len(out))
    out["pct_all"] = out["n"] / len(nodes)
    return out


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
    panel_label(ax, "A")


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
    panel_label(ax, "C")


def draw_component_size_distribution(ax, components: pd.DataFrame) -> None:
    comp_sizes = components["n_nodes"].dropna().sort_values().to_numpy()
    if comp_sizes.size:
        ccdf = 1.0 - np.arange(len(comp_sizes)) / len(comp_sizes)
        ax.plot(comp_sizes, ccdf, color="#1f4e79", lw=1.4)
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.set_xlabel("Weak-component size (nodes)")
    ax.set_ylabel("Pr(component size >= x)")
    panel_label(ax, "B")


def draw_selection_funnel(ax, table: pd.DataFrame) -> None:
    colors = [
        "#4C78A8",
        "#72A0C1",
        "#B8B8B8",
        "#E6AB02",
        "#D55E00",
        "#0072B2",
        "#7B3294",
    ]
    y = np.arange(len(table))[::-1]
    bars = ax.barh(y, table["n"], color=colors)
    ax.set_xscale("log")
    ax.set_yticks(y, [fill(stage, width=18) for stage in table["stage"]])
    ax.set_xlabel("Clusters (log scale)")
    for bar, value in zip(bars, table["n"]):
        label_x = max(value, 0.8) * 1.12 if value > 0 else 35
        ax.text(
            label_x,
            bar.get_y() + bar.get_height() / 2,
            f"{int(value):,}",
            va="center",
        )
    panel_label(ax, "D")


def build(paths: Paths) -> dict[str, object]:
    transition_window = read_table(paths, "transition_window_summary")
    transition_nodes = read_table(paths, "transition_node_table")
    components = read_table(paths, "transition_component_summary")
    clusters = read_table(paths, "cluster_table")

    table = build_selection_funnel(clusters)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")
    table = table.sort_values("order")

    window, policy_context, role_pivot = build_transition_context(
        transition_window, transition_nodes
    )

    fig, axes = styled_new_figure(
        nrows=2,
        ncols=2,
        width="double", 
        height_in=6, 
        constrained_layout=True
    )

    draw_transition_counts(axes[0, 0], window, policy_context)
    draw_component_size_distribution(axes[0, 1], components)
    draw_role_distribution(axes[1, 0], role_pivot, policy_context)
    draw_selection_funnel(axes[1, 1], table)

    outputs = styled_save_figure(fig, paths, f"fig_{FILE_NAME}", tight=False)
    return {"figure": fig, "outputs": outputs}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    paths = paths_from_args(parser.parse_args())
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
