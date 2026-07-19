"""Build Chapter 5 Figure 7: representative candidate trajectories."""

from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from .common import (
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)
from ..sse.io import write_table

FILE_NAME = "ch5_candidate_exemplars"
COLORS = {
    "background_or_low_information": "#B8B8B8",
    "possible_review": "#E6AB02",
    "high_priority_burst": "#D55E00",
    "high_priority_burden": "#0072B2",
    "high_priority_both_axes": "#7B3294",
}


def build_exemplars(
    nodes: pd.DataFrame, edges: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tested = nodes.loc[nodes["sse_tested"].fillna(False)].copy()
    selections: list[tuple[str, pd.Series]] = []
    specs = [
        (
            "Local burst",
            tested["candidate_tier"].isin(
                {"high_priority_burst", "high_priority_both_axes"}
            ),
            "burst_score_null_z",
        ),
        (
            "Onward burden",
            tested["candidate_tier"].isin(
                {"high_priority_burden", "high_priority_both_axes"}
            ),
            "burden_score_null_z",
        ),
        (
            "Possible review",
            tested["candidate_tier"].eq("possible_review"),
            "max_axis_null_z",
        ),
    ]
    for label, mask, score in specs:
        pool = tested.loc[mask].sort_values([score, "cluster_size"], ascending=False)
        if not pool.empty:
            selections.append((label, pool.iloc[0]))
    anchor = selections[0][1]
    background = tested.loc[
        tested["candidate_tier"].eq("background_or_low_information")
    ].copy()
    background["match_distance"] = (
        (background["window_idx"] - anchor["window_idx"]).abs()
        + (
            np.log1p(background["cluster_size"]) - np.log1p(anchor["cluster_size"])
        ).abs()
        + (~background["clade"].eq(anchor["clade"])).astype(float)
    )
    selections.append(
        ("Matched background", background.sort_values("match_distance").iloc[0])
    )

    node_parts, edge_parts = [], []
    for label, focal in selections:
        focal_id = focal["cluster_id"]
        related_edges = edges.loc[
            (edges["source"].eq(focal_id)) | (edges["target"].eq(focal_id))
        ].copy()
        ids = {focal_id} | set(related_edges["source"]) | set(related_edges["target"])
        subnodes = nodes.loc[
            nodes["cluster_id"].isin(ids),
            [
                "cluster_id",
                "window_idx",
                "window_id",
                "cluster_size",
                "candidate_tier",
                "burst_score_null_z",
                "burden_score_null_z",
                "policy_era",
                "policy_period",
                "who_voc",
            ],
        ].copy()
        subnodes["example"] = label
        subnodes["is_focal"] = subnodes["cluster_id"].eq(focal_id)
        related_edges["example"] = label
        node_parts.append(subnodes)
        edge_parts.append(related_edges)
    return pd.concat(node_parts, ignore_index=True), pd.concat(
        edge_parts, ignore_index=True
    )


def _draw(ax: Axes, nodes: pd.DataFrame, edges: pd.DataFrame, title: str) -> None:
    windows = sorted(nodes["window_idx"].unique())
    x_lookup = {w: i for i, w in enumerate(windows)}
    positions = {}
    for window, group in nodes.groupby("window_idx"):
        ordered = group.sort_values(
            ["is_focal", "cluster_size"], ascending=[False, False]
        )
        ys = np.arange(len(ordered)) - (len(ordered) - 1) / 2
        for (_, row), y in zip(ordered.iterrows(), ys):
            positions[row["cluster_id"]] = (x_lookup[window], float(y))
    for edge in edges.itertuples(index=False):
        if edge.source in positions and edge.target in positions:
            ax.annotate(
                "",
                xy=positions[edge.target],
                xytext=positions[edge.source],
                arrowprops={
                    "arrowstyle": "->",
                    "color": "#777777",
                    "lw": 0.7 + np.log1p(edge.n_shared_sequences) / 3, # type: ignore
                },
                zorder=1,
            )
    for row in nodes.itertuples(index=False):
        x, y = positions[row.cluster_id]
        color = COLORS.get(row.candidate_tier, "#B8B8B8") # type: ignore
        ax.scatter(
            x,
            y,
            s=np.clip(np.sqrt(row.cluster_size) * 22, 45, 250), # type: ignore
            color=color,
            edgecolor="black" if row.is_focal else "white",
            linewidth=1 if row.is_focal else 0.4,
            zorder=2,
        )
        if row.is_focal:
            ax.text(x, y - 0.42, f"n={int(row.cluster_size)}", ha="center", va="top") # type: ignore
    focal = nodes.loc[nodes["is_focal"]].iloc[0]
    ax.set_title(f"{title}\n{focal['policy_period']} · {focal['who_voc']}")
    ax.set_xticks(range(len(windows)), [f"W{int(w):03d}" for w in windows])
    ax.set_yticks([])
    ax.set_ylim(-max(1.4, len(nodes) / 2), max(1.4, len(nodes) / 2))
    ax.set_xlabel("Transition window")


def build(paths: Paths) -> dict[str, object]:
    all_nodes = read_table(paths, "cluster_table")
    all_edges = read_table(paths, "edge_table")
    nodes, edges = build_exemplars(all_nodes, all_edges)
    write_table(nodes, paths.result_table_dir, f"tab_{FILE_NAME}_nodes")
    write_table(edges, paths.result_table_dir, f"tab_{FILE_NAME}_edges")
    order = ["Local burst", "Onward burden", "Possible review", "Matched background"]
    fig, axes = styled_new_figure(
        width="double", height_in=6.2, nrows=2, ncols=2, constrained_layout=True
    )
    for ax, example, label in zip(axes.flat, order, "ABCD"):
        _draw(
            ax,
            nodes.loc[nodes["example"].eq(example)],
            edges.loc[edges["example"].eq(example)],
            example,
        )
        panel_label(ax, label)
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
