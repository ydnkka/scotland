"""Temporal cluster-transition graph construction and descriptive summaries."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import networkx as nx
import numpy as np
import pandas as pd

from .clusters import build_cluster_table


@dataclass(frozen=True)
class TransitionOutputs:
    """Container for Chapter 4 temporal cluster-transition outputs."""

    cluster_table: pd.DataFrame
    edge_table: pd.DataFrame
    node_table: pd.DataFrame
    graph_summary: pd.DataFrame
    window_summary: pd.DataFrame
    component_summary: pd.DataFrame


def build_transition_edges(df: pd.DataFrame) -> pd.DataFrame:
    """Build adjacent-window directed edges from shared sequence membership."""
    required = {"sequence_id", "cluster_id", "window_id", "window_idx"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing transition edge columns: {sorted(missing)}")

    seq_nodes = (
        df[["sequence_id", "cluster_id", "window_id", "window_idx"]]
        .drop_duplicates()
        .sort_values(["sequence_id", "window_idx"])
    )

    edge_map: dict[tuple[object, object, object, object, int, int], set[str]] = (
        defaultdict(set)
    )
    for sequence_id, group in seq_nodes.groupby("sequence_id", sort=False):
        records = list(
            group[["cluster_id", "window_id", "window_idx"]].itertuples(
                index=False,
                name=None,
            )
        )
        for i, (source_node, source_window, source_idx) in enumerate(records):
            for target_node, target_window, target_idx in records[i + 1 :]:
                delta = int(target_idx) - int(source_idx)
                if delta == 1 and source_node != target_node:
                    edge_map[
                        (
                            source_node,
                            target_node,
                            source_window,
                            target_window,
                            int(source_idx),
                            int(target_idx),
                        )
                    ].add(str(sequence_id))
                elif delta > 1:
                    break

    rows = [
        {
            "source": source,
            "target": target,
            "source_window_id": source_window,
            "target_window_id": target_window,
            "source_window_idx": source_idx,
            "target_window_idx": target_idx,
            "n_shared_sequences": len(shared),
        }
        for (
            source,
            target,
            source_window,
            target_window,
            source_idx,
            target_idx,
        ), shared in edge_map.items()
    ]

    out = pd.DataFrame(
        rows,
        columns=[
            "source",
            "target",
            "source_window_id",
            "target_window_id",
            "source_window_idx",
            "target_window_idx",
            "n_shared_sequences",
        ],
    )
    if out.empty:
        return out
    return out.sort_values(
        [
            "source_window_idx",
            "target_window_idx",
            "n_shared_sequences",
            "source",
            "target",
        ],
        ascending=[True, True, False, True, True],
        ignore_index=True,
    )


def _component_map(nodes: pd.Series, edge_table: pd.DataFrame) -> tuple[dict, dict]:
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes.dropna().unique())
    if not edge_table.empty:
        graph.add_weighted_edges_from(
            edge_table[["source", "target", "n_shared_sequences"]].itertuples(
                index=False,
                name=None,
            )
        )

    components = sorted(
        nx.weakly_connected_components(graph),
        key=lambda items: (-len(items), sorted(items)[0]),
    )
    node_to_component = {
        node: f"CC{idx:05d}"
        for idx, component in enumerate(components, start=1)
        for node in component
    }
    component_size = {
        f"CC{idx:05d}": len(component)
        for idx, component in enumerate(components, start=1)
    }
    return node_to_component, component_size


def _add_degree_metrics(
    node_table: pd.DataFrame, edge_table: pd.DataFrame
) -> pd.DataFrame:
    out = node_table.copy()
    if edge_table.empty:
        for col in ("out_degree", "out_strength", "in_degree", "in_strength"):
            out[col] = 0
        return out

    out_metrics = (
        edge_table.groupby("source", dropna=False)
        .agg(
            out_degree=("target", "nunique"),
            out_strength=("n_shared_sequences", "sum"),
        )
        .rename_axis("cluster_id")
        .reset_index()
    )
    in_metrics = (
        edge_table.groupby("target", dropna=False)
        .agg(
            in_degree=("source", "nunique"),
            in_strength=("n_shared_sequences", "sum"),
        )
        .rename_axis("cluster_id")
        .reset_index()
    )
    out = out.merge(out_metrics, on="cluster_id", how="left").merge(
        in_metrics,
        on="cluster_id",
        how="left",
    )
    for col in ("out_degree", "out_strength", "in_degree", "in_strength"):
        out[col] = out[col].fillna(0).astype(int)
    return out


def _add_downstream_burden(
    node_table: pd.DataFrame, edge_table: pd.DataFrame
) -> pd.DataFrame:
    out = node_table.copy()
    out["downstream_cluster_burden"] = 0
    out["mean_successor_cluster_size"] = np.nan
    if edge_table.empty:
        return out

    target_sizes = out[["cluster_id", "cluster_size"]].rename(
        columns={"cluster_id": "target", "cluster_size": "_target_cluster_size"}
    )
    downstream = (
        edge_table.merge(target_sizes, on="target", how="left")
        .groupby("source", dropna=False)
        .agg(
            downstream_cluster_burden=("_target_cluster_size", "sum"),
            mean_successor_cluster_size=("_target_cluster_size", "mean"),
        )
        .rename_axis("cluster_id")
        .reset_index()
    )
    out = out.drop(columns=["downstream_cluster_burden", "mean_successor_cluster_size"])
    out = out.merge(downstream, on="cluster_id", how="left")
    out["downstream_cluster_burden"] = out["downstream_cluster_burden"].fillna(0)
    return out


def _assign_graph_roles(node_table: pd.DataFrame) -> pd.DataFrame:
    out = node_table.copy()
    in_degree = out["in_degree"].fillna(0)
    out_degree = out["out_degree"].fillna(0)
    out["has_incoming"] = in_degree.gt(0)
    out["has_outgoing"] = out_degree.gt(0)
    out["has_branching"] = out_degree.gt(1)
    out["has_merging"] = in_degree.gt(1)
    out["primary_graph_role"] = "other"
    out.loc[in_degree.eq(0) & out_degree.eq(0), "primary_graph_role"] = "isolated"
    out.loc[
        in_degree.eq(0) & out_degree.eq(1),
        "primary_graph_role",
    ] = "single_outgoing_source"
    out.loc[in_degree.eq(0) & out_degree.gt(1), "primary_graph_role"] = (
        "source_branching"
    )
    out.loc[in_degree.eq(1) & out_degree.eq(0), "primary_graph_role"] = (
        "single_incoming_sink"
    )
    out.loc[in_degree.gt(1) & out_degree.eq(0), "primary_graph_role"] = "merging_sink"
    out.loc[in_degree.eq(1) & out_degree.eq(1), "primary_graph_role"] = (
        "linear_continuation"
    )
    out.loc[in_degree.eq(1) & out_degree.gt(1), "primary_graph_role"] = (
        "internal_branching"
    )
    out.loc[in_degree.gt(1) & out_degree.eq(1), "primary_graph_role"] = (
        "internal_merging"
    )
    out.loc[in_degree.gt(1) & out_degree.gt(1), "primary_graph_role"] = (
        "merge_and_branch"
    )
    return out


def build_transition_node_table(
    cluster_table: pd.DataFrame,
    edge_table: pd.DataFrame,
) -> pd.DataFrame:
    """Attach graph degree, component, role, and downstream baseline metrics."""
    out = _add_degree_metrics(cluster_table, edge_table)
    out = _add_downstream_burden(out, edge_table)
    out = _assign_graph_roles(out)

    node_to_component, component_size = _component_map(out["cluster_id"], edge_table)
    out["weak_component_id"] = out["cluster_id"].map(node_to_component)
    out["weak_component_size"] = out["weak_component_id"].map(component_size)
    return out


def build_graph_summary(
    node_table: pd.DataFrame, edge_table: pd.DataFrame
) -> pd.DataFrame:
    """Return scalar descriptive summaries for the transition graph."""
    rows = [
        ("nodes", node_table["cluster_id"].nunique()),
        ("edges", len(edge_table)),
        ("weak_components", node_table["weak_component_id"].nunique()),
        ("isolated_nodes", int(node_table["primary_graph_role"].eq("isolated").sum())),
        ("branching_nodes", int(node_table["has_branching"].sum())),
        ("merging_nodes", int(node_table["has_merging"].sum())),
        ("max_weak_component_size", node_table["weak_component_size"].max()),
        ("mean_out_degree", node_table["out_degree"].mean()),
        ("mean_in_degree", node_table["in_degree"].mean()),
        ("total_shared_sequence_weight", edge_table["n_shared_sequences"].sum()),
    ]
    return pd.DataFrame(
        {
            "metric": [metric for metric, _ in rows],
            "value": [value for _, value in rows],
        }
    )


def build_transition_window_summary(
    node_table: pd.DataFrame,
    edge_table: pd.DataFrame,
) -> pd.DataFrame:
    """Summarise transition graph nodes and outgoing edges by source window."""
    nodes = (
        node_table.groupby(["window_id", "window_idx"], dropna=False)
        .agg(
            n_nodes=("cluster_id", "nunique"),
            n_isolates=("primary_graph_role", lambda x: (x == "isolated").sum()),
            n_branching=("has_branching", "sum"),
            n_merging=("has_merging", "sum"),
            median_out_degree=("out_degree", "median"),
            max_out_degree=("out_degree", "max"),
        )
        .reset_index()
    )
    if edge_table.empty:
        nodes["n_out_edges"] = 0
        nodes["out_edge_weight"] = 0
        return nodes.sort_values("window_idx")

    edges = (
        edge_table.groupby(["source_window_id", "source_window_idx"], dropna=False)
        .agg(
            n_out_edges=("target", "size"),
            out_edge_weight=("n_shared_sequences", "sum"),
        )
        .reset_index()
        .rename(
            columns={
                "source_window_id": "window_id",
                "source_window_idx": "window_idx",
            }
        )
    )
    out = nodes.merge(edges, on=["window_id", "window_idx"], how="left")
    out[["n_out_edges", "out_edge_weight"]] = out[
        ["n_out_edges", "out_edge_weight"]
    ].fillna(0)
    return out.sort_values("window_idx")


def build_component_summary(node_table: pd.DataFrame) -> pd.DataFrame:
    """Summarise weakly connected components in the transition graph."""
    return (
        node_table.groupby("weak_component_id", dropna=False)
        .agg(
            n_nodes=("cluster_id", "nunique"),
            first_window_idx=("window_idx", "min"),
            last_window_idx=("window_idx", "max"),
            n_windows=("window_id", "nunique"),
            total_cluster_size=("cluster_size", "sum"),
            max_cluster_size=("cluster_size", "max"),
            n_health_boards=("modal_health_board", "nunique")
            if "modal_health_board" in node_table.columns
            else ("cluster_id", "size"),
        )
        .reset_index()
        .sort_values(["n_nodes", "total_cluster_size"], ascending=[False, False])
    )


def build_transition_outputs(sequence_df: pd.DataFrame) -> TransitionOutputs:
    """Build all Chapter 4 transition-graph baseline outputs."""
    cluster_table = build_cluster_table(sequence_df)
    edge_table = build_transition_edges(sequence_df)
    node_table = build_transition_node_table(cluster_table, edge_table)
    return TransitionOutputs(
        cluster_table=cluster_table,
        edge_table=edge_table,
        node_table=node_table,
        graph_summary=build_graph_summary(node_table, edge_table),
        window_summary=build_transition_window_summary(node_table, edge_table),
        component_summary=build_component_summary(node_table),
    )
