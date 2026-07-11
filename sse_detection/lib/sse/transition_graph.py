"""Temporal cluster-transition graph construction and summaries for SSE outputs."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import networkx as nx
import numpy as np
import pandas as pd

from .entropy import onward_edge_entropy


@dataclass(frozen=True)
class TransitionSummaryOutputs:
    """Descriptive transition-graph tables written by the SSE pipeline."""

    node_table: pd.DataFrame
    graph_summary: pd.DataFrame
    window_summary: pd.DataFrame
    component_summary: pd.DataFrame

    def tables(self) -> dict[str, pd.DataFrame]:
        """Return the standard output-table mapping."""
        return {
            "transition_node_table": self.node_table,
            "transition_graph_summary": self.graph_summary,
            "transition_window_summary": self.window_summary,
            "transition_component_summary": self.component_summary,
        }


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


def build_transition_network(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, nx.DiGraph, dict[object, str]]:
    """Build the directed adjacent-window cluster-transition graph."""
    edge_table = build_transition_edges(df)
    nodes = df["cluster_id"].dropna().unique()

    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    if not edge_table.empty:
        graph.add_weighted_edges_from(
            edge_table[["source", "target", "n_shared_sequences"]].itertuples(
                index=False,
                name=None,
            )
        )

    node_to_component, _ = component_maps(graph)
    return edge_table, graph, node_to_component


def component_maps(graph: nx.DiGraph) -> tuple[dict[object, str], dict[str, int]]:
    """Return deterministic weak-component ids and sizes for graph nodes."""
    components = sorted(
        nx.weakly_connected_components(graph),
        key=lambda nodes: (-len(nodes), sorted(nodes)[0]),
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


def _graph_from_node_edge_tables(
    nodes: pd.Series,
    edge_table: pd.DataFrame,
) -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes.dropna().unique())
    if not edge_table.empty:
        graph.add_weighted_edges_from(
            edge_table[["source", "target", "n_shared_sequences"]].itertuples(
                index=False,
                name=None,
            )
        )
    return graph


def build_edge_flow_metrics(edge_table: pd.DataFrame) -> pd.DataFrame:
    """Build per-node incoming and outgoing transition-flow metrics."""
    edge_stats = onward_edge_entropy(
        edge_table,
        source_col="source",
        weight_col="n_shared_sequences",
    )

    in_cols = ["cluster_id", "in_degree", "in_strength"]
    if edge_table.empty:
        in_metrics = pd.DataFrame(columns=in_cols)
    else:
        in_metrics = (
            edge_table.groupby("target", dropna=False)
            .agg(
                in_degree=("source", "nunique"),
                in_strength=("n_shared_sequences", "sum"),
            )
            .rename_axis("cluster_id")
            .reset_index()
        )

    return edge_stats.merge(in_metrics, on="cluster_id", how="outer")


def add_edge_flow_metrics(
    node_table: pd.DataFrame,
    edge_table: pd.DataFrame,
) -> pd.DataFrame:
    """Attach per-node incoming/outgoing transition-flow metrics."""
    flow = build_edge_flow_metrics(edge_table)
    if flow.empty:
        return node_table.copy()

    metric_cols = [col for col in flow.columns if col != "cluster_id"]
    out = node_table.drop(columns=metric_cols, errors="ignore").copy()
    return out.merge(flow, on="cluster_id", how="left")


def add_downstream_burden(
    node_table: pd.DataFrame,
    edge_table: pd.DataFrame,
) -> pd.DataFrame:
    out = node_table.drop(
        columns=["downstream_cluster_burden", "mean_successor_cluster_size"],
        errors="ignore",
    ).copy()
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


def build_new_downstream_metrics(
    sequence_df: pd.DataFrame,
    edge_table: pd.DataFrame,
    *,
    min_shared_sequences: int = 2,
) -> pd.DataFrame:
    """Compute overlap-adjusted downstream burden for each source node.

    ``new_downstream_burden`` counts unique sequences in all adjacent-window
    child clusters that are not already present in the source cluster. The
    supported version applies the same calculation after excluding weak edges
    with fewer than ``min_shared_sequences`` shared sequences.
    """
    if min_shared_sequences < 1:
        raise ValueError("min_shared_sequences must be at least 1.")

    stat_cols = [
        "new_downstream_burden",
        "supported_new_downstream_burden",
        "new_downstream_children",
        "supported_new_downstream_children",
        "mean_successor_new_sequences",
    ]

    required_edge_cols = {"source", "target", "n_shared_sequences"}
    missing_edge_cols = sorted(required_edge_cols - set(edge_table.columns))
    if missing_edge_cols:
        raise ValueError(f"Missing edge columns: {missing_edge_cols}")

    required_sequence_cols = {"cluster_id", "sequence_id"}
    missing_sequence_cols = sorted(required_sequence_cols - set(sequence_df.columns))
    if missing_sequence_cols:
        raise ValueError(f"Missing sequence columns: {missing_sequence_cols}")

    if edge_table.empty:
        return pd.DataFrame(columns=["cluster_id", *stat_cols])

    membership = sequence_df[["cluster_id", "sequence_id"]].dropna().drop_duplicates()
    cluster_sequences = (
        membership.groupby("cluster_id", sort=False)["sequence_id"]
        .agg(lambda values: frozenset(values))
        .to_dict()
    )

    rows = []
    for source, edges in edge_table.groupby("source", sort=False, dropna=False):
        source_sequences = cluster_sequences.get(source, frozenset())
        new_sequences: set = set()
        supported_new_sequences: set = set()
        child_new_counts: list[int] = []
        new_child_count = 0
        supported_new_child_count = 0

        for edge in edges.itertuples(index=False):
            target_sequences = cluster_sequences.get(edge.target, frozenset())
            child_new_sequences = target_sequences.difference(source_sequences)
            child_new_count = len(child_new_sequences)
            child_new_counts.append(child_new_count)

            if child_new_count > 0:
                new_child_count += 1
                new_sequences.update(child_new_sequences)

            if edge.n_shared_sequences >= min_shared_sequences:  # type: ignore
                if child_new_count > 0:
                    supported_new_child_count += 1
                    supported_new_sequences.update(child_new_sequences)

        rows.append(
            {
                "cluster_id": source,
                "new_downstream_burden": len(new_sequences),
                "supported_new_downstream_burden": len(supported_new_sequences),
                "new_downstream_children": new_child_count,
                "supported_new_downstream_children": supported_new_child_count,
                "mean_successor_new_sequences": (
                    float(np.mean(child_new_counts)) if child_new_counts else np.nan
                ),
            }
        )

    return pd.DataFrame(rows, columns=["cluster_id", *stat_cols])


def add_new_downstream_metrics(
    node_table: pd.DataFrame,
    sequence_df: pd.DataFrame,
    edge_table: pd.DataFrame,
    *,
    min_shared_sequences: int = 2,
) -> pd.DataFrame:
    """Attach overlap-adjusted downstream burden metrics to a node table."""
    new_downstream = build_new_downstream_metrics(
        sequence_df,
        edge_table,
        min_shared_sequences=min_shared_sequences,
    )
    out = node_table.merge(new_downstream, on="cluster_id", how="left")
    fill_zero_cols = [
        "new_downstream_burden",
        "supported_new_downstream_burden",
        "new_downstream_children",
        "supported_new_downstream_children",
    ]
    for col in fill_zero_cols:
        out[col] = out[col].fillna(0)
    return out


def add_graph_role_indicators(
    node_table: pd.DataFrame,
    *,
    window_col: str = "window_idx",
    add_censoring: bool = True,
) -> pd.DataFrame:
    """Attach graph-role flags and optional boundary-censoring flags."""
    out = node_table.copy()
    in_degree = out["in_degree"].fillna(0)
    out_degree = out["out_degree"].fillna(0)
    out["has_incoming"] = in_degree.gt(0)
    out["has_outgoing"] = out_degree.gt(0)
    out["is_source_boundary"] = in_degree.eq(0)
    out["is_sink_boundary"] = out_degree.eq(0)
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

    if add_censoring and window_col in out.columns:
        first_window = out[window_col].min()
        last_window = out[window_col].max()
        out["left_censored"] = out["is_source_boundary"] & out[window_col].eq(
            first_window
        )
        out["right_censored"] = out["is_sink_boundary"] & out[window_col].eq(
            last_window
        )
        out["boundary_censored"] = out["left_censored"] | out["right_censored"]
    return out


def build_transition_node_table(
    cluster_table: pd.DataFrame,
    edge_table: pd.DataFrame,
) -> pd.DataFrame:
    """Attach transition role/component summaries to an SSE cluster table copy."""
    required = {"cluster_id", "window_id", "window_idx", "cluster_size"}
    missing = required - set(cluster_table.columns)
    if missing:
        raise KeyError(f"Missing transition node columns: {sorted(missing)}")

    out = cluster_table.copy()
    flow_cols = {
        "out_degree",
        "out_strength",
        "onward_entropy",
        "onward_entropy_norm",
        "effective_successors",
        "dominant_successor_frac",
        "in_degree",
        "in_strength",
    }
    if flow_cols - set(out.columns):
        out = add_edge_flow_metrics(out, edge_table)
    if {"downstream_cluster_burden", "mean_successor_cluster_size"} - set(out.columns):
        out = add_downstream_burden(out, edge_table)
    for col in ("out_degree", "out_strength", "in_degree", "in_strength"):
        if col not in out.columns:
            out[col] = 0
        else:
            out[col] = out[col].fillna(0)
    out = add_graph_role_indicators(out)

    graph = _graph_from_node_edge_tables(out["cluster_id"], edge_table)
    node_to_component, component_size = component_maps(graph)
    out["weak_component_id"] = out["cluster_id"].map(node_to_component)
    out["weak_component_size"] = out["weak_component_id"].map(component_size)
    return out


def build_graph_summary(
    node_table: pd.DataFrame,
    edge_table: pd.DataFrame,
) -> pd.DataFrame:
    """Return scalar descriptive summaries for the transition graph."""
    total_weight = (
        edge_table["n_shared_sequences"].sum()
        if "n_shared_sequences" in edge_table.columns
        else 0
    )
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
        ("total_shared_sequence_weight", total_weight),
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
    agg = {
        "n_nodes": ("cluster_id", "nunique"),
        "first_window_idx": ("window_idx", "min"),
        "last_window_idx": ("window_idx", "max"),
        "n_windows": ("window_id", "nunique"),
        "total_cluster_size": ("cluster_size", "sum"),
        "max_cluster_size": ("cluster_size", "max"),
    }
    if "modal_health_board" in node_table.columns:
        agg["n_health_boards"] = ("modal_health_board", "nunique")
    return (
        node_table.groupby("weak_component_id", dropna=False)
        .agg(**agg)
        .reset_index()
        .sort_values(["n_nodes", "total_cluster_size"], ascending=[False, False])
    )


def build_transition_summary_outputs(
    cluster_table: pd.DataFrame,
    edge_table: pd.DataFrame,
) -> TransitionSummaryOutputs:
    """Build descriptive transition graph outputs from SSE detector tables."""
    node_table = build_transition_node_table(cluster_table, edge_table)
    return TransitionSummaryOutputs(
        node_table=node_table,
        graph_summary=build_graph_summary(node_table, edge_table),
        window_summary=build_transition_window_summary(node_table, edge_table),
        component_summary=build_component_summary(node_table),
    )
