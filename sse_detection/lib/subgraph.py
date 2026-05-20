"""Plot the induced subgraph of a single meta-cluster in a dot-style layout.

Layout
------
The graph is a DAG whose edges only span adjacent windows by construction
(see ``sse_detection.ipynb``). That makes it a natural fit for the layered
("Sugiyama") layout that Graphviz's ``dot`` produces:

* Each node's *rank* is its ``window_idx``.
* Within a rank, nodes are ordered to minimise edge crossings using a
  bidirectional barycentre sweep over the bipartite layers.
* Edges are drawn as curved arrows running left-to-right (or top-to-bottom).

If ``pygraphviz`` is importable in the active environment, Graphviz's ``dot``
engine is used through it. If not, the ``dot`` binary / ``pydot`` are tried
before falling back to a pure-Python Sugiyama layout. Both produce visually
similar layouts.

Encodings
---------
* Node colour: ``sse_category`` (the role-derived palette in
  :mod:`palettes`). Also accepts ``"sse_role"`` and ``"sse_onward_dynamic"``.
* Node size: proportional to ``log1p(cluster_size)``.
* Edge width: proportional to ``log1p(n_shared_sequences)``.
"""

from __future__ import annotations

import shutil
import subprocess
from typing import Iterable

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch

from utils import style

from .palettes import (
    DYNAMIC_ORDER,
    DYNAMIC_PALETTE,
    NOT_SSE_COLOR,
    ROLE_ORDER,
    ROLE_PALETTE,
    sse_category_palette_from,
)


# ---------------------------------------------------------------------------
# Subgraph selection
# ---------------------------------------------------------------------------


def _select_subgraph(
    node_stats: pd.DataFrame,
    edge_table: pd.DataFrame,
    meta_cluster_id: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    _require_columns(
        node_stats,
        {
            "cluster_id",
            "meta_cluster_id",
            "window_idx",
            "wn_mid_date",
            "cluster_size",
        },
        "node_stats",
    )
    _require_columns(edge_table, {"source", "target", "n_shared_sequences"}, "edge_table")

    nodes = node_stats.loc[node_stats["meta_cluster_id"] == meta_cluster_id].copy()
    if nodes.empty:
        raise ValueError(f"No nodes found for meta_cluster_id={meta_cluster_id!r}")
    if nodes["window_idx"].isna().any():
        raise ValueError("window_idx values must be non-null for subgraph layout")
    if nodes["cluster_id"].duplicated().any():
        duplicated = nodes.loc[nodes["cluster_id"].duplicated(), "cluster_id"].iloc[0]
        raise ValueError(f"cluster_id values must be unique; found duplicate {duplicated!r}")
    node_ids = set(nodes["cluster_id"])
    edges = edge_table.loc[
        edge_table["source"].isin(node_ids) & edge_table["target"].isin(node_ids)
    ].copy()
    return nodes, edges


def _require_columns(
    df: pd.DataFrame,
    required: set[str],
    name: str,
) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        cols = ", ".join(repr(c) for c in missing)
        raise KeyError(f"{name} must include required column(s): {cols}")


# ---------------------------------------------------------------------------
# Layered layout (dot-style)
# ---------------------------------------------------------------------------


def _initial_order(
    nodes: pd.DataFrame,
    rank_col: str = "window_idx",
) -> dict[int, list[str]]:
    """Deterministic initial within-rank order: largest cluster_size first."""
    n = nodes.sort_values(
        [rank_col, "cluster_size", "cluster_id"],
        ascending=[True, False, True],
    )
    return {
        int(rank): sub["cluster_id"].tolist()
        for rank, sub in n.groupby(rank_col)
    }


def _barycentre_sweep(
    order_by_rank: dict[int, list[str]],
    successors: dict[str, list[str]],
    predecessors: dict[str, list[str]],
    *,
    n_sweeps: int = 24,
) -> dict[int, list[str]]:
    """Reduce edge crossings via a bidirectional barycentre heuristic."""
    ranks = sorted(order_by_rank.keys())
    if len(ranks) < 2:
        return order_by_rank

    def _sort_layer(layer_nodes, neighbour_pos, neighbour_map):
        def key(n):
            ns = neighbour_map.get(n, [])
            positions = [neighbour_pos[m] for m in ns if m in neighbour_pos]
            if not positions:
                # Keep node at its current position by returning a neutral key.
                return float("inf")
            return sum(positions) / len(positions)
        return sorted(layer_nodes, key=key)

    for _ in range(n_sweeps):
        # Forward sweep: order by the previous present rank, so skipped
        # window_idx values do not break ordering.
        for prev_rank, r in zip(ranks[:-1], ranks[1:]):
            prev_pos = {n: i for i, n in enumerate(order_by_rank[prev_rank])}
            order_by_rank[r] = _sort_layer(order_by_rank[r], prev_pos, predecessors)
        # Backward sweep: order by the next present rank.
        for r, next_rank in zip(reversed(ranks[:-1]), reversed(ranks[1:])):
            next_pos = {n: i for i, n in enumerate(order_by_rank[next_rank])}
            order_by_rank[r] = _sort_layer(order_by_rank[r], next_pos, successors)
    return order_by_rank


def _sugiyama_positions(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    rank_col: str = "window_idx",
    n_sweeps: int = 24,
) -> dict[str, tuple[float, float]]:
    """Layered layout. Returns rank-space positions (x = layer, y = order)."""
    order_by_rank = _initial_order(nodes, rank_col=rank_col)

    successors: dict[str, list[str]] = {nid: [] for nid in nodes["cluster_id"]}
    predecessors: dict[str, list[str]] = {nid: [] for nid in nodes["cluster_id"]}
    for s, t in zip(edges["source"], edges["target"]):
        if s in successors and t in predecessors:
            successors[s].append(t)
            predecessors[t].append(s)

    order_by_rank = _barycentre_sweep(order_by_rank, successors, predecessors, n_sweeps=n_sweeps)

    # Compress ranks to consecutive indices so x-axis is dot-style equi-spaced.
    sorted_ranks = sorted(order_by_rank.keys())
    rank_to_x = {r: i for i, r in enumerate(sorted_ranks)}

    pos: dict[str, tuple[float, float]] = {}
    for rank, layer in order_by_rank.items():
        n = len(layer)
        # Centre the layer vertically around 0; use unit spacing.
        for i, nid in enumerate(layer):
            y = (i - (n - 1) / 2.0)
            pos[nid] = (rank_to_x[rank], y)
    return pos


def _try_graphviz_positions(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    rankdir: str,
    rank_col: str = "window_idx",
) -> dict[str, tuple[float, float]] | None:
    """Use Graphviz's ``dot`` engine when the binary or Python wrappers exist."""
    id_by_cluster, anchor_by_rank = _graphviz_ids(nodes, rank_col=rank_col)

    pos = _try_pygraphviz_positions(nodes, edges, id_by_cluster, anchor_by_rank, rankdir)
    if pos is not None:
        return pos

    pos = _try_graphviz_binary_positions(nodes, edges, id_by_cluster, anchor_by_rank, rankdir)
    if pos is not None:
        return pos

    return _try_pydot_positions(nodes, edges, id_by_cluster, anchor_by_rank, rankdir)


def _graphviz_ids(
    nodes: pd.DataFrame,
    *,
    rank_col: str,
) -> tuple[dict[str, str], dict[int, str]]:
    id_by_cluster = {
        cluster_id: f"n{i}"
        for i, cluster_id in enumerate(nodes["cluster_id"])
    }
    ranks = sorted(int(rank) for rank in nodes[rank_col].dropna().unique())
    anchor_by_rank = {rank: f"rank_anchor_{i}" for i, rank in enumerate(ranks)}
    return id_by_cluster, anchor_by_rank


def _build_dot_source(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    id_by_cluster: dict[str, str],
    anchor_by_rank: dict[int, str],
    rankdir: str,
    *,
    rank_col: str = "window_idx",
) -> str:
    lines = [
        "digraph G {",
        f'  graph [rankdir="{rankdir}", nodesep="0.3", ranksep="0.6"];',
        '  node [shape=point, width="0.08", height="0.08", label=""];',
        '  edge [arrowsize="0.7"];',
    ]

    for graph_id in id_by_cluster.values():
        lines.append(f"  {graph_id};")
    for anchor_id in anchor_by_rank.values():
        lines.append(
            f'  {anchor_id} [style="invis", width="0.01", height="0.01", label=""];'
        )

    for i, (rank, sub) in enumerate(nodes.groupby(rank_col, sort=True)):
        rank_id = int(rank)
        members = [anchor_by_rank[rank_id]]
        members.extend(id_by_cluster[c] for c in sub["cluster_id"])
        lines.append(f"  subgraph rank_{i} {{")
        lines.append("    rank=same;")
        lines.extend(f"    {member};" for member in members)
        lines.append("  }")

    anchors = [anchor_by_rank[rank] for rank in sorted(anchor_by_rank)]
    for left, right in zip(anchors[:-1], anchors[1:]):
        lines.append(f'  {left} -> {right} [style="invis", weight="100"];')

    for src, tgt in zip(edges["source"], edges["target"]):
        if src in id_by_cluster and tgt in id_by_cluster:
            lines.append(f"  {id_by_cluster[src]} -> {id_by_cluster[tgt]};")

    lines.append("}")
    return "\n".join(lines)


def _try_graphviz_binary_positions(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    id_by_cluster: dict[str, str],
    anchor_by_rank: dict[int, str],
    rankdir: str,
) -> dict[str, tuple[float, float]] | None:
    dot = shutil.which("dot")
    if dot is None:
        return None

    dot_source = _build_dot_source(nodes, edges, id_by_cluster, anchor_by_rank, rankdir)
    try:
        result = subprocess.run(
            [dot, "-Tplain"],
            input=dot_source,
            capture_output=True,
            check=False,
            text=True,
        )
    except OSError:
        return None

    if result.returncode != 0:
        return None
    return _parse_dot_plain_positions(result.stdout, id_by_cluster)


def _try_pygraphviz_positions(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    id_by_cluster: dict[str, str],
    anchor_by_rank: dict[int, str],
    rankdir: str,
    *,
    rank_col: str = "window_idx",
) -> dict[str, tuple[float, float]] | None:
    try:
        import pygraphviz as pgv
    except ImportError:
        return None

    try:
        graph = pgv.AGraph(
            directed=True,
            strict=False,
            rankdir=rankdir,
            nodesep="0.3",
            ranksep="0.6",
        )
        for graph_id in id_by_cluster.values():
            graph.add_node(graph_id, shape="point", width="0.08", height="0.08", label="")
        for anchor_id in anchor_by_rank.values():
            graph.add_node(
                anchor_id,
                style="invis",
                width="0.01",
                height="0.01",
                label="",
            )
        for _, sub in nodes.groupby(rank_col, sort=True):
            rank_id = int(sub[rank_col].iloc[0])
            members = [anchor_by_rank[rank_id]]
            members.extend(id_by_cluster[c] for c in sub["cluster_id"])
            graph.add_subgraph(members, rank="same")
        anchors = [anchor_by_rank[rank] for rank in sorted(anchor_by_rank)]
        for left, right in zip(anchors[:-1], anchors[1:]):
            graph.add_edge(left, right, style="invis", weight="100")
        for src, tgt in zip(edges["source"], edges["target"]):
            if src in id_by_cluster and tgt in id_by_cluster:
                graph.add_edge(id_by_cluster[src], id_by_cluster[tgt])
        graph.layout(prog="dot")
        return _parse_pygraphviz_positions(graph, id_by_cluster)
    except Exception:
        return None


def _parse_pygraphviz_positions(
    graph,
    id_by_cluster: dict[str, str],
) -> dict[str, tuple[float, float]] | None:
    pos: dict[str, tuple[float, float]] = {}
    for cluster_id, graph_id in id_by_cluster.items():
        raw = graph.get_node(graph_id).attr.get("pos")
        if not raw:
            return None
        x, y = raw.split(",", 1)
        pos[cluster_id] = (float(x), float(y))
    return pos


def _try_pydot_positions(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    id_by_cluster: dict[str, str],
    anchor_by_rank: dict[int, str],
    rankdir: str,
    *,
    rank_col: str = "window_idx",
) -> dict[str, tuple[float, float]] | None:
    try:
        import pydot
    except ImportError:
        return None

    try:
        graph = pydot.Dot(
            graph_type="digraph",
            rankdir=rankdir,
            nodesep="0.3",
            ranksep="0.6",
        )
        for graph_id in id_by_cluster.values():
            graph.add_node(
                pydot.Node(graph_id, shape="point", width="0.08", height="0.08", label="")
            )
        for anchor_id in anchor_by_rank.values():
            graph.add_node(
                pydot.Node(
                    anchor_id,
                    style="invis",
                    width="0.01",
                    height="0.01",
                    label="",
                )
            )
        for i, (_, sub) in enumerate(nodes.groupby(rank_col, sort=True)):
            rank_id = int(sub[rank_col].iloc[0])
            subgraph = pydot.Subgraph(graph_name=f"rank_{i}", rank="same")
            subgraph.add_node(pydot.Node(anchor_by_rank[rank_id]))
            for cluster_id in sub["cluster_id"]:
                subgraph.add_node(pydot.Node(id_by_cluster[cluster_id]))
            graph.add_subgraph(subgraph)
        anchors = [anchor_by_rank[rank] for rank in sorted(anchor_by_rank)]
        for left, right in zip(anchors[:-1], anchors[1:]):
            graph.add_edge(pydot.Edge(left, right, style="invis", weight="100"))
        for src, tgt in zip(edges["source"], edges["target"]):
            if src in id_by_cluster and tgt in id_by_cluster:
                graph.add_edge(pydot.Edge(id_by_cluster[src], id_by_cluster[tgt]))

        plain = graph.create_plain(prog="dot")
        if isinstance(plain, bytes):
            plain = plain.decode("utf-8")
        return _parse_dot_plain_positions(plain, id_by_cluster)
    except Exception:
        return None


def _parse_dot_plain_positions(
    plain: str,
    id_by_cluster: dict[str, str],
) -> dict[str, tuple[float, float]] | None:
    cluster_by_id = {graph_id: cluster_id for cluster_id, graph_id in id_by_cluster.items()}
    pos: dict[str, tuple[float, float]] = {}
    for line in plain.splitlines():
        parts = line.split()
        if len(parts) < 4 or parts[0] != "node":
            continue
        graph_id = parts[1]
        if graph_id in cluster_by_id:
            pos[cluster_by_id[graph_id]] = (float(parts[2]), float(parts[3]))
    if len(pos) != len(id_by_cluster):
        return None
    return pos


def _dot_layout(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    layout: str = "auto",
    rankdir: str = "LR",
    n_sweeps: int = 24,
) -> tuple[dict[str, tuple[float, float]], str]:
    """Compute dot-style positions for the subgraph.

    Returns
    -------
    pos : dict
        Mapping ``cluster_id -> (x, y)``.
    engine : str
        Which engine was used (``"graphviz"`` or ``"sugiyama"``).
    """
    if layout not in {"auto", "dot", "sugiyama"}:
        raise ValueError("layout must be one of {'auto','dot','sugiyama'}")
    if rankdir not in {"LR", "TB"}:
        raise ValueError("rankdir must be one of {'LR','TB'}")

    if layout in {"auto", "dot"}:
        pos = _try_graphviz_positions(nodes, edges, rankdir=rankdir)
        if pos is not None:
            return pos, "graphviz"
        if layout == "dot":
            raise RuntimeError(
                "layout='dot' requires the Graphviz dot binary, pygraphviz, "
                "or pydot to be available in the active environment."
            )

    pos = _sugiyama_positions(nodes, edges, n_sweeps=n_sweeps)
    if rankdir == "TB":
        # Swap axes: (x, y) -> (y, -x) so layers stack top-to-bottom.
        pos = {nid: (y, -x) for nid, (x, y) in pos.items()}
    return pos, "sugiyama"


# ---------------------------------------------------------------------------
# Visual encoding helpers
# ---------------------------------------------------------------------------


def _node_sizes(
    cluster_sizes: np.ndarray,
    *,
    size_min: float = 30,
    size_max: float = 800,
) -> np.ndarray:
    log_sizes = np.log1p(np.asarray(cluster_sizes, dtype=float).clip(min=0))
    if log_sizes.max() == log_sizes.min():
        return np.full_like(log_sizes, (size_min + size_max) / 2.0)
    norm = (log_sizes - log_sizes.min()) / (log_sizes.max() - log_sizes.min())
    return size_min + norm * (size_max - size_min)


def _edge_widths(
    weights: np.ndarray,
    *,
    width_min: float = 0.4,
    width_max: float = 4.0,
) -> np.ndarray:
    log_w = np.log1p(np.asarray(weights, dtype=float).clip(min=0))
    if log_w.max() == 0:
        return np.full_like(log_w, width_min)
    return width_min + (log_w / log_w.max()) * (width_max - width_min)


def _candidate_mask(layout: pd.DataFrame) -> pd.Series:
    return layout["sse_candidate"].fillna(False).astype(bool)


def _ordered_values(values: Iterable[str], order: Iterable[str]) -> list[str]:
    present = set(values)
    ordered = [v for v in order if v in present]
    ordered.extend(sorted(present - set(ordered)))
    return ordered


def _ordered_categories(values: Iterable[str]) -> list[str]:
    role_index = {role: i for i, role in enumerate(ROLE_ORDER)}
    dynamic_index = {dyn: i for i, dyn in enumerate(DYNAMIC_ORDER)}

    def key(category: str) -> tuple[int, int, str]:
        if category == "not_sse_like":
            return (len(role_index), len(dynamic_index), category)
        if "__" not in category:
            return (len(role_index) + 1, len(dynamic_index) + 1, category)
        role, dynamic = category.split("__", 1)
        return (
            role_index.get(role, len(role_index) + 1),
            dynamic_index.get(dynamic, len(dynamic_index) + 1),
            category,
        )

    return sorted(set(values), key=key)


def _resolve_colours(
    layout: pd.DataFrame,
    color_by: str,
) -> tuple[pd.Series, dict]:
    if color_by == "sse_category":
        _require_columns(layout, {"sse_category"}, "node_stats")
        cats = _ordered_categories(
            set(layout["sse_category"].dropna().unique()) | {"not_sse_like"}
        )
        palette = sse_category_palette_from(cats)
        colours = layout["sse_category"].map(palette).fillna(NOT_SSE_COLOR)
    elif color_by == "sse_role":
        _require_columns(layout, {"sse_role", "sse_candidate"}, "node_stats")
        palette = {**ROLE_PALETTE, "not_sse_like": NOT_SSE_COLOR}
        role_or_neutral = layout["sse_role"].where(_candidate_mask(layout), "not_sse_like")
        colours = role_or_neutral.map(palette).fillna(NOT_SSE_COLOR)
    elif color_by == "sse_onward_dynamic":
        _require_columns(layout, {"sse_onward_dynamic", "sse_candidate"}, "node_stats")
        uniq = _ordered_values(layout["sse_onward_dynamic"].dropna().unique(), DYNAMIC_ORDER)
        palette = {v: DYNAMIC_PALETTE[v] for v in uniq if v in DYNAMIC_PALETTE}
        unknown = [v for v in uniq if v not in palette]
        if unknown:
            cmap = plt.get_cmap("tab20")
            palette.update({v: cmap(i % 20) for i, v in enumerate(unknown)})
        palette["not_sse_like"] = NOT_SSE_COLOR
        dyn = layout["sse_onward_dynamic"].where(_candidate_mask(layout), "not_sse_like")
        colours = dyn.map(palette).fillna(NOT_SSE_COLOR)
    else:
        raise ValueError(
            "color_by must be one of {'sse_category','sse_role','sse_onward_dynamic'}"
        )
    return colours, palette


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def plot_meta_cluster_subgraph(
    node_stats: pd.DataFrame,
    edge_table: pd.DataFrame,
    meta_cluster_id: str,
    *,
    color_by: str = "sse_category",
    annotate_top_n: int = 6,
    annotate_col: str = "cluster_size",
    show_legend: bool = True,
    layout: str = "auto",
    rankdir: str = "LR",
    curvature: float = 0.18,
    show_arrows: bool = True,
    show_window_axis: bool = True,
    title: str | None = None,
    width: str = "slide",
    height_in: float | None = None,
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """Draw the induced subgraph of a single meta-cluster in a dot-tree layout.

    Parameters
    ----------
    node_stats
        Node-level dataframe produced by ``sse_detection.ipynb``. Must
        include ``cluster_id``, ``meta_cluster_id``, ``window_idx``,
        ``wn_mid_date``, ``cluster_size``, and the ``sse_*`` label columns.
    edge_table
        Edge dataframe with ``source``, ``target``, ``n_shared_sequences``.
    meta_cluster_id
        ID of the meta-cluster to plot (e.g. ``"AM00001"``). Pass a single
        ID; the function intentionally restricts itself to one cluster at a
        time for a clean tree layout. Use a loop in the caller to render
        several.
    color_by
        Node attribute used to colour points. ``"sse_category"`` (default),
        ``"sse_role"``, or ``"sse_onward_dynamic"``.
    annotate_top_n
        Number of largest nodes to label with their ``cluster_id`` suffix.
    annotate_col
        Column used to pick which nodes to annotate (default
        ``"cluster_size"``).
    show_legend
        Add a colour legend below the axes.
    layout
        ``"auto"`` (try Graphviz dot, fall back to Sugiyama), ``"dot"``
        (require Graphviz dot), or ``"sugiyama"`` (force pure-Python).
    rankdir
        ``"LR"`` (left to right; default) or ``"TB"`` (top to bottom).
        Mirrors Graphviz's ``rankdir`` attribute.
    curvature
        Bezier curvature for edges (matplotlib ``arc3,rad=`` value). Pass 0
        for straight lines.
    show_arrows
        Draw arrowheads on edges (default True).
    show_window_axis
        If True (default and ``rankdir='LR'``), the x-axis is labelled with
        ``window_idx`` values and ``wn_mid_date`` underneath.
    title
        Optional figure title.
    width, height_in
        Forwarded to ``utils.style.new_figure``. ``height_in`` defaults to
        a value scaled by the widest layer when omitted.
    ax
        Pre-existing axis to draw onto. If provided, ``width``/``height_in``
        are ignored.

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the subgraph.
    """
    if not isinstance(meta_cluster_id, str):
        raise TypeError(
            "meta_cluster_id must be a single string; the plotter intentionally "
            "handles one meta-cluster per call. Loop in the caller to produce "
            "multiple figures."
        )

    nodes, edges = _select_subgraph(node_stats, edge_table, meta_cluster_id)

    pos, engine = _dot_layout(nodes, edges, layout=layout, rankdir=rankdir)

    # Decorate the nodes dataframe with positions for downstream lookups.
    nodes = nodes.copy()
    nodes["x_pos"] = nodes["cluster_id"].map(lambda c: pos[c][0])
    nodes["y_pos"] = nodes["cluster_id"].map(lambda c: pos[c][1])

    colours, palette = _resolve_colours(nodes, color_by)
    sizes = pd.Series(_node_sizes(nodes["cluster_size"].to_numpy()), index=nodes.index)

    if ax is None:
        if height_in is None:
            # Scale height by the widest layer so dense layers stay legible.
            widest = nodes.groupby("window_idx")["cluster_id"].count().max()
            height_in = float(np.clip(3.5 + 0.18 * widest, 4.0, 9.0))
        fig, ax = style.new_figure(
            width=width,
            height_in=height_in,
            context="talk",
            font_scale=0.78,
        )
    else:
        fig = ax.figure

    # ----- Edges -----
    if not edges.empty:
        widths = _edge_widths(edges["n_shared_sequences"].to_numpy())
        arrowstyle = "-|>" if show_arrows else "-"
        rad = float(curvature)
        for src, tgt, w in zip(edges["source"], edges["target"], widths):
            if src not in pos or tgt not in pos:
                continue
            arrow = FancyArrowPatch(
                posA=pos[src],
                posB=pos[tgt],
                connectionstyle=f"arc3,rad={rad}",
                arrowstyle=arrowstyle,
                mutation_scale=10 if show_arrows else 0,
                linewidth=float(w),
                color="#9c9c9c",
                alpha=0.7,
                shrinkA=6,
                shrinkB=6,
                zorder=1,
            )
            ax.add_patch(arrow)

    # ----- Nodes -----
    ax.scatter(
        nodes["x_pos"], nodes["y_pos"],
        s=sizes.to_numpy(),
        c=list(colours),
        alpha=0.94,
        edgecolor="white",
        linewidth=0.7,
        zorder=2,
    )

    # Black outline for SSE candidates so the colour-as-category encoding
    # remains readable on the small-size end of the distribution.
    if "sse_candidate" in nodes.columns:
        cand = nodes.loc[_candidate_mask(nodes)]
        if not cand.empty:
            ax.scatter(
                cand["x_pos"], cand["y_pos"],
                s=sizes.loc[cand.index].to_numpy(),
                facecolor="none",
                edgecolor="#14151F",
                linewidth=0.8,
                zorder=3,
            )

    # ----- Annotations -----
    if annotate_top_n > 0 and annotate_col in nodes.columns:
        top = nodes.nlargest(annotate_top_n, annotate_col)
        for row in top.itertuples(index=False):
            label = row.cluster_id.split("|")[-1]
            ax.annotate(
                label,
                xy=(row.x_pos, row.y_pos),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=7,
                color="#14151F",
                alpha=0.9,
            )

    # ----- Axes / chrome -----
    _format_axes(ax, nodes, rankdir, show_window_axis)
    ax.set_title(title or _default_title(meta_cluster_id, color_by, engine))

    if show_legend:
        _add_color_legend(ax, color_by, nodes, palette)

    return fig


# ---------------------------------------------------------------------------
# Labels and legend
# ---------------------------------------------------------------------------


def _format_axes(
    ax: plt.Axes,
    nodes: pd.DataFrame,
    rankdir: str,
    show_window_axis: bool,
) -> None:
    ax.set_yticks([])
    ax.set_xticks([])

    # Add a strip of date labels for the layered axis (window_idx ranks).
    if show_window_axis and rankdir == "LR":
        per_rank = (
            nodes.groupby("window_idx")
            .agg(_x=("x_pos", "first"), wn_mid_date=("wn_mid_date", "first"))
            .reset_index()
            .sort_values("_x")
        )
        ax.set_xticks(per_rank["_x"].to_numpy())
        many_windows = len(per_rank) > 10

        def tick_label(window_idx, date) -> str:
            prefix = f"W{int(window_idx)}"
            if pd.isna(date):
                return prefix
            date_label = pd.to_datetime(date).strftime("%Y-%m-%d")
            return f"{prefix} {date_label}" if many_windows else f"{prefix}\n{date_label}"

        labels = [
            tick_label(w, d)
            for w, d in zip(per_rank["window_idx"], per_rank["wn_mid_date"])
        ]
        ax.set_xticklabels(
            labels,
            fontsize=7,
            rotation=90 if many_windows else 0,
            ha="center",
            va="top",
        )
        ax.set_xlabel("window index / midpoint")
    elif show_window_axis and rankdir == "TB":
        # Y-axis becomes the rank axis under TB.
        per_rank = (
            nodes.groupby("window_idx")
            .agg(_y=("y_pos", "first"), wn_mid_date=("wn_mid_date", "first"))
            .reset_index()
            .sort_values("_y", ascending=False)
        )
        ax.set_yticks(per_rank["_y"].to_numpy())
        labels = [
            (
                f"W{int(w)} {pd.to_datetime(d).strftime('%Y-%m-%d')}"
                if pd.notna(d)
                else f"W{int(w)}"
            )
            for w, d in zip(per_rank["window_idx"], per_rank["wn_mid_date"])
        ]
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_ylabel("window index / midpoint")
    else:
        ax.set_xlabel("")
        ax.set_ylabel("")

    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(left=False, bottom=False)
    # Keep a little padding around the layout so nodes aren't clipped.
    xs = nodes["x_pos"].to_numpy()
    ys = nodes["y_pos"].to_numpy()
    if len(xs):
        pad_x = max(0.5, (xs.max() - xs.min()) * 0.04 + 0.5)
        pad_y = max(0.5, (ys.max() - ys.min()) * 0.06 + 0.5)
        ax.set_xlim(xs.min() - pad_x, xs.max() + pad_x)
        ax.set_ylim(ys.min() - pad_y, ys.max() + pad_y)


def _default_title(meta_cluster_id: str, color_by: str, engine: str) -> str:
    return f"Meta-cluster {meta_cluster_id} — dot tree ({engine}; colour: {color_by})"


def _add_color_legend(
    ax: plt.Axes,
    color_by: str,
    layout: pd.DataFrame,
    palette: dict,
) -> None:
    if color_by == "sse_category":
        cats: Iterable[str] = _ordered_categories(layout["sse_category"].dropna().unique())
        cats = list(cats)[:14]  # cap legend size
    elif color_by == "sse_role":
        cats = [
            r
            for r in ROLE_ORDER
            if r in set(layout["sse_role"].dropna().unique())
        ]
        if (~_candidate_mask(layout)).any():
            cats.append("not_sse_like")
    else:
        cats = _ordered_values(
            layout["sse_onward_dynamic"].dropna().unique(),
            DYNAMIC_ORDER,
        )
        if (~_candidate_mask(layout)).any():
            cats.append("not_sse_like")

    handles = [
        mpatches.Patch(color=palette.get(c, NOT_SSE_COLOR), label=c)
        for c in cats
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=min(4, max(1, len(handles))),
        frameon=False,
        fontsize=7,
    )
