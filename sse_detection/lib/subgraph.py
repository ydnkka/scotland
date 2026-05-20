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

If ``pygraphviz`` or ``pydot`` is importable, Graphviz's ``dot`` engine is
used directly via ``networkx.drawing.nx_agraph`` / ``nx_pydot``; otherwise a
pure-Python Sugiyama fallback runs. Both produce visually similar layouts.

Encodings
---------
* Node colour: ``sse_category`` (the role-derived palette in
  :mod:`palettes`). Also accepts ``"sse_role"`` and ``"sse_onward_dynamic"``.
* Node size: proportional to ``log1p(cluster_size)``.
* Edge width: proportional to ``log1p(n_shared_sequences)``.
"""

from __future__ import annotations

from typing import Iterable

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch

from utils import style

from .palettes import (
    NOT_SSE_COLOR,
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
    if "meta_cluster_id" not in node_stats.columns:
        raise KeyError(
            "node_stats must include 'meta_cluster_id'. Run sse_detection.ipynb "
            "through the meta-cluster step before plotting."
        )
    nodes = node_stats.loc[node_stats["meta_cluster_id"] == meta_cluster_id].copy()
    if nodes.empty:
        raise ValueError(f"No nodes found for meta_cluster_id={meta_cluster_id!r}")
    node_ids = set(nodes["cluster_id"])
    edges = edge_table.loc[
        edge_table["source"].isin(node_ids) & edge_table["target"].isin(node_ids)
    ].copy()
    return nodes, edges


# ---------------------------------------------------------------------------
# Layered layout (dot-style)
# ---------------------------------------------------------------------------


def _initial_order(
    nodes: pd.DataFrame,
    rank_col: str = "window_idx",
) -> dict[int, list[str]]:
    """Deterministic initial within-rank order: largest cluster_size first."""
    n = nodes.sort_values([rank_col, "cluster_size", "cluster_id"], ascending=[True, False, True])
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
            if not ns:
                # Keep node at its current position by returning a neutral key.
                return float("inf")
            return sum(neighbour_pos.get(m, 0) for m in ns) / len(ns)
        return sorted(layer_nodes, key=key)

    for _ in range(n_sweeps):
        # Forward sweep: order rank r by mean position in rank r-1.
        for r in ranks[1:]:
            prev_pos = {n: i for i, n in enumerate(order_by_rank[r - 1] if (r - 1) in order_by_rank else [])}
            order_by_rank[r] = _sort_layer(order_by_rank[r], prev_pos, predecessors)
        # Backward sweep: order rank r by mean position in rank r+1.
        for r in reversed(ranks[:-1]):
            next_pos = {n: i for i, n in enumerate(order_by_rank[r + 1] if (r + 1) in order_by_rank else [])}
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
) -> dict[str, tuple[float, float]] | None:
    """Use Graphviz's ``dot`` engine if pygraphviz or pydot is available."""
    try:
        import networkx as nx  # noqa: F401
    except ImportError:
        return None

    try:
        import networkx as nx
        G = nx.DiGraph()
        G.add_nodes_from(nodes["cluster_id"])
        G.add_edges_from(zip(edges["source"], edges["target"]))
    except ImportError:
        return None

    # Encode rank constraints so dot respects window_idx ordering.
    # Graphviz reads node attribute 'rank' on subgraphs, not nodes; for plain
    # graphviz_layout the rank is inferred from edge direction, which is what
    # we want.
    args = f"-Grankdir={rankdir} -Gnodesep=0.3 -Granksep=0.6"

    # Try pygraphviz first.
    try:
        from networkx.drawing.nx_agraph import graphviz_layout
        return graphviz_layout(G, prog="dot", args=args)
    except (ImportError, Exception):
        pass

    # Fall back to pydot.
    try:
        from networkx.drawing.nx_pydot import graphviz_layout
        return graphviz_layout(G, prog="dot")
    except (ImportError, Exception):
        return None


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

    if layout in {"auto", "dot"}:
        pos = _try_graphviz_positions(nodes, edges, rankdir=rankdir)
        if pos is not None:
            return pos, "graphviz"
        if layout == "dot":
            raise RuntimeError(
                "layout='dot' requires pygraphviz or pydot + the graphviz "
                "binary to be installed."
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


def _resolve_colours(
    layout: pd.DataFrame,
    color_by: str,
) -> tuple[pd.Series, dict]:
    if color_by == "sse_category":
        cats = sorted(set(layout["sse_category"].dropna().unique()) | {"not_sse_like"})
        palette = sse_category_palette_from(cats)
        colours = layout["sse_category"].map(palette).fillna(NOT_SSE_COLOR)
    elif color_by == "sse_role":
        palette = {**ROLE_PALETTE, "not_sse_like": NOT_SSE_COLOR}
        role_or_neutral = layout["sse_role"].where(layout["sse_candidate"], "not_sse_like")
        colours = role_or_neutral.map(palette).fillna(NOT_SSE_COLOR)
    elif color_by == "sse_onward_dynamic":
        uniq = [v for v in layout["sse_onward_dynamic"].dropna().unique()]
        cmap = plt.get_cmap("tab20")
        palette = {v: cmap(i % 20) for i, v in enumerate(uniq)}
        palette["not_sse_like"] = NOT_SSE_COLOR
        dyn = layout["sse_onward_dynamic"].where(layout["sse_candidate"], "not_sse_like")
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
    sizes = _node_sizes(nodes["cluster_size"].to_numpy())

    if ax is None:
        if height_in is None:
            # Scale height by the widest layer so dense layers stay legible.
            widest = nodes.groupby("window_idx")["cluster_id"].count().max()
            height_in = float(np.clip(3.5 + 0.18 * widest, 4.0, 9.0))
        fig, ax = style.new_figure(width=width, height_in=height_in, context="talk", font_scale=0.78)
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
        s=sizes,
        c=list(colours),
        alpha=0.94,
        edgecolor="white",
        linewidth=0.7,
        zorder=2,
    )

    # Black outline for SSE candidates so the colour-as-category encoding
    # remains readable on the small-size end of the distribution.
    if "sse_candidate" in nodes.columns:
        cand = nodes.loc[nodes["sse_candidate"]]
        if not cand.empty:
            cand_sizes = _node_sizes(cand["cluster_size"].to_numpy())
            ax.scatter(
                cand["x_pos"], cand["y_pos"],
                s=cand_sizes,
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
    _format_axes(ax, nodes, engine, rankdir, show_window_axis)
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
    engine: str,
    rankdir: str,
    show_window_axis: bool,
) -> None:
    ax.set_yticks([])
    ax.set_xticks([])

    # Add a strip of date labels for the layered axis (window_idx ranks).
    if show_window_axis and engine == "sugiyama" and rankdir == "LR":
        per_rank = (
            nodes.groupby("window_idx")
            .agg(_x=("x_pos", "first"), wn_mid_date=("wn_mid_date", "first"))
            .reset_index()
            .sort_values("x_pos")
        )
        ax.set_xticks(per_rank["x_pos"].to_numpy())
        labels = [
            f"W{int(w)}\n{pd.to_datetime(d).strftime('%Y-%m-%d')}" if pd.notna(d) else f"W{int(w)}"
            for w, d in zip(per_rank["window_idx"], per_rank["wn_mid_date"])
        ]
        ax.set_xticklabels(labels, fontsize=7, rotation=0)
        ax.set_xlabel("window index / midpoint")
    elif show_window_axis and engine == "sugiyama" and rankdir == "TB":
        # Y-axis becomes the rank axis under TB.
        per_rank = (
            nodes.groupby("window_idx")
            .agg(_y=("y_pos", "first"), wn_mid_date=("wn_mid_date", "first"))
            .reset_index()
            .sort_values("y_pos", ascending=False)
        )
        ax.set_yticks(per_rank["y_pos"].to_numpy())
        labels = [
            f"W{int(w)} {pd.to_datetime(d).strftime('%Y-%m-%d')}" if pd.notna(d) else f"W{int(w)}"
            for w, d in zip(per_rank["window_idx"], per_rank["wn_mid_date"])
        ]
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_ylabel("window index / midpoint")
    else:
        # Graphviz coordinates are arbitrary; just drop the chrome.
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
        cats: Iterable[str] = sorted(layout["sse_category"].dropna().unique())
        cats = list(cats)[:14]  # cap legend size
    elif color_by == "sse_role":
        cats = [r for r in ROLE_PALETTE if r in set(layout["sse_role"].dropna().unique())]
        if (~layout["sse_candidate"]).any():
            cats.append("not_sse_like")
    else:
        cats = sorted([v for v in layout["sse_onward_dynamic"].dropna().unique()])
        if (~layout["sse_candidate"]).any():
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
