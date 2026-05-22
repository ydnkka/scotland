"""Plot the induced subgraph of a single meta-cluster in a dot-style layout.

Layout
------
The graph is a DAG whose edges only span adjacent windows by construction
(see ``sse_detection.ipynb``). That makes it a natural fit for Graphviz's
layered ``dot`` layout:

* Each node's *rank* is its ``window_idx``.
* Edges are drawn as curved arrows running left-to-right (or top-to-bottom).

Graphviz is used through ``pygraphviz``. The public ``layout="auto"`` option
is retained as a compatibility alias for ``layout="dot"``.

Encodings
---------
* Node colour: ``sse_category`` (the role-derived palette in
:mod:`palettes`). Also accepts ``"sse_role"`` and ``"sse_onward_dynamic"``.
* Node size: proportional to ``log1p(cluster_size)``.
* Edge width: proportional to ``log1p(n_shared_sequences)``.
"""

from __future__ import annotations

import textwrap
from typing import Iterable

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch

from utils.style import (
    new_figure,
    WIDTHS,
    CONTEXTS,
)

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
# Graphviz layout
# ---------------------------------------------------------------------------


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


def _pygraphviz_positions(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    id_by_cluster: dict[str, str],
    anchor_by_rank: dict[int, str],
    rankdir: str,
    *,
    rank_col: str = "window_idx",
) -> dict[str, tuple[float, float]]:
    try:
        import pygraphviz as pgv
    except ImportError as exc:
        raise RuntimeError(
            "pygraphviz is required for plot_meta_cluster_subgraph(). "
            "Install pygraphviz in the active Python environment."
        ) from exc

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


def _parse_pygraphviz_positions(
    graph,
    id_by_cluster: dict[str, str],
) -> dict[str, tuple[float, float]]:
    pos: dict[str, tuple[float, float]] = {}
    for cluster_id, graph_id in id_by_cluster.items():
        raw = graph.get_node(graph_id).attr.get("pos")
        if not raw:
            raise RuntimeError(f"Graphviz did not return a position for node {cluster_id!r}")
        x, y = str(raw).strip('"').split(",", 1)
        pos[cluster_id] = (float(x), float(y.rstrip("!")))
    return pos


def _dot_layout(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    rankdir: str = "LR",
) -> dict[str, tuple[float, float]]:
    """Compute dot-style positions for the subgraph.

    Returns
    -------
    pos : dict
        Mapping ``cluster_id -> (x, y)``.
    """
    if rankdir not in {"LR", "TB"}:
        raise ValueError("rankdir must be one of {'LR','TB'}")

    id_by_cluster, anchor_by_rank = _graphviz_ids(nodes, rank_col="window_idx")
    pos = _pygraphviz_positions(nodes, edges, id_by_cluster, anchor_by_rank, rankdir)
    return pos


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
    width_max: float = 5.0,
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
        palette  = sse_category_palette_from(cats)
        colours = layout["sse_category"].map(palette).fillna(NOT_SSE_COLOR)
    elif color_by == "sse_role":
        _require_columns(layout, {"sse_role", "sse_candidate"}, "node_stats")
        palette = {**ROLE_PALETTE, "not_sse_like": NOT_SSE_COLOR}
        role_or_neutral = layout["sse_role"].where(_candidate_mask(layout), "not_sse_like")
        colours = role_or_neutral.map(palette).fillna(NOT_SSE_COLOR)
    elif color_by == "sse_onward_dynamic":
        _require_columns(layout, {"sse_onward_dynamic", "sse_candidate"}, "node_stats")
        uniq = _ordered_values(layout["sse_onward_dynamic"].dropna().unique(), DYNAMIC_ORDER)
        palette: dict[str, object] = {v: DYNAMIC_PALETTE[v] for v in uniq if v in DYNAMIC_PALETTE}
        unknown = [v for v in uniq if v not in palette]
        if unknown:
            cmap = plt.get_cmap("tab20")
            for i, v in enumerate(unknown):
                palette[v] = cmap(i % 20)
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
    rankdir: str = "LR",
    curvature: float = 0.18,
    show_arrows: bool = True,
    show_window_axis: bool = True,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 5,
    context: CONTEXTS = "paper",
    font_scale: float = 1,
) -> Figure:
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
    width, width_in, height_in, context, font_scale
        Forwarded to ``utils.style.new_figure``.
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

    pos = _dot_layout(nodes, edges, rankdir=rankdir)

    # Decorate the nodes dataframe with positions for downstream lookups.
    nodes = nodes.copy()
    nodes["x_pos"] = nodes["cluster_id"].map(lambda c: pos[c][0])
    nodes["y_pos"] = nodes["cluster_id"].map(lambda c: pos[c][1])

    colours, palette = _resolve_colours(nodes, color_by)
    sizes = pd.Series(_node_sizes(nodes["cluster_size"].to_numpy()), index=nodes.index)

    # Scale height by the widest layer so dense layers stay legible.
    widest = nodes.groupby("window_idx")["cluster_id"].count().max()
    height_in = float(np.clip(3.5 + 0.18 * widest, 4.0, 9.0)) if height_in is None else height_in
    fig, ax = new_figure(
    width=width,
    width_in=width_in,
    height_in=height_in,
    context=context,
    font_scale=font_scale
    )

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
            label = str(row.cluster_id).split("|")[-1]
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

    if show_legend:
        _add_color_legend(
            ax,
            color_by,
            nodes,
            palette,
            rankdir=rankdir,
            show_window_axis=show_window_axis,
        )

    plt.close(fig)
    return fig


# ---------------------------------------------------------------------------
# Labels and legend
# ---------------------------------------------------------------------------


def _format_axes(
    ax: Axes,
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


def _legend_text_parts(category: str) -> list[str]:
    return [
        part.replace("_", " ").capitalize()
        for part in str(category).split("__", 1)
    ]


def _legend_label(category: str, *, line_width: int) -> str:
    lines: list[str] = []
    for part in _legend_text_parts(category):
        lines.extend(
            textwrap.wrap(part, width=line_width, break_long_words=False)
            or [part]
        )
    return "\n".join(lines)


def _legend_column_count(categories: list[str], fig_width_in: float) -> int:
    if not categories:
        return 1
    longest_line = max(
        len(part)
        for category in categories
        for part in _legend_text_parts(category)
    )
    target_col_width = min(2.0, max(1.3, 0.06 * longest_line + 0.35))
    width_limited_cols = max(1, int(fig_width_in // target_col_width))
    max_cols = 5 if fig_width_in >= 7 else 4
    max_rows = 4 if fig_width_in >= 7 else 5
    row_limited_cols = (len(categories) + max_rows - 1) // max_rows
    return min(
        len(categories),
        max(1, min(max_cols, max(width_limited_cols, row_limited_cols))),
    )


def _legend_line_width(fig_width_in: float, ncol: int) -> int:
    column_width = fig_width_in / max(1, ncol)
    return max(10, min(28, int((column_width - 0.35) / 0.06)))


def _legend_margin(
    labels: list[str],
    *,
    ncol: int,
    fig_height_in: float,
    has_bottom_axis: bool,
    many_windows: bool,
    font_size: float,
) -> tuple[float, float]:
    rows = (len(labels) + ncol - 1) // ncol
    max_lines = max(label.count("\n") + 1 for label in labels)
    line_height_in = font_size * 1.25 / 72.0
    legend_height_in = (
        rows * max_lines * line_height_in
        + 0.05 * max(0, rows - 1)
        + 0.12
    )
    legend_bottom_in = 0.05
    axis_space_in = 0.3
    if has_bottom_axis:
        legend_bottom_in = 0.08
        axis_space_in = 1.18 if many_windows else 0.72

    bottom = (
        legend_bottom_in
        + legend_height_in
        + axis_space_in
    ) / max(fig_height_in, 1.0)
    bottom = min(max(bottom, 0.2), 0.62)
    legend_top = max(
        0.02,
        (legend_bottom_in + legend_height_in) / max(fig_height_in, 1.0),
    )
    return bottom, legend_top


def _add_color_legend(
    ax: Axes,
    color_by: str,
    layout: pd.DataFrame,
    palette: dict,
    *,
    rankdir: str,
    show_window_axis: bool,
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

    cats = list(cats)
    fig_width_in, fig_height_in = ax.figure.get_size_inches() # type: ignore
    ncol = _legend_column_count(cats, fig_width_in)
    line_width = _legend_line_width(fig_width_in, ncol)
    labels = [_legend_label(c, line_width=line_width) for c in cats]
    handles = [
        mpatches.Patch(color=palette.get(c, NOT_SSE_COLOR), label=label)
        for c, label in zip(cats, labels)
    ]
    if not handles:
        return

    has_bottom_axis = show_window_axis and rankdir == "LR"
    many_windows = has_bottom_axis and layout["window_idx"].nunique() > 10
    bottom, legend_top = _legend_margin(
        labels,
        ncol=ncol,
        fig_height_in=fig_height_in,
        has_bottom_axis=has_bottom_axis,
        many_windows=many_windows,
        font_size=7,
    )

    if len(ax.figure.axes) == 1:
        ax.figure.subplots_adjust(bottom=bottom)
        bbox_to_anchor = (0.5, legend_top)
        bbox_transform = ax.figure.transFigure
    else:
        bbox_to_anchor = (0.5, -0.12)
        bbox_transform = ax.transAxes

    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=bbox_to_anchor,
        bbox_transform=bbox_transform,
        ncol=ncol,
        frameon=False,
        fontsize=7,
        columnspacing=0.9,
        handlelength=1.0,
        handletextpad=0.35,
        labelspacing=0.55,
        borderaxespad=0.0,
    )
