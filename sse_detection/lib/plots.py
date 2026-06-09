"""Plotting and table-formatting helpers for SSE detection outputs.

The functions in this module expect already-prepared analysis tables. They do
not load pipeline outputs directly; callers should use ``sse_detection.lib.io``
or notebook-side data preparation before plotting.
"""

from __future__ import annotations

import ast
import re
from typing import Any, Iterable, Mapping, Sequence

from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize, TwoSlopeNorm, to_rgb
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from shapely import wkb
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union

from utils import (
    set_theme,
    WIDTHS,
    CONTEXTS,
    new_figure,
    add_panel_labels,
    load_analysis_columns,
    CLADES,
    CLADE_PALETTE,
)

from .palettes import (
    BACKGROUND_COLOR,
    BACKGROUND_DARK,
    CANDIDATE_COLOR,
    CANDIDATE_DARK,
)
from .io import HIGH_PRIORITY_CANDIDATE_TIERS

BORDER = "#D0D7DE"
GRAY = "#808080"
GRAY_LIGHT = "#6C757D"
GRID = "#E9ECEF"
INK = "#212529"
INK_SOFT = "#495057"
ORANGE_DARK = "#D95F02"
OR_DIVERGING = plt.get_cmap("RdBu_r")
REFERENCE_COLOR = "#FF0000"
TEAL_DARK = "#007C89"
WARM_SEQ = plt.get_cmap("YlOrBr")

SSE_SIGNATURE_ORDER = ["burst", "burden", "burst+burden"]
SSE_SIGNATURE_PALETTE = {
    "burst": CANDIDATE_COLOR,
    "burden": TEAL_DARK,
    "burst+burden": ORANGE_DARK,
}


def _missing_optional_helper(name: str):
    def _missing(*args, **kwargs):
        raise ImportError(
            f"{name} is not available in this checkout; the optional report "
            "figure module that defines it is not present."
        )

    return _missing


try:
    from .sensitivity_figures import (
        make_clade_association_figures,
        make_clade_association_outputs,
        make_clade_association_summary_tables,
        make_sensitivity_analysis_figures,
        make_sensitivity_analysis_outputs,
        make_sensitivity_analysis_summary_tables,
    )
except ImportError:
    make_clade_association_figures = _missing_optional_helper(
        "make_clade_association_figures"
    )
    make_clade_association_outputs = _missing_optional_helper(
        "make_clade_association_outputs"
    )
    make_clade_association_summary_tables = _missing_optional_helper(
        "make_clade_association_summary_tables"
    )
    make_sensitivity_analysis_figures = _missing_optional_helper(
        "make_sensitivity_analysis_figures"
    )
    make_sensitivity_analysis_outputs = _missing_optional_helper(
        "make_sensitivity_analysis_outputs"
    )
    make_sensitivity_analysis_summary_tables = _missing_optional_helper(
        "make_sensitivity_analysis_summary_tables"
    )
try:
    from .policy_figures import make_policy_figures, plot_policy_report
except ImportError:
    make_policy_figures = _missing_optional_helper("make_policy_figures")
    plot_policy_report = _missing_optional_helper("plot_policy_report")

try:
    from .vaccination_figures import make_vaccination_figures, plot_vaccination_report
except ImportError:
    make_vaccination_figures = _missing_optional_helper("make_vaccination_figures")
    plot_vaccination_report = _missing_optional_helper("plot_vaccination_report")


__all__ = [
    "plot_cluster_size_distribution",
    "plot_meta_cluster_subgraph",
    "plot_role_dynamic_heatmap",
    "plot_candidate_rate_over_time",
    "plot_core_metric_space",
    "plot_composite_distributions",
    "plot_individual_categorical_distribution_bars",
    "plot_socio_demo_breakdown",
    "make_clade_association_figures",
    "make_clade_association_outputs",
    "make_clade_association_summary_tables",
    "make_sensitivity_analysis_figures",
    "make_sensitivity_analysis_outputs",
    "make_sensitivity_analysis_summary_tables",
    "make_policy_figures",
    "make_vaccination_figures",
    "plot_policy_report",
    "plot_vaccination_report",
    "plot_regression_wald_heatmap",
    "make_regression_wald_table",
    "make_regression_odds_ratio_table",
    "make_regression_fit_table",
    "load_health_board_geometries",
    "plot_health_board_enrichment_map",
    "plot_age_sex_simd_forest",
    "collect_sensitivity_matrix_results",
    "plot_sensitivity_matrix",
]


# ---------------------------------------------------------------------------
# Basic descriptive figures for candidate clusters
# ---------------------------------------------------------------------------


_ROLE_DYNAMIC_LABELS = {
    "putative_birth": "Putative birth",
    "relay_amplifier": "Relay amplifier",
    "merged_relay": "Merged relay",
    "terminal_sink": "Terminal sink",
    "isolated_burst": "Isolated burst",
    "unclear_origin": "Unclear origin",
    "no_observed_onward_spread": "No observed onward spread",
    "contained_burst": "Contained burst",
    "single_successor_chain": "Single-successor chain",
    "dominant_branch": "Dominant branch",
    "high_volume_onward_spread": "High-volume onward spread",
    "multi_branch_seeder": "Multi-branch seeder",
    "multi_branch_expander": "Multi-branch expander",
    "diverse_population_broadcaster": "Diverse population broadcaster",
    "weak_or_ambiguous_onward_spread": "Weak/ambiguous onward spread",
    "not_sse_like": "Not SSE-like",
    "mixed_population_dissemination": "Mixed-population dissemination",
    "putative_introduction_burst": "Putative introduction burst",
    "secondary_relay_amplification": "Secondary relay amplification",
    "diffuse_branching_transmission": "Diffuse branching transmission",
    "focused_branching_transmission": "Focused branching transmission",
    "sustained_single_chain": "Sustained single chain",
    "contained_local_burst": "Contained local burst",
    "high_volume_onward_transmission": "High-volume onward transmission",
    "ambiguous_amplification_signal": "Ambiguous amplification signal",
}


def _pretty_role_dynamic(value: Any, label_map: Mapping[str, str] | None = None) -> str:
    if pd.isna(value):
        return ""
    text = str(value)
    if label_map and text in label_map:
        return label_map[text]
    if text in _ROLE_DYNAMIC_LABELS:
        return _ROLE_DYNAMIC_LABELS[text]
    return text.replace("_", " ").strip().capitalize()


def _candidate_mask(
    df: pd.DataFrame,
    *,
    sse_col: str = "sse_candidate",
) -> pd.Series:
    """Return a boolean candidate mask from current or compatibility columns."""
    if sse_col in df.columns:
        return df[sse_col].fillna(False).astype(bool)
    if "candidate_tier" in df.columns:
        return df["candidate_tier"].isin(HIGH_PRIORITY_CANDIDATE_TIERS)
    raise KeyError(f"data needs '{sse_col}' or 'candidate_tier'")


def _with_candidate_status(
    df: pd.DataFrame,
    *,
    sse_col: str = "sse_candidate",
) -> pd.DataFrame:
    out = df.copy()
    mask = _candidate_mask(out, sse_col=sse_col)
    out["_sse_candidate_mask"] = mask
    out["_sse_status"] = np.where(mask, "Candidate", "Background")
    return out


def _signature_series(df: pd.DataFrame) -> pd.Series:
    if "sse_signature" in df.columns:
        signature = df["sse_signature"].astype("string").fillna("none")
    elif "axes_fired" in df.columns:
        signature = df["axes_fired"].astype("string").fillna("none")
    elif "sse_category" in df.columns:
        signature = df["sse_category"].astype("string").fillna("none")
    else:
        signature = pd.Series("none", index=df.index, dtype="string")
    return signature.mask(~_candidate_mask(df), "none")


def _signature_order(values: Iterable[Any]) -> list[str]:
    observed = [str(v) for v in pd.Series(list(values)).dropna().unique()]
    ordered = [v for v in SSE_SIGNATURE_ORDER if v in observed]
    ordered.extend(sorted(set(observed) - set(ordered) - {"none"}))
    return ordered


def _signature_color(value: Any) -> str:
    return SSE_SIGNATURE_PALETTE.get(str(value), CANDIDATE_COLOR)


def plot_meta_cluster_subgraph(
    node_stats: pd.DataFrame,
    edge_table: pd.DataFrame,
    meta_cluster_id: Any,
    *,
    meta_col: str = "meta_cluster_id",
    sse_col: str = "sse_candidate",
    max_nodes: int | None = 250,
    annotate_top_n: int = 0,
    rankdir: str = "LR",
    layout_method: str = "sugiyama",
    jitter: float | tuple[float, float] = 0.08,
    random_state: int | None = 42,
    edge_weight_col: str = "n_shared_sequences",
    scale_edges_by_weight: bool = True,
    edge_curve: float = 0.16,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 5.0,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Draw one connected-component subgraph with igraph and binary SSE status.

    Nodes are coloured only by candidate/background status. ``layout_method``
    can be ``"sugiyama"`` for igraph's layered DAG layout or ``"tree"`` /
    ``"dot"`` for Graphviz's tree-like hierarchical layout via pygraphviz.
    ``jitter`` applies a deterministic two-dimensional layout perturbation to
    reduce overlap. Pass a scalar to use the same relative spread on x and y, or
    ``(x_jitter, y_jitter)`` to tune horizontal and vertical spread separately.
    Edge weight defaults to ``n_shared_sequences``, the count of sequence IDs
    shared across adjacent-window clusters; it is overlap support, not observed
    transmission. When ``max_nodes`` is set and the component is larger,
    candidates are retained first and the remaining slots are filled by largest
    background nodes.
    """
    if "cluster_id" not in node_stats.columns:
        raise KeyError("node_stats needs 'cluster_id'")
    if meta_col not in node_stats.columns:
        if "connected_components" in node_stats.columns:
            meta_col = "connected_components"
        else:
            raise KeyError(f"node_stats needs '{meta_col}' or 'connected_components'")
    if not {"source", "target"}.issubset(edge_table.columns):
        raise KeyError("edge_table needs 'source' and 'target'")

    nodes = node_stats.loc[node_stats[meta_col].eq(meta_cluster_id)].copy()
    if nodes.empty:
        raise ValueError(f"No nodes found for {meta_col}={meta_cluster_id!r}")

    nodes = _with_candidate_status(nodes, sse_col=sse_col)
    if max_nodes is not None and len(nodes) > max_nodes:
        sort_cols = ["_sse_candidate_mask"]
        ascending = [False]
        if "cluster_size" in nodes.columns:
            sort_cols.append("cluster_size")
            ascending.append(False)
        if "window_idx" in nodes.columns:
            sort_cols.append("window_idx")
            ascending.append(True)
        nodes = nodes.sort_values(sort_cols, ascending=ascending).head(max_nodes)

    node_ids = nodes["cluster_id"].astype(str).tolist()
    node_id_set = set(node_ids)
    edges = edge_table.loc[
        edge_table["source"].astype(str).isin(node_id_set)
        & edge_table["target"].astype(str).isin(node_id_set)
    ].copy()
    edges["source"] = edges["source"].astype(str)
    edges["target"] = edges["target"].astype(str)

    import igraph as ig

    graph = ig.Graph(directed=True)
    graph.add_vertices(node_ids)
    edge_pairs = list(edges[["source", "target"]].itertuples(index=False, name=None))
    if edge_pairs:
        graph.add_edges(edge_pairs)
        if edge_weight_col in edges.columns:
            graph.es["weight"] = pd.to_numeric(
                edges[edge_weight_col],
                errors="coerce",
            ).fillna(1).clip(lower=1).to_numpy(dtype=float)
        else:
            graph.es["weight"] = np.ones(len(edge_pairs), dtype=float)

    node_lookup = nodes.assign(cluster_id=nodes["cluster_id"].astype(str)).set_index(
        "cluster_id"
    )
    graph.vs["candidate"] = [
        bool(node_lookup.loc[name, "_sse_candidate_mask"]) for name in graph.vs["name"]
    ]

    if "window_idx" in node_lookup.columns:
        windows = pd.Series( pd.to_numeric(
            node_lookup.loc[graph.vs["name"], "window_idx"], errors="coerce"
        ))
        layers = pd.factorize(windows.fillna(windows.min()).to_numpy(), sort=True)[0]
    else:
        layers = None

    layout_method = layout_method.lower()
    if layout_method not in {"sugiyama", "tree", "dot", "graphviz", "circle"}:
        raise ValueError(
            "layout_method must be one of 'sugiyama', 'tree', 'dot', "
            "'graphviz', or 'circle'."
        )

    if graph.vcount() == 1:
        coords = np.array([[0.0, 0.0]])
    elif layout_method in {"tree", "dot", "graphviz"}:
        import pygraphviz as pgv

        ag = pgv.AGraph(directed=True, strict=False, rankdir=rankdir.upper())
        ag.graph_attr.update(splines="true", overlap="false")
        for name in node_ids:
            ag.add_node(name)
        for source, target in edge_pairs:
            ag.add_edge(source, target)
        ag.layout(prog="dot")
        coords = np.asarray(
            [
                [
                    float(value)
                    for value in str(getattr(ag.get_node(name), "attr", {}).get("pos")).split(",")[:2]
                ]
                for name in graph.vs["name"]
            ],
            dtype=float,
        )
    elif layout_method == "circle":
        layout = graph.layout_circle()
        coords = np.asarray(layout.coords, dtype=float)
    elif graph.ecount() > 0:
        layout = graph.layout_sugiyama(layers=layers)
        coords = np.asarray(layout.coords, dtype=float)
    else:
        layout = graph.layout_circle()
        coords = np.asarray(layout.coords, dtype=float)

    horizontal = rankdir.upper() in {"LR", "RL"}
    if layout_method == "sugiyama" and horizontal:
        coords = coords[:, [1, 0]]
    coords = coords - coords.mean(axis=0, keepdims=True)

    if isinstance(jitter, tuple):
        if len(jitter) != 2:
            raise ValueError("jitter tuple must be (x_jitter, y_jitter).")
        jitter_xy = np.asarray(jitter, dtype=float)
    else:
        jitter_xy = np.asarray([float(jitter), float(jitter)], dtype=float)
    if np.any(jitter_xy < 0):
        raise ValueError("jitter values must be non-negative.")

    if np.any(jitter_xy > 0) and graph.vcount() > 1:
        rng = np.random.default_rng(random_state)
        span = np.ptp(coords, axis=0)
        fallback_span = max(float(np.nanmax(span)), 1.0)
        axis_scale = np.where(span > 0, span, fallback_span)
        coords += rng.normal(
            0.0,
            jitter_xy * axis_scale,
            size=(graph.vcount(), 2),
        )

    fig, ax = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
    )

    for edge in graph.es:
        source, target = edge.tuple
        x0, y0 = coords[source]
        x1, y1 = coords[target]
        weight = float(edge["weight"]) if "weight" in edge.attributes() else 1.0
        if scale_edges_by_weight and graph.ecount() > 0:
            weights = np.asarray(graph.es["weight"], dtype=float)
            max_weight = max(float(np.nanmax(weights)), 1.0)
            width_scale = np.sqrt(weight / max_weight)
            edge_lw = 0.45 + 1.8 * width_scale
            edge_alpha = 0.18 + 0.34 * width_scale
        else:
            edge_lw = 0.7
            edge_alpha = 0.32
        curve_sign = -1 if edge.index % 2 else 1
        connectionstyle = f"arc3,rad={curve_sign * edge_curve:.3f}"
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops={
                "arrowstyle": "-|>",
                "color": GRAY_LIGHT,
                "lw": edge_lw,
                "alpha": edge_alpha,
                "shrinkA": 4,
                "shrinkB": 4,
                "connectionstyle": connectionstyle,
            },
            zorder=1,
        )

    if "cluster_size" in node_lookup.columns:
        sizes = pd.Series(pd.to_numeric(
            node_lookup.loc[graph.vs["name"], "cluster_size"], errors="coerce"
        )).fillna(1)
        marker_sizes = 26 + 120 * np.sqrt(sizes / max(float(sizes.max()), 1.0))
    else:
        marker_sizes = np.full(graph.vcount(), 58.0)

    candidate = np.asarray(graph.vs["candidate"], dtype=bool)
    for mask, label, color, edge_color, zorder in [
        (~candidate, "Background", BACKGROUND_COLOR, BACKGROUND_DARK, 2),
        (candidate, "Candidate", CANDIDATE_COLOR, CANDIDATE_DARK, 3),
    ]:
        if not mask.any():
            continue
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=np.asarray(marker_sizes)[mask],
            c=color,
            edgecolors=edge_color,
            linewidths=0.45,
            alpha=0.88,
            label=f"{label} (n={int(mask.sum()):,})",
            zorder=zorder,
        )

    if annotate_top_n > 0 and "cluster_size" in node_lookup.columns:
        top_names = (
            node_lookup["cluster_size"].sort_values(ascending=False).head(annotate_top_n)
        ).index
        coord_lookup = dict(zip(graph.vs["name"], coords))
        for name in top_names:
            if name not in coord_lookup:
                continue
            x, y = coord_lookup[name]
            ax.text(x, y, str(name), fontsize="x-small", color=INK, zorder=4)

    ax.set_title(f"{meta_col}: {meta_cluster_id} ({graph.vcount():,} nodes)")
    ax.set_axis_off()
    legend_handles, legend_labels = ax.get_legend_handles_labels()
    if scale_edges_by_weight and graph.ecount() > 0:
        weights = np.asarray(graph.es["weight"], dtype=float)
        finite = weights[np.isfinite(weights)]
        if finite.size:
            for value in sorted(set([float(np.nanmin(finite)), float(np.nanmax(finite))])):
                max_weight = max(float(np.nanmax(finite)), 1.0)
                width_scale = np.sqrt(value / max_weight)
                legend_handles.append(
                    Line2D(
                        [0],
                        [0],
                        color=GRAY_LIGHT,
                        lw=0.45 + 1.8 * width_scale,
                        alpha=0.55,
                    )
                )
                legend_labels.append(f"{edge_weight_col}={value:g}")
    ax.legend(legend_handles, legend_labels, loc="upper left", frameon=False)
    plt.close(fig)
    return fig


def plot_cluster_size_distribution(
    df: pd.DataFrame,
    size_col: str = "cluster_size",
    *,
    by: str | None = "clade",
    sse_col: str = "sse_candidate",
    sse_labels: dict[int, str] | None = None,
    min_size: int = 1,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 6.5,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    complementary: bool = True,
) -> Figure:
    """Four-panel plot of cluster-size distributions, split by `sse_col`.

    Top row / bottom row correspond to the two values of the binary
    `sse_col`. Columns are ECDF (left) and violin (right), matching the
    two-panel layout. Size scales are shared across rows for comparability.
    """

    df = df.copy()
    transformed_min = np.log1p(min_size)
    df["_plot_size"] = np.log1p(df[size_col])
    df = df.loc[df["_plot_size"] >= transformed_min]

    # --- grouping (computed once, on full data, so rows stay consistent) ---
    if by is None or by not in df.columns:
        df["_plot_group"] = "All clusters"
        group_order = ["All clusters"]
        palette = {"All clusters": "#000000"}
    else:
        df["_plot_group"] = df[by].map(CLADES).fillna("Other")
        observed_groups = set(df["_plot_group"])
        group_order = [g for g in CLADES.values() if g in observed_groups]
        if "Other" in observed_groups:
            group_order.append("Other")
        palette = {g: CLADE_PALETTE.get(g, GRAY) for g in group_order}

    # --- row split on the binary column ---
    if sse_labels is None:
        sse_labels = {1: "Candidate", 0: "Background"}
    # top = 1, bottom = 0 (only rows present in the data)
    row_values = [v for v in (1, 0) if v in set(df[sse_col].dropna())]

    # --- shared size-axis ticks, from full filtered data ---
    smin = np.floor(df["_plot_size"].min())
    smax = np.ceil(df["_plot_size"].max())
    log_ticks = np.arange(smin, smax + 1)
    size_ticks = np.expm1(log_ticks)
    size_tick_labels = [f"{int(round(t))}" for t in size_ticks]

    clade_label_to_key = {label: key for key, label in CLADES.items()}

    fig, axes = new_figure(
        nrows=2,
        ncols=2,
        layout="constrained",
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
    )
    axes = axes

    legend_handles: list = []
    legend_labels: list = []

    def _draw_row(ax_ecdf, ax_violin, sub_df, row_label, *, is_bottom):
        nonlocal legend_handles, legend_labels

        if sub_df.empty:
            for ax in (ax_ecdf, ax_violin):
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes)
            return

        # Left: ECDF
        sns.ecdfplot(
            data=sub_df,
            x="_plot_size",
            hue="_plot_group",
            hue_order=group_order,
            palette=palette,
            complementary=complementary,
            stat="proportion",
            linewidth=1.5,
            ax=ax_ecdf,
        )
        ax_ecdf.set_ylabel("P(X ≥ Cluster size)")
        ax_ecdf.set_xticks(log_ticks)
        ax_ecdf.set_xlim(smin, smax)

        leg = ax_ecdf.get_legend()
        if leg:
            if not legend_handles:  # capture once, from the first populated row
                legend_handles = list(leg.legend_handles)
                legend_labels = [t.get_text() for t in leg.get_texts()]
            leg.remove()

        # Right: Violin
        sns.violinplot(
            data=sub_df,
            x="_plot_group",
            y="_plot_size",
            hue="_plot_group",
            order=group_order,
            hue_order=group_order,
            palette=palette,
            cut=0,
            inner="quartile",
            linewidth=0.8,
            ax=ax_violin,
        )
        ax_violin.set_ylabel("Cluster size")
        ax_violin.set_yticks(log_ticks)
        ax_violin.set_yticklabels(size_tick_labels)
        ax_violin.set_ylim(smin, smax)
        if ax_violin.get_legend():
            ax_violin.get_legend().remove()

        # X handling differs top vs bottom
        if is_bottom:
            ax_ecdf.set_xticklabels(size_tick_labels)
            ax_ecdf.set_xlabel("Cluster size")
            if by == "clade":
                x_tick_labels = [clade_label_to_key.get(g, g) for g in group_order]
                ax_violin.set_xticks(np.arange(len(group_order)))
                ax_violin.set_xticklabels(x_tick_labels, rotation=35)
                ax_violin.set_xlabel("SARS-CoV-2 clade")
            else:
                ax_violin.tick_params(axis="x")
                ax_violin.set_xlabel("")
        else:
            ax_ecdf.tick_params(labelbottom=False)
            ax_ecdf.set_xlabel("")
            ax_violin.tick_params(labelbottom=False)
            ax_violin.set_xlabel("")

        # Row label on the far left, rotated (won't collide with panel letters)
        ax_ecdf.annotate(
            row_label,
            xy=(0, 0.5), xycoords="axes fraction",
            xytext=(-ax_ecdf.yaxis.labelpad - 28, 0), textcoords="offset points",
            ha="right", va="center", rotation=90, fontweight="bold",
        )

    panel_axes = []
    for r, val in enumerate(row_values):
        ax_ecdf, ax_violin = axes[r, 0], axes[r, 1]
        sub = df.loc[df[sse_col] == val]
        _draw_row(
            ax_ecdf, ax_violin, sub,
            row_label=sse_labels.get(val, str(val)),
            is_bottom=(r == len(row_values) - 1),
        )
        panel_axes.extend([ax_ecdf, ax_violin])

    fig.legend(
        handles=legend_handles,
        labels=legend_labels,
        title="",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
    )
    add_panel_labels(panel_axes)
    plt.close(fig)

    return fig


def plot_role_dynamic_heatmap(
    candidates: pd.DataFrame,
    *,
    label_map: Mapping[str, str] | None = None,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 5.0,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Counts of candidate signatures by review tier or burden status."""
    d = candidates.copy()
    d = d.loc[_candidate_mask(d)].copy()
    if d.empty:
        raise ValueError("No candidate rows available for signature heatmap.")

    d["sse_signature"] = _signature_series(d)
    row_col = "candidate_tier" if "candidate_tier" in d.columns else "burden_status"
    if row_col not in d.columns:
        d[row_col] = "candidate"

    row_order = [v for v in d[row_col].dropna().astype(str).unique()]
    sig_order = _signature_order(d["sse_signature"])
    heat = pd.crosstab(d[row_col].astype(str), d["sse_signature"].astype(str)).reindex(
        index=row_order,
        columns=sig_order,
        fill_value=0,
    )
    heat_plot = np.log10(heat + 1)
    heat_plot.index = [_pretty_role_dynamic(v, label_map) for v in heat_plot.index]  # type: ignore
    heat_plot.columns = [_pretty_role_dynamic(v, label_map) for v in heat_plot.columns]  # type: ignore
    annot = heat.copy()
    annot.index = heat_plot.index  # type: ignore
    annot.columns = heat_plot.columns  # type: ignore

    fig, ax = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
    )
    sns.heatmap(
        heat_plot,
        annot=annot,
        fmt="d",
        cmap="YlGnBu",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Number of candidates [log10(n + 1)]"},
        ax=ax,
    )
    ax.set_ylabel(_pretty_role_dynamic(row_col))
    ax.set_xlabel("Detection signature")
    ax.tick_params(axis="x", rotation=35)
    ax.tick_params(axis="y", rotation=0)

    plt.close(fig)
    return fig


def plot_candidate_rate_over_time(
    node_stats: pd.DataFrame,
    *,
    sequence_df: pd.DataFrame | None = None,
    window_stride: int | None = 2,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 5.0,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Per-window candidate rate, with signature composition stacked underneath.

    Top panel: stacked % of unique sequences present in each retained window,
    where a sequence is candidate-associated if it appears in at least one
    candidate SSE node across the overlapping windows. Bottom panel: stacked
    raw counts of candidate nodes by ``sse_signature``.
    """
    required = {
        "cluster_id",
        "window_idx",
        "wn_mid_date",
    }
    missing = required.difference(node_stats.columns)
    if missing:
        raise KeyError(f"node_stats needs {sorted(missing)}")

    plot_df = _with_candidate_status(node_stats)
    plot_df["sse_signature"] = _signature_series(plot_df)
    candidate_mask = plot_df["_sse_candidate_mask"]
    candidate_nodes = set(plot_df.loc[candidate_mask, "cluster_id"].dropna())

    if sequence_df is None:
        sequence_df = load_analysis_columns(
            ["window_id", "window_idx", "wn_mid_date", "cluster_id", "sequence_id"],
            add_policy=False,
            window_stride=window_stride,
        )

    sequence_required = {"window_idx", "cluster_id", "sequence_id"}
    missing_sequence = sequence_required.difference(sequence_df.columns)
    if missing_sequence:
        raise KeyError(f"sequence_df needs {sorted(missing_sequence)}")

    seq = sequence_df.copy()
    seq = seq.loc[
        seq["window_idx"].isin(plot_df["window_idx"].dropna().unique())
    ].copy()
    if "wn_mid_date" not in seq.columns:
        window_dates = plot_df[["window_idx", "wn_mid_date"]].drop_duplicates(
            "window_idx"
        )
        seq = seq.merge(window_dates, on="window_idx", how="left")
    if "wn_mid_date" not in seq.columns:
        raise KeyError("sequence_df needs 'wn_mid_date' or node_stats must provide it")

    seq["wn_mid_date"] = pd.to_datetime(seq["wn_mid_date"], errors="coerce")
    seq["_in_candidate_node"] = seq["cluster_id"].isin(candidate_nodes)
    candidate_sequences = set(
        seq.loc[seq["_in_candidate_node"], "sequence_id"].dropna()
    )
    seq["_candidate_sequence"] = seq["sequence_id"].isin(candidate_sequences)

    seq_windows = (
        seq[["window_idx", "wn_mid_date", "sequence_id", "_candidate_sequence"]]
        .dropna(subset=["window_idx", "wn_mid_date", "sequence_id"])
        .drop_duplicates(["window_idx", "sequence_id"])
    )

    summary = seq_windows.groupby(["window_idx", "wn_mid_date"], as_index=False).agg(
        n_sequences=("sequence_id", "nunique"),
        n_candidate_sequences=("_candidate_sequence", "sum"),
    )
    summary = summary.sort_values("window_idx")
    summary["n_background_sequences"] = (
        summary["n_sequences"] - summary["n_candidate_sequences"]
    )
    summary["candidate_sequence_pct"] = np.where(
        summary["n_sequences"].gt(0),
        100 * summary["n_candidate_sequences"] / summary["n_sequences"],
        0.0,
    )
    summary["background_sequence_pct"] = np.where(
        summary["n_sequences"].gt(0),
        100 * summary["n_background_sequences"] / summary["n_sequences"],
        0.0,
    )
    window_candidate_rate = (
        plot_df.assign(_candidate_node=candidate_mask)
        .groupby(["window_idx", "wn_mid_date"], as_index=False)
        .agg(candidate_rate=("_candidate_node", lambda s: 100 * s.fillna(False).mean()))
        .sort_values("window_idx")
    )
    window_candidate_rate["wn_mid_date"] = pd.to_datetime(
        window_candidate_rate["wn_mid_date"],
        errors="coerce",
    )

    signature_counts = (
        plot_df.loc[candidate_mask]
        .groupby(["wn_mid_date", "sse_signature"], as_index=False)
        .size()
        .rename(columns={"size": "n"})
    )
    signature_pivot = (
        signature_counts.pivot(
            index="wn_mid_date",
            columns="sse_signature",
            values="n",
        )
        .fillna(0)
        .sort_index()
    )
    signature_pivot = signature_pivot.reindex(
        columns=_signature_order(signature_pivot.columns),
        fill_value=0,
    )

    fig, axes = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        nrows=2,
        sharex=True,
        context=context,
        font_scale=font_scale,
        layout="constrained",
    )
    ax = axes[0]
    ax.bar(
        summary["wn_mid_date"],
        summary["background_sequence_pct"],
        width=5,
        color=BACKGROUND_COLOR,
        edgecolor=BACKGROUND_DARK,
        alpha=0.86,
        label="background sequences",
    )
    ax.bar(
        summary["wn_mid_date"],
        summary["candidate_sequence_pct"],
        width=5,
        bottom=summary["background_sequence_pct"],
        color=CANDIDATE_COLOR,
        edgecolor=CANDIDATE_DARK,
        alpha=0.86,
        label="candidate-associated sequences",
    )
    ax.plot(
        window_candidate_rate["wn_mid_date"],
        window_candidate_rate["candidate_rate"],
        color="red",
        linestyle="--",
        marker="o",
        markersize=3,
        linewidth=1.2,
        label="Candidate rate by window",
    )
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.0f}%"))
    ax.set_ylabel("Percent")
    handles1, labels1 = ax.get_legend_handles_labels()
    ax.legend(
        handles1,
        labels1,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.2),
        ncol=3,
        frameon=False,
    )

    ax = axes[1]
    if len(signature_pivot.columns):
        colors = [_signature_color(signature) for signature in signature_pivot.columns]
        ax.stackplot(
            signature_pivot.index,
            [signature_pivot[c].to_numpy() for c in signature_pivot.columns],
            labels=[_pretty_role_dynamic(c) for c in signature_pivot.columns],
            colors=colors,
            alpha=0.86,
        )
        ax.legend(loc="upper left", ncol=2, frameon=False)
    ax.set_ylabel("Candidate nodes")
    ax.set_xlabel("Window midpoint date")
    fig.autofmt_xdate()

    add_panel_labels(list(axes))

    plt.close(fig)
    return fig


def plot_core_metric_space(
    node_stats: pd.DataFrame,
    *,
    x_col: str | None = None,
    y_col: str | None = None,
    height_in: float = 4,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    min_size: int = 1,
) -> Figure:
    """Scatter of local-burst vs onward-burden metrics by SSE status."""
    set_theme(
        context=context,
        font_scale=font_scale,
    )

    plot_df = _with_candidate_status(
        node_stats.loc[node_stats["cluster_size"].ge(min_size)].copy()
    )
    x_candidates = [
        "burst_score_null_z",
        "burst_score",
        "sampling_adjusted_excess_size",
        "log_cluster_size",
    ]
    y_candidates = [
        "burden_score_null_z",
        "burden_score",
        "log_new_downstream_burden_ratio",
        "log_supported_new_downstream_burden_ratio",
    ]
    x_col = x_col or next((col for col in x_candidates if col in plot_df.columns), None)
    y_col = y_col or next(
        (
            col
            for col in y_candidates
            if col in plot_df.columns and plot_df[col].notna().sum() >= 5
        ),
        None,
    )
    if x_col is None or y_col is None:
        raise ValueError("No usable burst/burden metric columns are present.")
    plot_df = plot_df.dropna(subset=[x_col, y_col]).copy()

    clip_hi = plot_df["cluster_size"].quantile(0.995)
    plot_df["marker_size"] = np.sqrt(
        plot_df["cluster_size"].clip(lower=1, upper=clip_hi)
    )

    fig, ax = new_figure(
        width="double",
        height_in=height_in,
        context=context,
        font_scale=font_scale,
    )
    for status, color, edgecolor, alpha, zorder in [
        ("Background", BACKGROUND_COLOR, BACKGROUND_DARK, 0.28, 1),
        ("Candidate", CANDIDATE_COLOR, CANDIDATE_DARK, 0.82, 2),
    ]:
        sub = plot_df.loc[plot_df["_sse_status"].eq(status)]
        if sub.empty:
            continue
        size_min = plot_df["marker_size"].min()
        size_max = plot_df["marker_size"].max()
        if size_min == size_max:
            point_sizes = np.full(len(sub), 60.0)
        else:
            point_sizes = np.interp(
                sub["marker_size"],
                (size_min, size_max),
                (12, 180),
            )
        ax.scatter(
            sub[x_col],
            sub[y_col],
            s=point_sizes,
            c=color,
            edgecolors=edgecolor,
            linewidths=0.25,
            alpha=alpha,
            label=f"{status} (n={len(sub):,})",
            zorder=zorder,
        )

    ax.set_xlabel(_pretty_role_dynamic(x_col))
    ax.set_ylabel(_pretty_role_dynamic(y_col))
    ax.legend(frameon=False, loc="best")
    plt.close(fig)
    return fig


def plot_composite_distributions(
    node_stats: pd.DataFrame,
    *,
    columns: Iterable[tuple[str, str]] = (
        ("log_cluster_size_pct_window", "Cluster size"),
        ("sampling_adjusted_excess_size_pct_window", "Context excess size"),
        ("log_excess_over_upstream_pct_window", "Excess over upstream"),
        ("log_new_downstream_burden_ratio_pct_window", "New onward burden ratio"),
        (
            "log_supported_new_downstream_burden_ratio_pct_window",
            "Supported onward burden ratio",
        ),
        ("burst_score", "Burst score"),
        ("burden_score", "Burden score"),
    ),
    nrows: int = 2,
    ncols: int = 3,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 4.5,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    min_size: int = 1,
) -> Figure:
    """Overlaid KDE of score components for candidates vs background.

    Defaults use the current two-axis detector components: local burst,
    upstream excess, and onward-burden ratios/scores.
    """
    node_stats = _with_candidate_status(node_stats)

    columns = [s for s in columns if s[0] in node_stats.columns]
    if not columns:
        raise ValueError("None of the requested score columns are present.")

    fig, axes = new_figure(
        nrows=nrows,
        ncols=ncols,
        layout="constrained",
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
        sharex=True,
    )

    background_label = (
        f"Background (n≥{min_size} sequences, "
        f"{len(node_stats[~node_stats['_sse_candidate_mask'] & node_stats['cluster_size'].ge(min_size)]):,} nodes)"
    )

    candidate_label = (
        f"Candidate (n≥{min_size} sequences, "
        f"{len(node_stats[node_stats['_sse_candidate_mask'] & node_stats['cluster_size'].ge(min_size)]):,} nodes)"
    )

    axes = axes.flatten()

    for ax, col in zip(axes, columns):
        for label, color, mask in [
            (
                "background",
                BACKGROUND_COLOR,
                (
                    ~node_stats["sse_candidate"]
                    if "sse_candidate" in node_stats.columns
                    else ~node_stats["_sse_candidate_mask"]
                ),
            ),
            ("candidate", CANDIDATE_COLOR, node_stats["_sse_candidate_mask"]),
        ]:
            if isinstance(mask, pd.Series):
                mask = mask & node_stats["cluster_size"].ge(min_size)
            values = node_stats[mask][col[0]].dropna().to_numpy()
            if len(values) < 5:
                continue
            sns.kdeplot(
                values,
                ax=ax,
                fill=True,
                color=color,
                alpha=0.35,
                linewidth=1.2,
                label=label,
                common_norm=False,
            )
        ax.set_title(col[1])
        ax.set_ylabel("Density")
        if col[0].endswith("_pct_window") or col[0].endswith("_pct_onward_window"):
            ax.set_xlim(0, 1)
    for ax in axes[len(columns) :]:
        ax.set_visible(False)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=BACKGROUND_COLOR,
            linewidth=1.8,
            label=background_label,
        ),
        Line2D(
            [0],
            [0],
            color=CANDIDATE_COLOR,
            linewidth=1.8,
            label=candidate_label,
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="outside upper center",
        bbox_to_anchor=(0.5, 1.1),
        ncol=2,
        frameon=False,
        handlelength=1.8,
        columnspacing=1.6,
    )
    fig.supxlabel("Within window percentile score (higher = more extreme)")
    add_panel_labels([ax for ax in axes[: len(columns)] if ax.get_visible()])
    plt.close(fig)
    return fig


_DEFAULT_INDIVIDUAL_CATEGORICAL_VARIABLES: list[tuple[str, str]] = [
    ("sex", "Sex"),
    ("age_band", "Age band"),
    ("dz_simd_quintile", "SIMD quintile"),
    ("dz_urban_rural_class", "Urban/rural class"),
    ("dz_health_board", "Health board"),
]

_SEX_ORDER = ["Male", "Female", "Other", "Unknown", "Missing"]
_SIMD_ORDER = ["1", "2", "3", "4", "5", "Unknown", "Missing"]
_URBAN_RURAL_ORDER = [
    "Large Urban Areas",
    "Other Urban Areas",
    "Accessible Small Towns",
    "Remote Small Towns",
    "Very Remote Small Towns",
    "Accessible Rural Areas",
    "Remote Rural Areas",
    "Very Remote Rural Areas",
    "Unknown",
    "Missing",
]


def _normalise_categorical_specs(
    variables: Sequence[tuple[str, str] | Mapping[str, Any]] | None,
) -> list[tuple[str, str]]:
    if variables is None:
        return list(_DEFAULT_INDIVIDUAL_CATEGORICAL_VARIABLES)

    specs: list[tuple[str, str]] = []
    for item in variables:
        if isinstance(item, Mapping):
            column = str(item["column"])
            label = str(item.get("label", _pretty_text(column)))
        else:
            column, label = item
            column = str(column)
            label = str(label)
        specs.append((column, label))
    return specs


def _age_band_sort_key(value: Any) -> tuple[int, int, str]:
    text = str(value).strip()
    if text.lower() in {"missing", "unknown", "nan", "none"}:
        return (10_000, 0, text)
    match = re.match(r"^<?\s*(\d+)\s*(?:[-–]\s*(\d+)|\+)?$", text)
    if match:
        lower = int(match.group(1))
        upper = int(match.group(2)) if match.group(2) else lower
        return (lower, upper, text)
    return (9_000, 0, text)


def _ordered_category_levels(
    summary: pd.DataFrame,
    *,
    column: str,
    max_levels: int | None,
) -> tuple[pd.DataFrame, list[str]]:
    """Return a plotting summary with stable level order for common fields."""
    summary = summary.copy()
    column_key = column.lower()
    observed = summary.groupby("_level", as_index=False)["n"].sum()
    observed["_level"] = observed["_level"].astype(str)
    observed_set = set(observed["_level"])

    if column_key in {"sex"}:
        order = [level for level in _SEX_ORDER if level in observed_set]
        order.extend(
            sorted(observed_set - set(order), key=lambda value: str(value).casefold())
        )
    elif "simd" in column_key and "quintile" in column_key:
        order = [level for level in _SIMD_ORDER if level in observed_set]
        order.extend(
            sorted(observed_set - set(order), key=lambda value: str(value).casefold())
        )
    elif "age" in column_key:
        order = sorted(observed["_level"], key=_age_band_sort_key)
    elif "urban" in column_key or "rural" in column_key:
        order = [level for level in _URBAN_RURAL_ORDER if level in observed_set]
        order.extend(
            sorted(observed_set - set(order), key=lambda value: str(value).casefold())
        )
    else:
        order = (
            observed.sort_values(["n", "_level"], ascending=[False, True])["_level"]
            .astype(str)
            .tolist()
        )

    if max_levels is not None and max_levels > 1 and len(order) > max_levels:
        keep = set(order[: max_levels - 1])
        other_label = f"Other ({len(order) - len(keep)} levels)"
        summary["_level"] = summary["_level"].where(
            summary["_level"].isin(keep), other_label
        )
        summary = summary.groupby(["_level", "_sse_status"], as_index=False)["n"].sum()
        order = order[: max_levels - 1] + [other_label]

    return summary, order


def _individual_category_summary(
    data: pd.DataFrame,
    *,
    column: str,
    status_col: str,
    sequence_col: str,
    unit: str,
) -> pd.DataFrame:
    if column not in data.columns:
        raise KeyError(f"{column!r} is not present in the input data.")

    df = _with_candidate_status(data, sse_col=status_col)
    keep_cols = [column, "_sse_status"]
    if unit == "sequences" and sequence_col in df.columns:
        keep_cols.append(sequence_col)
    d = df.loc[:, keep_cols].copy()
    d["_level"] = (
        d[column]
        .astype("string")
        .fillna("Missing")
        .str.strip()
        .replace({"": "Missing", "<NA>": "Missing", "nan": "Missing"})
    )

    if unit == "sequences" and sequence_col in d.columns:
        d = d.dropna(subset=[sequence_col])
        d = d.drop_duplicates([sequence_col, "_sse_status", "_level"])
    elif unit != "rows":
        raise ValueError("unit must be 'sequences' or 'rows'.")

    return d.groupby(["_level", "_sse_status"], as_index=False).size().rename(
        columns={"size": "n"}
    )


def plot_individual_categorical_distribution_bars(
    data: pd.DataFrame,
    *,
    variables: Sequence[tuple[str, str] | Mapping[str, Any]] | None = None,
    status_col: str = "candidate",
    sequence_col: str = "sequence_id",
    unit: str = "sequences",
    max_levels: int | None = None,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float | None = None,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Plot within-status categorical distributions for individual-level rows.

    The default input is ``AssociationFrames.composition_base`` from
    ``load_association_frames(run_composition=True)``. Percentages are computed
    separately within candidate-associated and background-associated records.
    """

    specs = _normalise_categorical_specs(variables)
    if not specs:
        raise ValueError("At least one categorical variable is required.")
    if unit not in {"sequences", "rows"}:
        raise ValueError("unit must be 'sequences' or 'rows'.")

    ncols = 1 if len(specs) == 1 else 2
    nrows = int(np.ceil(len(specs) / ncols))
    if height_in is None:
        height_in = max(3.0, 2.4 * nrows)

    fig, axes = new_figure(
        nrows=nrows,
        ncols=ncols,
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
        layout="constrained",
        squeeze=False,
    )
    axes_flat = np.asarray(axes).reshape(-1)

    bar_height = 0.34
    statuses = [
        ("Background", BACKGROUND_COLOR, BACKGROUND_DARK, -bar_height / 2),
        ("Candidate", CANDIDATE_COLOR, CANDIDATE_DARK, bar_height / 2),
    ]

    for ax, (column, label) in zip(axes_flat, specs):
        summary = _individual_category_summary(
            data,
            column=column,
            status_col=status_col,
            sequence_col=sequence_col,
            unit=unit,
        )
        summary, order = _ordered_category_levels(
            summary,
            column=column,
            max_levels=max_levels,
        )

        totals = summary.groupby("_sse_status")["n"].transform("sum")
        summary["pct"] = np.where(totals.gt(0), 100 * summary["n"] / totals, 0.0)
        pivot = (
            summary.pivot_table(
                index="_level",
                columns="_sse_status",
                values="pct",
                aggfunc="sum",
                fill_value=0,
            )
            .reindex(order, fill_value=0)
            .reindex(columns=["Background", "Candidate"], fill_value=0)
        )

        y = np.arange(len(order))
        for status, color, edgecolor, offset in statuses:
            ax.barh(
                y + offset,
                pivot[status].to_numpy(),
                height=bar_height,
                color=color,
                edgecolor=edgecolor,
                linewidth=0.4,
                alpha=0.9,
                label=status,
            )

        ax.set_yticks(y)
        ax.set_yticklabels(order)
        ax.invert_yaxis()
        ax.set_title(label)
        ax.set_xlabel("Within-status share (%)")
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.0f}%"))
        ax.grid(axis="x", color=GRID, linewidth=0.8)
        ax.grid(axis="y", visible=False)
        xmax = float(pivot.to_numpy().max()) if not pivot.empty else 0.0
        ax.set_xlim(0, max(5.0, np.ceil(xmax / 5) * 5))

    for ax in axes_flat[len(specs) :]:
        ax.set_visible(False)

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=color, edgecolor=edgecolor, label=status)
        for status, color, edgecolor, _ in statuses
    ]
    fig.legend(
        handles=handles,
        loc="outside upper center",
        bbox_to_anchor=(0.5, 1.04),
        ncol=2,
        frameon=False,
        columnspacing=1.5,
    )

    add_panel_labels([ax for ax in axes_flat[: len(specs)] if ax.get_visible()])
    plt.close(fig)
    return fig


def plot_socio_demo_breakdown(
    node_stats: pd.DataFrame,
    col: str = "top_simd_quintiles",
    score: str = "simd_entropy_obs",
    labels: tuple[str, str] = (
        "Socioeconomic mixing score",
        "SIMD Quintile (1 = most deprived)",
    ),
    *,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 4.0,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    min_size: int = 1,
) -> Figure:
    """
    Plot the distribution of mixing scores for candidate vs background nodes,
    and the distribution of class-label frequencies for candidate vs background nodes.
    """

    node_stats = _with_candidate_status(
        node_stats.loc[node_stats["cluster_size"].ge(min_size)].copy()
    )

    fig, axes = new_figure(
        ncols=2,
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
        layout="constrained",
    )

    ax = axes[0]

    for label, color, mask in [
        ("background", BACKGROUND_COLOR, ~node_stats["_sse_candidate_mask"]),
        ("candidate", CANDIDATE_COLOR, node_stats["_sse_candidate_mask"]),
    ]:
        values = node_stats[mask][score].dropna().to_numpy()

        if len(values) < 5:
            continue

        sns.kdeplot(
            values,
            ax=ax,
            fill=True,
            color=color,
            alpha=0.35,
            linewidth=1.2,
            label=label,
            common_norm=False,
        )

    ax.set_xlabel(labels[0])
    ax.set_ylabel("Density")
    ax.legend(loc="best", frameon=False)

    ax = axes[1]
    if col not in node_stats.columns:
        ax.set_visible(False)
        add_panel_labels([axes[0]])
        plt.close(fig)
        return fig

    def _parse_freq_counts(s: str):
        """
        Parse strings like:
        'class label one (3); class label two (10); class label three (1)'

        Returns:
        [('class label one', 3), ('class label two', 10), ...]
        """
        if not isinstance(s, str) or not s.strip():
            return []

        out = []

        for part in s.split(";"):
            part = part.strip()

            if not part:
                continue

            match = re.match(r"^(.*?)\s*\((\d+)\)\s*$", part)

            if match is None:
                continue

            label = match.group(1).strip()
            count = int(match.group(2))

            if label and count > 0:
                out.append((label, count))

        return out

    records = []

    for _, row in node_stats.iterrows():
        candidate = row["_sse_candidate_mask"]

        for q, n in _parse_freq_counts(row[col]):
            records.append(
                {
                    "q": q,
                    "candidate": candidate,
                    "n": n,
                }
            )

    if records:
        share = (
            pd.DataFrame(records).groupby(["q", "candidate"], as_index=False)["n"].sum()
        )

        denom = share.groupby("candidate")["n"].transform("sum")
        share["frac"] = share["n"] / denom

        order = sorted(share["q"].unique())

        bar_h = 0.4
        y = np.arange(len(order))

        for off, candidate_value, lbl, color in [
            (-bar_h / 2, False, "background", BACKGROUND_COLOR),
            (bar_h / 2, True, "candidate", CANDIDATE_COLOR),
        ]:
            sub = (
                share.loc[share["candidate"] == candidate_value]
                .set_index("q")
                .reindex(order, fill_value=0)
            )

            ax.barh(
                y + off,
                sub["frac"].to_numpy(),
                height=bar_h,
                color=color,
                label=lbl,
            )

        ax.set_yticks(y)
        ax.set_yticklabels(order)

    ax.set_ylabel(labels[1])
    ax.set_xlabel("Fraction of sequence counts")

    add_panel_labels(axes)

    plt.close(fig)
    return fig


# ---------------------------------------------------------------------------
# Regression output figures and manuscript tables
# ---------------------------------------------------------------------------


# Shared regression labels/table helpers -------------------------------------


_COMPOSITION_ORDER = [
    "sex",
    "age_band",
    "simd_quintile",
    "urban_rural_class",
    "health_board",
]

_MIXING_ORDER = [
    "sex_entropy_z",
    "age_entropy_z",
    "simd_entropy_z",
    "urban_rural_entropy_z",
    "health_board_entropy_z",
    "all_mixing",
    "sex_entropy_obs",
    "age_entropy_obs",
    "simd_entropy_obs",
    "urban_rural_entropy_obs",
    "health_board_entropy_obs",
]

_PRETTY_LABELS = {
    "sex": "Sex",
    "age_band": "Age band",
    "simd_quintile": "SIMD quintile",
    "urban_rural_class": "Urban/rural class",
    "health_board": "Health board",
    "dz_simd_quintile": "SIMD quintile",
    "dz_urban_rural_class": "Urban/rural class",
    "dz_health_board": "Health board",
    "sex_entropy_obs": "Sex entropy",
    "age_entropy_obs": "Age entropy",
    "simd_entropy_obs": "SIMD entropy",
    "urban_rural_entropy_obs": "Urban/rural entropy",
    "health_board_entropy_obs": "Health-board entropy",
    "sex_entropy_z": "Sex entropy z-score",
    "age_entropy_z": "Age entropy z-score",
    "simd_entropy_z": "SIMD entropy z-score",
    "urban_rural_entropy_z": "Urban/rural entropy z-score",
    "health_board_entropy_z": "Health-board entropy z-score",
    "all_mixing": "All mixing predictors",
    "all_mixing_predictors": "All mixing predictors",
    "primary": "Primary",
    "expanded": "Expanded",
    "single": "Single predictor",
    "joint": "Joint",
    "composition": "Composition",
    "node_mixing": "Node mixing",
}


def _pretty_text(value: Any, label_map: Mapping[str, str] | None = None) -> str:
    """Human-readable label with conservative project-specific replacements."""
    if pd.isna(value):
        return ""
    text = str(value)
    if label_map and text in label_map:
        return label_map[text]
    if text in _PRETTY_LABELS:
        return _PRETTY_LABELS[text]
    return text.replace("_", " ").strip().capitalize()


def _term_level(term: Any) -> str:
    """Extract the displayed contrast level from a simple patsy term."""
    if pd.isna(term):
        return ""
    term = str(term)
    match = re.search(r"\[T\.(.*)\]$", term)
    if match:
        return match.group(1)
    return _pretty_text(term)


def _filter_regression_table(
    table: pd.DataFrame,
    *,
    domain: str | Iterable[str] | None = None,
    model_set: str | Iterable[str] | None = None,
    predictor_set: str | Iterable[str] | None = None,
    predictors: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Filter a regression output table without mutating the caller's data."""
    df = table.copy()
    df.columns = [str(col).strip() for col in df.columns]
    string_cols = df.select_dtypes(include=["object", "string"]).columns
    for col in string_cols:
        df[col] = df[col].map(lambda x: x.strip() if isinstance(x, str) else x)

    for col, value in [
        ("domain", domain),
        ("model_set", model_set),
        ("predictor_set", predictor_set),
    ]:
        if value is None or col not in df.columns:
            continue
        values = [value] if isinstance(value, str) else list(value)
        df = df.loc[df[col].isin(values)]

    if predictors is not None and "predictor" in df.columns:
        df = df.loc[df["predictor"].isin(predictors)]

    return df.copy()


def _with_regression_display_labels(
    df: pd.DataFrame,
    *,
    label_map: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Add stable display labels used by both figures and tables."""
    out = df.copy()

    if "label" in out.columns:
        out["_predictor_label"] = out["label"].map(lambda x: _pretty_text(x, label_map))
    elif "predictor" in out.columns:
        out["_predictor_label"] = out["predictor"].map(
            lambda x: _pretty_text(x, label_map)
        )
    else:
        out["_predictor_label"] = ""

    if "predictor" in out.columns:
        predictor_fallback = out["predictor"].map(lambda x: _pretty_text(x, label_map))
        all_mixing = out["_predictor_label"].str.lower().eq("all mixing predictors")
        out.loc[all_mixing, "_predictor_label"] = predictor_fallback.loc[all_mixing]
        if "term" in out.columns:
            term_fallback = out["term"].map(lambda x: _pretty_text(x, label_map))
            out.loc[all_mixing, "_predictor_label"] = term_fallback.loc[all_mixing]

    if "term" in out.columns:
        term_label = out["term"].map(_term_level)
    else:
        term_label = pd.Series("", index=out.index)
    out["_term_label"] = term_label

    if "reference" in out.columns:
        out["_contrast_label"] = np.where(
            out["_term_label"].eq(out["_predictor_label"]) | out["_term_label"].eq(""),
            out["_predictor_label"],
            out["_term_label"] + " vs " + out["reference"].astype(str),
        )
    else:
        out["_contrast_label"] = np.where(
            out["_term_label"].eq(out["_predictor_label"]) | out["_term_label"].eq(""),
            out["_predictor_label"],
            out["_term_label"],
        )

    if {"model_set", "predictor_set"}.issubset(out.columns):
        out["_model_label"] = (
            out["model_set"].map(lambda x: _pretty_text(x, label_map))
            + "\n"
            + out["predictor_set"].map(lambda x: _pretty_text(x, label_map))
        )
    elif "model_set" in out.columns:
        out["_model_label"] = out["model_set"].map(lambda x: _pretty_text(x, label_map))
    elif "predictor_set" in out.columns:
        out["_model_label"] = out["predictor_set"].map(
            lambda x: _pretty_text(x, label_map)
        )
    else:
        out["_model_label"] = ""

    return out


def _predictor_sort_key(values: pd.Series) -> pd.Series:
    order = {name: i for i, name in enumerate(_COMPOSITION_ORDER + _MIXING_ORDER)}
    return values.map(lambda x: order.get(str(x), len(order)))


def _regression_sort_source(df: pd.DataFrame) -> pd.Series:
    """Predictor-like series for stable display ordering."""
    if "predictor" not in df.columns:
        return pd.Series("", index=df.index)
    source = df["predictor"].astype(str).copy()
    if "term" in df.columns:
        use_term = source.eq("all_mixing") | source.eq("all_mixing_predictors")
        source.loc[use_term] = df.loc[use_term, "term"].astype(str)
    return source


def _format_p_value(
    value: Any,
    *,
    threshold: float = 0.001,
    digits: int = 3,
) -> str:
    if pd.isna(value):
        return ""
    value = float(value)
    if value < threshold:
        return f"<{threshold:.{digits}f}"
    return f"{value:.{digits}f}"


def _format_number(value: Any, *, digits: int = 2) -> str:
    if pd.isna(value):
        return ""
    return f"{float(value):.{digits}f}"


def _format_int(value: Any) -> str:
    if pd.isna(value):
        return ""
    return f"{int(round(float(value))):,}"


def _format_or_ci(row: pd.Series, *, digits: int = 2) -> str:
    if row[["odds_ratio", "or_low", "or_high"]].isna().any():
        return ""
    return (
        f"{row['odds_ratio']:.{digits}f} "
        f"({row['or_low']:.{digits}f}-{row['or_high']:.{digits}f})"
    )


# Manuscript table builders --------------------------------------------------


def make_regression_wald_table(
    wald_df: pd.DataFrame,
    *,
    domain: str | Iterable[str] | None = None,
    model_set: str | Iterable[str] | None = None,
    predictor_set: str | Iterable[str] | None = None,
    predictors: Iterable[str] | None = None,
    p_col: str = "P>chi2",
    p_adj_col: str = "p_adj_bh",
    label_map: Mapping[str, str] | None = None,
    digits: int = 2,
) -> pd.DataFrame:
    """Return a manuscript-facing omnibus Wald table.

    The table keeps the publication-facing model identifiers, test statistic,
    adjusted p-value, and interpretable sample sizes.
    """
    df = _filter_regression_table(
        wald_df,
        domain=domain,
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
    )
    if df.empty:
        return pd.DataFrame()

    df = _with_regression_display_labels(df, label_map=label_map)
    df["_predictor_order"] = _predictor_sort_key(_regression_sort_source(df))
    df = df.sort_values(
        [c for c in ["domain", "model_set", "predictor_set"] if c in df.columns]
        + (["_predictor_order"] if "predictor" in df.columns else []),
    )

    out = pd.DataFrame(index=df.index)
    if "model_set" in df.columns:
        out["Model set"] = df["model_set"].map(lambda x: _pretty_text(x, label_map))
    if "predictor_set" in df.columns:
        out["Specification"] = df["predictor_set"].map(
            lambda x: _pretty_text(x, label_map)
        )
    out["Predictor"] = df["_predictor_label"]
    if "reference" in df.columns:
        out["Reference/scale"] = df["reference"].astype(str)
    if "df" in df.columns:
        out["df"] = df["df"].map(_format_int)
    if "chi2" in df.columns:
        out["Wald chi-square"] = df["chi2"].map(
            lambda x: _format_number(x, digits=digits)
        )
    if p_adj_col in df.columns:
        out["FDR-adjusted P value"] = df[p_adj_col].map(_format_p_value)
    elif p_col in df.columns:
        out["P value"] = df[p_col].map(_format_p_value)
    for source, target in [
        ("n_sequences", "Sequences"),
        ("n_nodes", "Nodes"),
    ]:
        if source in df.columns:
            out[target] = df[source].map(_format_int)

    return out.reset_index(drop=True)


def make_regression_odds_ratio_table(
    odds_df: pd.DataFrame,
    *,
    domain: str | Iterable[str] | None = None,
    model_set: str | Iterable[str] | None = None,
    predictor_set: str | Iterable[str] | None = None,
    predictors: Iterable[str] | None = None,
    terms: Iterable[str] | None = None,
    p_col: str = "p_value",
    label_map: Mapping[str, str] | None = None,
    digits: int = 2,
) -> pd.DataFrame:
    """Return a manuscript-facing coefficient odds-ratio table."""
    df = _filter_regression_table(
        odds_df,
        domain=domain,
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
    )
    if terms is not None and "term" in df.columns:
        df = df.loc[df["term"].isin(terms)].copy()
    if df.empty:
        return pd.DataFrame()

    df = _with_regression_display_labels(df, label_map=label_map)
    if "predictor" in df.columns:
        df["_predictor_order"] = _predictor_sort_key(_regression_sort_source(df))
        df = df.sort_values(
            [
                c
                for c in [
                    "domain",
                    "model_set",
                    "predictor_set",
                    "_predictor_order",
                    "term",
                ]
                if c in df.columns
            ]
        )

    out = pd.DataFrame(index=df.index)
    if "model_set" in df.columns:
        out["Model set"] = df["model_set"].map(lambda x: _pretty_text(x, label_map))
    if "predictor_set" in df.columns:
        out["Specification"] = df["predictor_set"].map(
            lambda x: _pretty_text(x, label_map)
        )
    out["Predictor"] = df["_predictor_label"]
    out["Contrast"] = df["_contrast_label"]
    if "reference" in df.columns:
        out["Reference/scale"] = df["reference"].astype(str)
    out["Odds ratio (95% CI)"] = df.apply(
        lambda row: _format_or_ci(row, digits=digits), axis=1
    )
    if p_col in df.columns:
        out["P value"] = df[p_col].map(_format_p_value)
    for source, target in [
        ("n_sequences", "Sequences"),
        ("n_nodes", "Nodes"),
    ]:
        if source in df.columns:
            out[target] = df[source].map(_format_int)

    return out.reset_index(drop=True)


def make_regression_fit_table(
    fit_df: pd.DataFrame,
    *,
    domain: str | Iterable[str] | None = None,
    model_set: str | Iterable[str] | None = None,
    predictor_set: str | Iterable[str] | None = None,
    predictors: Iterable[str] | None = None,
    label_map: Mapping[str, str] | None = None,
    digits: int = 3,
) -> pd.DataFrame:
    """Return a manuscript-facing model fit comparison table."""
    df = _filter_regression_table(
        fit_df,
        domain=domain,
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
    )
    if df.empty:
        return pd.DataFrame()

    df = _with_regression_display_labels(df, label_map=label_map)
    if "predictor" in df.columns:
        df["_predictor_order"] = _predictor_sort_key(_regression_sort_source(df))
        df = df.sort_values(
            [
                c
                for c in ["domain", "model_set", "predictor_set", "_predictor_order"]
                if c in df.columns
            ]
        )

    out = pd.DataFrame(index=df.index)
    if "model_set" in df.columns:
        out["Model set"] = df["model_set"].map(lambda x: _pretty_text(x, label_map))
    if "predictor_set" in df.columns:
        out["Specification"] = df["predictor_set"].map(
            lambda x: _pretty_text(x, label_map)
        )
    out["Predictor"] = df["_predictor_label"]
    if "r2_mcfadden" in df.columns:
        out["McFadden pseudo-R2"] = df["r2_mcfadden"].map(
            lambda x: _format_number(x, digits=digits)
        )
    for source, target in [
        ("aic", "AIC"),
        ("bic_llf", "BIC"),
        ("log_likelihood", "Log likelihood"),
    ]:
        if source in df.columns:
            out[target] = df[source].map(lambda x: _format_number(x, digits=1))
    for source, target in [
        ("n_sequences", "Sequences"),
        ("n_nodes", "Nodes"),
    ]:
        if source in df.columns:
            out[target] = df[source].map(_format_int)
    if "converged" in df.columns:
        out["Converged"] = df["converged"].map(
            lambda x: "" if pd.isna(x) else str(bool(x))
        )

    return out.reset_index(drop=True)


# Wald heatmap ---------------------------------------------------------------


def _prepare_regression_wald_heatmap(
    wald_df: pd.DataFrame,
    *,
    domain: str | Iterable[str] | None = None,
    model_set: str | Iterable[str] | None = None,
    predictor_set: str | Iterable[str] | None = None,
    predictors: Iterable[str] | None = None,
    p_col: str = "p_adj_bh",
    row_col: str | None = None,
    label_map: Mapping[str, str] | None = None,
    cap_neg_log10_p: float = 20.0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """Prepare pivoted Wald evidence and p-value annotation tables."""
    df = _filter_regression_table(
        wald_df,
        domain=domain,
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
    )
    if df.empty:
        raise ValueError("No Wald rows remain after filtering.")

    if p_col not in df.columns:
        fallback = "P>chi2"
        if fallback not in df.columns:
            raise ValueError(f"`wald_df` must contain `{p_col}` or `{fallback}`.")
        p_col = fallback

    df = _with_regression_display_labels(df, label_map=label_map)
    if row_col is None:
        df["_row_id"] = _regression_sort_source(df)
    elif row_col not in df.columns:
        raise ValueError(f"`wald_df` does not contain row column `{row_col}`.")
    else:
        df["_row_id"] = df[row_col].astype(str)

    row_labels = (
        df.sort_values("_row_id", key=_predictor_sort_key)
        .drop_duplicates("_row_id")
        .set_index("_row_id")["_predictor_label"]
    )
    missing_label = row_labels.str.lower().eq("all mixing predictors")
    if missing_label.any() and "term" in df.columns:
        replacements = (
            df.drop_duplicates("_row_id").set_index("_row_id")["term"].map(_pretty_text)
        )
        row_labels.loc[missing_label] = replacements.loc[missing_label]

    df["_p_for_plot"] = pd.to_numeric(df[p_col], errors="coerce")
    positive = df.loc[df["_p_for_plot"].gt(0), "_p_for_plot"]
    min_positive = positive.min() if not positive.empty else 10 ** (-cap_neg_log10_p)
    zero_floor = min(min_positive, 10 ** (-cap_neg_log10_p))
    df["_neg_log10_p"] = -np.log10(df["_p_for_plot"].replace(0, zero_floor))
    df["_neg_log10_p"] = df["_neg_log10_p"].clip(upper=cap_neg_log10_p)

    pivot = df.pivot_table(
        index="_row_id",
        columns="_model_label",
        values="_neg_log10_p",
        aggfunc="first",
    )
    row_order = [row for row in row_labels.index if row in pivot.index]
    pivot = pivot.reindex(row_order)

    if {"model_set", "predictor_set"}.issubset(df.columns):
        col_order_df = (
            df[["model_set", "predictor_set", "_model_label"]]
            .drop_duplicates()
            .assign(
                _model_set_order=lambda x: (
                    x["model_set"].map({"primary": 0, "expanded": 1}).fillna(99)
                ),
                _predictor_set_order=lambda x: (
                    x["predictor_set"].map({"single": 0, "joint": 1}).fillna(99)
                ),
            )
            .sort_values(
                [
                    "_model_set_order",
                    "_predictor_set_order",
                    "model_set",
                    "predictor_set",
                ]
            )
        )
        pivot = pivot.reindex(columns=col_order_df["_model_label"])

    p_pivot = df.pivot_table(
        index="_row_id",
        columns="_model_label",
        values=p_col,
        aggfunc="first",
    ).reindex(index=pivot.index, columns=pivot.columns)
    annot = p_pivot.apply(
        lambda col: col.map(
            lambda x: (
                ""
                if pd.isna(x)
                else (
                    f"<1e-{int(cap_neg_log10_p)}"
                    if float(x) == 0
                    else _format_p_value(x)
                )
            )
        )
    )
    return pivot, annot, row_labels


def _ordered_column_union(*columns: Iterable[Any]) -> list[Any]:
    ordered: list[Any] = []
    for values in columns:
        for value in values:
            if value not in ordered:
                ordered.append(value)
    return ordered


def _hatch_missing_heatmap_cells(ax: Any, values: pd.DataFrame) -> None:
    for row, col in np.argwhere(values.isna().to_numpy()):
        ax.add_patch(
            Rectangle(
                (col, row),
                1,
                1,
                facecolor="#F1F3F5",
                edgecolor=BORDER,
                hatch="////",
                linewidth=0.6,
                zorder=2,
            )
        )
        ax.text(
            col + 0.5,
            row + 0.5,
            "n/a",
            ha="center",
            va="center",
            fontsize="small",
            color=GRAY_LIGHT,
            zorder=3,
        )


def plot_regression_wald_heatmap(
    composition_wald_df: pd.DataFrame,
    mixing_wald_df: pd.DataFrame | None = None,
    *,
    model_set: str | Iterable[str] | None = None,
    predictor_set: str | Iterable[str] | None = None,
    predictors: Iterable[str] | None = None,
    p_col: str = "p_adj_bh",
    row_col: str | None = None,
    label_map: Mapping[str, str] | None = None,
    cap_neg_log10_p: float = 20.0,
    annotate_p: bool = True,
    title: str | None = None,
    panel_titles: tuple[str, str] = ("Composition", "Mixing"),
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float | None = None,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Heatmap of Wald evidence across regression specifications.

    Passing both composition and mixing Wald tables draws a two-panel figure
    with composition on top and mixing below. Passing only the first table keeps
    the historical single-panel behaviour. The colour scale is ``-log10(p)``
    using the BH-adjusted column when available. Exact zero p-values are capped
    for display and annotated as ``<1e-k`` rather than being dropped.
    """
    if mixing_wald_df is None:
        pivot, annot, row_labels = _prepare_regression_wald_heatmap(
            composition_wald_df,
            model_set=model_set,
            predictor_set=predictor_set,
            predictors=predictors,
            p_col=p_col,
            row_col=row_col,
            label_map=label_map,
            cap_neg_log10_p=cap_neg_log10_p,
        )
        if height_in is None:
            height_in = max(2.5, 0.35 * len(pivot.index) + 1.2)

        fig, ax = new_figure(
            width=width,
            width_in=width_in,
            height_in=height_in,
            context=context,
            font_scale=font_scale,
        )
        sns.heatmap(
            pivot,
            ax=ax,
            cmap="YlOrRd",
            vmin=0,
            vmax=cap_neg_log10_p,
            linewidths=0.5,
            linecolor="white",
            cbar_kws={"label": "Evidence strength (-log10 p-value)"},
            annot=annot if annotate_p else False,
            fmt="" if annotate_p else ".2f",
        )
        _hatch_missing_heatmap_cells(ax, pivot)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_yticklabels(
            [row_labels.get(row, row) for row in pivot.index],
            rotation=0,
        )
        ax.tick_params(axis="x", rotation=0)
        if title is not None:
            ax.set_title(title)

        plt.close(fig)
        return fig

    composition_pivot, composition_annot, composition_labels = (
        _prepare_regression_wald_heatmap(
            composition_wald_df,
            domain="composition",
            model_set=model_set,
            predictor_set=predictor_set,
            predictors=predictors,
            p_col=p_col,
            row_col=row_col,
            label_map=label_map,
            cap_neg_log10_p=cap_neg_log10_p,
        )
    )
    mixing_pivot, mixing_annot, mixing_labels = _prepare_regression_wald_heatmap(
        mixing_wald_df,
        domain="node_mixing",
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
        p_col=p_col,
        row_col=row_col,
        label_map=label_map,
        cap_neg_log10_p=cap_neg_log10_p,
    )

    col_order = _ordered_column_union(composition_pivot.columns, mixing_pivot.columns)
    composition_pivot = composition_pivot.reindex(columns=col_order)
    mixing_pivot = mixing_pivot.reindex(columns=col_order)
    composition_annot = composition_annot.reindex(
        index=composition_pivot.index,
        columns=col_order,
    ).fillna("")
    mixing_annot = mixing_annot.reindex(
        index=mixing_pivot.index,
        columns=col_order,
    ).fillna("")

    if height_in is None:
        n_rows = len(composition_pivot.index) + len(mixing_pivot.index)
        height_in = max(3.6, 0.35 * n_rows + 1.6)

    fig, axes = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        nrows=2,
        ncols=1,
        context=context,
        font_scale=font_scale,
        gridspec_kw={
            "height_ratios": [
                max(1, len(composition_pivot.index)),
                max(1, len(mixing_pivot.index)),
            ],
            "hspace": 0.10,
        },
        layout="constrained",
    )

    axes = axes.ravel()
    cmap = plt.get_cmap("YlOrRd")
    for ax, pivot, annot, row_labels, panel_title in [
        (
            axes[0],
            composition_pivot,
            composition_annot,
            composition_labels,
            panel_titles[0],
        ),
        (axes[1], mixing_pivot, mixing_annot, mixing_labels, panel_titles[1]),
    ]:
        sns.heatmap(
            pivot,
            ax=ax,
            cmap=cmap,
            vmin=0,
            vmax=cap_neg_log10_p,
            linewidths=0.5,
            linecolor="white",
            cbar=False,
            annot=annot if annotate_p else False,
            fmt="" if annotate_p else ".2f",
        )
        _hatch_missing_heatmap_cells(ax, pivot)
        ax.set_title(panel_title)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_yticklabels(
            [row_labels.get(row, row) for row in pivot.index],
            rotation=0,
        )
        ax.tick_params(axis="x", rotation=0)

    axes[0].tick_params(axis="x", labelbottom=False)
    sm = ScalarMappable(norm=Normalize(0, cap_neg_log10_p), cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.86, pad=0.02)
    cbar.set_label("Evidence strength (-log10 p-value)")

    if title is not None:
        fig.suptitle(title)

    add_panel_labels(axes)

    plt.close(fig)
    return fig


# Focused odds-ratio displays ------------------------------------------------


def _read_table(table: pd.DataFrame | str | Any) -> pd.DataFrame:
    """Return a dataframe from an in-memory table or a CSV/parquet path."""
    if isinstance(table, pd.DataFrame):
        return table.copy()
    path = str(table)
    if path.endswith(".parquet"):
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _coerce_shapely_geometry(value: Any):
    if isinstance(value, (bytes, bytearray, memoryview)):
        return wkb.loads(bytes(value))
    return value


def load_health_board_geometries(
    geography: pd.DataFrame | str | Any,
    *,
    board_col: str = "dz_health_board",
    geometry_col: str = "geometry",
) -> dict[str, Any]:
    """Dissolve data-zone geometries into health-board geometries.

    ``geography`` may be a dataframe or a parquet path. Geometry values may
    already be shapely objects or WKB bytes, as in the processed geography
    parquet used by the project.
    """
    if isinstance(geography, pd.DataFrame):
        df = geography[[geometry_col, board_col]].copy()
    else:
        df = pd.read_parquet(geography, columns=[geometry_col, board_col])

    df = df.dropna(subset=[geometry_col, board_col])
    geoms: dict[str, Any] = {}
    for board, sub in df.groupby(board_col):
        geoms[str(board)] = unary_union(
            [_coerce_shapely_geometry(value) for value in sub[geometry_col]]
        )
    return geoms


def _polygon_rings(geom: Any) -> list[np.ndarray]:
    members = (
        [geom]
        if isinstance(geom, Polygon)
        else list(geom.geoms)
        if isinstance(geom, MultiPolygon)
        else list(getattr(geom, "geoms", [geom]))
    )
    return [
        np.asarray(poly.exterior.coords)
        for poly in members
        if isinstance(poly, Polygon)
    ]


def _text_on_fill(color: Any) -> str:
    r, g, b = to_rgb(color)
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "#FFFFFF" if luminance < 0.55 else INK


def _or_effect_color(
    odds: Any,
    low: Any = np.nan,
    high: Any = np.nan,
    *,
    is_reference: bool = False,
) -> str:
    """Color ORs by direction when the CI excludes one."""
    if is_reference:
        return REFERENCE_COLOR

    odds_value = pd.to_numeric(pd.Series([odds]), errors="coerce").iloc[0]
    low_value = pd.to_numeric(pd.Series([low]), errors="coerce").iloc[0]
    high_value = pd.to_numeric(pd.Series([high]), errors="coerce").iloc[0]
    if pd.isna(odds_value):
        return GRAY

    has_ci = pd.notna(low_value) and pd.notna(high_value)
    if has_ci and min(low_value, high_value) <= 1.0 <= max(low_value, high_value):
        return GRAY
    return ORANGE_DARK if odds_value > 1 else TEAL_DARK


def _odds_ratio_lookup(
    odds_df: pd.DataFrame,
    term_token: str,
    *,
    reference: str,
) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    subset = odds_df.loc[
        odds_df["term"].astype(str).str.contains(term_token, regex=False)
    ]
    odds = dict(zip(subset["term"].map(_term_level), subset["odds_ratio"]))
    low = dict(zip(subset["term"].map(_term_level), subset["or_low"]))
    high = dict(zip(subset["term"].map(_term_level), subset["or_high"]))
    odds[reference] = 1.0
    low[reference] = 1.0
    high[reference] = 1.0
    return odds, low, high


def _reference_from_table(
    odds_df: pd.DataFrame,
    predictor: str,
    *,
    default: str,
    aliases: Iterable[str] = (),
) -> str:
    """Best-effort reference lookup for single and joint composition tables."""
    if "reference" not in odds_df.columns:
        return default

    keys = [predictor, *aliases]
    ref_source = odds_df
    if "predictor" in odds_df.columns:
        subset = odds_df.loc[odds_df["predictor"].astype(str).eq(predictor)]
        if not subset.empty:
            ref_source = subset

    for value in ref_source["reference"]:
        if value is None:
            continue
        if hasattr(value, "get"):
            for key in keys:
                ref_value = value.get(key)
                if ref_value is not None:
                    return str(ref_value)
            continue
        if pd.isna(value):
            continue

        text = str(value).strip()
        if not text:
            continue
        if text.startswith("{") and text.endswith("}"):
            try:
                parsed = ast.literal_eval(text)
            except (SyntaxError, ValueError):
                parsed = None
            if parsed is not None and hasattr(parsed, "get"):
                for key in keys:
                    ref_value = parsed.get(key)
                    if ref_value is not None:
                        return str(ref_value)
                continue
        return text

    return default


def _categorical_or_frame(
    odds_df: pd.DataFrame,
    term_token: str,
    *,
    reference: str,
    order: Iterable[str] | None = None,
    label_map: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Return ordered OR rows including the reference level."""
    reference = str(reference)
    odds, low, high = _odds_ratio_lookup(
        odds_df,
        term_token,
        reference=reference,
    )
    if order is None:
        term_levels = (
            odds_df.loc[
                odds_df["term"].astype(str).str.contains(term_token, regex=False),
                "term",
            ]
            .map(_term_level)
            .astype(str)
            .tolist()
        )
        order_values = [reference, *term_levels]
    else:
        order_values = [str(value) for value in order]
        if reference not in order_values:
            order_values.insert(0, reference)

    ordered = list(dict.fromkeys(order_values))
    ordered.extend(level for level in odds if level not in ordered)

    records = []
    for level in ordered:
        label = label_map.get(level, level) if label_map else level
        records.append(
            {
                "level": level,
                "label": f"{label} (ref)" if level == reference else label,
                "odds_ratio": odds.get(level, np.nan),
                "or_low": low.get(level, np.nan),
                "or_high": high.get(level, np.nan),
                "is_reference": level == reference,
            }
        )

    out = pd.DataFrame(records)
    for col in ["odds_ratio", "or_low", "or_high"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["odds_ratio"])
    return out


_DEFAULT_MIXING_OR_ORDER = (
    "sex_entropy_z",
    "age_entropy_z",
    "simd_entropy_z",
    "urban_rural_entropy_z",
    "health_board_entropy_z",
    "sex_entropy_obs_x10",
    "age_entropy_obs_x10",
    "simd_entropy_obs_x10",
    "urban_rural_entropy_obs_x10",
    "health_board_entropy_obs_x10",
    "sex_entropy_obs",
    "age_entropy_obs",
    "simd_entropy_obs",
    "urban_rural_entropy_obs",
    "health_board_entropy_obs",
)


def _mixing_or_label(value: Any, label_map: Mapping[str, str] | None = None) -> str:
    text = str(value)
    if label_map and text in label_map:
        return label_map[text]
    canonical = (
        text.replace("_entropy_z", "")
        .replace("_entropy_obs_x10", "")
        .replace("_entropy_obs", "")
    )
    labels = {
        "sex": "Sex entropy",
        "age": "Age entropy",
        "simd": "SIMD entropy",
        "urban_rural": "Urban/rural entropy",
        "health_board": "Health-board entropy",
    }
    return labels.get(canonical, _pretty_text(text, label_map))


def _mixing_or_frame(
    odds_df: pd.DataFrame,
    *,
    model_set: str,
    predictor_set: str,
    order: Iterable[str] | None = None,
    label_map: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Return ordered continuous mixing OR rows."""
    data = _filter_regression_table(
        odds_df,
        domain="node_mixing",
        model_set=model_set,
        predictor_set=predictor_set,
    )
    if predictor_set == "joint" and "predictor" in data.columns:
        joint = data.loc[data["predictor"].eq("all_mixing")]
        if not joint.empty:
            data = joint
    if data.empty:
        return pd.DataFrame()

    data = data.copy()
    for col in ["odds_ratio", "or_low", "or_high"]:
        data[col] = pd.to_numeric(data[col], errors="coerce")

    observed_terms = list(dict.fromkeys(data["term"].astype(str).tolist()))
    if order is None:
        order = [term for term in _DEFAULT_MIXING_OR_ORDER if term in observed_terms]
    ordered_terms = list(dict.fromkeys([str(term) for term in order]))
    ordered_terms.extend(term for term in observed_terms if term not in ordered_terms)

    rows = []
    data = data.set_index(data["term"].astype(str), drop=False)
    for term in ordered_terms:
        if term not in data.index:
            continue
        row = data.loc[term]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        rows.append(
            {
                "level": term,
                "label": _mixing_or_label(term, label_map),
                "odds_ratio": row["odds_ratio"],
                "or_low": row["or_low"],
                "or_high": row["or_high"],
                "is_reference": False,
            }
        )

    out = pd.DataFrame(rows)
    return out.replace([np.inf, -np.inf], np.nan).dropna(subset=["odds_ratio"])


def _plot_forest_panel(
    ax: Any,
    panel: pd.DataFrame,
    *,
    title: str,
    xlabel: str | None = None,
    annotate_values: bool = True,
) -> None:
    y_positions = np.arange(len(panel))[::-1]
    for y, row in zip(y_positions, panel.itertuples(index=False)):
        odds = row.odds_ratio
        has_ci = not (pd.isna(row.or_low) or pd.isna(row.or_high))
        low = row.or_low if has_ci else odds
        high = row.or_high if has_ci else odds
        ci_overlaps_one = (
            True if not has_ci else min(low, high) <= 1.0 <= max(low, high)  # type: ignore
        )
        color = _or_effect_color(
            odds,
            low,
            high,
            is_reference=bool(row.is_reference),
        )
        if row.is_reference:
            ax.plot(1.0, y, marker="o", ms=7.5, color=color, zorder=3)
        else:
            ax.plot([low, high], [y, y], color=color, lw=2.2, alpha=0.40, zorder=2)
            ax.plot([1.0, odds], [y, y], color=color, lw=0.9, ls=":", alpha=0.6)
            ax.plot(odds, y, marker="o", ms=7.5, color=color, zorder=3)
        if annotate_values:
            ax.annotate(
                f"{odds:.2f}",
                (odds, y),
                textcoords="offset points",
                xytext=(0, 9),
                ha="center",
                va="center",
                fontsize="small",
                color=color,
                fontweight=(
                    "bold" if not row.is_reference and not ci_overlaps_one
                    else "normal"
                ),
            )

    ax.axvline(1.0, color=GRAY_LIGHT, lw=1.0, ls="--", zorder=0)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(panel["label"].tolist())
    ax.set_ylim(-0.6, len(panel) - 0.4)
    ax.set_xlabel(xlabel or "")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.set_axisbelow(True)


def plot_health_board_enrichment_map(
    geography: pd.DataFrame | str | Any,
    odds_ratios: pd.DataFrame | str | Any,
    *,
    board_col: str = "dz_health_board",
    geometry_col: str = "geometry",
    model_set: str = "primary",
    predictor_set: str = "joint",
    reference_board: str = "Greater Glasgow and Clyde",
    reference_urban_rural: str = "Large Urban Areas",
    board_or_limits: tuple[float, float, float] | None = None,
    annotate_top_n: int = 3,
    annotate_all_boards: bool = False,
    show_annotation_ci: bool = True,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 4.2,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Manuscript map of health-board ORs with an urban/rural companion panel.

    Health-board annotations are selected from the largest ORs above and below
    one, excluding the reference board. Confidence intervals are shown in the
    annotation text when available.
    """
    geoms = load_health_board_geometries(
        geography,
        board_col=board_col,
        geometry_col=geometry_col,
    )
    odds_df = _read_table(odds_ratios)
    odds_df = _filter_regression_table(
        odds_df,
        domain="composition",
        model_set=model_set,
        predictor_set=predictor_set,
    )
    if odds_df.empty:
        raise ValueError("No matching composition odds-ratio rows remain.")

    if predictor_set == "joint" and "predictor" in odds_df.columns:
        joint = odds_df.loc[odds_df["predictor"].eq("all_composition")]
        if not joint.empty:
            odds_df = joint

    hb_or, hb_low, hb_high = _odds_ratio_lookup(
        odds_df,
        "health_board",
        reference=reference_board,
    )
    ur_or, ur_low, ur_high = _odds_ratio_lookup(
        odds_df,
        "urban_rural",
        reference=reference_urban_rural,
    )

    if board_or_limits is None:
        finite_board_or = np.asarray(
            [value for value in hb_or.values() if np.isfinite(value)],
            dtype=float,
        )
        if finite_board_or.size == 0:
            finite_board_or = np.asarray([1.0])
        vmin = min(float(finite_board_or.min()), 1.0)
        vmax = max(float(finite_board_or.max()), 1.0)
        pad = max((vmax - vmin) * 0.05, 0.02)
        board_or_limits = (max(0.01, vmin - pad), 1.0, vmax + pad)

    norm = TwoSlopeNorm(
        vmin=board_or_limits[0],
        vcenter=board_or_limits[1],
        vmax=board_or_limits[2],
    )
    cmap = OR_DIVERGING

    fig, axes = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        nrows=1,
        ncols=2,
        context=context,
        font_scale=font_scale,
        gridspec_kw={"width_ratios": [1.12, 0.88], "wspace": 0.2},
        layout="constrained",
    )

    ax_map, ax_or = axes.ravel()
    ax_map.set_aspect("equal")
    ax_map.axis("off")

    reps = {}
    reference_fill = REFERENCE_COLOR
    for board, geom in geoms.items():
        fill = (
            reference_fill
            if board == reference_board
            else cmap(norm(hb_or.get(board, 1.0)))
        )
        ax_map.add_collection(
            PolyCollection(
                _polygon_rings(geom),
                facecolors=[fill],
                edgecolors="#FFFFFF",
                linewidths=0.6,
                zorder=2,
            )
        )
        reps[board] = geom.representative_point()

    all_x = np.concatenate(
        [ring[:, 0] for geom in geoms.values() for ring in _polygon_rings(geom)]
    )
    all_y = np.concatenate(
        [ring[:, 1] for geom in geoms.values() for ring in _polygon_rings(geom)]
    )
    span_x = all_x.max() - all_x.min()
    span_y = all_y.max() - all_y.min()
    label_margin = 0.30 if annotate_all_boards or annotate_top_n > 0 else 0.05
    ax_map.set_xlim(
        all_x.min() - label_margin * span_x,
        all_x.max() + label_margin * span_x,
    )
    ax_map.set_ylim(all_y.min() - 0.03 * span_y, all_y.max() + 0.03 * span_y)

    if reference_board in reps:
        point = reps[reference_board]
        ax_map.text(
            point.x,
            point.y,
            "ref",
            ha="center",
            va="center",
            fontsize="small",
            fontweight="bold",
            color=_text_on_fill(reference_fill),
            zorder=4,
        )

    def _format_board_annotation(row: pd.Series) -> str:
        text = f"{row['board']}\nOR {row['odds_ratio']:.2f}"
        if (
            show_annotation_ci
            and pd.notna(row["or_low"])
            and pd.notna(row["or_high"])
        ):
            text += f"\n[{row['or_low']:.2f}--{row['or_high']:.2f}]"
        return text

    hb_records = []
    for board in geoms:
        if board == reference_board or board not in reps:
            continue
        odds = hb_or.get(board, np.nan)
        low = hb_low.get(board, np.nan)
        high = hb_high.get(board, np.nan)
        if not np.isfinite(odds):
            continue
        point = reps[board]
        hb_records.append(
            {
                "board": board,
                "odds_ratio": float(odds),
                "or_low": low,
                "or_high": high,
                "point_x": point.x,
                "point_y": point.y,
                "color": _or_effect_color(odds, low, high),
            }
        )
    hb_annotations = pd.DataFrame(hb_records)
    if annotate_all_boards:
        label_rows = hb_annotations.sort_values("point_y", ascending=False)
    elif annotate_top_n > 0 and not hb_annotations.empty:
        high_rows = (
            hb_annotations.loc[hb_annotations["odds_ratio"].gt(1.0)]
            .sort_values("odds_ratio", ascending=False)
            .head(annotate_top_n)
        )
        low_rows = (
            hb_annotations.loc[hb_annotations["odds_ratio"].lt(1.0)]
            .sort_values("odds_ratio", ascending=True)
            .head(annotate_top_n)
        )
        label_rows = (
            pd.concat([high_rows, low_rows], ignore_index=True)
            .drop_duplicates("board")
            .sort_values(["odds_ratio", "board"], ascending=[False, True])
        )
    else:
        label_rows = pd.DataFrame()

    if not label_rows.empty:
        x_left = all_x.min() - 0.24 * span_x
        x_right = all_x.max() + 0.24 * span_x
        y_top = all_y.max() - 0.03 * span_y
        y_bottom = all_y.min() + 0.07 * span_y
        for side, side_rows in [
            ("right", label_rows.loc[label_rows["odds_ratio"].ge(1.0)]),
            ("left", label_rows.loc[label_rows["odds_ratio"].lt(1.0)]),
        ]:
            if side_rows.empty:
                continue
            side_rows = side_rows.sort_values("point_y", ascending=False)
            y_targets = np.linspace(y_top, y_bottom, len(side_rows))
            for label_y, (_, row) in zip(y_targets, side_rows.iterrows()):
                point = reps[row["board"]]
                label_x = x_right if side == "right" else x_left
                ax_map.annotate(
                    _format_board_annotation(row),
                    xy=(point.x, point.y),
                    xytext=(label_x, label_y),
                    fontsize="small",
                    color=row["color"],
                    fontweight=(
                        "bold" if row["color"] in {ORANGE_DARK, TEAL_DARK}
                        else "normal"
                    ),
                    va="center",
                    ha="left" if side == "right" else "right",
                    linespacing=1.0,
                    arrowprops=dict(
                        arrowstyle="-",
                        color=GRAY_LIGHT,
                        lw=0.9,
                        shrinkA=1,
                        shrinkB=2,
                    ),
                    zorder=5,
                )

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = ax_map.inset_axes((0.5, 0.0, 0.42, 0.03))
    upper_tick = board_or_limits[2]
    board_ticks = [
        board_or_limits[0],
        1.0,
        1.0 + (upper_tick - 1.0) / 2,
        upper_tick,
    ]
    cb = fig.colorbar(
        sm,
        cax=cax,
        orientation="horizontal",
        ticks=board_ticks,
    )
    cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.2g}"))
    cb.ax.tick_params(length=2)
    cb.outline.set_visible(False)  # type: ignore
    cb.set_label(
        f"Candidate odds ratio vs {reference_board}",
        labelpad=3,
    )

    urban_order = [
        "Large Urban Areas",
        "Other Urban Areas",
        "Accessible Small Towns",
        "Remote Small Towns",
        "Accessible Rural",
        "Remote Rural",
    ]
    urban_short = {
        "Large Urban Areas": "Large urban (ref)",
        "Other Urban Areas": "Other urban",
        "Accessible Small Towns": "Accessible towns",
        "Remote Small Towns": "Remote towns",
        "Accessible Rural": "Accessible rural",
        "Remote Rural": "Remote rural",
    }
    y_positions = np.arange(len(urban_order))[::-1]
    for y, key in zip(y_positions, urban_order):
        odds = ur_or.get(key, np.nan)
        if key == reference_urban_rural:
            ax_or.plot(1, y, marker="o", ms=9, color=REFERENCE_COLOR, zorder=3)
            ax_or.annotate(
                "1.00",
                (1, y),
                textcoords="offset points",
                xytext=(0, 11),
                ha="center",
                fontsize="small",
                color=REFERENCE_COLOR,
            )
            continue
        low = ur_low.get(key, odds)
        high = ur_high.get(key, odds)
        color = _or_effect_color(odds, low, high)
        ax_or.plot([low, high], [y, y], color=color, lw=2.6, alpha=0.40, zorder=2)
        ax_or.plot([1, odds], [y, y], color=color, lw=1.0, ls=":", alpha=0.6, zorder=1)
        ax_or.plot(odds, y, marker="o", ms=9, color=color, zorder=3)
        ax_or.annotate(
            f"{odds:.2f}",
            (odds, y),
            textcoords="offset points",
            xytext=(0, 11),
            ha="center",
            fontsize="small",
            color=color,
            fontweight="bold" if color in {ORANGE_DARK, TEAL_DARK} else "normal",
        )

    ax_or.axvline(1.0, color=GRAY_LIGHT, lw=1.0, ls="--", zorder=0)
    ax_or.set_yticks(y_positions)
    ax_or.set_yticklabels([urban_short[key] for key in urban_order])
    finite_urban_or = [
        value
        for mapping in (ur_low, ur_or, ur_high)
        for value in mapping.values()
        if np.isfinite(value)
    ]
    x_min = min(0.98, min(finite_urban_or, default=1.0) - 0.01)
    x_max = max(1.12, max(finite_urban_or, default=1.0) + 0.02)
    ax_or.set_xlim(x_min, x_max)
    ax_or.set_ylim(-0.6, len(urban_order) - 0.25)
    ax_or.set_xlabel("Candidate odds ratio vs large urban areas")
    ax_or.spines["left"].set_visible(False)
    ax_or.tick_params(axis="y", length=0)
    ax_or.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax_or.set_axisbelow(True)

    add_panel_labels([ax_map, ax_or])
    plt.close(fig)
    return fig


def plot_age_sex_simd_forest(
    composition_odds_ratios: pd.DataFrame | str | Any,
    mixing_odds_ratios: pd.DataFrame | str | Any,
    *,
    model_set: str = "primary",
    predictor_set: str = "joint",
    mixing_predictor_set: str | None = None,
    reference_age: str | None = None,
    reference_sex: str | None = None,
    reference_simd: str | None = None,
    age_order: Iterable[str] | None = None,
    sex_order: Iterable[str] | None = None,
    simd_order: Iterable[str] | None = None,
    mixing_order: Iterable[str] | None = None,
    simd_label_map: Mapping[str, str] | None = None,
    mixing_label_map: Mapping[str, str] | None = None,
    width_ratios: tuple[float, float] = (1.12, 0.88),
    height_ratios: tuple[float, ...] | None = None,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 4.8,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    annotate_values: bool = True,
) -> Figure:
    """Age, sex, SIMD, and mixing odds-ratio forest figure.

    The layout is a 3x2 grid where age spans the first column and sex, SIMD,
    and mixing occupy the right column.
    """
    composition_df = _read_table(composition_odds_ratios)
    missing = {"term", "odds_ratio", "or_low", "or_high"} - set(composition_df.columns)
    if missing:
        raise ValueError(
            f"`composition_odds_ratios` is missing required columns: {sorted(missing)}"
        )

    odds_df = _filter_regression_table(
        composition_df,
        domain="composition",
        model_set=model_set,
        predictor_set=predictor_set,
    )
    if odds_df.empty:
        raise ValueError("No matching composition odds-ratio rows remain.")

    if predictor_set == "joint" and "predictor" in odds_df.columns:
        joint = odds_df.loc[odds_df["predictor"].eq("all_composition")]
        if not joint.empty:
            odds_df = joint

    reference_age = reference_age or _reference_from_table(
        odds_df,
        "age_band",
        default="20-24",
    )
    reference_sex = reference_sex or _reference_from_table(
        odds_df,
        "sex",
        default="Male",
    )
    reference_simd = reference_simd or _reference_from_table(
        odds_df,
        "simd_quintile",
        default="1",
        aliases=("dz_simd_quintile",),
    )

    if age_order is None:
        age_order = [
            "00-04",
            "05-09",
            "10-14",
            "15-19",
            "20-24",
            "25-29",
            "30-34",
            "35-39",
            "40-44",
            "45-49",
            "50-54",
            "55-59",
            "60-64",
            "65-69",
            "70-74",
            "75+",
        ]
    if sex_order is None:
        sex_order = ["Male", "Female"]
    if simd_order is None:
        simd_order = ["1", "2", "3", "4", "5"]
    if simd_label_map is None:
        simd_label_map = {
            "1": "1 most deprived",
            "2": "2",
            "3": "3",
            "4": "4",
            "5": "5 least deprived",
        }

    age_panel = _categorical_or_frame(
        odds_df,
        "age_band",
        reference=reference_age,
        order=age_order,
    )
    sex_panel = _categorical_or_frame(
        odds_df,
        "sex",
        reference=reference_sex,
        order=sex_order,
    )
    simd_panel = _categorical_or_frame(
        odds_df,
        "simd_quintile",
        reference=reference_simd,
        order=simd_order,
        label_map=simd_label_map,
    )
    mixing_df = _read_table(mixing_odds_ratios)
    mixing_panel = _mixing_or_frame(
        mixing_df,
        model_set=model_set,
        predictor_set=mixing_predictor_set or predictor_set,
        order=mixing_order,
        label_map=mixing_label_map,
    )
    if mixing_panel.empty:
        raise ValueError("No matching mixing odds-ratio rows remain.")

    if age_panel.empty or sex_panel.empty or simd_panel.empty:
        raise ValueError(
            "Age, sex, and SIMD panels must each contain at least one row."
        )

    nrows = 3
    if height_ratios is None:
        height_ratios = (0.55, 1.0, 1.0)
    if len(height_ratios) != nrows:
        raise ValueError(f"`height_ratios` must contain {nrows} values.")

    fig, axes = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        nrows=nrows,
        ncols=2,
        sharex=True,
        context=context,
        font_scale=font_scale,
        gridspec_kw={
            "width_ratios": list(width_ratios),
            "height_ratios": list(height_ratios),
            "wspace": 0.10,
        },
        layout="constrained",
    )
    grid = axes[0, 0].get_gridspec()
    for row_idx in range(nrows):
        axes[row_idx, 0].remove()
    ax_age = fig.add_subplot(grid[:, 0])
    ax_sex = axes[0, 1]
    ax_simd = axes[1, 1]
    ax_mixing = axes[2, 1]

    xlabel = "Candidate odds ratio vs reference level"
    _plot_forest_panel(
        ax_age,
        age_panel,
        title="",
        xlabel=xlabel,
        annotate_values=annotate_values,
    )
    _plot_forest_panel(
        ax_sex,
        sex_panel,
        title="",
        annotate_values=annotate_values,
    )
    _plot_forest_panel(
        ax_simd,
        simd_panel,
        title="",
        xlabel="",
        annotate_values=annotate_values,
    )
    ax_sex.tick_params(axis="x", labelbottom=False)
    _plot_forest_panel(
        ax_mixing,
        mixing_panel,
        title="",
        xlabel=xlabel,
        annotate_values=annotate_values,
    )
    ax_simd.tick_params(axis="x", labelbottom=False)

    add_panel_labels([ax_age, ax_sex, ax_simd, ax_mixing])
    plt.close(fig)
    return fig


# ---------------------------------------------------------------------------
# Sensitivity matrix figure
# ---------------------------------------------------------------------------


_DEFAULT_SENSITIVITY_FAMILIES = [
    ("Main", "association_outputs"),
    ("Clade\ngrouped", "sensitivity_clade"),
    ("Window\nadjusters", "sensitivity_window"),
    ("Observed\nentropy", "sensitivity_observed_entropy"),
]

_SENSITIVITY_COMPOSITION_ROWS = [
    ("health_board", "Health board"),
    ("urban_rural_class", "Urban/rural class"),
    ("simd_quintile", "SIMD quintile"),
    ("age_band", "Age band"),
    ("sex", "Sex"),
]

_SENSITIVITY_MIXING_ROWS = [
    ("age", "Age entropy"),
    ("health_board", "Health-board entropy"),
    ("simd", "SIMD entropy"),
    ("urban_rural", "Urban/rural entropy"),
    ("sex", "Sex entropy"),
    ("all_mixing", "All mixing predictors"),
]


def _canonical_mixing_predictor(value: Any) -> str:
    return (
        str(value)
        .replace("_entropy_z", "")
        .replace("_entropy_obs_x10", "")
        .replace("_entropy_obs", "")
    )


def collect_sensitivity_matrix_results(
    results_root: str | Any,
    *,
    families: list[tuple[str, str]] = _DEFAULT_SENSITIVITY_FAMILIES,
    p_col: str = "p_adj_bh",
    alpha: float = 0.05,
) -> tuple[
    dict[str, dict[str, tuple[int, int]] | None],
    dict[str, dict[str, tuple[int, int]] | None],
]:
    """Collect significant-predictor counts for the sensitivity matrix."""
    root = pd.io.common.stringify_path(results_root)  # type: ignore
    composition: dict[str, dict[str, tuple[int, int]] | None] = {}
    mixing: dict[str, dict[str, tuple[int, int]] | None] = {}

    for family, subdir in families:
        path = f"{root}/{subdir}"
        comp_path = f"{path}/composition_wald.csv"
        if pd.io.common.file_exists(comp_path):  # type: ignore
            comp = pd.read_csv(comp_path)
            comp = comp.loc[comp["predictor"].ne("all_composition")].copy()
            composition[family] = {  # type: ignore
                predictor: (int((group[p_col] < alpha).sum()), len(group))
                for predictor, group in comp.groupby("predictor")
            }
        else:
            composition[family] = None

        mix_path = f"{path}/mixing_wald.csv"
        if pd.io.common.file_exists(mix_path):  # type: ignore
            mix = pd.read_csv(mix_path)
            source = mix["predictor"].astype(str)
            if "term" in mix.columns:
                source = source.where(~source.eq("all_mixing"), mix["term"].astype(str))
            mix["_canonical_predictor"] = source.map(_canonical_mixing_predictor)
            mixing[family] = {  # type: ignore
                predictor: (int((group[p_col] < alpha).sum()), len(group))
                for predictor, group in mix.groupby("_canonical_predictor")
            }
        else:
            mixing[family] = None

    return composition, mixing


def _draw_sensitivity_block(
    ax,
    rows: list[tuple[str, str]],
    data: Mapping[str, Mapping[str, tuple[int, int]] | None],
    *,
    families: list[tuple[str, str]],
    y_top: float,
    row_h: float,
    col_x: list[float],
    col_w: float,
    cmap,
    norm,
    cell_size: float | str,
    percent_size: float | str,
) -> None:
    for row_idx, (key, label) in enumerate(rows):
        y_center = y_top - (row_idx + 0.5) * row_h
        ax.text(
            col_x[0] - 0.4,
            y_center,
            label,
            ha="right",
            va="center",
            color=INK_SOFT,
        )
        for col_idx, (family, _) in enumerate(families):
            x = col_x[col_idx]
            cell = data.get(family)
            if cell is None or key not in cell:
                ax.add_patch(
                    Rectangle(
                        (x, y_center - row_h / 2 + 0.06),
                        col_w,
                        row_h - 0.12,
                        fc="#F1F3F5",
                        ec=BORDER,
                        lw=0.6,
                        hatch="////",
                        zorder=2,
                    )
                )
                ax.text(
                    x + col_w / 2,
                    y_center,
                    "n/a",
                    ha="center",
                    va="center",
                    fontsize=percent_size,
                    color=GRAY_LIGHT,
                    zorder=3,
                )
                continue

            n_sig, n_total = cell[key]
            share = n_sig / n_total if n_total else 0
            face = cmap(norm(share))
            ax.add_patch(
                Rectangle(
                    (x, y_center - row_h / 2 + 0.06),
                    col_w,
                    row_h - 0.12,
                    fc=face,
                    ec="#FFFFFF",
                    lw=1.2,
                    zorder=2,
                )
            )
            text_color = "#FFFFFF" if share > 0.5 else (GRAY if share == 0 else INK)
            ax.text(
                x + col_w / 2,
                y_center + 0.10,
                f"{n_sig}/{n_total}",
                ha="center",
                va="center",
                fontsize=cell_size,
                fontweight="bold",
                color=text_color,
                zorder=3,
            )
            ax.text(
                x + col_w / 2,
                y_center - row_h * 0.28,
                f"{share * 100:.0f}%",
                ha="center",
                va="center",
                fontsize=percent_size,
                color=text_color,
                alpha=0.9,
                zorder=3,
            )


def plot_sensitivity_matrix(
    results_root: str | Any,
    *,
    families: list[tuple[str, str]] = _DEFAULT_SENSITIVITY_FAMILIES,
    composition_rows: list[tuple[str, str]] = _SENSITIVITY_COMPOSITION_ROWS,
    mixing_rows: list[tuple[str, str]] = _SENSITIVITY_MIXING_ROWS,
    p_col: str = "p_adj_bh",
    alpha: float = 0.05,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 4.8,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Supplementary manuscript matrix of sensitivity-analysis recurrence."""
    composition, mixing = collect_sensitivity_matrix_results(
        results_root,
        families=families,
        p_col=p_col,
        alpha=alpha,
    )
    cmap = WARM_SEQ
    norm = Normalize(0, 1)

    fig, ax = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
        layout="constrained",
    )
    ax.set_xlim(2.0, 12)
    ax.set_ylim(0, 12)
    ax.axis("off")

    left = 4.6
    col_w = (12 - left - 0.2) / len(families)
    col_x = [left + i * col_w for i in range(len(families))]
    row_h = 0.92

    for col_idx, (family, _) in enumerate(families):
        ax.text(
            col_x[col_idx] + col_w / 2,
            11.35,
            family,
            ha="center",
            va="center",
            fontweight="bold",
            color=INK,
            linespacing=1.0,
        )

    comp_top = 10.7
    # ax.text(
    #     left - 0.4,
    #     comp_top + 0.28,
    #     "COMPOSITION",
    #     ha="right",
    #     va="bottom",
    #     fontweight="bold",
    #     color=TEAL_DARK,
    # )
    _draw_sensitivity_block(
        ax,
        composition_rows,
        composition,
        families=families,
        y_top=comp_top,
        row_h=row_h,
        col_x=col_x,
        col_w=col_w,
        cmap=cmap,
        norm=norm,
        cell_size="medium",
        percent_size="small",
    )

    mix_top = comp_top - len(composition_rows) * row_h - 0.55
    ax.plot([left - 0.1, 12], [mix_top + 0.30, mix_top + 0.30], color=BORDER, lw=1.0)
    # ax.text(
    #     left - 0.4,
    #     mix_top + 0.30,
    #     "MIXING",
    #     ha="right",
    #     va="bottom",
    #     fontweight="bold",
    #     color=TEAL_DARK,
    # )
    _draw_sensitivity_block(
        ax,
        mixing_rows,
        mixing,
        families=families,
        y_top=mix_top,
        row_h=row_h,
        col_x=col_x,
        col_w=col_w,
        cmap=cmap,
        norm=norm,
        cell_size="medium",
        percent_size="small",
    )

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = ax.inset_axes((0.5, -0.05, (12 - left) / 12 * 0.5, 0.022))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[0, 0.5, 1.0])
    cb.ax.set_xticklabels(["0%", "50%", "100%"])
    cb.outline.set_visible(False)  # type: ignore
    cb.set_label(
        f"Share of sensitivity models with FDR-adjusted p < {alpha:g}",
        labelpad=4,
    )

    plt.close(fig)
    return fig
