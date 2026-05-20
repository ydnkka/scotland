"""Utilities for the SSE-detection pipeline and its output notebook.

The module is organised so each figure is a single function that takes
ready-to-plot dataframes and returns a matplotlib ``Figure``. The notebook
``sse_detection_output_plots.ipynb`` is intended to be a thin orchestrator
on top of these functions.

Sub-modules
-----------
stats
    Statistical machinery for the SSE pipeline (entropy z-scores, downstream
    entropy, node-metric assembly, candidate categorisation, weekly growth
    flag). These used to live in ``utils/stats.py``.
io
    Loaders for the parquet outputs produced by ``sse_detection.ipynb``.
palettes
    Stable categorical colour maps for the SSE label space.
plots
    Figure functions matching the suggested figure list (overview,
    Layer 1, Layer 2, Layer 1 x Layer 2 concordance, spatial /
    demographic, sensitivity).
subgraph
    Meta-cluster subgraph plotter (colour = ``sse_category``,
    size proportional to ``log1p(cluster_size)``).
"""

from .stats import (
    add_sse_node_metrics,
    attach_entropy_zscore,
    categorise_sse_nodes,
    downstream_entropy_fast,
    downstream_spread_entropy,
    flag_sse,
)
from .io import SseOutputs, load_sse_outputs, load_weekly_growth
from .palettes import (
    DYNAMIC_PALETTE,
    DYNAMIC_ORDER,
    LIFECYCLE_PALETTE,
    ROLE_ORDER,
    ROLE_PALETTE,
    SSE_CATEGORY_PALETTE,
    sse_category_palette_from,
)
from .plots import (
    plot_candidate_rate_over_time,
    plot_cluster_size_ccdf,
    plot_composite_score_distributions,
    plot_layer_concordance,
    plot_metric_space_scatter,
    plot_meta_cluster_trajectories,
    plot_norm_change_histogram,
    plot_null_comparison,
    plot_role_dynamic_heatmap,
    plot_sequence_volume_timeline,
    plot_simd_breakdown,
    plot_threshold_sensitivity,
)
from .subgraph import plot_meta_cluster_subgraph

__all__ = [
    "add_sse_node_metrics",
    "attach_entropy_zscore",
    "categorise_sse_nodes",
    "downstream_entropy_fast",
    "downstream_spread_entropy",
    "flag_sse",
    "SseOutputs",
    "load_sse_outputs",
    "load_weekly_growth",
    "DYNAMIC_PALETTE",
    "DYNAMIC_ORDER",
    "LIFECYCLE_PALETTE",
    "ROLE_ORDER",
    "ROLE_PALETTE",
    "SSE_CATEGORY_PALETTE",
    "sse_category_palette_from",
    "plot_candidate_rate_over_time",
    "plot_cluster_size_ccdf",
    "plot_composite_score_distributions",
    "plot_layer_concordance",
    "plot_meta_cluster_subgraph",
    "plot_meta_cluster_trajectories",
    "plot_metric_space_scatter",
    "plot_norm_change_histogram",
    "plot_null_comparison",
    "plot_role_dynamic_heatmap",
    "plot_sequence_volume_timeline",
    "plot_simd_breakdown",
    "plot_threshold_sensitivity",
]
