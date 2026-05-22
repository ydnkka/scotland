"""Utilities for the SSE-detection pipeline and its output notebook.

The module is organised so each figure is a single function that takes
ready-to-plot dataframes and returns a matplotlib ``Figure``. The notebook
``sse_detection_plots.ipynb`` is intended to be a thin orchestrator
on top of these functions.

Sub-modules
-----------
stats
    Statistical machinery for the SSE pipeline (entropy z-scores, downstream
    entropy, node-metric assembly, candidate categorisation, weekly growth
    flag).
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
regression
    Reusable regression helpers for association notebooks. These functions
    expect already-prepared model frames and do not perform data cleaning.
"""

from .stats import (
    add_sse_node_metrics,
    categorise_sse_nodes,
    test_category_distribution,
    flag_sse,
    frequencies,
    safe_mode,
)

from .entropy import (
    max_entropy, 
    shannon_entropy, 
    shannon_entropy_grouped,
    cluster_socio_demo_entropy,
    downstream_edge_entropy,
    )

from .io import SseOutputs, load_sse_outputs, load_weekly_growth
from .regression import (
    AssociationModel,
    bh_adjust,
    bounded_exp,
    categorical_term,
    fit_binomial_glm,
    fit_exposure_model,
    make_formula,
    model_fit_stats,
    model_variables_from_terms,
    robust_wald_for_params,
    robust_wald_for_prefix,
    tidy_odds_ratios,
    tidy_single_parameter_wald,
)
from .palettes import (
    DYNAMIC_ORDER,
    DYNAMIC_PALETTE,
    LIFECYCLE_PALETTE,
    ROLE_ORDER,
    ROLE_PALETTE,
    SSE_CATEGORY_PALETTE,
    sse_category_palette_from,
)
from .plots import (
    plot_candidate_rate_over_time,
    plot_cluster_size_distribution,
    plot_composite_distributions,
    plot_core_metric_space,
    plot_role_dynamic_heatmap,
    plot_socio_demo_breakdown,
    plot_socio_demo_candidate_background_diff,
)
from .subgraph import plot_meta_cluster_subgraph

__all__ = [
    "add_sse_node_metrics",
    "categorise_sse_nodes",
    "test_category_distribution",
    "flag_sse",
    "frequencies",
    "safe_mode",
    "max_entropy",
    "shannon_entropy",
    "shannon_entropy_grouped",
    "cluster_socio_demo_entropy",
    "downstream_edge_entropy",
    "SseOutputs",
    "load_sse_outputs",
    "load_weekly_growth",
    "AssociationModel",
    "bh_adjust",
    "bounded_exp",
    "categorical_term",
    "fit_binomial_glm",
    "fit_exposure_model",
    "make_formula",
    "model_fit_stats",
    "model_variables_from_terms",
    "robust_wald_for_params",
    "robust_wald_for_prefix",
    "tidy_odds_ratios",
    "tidy_single_parameter_wald",
    "DYNAMIC_PALETTE",
    "DYNAMIC_ORDER",
    "LIFECYCLE_PALETTE",
    "ROLE_ORDER",
    "ROLE_PALETTE",
    "SSE_CATEGORY_PALETTE",
    "sse_category_palette_from",
    "plot_candidate_rate_over_time",
    "plot_cluster_size_distribution",
    "plot_composite_distributions",
    "plot_meta_cluster_subgraph",
    "plot_core_metric_space",
    "plot_role_dynamic_heatmap",
    "plot_socio_demo_breakdown",
    "plot_socio_demo_candidate_background_diff",
]
