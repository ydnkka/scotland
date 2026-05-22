"""Utilities for the SSE-detection pipeline and its output notebook.
"""

from .stats import (
    add_sse_node_metrics,
    categorise_sse_nodes,
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
    cluster_se_diagnostics,
    fit_binomial_glm,
    fit_conditional_logit,
    fit_exposure_model,
    fit_firth_logit,
    make_formula,
    model_fit_stats,
    model_variables_from_terms,
    parameter_names_for_term,
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
    make_regression_fit_table,
    make_regression_odds_ratio_table,
    make_regression_wald_table,
    plot_candidate_rate_over_time,
    plot_cluster_size_distribution,
    plot_composite_distributions,
    plot_core_metric_space,
    plot_regression_odds_ratio_forest,
    plot_regression_wald_heatmap,
    plot_role_dynamic_heatmap,
    plot_socio_demo_breakdown,
) 

from .subgraph import plot_meta_cluster_subgraph

__all__ = [
    "add_sse_node_metrics",
    "categorise_sse_nodes",
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
    "cluster_se_diagnostics",
    "fit_binomial_glm",
    "fit_conditional_logit",
    "fit_exposure_model",
    "fit_firth_logit",
    "make_formula",
    "model_fit_stats",
    "model_variables_from_terms",
    "parameter_names_for_term",
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
    "make_regression_fit_table",
    "make_regression_odds_ratio_table",
    "make_regression_wald_table",
    "plot_candidate_rate_over_time",
    "plot_cluster_size_distribution",
    "plot_composite_distributions",
    "plot_meta_cluster_subgraph",
    "plot_core_metric_space",
    "plot_regression_odds_ratio_forest",
    "plot_regression_wald_heatmap",
    "plot_role_dynamic_heatmap",
    "plot_socio_demo_breakdown",
]
