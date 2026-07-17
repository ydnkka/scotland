"""Shared utilities for Scotland clustering chapters.

Submodules
----------
data
    Data loading and preprocessing helpers that use `config.yaml` paths and
    read the processed parquet file in a memory-conscious way.
style
    Publication-grade matplotlib theme, figure size presets, and consistent
    categorical colour palettes for SIMD quintiles, WHO variants, Leiden
    resolutions, and Nextclade clades.
stats
    Light statistical helpers: negative-binomial cluster-size regression,
    bootstrapped confidence intervals, singleton odds ratios, and tidy
    coefficient tables for forest plots.
"""

from . import style, data

from .data import (
    QCStatus,
    VALID_QC_STATUSES,
    load_analysis_columns,
    load_policy_data,
    load_datazone_info,
    load_pairwise_edges,
    pango_lineages_for_clades,
    CLADES,
    CLADE_PALETTE,
    PRIMARY_RESOLUTION,
)
from .style import (
    set_theme,
    save_figure,
    add_panel_labels,
    new_figure,
    lighten,
    FIG_WIDTHS_IN,
    WIDTHS,
    CONTEXTS,
    DEFAULT_HEIGHT_IN,
)

__all__ = [
    "data",
    "style",
    "QCStatus",
    "VALID_QC_STATUSES",
    "load_analysis_columns",
    "load_policy_data",
    "load_datazone_info",
    "load_pairwise_edges",
    "pango_lineages_for_clades",
    "PRIMARY_RESOLUTION",
    "CLADES",
    "CLADE_PALETTE",
    "set_theme",
    "save_figure",
    "add_panel_labels",
    "new_figure",
    "lighten",
    "FIG_WIDTHS_IN",
    "WIDTHS",
    "CONTEXTS",
    "DEFAULT_HEIGHT_IN",
]
