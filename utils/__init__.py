"""Shared utilities for Scotland clustering chapters.
"""

from . import style, data

from .data import (
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
