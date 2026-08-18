"""Shared builder registry for top-level publication results."""

from .config import FIGURES_DIR, PROJECT_ROOT, RESULTS_DIR, TABLES_DIR
from .registry import (
    ArtifactBuilder,
    BuildContext,
    build_figures,
    build_tables,
    figure_builders,
    table_builders,
)

__all__ = [
    "ArtifactBuilder",
    "BuildContext",
    "FIGURES_DIR",
    "PROJECT_ROOT",
    "RESULTS_DIR",
    "TABLES_DIR",
    "build_figures",
    "build_tables",
    "figure_builders",
    "table_builders",
]

