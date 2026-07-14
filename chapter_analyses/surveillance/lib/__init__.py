"""Reusable configuration for the surveillance analysis package."""

from .config import FIGURES_DIR, PACKAGE_DIR, PROJECT_ROOT, RESULTS_DIR, TABLES_DIR
from .io import ensure_results_dirs, write_table

__all__ = [
    "FIGURES_DIR",
    "PACKAGE_DIR",
    "PROJECT_ROOT",
    "RESULTS_DIR",
    "TABLES_DIR",
    "ensure_results_dirs",
    "write_table",
]
