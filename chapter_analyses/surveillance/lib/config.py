"""Paths and reproducibility constants for surveillance analyses."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
PACKAGE_DIR = PROJECT_ROOT / "chapter_analyses/surveillance"
RESULTS_DIR = PACKAGE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"

DAILY_SMOOTH_WINDOW = 7
SEQUENCE_WINDOW_STRIDE = 3
FIGURE_NAME = "policy_sequences_over_time"
POLICY_INDEX_FIGURE_NAME = "policy_index_comparison"
