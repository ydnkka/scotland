"""Configuration constants for SSE detection analyses."""

from __future__ import annotations

from pathlib import Path


# config.py lives at <project>/chapter_analyses/sse_detection/lib/sse/config.py.
PROJECT_ROOT = Path(__file__).resolve().parents[4]
PACKAGE_DIR = PROJECT_ROOT / "chapter_analyses/sse_detection"
RESULTS_DIR = PACKAGE_DIR / "results"
SSE_OUTPUT_DIR = RESULTS_DIR / "sse_outputs"
BAYESIAN_OUTPUT_DIR = RESULTS_DIR / "bayesian_outputs"
FIGURE_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"

TRANSITION_WINDOW_STRIDE = 2
DETECTION_RANDOM_SEED = 42
N_ENTROPY_DRAWS = 1000
N_PERMUTATIONS = 1000
MIN_CLUSTER_SIZE = 6

ANALYSIS_COLUMNS: tuple[str, ...] = (
    "window_id",
    "window_idx",
    "wn_start_date",
    "wn_mid_date",
    "wn_end_date",
    "wn_positive_tests",
    "wn_no_sequences",
    "wn_prop_sequenced",
    "cluster_id",
    "cluster_size",
    "cluster_n_datazones",
    "cluster_duration_days",
    "age_group",
    "sex",
    "is_female",
    "who_voc",
    "clade",
    "test_reason",
    "is_reinfection",
    "datazone",
    "dz_care_home_tests",
    "dz_urban_rural_class",
    "dz_health_board",
    "dz_local_authority",
    "dz_simd_quintile",
    "dz_7d_test_positivity",
    "dz_cum_sequences",
    "dz_cum_positive_tests",
    "dz_cum_prop_sequenced",
    "dz_cum_incidence_per_capita",
)
