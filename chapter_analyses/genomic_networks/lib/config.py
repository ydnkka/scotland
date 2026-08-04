"""Configuration constants for Chapter 4 observation/network analyses."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CHAPTER_DIR = PROJECT_ROOT / "chapter_analyses/genomic_networks"
RESULTS_DIR = CHAPTER_DIR / "results"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"
INTERMEDIATE_DIR = RESULTS_DIR / "intermediate"

ANALYSIS_RESOLUTION = 0.3
SPARSIFICATION_THRESHOLD = 0.001
TRANSITION_WINDOW_STRIDE = 2
DISCLOSURE_MIN_CELL = 5


@dataclass(frozen=True)
class AttributeSpec:
    """Description of a categorical variable used in Chapter 4 mixing tables."""

    name: str
    column: str
    label: str
    ordered: bool = False


ASSORTATIVITY_ATTRIBUTES: tuple[AttributeSpec, ...] = (
    AttributeSpec("sex", "sex", "Sex"),
    AttributeSpec("age_band", "age_band", "Age band", ordered=True),
    AttributeSpec("age_group", "age_group", "Age group", ordered=True),
    AttributeSpec("simd_quintile", "dz_simd_quintile", "SIMD quintile", ordered=True),
    AttributeSpec("urban_rural", "dz_urban_rural_class", "Urban/rural class"),
    AttributeSpec("local_authority", "dz_local_authority", "Local authority"),
    AttributeSpec("health_board", "dz_health_board", "Health board"),
)

CLUSTER_COLUMNS: tuple[str, ...] = (
    "window_id",
    "window_idx",
    "wn_start_date",
    "wn_mid_date",
    "wn_end_date",
    "wn_no_sequences",
    "wn_positive_tests",
    "wn_prop_sequenced",
    "pango_lineage",
    "cluster_id",
    "cluster_size",
    "cluster_n_datazones",
    "cluster_start_date",
    "cluster_end_date",
    "cluster_duration_days",
)

ANALYSIS_COLUMNS: tuple[str, ...] = (
    *CLUSTER_COLUMNS,
    "sequence_id",
    "patient_id",
    "collection_date",
    "sex",
    "is_female",
    "age_band",
    "age_group",
    "is_vaccinated",
    "vacc_dose_number",
    "vacc_booster",
    "days_since_vaccination",
    "test_reason",
    "is_reinfection",
    "clade",
    "who_voc",
    "policy_era",
    "policy_period",
    "policy_period_label",
    "datazone",
    "dz_xcoord",
    "dz_ycoord",
    "dz_population",
    "dz_population_density",
    "dz_simd_quintile",
    "dz_urban_rural_class",
    "dz_local_authority",
    "dz_health_board",
    "dz_cum_sequences",
    "dz_cum_positive_tests",
    "dz_cum_prop_sequenced",
    "dz_cum_incidence_per_capita",
)
