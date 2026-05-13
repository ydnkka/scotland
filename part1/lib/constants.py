"""Shared constants for the Part 1 analysis.

All labels, term sets, and model specifications used by the data-preparation,
fitting, and plotting modules live here so that they have a single source of
truth.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


# ---------------------------------------------------------------------------
# Run defaults
# ---------------------------------------------------------------------------

QC_DEFAULT = "good"
PRIMARY_RESOLUTION = 0.3
LINEAGE_MIN_CLUSTERS = 50
CALENDAR_SPLINE_DF = 8


# ---------------------------------------------------------------------------
# Columns pulled from the analysis dataset
# ---------------------------------------------------------------------------

SEQUENCE_COLUMNS: Iterable[str] = [
    "cluster_id",
    "sequence_id",
    "resolution",
    "window_id",
    "window_idx",
    "wn_mid_date",
    "wn_prop_sequenced",
    "collection_date",
    "datazone",
    "pango_lineage",
    "nextclade_qc",
    "age_band",
    "sex",
    "dz_simd_rank",
    "dz_simd_quintile",
    "dz_cum_incidence_per_capita",
    "dz_cum_prop_sequenced",
    "dz_7d_test_positivity",
    "wn_no_sequences",
    "dz_health_board_code",
]

# Domain/wave analyses use the same base columns plus the SIMD subdomain ranks.
BASE_SEQUENCE_COLUMNS: Iterable[str] = [
    "cluster_id",
    "sequence_id",
    "resolution",
    "window_id",
    "window_idx",
    "wn_mid_date",
    "wn_prop_sequenced",
    "collection_date",
    "datazone",
    "pango_lineage",
    "nextclade_qc",
    "age_band",
    "sex",
    "dz_simd_quintile",
    "dz_cum_incidence_per_capita",
    "dz_cum_prop_sequenced",
    "dz_7d_test_positivity",
]


# ---------------------------------------------------------------------------
# Primary covariate sets (Line 1 — deprivation)
# ---------------------------------------------------------------------------

PRIMARY_TERMS: Iterable[str] = [
    "deprivation_z",
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
]


# Labels used in tables and plots.  The mixing-predictor labels are appended
# after ``MIXING_VARIABLES`` is defined below.
TERM_LABELS: dict[str, str] = {
    "deprivation_z": "Mean SIMD deprivation",
    "index_deprivation_z": "Index-case SIMD deprivation",
    "local_incidence_z": "Local cumulative incidence",
    "local_seq_fraction_z": "Local sequencing fraction",
    "window_seq_fraction_z": "Window sequencing proportion",
    "test_positivity_z": "Local test positivity",
    "log_cluster_size_z": "Cluster size",
}


# ---------------------------------------------------------------------------
# Mixing variables and predictors (Line 1 outcomes / Line 2 exposures)
# ---------------------------------------------------------------------------

MIXING_VARIABLES: dict[str, dict[str, str]] = {
    "simd": {
        "column": "dz_simd_quintile",
        "label": "SIMD quintile mixing",
        "short_label": "SIMD",
    },
    "age": {
        "column": "age_band",
        "label": "Age-band mixing",
        "short_label": "Age",
    },
    "sex": {
        "column": "sex",
        "label": "Sex mixing",
        "short_label": "Sex",
    },
    "profile": {
        "column": "socio_demographic_profile",
        "label": "Joint SIMD-age-sex profile mixing",
        "short_label": "Joint profile",
    },
}

MIXING_PREDICTOR_TERMS: Iterable[str] = [
    f"{prefix}_excess_mixing_z" for prefix in MIXING_VARIABLES
]

TERM_LABELS.update(
    {
        f"{prefix}_excess_mixing_z": f"{spec['short_label']} excess mixing"
        for prefix, spec in MIXING_VARIABLES.items()
    }
)


# ---------------------------------------------------------------------------
# Count model specifications (Line 1 + Line 2 outcomes)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CountModelSpec:
    name: str
    label: str
    raw_outcome: str
    binary_col: str
    positive_col: str
    positive_label: str
    include_size: bool = False


# Duration is excluded from the primary count outcomes because the fixed
# three-week clustering windows mechanically constrain the observed span.
COUNT_MODEL_SPECS: Iterable[CountModelSpec] = [
    CountModelSpec(
        name="cluster_size",
        label="Cluster size",
        raw_outcome="cluster_size",
        binary_col="cluster_size_gt1",
        positive_col="cluster_size_excess",
        positive_label="Additional sequences among non-singleton clusters",
    ),
    CountModelSpec(
        name="geographic_dispersion",
        label="Geographic dispersion",
        raw_outcome="cluster_n_datazones",
        binary_col="datazones_gt1",
        positive_col="datazones_excess",
        positive_label="Additional datazones among multi-datazone clusters",
    ),
    CountModelSpec(
        name="geographic_dispersion_size_adjusted",
        label="Geographic dispersion, size-adjusted",
        raw_outcome="cluster_n_datazones",
        binary_col="datazones_gt1",
        positive_col="datazones_excess",
        positive_label="Additional datazones among multi-datazone clusters",
        include_size=True,
    ),
]


# ---------------------------------------------------------------------------
# SIMD subdomain configuration (supplementary domain analyses)
# ---------------------------------------------------------------------------

DOMAINS: dict[str, dict[str, str]] = {
    "overall": {
        "label": "Overall",
        "rank_col": "dz_simd_rank",
        "quintile_col": "dz_simd_quintile",
    },
    "income":     {"label": "Income",     "rank_col": "dz_simd_income_rank"},
    "employment": {"label": "Employment", "rank_col": "dz_simd_employment_rank"},
    "education":  {"label": "Education",  "rank_col": "dz_simd_education_rank"},
    "health":     {"label": "Health",     "rank_col": "dz_simd_health_rank"},
    "access":     {"label": "Access",     "rank_col": "dz_simd_access_rank"},
    "crime":      {"label": "Crime",      "rank_col": "dz_simd_crime_rank"},
    "housing":    {"label": "Housing",    "rank_col": "dz_simd_housing_rank"},
}

SHARED_COUNT_TERMS: list[str] = [
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
]

SHARED_MIXING_TERMS: list[str] = SHARED_COUNT_TERMS + ["log_cluster_size_z"]

DEMOGRAPHIC_MIXING: dict[str, dict[str, str]] = {
    "age":     {"column": "age_band",        "label": "Age-band mixing",            "short_label": "Age"},
    "sex":     {"column": "sex",             "label": "Sex mixing",                 "short_label": "Sex"},
    "age_sex": {"column": "age_sex_profile", "label": "Joint age-sex profile mixing", "short_label": "Age-sex"},
}

DEMOGRAPHIC_MIXING_PREDICTOR_TERMS: Iterable[str] = [
    f"{prefix}_excess_mixing_z" for prefix in DEMOGRAPHIC_MIXING
]


# ---------------------------------------------------------------------------
# Wave configuration
# ---------------------------------------------------------------------------

WAVE_ORDER: Iterable[str] = [
    "B.1.177",
    "Alpha",
    "Delta",
    "BA.1",
    "BA.2",
    "BA.4",
    "BA.5",
    "BQ.1",
    "XBB",
]

WAVE_LABELS: dict[str, str] = {wave: wave for wave in WAVE_ORDER}
WAVE_LABELS["Other"] = "Other"


# ---------------------------------------------------------------------------
# Observed-vs-expected pair-probability matrices
# ---------------------------------------------------------------------------

MATRIX_VARIABLES: dict[str, dict[str, object]] = {
    "simd": {
        "column": "dz_simd_quintile",
        "label": "SIMD quintile",
        "levels": [1, 2, 3, 4, 5],
    },
    "age": {
        "column": "age_band",
        "label": "Age band",
        "levels": [
            "00-04", "05-09", "10-14", "15-19", "20-24", "25-29", "30-34",
            "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
            "70-74", "75+",
        ],
    },
}
