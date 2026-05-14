"""Shared constants for the Chapter 1 analysis.

Chapter 1 framing
-----------------
*Outcome* — cluster scale: ``cluster_size`` and ``cluster_n_datazones``,
modelled with a zero-truncated negative binomial (ZTNB) model.

*Exposure* — excess sociodemographic mixing: observed minus expected
pairwise discordance within lineage and analysis window, separately for
age band, sex, and SIMD area deprivation (and a joint profile as a
supplementary predictor).

*Question* — do clusters with more boundary-crossing mixing tend to be
larger or more geographically dispersed, after adjustment for lineage,
calendar time, surveillance intensity, and local incidence?
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


# ---------------------------------------------------------------------------
# Run defaults
# ---------------------------------------------------------------------------

QC_DEFAULT = "good"
PRIMARY_RESOLUTION = 0.3
LINEAGE_MIN_CLUSTERS = 30
CALENDAR_SPLINE_DF = 8


# ---------------------------------------------------------------------------
# Sequence-row columns pulled from the analysis dataset
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
    "dz_simd_decile",
    "dz_cum_incidence_per_capita",
    "dz_cum_prop_sequenced",
    "dz_7d_test_positivity",
    "wn_no_sequences",
    "dz_health_board_code",
]

# Subdomain rank columns used by the domain analysis.
DOMAINS: dict[str, dict[str, str]] = {
    "overall":    {"label": "Overall SIMD",  "rank_col": "dz_simd_rank"},
    "income":     {"label": "Income",        "rank_col": "dz_simd_income_rank"},
    "employment": {"label": "Employment",    "rank_col": "dz_simd_employment_rank"},
    "education":  {"label": "Education",     "rank_col": "dz_simd_education_rank"},
    "health":     {"label": "Health",        "rank_col": "dz_simd_health_rank"},
    "access":     {"label": "Access",        "rank_col": "dz_simd_access_rank"},
    "crime":      {"label": "Crime",         "rank_col": "dz_simd_crime_rank"},
    "housing":    {"label": "Housing",       "rank_col": "dz_simd_housing_rank"},
}


# ---------------------------------------------------------------------------
# Mixing predictors
# ---------------------------------------------------------------------------

MIXING_VARIABLES: dict[str, dict[str, str]] = {
    "age":                       {"column": "age_band",                   "short_label": "Age"},
    "sex":                       {"column": "sex",                        "short_label": "Sex"},
    "simd":                      {"column": "dz_simd_quintile",           "short_label": "SIMD"},
    "simd_decile":               {"column": "dz_simd_decile",             "short_label": "SIMD decile"},
    "demographic_profile":       {"column": "demographic_profile",        "short_label": "Demographic"},
    "socio_demographic_profile": {"column": "socio_demographic_profile",  "short_label": "Sociodemographic"},
}

MATRIX_VARIABLES: dict[str, dict[str, object]] = {
    "simd_quintile": {
        "column": "dz_simd_quintile",
        "label": "SIMD quintile",
        "levels": [1, 2, 3, 4, 5],
    },
    "simd_decile": {
        "column": "dz_simd_decile",
        "label": "SIMD decile",
        "levels": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
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

# Core predictors entered together in the main model.  SIMD-decile mixing
# is fit only as a sensitivity (see :func:`fit_models.fit_simd_decile_sensitivity`).
CORE_MIXING_PREDICTORS: Iterable[str] = ("age", "sex", "simd")

# Joint-profile predictors fit as separate single-predictor sensitivities:
#   * ``demographic_profile``       = age × sex
#   * ``socio_demographic_profile`` = age × sex × SIMD quintile
PROFILE_PREDICTORS: Iterable[str] = (
    "demographic_profile",
    "socio_demographic_profile",
)

EXCESS_MIXING_TERMS: Iterable[str] = [
    f"{prefix}_excess_mixing_z" for prefix in CORE_MIXING_PREDICTORS
]


# ---------------------------------------------------------------------------
# Adjustment covariates (z-scored)
# ---------------------------------------------------------------------------

ADJUSTMENT_TERMS: Iterable[str] = (
    "deprivation_z",
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
)


# ---------------------------------------------------------------------------
# Labels for output tables
# ---------------------------------------------------------------------------

TERM_LABELS: dict[str, str] = {
    "deprivation_z":          "Mean SIMD deprivation",
    "local_incidence_z":      "Local cumulative incidence",
    "local_seq_fraction_z":   "Local sequencing fraction",
    "window_seq_fraction_z":  "Window sequencing proportion",
    "test_positivity_z":      "Local test positivity",
    "log_cluster_size_z":     "log(cluster size)",
}
TERM_LABELS.update(
    {
        f"{prefix}_excess_mixing_z": f"{spec['short_label']} excess mixing"
        for prefix, spec in MIXING_VARIABLES.items()
    }
)

# Friendly per-table label for the SIMD-decile sensitivity (overrides the
# default ``SIMD decile excess mixing`` so it reads consistently in plots).
TERM_LABELS["simd_decile_excess_mixing_z"] = "SIMD-decile excess mixing"


# ---------------------------------------------------------------------------
# Outcome / model specifications
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class OutcomeSpec:
    """A single ZTNB outcome.

    Chapter 1 fits ZTNB only — both cluster size and geographic spread are
    modelled on the non-singleton positive sub-population (``raw > 1``).
    ``positive_col`` is the excess count (``raw - 1``).  ``hurdle_col`` is
    kept for descriptive use (e.g. dataset summaries) but no binary
    component is fit.
    """

    name: str
    label: str
    raw_col: str
    hurdle_col: str
    positive_col: str


OUTCOMES: Iterable[OutcomeSpec] = [
    OutcomeSpec(
        name="cluster_size",
        label="Cluster size",
        raw_col="cluster_size",
        hurdle_col="cluster_size_gt1",
        positive_col="cluster_size_excess",
    ),
    OutcomeSpec(
        name="geographic_spread",
        label="Geographic spread",
        raw_col="cluster_n_datazones",
        hurdle_col="datazones_gt1",
        positive_col="geographic_spread",
    ),
]


# ---------------------------------------------------------------------------
# Epidemic waves (used for the wave-interaction model + wave stratification)
# ---------------------------------------------------------------------------

WAVE_ORDER: Iterable[str] = (
    "B.1.177", "Alpha", "Delta",
    "BA.1", "BA.2", "BA.4", "BA.5", "BQ.1", "XBB",
)
WAVE_LABELS: dict[str, str] = {wave: wave for wave in WAVE_ORDER}
WAVE_LABELS["Other"] = "Other"

# Reference wave for the interaction model.  Delta is chosen because it is the
# largest single wave in the dataset and so the reference effect is the most
# precisely estimated.
WAVE_REFERENCE: str = "Delta"
