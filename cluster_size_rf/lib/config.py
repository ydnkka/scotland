from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = PROJECT_ROOT / "data" / "processed" / "scotland_clustering_analysis_dataset.parquet"
DATAZONE_SIMD_PATH = PROJECT_ROOT / "data" / "processed" / "scotland_datazone_simd_data.parquet"
OUT_DIR = PROJECT_ROOT / "cluster_size_rf" / "outputs"


RUN_MODE_PRESETS = {
    "quick": {
        "n_estimators": 250,
        "min_samples_leaf": 50,
        "max_train_rows": 120_000,
        "max_perm_rows": 30_000,
        "perm_repeats": 5,
    },
    "final": {
        "n_estimators": 1000,
        "min_samples_leaf": 20,
        "max_train_rows": None,
        "max_perm_rows": 75_000,
        "perm_repeats": 20,
    },
}


@dataclass(slots=True)
class AnalysisConfig:
    run_mode: str = "quick"
    primary_resolution: float = 0.3
    primary_large_min: int = 13
    random_state: int = 42
    window_dedup_strategy: str = "largest_cluster"
    simd_overall_mode: str = "quintile"
    simd_domain_mode: str = "quintile"
    simd_band_weighting: str = "population"
    use_sensitivity_context_controls: bool = False
    fit_singleton_binary: bool = True
    fit_simd_domain_decomposition: bool = True
    fit_secondary_multiclass: bool = True
    run_threshold_sensitivity: bool = False
    run_resolution_sensitivity: bool = False
    sensitivity_large_mins: tuple[int, ...] = (10, 20)
    sensitivity_resolutions: tuple[float, ...] = (0.1, 0.2, 0.4, 0.5)
    save_models: bool = True
    n_estimators: int = 250
    min_samples_leaf: int = 50
    max_train_rows: int | None = 120_000
    max_perm_rows: int | None = 30_000
    perm_repeats: int = 5

    @classmethod
    def from_run_mode(cls, run_mode: str = "quick", **overrides) -> "AnalysisConfig":
        if run_mode not in RUN_MODE_PRESETS:
            raise ValueError(f"Unknown run_mode {run_mode!r}; choose one of {sorted(RUN_MODE_PRESETS)}")
        values = {"run_mode": run_mode, **RUN_MODE_PRESETS[run_mode], **overrides}
        return cls(**values)

    def to_dict(self) -> dict:
        return asdict(self)


ID_COLUMNS = [
    "sequence_id",
    "patient_id",
    "cluster_id",
    "datazone",
    "window_idx",
    "window_id",
    "resolution",
]

DATE_COLUMNS = [
    "collection_date",
    "wn_mid_date",
]

OUTCOME_COLUMNS = [
    "cluster_size",
]

QC_COLUMNS = [
    "nextclade_qc",
]

FOCAL_SOCIODEMOGRAPHIC_COLUMNS = [
    "age_band",
    "age_midpoint",
    "sex",
    "is_female",
    "dz_simd_rank",
    "dz_simd_quintile",
    "dz_simd_decile",
    "dz_simd_income_rank",
    "dz_simd_employment_rank",
    "dz_simd_education_rank",
    "dz_simd_health_rank",
    "dz_simd_access_rank",
    "dz_simd_crime_rank",
    "dz_simd_housing_rank",
    "dz_urban_rural_class",
    "dz_population_density",
]

FOCAL_VACCINATION_COLUMNS = [
    "is_vaccinated",
    "vacc_dose_number",
    "days_since_vaccination",
    "vacc_booster",
    "vacc_product_name",
    "is_reinfection",
]

SURVEILLANCE_INCIDENCE_COLUMNS = [
    "dz_cum_prop_sequenced",
    "dz_cum_incidence_per_capita",
    "dz_7d_test_positivity",
    "wn_prop_sequenced",
    "wn_no_sequences",
]

CONTEXT_COLUMNS = [
    "pango_lineage",
    "who_voc",
    "dz_health_board",
    "test_reason",
    "test_type",
]

SENSITIVITY_CONTEXT_COLUMNS = [
    "hb_hospital_occupancy",
    "hb_reinfection_rate",
    "dz_cum_prop_vaccinated",
]

READ_COLUMNS = sorted(
    set(
        ID_COLUMNS
        + DATE_COLUMNS
        + OUTCOME_COLUMNS
        + QC_COLUMNS
        + FOCAL_SOCIODEMOGRAPHIC_COLUMNS
        + FOCAL_VACCINATION_COLUMNS
        + SURVEILLANCE_INCIDENCE_COLUMNS
        + CONTEXT_COLUMNS
        + SENSITIVITY_CONTEXT_COLUMNS
    )
)

SIMD_OVERALL_FEATURE_BY_MODE = {
    "rank": "dz_simd_rank",
    "quintile": "dz_simd_quintile",
    "decile": "dz_simd_decile",
}

SIMD_DOMAIN_FEATURES = [
    "dz_simd_income_rank",
    "dz_simd_employment_rank",
    "dz_simd_education_rank",
    "dz_simd_health_rank",
    "dz_simd_access_rank",
    "dz_simd_crime_rank",
    "dz_simd_housing_rank",
]

SIMD_DOMAIN_RANK_FEATURES = list(SIMD_DOMAIN_FEATURES)

SIMD_DOMAIN_BASES = [
    "income",
    "employment",
    "education",
    "health",
    "access",
    "crime",
    "housing",
]

SIMD_DOMAIN_FEATURES_BY_MODE = {
    "rank": [f"dz_simd_{base}_rank" for base in SIMD_DOMAIN_BASES],
    "quintile": [f"dz_simd_{base}_quintile" for base in SIMD_DOMAIN_BASES],
    "decile": [f"dz_simd_{base}_decile" for base in SIMD_DOMAIN_BASES],
}

SIMD_DOMAIN_POP_FEATURES_BY_MODE = {
    "rank": [f"dz_simd_{base}_rank" for base in SIMD_DOMAIN_BASES],
    "quintile": [f"dz_simd_{base}_pop_quintile" for base in SIMD_DOMAIN_BASES],
    "decile": [f"dz_simd_{base}_pop_decile" for base in SIMD_DOMAIN_BASES],
}

PRIMARY_CATEGORICAL_FEATURES = [
    "who_voc_model",
    "dz_health_board_model",
    "dz_urban_rural_class_model",
    "test_reason_model",
    "test_type_model",
    "vacc_product_name_model",
]

SENSITIVITY_NUMERIC_FEATURES = [
    "log_hb_hospital_occupancy",
    "hb_reinfection_rate",
    "dz_cum_prop_vaccinated",
]

TEST_REASON_MAP = {
    "symptomatic-citizen": "symptomatic_citizen",
    "I have coronavirus symptoms": "symptomatic_citizen",
    "I live~ work or study in a lockdown area with a coronavirus outbreak": "symptomatic_citizen",
    "symptomatic-essential-worker": "symptomatic_essential_worker",
    "Im an essential worker": "symptomatic_essential_worker",
    "scotland-wales-keyworker": "symptomatic_essential_worker",
    "wales-keyworker": "symptomatic_essential_worker",
    "test-for-contact-tracing": "contact_tracing",
    "test-for-contact-tracing-app": "contact_tracing",
    "test-for-contact-self-referral": "contact_tracing",
    "for-symptomatic-household-member": "contact_tracing",
    "Ive been in contact with a person who has tested positive for coronavirus and Ive been asked to take a test by a contact tracer (Northern Ireland and Scotland)": "contact_tracing",
    "Ive been in contact with a person who has tested positive for coronavirus and have since developed symptoms": "contact_tracing",
    "confirmatory-positive-test": "confirmatory",
    "confirmatory-other-reason": "confirmatory",
    "confirmatory-test-unclear": "confirmatory",
    "confirmatory-test-borders": "confirmatory",
    "told-to-order-repeat-test": "confirmatory",
    "self-isolation-support-grant": "isolation_scheme",
    "isolation-testing-home": "isolation_scheme",
    "isolation-testing-facility": "isolation_scheme",
    "gp-healthcare-request": "clinical",
    "antiviral-order": "clinical",
    "dental-patient-testing": "clinical",
    "I have been told to have a test before I go into hospital~ for example~ for surgery": "clinical",
    "zoe-symptom-study": "surveillance_research",
    "contact-testing-study": "surveillance_research",
    "events-research-programme": "surveillance_research",
    "serial-testing": "surveillance_research",
    "ntrg-member": "surveillance_research",
    "local-council-request": "local_outbreak",
    "attended-outbreak-venue": "local_outbreak",
    "community-testing": "local_outbreak",
    "scotland-university": "local_outbreak",
    "wales-university": "local_outbreak",
    "green-traveller": "travel",
    "other": "other",
    "Other": "other",
    "none": "other",
    "do-not-know": "other",
    "general-cta-referral": "other",
    "personal-assistant": "other",
    "Im a visiting professional": "other",
    "asymptomatic-home-order": "other",
}
