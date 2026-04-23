from utils import bambi

def _standardise(values):
    return (values - values.mean()) / values.std()


# ---------------------------------------------------------------------------
# Cluster-level formulas  (fit_cluster_model)
# response = n_sequences, offset = log_seq_prop, random = (1 | window_id)
# ---------------------------------------------------------------------------

CLUSTER_FORMULAS = {

    # Q1 — Does deprivation predict transmission cluster size?
    # Baseline: SIMD quintile only, reference = quintile 3 (middle)
    "deprivation_main": dict(
        fixed_effects=["C(simd_quintile_mode, Treatment(3))"],
    ),

    # Q2 — Does the deprivation effect vary across VOC epochs?
    # Key question: did Alpha/Delta/Omicron amplify or attenuate inequality?
    "deprivation_x_epoch": dict(
        fixed_effects=["C(simd_quintile_mode, Treatment(3))", "C(epoch)"],
        interaction_effects=["C(simd_quintile_mode, Treatment(3)):C(epoch)"],
    ),

    # Q3 — Is within-cluster socioeconomic mixing associated with cluster size?
    # simd_quintile_std captures heterogeneity of deprivation within a cluster.
    # High std = mixing across deprivation groups; low = homogeneous clusters.
    "within_cluster_mixing": dict(
        fixed_effects=[
            "C(simd_quintile_mode, Treatment(3))",
            "simd_quintile_std",
        ],
    ),

    # Q4 — Does vaccination coverage moderate cluster size?
    # Controls for deprivation since vaccination uptake is deprivation-patterned.
    "vaccination_deprivation": dict(
        fixed_effects=[
            "C(simd_quintile_mode, Treatment(3))",
            "frac_vaccinated",
            "C(epoch)",
        ],
    ),

    # Q5 — Multidimensional deprivation: is the overall SIMD rank effect
    # driven by a specific domain (income, employment, housing…)?
    # Fit separately per domain to decompose the deprivation signal.
    # Generate programmatically:
    **{
        f"domain_{dom}": dict(
            fixed_effects=[
                f"{dom}_zscore",
                "C(age_band)",
                "is_female",
                "is_vaccinated",
            ],
        )
        for dom in [
            "income", "employment", "education",
            "health", "access", "crime", "housing",
        ]
    },

    # Q6 — Age and sex structure as transmission drivers
    # median_age and age_diversity together capture both the central tendency
    # and spread of ages in the cluster.
    "age_sex": dict(
        fixed_effects=[
            "median_age",
            "age_diversity",
            "frac_female",
            "C(simd_quintile_mode, Treatment(3))",
            "C(epoch)",
        ],
    ),

    # Q7 — Full adjusted model
    # Combines deprivation, vaccination, age/sex, and epoch.
    # Useful for decomposing independent contributions.
    "fully_adjusted": dict(
        fixed_effects=[
            "C(simd_quintile_mode, Treatment(3))",
            "C(epoch)",
            "frac_vaccinated",
            "median_age",
            "age_diversity",
            "frac_female",
            "simd_quintile_std",
        ],
    ),
}


# ---------------------------------------------------------------------------
# Individual-level formulas  (fit_individual_model)
# response = non_singleton_k / non_singleton_n
# ---------------------------------------------------------------------------

INDIVIDUAL_FORMULAS = {

    # Q1 — Does deprivation predict transmission chain membership?
    # non_singleton_k/n = fraction of windows where patient was in a cluster >1
    # i.e. evidence of onward transmission vs isolated case
    "deprivation_main": dict(
        fixed_effects=["C(dz_simd_quintile, Treatment(3))"],
    ),

    # Q2 — Does the deprivation effect hold after adjusting for age, sex,
    # and vaccination? Core confounders for individual-level analyses.
    "deprivation_adjusted": dict(
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(age_band)",
            "is_female",
            "is_vaccinated",
        ],
    ),

    # Q3 — Does vaccination protect against cluster membership,
    # and does this vary by deprivation?
    # Addresses vaccine equity — did vaccination close the deprivation gap?
    "vaccination_x_deprivation": dict(
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "is_vaccinated",
        ],
        interaction_effects=["C(dz_simd_quintile, Treatment(3)):is_vaccinated"],
    ),

    # Q4 — Multidimensional deprivation: is the overall SIMD rank effect
    # driven by a specific domain (income, employment, housing…)?
    # Fit separately per domain to decompose the deprivation signal.
    # Generate programmatically:
    **{
        f"domain_{dom}": dict(
            fixed_effects=[
                f"{dom}_zscore",
                "C(age_band)",
                "is_female",
                "is_vaccinated",
            ],
        )
        for dom in [
            "income", "employment", "education",
            "health", "access", "crime", "housing",
        ]
    },

    # Q5 — Age interaction with deprivation
    # Are older people in deprived areas disproportionately likely to appear
    # in transmission clusters?
    "age_x_deprivation": dict(
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(age_band)",
            "is_female",
            "is_vaccinated",
        ],
        interaction_effects=["C(dz_simd_quintile, Treatment(3)):C(age_band)"],
    ),

    # Q6 — Full adjusted model
    "fully_adjusted": dict(
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(age_band)",
            "is_female",
            "is_vaccinated",
        ],
    ),
}

def run_all(cluster_df, individual_df, save_dir="bambi_outputs"):

    cluster_df["overall_zscore"] = _standardise(cluster_df["simd_overall_mean"])
    cluster_df["income_zscore"] = _standardise(cluster_df["simd_income_mean"])
    cluster_df["employment_zscore"] = _standardise(cluster_df["simd_employment_mean"])
    cluster_df["education_zscore"] = _standardise(cluster_df["simd_education_mean"])
    cluster_df["health_zscore"] = _standardise(cluster_df["simd_health_mean"])
    cluster_df["access_zscore"] = _standardise(cluster_df["simd_access_mean"])
    cluster_df["crime_zscore"] = _standardise(cluster_df["simd_crime_mean"])
    cluster_df["housing_zscore"] = _standardise(cluster_df["simd_housing_mean"])

    individual_df["income_zscore"] = _standardise(individual_df["dz_simd_income_rank"])
    individual_df["employment_zscore"] = _standardise(individual_df["dz_simd_employment_rank"])
    individual_df["education_zscore"] = _standardise(individual_df["dz_simd_education_rank"])
    individual_df["health_zscore"] = _standardise(individual_df["dz_simd_health_rank"])
    individual_df["access_zscore"] = _standardise(individual_df["dz_simd_access_rank"])
    individual_df["crime_zscore"] = _standardise(individual_df["dz_simd_crime_rank"])
    individual_df["housing_zscore"] = _standardise(individual_df["dz_simd_housing_rank"])
    cluster_traces, individual_traces = {}, {}

    for run_id, kwargs in CLUSTER_FORMULAS.items():
        cluster_traces[run_id] = bambi.fit_cluster_model(
            cluster_df, run_id=run_id, save_dir=save_dir, **kwargs
        )

    for run_id, kwargs in INDIVIDUAL_FORMULAS.items():
        individual_traces[run_id] = bambi.fit_individual_model(
            individual_df, run_id=run_id, save_dir=save_dir, **kwargs
        )

    return cluster_traces, individual_traces