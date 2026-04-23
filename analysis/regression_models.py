"""Model registry and fitting utilities for the Scotland clustering analysis.

Public API
----------
get_model_registry()
    Return a DataFrame describing every registered model.

run_models(cluster_df, individual_df, ...)
    Fit a subset (or all) models and return posterior traces.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import polars as pl

from utils import bambi


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _standardise(values):
    """Negate z-score so higher values = greater deprivation."""
    return -(values - values.mean()) / values.std()


def _add_domain_zscores(
    cluster_df: pl.DataFrame,
    individual_df: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Append standardised domain columns to both DataFrames in-place copies."""
    domains = ["income", "employment", "education", "health", "access", "crime", "housing"]
    for dom in domains:
        cluster_df    = cluster_df.with_columns(
            _standardise(cluster_df[f"simd_{dom}_mean"]).alias(f"{dom}_zscore")
        )
        individual_df = individual_df.with_columns(
            _standardise(individual_df[f"dz_simd_{dom}_rank"]).alias(f"{dom}_zscore")
        )
    return cluster_df, individual_df


# ---------------------------------------------------------------------------
# Model spec dataclass
# ---------------------------------------------------------------------------

Level = Literal["cluster", "individual"]

@dataclass(frozen=True)
class ModelSpec:
    model:       str
    level:       Level
    description: str
    fixed_effects:       list[str] = field(default_factory=list)
    interaction_effects: list[str] = field(default_factory=list)

    def fit_kwargs(self) -> dict:
        kw: dict = {"fixed_effects": self.fixed_effects}
        if self.interaction_effects:
            kw["interaction_effects"] = self.interaction_effects
        return kw


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_DOMAINS = ["income", "employment", "education", "health", "access", "crime", "housing"]

_REGISTRY: list[ModelSpec] = [

    # ══════════════════════════════════════════════════════════════════════
    # CLUSTER-LEVEL
    # response   : n_sequences - 1  (secondary cases, singletons anchored at 0)
    # offset     : log_seq_prop
    # random     : (1 | window_id)
    # likelihood : negative-binomial
    # ══════════════════════════════════════════════════════════════════════

    # ── C1. Baseline deprivation ───────────────────────────────────────────
    # Establishes the marginal SIMD gradient before any controls.
    # No epoch: intentional — epoch effects absorbed into residual here
    # to keep the estimand clean for comparison with C2.
    ModelSpec(
        model="deprivation_main",
        level="cluster",
        description="Does area deprivation predict transmission cluster size?",
        fixed_effects=["C(simd_quintile_mode, Treatment(3))"],
    ),

    # ── C2. Deprivation × epoch ────────────────────────────────────────────
    # Core moderation question: did the SIMD gradient on cluster size
    # amplify or attenuate across VOC epochs?
    # Both main effects required by marginality before adding interaction.
    ModelSpec(
        model="deprivation_x_epoch",
        level="cluster",
        description=(
            "Does the deprivation effect on cluster size vary across VOC epochs? "
            "(did Alpha / Delta / Omicron amplify or attenuate the SIMD gradient?)"
        ),
        fixed_effects=[
            "C(simd_quintile_mode, Treatment(3))",
            "C(epoch, Treatment('Delta'))",
        ],
        interaction_effects=["C(simd_quintile_mode, Treatment(3)):C(epoch, Treatment('Delta'))"],
    ),

    # ── C3. Within-cluster socioeconomic mixing ───────────────────────────
    # simd_quintile_std captures deprivation heterogeneity within a cluster.
    # Epoch included because mixing patterns differ structurally across variants
    # (Omicron spread more broadly across socioeconomic groups than Alpha/Delta).
    # Omitting epoch would confound epoch-level variation in mixing with the
    # within-cluster mixing coefficient.
    ModelSpec(
        model="within_cluster_mixing",
        level="cluster",
        description=(
            "Does socioeconomic mixing within clusters predict size, "
            "adjusting for modal deprivation and epoch?"
        ),
        fixed_effects=[
            "C(simd_quintile_mode, Treatment(3))",
            "simd_quintile_std",
            "C(epoch, Treatment('Delta'))",
        ],
    ),

    # ── C4. Vaccination coverage × deprivation ────────────────────────────
    # frac_vaccinated is cluster-level aggregate coverage.
    # Interaction tests whether vaccination coverage moderates the SIMD
    # gradient: did higher vaccination in a cluster reduce the deprivation
    # penalty on cluster size, and did this differ by quintile?
    # Epoch controls for the confound that vaccination uptake is time-varying.
    ModelSpec(
        model="vaccination_x_deprivation",
        level="cluster",
        description=(
            "Does cluster-level vaccination coverage moderate the deprivation "
            "effect on cluster size, adjusting for epoch?"
        ),
        fixed_effects=[
            "C(simd_quintile_mode, Treatment(3))",
            "frac_vaccinated",
            "C(epoch, Treatment('Delta'))",
        ],
        interaction_effects=["C(simd_quintile_mode, Treatment(3)):frac_vaccinated"],
    ),

    # ── C5. SIMD domain decomposition ─────────────────────────────────────
    # One model per domain; each replaces SIMD quintile with a standardised
    # domain z-score (negated so higher = more deprived).
    # Age and sex structure included to isolate the domain signal from
    # compositional differences across clusters.
    # Epoch deliberately omitted: models are focused decompositions run
    # in isolation for comparability across domains. Epoch confounding is
    # minimal because domain ranks are stable over the study period.
    *[
        ModelSpec(
            model=f"domain_{dom}",
            level="cluster",
            description=(
                f"Does the {dom} deprivation domain drive cluster size, "
                f"adjusting for age and sex structure?"
            ),
            fixed_effects=[
                f"{dom}_zscore",
                "median_age",
                "age_diversity",
                "frac_female",
            ],
        )
        for dom in _DOMAINS
    ],

    # ── C6. Age and sex structure ──────────────────────────────────────────
    # Tests whether demographic composition of clusters predicts size
    # independently of area deprivation.
    # Epoch included because age structure of infected cohorts shifted
    # substantially across variants (e.g. younger Delta wave).
    ModelSpec(
        model="age_sex",
        level="cluster",
        description=(
            "Do age and sex structure predict cluster size, "
            "independent of area deprivation and epoch?"
        ),
        fixed_effects=[
            "C(simd_quintile_mode, Treatment(3))",
            "C(epoch, Treatment('Delta'))",
            "median_age",
            "age_diversity",
            "frac_female",
        ],
    ),

    # ── C7. Fully adjusted ────────────────────────────────────────────────
    # Combines all predictors. SIMD:epoch interaction retained because the
    # core narrative includes epoch moderation; without it the fully adjusted
    # model would understate the main finding.
    # Reports independent contributions of each predictor block.
    ModelSpec(
        model="fully_adjusted",
        level="cluster",
        description=(
            "Independent contributions of deprivation, epoch, age, sex, "
            "and within-cluster mixing, with deprivation–epoch moderation."
        ),
        fixed_effects=[
            "C(simd_quintile_mode, Treatment(3))",
            "C(epoch, Treatment('Delta'))",
            "median_age",
            "age_diversity",
            "frac_female",
            "simd_quintile_std",
        ],
        interaction_effects=["C(simd_quintile_mode, Treatment(3)):C(epoch, Treatment('Delta'))"],
    ),


    # ══════════════════════════════════════════════════════════════════════
    # INDIVIDUAL-LEVEL
    # response   : non_singleton_k / non_singleton_n  (binomial)
    # covariate  : log_seq_prop  (adjusts for surveillance density)
    # likelihood : binomial (logit link)
    # reference  : SIMD quintile 3 (middle deprivation)
    # ══════════════════════════════════════════════════════════════════════

    # ── I1. Baseline deprivation ───────────────────────────────────────────
    # Marginal SIMD gradient before any controls.
    # log_seq_prop is always included at individual level: higher local
    # sequencing effort mechanically inflates detected transmission links.
    ModelSpec(
        model="deprivation_main",
        level="individual",
        description=(
            "Does area deprivation predict individual transmission chain membership?"
        ),
        fixed_effects=["C(dz_simd_quintile, Treatment(3))", "log_seq_prop"],
    ),

    # ── I2. Deprivation × epoch ────────────────────────────────────────────
    # Did the SIMD gradient on clustering consistency vary across VOC epochs?
    # Both main effects required by marginality.
    ModelSpec(
        model="deprivation_x_epoch",
        level="individual",
        description=(
            "Does the deprivation gradient in cluster membership vary "
            "across VOC epochs?"
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(epoch, Treatment('Delta'))",
            "log_seq_prop",
        ],
        interaction_effects=["C(dz_simd_quintile, Treatment(3)):C(epoch, Treatment('Delta'))"],
    ),

    # ── I3. Deprivation adjusted (confounding check) ───────────────────────
    # Additive model — no interaction. Estimand: does the SIMD gradient
    # survive adjustment for the main confounders?
    # Distinct from I4: this is a confounding question, not a moderation question.
    # Including the interaction here would prevent a clean confounding answer.
    ModelSpec(
        model="deprivation_adjusted",
        level="individual",
        description=(
            "Does the deprivation gradient persist after adjusting "
            "for epoch, age, sex, and vaccination status?"
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(epoch, Treatment('Delta'))",
            "C(age_band, Treatment('40-44'))",
            "is_female",
            "is_vaccinated",
            "log_seq_prop",
        ],
    ),

    # ── I4. Deprivation × epoch, adjusted ─────────────────────────────────
    # Does the epoch moderation of the SIMD gradient survive demographic
    # adjustment? Extends I2 by adding age, sex, vaccination.
    # Separating I3 and I4 keeps confounding and moderation as distinct
    # estimands with clean interpretations.
    ModelSpec(
        model="deprivation_x_epoch_adjusted",
        level="individual",
        description=(
            "Does the deprivation–epoch interaction in cluster membership "
            "persist after adjusting for age, sex, and vaccination?"
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(epoch, Treatment('Delta'))",
            "C(age_band, Treatment('40-44'))",
            "is_female",
            "is_vaccinated",
            "log_seq_prop",
        ],
        interaction_effects=["C(dz_simd_quintile, Treatment(3)):C(epoch, Treatment('Delta'))"],
    ),

    # ── I5. Vaccination × deprivation ─────────────────────────────────────
    # Did vaccination close the deprivation gap in cluster membership?
    # No epoch: isolates the vaccine equity question without epoch confounding
    # the interaction cells. Compare with I6 to see if the pattern is
    # epoch-varying.
    ModelSpec(
        model="vaccination_x_deprivation",
        level="individual",
        description=(
            "Did vaccination close the deprivation gap in cluster membership?"
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "is_vaccinated",
            "log_seq_prop",
        ],
        interaction_effects=["C(dz_simd_quintile, Treatment(3)):is_vaccinated"],
    ),

    # ── I6. Vaccination × epoch ────────────────────────────────────────────
    # Did vaccination reduce cluster membership more in later epochs, when
    # population coverage was higher and immunity more established?
    # Interaction is is_vaccinated:C(epoch, Treatment('Delta')) — NOT SIMD:epoch.
    # SIMD included as main effect only (no SIMD:epoch here) to keep
    # this focused on the vaccination–time question.
    ModelSpec(
        model="vaccination_x_epoch",
        level="individual",
        description=(
            "Did vaccination reduce cluster membership more in later epochs, "
            "when population coverage was higher?"
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(epoch, Treatment('Delta'))",
            "is_vaccinated",
            "log_seq_prop",
        ],
        interaction_effects=["is_vaccinated:C(epoch, Treatment('Delta'))"],
    ),

    # ── I7. SIMD domain decomposition ─────────────────────────────────────
    # One model per domain; each replaces SIMD quintile with the domain
    # z-score. Age, sex, vaccination included to isolate the domain signal.
    # Epoch omitted deliberately for the same reason as cluster-level C5:
    # focused decomposition, cross-domain comparability.
    *[
        ModelSpec(
            model=f"domain_{dom}",
            level="individual",
            description=(
                f"Does the {dom} deprivation domain drive individual cluster "
                f"membership, adjusting for age, sex, and vaccination?"
            ),
            fixed_effects=[
                f"{dom}_zscore",
                "C(age_band, Treatment('40-44'))",
                "is_female",
                "is_vaccinated",
                "log_seq_prop",
            ],
        )
        for dom in _DOMAINS
    ],

    # ── I8. Age × deprivation ──────────────────────────────────────────────
    # Are older individuals in deprived areas disproportionately linked
    # into transmission clusters? Tests a specific mechanistic hypothesis
    # about age-deprivation synergy in transmission risk.
    ModelSpec(
        model="age_x_deprivation",
        level="individual",
        description=(
            "Are older people in deprived areas disproportionately "
            "linked into transmission clusters?"
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(age_band, Treatment('40-44'))",
            "is_female",
            "is_vaccinated",
            "log_seq_prop",
        ],
        interaction_effects=["C(dz_simd_quintile, Treatment(3)):C(age_band, Treatment('40-44'))"],
    ),

    # ── I9. Fully adjusted ────────────────────────────────────────────────
    # All predictors combined. SIMD:epoch interaction retained as the primary
    # moderation term. Reports independent contributions of each predictor
    # block after mutual adjustment.
    ModelSpec(
        model="fully_adjusted",
        level="individual",
        description=(
            "Independent contributions of all predictors, "
            "with deprivation–epoch moderation."
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(epoch, Treatment('Delta'))",
            "C(age_band, Treatment('40-44'))",
            "is_female",
            "is_vaccinated",
            "log_seq_prop",
        ],
        interaction_effects=["C(dz_simd_quintile, Treatment(3)):C(epoch, Treatment('Delta'))"],
    ),
]


# ---------------------------------------------------------------------------
# Public: registry inspection
# ---------------------------------------------------------------------------

def get_model_registry() -> pl.DataFrame:
    """Return a tidy DataFrame describing every registered model.

    Returns
    -------
    pl.DataFrame with columns:
        model        - model name (str)
        level        - "cluster" or "individual"
        description  - plain-English research question (str)
        fixed_effects       - list[str]
        interaction_effects - list[str]
    """
    return pl.DataFrame(
        {
            "model":               [s.model       for s in _REGISTRY],
            "level":               [s.level       for s in _REGISTRY],
            "description":         [s.description for s in _REGISTRY],
            "fixed_effects":       [s.fixed_effects       for s in _REGISTRY],
            "interaction_effects": [s.interaction_effects for s in _REGISTRY],
        }
    )


# ---------------------------------------------------------------------------
# Public: fitting
# ---------------------------------------------------------------------------

def run_models(
    cluster_df:    pl.DataFrame,
    individual_df: pl.DataFrame,
    *,
    models: list[str] | None = None,
    levels: list[Level] | None = None,
    run_all: bool = False,
    save_dir: str = "bambi_outputs",
) -> dict[str, dict[str, object]]:
    """Fit a selection of registered models and return posterior traces.

    Parameters
    ----------
    cluster_df, individual_df
        Input DataFrames as returned by the data-loading helpers.
    models
        Model names to run, e.g. ``["deprivation_main", "deprivation_adjusted"]``.
        Matched against ``ModelSpec.model``; the same name may appear at both
        levels — use ``levels`` to restrict if needed.
        Ignored when ``run_all=True``.
    levels
        Restrict fitting to ``["cluster"]``, ``["individual"]``, or both
        (default). Applied after ``models`` filtering.
    run_all
        If ``True``, fit every registered model regardless of ``models``.
    save_dir
        Directory passed through to the Bambi fitting helpers.

    Returns
    -------
    dict with keys ``"cluster"`` and ``"individual"``, each mapping
    ``model_name → posterior_trace``.

    Examples
    --------
    # Fit just the three talk models (individual level only)
    run_models(
        cluster_df, individual_df,
        models=["deprivation_main", "deprivation_adjusted", "vaccination_x_deprivation"],
        levels=["individual"],
    )

    # Fit every registered model
    run_models(cluster_df, individual_df, run_all=True)
    """
    if not run_all and models is None:
        raise ValueError("Provide a list of model names via `models` or set `run_all=True`.")

    # ── Resolve which specs to fit ─────────────────────────────────────────
    specs = _REGISTRY if run_all else [s for s in _REGISTRY if s.model in models]

    if levels is not None:
        specs = [s for s in specs if s.level in levels]

    if not specs:
        raise ValueError(
            f"No registered models matched models={models!r}, levels={levels!r}. "
            f"Call get_model_registry() to see available options."
        )

    # ── Pre-processing ─────────────────────────────────────────────────────
    cluster_df    = cluster_df.clone()
    individual_df = individual_df.clone()

    cluster_df    = cluster_df.with_columns(
        (pl.col("n_sequences") - 1).alias("n_sequences")
    )
    cluster_df, individual_df = _add_domain_zscores(cluster_df, individual_df)

    # ── Fit ────────────────────────────────────────────────────────────────
    traces: dict[str, dict[str, object]] = {"cluster": {}, "individual": {}}

    cluster_specs    = [s for s in specs if s.level == "cluster"]
    individual_specs = [s for s in specs if s.level == "individual"]

    for spec in cluster_specs:
        traces["cluster"][spec.model] = bambi.fit_cluster_model(
            cluster_df, run_id=spec.model, save_dir=save_dir, **spec.fit_kwargs()
        )

    for spec in individual_specs:
        traces["individual"][spec.model] = bambi.fit_individual_model(
            individual_df, run_id=spec.model, save_dir=save_dir, **spec.fit_kwargs()
        )

    return traces