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
    random_effects:  list[str] = field(default_factory=list)

    def fit_kwargs(self) -> dict:
        kw: dict = {"fixed_effects": self.fixed_effects}
        if self.interaction_effects:
            kw["interaction_effects"] = self.interaction_effects
        if self.random_effects:
            kw["random_effects"] = self.random_effects
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
    # random     : (1 | window_id) or (1 | pango_lineage)
    # likelihood : negative-binomial
    # ══════════════════════════════════════════════════════════════════════

    ModelSpec(
        model="deprivation_main",
        level="cluster",
        description="Does area deprivation predict transmission cluster size?",
        fixed_effects=["C(simd_quintile_mode, Treatment(3))"],
    ),

    # ══════════════════════════════════════════════════════════════════════
    # INDIVIDUAL-LEVEL
    # response   : non_singleton_k / non_singleton_n  (binomial)
    # covariate  : log_seq_prop  (adjusts for surveillance density)
    # random     : (1 | window_id) or (1 | pango_lineage)
    # likelihood : binomial (logit link)
    # reference  : SIMD quintile 3 (middle deprivation)
    # reference  : AY.4 pango lineage (Delta most common in dataset)
    # ══════════════════════════════════════════════════════════════════════

    ModelSpec(
        model="deprivation_main",
        level="individual",
        description=(
            "Does area deprivation predict individual transmission chain membership?"
        ),
        fixed_effects=["C(dz_simd_quintile, Treatment(3))", "log_seq_prop"],
    ),

    ModelSpec(
        model="deprivation_main_wn",
        level="individual",
        description=(
            "Does area deprivation predict individual transmission chain membership, "
            "with a window-level random intercept to account for clustering?"
        ),
        fixed_effects=["C(dz_simd_quintile, Treatment(3))", "log_seq_prop"],
        random_effects=["(1 | window_id)"]
    ),

    ModelSpec(
        model="deprivation_demographic",
        level="individual",
        description=(
            "Does the deprivation gradient in cluster membership persist after adjusting for "
            "age, sex, and vaccination status?"
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))", "C(age_band, Treatment('40-44'))",
            "is_female", "is_vaccinated", "log_seq_prop"],
        random_effects=["(1 | window_id)"]
    ),

    ModelSpec(
        model="deprivation_x_demographic",
        level="individual",
        description=(
            "Does area deprivation predict individual transmission chain membership, "
            "with a window-level random intercept to account for clustering?"
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))", "C(age_band, Treatment('40-44'))",
            "is_female", "is_vaccinated", "log_seq_prop"],
        random_effects=["(1 | window_id)"],
        interaction_effects=["C(dz_simd_quintile, Treatment(3)):C(age_band, Treatment('40-44'))"],
    ),

    ModelSpec(
        model="deprivation_x_pango_lineage",
        level="individual",
        description=(
            "Does the deprivation gradient in cluster membership vary "
            "across pango_lineages?"
        ),
        fixed_effects=[
            "C(dz_simd_quintile, Treatment(3))",
            "C(pango_lineage, Treatment('AY.4'))",
            "log_seq_prop",
        ],
        random_effects=["(1 | window_id)"],
        interaction_effects=["C(dz_simd_quintile, Treatment(3)):C(pango_lineage, Treatment('AY.4'))"],
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