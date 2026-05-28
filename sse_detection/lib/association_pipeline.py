"""Parameterized association-regression runner for SSE sensitivity notebooks."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from utils import CLADES, load_analysis_columns
from utils import PERIOD_INTENSITY

from .entropy import cluster_age_conditional_binary_entropy
from .io import load_sse_outputs
from .regression import (
    AssociationModel,
    bh_adjust,
    categorical_term,
    cluster_se_diagnostics,
    fit_binomial_glm,
    fit_conditional_logit,
    fit_firth_logit,
    make_formula,
    model_fit_stats,
    model_variables_from_terms,
    parameter_names_for_term,
    robust_wald_for_params,
    robust_wald_for_prefix,
    tidy_odds_ratios,
    tidy_single_parameter_wald,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]


__all__ = [
    "COMPOSITION_SPECS",
    "DEFAULT_MIXING_FEATURES",
    "OBSERVED_MIXING_FEATURES",
    "OBSERVED_MIXING_FEATURES_X10",
    "OBSERVED_MIXING_REFERENCE_X10",
    "POLICY_ERA_BY_PERIOD",
    "POLICY_ERA_ORDER",
    "VACCINATION_COMPOSITION_EXPANDED_ADJUSTERS",
    "VACCINATION_COMPOSITION_JOINT_GROUPS",
    "VACCINATION_COMPOSITION_SPECS",
    "VACCINATION_NODE_EXPANDED_ADJUSTERS",
    "VACCINATION_NODE_FEATURES",
    "VACCINATION_NODE_JOINT_GROUPS",
    "VACCINATION_MIXING_TERTILE_ORDER",
    "WINDOW_SURVEILLANCE_ADJUSTERS",
    "EXPANDED_CONTEXT_ADJUSTERS",
    "TABLE_DISPLAY_COLUMNS",
    "AssociationFrames",
    "add_clade_group",
    "add_policy_era",
    "add_standardised_adjusters",
    "add_vaccination_composition_features",
    "add_vaccination_mixing_features",
    "add_vaccination_node_features",
    "default_model_sets",
    "load_association_frames",
    "run_association_pipeline",
    "run_main_association_analysis",
    "run_policy_analysis",
    "run_vaccination_analysis",
    "select_table_columns",
]


WINDOW_SURVEILLANCE_ADJUSTERS = [
    "z_wn_prop_sequenced",
    "z_log1p_wn_positive_tests",
]

EXPANDED_CONTEXT_ADJUSTERS = [
    "z_dz_cum_prop_sequenced",
    "z_dz_cum_incidence_per_capita",
    "z_dz_7d_test_positivity",
    "z_log1p_dz_cum_positive_tests",
]

POLICY_ERA_BY_PERIOD = {
    "E0": "early_restriction_easing",
    "L1": "early_restriction_easing",
    "P1": "early_restriction_easing",
    "P2": "early_restriction_easing",
    "P3": "early_restriction_easing",
    "T1": "autumn_winter_restrictions",
    "F5": "autumn_winter_restrictions",
    "L2": "autumn_winter_restrictions",
    "SL": "spring_summer_2021_easing",
    "L3": "spring_summer_2021_easing",
    "L21": "spring_summer_2021_easing",
    "L0": "spring_summer_2021_easing",
    "NN": "near_normal_delta",
    "OM": "omicron_response",
    "FE": "omicron_response",
    "PR": "post_restriction",
}

POLICY_ERA_ORDER = [
    "early_restriction_easing",
    "autumn_winter_restrictions",
    "spring_summer_2021_easing",
    "near_normal_delta",
    "omicron_response",
    "post_restriction",
]

COMPOSITION_SPECS = [
    {
        "name": "sex",
        "column": "sex",
        "reference": "Male",
        "label": "Sex",
    },
    {
        "name": "age_band",
        "column": "age_band",
        "reference": "20-24",
        "fallback_references": ["25-29", "30-34"],
        "label": "Age band",
    },
    {
        "name": "simd_quintile",
        "column": "dz_simd_quintile",
        "reference": "1",
        "label": "SIMD quintile",
    },
    {
        "name": "urban_rural_class",
        "column": "dz_urban_rural_class",
        "reference": "Large Urban Areas",
        "label": "Urban/rural class",
    },
    {
        "name": "health_board",
        "column": "dz_health_board",
        "reference": "Greater Glasgow and Clyde",
        "label": "Health board",
    },
]

DEFAULT_MIXING_FEATURES = [
    "sex_entropy_z",
    "age_entropy_z",
    "simd_entropy_z",
    "urban_rural_entropy_z",
    "health_board_entropy_z",
]

OBSERVED_MIXING_FEATURES = [
    "sex_entropy_obs",
    "age_entropy_obs",
    "simd_entropy_obs",
    "urban_rural_entropy_obs",
    "health_board_entropy_obs",
]
OBSERVED_MIXING_FEATURES_X10 = [
    f"{feature}_x10" for feature in OBSERVED_MIXING_FEATURES
]
OBSERVED_MIXING_REFERENCE_X10 = "per 0.1 increase in observed normalised entropy"

VACCINATION_COMPOSITION_SPECS = [
    {
        "name": "is_vaccinated",
        "column": "is_vaccinated",
        "label": "Individual vaccinated",
        "reference": "vaccinated vs unvaccinated",
        "categorical": False,
    },
    {
        "name": "vacc_dose_cat",
        "column": "vacc_dose_cat",
        "reference": "0",
        "label": "Prior vaccination dose category",
        "categorical": True,
    },
    {
        "name": "vacc_booster_status",
        "column": "vacc_booster_status",
        "reference": "unvaccinated",
        "label": "Prior booster status",
        "categorical": True,
    },
    {
        "name": "days_since_vaccination_cat",
        "column": "days_since_vaccination_cat",
        "reference": "unvaccinated",
        "label": "Days since prior vaccination",
        "categorical": True,
    },
    {
        "name": "dz_cum_prop_vaccinated",
        "column": "z_dz_cum_prop_vaccinated",
        "label": "Datazone cumulative vaccination events per capita",
        "reference": "per 1 SD increase in vaccination-event coverage",
        "categorical": False,
    },
]

VACCINATION_NODE_FEATURES = [
    {
        "name": "node_prop_vaccinated",
        "column": "z_node_prop_vaccinated",
        "label": "Node proportion vaccinated",
        "reference": "per 1 SD increase",
    },
    {
        "name": "node_prop_boosted",
        "column": "z_node_prop_boosted",
        "label": "Node proportion boosted",
        "reference": "per 1 SD increase",
    },
    {
        "name": "node_mean_vacc_dose",
        "column": "z_node_mean_vacc_dose",
        "label": "Node mean prior dose count",
        "reference": "per 1 SD increase",
    },
    {
        "name": "node_median_days_since_vaccination",
        "column": "z_node_median_days_since_vaccination",
        "label": "Node median days since prior vaccination",
        "reference": "per 1 SD increase among nodes with vaccinated sequences",
    },
    {
        "name": "node_mean_dz_cum_prop_vaccinated",
        "column": "z_node_mean_dz_cum_prop_vaccinated",
        "label": "Node mean datazone vaccination-event coverage",
        "reference": "per 1 SD increase",
    },
]

VACCINATION_COMPOSITION_EXPANDED_ADJUSTERS = [
    "C(age_band)",
    "C(sex)",
    "C(dz_simd_quintile)",
    "C(dz_urban_rural_class)",
    "C(dz_health_board)",
]

VACCINATION_NODE_EXPANDED_ADJUSTERS = [
    *EXPANDED_CONTEXT_ADJUSTERS,
]

VACCINATION_COMPOSITION_JOINT_GROUPS = [
    {
        "name": "joint_status_area",
        "label": "Vaccination status plus datazone vaccination-event coverage",
        "specs": ["is_vaccinated", "dz_cum_prop_vaccinated"],
    },
    {
        "name": "joint_dose_area",
        "label": "Dose category plus datazone vaccination-event coverage",
        "specs": ["vacc_dose_cat", "dz_cum_prop_vaccinated"],
    },
    {
        "name": "joint_booster_area",
        "label": "Booster status plus datazone vaccination-event coverage",
        "specs": ["vacc_booster_status", "dz_cum_prop_vaccinated"],
    },
    {
        "name": "joint_recency_area",
        "label": "Vaccination recency plus datazone vaccination-event coverage",
        "specs": ["days_since_vaccination_cat", "dz_cum_prop_vaccinated"],
    },
]

VACCINATION_NODE_JOINT_GROUPS = [
    {
        "name": "joint_status_area",
        "label": "Node vaccination status plus datazone vaccination-event coverage",
        "specs": ["node_prop_vaccinated", "node_mean_dz_cum_prop_vaccinated"],
    },
]

VACCINATION_MIXING_TERTILE_ORDER = [
    "more_homogeneous",
    "as_expected",
    "more_mixed",
]

STANDARDISE_SPECS = {
    "z_wn_prop_sequenced": "wn_prop_sequenced",
    "z_log1p_wn_positive_tests": "log1p_wn_positive_tests",
    "z_dz_cum_prop_sequenced": "dz_cum_prop_sequenced",
    "z_dz_cum_incidence_per_capita": "dz_cum_incidence_per_capita",
    "z_dz_7d_test_positivity": "dz_7d_test_positivity",
    "z_log1p_dz_cum_positive_tests": "log1p_dz_cum_positive_tests",
    "z_dz_cum_prop_vaccinated": "dz_cum_prop_vaccinated",
}

TABLE_DISPLAY_COLUMNS = {
    "wald": [
        "domain",
        "model_set",
        "predictor_set",
        "joint_model",
        "predictor",
        "label",
        "reference",
        "term",
        "chi2",
        "df",
        "P>chi2",
        "p_adj_bh",
        "n_model_rows",
        "n_sequences",
        "n_nodes",
        "dropped_nonvarying_rows",
        "dropped_nonvarying_strata",
        "dropped_nonvarying_detail",
    ],
    "odds_ratios": [
        "domain",
        "model_set",
        "predictor_set",
        "joint_model",
        "predictor",
        "label",
        "reference",
        "term",
        "estimate",
        "std_error",
        "p_value",
        "odds_ratio",
        "or_low",
        "or_high",
    ],
    "fit_stats": [
        "domain",
        "model_set",
        "predictor_set",
        "joint_model",
        "predictor",
        "r2_mcfadden",
        "converged",
        "aic",
        "bic_llf",
        "log_likelihood",
        "ll_null",
        "n_model_rows",
        "n_sequences",
        "n_nodes",
    ],
}


@dataclass
class AssociationFrames:
    """Prepared node and sequence-level model frames."""

    node_stats: pd.DataFrame
    eligible_nodes: pd.DataFrame
    node_model_base: pd.DataFrame
    composition_base: pd.DataFrame
    cluster_diagnostics: pd.DataFrame
    min_candidate_size: int


def _dedupe(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))


def _slug(value: object) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_").lower()
    return slug or "all"


def _concat_or_empty(tables: list[pd.DataFrame]) -> pd.DataFrame:
    tables = [table for table in tables if table is not None and not table.empty]
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def _has_categorical_adjuster(adjusters: Sequence[str], variable: str) -> bool:
    pattern = re.compile(rf"^C\(\s*{re.escape(variable)}\s*(?:,|\))")
    return any(pattern.search(term) for term in adjusters)


def _model_name(base: str, group_label: object | None) -> str:
    if group_label is None:
        return base
    return f"{base}__{_slug(group_label)}"


def _group_metadata(
    group_col: str | None, group_label: object | None
) -> dict[str, object]:
    if group_col is None:
        return {}
    return {"analysis_group_col": group_col, group_col: group_label}


def _iter_groups(
    data: pd.DataFrame,
    group_col: str | None,
    group_values: Sequence[object] | None,
    *,
    outcome: str = "candidate",
) -> Iterable[tuple[object | None, pd.DataFrame]]:
    if group_col is None:
        if data[outcome].nunique(dropna=True) >= 2:
            yield None, data
        return

    observed = set(data[group_col].dropna().unique())
    ordered = [value for value in (group_values or []) if value in observed]
    extras = sorted(observed - set(ordered), key=str)

    for value in [*ordered, *extras]:
        d = data.loc[data[group_col].eq(value)].copy()
        if d.empty:
            continue
        if d[outcome].nunique(dropna=True) < 2:
            print(f"Skipping {group_col}={value!r}: outcome does not vary", flush=True)
            continue
        yield value, d


def add_clade_group(
    data: pd.DataFrame,
    *,
    source_col: str = "clade",
    target_col: str = "clade_group",
) -> pd.DataFrame:
    """Map Nextclade clades onto the curated top-clade labels plus Other."""
    if source_col not in data.columns:
        raise KeyError(f"{source_col!r} is required to build clade groups.")
    out = data.copy()
    out[target_col] = out[source_col].map(CLADES).fillna("Other")
    return out


def add_policy_era(
    data: pd.DataFrame,
    *,
    source_col: str = "policy_period",
    target_col: str = "policy_era",
) -> pd.DataFrame:
    """Attach grouped policy-era labels and policy intensity scores."""
    if source_col not in data.columns:
        raise KeyError(f"{source_col!r} is required to build policy eras.")

    observed = set(data[source_col].dropna().astype(str).unique())
    unmapped = sorted(observed - set(POLICY_ERA_BY_PERIOD))
    if unmapped:
        raise ValueError(
            "Policy period values are not mapped to policy eras: "
            + ", ".join(unmapped)
        )

    out = data.copy()
    period = out[source_col].astype("string")
    out[target_col] = pd.Categorical(
        period.map(POLICY_ERA_BY_PERIOD),
        categories=POLICY_ERA_ORDER,
        ordered=True,
    )
    if "policy_intensity" not in out.columns:
        out["policy_intensity"] = period.map(PERIOD_INTENSITY).astype(float)
    return out


def _add_standardised_columns(
    data: pd.DataFrame,
    mapping: Mapping[str, str],
) -> pd.DataFrame:
    """Add z-scored columns from ``target -> source`` mappings."""
    out = data.copy()
    for target, source in mapping.items():
        if source not in out.columns:
            continue
        values = pd.to_numeric(out[source], errors="coerce")
        sd = values.std(skipna=True)
        if pd.isna(sd) or sd == 0:
            out[target] = np.nan
        else:
            out[target] = (values - values.mean(skipna=True)) / sd
    return out


def add_standardised_adjusters(data: pd.DataFrame) -> pd.DataFrame:
    """Add standardised surveillance and context adjusters used by notebooks."""
    out = data.copy()
    for feature, scaled_feature in zip(
        OBSERVED_MIXING_FEATURES,
        OBSERVED_MIXING_FEATURES_X10,
    ):
        if feature in out.columns:
            out[scaled_feature] = out[feature].astype(float) * 10
    if "wn_positive_tests" in out.columns:
        out["log1p_wn_positive_tests"] = np.log1p(out["wn_positive_tests"])
    if "dz_cum_positive_tests" in out.columns:
        out["log1p_dz_cum_positive_tests"] = np.log1p(out["dz_cum_positive_tests"])
    for target, source in STANDARDISE_SPECS.items():
        if source not in out.columns:
            continue
        values = out[source].astype(float)
        sd = values.std(skipna=True)
        if pd.isna(sd) or sd == 0:
            out[target] = np.nan
        else:
            out[target] = (values - values.mean(skipna=True)) / sd
    return out


def add_vaccination_composition_features(data: pd.DataFrame) -> pd.DataFrame:
    """Add sequence-level vaccination predictors for composition models."""
    required = {
        "is_vaccinated",
        "vacc_dose_number",
        "vacc_booster",
        "days_since_vaccination",
    }
    missing = sorted(required - set(data.columns))
    if missing:
        raise KeyError(f"Missing vaccination composition columns: {missing}")

    out = data.copy()
    dose = pd.to_numeric(out["vacc_dose_number"], errors="coerce").fillna(0)
    out["vacc_dose_cat"] = pd.Categorical(
        np.select(
            [dose.le(0), dose.eq(1), dose.eq(2), dose.ge(3)],
            ["0", "1", "2", "3+"],
            default=None,
        ),
        categories=["0", "1", "2", "3+"],
        ordered=True,
    )

    vaccinated = pd.to_numeric(out["is_vaccinated"], errors="coerce").fillna(0).gt(0)
    booster = pd.to_numeric(out["vacc_booster"], errors="coerce")
    out["vacc_booster_status"] = pd.Categorical(
        np.select(
            [
                ~vaccinated,
                vaccinated & booster.eq(0),
                vaccinated & booster.eq(1),
            ],
            ["unvaccinated", "not_booster", "booster"],
            default=None,
        ),
        categories=["unvaccinated", "not_booster", "booster"],
        ordered=True,
    )

    days = pd.to_numeric(out["days_since_vaccination"], errors="coerce")
    out["days_since_vaccination_cat"] = pd.Categorical(
        np.select(
            [
                ~vaccinated,
                vaccinated & days.between(0, 13, inclusive="both"),
                vaccinated & days.between(14, 89, inclusive="both"),
                vaccinated & days.between(90, 179, inclusive="both"),
                vaccinated & days.ge(180),
            ],
            ["unvaccinated", "0-13", "14-89", "90-179", "180+"],
            default=None,
        ),
        categories=["unvaccinated", "0-13", "14-89", "90-179", "180+"],
        ordered=True,
    )

    return add_standardised_adjusters(out)


def add_vaccination_node_features(
    node_data: pd.DataFrame,
    sequence_data: pd.DataFrame,
    *,
    cluster_col: str = "cluster_id",
) -> pd.DataFrame:
    """Aggregate sequence-level vaccination context onto eligible nodes."""
    required = {
        cluster_col,
        "is_vaccinated",
        "vacc_dose_number",
        "vacc_booster",
        "days_since_vaccination",
        "dz_cum_prop_vaccinated",
    }
    missing = sorted(required - set(sequence_data.columns))
    if missing:
        raise KeyError(f"Missing vaccination node columns: {missing}")
    if cluster_col not in node_data.columns:
        raise KeyError(f"{cluster_col!r} is required in node_data.")

    seq = sequence_data.copy()
    seq["is_vaccinated"] = pd.to_numeric(seq["is_vaccinated"], errors="coerce")
    seq["vacc_dose_number"] = pd.to_numeric(seq["vacc_dose_number"], errors="coerce")
    seq["vacc_booster_filled"] = pd.to_numeric(
        seq["vacc_booster"], errors="coerce"
    ).fillna(0)
    seq["days_since_vaccination"] = pd.to_numeric(
        seq["days_since_vaccination"], errors="coerce"
    )
    seq["dz_cum_prop_vaccinated"] = pd.to_numeric(
        seq["dz_cum_prop_vaccinated"], errors="coerce"
    )

    agg = (
        seq.groupby(cluster_col, dropna=False)
        .agg(
            node_prop_vaccinated=("is_vaccinated", "mean"),
            node_prop_boosted=("vacc_booster_filled", "mean"),
            node_mean_vacc_dose=("vacc_dose_number", "mean"),
            node_median_days_since_vaccination=(
                "days_since_vaccination",
                "median",
            ),
            node_mean_dz_cum_prop_vaccinated=("dz_cum_prop_vaccinated", "mean"),
        )
        .reset_index()
    )
    out = node_data.merge(agg, on=cluster_col, how="left")
    return _add_standardised_columns(
        out,
        {
            "z_node_prop_vaccinated": "node_prop_vaccinated",
            "z_node_prop_boosted": "node_prop_boosted",
            "z_node_mean_vacc_dose": "node_mean_vacc_dose",
            "z_node_median_days_since_vaccination": (
                "node_median_days_since_vaccination"
            ),
            "z_node_mean_dz_cum_prop_vaccinated": (
                "node_mean_dz_cum_prop_vaccinated"
            ),
        },
    )


def _ordered_tertile(values: pd.Series, *, categories: list[str]) -> pd.Categorical:
    """Return rank-based ordered tertiles for a numeric series."""
    out = pd.Series(pd.NA, index=values.index, dtype="object")
    present = values.notna()
    if present.sum() >= 3:
        ranks = values.loc[present].rank(method="first")
        bins = pd.qcut(
            ranks,
            q=3,
            labels=categories,
        )
        out.loc[present] = bins.astype(str).to_numpy()
    return pd.Categorical(out, categories=categories, ordered=True)


def add_vaccination_mixing_features(
    node_data: pd.DataFrame,
    sequence_data: pd.DataFrame,
    *,
    cluster_col: str = "cluster_id",
    vaccination_col: str = "is_vaccinated",
    window_col: str = "window_idx",
    age_col: str = "age_band",
    n_random: int = 1000,
    random_state: int = 42,
) -> pd.DataFrame:
    """Attach age-conditional vaccination-mixing entropy features to nodes."""
    required = {cluster_col, vaccination_col, window_col, age_col}
    missing = sorted(required - set(sequence_data.columns))
    if missing:
        raise KeyError(f"Missing vaccination mixing columns: {missing}")
    if cluster_col not in node_data.columns:
        raise KeyError(f"{cluster_col!r} is required in node_data.")

    seq = sequence_data.copy()
    seq["_vaccination_mix_positive"] = (
        pd.to_numeric(seq[vaccination_col], errors="coerce").fillna(0).gt(0).astype(int)
    )
    seq[age_col] = seq[age_col].astype("string").fillna("Missing")
    with_entropy = cluster_age_conditional_binary_entropy(
        seq,
        cluster_col=cluster_col,
        binary_col="_vaccination_mix_positive",
        window_col=window_col,
        age_col=age_col,
        n_random=n_random,
        random_state=random_state,
        prefix="vaccination_mix",
    )
    feature_cols = [
        "vaccination_mix_n",
        "vaccination_mix_prop_positive",
        "vaccination_mix_entropy_obs",
        "vaccination_mix_entropy_null_mean",
        "vaccination_mix_entropy_null_sd",
        "vaccination_mix_entropy_z",
    ]
    features = with_entropy[[cluster_col, *feature_cols]].drop_duplicates(cluster_col)
    out = node_data.merge(features, on=cluster_col, how="left")
    out["vaccination_mix_tertile"] = _ordered_tertile(
        out["vaccination_mix_entropy_z"],
        categories=VACCINATION_MIXING_TERTILE_ORDER,
    )
    return out


def default_model_sets(
    *,
    variant_adjuster: str | None = "clade",
    window_adjustment: str = "fixed_effects",
) -> dict[str, list[str]]:
    """Build primary and expanded adjustment sets for association models."""
    if window_adjustment == "fixed_effects":
        window_terms = ["C(window_idx)"]
    elif window_adjustment == "surveillance":
        window_terms = list(WINDOW_SURVEILLANCE_ADJUSTERS)
    else:
        raise ValueError("window_adjustment must be 'fixed_effects' or 'surveillance'.")

    variant_terms = [f"C({variant_adjuster})"] if variant_adjuster else []
    primary = [*window_terms, *variant_terms]
    return {
        "primary": primary,
        "expanded": [*primary, *EXPANDED_CONTEXT_ADJUSTERS],
    }


def load_association_frames(
    *,
    output_dir: Path | str | None = None,
    cluster_se: str = "cluster_id",
    variant_adjuster: str | None = "clade",
    group_by_clade: bool = False,
    clade_group_col: str = "clade_group",
    window_stride: int = 2,
    run_composition: bool = True,
    extra_sequence_columns: Sequence[str] | None = None,
) -> AssociationFrames:
    """Load SSE outputs and construct complete model base frames."""
    output_path = (
        Path(output_dir)
        if output_dir
        else PROJECT_ROOT / "sse_detection" / "results" / "sse_outputs"
    )
    outs = load_sse_outputs(output_path)
    node_stats = outs.node_stats.copy()

    if group_by_clade:
        node_stats = add_clade_group(node_stats, target_col=clade_group_col)
        if variant_adjuster == "clade":
            variant_adjuster = None

    min_candidate_size = int(
        node_stats.loc[node_stats["sse_candidate"], "cluster_size"].min()
    )
    eligible_nodes = node_stats.loc[
        node_stats["cluster_size"].ge(min_candidate_size)
    ].copy()

    if cluster_se not in eligible_nodes.columns:
        raise KeyError(f"{cluster_se!r} is not present in node_stats.")
    if variant_adjuster and variant_adjuster not in eligible_nodes.columns:
        raise KeyError(f"{variant_adjuster!r} is not present in node_stats.")

    node_model_base = eligible_nodes.copy()
    node_model_base["candidate"] = node_model_base["sse_candidate"].astype(int)
    if variant_adjuster:
        node_model_base[variant_adjuster] = (
            node_model_base[variant_adjuster].fillna("Missing").astype(str)
        )
    node_model_base[cluster_se] = (
        node_model_base[cluster_se].fillna(node_model_base["cluster_id"]).astype(str)
    )
    node_model_base = add_standardised_adjusters(node_model_base)

    composition_base = pd.DataFrame()
    diagnostics = [
        cluster_se_diagnostics(
            node_model_base,
            cluster_se,
            outcome="candidate",
        ).assign(analysis_frame="node_mixing")
    ]

    if run_composition:
        node_key_cols = [
            "cluster_id",
            "meta_cluster_id",
            "sse_candidate",
            cluster_se,
            "cluster_size",
        ]
        if group_by_clade:
            node_key_cols.append(clade_group_col)
        node_key = eligible_nodes[_dedupe(node_key_cols)].drop_duplicates("cluster_id")

        sequence_columns: list[str] = list(
            {
                "window_id",
                "window_idx",
                "cluster_id",
                "sequence_id",
                "wn_prop_sequenced",
                "wn_positive_tests",
                "dz_cum_prop_sequenced",
                "dz_cum_incidence_per_capita",
                "dz_7d_test_positivity",
                "dz_cum_positive_tests",
                *(spec["column"] for spec in COMPOSITION_SPECS),
                *([variant_adjuster] if variant_adjuster else []),
                *(extra_sequence_columns or []),
            }
        )

        sequence_raw = load_analysis_columns(
            sequence_columns,
            add_policy=False,
            window_stride=window_stride,
        )
        composition_base = sequence_raw.merge(node_key, on="cluster_id", how="inner")
        composition_base["candidate"] = composition_base["sse_candidate"].astype(int)
        if variant_adjuster:
            composition_base[variant_adjuster] = (
                composition_base[variant_adjuster].fillna("Missing").astype(str)
            )
        composition_base[cluster_se] = (
            composition_base[cluster_se]
            .fillna(composition_base["cluster_id"])
            .astype(str)
        )
        composition_base = add_standardised_adjusters(composition_base)
        diagnostics.insert(
            0,
            cluster_se_diagnostics(
                composition_base,
                cluster_se,
                outcome="candidate",
            ).assign(analysis_frame="composition"),
        )

    return AssociationFrames(
        node_stats=node_stats,
        eligible_nodes=eligible_nodes,
        node_model_base=node_model_base,
        composition_base=composition_base,
        cluster_diagnostics=pd.concat(diagnostics, ignore_index=True),
        min_candidate_size=min_candidate_size,
    )


def fit_association_result(
    data: pd.DataFrame,
    formula: str,
    *,
    model_method: str,
    cluster_se: str,
    window_strata: str,
):
    if model_method == "conditional_logit_by_window":
        return fit_conditional_logit(data, formula, strata_col=window_strata)
    if model_method == "firth_glm":
        return fit_firth_logit(data, formula)
    if model_method == "glm_clustered":
        return fit_binomial_glm(data, formula, cluster_col=cluster_se)
    raise ValueError(f"Unknown model_method={model_method!r}")


def fit_exposure_association(
    data: pd.DataFrame,
    *,
    outcome: str,
    exposure: str,
    adjusters: list[str],
    model_name: str,
    model_method: str,
    cluster_se: str,
    window_strata: str,
    reference=None,
    categorical: bool = True,
) -> AssociationModel:
    exposure_term = (
        categorical_term(exposure, reference=reference) if categorical else exposure
    )
    formula = make_formula(outcome, exposure_term, adjusters)
    result = fit_association_result(
        data,
        formula,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
    )
    odds = tidy_odds_ratios(
        result,
        model_name=model_name,
        term_filter=exposure_term,
    )
    if categorical:
        wald = robust_wald_for_prefix(
            result,
            exposure_term,
            model_name=model_name,
            term=exposure,
        )
    else:
        wald = tidy_single_parameter_wald(
            result,
            [exposure],
            model_name=model_name,
        )
    wald["formula"] = formula
    return AssociationModel(result=result, odds_ratios=odds, wald=wald, formula=formula)


def resolve_reference(data: pd.DataFrame, column: str, preferred, fallbacks=None):
    levels = set(data[column].dropna().astype(str))
    candidates = [preferred, *(fallbacks or [])]
    for ref in candidates:
        if ref is None:
            continue
        ref_str = str(ref)
        if ref_str in levels:
            return ref_str
        for level in levels:
            try:
                if float(level) == float(ref_str):
                    return level
            except ValueError:
                pass
    counts = data[column].dropna().astype(str).value_counts()
    if counts.empty:
        raise ValueError(f"No observed levels for {column!r}.")
    return str(counts.index[0])


def complete_case(data: pd.DataFrame, required: list[str]) -> pd.DataFrame:
    missing = [col for col in required if col not in data.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")
    return data.dropna(subset=required).copy()


def drop_nonvarying_levels(
    data: pd.DataFrame,
    columns: list[str],
    *,
    outcome: str = "candidate",
) -> tuple[pd.DataFrame, dict[str, int], dict[str, int]]:
    d = data.copy()
    dropped_rows: dict[str, int] = {}
    dropped_strata: dict[str, int] = {}
    changed = True
    while changed:
        changed = False
        for col in columns:
            if col not in d.columns:
                continue
            strata_nunique = d.groupby(col, dropna=False)[outcome].transform("nunique")
            varies = strata_nunique.gt(1)
            n_drop = int((~varies).sum())
            if n_drop:
                dropped_rows[col] = dropped_rows.get(col, 0) + n_drop
                dropped_strata[col] = dropped_strata.get(col, 0) + int(
                    d.loc[~varies, col].nunique(dropna=False)
                )
                d = d.loc[varies].copy()
                changed = True
    return d, dropped_rows, dropped_strata


def dropped_metadata(
    dropped_rows: dict[str, int], dropped_strata: dict[str, int]
) -> dict[str, object]:
    return {
        "dropped_nonvarying_rows": sum(dropped_rows.values()),
        "dropped_nonvarying_strata": sum(dropped_strata.values()),
        "dropped_nonvarying_detail": repr(
            {
                "rows": dropped_rows,
                "strata": dropped_strata,
            }
        ),
    }


def add_model_metadata(table: pd.DataFrame, **metadata) -> pd.DataFrame:
    out = table.copy()
    for key, value in reversed(list(metadata.items())):
        out.insert(0, key, value)
    return out


def add_fit_metadata(fit_stats: pd.DataFrame, **metadata) -> pd.DataFrame:
    out = fit_stats.copy()
    for key, value in reversed(list(metadata.items())):
        out.insert(0, key, value)
    return out


def bh_adjust_by(
    table: pd.DataFrame,
    group_cols: list[str],
    p_col: str = "P>chi2",
) -> pd.DataFrame:
    if table.empty:
        return table.copy()
    out = table.copy()
    out["p_adj_bh"] = np.nan
    for _, idx in out.groupby(group_cols, dropna=False).groups.items():
        adjusted = bh_adjust(out.loc[idx], p_col=p_col)
        out.loc[idx, "p_adj_bh"] = adjusted["p_adj_bh"].to_numpy()
    return out


def _prepare_model_frame(
    source: pd.DataFrame,
    *,
    predictors: list[str],
    adjusters: list[str],
    required_base: list[str],
    window_strata: str,
    drop_window_nonvarying: bool,
    categorical_predictors: bool = False,
) -> tuple[pd.DataFrame, dict[str, int], dict[str, int]]:
    required = required_base + predictors + model_variables_from_terms(adjusters)
    required = _dedupe(required)
    d = complete_case(source, required)
    if categorical_predictors:
        for col in predictors:
            d[col] = d[col].astype(str)
    strata = [window_strata] if drop_window_nonvarying else []
    d, dropped_rows, dropped_strata = drop_nonvarying_levels(d, strata)
    if d.empty:
        raise ValueError("No complete-case rows remain after filtering.")
    if d["candidate"].nunique(dropna=True) < 2:
        raise ValueError("Outcome does not vary after filtering.")
    return d, dropped_rows, dropped_strata


def _record_failure(
    failures: list[dict[str, object]],
    *,
    domain: str,
    model_set: str,
    predictor_set: str,
    predictor: str,
    group_col: str | None,
    group_label: object | None,
    error: Exception,
) -> None:
    row = {
        "domain": domain,
        "model_set": model_set,
        "predictor_set": predictor_set,
        "predictor": predictor,
        "error_type": type(error).__name__,
        "error": str(error),
    }
    row |= _group_metadata(group_col, group_label)
    failures.append(row)
    where = f", {group_col}={group_label!r}" if group_col else ""
    print(
        f"Failed {domain} {model_set} {predictor_set} {predictor}{where}: {error}",
        flush=True,
    )


def fit_single_composition_models(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    wald_tables = []
    or_tables = []
    fit_tables = []
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )

    for spec in COMPOSITION_SPECS:
        predictor = spec["column"]
        try:
            d, dropped_rows, dropped_strata = _prepare_model_frame(
                source,
                predictors=[predictor],
                adjusters=adjusters,
                required_base=[
                    "candidate",
                    "cluster_id",
                    "sequence_id",
                    cluster_se,
                    window_strata,
                ],
                window_strata=window_strata,
                drop_window_nonvarying=drop_window,
                categorical_predictors=True,
            )
            reference = resolve_reference(
                d,
                predictor,
                spec.get("reference"),
                spec.get("fallback_references"),
            )
            base_name = f"composition__{model_set}__single__{spec['name']}"
            model_name = _model_name(base_name, group_label)
            fit = fit_exposure_association(
                d,
                outcome="candidate",
                exposure=predictor,
                adjusters=adjusters,
                model_name=model_name,
                model_method=model_method,
                cluster_se=cluster_se,
                window_strata=window_strata,
                reference=reference,
                categorical=True,
            )
            meta = {
                "domain": "composition",
                "model_set": model_set,
                "predictor_set": "single",
                "predictor": spec["name"],
                "label": spec["label"],
                "reference": reference,
                "n_model_rows": len(d),
                "n_sequences": d["sequence_id"].nunique(),
                "n_nodes": d["cluster_id"].nunique(),
                **_group_metadata(group_col, group_label),
                **dropped_metadata(dropped_rows, dropped_strata),
            }
            wald_tables.append(add_model_metadata(fit.wald, **meta))
            or_tables.append(add_model_metadata(fit.odds_ratios, **meta))
            fit_tables.append(
                add_fit_metadata(
                    model_fit_stats(
                        fit.result, model_name=model_name, formula=fit.formula
                    ),
                    **meta,
                )
            )
            print(f"Fitted {model_name}: {len(d):,} rows", flush=True)
        except Exception as exc:  # keep other clade/predictor fits moving
            _record_failure(
                failures,
                domain="composition",
                model_set=model_set,
                predictor_set="single",
                predictor=spec["name"],
                group_col=group_col,
                group_label=group_label,
                error=exc,
            )

    return {
        "wald": _concat_or_empty(wald_tables),
        "odds": _concat_or_empty(or_tables),
        "fit_stats": _concat_or_empty(fit_tables),
    }


def fit_joint_composition_model(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )
    predictors = [spec["column"] for spec in COMPOSITION_SPECS]
    try:
        d, dropped_rows, dropped_strata = _prepare_model_frame(
            source,
            predictors=predictors,
            adjusters=adjusters,
            required_base=[
                "candidate",
                "cluster_id",
                "sequence_id",
                cluster_se,
                window_strata,
            ],
            window_strata=window_strata,
            drop_window_nonvarying=drop_window,
            categorical_predictors=True,
        )

        terms = []
        references = {}
        for spec in COMPOSITION_SPECS:
            reference = resolve_reference(
                d,
                spec["column"],
                spec.get("reference"),
                spec.get("fallback_references"),
            )
            references[spec["name"]] = reference
            terms.append(categorical_term(spec["column"], reference))

        formula = "candidate ~ " + " + ".join(terms + adjusters)
        model_name = _model_name(f"composition__{model_set}__joint", group_label)
        result = fit_association_result(
            d,
            formula,
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
        )

        wald_tables = []
        for spec, term in zip(COMPOSITION_SPECS, terms):
            wald = robust_wald_for_prefix(
                result,
                term,
                model_name=model_name,
                term=spec["name"],
            )
            meta = {
                "domain": "composition",
                "model_set": model_set,
                "predictor_set": "joint",
                "predictor": spec["name"],
                "label": spec["label"],
                "reference": references[spec["name"]],
                "n_model_rows": len(d),
                "n_sequences": d["sequence_id"].nunique(),
                "n_nodes": d["cluster_id"].nunique(),
                **_group_metadata(group_col, group_label),
                **dropped_metadata(dropped_rows, dropped_strata),
            }
            wald_tables.append(add_model_metadata(wald, **meta))

        term_names = [
            name for term in terms for name in parameter_names_for_term(result, term)
        ]
        odds = tidy_odds_ratios(result, model_name=model_name)
        odds = odds.loc[odds["term"].isin(term_names)].copy()
        common_meta = {
            "domain": "composition",
            "model_set": model_set,
            "predictor_set": "joint",
            "predictor": "all_composition",
            "label": "All composition predictors",
            "reference": repr(references),
            "n_model_rows": len(d),
            "n_sequences": d["sequence_id"].nunique(),
            "n_nodes": d["cluster_id"].nunique(),
            **_group_metadata(group_col, group_label),
            **dropped_metadata(dropped_rows, dropped_strata),
        }
        fit_stats = model_fit_stats(result, model_name=model_name, formula=formula)
        print(f"Fitted {model_name}: {len(d):,} rows", flush=True)
        return {
            "wald": _concat_or_empty(wald_tables),
            "odds": add_model_metadata(odds, **common_meta),
            "fit_stats": add_fit_metadata(fit_stats, **common_meta),
        }
    except Exception as exc:
        _record_failure(
            failures,
            domain="composition",
            model_set=model_set,
            predictor_set="joint",
            predictor="all_composition",
            group_col=group_col,
            group_label=group_label,
            error=exc,
        )
        return {
            "wald": pd.DataFrame(),
            "odds": pd.DataFrame(),
            "fit_stats": pd.DataFrame(),
        }


def run_composition_model_set(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    single = fit_single_composition_models(
        source=source,
        model_set=model_set,
        adjusters=adjusters,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
        group_col=group_col,
        group_label=group_label,
        failures=failures,
    )
    joint = fit_joint_composition_model(
        source=source,
        model_set=model_set,
        adjusters=adjusters,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
        group_col=group_col,
        group_label=group_label,
        failures=failures,
    )
    return {
        "wald": _concat_or_empty([single["wald"], joint["wald"]]),
        "odds": _concat_or_empty([single["odds"], joint["odds"]]),
        "fit_stats": _concat_or_empty([single["fit_stats"], joint["fit_stats"]]),
    }


def fit_single_mixing_models(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    mixing_features: Sequence[str],
    mixing_reference: str,
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    wald_tables = []
    or_tables = []
    fit_tables = []
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )

    for feature in mixing_features:
        if feature not in source.columns:
            print(f"Skipping {feature}: not found", flush=True)
            continue
        try:
            d, dropped_rows, dropped_strata = _prepare_model_frame(
                source,
                predictors=[feature],
                adjusters=adjusters,
                required_base=["candidate", "cluster_id", cluster_se, window_strata],
                window_strata=window_strata,
                drop_window_nonvarying=drop_window,
            )
            base_name = f"mixing__{model_set}__single__{feature}"
            model_name = _model_name(base_name, group_label)
            fit = fit_exposure_association(
                d,
                outcome="candidate",
                exposure=feature,
                adjusters=adjusters,
                model_name=model_name,
                model_method=model_method,
                cluster_se=cluster_se,
                window_strata=window_strata,
                categorical=False,
            )
            meta = {
                "domain": "node_mixing",
                "model_set": model_set,
                "predictor_set": "single",
                "predictor": feature,
                "label": feature.replace("_", " "),
                "reference": mixing_reference,
                "n_model_rows": len(d),
                "n_nodes": d["cluster_id"].nunique(),
                **_group_metadata(group_col, group_label),
                **dropped_metadata(dropped_rows, dropped_strata),
            }
            wald_tables.append(add_model_metadata(fit.wald, **meta))
            or_tables.append(
                add_model_metadata(
                    fit.odds_ratios.loc[fit.odds_ratios["term"].eq(feature)].copy(),
                    **meta,
                )
            )
            fit_tables.append(
                add_fit_metadata(
                    model_fit_stats(
                        fit.result, model_name=model_name, formula=fit.formula
                    ),
                    **meta,
                )
            )
            print(f"Fitted {model_name}: {len(d):,} nodes", flush=True)
        except Exception as exc:
            _record_failure(
                failures,
                domain="node_mixing",
                model_set=model_set,
                predictor_set="single",
                predictor=feature,
                group_col=group_col,
                group_label=group_label,
                error=exc,
            )

    return {
        "wald": _concat_or_empty(wald_tables),
        "odds": _concat_or_empty(or_tables),
        "fit_stats": _concat_or_empty(fit_tables),
    }


def fit_joint_mixing_model(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    mixing_features: Sequence[str],
    mixing_reference: str,
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    features = [feature for feature in mixing_features if feature in source.columns]
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )
    try:
        d, dropped_rows, dropped_strata = _prepare_model_frame(
            source,
            predictors=features,
            adjusters=adjusters,
            required_base=["candidate", "cluster_id", cluster_se, window_strata],
            window_strata=window_strata,
            drop_window_nonvarying=drop_window,
        )
        formula = "candidate ~ " + " + ".join(features + adjusters)
        model_name = _model_name(f"mixing__{model_set}__joint", group_label)
        result = fit_association_result(
            d,
            formula,
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
        )

        common_meta = {
            "domain": "node_mixing",
            "model_set": model_set,
            "predictor_set": "joint",
            "predictor": "all_mixing",
            "label": "All mixing predictors",
            "reference": mixing_reference,
            "n_model_rows": len(d),
            "n_nodes": d["cluster_id"].nunique(),
            **_group_metadata(group_col, group_label),
            **dropped_metadata(dropped_rows, dropped_strata),
        }
        wald = robust_wald_for_params(
            result,
            features,
            model_name=model_name,
            term="all_mixing",
        )
        odds = tidy_odds_ratios(result, model_name=model_name)
        odds = odds.loc[odds["term"].isin(features)].copy()
        fit_stats = model_fit_stats(result, model_name=model_name, formula=formula)
        print(f"Fitted {model_name}: {len(d):,} nodes", flush=True)
        return {
            "wald": add_model_metadata(wald, **common_meta),
            "odds": add_model_metadata(odds, **common_meta),
            "fit_stats": add_fit_metadata(fit_stats, **common_meta),
        }
    except Exception as exc:
        _record_failure(
            failures,
            domain="node_mixing",
            model_set=model_set,
            predictor_set="joint",
            predictor="all_mixing",
            group_col=group_col,
            group_label=group_label,
            error=exc,
        )
        return {
            "wald": pd.DataFrame(),
            "odds": pd.DataFrame(),
            "fit_stats": pd.DataFrame(),
        }


def run_mixing_model_set(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    mixing_features: Sequence[str],
    mixing_reference: str,
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    single = fit_single_mixing_models(
        source=source,
        model_set=model_set,
        adjusters=adjusters,
        mixing_features=mixing_features,
        mixing_reference=mixing_reference,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
        group_col=group_col,
        group_label=group_label,
        failures=failures,
    )
    joint = fit_joint_mixing_model(
        source=source,
        model_set=model_set,
        adjusters=adjusters,
        mixing_features=mixing_features,
        mixing_reference=mixing_reference,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
        group_col=group_col,
        group_label=group_label,
        failures=failures,
    )
    return {
        "wald": _concat_or_empty([single["wald"], joint["wald"]]),
        "odds": _concat_or_empty([single["odds"], joint["odds"]]),
        "fit_stats": _concat_or_empty([single["fit_stats"], joint["fit_stats"]]),
    }


def fit_single_exposure_specs(
    *,
    source: pd.DataFrame,
    specs: Sequence[Mapping[str, object]],
    domain: str,
    model_set: str,
    adjusters: list[str],
    model_method: str,
    cluster_se: str,
    window_strata: str,
    required_base: list[str],
    failures: list[dict[str, object]],
    n_label: str = "rows",
) -> dict[str, pd.DataFrame]:
    """Fit one single-predictor model for each exposure spec."""
    wald_tables = []
    or_tables = []
    fit_tables = []
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )

    for spec in specs:
        predictor = str(spec["column"])
        predictor_name = str(spec.get("name", predictor))
        categorical = bool(spec.get("categorical", False))
        try:
            d, dropped_rows, dropped_strata = _prepare_model_frame(
                source,
                predictors=[predictor],
                adjusters=adjusters,
                required_base=required_base,
                window_strata=window_strata,
                drop_window_nonvarying=drop_window,
                categorical_predictors=categorical,
            )
            if categorical:
                reference = resolve_reference(
                    d,
                    predictor,
                    spec.get("reference"),
                    spec.get("fallback_references"),
                )
            else:
                reference = spec.get("reference")

            model_name = f"{domain}__{model_set}__single__{predictor_name}"
            fit = fit_exposure_association(
                d,
                outcome="candidate",
                exposure=predictor,
                adjusters=adjusters,
                model_name=model_name,
                model_method=model_method,
                cluster_se=cluster_se,
                window_strata=window_strata,
                reference=reference,
                categorical=categorical,
            )
            meta = {
                "domain": domain,
                "model_set": model_set,
                "predictor_set": "single",
                "predictor": predictor_name,
                "label": spec.get("label", predictor_name.replace("_", " ")),
                "reference": reference,
                "n_model_rows": len(d),
                "n_nodes": d["cluster_id"].nunique()
                if "cluster_id" in d.columns
                else np.nan,
                **dropped_metadata(dropped_rows, dropped_strata),
            }
            if "sequence_id" in d.columns:
                meta["n_sequences"] = d["sequence_id"].nunique()

            wald_tables.append(add_model_metadata(fit.wald, **meta))
            or_tables.append(add_model_metadata(fit.odds_ratios, **meta))
            fit_tables.append(
                add_fit_metadata(
                    model_fit_stats(
                        fit.result,
                        model_name=model_name,
                        formula=fit.formula,
                    ),
                    **meta,
                )
            )
            print(f"Fitted {model_name}: {len(d):,} {n_label}", flush=True)
        except Exception as exc:
            _record_failure(
                failures,
                domain=domain,
                model_set=model_set,
                predictor_set="single",
                predictor=predictor_name,
                group_col=None,
                group_label=None,
                error=exc,
            )

    return {
        "wald": _concat_or_empty(wald_tables),
        "odds": _concat_or_empty(or_tables),
        "fit_stats": _concat_or_empty(fit_tables),
    }


def fit_joint_exposure_specs(
    *,
    source: pd.DataFrame,
    specs: Sequence[Mapping[str, object]],
    joint_groups: Sequence[Mapping[str, object]],
    domain: str,
    model_set: str,
    adjusters: list[str],
    model_method: str,
    cluster_se: str,
    window_strata: str,
    required_base: list[str],
    failures: list[dict[str, object]],
    n_label: str = "rows",
) -> dict[str, pd.DataFrame]:
    """Fit selected joint models over exposure specs.

    Joint groups intentionally remain small because many vaccination fields are
    structurally correlated. Each group is fitted once, then Wald and OR rows
    are emitted for the member predictors.
    """
    spec_by_name = {str(spec.get("name", spec["column"])): spec for spec in specs}
    wald_tables = []
    or_tables = []
    fit_tables = []
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )

    for group in joint_groups:
        joint_name = str(group["name"])
        try:
            group_specs = [spec_by_name[str(name)] for name in group.get("specs", [])]
            predictors = [str(spec["column"]) for spec in group_specs]
            d, dropped_rows, dropped_strata = _prepare_model_frame(
                source,
                predictors=predictors,
                adjusters=adjusters,
                required_base=required_base,
                window_strata=window_strata,
                drop_window_nonvarying=drop_window,
                categorical_predictors=False,
            )

            terms = []
            references = {}
            for spec in group_specs:
                predictor = str(spec["column"])
                predictor_name = str(spec.get("name", predictor))
                if bool(spec.get("categorical", False)):
                    reference = resolve_reference(
                        d,
                        predictor,
                        spec.get("reference"),
                        spec.get("fallback_references"),
                    )
                    term = categorical_term(predictor, reference)
                else:
                    reference = spec.get("reference")
                    term = predictor
                terms.append(term)
                references[predictor_name] = reference

            formula = "candidate ~ " + " + ".join(terms + adjusters)
            model_name = f"{domain}__{model_set}__joint__{joint_name}"
            result = fit_association_result(
                d,
                formula,
                model_method=model_method,
                cluster_se=cluster_se,
                window_strata=window_strata,
            )
            all_odds = tidy_odds_ratios(result, model_name=model_name)

            for spec, term in zip(group_specs, terms):
                predictor = str(spec["column"])
                predictor_name = str(spec.get("name", predictor))
                if bool(spec.get("categorical", False)):
                    wald = robust_wald_for_prefix(
                        result,
                        term,
                        model_name=model_name,
                        term=predictor_name,
                    )
                else:
                    wald = robust_wald_for_params(
                        result,
                        [predictor],
                        model_name=model_name,
                        term=predictor_name,
                    )

                meta = {
                    "domain": domain,
                    "model_set": model_set,
                    "predictor_set": "joint",
                    "joint_model": joint_name,
                    "predictor": predictor_name,
                    "label": spec.get("label", predictor_name.replace("_", " ")),
                    "reference": references[predictor_name],
                    "n_model_rows": len(d),
                    "n_nodes": d["cluster_id"].nunique()
                    if "cluster_id" in d.columns
                    else np.nan,
                    **dropped_metadata(dropped_rows, dropped_strata),
                }
                if "sequence_id" in d.columns:
                    meta["n_sequences"] = d["sequence_id"].nunique()

                term_names = parameter_names_for_term(result, term)
                odds = all_odds.loc[all_odds["term"].isin(term_names)].copy()
                wald_tables.append(add_model_metadata(wald, **meta))
                or_tables.append(add_model_metadata(odds, **meta))

            fit_meta = {
                "domain": domain,
                "model_set": model_set,
                "predictor_set": "joint",
                "joint_model": joint_name,
                "predictor": joint_name,
                "label": group.get("label", joint_name.replace("_", " ")),
                "reference": repr(references),
                "n_model_rows": len(d),
                "n_nodes": d["cluster_id"].nunique()
                if "cluster_id" in d.columns
                else np.nan,
                **dropped_metadata(dropped_rows, dropped_strata),
            }
            if "sequence_id" in d.columns:
                fit_meta["n_sequences"] = d["sequence_id"].nunique()
            fit_tables.append(
                add_fit_metadata(
                    model_fit_stats(
                        result,
                        model_name=model_name,
                        formula=formula,
                    ),
                    **fit_meta,
                )
            )
            print(f"Fitted {model_name}: {len(d):,} {n_label}", flush=True)
        except Exception as exc:
            _record_failure(
                failures,
                domain=domain,
                model_set=model_set,
                predictor_set="joint",
                predictor=joint_name,
                group_col=None,
                group_label=None,
                error=exc,
            )
            if failures:
                failures[-1]["joint_model"] = joint_name

    return {
        "wald": _concat_or_empty(wald_tables),
        "odds": _concat_or_empty(or_tables),
        "fit_stats": _concat_or_empty(fit_tables),
    }


def make_policy_era_candidate_summary(source: pd.DataFrame) -> pd.DataFrame:
    """Summarise eligible-node candidate frequency by grouped policy era."""
    d = source.copy()
    grouped = d.groupby("policy_era", observed=False, dropna=False)
    out = grouped.agg(
        n_nodes=("cluster_id", "nunique"),
        n_candidates=("candidate", "sum"),
        n_sequences=("cluster_size", "sum"),
        mean_cluster_size=("cluster_size", "mean"),
        median_cluster_size=("cluster_size", "median"),
        mean_policy_intensity=("policy_intensity", "mean"),
        policy_periods=(
            "policy_period",
            lambda s: "; ".join(sorted(s.dropna().astype(str).unique())),
        ),
    ).reset_index()
    out["candidate_rate"] = np.where(
        out["n_nodes"].gt(0),
        out["n_candidates"] / out["n_nodes"],
        np.nan,
    )
    return out


def make_policy_era_category_summary(source: pd.DataFrame) -> pd.DataFrame:
    """Summarise candidate category mix by grouped policy era."""
    d = source.loc[source["candidate"].eq(1)].copy()
    if d.empty:
        return pd.DataFrame()
    out = (
        d.groupby(["policy_era", "sse_category"], observed=False, dropna=False)
        .agg(n_candidates=("cluster_id", "nunique"))
        .reset_index()
    )
    out = out.loc[out["n_candidates"].gt(0)].copy()
    totals = out.groupby("policy_era", observed=False)["n_candidates"].transform("sum")
    out["candidate_category_share"] = out["n_candidates"] / totals
    return out


def make_vaccination_candidate_summary(
    sequence_source: pd.DataFrame,
    node_source: pd.DataFrame,
) -> pd.DataFrame:
    """Summarise vaccination context by candidate status at sequence and node levels."""
    seq = sequence_source.copy()
    seq["vacc_booster_filled"] = pd.to_numeric(
        seq["vacc_booster"], errors="coerce"
    ).fillna(0)
    sequence_summary = (
        seq.groupby("candidate", dropna=False)
        .agg(
            n_rows=("sequence_id", "size"),
            n_sequences=("sequence_id", "nunique"),
            n_nodes=("cluster_id", "nunique"),
            prop_vaccinated=("is_vaccinated", "mean"),
            mean_vacc_dose=("vacc_dose_number", "mean"),
            prop_boosted=("vacc_booster_filled", "mean"),
            median_days_since_vaccination=("days_since_vaccination", "median"),
            mean_dz_vaccination_events_per_capita=(
                "dz_cum_prop_vaccinated",
                "mean",
            ),
        )
        .reset_index()
    )
    sequence_summary.insert(0, "analysis_level", "sequence_composition")

    node_summary = (
        node_source.groupby("candidate", dropna=False)
        .agg(
            n_rows=("cluster_id", "size"),
            n_sequences=("cluster_size", "sum"),
            n_nodes=("cluster_id", "nunique"),
            prop_vaccinated=("node_prop_vaccinated", "mean"),
            mean_vacc_dose=("node_mean_vacc_dose", "mean"),
            prop_boosted=("node_prop_boosted", "mean"),
            median_days_since_vaccination=(
                "node_median_days_since_vaccination",
                "median",
            ),
            mean_dz_vaccination_events_per_capita=(
                "node_mean_dz_cum_prop_vaccinated",
                "mean",
            ),
        )
        .reset_index()
    )
    node_summary.insert(0, "analysis_level", "node_context")
    return pd.concat([sequence_summary, node_summary], ignore_index=True)


def make_vaccination_mixing_age_conditional_summary(
    source: pd.DataFrame,
) -> pd.DataFrame:
    """Summarise age-conditional vaccination mixing by tertile."""
    d = source.dropna(subset=["vaccination_mix_tertile"]).copy()
    if d.empty:
        return pd.DataFrame(
            columns=[
                "vaccination_mix_tertile",
                "n_nodes",
                "n_candidates",
                "candidate_rate",
                "mean_vaccination_mix_entropy_z",
                "median_vaccination_mix_entropy_z",
                "mean_vaccination_mix_entropy_obs",
                "mean_vaccination_mix_null_mean",
                "mean_vaccination_mix_prop_positive",
            ]
        )
    grouped = d.groupby("vaccination_mix_tertile", observed=False, dropna=False)
    out = (
        grouped.agg(
            n_nodes=("cluster_id", "nunique"),
            n_candidates=("candidate", "sum"),
            mean_vaccination_mix_entropy_z=("vaccination_mix_entropy_z", "mean"),
            median_vaccination_mix_entropy_z=("vaccination_mix_entropy_z", "median"),
            mean_vaccination_mix_entropy_obs=("vaccination_mix_entropy_obs", "mean"),
            mean_vaccination_mix_null_mean=(
                "vaccination_mix_entropy_null_mean",
                "mean",
            ),
            mean_vaccination_mix_prop_positive=(
                "vaccination_mix_prop_positive",
                "mean",
            ),
        )
        .reset_index()
    )
    out["candidate_rate"] = out["n_candidates"] / out["n_nodes"]
    return out


def make_vaccination_mixing_age_conditional_category_summary(
    source: pd.DataFrame,
) -> pd.DataFrame:
    """Summarise candidate phenotype mix by vaccination-mixing tertile."""
    d = source.dropna(subset=["vaccination_mix_tertile", "sse_category"]).copy()
    d = d.loc[d["candidate"].astype(bool)].copy()
    if d.empty:
        return pd.DataFrame(
            columns=[
                "vaccination_mix_tertile",
                "sse_category",
                "n_candidates",
                "candidate_category_share",
            ]
        )
    out = (
        d.groupby(
            ["vaccination_mix_tertile", "sse_category"],
            observed=False,
            dropna=False,
        )
        .agg(n_candidates=("cluster_id", "nunique"))
        .reset_index()
    )
    out = out.loc[out["n_candidates"].gt(0)].copy()
    totals = out.groupby("vaccination_mix_tertile", observed=False)[
        "n_candidates"
    ].transform("sum")
    out["candidate_category_share"] = out["n_candidates"] / totals
    return out


def vaccination_mixing_node_feature_table(source: pd.DataFrame) -> pd.DataFrame:
    """Return the node-level columns needed to audit vaccination-mixing features."""
    columns = [
        "cluster_id",
        "window_idx",
        "clade",
        "candidate",
        "sse_candidate",
        "sse_category",
        "cluster_size",
        "vaccination_mix_n",
        "vaccination_mix_prop_positive",
        "vaccination_mix_entropy_obs",
        "vaccination_mix_entropy_null_mean",
        "vaccination_mix_entropy_null_sd",
        "vaccination_mix_entropy_z",
        "vaccination_mix_tertile",
        "age_entropy_z",
    ]
    return source[[col for col in columns if col in source.columns]].copy()


def _export_summary_tables(
    summary_tables: Mapping[str, pd.DataFrame],
    result_path: Path,
) -> None:
    for filename, table in summary_tables.items():
        path = result_path / filename
        clean_export_table(table).to_csv(path, index=False)
        print(f"saved {filename}: {len(table):,} rows", flush=True)


def clean_export_table(table: pd.DataFrame) -> pd.DataFrame:
    out = table.copy()
    out.columns = [str(col).strip() for col in out.columns]
    for col in out.select_dtypes(include=["object", "string"]).columns:
        present = out[col].notna()
        out.loc[present, col] = out.loc[present, col].astype(str).str.strip()
    return out


def select_table_columns(
    table: pd.DataFrame,
    kind: str,
    *,
    sort_by: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Return a notebook-display view with available standard columns only."""
    if kind not in TABLE_DISPLAY_COLUMNS:
        raise KeyError(
            f"Unknown table kind {kind!r}. "
            f"Expected one of: {sorted(TABLE_DISPLAY_COLUMNS)}."
        )
    columns = [col for col in TABLE_DISPLAY_COLUMNS[kind] if col in table.columns]
    out = table.loc[:, columns].copy()
    if sort_by:
        present_sort = [col for col in sort_by if col in out.columns]
        if present_sort:
            out = out.sort_values(present_sort)
    return out


def _summarise_tables(
    *,
    composition_results: list[dict[str, pd.DataFrame]],
    mixing_results: list[dict[str, pd.DataFrame]],
) -> dict[str, pd.DataFrame]:
    summary: dict[str, pd.DataFrame] = {}

    if composition_results:
        composition_wald = _concat_or_empty(
            [result["wald"] for result in composition_results]
        )
        composition_or = _concat_or_empty(
            [result["odds"] for result in composition_results]
        )
        composition_fit_stats = _concat_or_empty(
            [result["fit_stats"] for result in composition_results]
        )
        if not composition_wald.empty:
            group_cols = ["domain", "model_set", "predictor_set"]
            if "analysis_group_col" in composition_wald.columns:
                group_cols.extend(
                    col
                    for col in composition_wald["analysis_group_col"]
                    .dropna()
                    .astype(str)
                    .unique()
                    if col in composition_wald.columns
                )
            composition_wald = bh_adjust_by(composition_wald, group_cols)
        summary.update(
            {
                "composition_wald.csv": composition_wald,
                "composition_odds_ratios.csv": composition_or,
                "composition_fit_stats.csv": composition_fit_stats,
            }
        )

    if mixing_results:
        mixing_wald = _concat_or_empty([result["wald"] for result in mixing_results])
        mixing_or = _concat_or_empty([result["odds"] for result in mixing_results])
        mixing_fit_stats = _concat_or_empty(
            [result["fit_stats"] for result in mixing_results]
        )
        if not mixing_wald.empty:
            group_cols = ["domain", "model_set", "predictor_set"]
            if "analysis_group_col" in mixing_wald.columns:
                group_cols.extend(
                    col
                    for col in mixing_wald["analysis_group_col"]
                    .dropna()
                    .astype(str)
                    .unique()
                    if col in mixing_wald.columns
                )
            mixing_wald = bh_adjust_by(mixing_wald, group_cols)
        summary.update(
            {
                "mixing_wald.csv": mixing_wald,
                "mixing_odds_ratios.csv": mixing_or,
                "mixing_fit_stats.csv": mixing_fit_stats,
            }
        )

    return summary


def run_association_pipeline(
    *,
    output_dir: Path | str | None = None,
    result_dir: Path | str | None = None,
    result_subdir: str = "association_outputs",
    cluster_se: str = "cluster_id",
    window_strata: str = "window_idx",
    model_method: str = "firth_glm",
    variant_adjuster: str | None = "clade",
    window_adjustment: str = "fixed_effects",
    composition_model_sets: Mapping[str, Sequence[str]] | None = None,
    mixing_model_sets: Mapping[str, Sequence[str]] | None = None,
    mixing_features: Sequence[str] = DEFAULT_MIXING_FEATURES,
    mixing_reference: str = "per 1 null-model SD increase in entropy",
    group_by_clade: bool = False,
    clade_group_col: str = "clade_group",
    clade_group_values: Sequence[object] | None = None,
    window_stride: int = 2,
    run_composition: bool = True,
    run_mixing: bool = True,
) -> dict[str, Any]:
    """Run composition and/or mixing association models and export CSV tables."""
    output_path = (
        Path(output_dir)
        if output_dir
        else PROJECT_ROOT / "sse_detection" / "results" / "sse_outputs"
    )
    result_path = (
        Path(result_dir)
        if result_dir
        else PROJECT_ROOT / "sse_detection" / "results" / result_subdir
    )
    result_path.mkdir(parents=True, exist_ok=True)

    if group_by_clade and variant_adjuster == "clade":
        variant_adjuster = None

    if composition_model_sets is None:
        composition_model_sets = default_model_sets(
            variant_adjuster=variant_adjuster,
            window_adjustment=window_adjustment,
        )
    if mixing_model_sets is None:
        mixing_model_sets = default_model_sets(
            variant_adjuster=variant_adjuster,
            window_adjustment=window_adjustment,
        )

    frames = load_association_frames(
        output_dir=output_path,
        cluster_se=cluster_se,
        variant_adjuster=variant_adjuster,
        group_by_clade=group_by_clade,
        clade_group_col=clade_group_col,
        window_stride=window_stride,
        run_composition=run_composition,
    )

    group_col = clade_group_col if group_by_clade else None
    if clade_group_values is None:
        clade_group_values = [*CLADES.values(), "Other"]

    failures: list[dict[str, object]] = []
    composition_results: list[dict[str, pd.DataFrame]] = []
    mixing_results: list[dict[str, pd.DataFrame]] = []

    if run_composition:
        for model_set, adjusters in composition_model_sets.items():
            for group_label, source in _iter_groups(
                frames.composition_base,
                group_col,
                clade_group_values,
            ):
                composition_results.append(
                    run_composition_model_set(
                        source=source,
                        model_set=model_set,
                        adjusters=list(adjusters),
                        model_method=model_method,
                        cluster_se=cluster_se,
                        window_strata=window_strata,
                        group_col=group_col,
                        group_label=group_label,
                        failures=failures,
                    )
                )

    if run_mixing:
        for model_set, adjusters in mixing_model_sets.items():
            for group_label, source in _iter_groups(
                frames.node_model_base,
                group_col,
                clade_group_values,
            ):
                mixing_results.append(
                    run_mixing_model_set(
                        source=source,
                        model_set=model_set,
                        adjusters=list(adjusters),
                        mixing_features=mixing_features,
                        mixing_reference=mixing_reference,
                        model_method=model_method,
                        cluster_se=cluster_se,
                        window_strata=window_strata,
                        group_col=group_col,
                        group_label=group_label,
                        failures=failures,
                    )
                )

    summary_tables = _summarise_tables(
        composition_results=composition_results,
        mixing_results=mixing_results,
    )

    for filename, table in summary_tables.items():
        path = result_path / filename
        clean_export_table(table).to_csv(path, index=False)
        print(f"saved {filename}: {len(table):,} rows", flush=True)

    failures_df = pd.DataFrame(failures)
    if not failures_df.empty:
        clean_export_table(failures_df).to_csv(
            result_path / "model_failures.csv", index=False
        )
        print(f"saved model_failures.csv: {len(failures_df):,} rows", flush=True)

    frames.cluster_diagnostics.to_csv(
        result_path / "cluster_diagnostics.csv", index=False
    )

    return {
        "result_dir": result_path,
        "frames": frames,
        "summary_tables": summary_tables,
        "failures": failures_df,
        "cluster_diagnostics": frames.cluster_diagnostics,
    }


def run_policy_analysis(
    *,
    output_dir: Path | str | None = None,
    result_dir: Path | str | None = None,
    result_subdir: str = "policy_outputs",
    cluster_se: str = "cluster_id",
    window_strata: str = "window_idx",
    model_method: str = "firth_glm",
    window_stride: int = 2,
) -> dict[str, Any]:
    """Run policy-era association analyses."""
    output_path = (
        Path(output_dir)
        if output_dir
        else PROJECT_ROOT / "sse_detection" / "results" / "sse_outputs"
    )
    result_path = (
        Path(result_dir)
        if result_dir
        else PROJECT_ROOT / "sse_detection" / "results" / result_subdir
    )
    result_path.mkdir(parents=True, exist_ok=True)

    frames = load_association_frames(
        output_dir=output_path,
        cluster_se=cluster_se,
        variant_adjuster="clade",
        window_stride=window_stride,
        run_composition=False,
    )
    node_base = add_policy_era(add_clade_group(frames.node_model_base))
    print(
        "Policy eligible nodes: "
        f"{len(node_base):,}; candidates: {int(node_base['candidate'].sum()):,}",
        flush=True,
    )

    failures: list[dict[str, object]] = []
    policy_specs = [
        {
            "name": "policy_era",
            "column": "policy_era",
            "label": "Policy era",
            "reference": "post_restriction",
            "categorical": True,
        },
        {
            "name": "policy_intensity",
            "column": "policy_intensity",
            "label": "Policy intensity sensitivity",
            "reference": "per 1 policy-intensity point",
            "categorical": False,
        },
    ]
    policy_model_sets = {
        "primary": ["C(clade_group)", *WINDOW_SURVEILLANCE_ADJUSTERS],
        "expanded": [
            "C(clade_group)",
            *WINDOW_SURVEILLANCE_ADJUSTERS,
            *EXPANDED_CONTEXT_ADJUSTERS,
        ],
    }
    policy_results = [
        fit_single_exposure_specs(
            source=node_base,
            specs=policy_specs,
            domain="policy",
            model_set=model_set,
            adjusters=list(adjusters),
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
            required_base=["candidate", "cluster_id", cluster_se, window_strata],
            failures=failures,
            n_label="nodes",
        )
        for model_set, adjusters in policy_model_sets.items()
    ]
    policy_wald = _concat_or_empty([result["wald"] for result in policy_results])
    if not policy_wald.empty:
        policy_wald = bh_adjust_by(
            policy_wald,
            ["domain", "model_set", "predictor_set"],
        )

    summary_tables = {
        "policy_era_candidate_summary.csv": make_policy_era_candidate_summary(
            node_base
        ),
        "policy_era_category_summary.csv": make_policy_era_category_summary(
            node_base
        ),
        "policy_wald.csv": policy_wald,
        "policy_odds_ratios.csv": _concat_or_empty(
            [result["odds"] for result in policy_results]
        ),
        "policy_fit_stats.csv": _concat_or_empty(
            [result["fit_stats"] for result in policy_results]
        ),
    }
    _export_summary_tables(summary_tables, result_path)

    failures_df = pd.DataFrame(failures)
    if not failures_df.empty:
        clean_export_table(failures_df).to_csv(
            result_path / "model_failures.csv",
            index=False,
        )
        print(f"saved model_failures.csv: {len(failures_df):,} rows", flush=True)
    frames.cluster_diagnostics.to_csv(
        result_path / "cluster_diagnostics.csv",
        index=False,
    )
    return {
        "result_dir": result_path,
        "frames": frames,
        "policy_node_base": node_base,
        "summary_tables": summary_tables,
        "failures": failures_df,
        "cluster_diagnostics": frames.cluster_diagnostics,
    }


def run_vaccination_analysis(
    *,
    output_dir: Path | str | None = None,
    result_dir: Path | str | None = None,
    result_subdir: str = "vaccination_outputs",
    cluster_se: str = "cluster_id",
    window_strata: str = "window_idx",
    model_method: str = "firth_glm",
    window_stride: int = 2,
) -> dict[str, Any]:
    """Run vaccination-context association analyses."""
    output_path = (
        Path(output_dir)
        if output_dir
        else PROJECT_ROOT / "sse_detection" / "results" / "sse_outputs"
    )
    result_path = (
        Path(result_dir)
        if result_dir
        else PROJECT_ROOT / "sse_detection" / "results" / result_subdir
    )
    result_path.mkdir(parents=True, exist_ok=True)

    vaccination_columns = [
        "is_vaccinated",
        "vacc_dose_number",
        "vacc_booster",
        "days_since_vaccination",
        "dz_cum_prop_vaccinated",
    ]
    frames = load_association_frames(
        output_dir=output_path,
        cluster_se=cluster_se,
        variant_adjuster="clade",
        window_stride=window_stride,
        run_composition=True,
        extra_sequence_columns=vaccination_columns,
    )
    node_base = add_clade_group(frames.node_model_base)
    composition_base = add_vaccination_composition_features(frames.composition_base)
    vaccination_node_base = add_vaccination_node_features(
        node_base,
        composition_base,
        cluster_col="cluster_id",
    )
    vaccination_mixing_base = add_vaccination_mixing_features(
        vaccination_node_base,
        composition_base,
        cluster_col="cluster_id",
        window_col=window_strata,
        age_col="age_band",
    )
    print(
        "Vaccination eligible nodes: "
        f"{len(node_base):,}; candidates: {int(node_base['candidate'].sum()):,}",
        flush=True,
    )

    failures: list[dict[str, object]] = []
    vaccination_primary_adjusters = ["C(window_idx)", "C(clade)"]
    vaccination_composition_model_sets = {
        "primary": vaccination_primary_adjusters,
        "expanded": [
            *vaccination_primary_adjusters,
            *VACCINATION_COMPOSITION_EXPANDED_ADJUSTERS,
        ],
    }
    vaccination_node_model_sets = {
        "primary": vaccination_primary_adjusters,
        "expanded": [
            *vaccination_primary_adjusters,
            *VACCINATION_NODE_EXPANDED_ADJUSTERS,
        ],
    }

    vaccination_composition_results = []
    for model_set, adjusters in vaccination_composition_model_sets.items():
        single = fit_single_exposure_specs(
            source=composition_base,
            specs=VACCINATION_COMPOSITION_SPECS,
            domain="vaccination_composition",
            model_set=model_set,
            adjusters=list(adjusters),
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
            required_base=[
                "candidate",
                "cluster_id",
                "sequence_id",
                cluster_se,
                window_strata,
            ],
            failures=failures,
            n_label="rows",
        )
        joint = fit_joint_exposure_specs(
            source=composition_base,
            specs=VACCINATION_COMPOSITION_SPECS,
            joint_groups=VACCINATION_COMPOSITION_JOINT_GROUPS,
            domain="vaccination_composition",
            model_set=model_set,
            adjusters=list(adjusters),
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
            required_base=[
                "candidate",
                "cluster_id",
                "sequence_id",
                cluster_se,
                window_strata,
            ],
            failures=failures,
            n_label="rows",
        )
        vaccination_composition_results.append(
            {
                "wald": _concat_or_empty([single["wald"], joint["wald"]]),
                "odds": _concat_or_empty([single["odds"], joint["odds"]]),
                "fit_stats": _concat_or_empty(
                    [single["fit_stats"], joint["fit_stats"]]
                ),
            }
        )

    vaccination_node_specs = [
        {**spec, "categorical": False} for spec in VACCINATION_NODE_FEATURES
    ]
    vaccination_node_results = []
    for model_set, adjusters in vaccination_node_model_sets.items():
        single = fit_single_exposure_specs(
            source=vaccination_node_base,
            specs=vaccination_node_specs,
            domain="vaccination_node",
            model_set=model_set,
            adjusters=list(adjusters),
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
            required_base=["candidate", "cluster_id", cluster_se, window_strata],
            failures=failures,
            n_label="nodes",
        )
        joint = fit_joint_exposure_specs(
            source=vaccination_node_base,
            specs=vaccination_node_specs,
            joint_groups=VACCINATION_NODE_JOINT_GROUPS,
            domain="vaccination_node",
            model_set=model_set,
            adjusters=list(adjusters),
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
            required_base=["candidate", "cluster_id", cluster_se, window_strata],
            failures=failures,
            n_label="nodes",
        )
        vaccination_node_results.append(
            {
                "wald": _concat_or_empty([single["wald"], joint["wald"]]),
                "odds": _concat_or_empty([single["odds"], joint["odds"]]),
                "fit_stats": _concat_or_empty(
                    [single["fit_stats"], joint["fit_stats"]]
                ),
            }
        )

    vaccination_mixing_continuous_spec = [
        {
            "name": "vaccination_mix_entropy_z",
            "column": "vaccination_mix_entropy_z",
            "label": "Age-conditional vaccination mixing",
            "reference": "per 1 age-conditional null SD increase",
            "categorical": False,
        }
    ]
    vaccination_mixing_tertile_spec = [
        {
            "name": "vaccination_mix_tertile",
            "column": "vaccination_mix_tertile",
            "label": "Age-conditional vaccination mixing tertile",
            "reference": "as_expected",
            "categorical": True,
        }
    ]
    vaccination_mixing_model_sets = {
        "primary": vaccination_primary_adjusters,
        "age_mixing": [
            *vaccination_primary_adjusters,
            "age_entropy_z",
        ],
        "expanded": [
            *vaccination_primary_adjusters,
            "age_entropy_z",
            *VACCINATION_NODE_EXPANDED_ADJUSTERS,
        ],
    }
    vaccination_mixing_results = []
    for model_set, adjusters in vaccination_mixing_model_sets.items():
        continuous = fit_single_exposure_specs(
            source=vaccination_mixing_base,
            specs=vaccination_mixing_continuous_spec,
            domain="vaccination_mixing_age_conditional",
            model_set=model_set,
            adjusters=list(adjusters),
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
            required_base=["candidate", "cluster_id", cluster_se, window_strata],
            failures=failures,
            n_label="nodes",
        )
        tertile = fit_single_exposure_specs(
            source=vaccination_mixing_base,
            specs=vaccination_mixing_tertile_spec,
            domain="vaccination_mixing_age_conditional",
            model_set=model_set,
            adjusters=list(adjusters),
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
            required_base=["candidate", "cluster_id", cluster_se, window_strata],
            failures=failures,
            n_label="nodes",
        )
        vaccination_mixing_results.append(
            {
                "wald": _concat_or_empty([continuous["wald"], tertile["wald"]]),
                "odds": _concat_or_empty([continuous["odds"], tertile["odds"]]),
                "fit_stats": _concat_or_empty(
                    [continuous["fit_stats"], tertile["fit_stats"]]
                ),
            }
        )

    vaccination_composition_wald = _concat_or_empty(
        [result["wald"] for result in vaccination_composition_results]
    )
    if not vaccination_composition_wald.empty:
        vaccination_composition_wald = bh_adjust_by(
            vaccination_composition_wald,
            ["domain", "model_set", "predictor_set", "joint_model"],
        )
    vaccination_node_wald = _concat_or_empty(
        [result["wald"] for result in vaccination_node_results]
    )
    if not vaccination_node_wald.empty:
        vaccination_node_wald = bh_adjust_by(
            vaccination_node_wald,
            ["domain", "model_set", "predictor_set", "joint_model"],
        )
    vaccination_mixing_wald = _concat_or_empty(
        [result["wald"] for result in vaccination_mixing_results]
    )
    if not vaccination_mixing_wald.empty:
        vaccination_mixing_wald = bh_adjust_by(
            vaccination_mixing_wald,
            ["domain", "model_set", "predictor_set"],
        )

    summary_tables = {
        "vaccination_candidate_summary.csv": make_vaccination_candidate_summary(
            composition_base,
            vaccination_node_base,
        ),
        "vaccination_mixing_age_conditional_node_features.csv": (
            vaccination_mixing_node_feature_table(vaccination_mixing_base)
        ),
        "vaccination_mixing_age_conditional_summary.csv": (
            make_vaccination_mixing_age_conditional_summary(vaccination_mixing_base)
        ),
        "vaccination_mixing_age_conditional_category_summary.csv": (
            make_vaccination_mixing_age_conditional_category_summary(
                vaccination_mixing_base
            )
        ),
        "vaccination_composition_wald.csv": vaccination_composition_wald,
        "vaccination_composition_odds_ratios.csv": _concat_or_empty(
            [result["odds"] for result in vaccination_composition_results]
        ),
        "vaccination_composition_fit_stats.csv": _concat_or_empty(
            [result["fit_stats"] for result in vaccination_composition_results]
        ),
        "vaccination_node_wald.csv": vaccination_node_wald,
        "vaccination_node_odds_ratios.csv": _concat_or_empty(
            [result["odds"] for result in vaccination_node_results]
        ),
        "vaccination_node_fit_stats.csv": _concat_or_empty(
            [result["fit_stats"] for result in vaccination_node_results]
        ),
        "vaccination_mixing_age_conditional_wald.csv": vaccination_mixing_wald,
        "vaccination_mixing_age_conditional_odds_ratios.csv": _concat_or_empty(
            [result["odds"] for result in vaccination_mixing_results]
        ),
        "vaccination_mixing_age_conditional_fit_stats.csv": _concat_or_empty(
            [result["fit_stats"] for result in vaccination_mixing_results]
        ),
    }
    _export_summary_tables(summary_tables, result_path)

    failures_df = pd.DataFrame(failures)
    if not failures_df.empty:
        clean_export_table(failures_df).to_csv(
            result_path / "model_failures.csv",
            index=False,
        )
        print(f"saved model_failures.csv: {len(failures_df):,} rows", flush=True)
    frames.cluster_diagnostics.to_csv(
        result_path / "cluster_diagnostics.csv",
        index=False,
    )
    return {
        "result_dir": result_path,
        "frames": frames,
        "vaccination_composition_base": composition_base,
        "vaccination_node_base": vaccination_node_base,
        "vaccination_mixing_base": vaccination_mixing_base,
        "summary_tables": summary_tables,
        "failures": failures_df,
        "cluster_diagnostics": frames.cluster_diagnostics,
    }


def run_main_association_analysis(
    *,
    result_dir: Path | str | None = None,
    model_method: str = "firth_glm",
    variant_adjuster: str | None = "clade",
    window_adjustment: str = "fixed_effects",
    run_composition: bool = True,
    run_mixing: bool = True,
    **kwargs,
) -> dict[str, Any]:
    """Run the primary overall socio-geodemographic association analysis.

    This is a thin preset around :func:`run_association_pipeline` for the main
    notebook: no clade stratification, primary and expanded model sets,
    sequence-level composition models, and node-level entropy z-score mixing
    models saved under ``sse_detection/results/association_outputs`` by default.
    """
    return run_association_pipeline(
        result_dir=result_dir,
        result_subdir="association_outputs",
        model_method=model_method,
        variant_adjuster=variant_adjuster,
        window_adjustment=window_adjustment,
        group_by_clade=False,
        run_composition=run_composition,
        run_mixing=run_mixing,
        **kwargs,
    )
