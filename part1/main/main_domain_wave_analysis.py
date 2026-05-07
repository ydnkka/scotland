"""Main-formulation SIMD-domain and wave-specific Part 1 analyses.

This script extends ``main_analysis.py`` without returning to the older
all-resolution log-linear models. It uses the same primary Leiden resolution,
QC filter, lineage pooling, calendar spline, and window-clustered standard
errors as the main Part 1 analysis.

Outputs are written under ``part1/main``:

* SIMD-domain hurdle/ZTNB count models
* SIMD-domain hurdle/ZTNB count models with mixing predictors
* SIMD-domain quintile mixing models
* SIMD-domain demographic mixing models
* wave-specific SIMD-domain demographic mixing models
* primary-resolution observed-vs-expected mixing matrices
"""

from __future__ import annotations

import argparse
import gc
import math
import os
from pathlib import Path
from typing import Iterable
import warnings

import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning


from main_analysis import (  # noqa: E402
    CALENDAR_SPLINE_DF,
    COUNT_MODEL_SPECS,
    LINEAGE_MIN_CLUSTERS,
    PRIMARY_RESOLUTION,
    QC_DEFAULT,
    build_exog,
    expected_stratum_discordance,
    fit_ztnb,
    lineage_levels,
    logit_clipped,
    load_analysis_columns_pandas,
    load_simd_columns_pandas,
    observed_cluster_discordance,
    repo_root,
    zscore,
)


DOMAINS = {
    "overall": {
        "label": "Overall",
        "rank_col": "dz_simd_rank",
        "quintile_col": "dz_simd_quintile",
    },
    "income": {"label": "Income", "rank_col": "dz_simd_income_rank"},
    "employment": {"label": "Employment", "rank_col": "dz_simd_employment_rank"},
    "education": {"label": "Education", "rank_col": "dz_simd_education_rank"},
    "health": {"label": "Health", "rank_col": "dz_simd_health_rank"},
    "access": {"label": "Access", "rank_col": "dz_simd_access_rank"},
    "crime": {"label": "Crime", "rank_col": "dz_simd_crime_rank"},
    "housing": {"label": "Housing", "rank_col": "dz_simd_housing_rank"},
}

SHARED_COUNT_TERMS = [
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
]

SHARED_MIXING_TERMS = SHARED_COUNT_TERMS + ["log_cluster_size_z"]

DEMOGRAPHIC_MIXING = {
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
    "age_sex": {
        "column": "age_sex_profile",
        "label": "Joint age-sex profile mixing",
        "short_label": "Age-sex",
    },
}

DEMOGRAPHIC_MIXING_PREDICTOR_TERMS = [
    f"{prefix}_excess_mixing_z"
    for prefix in DEMOGRAPHIC_MIXING
]

WAVE_ORDER = [
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

WAVE_LABELS = {
    "B.1.177": "B.1.177",
    "Alpha": "Alpha",
    "Delta": "Delta",
    "BA.1": "BA.1",
    "BA.2": "BA.2",
    "BA.4": "BA.4",
    "BA.5": "BA.5",
    "BQ.1": "BQ.1",
    "XBB": "XBB",
    "Other": "Other",
}

MATRIX_VARIABLES = {
    "simd": {
        "column": "dz_simd_quintile",
        "label": "SIMD quintile",
        "levels": [1, 2, 3, 4, 5],
    },
    "age": {
        "column": "age_band",
        "label": "Age band",
        "levels": [
            "00-04",
            "05-09",
            "10-14",
            "15-19",
            "20-24",
            "25-29",
            "30-34",
            "35-39",
            "40-44",
            "45-49",
            "50-54",
            "55-59",
            "60-64",
            "65-69",
            "70-74",
            "75+",
        ],
    },
}

BASE_SEQUENCE_COLUMNS = [
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


def domain_rank_maxima() -> dict[str, float]:
    cols = [spec["rank_col"] for spec in DOMAINS.values()]
    simd = load_simd_columns_pandas(columns=cols)
    return {
        domain: float(simd[spec["rank_col"]].max())
        for domain, spec in DOMAINS.items()
    }


def rank_to_quintile(rank: pd.Series, max_rank: float) -> pd.Series:
    quintile = np.ceil(rank.astype(float) / (max_rank / 5.0))
    return quintile.clip(1, 5).astype("Int64").astype("category")


def assign_wave(lineage: str) -> str:
    if not isinstance(lineage, str):
        return "Other"
    if lineage.startswith("B.1.177"):
        return "B.1.177"
    if lineage == "B.1.1.7" or lineage.startswith("B.1.1.7."):
        return "Alpha"
    if lineage.startswith("AY.") or lineage == "B.1.617.2":
        return "Delta"
    if lineage.startswith("BA.1"):
        return "BA.1"
    if lineage.startswith("BA.2"):
        return "BA.2"
    if lineage.startswith("BA.4"):
        return "BA.4"
    if lineage.startswith("BA.5") or lineage.startswith("BE."):
        return "BA.5"
    if lineage.startswith("BQ."):
        return "BQ.1"
    if lineage.startswith("XBB"):
        return "XBB"
    return "Other"


def read_sequence_rows(
    _: Path,
    qc: str | None,
    primary_resolution: float,
) -> pd.DataFrame:
    rank_cols = [spec["rank_col"] for spec in DOMAINS.values()]
    columns = list(dict.fromkeys([*BASE_SEQUENCE_COLUMNS, *rank_cols]))
    seq = load_analysis_columns_pandas(
        columns=columns,
        resolution=primary_resolution,
        qc=qc,
    )

    maxima = domain_rank_maxima()
    for domain, spec in DOMAINS.items():
        q_col = f"{domain}_domain_quintile"
        if domain == "overall":
            seq[q_col] = seq["dz_simd_quintile"].astype("category")
        else:
            seq[q_col] = rank_to_quintile(seq[spec["rank_col"]], maxima[domain])

    categorical = [
        "cluster_id",
        "sequence_id",
        "window_id",
        "datazone",
        "pango_lineage",
        "nextclade_qc",
        "age_band",
        "sex",
        "dz_simd_quintile",
    ]
    for col in categorical:
        seq[col] = seq[col].astype("category")

    seq["collection_date"] = pd.to_datetime(seq["collection_date"])
    seq["wn_mid_date"] = pd.to_datetime(seq["wn_mid_date"])
    seq["wave_group"] = seq["pango_lineage"].astype(str).map(assign_wave).astype("category")
    complete_age_sex = seq[["age_band", "sex"]].notna().all(axis=1)
    seq["age_sex_profile"] = pd.NA
    seq.loc[complete_age_sex, "age_sex_profile"] = (
        seq.loc[complete_age_sex, "age_band"].astype(str)
        + "|"
        + seq.loc[complete_age_sex, "sex"].astype(str)
    )
    seq["age_sex_profile"] = seq["age_sex_profile"].astype("category")
    return seq


def build_cluster_table(
    seq: pd.DataFrame,
    lineage_min_clusters: int,
    calendar_spline_df: int,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    rank_cols = [spec["rank_col"] for spec in DOMAINS.values()]
    required = [
        "cluster_id",
        "sequence_id",
        "window_id",
        "window_idx",
        "collection_date",
        "datazone",
        "pango_lineage",
        "dz_cum_incidence_per_capita",
        "dz_cum_prop_sequenced",
        "wn_prop_sequenced",
        "dz_7d_test_positivity",
        *rank_cols,
    ]
    before = len(seq)
    seq = seq.dropna(subset=required).copy()
    dropped = before - len(seq)

    agg = {
        "cluster_size": ("sequence_id", "nunique"),
        "cluster_n_datazones": ("datazone", "nunique"),
        "cluster_start_date": ("collection_date", "min"),
        "cluster_end_date": ("collection_date", "max"),
        "resolution": ("resolution", "first"),
        "window_id": ("window_id", "first"),
        "window_idx": ("window_idx", "first"),
        "wn_mid_date": ("wn_mid_date", "first"),
        "pango_lineage": ("pango_lineage", "first"),
        "wave_group": ("wave_group", "first"),
        "mean_local_incidence_per_capita": ("dz_cum_incidence_per_capita", "mean"),
        "mean_local_seq_fraction": ("dz_cum_prop_sequenced", "mean"),
        "mean_window_seq_fraction": ("wn_prop_sequenced", "mean"),
        "mean_test_positivity": ("dz_7d_test_positivity", "mean"),
    }
    for domain, spec in DOMAINS.items():
        agg[f"{domain}_mean_rank"] = (spec["rank_col"], "mean")

    clusters = seq.groupby("cluster_id", observed=True, sort=False).agg(**agg).reset_index()
    clusters["duration_days"] = (
        clusters["cluster_end_date"] - clusters["cluster_start_date"]
    ).dt.days.astype(int)
    clusters["cluster_size_gt1"] = (clusters["cluster_size"] > 1).astype(int)
    clusters["duration_gt0"] = (clusters["duration_days"] > 0).astype(int)
    clusters["datazones_gt1"] = (clusters["cluster_n_datazones"] > 1).astype(int)
    clusters["cluster_size_excess"] = clusters["cluster_size"] - 1
    clusters["duration_positive_days"] = clusters["duration_days"]
    clusters["datazones_excess"] = clusters["cluster_n_datazones"] - 1

    for domain in DOMAINS:
        variable = f"{domain}_domain_quintile"
        prefix = f"{domain}_domain"
        obs = observed_cluster_discordance(seq, variable, prefix)
        exp = expected_stratum_discordance(seq, variable, prefix)
        clusters = clusters.merge(obs, on="cluster_id", how="left")
        clusters = clusters.merge(exp, on=["window_id", "pango_lineage"], how="left")
        clusters[f"{prefix}_excess_discordance"] = (
            clusters[f"{prefix}_discordance"] - clusters[f"{prefix}_expected_discordance"]
        )

    for prefix, spec in DEMOGRAPHIC_MIXING.items():
        obs = observed_cluster_discordance(seq, spec["column"], prefix)
        exp = expected_stratum_discordance(seq, spec["column"], prefix)
        clusters = clusters.merge(obs, on="cluster_id", how="left")
        clusters = clusters.merge(exp, on=["window_id", "pango_lineage"], how="left")
        clusters[f"{prefix}_excess_discordance"] = (
            clusters[f"{prefix}_discordance"] - clusters[f"{prefix}_expected_discordance"]
        )

    derived_cols = pd.DataFrame(
        {
            "local_incidence_log": np.log1p(
                clusters["mean_local_incidence_per_capita"].clip(lower=0) * 1000
            ),
            "local_seq_fraction_logit": logit_clipped(clusters["mean_local_seq_fraction"]),
            "window_seq_fraction_logit": logit_clipped(clusters["mean_window_seq_fraction"]),
            "test_positivity_logit": logit_clipped(
                clusters["mean_test_positivity"].fillna(0)
            ),
            "log_cluster_size": np.log(clusters["cluster_size"]),
        },
        index=clusters.index,
    )

    scaling_rows = []
    transforms = {
        "local_incidence_z": "local_incidence_log",
        "local_seq_fraction_z": "local_seq_fraction_logit",
        "window_seq_fraction_z": "window_seq_fraction_logit",
        "test_positivity_z": "test_positivity_logit",
        "log_cluster_size_z": "log_cluster_size",
    }
    transforms.update(
        {
            f"{domain}_domain_excess_mixing_z": f"{domain}_domain_excess_discordance"
            for domain in DOMAINS
        }
    )
    transforms.update(
        {
            f"{prefix}_excess_mixing_z": f"{prefix}_excess_discordance"
            for prefix in DEMOGRAPHIC_MIXING
        }
    )

    standardised_cols: dict[str, pd.Series] = {}
    for z_col, raw_col in transforms.items():
        source = derived_cols[raw_col] if raw_col in derived_cols else clusters[raw_col]
        standardised_cols[z_col], mean, sd = zscore(source)
        scaling_rows.append(
            {
                "standardised_column": z_col,
                "source_column": raw_col,
                "source_mean": mean,
                "source_sd": sd,
            }
        )

    for domain in DOMAINS:
        raw_col = f"{domain}_deprivation_raw"
        z_col = f"{domain}_deprivation_z"
        raw_values = -clusters[f"{domain}_mean_rank"]
        derived_cols[raw_col] = raw_values
        standardised_cols[z_col], mean, sd = zscore(raw_values)
        scaling_rows.append(
            {
                "standardised_column": z_col,
                "source_column": raw_col,
                "source_mean": mean,
                "source_sd": sd,
            }
        )

    lineage_counts = clusters["pango_lineage"].astype(str).value_counts()
    common_lineages = set(lineage_counts[lineage_counts >= lineage_min_clusters].index)
    lineage = clusters["pango_lineage"].astype(str)
    derived_cols["lineage_model"] = np.where(
        lineage.isin(common_lineages),
        lineage,
        "Other rare lineages",
    )

    clusters = pd.concat(
        [
            clusters.reset_index(drop=True),
            derived_cols.reset_index(drop=True),
            pd.DataFrame(standardised_cols).reset_index(drop=True),
        ],
        axis=1,
    ).copy()

    calendar = dmatrix(
        f"bs(window_idx, df={calendar_spline_df}, degree=3, include_intercept=False) - 1",
        clusters,
        return_type="dataframe",
    )
    calendar.columns = [f"calendar_spline_{i + 1}" for i in range(calendar.shape[1])]
    clusters = pd.concat([clusters.reset_index(drop=True), calendar.reset_index(drop=True)], axis=1)

    scaling = pd.DataFrame(scaling_rows)
    scaling.attrs["dropped_sequence_rows_missing_model_fields"] = dropped
    scaling.attrs["lineages_total"] = len(lineage_counts)
    scaling.attrs["lineages_modelled"] = len(common_lineages) + int(
        len(common_lineages) < len(lineage_counts)
    )
    scaling.attrs["lineage_min_clusters"] = lineage_min_clusters
    return clusters, scaling, dropped


def domain_term_label(domain: str) -> str:
    return f"{DOMAINS[domain]['label']} deprivation"


def domain_mixing_predictor_terms(domain: str) -> list[str]:
    return [f"{domain}_domain_excess_mixing_z", *DEMOGRAPHIC_MIXING_PREDICTOR_TERMS]


def term_label(domain: str, term: str) -> str:
    if term == f"{domain}_deprivation_z":
        return domain_term_label(domain)
    if term == f"{domain}_domain_excess_mixing_z":
        return f"{DOMAINS[domain]['label']} domain-quintile excess mixing"
    for prefix, spec in DEMOGRAPHIC_MIXING.items():
        if term == f"{prefix}_excess_mixing_z":
            return f"{spec['short_label']} excess mixing"
    return term


def extract_ratio_rows(
    *,
    params: np.ndarray,
    bse: np.ndarray,
    pvalues: np.ndarray,
    exog_names: list[str],
    terms: list[str],
    domain: str,
    outcome: str,
    outcome_label: str,
    component: str,
    component_label: str,
    model_family: str,
    response: str,
    n_observations: int,
    n_events: int | None,
) -> pd.DataFrame:
    idx = {name: i for i, name in enumerate(exog_names)}
    rows = []
    for term in terms:
        if term not in idx:
            continue
        i = idx[term]
        coef = float(params[i])
        stderr = float(bse[i])
        rows.append(
            {
                "domain": domain,
                "domain_label": DOMAINS[domain]["label"],
                "outcome": outcome,
                "outcome_label": outcome_label,
                "component": component,
                "component_label": component_label,
                "model_family": model_family,
                "response": response,
                "term": term,
                "term_label": term_label(domain, term),
                "coefficient": coef,
                "std_error_clustered_by_window": stderr,
                "z": coef / stderr if stderr > 0 else np.nan,
                "p_value": float(pvalues[i]),
                "ratio": float(np.exp(coef)),
                "ratio_ci_low": float(np.exp(coef - 1.96 * stderr)),
                "ratio_ci_high": float(np.exp(coef + 1.96 * stderr)),
                "n_observations": n_observations,
                "n_events": n_events,
            }
        )
    return pd.DataFrame(rows)


def fit_domain_binary_component(
    clusters: pd.DataFrame,
    spec,
    domain: str,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
    extra_terms: list[str] | None = None,
) -> tuple[pd.DataFrame, dict]:
    terms = [f"{domain}_deprivation_z", *SHARED_COUNT_TERMS]
    if extra_terms:
        terms.extend(extra_terms)
    use = clusters.dropna(subset=[spec.binary_col, *terms, *calendar_cols, "lineage_model"]).copy()
    y = use[spec.binary_col].astype(int)
    x = build_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use["window_id"].astype(str).to_numpy()
    model = sm.GLM(y, x, family=sm.families.Binomial())

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        result = model.fit(maxiter=maxiter, cov_type="cluster", cov_kwds={"groups": groups})

    rows = extract_ratio_rows(
        params=np.asarray(result.params, dtype=float),
        bse=np.asarray(result.bse, dtype=float),
        pvalues=np.asarray(result.pvalues, dtype=float),
        exog_names=list(result.model.exog_names),
        terms=terms,
        domain=domain,
        outcome=spec.name,
        outcome_label=spec.label,
        component="hurdle_binary",
        component_label="Probability of exceeding structural minimum",
        model_family="Binomial GLM with logit link",
        response=spec.binary_col,
        n_observations=int(len(use)),
        n_events=int(y.sum()),
    )
    diag = {
        "domain": domain,
        "domain_label": DOMAINS[domain]["label"],
        "outcome": spec.name,
        "component": "hurdle_binary",
        "model_family": "Binomial GLM with logit link",
        "response": spec.binary_col,
        "n_observations": int(len(use)),
        "n_events": int(y.sum()),
        "event_fraction": float(y.mean()),
        "n_features": int(x.shape[1]),
        "n_lineage_levels_available": int(len(lineage_levels_all)),
        "n_lineage_terms_used": int(sum(col.startswith("lineage_") for col in x.columns)),
        "n_windows": int(use["window_id"].nunique()),
        "converged": bool(getattr(result, "converged", False)),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "warnings": "; ".join(str(w.message) for w in caught),
    }
    return rows, diag


def fit_domain_positive_component(
    clusters: pd.DataFrame,
    spec,
    domain: str,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
    extra_terms: list[str] | None = None,
) -> tuple[pd.DataFrame, dict]:
    terms = [f"{domain}_deprivation_z", *SHARED_COUNT_TERMS]
    if extra_terms:
        terms.extend(extra_terms)
    use = clusters.loc[clusters[spec.positive_col] > 0].dropna(
        subset=[spec.positive_col, *terms, *calendar_cols, "lineage_model"]
    )
    use = use.copy()
    y = use[spec.positive_col].astype(int).to_numpy()
    x = build_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use["window_id"].astype(str).to_numpy()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fit_ztnb(y, x, groups, maxiter=maxiter)

    idx = {name: i for i, name in enumerate(result.exog_names)}
    pvalues = np.asarray(result.pvalues, dtype=float)
    rows = extract_ratio_rows(
        params=np.asarray(result.params, dtype=float),
        bse=np.asarray(result.bse, dtype=float),
        pvalues=pvalues,
        exog_names=result.exog_names,
        terms=terms,
        domain=domain,
        outcome=spec.name,
        outcome_label=spec.label,
        component="positive_zero_truncated_count",
        component_label=spec.positive_label,
        model_family="Zero-truncated negative binomial",
        response=spec.positive_col,
        n_observations=int(len(use)),
        n_events=None,
    )
    diag = {
        "domain": domain,
        "domain_label": DOMAINS[domain]["label"],
        "outcome": spec.name,
        "component": "positive_zero_truncated_count",
        "model_family": "Zero-truncated negative binomial",
        "response": spec.positive_col,
        "n_observations": int(len(use)),
        "n_events": None,
        "event_fraction": None,
        "mean_response": float(np.mean(y)),
        "max_response": int(np.max(y)),
        "n_features": int(x.shape[1]),
        "n_lineage_levels_available": int(len(lineage_levels_all)),
        "n_lineage_terms_used": int(sum(col.startswith("lineage_") for col in x.columns)),
        "n_windows": int(use["window_id"].nunique()),
        "converged": bool(result.converged),
        "iterations": int(result.nit),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "alpha": float(result.alpha),
        "alpha_at_upper_bound": bool(np.isclose(result.alpha, math.exp(8.0))),
        "optimizer_message": result.message,
        "warnings": "; ".join(str(w.message) for w in caught),
    }
    if f"{domain}_deprivation_z" in idx:
        diag["domain_deprivation_ratio"] = float(np.exp(result.params[idx[f"{domain}_deprivation_z"]]))
    return rows, diag


def fit_domain_count_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
    include_mixing_predictors: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary_specs = [spec for spec in COUNT_MODEL_SPECS if not spec.include_size]
    frames = []
    diagnostics = []
    predictor_set = "domain_primary_plus_mixing"
    for domain in DOMAINS:
        extra_terms = domain_mixing_predictor_terms(domain) if include_mixing_predictors else []
        for spec in primary_specs:
            if include_mixing_predictors and spec.name == "cluster_size":
                diagnostics.append(
                    {
                        "domain": domain,
                        "domain_label": DOMAINS[domain]["label"],
                        "outcome": spec.name,
                        "component": "hurdle_binary",
                        "model_family": None,
                        "response": spec.binary_col,
                        "skipped": True,
                        "reason": (
                            "mixing predictors require at least two valid cases, "
                            "so the cluster-size hurdle has no singleton comparison group"
                        ),
                        "predictor_set": predictor_set,
                        "extra_predictor_terms": ";".join(extra_terms),
                    }
                )
            else:
                suffix = " with mixing predictors" if include_mixing_predictors else ""
                print(f"  - {domain} {spec.name}: hurdle binary{suffix}", flush=True)
                rows, diag = fit_domain_binary_component(
                    clusters,
                    spec,
                    domain,
                    lineage_levels_all,
                    calendar_cols,
                    maxiter,
                    extra_terms=extra_terms,
                )
                if include_mixing_predictors:
                    rows = rows.copy()
                    rows["predictor_set"] = predictor_set
                    diag["predictor_set"] = predictor_set
                    diag["extra_predictor_terms"] = ";".join(extra_terms)
                frames.append(rows)
                diagnostics.append(diag)

            suffix = " with mixing predictors" if include_mixing_predictors else ""
            print(f"  - {domain} {spec.name}: zero-truncated NB positive count{suffix}", flush=True)
            rows, diag = fit_domain_positive_component(
                clusters,
                spec,
                domain,
                lineage_levels_all,
                calendar_cols,
                maxiter,
                extra_terms=extra_terms,
            )
            if include_mixing_predictors:
                rows = rows.copy()
                rows["predictor_set"] = predictor_set
                diag["predictor_set"] = predictor_set
                diag["extra_predictor_terms"] = ";".join(extra_terms)
            frames.append(rows)
            diagnostics.append(diag)
            gc.collect()
    return pd.concat(frames, ignore_index=True), pd.DataFrame(diagnostics)


def fit_linear_model(
    df: pd.DataFrame,
    *,
    domain: str,
    outcome: str,
    outcome_label: str,
    terms: list[str],
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    extra_result_fields: dict[str, object] | None = None,
    extra_diag_fields: dict[str, object] | None = None,
) -> tuple[pd.DataFrame, dict]:
    use = df.dropna(subset=[outcome, *terms, *calendar_cols, "lineage_model"]).copy()
    y = use[outcome].astype(float)
    x = build_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use["window_id"].astype(str).to_numpy()
    result = sm.OLS(y, x).fit(cov_type="cluster", cov_kwds={"groups": groups})

    names = list(result.model.exog_names)
    params = np.asarray(result.params, dtype=float)
    cov = np.asarray(result.cov_params(), dtype=float)
    variances = np.diag(cov)
    negative_variance = variances < -1e-12
    bse = np.full_like(variances, np.nan, dtype=float)
    valid_variance = ~negative_variance & np.isfinite(variances)
    bse[valid_variance] = np.sqrt(np.clip(variances[valid_variance], 0, None))
    idx = {name: i for i, name in enumerate(names)}
    negative_variance_terms = [
        name for name, invalid in zip(names, negative_variance) if invalid
    ]

    rows = []
    for term in terms:
        if term not in idx:
            continue
        i = idx[term]
        coef = float(params[i])
        stderr = float(bse[i])
        z_value = coef / stderr if np.isfinite(stderr) and stderr > 0 else np.nan
        p_value = float(2 * norm.sf(abs(z_value))) if np.isfinite(z_value) else np.nan
        row = {
            "domain": domain,
            "domain_label": DOMAINS[domain]["label"],
            "outcome": outcome,
            "outcome_label": outcome_label,
            "response": outcome,
            "term": term,
            "term_label": term_label(domain, term),
            "coefficient_excess_discordance": coef,
            "coefficient_percentage_points": coef * 100,
            "std_error_clustered_by_window": stderr,
            "std_error_note": (
                "negative clustered covariance diagonal"
                if term in negative_variance_terms
                else ""
            ),
            "z": z_value,
            "p_value": p_value,
            "ci_low": coef - 1.96 * stderr,
            "ci_high": coef + 1.96 * stderr,
            "ci_low_percentage_points": (coef - 1.96 * stderr) * 100,
            "ci_high_percentage_points": (coef + 1.96 * stderr) * 100,
            "n_observations": int(len(use)),
        }
        if extra_result_fields:
            row.update(extra_result_fields)
        rows.append(row)

    diag = {
        "domain": domain,
        "domain_label": DOMAINS[domain]["label"],
        "outcome": outcome,
        "outcome_label": outcome_label,
        "model_family": "Linear model for excess pairwise discordance",
        "response": outcome,
        "n_observations": int(len(use)),
        "n_features": int(x.shape[1]),
        "n_lineage_levels_available": int(len(lineage_levels_all)),
        "n_lineage_terms_used": int(sum(col.startswith("lineage_") for col in x.columns)),
        "n_windows": int(use["window_id"].nunique()),
        "mean_response": float(y.mean()),
        "sd_response": float(y.std(ddof=0)),
        "r2": float(result.rsquared),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "n_negative_cluster_variance_terms": int(negative_variance.sum()),
        "negative_cluster_variance_terms": ";".join(negative_variance_terms),
    }
    if extra_diag_fields:
        diag.update(extra_diag_fields)
    return pd.DataFrame(rows), diag


def fit_domain_quintile_mixing_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    diagnostics = []
    multi = clusters.loc[clusters["cluster_size"] >= 2].copy()
    for domain in DOMAINS:
        outcome = f"{domain}_domain_excess_discordance"
        terms = [f"{domain}_deprivation_z", *SHARED_MIXING_TERMS]
        rows, diag = fit_linear_model(
            multi,
            domain=domain,
            outcome=outcome,
            outcome_label=f"{DOMAINS[domain]['label']} quintile mixing",
            terms=terms,
            lineage_levels_all=lineage_levels_all,
            calendar_cols=calendar_cols,
            extra_result_fields={"mixing": "domain_quintile", "mixing_label": "Domain quintile mixing"},
            extra_diag_fields={"mixing": "domain_quintile", "mixing_label": "Domain quintile mixing"},
        )
        frames.append(rows)
        diagnostics.append(diag)
    return pd.concat(frames, ignore_index=True), pd.DataFrame(diagnostics)


def fit_domain_demographic_mixing_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    diagnostics = []
    multi = clusters.loc[clusters["cluster_size"] >= 2].copy()
    for domain in DOMAINS:
        for mixing, mixing_spec in DEMOGRAPHIC_MIXING.items():
            outcome = f"{mixing}_excess_discordance"
            terms = [f"{domain}_deprivation_z", *SHARED_MIXING_TERMS]
            rows, diag = fit_linear_model(
                multi,
                domain=domain,
                outcome=outcome,
                outcome_label=mixing_spec["label"],
                terms=terms,
                lineage_levels_all=lineage_levels_all,
                calendar_cols=calendar_cols,
                extra_result_fields={"mixing": mixing, "mixing_label": mixing_spec["label"]},
                extra_diag_fields={"mixing": mixing, "mixing_label": mixing_spec["label"]},
            )
            frames.append(rows)
            diagnostics.append(diag)
    return pd.concat(frames, ignore_index=True), pd.DataFrame(diagnostics)


def fit_wave_domain_demographic_mixing_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    *,
    min_clusters: int,
    min_windows: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    diagnostics = []
    multi = clusters.loc[
        (clusters["cluster_size"] >= 2) & (clusters["wave_group"].isin(WAVE_ORDER))
    ].copy()
    for wave in WAVE_ORDER:
        wave_df = multi.loc[multi["wave_group"] == wave].copy()
        n_windows = int(wave_df["window_id"].nunique())
        if len(wave_df) < min_clusters or n_windows < min_windows:
            diagnostics.append(
                {
                    "wave_group": wave,
                    "wave_label": WAVE_LABELS[wave],
                    "skipped": True,
                    "reason": "below minimum clusters/windows",
                    "n_clusters": int(len(wave_df)),
                    "n_windows": n_windows,
                }
            )
            continue

        for domain in DOMAINS:
            for mixing, mixing_spec in DEMOGRAPHIC_MIXING.items():
                outcome = f"{mixing}_excess_discordance"
                terms = [f"{domain}_deprivation_z", *SHARED_MIXING_TERMS]
                rows, diag = fit_linear_model(
                    wave_df,
                    domain=domain,
                    outcome=outcome,
                    outcome_label=mixing_spec["label"],
                    terms=terms,
                    lineage_levels_all=lineage_levels_all,
                    calendar_cols=calendar_cols,
                    extra_result_fields={
                        "wave_group": wave,
                        "wave_label": WAVE_LABELS[wave],
                        "mixing": mixing,
                        "mixing_label": mixing_spec["label"],
                    },
                    extra_diag_fields={
                        "wave_group": wave,
                        "wave_label": WAVE_LABELS[wave],
                        "mixing": mixing,
                        "mixing_label": mixing_spec["label"],
                        "skipped": False,
                        "reason": "",
                    },
                )
                frames.append(rows)
                diagnostics.append(diag)
    if not frames:
        return pd.DataFrame(), pd.DataFrame(diagnostics)
    return pd.concat(frames, ignore_index=True), pd.DataFrame(diagnostics)


def observed_ordered_pairs(cluster_counts: pd.DataFrame, levels: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide = (
        cluster_counts.pivot_table(
            index=["cluster_id", "wave_group", "window_id", "pango_lineage"],
            columns="category",
            values="n",
            aggfunc="sum",
            fill_value=0,
            observed=True,
        )
        .reindex(columns=levels, fill_value=0)
        .astype(np.int32)
        .reset_index()
    )
    wide["n_valid"] = wide[levels].sum(axis=1)
    wide = wide[wide["n_valid"] >= 2].copy()

    rows = []
    for left in levels:
        for right in levels:
            values = wide[left].astype(np.int64) * wide[right].astype(np.int64)
            if left == right:
                values = wide[left].astype(np.int64) * (wide[left].astype(np.int64) - 1)
            by_wave = values.groupby(wide["wave_group"], observed=True).sum()
            for wave, n_pairs in by_wave.items():
                rows.append(
                    {
                        "wave_group": wave,
                        "category_i": left,
                        "category_j": right,
                        "observed_pairs": float(n_pairs),
                    }
                )
    return pd.DataFrame(rows), wide


def expected_ordered_pairs(cluster_wide: pd.DataFrame, stratum_counts: pd.DataFrame, levels: list) -> pd.DataFrame:
    cluster_wide = cluster_wide.copy()
    cluster_wide["ordered_pairs"] = cluster_wide["n_valid"] * (cluster_wide["n_valid"] - 1)
    stratum_cols = ["wave_group", "window_id", "pango_lineage"]
    stratum_pair_totals = (
        cluster_wide.groupby(stratum_cols, observed=True)["ordered_pairs"]
        .sum()
        .rename("cluster_ordered_pairs")
        .reset_index()
    )

    stratum_wide = (
        stratum_counts.pivot_table(
            index=stratum_cols,
            columns="category",
            values="n",
            aggfunc="sum",
            fill_value=0,
            observed=True,
        )
        .reindex(columns=levels, fill_value=0)
        .astype(np.int64)
        .reset_index()
    )
    stratum_wide["stratum_n"] = stratum_wide[levels].sum(axis=1)
    stratum_wide = stratum_wide.merge(stratum_pair_totals, on=stratum_cols, how="inner")
    denom = stratum_wide["stratum_n"] * (stratum_wide["stratum_n"] - 1)

    rows = []
    for left in levels:
        for right in levels:
            numerator = stratum_wide[left].astype(np.float64) * stratum_wide[right].astype(np.float64)
            if left == right:
                numerator = stratum_wide[left].astype(np.float64) * (
                    stratum_wide[left].astype(np.float64) - 1
                )
            expected = stratum_wide["cluster_ordered_pairs"] * numerator / denom
            expected = expected.replace([np.inf, -np.inf], np.nan).fillna(0)
            by_wave = expected.groupby(stratum_wide["wave_group"], observed=True).sum()
            for wave, n_pairs in by_wave.items():
                rows.append(
                    {
                        "wave_group": wave,
                        "category_i": left,
                        "category_j": right,
                        "expected_pairs": float(n_pairs),
                    }
                )
    return pd.DataFrame(rows)


def build_matrix_for_variable(seq: pd.DataFrame, variable: str) -> pd.DataFrame:
    spec = MATRIX_VARIABLES[variable]
    levels = spec["levels"]
    work = seq.dropna(subset=[spec["column"]]).copy()
    work = work[work["wave_group"].isin(WAVE_ORDER)].copy()
    work["category"] = work[spec["column"]]

    cluster_counts = (
        work.groupby(
            ["cluster_id", "wave_group", "window_id", "pango_lineage", "category"],
            observed=True,
        )
        .size()
        .rename("n")
        .reset_index()
    )
    observed, cluster_wide = observed_ordered_pairs(cluster_counts, levels)

    stratum_counts = (
        work.groupby(["wave_group", "window_id", "pango_lineage", "category"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    expected = expected_ordered_pairs(cluster_wide, stratum_counts, levels)

    matrix = observed.merge(
        expected,
        on=["wave_group", "category_i", "category_j"],
        how="outer",
    ).fillna({"observed_pairs": 0, "expected_pairs": 0})
    matrix["variable"] = variable
    matrix["variable_label"] = spec["label"]

    overall = (
        matrix.groupby(["variable", "variable_label", "category_i", "category_j"], observed=True)[
            ["observed_pairs", "expected_pairs"]
        ]
        .sum()
        .reset_index()
    )
    overall["wave_group"] = "Overall"
    matrix = pd.concat([matrix, overall], ignore_index=True)

    totals = (
        matrix.groupby(["variable", "wave_group"], observed=True)[["observed_pairs", "expected_pairs"]]
        .sum()
        .rename(columns={"observed_pairs": "total_observed_pairs", "expected_pairs": "total_expected_pairs"})
        .reset_index()
    )
    matrix = matrix.merge(totals, on=["variable", "wave_group"], how="left")
    matrix["observed_probability"] = matrix["observed_pairs"] / matrix["total_observed_pairs"]
    matrix["expected_probability"] = matrix["expected_pairs"] / matrix["total_expected_pairs"]
    matrix["excess_probability"] = matrix["observed_probability"] - matrix["expected_probability"]
    matrix["excess_percentage_points"] = matrix["excess_probability"] * 100
    matrix["observed_expected_ratio"] = matrix["observed_probability"] / matrix["expected_probability"]
    matrix["wave_label"] = matrix["wave_group"].map(
        lambda w: "Overall" if w == "Overall" else WAVE_LABELS.get(w, w)
    )
    return matrix[
        [
            "variable",
            "variable_label",
            "wave_group",
            "wave_label",
            "category_i",
            "category_j",
            "observed_pairs",
            "expected_pairs",
            "observed_probability",
            "expected_probability",
            "excess_probability",
            "excess_percentage_points",
            "observed_expected_ratio",
        ]
    ]


def summarise_domain_wave_dataset(
    seq: pd.DataFrame,
    clusters: pd.DataFrame,
    qc: str | None,
    primary_resolution: float,
    dropped: int,
) -> pd.DataFrame:
    rows = [
        {"measure": "sequence_rows_used", "statistic": "count", "value": len(seq)},
        {
            "measure": "sequence_rows_dropped_missing_model_fields",
            "statistic": "count",
            "value": dropped,
        },
        {"measure": "clusters", "statistic": "count", "value": len(clusters)},
        {"measure": "clusters_size_ge_2", "statistic": "count", "value": int((clusters["cluster_size"] >= 2).sum())},
        {"measure": "primary_leiden_resolution", "statistic": "value", "value": primary_resolution},
        {"measure": "qc_filter", "statistic": "value", "value": qc or "none"},
    ]
    wave_counts = (
        clusters.loc[clusters["cluster_size"] >= 2, "wave_group"]
        .astype(str)
        .value_counts()
        .sort_index()
    )
    rows.extend(
        {"measure": f"non_singleton_clusters_wave_{wave}", "statistic": "count", "value": int(value)}
        for wave, value in wave_counts.items()
    )
    return pd.DataFrame(rows)


def run(
    root: Path,
    qc: str | None,
    primary_resolution: float,
    lineage_min_clusters: int,
    calendar_spline_df: int,
    maxiter: int,
    min_wave_clusters: int,
    min_wave_windows: int,
) -> None:
    out_dir = root / "part1" / "main"
    tables_dir = out_dir / "tables"
    cache_dir = out_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    print("Reading main-formulation domain/wave sequence rows", flush=True)
    seq = read_sequence_rows(root, qc=qc, primary_resolution=primary_resolution)
    print(f"Building domain/wave cluster table from {len(seq):,} sequence rows", flush=True)
    clusters, scaling, dropped = build_cluster_table(
        seq,
        lineage_min_clusters=lineage_min_clusters,
        calendar_spline_df=calendar_spline_df,
    )
    calendar_cols = [col for col in clusters.columns if col.startswith("calendar_spline_")]
    lineage_levels_all = lineage_levels(clusters)

    clusters.to_parquet(cache_dir / "main_domain_wave_cluster_table.parquet", index=False)
    scaling.to_csv(tables_dir / "main_domain_wave_covariate_scaling.csv", index=False)
    summarise_domain_wave_dataset(seq, clusters, qc, primary_resolution, dropped).to_csv(
        tables_dir / "main_domain_wave_dataset_descriptives.csv",
        index=False,
    )

    print("Fitting SIMD-domain hurdle/ZTNB count models", flush=True)
    count_results, count_diagnostics = fit_domain_count_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
        maxiter=maxiter,
    )
    count_results.to_csv(
        tables_dir / "main_simd_domain_hurdle_count_model_results.csv",
        index=False,
    )
    count_diagnostics.to_csv(
        tables_dir / "main_simd_domain_hurdle_count_model_diagnostics.csv",
        index=False,
    )

    print("Fitting SIMD-domain hurdle/ZTNB count models with mixing predictors", flush=True)
    mixing_predictor_count_results, mixing_predictor_count_diagnostics = fit_domain_count_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
        maxiter=maxiter,
        include_mixing_predictors=True,
    )
    mixing_predictor_count_results.to_csv(
        tables_dir / "main_simd_domain_mixing_predictor_hurdle_count_model_results.csv",
        index=False,
    )
    mixing_predictor_count_diagnostics.to_csv(
        tables_dir / "main_simd_domain_mixing_predictor_hurdle_count_model_diagnostics.csv",
        index=False,
    )

    print("Fitting SIMD-domain quintile mixing models", flush=True)
    domain_mixing, domain_mixing_diag = fit_domain_quintile_mixing_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
    )
    domain_mixing.to_csv(
        tables_dir / "main_simd_domain_quintile_mixing_model_results.csv",
        index=False,
    )
    domain_mixing_diag.to_csv(
        tables_dir / "main_simd_domain_quintile_mixing_model_diagnostics.csv",
        index=False,
    )

    print("Fitting SIMD-domain demographic mixing models", flush=True)
    domain_demo, domain_demo_diag = fit_domain_demographic_mixing_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
    )
    domain_demo.to_csv(
        tables_dir / "main_simd_domain_demographic_mixing_model_results.csv",
        index=False,
    )
    domain_demo_diag.to_csv(
        tables_dir / "main_simd_domain_demographic_mixing_model_diagnostics.csv",
        index=False,
    )

    print("Fitting wave-specific SIMD-domain demographic mixing models", flush=True)
    wave_demo, wave_demo_diag = fit_wave_domain_demographic_mixing_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
        min_clusters=min_wave_clusters,
        min_windows=min_wave_windows,
    )
    wave_demo.to_csv(
        tables_dir / "main_wave_specific_domain_demographic_mixing_model_results.csv",
        index=False,
    )
    wave_demo_diag.to_csv(
        tables_dir / "main_wave_specific_domain_demographic_mixing_model_diagnostics.csv",
        index=False,
    )

    print("Building primary-resolution observed-vs-expected mixing matrices", flush=True)
    matrices = pd.concat(
        [build_matrix_for_variable(seq, variable) for variable in MATRIX_VARIABLES],
        ignore_index=True,
    )
    matrices.to_csv(tables_dir / "main_observed_expected_mixing_matrices.csv", index=False)

    print(f"Wrote main domain/wave tables to {tables_dir}", flush=True)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument("--qc", default=QC_DEFAULT)
    parser.add_argument("--primary-resolution", type=float, default=PRIMARY_RESOLUTION)
    parser.add_argument("--lineage-min-clusters", type=int, default=LINEAGE_MIN_CLUSTERS)
    parser.add_argument("--calendar-spline-df", type=int, default=CALENDAR_SPLINE_DF)
    parser.add_argument("--maxiter", type=int, default=1000)
    parser.add_argument("--min-wave-clusters", type=int, default=1000)
    parser.add_argument("--min-wave-windows", type=int, default=4)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    qc = None if str(args.qc).lower() == "none" else str(args.qc)
    run(
        root=args.root.resolve(),
        qc=qc,
        primary_resolution=args.primary_resolution,
        lineage_min_clusters=args.lineage_min_clusters,
        calendar_spline_df=args.calendar_spline_df,
        maxiter=args.maxiter,
        min_wave_clusters=args.min_wave_clusters,
        min_wave_windows=args.min_wave_windows,
    )


if __name__ == "__main__":
    main()
