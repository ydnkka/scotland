"""High-level fit routines for the Part 1 analysis.

Each function in this module orchestrates one family of fits: it takes a
prepared cluster table from :mod:`lib.data_prep`, calls the appropriate
low-level estimator from :mod:`lib.estimators`, and returns tidy
``(results, diagnostics)`` data frames.

The two lines of inquiry described in the manuscript map onto these
high-level functions as follows:

* **Line 1 — deprivation as exposure.**
    :func:`fit_count_models`, :func:`fit_mixing_models`,
    :func:`fit_domain_count_models`, :func:`fit_domain_quintile_mixing_models`,
    :func:`fit_domain_demographic_mixing_models`,
    :func:`fit_wave_domain_demographic_mixing_models`,
    :func:`fit_wave_outcome_models` (with ``extra_terms=None``),
    :func:`fit_loglinear_models` (with ``extra_terms=None``).

* **Line 2 — excess mixing as predictor.**
    :func:`fit_mixing_predictor_count_models`,
    :func:`fit_domain_count_models` with ``include_mixing_predictors=True``,
    :func:`fit_wave_outcome_models` with ``extra_terms``,
    :func:`fit_loglinear_models` with ``extra_terms``.
"""

from __future__ import annotations

import gc
import math
import warnings
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from .constants import (
    COUNT_MODEL_SPECS,
    CountModelSpec,
    DEMOGRAPHIC_MIXING,
    DEMOGRAPHIC_MIXING_PREDICTOR_TERMS,
    DOMAINS,
    MIXING_PREDICTOR_TERMS,
    MIXING_VARIABLES,
    PRIMARY_TERMS,
    SHARED_COUNT_TERMS,
    SHARED_MIXING_TERMS,
    TERM_LABELS,
    WAVE_LABELS,
    WAVE_ORDER,
)
from .estimators import (
    build_exog,
    build_wave_exog,
    clustered_logit_inference,
    fit_ztnb,
    lineage_levels,
    stable_binomial_fit_stats,
)


# ---------------------------------------------------------------------------
# Term and label helpers
# ---------------------------------------------------------------------------


def _model_terms(
    spec: CountModelSpec,
    primary_terms: Iterable[str] = PRIMARY_TERMS,
    extra_terms: Iterable[str] | None = None,
) -> list[str]:
    terms = list(primary_terms)
    if extra_terms:
        terms.extend(extra_terms)
    if spec.include_size:
        terms.append("log_cluster_size_z")
    return terms


def _domain_term_label(domain: str, term: str) -> str:
    if term == f"{domain}_deprivation_z":
        return f"{DOMAINS[domain]['label']} deprivation"
    if term == f"{domain}_domain_excess_mixing_z":
        return f"{DOMAINS[domain]['label']} domain-quintile excess mixing"
    for prefix, spec in DEMOGRAPHIC_MIXING.items():
        if term == f"{prefix}_excess_mixing_z":
            return f"{spec['short_label']} excess mixing"
    return TERM_LABELS.get(term, term)


def _domain_mixing_predictor_terms(domain: str) -> list[str]:
    return [f"{domain}_domain_excess_mixing_z", *DEMOGRAPHIC_MIXING_PREDICTOR_TERMS]


def _extract_ratio_rows_from_arrays(
    *,
    params: np.ndarray,
    bse: np.ndarray,
    pvalues: np.ndarray,
    exog_names: list[str],
    terms: list[str],
    base_fields: dict[str, object],
    label_resolver=lambda term: TERM_LABELS[term],
) -> pd.DataFrame:
    idx = {name: i for i, name in enumerate(exog_names)}
    rows: list[dict[str, object]] = []
    for term in terms:
        if term not in idx:
            continue
        i = idx[term]
        coef = float(params[i])
        stderr = float(bse[i])
        row = dict(base_fields)
        row.update(
            {
                "term": term,
                "term_label": label_resolver(term),
                "coefficient": coef,
                "std_error_clustered_by_window": stderr,
                "z": coef / stderr if stderr > 0 else np.nan,
                "p_value": float(pvalues[i]),
                "ratio": float(np.exp(coef)),
                "ratio_ci_low": float(np.exp(coef - 1.96 * stderr)),
                "ratio_ci_high": float(np.exp(coef + 1.96 * stderr)),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Hurdle (binary) and ZTNB (positive count) components
# ---------------------------------------------------------------------------


def fit_binary_component(
    clusters: pd.DataFrame,
    spec: CountModelSpec,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
    cluster_by: str = "window_id",
    primary_terms: Iterable[str] = PRIMARY_TERMS,
    extra_terms: list[str] | None = None,
    analysis_population_label: str | None = None,
) -> tuple[pd.DataFrame, dict]:
    terms = _model_terms(spec, primary_terms=primary_terms, extra_terms=extra_terms)
    if spec.include_size:
        use = clusters.loc[clusters["cluster_size"] > 1].copy()
    else:
        use = clusters.copy()
    use = use.dropna(subset=[spec.binary_col, *terms, *calendar_cols, "lineage_model"]).copy()
    y = use[spec.binary_col].astype(int)
    x = build_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use[cluster_by].astype(str).to_numpy()

    model = sm.GLM(y, x, family=sm.families.Binomial())
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        result = model.fit(
            maxiter=maxiter,
            cov_type="cluster",
            cov_kwds={"groups": groups},
        )

    rows = _extract_ratio_rows_from_arrays(
        params=np.asarray(result.params, dtype=float),
        bse=np.asarray(result.bse, dtype=float),
        pvalues=np.asarray(result.pvalues, dtype=float),
        exog_names=list(result.model.exog_names),
        terms=terms,
        base_fields={
            "outcome": spec.name,
            "outcome_label": spec.label,
            "component": "hurdle_binary",
            "component_label": "Probability of exceeding structural minimum",
            "model_family": "Binomial GLM with logit link",
            "response": spec.binary_col,
            "n_observations": len(use),
            "n_events": int(y.sum()),
        },
    )
    diag = {
        "outcome": spec.name,
        "component": "hurdle_binary",
        "model_family": "Binomial GLM with logit link",
        "response": spec.binary_col,
        "n_observations": int(len(use)),
        "n_events": int(y.sum()),
        "event_fraction": float(y.mean()),
        "analysis_population": (
            analysis_population_label
            or ("non-singleton clusters" if spec.include_size else "all primary-resolution clusters")
        ),
        "n_features": int(x.shape[1]),
        "n_lineage_levels_available": int(len(lineage_levels_all)),
        "n_lineage_terms_used": int(sum(col.startswith("lineage_") for col in x.columns)),
        "n_windows": int(use["window_id"].nunique()),
        "cluster_by": cluster_by,
        "converged": bool(getattr(result, "converged", False)),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "warnings": "; ".join(str(w.message) for w in caught),
    }
    return rows, diag


def fit_positive_component(
    clusters: pd.DataFrame,
    spec: CountModelSpec,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
    cluster_by: str = "window_id",
    use_size_offset: bool = False,
    winsorise_quantile: float = 0.0,
    primary_terms: Iterable[str] = PRIMARY_TERMS,
    extra_terms: list[str] | None = None,
) -> tuple[pd.DataFrame, dict]:
    terms = _model_terms(spec, primary_terms=primary_terms, extra_terms=extra_terms)
    use = clusters.loc[clusters[spec.positive_col] > 0].dropna(
        subset=[spec.positive_col, *terms, *calendar_cols, "lineage_model"]
    )
    use = use.copy()
    y = use[spec.positive_col].astype(int).to_numpy()

    # Tail winsorisation sensitivity: cap extreme positive counts at the
    # given quantile to check sensitivity to heavy-right-tail clusters.
    winsorised = False
    winsorise_cap: int | None = None
    if winsorise_quantile > 0.0:
        winsorise_cap = int(np.quantile(y, winsorise_quantile))
        y = np.minimum(y, winsorise_cap)
        winsorised = True

    # Window-pool size offset for cluster-size model (SAP §6.1).
    offset: np.ndarray | None = None
    if use_size_offset and spec.name == "cluster_size":
        wn_seq = use["wn_no_sequences"].to_numpy(dtype=float)
        offset = np.log(np.clip(wn_seq, 1.0, None))

    x = build_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use[cluster_by].astype(str).to_numpy()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fit_ztnb(y, x, groups, maxiter=maxiter, offset=offset)

    rows = _extract_ratio_rows_from_arrays(
        params=np.asarray(result.params, dtype=float),
        bse=np.asarray(result.bse, dtype=float),
        pvalues=np.asarray(result.pvalues, dtype=float),
        exog_names=result.exog_names,
        terms=terms,
        base_fields={
            "outcome": spec.name,
            "outcome_label": spec.label,
            "component": "positive_zero_truncated_count",
            "component_label": spec.positive_label,
            "model_family": "Zero-truncated negative binomial",
            "response": spec.positive_col,
            "n_observations": len(use),
            "n_events": None,
        },
    )
    diag = {
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
        "cluster_by": cluster_by,
        "size_offset_used": bool(offset is not None),
        "winsorised": winsorised,
        "winsorise_quantile": winsorise_quantile if winsorised else None,
        "winsorise_cap": winsorise_cap,
        "converged": bool(result.converged),
        "iterations": int(result.nit),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "alpha": float(result.alpha),
        "alpha_at_upper_bound": bool(np.isclose(result.alpha, math.exp(8.0))),
        "optimizer_message": result.message,
        "warnings": "; ".join(str(w.message) for w in caught),
    }
    return rows, diag


# ---------------------------------------------------------------------------
# Hurdle / ZTNB count models (Line 1: primary; Line 2: with mixing predictors)
# ---------------------------------------------------------------------------


def fit_count_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
    cluster_by: str = "window_id",
    use_size_offset: bool = False,
    winsorise_quantile: float = 0.0,
    primary_terms: Iterable[str] = PRIMARY_TERMS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Line 1 primary count models (hurdle + positive ZTNB) per spec."""
    result_frames: list[pd.DataFrame] = []
    diagnostics: list[dict] = []
    for spec in COUNT_MODEL_SPECS:
        if not spec.include_size:
            print(f"  - {spec.name}: hurdle binary", flush=True)
            rows, diag = fit_binary_component(
                clusters, spec, lineage_levels_all, calendar_cols, maxiter,
                cluster_by=cluster_by, primary_terms=primary_terms,
            )
            result_frames.append(rows)
            diagnostics.append(diag)

        print(f"  - {spec.name}: zero-truncated NB positive count", flush=True)
        rows, diag = fit_positive_component(
            clusters, spec, lineage_levels_all, calendar_cols, maxiter,
            cluster_by=cluster_by,
            use_size_offset=use_size_offset,
            winsorise_quantile=winsorise_quantile,
            primary_terms=primary_terms,
        )
        result_frames.append(rows)
        diagnostics.append(diag)
        gc.collect()

    return pd.concat(result_frames, ignore_index=True), pd.DataFrame(diagnostics)


def _skipped_mixing_predictor_diag(
    spec: CountModelSpec,
    component: str,
    response: str,
    reason: str,
) -> dict:
    return {
        "outcome": spec.name,
        "component": component,
        "model_family": None,
        "response": response,
        "skipped": True,
        "reason": reason,
        "analysis_population": "clusters with non-missing mixing predictors",
        "predictor_set": "primary_plus_mixing",
        "extra_predictor_terms": ";".join(MIXING_PREDICTOR_TERMS),
    }


def _add_predictor_set_metadata(
    rows: pd.DataFrame,
    diag: dict,
    *,
    predictor_set: str,
    extra_terms: list[str],
) -> tuple[pd.DataFrame, dict]:
    rows = rows.copy()
    if not rows.empty:
        rows["predictor_set"] = predictor_set
    diag = dict(diag)
    diag["predictor_set"] = predictor_set
    diag["extra_predictor_terms"] = ";".join(extra_terms)
    return rows, diag


def fit_mixing_predictor_count_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
    cluster_by: str = "window_id",
    use_size_offset: bool = False,
    winsorise_quantile: float = 0.0,
    primary_terms: Iterable[str] = PRIMARY_TERMS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Line 2 count models: excess-mixing predictors added to the primary set."""
    result_frames: list[pd.DataFrame] = []
    diagnostics: list[dict] = []
    predictor_set = "primary_plus_mixing"
    analysis_population = "clusters with non-missing mixing predictors"

    for spec in COUNT_MODEL_SPECS:
        if not spec.include_size:
            if spec.name == "cluster_size":
                diagnostics.append(
                    _skipped_mixing_predictor_diag(
                        spec,
                        "hurdle_binary",
                        spec.binary_col,
                        (
                            "mixing predictors require at least two valid cases, "
                            "so the cluster-size hurdle has no singleton comparison group"
                        ),
                    )
                )
            else:
                print(f"  - {spec.name}: hurdle binary with mixing predictors", flush=True)
                rows, diag = fit_binary_component(
                    clusters,
                    spec,
                    lineage_levels_all,
                    calendar_cols,
                    maxiter,
                    cluster_by=cluster_by,
                    primary_terms=primary_terms,
                    extra_terms=MIXING_PREDICTOR_TERMS,
                    analysis_population_label=analysis_population,
                )
                rows, diag = _add_predictor_set_metadata(
                    rows, diag,
                    predictor_set=predictor_set,
                    extra_terms=MIXING_PREDICTOR_TERMS,
                )
                result_frames.append(rows)
                diagnostics.append(diag)

        print(f"  - {spec.name}: zero-truncated NB positive count with mixing predictors", flush=True)
        rows, diag = fit_positive_component(
            clusters,
            spec,
            lineage_levels_all,
            calendar_cols,
            maxiter,
            cluster_by=cluster_by,
            use_size_offset=use_size_offset,
            winsorise_quantile=winsorise_quantile,
            primary_terms=primary_terms,
            extra_terms=MIXING_PREDICTOR_TERMS,
        )
        rows, diag = _add_predictor_set_metadata(
            rows, diag,
            predictor_set=predictor_set,
            extra_terms=MIXING_PREDICTOR_TERMS,
        )
        result_frames.append(rows)
        diagnostics.append(diag)
        gc.collect()

    return pd.concat(result_frames, ignore_index=True), pd.DataFrame(diagnostics)


# ---------------------------------------------------------------------------
# Linear models for excess pairwise discordance (Line 1 mixing outcomes)
# ---------------------------------------------------------------------------


def fit_mixing_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    cluster_by: str = "window_id",
    primary_terms: Iterable[str] = PRIMARY_TERMS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Line 1: mixing outcomes (SIMD, age, sex, joint profile) ~ deprivation + …"""
    result_frames: list[pd.DataFrame] = []
    diagnostics: list[dict] = []
    terms = list(primary_terms) + ["log_cluster_size_z"]
    for prefix, spec in MIXING_VARIABLES.items():
        outcome = f"{prefix}_excess_discordance"
        print(f"  - {prefix} mixing", flush=True)
        use = clusters.loc[clusters["cluster_size"] >= 2].dropna(
            subset=[outcome, *terms, *calendar_cols, "lineage_model"]
        )
        use = use.copy()
        y = use[outcome].astype(float)
        x = build_exog(use, terms, calendar_cols, lineage_levels_all)
        groups = use[cluster_by].astype(str).to_numpy()

        model = sm.OLS(y, x)
        result = model.fit(cov_type="cluster", cov_kwds={"groups": groups})
        names = list(result.model.exog_names)
        params = np.asarray(result.params, dtype=float)
        bse = np.asarray(result.bse, dtype=float)
        pvalues = np.asarray(result.pvalues, dtype=float)
        idx = {name: i for i, name in enumerate(names)}

        rows: list[dict[str, object]] = []
        for term in terms:
            i = idx[term]
            coef = float(params[i])
            stderr = float(bse[i])
            rows.append(
                {
                    "outcome": prefix,
                    "outcome_label": spec["label"],
                    "response": outcome,
                    "term": term,
                    "term_label": TERM_LABELS[term],
                    "coefficient_excess_discordance": coef,
                    "coefficient_percentage_points": coef * 100,
                    "std_error_clustered_by_window": stderr,
                    "z": coef / stderr if stderr > 0 else np.nan,
                    "p_value": float(pvalues[i]),
                    "ci_low": coef - 1.96 * stderr,
                    "ci_high": coef + 1.96 * stderr,
                    "ci_low_percentage_points": (coef - 1.96 * stderr) * 100,
                    "ci_high_percentage_points": (coef + 1.96 * stderr) * 100,
                    "n_observations": int(len(use)),
                }
            )
        result_frames.append(pd.DataFrame(rows))
        diagnostics.append(
            {
                "outcome": prefix,
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
            }
        )

    return pd.concat(result_frames, ignore_index=True), pd.DataFrame(diagnostics)


# ---------------------------------------------------------------------------
# Log-linear sensitivity
# ---------------------------------------------------------------------------


_LOGLINEAR_OUTCOMES: dict[str, dict[str, object]] = {
    "cluster_size": {
        "label": "Cluster size",
        "source": "cluster_size",
        "log_plus": 0,
    },
    "geographic_dispersion": {
        "label": "Geographic dispersion",
        "source": "cluster_n_datazones",
        "log_plus": 0,
    },
}


def fit_loglinear_models(
    clusters: pd.DataFrame,
    *,
    extra_terms: list[str] | None = None,
    predictor_set: str | None = None,
    primary_terms: Iterable[str] = PRIMARY_TERMS,
) -> pd.DataFrame:
    """Log-linear (single-component) sensitivity vs the hurdle/ZTNB split."""
    calendar_cols = [col for col in clusters.columns if col.startswith("calendar_spline_")]
    lineage_levels_all = lineage_levels(clusters)
    frames: list[pd.DataFrame] = []
    for outcome, spec in _LOGLINEAR_OUTCOMES.items():
        terms = list(primary_terms)
        if extra_terms:
            terms.extend(extra_terms)
        use = clusters.dropna(
            subset=[spec["source"], *terms, *calendar_cols, "lineage_model"]
        ).copy()
        y_raw = use[spec["source"]].astype(float) + float(spec["log_plus"])
        y = np.log(y_raw)
        x = build_exog(use, terms, calendar_cols, lineage_levels_all)
        groups = use["window_id"].astype(str).to_numpy()
        result = sm.OLS(y, x).fit(cov_type="cluster", cov_kwds={"groups": groups})

        names = list(result.model.exog_names)
        params = np.asarray(result.params, dtype=float)
        bse = np.asarray(result.bse, dtype=float)
        pvalues = np.asarray(result.pvalues, dtype=float)
        idx = {name: i for i, name in enumerate(names)}

        rows: list[dict[str, object]] = []
        for term in terms:
            i = idx[term]
            coef = float(params[i])
            stderr = float(bse[i])
            row = {
                "model": outcome,
                "model_label": spec["label"],
                "outcome": spec["source"],
                "term": term,
                "term_label": TERM_LABELS[term],
                "coefficient_log_ratio": coef,
                "std_error_clustered_by_window": stderr,
                "z": coef / stderr if stderr > 0 else np.nan,
                "p_value": float(pvalues[i]),
                "geometric_mean_ratio": float(np.exp(coef)),
                "ci_low": float(np.exp(coef - 1.96 * stderr)),
                "ci_high": float(np.exp(coef + 1.96 * stderr)),
                "n_observations": int(len(use)),
                "r2": float(result.rsquared),
            }
            if predictor_set is not None:
                row["predictor_set"] = predictor_set
                row["extra_predictor_terms"] = ";".join(extra_terms or [])
            rows.append(row)
        frames.append(pd.DataFrame(rows))
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# SIMD-domain count and mixing models
# ---------------------------------------------------------------------------


def _domain_extract_rows(
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
    base_fields = {
        "domain": domain,
        "domain_label": DOMAINS[domain]["label"],
        "outcome": outcome,
        "outcome_label": outcome_label,
        "component": component,
        "component_label": component_label,
        "model_family": model_family,
        "response": response,
        "n_observations": n_observations,
        "n_events": n_events,
    }
    return _extract_ratio_rows_from_arrays(
        params=params,
        bse=bse,
        pvalues=pvalues,
        exog_names=exog_names,
        terms=terms,
        base_fields=base_fields,
        label_resolver=lambda term: _domain_term_label(domain, term),
    )


def fit_domain_binary_component(
    clusters: pd.DataFrame,
    spec: CountModelSpec,
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

    rows = _domain_extract_rows(
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
    spec: CountModelSpec,
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
    rows = _domain_extract_rows(
        params=np.asarray(result.params, dtype=float),
        bse=np.asarray(result.bse, dtype=float),
        pvalues=np.asarray(result.pvalues, dtype=float),
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
        diag["domain_deprivation_ratio"] = float(
            np.exp(result.params[idx[f"{domain}_deprivation_z"]])
        )
    return rows, diag


def fit_domain_count_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
    include_mixing_predictors: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-SIMD-domain hurdle/ZTNB count models for the supplements."""
    primary_specs = [spec for spec in COUNT_MODEL_SPECS if not spec.include_size]
    frames: list[pd.DataFrame] = []
    diagnostics: list[dict] = []
    predictor_set = "domain_primary_plus_mixing"
    for domain in DOMAINS:
        extra_terms = _domain_mixing_predictor_terms(domain) if include_mixing_predictors else []
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
                    clusters, spec, domain,
                    lineage_levels_all, calendar_cols, maxiter,
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
                clusters, spec, domain,
                lineage_levels_all, calendar_cols, maxiter,
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


# ---------------------------------------------------------------------------
# Linear mixing models for the SIMD-domain supplements
# ---------------------------------------------------------------------------


def _fit_linear_model(
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

    rows: list[dict[str, object]] = []
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
            "term_label": _domain_term_label(domain, term),
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
    frames: list[pd.DataFrame] = []
    diagnostics: list[dict] = []
    multi = clusters.loc[clusters["cluster_size"] >= 2].copy()
    for domain in DOMAINS:
        outcome = f"{domain}_domain_excess_discordance"
        terms = [f"{domain}_deprivation_z", *SHARED_MIXING_TERMS]
        rows, diag = _fit_linear_model(
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
    frames: list[pd.DataFrame] = []
    diagnostics: list[dict] = []
    multi = clusters.loc[clusters["cluster_size"] >= 2].copy()
    for domain in DOMAINS:
        for mixing, mixing_spec in DEMOGRAPHIC_MIXING.items():
            outcome = f"{mixing}_excess_discordance"
            terms = [f"{domain}_deprivation_z", *SHARED_MIXING_TERMS]
            rows, diag = _fit_linear_model(
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
    frames: list[pd.DataFrame] = []
    diagnostics: list[dict] = []
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
                rows, diag = _fit_linear_model(
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


# ---------------------------------------------------------------------------
# Wave-stratified hurdle / ZTNB count models
# ---------------------------------------------------------------------------


def _wave_extract_rows(
    *,
    params: np.ndarray,
    bse: np.ndarray,
    pvalues: np.ndarray,
    exog_names: list[str],
    terms: list[str],
    wave: str,
    spec: CountModelSpec,
    component: str,
    component_label: str,
    model_family: str,
    response: str,
    n_observations: int,
    n_events: int | None,
) -> pd.DataFrame:
    base_fields = {
        "wave_group": wave,
        "wave_label": WAVE_LABELS.get(wave, wave),
        "outcome": spec.name,
        "outcome_label": spec.label,
        "component": component,
        "component_label": component_label,
        "model_family": model_family,
        "response": response,
        "n_observations": n_observations,
        "n_events": n_events,
    }
    return _extract_ratio_rows_from_arrays(
        params=params,
        bse=bse,
        pvalues=pvalues,
        exog_names=exog_names,
        terms=terms,
        base_fields=base_fields,
        label_resolver=lambda term: TERM_LABELS[term],
    )


def _skipped_wave_diag(
    wave: str,
    spec: CountModelSpec,
    component: str,
    reason: str,
    n_observations: int,
    n_windows: int,
) -> dict:
    return {
        "wave_group": wave,
        "wave_label": WAVE_LABELS.get(wave, wave),
        "outcome": spec.name,
        "component": component,
        "skipped": True,
        "reason": reason,
        "n_observations": int(n_observations),
        "n_windows": int(n_windows),
    }


def fit_wave_binary_component(
    wave_df: pd.DataFrame,
    wave: str,
    spec: CountModelSpec,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    *,
    maxiter: int,
    min_events: int,
    extra_terms: list[str] | None = None,
) -> tuple[pd.DataFrame, dict]:
    terms = list(PRIMARY_TERMS)
    if extra_terms:
        terms.extend(extra_terms)
    use = wave_df.dropna(subset=[spec.binary_col, *terms, *calendar_cols, "lineage_model"]).copy()
    y = use[spec.binary_col].astype(int)
    n_events = int(y.sum())
    n_nonevents = int(len(y) - n_events)
    n_windows = int(use["window_id"].nunique())
    if n_events < min_events or n_nonevents < min_events:
        return pd.DataFrame(), _skipped_wave_diag(
            wave, spec, "hurdle_binary",
            "below minimum events/non-events",
            len(use), n_windows,
        )

    x = build_wave_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use["window_id"].astype(str).to_numpy()
    model = sm.GLM(y, x, family=sm.families.Binomial())
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        try:
            result = model.fit(maxiter=maxiter)
        except Exception as exc:
            return pd.DataFrame(), _skipped_wave_diag(
                wave, spec, "hurdle_binary",
                f"fit failed: {exc}",
                len(use), n_windows,
            )

    bse, pvalues = clustered_logit_inference(result, y, x, groups)
    log_likelihood, aic = stable_binomial_fit_stats(result, y, x)

    rows = _wave_extract_rows(
        params=np.asarray(result.params, dtype=float),
        bse=bse,
        pvalues=pvalues,
        exog_names=list(result.model.exog_names),
        terms=terms,
        wave=wave,
        spec=spec,
        component="hurdle_binary",
        component_label="Probability of exceeding structural minimum",
        model_family="Binomial GLM with logit link",
        response=spec.binary_col,
        n_observations=int(len(use)),
        n_events=n_events,
    )
    diag = {
        "wave_group": wave,
        "wave_label": WAVE_LABELS.get(wave, wave),
        "outcome": spec.name,
        "component": "hurdle_binary",
        "skipped": False,
        "reason": "",
        "model_family": "Binomial GLM with logit link",
        "response": spec.binary_col,
        "n_observations": int(len(use)),
        "n_events": n_events,
        "event_fraction": float(y.mean()),
        "n_features": int(x.shape[1]),
        "n_lineage_levels_available": int(len(lineage_levels_all)),
        "n_lineage_terms_used": int(sum(col.startswith("lineage_") for col in x.columns)),
        "lineage_adjustment": "wave-stratified; lineage dummies included, rank-dropped if collinear",
        "n_windows": n_windows,
        "covariance_method": "window-clustered sandwich with pseudo-inverse bread",
        "converged": bool(getattr(result, "converged", False)),
        "log_likelihood": log_likelihood,
        "aic": aic,
        "warnings": "; ".join(str(w.message) for w in caught),
    }
    return rows, diag


def fit_wave_positive_component(
    wave_df: pd.DataFrame,
    wave: str,
    spec: CountModelSpec,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    *,
    maxiter: int,
    min_positive: int,
    min_windows: int,
    extra_terms: list[str] | None = None,
) -> tuple[pd.DataFrame, dict]:
    terms = list(PRIMARY_TERMS)
    if extra_terms:
        terms.extend(extra_terms)
    use = wave_df.loc[wave_df[spec.positive_col] > 0].dropna(
        subset=[spec.positive_col, *terms, *calendar_cols, "lineage_model"]
    )
    use = use.copy()
    n_windows = int(use["window_id"].nunique())
    if len(use) < min_positive or n_windows < min_windows:
        return pd.DataFrame(), _skipped_wave_diag(
            wave, spec, "positive_zero_truncated_count",
            "below minimum positive clusters/windows",
            len(use), n_windows,
        )

    y = use[spec.positive_col].astype(int).to_numpy()
    x = build_wave_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use["window_id"].astype(str).to_numpy()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fit_ztnb(y, x, groups, maxiter=maxiter)

    rows = _wave_extract_rows(
        params=np.asarray(result.params, dtype=float),
        bse=np.asarray(result.bse, dtype=float),
        pvalues=np.asarray(result.pvalues, dtype=float),
        exog_names=result.exog_names,
        terms=terms,
        wave=wave,
        spec=spec,
        component="positive_zero_truncated_count",
        component_label=spec.positive_label,
        model_family="Zero-truncated negative binomial",
        response=spec.positive_col,
        n_observations=int(len(use)),
        n_events=None,
    )
    diag = {
        "wave_group": wave,
        "wave_label": WAVE_LABELS.get(wave, wave),
        "outcome": spec.name,
        "component": "positive_zero_truncated_count",
        "skipped": False,
        "reason": "",
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
        "lineage_adjustment": "wave-stratified; lineage dummies included, rank-dropped if collinear",
        "n_windows": n_windows,
        "converged": bool(result.converged),
        "iterations": int(result.nit),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "alpha": float(result.alpha),
        "alpha_at_upper_bound": bool(np.isclose(result.alpha, math.exp(8.0))),
        "optimizer_message": result.message,
        "warnings": "; ".join(str(w.message) for w in caught),
    }
    return rows, diag


def fit_wave_outcome_models(
    clusters: pd.DataFrame,
    *,
    maxiter: int,
    min_clusters: int,
    min_windows: int,
    min_positive: int,
    min_events: int,
    extra_terms: list[str] | None = None,
    predictor_set: str | None = None,
    skip_cluster_size_binary: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Wave-stratified hurdle/ZTNB count models."""
    primary_specs = [spec for spec in COUNT_MODEL_SPECS if not spec.include_size]
    calendar_cols = [col for col in clusters.columns if col.startswith("calendar_spline_")]
    frames: list[pd.DataFrame] = []
    diagnostics: list[dict] = []

    for wave in WAVE_ORDER:
        wave_df = clusters.loc[clusters["wave_group"] == wave].copy()
        n_windows = int(wave_df["window_id"].nunique())
        if len(wave_df) < min_clusters or n_windows < min_windows:
            for spec in primary_specs:
                diag = _skipped_wave_diag(
                    wave, spec, "all_components",
                    "below minimum clusters/windows",
                    len(wave_df), n_windows,
                )
                if predictor_set is not None:
                    diag["predictor_set"] = predictor_set
                    diag["extra_predictor_terms"] = ";".join(extra_terms or [])
                diagnostics.append(diag)
            continue

        lineage_levels_wave = lineage_levels(wave_df)
        for spec in primary_specs:
            if skip_cluster_size_binary and spec.name == "cluster_size":
                diag = _skipped_wave_diag(
                    wave, spec, "hurdle_binary",
                    (
                        "mixing predictors require at least two valid cases, "
                        "so the cluster-size hurdle has no singleton comparison group"
                    ),
                    len(wave_df.dropna(subset=[*MIXING_PREDICTOR_TERMS])),
                    n_windows,
                )
                if predictor_set is not None:
                    diag["predictor_set"] = predictor_set
                    diag["extra_predictor_terms"] = ";".join(extra_terms or [])
                diagnostics.append(diag)
            else:
                print(f"  - {wave} {spec.name}: hurdle binary", flush=True)
                rows, diag = fit_wave_binary_component(
                    wave_df, wave, spec,
                    lineage_levels_wave, calendar_cols,
                    maxiter=maxiter,
                    min_events=min_events,
                    extra_terms=extra_terms,
                )
                if predictor_set is not None:
                    if not rows.empty:
                        rows = rows.copy()
                        rows["predictor_set"] = predictor_set
                    diag["predictor_set"] = predictor_set
                    diag["extra_predictor_terms"] = ";".join(extra_terms or [])
                if not rows.empty:
                    frames.append(rows)
                diagnostics.append(diag)

            print(f"  - {wave} {spec.name}: zero-truncated NB positive count", flush=True)
            rows, diag = fit_wave_positive_component(
                wave_df, wave, spec,
                lineage_levels_wave, calendar_cols,
                maxiter=maxiter,
                min_positive=min_positive,
                min_windows=min_windows,
                extra_terms=extra_terms,
            )
            if predictor_set is not None:
                if not rows.empty:
                    rows = rows.copy()
                    rows["predictor_set"] = predictor_set
                diag["predictor_set"] = predictor_set
                diag["extra_predictor_terms"] = ";".join(extra_terms or [])
            if not rows.empty:
                frames.append(rows)
            diagnostics.append(diag)
            gc.collect()

    results = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return results, pd.DataFrame(diagnostics)


def summarise_wave_outcomes(clusters: pd.DataFrame) -> pd.DataFrame:
    """Per-wave descriptive outcome statistics."""
    rows: list[dict[str, object]] = []
    outcome_specs = [
        ("cluster_size", "Cluster size", "cluster_size_gt1", "cluster_size_excess"),
        ("duration_days", "Duration", "duration_gt0", "duration_positive_days"),
        ("cluster_n_datazones", "Geographic spread", "datazones_gt1", "datazones_excess"),
    ]
    for wave in WAVE_ORDER:
        sub = clusters.loc[clusters["wave_group"] == wave]
        if sub.empty:
            continue
        for raw_col, label, binary_col, positive_col in outcome_specs:
            values = sub[raw_col].dropna()
            pos = sub.loc[sub[positive_col] > 0, positive_col].dropna()
            rows.append(
                {
                    "wave_group": wave,
                    "wave_label": WAVE_LABELS.get(wave, wave),
                    "outcome": raw_col,
                    "outcome_label": label,
                    "n_clusters": int(len(sub)),
                    "n_windows": int(sub["window_id"].nunique()),
                    "structural_minimum_fraction": float(1 - sub[binary_col].mean()),
                    "exceeds_minimum_fraction": float(sub[binary_col].mean()),
                    "mean": float(values.mean()),
                    "median": float(values.median()),
                    "p75": float(values.quantile(0.75)),
                    "p90": float(values.quantile(0.90)),
                    "p95": float(values.quantile(0.95)),
                    "positive_n": int(len(pos)),
                    "positive_mean": float(pos.mean()) if len(pos) else np.nan,
                    "positive_median": float(pos.median()) if len(pos) else np.nan,
                    "positive_p90": float(pos.quantile(0.90)) if len(pos) else np.nan,
                }
            )
    return pd.DataFrame(rows)
