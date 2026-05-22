"""Reusable regression helpers for SSE association analyses.

This module deliberately does not load data, filter rows, fill missing values,
or otherwise prepare analysis frames. Notebook code should do all data
preparation explicitly before calling these functions. Model fitting uses
``missing="raise"`` so accidental missing values fail loudly instead of being
silently dropped by patsy/statsmodels.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
import re
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import patsy
from scipy.special import expit
from scipy.stats import chi2, norm
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.discrete.conditional_models import ConditionalLogit


__all__ = [
    "AssociationModel",
    "bh_adjust",
    "bounded_exp",
    "categorical_term",
    "cluster_se_diagnostics",
    "fit_binomial_glm",
    "fit_conditional_logit",
    "fit_exposure_model",
    "fit_firth_logit",
    "make_formula",
    "model_fit_stats",
    "model_variables_from_terms",
    "parameter_names_for_term",
    "robust_wald_for_params",
    "robust_wald_for_prefix",
    "tidy_odds_ratios",
    "tidy_single_parameter_wald",
]


@dataclass
class AssociationModel:
    """Container for one fitted association model and its tidy summaries."""

    result: object
    odds_ratios: pd.DataFrame
    wald: pd.DataFrame
    formula: str


@dataclass
class FirthLogitResult:
    """Minimal statsmodels-like result for Firth-penalised logistic fits."""

    params: pd.Series
    bse: pd.Series
    pvalues: pd.Series
    cov: pd.DataFrame
    model: object
    nobs: int
    df_model: int
    df_resid: int
    llf: float
    llnull: float
    aic: float
    bic_llf: float
    converged: bool
    fit_history: dict[str, object]

    def cov_params(self) -> pd.DataFrame:
        return self.cov

    def conf_int(self, alpha: float = 0.05) -> pd.DataFrame:
        z = norm.ppf(1 - alpha / 2)
        return pd.DataFrame(
            {
                0: self.params - z * self.bse,
                1: self.params + z * self.bse,
            },
            index=self.params.index,
        )


def categorical_term(variable: str, reference: str | int | float | None = None) -> str:
    """Return a patsy categorical term, optionally with a treatment reference."""
    if reference is None:
        return f"C({variable})"
    return f"C({variable}, Treatment(reference={reference!r}))"


def make_formula(
    outcome: str,
    exposure_term: str,
    adjusters: Sequence[str] | None = None,
) -> str:
    """Build a formula from an outcome, exposure term, and optional adjusters."""
    rhs = [exposure_term, *(adjusters or [])]
    return f"{outcome} ~ " + " + ".join(rhs)


def model_variables_from_terms(terms: Iterable[str]) -> list[str]:
    """Best-effort extraction of raw column names from simple formula terms.

    This is intended for notebook-side diagnostics and validation. It supports
    plain terms, ``C(col)``, ``C(col, ...)``, and ``np.log1p(col)``.
    """
    variables: set[str] = set()
    for term in terms:
        variables.update(re.findall(r"C\(([^,\)]+)", term))
        variables.update(re.findall(r"np\.log1p\(([^)]+)\)", term))
        if "(" not in term:
            variables.add(term)
    return sorted(variables)


def fit_binomial_glm(
    data: pd.DataFrame,
    formula: str,
    *,
    cluster_col: str | None = None,
):
    """Fit a binomial GLM, optionally with cluster-robust standard errors.

    Missing data are not handled here. ``missing="raise"`` forces the caller to
    prepare a complete model frame before fitting.
    """
    model = smf.glm(
        formula,
        data=data,
        family=sm.families.Binomial(),
        missing="raise",
    )

    if cluster_col is None:
        return model.fit()

    return model.fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data[cluster_col],
            "use_correction": True,
            "df_correction": True,
        },
    )


def _design_term_slices(result) -> dict[str, slice]:
    """Return Patsy term slices for supported statsmodels result objects."""
    custom = getattr(result, "_patsy_term_name_slices", None)
    if custom is not None:
        return dict(custom)

    model = getattr(result, "model", None)
    data = getattr(model, "data", None)
    design_info = getattr(data, "design_info", None)
    if design_info is None:
        return {}
    return dict(design_info.term_name_slices)


def _param_index(result) -> pd.Index:
    params = getattr(result, "params")
    if isinstance(params, pd.Series):
        return params.index
    names = getattr(getattr(result, "model", None), "exog_names", None)
    if names is None:
        names = [f"x{i}" for i in range(len(params))]
    return pd.Index(names)


def _matching_design_terms(term_name_slices: dict[str, slice], term: str) -> list[str]:
    if term in term_name_slices:
        return [term]

    # Backwards-compatible handling for older callers that passed "C(var"
    # after splitting a treatment-coded term at the comma. Match only the Patsy
    # factor token boundary so C(age, ...) cannot also collect C(age_band, ...).
    if term.startswith("C(") and not term.endswith(")"):
        return [
            name for name in term_name_slices
            if name.startswith(f"{term},") or name == f"{term})"
        ]

    return []


def parameter_names_for_term(result, term: str) -> list[str]:
    """Return fitted parameter names belonging to one full Patsy term.

    Prefer this over prefix or substring matching. It uses Patsy's
    ``term_name_slices`` when available, then falls back to exact parameter
    names and categorical-token boundaries.
    """
    params = _param_index(result)
    term_name_slices = _design_term_slices(result)
    matched_terms = _matching_design_terms(term_name_slices, term)

    if matched_terms:
        names: list[str] = []
        for matched in matched_terms:
            term_slice = term_name_slices[matched]
            names.extend(params[term_slice].tolist())
        return [name for name in names if name in params]

    if term in params:
        return [term]

    categorical_token = f"{term}["
    return [name for name in params if name.startswith(categorical_token)]


def fit_conditional_logit(
    data: pd.DataFrame,
    formula: str,
    *,
    strata_col: str,
    maxiter: int = 100,
    disp: bool = False,
):
    """Fit a conditional logistic model stratified by ``strata_col``.

    This is intended as the association-analysis primary model when window
    strata are nuisance parameters. Strata without outcome variation are
    dropped before fitting and reported on the result as
    ``_sse_diagnostics``.
    """
    y, x = patsy.dmatrices(
        formula,
        data=data,
        return_type="dataframe",
        NA_action="raise",
    )
    design_info = x.design_info
    original_columns = list(x.columns)
    if "Intercept" in x.columns:
        x = x.drop(columns="Intercept")
    term_name_slices = {}
    for term_name, term_slice in design_info.term_name_slices.items():
        columns = [
            col for col in original_columns[term_slice]
            if col in x.columns
        ]
        if not columns:
            continue
        positions = [x.columns.get_loc(col) for col in columns]
        term_name_slices[term_name] = slice(min(positions), max(positions) + 1)

    groups = data.loc[x.index, strata_col]
    outcome = y.iloc[:, 0]
    varying = outcome.groupby(groups, dropna=False).transform("nunique").gt(1)
    dropped_rows = int((~varying).sum())
    dropped_strata = int(groups.loc[~varying].nunique(dropna=False))

    y_fit = outcome.loc[varying]
    x_fit = x.loc[varying]
    groups_fit = groups.loc[varying]

    model = ConditionalLogit(
        y_fit,
        x_fit,
        groups=groups_fit,
        missing="raise",
    )
    result = model.fit(maxiter=maxiter, disp=disp)
    result._patsy_term_name_slices = term_name_slices
    result._sse_diagnostics = {
        "strata_col": strata_col,
        "input_rows": int(len(data)),
        "model_rows": int(len(y_fit)),
        "dropped_nonvarying_rows": dropped_rows,
        "dropped_nonvarying_strata": dropped_strata,
        "n_strata": int(groups_fit.nunique(dropna=False)),
    }
    return result


def _logistic_loglike(y: np.ndarray, x: np.ndarray, beta: np.ndarray) -> float:
    eta = np.clip(x @ beta, -35, 35)
    return float(np.sum(y * eta - np.logaddexp(0, eta)))


def _firth_penalized_loglike(
    y: np.ndarray,
    x: np.ndarray,
    beta: np.ndarray,
) -> float:
    eta = np.clip(x @ beta, -35, 35)
    mu = expit(eta)
    w = np.clip(mu * (1 - mu), 1e-9, None)
    info = (x.T * w) @ x
    sign, logdet = np.linalg.slogdet(info)
    if sign <= 0:
        return -np.inf
    return _logistic_loglike(y, x, beta) + 0.5 * float(logdet)


def fit_firth_logit(
    data: pd.DataFrame,
    formula: str,
    *,
    maxiter: int = 100,
    tol: float = 1e-6,
    max_step_halving: int = 25,
) -> FirthLogitResult:
    """Fit Firth's bias-reduced logistic regression for sparse binary data.

    The implementation follows the standard adjusted-score iteration for
    logistic regression. It is intended for the association notebook's sparse
    fixed-effect models when exact conditional logistic regression is
    computationally infeasible for very large window strata.
    """
    y_df, x_df = patsy.dmatrices(
        formula,
        data=data,
        return_type="dataframe",
        NA_action="raise",
    )
    y = y_df.iloc[:, 0].to_numpy(dtype=float)
    x = x_df.to_numpy(dtype=float)
    param_names = pd.Index(x_df.columns)
    beta = np.zeros(x.shape[1], dtype=float)
    converged = False
    pll = _firth_penalized_loglike(y, x, beta)

    for iteration in range(1, maxiter + 1):
        eta = np.clip(x @ beta, -35, 35)
        mu = expit(eta)
        w = np.clip(mu * (1 - mu), 1e-9, None)
        info = (x.T * w) @ x
        info_inv = np.linalg.pinv(info)
        leverage = w * np.einsum("ij,jk,ik->i", x, info_inv, x, optimize=True)
        adjusted_residual = y - mu + leverage * (0.5 - mu)
        score = x.T @ adjusted_residual
        step = info_inv @ score

        if not np.all(np.isfinite(step)):
            break

        step_scale = 1.0
        accepted = False
        for _ in range(max_step_halving + 1):
            candidate = beta + step_scale * step
            candidate_pll = _firth_penalized_loglike(y, x, candidate)
            if np.isfinite(candidate_pll) and candidate_pll >= pll:
                accepted = True
                break
            step_scale *= 0.5

        if not accepted:
            break

        delta = np.max(np.abs(candidate - beta))
        beta = candidate
        pll = candidate_pll
        if delta < tol:
            converged = True
            break

    eta = np.clip(x @ beta, -35, 35)
    mu = expit(eta)
    w = np.clip(mu * (1 - mu), 1e-9, None)
    info = (x.T * w) @ x
    cov = pd.DataFrame(np.linalg.pinv(info), index=param_names, columns=param_names)
    params = pd.Series(beta, index=param_names)
    bse = pd.Series(np.sqrt(np.diag(cov)), index=param_names)
    z_values = params / bse.replace(0, np.nan)
    pvalues = pd.Series(2 * norm.sf(np.abs(z_values)), index=param_names)
    llf = _logistic_loglike(y, x, beta)
    df_model = int(x.shape[1] - 1)
    nobs = int(x.shape[0])

    model = SimpleNamespace(
        data=SimpleNamespace(design_info=x_df.design_info),
        exog_names=list(param_names),
    )
    return FirthLogitResult(
        params=params,
        bse=bse,
        pvalues=pvalues,
        cov=cov,
        model=model,
        nobs=nobs,
        df_model=df_model,
        df_resid=int(nobs - x.shape[1]),
        llf=llf,
        llnull=np.nan,
        aic=float(-2 * llf + 2 * x.shape[1]),
        bic_llf=np.nan,
        converged=converged,
        fit_history={
            "iterations": iteration if "iteration" in locals() else 0,
            "penalized_loglike": pll,
        },
    )


def cluster_se_diagnostics(
    data: pd.DataFrame,
    cluster_col: str,
    *,
    outcome: str | None = None,
) -> pd.DataFrame:
    """Summarise cluster-robust SE support for an analysis frame."""
    if cluster_col not in data.columns:
        raise KeyError(f"{cluster_col!r} is not present in the analysis frame")

    clusters = data[cluster_col].dropna()
    row = {
        "cluster_col": cluster_col,
        "n_rows": int(len(data)),
        "n_clusters": int(clusters.nunique()),
        "min_rows_per_cluster": int(clusters.value_counts().min()) if not clusters.empty else 0,
        "median_rows_per_cluster": (
            float(clusters.value_counts().median()) if not clusters.empty else np.nan
        ),
    }

    if outcome is not None:
        if outcome not in data.columns:
            raise KeyError(f"{outcome!r} is not present in the analysis frame")
        by_cluster = data.groupby(cluster_col, dropna=False)[outcome]
        row["outcome_positive_clusters"] = int(by_cluster.max().fillna(0).gt(0).sum())
        row["outcome_varying_clusters"] = int(by_cluster.nunique(dropna=True).gt(1).sum())

    return pd.DataFrame([row])


def bounded_exp(values) -> np.ndarray:
    """Exponentiate log-odds while avoiding display-table overflow warnings."""
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    out[finite] = np.exp(np.clip(values[finite], -50, 50))
    return out


def tidy_odds_ratios(
    result,
    *,
    model_name: str | None = None,
    term_filter: str | None = None,
) -> pd.DataFrame:
    """Return coefficient, robust SE, p-value, and exponentiated OR columns."""
    ci = result.conf_int()
    out = pd.DataFrame(
        {
            "term": result.params.index,
            "estimate": result.params.to_numpy(),
            "std_error": result.bse.to_numpy(),
            "p_value": result.pvalues.to_numpy(),
            "conf_low": ci.iloc[:, 0].to_numpy(),
            "conf_high": ci.iloc[:, 1].to_numpy(),
        }
    )
    if model_name is not None:
        out.insert(0, "model", model_name)

    out["odds_ratio"] = bounded_exp(out["estimate"])
    out["or_low"] = bounded_exp(out["conf_low"])
    out["or_high"] = bounded_exp(out["conf_high"])

    if term_filter is not None:
        term_names = parameter_names_for_term(result, term_filter)
        out = out.loc[out["term"].isin(term_names)].copy()

    return out


def robust_wald_for_params(
    result,
    param_names: Sequence[str],
    *,
    model_name: str | None = None,
    term: str = "term",
) -> pd.DataFrame:
    """Manual Wald test for a selected set of fitted parameters.

    Uses a pseudoinverse and the covariance matrix rank, which makes the test
    more tolerant of robust covariance matrices that are singular or nearly
    singular. Non-finite estimates/covariances return ``NaN`` instead of
    raising inside statsmodels.
    """
    param_names = [p for p in param_names if p in result.params.index]
    if not param_names:
        row = {"term": term, "chi2": np.nan, "df": 0, "P>chi2": np.nan}
        if model_name is not None:
            row = {"model": model_name, **row}
        return pd.DataFrame([row])

    beta = result.params.loc[param_names].to_numpy(dtype=float)
    cov = result.cov_params().loc[param_names, param_names].to_numpy(dtype=float)
    finite = np.isfinite(beta) & np.isfinite(cov).all(axis=0) & np.isfinite(cov).all(axis=1)
    beta = beta[finite]
    cov = cov[np.ix_(finite, finite)]

    if beta.size == 0:
        chi2_stat = np.nan
        df = 0
        p_value = np.nan
    else:
        df = int(np.linalg.matrix_rank(cov))
        if df == 0:
            chi2_stat = np.nan
            p_value = np.nan
        else:
            chi2_stat = float(beta.T @ np.linalg.pinv(cov) @ beta)
            p_value = float(chi2.sf(chi2_stat, df))

    row = {"term": term, "chi2": chi2_stat, "df": df, "P>chi2": p_value}
    if model_name is not None:
        row = {"model": model_name, **row}
    return pd.DataFrame([row])


def robust_wald_for_prefix(
    result,
    term_prefix: str,
    *,
    model_name: str | None = None,
    term: str | None = None,
) -> pd.DataFrame:
    """Manual Wald test for all fitted parameters belonging to one term.

    The historical name is retained for notebook compatibility, but matching
    now uses Patsy factor tokens rather than raw string prefixes.
    """
    param_names = parameter_names_for_term(result, term_prefix)
    return robust_wald_for_params(
        result,
        param_names,
        model_name=model_name,
        term=term or term_prefix,
    )


def tidy_single_parameter_wald(
    result,
    terms: Sequence[str],
    *,
    model_name: str | None = None,
) -> pd.DataFrame:
    """Return one-row Wald tests for selected single-parameter terms."""
    rows = [
        robust_wald_for_params(
            result,
            [term],
            model_name=model_name,
            term=term,
        )
        for term in terms
    ]
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def fit_exposure_model(
    data: pd.DataFrame,
    *,
    outcome: str,
    exposure: str,
    adjusters: Sequence[str] | None = None,
    model_name: str | None = None,
    reference: str | int | float | None = None,
    cluster_col: str | None = None,
    categorical: bool = True,
) -> AssociationModel:
    """Fit one exposure model and return tidy OR and omnibus Wald summaries.

    The caller is responsible for complete-case filtering, type conversion,
    reference-level choice, and removal of non-identifiable strata.
    """
    exposure_term = (
        categorical_term(exposure, reference=reference)
        if categorical
        else exposure
    )
    formula = make_formula(outcome, exposure_term, adjusters)
    result = fit_binomial_glm(data, formula, cluster_col=cluster_col)

    odds = tidy_odds_ratios(
        result,
        model_name=model_name,
        term_filter=exposure_term.split(",")[0],
    )

    if categorical:
        wald = robust_wald_for_prefix(
            result,
            exposure_term.split(",")[0],
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

    return AssociationModel(
        result=result,
        odds_ratios=odds,
        wald=wald,
        formula=formula,
    )


def model_fit_stats(
    result,
    *,
    model_name: str | None = None,
    formula: str | None = None,
) -> pd.DataFrame:
    """Return available fit statistics for a fitted statsmodels result.

    Includes McFadden pseudo-R2 when both model and null log-likelihoods are
    available. Values unavailable on a particular statsmodels result are
    returned as ``NaN``.
    """
    llf = getattr(result, "llf", np.nan)
    llnull = getattr(result, "llnull", np.nan)
    if pd.notna(llf) and pd.notna(llnull) and llnull != 0:
        r2_mcfadden = 1 - (llf / llnull)
    else:
        r2_mcfadden = np.nan

    bic_llf = getattr(result, "bic_llf", np.nan)
    row = {
        "nobs": getattr(result, "nobs", np.nan),
        "df_model": getattr(result, "df_model", np.nan),
        "df_resid": getattr(result, "df_resid", np.nan),
        "log_likelihood": llf,
        "ll_null": llnull,
        "r2_mcfadden": r2_mcfadden,
        "aic": getattr(result, "aic", np.nan),
        "bic": bic_llf,
        "bic_llf": bic_llf,
        "converged": getattr(result, "converged", np.nan),
    }
    if model_name is not None:
        row = {"model": model_name, **row}
    if formula is not None:
        row["formula"] = formula
    return pd.DataFrame([row])


def bh_adjust(
    table: pd.DataFrame,
    p_col: str = "P>chi2",
    out_col: str = "p_adj_bh",
) -> pd.DataFrame:
    """Return a copy of ``table`` with Benjamini-Hochberg adjusted p-values.

    Missing p-values are preserved as missing and are not included in the
    multiple-testing correction.
    """
    from statsmodels.stats.multitest import multipletests

    out = table.copy()
    out[out_col] = np.nan
    valid = out[p_col].notna()
    if not valid.any():
        return out
    out.loc[valid, out_col] = multipletests(out.loc[valid, p_col], method="fdr_bh")[1]
    return out
