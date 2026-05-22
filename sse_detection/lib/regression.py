"""Reusable regression helpers for SSE association analyses.

This module deliberately does not load data, filter rows, fill missing values,
or otherwise prepare analysis frames. Notebook code should do all data
preparation explicitly before calling these functions. Model fitting uses
``missing="raise"`` so accidental missing values fail loudly instead of being
silently dropped by patsy/statsmodels.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import chi2
import statsmodels.api as sm
import statsmodels.formula.api as smf


__all__ = [
    "AssociationModel",
    "bh_adjust",
    "bounded_exp",
    "categorical_term",
    "fit_binomial_glm",
    "fit_exposure_model",
    "make_formula",
    "model_fit_stats",
    "model_variables_from_terms",
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
        cov_kwds={"groups": data[cluster_col]},
    )


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
        out = out.loc[out["term"].str.contains(term_filter, regex=False)].copy()

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
    """Manual Wald test for all fitted parameters with a shared name prefix."""
    param_names = [p for p in result.params.index if p.startswith(term_prefix)]
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

    row = {
        "nobs": getattr(result, "nobs", np.nan),
        "df_model": getattr(result, "df_model", np.nan),
        "df_resid": getattr(result, "df_resid", np.nan),
        "log_likelihood": llf,
        "ll_null": llnull,
        "r2_mcfadden": r2_mcfadden,
        "aic": getattr(result, "aic", np.nan),
        "bic": getattr(result, "bic", np.nan),
        "bic_llf": getattr(result, "bic_llf", np.nan),
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
