"""Low-level fitting primitives for Chapter 1.

* :func:`zscore`, :func:`logit_clipped` — scaling helpers.
* :func:`build_exog` — design-matrix assembly with QR-based redundant-column
  pruning.  Optionally takes a ``categorical`` block (for a wave or other
  factor) and an ``interaction_pairs`` list (e.g. ``[("age_excess_mixing_z",
  "wave_group"), ...]``) which generates the corresponding interaction
  columns.
* :func:`fit_ztnb` — zero-truncated negative binomial with analytical
  gradient and cluster-robust sandwich SE.
* :func:`clustered_logit_inference`, :func:`stable_binomial_fit_stats` —
  cluster-robust inference and stable log-likelihood for binomial GLMs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.linalg import pinvh, qr
from scipy.optimize import minimize
from scipy.special import digamma, expit, gammaln
from scipy.stats import norm
import statsmodels.api as sm


# ---------------------------------------------------------------------------
# Scaling helpers
# ---------------------------------------------------------------------------


def zscore(values: pd.Series) -> tuple[pd.Series, float, float]:
    mean = float(values.mean())
    sd = float(values.std(ddof=0))
    if not math.isfinite(sd) or sd == 0:
        raise ValueError(
            f"Cannot standardise {values.name!r}: zero or invalid SD."
        )
    return (values - mean) / sd, mean, sd


def logit_clipped(values: pd.Series, eps: float = 1e-5) -> pd.Series:
    clipped = values.clip(lower=eps, upper=1 - eps)
    return np.log(clipped / (1 - clipped))


# ---------------------------------------------------------------------------
# Design-matrix assembly
# ---------------------------------------------------------------------------


def drop_redundant_columns_qr(
    x: pd.DataFrame, tol: float = 1e-8, *, check_finite: bool = True,
) -> pd.DataFrame:
    """Drop numerically redundant columns via column-pivoted QR."""
    matrix = x.to_numpy(dtype=float)
    if check_finite and not np.all(np.isfinite(matrix)):
        for col in x.columns:
            if not np.all(np.isfinite(x[col].to_numpy(dtype=float))):
                raise ValueError(
                    f"Column {col!r} contains non-finite values."
                )
    if matrix.shape[0] == 0:
        return x
    R, P = qr(matrix, pivoting=True, mode="r")
    diag = np.abs(np.diag(R))
    if diag.size == 0:
        return x
    rank = int(np.sum(diag > tol * diag[0]))
    keep = np.sort(P[:rank])
    return x.iloc[:, keep]


def _categorical_dummies(
    series: pd.Series, prefix: str, reference: str | None = None,
) -> pd.DataFrame:
    """Reference-coded dummies for a categorical Series."""
    s = series.astype("category")
    levels = list(s.cat.categories)
    if reference is None:
        reference = levels[0]
    levels = [reference] + [lv for lv in levels if lv != reference]
    s = s.cat.reorder_categories(levels, ordered=False)
    dummies = pd.get_dummies(s, prefix=prefix, drop_first=True, dtype=float)
    dummies.index = series.index
    return dummies


def build_exog(
    df: pd.DataFrame,
    numeric_terms: Iterable[str],
    calendar_cols: Iterable[str] = (),
    lineage_levels: Iterable[str] | None = None,
    wave_reference: str | None = None,
    interaction_with_wave: Iterable[str] = (),
    *,
    tol: float = 1e-8,
) -> pd.DataFrame:
    """Assemble the design matrix.

    Layout
    ------
    1. Intercept.
    2. ``numeric_terms`` (mixing predictors + adjustments + log_size if used).
    3. Calendar spline columns.
    4. Either lineage dummies *or* wave-group dummies (mutually exclusive,
       since wave is a coarse grouping of lineage; the caller picks which).
    5. ``interaction_with_wave`` × wave dummies (only when ``wave_reference``
       is non-None and the corresponding numeric terms are listed).
    """
    numeric_terms = list(numeric_terms)
    calendar_cols = list(calendar_cols)
    interaction_with_wave = list(interaction_with_wave)

    parts: list[pd.DataFrame] = [
        pd.DataFrame({"const": np.ones(len(df), dtype=float)}, index=df.index),
        df[numeric_terms].astype(float),
    ]
    if calendar_cols:
        parts.append(df[calendar_cols].astype(float))

    wave_dummies: pd.DataFrame | None = None
    if wave_reference is not None:
        wave_dummies = _categorical_dummies(
            df["wave_group"], prefix="wave", reference=wave_reference,
        )
        parts.append(wave_dummies)
    elif lineage_levels is not None:
        lineage = pd.Categorical(
            df["lineage_model"].astype(str),
            categories=list(lineage_levels),
            ordered=False,
        )
        lineage_dummies = pd.get_dummies(
            pd.Series(lineage, index=df.index, name="lineage_model"),
            prefix="lineage",
            drop_first=True,
            dtype=float,
        )
        parts.append(lineage_dummies)

    if interaction_with_wave and wave_dummies is not None:
        for term in interaction_with_wave:
            if term not in df.columns:
                continue
            base = df[term].astype(float).values
            cols = {}
            for dummy_col in wave_dummies.columns:
                cols[f"{term}__x__{dummy_col}"] = (
                    base * wave_dummies[dummy_col].values
                )
            parts.append(pd.DataFrame(cols, index=df.index))

    x = pd.concat(parts, axis=1)
    return drop_redundant_columns_qr(x, tol=tol, check_finite=True)


def lineage_levels(clusters: pd.DataFrame) -> list[str]:
    """Lineage levels in descending count order."""
    counts = clusters["lineage_model"].astype(str).value_counts()
    return counts.index.tolist()


# ---------------------------------------------------------------------------
# Zero-truncated negative binomial
# ---------------------------------------------------------------------------


@dataclass
class ZTNBResult:
    params: np.ndarray
    bse: np.ndarray
    pvalues: np.ndarray
    exog_names: list[str]
    converged: bool
    message: str
    nit: int
    llf: float
    aic: float
    alpha: float


def ztnb_loglike_score(
    params: np.ndarray,
    y: np.ndarray,
    x: np.ndarray,
    offset: np.ndarray | None = None,
) -> tuple[float, np.ndarray, np.ndarray]:
    beta = params[:-1]
    log_alpha = float(np.clip(params[-1], -10.0, 8.0))
    alpha = float(np.exp(log_alpha))
    r = 1.0 / alpha

    eta = x @ beta
    if offset is not None:
        eta = eta + offset
    if not np.all(np.isfinite(eta)):
        return -np.inf, np.full_like(params, np.nan), np.full((len(y), len(params)), np.nan)
    eta = np.clip(eta, -30.0, 30.0)
    mu = np.exp(eta)

    log_r = math.log(r)
    log_r_mu = np.log(r + mu)
    log_p0 = r * (log_r - log_r_mu)
    p0 = np.exp(np.clip(log_p0, -745.0, -1e-12))
    one_minus_p0 = np.clip(1.0 - p0, 1e-12, 1.0)
    log_one_minus_p0 = np.log(one_minus_p0)

    logpmf = (
        gammaln(y + r) - gammaln(r) - gammaln(y + 1)
        + r * (log_r - log_r_mu)
        + y * (eta - log_r_mu)
    )
    ll_obs = logpmf - log_one_minus_p0
    llf = float(np.sum(ll_obs))

    p0_ratio = p0 / one_minus_p0
    score_eta = r * (y - mu) / (r + mu) - mu * r * p0_ratio / (r + mu)
    score_beta = x * score_eta[:, None]

    dlogpmf_dr = (
        digamma(y + r) - digamma(r) + log_r + 1.0
        - log_r_mu - (y + r) / (r + mu)
    )
    dlogp0_dr = log_r - log_r_mu + 1.0 - r / (r + mu)
    score_r = dlogpmf_dr + p0_ratio * dlogp0_dr
    score_log_alpha = -r * score_r

    score_obs = np.column_stack([score_beta, score_log_alpha])
    return llf, score_obs.sum(axis=0), score_obs


def _ztnb_objective(
    params: np.ndarray,
    y: np.ndarray,
    x: np.ndarray,
    offset: np.ndarray | None = None,
) -> tuple[float, np.ndarray]:
    llf, score, _ = ztnb_loglike_score(params, y, x, offset)
    if not np.isfinite(llf) or not np.all(np.isfinite(score)):
        return 1e100, np.zeros_like(params)
    return -llf, -score


def _ztnb_numerical_hessian(
    params: np.ndarray,
    y: np.ndarray,
    x: np.ndarray,
    offset: np.ndarray | None = None,
    eps: float = 1e-4,
) -> np.ndarray:
    n = len(params)
    _, grad0, _ = ztnb_loglike_score(params, y, x, offset)
    H = np.zeros((n, n))
    for i in range(n):
        p_fwd = params.copy()
        p_fwd[i] += eps
        _, grad_fwd, _ = ztnb_loglike_score(p_fwd, y, x, offset)
        H[i] = (grad_fwd - grad0) / eps
    return (H + H.T) / 2.0


def _ztnb_start_params(y: np.ndarray, x: pd.DataFrame) -> np.ndarray:
    beta = np.zeros(x.shape[1], dtype=float)
    beta[0] = float(np.log(np.mean(y)))
    try:
        poisson = sm.GLM(y, x, family=sm.families.Poisson()).fit(maxiter=50, disp=0)
        beta = np.asarray(poisson.params, dtype=float)
    except Exception:
        pass
    mean = float(np.mean(y))
    var = float(np.var(y, ddof=1))
    alpha = max((var - mean) / (mean * mean), 0.05)
    return np.r_[beta, math.log(alpha)]


def fit_ztnb(
    y: np.ndarray,
    x: pd.DataFrame,
    groups: np.ndarray,
    maxiter: int = 1000,
    offset: np.ndarray | None = None,
) -> ZTNBResult:
    """Fit ZTNB with cluster-robust SE clustered by ``groups``."""
    x_array = np.asarray(x, dtype=float)
    start = _ztnb_start_params(y, x)
    opt = minimize(
        _ztnb_objective,
        start,
        args=(y, x_array, offset),
        method="L-BFGS-B",
        jac=True,
        bounds=[(None, None)] * x_array.shape[1] + [(-10.0, 8.0)],
        options={"maxiter": maxiter, "ftol": 1e-8, "gtol": 1e-5, "maxls": 50},
    )
    params = opt.x.copy()
    params[-1] = float(np.clip(params[-1], -10.0, 8.0))
    llf, _, score_obs = ztnb_loglike_score(params, y, x_array, offset)

    H = _ztnb_numerical_hessian(params, y, x_array, offset)
    info = -H
    bread_inv = pinvh(info, rtol=1e-10)

    group_codes, inverse = np.unique(groups, return_inverse=True)
    cluster_scores = np.zeros((len(group_codes), len(params)), dtype=float)
    for g in range(len(group_codes)):
        cluster_scores[g, :] = score_obs[inverse == g].sum(axis=0)

    meat = cluster_scores.T @ cluster_scores
    cov = bread_inv @ meat @ bread_inv
    n = x_array.shape[0]
    k = len(params)
    if len(group_codes) > 1 and n > k:
        correction = (
            (len(group_codes) / (len(group_codes) - 1))
            * ((n - 1) / (n - k))
        )
        cov *= correction
    bse = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
    z = np.divide(params, bse, out=np.full_like(params, np.nan), where=bse > 0)
    pvalues = np.asarray(2 * norm.sf(np.abs(z)))

    return ZTNBResult(
        params=params,
        bse=bse,
        pvalues=pvalues,
        exog_names=list(x.columns) + ["log_alpha"],
        converged=bool(opt.success),
        message=str(opt.message),
        nit=int(opt.nit),
        llf=float(llf),
        aic=float(-2 * llf + 2 * len(params)),
        alpha=float(np.exp(params[-1])),
    )


# ---------------------------------------------------------------------------
# Cluster-robust binomial inference
# ---------------------------------------------------------------------------


def clustered_logit_inference(
    result, y: pd.Series, x: pd.DataFrame, groups: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    params = np.asarray(result.params, dtype=float)
    mu = expit(x_array @ params)
    mu = np.clip(mu, 1e-9, 1.0 - 1e-9)

    weights = mu * (1.0 - mu)
    info = (x_array * weights[:, None]).T @ x_array
    bread_inv = pinvh(info, rtol=1e-10)

    score_obs = x_array * (y_array - mu)[:, None]
    group_codes, inverse = np.unique(groups, return_inverse=True)
    cluster_scores = np.zeros((len(group_codes), x_array.shape[1]), dtype=float)
    for g in range(len(group_codes)):
        cluster_scores[g, :] = score_obs[inverse == g].sum(axis=0)

    meat = cluster_scores.T @ cluster_scores
    cov = bread_inv @ meat @ bread_inv
    n, p = x_array.shape
    if len(group_codes) > 1 and n > p:
        correction = (
            (len(group_codes) / (len(group_codes) - 1))
            * ((n - 1) / (n - p))
        )
        cov *= correction
    bse = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
    z = np.divide(params, bse, out=np.full_like(params, np.nan), where=bse > 0)
    pvalues = 2 * norm.sf(np.abs(z))
    return bse, pvalues


def stable_binomial_fit_stats(result, y: pd.Series, x: pd.DataFrame) -> tuple[float, float]:
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    params = np.asarray(result.params, dtype=float)
    mu = expit(x_array @ params)
    mu = np.clip(mu, 1e-12, 1.0 - 1e-12)
    llf = float(np.sum(y_array * np.log(mu) + (1.0 - y_array) * np.log1p(-mu)))
    aic = float(-2.0 * llf + 2.0 * len(params))
    return llf, aic
