"""Low-level fitting primitives for the Part 1 analysis.

These functions are the numerical core of the analysis: design-matrix
construction, the zero-truncated negative binomial (ZTNB) implementation
(with analytical gradient, numerical Hessian, and a cluster-robust sandwich
covariance), and replacement inference for cluster-robust binomial GLMs in
the wave-stratified models.

The higher-level orchestration that decides *which* models to fit, on *which*
subset of the cluster table, and *which* coefficients to extract lives in
``fit_models``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.linalg import pinvh
from scipy.optimize import minimize
from scipy.special import digamma, expit, gammaln
from scipy.stats import norm
import statsmodels.api as sm


# ---------------------------------------------------------------------------
# Scaling helpers
# ---------------------------------------------------------------------------


def zscore(values: pd.Series) -> tuple[pd.Series, float, float]:
    """Return ``(z, mean, sd)`` for ``values``; SD is the population SD."""
    mean = float(values.mean())
    sd = float(values.std(ddof=0))
    if not math.isfinite(sd) or sd == 0:
        raise ValueError(f"Cannot standardise {values.name!r}: zero or invalid SD.")
    return (values - mean) / sd, mean, sd


def logit_clipped(values: pd.Series, eps: float = 1e-5) -> pd.Series:
    clipped = values.clip(lower=eps, upper=1 - eps)
    return np.log(clipped / (1 - clipped))


# ---------------------------------------------------------------------------
# Design-matrix construction
# ---------------------------------------------------------------------------


def build_exog(
    df: pd.DataFrame,
    numeric_terms: Iterable[str],
    calendar_cols: Iterable[str],
    all_lineage_levels: Iterable[str],
) -> pd.DataFrame:
    """Assemble the model design matrix.

    Columns are ordered: intercept → numeric covariates → calendar spline →
    lineage dummies (one-hot encoded with the reference level dropped).
    Empty lineage columns (i.e. levels with no observations in ``df``) are
    pruned so the matrix stays full rank.
    """
    numeric_terms = list(numeric_terms)
    calendar_cols = list(calendar_cols)
    all_lineage_levels = list(all_lineage_levels)

    parts: list[pd.DataFrame] = [
        pd.DataFrame({"const": np.ones(len(df), dtype=float)}, index=df.index),
        df[numeric_terms].astype(float),
        df[calendar_cols].astype(float),
    ]

    lineages = pd.Categorical(
        df["lineage_model"].astype(str),
        categories=all_lineage_levels,
        ordered=False,
    )
    lineage_dummies = pd.get_dummies(
        pd.Series(lineages, index=df.index, name="lineage_model"),
        prefix="lineage",
        drop_first=True,
        dtype=float,
    )
    parts.append(lineage_dummies)

    x = pd.concat(parts, axis=1)
    zero_columns = [
        col
        for col in x.columns
        if col not in {"const", *numeric_terms} and float(x[col].abs().sum()) == 0
    ]
    if zero_columns:
        x = x.drop(columns=zero_columns)
    return x


def drop_redundant_columns(x: pd.DataFrame, tol: float = 1e-8) -> pd.DataFrame:
    """Keep a full-rank column set, preserving earlier columns first.

    The intended use is the wave-stratified design matrix, where sparse
    lineages may collide with the intercept or with calendar splines.  The
    column ordering produced by :func:`build_wave_exog` is intercept →
    primary terms → calendar spline → lineage dummies, so the sequential
    orthogonalisation always preserves the substantive covariates unless
    they are themselves redundant within a wave.
    """
    keep: list[str] = []
    basis: list[np.ndarray] = []
    n = max(len(x), 1)
    base_scale = math.sqrt(n)
    for col in x.columns:
        values = np.asarray(x[col], dtype=float)
        if not np.all(np.isfinite(values)):
            continue
        col_norm = float(np.linalg.norm(values))
        if col_norm <= tol:
            continue
        residual = values.copy()
        for q in basis:
            residual -= float(np.dot(q, residual)) * q
        residual_norm = float(np.linalg.norm(residual))
        if residual_norm > tol * max(col_norm, base_scale):
            keep.append(col)
            basis.append(residual / residual_norm)
    return x[keep]


def build_wave_exog(
    df: pd.DataFrame,
    terms: Iterable[str],
    calendar_cols: Iterable[str],
    lineage_levels_wave: Iterable[str],
) -> pd.DataFrame:
    """Wave-stratified design matrix.

    Lineage dummies are included when ``lineage_levels_wave`` is non-empty.
    ``drop_redundant_columns`` handles rank deficiency caused by sparse
    lineages within a wave, so no manual lineage pre-filtering is needed.
    """
    terms = list(terms)
    calendar_cols = list(calendar_cols)
    lineage_levels_wave = list(lineage_levels_wave)

    parts = [
        pd.DataFrame({"const": np.ones(len(df), dtype=float)}, index=df.index),
        df[terms].astype(float),
        df[calendar_cols].astype(float),
    ]
    if lineage_levels_wave:
        lineages = pd.Categorical(
            df["lineage_model"].astype(str),
            categories=lineage_levels_wave,
            ordered=False,
        )
        lineage_dummies = pd.get_dummies(
            pd.Series(lineages, index=df.index, name="lineage_model"),
            prefix="lineage",
            drop_first=True,
            dtype=float,
        )
        parts.append(lineage_dummies)
    return drop_redundant_columns(pd.concat(parts, axis=1))


def lineage_levels(clusters: pd.DataFrame) -> list[str]:
    """Lineage levels in the order they appear (descending count) in ``clusters``."""
    counts = clusters["lineage_model"].astype(str).value_counts()
    return counts.index.tolist()


# ---------------------------------------------------------------------------
# Zero-truncated negative binomial implementation
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
    """ZTNB log-likelihood, summed score, and observation-level scores."""
    beta = params[:-1]
    log_alpha = float(np.clip(params[-1], -10.0, 8.0))
    alpha = float(np.exp(log_alpha))
    r = 1.0 / alpha

    eta = x @ beta
    if offset is not None:
        eta = eta + offset
    if not np.all(np.isfinite(eta)):
        bad = np.full_like(params, np.nan)
        return -np.inf, bad, np.full((len(y), len(params)), np.nan)
    eta = np.clip(eta, -30.0, 30.0)
    mu = np.exp(eta)

    log_r = math.log(r)
    log_r_mu = np.log(r + mu)
    log_p0 = r * (log_r - log_r_mu)
    p0 = np.exp(np.clip(log_p0, -745.0, -1e-12))
    one_minus_p0 = np.clip(1.0 - p0, 1e-12, 1.0)
    log_one_minus_p0 = np.log(one_minus_p0)

    logpmf = (
        gammaln(y + r)
        - gammaln(r)
        - gammaln(y + 1)
        + r * (log_r - log_r_mu)
        + y * (eta - log_r_mu)
    )
    ll_obs = logpmf - log_one_minus_p0
    llf = float(np.sum(ll_obs))

    p0_ratio = p0 / one_minus_p0
    score_eta = (
        r * (y - mu) / (r + mu)
        - mu * r * p0_ratio / (r + mu)
    )
    score_beta = x * score_eta[:, None]

    dlogpmf_dr = (
        digamma(y + r)
        - digamma(r)
        + log_r
        + 1.0
        - log_r_mu
        - (y + r) / (r + mu)
    )
    dlogp0_dr = log_r - log_r_mu + 1.0 - r / (r + mu)
    score_r = dlogpmf_dr + p0_ratio * dlogp0_dr
    score_log_alpha = -r * score_r

    score_obs = np.column_stack([score_beta, score_log_alpha])
    score = np.sum(score_obs, axis=0)
    return llf, score, score_obs


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
    """Numerical Hessian of the ZTNB log-likelihood via first-order finite
    differences on the analytical gradient.

    Returns ``H`` such that ``H[i, j] = ∂²ℓ / ∂θᵢ ∂θⱼ`` at ``params``. The
    Fisher information matrix (bread of the sandwich) is ``-H``.
    """
    n = len(params)
    _, grad0, _ = ztnb_loglike_score(params, y, x, offset)
    H = np.zeros((n, n))
    for i in range(n):
        p_fwd = params.copy()
        p_fwd[i] += eps
        _, grad_fwd, _ = ztnb_loglike_score(p_fwd, y, x, offset)
        H[i] = (grad_fwd - grad0) / eps
    return (H + H.T) / 2.0  # symmetrise to correct for numerical asymmetry


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
    maxiter: int,
    offset: np.ndarray | None = None,
) -> ZTNBResult:
    """Fit a zero-truncated negative binomial with cluster-robust SEs.

    Uses the numerical Hessian as the bread of the sandwich estimator, which
    is more robust than the OPG under model misspecification (precisely the
    setting where cluster-robust inference is most needed).
    """
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
    for group_idx in range(len(group_codes)):
        cluster_scores[group_idx, :] = score_obs[inverse == group_idx].sum(axis=0)

    meat = cluster_scores.T @ cluster_scores
    cov = bread_inv @ meat @ bread_inv
    n, p = x_array.shape
    if len(group_codes) > 1 and n > p:
        correction = (len(group_codes) / (len(group_codes) - 1)) * ((n - 1) / (n - p))
        cov *= correction
    bse = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
    z_values = np.divide(params, bse, out=np.full_like(params, np.nan), where=bse > 0)
    pvalues = 2 * norm.sf(np.abs(z_values))

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
# Cluster-robust inference for binomial GLMs (wave-stratified path)
# ---------------------------------------------------------------------------


def clustered_logit_inference(
    result,
    y: pd.Series,
    x: pd.DataFrame,
    groups: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Cluster-robust logit standard errors using a pseudo-inverse bread.

    Statsmodels' clustered covariance path can fail with a singular Hessian
    in wave-stratified models. The fitted coefficients are still usable, so
    compute the sandwich covariance directly with a generalised inverse.
    """
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    params = np.asarray(result.params, dtype=float)
    mu = expit(x_array @ params)
    mu = np.clip(mu, 1e-9, 1.0 - 1e-9)

    weights = mu * (1.0 - mu)
    information = (x_array * weights[:, None]).T @ x_array
    bread_inv = pinvh(information, rtol=1e-10)

    score_obs = x_array * (y_array - mu)[:, None]
    group_codes, inverse = np.unique(groups, return_inverse=True)
    cluster_scores = np.zeros((len(group_codes), x_array.shape[1]), dtype=float)
    for group_idx in range(len(group_codes)):
        cluster_scores[group_idx, :] = score_obs[inverse == group_idx].sum(axis=0)

    meat = cluster_scores.T @ cluster_scores
    cov = bread_inv @ meat @ bread_inv
    n, p = x_array.shape
    if len(group_codes) > 1 and n > p:
        correction = (len(group_codes) / (len(group_codes) - 1)) * ((n - 1) / (n - p))
        cov *= correction

    bse = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
    z_values = np.divide(params, bse, out=np.full_like(params, np.nan), where=bse > 0)
    pvalues = 2 * norm.sf(np.abs(z_values))
    return bse, pvalues


def stable_binomial_fit_stats(result, y: pd.Series, x: pd.DataFrame) -> tuple[float, float]:
    """Binomial log-likelihood and AIC with clipped fitted probabilities.

    In some wave-stratified hurdle models, fitted probabilities can saturate
    at exactly 0 or 1 for a small number of rows. Calling ``result.llf``
    then asks statsmodels to evaluate ``log(0)``, producing noisy runtime
    warnings and NaN diagnostics even though the fitted coefficients and
    robust inference are usable.
    """
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    params = np.asarray(result.params, dtype=float)
    mu = expit(x_array @ params)
    mu = np.clip(mu, 1e-12, 1.0 - 1e-12)
    llf = float(np.sum(y_array * np.log(mu) + (1.0 - y_array) * np.log1p(-mu)))
    aic = float(-2.0 * llf + 2.0 * len(params))
    return llf, aic
