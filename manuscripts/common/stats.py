"""Statistical helpers shared by the three manuscripts.

Kept deliberately small: pandas/statsmodels/scipy only. None of the figures
need Bayesian inference in the primary analysis; sensitivity analyses can
be added per paper.

Functions
---------
negbin_cluster_size
    Fit a negative-binomial GLM of cluster size on a design matrix, with
    optional offset for sequencing proportion.
logit_singleton
    Fit a logistic GLM for Pr(singleton) vs. covariates.
tidy_glm
    Convert a statsmodels GLM result to a tidy DataFrame with IRR/OR and 95% CI.
bootstrap_ci
    Non-parametric bootstrap confidence intervals for a scalar statistic.
"""

from __future__ import annotations

import contextlib
import warnings
from typing import Callable, Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm


@contextlib.contextmanager
def _suppress_glm_runtime_warnings():
    """Silence harmless numerical-stability warnings from statsmodels.

    During IRLS the logistic link computes ``exp(-z)``; whenever a single
    iterate of the linear predictor becomes very negative this overflows
    (saturating correctly at inf → sigmoid = 0), which does not affect the
    fitted coefficients but clutters the output. Ditto rare
    `log(0)` warnings in the deviance computation on perfectly-fitted rows.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=RuntimeWarning, message=".*overflow encountered in exp.*",
        )
        warnings.filterwarnings(
            "ignore", category=RuntimeWarning, message=".*divide by zero encountered in log.*",
        )
        yield


# ---------------------------------------------------------------------------
# GLMs
# ---------------------------------------------------------------------------


def negbin_cluster_size(
    df: pd.DataFrame,
    outcome: str,
    predictors: Iterable[str],
    *,
    offset: str | None = None,
    alpha: float = 1.0,
):
    """Fit a negative-binomial GLM.

    Parameters
    ----------
    df : DataFrame
        One row per cluster.
    outcome : str
        Column name for the count outcome (typically `n_sequences`).
    predictors : iterable of str
        Design matrix column names. Categorical dummies should already be
        materialised (e.g. via `pd.get_dummies`).
    offset : str, optional
        Column name for an offset on the log scale (e.g. log of sequencing
        proportion, to adjust for surveillance intensity).
    alpha : float
        Negative-binomial dispersion parameter passed to statsmodels.
    """
    X = sm.add_constant(df[list(predictors)].astype(float), has_constant="add")
    y = df[outcome].astype(float)
    kwargs = {}
    if offset is not None:
        kwargs["offset"] = np.log(df[offset].clip(lower=1e-6).astype(float))
    model = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha), **kwargs)
    with _suppress_glm_runtime_warnings():
        return model.fit()


def logit_singleton(
    df: pd.DataFrame,
    predictors: Iterable[str],
    *,
    outcome: str = "is_singleton",
):
    """Fit a logistic GLM for Pr(singleton)."""
    X = sm.add_constant(df[list(predictors)].astype(float), has_constant="add")
    y = df[outcome].astype(float)
    model = sm.GLM(y, X, family=sm.families.Binomial())
    with _suppress_glm_runtime_warnings():
        return model.fit()


# ---------------------------------------------------------------------------
# Tidy outputs
# ---------------------------------------------------------------------------


def tidy_glm(fit, *, exponentiate: bool = True, alpha: float = 0.05) -> pd.DataFrame:
    """Convert a statsmodels fit to a tidy data frame.

    With `exponentiate=True`, the estimate column becomes an IRR (for NB)
    or OR (for logistic).
    """
    params = fit.params
    se = fit.bse
    z = params / se
    # Two-sided normal-approx p-value.
    from scipy.stats import norm
    # p = 2 * (1 - norm.cdf(np.abs(z)))
    p = 2 * norm.sf(np.abs(z))
    q = norm.ppf(1 - alpha / 2)
    # low, high = params - q * se, params + q * se
    low, high = norm.interval(1 - alpha, loc=params, scale=se)
    # norm.interval may return a numpy ndarray (newer scipy) or a pandas
    # Series (older scipy) depending on the type of loc/scale.  np.asarray
    # normalises both cases so the DataFrame constructor is always happy.
    tidy = pd.DataFrame(
        {
            "term": params.index,
            "estimate": np.asarray(params),
            "std_error": np.asarray(se),
            "conf_low": np.asarray(low),
            "conf_high": np.asarray(high),
            "z": np.asarray(z),
            "p_value": np.asarray(p),
        }
    )
    if exponentiate:
        # Very large CI endpoints (near-separation) can overflow `exp`; the
        # resulting `inf` is the correct OR/IRR, so silence the warning.
        with _suppress_glm_runtime_warnings():
            for c in ("estimate", "conf_low", "conf_high"):
                tidy[c] = np.exp(tidy[c])
    return tidy


# ---------------------------------------------------------------------------
# Bootstrap CI
# ---------------------------------------------------------------------------


def bootstrap_ci(
    values: np.ndarray,
    statistic: Callable[[np.ndarray], float] = np.mean,
    *,
    n_boot: int = 2000,
    alpha: float = 0.05,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Return (point, low, high) via the percentile bootstrap."""
    rng = np.random.default_rng(seed)
    values = np.asarray(values)
    n = len(values)
    if n == 0:
        return (np.nan, np.nan, np.nan)
    idx = rng.integers(0, n, size=(n_boot, n))
    boots = np.fromiter(
        (statistic(values[i]) for i in idx), dtype=float, count=n_boot
    )
    return (
        float(statistic(values)),
        float(np.quantile(boots, alpha / 2)),
        float(np.quantile(boots, 1 - alpha / 2)),
    )


# ---------------------------------------------------------------------------
# Convenience: one-hot with reference level control
# ---------------------------------------------------------------------------


def one_hot(
    series: pd.Series,
    *,
    reference: str | int | None = None,
    prefix: str | None = None,
) -> pd.DataFrame:
    """Dummy-encode `series`, dropping the reference level explicitly."""
    dummies = pd.get_dummies(series, prefix=prefix or series.name, drop_first=False)
    if reference is not None:
        col = f"{prefix or series.name}_{reference}"
        if col in dummies.columns:
            dummies = dummies.drop(columns=col)
    return dummies.astype(float)
