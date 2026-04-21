"""Regression models for Chapter 1 (socioeconomic deprivation).

All cluster-level regressions used by the Chapter 1 scripts live here so
every figure pulls from one source of truth. The four public model
entry points are:

- ``cluster_size_model``             — headline NB GLM of cluster size on
                                        deprivation quintile (or a
                                        standardised rank), VOC and a
                                        natural-cubic-spline on
                                        ``wn_mid_date``.
- ``build_domain_forest_table``      — per-SIMD-domain NB GLMs (one per
                                        domain, ``estimate`` is an IRR per
                                        1-SD of deprivation).
- ``build_domain_decomposition_table`` — single joint NB GLM with all
                                        seven domain ranks plus spline /
                                        VOC; returns Shapley shares of the
                                        domain-attributable log-likelihood
                                        gain.
- ``build_singleton_epoch_table``    — per-epoch logistic GLMs of
                                        singleton status on quintile
                                        dummies.

The module was previously named ``cluster_size_simd``; it was renamed to
reflect that it covers both the cluster-size NB and the singleton
logistic models.
"""

from __future__ import annotations

from itertools import combinations
from functools import lru_cache
from typing import Iterable
from dataclasses import dataclass, asdict
from collections import defaultdict
from math import factorial

import numpy as np
import pandas as pd
import patsy
from scipy.special import gammaln
from sklearn.model_selection import KFold

from analysis.utils import data, stats


# ---------------------------------------------------------------------------
# Column maps
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DOMAINS:
    overall: str = "dz_simd_rank"
    income: str = "dz_simd_income_rank"
    employment: str = "dz_simd_employment_rank"
    education: str = "dz_simd_education_rank"
    health: str = "dz_simd_health_rank"
    access: str = "dz_simd_access_rank"
    crime: str = "dz_simd_crime_rank"
    housing: str = "dz_simd_housing_rank"

def _domain_only(domains: DOMAINS | None = None) -> dict[str, str]:
    """Return only the seven raw SIMD domain columns, excluding overall rank."""
    return {
        name: col
        for name, col in asdict(domains or DOMAINS()).items()
        if name != "overall"
    }

def _significant_domains(domains: DOMAINS | None = None) -> dict[str, str]:
    """Return only the five significant SIMD domain columns.

    Initial univariate screening of the seven SIMD domains revealed
    that Geographic Access and Housing were not significantly associated with
    [n_sequences] ($p=0.92$ and $p=0.52$ respectively). Consequently, these
    domains were excluded from the Shapley log-likelihood decomposition
    to focus on the relative contributions of the primary socio-economic drivers.
    """
    return {
        name: col
        for name, col in asdict(domains or DOMAINS()).items()
        if name not in {"access", "housing", "overall"}
    }

# Spline degrees of freedom. VOC dummies already absorb the coarse wave
# structure, so the cr() spline only needs to capture within-VOC drift.
# df=4 is plenty for the 2.5-year series and keeps the IRLS well-
# conditioned alongside the VOC dummies (a higher df fights the VOC
# columns for the same variance and made NB IRLS blow up memory).
TIME_SPLINE_DF = 4
TIME_SPLINE_DF_EPOCH = 3
MODEL_CACHE_VERSION = 2


# ---------------------------------------------------------------------------
# Regression frame
# ---------------------------------------------------------------------------


@lru_cache(maxsize=4)
def build_cluster_regression_frame(
    resolution: float = data.PRIMARY_RESOLUTION,
    min_size: int = 1,
) -> pd.DataFrame:
    """Build the one-row-per-cluster regression frame used by Chapter 1.

    Columns
    -------
    cluster_id, window_id, wn_mid_date, wn_prop_sequenced, who_voc, pango_lineage,
    n_sequences, is_singleton, qc_frac_mediocre, qc_frac_bad, simd_quintile_mode, simd_rank_mean,
    dz_simd_income_rank, dz_simd_employment_rank, dz_simd_education_rank,
    dz_simd_health_rank, dz_simd_access_rank, dz_simd_crime_rank, dz_simd_housing_rank.
    """
    cols = [
        "window_id", "cluster_id", "resolution", "sequence_id",
        "wn_mid_date", "wn_prop_sequenced", "who_voc", "pango_lineage",
        "datazone", "dz_simd_rank", "dz_simd_quintile",
        "dz_simd_income_rank", "dz_simd_employment_rank", "dz_simd_education_rank",
        "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank",
        "dz_simd_housing_rank", "nextclade_qc",
    ]
    df = data.load_analysis_columns(cols, resolution=resolution, qc=None)
    df["_is_mediocre"] = (df["nextclade_qc"] == "mediocre").astype(float)
    df["_is_bad"]      = (df["nextclade_qc"] == "bad").astype(float)

    grp = df.groupby(["window_id", "cluster_id"], observed=True)

    def _mode(s: pd.Series):
        m = s.mode()
        return m.iloc[0] if len(m) > 0 else np.nan

    agg: dict[str, tuple] = {
        "n_sequences":         ("sequence_id", "nunique"),
        "wn_mid_date":         ("wn_mid_date", "first"),
        "wn_prop_sequenced":   ("wn_prop_sequenced", "first"),
        "who_voc":             ("who_voc", _mode),
        "pango_lineage":       ("pango_lineage", _mode),
        "qc_frac_mediocre":    ("_is_mediocre", "mean"),
        "qc_frac_bad":         ("_is_bad", "mean"),
        "simd_quintile_mode":  ("dz_simd_quintile", _mode),
        "dz_simd_rank":      ("dz_simd_rank", "mean"),
    }
    for col in _domain_only().values():
        agg[col] = (col, "mean")

    out = grp.agg(**agg).reset_index()
    out["is_singleton"] = (out["n_sequences"] == 1).astype(int)
    if min_size > 1:
        out = out[out["n_sequences"] >= min_size].reset_index(drop=True)
    return out


# ---------------------------------------------------------------------------
# Shared design-matrix helpers
# ---------------------------------------------------------------------------


def _time_spline(
    dates: pd.Series,
    *,
    df: int = TIME_SPLINE_DF,
    origin: pd.Timestamp = None,
) -> pd.DataFrame:
    """Natural-cubic-spline design matrix on ``dates`` (no intercept column).

    Returns a DataFrame with ``df`` basis columns named ``t_spline_1…_df``
    and an index that matches the input. The intercept is dropped (``-1``
    in the patsy formula) because the downstream GLM helper calls
    ``sm.add_constant``.
    """
    t = pd.to_datetime(dates)
    if origin is None:
        origin = t.min()
    days = (t - origin).dt.days.astype(float).to_numpy()
    # Normalise to [0, 1]. The cr() basis is invariant to monotone scaling,
    # but keeping the predictor in [0, 1] markedly improves IRLS
    # conditioning alongside unit-SD standardised deprivation and 0/1 VOC
    # dummies — an unnormalised ``days`` column (values up to ~930) made
    # the NB fits oscillate and blow up memory at df ≥ 5.
    span = float(days.max() - days.min()) or 1.0
    x = (days - days.min()) / span
    # ``cr`` is the natural (restricted) cubic regression spline. The
    # ``- 1`` drops the intercept column so statsmodels' own
    # ``add_constant`` doesn't create a duplicate.
    basis = patsy.dmatrix(
        f"cr(x, df={df}) - 1", {"x": x}, return_type="dataframe"
    )
    basis.columns = [f"t_spline_{i + 1}" for i in range(basis.shape[1])]
    basis.index = dates.index
    return basis


def _voc_dummies(voc_col: pd.Series) -> pd.DataFrame:
    return stats.one_hot(voc_col.fillna("None"), reference="Omicron", prefix="voc")


def _standardise(s: pd.Series, *, flip_sign: bool = False) -> pd.Series:
    z = (s - s.mean()) / s.std()
    return -z if flip_sign else z


def _good_qc_only(frame: pd.DataFrame) -> pd.DataFrame:
    """Keep only clusters composed entirely of Nextclade-good genomes."""
    return frame[
        (frame["qc_frac_mediocre"] <= 0) & (frame["qc_frac_bad"] <= 0)
    ].copy()


def _append_qc_covariates(
    X: pd.DataFrame,
    frame: pd.DataFrame,
    *,
    qc_sensitivity: bool,
) -> pd.DataFrame:
    if not qc_sensitivity:
        return X
    return pd.concat([X, frame[["qc_frac_mediocre", "qc_frac_bad"]]], axis=1)


def is_current_model_output(tab: pd.DataFrame) -> bool:
    """True when a cached table was written by the current model version."""
    if tab.empty or "model_version" not in tab.columns:
        return False
    versions = pd.to_numeric(tab["model_version"], errors="coerce")
    if versions.isna().any():
        return False
    return versions.nunique() == 1 and int(versions.iloc[0]) == MODEL_CACHE_VERSION


def tag_model_output(tab: pd.DataFrame) -> pd.DataFrame:
    """Stamp a cached output table with the current model version."""
    out = tab.copy()
    out["model_version"] = MODEL_CACHE_VERSION
    return out


# ---------------------------------------------------------------------------
# Cluster-size NB model
# ---------------------------------------------------------------------------

def cluster_size_model(
    frame: pd.DataFrame,
    *,
    deprivation_measure: str = "simd_quintile_mode",
    include_voc: bool = True,
    include_time_spline: bool = True,
    time_spline_df: int = TIME_SPLINE_DF,
    qc_sensitivity: bool = False,
):
    """NB fit: cluster size ~ deprivation + VOC + cr(wn_mid_date) (offset=wn_prop_sequenced).

    The `cr` spline absorbs residual temporal trend within each VOC epoch
     (e.g. the long Delta tail in some regions; Alpha's surveillance ramp-up).
     Return model data and statsmodels GLMResults object.
    """
    required = [
        deprivation_measure, "n_sequences",
        "wn_prop_sequenced", "wn_mid_date", "who_voc",
        "qc_frac_mediocre", "qc_frac_bad",
    ]
    df = frame.dropna(subset=required).copy()
    if df.empty:
        return None, df

    if not qc_sensitivity:
        df = _good_qc_only(df)
    if df.empty:
        return None, df

    if deprivation_measure == "simd_quintile_mode":
        X = stats.one_hot(df[deprivation_measure].astype(int), reference=5, prefix="q")
    else:
        X = pd.DataFrame(
            {"deprivation_sd": _standardise(df[deprivation_measure], flip_sign=True)},
            index=df.index,
        )

    if include_voc:
        X = pd.concat([X, _voc_dummies(df["who_voc"])], axis=1)

    if include_time_spline:
        X = pd.concat([X, _time_spline(df["wn_mid_date"], df=time_spline_df)], axis=1)
    X = _append_qc_covariates(X, df, qc_sensitivity=qc_sensitivity)

    df_model = pd.concat([df[["n_sequences", "wn_prop_sequenced"]], X], axis=1)
    fit = stats.negbin_cluster_size(
        df_model,
        outcome="n_sequences",
        predictors=X.columns.tolist(),
        offset="wn_prop_sequenced",
    )
    return fit, df_model


# ---------------------------------------------------------------------------
# Per-domain forest
# ---------------------------------------------------------------------------

def build_domain_forest_table(
    frame: pd.DataFrame,
    domains: DOMAINS = None,
) -> pd.DataFrame:
    """Tidy table of IRRs for a 1-SD increase in deprivation, one row per domain."""
    domains: DOMAINS = domains or DOMAINS()
    rows = []
    # Use asdict to allow iteration over field names and values
    for name, col in asdict(domains).items():
        fit, df_model = cluster_size_model(frame, deprivation_measure=col)

        if fit is None:
            rows.append({
                "estimate": np.nan, "conf_low": np.nan, "conf_high": np.nan,
                "p_value": np.nan, "n": 0, "domain": name,
            })
            continue

        tidy = stats.tidy_glm(fit)

        # Ensure the term exists before accessing
        mask = tidy["term"] == "deprivation_sd"
        if not mask.any():
            continue

        r = tidy.loc[mask].iloc[0]
        rows.append({
            "estimate": r["estimate"],
            "conf_low": r["conf_low"],
            "conf_high": r["conf_high"],
            "p_value": r["p_value"],
            "n": int(len(df_model)),
            "domain": name,
        })

    return pd.DataFrame(rows).set_index("domain")


# ---------------------------------------------------------------------------
# Per-epoch singleton odds
# ---------------------------------------------------------------------------

def singleton_epoch_model(
    frame: pd.DataFrame,
    *,
    min_prop_seq: float = 1e-3,
    include_time_spline: bool = True,
    time_spline_df: int = TIME_SPLINE_DF_EPOCH,
    qc_sensitivity: bool = False,
):
    """Logistic GLM: is_singleton ~ Q1..Q4 (vs Q5) + VOC + cr(wn_mid_date) (offset=wn_prop_sequenced)

    ``frame`` must already be restricted to a single epoch; within-epoch
    the VOC dummies collapse, so only the time spline captures residual
    temporal trend. Returns ``(fit, df_model)`` or ``(None, df)`` if the
    frame is too sparse to fit.
    """
    required = [
        "simd_quintile_mode", "is_singleton",
        "wn_prop_sequenced", "wn_mid_date",
        "qc_frac_mediocre", "qc_frac_bad",
    ]
    df = frame.dropna(subset=required).copy()
    df = df[df["wn_prop_sequenced"] >= min_prop_seq]

    if not qc_sensitivity:
        df = _good_qc_only(df)
    if len(df) < 50 or df["simd_quintile_mode"].nunique() < 2:
        return None, df

    q_dummies = stats.one_hot(
        df["simd_quintile_mode"].astype(int), reference=5, prefix="q"
    )
    log_prop = np.log(df["wn_prop_sequenced"].astype(float))
    log_prop_z = (log_prop - log_prop.mean()) / log_prop.std()

    X = q_dummies.copy()
    X["log_prop_seq_z"] = log_prop_z.values
    if include_time_spline:
        X = pd.concat([X, _time_spline(df["wn_mid_date"], df=time_spline_df)], axis=1)
    # QC fraction covariates are already on [0, 1]; no further scaling needed.
    X = _append_qc_covariates(X, df, qc_sensitivity=qc_sensitivity)

    df_model = pd.concat(
        [df[["is_singleton"]].reset_index(drop=True), X.reset_index(drop=True)],
        axis=1,
    )
    fit = stats.logit_singleton(df_model, predictors=X.columns.tolist())
    return fit, df_model

def build_singleton_epoch_table(
    frame: pd.DataFrame,
    *,
    epochs: Iterable[str] = None,
) -> pd.DataFrame:
    """Loop over epochs; return a tidy OR table for quintile dummies.

    Columns: ``estimate`` (OR vs Q5), ``conf_low``/``conf_high``,
    ``p_value``, ``z``, ``std_error``, ``quintile`` (1..4), ``epoch``,
    ``model``.
    """
    frame = frame.copy()
    frame["epoch"] = data.assign_epoch(frame["wn_mid_date"])

    if epochs is None:
        epochs = [lbl for lbl, *_ in data.VOC_EPOCHS]

    out: dict[str, list[pd.DataFrame]] = defaultdict(list)

    def _get_tidy(tidy_df, model):
        tidy_df = tidy_df[tidy_df["term"].str.startswith("q_")].copy()
        tidy_df["quintile"] = (
            tidy_df["term"]
            .str.extract(r"^q_(\d+)$", expand=False)
            .astype(int)
        )
        tidy_df["epoch"] = epoch
        tidy_df["model"] = model
        return tidy_df

    for epoch in epochs:
        sub = frame[frame["epoch"] == epoch]
        if sub.empty:
            continue

        fit, _ = singleton_epoch_model(sub, qc_sensitivity=False)
        if fit is not None:
            tidy = stats.tidy_glm(fit)
            out["primary"].append(_get_tidy(tidy, "primary"))

        fit, _ = singleton_epoch_model(sub, qc_sensitivity=True)
        if fit is not None:
            tidy = stats.tidy_glm(fit)
            out["qc_adjusted"].append(_get_tidy(tidy, "qc_adjusted"))

    frames = [df for dfs in out.values() for df in dfs]

    if not frames:
        return pd.DataFrame(
            columns=[
                "term", "estimate", "std_error", "conf_low", "conf_high",
                "z", "p_value", "quintile", "epoch", "model"
            ]
        )

    return pd.concat(frames, ignore_index=True)[
        ["term", "estimate", "std_error", "conf_low", "conf_high",
         "z", "p_value", "quintile", "epoch", "model"]
    ]


# ---------------------------------------------------------------------------
# Cross-validated domain Shapley decomposition
# ---------------------------------------------------------------------------

def _fit_nb_design(
    frame: pd.DataFrame,
    X: pd.DataFrame,
):
    df_model = pd.concat([frame[["n_sequences", "wn_prop_sequenced"]], X], axis=1)
    fit = stats.negbin_cluster_size(
        df_model,
        outcome="n_sequences",
        predictors=X.columns.tolist(),
        offset="wn_prop_sequenced",
    )
    return fit, df_model


def _prepare_domain_decomposition_design(
    frame: pd.DataFrame,
    *,
    domains: DOMAINS,
    include_time_spline: bool,
    qc_sensitivity: bool,
) -> tuple[pd.DataFrame, dict[str, str], pd.DataFrame]:
    required = [
        "n_sequences", "wn_prop_sequenced", "wn_mid_date",
        "who_voc", "qc_frac_mediocre", "qc_frac_bad",
    ]
    domain_only = _significant_domains(domains)
    required.extend(domain_only.values())

    df = frame.dropna(subset=required).copy()
    if not qc_sensitivity:
        df = _good_qc_only(df)
    if df.empty:
        raise RuntimeError("No complete cases for joint-domain model.")

    dep_name_to_col: dict[str, str] = {}
    for name, col in domain_only.items():
        dep_col = f"dep_{name}"
        df[dep_col] = _standardise(df[col], flip_sign=True)
        dep_name_to_col[name] = dep_col

    base_blocks = [_voc_dummies(df["who_voc"])]
    if include_time_spline:
        base_blocks.append(_time_spline(df["wn_mid_date"]))
    if qc_sensitivity:
        base_blocks.append(df[["qc_frac_mediocre", "qc_frac_bad"]])
    X_base = pd.concat(base_blocks, axis=1)
    return df, dep_name_to_col, X_base

def _align_test_design_to_train(
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
) -> pd.DataFrame:
    """
    Align test design matrix columns to the training matrix.
    Missing columns are filled with 0; extra columns are dropped.
    """
    X_test = X_test.copy()
    for col in X_train.columns:
        if col not in X_test.columns:
            X_test[col] = 0.0
    return X_test[X_train.columns]


def _build_base_design_for_split(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    *,
    include_time_spline: bool,
    qc_sensitivity: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build train/test versions of the non-domain covariate design matrix.
    Ensures dummy/spline columns are aligned between train and test.
    """
    train_blocks = [_voc_dummies(train_df["who_voc"])]
    test_blocks = [_voc_dummies(test_df["who_voc"])]

    if include_time_spline:
        X_time_train = _time_spline(train_df["wn_mid_date"])
        X_time_test = _time_spline(test_df["wn_mid_date"])
        test_blocks[0] = _align_test_design_to_train(train_blocks[0], test_blocks[0])
        X_time_test = _align_test_design_to_train(X_time_train, X_time_test)
        train_blocks.append(X_time_train)
        test_blocks.append(X_time_test)
    else:
        test_blocks[0] = _align_test_design_to_train(train_blocks[0], test_blocks[0])

    if qc_sensitivity:
        train_blocks.append(train_df[["qc_frac_mediocre", "qc_frac_bad"]])
        test_blocks.append(test_df[["qc_frac_mediocre", "qc_frac_bad"]])

    X_base_train = pd.concat(train_blocks, axis=1)
    X_base_test = pd.concat(test_blocks, axis=1)

    X_base_test = _align_test_design_to_train(X_base_train, X_base_test)
    return X_base_train, X_base_test


def _prepare_domain_decomposition_train_test(
    train_frame: pd.DataFrame,
    test_frame: pd.DataFrame,
    *,
    domains: DOMAINS,
    include_time_spline: bool,
    qc_sensitivity: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str], pd.DataFrame, pd.DataFrame]:
    """
    Prepare consistent train/test data for CV domain decomposition.
    Standardisation for domain variables is learned on the training set only.
    """
    required = [
        "n_sequences", "wn_prop_sequenced", "wn_mid_date",
        "who_voc", "qc_frac_mediocre", "qc_frac_bad",
    ]
    domain_only = _significant_domains(domains)
    required.extend(domain_only.values())

    train_df = train_frame.dropna(subset=required).copy()
    test_df = test_frame.dropna(subset=required).copy()

    if not qc_sensitivity:
        train_df = _good_qc_only(train_df)
        test_df = _good_qc_only(test_df)

    if train_df.empty or test_df.empty:
        raise RuntimeError("Train/test split produced no complete cases for CV domain decomposition.")

    dep_name_to_col: dict[str, str] = {}

    # Standardise using training statistics only
    for name, col in domain_only.items():
        mu = train_df[col].mean()
        sd = train_df[col].std(ddof=0)

        if not np.isfinite(sd) or sd <= 0:
            train_df[f"dep_{name}"] = 0.0
            test_df[f"dep_{name}"] = 0.0
        else:
            train_df[f"dep_{name}"] = -(train_df[col] - mu) / sd
            test_df[f"dep_{name}"] = -(test_df[col] - mu) / sd

        dep_name_to_col[name] = f"dep_{name}"

    X_base_train, X_base_test = _build_base_design_for_split(
        train_df,
        test_df,
        include_time_spline=include_time_spline,
        qc_sensitivity=qc_sensitivity,
    )

    return train_df, test_df, dep_name_to_col, X_base_train, X_base_test


def _extract_alpha_from_fit(fit) -> float:
    """
    Best-effort extraction of NB dispersion parameter alpha.
    Assumes NB2 variance: Var(Y) = mu + alpha * mu^2.
    """
    # Common case: alpha estimated as a named parameter
    if hasattr(fit, "params"):
        params = fit.params
        if isinstance(params, pd.Series):
            if "alpha" in params.index:
                return float(params["alpha"])
        else:
            # statsmodels may return ndarray; sometimes alpha is last parameter,
            # but we avoid guessing unless there is metadata.
            pass

    # Some wrappers keep alpha directly
    for attr in ("alpha", "_alpha"):
        if hasattr(fit, attr):
            alpha = getattr(fit, attr)
            if np.isscalar(alpha) and np.isfinite(alpha) and alpha > 0:
                return float(alpha)

    # Sometimes stored on the model
    if hasattr(fit, "model"):
        model = fit.model
        for attr in ("alpha", "_alpha"):
            if hasattr(model, attr):
                alpha = getattr(model, attr)
                if np.isscalar(alpha) and np.isfinite(alpha) and alpha > 0:
                    return float(alpha)

    raise RuntimeError(
        "Could not extract NB dispersion parameter 'alpha' from fitted model. "
        "Inspect the object returned by stats.negbin_cluster_size(...) and "
        "adapt _extract_alpha_from_fit(...) accordingly."
    )


def _predict_mean_from_fit(fit, df_model: pd.DataFrame) -> np.ndarray:
    """
    Best-effort mean prediction from fitted model on a full model frame.
    """
    # Most statsmodels results wrappers
    try:
        mu = fit.predict(df_model)
        mu = np.asarray(mu, dtype=float)
        if mu.ndim == 1 and len(mu) == len(df_model):
            return np.clip(mu, 1e-12, np.inf)
    except Exception:
        pass

    # Sometimes predict needs exog only
    try:
        exog_cols = [c for c in df_model.columns if c != "n_sequences"]
        mu = fit.predict(df_model[exog_cols])
        mu = np.asarray(mu, dtype=float)
        if mu.ndim == 1 and len(mu) == len(df_model):
            return np.clip(mu, 1e-12, np.inf)
    except Exception:
        pass

    raise RuntimeError(
        "Could not obtain predictions from fitted model. "
        "Inspect the object returned by stats.negbin_cluster_size(...) and "
        "adapt _predict_mean_from_fit(...) accordingly."
    )


def _nb2_logpmf(
    y: np.ndarray,
    mu: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """
    Negative binomial NB2 log PMF:
        Var(Y) = mu + alpha * mu^2
    Equivalent shape parameter r = 1 / alpha.
    """
    y = np.asarray(y, dtype=float)
    mu = np.asarray(mu, dtype=float)

    mu = np.clip(mu, 1e-12, np.inf)
    alpha = float(alpha)
    if not np.isfinite(alpha) or alpha <= 0:
        raise RuntimeError(f"Invalid alpha for NB2 predictive log-likelihood: {alpha!r}")

    r = 1.0 / alpha
    return (
        gammaln(y + r)
        - gammaln(r)
        - gammaln(y + 1.0)
        + r * np.log(r / (r + mu))
        + y * np.log(mu / (r + mu))
    )


def _fit_domain_subset_predictive_loglik(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    X_base_train: pd.DataFrame,
    X_base_test: pd.DataFrame,
    dep_name_to_col: dict[str, str],
    subset: tuple[str, ...],
) -> float:
    """
    Fit the subset model on the training split and evaluate predictive
    NB log-likelihood on the test split.
    """
    train_blocks = []
    test_blocks = []

    dep_cols = [dep_name_to_col[name] for name in subset]
    if dep_cols:
        train_blocks.append(train_df[dep_cols])
        test_blocks.append(test_df[dep_cols])

    if X_base_train.shape[1] > 0:
        train_blocks.append(X_base_train)
        test_blocks.append(X_base_test)

    X_train = pd.concat(train_blocks, axis=1) if train_blocks else pd.DataFrame(index=train_df.index)
    X_test = pd.concat(test_blocks, axis=1) if test_blocks else pd.DataFrame(index=test_df.index)

    fit, _ = _fit_nb_design(train_df, X_train)

    fit, df_model = _fit_nb_design(df, X_full)

    # ── introspect the fit object ──────────────────────────────────────────────
    import pprint

    print(type(fit))
    print(dir(fit))

    # If it has .params (a Series or dict), print them
    if hasattr(fit, "params"):
        print("\n── params ──")
        print(fit.params)

    # If it has a .model attribute, inspect that too
    if hasattr(fit, "model"):
        print("\n── fit.model attrs ──")
        print(dir(fit.model))
        if hasattr(fit.model, "params"):
            print(fit.model.params)

    # Print all scalar numeric attributes on fit itself
    print("\n── scalar numeric attrs on fit ──")
    for attr in dir(fit):
        if attr.startswith("_"):
            continue
        try:
            val = getattr(fit, attr)
            if np.isscalar(val) and np.isfinite(float(val)):
                print(f"  {attr}: {val}")
        except Exception:
            pass

    # Build the test model frame expected by the predictor
    df_model_test = pd.concat(
        [test_df[["n_sequences", "wn_prop_sequenced"]], X_test],
        axis=1,
    )

    mu_test = _predict_mean_from_fit(fit, df_model_test)
    alpha = _extract_alpha_from_fit(fit)

    y_test = test_df["n_sequences"].to_numpy(dtype=float)
    ll = _nb2_logpmf(y_test, mu_test, alpha).sum()
    return float(ll)


def _domain_cv_shapley_loglik_gain(
    frame: pd.DataFrame,
    domains: DOMAINS,
    *,
    include_time_spline: bool,
    qc_sensitivity: bool,
    n_splits: int = 5,
    random_state: int = 42,
) -> pd.Series:
    """
    Exact Shapley decomposition of cross-validated predictive log-likelihood gain
    attributable to the domain blocks.
    """
    # Use the same complete-case / QC-filtered pool before splitting
    required = [
        "n_sequences", "wn_prop_sequenced", "wn_mid_date",
        "who_voc", "qc_frac_mediocre", "qc_frac_bad",
    ]
    domain_only = _significant_domains(domains)
    required.extend(domain_only.values())

    df_all = frame.dropna(subset=required).copy()
    if not qc_sensitivity:
        df_all = _good_qc_only(df_all)

    if df_all.empty:
        raise RuntimeError("No complete cases available for CV domain decomposition.")

    if len(df_all) < n_splits:
        raise RuntimeError(
            f"Not enough rows ({len(df_all)}) for {n_splits}-fold CV "
            "in domain decomposition."
        )

    domain_names = list(domain_only.keys())
    subset_scores = {frozenset(s): 0.0 for r in range(len(domain_names) + 1)
                     for s in combinations(domain_names, r)}

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    for train_idx, test_idx in kf.split(df_all):
        train_raw = df_all.iloc[train_idx].copy()
        test_raw = df_all.iloc[test_idx].copy()

        train_df, test_df, dep_name_to_col, X_base_train, X_base_test = (
            _prepare_domain_decomposition_train_test(
                train_raw,
                test_raw,
                domains=domains,
                include_time_spline=include_time_spline,
                qc_sensitivity=qc_sensitivity,
            )
        )

        for r in range(len(domain_names) + 1):
            for subset in combinations(domain_names, r):
                subset_scores[frozenset(subset)] += _fit_domain_subset_predictive_loglik(
                    train_df=train_df,
                    test_df=test_df,
                    X_base_train=X_base_train,
                    X_base_test=X_base_test,
                    dep_name_to_col=dep_name_to_col,
                    subset=subset,
                )

    n_domains = len(domain_names)
    weight_denom = factorial(n_domains)
    gain_by_domain: dict[str, float] = {}

    for name in domain_names:
        others = [other for other in domain_names if other != name]
        gain = 0.0
        for r in range(len(others) + 1):
            weight = factorial(r) * factorial(n_domains - r - 1) / weight_denom
            for subset in combinations(others, r):
                s = frozenset(subset)
                marginal = subset_scores[s | {name}] - subset_scores[s]
                gain += weight * marginal
        gain_by_domain[name] = gain

    return pd.Series(gain_by_domain, dtype=float)


def build_domain_decomposition_table_cv_shapley(
    frame: pd.DataFrame,
    domains: DOMAINS = None,
    *,
    include_time_spline: bool = True,
    qc_sensitivity: bool = False,
    n_splits: int = 5,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Fit the full joint NB model for coefficient estimates, but quantify
    domain contribution using exact Shapley decomposition of cross-validated
    predictive log-likelihood.

    Returns columns:
    domain, estimate, std_error, conf_low, conf_high, p_value,
    cv_loglik_gain, share

    Notes
    -----
    - estimate / std_error / CIs / p-values come from the full-sample joint NB fit
    - cv_loglik_gain / share come from cross-validated predictive performance
    """
    domains = domains or DOMAINS()

    # Full-data joint fit for coefficient summaries
    df, dep_name_to_col, X_base = _prepare_domain_decomposition_design(
        frame,
        domains=domains,
        include_time_spline=include_time_spline,
        qc_sensitivity=qc_sensitivity,
    )

    dep_cols = [dep_name_to_col[name] for name in dep_name_to_col]
    blocks = [df[dep_cols]]
    if X_base.shape[1] > 0:
        blocks.append(X_base)
    X_full = pd.concat(blocks, axis=1)

    fit, _ = _fit_nb_design(df, X_full)
    tidy = stats.tidy_glm(fit, exponentiate=False).copy()
    tidy = tidy[tidy["term"].str.startswith("dep_")].copy()
    tidy["domain"] = tidy["term"].str.replace("dep_", "", regex=False)

    # CV Shapley gains
    cv_gain = _domain_cv_shapley_loglik_gain(
        frame=frame,
        domains=domains,
        include_time_spline=include_time_spline,
        qc_sensitivity=qc_sensitivity,
        n_splits=n_splits,
        random_state=random_state,
    )

    tidy["cv_loglik_gain"] = tidy["domain"].map(cv_gain)

    total_gain = float(tidy["cv_loglik_gain"].sum())
    if not np.isfinite(total_gain) or total_gain == 0:
        raise RuntimeError(
            "Cross-validated Shapley decomposition produced zero or non-finite "
            "total predictive log-likelihood gain."
        )

    tidy["share"] = tidy["cv_loglik_gain"] / total_gain
    return tidy.sort_values("share", ascending=True).reset_index(drop=True)
