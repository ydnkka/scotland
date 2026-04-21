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
                                        VOC; returns shares of total
                                        |standardised coef|.
- ``build_singleton_epoch_table``    — per-epoch logistic GLMs of
                                        singleton status on quintile
                                        dummies.

The module was previously named ``cluster_size_simd``; it was renamed to
reflect that it covers both the cluster-size NB and the singleton
logistic models.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Iterable
from dataclasses import dataclass, asdict
from collections import defaultdict

import numpy as np
import pandas as pd
import patsy

from manuscripts.common import data, stats


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

# Spline degrees of freedom. VOC dummies already absorb the coarse wave
# structure, so the cr() spline only needs to capture within-VOC drift.
# df=4 is plenty for the 2.5-year series and keeps the IRLS well-
# conditioned alongside the VOC dummies (a higher df fights the VOC
# columns for the same variance and made NB IRLS blow up memory).
TIME_SPLINE_DF = 4
TIME_SPLINE_DF_EPOCH = 3


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
        df = df[(df["qc_frac_mediocre"] <= 0) | (df["qc_frac_bad"] <= 0)]

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

    if len(df) < 50 or df["simd_quintile_mode"].nunique() < 2:
        return None, df

    if not qc_sensitivity:
        df = df[(df["qc_frac_mediocre"] <= 0) | (df["qc_frac_bad"] <= 0)]

    q_dummies = stats.one_hot(
        df["simd_quintile_mode"].astype(int), reference=5, prefix="q"
    )
    log_prop = np.log(df["wn_prop_sequenced"].astype(float))
    log_prop_z = (log_prop - log_prop.mean()) / log_prop.std()

    X = q_dummies.copy()
    X["log_prop_seq_z"] = log_prop_z.values
    if include_time_spline:
        X = pd.concat([X, _time_spline(df["wn_mid_date"], df=time_spline_df)], axis=1)

    if not qc_sensitivity:
        # QC fraction covariates — already on [0, 1];
        # no further standardisation needed,
        X = pd.concat([X, df[["qc_frac_mediocre", "qc_frac_bad"]]], axis=1)

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
        tidy_df["quintile"] = tidy["term"].str.replace("q_", "", regex=False).astype(int)
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
# Joint domain decomposition
# ---------------------------------------------------------------------------

def build_domain_decomposition_table(
    frame: pd.DataFrame,
    domains: DOMAINS = None,
    *,
    include_time_spline: bool = True,
    qc_sensitivity: bool = False,
) -> pd.DataFrame:
    """Joint NB fit across all seven domain ranks; return tidy share table.

    Columns
    -------
    domain, estimate (standardised log-IRR — NOT exponentiated),
    std_error, conf_low, conf_high, p_value, abs_std_coef, share.
    Shares sum to 1 across rows; rows are sorted ascending on ``share``.
    """
    domains: DOMAINS = domains or DOMAINS()
    required = [
        "n_sequences", "wn_prop_sequenced", "wn_mid_date",
        "who_voc", "qc_frac_mediocre", "qc_frac_bad"
    ]
    domain_only = _domain_only(domains)
    required.extend(domain_only.values())

    df = frame.dropna(subset=required).copy()
    if df.empty:
        raise RuntimeError("No complete cases for joint-domain model.")


    if not qc_sensitivity:
        df = df[(df["qc_frac_mediocre"] <= 0) | (df["qc_frac_bad"] <= 0)]

    dep_cols: list[str] = []
    for name, col in domain_only.items():
        c = f"dep_{name}"
        df[c] = _standardise(df[col], flip_sign=True)
        dep_cols.append(c)

    X = df[dep_cols].copy()
    X = pd.concat([X, _voc_dummies(df["who_voc"])], axis=1)
    if include_time_spline:
        X = pd.concat([X, _time_spline(df["wn_mid_date"])], axis=1)

    df_model = pd.concat([df[["n_sequences", "wn_prop_sequenced"]], X], axis=1)
    fit = stats.negbin_cluster_size(
        df_model,
        outcome="n_sequences",
        predictors=X.columns.tolist(),
        offset="wn_prop_sequenced",
    )
    tidy = stats.tidy_glm(fit, exponentiate=False)  # keep on log-IRR scale
    tidy = tidy[tidy["term"].str.startswith("dep_")].copy()
    tidy["domain"] = tidy["term"].str.replace("dep_", "")
    tidy["abs_std_coef"] = tidy["estimate"].abs()
    tidy["share"] = tidy["abs_std_coef"] / tidy["abs_std_coef"].sum()
    return tidy.sort_values("share", ascending=True).reset_index(drop=True)
