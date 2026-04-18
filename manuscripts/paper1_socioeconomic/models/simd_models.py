"""Regression models for Paper 1 (socioeconomic deprivation).

All cluster-level regressions used by the Paper 1 figures live here so
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

import numpy as np
import pandas as pd
import patsy

from manuscripts.common import data, stats


# ---------------------------------------------------------------------------
# Column maps
# ---------------------------------------------------------------------------

DOMAIN_COLUMNS: dict[str, str] = {
    "overall":    "simd_rank_mean",
    "income":     "dz_simd_income_rank",
    "employment": "dz_simd_employment_rank",
    "education":  "dz_simd_education_rank",
    "health":     "dz_simd_health_rank",
    "access":     "dz_simd_access_rank",
    "crime":      "dz_simd_crime_rank",
    "housing":    "dz_simd_housing_rank",
}

# "Overall" is not a per-domain column and must be excluded from any
# decomposition that fits all seven domain ranks jointly.
DOMAIN_ONLY = {k: v for k, v in DOMAIN_COLUMNS.items() if k != "overall"}

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
    """Build the one-row-per-cluster regression frame used by Paper 1.

    Columns
    -------
    cluster_id, window_id, wn_mid_date, wn_prop_sequenced, who_voc, pango_lineage,
    n_sequences, is_singleton,
    simd_quintile_mode, simd_rank_mean,
    dz_simd_income_rank, dz_simd_employment_rank, dz_simd_education_rank,
    dz_simd_health_rank, dz_simd_access_rank, dz_simd_crime_rank, dz_simd_housing_rank.
    """
    cols = [
        "window_id", "cluster_id", "resolution", "sequence_id",
        "wn_mid_date", "wn_prop_sequenced", "who_voc", "pango_lineage",
        "datazone", "dz_simd_rank", "dz_simd_quintile",
        "dz_simd_income_rank", "dz_simd_employment_rank", "dz_simd_education_rank",
        "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank",
        "dz_simd_housing_rank",
    ]
    df = data.load_analysis_columns(cols, resolution=resolution)
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
        "simd_quintile_mode":  ("dz_simd_quintile", _mode),
        "simd_rank_mean":      ("dz_simd_rank", "mean"),
    }
    for col in DOMAIN_ONLY.values():
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
    origin: pd.Timestamp | None = None,
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
# 1. Headline cluster-size NB model
# ---------------------------------------------------------------------------


def cluster_size_model(
    frame: pd.DataFrame,
    *,
    deprivation_measure: str = "simd_quintile_mode",
    include_voc: bool = True,
    include_time_spline: bool = True,
    time_spline_df: int = TIME_SPLINE_DF,
):
    """Headline NB fit: cluster size ~ deprivation + VOC + cr(wn_mid_date).

    The spline term matches what the README describes and absorbs
    residual temporal trend within each VOC epoch (e.g. the long Delta
    tail in some regions; Alpha's surveillance ramp-up). Returns a
    statsmodels GLMResults object.
    """
    required = [deprivation_measure, "n_sequences", "wn_prop_sequenced", "wn_mid_date"]
    df = frame.dropna(subset=required).copy()

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
    return stats.negbin_cluster_size(
        df_model,
        outcome="n_sequences",
        predictors=X.columns.tolist(),
        offset="wn_prop_sequenced",
    )


# ---------------------------------------------------------------------------
# 2. Per-domain forest (Fig. 3)
# ---------------------------------------------------------------------------


def domain_irr_model(
    frame: pd.DataFrame,
    *,
    domain_col: str,
    include_time_spline: bool = True,
):
    """One-domain NB fit: n_sequences ~ standardised-deprivation + VOC + spline."""
    required = [domain_col, "n_sequences", "wn_prop_sequenced", "wn_mid_date", "who_voc"]
    df = frame.dropna(subset=required).copy()
    if df.empty:
        return None, df
    X = pd.DataFrame(
        {"deprivation_sd": _standardise(df[domain_col], flip_sign=True)},
        index=df.index,
    )
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
    return fit, df_model


def build_domain_forest_table(
    frame: pd.DataFrame,
    domains: dict[str, str] = DOMAIN_COLUMNS,
) -> pd.DataFrame:
    """Tidy table of IRRs for a 1-SD increase in deprivation, one row per domain.

    Index: domain name. Columns: ``estimate``, ``conf_low``, ``conf_high``,
    ``p_value``, ``n``.
    """
    rows = []
    for name, col in domains.items():
        fit, df_model = domain_irr_model(frame, domain_col=col)
        if fit is None:
            rows.append(pd.Series({
                "estimate": np.nan, "conf_low": np.nan, "conf_high": np.nan,
                "p_value": np.nan, "n": 0, "domain": name,
            }))
            continue
        tidy = stats.tidy_glm(fit)
        r = tidy.loc[tidy["term"] == "deprivation_sd"].iloc[0]
        rows.append(pd.Series({
            "estimate": r["estimate"],
            "conf_low": r["conf_low"],
            "conf_high": r["conf_high"],
            "p_value": r["p_value"],
            "n": int(len(df_model)),
            "domain": name,
        }))
    return pd.DataFrame(rows).set_index("domain")


# ---------------------------------------------------------------------------
# 4. Per-epoch singleton odds (Fig. 4)
# ---------------------------------------------------------------------------


def singleton_epoch_model(
    frame: pd.DataFrame,
    *,
    min_prop_seq: float = 1e-3,
    include_time_spline: bool = True,
    time_spline_df: int = TIME_SPLINE_DF_EPOCH,
):
    """Logistic GLM: is_singleton ~ Q1..Q4 (vs Q5) + log-prop-seq + cr(time).

    ``frame`` must already be restricted to a single epoch; within-epoch
    the VOC dummies collapse, so only the time spline captures residual
    temporal trend. Returns ``(fit, df_model)`` or ``(None, df)`` if the
    frame is too sparse to fit.
    """
    df = frame.dropna(
        subset=["simd_quintile_mode", "is_singleton", "wn_prop_sequenced", "wn_mid_date"]
    ).copy()
    df = df[df["wn_prop_sequenced"] >= min_prop_seq]
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

    df_model = pd.concat(
        [df[["is_singleton"]].reset_index(drop=True), X.reset_index(drop=True)],
        axis=1,
    )
    fit = stats.logit_singleton(df_model, predictors=X.columns.tolist())
    return fit, df_model


def build_singleton_epoch_table(
    frame: pd.DataFrame,
    *,
    epochs: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Loop over epochs; return a tidy OR table for quintile dummies.

    Columns: ``estimate`` (OR vs Q5), ``conf_low``/``conf_high``,
    ``p_value``, ``z``, ``std_error``, ``quintile`` (1..4), ``epoch``.
    """
    frame = frame.copy()
    frame["epoch"] = data.assign_epoch(frame["wn_mid_date"])
    if epochs is None:
        epochs = [lbl for lbl, *_ in data.VOC_EPOCHS]

    out: list[pd.DataFrame] = []
    for epoch in epochs:
        sub = frame[frame["epoch"] == epoch]
        fit, _ = singleton_epoch_model(sub)
        if fit is None:
            continue
        tidy = stats.tidy_glm(fit)
        tidy = tidy[tidy["term"].str.startswith("q_")].copy()
        tidy["quintile"] = tidy["term"].str.replace("q_", "").astype(int)
        tidy["epoch"] = epoch
        out.append(tidy)
    return (
        pd.concat(out, ignore_index=True)
        if out
        else pd.DataFrame(
            columns=["term", "estimate", "std_error", "conf_low", "conf_high",
                     "z", "p_value", "quintile", "epoch"]
        )
    )



# ---------------------------------------------------------------------------
# 3. Joint domain decomposition (Fig. 6)
# ---------------------------------------------------------------------------


def build_domain_decomposition_table(
    frame: pd.DataFrame,
    domains: dict[str, str] = DOMAIN_ONLY,
    *,
    include_time_spline: bool = True,
) -> pd.DataFrame:
    """Joint NB fit across all seven domain ranks; return tidy share table.

    Columns
    -------
    domain, estimate (standardised log-IRR — NOT exponentiated),
    std_error, conf_low, conf_high, p_value, abs_std_coef, share.
    Shares sum to 1 across rows; rows are sorted ascending on ``share``.
    """
    required = list(domains.values()) + [
        "n_sequences", "wn_prop_sequenced", "wn_mid_date", "who_voc"
    ]
    df = frame.dropna(subset=required).copy()
    if df.empty:
        raise RuntimeError("No complete cases for joint-domain model.")

    dep_cols: list[str] = []
    for name, col in domains.items():
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