"""Main Part 1 models with hurdle and zero-truncated count components.

This script is the primary modelling pass for Part 1. It uses one Leiden
resolution to avoid treating repeated cluster resolutions as independent, then
fits:

* hurdle models for cluster size, duration, and geographic dispersion
* zero-truncated negative-binomial models for the positive count components
* linear models for observed-minus-expected within-cluster mixing

Run from the repository root with:

    conda run -n PhD python part1/main/main_analysis.py
"""

from __future__ import annotations

import argparse
import gc
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import warnings

import numpy as np
import pandas as pd
import yaml
from patsy import dmatrix
from scipy.linalg import pinvh
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.special import digamma, gammaln
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning


QC_DEFAULT = "good"
PRIMARY_RESOLUTION = 0.3
LINEAGE_MIN_CLUSTERS = 50
CALENDAR_SPLINE_DF = 8

SEQUENCE_COLUMNS = [
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
    "dz_simd_rank",
    "dz_simd_quintile",
    "dz_cum_incidence_per_capita",
    "dz_cum_prop_sequenced",
    "dz_7d_test_positivity",
]

PRIMARY_TERMS = [
    "deprivation_z",
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
]

TERM_LABELS = {
    "deprivation_z": "Mean SIMD deprivation",
    "local_incidence_z": "Local cumulative incidence",
    "local_seq_fraction_z": "Local sequencing fraction",
    "window_seq_fraction_z": "Window sequencing proportion",
    "test_positivity_z": "Local test positivity",
    "log_cluster_size_z": "Cluster size",
}

MIXING_VARIABLES = {
    "simd": {
        "column": "dz_simd_quintile",
        "label": "SIMD quintile mixing",
        "short_label": "SIMD",
    },
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
    "profile": {
        "column": "socio_demographic_profile",
        "label": "Joint SIMD-age-sex profile mixing",
        "short_label": "Joint profile",
    },
}


@dataclass(frozen=True)
class CountModelSpec:
    name: str
    label: str
    raw_outcome: str
    binary_col: str
    positive_col: str
    positive_label: str
    include_size: bool = False


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


COUNT_MODEL_SPECS = [
    CountModelSpec(
        name="cluster_size",
        label="Cluster size",
        raw_outcome="cluster_size",
        binary_col="cluster_size_gt1",
        positive_col="cluster_size_excess",
        positive_label="Additional sequences among non-singleton clusters",
    ),
    CountModelSpec(
        name="duration",
        label="Duration",
        raw_outcome="duration_days",
        binary_col="duration_gt0",
        positive_col="duration_positive_days",
        positive_label="Days among clusters lasting more than one day",
    ),
    CountModelSpec(
        name="geographic_dispersion",
        label="Geographic dispersion",
        raw_outcome="cluster_n_datazones",
        binary_col="datazones_gt1",
        positive_col="datazones_excess",
        positive_label="Additional datazones among multi-datazone clusters",
    ),
    CountModelSpec(
        name="duration_size_adjusted",
        label="Duration, size-adjusted",
        raw_outcome="duration_days",
        binary_col="duration_gt0",
        positive_col="duration_positive_days",
        positive_label="Days among clusters lasting more than one day",
        include_size=True,
    ),
    CountModelSpec(
        name="geographic_dispersion_size_adjusted",
        label="Geographic dispersion, size-adjusted",
        raw_outcome="cluster_n_datazones",
        binary_col="datazones_gt1",
        positive_col="datazones_excess",
        positive_label="Additional datazones among multi-datazone clusters",
        include_size=True,
    ),
]


def repo_root(start: Path | None = None) -> Path:
    p = (start or Path(__file__)).resolve()
    for candidate in [p, *p.parents]:
        if (candidate / "config.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate config.yaml.")


def analysis_dataset_path(root: Path) -> Path:
    with open(root / "config.yaml", "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return root / cfg["data"]["processed"]["analysis_dataset"]


def zscore(values: pd.Series) -> tuple[pd.Series, float, float]:
    mean = float(values.mean())
    sd = float(values.std(ddof=0))
    if not math.isfinite(sd) or sd == 0:
        raise ValueError(f"Cannot standardise {values.name!r}: zero or invalid SD.")
    return (values - mean) / sd, mean, sd


def logit_clipped(values: pd.Series, eps: float = 1e-5) -> pd.Series:
    clipped = values.clip(lower=eps, upper=1 - eps)
    return np.log(clipped / (1 - clipped))


def read_sequence_rows(
    path: Path,
    qc: str | None,
    primary_resolution: float,
) -> pd.DataFrame:
    filters: list[tuple[str, str, object]] = [("resolution", "==", primary_resolution)]
    if qc is not None:
        filters.append(("nextclade_qc", "==", qc))

    seq = pd.read_parquet(
        path,
        columns=SEQUENCE_COLUMNS,
        filters=filters,
        engine="pyarrow",
    )

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

    complete_profile = seq[["dz_simd_quintile", "age_band", "sex"]].notna().all(axis=1)
    seq["socio_demographic_profile"] = pd.NA
    seq.loc[complete_profile, "socio_demographic_profile"] = (
        seq.loc[complete_profile, "dz_simd_quintile"].astype(str)
        + "|"
        + seq.loc[complete_profile, "age_band"].astype(str)
        + "|"
        + seq.loc[complete_profile, "sex"].astype(str)
    )
    seq["socio_demographic_profile"] = seq["socio_demographic_profile"].astype("category")
    return seq


def pairwise_discordance_from_counts(
    counts: pd.DataFrame,
    group_cols: list[str],
) -> pd.DataFrame:
    totals = counts.groupby(group_cols, observed=True)["n"].sum().rename("n_valid")
    same_pairs = (
        counts.assign(same_pairs=counts["n"] * (counts["n"] - 1))
        .groupby(group_cols, observed=True)["same_pairs"]
        .sum()
    )
    out = pd.concat([totals, same_pairs], axis=1).reset_index()
    denom = out["n_valid"] * (out["n_valid"] - 1)
    out["discordance"] = np.where(denom > 0, 1 - out["same_pairs"] / denom, np.nan)
    return out.drop(columns=["same_pairs"])


def observed_cluster_discordance(
    seq: pd.DataFrame,
    variable: str,
    prefix: str,
) -> pd.DataFrame:
    counts = (
        seq.dropna(subset=[variable])
        .groupby(["cluster_id", variable], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = pairwise_discordance_from_counts(counts, ["cluster_id"])
    return out.rename(
        columns={
            "n_valid": f"{prefix}_n_valid",
            "discordance": f"{prefix}_discordance",
        }
    )


def expected_stratum_discordance(
    seq: pd.DataFrame,
    variable: str,
    prefix: str,
) -> pd.DataFrame:
    strata = ["window_id", "pango_lineage"]
    counts = (
        seq.dropna(subset=[variable])
        .groupby(strata + [variable], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = pairwise_discordance_from_counts(counts, strata)
    return out.rename(
        columns={
            "n_valid": f"{prefix}_stratum_n_valid",
            "discordance": f"{prefix}_expected_discordance",
        }
    )


def build_cluster_table(
    seq: pd.DataFrame,
    lineage_min_clusters: int,
    calendar_spline_df: int,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    required = [
        "cluster_id",
        "sequence_id",
        "window_id",
        "window_idx",
        "collection_date",
        "datazone",
        "pango_lineage",
        "dz_simd_rank",
        "dz_cum_incidence_per_capita",
        "dz_cum_prop_sequenced",
        "wn_prop_sequenced",
        "dz_7d_test_positivity",
    ]
    before = len(seq)
    seq = seq.dropna(subset=required).copy()
    dropped = before - len(seq)

    clusters = (
        seq.groupby("cluster_id", observed=True, sort=False)
        .agg(
            cluster_size=("sequence_id", "nunique"),
            cluster_n_datazones=("datazone", "nunique"),
            cluster_start_date=("collection_date", "min"),
            cluster_end_date=("collection_date", "max"),
            resolution=("resolution", "first"),
            window_id=("window_id", "first"),
            window_idx=("window_idx", "first"),
            wn_mid_date=("wn_mid_date", "first"),
            pango_lineage=("pango_lineage", "first"),
            mean_simd_rank=("dz_simd_rank", "mean"),
            mean_local_incidence_per_capita=("dz_cum_incidence_per_capita", "mean"),
            mean_local_seq_fraction=("dz_cum_prop_sequenced", "mean"),
            mean_window_seq_fraction=("wn_prop_sequenced", "mean"),
            mean_test_positivity=("dz_7d_test_positivity", "mean"),
        )
        .reset_index()
    )

    clusters["duration_days"] = (
        clusters["cluster_end_date"] - clusters["cluster_start_date"]
    ).dt.days.astype(int)
    clusters["cluster_size_gt1"] = (clusters["cluster_size"] > 1).astype(int)
    clusters["duration_gt0"] = (clusters["duration_days"] > 0).astype(int)
    clusters["datazones_gt1"] = (clusters["cluster_n_datazones"] > 1).astype(int)
    clusters["cluster_size_excess"] = clusters["cluster_size"] - 1
    clusters["duration_positive_days"] = clusters["duration_days"]
    clusters["datazones_excess"] = clusters["cluster_n_datazones"] - 1

    for prefix, spec in MIXING_VARIABLES.items():
        obs = observed_cluster_discordance(seq, spec["column"], prefix)
        exp = expected_stratum_discordance(seq, spec["column"], prefix)
        clusters = clusters.merge(obs, on="cluster_id", how="left")
        clusters = clusters.merge(exp, on=["window_id", "pango_lineage"], how="left")
        clusters[f"{prefix}_excess_discordance"] = (
            clusters[f"{prefix}_discordance"] - clusters[f"{prefix}_expected_discordance"]
        )

    clusters["deprivation_raw"] = -clusters["mean_simd_rank"]
    clusters["local_incidence_log"] = np.log1p(
        clusters["mean_local_incidence_per_capita"].clip(lower=0) * 1000
    )
    clusters["local_seq_fraction_logit"] = logit_clipped(clusters["mean_local_seq_fraction"])
    clusters["window_seq_fraction_logit"] = logit_clipped(clusters["mean_window_seq_fraction"])
    clusters["test_positivity_logit"] = logit_clipped(clusters["mean_test_positivity"].fillna(0))
    clusters["log_cluster_size"] = np.log(clusters["cluster_size"])

    scaling_rows = []
    transforms = {
        "deprivation_z": "deprivation_raw",
        "local_incidence_z": "local_incidence_log",
        "local_seq_fraction_z": "local_seq_fraction_logit",
        "window_seq_fraction_z": "window_seq_fraction_logit",
        "test_positivity_z": "test_positivity_logit",
        "log_cluster_size_z": "log_cluster_size",
    }
    for z_col, raw_col in transforms.items():
        clusters[z_col], mean, sd = zscore(clusters[raw_col])
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
    clusters["lineage_model"] = np.where(
        clusters["pango_lineage"].astype(str).isin(common_lineages),
        clusters["pango_lineage"].astype(str),
        "Other rare lineages",
    )

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


def lineage_levels(clusters: pd.DataFrame) -> list[str]:
    counts = clusters["lineage_model"].astype(str).value_counts()
    return counts.index.tolist()


def build_exog(
    df: pd.DataFrame,
    numeric_terms: list[str],
    calendar_cols: list[str],
    all_lineage_levels: list[str],
) -> pd.DataFrame:
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


def model_terms(spec: CountModelSpec) -> list[str]:
    terms = PRIMARY_TERMS.copy()
    if spec.include_size:
        terms.append("log_cluster_size_z")
    return terms


def extract_ratio_rows(
    result,
    term_names: list[str],
    *,
    outcome: str,
    outcome_label: str,
    component: str,
    component_label: str,
    model_family: str,
    response: str,
    n_observations: int,
    n_events: int | None,
) -> pd.DataFrame:
    names = list(result.model.exog_names)
    params = np.asarray(result.params, dtype=float)
    bse = np.asarray(result.bse, dtype=float)
    pvalues = np.asarray(result.pvalues, dtype=float)
    idx = {name: i for i, name in enumerate(names)}

    rows = []
    for term in term_names:
        if term not in idx:
            raise KeyError(f"Term {term!r} was not present in fitted model for {outcome}.")
        i = idx[term]
        coef = float(params[i])
        stderr = float(bse[i])
        rows.append(
            {
                "outcome": outcome,
                "outcome_label": outcome_label,
                "component": component,
                "component_label": component_label,
                "model_family": model_family,
                "response": response,
                "term": term,
                "term_label": TERM_LABELS[term],
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


def fit_binary_component(
    clusters: pd.DataFrame,
    spec: CountModelSpec,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
) -> tuple[pd.DataFrame, dict]:
    terms = model_terms(spec)
    if spec.include_size:
        use = clusters.loc[clusters["cluster_size"] > 1].copy()
    else:
        use = clusters.copy()
    use = use.dropna(subset=[spec.binary_col, *terms, *calendar_cols, "lineage_model"]).copy()
    y = use[spec.binary_col].astype(int)
    x = build_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use["window_id"].astype(str).to_numpy()

    model = sm.GLM(y, x, family=sm.families.Binomial())
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        result = model.fit(
            maxiter=maxiter,
            cov_type="cluster",
            cov_kwds={"groups": groups},
        )

    rows = extract_ratio_rows(
        result,
        terms,
        outcome=spec.name,
        outcome_label=spec.label,
        component="hurdle_binary",
        component_label="Probability of exceeding structural minimum",
        model_family="Binomial GLM with logit link",
        response=spec.binary_col,
        n_observations=len(use),
        n_events=int(y.sum()),
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
            "non-singleton clusters" if spec.include_size else "all primary-resolution clusters"
        ),
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


def ztnb_loglike_score(
    params: np.ndarray,
    y: np.ndarray,
    x: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Return ZTNB log likelihood, summed score, and observation-level scores."""
    beta = params[:-1]
    log_alpha = float(np.clip(params[-1], -10.0, 8.0))
    alpha = float(np.exp(log_alpha))
    r = 1.0 / alpha

    eta = x @ beta
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


def ztnb_objective(
    params: np.ndarray,
    y: np.ndarray,
    x: np.ndarray,
) -> tuple[float, np.ndarray]:
    llf, score, _ = ztnb_loglike_score(params, y, x)
    if not np.isfinite(llf) or not np.all(np.isfinite(score)):
        return 1e100, np.zeros_like(params)
    return -llf, -score


def ztnb_start_params(y: np.ndarray, x: pd.DataFrame) -> np.ndarray:
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
) -> ZTNBResult:
    x_array = np.asarray(x, dtype=float)
    start = ztnb_start_params(y, x)
    opt = minimize(
        ztnb_objective,
        start,
        args=(y, x_array),
        method="L-BFGS-B",
        jac=True,
        bounds=[(None, None)] * x_array.shape[1] + [(-10.0, 8.0)],
        options={"maxiter": maxiter, "ftol": 1e-8, "gtol": 1e-5, "maxls": 50},
    )
    params = opt.x.copy()
    params[-1] = float(np.clip(params[-1], -10.0, 8.0))
    llf, _, score_obs = ztnb_loglike_score(params, y, x_array)

    opg = score_obs.T @ score_obs
    bread_inv = pinvh(opg, rtol=1e-10)

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


def fit_positive_component(
    clusters: pd.DataFrame,
    spec: CountModelSpec,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
) -> tuple[pd.DataFrame, dict]:
    terms = model_terms(spec)
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
    rows = []
    for term in terms:
        i = idx[term]
        coef = float(result.params[i])
        stderr = float(result.bse[i])
        rows.append(
            {
                "outcome": spec.name,
                "outcome_label": spec.label,
                "component": "positive_zero_truncated_count",
                "component_label": spec.positive_label,
                "model_family": "Zero-truncated negative binomial",
                "response": spec.positive_col,
                "term": term,
                "term_label": TERM_LABELS[term],
                "coefficient": coef,
                "std_error_clustered_by_window": stderr,
                "z": coef / stderr if stderr > 0 else np.nan,
                "p_value": float(result.pvalues[i]),
                "ratio": float(np.exp(coef)),
                "ratio_ci_low": float(np.exp(coef - 1.96 * stderr)),
                "ratio_ci_high": float(np.exp(coef + 1.96 * stderr)),
                "n_observations": len(use),
                "n_events": None,
            }
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
        "converged": bool(result.converged),
        "iterations": int(result.nit),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "alpha": float(result.alpha),
        "alpha_at_upper_bound": bool(np.isclose(result.alpha, math.exp(8.0))),
        "optimizer_message": result.message,
        "warnings": "; ".join(str(w.message) for w in caught),
    }
    return pd.DataFrame(rows), diag


def fit_count_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    maxiter: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    result_frames = []
    diagnostics = []
    for spec in COUNT_MODEL_SPECS:
        if not spec.include_size:
            print(f"  - {spec.name}: hurdle binary", flush=True)
            rows, diag = fit_binary_component(
                clusters, spec, lineage_levels_all, calendar_cols, maxiter
            )
            result_frames.append(rows)
            diagnostics.append(diag)

        print(f"  - {spec.name}: zero-truncated NB positive count", flush=True)
        rows, diag = fit_positive_component(
            clusters, spec, lineage_levels_all, calendar_cols, maxiter
        )
        result_frames.append(rows)
        diagnostics.append(diag)
        gc.collect()

    return pd.concat(result_frames, ignore_index=True), pd.DataFrame(diagnostics)


def fit_mixing_models(
    clusters: pd.DataFrame,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    result_frames = []
    diagnostics = []
    terms = PRIMARY_TERMS + ["log_cluster_size_z"]
    for prefix, spec in MIXING_VARIABLES.items():
        outcome = f"{prefix}_excess_discordance"
        print(f"  - {prefix} mixing", flush=True)
        use = clusters.loc[clusters["cluster_size"] >= 2].dropna(
            subset=[outcome, *terms, *calendar_cols, "lineage_model"]
        )
        use = use.copy()
        y = use[outcome].astype(float)
        x = build_exog(use, terms, calendar_cols, lineage_levels_all)
        groups = use["window_id"].astype(str).to_numpy()

        model = sm.OLS(y, x)
        result = model.fit(cov_type="cluster", cov_kwds={"groups": groups})
        names = list(result.model.exog_names)
        params = np.asarray(result.params, dtype=float)
        bse = np.asarray(result.bse, dtype=float)
        pvalues = np.asarray(result.pvalues, dtype=float)
        idx = {name: i for i, name in enumerate(names)}

        rows = []
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


def summarise_dataset(
    seq: pd.DataFrame,
    clusters: pd.DataFrame,
    qc: str | None,
    primary_resolution: float,
    dropped: int,
    scaling: pd.DataFrame,
) -> pd.DataFrame:
    rows = [
        {"measure": "sequence_rows_used", "statistic": "count", "value": len(seq)},
        {
            "measure": "sequence_rows_dropped_missing_model_fields",
            "statistic": "count",
            "value": dropped,
        },
        {"measure": "clusters", "statistic": "count", "value": len(clusters)},
        {"measure": "primary_leiden_resolution", "statistic": "value", "value": primary_resolution},
        {"measure": "windows", "statistic": "count", "value": clusters["window_id"].nunique()},
        {
            "measure": "pango_lineages_raw",
            "statistic": "count",
            "value": clusters["pango_lineage"].nunique(),
        },
        {
            "measure": "pango_lineage_model_levels",
            "statistic": "count",
            "value": clusters["lineage_model"].nunique(),
        },
        {
            "measure": "lineage_min_clusters",
            "statistic": "value",
            "value": int(scaling.attrs["lineage_min_clusters"]),
        },
        {"measure": "qc_filter", "statistic": "value", "value": qc or "none"},
    ]

    outcomes = [
        ("cluster_size", "Cluster size", "singleton_fraction", clusters["cluster_size"] == 1),
        ("duration_days", "Duration days", "zero_fraction", clusters["duration_days"] == 0),
        (
            "cluster_n_datazones",
            "Distinct datazones",
            "single_datazone_fraction",
            clusters["cluster_n_datazones"] == 1,
        ),
    ]
    percentiles = [0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
    for col, label, structural_name, structural_mask in outcomes:
        desc = clusters[col].describe(percentiles=percentiles)
        rows.extend(
            {
                "measure": label,
                "statistic": str(stat),
                "value": float(value),
            }
            for stat, value in desc.items()
        )
        rows.append(
            {
                "measure": label,
                "statistic": structural_name,
                "value": float(structural_mask.mean()),
            }
        )

    for prefix, spec in MIXING_VARIABLES.items():
        outcome = f"{prefix}_excess_discordance"
        values = clusters.loc[clusters["cluster_size"] >= 2, outcome].dropna()
        desc = values.describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9])
        rows.extend(
            {
                "measure": f"{spec['short_label']} excess mixing",
                "statistic": str(stat),
                "value": float(value),
            }
            for stat, value in desc.items()
        )

    return pd.DataFrame(rows)


def plot_count_effects(results: pd.DataFrame, out_base: Path) -> None:
    setup_matplotlib_cache()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    primary = results[
        results["outcome"].isin(["cluster_size", "duration", "geographic_dispersion"])
        & results["term"].isin(PRIMARY_TERMS)
    ].copy()
    outcomes = ["cluster_size", "duration", "geographic_dispersion"]
    components = ["hurdle_binary", "positive_zero_truncated_count"]
    colours = {
        "deprivation_z": "#2b2b2b",
        "local_incidence_z": "#4e79a7",
        "local_seq_fraction_z": "#59a14f",
        "window_seq_fraction_z": "#f28e2b",
        "test_positivity_z": "#b07aa1",
    }

    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 8,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(3, 2, figsize=(8.0, 6.2), sharex=False)
    term_offsets = np.linspace(-0.3, 0.3, len(PRIMARY_TERMS))
    term_positions = dict(zip(PRIMARY_TERMS, term_offsets))

    for i, outcome in enumerate(outcomes):
        for j, component in enumerate(components):
            ax = axes[i, j]
            sub = primary[(primary["outcome"] == outcome) & (primary["component"] == component)]
            for _, row in sub.iterrows():
                y = term_positions[row["term"]]
                ax.plot(
                    [row["ratio_ci_low"], row["ratio_ci_high"]],
                    [y, y],
                    color=colours[row["term"]],
                    linewidth=1.1,
                    solid_capstyle="round",
                )
                ax.scatter(
                    row["ratio"],
                    y,
                    color=colours[row["term"]],
                    s=18,
                    zorder=3,
                    label=TERM_LABELS[row["term"]],
                )
            ax.axvline(1, color="#666666", linewidth=0.8, linestyle="--")
            ax.set_xscale("log")
            ax.set_yticks([])
            ax.grid(axis="x", color="#dddddd", linewidth=0.6)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if i == 0:
                ax.set_title(
                    "Hurdle: any excess" if component == "hurdle_binary" else "Positive: ZTNB",
                    fontsize=8,
                )
            if j == 0:
                ax.set_ylabel(
                    {"cluster_size": "Size", "duration": "Duration", "geographic_dispersion": "Datazones"}[
                        outcome
                    ],
                    rotation=0,
                    ha="right",
                    va="center",
                    labelpad=32,
                )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    fig.legend(
        unique.values(),
        unique.keys(),
        loc="lower center",
        ncol=2,
        frameon=False,
        columnspacing=1.4,
        handlelength=1.2,
    )
    fig.supxlabel("Adjusted ratio per 1 SD higher cluster-level covariate", y=0.08, fontsize=8)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.93, bottom=0.24, hspace=0.42, wspace=0.22)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_mixing_effects(results: pd.DataFrame, out_base: Path) -> None:
    setup_matplotlib_cache()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 8,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    models = ["simd", "age", "sex", "profile"]
    model_positions = {model: i for i, model in enumerate(models)}
    terms = PRIMARY_TERMS + ["log_cluster_size_z"]
    term_offsets = np.linspace(-0.32, 0.32, len(terms))
    term_positions = dict(zip(terms, term_offsets))
    colours = {
        "deprivation_z": "#2b2b2b",
        "local_incidence_z": "#4e79a7",
        "local_seq_fraction_z": "#59a14f",
        "window_seq_fraction_z": "#f28e2b",
        "test_positivity_z": "#b07aa1",
        "log_cluster_size_z": "#7f7f7f",
    }

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    for _, row in results.iterrows():
        y = model_positions[row["outcome"]] + term_positions[row["term"]]
        ax.plot(
            [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
            [y, y],
            color=colours[row["term"]],
            linewidth=1.2,
            solid_capstyle="round",
        )
        ax.scatter(
            row["coefficient_percentage_points"],
            y,
            color=colours[row["term"]],
            s=18,
            zorder=3,
            label=TERM_LABELS[row["term"]],
        )

    ax.axvline(0, color="#666666", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Change in excess pairwise discordance, percentage points per 1 SD higher covariate")
    ax.set_yticks(list(model_positions.values()))
    ax.set_yticklabels([MIXING_VARIABLES[m]["short_label"] for m in models])
    ax.set_ylim(-0.6, len(models) - 0.4)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(
        unique.values(),
        unique.keys(),
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=2,
        frameon=False,
        columnspacing=1.4,
        handlelength=1.2,
    )
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(bottom=0.34, left=0.2, right=0.98)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def setup_matplotlib_cache() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)


def run(
    root: Path,
    qc: str | None,
    primary_resolution: float,
    lineage_min_clusters: int,
    calendar_spline_df: int,
    maxiter: int,
) -> None:
    out_dir = root / "part1" / "main"
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    cache_dir = out_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    path = analysis_dataset_path(root)
    print(f"Reading primary-resolution sequence rows from {path}", flush=True)
    seq = read_sequence_rows(path, qc=qc, primary_resolution=primary_resolution)
    print(f"Building cluster table from {len(seq):,} sequence rows", flush=True)
    clusters, scaling, dropped = build_cluster_table(
        seq,
        lineage_min_clusters=lineage_min_clusters,
        calendar_spline_df=calendar_spline_df,
    )
    calendar_cols = [col for col in clusters.columns if col.startswith("calendar_spline_")]
    lineage_levels_all = lineage_levels(clusters)

    clusters.to_parquet(cache_dir / "main_cluster_table.parquet", index=False)
    scaling.to_csv(tables_dir / "main_covariate_scaling.csv", index=False)
    descriptives = summarise_dataset(seq, clusters, qc, primary_resolution, dropped, scaling)
    descriptives.to_csv(tables_dir / "main_dataset_descriptives.csv", index=False)
    print(
        f"Fitting count models for {len(clusters):,} clusters, "
        f"{len(lineage_levels_all)} lineage model levels, {len(calendar_cols)} calendar spline terms",
        flush=True,
    )
    count_results, count_diagnostics = fit_count_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
        maxiter=maxiter,
    )
    count_results.to_csv(tables_dir / "main_hurdle_count_model_results.csv", index=False)
    count_diagnostics.to_csv(tables_dir / "main_hurdle_count_model_diagnostics.csv", index=False)

    print("Fitting mixing models", flush=True)
    mixing_results, mixing_diagnostics = fit_mixing_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
    )
    mixing_results.to_csv(tables_dir / "main_mixing_model_results.csv", index=False)
    mixing_diagnostics.to_csv(tables_dir / "main_mixing_model_diagnostics.csv", index=False)

    plot_count_effects(count_results, figures_dir / "main_hurdle_count_effects")
    plot_mixing_effects(mixing_results, figures_dir / "main_mixing_effects")

    print(f"Wrote {tables_dir / 'main_hurdle_count_model_results.csv'}", flush=True)
    print(f"Wrote {tables_dir / 'main_mixing_model_results.csv'}", flush=True)
    print(f"Wrote {figures_dir / 'main_hurdle_count_effects.png'}", flush=True)
    print(f"Wrote {figures_dir / 'main_mixing_effects.png'}", flush=True)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument(
        "--qc",
        default=QC_DEFAULT,
        help="Nextclade QC status to retain. Use 'none' to disable QC filtering.",
    )
    parser.add_argument("--primary-resolution", type=float, default=PRIMARY_RESOLUTION)
    parser.add_argument("--lineage-min-clusters", type=int, default=LINEAGE_MIN_CLUSTERS)
    parser.add_argument("--calendar-spline-df", type=int, default=CALENDAR_SPLINE_DF)
    parser.add_argument("--maxiter", type=int, default=1000)
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
    )


if __name__ == "__main__":
    main()
