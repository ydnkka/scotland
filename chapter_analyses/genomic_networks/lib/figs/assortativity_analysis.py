"""Shared assortativity analysis helpers for Chapter 4 figures and tables."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from common import ATTRIBUTE_ORDER, Paths, read_table  # noqa: E402
from utils.data import Paths as DataPaths  # noqa: E402
from utils.policy import window_policy_lookup  # noqa: E402


MIN_POOL_EDGE_CONTRIBUTIONS = 35
VARIANCE_WINDSORISE_VALUES = tuple(list(range(0, 100, 5)) + [99])
VARIANCE_REFERENCE_WINDSORISE = 90


def _weighted_mean(y: np.ndarray, v: np.ndarray, tau2: float) -> tuple[float, np.ndarray]:
    w = 1.0 / (v + tau2)
    return float(np.sum(w * y) / np.sum(w)), w


def _tau2_dl(y: np.ndarray, v: np.ndarray) -> float:
    """DerSimonian-Laird (closed form)."""
    w = 1.0 / v
    sw = w.sum()
    mu = np.sum(w * y) / sw
    q = np.sum(w * (y - mu) ** 2)
    df = len(y) - 1
    c = sw - np.sum(w**2) / sw
    return max(0.0, (q - df) / c) if c > 0 else 0.0


def _tau2_pm(y: np.ndarray, v: np.ndarray, tol: float = 1e-8, max_iter: int = 200) -> float:
    """
    Paule-Mandel: find tau2 s.t. sum w_i (y_i - mu)^2 = k - 1,
    with w_i = 1/(v_i + tau2). Solved by bisection.
    """
    df = len(y) - 1

    def F(tau2: float) -> float:
        mu, w = _weighted_mean(y, v, tau2)
        return float(np.sum(w * (y - mu) ** 2) - df)

    if F(0.0) <= 0:
        return 0.0

    lo, hi = 0.0, max(float(v.max()), 1.0)
    it = 0
    while F(hi) > 0 and it < 60:
        hi *= 2.0
        it += 1

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f = F(mid)
        if abs(f) < tol or (hi - lo) < tol:
            return max(0.0, mid)
        if f > 0:
            lo = mid
        else:
            hi = mid
    return max(0.0, 0.5 * (lo + hi))


def _tau2_reml(y: np.ndarray, v: np.ndarray, tol: float = 1e-8, max_iter: int = 500) -> float:
    """REML via fixed-point iteration (Viechtbauer 2005)."""
    tau2 = max(_tau2_dl(y, v), 0.0)
    for _ in range(max_iter):
        w = 1.0 / (v + tau2)
        sw = w.sum()
        mu = np.sum(w * y) / sw
        num = np.sum(w**2 * ((y - mu) ** 2 - v)) + (1.0 / sw)
        den = np.sum(w**2)
        tau2_new = max(0.0, num / den)
        if abs(tau2_new - tau2) < tol:
            tau2 = tau2_new
            break
        tau2 = tau2_new
    return tau2


_TAU2_ESTIMATORS = {"DL": _tau2_dl, "PM": _tau2_pm, "REML": _tau2_reml}


def _weighted_mean_ci_from_se(
    values: pd.Series,
    weights: pd.Series,
    standard_errors: pd.Series,
) -> dict[str, float]:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return {
            "weighted_mean": np.nan,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "ci_weight_share": np.nan,
        }

    values = values.loc[mask].astype(float)
    weights = weights.loc[mask].astype(float)
    weighted_mean_value = float(np.average(values, weights=weights))

    se_mask = standard_errors.loc[mask].notna()
    if not se_mask.any():
        return {
            "weighted_mean": weighted_mean_value,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "ci_weight_share": np.nan,
        }

    ci_weights = weights.loc[se_mask]
    ci_standard_errors = standard_errors.loc[mask].loc[se_mask].astype(float)
    normalized = ci_weights / weights.sum()
    combined_se = float(np.sqrt(np.sum((normalized * ci_standard_errors) ** 2)))
    return {
        "weighted_mean": weighted_mean_value,
        "combined_se": combined_se,
        "ci_low": weighted_mean_value - 1.96 * combined_se,
        "ci_high": weighted_mean_value + 1.96 * combined_se,
        "ci_weight_share": float(ci_weights.sum() / weights.sum()),
    }


def compatibility_window_lookup(windows: object | None = None) -> pd.DataFrame:
    """Return window metadata with policy labels and analysis-window bounds."""
    lookup = window_policy_lookup(windows)
    data_paths = DataPaths.from_config()
    bounds = pd.read_parquet(
        data_paths.analysis_dataset,
        columns=["window_id", "window_idx", "wn_start_date", "wn_end_date"],
    )
    bounds["window_id"] = bounds["window_id"].astype(str)
    bounds["wn_start_date"] = pd.to_datetime(bounds["wn_start_date"], errors="coerce")
    bounds["wn_end_date"] = pd.to_datetime(bounds["wn_end_date"], errors="coerce")
    bounds = bounds.drop_duplicates(["window_id", "window_idx"])
    lookup = lookup.merge(bounds, on=["window_id", "window_idx"], how="left")
    return lookup.sort_values("window_idx").reset_index(drop=True)


def _overlap_based_covariance(
    window_lookup: pd.DataFrame,
    standard_errors: pd.Series,
) -> np.ndarray:
    """Build a working covariance matrix from overlap between rolling windows."""
    work = window_lookup[["wn_start_date", "wn_end_date"]].copy()
    starts = pd.to_datetime(work["wn_start_date"], errors="coerce")
    ends = pd.to_datetime(work["wn_end_date"], errors="coerce")
    if starts.isna().any() or ends.isna().any():
        raise ValueError("Window start/end dates are required for correlated GLS.")

    starts_ns = starts.to_numpy(dtype="datetime64[ns]")
    ends_ns = ends.to_numpy(dtype="datetime64[ns]")
    widths = ((ends_ns - starts_ns) / np.timedelta64(1, "D")) + 1.0

    overlap = np.minimum(ends_ns[:, None], ends_ns[None, :]) - np.maximum(
        starts_ns[:, None], starts_ns[None, :]
    )
    overlap = (overlap / np.timedelta64(1, "D")) + 1.0
    overlap = np.clip(np.asarray(overlap, dtype=float), 0.0, None)

    denom = np.sqrt(np.outer(widths, widths))
    rho = np.divide(overlap, denom, out=np.zeros_like(overlap), where=denom > 0)
    np.fill_diagonal(rho, 1.0)

    se = standard_errors.to_numpy(dtype=float)
    cov = rho * np.outer(se, se)
    diag_mean = float(np.nanmean(np.diag(cov))) if cov.size else 0.0
    ridge = max(1e-12, 1e-8 * diag_mean) if np.isfinite(diag_mean) and diag_mean > 0 else 1e-12
    return cov + np.eye(len(cov)) * ridge


def _correlated_gls_summary(
    work: pd.DataFrame,
    *,
    conf_level: float = 0.95,
) -> dict[str, float]:
    """Estimate a correlated GLS intercept using overlap-based covariance."""
    cols = ["window_idx", "pooled_mean", "pooled_se", "wn_start_date", "wn_end_date"]
    clean = work[cols].replace([np.inf, -np.inf], np.nan).dropna().copy()
    clean = clean.loc[clean["pooled_se"].gt(0)]
    if clean.empty:
        return {
            "gls_mean": np.nan,
            "gls_se": np.nan,
            "gls_ci_low": np.nan,
            "gls_ci_high": np.nan,
        }

    y = clean["pooled_mean"].astype(float).to_numpy()
    se = clean["pooled_se"].astype(float)
    if len(clean) == 1:
        est = float(y[0])
        se_est = float(se.iloc[0])
    else:
        cov = _overlap_based_covariance(clean, se)
        ones = np.ones(len(clean), dtype=float)
        try:
            v_1 = np.linalg.solve(cov, ones)
            v_y = np.linalg.solve(cov, y)
        except np.linalg.LinAlgError:
            cov_inv = np.linalg.pinv(cov, rcond=1e-10)
            v_1 = cov_inv @ ones
            v_y = cov_inv @ y
        denom = float(ones @ v_1)
        if not np.isfinite(denom) or denom <= 0:
            weights = 1.0 / se.to_numpy(dtype=float) ** 2
            est = float(np.average(y, weights=weights))
            se_est = float(np.sqrt(1.0 / weights.sum()))
        else:
            est = float((ones @ v_y) / denom)
            se_est = float(np.sqrt(1.0 / denom))

    crit = float(stats.norm.ppf(0.5 + conf_level / 2.0))
    return {
        "gls_mean": est,
        "gls_se": se_est,
        "gls_ci_low": est - crit * se_est,
        "gls_ci_high": est + crit * se_est,
    }


def random_effects_pool(
    data: pd.DataFrame,
    est_col: str = "assortativity",
    se_col: str = "assortativity_se",
    method: str = "REML",  # "DL", "PM", or "REML"
    bounds: tuple[float, float] | None = (-1.0, 1.0),
    use_hksj: bool = True,
    conf_level: float = 0.95,
) -> pd.DataFrame:
    """
    Random-effects meta-analysis with selectable tau^2 estimator,
    HKSJ small-sample CIs, and a prediction interval.
    Returns a one-row DataFrame so grouping stays stable across pandas versions.
    """
    if method not in _TAU2_ESTIMATORS:
        raise ValueError(f"method must be one of {list(_TAU2_ESTIMATORS)}")

    cols = [
        "n_lineages",
        "method",
        "pooled_mean",
        "pooled_se",
        "pooled_ci_low",
        "pooled_ci_high",
        "tau2",
        "i2",
        "q",
        "q_pvalue",
        "pred_int_low",
        "pred_int_high",
        "q05",
        "q25",
        "median",
        "q75",
        "q95",
    ]

    def empty(k: int) -> pd.DataFrame:
        out = {c: np.nan for c in cols}
        out["n_lineages"] = k
        out["method"] = method # type: ignore
        return pd.DataFrame([out])[cols]

    d = data[[est_col, se_col]].copy()
    d = d.replace([np.inf, -np.inf], np.nan).dropna()
    d = d[d[se_col] > 0]

    k = len(d)
    if k < 2:
        return empty(k)

    y = d[est_col].to_numpy(dtype=float)
    se = d[se_col].to_numpy(dtype=float)
    v = se**2
    alpha = 1.0 - conf_level

    w_fe = 1.0 / v
    mu_fixed = np.sum(w_fe * y) / w_fe.sum()
    q = np.sum(w_fe * (y - mu_fixed) ** 2)
    df_q = k - 1
    q_pvalue = stats.chi2.sf(q, df_q)
    i2 = max(0.0, (q - df_q) / q) if q > 0 else 0.0

    tau2 = _TAU2_ESTIMATORS[method](y, v)

    w_re = 1.0 / (v + tau2)
    sw_re = w_re.sum()
    pooled_mean = np.sum(w_re * y) / sw_re
    var_wald = 1.0 / sw_re

    if use_hksj and k >= 2:
        q_hk = np.sum(w_re * (y - pooled_mean) ** 2) / df_q
        var_ci = max(q_hk / sw_re, var_wald)
        crit = stats.t.ppf(1 - alpha / 2, df=df_q)
    else:
        var_ci = var_wald
        crit = stats.norm.ppf(1 - alpha / 2)

    pooled_se = np.sqrt(var_ci)
    ci_low = pooled_mean - crit * pooled_se
    ci_high = pooled_mean + crit * pooled_se

    if k >= 3:
        t_pred = stats.t.ppf(1 - alpha / 2, df=k - 2)
        pred_se = np.sqrt(tau2 + var_wald)
        pred_low = pooled_mean - t_pred * pred_se
        pred_high = pooled_mean + t_pred * pred_se
    else:
        pred_low = pred_high = np.nan

    if bounds is not None:
        lo, hi = bounds
        ci_low, ci_high = max(lo, ci_low), min(hi, ci_high)
        if np.isfinite(pred_low):
            pred_low = max(lo, pred_low)
        if np.isfinite(pred_high):
            pred_high = min(hi, pred_high)

    out = {
        "n_lineages": k,
        "method": method,
        "pooled_mean": pooled_mean,
        "pooled_se": pooled_se,
        "pooled_ci_low": ci_low,
        "pooled_ci_high": ci_high,
        "tau2": tau2,
        "i2": i2,
        "q": q,
        "q_pvalue": q_pvalue,
        "pred_int_low": pred_low,
        "pred_int_high": pred_high,
        "q05": np.quantile(y, 0.05),
        "q25": np.quantile(y, 0.25),
        "median": np.median(y),
        "q75": np.quantile(y, 0.75),
        "q95": np.quantile(y, 0.95),
    }
    return pd.DataFrame([out])[cols]


def compatibility_pooled_window_inputs(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and filter the bootstrap assortativity table used in the notebook."""
    ass_raw = read_table(paths, "compatibility_assortativity_bootstrap")
    windows = compatibility_window_lookup(ass_raw["window_id"].unique())
    ass_raw = ass_raw.merge(windows, on="window_id", how="left")
    ass = ass_raw.dropna(subset=["assortativity"])
    ass = ass.loc[ass["n_edges_used"].ge(MIN_POOL_EDGE_CONTRIBUTIONS)].copy()
    return ass, windows.sort_values("window_idx").copy()


def compatibility_window_pooled_meta(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the pooled window meta-analysis table and the window lookup."""
    ass, windows = compatibility_pooled_window_inputs(paths)
    window_meta = (
        ass.groupby(["attribute", "attribute_label", "window_idx"], observed=True)
        .apply(random_effects_pool)
        .reset_index(level=-1, drop=True)
        .reset_index()
    )
    return window_meta, windows


def pooled_window_attribute_summary(window_meta: pd.DataFrame) -> pd.DataFrame:
    """Summarise pooled window estimates by attribute using correlated GLS."""
    window_lookup = compatibility_window_lookup(window_meta["window_idx"].dropna().unique())
    rows: list[dict[str, object]] = []
    attribute_labels = [
        label
        for label in ATTRIBUTE_ORDER
        if label in window_meta["attribute_label"].dropna().unique()
    ]
    for label in attribute_labels:
        group = window_meta.loc[window_meta["attribute_label"].eq(label)].copy()
        if group.empty:
            continue

        merged = group.merge(
            window_lookup[["window_idx", "wn_start_date", "wn_end_date"]],
            on="window_idx",
            how="left",
        )
        gls = _correlated_gls_summary(merged)
        pooled_mean = pd.to_numeric(group["pooled_mean"], errors="coerce").replace(
            [np.inf, -np.inf], np.nan
        )
        clean = pooled_mean.dropna()
        rows.append(
            {
                "attribute": group["attribute"].iloc[0],
                "attribute_label": label,
                "method": group["method"].iloc[0] if "method" in group.columns else "REML",
                "n_windows": int(group["window_idx"].nunique()),
                "n_estimated_windows": int(pooled_mean.notna().sum()),
                "median_n_lineages": float(pd.to_numeric(group["n_lineages"], errors="coerce").median()),
                "gls_mean": gls["gls_mean"],
                "gls_se": gls["gls_se"],
                "gls_ci_low": gls["gls_ci_low"],
                "gls_ci_high": gls["gls_ci_high"],
                "window_median": float(clean.median()) if not clean.empty else np.nan,
                "window_q10": float(clean.quantile(0.10)) if not clean.empty else np.nan,
                "window_q90": float(clean.quantile(0.90)) if not clean.empty else np.nan,
                "median_tau2": float(pd.to_numeric(group["tau2"], errors="coerce").median()),
                "median_i2": float(pd.to_numeric(group["i2"], errors="coerce").median()),
            }
        )

    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary

    summary["attribute_label"] = pd.Categorical(
        summary["attribute_label"], categories=ATTRIBUTE_ORDER, ordered=True
    )
    return summary.sort_values("attribute_label").reset_index(drop=True)


def variance_decomposition_additive_model(
    data: pd.DataFrame,
    min_rows: int = 5,
    group_keys: tuple[str, str] = ("window_id", "pango_lineage"),
    winsorize_weights: int | None = None,
) -> pd.Series:
    """
    Weighted (inverse-variance) variance decomposition of `assortativity`
    across `group_keys`.
    """
    if len(group_keys) != 2:
        raise ValueError(f"group_keys must contain exactly 2 columns, got {group_keys}")

    k0, k1 = group_keys
    ci_width: pd.Series | None = None
    if "ci_width" in data.columns:
        ci_width = pd.to_numeric(data["ci_width"], errors="coerce")
    elif {"assortativity_ci_low", "assortativity_ci_high"}.issubset(data.columns):
        ci_width = pd.to_numeric(data["assortativity_ci_high"], errors="coerce") - pd.to_numeric(
            data["assortativity_ci_low"], errors="coerce"
        )

    base_cols = list(group_keys) + ["assortativity", "assortativity_se"]
    missing = [c for c in base_cols if c not in data.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    d = data[base_cols].copy()
    if ci_width is not None:
        d["ci_width"] = ci_width

    d = d.dropna(subset=base_cols)
    d["assortativity"] = d["assortativity"].astype(float)
    d["assortativity_se"] = d["assortativity_se"].astype(float)
    d = d[
        np.isfinite(d["assortativity"])
        & np.isfinite(d["assortativity_se"])
        & (d["assortativity_se"] > 0)
    ]

    if len(d) < min_rows:
        return pd.Series(dtype=float)

    d[k0] = d[k0].astype("category")
    d[k1] = d[k1].astype("category")

    n = len(d)
    n_windows = d[k0].nunique()
    n_lineages = d[k1].nunique()
    if n_windows < 2 or n_lineages < 2:
        return pd.Series(dtype=float)

    w = 1.0 / d["assortativity_se"] ** 2
    weight_cap = np.nan
    n_weights_capped = 0
    if winsorize_weights is not None:
        weight_cap = np.nanpercentile(w, winsorize_weights)
        n_weights_capped = int((w > weight_cap).sum())
        w = np.minimum(w, weight_cap)
    d["w"] = w

    def _fit(formula: str):
        return smf.wls(formula, data=d, weights=d["w"]).fit()

    m0 = _fit("assortativity ~ 1")
    m_window = _fit(f"assortativity ~ C({k0})")
    m_lineage = _fit(f"assortativity ~ C({k1})")
    m_both = _fit(f"assortativity ~ C({k0}) + C({k1})")

    ss_total = m0.ssr
    if not np.isfinite(ss_total) or ss_total <= 0:
        return pd.Series(dtype=float)

    window_alone = m0.ssr - m_window.ssr
    lineage_alone = m0.ssr - m_lineage.ssr
    additive_explained = m0.ssr - m_both.ssr

    lineage_given_window = m_window.ssr - m_both.ssr
    window_given_lineage = m_lineage.ssr - m_both.ssr
    residual = m_both.ssr

    df_model_both = m_both.df_model
    df_resid_both = m_both.df_resid
    if df_resid_both > 0:
        adj_additive_fraction = 1.0 - (residual / df_resid_both) / (ss_total / (n - 1))
    else:
        adj_additive_fraction = np.nan

    return pd.Series(
        {
            "n": n,
            "n_windows": n_windows,
            "n_lineages": n_lineages,
            "df_model_additive": df_model_both,
            "df_resid_additive": df_resid_both,
            "weighted_total_ss": ss_total,
            "window_alone_fraction": window_alone / ss_total,
            "lineage_alone_fraction": lineage_alone / ss_total,
            "additive_model_fraction": additive_explained / ss_total,
            "adj_additive_model_fraction": adj_additive_fraction,
            "residual_fraction": residual / ss_total,
            "lineage_given_window_fraction": lineage_given_window / ss_total,
            "window_given_lineage_fraction": window_given_lineage / ss_total,
            "weight_cap": weight_cap,
            "n_weights_capped": n_weights_capped,
            "median_boot_se": d["assortativity_se"].median(),
            "median_ci_width": d["ci_width"].median() if "ci_width" in d.columns else np.nan,
        }
    )


def compatibility_variance_decomposition_long(
    paths: Paths,
    *,
    winsorize_values: tuple[int, ...] = VARIANCE_WINDSORISE_VALUES,
) -> pd.DataFrame:
    """Compute the full winsorisation sensitivity table from the notebook."""
    ass, _ = compatibility_pooled_window_inputs(paths)
    vd_results: dict[int, pd.DataFrame] = {}
    for winsorize in winsorize_values:
        vd = (
            ass.groupby(["attribute", "attribute_label"], observed=True)
            .apply(variance_decomposition_additive_model, winsorize_weights=winsorize)
            .reset_index()
        )
        vd["winsorize"] = winsorize
        vd_results[winsorize] = vd
    if not vd_results:
        return pd.DataFrame()
    return pd.concat(vd_results, ignore_index=True)


def variance_decomposition_summary(
    vd_long: pd.DataFrame,
    *,
    winsorize: int = VARIANCE_REFERENCE_WINDSORISE,
) -> pd.DataFrame:
    """Summarise a single winsorisation level for LaTeX tables."""
    work = vd_long.loc[vd_long["winsorize"].eq(winsorize)].copy()
    if work.empty:
        raise ValueError(f"No variance decomposition rows found for winsorize={winsorize}.")

    rows: list[dict[str, object]] = []
    attribute_labels = [label for label in ATTRIBUTE_ORDER if label in work["attribute_label"].dropna().unique()]
    for label in attribute_labels:
        group = work.loc[work["attribute_label"].eq(label)].copy()
        if group.empty:
            continue
        row = group.iloc[0]
        rows.append(
            {
                "attribute": row["attribute"],
                "attribute_label": label,
                "n": row["n"],
                "n_windows": row["n_windows"],
                "n_lineages": row["n_lineages"],
                "df_model_additive": row["df_model_additive"],
                "df_resid_additive": row["df_resid_additive"],
                "weighted_total_ss": row["weighted_total_ss"],
                "window_alone_fraction": row["window_alone_fraction"],
                "lineage_alone_fraction": row["lineage_alone_fraction"],
                "additive_model_fraction": row["additive_model_fraction"],
                "adj_additive_model_fraction": row["adj_additive_model_fraction"],
                "residual_fraction": row["residual_fraction"],
                "lineage_given_window_fraction": row["lineage_given_window_fraction"],
                "window_given_lineage_fraction": row["window_given_lineage_fraction"],
                "weight_cap": row["weight_cap"],
                "n_weights_capped": row["n_weights_capped"],
                "median_boot_se": row["median_boot_se"],
                "median_ci_width": row["median_ci_width"],
            }
        )

    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary

    summary["attribute_label"] = pd.Categorical(
        summary["attribute_label"], categories=ATTRIBUTE_ORDER, ordered=True
    )
    return summary.sort_values("attribute_label").reset_index(drop=True)
