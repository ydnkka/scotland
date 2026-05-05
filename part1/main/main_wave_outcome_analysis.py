"""Wave-specific main-formulation cluster outcome models.

Fits the main hurdle/ZTNB count models separately within epidemic wave groups
for cluster size, duration, and geographic spread. The figure generated from
these tables focuses on the wave-specific SIMD-deprivation coefficient; the
tables retain all main covariates.
"""

from __future__ import annotations

import argparse
import gc
import math
from pathlib import Path
from typing import Iterable
import warnings

import numpy as np
import pandas as pd
from scipy.linalg import pinvh
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning

try:
    from .main_analysis import (
        COUNT_MODEL_SPECS,
        PRIMARY_TERMS,
        fit_ztnb,
        lineage_levels,
        repo_root,
    )
    from .main_domain_wave_analysis import WAVE_LABELS, WAVE_ORDER, assign_wave
except ImportError:
    from main_analysis import (
        COUNT_MODEL_SPECS,
        PRIMARY_TERMS,
        fit_ztnb,
        lineage_levels,
        repo_root,
    )
    from main_domain_wave_analysis import WAVE_LABELS, WAVE_ORDER, assign_wave


TERM_LABELS = {
    "deprivation_z": "Mean SIMD deprivation",
    "local_incidence_z": "Local cumulative incidence",
    "local_seq_fraction_z": "Local sequencing fraction",
    "window_seq_fraction_z": "Window sequencing proportion",
    "test_positivity_z": "Local test positivity",
}


def drop_redundant_columns(x: pd.DataFrame, tol: float = 1e-8) -> pd.DataFrame:
    """Keep a full-rank column set, preserving earlier columns first.

    ``build_exog`` orders columns as intercept, primary terms, calendar spline,
    then lineage dummies. A sequential orthogonalisation therefore preserves the
    substantive covariates unless they are themselves redundant within a wave.
    """
    keep: list[str] = []
    basis: list[np.ndarray] = []
    n = max(len(x), 1)
    base_scale = math.sqrt(n)
    for col in x.columns:
        values = np.asarray(x[col], dtype=float)
        if not np.all(np.isfinite(values)):
            continue
        norm = float(np.linalg.norm(values))
        if norm <= tol:
            continue
        residual = values.copy()
        for q in basis:
            residual -= float(np.dot(q, residual)) * q
        residual_norm = float(np.linalg.norm(residual))
        if residual_norm > tol * max(norm, base_scale):
            keep.append(col)
            basis.append(residual / residual_norm)
    return x[keep]


def build_wave_exog(
    df: pd.DataFrame,
    terms: list[str],
    calendar_cols: list[str],
    lineage_levels_wave: list[str],
) -> pd.DataFrame:
    """Build the design matrix for a wave-specific model.

    Lineage dummies are included when ``lineage_levels_wave`` is non-empty.
    ``drop_redundant_columns`` handles rank deficiency caused by sparse
    lineages within a wave, so no manual lineage pre-filtering is needed.
    The column ordering (intercept → covariates → calendar → lineage) means
    ``drop_redundant_columns`` always preserves the substantive terms first.
    """
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


def extract_ratio_rows(
    *,
    params: np.ndarray,
    bse: np.ndarray,
    pvalues: np.ndarray,
    exog_names: list[str],
    terms: list[str],
    wave: str,
    spec,
    component: str,
    component_label: str,
    model_family: str,
    response: str,
    n_observations: int,
    n_events: int | None,
) -> pd.DataFrame:
    idx = {name: i for i, name in enumerate(exog_names)}
    rows = []
    for term in terms:
        if term not in idx:
            continue
        i = idx[term]
        coef = float(params[i])
        stderr = float(bse[i])
        rows.append(
            {
                "wave_group": wave,
                "wave_label": WAVE_LABELS.get(wave, wave),
                "outcome": spec.name,
                "outcome_label": spec.label,
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


def skipped_diag(wave: str, spec, component: str, reason: str, n_observations: int, n_windows: int) -> dict:
    return {
        "wave_group": wave,
        "wave_label": WAVE_LABELS.get(wave, wave),
        "outcome": spec.name,
        "component": component,
        "skipped": True,
        "reason": reason,
        "n_observations": int(n_observations),
        "n_windows": int(n_windows),
    }


def clustered_logit_inference(
    result,
    y: pd.Series,
    x: pd.DataFrame,
    groups: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Cluster-robust logit standard errors using a pseudo-inverse bread.

    Statsmodels' clustered covariance path can fail with a singular Hessian in
    wave-stratified models. The fitted coefficients are still usable, so compute
    the sandwich covariance directly with a generalized inverse.
    """
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    mu = np.asarray(result.fittedvalues, dtype=float)
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
    params = np.asarray(result.params, dtype=float)
    z_values = np.divide(params, bse, out=np.full_like(params, np.nan), where=bse > 0)
    pvalues = 2 * norm.sf(np.abs(z_values))
    return bse, pvalues


def fit_wave_binary_component(
    wave_df: pd.DataFrame,
    wave: str,
    spec,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    *,
    maxiter: int,
    min_events: int,
) -> tuple[pd.DataFrame, dict]:
    terms = PRIMARY_TERMS.copy()
    use = wave_df.dropna(subset=[spec.binary_col, *terms, *calendar_cols, "lineage_model"]).copy()
    y = use[spec.binary_col].astype(int)
    n_events = int(y.sum())
    n_nonevents = int(len(y) - n_events)
    n_windows = int(use["window_id"].nunique())
    if n_events < min_events or n_nonevents < min_events:
        return pd.DataFrame(), skipped_diag(
            wave,
            spec,
            "hurdle_binary",
            "below minimum events/non-events",
            len(use),
            n_windows,
        )

    x = build_wave_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use["window_id"].astype(str).to_numpy()
    model = sm.GLM(y, x, family=sm.families.Binomial())
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        try:
            result = model.fit(maxiter=maxiter)
        except Exception as exc:
            return pd.DataFrame(), skipped_diag(
                wave,
                spec,
                "hurdle_binary",
                f"fit failed: {exc}",
                len(use),
                n_windows,
            )

    bse, pvalues = clustered_logit_inference(result, y, x, groups)

    rows = extract_ratio_rows(
        params=np.asarray(result.params, dtype=float),
        bse=bse,
        pvalues=pvalues,
        exog_names=list(result.model.exog_names),
        terms=terms,
        wave=wave,
        spec=spec,
        component="hurdle_binary",
        component_label="Probability of exceeding structural minimum",
        model_family="Binomial GLM with logit link",
        response=spec.binary_col,
        n_observations=int(len(use)),
        n_events=n_events,
    )
    diag = {
        "wave_group": wave,
        "wave_label": WAVE_LABELS.get(wave, wave),
        "outcome": spec.name,
        "component": "hurdle_binary",
        "skipped": False,
        "reason": "",
        "model_family": "Binomial GLM with logit link",
        "response": spec.binary_col,
        "n_observations": int(len(use)),
        "n_events": n_events,
        "event_fraction": float(y.mean()),
        "n_features": int(x.shape[1]),
        "n_lineage_levels_available": int(len(lineage_levels_all)),
        "n_lineage_terms_used": int(sum(col.startswith("lineage_") for col in x.columns)),
        "lineage_adjustment": "wave-stratified; lineage dummies included, rank-dropped if collinear",
        "n_windows": n_windows,
        "covariance_method": "window-clustered sandwich with pseudo-inverse bread",
        "converged": bool(getattr(result, "converged", False)),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "warnings": "; ".join(str(w.message) for w in caught),
    }
    return rows, diag


def fit_wave_positive_component(
    wave_df: pd.DataFrame,
    wave: str,
    spec,
    lineage_levels_all: list[str],
    calendar_cols: list[str],
    *,
    maxiter: int,
    min_positive: int,
    min_windows: int,
) -> tuple[pd.DataFrame, dict]:
    terms = PRIMARY_TERMS.copy()
    use = wave_df.loc[wave_df[spec.positive_col] > 0].dropna(
        subset=[spec.positive_col, *terms, *calendar_cols, "lineage_model"]
    )
    use = use.copy()
    n_windows = int(use["window_id"].nunique())
    if len(use) < min_positive or n_windows < min_windows:
        return pd.DataFrame(), skipped_diag(
            wave,
            spec,
            "positive_zero_truncated_count",
            "below minimum positive clusters/windows",
            len(use),
            n_windows,
        )

    y = use[spec.positive_col].astype(int).to_numpy()
    x = build_wave_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use["window_id"].astype(str).to_numpy()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fit_ztnb(y, x, groups, maxiter=maxiter)

    rows = extract_ratio_rows(
        params=np.asarray(result.params, dtype=float),
        bse=np.asarray(result.bse, dtype=float),
        pvalues=np.asarray(result.pvalues, dtype=float),
        exog_names=result.exog_names,
        terms=terms,
        wave=wave,
        spec=spec,
        component="positive_zero_truncated_count",
        component_label=spec.positive_label,
        model_family="Zero-truncated negative binomial",
        response=spec.positive_col,
        n_observations=int(len(use)),
        n_events=None,
    )
    diag = {
        "wave_group": wave,
        "wave_label": WAVE_LABELS.get(wave, wave),
        "outcome": spec.name,
        "component": "positive_zero_truncated_count",
        "skipped": False,
        "reason": "",
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
        "lineage_adjustment": "wave-stratified; lineage dummies included, rank-dropped if collinear",
        "n_windows": n_windows,
        "converged": bool(result.converged),
        "iterations": int(result.nit),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "alpha": float(result.alpha),
        "alpha_at_upper_bound": bool(np.isclose(result.alpha, math.exp(8.0))),
        "optimizer_message": result.message,
        "warnings": "; ".join(str(w.message) for w in caught),
    }
    return rows, diag


def fit_wave_outcome_models(
    clusters: pd.DataFrame,
    *,
    maxiter: int,
    min_clusters: int,
    min_windows: int,
    min_positive: int,
    min_events: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary_specs = [spec for spec in COUNT_MODEL_SPECS if not spec.include_size]
    calendar_cols = [col for col in clusters.columns if col.startswith("calendar_spline_")]
    frames = []
    diagnostics = []

    for wave in WAVE_ORDER:
        wave_df = clusters.loc[clusters["wave_group"] == wave].copy()
        n_windows = int(wave_df["window_id"].nunique())
        if len(wave_df) < min_clusters or n_windows < min_windows:
            for spec in primary_specs:
                diagnostics.append(
                    skipped_diag(
                        wave,
                        spec,
                        "all_components",
                        "below minimum clusters/windows",
                        len(wave_df),
                        n_windows,
                    )
                )
            continue

        lineage_levels_wave = lineage_levels(wave_df)
        for spec in primary_specs:
            print(f"  - {wave} {spec.name}: hurdle binary", flush=True)
            rows, diag = fit_wave_binary_component(
                wave_df,
                wave,
                spec,
                lineage_levels_wave,
                calendar_cols,
                maxiter=maxiter,
                min_events=min_events,
            )
            if not rows.empty:
                frames.append(rows)
            diagnostics.append(diag)

            print(f"  - {wave} {spec.name}: zero-truncated NB positive count", flush=True)
            rows, diag = fit_wave_positive_component(
                wave_df,
                wave,
                spec,
                lineage_levels_wave,
                calendar_cols,
                maxiter=maxiter,
                min_positive=min_positive,
                min_windows=min_windows,
            )
            if not rows.empty:
                frames.append(rows)
            diagnostics.append(diag)
            gc.collect()

    return pd.concat(frames, ignore_index=True), pd.DataFrame(diagnostics)


def summarise_wave_outcomes(clusters: pd.DataFrame) -> pd.DataFrame:
    rows = []
    outcome_specs = [
        ("cluster_size", "Cluster size", "cluster_size_gt1", "cluster_size_excess"),
        ("duration_days", "Duration", "duration_gt0", "duration_positive_days"),
        ("cluster_n_datazones", "Geographic spread", "datazones_gt1", "datazones_excess"),
    ]
    for wave in WAVE_ORDER:
        sub = clusters.loc[clusters["wave_group"] == wave]
        if sub.empty:
            continue
        for raw_col, label, binary_col, positive_col in outcome_specs:
            values = sub[raw_col].dropna()
            pos = sub.loc[sub[positive_col] > 0, positive_col].dropna()
            rows.append(
                {
                    "wave_group": wave,
                    "wave_label": WAVE_LABELS.get(wave, wave),
                    "outcome": raw_col,
                    "outcome_label": label,
                    "n_clusters": int(len(sub)),
                    "n_windows": int(sub["window_id"].nunique()),
                    "structural_minimum_fraction": float(1 - sub[binary_col].mean()),
                    "exceeds_minimum_fraction": float(sub[binary_col].mean()),
                    "mean": float(values.mean()),
                    "median": float(values.median()),
                    "p75": float(values.quantile(0.75)),
                    "p90": float(values.quantile(0.90)),
                    "p95": float(values.quantile(0.95)),
                    "positive_n": int(len(pos)),
                    "positive_mean": float(pos.mean()) if len(pos) else np.nan,
                    "positive_median": float(pos.median()) if len(pos) else np.nan,
                    "positive_p90": float(pos.quantile(0.90)) if len(pos) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def run(
    root: Path,
    *,
    maxiter: int,
    min_clusters: int,
    min_windows: int,
    min_positive: int,
    min_events: int,
    tables_dir: Path | None = None,
    cache_dir: Path | None = None,
) -> None:
    main_dir = root / "part1" / "main"
    if tables_dir is None:
        tables_dir = main_dir / "tables"
    if cache_dir is None:
        cache_dir = main_dir / "cache"
    clusters = pd.read_parquet(cache_dir / "main_cluster_table.parquet")
    clusters["wave_group"] = clusters["pango_lineage"].astype(str).map(assign_wave)

    print("Fitting wave-specific cluster outcome models", flush=True)
    results, diagnostics = fit_wave_outcome_models(
        clusters,
        maxiter=maxiter,
        min_clusters=min_clusters,
        min_windows=min_windows,
        min_positive=min_positive,
        min_events=min_events,
    )
    results.to_csv(tables_dir / "main_wave_specific_hurdle_count_model_results.csv", index=False)
    diagnostics.to_csv(tables_dir / "main_wave_specific_hurdle_count_model_diagnostics.csv", index=False)
    summarise_wave_outcomes(clusters).to_csv(
        tables_dir / "main_wave_cluster_outcome_descriptives.csv",
        index=False,
    )
    print(f"Wrote {tables_dir / 'main_wave_specific_hurdle_count_model_results.csv'}", flush=True)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument("--maxiter", type=int, default=1000)
    parser.add_argument("--min-clusters", type=int, default=1000)
    parser.add_argument("--min-windows", type=int, default=4)
    parser.add_argument("--min-positive", type=int, default=500)
    parser.add_argument("--min-events", type=int, default=50)
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=None,
        help=(
            "Directory for output CSV tables. Defaults to part1/main/tables. "
            "Set to match the --tables-dir used by main_analysis.py for the "
            "same sensitivity run."
        ),
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing the main_cluster_table.parquet cache. "
            "Defaults to part1/main/cache. Must match the --cache-dir used "
            "by main_analysis.py for the same sensitivity run."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    root = args.root.resolve()
    run(
        root=root,
        maxiter=args.maxiter,
        min_clusters=args.min_clusters,
        min_windows=args.min_windows,
        min_positive=args.min_positive,
        min_events=args.min_events,
        tables_dir=args.tables_dir.resolve() if args.tables_dir else None,
        cache_dir=args.cache_dir.resolve() if args.cache_dir else None,
    )


if __name__ == "__main__":
    main()
