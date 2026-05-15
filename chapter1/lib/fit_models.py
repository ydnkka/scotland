"""High-level fits for Chapter 1.

Each public function in this module produces a tidy
``(results, diagnostics)`` pair of DataFrames.

Model families
--------------
* :func:`fit_main_effects` — outcomes ~ three excess-mixing predictors +
  adjustments; lineage-adjusted; one row per (outcome, component, term).
* :func:`fit_wave_interactions` — adds excess_mixing × wave interactions;
  lineage adjustment is replaced by wave dummies because the two are
  collinear at the broad-group level.
* :func:`fit_size_spline_sensitivity` — refits the size-adjusted
  geographic-spread ZTNB with a natural spline on ``log_cluster_size``
  in place of the linear z-score.
* :func:`fit_finite_sample_mixing_sensitivity` — replaces raw
  observed-minus-expected excess mixing with finite-sample standardised
  excess mixing.
* :func:`fit_joint_profile_adjusted_sensitivity` — adds the joint
  age × sex × SIMD profile predictor to the main age, sex, and SIMD
  predictor set.
* :func:`fit_null_residual_sensitivity` — refits the main model after
  replacing each observed-minus-expected excess-mixing predictor with the
  residual from a per-mixing-dimension null regression of observed
  discordance on log_size + composition + lineage + window.
* :func:`fit_domain_main_effects` — refits the main model with each SIMD
  domain's quintile mixing in place of overall SIMD mixing (one model per
  domain), keeping age and sex predictors fixed.
* :func:`fit_wave_stratified` — refits the main model within each wave.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from patsy import dmatrix

from .constants import (
    ADJUSTMENT_TERMS,
    CORE_MIXING_PREDICTORS,
    DOMAINS,
    EXCESS_MIXING_TERMS,
    OUTCOMES,
    OutcomeSpec,
    PROFILE_PREDICTORS,
    TERM_LABELS,
    WAVE_ORDER,
    WAVE_REFERENCE,
)
from .estimators import (
    build_exog,
    fit_ztnb,
    lineage_levels,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _calendar_cols(clusters: pd.DataFrame) -> list[str]:
    return [c for c in clusters.columns if c.startswith("calendar_spline_")]


def _label_for(term: str) -> str:
    """Human label for a coefficient name; handles interaction columns."""
    if "__x__" in term:
        left, right = term.split("__x__", 1)
        return f"{TERM_LABELS.get(left, left)} × {right}"
    return TERM_LABELS.get(term, term)


def _extract_rows(
    *,
    outcome: str,
    component: str,
    model_label: str,
    params: np.ndarray,
    bse: np.ndarray,
    pvalues: np.ndarray,
    exog_names: Sequence[str],
    terms_of_interest: Iterable[str],
) -> pd.DataFrame:
    """Pull rows for the named terms (plus any interaction column starting
    with one of them) from a fitted parameter vector."""
    terms_of_interest = list(terms_of_interest)
    rows: list[dict[str, object]] = []
    name_to_idx = {name: i for i, name in enumerate(exog_names)}
    for term in terms_of_interest:
        for name in exog_names:
            if name != term and not name.startswith(term + "__x__"):
                continue
            idx = name_to_idx[name]
            est = float(params[idx])
            se = float(bse[idx])
            p = float(pvalues[idx])
            rows.append({
                "outcome": outcome,
                "component": component,
                "model": model_label,
                "term": name,
                "term_label": _label_for(name),
                "estimate": est,
                "std_error": se,
                "p_value": p,
                "ratio": float(np.exp(est)) if component != "linear" else np.nan,
                "ratio_lower": (
                    float(np.exp(est - 1.96 * se))
                    if component != "linear" else np.nan
                ),
                "ratio_upper": (
                    float(np.exp(est + 1.96 * se))
                    if component != "linear" else np.nan
                ),
                "estimate_lower": float(est - 1.96 * se),
                "estimate_upper": float(est + 1.96 * se),
            })
    return pd.DataFrame(rows)


def _diagnostic_row(
    *,
    outcome: str,
    component: str,
    model_label: str,
    n_obs: int,
    converged: bool,
    llf: float,
    aic: float,
    alpha: float | None = None,
    extra: dict[str, object] | None = None,
) -> dict[str, object]:
    row = {
        "outcome": outcome,
        "component": component,
        "model": model_label,
        "n_obs": n_obs,
        "converged": converged,
        "llf": llf,
        "aic": aic,
        "alpha": alpha,
    }
    if extra:
        row.update(extra)
    return row


def _tail_sensitivity_extra(
    y: np.ndarray,
    *,
    winsorise_quantile: float,
    winsorise_cap: int | None,
    exclude_tail_quantile: float,
    exclude_tail_cap: float | None,
    n_before_tail_filter: int,
    n_after_tail_filter: int,
) -> dict[str, object]:
    """Diagnostic metadata for tail-influence sensitivities."""
    return {
        "n_obs_before_tail_filter": n_before_tail_filter,
        "n_tail_excluded": n_before_tail_filter - n_after_tail_filter,
        "tail_excluded": bool(exclude_tail_quantile > 0.0),
        "tail_exclude_quantile": (
            exclude_tail_quantile if exclude_tail_quantile > 0.0 else None
        ),
        "tail_exclude_cap": exclude_tail_cap,
        "winsorised": bool(winsorise_quantile > 0.0),
        "winsorise_quantile": (
            winsorise_quantile if winsorise_quantile > 0.0 else None
        ),
        "winsorise_cap": winsorise_cap,
        "mean_response": float(np.mean(y)) if len(y) else np.nan,
        "max_response": float(np.max(y)) if len(y) else np.nan,
    }


# ---------------------------------------------------------------------------
# Single-outcome fitter
# ---------------------------------------------------------------------------


def _fit_outcome(
    *,
    spec: OutcomeSpec,
    clusters: pd.DataFrame,
    numeric_terms: Sequence[str],
    calendar_cols: Sequence[str],
    lineage_levels_all: Sequence[str] | None,
    wave_reference: str | None,
    interaction_with_wave: Sequence[str],
    cluster_by: str,
    model_label: str,
    maxiter: int,
    terms_of_interest: Sequence[str],
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    """Fit a ZTNB on the non-singleton sub-population.

    Chapter 1 is ZTNB-only: the population is non-singleton clusters
    """
    full = clusters.dropna(subset=list(numeric_terms)).copy()
    df_pos = full.loc[full["cluster_size"] >= 2].copy()

    if len(df_pos) == 0:
        return pd.DataFrame(), [_diagnostic_row(
            outcome=spec.name, component="positive_count",
            model_label=model_label,
            n_obs=0, converged=False, llf=np.nan, aic=np.nan,
            extra={"error": "no_non_singleton_rows"},
        )]

    n_before_tail_filter = len(df_pos)
    y_pos = df_pos[spec.positive_col].to_numpy(dtype=float)
    exclude_tail_cap: float | None = None
    if exclude_tail_quantile > 0.0:
        exclude_tail_cap = float(np.quantile(y_pos, exclude_tail_quantile))
        keep = y_pos <= exclude_tail_cap
        df_pos = df_pos.loc[keep].copy()
        y_pos = y_pos[keep]

    if len(df_pos) == 0:
        return pd.DataFrame(), [_diagnostic_row(
            outcome=spec.name, component="positive_count",
            model_label=model_label,
            n_obs=0, converged=False, llf=np.nan, aic=np.nan,
            extra={
                "error": "no_rows_after_tail_filter",
                "n_obs_before_tail_filter": n_before_tail_filter,
                "tail_excluded": True,
                "tail_exclude_quantile": exclude_tail_quantile,
                "tail_exclude_cap": exclude_tail_cap,
            },
        )]

    winsorise_cap: int | None = None
    if winsorise_quantile > 0.0:
        winsorise_cap = int(np.quantile(y_pos, winsorise_quantile))
        y_pos = np.minimum(y_pos, winsorise_cap)

    x_pos = build_exog(
        df_pos,
        numeric_terms=numeric_terms,
        calendar_cols=calendar_cols,
        lineage_levels=lineage_levels_all,
        wave_reference=wave_reference,
        interaction_with_wave=interaction_with_wave,
    )
    groups_pos = df_pos[cluster_by].to_numpy()

    try:
        ztnb = fit_ztnb(
            y=y_pos, x=x_pos, groups=groups_pos, maxiter=maxiter,
        )
    except Exception as exc:  # noqa: BLE001
        return pd.DataFrame(), [_diagnostic_row(
            outcome=spec.name, component="positive_count",
            model_label=model_label,
            n_obs=len(df_pos), converged=False, llf=np.nan, aic=np.nan,
            extra={"error": str(exc)},
        )]

    results = _extract_rows(
        outcome=spec.name,
        component="positive_count",
        model_label=model_label,
        params=ztnb.params[:-1],
        bse=ztnb.bse[:-1],
        pvalues=ztnb.pvalues[:-1],
        exog_names=list(x_pos.columns),
        terms_of_interest=terms_of_interest,
    )
    diagnostics = [_diagnostic_row(
        outcome=spec.name, component="positive_count",
        model_label=model_label,
        n_obs=len(y_pos), converged=ztnb.converged,
        llf=ztnb.llf, aic=ztnb.aic, alpha=ztnb.alpha,
        extra=_tail_sensitivity_extra(
            y_pos,
            winsorise_quantile=winsorise_quantile,
            winsorise_cap=winsorise_cap,
            exclude_tail_quantile=exclude_tail_quantile,
            exclude_tail_cap=exclude_tail_cap,
            n_before_tail_filter=n_before_tail_filter,
            n_after_tail_filter=len(df_pos),
        ),
    )]
    return results, diagnostics


# ---------------------------------------------------------------------------
# Public fits
# ---------------------------------------------------------------------------


def fit_main_effects(
    clusters: pd.DataFrame,
    *,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    include_log_size: bool = False,
    extra_mixing_terms: Iterable[str] = (),
    model_label: str = "main",
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Main-effects model: outcomes ~ excess mixing + adjustments + lineage.

    ``include_log_size=True`` is used for the size-adjusted geographic-spread
    model.  ``extra_mixing_terms`` lets the caller drop in alternative
    mixing variables (e.g. profile-only or null-residual).
    """
    calendar_cols = _calendar_cols(clusters)
    lineage_all = lineage_levels(clusters)

    mixing_terms = list(extra_mixing_terms) or list(EXCESS_MIXING_TERMS)
    adjustment = list(ADJUSTMENT_TERMS)
    if include_log_size:
        adjustment = adjustment + ["log_cluster_size_z"]

    numeric_terms = mixing_terms + adjustment
    terms_of_interest = mixing_terms + adjustment

    all_results: list[pd.DataFrame] = []
    all_diag: list[dict[str, object]] = []
    for spec in OUTCOMES:
        # For the size-adjusted variant, only fit geographic_spread.
        if include_log_size and spec.name != "geographic_spread":
            continue
        results, diag = _fit_outcome(
            spec=spec,
            clusters=clusters,
            numeric_terms=numeric_terms,
            calendar_cols=calendar_cols,
            lineage_levels_all=lineage_all,
            wave_reference=None,
            interaction_with_wave=(),
            cluster_by=cluster_by,
            model_label=(
                f"{model_label}_size_adjusted" if include_log_size
                else model_label
            ),
            maxiter=maxiter,
            terms_of_interest=terms_of_interest,
            winsorise_quantile=winsorise_quantile,
            exclude_tail_quantile=exclude_tail_quantile,
        )
        all_results.append(results)
        all_diag.extend(diag)

    return (
        pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame(),
        pd.DataFrame(all_diag),
    )


def fit_wave_interactions(
    clusters: pd.DataFrame,
    *,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    wave_reference: str = WAVE_REFERENCE,
    model_label: str = "wave_interaction",
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Interaction model: mixing × wave terms.

    Lineage dummies are *replaced* by wave dummies because including both
    would be over-parameterised — wave is essentially a coarse pooling of
    lineage.  The wave reference category is dropped from the design.
    """
    calendar_cols = _calendar_cols(clusters)
    mixing_terms = list(EXCESS_MIXING_TERMS)
    numeric_terms = mixing_terms + list(ADJUSTMENT_TERMS)
    terms_of_interest = numeric_terms + ["wave"]

    # Restrict to the named waves with enough non-singleton clusters
    # available; the wave_group "Other" category is dropped.
    in_panel = clusters["wave_group"].astype(str).isin(WAVE_ORDER)
    panel = clusters.loc[in_panel].copy()

    all_results: list[pd.DataFrame] = []
    all_diag: list[dict[str, object]] = []
    for spec in OUTCOMES:
        results, diag = _fit_outcome(
            spec=spec,
            clusters=panel,
            numeric_terms=numeric_terms,
            calendar_cols=calendar_cols,
            lineage_levels_all=None,
            wave_reference=wave_reference,
            interaction_with_wave=mixing_terms,
            cluster_by=cluster_by,
            model_label=model_label,
            maxiter=maxiter,
            terms_of_interest=terms_of_interest,
            winsorise_quantile=winsorise_quantile,
            exclude_tail_quantile=exclude_tail_quantile,
        )
        all_results.append(results)
        all_diag.extend(diag)

    return (
        pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame(),
        pd.DataFrame(all_diag),
    )


# ---------------------------------------------------------------------------
# Sensitivity 1 — non-linear log(size) in size-adjusted spread
# ---------------------------------------------------------------------------


def fit_size_spline_sensitivity(
    clusters: pd.DataFrame,
    *,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    spline_df: int = 4,
    model_label: str = "size_spline",
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refit the size-adjusted spread ZTNB with a natural spline on log_size."""
    calendar_cols = _calendar_cols(clusters)
    lineage_all = lineage_levels(clusters)

    panel = clusters.dropna(
        subset=list(EXCESS_MIXING_TERMS) + list(ADJUSTMENT_TERMS)
    ).copy()
    panel["log_size"] = np.log(panel["cluster_size"].astype(float))
    spline = dmatrix(
        f"bs(log_size, df={spline_df}, degree=3, include_intercept=False) - 1",
        panel,
        return_type="dataframe",
    )
    spline.columns = [f"log_size_spline_{i + 1}" for i in range(spline.shape[1])]
    panel = pd.concat(
        [panel.reset_index(drop=True), spline.reset_index(drop=True)],
        axis=1,
    )

    spline_cols = list(spline.columns)
    numeric_terms = (
        list(EXCESS_MIXING_TERMS) + list(ADJUSTMENT_TERMS) + spline_cols
    )
    terms_of_interest = list(EXCESS_MIXING_TERMS) + spline_cols

    spec = next(o for o in OUTCOMES if o.name == "geographic_spread")
    results, diag = _fit_outcome(
        spec=spec,
        clusters=panel,
        numeric_terms=numeric_terms,
        calendar_cols=calendar_cols,
        lineage_levels_all=lineage_all,
        wave_reference=None,
        interaction_with_wave=(),
        cluster_by=cluster_by,
        model_label=model_label,
        maxiter=maxiter,
        terms_of_interest=terms_of_interest,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )
    return results, pd.DataFrame(diag)


# ---------------------------------------------------------------------------
# Sensitivity 2 — null-regression residual excess mixing
# ---------------------------------------------------------------------------


def build_null_residual_mixing(
    clusters: pd.DataFrame,
    *,
    mixing_prefixes: Iterable[str] = CORE_MIXING_PREDICTORS,
) -> pd.DataFrame:
    """Add ``{prefix}_excess_mixing_null_z`` columns to ``clusters``.

    The null model regresses *observed* discordance on log_size + lineage +
    window-spline + the cluster's own marginal entropy.  The residual is
    standardised and stored alongside the observed-minus-expected predictor.
    Population: clusters with non-missing observed discordance for the
    variable.
    """
    out = clusters.copy()
    calendar_cols = _calendar_cols(out)
    lineage_all = lineage_levels(out)

    for prefix in mixing_prefixes:
        obs_col = f"{prefix}_discordance"
        ent_col = f"{prefix}_entropy"
        if obs_col not in out.columns or ent_col not in out.columns:
            continue
        mask = out[obs_col].notna() & out[ent_col].notna()
        df = out.loc[mask].copy()
        df["log_size"] = np.log(df["cluster_size"].astype(float))

        x_cols = ["log_size", ent_col]
        design = build_exog(
            df,
            numeric_terms=x_cols,
            calendar_cols=calendar_cols,
            lineage_levels=lineage_all,
            wave_reference=None,
            interaction_with_wave=(),
        )
        y = df[obs_col].astype(float).to_numpy()
        try:
            beta, *_ = np.linalg.lstsq(design.values, y, rcond=None)
            fitted = design.values @ beta
            residuals = y - fitted
            null_col = f"{prefix}_excess_mixing_null_raw"
            z_col = f"{prefix}_excess_mixing_null_z"
            out.loc[mask, null_col] = residuals
            mean = float(np.nanmean(residuals))
            sd = float(np.nanstd(residuals, ddof=0))
            if sd > 0:
                out.loc[mask, z_col] = (residuals - mean) / sd
        except Exception:  # noqa: BLE001
            continue
    return out


def fit_null_residual_sensitivity(
    clusters_with_null: pd.DataFrame,
    *,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    model_label: str = "null_residual",
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    null_terms = [f"{p}_excess_mixing_null_z" for p in CORE_MIXING_PREDICTORS]
    return fit_main_effects(
        clusters_with_null,
        cluster_by=cluster_by,
        maxiter=maxiter,
        include_log_size=False,
        extra_mixing_terms=null_terms,
        model_label=model_label,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )


# ---------------------------------------------------------------------------
# Joint-profile predictor
# ---------------------------------------------------------------------------


def fit_simd_decile_sensitivity(
    clusters: pd.DataFrame,
    *,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    model_label: str = "simd_decile",
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refit the main model swapping SIMD-quintile mixing for SIMD-decile.

    Age and sex excess mixing are unchanged.  The SIMD channel is
    ``simd_decile_excess_mixing_z`` (observed minus expected pairwise
    decile-discordance, z-scored) instead of
    ``simd_excess_mixing_z`` (the quintile version).

    The purpose is to check that the SIMD-mixing → cluster-scale finding
    is not an artefact of where the quintile cutpoints fall.  Decile
    resolution is more sensitive to within-quintile gradient but exposes
    the expected-discordance calculation to more small-cell noise within
    window × lineage strata.
    """
    if "simd_decile_excess_mixing_z" not in clusters.columns:
        raise KeyError(
            "simd_decile_excess_mixing_z is not in the cluster table — "
            "make sure the cluster table was built with dz_simd_decile in "
            "SEQUENCE_COLUMNS and simd_decile registered in MIXING_VARIABLES."
        )
    mixing_terms = [
        "age_excess_mixing_z",
        "sex_excess_mixing_z",
        "simd_decile_excess_mixing_z",
    ]
    return fit_main_effects(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        include_log_size=False,
        extra_mixing_terms=mixing_terms,
        model_label=model_label,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )


def fit_finite_sample_mixing_sensitivity(
    clusters: pd.DataFrame,
    *,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    model_label: str = "finite_sample_mixing",
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refit main effects using finite-sample standardised excess mixing.

    The primary observed-minus-expected excess-discordance score is scaled
    only across clusters.  This sensitivity first divides each cluster's
    excess discordance by the approximate binomial pair-sampling standard
    error under the lineage-window null, then z-scores the resulting value.
    """
    finite_terms = [
        f"{prefix}_finite_sample_mixing_z"
        for prefix in CORE_MIXING_PREDICTORS
    ]
    missing = [term for term in finite_terms if term not in clusters.columns]
    if missing:
        raise KeyError(
            "Missing finite-sample mixing columns: " + ", ".join(missing)
        )
    return fit_main_effects(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        include_log_size=False,
        extra_mixing_terms=finite_terms,
        model_label=model_label,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )


def fit_joint_profile_adjusted_sensitivity(
    clusters: pd.DataFrame,
    *,
    profile_name: str = "socio_demographic_profile",
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    model_label: str = "joint_profile_adjusted",
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refit main effects with age, sex, SIMD, and a joint-profile term."""
    if profile_name not in PROFILE_PREDICTORS:
        raise ValueError(
            f"profile_name must be one of {tuple(PROFILE_PREDICTORS)}; "
            f"got {profile_name!r}"
        )
    profile_term = f"{profile_name}_excess_mixing_z"
    return fit_main_effects(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        include_log_size=False,
        extra_mixing_terms=[*EXCESS_MIXING_TERMS, profile_term],
        model_label=model_label,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )


def fit_profile_predictor(
    clusters: pd.DataFrame,
    *,
    profile_name: str = "socio_demographic_profile",
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    model_label: str | None = None,
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refit main model using a single joint-profile excess-mixing predictor.

    ``profile_name`` is one of :data:`constants.PROFILE_PREDICTORS`:

    * ``"demographic_profile"`` — age × sex.
    * ``"socio_demographic_profile"`` — age × sex × SIMD quintile.
    """
    if profile_name not in PROFILE_PREDICTORS:
        raise ValueError(
            f"profile_name must be one of {tuple(PROFILE_PREDICTORS)}; "
            f"got {profile_name!r}"
        )
    profile_term = f"{profile_name}_excess_mixing_z"
    return fit_main_effects(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        include_log_size=False,
        extra_mixing_terms=[profile_term],
        model_label=model_label or f"profile_{profile_name}",
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )


# ---------------------------------------------------------------------------
# Domain stratification
# ---------------------------------------------------------------------------


def fit_domain_main_effects(
    clusters: pd.DataFrame,
    *,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refit the main model swapping overall SIMD mixing for each domain.

    Age and sex excess mixing remain in every model so the SIMD-channel
    estimate is comparable across domains.
    """
    all_results: list[pd.DataFrame] = []
    all_diag: list[pd.DataFrame] = []
    for domain in DOMAINS:
        domain_term = f"{domain}_domain_excess_mixing_z"
        if domain_term not in clusters.columns:
            continue
        mixing_terms = [
            "age_excess_mixing_z",
            "sex_excess_mixing_z",
            domain_term,
        ]
        results, diag = fit_main_effects(
            clusters,
            cluster_by=cluster_by,
            maxiter=maxiter,
            include_log_size=False,
            extra_mixing_terms=mixing_terms,
            model_label=f"domain_{domain}",
            winsorise_quantile=winsorise_quantile,
            exclude_tail_quantile=exclude_tail_quantile,
        )
        if not results.empty:
            results["domain"] = domain
        diag["domain"] = domain
        all_results.append(results)
        all_diag.append(diag)
    return (
        pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame(),
        pd.concat(all_diag, ignore_index=True) if all_diag else pd.DataFrame(),
    )


# ---------------------------------------------------------------------------
# Wave stratification
# ---------------------------------------------------------------------------


def fit_wave_stratified(
    clusters: pd.DataFrame,
    *,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    min_clusters_per_wave: int = 50,
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refit the main model separately within each wave.

    Within-wave fits use a wave-restricted lineage list, drop the calendar
    spline columns that become collinear, and rely on the QR pruning to
    handle any rank-deficient lineage columns.
    """
    calendar_cols = _calendar_cols(clusters)
    mixing_terms = list(EXCESS_MIXING_TERMS)
    numeric_terms = mixing_terms + list(ADJUSTMENT_TERMS)

    all_results: list[pd.DataFrame] = []
    all_diag: list[dict[str, object]] = []
    for wave in WAVE_ORDER:
        panel = clusters.loc[clusters["wave_group"].astype(str) == wave].copy()
        n_non_singleton = int((panel["cluster_size"] >= 2).sum())
        if n_non_singleton < min_clusters_per_wave:
            all_diag.append({
                "outcome": "ALL",
                "component": "ALL",
                "model": f"wave_{wave}",
                "n_obs": n_non_singleton,
                "converged": False,
                "llf": np.nan,
                "aic": np.nan,
                "alpha": None,
                "wave": wave,
                "error": "below_min_threshold",
            })
            continue

        wave_lineages = lineage_levels(panel)

        for spec in OUTCOMES:
            results, diag = _fit_outcome(
                spec=spec,
                clusters=panel,
                numeric_terms=numeric_terms,
                calendar_cols=calendar_cols,
                lineage_levels_all=wave_lineages,
                wave_reference=None,
                interaction_with_wave=(),
                cluster_by=cluster_by,
                model_label=f"wave_{wave}",
                maxiter=maxiter,
                terms_of_interest=numeric_terms,
                winsorise_quantile=winsorise_quantile,
                exclude_tail_quantile=exclude_tail_quantile,
            )
            if not results.empty:
                results["wave"] = wave
            for d in diag:
                d["wave"] = wave
            all_results.append(results)
            all_diag.extend(diag)

    return (
        pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame(),
        pd.DataFrame(all_diag),
    )
