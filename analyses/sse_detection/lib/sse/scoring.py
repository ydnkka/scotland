"""SSE node scoring, null calibration, and candidate tiering."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from .config import DETECTION_RANDOM_SEED, MIN_CLUSTER_SIZE, N_PERMUTATIONS
from .transition_graph import add_graph_role_indicators


def safe_divide(numerator, denominator) -> pd.Series:
    """Divide safely, returning NaN where the denominator is zero or missing."""
    denom = pd.Series(denominator).replace(0, np.nan)
    out = pd.Series(numerator) / denom
    return out.replace([np.inf, -np.inf], np.nan)


def log1p_ratio(numerator, denominator) -> pd.Series:
    """Return ``log((1 + numerator) / (1 + denominator))``."""
    return np.log1p(numerator) - np.log1p(denominator)


def mean_with_count(
    df: pd.DataFrame,
    cols: Sequence[str],
) -> tuple[pd.Series, pd.Series]:
    """Row-wise mean over existing columns plus number of valid components."""
    existing = [col for col in cols if col in df.columns]
    if not existing:
        return (
            pd.Series(np.nan, index=df.index, dtype=float),
            pd.Series(0, index=df.index, dtype=int),
        )
    sub = df[existing]
    return sub.mean(axis=1, skipna=True), sub.notna().sum(axis=1)


def empirical_upper_tail_p(
    obs: np.ndarray,
    null_scores: np.ndarray,
    *,
    tie_method: str = "conservative",
    random_state: int | np.random.Generator | None = None,
) -> np.ndarray:
    """Compute smoothed upper-tail empirical p-values from null scores.

    ``conservative`` counts every null tie in the upper tail. ``randomized``
    distributes ties reproducibly across their attainable p-value interval;
    this is useful when discrete composite scores otherwise accumulate at 1.
    """
    obs = np.asarray(obs, dtype=float)
    null_scores = np.asarray(null_scores, dtype=float)

    if null_scores.ndim != 2:
        raise ValueError("null_scores must be a 2D array with shape (n, B).")
    if obs.ndim != 1:
        raise ValueError("obs must be a 1D array with shape (n,).")
    if null_scores.shape[0] != obs.shape[0]:
        raise ValueError("obs and null_scores must have the same number of rows.")
    if tie_method not in {"conservative", "randomized"}:
        raise ValueError("tie_method must be 'conservative' or 'randomized'.")

    valid_null = ~np.isnan(null_scores)
    n_valid = valid_null.sum(axis=1)
    greater = np.sum((null_scores > obs[:, None]) & valid_null, axis=1)
    ties = np.sum((null_scores == obs[:, None]) & valid_null, axis=1)
    if tie_method == "conservative":
        upper_count = greater + ties
    else:
        rng = (
            random_state
            if isinstance(random_state, np.random.Generator)
            else np.random.default_rng(random_state)
        )
        upper_count = greater + rng.random(obs.shape[0]) * ties
    p = (1 + upper_count) / (1 + n_valid)
    p[np.isnan(obs) | (n_valid == 0)] = np.nan
    return p


def choose_permutation_strata(
    df: pd.DataFrame,
    *,
    strata_cols: Sequence[str],
    min_stratum_n: int = 20,
) -> pd.Series:
    """Construct adaptive permutation-stratum labels."""
    if min_stratum_n < 1:
        raise ValueError("min_stratum_n must be at least 1.")

    available = [col for col in strata_cols if col in df.columns]
    if not available:
        return pd.Series("all", index=df.index, dtype="string")

    best_key = df[available[:1]].astype("string").agg("|".join, axis=1)
    for n_cols in range(1, len(available) + 1):
        candidate_key = df[available[:n_cols]].astype("string").agg("|".join, axis=1)
        candidate_counts = candidate_key.map(candidate_key.value_counts())
        best_key = pd.Series(
            np.where(candidate_counts >= min_stratum_n, candidate_key, best_key),
            index=df.index,
            dtype="string",
        )
    return pd.Series(best_key, index=df.index, dtype="string")


def _strata_nested_in_rank_groups(
    strata_codes: np.ndarray,
    rank_codes: np.ndarray,
) -> bool:
    """Return True when every permutation stratum lies inside one rank group."""
    strata_codes = np.asarray(strata_codes)
    rank_codes = np.asarray(rank_codes)

    if strata_codes.shape[0] != rank_codes.shape[0]:
        raise ValueError("strata_codes and rank_codes must have the same length.")

    order = np.argsort(strata_codes, kind="stable")
    s = strata_codes[order]
    r = rank_codes[order]
    is_boundary = np.r_[True, s[1:] != s[:-1]]
    within_run = ~is_boundary
    previous_positions = np.flatnonzero(within_run) - 1
    return bool(np.all(r[within_run] == r[previous_positions]))


def add_composite_null_scores(
    df: pd.DataFrame,
    *,
    score_col: str,
    raw_component_cols: Sequence[str],
    rank_within: str,
    strata_cols: Sequence[str],
    min_stratum_n: int = 20,
    n_permutations: int = 500,
    random_state: int = 42,
    null_mode: str = "profile",
    p_value_mode: str = "conservative",
) -> pd.DataFrame:
    """Add permutation-based null statistics for a composite node score."""
    if n_permutations < 2:
        raise ValueError("n_permutations must be at least 2.")
    if null_mode not in {"profile", "independent"}:
        raise ValueError("null_mode must be either 'profile' or 'independent'.")
    if p_value_mode not in {"conservative", "randomized"}:
        raise ValueError("p_value_mode must be 'conservative' or 'randomized'.")

    required_cols = [score_col, rank_within]
    missing_required = [col for col in required_cols if col not in df.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {missing_required}")

    out = df.copy()
    null_cols = [
        f"{score_col}_null_mean",
        f"{score_col}_null_sd",
        f"{score_col}_null_z",
        f"{score_col}_upper_p_conservative",
        f"{score_col}_upper_p_randomized",
        f"{score_col}_upper_p",
    ]

    components = [col for col in raw_component_cols if col in out.columns]
    if not components:
        for col in null_cols:
            out[col] = np.nan
        return out

    rng = np.random.default_rng(random_state)
    strata = choose_permutation_strata(
        out,
        strata_cols=strata_cols,
        min_stratum_n=min_stratum_n,
    )
    strata_values = strata.to_numpy()
    strata_groups = [
        np.flatnonzero(strata_values == stratum) for stratum in pd.unique(strata_values)
    ]

    rank_codes = out[rank_within].astype("category").cat.codes.to_numpy()
    unique_rank_codes = np.unique(rank_codes)
    raw = out[components].to_numpy(dtype=float)
    n_rows, n_components = raw.shape
    null_scores = np.empty((n_rows, n_permutations), dtype=float)

    def pct_rank_within(values_2d: np.ndarray) -> np.ndarray:
        ranked = np.empty_like(values_2d, dtype=float)
        for group_code in unique_rank_codes:
            mask = rank_codes == group_code
            block = values_2d[mask]
            ranked[mask] = (
                pd.DataFrame(block)
                .rank(pct=True, method="average")
                .to_numpy(dtype=float)
            )
        return ranked

    def permute_matrix_within_strata(values_2d: np.ndarray) -> np.ndarray:
        permuted = values_2d.copy()
        for idx in strata_groups:
            if idx.size <= 1:
                continue
            if null_mode == "profile":
                permuted[idx, :] = values_2d[rng.permutation(idx), :]
            else:
                for j in range(n_components):
                    permuted[idx, j] = rng.permutation(values_2d[idx, j])
        return permuted

    strata_codes = pd.factorize(strata, sort=False)[0]
    nested = _strata_nested_in_rank_groups(strata_codes, rank_codes)

    if nested:
        base_pct = pct_rank_within(raw)
        for b in range(n_permutations):
            permuted_pct = permute_matrix_within_strata(base_pct)
            with np.errstate(invalid="ignore"):
                null_scores[:, b] = np.nanmean(permuted_pct, axis=1)
    else:
        for b in range(n_permutations):
            permuted_raw = permute_matrix_within_strata(raw)
            permuted_pct = pct_rank_within(permuted_raw)
            with np.errstate(invalid="ignore"):
                null_scores[:, b] = np.nanmean(permuted_pct, axis=1)

    obs = out[score_col].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        null_mean = np.nanmean(null_scores, axis=1)
        null_sd = np.nanstd(null_scores, axis=1, ddof=1)

    null_z = np.full(n_rows, np.nan, dtype=float)
    valid_z = ~np.isnan(obs) & ~np.isnan(null_mean) & ~np.isnan(null_sd) & (null_sd > 0)
    null_z[valid_z] = (obs[valid_z] - null_mean[valid_z]) / null_sd[valid_z]

    out[f"{score_col}_null_mean"] = null_mean
    out[f"{score_col}_null_sd"] = null_sd
    out[f"{score_col}_null_z"] = null_z
    conservative_p = empirical_upper_tail_p(obs, null_scores)
    randomized_p = empirical_upper_tail_p(
        obs,
        null_scores,
        tie_method="randomized",
        random_state=random_state + 1_000_003,
    )
    out[f"{score_col}_upper_p_conservative"] = conservative_p
    out[f"{score_col}_upper_p_randomized"] = randomized_p
    out[f"{score_col}_upper_p"] = (
        conservative_p if p_value_mode == "conservative" else randomized_p
    )
    return out


def add_sse_node_metrics(
    df: pd.DataFrame,
    *,
    min_stratum_n: int = 20,
    n_permutations: int = N_PERMUTATIONS,
    random_state: int = DETECTION_RANDOM_SEED,
    min_cluster_size: int = MIN_CLUSTER_SIZE,
    null_mode: str = "profile",
    p_value_mode: str = "randomized",
) -> pd.DataFrame:
    """Add calibrated SSE candidate scores and candidate tiers to cluster nodes.

    Randomized tie handling is the operational default because the composite
    scores are discrete. Conservative p-values are retained alongside it for
    audit and sensitivity analysis.
    """
    if null_mode not in {"profile", "independent"}:
        raise ValueError("null_mode must be either 'profile' or 'independent'.")
    if p_value_mode not in {"conservative", "randomized"}:
        raise ValueError("p_value_mode must be 'conservative' or 'randomized'.")

    required_cols = ["cluster_size", "window_idx"]
    missing_required = [col for col in required_cols if col not in df.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {missing_required}")

    out = df.copy()
    out["sse_tested"] = out["cluster_size"].ge(min_cluster_size)
    out["candidate_tier"] = "size_ineligible"

    basic_edge_default_cols = {
        "in_degree": 0,
        "out_degree": 0,
        "in_strength": 0,
        "out_strength": 0,
        "downstream_cluster_burden": 0,
        "new_downstream_burden": 0,
        "supported_new_downstream_burden": 0,
        "source_attributable_new_downstream_burden": 0,
        "cumulative_unique_new_sequences": 0,
        "upstream_novelty_eligible": False,
        "unique_local_new_sequences_ratio": np.nan,
    }
    for col, default in basic_edge_default_cols.items():
        if col not in out.columns:
            out[col] = default
        else:
            out[col] = out[col].fillna(default)

    tested = out.loc[out["sse_tested"]].copy()
    if tested.empty:
        return out

    tested["log_cluster_size"] = np.log1p(tested["cluster_size"])
    context_cols = [col for col in ["window_idx", "clade"] if col in tested.columns]

    expected_log = pd.Series(
        tested["log_cluster_size"].mean(),
        index=tested.index,
        dtype=float,
    )
    scale_log = pd.Series(
        max(tested["log_cluster_size"].std(ddof=1), 1e-6),
        index=tested.index,
        dtype=float,
    )

    for n_cols in range(1, len(context_cols) + 1):
        cols = context_cols[:n_cols]
        grouped = tested.groupby(cols, dropna=False)["log_cluster_size"]
        group_mean = grouped.transform("mean")
        group_std = grouped.transform("std")
        group_n = grouped.transform("size")
        use_group = (group_n >= min_stratum_n) & group_std.gt(0) & group_std.notna()
        expected_log.loc[use_group] = group_mean.loc[use_group]
        scale_log.loc[use_group] = group_std.loc[use_group]

    scale_log = scale_log.clip(lower=1e-6)
    tested["expected_cluster_size"] = np.maximum(np.expm1(expected_log), 0.5)
    tested["sampling_adjusted_excess_size"] = (
        tested["log_cluster_size"] - expected_log
    ) / scale_log

    tested = add_graph_role_indicators(tested)

    tested["log_new_downstream_burden"] = np.log1p(tested["new_downstream_burden"])
    tested["log_excess_over_upstream"] = log1p_ratio(
        tested["cluster_size"],
        tested["in_strength"],
    )
    tested["log_new_downstream_burden_ratio"] = log1p_ratio(
        tested["new_downstream_burden"],
        tested["cluster_size"],
    )
    tested["log_supported_new_downstream_burden_ratio"] = log1p_ratio(
        tested["supported_new_downstream_burden"],
        tested["cluster_size"],
    )
    tested["log_source_attributable_new_downstream_burden_ratio"] = log1p_ratio(
        tested["source_attributable_new_downstream_burden"],
        tested["cluster_size"],
    )
    tested["log_cumulative_unique_new_sequences_ratio"] = log1p_ratio(
        tested["cumulative_unique_new_sequences"],
        tested["cluster_size"],
    )

    rank_cols = [
        "sampling_adjusted_excess_size",
        "unique_local_new_sequences_ratio",
        "log_source_attributable_new_downstream_burden_ratio",
        "log_cumulative_unique_new_sequences_ratio",
    ]
    for col in rank_cols:
        if col not in tested.columns:
            continue
        tested[f"{col}_pct_window"] = tested.groupby("window_idx", dropna=False)[
            col
        ].rank(pct=True, method="average")

    burst_components = [
        "sampling_adjusted_excess_size",
        "unique_local_new_sequences_ratio",
    ]
    burden_components = [
        "log_source_attributable_new_downstream_burden_ratio",
        "log_cumulative_unique_new_sequences_ratio",
    ]
    axis_strata_cols = [c for c in ["window_idx", "clade"] if c in tested.columns]
    # Keep component availability in every adaptive stratum so parentless
    # one-component profiles are never permuted against two-component profiles.
    burst_strata_cols = ["upstream_novelty_eligible", *axis_strata_cols]

    burst_pct = [f"{c}_pct_window" for c in burst_components]
    tested["burst_score"], tested["burst_score_n"] = mean_with_count(tested, burst_pct)
    tested = add_composite_null_scores(
        tested,
        score_col="burst_score",
        raw_component_cols=burst_components,
        rank_within="window_idx",
        strata_cols=burst_strata_cols,
        min_stratum_n=min_stratum_n,
        n_permutations=n_permutations,
        random_state=random_state,
        null_mode=null_mode,
        p_value_mode=p_value_mode,
    )

    burden_present = [c for c in burden_components if c in tested.columns]
    tested["burden_eligible"] = (
        tested["source_attributable_new_downstream_burden"].fillna(0).gt(0)
        | tested["cumulative_unique_new_sequences"].fillna(0).gt(0)
    )

    burden_idx = tested.index[tested["burden_eligible"]]
    for col in (
        "burden_score",
        "burden_score_n",
        "burden_score_null_mean",
        "burden_score_null_sd",
        "burden_score_null_z",
        "burden_score_upper_p_conservative",
        "burden_score_upper_p_randomized",
        "burden_score_upper_p",
    ):
        tested[col] = np.nan

    if len(burden_idx) >= 2:
        burden_sub = tested.loc[burden_idx].copy()
        for c in burden_present:
            burden_sub[f"{c}_burden_pct"] = burden_sub.groupby(
                "window_idx",
                dropna=False,
            )[c].rank(pct=True, method="average")
        burden_sub_pct = [f"{c}_burden_pct" for c in burden_present]
        burden_sub["burden_score"], burden_sub["burden_score_n"] = mean_with_count(
            burden_sub,
            burden_sub_pct,
        )
        burden_strata_cols = [c for c in ["clade"] if c in burden_sub.columns]
        burden_sub = add_composite_null_scores(
            burden_sub,
            score_col="burden_score",
            raw_component_cols=burden_present,
            rank_within="window_idx",
            strata_cols=burden_strata_cols,
            min_stratum_n=min_stratum_n,
            n_permutations=n_permutations,
            random_state=random_state,
            null_mode=null_mode,
            p_value_mode=p_value_mode,
        )
        for col in (
            "burden_score",
            "burden_score_n",
            "burden_score_null_mean",
            "burden_score_null_sd",
            "burden_score_null_z",
            "burden_score_upper_p_conservative",
            "burden_score_upper_p_randomized",
            "burden_score_upper_p",
        ):
            tested.loc[burden_idx, col] = burden_sub[col]

    axis_z_cols = ["burst_score_null_z", "burden_score_null_z"]
    tested["max_axis_null_z"] = tested[axis_z_cols].max(axis=1)
    tested["candidate_rank"] = tested["max_axis_null_z"].rank(
        ascending=False,
        method="min",
        na_option="bottom",
    )

    enough_sequences = tested["cluster_size"].ge(min_cluster_size)
    burden_eligible = tested["burden_eligible"].fillna(False)
    burst_sig = tested["burst_score_upper_p"].le(0.05) & enough_sequences
    burden_sig = (
        tested["burden_score_upper_p"].le(0.05) & enough_sequences & burden_eligible
    )
    any_sig = burst_sig | burden_sig

    def _fired(row_burst, row_burden):
        fired = []
        if row_burst:
            fired.append("burst")
        if row_burden:
            fired.append("burden")
        return "+".join(fired) if fired else "none"

    tested["axes_fired"] = [_fired(b, u) for b, u in zip(burst_sig, burden_sig)]

    burst_possible = tested["burst_score_upper_p"].le(0.10)
    burden_possible = tested["burden_score_upper_p"].le(0.10) & burden_eligible
    any_possible = (burst_possible | burden_possible) & enough_sequences & ~any_sig

    uncalibrated = (
        (
            (tested["burst_score"].ge(0.90) & tested["burst_score_upper_p"].isna())
            | (tested["burden_score"].ge(0.90) & tested["burden_score_upper_p"].isna())
        )
        & enough_sequences
        & ~any_sig
    )

    both_sig = burst_sig & burden_sig
    tested["candidate_tier"] = np.select(
        [
            both_sig,
            burst_sig,
            burden_sig,
            any_possible,
            uncalibrated,
        ],
        [
            "high_priority_both_axes",
            "high_priority_burst",
            "high_priority_burden",
            "possible_review",
            "high_score_uncalibrated",
        ],
        default="background_or_low_information",
    )

    tested["burden_status"] = np.where(
        burden_eligible,
        "burden_tested",
        "burden_not_applicable",
    )
    tested["burst_vs_burden_contrast"] = np.where(
        burden_eligible,
        tested["burst_score_null_z"] - tested["burden_score_null_z"],
        np.nan,
    )

    for col in tested.columns:
        out.loc[tested.index, col] = tested[col]
    return out
