"""
Entropy utilities for cluster-category and edge-weight analyses.

Public API
----------
shannon_entropy, shannon_entropy_grouped, max_entropy
    Core entropy primitives.
empirical_tail_p_values
    Empirical lower/upper/two-sided tail probabilities.
cluster_socio_demo_entropy
    Cluster category-entropy vs. within-window multinomial null.
cluster_age_conditional_binary_entropy
    Cluster binary-label entropy vs. window-by-age conditional null.
onward_edge_entropy
    Per-source entropy of outgoing edge weights.
observed_mixing_entropy_scales
    Cluster-level scaled observed-entropy columns for mixing models.
vaccination_mixing_features
    Cluster-level age-conditional vaccination-mixing entropy features.
attach_cluster_stats
    Thin left-merge helper to broadcast cluster stats onto row-level data.

All entropies use log base 2 by default (units: bits).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


DEFAULT_MIXING_FEATURES = [
    "sex_entropy_z",
    "age_entropy_z",
    "simd_entropy_z",
    "datazone_entropy_z",
    "local_authority_entropy_z",
    "urban_rural_entropy_z",
    "health_board_entropy_z",
    "vaccination_entropy_z",
]

OBSERVED_MIXING_FEATURES = [
    "sex_entropy_obs",
    "age_entropy_obs",
    "simd_entropy_obs",
    "datazone_entropy_obs",
    "local_authority_entropy_obs",
    "urban_rural_entropy_obs",
    "health_board_entropy_obs",
    "vaccination_entropy_obs",
]
OBSERVED_MIXING_FEATURES_X10 = [f"{feature}_x10" for feature in OBSERVED_MIXING_FEATURES]

MIXING_TERTILE_FEATURES = DEFAULT_MIXING_FEATURES + OBSERVED_MIXING_FEATURES

MIXING_TERTILE_ORDER = [
    "more_homogeneous",
    "as_expected",
    "more_mixed",
]
VACCINATION_MIXING_TERTILE_ORDER = MIXING_TERTILE_ORDER


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _check_base(base: float) -> None:
    if base <= 0 or base == 1:
        raise ValueError("base must be positive and not equal to 1.")


def _max_entropy(k, base: float = 2):
    """Maximum Shannon entropy for ``k`` equally likely categories: 
    ``log_base(k)``.

    Vectorises over array ``k``.
    """
    return np.log(k) / np.log(base)


def _shannon_entropy(
    counts: np.ndarray,
    base: float = 2,
    axis: int = -1,
) -> np.ndarray:
    """Shannon entropy of one or more discrete count distributions.

    Parameters
    ----------
    counts : np.ndarray
        Non-negative counts. Entropy is computed along ``axis``.
    base : float, default 2
        Logarithm base. 2 -> bits, ``np.e`` -> nats.
    axis : int, default -1
        Axis along which to reduce.

    Returns
    -------
    np.ndarray
        Entropy per slice. Slices summing to zero return 0. Follows the
        ``0 * log(0) = 0`` convention.
    """
    _check_base(base)

    counts = np.asarray(counts, dtype=float)
    totals = counts.sum(axis=axis, keepdims=True)

    probs = np.divide(
        counts,
        totals,
        out=np.zeros_like(counts),
        where=totals > 0,
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        log_probs = np.where(probs > 0, np.log(probs), 0.0)

    return -(probs * log_probs).sum(axis=axis) / np.log(base)


def _shannon_entropy_grouped(
    values: np.ndarray,
    group_codes: np.ndarray,
    base: float = 2,
) -> dict:
    """Per-group Shannon entropy from ragged positive-valued data.

    Uses the identity ``H = log(S) - sum(w log w) / S`` (with ``S = sum w``)
    to avoid building a dense ``(n_groups, max_group_size)`` matrix.

    Parameters
    ----------
    values : np.ndarray
        Positive weights. Non-positive values should be filtered upstream.
    group_codes : np.ndarray
        Integer group labels of the same shape as ``values``.
    base : float, default 2
        Logarithm base.

    Returns
    -------
    dict
        ``unique_codes`` (ascending), ``size``, ``total``, ``entropy``
        (NaN where size <= 1), ``max_value``.
    """
    _check_base(base)
    if values.shape != group_codes.shape:
        raise ValueError(
            f"values and group_codes shape mismatch: "
            f"{values.shape} vs {group_codes.shape}"
        )
    if values.ndim != 1:
        raise ValueError("values and group_codes must be 1-D.")

    order = np.argsort(group_codes, kind="mergesort")
    g_sorted = group_codes[order]
    w_sorted = values[order].astype(float, copy=False)

    starts = np.r_[0, np.flatnonzero(np.diff(g_sorted)) + 1]
    sizes = np.diff(np.r_[starts, len(w_sorted)])
    unique_codes = g_sorted[starts]

    total = np.add.reduceat(w_sorted, starts)
    wlogw = np.add.reduceat(w_sorted * np.log(w_sorted), starts)
    max_value = np.maximum.reduceat(w_sorted, starts)

    with np.errstate(divide="ignore", invalid="ignore"):
        entropy = (np.log(total) - wlogw / total) / np.log(base)
    entropy[sizes <= 1] = np.nan

    return {
        "unique_codes": unique_codes,
        "size": sizes,
        "total": total,
        "entropy": entropy,
        "max_value": max_value,
    }


def _empirical_tail_p_values(
    observed: np.ndarray | float,
    null_values: np.ndarray,
    *,
    axis: int = -1,
    sorted_1d: bool = False,
) -> dict[str, np.ndarray]:
    """Empirical lower-, upper-, and two-sided tail probabilities.

    Uses +1 smoothing, so p-values are never exactly zero:
    upper ``(1 + #{null >= obs}) / (n_null + 1)``,
    lower ``(1 + #{null <= obs}) / (n_null + 1)``.

    Parameters
    ----------
    observed : array-like or float
        Observed statistic(s). Must broadcast against ``null_values`` with
        ``axis`` removed (broadcasting path), or be a 1-D array of
        observations (``sorted_1d`` path).
    null_values : np.ndarray
        Null draws.
    axis : int, default -1
        Axis containing null draws (ignored when ``sorted_1d=True``).
    sorted_1d : bool, default False
        If True, ``null_values`` is a single 1-D null sample shared by every
        observation; uses a sorted-search O(n log n) path instead of the
        O(n_obs * n_null) broadcast. Equivalent results.

    Returns
    -------
    dict
        ``lower_p``, ``upper_p``, ``two_sided_p`` arrays.
    """
    null_values = np.asarray(null_values, dtype=float)
    if null_values.size == 0:
        raise ValueError("null_values must contain at least one draw.")
    observed = np.asarray(observed, dtype=float)

    if sorted_1d:
        if null_values.ndim != 1:
            raise ValueError("sorted_1d requires a 1-D null_values array.")
        nulls = np.sort(null_values)
        n_null = nulls.size
        lower_count = np.searchsorted(nulls, observed, side="right")
        upper_count = n_null - np.searchsorted(nulls, observed, side="left")
        lower_p = (1 + lower_count) / (n_null + 1)
        upper_p = (1 + upper_count) / (n_null + 1)
    else:
        nulls = np.moveaxis(null_values, axis, -1)
        n_null = nulls.shape[-1]
        observed_expanded = np.expand_dims(observed, axis=-1)
        upper_p = (1 + np.sum(nulls >= observed_expanded, axis=-1)) / (n_null + 1)
        lower_p = (1 + np.sum(nulls <= observed_expanded, axis=-1)) / (n_null + 1)

    two_sided_p = np.minimum(1.0, 2 * np.minimum(lower_p, upper_p))
    return {"lower_p": lower_p, "upper_p": upper_p, "two_sided_p": two_sided_p}


def _zscore(observed: np.ndarray, null_mean: np.ndarray, null_sd: np.ndarray) -> np.ndarray:
    """Standardise observed against null moments; NaN where SD is 0/NaN."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(
            (null_sd > 0) & ~np.isnan(null_sd),
            (np.asarray(observed, dtype=float) - null_mean) / null_sd,
            np.nan,
        )


# --------------------------------------------------------------------------- #
# Cluster category entropy (sequence-level input -> cluster-level output)
# --------------------------------------------------------------------------- #
def cluster_socio_demo_entropy(
    df: pd.DataFrame,
    cluster_col: str,
    category_col: str,
    window_col: str,
    *,
    n_random: int = 1000,
    normalise: bool = True,
    base: float = 2,
    random_state: int = 42,
    prefix: str = "cluster",
) -> pd.DataFrame:
    """Cluster category-entropy statistics with a within-window null model.

    Accepts sequence-level data (one row per sequence). The null model draws
    cluster-sized multinomial samples from the window-pooled category
    frequencies; null draws are shared across clusters of the same
    ``(window, size)``.

    Singletons retain ``entropy = 0``, ``null_mean = 0``, ``null_sd = 0``,
    ``z = NaN`` (degenerate null distribution).

    Returns
    -------
    pd.DataFrame
        One row per cluster: ``[cluster_col, window_col]`` plus
        ``{prefix}_n``, ``{prefix}_entropy_obs``,
        ``{prefix}_entropy_null_mean``, ``{prefix}_entropy_null_sd``,
        ``{prefix}_entropy_z``, and the empirical tail probabilities
        ``{prefix}_entropy_lower_p``, ``{prefix}_entropy_upper_p``,
        ``{prefix}_entropy_two_sided_p``.

    Notes
    -----
    ``{prefix}_n`` counts non-missing category observations per cluster.
    The window column is retained because the null is window-stratified and a
    cluster could in principle span windows.
    """
    required_cols = [cluster_col, category_col, window_col]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in dataframe: {missing}")
    if n_random < 2:
        raise ValueError("n_random must be >= 2 to estimate null SD.")
    _check_base(base)

    rng = np.random.default_rng(random_state)

    codes, categories = pd.factorize(df[category_col], sort=False)
    k_global = len(categories)
    if k_global <= 1:
        raise ValueError("Cannot compute entropy for a single category level.")
    norm_factor = _max_entropy(k_global, base=base) if normalise else 1.0

    stat_cols = [
        f"{prefix}_n",
        f"{prefix}_entropy_obs",
        f"{prefix}_entropy_null_mean",
        f"{prefix}_entropy_null_sd",
        f"{prefix}_entropy_z",
        f"{prefix}_entropy_lower_p",
        f"{prefix}_entropy_upper_p",
        f"{prefix}_entropy_two_sided_p",
    ]

    valid = codes >= 0
    work = pd.DataFrame(
        {
            window_col: df.loc[valid, window_col].to_numpy(),
            cluster_col: df.loc[valid, cluster_col].to_numpy(),
            "_cat_code": codes[valid],
        }
    )
    if work.empty:
        return pd.DataFrame(columns=[cluster_col, window_col, *stat_cols])

    cluster_counts = (
        work.groupby([window_col, cluster_col, "_cat_code"], dropna=False, observed=True)
        .size()
        .unstack("_cat_code", fill_value=0)
        .reindex(columns=np.arange(k_global), fill_value=0)
    )

    count_matrix = cluster_counts.to_numpy()
    cluster_n = count_matrix.sum(axis=1).astype(int)
    observed_entropy = _shannon_entropy(count_matrix, base=base) / norm_factor

    n_rows = len(cluster_n)
    null_mean = np.full(n_rows, np.nan)
    null_sd = np.full(n_rows, np.nan)
    lower_p = np.full(n_rows, np.nan)
    upper_p = np.full(n_rows, np.nan)
    two_sided_p = np.full(n_rows, np.nan)

    index_to_row = {idx: i for i, idx in enumerate(cluster_counts.index)}

    for _, window_block in cluster_counts.groupby(level=0, dropna=False, sort=False):
        window_cat_counts = window_block.sum(axis=0).to_numpy(dtype=float)
        window_total = window_cat_counts.sum()
        if window_total == 0:
            continue
        window_probs = window_cat_counts / window_total

        window_rows = np.fromiter(
            (index_to_row[idx] for idx in window_block.index),
            dtype=int,
            count=len(window_block),
        )
        window_sizes = count_matrix[window_rows].sum(axis=1).astype(int)

        for size in np.unique(window_sizes):
            if size <= 0:
                continue
            target_rows = window_rows[window_sizes == size]

            if size == 1:
                null_mean[target_rows] = 0.0
                null_sd[target_rows] = 0.0
                lower_p[target_rows] = 1.0
                upper_p[target_rows] = 1.0
                two_sided_p[target_rows] = 1.0
                continue

            sampled = rng.multinomial(n=size, pvals=window_probs, size=n_random)
            null_H = _shannon_entropy(sampled, base=base) / norm_factor
            null_mean[target_rows] = null_H.mean()
            null_sd[target_rows] = null_H.std(ddof=1)
            tails = _empirical_tail_p_values(
                observed_entropy[target_rows], null_H, sorted_1d=True
            )
            lower_p[target_rows] = tails["lower_p"]
            upper_p[target_rows] = tails["upper_p"]
            two_sided_p[target_rows] = tails["two_sided_p"]

    stats_dict = {
        f"{prefix}_n": cluster_n,
        f"{prefix}_entropy_obs": observed_entropy,
        f"{prefix}_entropy_null_mean": null_mean,
        f"{prefix}_entropy_null_sd": null_sd,
        f"{prefix}_entropy_z": _zscore(observed_entropy, null_mean, null_sd),
        f"{prefix}_entropy_lower_p": lower_p,
        f"{prefix}_entropy_upper_p": upper_p,
        f"{prefix}_entropy_two_sided_p": two_sided_p,
    }

    return (
        pd.DataFrame(stats_dict, index=cluster_counts.index)
        .reset_index()
        .loc[:, [cluster_col, window_col, *stat_cols]]
    )


def cluster_age_conditional_binary_entropy(
    df: pd.DataFrame,
    cluster_col: str,
    binary_col: str,
    window_col: str,
    age_col: str,
    *,
    n_random: int = 1000,
    base: float = 2,
    random_state: int = 42,
    prefix: str = "cluster",
) -> pd.DataFrame:
    """Cluster entropy against a window-by-age conditional null.

    Accepts sequence-level data. The null preserves each cluster's observed
    ``window_col`` x ``age_col`` composition and redraws binary-positive
    counts from the pooled stratum positive probability.

    Returns
    -------
    pd.DataFrame
        One row per cluster: ``[cluster_col]`` plus ``{prefix}_n``,
        ``{prefix}_prop_positive``, ``{prefix}_entropy_obs``,
        ``{prefix}_entropy_null_mean``, ``{prefix}_entropy_null_sd``,
        ``{prefix}_entropy_z``, and the empirical tail probabilities
        ``{prefix}_entropy_lower_p``, ``{prefix}_entropy_upper_p``,
        ``{prefix}_entropy_two_sided_p``.
    """
    required_cols = [cluster_col, binary_col, window_col, age_col]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in dataframe: {missing}")
    if n_random < 2:
        raise ValueError("n_random must be >= 2 to estimate null SD.")
    _check_base(base)

    stat_cols = [
        f"{prefix}_n",
        f"{prefix}_prop_positive",
        f"{prefix}_entropy_obs",
        f"{prefix}_entropy_null_mean",
        f"{prefix}_entropy_null_sd",
        f"{prefix}_entropy_z",
        f"{prefix}_entropy_lower_p",
        f"{prefix}_entropy_upper_p",
        f"{prefix}_entropy_two_sided_p",
    ]

    binary = pd.to_numeric(df[binary_col], errors="coerce")
    valid = df[cluster_col].notna() & binary.notna()
    work = df.loc[valid, [cluster_col, window_col, age_col]].copy()
    work["_positive"] = binary.loc[valid].gt(0).astype(int).to_numpy()

    if work.empty:
        return pd.DataFrame(columns=[cluster_col, *stat_cols])

    stratum_prob = (
        work.groupby([window_col, age_col], dropna=False, observed=True)["_positive"]
        .mean()
        .rename("_stratum_positive_prob")
        .reset_index()
    )
    work = work.merge(stratum_prob, on=[window_col, age_col], how="left")

    cluster_totals = (
        work.groupby(cluster_col, dropna=False, observed=True)
        .agg(
            **{
                f"{prefix}_n": ("_positive", "size"),
                "_observed_positive": ("_positive", "sum"),
            }
        )
        .reset_index()
    )
    counts = np.column_stack(
        [
            cluster_totals["_observed_positive"].to_numpy(dtype=float),
            (
                cluster_totals[f"{prefix}_n"] - cluster_totals["_observed_positive"]
            ).to_numpy(dtype=float),
        ]
    )
    cluster_totals[f"{prefix}_prop_positive"] = (
        cluster_totals["_observed_positive"] / cluster_totals[f"{prefix}_n"]
    )
    cluster_totals[f"{prefix}_entropy_obs"] = _shannon_entropy(counts, base=base)

    rng = np.random.default_rng(random_state)
    n_rows = len(cluster_totals)
    null_mean = np.full(n_rows, np.nan)
    null_sd = np.full(n_rows, np.nan)
    lower_p = np.full(n_rows, np.nan)
    upper_p = np.full(n_rows, np.nan)
    two_sided_p = np.full(n_rows, np.nan)
    cluster_row = {c: i for i, c in enumerate(cluster_totals[cluster_col].tolist())}

    cluster_strata = (
        work.groupby([cluster_col, window_col, age_col], dropna=False, observed=True)
        .agg(n=("_positive", "size"), p=("_stratum_positive_prob", "first"))
        .reset_index()
    )

    for cluster, block in cluster_strata.groupby(cluster_col, dropna=False, sort=False):
        row_idx = cluster_row[cluster]
        n_by_stratum = block["n"].to_numpy(dtype=int)
        p_by_stratum = block["p"].to_numpy(dtype=float)
        total = int(n_by_stratum.sum())
        if total <= 0:
            continue

        sampled_positive = rng.binomial(
            n=n_by_stratum, p=p_by_stratum, size=(n_random, len(n_by_stratum))
        ).sum(axis=1)
        sampled_counts = np.column_stack([sampled_positive, total - sampled_positive])
        null_H = _shannon_entropy(sampled_counts, base=base)
        null_mean[row_idx] = null_H.mean()
        null_sd[row_idx] = null_H.std(ddof=1)
        observed = np.asarray(cluster_totals.loc[row_idx, f"{prefix}_entropy_obs"])
        tails = _empirical_tail_p_values(observed, null_H)
        lower_p[row_idx] = tails["lower_p"]
        upper_p[row_idx] = tails["upper_p"]
        two_sided_p[row_idx] = tails["two_sided_p"]

    cluster_totals[f"{prefix}_entropy_null_mean"] = null_mean
    cluster_totals[f"{prefix}_entropy_null_sd"] = null_sd
    cluster_totals[f"{prefix}_entropy_z"] = _zscore(
        cluster_totals[f"{prefix}_entropy_obs"].to_numpy(dtype=float),
        null_mean,
        null_sd,
    )
    cluster_totals[f"{prefix}_entropy_lower_p"] = lower_p
    cluster_totals[f"{prefix}_entropy_upper_p"] = upper_p
    cluster_totals[f"{prefix}_entropy_two_sided_p"] = two_sided_p

    return cluster_totals.loc[:, [cluster_col, *stat_cols]]


def vaccination_mixing_features(
    df: pd.DataFrame,
    *,
    cluster_col: str = "cluster_id",
    vaccination_col: str = "is_vaccinated",
    window_col: str = "window_idx",
    age_col: str = "age_band",
    n_random: int = 1000,
    random_state: int = 42,
) -> pd.DataFrame:
    """Cluster-level age-conditional vaccination-mixing entropy features.

    Accepts sequence-level data and returns one row per cluster: ``[cluster_col]``
    plus the vaccination entropy stats (``vaccination_entropy_obs``,
    ``vaccination_entropy_z``, etc.) and an ordered ``vaccination_entropy_z_tertile``.
    """
    required = {cluster_col, vaccination_col, window_col, age_col}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing vaccination mixing columns: {missing}")

    seq = df.copy()
    seq["_vaccination_positive"] = (
        pd.to_numeric(seq[vaccination_col], errors="coerce").fillna(0).gt(0).astype(int)
    )
    seq[age_col] = seq[age_col].astype("string").fillna("Missing")

    stats = cluster_age_conditional_binary_entropy(
        seq,
        cluster_col=cluster_col,
        binary_col="_vaccination_positive",
        window_col=window_col,
        age_col=age_col,
        n_random=n_random,
        random_state=random_state,
        prefix="vaccination",
    )
    feature_cols = [
        "vaccination_n",
        "vaccination_prop_positive",
        "vaccination_entropy_obs",
        "vaccination_entropy_null_mean",
        "vaccination_entropy_null_sd",
        "vaccination_entropy_z",
    ]
    out = stats.loc[:, [cluster_col, *feature_cols]].drop_duplicates(cluster_col)
    out["vaccination_entropy_z_tertile"] = _ordered_tertile(
        out["vaccination_entropy_z"], categories=MIXING_TERTILE_ORDER
    )
    return out


# --------------------------------------------------------------------------- #
# Onward edge entropy (edge-level input -> source-level output)
# --------------------------------------------------------------------------- #
def onward_edge_entropy(
    edge_df: pd.DataFrame,
    source_col: str,
    weight_col: str,
    *,
    base: float = 2,
    cluster_col: str = "cluster_id",
) -> pd.DataFrame:
    """Per-source entropy of outgoing edge weights, from edge-level data.

    Measures how evenly each source node distributes weight across its
    successors. Sources with a single successor have entropy NaN.

    Returns
    -------
    pd.DataFrame
        One row per source, keyed by ``cluster_col`` (the source id), with
        ``out_degree``, ``out_strength``, ``onward_entropy``,
        ``onward_entropy_norm``, ``effective_successors``,
        ``dominant_successor_frac``.

    Notes
    -----
    The source identifier is returned in the column named by ``cluster_col``
    to match the shared ``cluster_id + [stats]`` contract; set it to
    ``source_col`` if you prefer to keep the original name.
    ``onward_entropy_norm`` divides by ``log_base(out_degree)`` (in [0, 1]
    for >= 2 successors). ``effective_successors`` is the Hill number of
    order 1, ``base ** entropy``.
    """
    required_cols = [source_col, weight_col]
    missing = [c for c in required_cols if c not in edge_df.columns]
    if missing:
        raise ValueError(f"Missing columns in edge dataframe: {missing}")
    _check_base(base)

    stat_cols = [
        "out_degree",
        "out_strength",
        "onward_entropy",
        "onward_entropy_norm",
        "effective_successors",
        "dominant_successor_frac",
    ]

    df = edge_df[[source_col, weight_col]].copy()
    df[weight_col] = pd.to_numeric(df[weight_col], errors="coerce")
    df = df.dropna()
    df = df[df[weight_col] > 0]
    if df.empty:
        return pd.DataFrame(columns=[cluster_col, *stat_cols])

    codes, nodes = pd.factorize(df[source_col].to_numpy(), sort=False)
    agg = _shannon_entropy_grouped(df[weight_col].to_numpy(dtype=float), codes, base=base)

    sizes = agg["size"]
    H = agg["entropy"]

    H_norm = np.full(len(sizes), np.nan)
    multi = sizes > 1
    H_norm[multi] = H[multi] / _max_entropy(sizes[multi], base=base)
    H_norm = np.clip(H_norm, 0.0, 1.0)  # tame FP overshoot

    effective = np.power(float(base), H)
    dominant = agg["max_value"] / agg["total"]
    node_labels = nodes[agg["unique_codes"]]

    return pd.DataFrame(
        {
            cluster_col: node_labels,
            "out_degree": sizes,
            "out_strength": agg["total"],
            "onward_entropy": H,
            "onward_entropy_norm": H_norm,
            "effective_successors": effective,
            "dominant_successor_frac": dominant,
        }
    )


# --------------------------------------------------------------------------- #
# Mixing-model feature helpers (cluster-level out)
# --------------------------------------------------------------------------- #
def observed_mixing_entropy_scales(
    cluster_stats: pd.DataFrame,
    *,
    cluster_col: str = "cluster_id",
    features: list[str] | None = None,
    scaled_features: list[str] | None = None,
    multiplier: float = 10,
) -> pd.DataFrame:
    """Cluster-level scaled observed-entropy columns for mixing models.

    Accepts a cluster-level stats frame and returns ``[cluster_col]`` plus the
    scaled ``*_x10`` columns (per-0.1 scale for odds-ratio reporting). Only
    features present in the input are produced.
    """
    features = list(OBSERVED_MIXING_FEATURES if features is None else features)
    if scaled_features is None:
        scaled_features = [f"{feature}_x10" for feature in features]
    else:
        scaled_features = list(scaled_features)
    if len(features) != len(scaled_features):
        raise ValueError("features and scaled_features must have the same length.")
    if cluster_col not in cluster_stats.columns:
        raise KeyError(f"{cluster_col!r} is required in cluster_stats.")

    out = cluster_stats[[cluster_col]].copy()
    for feature, scaled_feature in zip(features, scaled_features):
        if feature in cluster_stats.columns:
            out[scaled_feature] = cluster_stats[feature].astype(float) * multiplier
    return out


def add_observed_mixing_entropy_scales(
    cluster_stats: pd.DataFrame,
    *,
    cluster_col: str = "cluster_id",
    features: list[str] | None = None,
    scaled_features: list[str] | None = None,
    multiplier: float = 10,
) -> pd.DataFrame:
    """Return ``cluster_stats`` with observed-entropy scale columns attached."""
    scaled = observed_mixing_entropy_scales(
        cluster_stats,
        cluster_col=cluster_col,
        features=features,
        scaled_features=scaled_features,
        multiplier=multiplier,
    )
    out = cluster_stats.copy()
    new_cols = [col for col in scaled.columns if col != cluster_col]
    for col in new_cols:
        out[col] = scaled[col]
    return out


def add_vaccination_mixing_features(
    node_data: pd.DataFrame,
    sequence_data: pd.DataFrame,
    *,
    cluster_col: str = "cluster_id",
    vaccination_col: str = "is_vaccinated",
    window_col: str = "window_idx",
    age_col: str = "age_band",
    n_random: int = 1000,
    random_state: int = 42,
) -> pd.DataFrame:
    """Attach age-conditional vaccination-mixing features to node rows."""
    required_node = {cluster_col}
    missing_node = sorted(required_node - set(node_data.columns))
    if missing_node:
        raise KeyError(f"Missing vaccination node columns: {missing_node}")

    required_seq = {cluster_col, vaccination_col, window_col, age_col}
    missing_seq = sorted(required_seq - set(sequence_data.columns))
    if missing_seq:
        raise KeyError(f"Missing vaccination mixing columns: {missing_seq}")

    seq = sequence_data.copy()
    seq["_vaccination_mix_positive"] = (
        pd.to_numeric(seq[vaccination_col], errors="coerce")
        .fillna(0)
        .gt(0)
        .astype(int)
    )
    seq[age_col] = seq[age_col].astype("string").fillna("Missing")

    stats = cluster_age_conditional_binary_entropy(
        seq,
        cluster_col=cluster_col,
        binary_col="_vaccination_mix_positive",
        window_col=window_col,
        age_col=age_col,
        n_random=n_random,
        random_state=random_state,
        prefix="vaccination_mix",
    )
    feature_cols = [
        "vaccination_mix_n",
        "vaccination_mix_prop_positive",
        "vaccination_mix_entropy_obs",
        "vaccination_mix_entropy_null_mean",
        "vaccination_mix_entropy_null_sd",
        "vaccination_mix_entropy_z",
    ]
    stats = stats.loc[:, [cluster_col, *feature_cols]].drop_duplicates(cluster_col)
    stats["vaccination_mix_tertile"] = _ordered_tertile(
        stats["vaccination_mix_entropy_z"],
        categories=VACCINATION_MIXING_TERTILE_ORDER,
    )
    return node_data.merge(stats, on=cluster_col, how="left")


def _ordered_tertile(values: pd.Series, *, categories: list[str]) -> pd.Categorical:
    """Rank-based ordered tertiles for a numeric series."""
    out = pd.Series(pd.NA, index=values.index, dtype="object")
    present = values.notna()
    if present.sum() >= 3:
        ranks = values.loc[present].rank(method="first")
        bins = pd.qcut(ranks, q=3, labels=categories)
        out.loc[present] = bins.astype(str).to_numpy()
    return pd.Categorical(out, categories=categories, ordered=True)


def add_mixing_tertiles(
    cluster_stats: pd.DataFrame,
    *,
    features: list[str] | None = None,
    categories: list[str] | None = None,
    suffix: str = "_tertile",
) -> pd.DataFrame:
    """Add rank-based ordered tertile columns for every present mixing feature.

    For each feature column in ``features`` that exists in ``cluster_stats``, a
    ``{feature}{suffix}`` column is added with low->high ordered labels
    (default ``more_homogeneous`` / ``as_expected`` / ``more_mixed``). Columns
    not present are skipped silently.
    """
    features = list(MIXING_TERTILE_FEATURES if features is None else features)
    categories = list(
        MIXING_TERTILE_ORDER if categories is None else categories
    )
    out = cluster_stats.copy()
    for feature in features:
        if feature in out.columns:
            out[f"{feature}{suffix}"] = _ordered_tertile(
                out[feature], categories=categories
            )
    return out
