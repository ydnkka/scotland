"""
Entropy utilities for cluster-category and edge-weight analyses.

Public API
----------
shannon_entropy
    Shannon entropy from dense count arrays.
shannon_entropy_grouped
    Per-group Shannon entropy from ragged positive-valued data.
max_entropy
    log_base(k); the maximum entropy of a k-category distribution.
cluster_socio_demo_entropy
    Cluster category-entropy with within-window null model and z-scores.
downstream_edge_entropy
    Per-source entropy of outgoing edge weights.

All entropies use log base 2 by default (units: bits).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "max_entropy",
    "shannon_entropy",
    "shannon_entropy_grouped",
    "cluster_socio_demo_entropy",
    "downstream_edge_entropy",
]


def max_entropy(k, base: float = 2):
    """Maximum Shannon entropy for ``k`` equally likely categories.

    Returns ``log_base(k)``. Vectorises over array ``k``.
    """
    return np.log(k) / np.log(base)


def shannon_entropy(
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
        Entropy per slice. Slices summing to zero return 0.
    
    Notes
    -----
    Follows the ``0 * log(0) = 0`` convention.
    """
    if base <= 0 or base == 1:
        raise ValueError("base must be positive and not equal to 1.")
    
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


def shannon_entropy_grouped(
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
        unique_codes : np.ndarray
            Group labels in ascending order.
        size : np.ndarray
            Number of elements per group.
        total : np.ndarray
            Sum of values per group.
        entropy : np.ndarray
            Shannon entropy per group, in given base. NaN where size <= 1.
        max_value : np.ndarray
            Largest value per group (for concentration measures).
    """
    if base <= 0 or base == 1:
        raise ValueError("base must be positive and not equal to 1.")
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
    """Attach cluster category-entropy statistics and within-window z-scores.

    The null model draws cluster-sized multinomial samples from the window-
    pooled category frequencies. Null draws are cached per ``(window, size)``,
    so clusters of the same size in the same window share null mean / SD.

    Singletons retain ``entropy = 0``, ``null_mean = 0``, ``null_sd = 0``,
    and ``z = NaN`` (degenerate null distribution).

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe.
    cluster_col : str
        Column identifying clusters.
    category_col : str
        Column containing the categorical variable.
    window_col : str
        Column identifying the window stratum for the null model.
    n_random : int, default 1000
        Multinomial draws per (window, size) null distribution.
    normalise : bool, default True
        If True, entropies are divided by ``max_entropy(k_global, base)``.
    base : float, default 2
        Logarithm base (default 2 -> bits).
    random_state : int, default 42
        RNG seed.
    prefix : str, default "cluster"
        Column prefix for outputs.

    Returns
    -------
    pd.DataFrame
        ``df`` with the following columns merged on (window_col, cluster_col):
        ``{prefix}_n``, ``{prefix}_entropy_obs``,
        ``{prefix}_entropy_null_mean``, ``{prefix}_entropy_null_sd``,
        ``{prefix}_entropy_z``.

    Notes
    -----
    ``{prefix}_n`` counts non-missing category observations per cluster, not
    the cluster's total row count. Clusters whose category values are entirely
    missing receive NaN for all statistic columns.
    """
    required_cols = [cluster_col, category_col, window_col]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in dataframe: {missing}")
    if n_random < 2:
        raise ValueError("n_random must be >= 2 to estimate null SD.")
    if base <= 0 or base == 1:
        raise ValueError("base must be positive and not equal to 1.")
    
    out = df.copy()
    rng = np.random.default_rng(random_state)

    codes, categories = pd.factorize(out[category_col], sort=False)
    k_global = len(categories)
    if k_global <= 1:
        raise ValueError("Cannot compute entropy for a single category level.")
    
    norm_factor = max_entropy(k_global, base=base) if normalise else 1.0

    stat_cols = [
        f"{prefix}_n",
        f"{prefix}_entropy_obs",
        f"{prefix}_entropy_null_mean",
        f"{prefix}_entropy_null_sd",
        f"{prefix}_entropy_z",
    ]

    valid = codes >= 0
    work = pd.DataFrame({
        window_col: out.loc[valid, window_col].to_numpy(),
        cluster_col: out.loc[valid, cluster_col].to_numpy(),
        "_cat_code": codes[valid],
    })

    if work.empty:
        for col in stat_cols:
            out[col] = np.nan
        return out
    
    # Observed counts: rows = (window, cluster), columns = category codes.
    cluster_counts = (
        work.groupby(
            [window_col, cluster_col, "_cat_code"],
            dropna=False,
            observed=True,
        )
        .size()
        .unstack("_cat_code", fill_value=0)
        .reindex(columns=np.arange(k_global), fill_value=0)
    )

    count_matrix = cluster_counts.to_numpy()
    cluster_n = count_matrix.sum(axis=1).astype(int)
    observed_entropy = shannon_entropy(count_matrix, base=base) / norm_factor

    null_mean = np.full(len(cluster_n), np.nan, dtype=float)
    null_sd = np.full(len(cluster_n), np.nan, dtype=float)

    # Position lookup so we can write into null arrays by integer row index.
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
                # Degenerate null. z structurally undefined.
                null_mean[target_rows] = 0.0
                null_sd[target_rows] = 0.0
                continue

            sampled = rng.multinomial(n=size, pvals=window_probs, size=n_random)
            null_H = shannon_entropy(sampled, base=base) / norm_factor
            null_mean[target_rows] = null_H.mean()
            null_sd[target_rows] = null_H.std(ddof=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(
            (null_sd > 0) & ~np.isnan(null_sd),
            (observed_entropy - null_mean) / null_sd,
            np.nan,
        )

    stats = pd.DataFrame(
        {
            f"{prefix}_n": cluster_n,
            f"{prefix}_entropy_obs": observed_entropy,
            f"{prefix}_entropy_null_mean": null_mean,
            f"{prefix}_entropy_null_sd": null_sd,
            f"{prefix}_entropy_z": z,
        },
        index=cluster_counts.index,
    ).reset_index()

    out = out.merge(stats, on=[window_col, cluster_col], how="left")

    return out


def downstream_edge_entropy(
    edge_df: pd.DataFrame,
    source_col: str,
    weight_col: str,
    *,
    base: float = 2,
) -> pd.DataFrame:
    """Per-source entropy of outgoing edge weights.

    Measures how evenly each source node distributes weight across its
    successors. Sources with a single successor have entropy NaN (no
    branching pattern to measure).

    Parameters
    ----------
    edge_df : pd.DataFrame
        Long-format edge dataframe.
    source_col : str
        Column identifying the source node.
    weight_col : str
        Numeric column of outgoing weights. Non-positive values are dropped.
    base : float, default 2
        Logarithm base (default 2 -> bits).

    Returns
    -------
    pd.DataFrame
        Indexed by source node, with columns:
        ``out_degree``, ``out_strength``, ``downstream_entropy``,
        ``downstream_entropy_norm``, ``effective_successors``,
        ``dominant_successor_frac``.

    Notes
    -----
    ``downstream_entropy_norm`` divides by ``log_base(out_degree)`` (a local
    reference, in [0, 1] for nodes with >= 2 successors).
    ``effective_successors`` is the Hill number of order 1, equal to
    ``base ** entropy`` regardless of base.
    """
    required_cols = [source_col, weight_col]
    missing = [c for c in required_cols if c not in edge_df.columns]
    if missing:
        raise ValueError(f"Missing columns in edge dataframe: {missing}")
    if base <= 0 or base == 1:
        raise ValueError("base must be positive and not equal to 1.")
    
    df = edge_df[[source_col, weight_col]].copy()
    df[weight_col] = pd.to_numeric(df[weight_col], errors="coerce")
    df = df.dropna()
    df = df[df[weight_col] > 0]

    out_cols = [
        "out_degree",
        "out_strength",
        "downstream_entropy",
        "downstream_entropy_norm",
        "effective_successors",
        "dominant_successor_frac",
    ]
    if df.empty:
        return pd.DataFrame(columns=out_cols, index=pd.Index([], name=source_col))
    
    codes, nodes = pd.factorize(df[source_col].to_numpy(), sort=False)

    agg = shannon_entropy_grouped(
        df[weight_col].to_numpy(dtype=float),
        codes,
        base=base,
    )

    sizes = agg["size"]
    H = agg["entropy"]

    # Local normalisation: divide by log_base(out_degree).
    H_norm = np.full(len(sizes), np.nan, dtype=float)
    multi = sizes > 1
    H_norm[multi] = H[multi] / max_entropy(sizes[multi], base=base)
    H_norm = np.clip(H_norm, 0.0, 1.0)  # tame FP overshoot

    # Effective successors = base ** H (Hill number of order 1).
    effective = np.power(float(base), H)

    dominant = agg["max_value"] / agg["total"]

    node_labels = nodes[agg["unique_codes"]]

    return pd.DataFrame(
        {
            "out_degree": sizes,
            "out_strength": agg["total"],
            "downstream_entropy": H,
            "downstream_entropy_norm": H_norm,
            "effective_successors": effective,
            "dominant_successor_frac": dominant,
        },
        index=pd.Index(node_labels, name=source_col),
    )
