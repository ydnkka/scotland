import numpy as np
import pandas as pd


def safe_divide(numerator, denominator):
    """Divide safely, returning NaN where the denominator is zero or missing."""
    denominator_safe = pd.Series(denominator).replace(0, np.nan)
    out = pd.Series(numerator) / denominator_safe
    return out.replace([np.inf, -np.inf], np.nan)


def log1p_ratio(numerator, denominator):
    """Return ``log((1 + numerator) / (1 + denominator))``."""
    return np.log1p(numerator) - np.log1p(denominator)


def add_window_percentile(
    df: pd.DataFrame,
    col: str,
    groupby_cols="window_idx",
) -> pd.Series:
    """Percentile-rank a column within temporal windows or window strata."""
    return df.groupby(groupby_cols)[col].rank(pct=True, method="average")


def add_window_zscore(
    df: pd.DataFrame,
    col: str,
    groupby_cols="window_idx",
) -> pd.Series:
    """Z-score a column within temporal windows or window strata."""
    grouped = df.groupby(groupby_cols)[col]
    mean = grouped.transform("mean")
    std = grouped.transform("std", ddof=0).replace(0, np.nan)
    return (df[col] - mean) / std


def mean_with_count(
    df: pd.DataFrame,
    cols: list[str],
    *,
    skipna: bool = False,
) -> tuple[pd.Series, pd.Series]:
    """Compute a row-wise mean and a companion count of valid components."""
    existing_cols = [col for col in cols if col in df.columns]
    sub = df[existing_cols]
    return sub.mean(axis=1, skipna=skipna), sub.notna().sum(axis=1)


def attach_entropy_zscore(
    df: pd.DataFrame,
    cluster_col: str,
    category_col: str,
    window_col: str,
    *,
    n_random: int = 1000,
    normalise: bool = True,
    entropy_base: float = 2,
    random_state: int = 42,
    prefix: str = "cluster",
) -> pd.DataFrame:
    """
    Attach entropy-based cluster statistics and z-scores to ``df``.

    The null distribution is generated within each window stratum and cached
    by cluster size, so clusters with the same size in the same window reuse the
    same null mean and standard deviation.

    Parameters
    ----------
    df:
        Input dataframe.
    cluster_col:
        Column identifying clusters.
    category_col:
        Column containing the categorical variable used to compute entropy.
    window_col:
        Column identifying the temporal/window stratum for the null model.
    n_random:
        Number of random draws for each window-size null distribution.
    normalise:
        If True, entropy is divided by the global maximum entropy.
    entropy_base:
        Logarithm base for entropy.
    random_state:
        Random seed.
    prefix:
        Prefix for output columns.

    Returns
    -------
    pd.DataFrame
        Original dataframe with entropy statistics merged back on.

    Notes
    -----
    Singletons are retained:
        - observed entropy = 0
        - null mean = 0
        - null SD = 0
        - z-score = NaN
    """

    required_cols = [cluster_col, category_col, window_col]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in dataframe: {missing_cols}")

    if n_random < 2:
        raise ValueError("n_random must be at least 2 to estimate a null standard deviation.")

    if entropy_base <= 0 or entropy_base == 1:
        raise ValueError("entropy_base must be positive and not equal to 1.")

    out = df.copy()
    rng = np.random.default_rng(random_state)

    # Encode categories globally. Missing category values receive code -1.
    codes, categories = pd.factorize(out[category_col], sort=False)
    k_global = len(categories)

    if k_global <= 1:
        raise ValueError("Cannot compute entropy for a single category level.")

    log_base = np.log(entropy_base)
    max_entropy = np.log(k_global) / log_base if normalise else 1.0

    valid = codes >= 0

    work = pd.DataFrame(
        {
            window_col: out.loc[valid, window_col].to_numpy(),
            cluster_col: out.loc[valid, cluster_col].to_numpy(),
            "_cat_code": codes[valid],
        }
    )

    # If all category values are missing, return empty statistic columns.
    if work.empty:
        for col in [
            "n",
            "entropy_obs",
            "entropy_null_mean",
            "entropy_null_sd",
            "entropy_z",
            "is_singleton",
        ]:
            out[f"{prefix}_{col}"] = np.nan
        return out

    def entropy_from_counts(counts: np.ndarray) -> np.ndarray:
        """
        Compute entropy row-wise from a count matrix.

        Parameters
        ----------
        counts:
            Array of shape (n_rows, k_global), where each row contains category counts.

        Returns
        -------
        np.ndarray
            Entropy value for each row.
        """
        counts = counts.astype(float, copy=False)
        totals = counts.sum(axis=1, keepdims=True)

        probs = np.divide(
            counts,
            totals,
            out=np.zeros_like(counts, dtype=float),
            where=totals > 0,
        )

        log_probs = np.zeros_like(probs, dtype=float)
        mask = probs > 0
        log_probs[mask] = np.log(probs[mask])

        entropy = -(probs * log_probs).sum(axis=1) / log_base

        if normalise:
            entropy = entropy / max_entropy

        return entropy

    # Observed category-count matrix for each window-cluster pair.
    cluster_counts = (
        work.groupby([window_col, cluster_col, "_cat_code"], dropna=False, observed=True)
        .size()
        .unstack("_cat_code", fill_value=0)
    )

    # Ensure every global category code exists as a column.
    cluster_counts = cluster_counts.reindex(columns=np.arange(k_global), fill_value=0)

    count_matrix = cluster_counts.to_numpy()
    cluster_n = count_matrix.sum(axis=1).astype(int)
    observed_entropy = entropy_from_counts(count_matrix)

    stats = pd.DataFrame(
        {
            f"{prefix}_n": cluster_n,
            f"{prefix}_entropy_obs": observed_entropy,
            f"{prefix}_entropy_null_mean": np.nan,
            f"{prefix}_entropy_null_sd": np.nan,
            f"{prefix}_entropy_z": np.nan,
            f"{prefix}_is_singleton": cluster_n == 1,
        },
        index=cluster_counts.index,
    )

    # Generate null distributions once per window and cluster size.
    for window_id, window_cluster_counts in cluster_counts.groupby(
        level=0,
        dropna=False,
        sort=False,
    ):
        window_category_counts = window_cluster_counts.sum(axis=0).to_numpy(dtype=float)
        window_total = window_category_counts.sum()

        if window_total == 0:
            continue

        window_probs = window_category_counts / window_total

        window_index = window_cluster_counts.index
        window_sizes = window_cluster_counts.sum(axis=1).to_numpy(dtype=int)

        for size in np.unique(window_sizes):
            matching_index = window_index[window_sizes == size]

            if size <= 0:
                continue

            # Singletons have entropy 0 and a degenerate null distribution.
            # Their z-score is structurally undefined.
            if size == 1:
                stats.loc[matching_index, f"{prefix}_entropy_null_mean"] = 0.0
                stats.loc[matching_index, f"{prefix}_entropy_null_sd"] = 0.0
                stats.loc[matching_index, f"{prefix}_entropy_z"] = np.nan
                continue

            sampled_counts = rng.multinomial(
                n=size,
                pvals=window_probs,
                size=n_random,
            )

            null_entropies = entropy_from_counts(sampled_counts)

            null_mean = np.nanmean(null_entropies)
            null_sd = np.nanstd(null_entropies, ddof=1)

            stats.loc[matching_index, f"{prefix}_entropy_null_mean"] = null_mean
            stats.loc[matching_index, f"{prefix}_entropy_null_sd"] = null_sd

    # Compute z-scores. Singletons remain NaN because their null SD is 0.
    null_sd = stats[f"{prefix}_entropy_null_sd"].to_numpy()
    obs = stats[f"{prefix}_entropy_obs"].to_numpy()
    null_mean = stats[f"{prefix}_entropy_null_mean"].to_numpy()

    z = np.divide(
        obs - null_mean,
        null_sd,
        out=np.full_like(obs, np.nan, dtype=float),
        where=(null_sd > 0) & ~np.isnan(null_sd),
    )

    stats[f"{prefix}_entropy_z"] = z

    stats = stats.reset_index()

    out = out.merge(
        stats,
        on=[window_col, cluster_col],
        how="left",
    )

    return out


def downstream_entropy_fast(
    edge_df: pd.DataFrame,
    source_col: str,
    weight_col: str,
) -> pd.DataFrame:
    """
    Summarise how evenly each source node distributes shared sequences downstream.

    The entropy is computed over outgoing edge weights. A node with one observed
    successor has undefined downstream entropy because there is no branching
    pattern to measure.
    """
    required_cols = [source_col, weight_col]
    missing_cols = [col for col in required_cols if col not in edge_df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in edge dataframe: {missing_cols}")

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

    codes, nodes = pd.factorize(df[source_col], sort=False)
    order = np.argsort(codes, kind="mergesort")

    codes = codes[order]
    weights = df[weight_col].to_numpy(float)[order]

    starts = np.r_[0, np.flatnonzero(np.diff(codes)) + 1]
    counts = np.diff(np.r_[starts, len(weights)])

    strength = np.add.reduceat(weights, starts)
    wlogw = np.add.reduceat(weights * np.log(weights), starts)

    entropy = np.log(strength) - (wlogw / strength)
    entropy[counts <= 1] = np.nan

    entropy_norm = np.full(len(counts), np.nan, dtype=float)
    multi_successor = counts > 1
    entropy_norm[multi_successor] = (
        entropy[multi_successor] / np.log(counts[multi_successor])
    )

    dominant_frac = np.maximum.reduceat(weights, starts) / strength

    return pd.DataFrame(
        {
            "out_degree": counts,
            "out_strength": strength,
            "downstream_entropy": entropy,
            "downstream_entropy_norm": entropy_norm,
            "effective_successors": np.exp(entropy),
            "dominant_successor_frac": dominant_frac,
        },
        index=pd.Index(nodes, name=source_col),
    )


def downstream_spread_entropy(
    edge_df: pd.DataFrame,
    source_col: str,
    weight_col: str,
    *,
    n_perm: int = 199,
    seed: int = 42,
    metric_col: str = "downstream_entropy_norm",
) -> pd.DataFrame:
    """
    Attach a permutation null for the downstream entropy metric.

    The null preserves the observed outgoing degree of each source node and the
    global set of outgoing edge weights, then repeatedly permutes the weights
    across source-edge slots. This asks whether a node's downstream distribution
    is unusually even or unusually concentrated relative to contemporaneous edge
    weight heterogeneity.
    """
    if n_perm < 1:
        raise ValueError("n_perm must be at least 1.")

    df = edge_df[[source_col, weight_col]].copy()
    df[weight_col] = pd.to_numeric(df[weight_col], errors="coerce")
    df = df.dropna()
    df = df[df[weight_col] > 0]

    observed = downstream_entropy_fast(df, source_col, weight_col)
    null_cols = [
        f"{metric_col}_null_mean",
        f"{metric_col}_null_sd",
        f"{metric_col}_z",
        f"{metric_col}_p_high",
        f"{metric_col}_p_low",
    ]
    if observed.empty:
        for col in null_cols:
            observed[col] = pd.Series(dtype=float)
        return observed

    if metric_col not in observed.columns:
        raise ValueError(f"Unknown metric_col {metric_col!r}.")

    obs = observed[metric_col].to_numpy()
    obs_finite = np.isfinite(obs)

    codes, nodes = pd.factorize(df[source_col], sort=False)
    order = np.argsort(codes, kind="mergesort")

    codes = codes[order]
    weights = df[weight_col].to_numpy(float)[order]

    starts = np.r_[0, np.flatnonzero(np.diff(codes)) + 1]
    counts = np.diff(np.r_[starts, len(weights)])

    rng = np.random.default_rng(seed)

    n = np.zeros(len(nodes), dtype=int)
    mean = np.zeros(len(nodes), dtype=float)
    m2 = np.zeros(len(nodes), dtype=float)
    ge_obs = np.zeros(len(nodes), dtype=int)
    le_obs = np.zeros(len(nodes), dtype=int)

    for _ in range(n_perm):
        permuted_weights = rng.permutation(weights)

        strength = np.add.reduceat(permuted_weights, starts)
        wlogw = np.add.reduceat(permuted_weights * np.log(permuted_weights), starts)

        entropy = np.log(strength) - (wlogw / strength)
        entropy[counts <= 1] = np.nan

        vals = np.full(len(counts), np.nan, dtype=float)
        multi_successor = counts > 1
        vals[multi_successor] = (
            entropy[multi_successor] / np.log(counts[multi_successor])
        )

        finite = np.isfinite(vals)

        n[finite] += 1
        delta = vals[finite] - mean[finite]
        mean[finite] += delta / n[finite]
        m2[finite] += delta * (vals[finite] - mean[finite])

        comparable = finite & obs_finite
        ge_obs[comparable] += vals[comparable] >= obs[comparable]
        le_obs[comparable] += vals[comparable] <= obs[comparable]

    sd = np.full(len(nodes), np.nan)
    ok = n > 1
    sd[ok] = np.sqrt(m2[ok] / (n[ok] - 1))

    z = np.divide(
        obs - mean,
        sd,
        out=np.full(len(nodes), np.nan, dtype=float),
        where=(sd > 0) & np.isfinite(sd),
    )

    p_high = np.full(len(nodes), np.nan)
    p_low = np.full(len(nodes), np.nan)
    ok = n > 0
    p_high[ok] = (ge_obs[ok] + 1) / (n[ok] + 1)
    p_low[ok] = (le_obs[ok] + 1) / (n[ok] + 1)

    null = pd.DataFrame(
        {
            f"{metric_col}_null_mean": mean,
            f"{metric_col}_null_sd": sd,
            f"{metric_col}_z": z,
            f"{metric_col}_p_high": p_high,
            f"{metric_col}_p_low": p_low,
        },
        index=pd.Index(nodes, name=source_col),
    )

    return observed.join(null)


def add_sse_node_metrics(
    node_df: pd.DataFrame,
    *,
    window_col: str = "window_idx",
    lineage: str = "clade",
) -> pd.DataFrame:
    """
    Add lifecycle, amplification, onward-spread, and mixing metrics to nodes.

    The resulting columns are intended for node-level superspreading signature
    detection in a temporal cluster-transition graph. Percentiles are computed
    both within windows and within window-lifecycle strata.
    """
    df = node_df.copy()

    for col in ["in_degree", "in_strength", "out_degree", "out_strength"]:
        if col not in df.columns:
            df[col] = 0
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    df["birth"] = df["in_degree"].eq(0)
    df["birth_like"] = df["in_strength"].le(1)
    df["death"] = df["out_degree"].eq(0)
    df["continuation"] = df["in_degree"].gt(0) & df["out_degree"].gt(0)
    df["branching"] = df["out_degree"].gt(1)
    df["merging"] = df["in_degree"].gt(1)

    df["lifecycle"] = np.select(
        [df["birth"], df["death"], df["continuation"]],
        ["birth", "death", "continuation"],
        default="unknown",
    )

    df["simple_chain"] = df["continuation"] & ~df["branching"] & ~df["merging"]
    df["birth_branching"] = df["birth"] & df["branching"]
    df["expansion"] = df["continuation"] & df["branching"] & ~df["merging"]
    df["hub"] = df["branching"] & df["merging"]
    df["sink"] = df["merging"] & df["death"]
    df["isolated"] = df["birth"] & df["death"]

    first_window = df[window_col].min()
    last_window = df[window_col].max()
    df["left_censored"] = df["birth"] & df[window_col].eq(first_window)
    df["right_censored"] = df["death"] & df[window_col].eq(last_window)

    if lineage in df.columns:
        first_per_group = df.groupby(lineage)[window_col].transform("min")
        last_per_gruop = df.groupby(lineage)[window_col].transform("max")
        df["epoch_left_censored"] = df["birth"] & df[window_col].eq(first_per_group)
        df["epoch_right_censored"] = df["death"] & df[window_col].eq(last_per_gruop)
    else:
        df["epoch_left_censored"] = False
        df["epoch_right_censored"] = False

    df["degree_imbalance"] = df["out_degree"] - df["in_degree"]
    df["strength_imbalance"] = df["out_strength"] - df["in_strength"]
    df["log_strength_ratio"] = log1p_ratio(df["out_strength"], df["in_strength"])
    df["strength_ratio"] = safe_divide(df["out_strength"], df["in_strength"])

    df["net_amplification"] = df["cluster_size"] - df["in_strength"]
    df["upstream_novelty_proxy"] = df["cluster_size"] / (1 + df["in_strength"])

    df["downstream_retention_ratio"] = safe_divide(
        df["out_strength"],
        df["cluster_size"],
    )
    df["log_downstream_retention_ratio"] = log1p_ratio(
        df["out_strength"],
        df["cluster_size"],
    )
    df["downstream_expansion_proxy"] = (
        df["downstream_retention_ratio"] * df["out_degree"]
    )

    for area_col, density_col, mean_col in [
        ("cluster_n_datazones", "datazone_density", "mean_sequences_per_datazone"),
        ("local_authorities", "local_authority_density", "mean_sequences_per_local_authority"),
        ("health_boards", "health_board_density", "mean_sequences_per_health_board"),
    ]:
        if area_col in df.columns:
            df[density_col] = safe_divide(df[area_col], df["cluster_size"])
            df[mean_col] = safe_divide(df["cluster_size"], df[area_col])

    percentile_cols = [
        "cluster_size",
        "out_degree",
        "out_strength",
        "degree_imbalance",
        "strength_imbalance",
        "log_strength_ratio",
        "net_amplification",
        "upstream_novelty_proxy",
        "downstream_retention_ratio",
        "downstream_expansion_proxy",
        "local_authorities",
        "local_authority_density",
    ]

    for col in percentile_cols:
        if col not in df.columns:
            continue
        df[f"{col}_pct_window"] = add_window_percentile(
            df,
            col,
            groupby_cols=window_col,
        )
        df[f"{col}_pct_window_lifecycle"] = add_window_percentile(
            df,
            col,
            groupby_cols=[window_col, "lifecycle"],
        )

    zscore_cols = [
        "cluster_size",
        "in_strength",
        "out_strength",
        "net_amplification",
        "upstream_novelty_proxy",
        "downstream_expansion_proxy",
        "local_authorities",
    ]

    for col in zscore_cols:
        if col not in df.columns:
            continue
        df[f"{col}_z_window"] = add_window_zscore(
            df,
            col,
            groupby_cols=window_col,
        )

    core_components = [
        "cluster_size_pct_window",
        "net_amplification_pct_window",
        "upstream_novelty_proxy_pct_window",
    ]
    df["core_amplification_score"], df["core_amplification_n"] = mean_with_count(
        df,
        core_components,
        skipna=False,
    )

    onward_components = [
        "out_degree_pct_window",
        "out_strength_pct_window",
        "downstream_expansion_proxy_pct_window",
        "downstream_retention_ratio_pct_window",
    ]
    (
        df["onward_dissemination_score"],
        df["onward_dissemination_n"],
    ) = mean_with_count(df, onward_components, skipna=False)

    mixing_components = [
        "simd_entropy_z",
        "health_board_entropy_z",
        "age_entropy_z",
        "sex_entropy_z",
    ]
    df["mixing_score"], df["mixing_score_n"] = mean_with_count(
        df,
        mixing_components,
        skipna=False,
    )

    return df


def categorise_sse_nodes(
    df: pd.DataFrame,
    *,
    high_q: float = 0.90,
    very_high_q: float = 0.95,
    entropy_high: float = 0.65,
    entropy_low: float = 0.35,
    dominant_high: float = 0.75,
) -> pd.DataFrame:
    """
    Assign interpretable superspreading-signature categories to node metrics.

    The labels are heuristic signatures, not confirmed epidemiological events.
    They combine local amplification, graph lifecycle, downstream branching,
    spatial breadth, and population-mixing evidence.
    """
    out = df.copy()

    def has(col):
        return col in out.columns

    def bool_col(col):
        if has(col):
            return out[col].fillna(False).astype(bool)
        return pd.Series(False, index=out.index)

    def num_col(col, default=np.nan):
        if has(col):
            return pd.to_numeric(out[col], errors="coerce")
        return pd.Series(default, index=out.index, dtype=float)

    def high_col(col, q=high_q):
        if not has(col):
            return pd.Series(False, index=out.index)

        x = num_col(col)
        if col.endswith("_pct_window") or col.endswith("_pct_window_lifecycle"):
            return x >= q

        cutoff = x.quantile(q)
        return x >= cutoff

    def low_col(col, q=None):
        if not has(col):
            return pd.Series(False, index=out.index)

        q = 1 - high_q if q is None else q
        x = num_col(col)
        if col.endswith("_pct_window") or col.endswith("_pct_window_lifecycle"):
            return x <= q

        cutoff = x.quantile(q)
        return x <= cutoff

    birth_like = bool_col("birth_like")
    death = bool_col("death")
    continuation = bool_col("continuation")
    merging = bool_col("merging")
    isolated = bool_col("isolated")
    left_censored = bool_col("left_censored") | bool_col("epoch_left_censored")
    right_censored = bool_col("right_censored") | bool_col("epoch_right_censored")

    out_degree = num_col("out_degree", default=0).fillna(0)
    in_degree = num_col("in_degree", default=0).fillna(0)
    out_strength = num_col("out_strength", default=0).fillna(0)
    in_strength = num_col("in_strength", default=0).fillna(0)

    entropy_norm = num_col("downstream_entropy_norm")
    entropy_z = num_col("downstream_entropy_norm_z")
    entropy_p_high = num_col("downstream_entropy_norm_p_high")
    entropy_p_low = num_col("downstream_entropy_norm_p_low")
    dominant_frac = num_col("dominant_successor_frac")

    high_entropy = (
        (entropy_norm >= entropy_high)
        | ((entropy_z >= 1.96) & (entropy_p_high <= 0.05))
    )
    low_entropy = (
        (entropy_norm <= entropy_low)
        | (dominant_frac >= dominant_high)
        | ((entropy_z <= -1.96) & (entropy_p_low <= 0.05))
    )

    high_core_amp = high_col("core_amplification_score", high_q)
    very_large = high_col("cluster_size_pct_window_lifecycle", very_high_q)
    high_novelty = high_col("upstream_novelty_proxy_pct_window_lifecycle", high_q)
    high_net_amp = high_col("net_amplification_pct_window_lifecycle", high_q)
    high_downstream_expansion = high_col(
        "downstream_expansion_proxy_pct_window_lifecycle",
        high_q,
    )
    high_out_strength = high_col("out_strength_pct_window_lifecycle", high_q)

    spatially_broad = (
        high_col("local_authorities_pct_window_lifecycle", high_q)
        | high_col("local_authority_density_pct_window_lifecycle", high_q)
        | high_col("mixing_score", high_q)
        | (num_col("health_board_entropy_z") >= 1.96)
    )

    out["sse_candidate"] = (
        high_core_amp
        | (very_large & (high_novelty | high_net_amp))
        | (high_novelty & high_downstream_expansion)
        | (
            num_col("cluster_size_z_window").ge(2)
            & num_col("net_amplification_z_window").ge(1)
        )
    )

    out["sse_role"] = np.select(
        [
            isolated,
            (in_strength > 0) & (out_strength == 0),
            merging & (in_degree >= 2) & (out_strength > 0),
            continuation & (in_strength > 0) & (out_strength > 0) & high_net_amp,
            birth_like & high_novelty,
        ],
        [
            "isolated_burst",
            "terminal_sink",
            "merged_relay",
            "relay_amplifier",
            "putative_birth",
        ],
        default="unclear_origin",
    )

    out["sse_onward_dynamic"] = np.select(
        [
            out_strength == 0,
            death | ((out_degree <= 1) & low_col("downstream_expansion_proxy_pct_window_lifecycle", 0.25)),
            (out_degree >= 1) & low_entropy,
            (out_degree >= 2) & high_entropy & high_downstream_expansion & spatially_broad,
            (out_degree >= 2) & high_entropy & high_downstream_expansion,
            (out_degree >= 2) & high_entropy,
            (out_degree >= 2) & high_out_strength,
        ],
        [
            "no_observed_onward_spread",
            "contained_burst",
            "single_dominant_chain",
            "diffuse_spatial_broadcaster",
            "multi_branch_expander",
            "multi_branch_seeder",
            "high_volume_onward_spread",
        ],
        default="weak_or_ambiguous_onward_spread",
    )

    out["sse_category"] = np.where(
        ~out["sse_candidate"],
        "not_sse_like",
        out["sse_role"] + "__" + out["sse_onward_dynamic"],
    )

    out["sse_censoring_note"] = np.select(
        [
            left_censored & right_censored,
            left_censored,
            right_censored,
        ],
        [
            "both_left_and_right_censored",
            "left_censored_origin_uncertain",
            "right_censored_onward_uncertain",
        ],
        default="not_censored",
    )

    return out


def flag_sse(
    df: pd.DataFrame,
    threshold: float = 9.0,
    *,
    drop_incomplete_first_week: bool = False,
    drop_incomplete_last_week: bool = True,
) -> pd.DataFrame:
    """
    Compute weekly cluster growth metrics.

    Weeks are Monday-starting epidemiological weeks, running Monday to Sunday.

    Parameters
    ----------
    df : pd.DataFrame
        Input metadata containing sequence_id, collection_date, meta_cluster_id, and clade.

    threshold : float, default 9.0
        Threshold used to define SSE weeks from normalised change.

    drop_incomplete_first_week : bool, default False
        Whether to remove the first week if the dataset starts after Monday.

    drop_incomplete_last_week : bool, default True
        Whether to remove the final week if the dataset ends before Sunday.

    Returns
    -------
    weekly : pd.DataFrame
        Weekly cluster-level growth metrics.
    """

    required_cols = ["sequence_id", "collection_date", "meta_cluster_id", "clade"]
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        raise ValueError(f"Missing columns in dataframe: {missing_cols}")

    # Monday-starting weeks: Monday to Sunday.
    df["week"] = (
        df["collection_date"]
        .dt.to_period("W-SUN")
        .dt.start_time
    )

    min_date = df["collection_date"].min().normalize()
    max_date = df["collection_date"].max().normalize()

    first_week_start = min_date.to_period("W-SUN").start_time
    last_week_start = max_date.to_period("W-SUN").start_time

    first_complete_week_start = first_week_start
    last_complete_week_start = last_week_start

    # Optionally drop first week if the dataset starts after Monday.
    if drop_incomplete_first_week:
        if min_date > first_week_start:
            first_complete_week_start = first_week_start + pd.Timedelta(weeks=1)

    # Optionally drop last week if the dataset ends before Sunday.
    if drop_incomplete_last_week:
        last_week_end = last_week_start + pd.Timedelta(days=6)

        if max_date < last_week_end:
            last_complete_week_start = last_week_start - pd.Timedelta(weeks=1)

    if first_complete_week_start > last_complete_week_start:
        weekly = pd.DataFrame(
            columns=[
                "week",
                "new_sequences",
                "meta_cluster_id",
                "clade",
                "cc_size",
                "cc_size_prev",
                "norm_change",
                "is_sse",
            ]
        )
        ecs = pd.DataFrame()
        return weekly, ecs

    all_weeks = pd.date_range(
        start=first_complete_week_start,
        end=last_complete_week_start,
        freq="W-MON"
    )

    # Keep only full weeks included in all_weeks.
    df = df[df["week"].isin(all_weeks)].copy()

    if df.empty:
        weekly = pd.DataFrame(
            columns=[
                "week",
                "new_sequences",
                "meta_cluster_id",
                "clade",
                "cc_size",
                "cc_size_prev",
                "norm_change",
                "is_sse",
            ]
        )
        return weekly

    weekly_observed = (
        df
        .groupby(["meta_cluster_id", "week"], as_index=False)
        .agg(
            new_sequences=("sequence_id", "nunique"),
            clade=("clade", "first"),
        )
        .sort_values(["meta_cluster_id", "week"])
    )

    def reindex_cluster(g: pd.DataFrame) -> pd.DataFrame:
        cluster_id = g.name
        clade_val = g["clade"].iloc[0]

        out = (
            g.set_index("week")[["new_sequences"]]
            .reindex(all_weeks, fill_value=0)
        )

        out["meta_cluster_id"] = cluster_id
        out["clade"] = clade_val

        return out

    weekly = (
        weekly_observed
        .groupby("meta_cluster_id", group_keys=False)
        .apply(reindex_cluster)
        .reset_index()
        .rename(columns={"index": "week"})
    )

    weekly = weekly.sort_values(["meta_cluster_id", "week"]).reset_index(drop=True)

    weekly["cc_size"] = (
        weekly
        .groupby("meta_cluster_id")["new_sequences"]
        .cumsum()
    )

    weekly["cc_size_prev"] = (
        weekly
        .groupby("meta_cluster_id")["cc_size"]
        .shift(1)
        .fillna(0)
    )

    weekly["norm_change"] = (
        weekly["new_sequences"] /
        np.sqrt(weekly["cc_size_prev"].clip(lower=1))
    )

    weekly["is_sse"] = weekly["norm_change"] > threshold

    return weekly
