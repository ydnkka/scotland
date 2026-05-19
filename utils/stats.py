import numpy as np
import pandas as pd


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
