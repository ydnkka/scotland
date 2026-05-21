"""Statistical machinery for the superspreading-signature detection pipeline.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


__all__ = [
    "add_sse_node_metrics",
    "categorise_sse_nodes",
    "flag_sse",
    "frequencies",
    "safe_mode",
    "test_category_distribution",
]


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
    groupby_cols: Any ="window_idx",
) -> pd.Series:
    """Percentile-rank a column within temporal windows"""
    return df.groupby(groupby_cols)[col].rank(pct=True, method="average")


def add_window_zscore(
    df: pd.DataFrame,
    col: str,
    groupby_cols: Any = "window_idx",
) -> pd.Series:
    """Z-score a column within temporal windows"""
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


def frequencies(values: pd.Series) -> str:
    counts = values.dropna().astype(str).value_counts()
    return "; ".join(f"{name} ({count})" for name, count in counts.items())


def safe_mode(values: pd.Series):
    values = values.dropna()
    if values.empty:
        return np.nan
    return values.mode().iloc[0]


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

    if "cluster_size" not in df.columns:
        raise ValueError(
            "add_sse_node_metrics requires a 'cluster_size' column on the node "
            "dataframe (used by net_amplification, retention ratios, and "
            "per-area densities)."
        )
    df["cluster_size"] = pd.to_numeric(df["cluster_size"], errors="coerce")

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
        last_per_group = df.groupby(lineage)[window_col].transform("max")
        df["epoch_left_censored"] = df["birth"] & df[window_col].eq(first_per_group)
        df["epoch_right_censored"] = df["death"] & df[window_col].eq(last_per_group)
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
        ("n_local_authorities", "local_authority_density", "mean_sequences_per_local_authority"),
        ("n_health_boards", "health_board_density", "mean_sequences_per_health_board"),
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
        "n_health_boards",
        "health_board_density",
        "n_local_authorities",
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
        "n_local_authorities",
        "n_health_boards",
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
    df["onward_dissemination_score"], df["onward_dissemination_n"] = mean_with_count(df, onward_components, skipna=False)


    mixing_components = [
        "sex_entropy_z",
        "age_entropy_z",
        "simd_entropy_z",
        "health_board_entropy_z",
        "urban_rural_entropy_z",
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
    dominant_frac = num_col("dominant_successor_frac")
    high_entropy = entropy_norm >= entropy_high
    low_entropy = (
        (entropy_norm <= entropy_low) | (dominant_frac >= dominant_high)
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
        high_col("n_local_authorities_pct_window_lifecycle", high_q)
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

    # Work on a copy so we never mutate the caller's frame.
    df = df.copy()

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

    empty_weekly = pd.DataFrame(
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

    if first_complete_week_start > last_complete_week_start:
        return empty_weekly

    all_weeks = pd.date_range(
        start=first_complete_week_start,
        end=last_complete_week_start,
        freq="W-MON"
    )

    # Keep only full weeks included in all_weeks.
    df = df[df["week"].isin(all_weeks)].copy()

    if df.empty:
        return empty_weekly

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


def test_category_distribution(df):
    from scipy.stats import chi2_contingency, fisher_exact
    from statsmodels.stats.multitest import multipletests

    df = df.copy()

    table = (
        df
        .pivot_table(
            index="q",
            columns="candidate",
            values="n",
            aggfunc="sum",
            fill_value=0,
        )
        .rename(columns={False: "background", True: "candidate"})
    )

    chi2, p, dof, expected = chi2_contingency(table)

    expected_df = pd.DataFrame(
        expected,
        index=table.index,
        columns=table.columns,
    )

    std_resid = (table - expected_df) / np.sqrt(expected_df)

    props = table.div(table.sum(axis=0), axis=1)
    props["diff_candidate_minus_background"] = (
        props["candidate"] - props["background"]
    )

    total_background = table["background"].sum()
    total_candidate = table["candidate"].sum()

    per_category = []

    for q, row in table.iterrows():
        candidate_in_q = row["candidate"]
        background_in_q = row["background"]

        test_table = np.array([
            [candidate_in_q, total_candidate - candidate_in_q],
            [background_in_q, total_background - background_in_q],
        ])

        odds_ratio, p_value = fisher_exact(test_table)

        per_category.append({
            "q": q,
            "candidate_n": candidate_in_q,
            "background_n": background_in_q,
            "candidate_frac": candidate_in_q / total_candidate,
            "background_frac": background_in_q / total_background,
            "diff_candidate_minus_background": (
                candidate_in_q / total_candidate
                - background_in_q / total_background
            ),
            "odds_ratio": odds_ratio,
            "p_value": p_value,
        })

    per_category = pd.DataFrame(per_category)

    per_category["p_adj_bh"] = multipletests(
        per_category["p_value"],
        method="fdr_bh",
    )[1]

    per_category = per_category.sort_values("p_adj_bh")

    return {
        "table": table,
        "chi2": chi2,
        "p_value": p,
        "dof": dof,
        "expected": expected_df,
        "standardized_residuals": std_resid,
        "proportions": props,
        "per_category_tests": per_category,
    }
