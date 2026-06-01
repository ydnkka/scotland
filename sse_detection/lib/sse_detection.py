#!/usr/bin/env python3
"""
Superspreading Signature Detection
==================================

Detect node-level superspreading signatures in a temporal cluster-transition
network. Each node is a genetically close, temporally proximate SARS-CoV-2
sequence cluster; directed edges link clusters in adjacent retained windows and
are weighted by shared sequence counts.

Inputs
------
- ``data/processed/scotland_clustering_analysis_dataset.parquet``

Usage
-----
    python sse_detection.py
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from typing import Sequence

import networkx as nx
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from utils import load_analysis_columns  # noqa: E402
import entropy as e  # noqa: E402

RANDOM_SEED = 42

N_ENTROPY_DRAWS = 1000

ANALYSIS_COLUMNS = [
    "window_id",
    "window_idx",
    "wn_start_date",
    "wn_mid_date",
    "wn_end_date",
    "wn_positive_tests",
    "wn_no_sequences",
    "wn_prop_sequenced",
    "cluster_id",
    "cluster_size",
    "cluster_n_datazones",
    "cluster_duration_days",
    "age_band",
    "sex",
    "is_female",
    "who_voc",
    "clade",
    "is_vaccinated",
    "vacc_dose_number",
    "vacc_booster",
    "days_since_vaccination",
    "test_reason",
    "is_reinfection",
    "datazone",
    "dz_care_home_tests",
    "dz_urban_rural_class",
    "dz_health_board",
    "dz_local_authority",
    "dz_simd_quintile",
    "dz_7d_test_positivity",
    "dz_cum_sequences",
    "dz_cum_positive_tests",
    "dz_cum_prop_sequenced",
    "dz_cum_incidence_per_capita",
    "dz_cum_vaccinated",
    "dz_cum_prop_vaccinated",
]

COHERENCE_OBS_COLS = [
    # "sex_entropy_obs",
    # "age_entropy_obs",
    # "urban_rural_entropy_obs",
    "simd_entropy_obs",
    "health_board_entropy_obs",
]


def load_sequence_data() -> pd.DataFrame:
    """Load sequence-window rows.

    Raw windows are overlapping 3-week windows stepped weekly; ``window_stride=2``
    keeps alternate windows by sorted-window position and reindexes them, giving
    a less redundant 2-week graph step while preserving overlap between adjacent
    retained windows.
    """
    df = load_analysis_columns(
        ANALYSIS_COLUMNS,
        add_policy=True,
        window_stride=2,
    )
    return df


def build_transition_network(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, nx.DiGraph, dict]:
    """Build the directed adjacent-window cluster-transition graph.

    Returns the edge table, the raw graph, and a node -> weakly-connected
    component map.
    """
    seq_nodes = (
        df[["sequence_id", "cluster_id", "window_id", "window_idx"]]
        .drop_duplicates()
        .sort_values(["sequence_id", "window_idx"])
    )

    edge_map: dict[tuple, set[str]] = defaultdict(set)
    for sequence_id, group in seq_nodes.groupby("sequence_id", sort=False):
        records = list(
            group[["cluster_id", "window_id", "window_idx"]].itertuples(
                index=False, name=None
            )
        )
        for i, (source_node, source_window, source_idx) in enumerate(records):
            for target_node, target_window, target_idx in records[i + 1 :]:
                delta = target_idx - source_idx
                if delta == 1 and source_node != target_node:
                    edge_key = (
                        source_node,
                        target_node,
                        source_window,
                        target_window,
                        source_idx,
                        target_idx,
                    )
                    edge_map[edge_key].add(str(sequence_id))
                elif delta > 1:
                    break

    edge_rows = [
        {
            "source": source,
            "target": target,
            "source_window_id": source_window,
            "target_window_id": target_window,
            "source_window_idx": source_idx,
            "target_window_idx": target_idx,
            "n_shared_sequences": len(shared),
        }
        for (
            source,
            target,
            source_window,
            target_window,
            source_idx,
            target_idx,
        ), shared in edge_map.items()
    ]

    edge_table = pd.DataFrame(
        edge_rows,
        columns=[
            "source",
            "target",
            "source_window_id",
            "target_window_id",
            "source_window_idx",
            "target_window_idx",
            "n_shared_sequences",
        ],
    )
    if not edge_table.empty:
        edge_table = edge_table.sort_values(
            [
                "source_window_idx",
                "target_window_idx",
                "n_shared_sequences",
                "source",
                "target",
            ],
            ascending=[True, True, False, True, True],
            ignore_index=True,
        )

    g_raw = nx.DiGraph()
    g_raw.add_nodes_from(seq_nodes["cluster_id"].unique())
    for row in edge_table.itertuples(index=False):
        g_raw.add_edge(row.source, row.target, weight=row.n_shared_sequences)

    components = sorted(
        nx.weakly_connected_components(g_raw),
        key=lambda nodes: (-len(nodes), sorted(nodes)[0]),
    )
    component_map = {
        node: f"CC{component_idx:05d}"
        for component_idx, nodes in enumerate(components, start=1)
        for node in nodes
    }

    return edge_table, g_raw, component_map


def build_cluster_stats(
    df: pd.DataFrame,
    edge_table: pd.DataFrame,
    *,
    n_entropy_draws: int = N_ENTROPY_DRAWS,
    random_state: int = RANDOM_SEED,
) -> pd.DataFrame:
    """Compute cluster-level mixing-entropy stats and edge-flow metrics.

    All entropy calls now pass ``n_random`` explicitly. Edge metrics are joined
    so that ALL clusters are retained (left joins); nodes with no in/out edges
    keep NaN here and are filled to 0 downstream.
    """
    sex_stats = e.cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "sex",
        "window_idx",
        prefix="sex",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    age_stats = e.cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "age_band",
        "window_idx",
        prefix="age",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    simd_stats = e.cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "dz_simd_quintile",
        "window_idx",
        prefix="simd",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    hb_stats = e.cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "dz_health_board",
        "window_idx",
        prefix="health_board",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    ur_stats = e.cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "dz_urban_rural_class",
        "window_idx",
        prefix="urban_rural",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    vacc_stats = e.vaccination_mixing_features(
        df,
        n_random=n_entropy_draws,
        random_state=random_state,
    )

    cluster_stats = (
        sex_stats[["cluster_id", "sex_entropy_z", "sex_entropy_obs"]]
        .merge(
            age_stats[["cluster_id", "age_entropy_z", "age_entropy_obs"]],
            on="cluster_id",
            how="outer",
        )
        .merge(
            simd_stats[["cluster_id", "simd_entropy_z", "simd_entropy_obs"]],
            on="cluster_id",
            how="outer",
        )
        .merge(
            hb_stats[
                ["cluster_id", "health_board_entropy_z", "health_board_entropy_obs"]
            ],
            on="cluster_id",
            how="outer",
        )
        .merge(
            ur_stats[
                ["cluster_id", "urban_rural_entropy_z", "urban_rural_entropy_obs"]
            ],
            on="cluster_id",
            how="outer",
        )
        .merge(
            vacc_stats[
                ["cluster_id", "vaccination_entropy_z", "vaccination_entropy_obs"]
            ],
            on="cluster_id",
            how="outer",
        )
    )

    cluster_stats = e.add_mixing_tertiles(cluster_stats)

    edge_stats = e.onward_edge_entropy(
        edge_table,
        source_col="source",
        weight_col="n_shared_sequences",
    )

    in_metrics = (
        edge_table.groupby("target", dropna=False)
        .agg(
            in_degree=("source", "nunique"),
            in_strength=("n_shared_sequences", "sum"),
        )
        .rename_axis("cluster_id")
        .reset_index()
    )

    # LEFT joins: keep every cluster, even pure sinks/sources with no edges.
    edge_stats = edge_stats.merge(in_metrics, on="cluster_id", how="outer")
    cluster_stats = cluster_stats.merge(edge_stats, on="cluster_id", how="left")
    return cluster_stats


def safe_mode(values: pd.Series):
    """Return the first modal non-missing value, or NaN when none exists."""
    values = values.dropna()
    if values.empty:
        return np.nan
    return values.mode().iloc[0]


def build_cluster_attributes(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse sequence-window rows to one row per (cluster, window)."""
    agg = {
        "wn_start_date": ("wn_start_date", "first"),
        "wn_mid_date": ("wn_mid_date", "first"),
        "wn_end_date": ("wn_end_date", "first"),
        "wn_positive_tests": ("wn_positive_tests", "first"),
        "wn_no_sequences": ("wn_no_sequences", "first"),
        "wn_prop_sequenced": ("wn_prop_sequenced", "first"),
        "cluster_size": ("cluster_size", "first"),
        "duration_days": ("cluster_duration_days", "first"),
        "first_collection_date": ("collection_date", "min"),
        "last_collection_date": ("collection_date", "max"),
        "who_voc": ("who_voc", "first"),
        "clade": ("clade", "first"),
        "pango_lineage": ("pango_lineage", "first"),
        "policy_period": ("policy_period", safe_mode),
        "policy_era": ("policy_era", safe_mode),
        "dz_7d_test_positivity": ("dz_7d_test_positivity", "mean"),
        "dz_cum_sequences": ("dz_cum_sequences", "mean"),
        "dz_cum_incidence_per_capita": ("dz_cum_incidence_per_capita", "mean"),
        "dz_cum_positive_tests": ("dz_cum_positive_tests", "mean"),
        "dz_cum_prop_sequenced": ("dz_cum_prop_sequenced", "mean"),
    }
    cluster_att = (
        df.groupby(["cluster_id", "window_id", "window_idx"], as_index=False)
        .agg(**agg)
        .sort_values(["window_idx", "cluster_size"], ascending=[True, False])
        .reset_index(drop=True)
    )
    return cluster_att


def build_cluster_table(
    cluster_att: pd.DataFrame,
    cluster_stats: pd.DataFrame,
    edge_table: pd.DataFrame,
) -> pd.DataFrame:
    """Join attributes + stats and attach downstream-burden metrics.

    LEFT joins throughout so the full node set survives to scoring.
    """
    cluster_table = cluster_att.merge(cluster_stats, on="cluster_id", how="left")

    target_sizes = cluster_table[["cluster_id", "cluster_size"]].rename(
        columns={"cluster_id": "target", "cluster_size": "_target_cluster_size"}
    )
    downstream = (
        edge_table.merge(target_sizes, on="target", how="left")
        .groupby("source", dropna=False)
        .agg(
            downstream_cluster_burden=("_target_cluster_size", "sum"),
            mean_successor_cluster_size=("_target_cluster_size", "mean"),
        )
        .rename_axis("cluster_id")
        .reset_index()
    )
    cluster_table = cluster_table.merge(downstream, on="cluster_id", how="left")
    cluster_table["downstream_cluster_burden"] = cluster_table[
        "downstream_cluster_burden"
    ].fillna(0)
    return cluster_table


def safe_divide(numerator, denominator) -> pd.Series:
    """Divide safely, returning NaN where the denominator is zero or missing."""
    denom = pd.Series(denominator).replace(0, np.nan)
    out = pd.Series(numerator) / denom
    return out.replace([np.inf, -np.inf], np.nan)


def log1p_ratio(numerator, denominator) -> pd.Series:
    """Return ``log((1 + numerator) / (1 + denominator))``."""
    return np.log1p(numerator) - np.log1p(denominator)


def mean_with_count(
    df: pd.DataFrame, cols: Sequence[str]
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
) -> np.ndarray:
    """
    Compute upper-tail empirical p-values from a simulated null distribution.

    For each observation i, this estimates:

        p_i = P(null_score_i >= observed_score_i)

    using the finite permutation distribution contained in ``null_scores``.

    A +1 smoothing correction is used:

        p_i = (1 + number of valid null scores >= obs_i) / (1 + number of valid null scores)

    This avoids returning p = 0 when none of the simulated null scores are as
    extreme as the observed score. Such a zero would be misleading because the
    permutation distribution is finite.

    Parameters
    ----------
    obs:
        One-dimensional array of observed scores, shape ``(n,)``.

    null_scores:
        Two-dimensional array of null scores, shape ``(n, B)``, where ``B`` is
        the number of permutations.

    Returns
    -------
    np.ndarray
        One empirical upper-tail p-value per observation, shape ``(n,)``.

        Returns ``NaN`` when the observed score is ``NaN`` or when no valid null
        scores are available for that row.

    Notes
    -----
    This is an upper-tail p-value. Small values mean that the observed score is
    unusually large compared with the null distribution.
    """
    obs = np.asarray(obs, dtype=float)
    null_scores = np.asarray(null_scores, dtype=float)

    if null_scores.ndim != 2:
        raise ValueError("null_scores must be a 2D array with shape (n, B).")

    if obs.ndim != 1:
        raise ValueError("obs must be a 1D array with shape (n,).")

    if null_scores.shape[0] != obs.shape[0]:
        raise ValueError("obs and null_scores must have the same number of rows.")

    # Only compare against valid null scores. Without this, NaNs would behave as
    # False in comparisons and could make p-values artificially small.
    valid_null = ~np.isnan(null_scores)
    n_valid = valid_null.sum(axis=1)

    ge = np.sum((null_scores >= obs[:, None]) & valid_null, axis=1)

    p = (1 + ge) / (1 + n_valid)

    # A p-value is undefined if the observed value is missing or if no valid
    # null scores exist for that observation.
    p[np.isnan(obs) | (n_valid == 0)] = np.nan

    return p


def choose_permutation_strata(
    df: pd.DataFrame,
    *,
    strata_cols: Sequence[str],
    min_stratum_n: int = 20,
) -> pd.Series:
    """
    Construct adaptive permutation-stratum labels.

    The aim is to define groups within which values can be permuted while
    preserving important contextual structure. The function starts with a broad
    grouping based on the first available stratum column, then progressively
    tries more detailed groupings by adding further columns.

    For example, with:

        strata_cols = ["window", "clade", "health_board"]

    candidate keys are tried in this order:

        window
        window|clade
        window|clade|health_board

    A row is moved to a more detailed key only if that candidate key has at
    least ``min_stratum_n`` rows. Otherwise, the row keeps its previous broader
    key.

    Parameters
    ----------
    df:
        Input dataframe.

    strata_cols:
        Ordered list of candidate columns used to build permutation strata.
        Earlier columns define broader context; later columns define finer
        context.

    min_stratum_n:
        Minimum size required before a candidate finer stratum is accepted.

    Returns
    -------
    pd.Series
        String-valued stratum key for each row in ``df``.

    Notes
    -----
    This function does not guarantee that every final stratum has at least
    ``min_stratum_n`` rows.

    Why? Some rows from a broad parent group may move into finer groups, leaving
    a small residual broad group behind. The function should therefore be read
    as:

        "Use the finest available key whose candidate group is large enough."

    rather than:

        "Guarantee all final groups are at least min_stratum_n."

    If no requested stratum columns are available, all rows are assigned to a
    single stratum called ``"all"``.
    """
    if min_stratum_n < 1:
        raise ValueError("min_stratum_n must be at least 1.")

    available = [col for col in strata_cols if col in df.columns]

    if not available:
        return pd.Series("all", index=df.index, dtype="string")

    # Start with the broadest available grouping.
    best_key = df[available[:1]].astype("string").agg("|".join, axis=1)

    # Try increasingly detailed keys. A row adopts the finer key only if that
    # finer candidate stratum is large enough.
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
    """
    Check whether every permutation stratum lies inside one ranking group.

    This is used to decide whether the faster permutation route is valid.

    The fast route is valid when each stratum maps to exactly one
    ``rank_within`` group. In that case, permuting values within strata can
    never move a value across ranking groups. Therefore, the set of values
    inside each ranking group is unchanged by permutation.

    This means percentile ranking and permutation commute:

        rank first, then permute

    gives the same null distribution as:

        permute first, then rank

    Parameters
    ----------
    strata_codes:
        Integer-coded permutation strata, shape ``(n,)``.

    rank_codes:
        Integer-coded ranking groups, shape ``(n,)``.

    Returns
    -------
    bool
        ``True`` if every stratum is nested within a single ranking group;
        otherwise ``False``.
    """
    strata_codes = np.asarray(strata_codes)
    rank_codes = np.asarray(rank_codes)

    if strata_codes.shape[0] != rank_codes.shape[0]:
        raise ValueError("strata_codes and rank_codes must have the same length.")

    # Sort rows by stratum so that each stratum appears as a contiguous run.
    order = np.argsort(strata_codes, kind="stable")

    s = strata_codes[order]
    r = rank_codes[order]

    # A boundary marks the first row of each new stratum.
    is_boundary = np.r_[True, s[1:] != s[:-1]]

    # For rows continuing an existing stratum run, check that the ranking group
    # is the same as the previous row's ranking group.
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
) -> pd.DataFrame:
    """
    Add permutation-based null statistics for a composite node score.

    The observed score is assumed to have been built from percentile-ranked raw
    components within ``rank_within`` groups. This function constructs a null
    score in the same way, but after permuting the raw component information
    within adaptive contextual strata.

    Null modes
    ----------
    ``null_mode="profile"`` is the recommended default for SSE prioritisation.
    It permutes the whole row-wise component profile within each stratum. This
    preserves the correlation structure among components such as cluster size,
    onward burden and downstream burden. It tests whether an intact multivariate
    SSE-like profile is unusually assigned to a node's context.

    ``null_mode="independent"`` permutes each component independently within
    each stratum. This preserves each component's within-stratum marginal
    distribution but breaks component co-occurrence. It is a more disruptive
    null and asks whether the observed combination of high component values is
    surprising if components were independently recombined.

    Parameters
    ----------
    df:
        Input dataframe containing the observed score, raw components, ranking
        groups and stratification variables.

    score_col:
        Name of the observed composite score column.

    raw_component_cols:
        Raw component columns used to generate the null composite score. Missing
        component columns are ignored.

    rank_within:
        Column defining the groups within which raw components are converted to
        percentile ranks.

    strata_cols:
        Ordered list of columns used to construct adaptive permutation strata.

    min_stratum_n:
        Minimum size required before a finer adaptive stratum is accepted.

    n_permutations:
        Number of null permutations.

    random_state:
        Seed used for reproducible permutations.

    null_mode:
        Either ``"profile"`` or ``"independent"``.

    Returns
    -------
    pd.DataFrame
        Copy of ``df`` with four additional columns:

        ``{score_col}_null_mean``
            Mean null score for each row.

        ``{score_col}_null_sd``
            Standard deviation of the null scores for each row.

        ``{score_col}_null_z``
            Z-score comparing the observed score with the null distribution.

        ``{score_col}_upper_p``
            Upper-tail empirical p-value.

    Notes
    -----
    The function uses a fast path when permutation strata are nested inside the
    ranking groups. In that case, percentile ranks are computed once and then
    permuted directly. Otherwise, raw values are permuted and re-ranked on every
    permutation.
    """
    if n_permutations < 2:
        raise ValueError("n_permutations must be at least 2.")

    if null_mode not in {"profile", "independent"}:
        raise ValueError("null_mode must be either 'profile' or 'independent'.")

    required_cols = [score_col, rank_within]
    missing_required = [col for col in required_cols if col not in df.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {missing_required}")

    out = df.copy()

    null_cols = [
        f"{score_col}_null_mean",
        f"{score_col}_null_sd",
        f"{score_col}_null_z",
        f"{score_col}_upper_p",
    ]

    components = [col for col in raw_component_cols if col in out.columns]
    if not components:
        for col in null_cols:
            out[col] = np.nan
        return out

    rng = np.random.default_rng(random_state)

    # ------------------------------------------------------------------
    # 1. Build adaptive permutation strata.
    # ------------------------------------------------------------------
    strata = choose_permutation_strata(
        out,
        strata_cols=strata_cols,
        min_stratum_n=min_stratum_n,
    )

    # Use positional indices, not dataframe index labels, because the matrices
    # below are NumPy arrays.
    strata_values = strata.to_numpy()
    strata_groups = [
        np.flatnonzero(strata_values == stratum) for stratum in pd.unique(strata_values)
    ]

    # ------------------------------------------------------------------
    # 2. Encode ranking groups and prepare raw component matrix.
    # ------------------------------------------------------------------
    rank_codes = out[rank_within].astype("category").cat.codes.to_numpy()
    unique_rank_codes = np.unique(rank_codes)

    raw = out[components].to_numpy(dtype=float)
    n_rows, n_components = raw.shape
    null_scores = np.empty((n_rows, n_permutations), dtype=float)

    def pct_rank_within(values_2d: np.ndarray) -> np.ndarray:
        """Percentile-rank each component within ``rank_within`` groups."""
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
        """
        Permute a component matrix within adaptive strata.

        In ``profile`` mode, rows are shuffled as intact multivariate profiles.
        In ``independent`` mode, each component column is shuffled separately.
        """
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

    # ------------------------------------------------------------------
    # 3. Generate null scores.
    # ------------------------------------------------------------------
    strata_codes = pd.factorize(strata, sort=False)[0]
    nested = _strata_nested_in_rank_groups(strata_codes, rank_codes)

    if nested:
        # Fast path: ranking groups retain the same value multiset after
        # within-stratum permutation, so ranking once is valid.
        base_pct = pct_rank_within(raw)

        for b in range(n_permutations):
            permuted_pct = permute_matrix_within_strata(base_pct)
            with np.errstate(invalid="ignore"):
                null_scores[:, b] = np.nanmean(permuted_pct, axis=1)
    else:
        # General path: permutation can change the value multiset inside ranking
        # groups, so raw values must be re-ranked after each permutation.
        for b in range(n_permutations):
            permuted_raw = permute_matrix_within_strata(raw)
            permuted_pct = pct_rank_within(permuted_raw)
            with np.errstate(invalid="ignore"):
                null_scores[:, b] = np.nanmean(permuted_pct, axis=1)

    # ------------------------------------------------------------------
    # 4. Summarise the null distribution.
    # ------------------------------------------------------------------
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
    out[f"{score_col}_upper_p"] = empirical_upper_tail_p(obs, null_scores)

    return out


def add_sse_node_metrics(
    df: pd.DataFrame,
    *,
    min_stratum_n: int = 20,
    n_permutations: int = 500,
    random_state: int = RANDOM_SEED,
    min_cluster_size: int = 6,
    null_mode: str = "profile",
) -> pd.DataFrame:
    """
    Add SSE-candidate metrics to a node-level cluster table.

    The SSE score is computed only for size-eligible nodes, defined by
    ``cluster_size >= min_cluster_size``. Smaller nodes are retained in the
    output but marked as ``size_ineligible`` and are not used for percentile
    ranking, null calibration or candidate tiering.

    The observed composite score is constructed by:

    1. deriving raw size, novelty, onward-flow and coherence components;
    2. percentile-ranking those components within ``window_idx``;
    3. averaging the available percentile-ranked components.

    A permutation-calibrated null is then added. By default, ``null_mode`` is
    ``"profile"``, which shuffles whole multivariate component profiles within
    adaptive strata. This preserves correlation among score components and is
    therefore more conservative than independently shuffling each component.

    Parameters
    ----------
    df:
        Node-level cluster table, usually produced by ``build_cluster_table``.

    min_stratum_n:
        Minimum candidate stratum size used for context-adjusted expected
        cluster size and adaptive permutation strata.

    n_permutations:
        Number of permutations used for composite-score null calibration.

    random_state:
        Random seed used for reproducible permutation results.

    min_cluster_size:
        Minimum cluster size required for a node to be tested as an SSE
        candidate.

    null_mode:
        Either ``"profile"`` or ``"independent"``. ``"profile"`` is the
        recommended default for calibrated prioritisation. ``"independent"`` is
        useful as a stronger co-occurrence-disruption sensitivity analysis.

    Returns
    -------
    pd.DataFrame
        Copy of ``df`` with additional SSE feature columns, sub-scores,
        permutation-null statistics and candidate review tiers.
    """
    if null_mode not in {"profile", "independent"}:
        raise ValueError("null_mode must be either 'profile' or 'independent'.")

    required_cols = ["cluster_size", "window_idx"]
    missing_required = [col for col in required_cols if col not in df.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {missing_required}")

    out = df.copy()
    out["sse_tested"] = out["cluster_size"].ge(min_cluster_size)
    out["candidate_tier"] = "size_ineligible"

    edge_default_cols = {
        "in_degree": 0,
        "out_degree": 0,
        "in_strength": 0,
        "out_strength": 0,
        "downstream_cluster_burden": 0,
    }
    for col, default in edge_default_cols.items():
        if col not in out.columns:
            out[col] = default
        else:
            out[col] = out[col].fillna(default)

    tested = out.loc[out["sse_tested"]].copy()
    if tested.empty:
        return out

    # ------------------------------------------------------------------
    # 1. Context-adjusted excess cluster size.
    # ------------------------------------------------------------------
    expected = pd.Series(
        tested["cluster_size"].mean(),
        index=tested.index,
        dtype=float,
    )

    context_cols = [col for col in ["window_idx", "clade"] if col in tested.columns]

    for n_cols in range(1, len(context_cols) + 1):
        cols = context_cols[:n_cols]
        grouped = tested.groupby(cols, dropna=False)["cluster_size"]
        group_mean = grouped.transform("mean")
        group_n = grouped.transform("size")
        use_group = group_n >= min_stratum_n
        expected.loc[use_group] = group_mean.loc[use_group]

    tested["expected_cluster_size"] = expected.clip(lower=0.5)
    tested["sampling_adjusted_excess_size"] = (
        tested["cluster_size"] - tested["expected_cluster_size"]
    ) / np.sqrt(tested["expected_cluster_size"])

    # ------------------------------------------------------------------
    # 2. Graph-position role indicators
    # ------------------------------------------------------------------
    in_degree = tested["in_degree"].fillna(0)
    out_degree = tested["out_degree"].fillna(0)

    # Basic directional flags. These are intentionally simple and may overlap.
    tested["has_incoming"] = in_degree.gt(0)
    tested["has_outgoing"] = out_degree.gt(0)
    tested["is_source_boundary"] = in_degree.eq(0)
    tested["is_sink_boundary"] = out_degree.eq(0)
    tested["has_branching"] = out_degree.gt(1)
    tested["has_merging"] = in_degree.gt(1)

    # Mutually exclusive primary graph role.
    # Each node receives exactly one role.
    tested["primary_graph_role"] = "other"

    tested.loc[
        in_degree.eq(0) & out_degree.eq(0),
        "primary_graph_role",
    ] = "isolated"

    tested.loc[
        in_degree.eq(0) & out_degree.eq(1),
        "primary_graph_role",
    ] = "single_outgoing_source"

    tested.loc[
        in_degree.eq(0) & out_degree.gt(1),
        "primary_graph_role",
    ] = "source_branching"

    tested.loc[
        in_degree.eq(1) & out_degree.eq(0),
        "primary_graph_role",
    ] = "single_incoming_sink"

    tested.loc[
        in_degree.gt(1) & out_degree.eq(0),
        "primary_graph_role",
    ] = "merging_sink"

    tested.loc[
        in_degree.eq(1) & out_degree.eq(1),
        "primary_graph_role",
    ] = "linear_continuation"

    tested.loc[
        in_degree.eq(1) & out_degree.gt(1),
        "primary_graph_role",
    ] = "internal_branching"

    tested.loc[
        in_degree.gt(1) & out_degree.eq(1),
        "primary_graph_role",
    ] = "internal_merging"

    tested.loc[
        in_degree.gt(1) & out_degree.gt(1),
        "primary_graph_role",
    ] = "merge_and_branch"

    first_window = tested["window_idx"].min()
    last_window = tested["window_idx"].max()

    tested["left_censored"] = tested["is_source_boundary"] & tested["window_idx"].eq(first_window)
    tested["right_censored"] = tested["is_sink_boundary"] & tested["window_idx"].eq(last_window)
    tested["boundary_censored"] = tested["left_censored"] | tested["right_censored"]

    # ------------------------------------------------------------------
    # 3. Size, novelty and onward-flow proxies.
    # ------------------------------------------------------------------
    tested["log_cluster_size"] = np.log1p(tested["cluster_size"])
    tested["log_excess_over_upstream"] = log1p_ratio(
        tested["cluster_size"],
        tested["in_strength"],
    )
    tested["novelty_fraction"] = safe_divide(
        tested["cluster_size"] - tested["in_strength"],
        tested["cluster_size"],
    ).clip(lower=0, upper=1)
    tested["onward_retention_ratio"] = safe_divide(
        tested["out_strength"],
        tested["cluster_size"],
    )
    tested["onward_expansion_proxy"] = (
        tested["onward_retention_ratio"] * tested["out_degree"]
    )
    tested["log_downstream_burden"] = np.log1p(tested["downstream_cluster_burden"])
    tested["log_downstream_burden_ratio"] = log1p_ratio(
        tested["downstream_cluster_burden"],
        tested["cluster_size"],
    )

    # ------------------------------------------------------------------
    # 4. Mixing / coherence scores.
    # ------------------------------------------------------------------
    obs_cols = [col for col in COHERENCE_OBS_COLS if col in tested.columns]
    if obs_cols:
        observed_coherence = 1 - tested[obs_cols]
        tested["coherence_score"] = observed_coherence.mean(axis=1, skipna=True)
        tested["coherence_score_n"] = observed_coherence.notna().sum(axis=1)
    else:
        tested["coherence_score"] = np.nan
        tested["coherence_score_n"] = 0

    # ------------------------------------------------------------------
    # 5. Within-window percentile ranks.
    # ------------------------------------------------------------------
    rank_cols = [
        "sampling_adjusted_excess_size",
        "log_cluster_size",
        "log_excess_over_upstream",
        "novelty_fraction",
        "out_strength",
        "out_degree",
        "effective_successors",
        "onward_entropy_norm",
        "onward_expansion_proxy",
        "log_downstream_burden",
        "coherence_score",
    ]

    for col in rank_cols:
        tested[f"{col}_pct_window"] = tested.groupby("window_idx", dropna=False)[
            col
        ].rank(pct=True, method="average")

    # ------------------------------------------------------------------
    # 6. Interpretable sub-scores.
    # ------------------------------------------------------------------
    focused_components = [
        "sampling_adjusted_excess_size_pct_window",
        "log_excess_over_upstream_pct_window",
        "novelty_fraction_pct_window",
        "coherence_score_pct_window",
    ]
    tested["focused_amplification_score"], tested["focused_amplification_n"] = (
        mean_with_count(tested, focused_components)
    )

    onward_components = [
        "out_strength_pct_window",
        "effective_successors_pct_window",
        "onward_entropy_norm_pct_window",
        "log_downstream_burden_pct_window",
    ]
    tested["onward_burden_score"], tested["onward_burden_n"] = mean_with_count(
        tested,
        onward_components,
    )

    # ------------------------------------------------------------------
    # 7. Composite SSE score and null calibration.
    # ------------------------------------------------------------------
    composite_raw_components = [
        "sampling_adjusted_excess_size",
        "log_cluster_size",
        "log_excess_over_upstream",
        "effective_successors",
        "onward_expansion_proxy",
        "log_downstream_burden",
        "coherence_score",
    ]
    composite_pct_components = [f"{col}_pct_window" for col in composite_raw_components]

    tested["composite_sse_score"], tested["composite_sse_score_n"] = mean_with_count(
        tested,
        composite_pct_components,
    )

    tested = add_composite_null_scores(
        tested,
        score_col="composite_sse_score",
        raw_component_cols=composite_raw_components,
        rank_within="window_idx",
        strata_cols=context_cols,
        min_stratum_n=min_stratum_n,
        n_permutations=n_permutations,
        random_state=random_state,
        null_mode=null_mode,
    )

    # ------------------------------------------------------------------
    # 8. Candidate ranking and review tiers.
    # ------------------------------------------------------------------
    tested["candidate_rank"] = tested["composite_sse_score_null_z"].rank(
        ascending=False,
        method="min",
        na_option="bottom",
    )

    enough_sequences = tested["cluster_size"].ge(min_cluster_size)
    high_priority = tested["composite_sse_score_upper_p"].le(0.05) & enough_sequences
    possible_review = tested["composite_sse_score_upper_p"].le(0.10) & enough_sequences
    high_score_uncalibrated = tested["composite_sse_score"].ge(0.90) & enough_sequences

    tested["candidate_tier"] = np.select(
        [
            high_priority,
            possible_review,
            high_score_uncalibrated,
        ],
        [
            "high_priority_review",
            "possible_review",
            "high_score_uncalibrated",
        ],
        default="background_or_low_information",
    )

    # Assign tested results back safely. This creates any new columns in the
    # full output while preserving size-ineligible rows.
    for col in tested.columns:
        out.loc[tested.index, col] = tested[col]

    return out


def main():
    print("...........loading data")
    df = load_sequence_data()

    print("...........building network")
    edge_table, _, component_map = build_transition_network(df)

    print("...........aggregating sequence data")
    cluster_stats = build_cluster_stats(df, edge_table)
    cluster_att = build_cluster_attributes(df)
    cluster_att["connected_components"] = cluster_att["cluster_id"].map(component_map)
    cluster_table = build_cluster_table(cluster_att, cluster_stats, edge_table)

    print("...........detecting SSEs")
    sse_df = add_sse_node_metrics(cluster_table)

    out_path = PROJECT_ROOT / "sse_detection/results/sse_outputs"
    out_path.mkdir(exist_ok=True)
    sse_df.to_parquet(out_path / "cluster_table")
    edge_table.to_parquet(out_path / "edge_table")

    print("Done.")


if __name__ == "__main__":
    main()
