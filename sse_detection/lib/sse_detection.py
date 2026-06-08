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
from entropy import (  # noqa: E402
    cluster_socio_demo_entropy,
    vaccination_mixing_features,
    add_mixing_tertiles,
    onward_edge_entropy,
)

RANDOM_SEED = 42
N_ENTROPY_DRAWS = 1000
N_PERMUTATIONS = 1000
MIN_CLUSTER_SIZE = 6

# Onward-spread SHAPE (dispersal: onward_entropy_norm and dominant_successor_frac)
# is NOT a detection axis. It was evaluated as one and demoted to
# characterisation: ~91% of the candidates it would uniquely contribute are
# out_degree==2 nodes, where the entropy is near-binary (balanced vs lopsided
# pair) and its rank-based p-value is poorly calibrated. The raw shape features
# are retained to characterise how detected candidates propagate (multiplicative
# fan-out vs concentrated hand-off), but they do not flag candidates.

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
    sex_stats = cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "sex",
        "window_idx",
        prefix="sex",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    age_stats = cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "age_band",
        "window_idx",
        prefix="age",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    simd_stats = cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "dz_simd_quintile",
        "window_idx",
        prefix="simd",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    dz_stats = cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "datazone",
        "window_idx",
        prefix="datazone",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    la_stats = cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "dz_local_authority",
        "window_idx",
        prefix="local_authority",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    hb_stats = cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "dz_health_board",
        "window_idx",
        prefix="health_board",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    ur_stats = cluster_socio_demo_entropy(
        df,
        "cluster_id",
        "dz_urban_rural_class",
        "window_idx",
        prefix="urban_rural",
        n_random=n_entropy_draws,
        random_state=random_state,
    )
    vacc_stats = vaccination_mixing_features(
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
            dz_stats[["cluster_id", "datazone_entropy_z", "datazone_entropy_obs"]],
            on="cluster_id",
            how="outer",
        )
        .merge(
            la_stats[
                [
                    "cluster_id",
                    "local_authority_entropy_z",
                    "local_authority_entropy_obs",
                ]
            ],
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

    cluster_stats = add_mixing_tertiles(cluster_stats)

    edge_stats = onward_edge_entropy(
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


def build_new_downstream_metrics(
    sequence_df: pd.DataFrame,
    edge_table: pd.DataFrame,
    *,
    min_shared_sequences: int = 2,
) -> pd.DataFrame:
    """Compute overlap-adjusted downstream burden for each source node.

    ``new_downstream_burden`` counts unique sequences in all adjacent-window
    child clusters that are not already present in the source cluster. The
    supported version applies the same calculation after excluding weak edges
    with fewer than ``min_shared_sequences`` shared sequences.
    """
    if min_shared_sequences < 1:
        raise ValueError("min_shared_sequences must be at least 1.")

    stat_cols = [
        "new_downstream_burden",
        "supported_new_downstream_burden",
        "new_downstream_children",
        "supported_new_downstream_children",
        "mean_successor_new_sequences",
    ]

    required_edge_cols = {"source", "target", "n_shared_sequences"}
    missing_edge_cols = sorted(required_edge_cols - set(edge_table.columns))
    if missing_edge_cols:
        raise ValueError(f"Missing edge columns: {missing_edge_cols}")

    required_sequence_cols = {"cluster_id", "sequence_id"}
    missing_sequence_cols = sorted(required_sequence_cols - set(sequence_df.columns))
    if missing_sequence_cols:
        raise ValueError(f"Missing sequence columns: {missing_sequence_cols}")

    if edge_table.empty:
        return pd.DataFrame(columns=["cluster_id", *stat_cols])

    membership = sequence_df[["cluster_id", "sequence_id"]].dropna().drop_duplicates()
    cluster_sequences = (
        membership.groupby("cluster_id", sort=False)["sequence_id"]
        .agg(lambda values: frozenset(values))
        .to_dict()
    )

    rows = []
    for source, edges in edge_table.groupby("source", sort=False, dropna=False):
        source_sequences = cluster_sequences.get(source, frozenset())
        new_sequences: set = set()
        supported_new_sequences: set = set()
        child_new_counts: list[int] = []
        new_child_count = 0
        supported_new_child_count = 0

        for edge in edges.itertuples(index=False):
            target_sequences = cluster_sequences.get(edge.target, frozenset())
            child_new_sequences = target_sequences.difference(source_sequences)
            child_new_count = len(child_new_sequences)
            child_new_counts.append(child_new_count)

            if child_new_count > 0:
                new_child_count += 1
                new_sequences.update(child_new_sequences)

            if edge.n_shared_sequences >= min_shared_sequences:  # type: ignore
                if child_new_count > 0:
                    supported_new_child_count += 1
                    supported_new_sequences.update(child_new_sequences)

        rows.append(
            {
                "cluster_id": source,
                "new_downstream_burden": len(new_sequences),
                "supported_new_downstream_burden": len(supported_new_sequences),
                "new_downstream_children": new_child_count,
                "supported_new_downstream_children": supported_new_child_count,
                "mean_successor_new_sequences": (
                    float(np.mean(child_new_counts)) if child_new_counts else np.nan
                ),
            }
        )

    return pd.DataFrame(rows, columns=["cluster_id", *stat_cols])


def build_cluster_table(
    cluster_att: pd.DataFrame,
    cluster_stats: pd.DataFrame,
    edge_table: pd.DataFrame,
    sequence_df: pd.DataFrame,
    *,
    min_supported_shared_sequences: int = 2,
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

    new_downstream = build_new_downstream_metrics(
        sequence_df,
        edge_table,
        min_shared_sequences=min_supported_shared_sequences,
    )
    cluster_table = cluster_table.merge(
        new_downstream,
        on="cluster_id",
        how="left",
    )
    fill_zero_cols = [
        "new_downstream_burden",
        "supported_new_downstream_burden",
        "new_downstream_children",
        "supported_new_downstream_children",
    ]
    for col in fill_zero_cols:
        cluster_table[col] = cluster_table[col].fillna(0)

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
        One-dimensional array of observed scores, dispersal ``(n,)``.

    null_scores:
        Two-dimensional array of null scores, dispersal ``(n, B)``, where ``B`` is
        the number of permutations.

    Returns
    -------
    np.ndarray
        One empirical upper-tail p-value per observation, dispersal ``(n,)``.

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
        raise ValueError("null_scores must be a 2D array with dispersal (n, B).")

    if obs.ndim != 1:
        raise ValueError("obs must be a 1D array with dispersal (n,).")

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
        Integer-coded permutation strata, dispersal ``(n,)``.

    rank_codes:
        Integer-coded ranking groups, dispersal ``(n,)``.

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
    n_permutations: int = N_PERMUTATIONS,
    random_state: int = RANDOM_SEED,
    min_cluster_size: int = MIN_CLUSTER_SIZE,
    null_mode: str = "profile",
) -> pd.DataFrame:
    """
    Add SSE-candidate metrics to a node-level cluster table.

    Detection rests on two demographic-blind, near-orthogonal calibrated axes:

    * **burst** — local accumulation: anomalously large / high-excess clusters,
        scored for every size-eligible node (``cluster_size >= min_cluster_size``);
    * **burden** — per-capita onward propagation: clusters that seed an
        unusually large downstream burden RELATIVE TO THEIR SIZE, scored only on
        the propagating subpopulation (clusters with non-zero onward burden), so
        the ~93% non-propagating mass does not dominate the ranking. Uses
        ratio/per-capita components only, so it is empirically orthogonal to
        burst rather than a restatement of size.

    The shape of onward spread (``onward_entropy_norm``,
    ``dominant_successor_frac``) was evaluated as a third "dispersal" axis but
    demoted to characterisation: it is dominated by two-successor nodes, where it
    is near-binary and poorly calibrated. It is retained raw to describe how a
    candidate propagates, not to flag candidates.

    Socio-demographic and geographic concentration features (the entropy-z
    coherence/mixing columns), raw onward counts (``effective_successors``,
    ``out_strength``), the onward-shape features and ``novelty_fraction`` are
    deliberately NOT used in any detection axis. They are retained on the output
    as descriptive characterisation features for the archetype classifier, so
    that the candidate set is not biased by demographic sampling.

    Smaller nodes are retained in the output but marked ``size_ineligible`` and
    are not used for percentile ranking, null calibration or tiering.

    Each observed axis score is the mean of its components' within-window
    percentile ranks; a profile-permutation null is then added per axis. By
    default ``null_mode="profile"`` shuffles whole multivariate component
    profiles within adaptive strata, preserving component correlation.

    Parameters
    ----------
    df:
        Node-level cluster table, usually produced by ``build_cluster_table``.

    min_stratum_n:
        Minimum candidate stratum size used for context-adjusted expected
        cluster size and adaptive permutation strata.

    n_permutations:
        Number of permutations used for axis-score null calibration.

    random_state:
        Random seed used for reproducible permutation results.

    min_cluster_size:
        Minimum cluster size required for a node to be tested at all.

    null_mode:
        Either ``"profile"`` or ``"independent"``. ``"profile"`` is the
        recommended default for calibrated prioritisation. ``"independent"`` is
        useful as a stronger co-occurrence-disruption sensitivity analysis.

    Returns
    -------
    pd.DataFrame
        Copy of ``df`` with additional SSE feature columns, sub-scores,
        per-axis permutation-null statistics, descriptive coherence features and
        candidate review tiers.
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

    basic_edge_default_cols = {
        "in_degree": 0,
        "out_degree": 0,
        "in_strength": 0,
        "out_strength": 0,
        "downstream_cluster_burden": 0,
    }
    for col, default in basic_edge_default_cols.items():
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
    # Cluster sizes are heavy-tailed (power-law-like), so a Poisson-style
    # standardisation (size - mean) / sqrt(mean) badly understates the spread
    # and is dominated by a few huge clusters. We instead work in log space,
    # where the bulk is approximately log-normal, and standardise each cluster
    # against the location and scale of its (window, clade) context. The mean is
    # also replaced implicitly by a log-space mean computed only over strata
    # large enough to estimate it, with a coarser-context backoff.
    tested["log_cluster_size"] = np.log1p(tested["cluster_size"])

    context_cols = [col for col in ["window_idx", "clade"] if col in tested.columns]

    # Global fallback location/scale (used until a large-enough stratum overrides).
    expected_log = pd.Series(
        tested["log_cluster_size"].mean(), index=tested.index, dtype=float
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
        # Require enough rows BOTH to trust the location and to estimate a SD
        # (std needs >= 2 and is noisy when thin, so reuse the stratum-size gate).
        use_group = (group_n >= min_stratum_n) & group_std.gt(0) & group_std.notna()
        expected_log.loc[use_group] = group_mean.loc[use_group]
        scale_log.loc[use_group] = group_std.loc[use_group]

    scale_log = scale_log.clip(lower=1e-6)

    # Report the back-transformed expectation for interpretability, and the
    # log-space robust z-score as the excess feature actually used downstream.
    tested["expected_cluster_size"] = np.maximum(np.expm1(expected_log), 0.5)
    tested["sampling_adjusted_excess_size"] = (
        tested["log_cluster_size"] - expected_log
    ) / scale_log

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

    tested["left_censored"] = tested["is_source_boundary"] & tested["window_idx"].eq(
        first_window
    )
    tested["right_censored"] = tested["is_sink_boundary"] & tested["window_idx"].eq(
        last_window
    )
    tested["boundary_censored"] = tested["left_censored"] | tested["right_censored"]

    # ------------------------------------------------------------------
    # 3. Size, novelty and onward-flow proxies.
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # 4. Within-window percentile ranks (detection components only).
    # ------------------------------------------------------------------
    # coherence_score is intentionally absent: it is a descriptive feature, not
    # a detection component.
    rank_cols = [
        "sampling_adjusted_excess_size",
        "log_cluster_size",
        "log_excess_over_upstream",
        "log_new_downstream_burden_ratio",
        "log_supported_new_downstream_burden_ratio",
    ]

    for col in rank_cols:
        if col not in tested.columns:
            continue
        tested[f"{col}_pct_window"] = tested.groupby("window_idx", dropna=False)[
            col
        ].rank(pct=True, method="average")

    # ------------------------------------------------------------------
    # 5. Detection axes and null calibration.
    # ------------------------------------------------------------------
    # Two near-orthogonal calibrated detection axes, each scored on the
    # population where it is actually defined:
    #   burst  (universal)                    — local accumulation / intensity
    #   burden (non-zero-burden subset, ~7%)  — per-capita onward volume
    # Pairwise (on the burden-eligible subset): burst~burden ~ -0.2 — distinct.
    # Burden is a SUBPOPULATION axis: the network is ~93% non-propagating, so
    # scoring burden across all nodes would rank a ~93% zero-mass and is
    # meaningless. It is scored and calibrated only on its eligible subset, with
    # the rest flagged "burden_not_applicable".

    # BURST: local accumulation, context-adjusted. Scored for every tested node.
    # Most size-dominated axis, so carries the small-cluster floor behaviour at
    # p ~ 1 (expected, conservative). sampling_adjusted_excess_size is kept (the
    # most epidemiologically faithful "large for its context" component); its
    # mild negative correlation with onward signals is acceptable.
    burst_components = [
        "log_cluster_size",
        "sampling_adjusted_excess_size",
        "log_excess_over_upstream",
    ]

    # BURDEN: per-capita onward volume on the propagating subpopulation. Ratio
    # forms only (raw counts would reintroduce size). Eligible = non-zero on
    # either burden ratio (a terminal cluster has no burden to be anomalous about).
    burden_components = [
        "log_new_downstream_burden_ratio",
        "log_supported_new_downstream_burden_ratio",
    ]

    axis_strata_cols = [c for c in ["window_idx", "clade"] if c in tested.columns]

    # --- BURST: universal axis, calibrated across all tested nodes. ---
    burst_pct = [f"{c}_pct_window" for c in burst_components]
    tested["burst_score"], tested["burst_score_n"] = mean_with_count(
        tested, burst_pct
    )
    tested = add_composite_null_scores(
        tested,
        score_col="burst_score",
        raw_component_cols=burst_components,
        rank_within="window_idx",
        strata_cols=axis_strata_cols,
        min_stratum_n=min_stratum_n,
        n_permutations=n_permutations,
        random_state=random_state,
        null_mode=null_mode,
    )

    # --- BURDEN: subpopulation axis on nodes with non-zero onward burden. ---
    burden_present = [c for c in burden_components if c in tested.columns]
    burden_nonzero = pd.Series(False, index=tested.index)
    for c in burden_present:
        burden_nonzero = burden_nonzero | (tested[c].fillna(0) > 0)
    tested["burden_eligible"] = burden_nonzero

    burden_idx = tested.index[tested["burden_eligible"]]
    for col in (
        "burden_score",
        "burden_score_n",
        "burden_score_null_mean",
        "burden_score_null_sd",
        "burden_score_null_z",
        "burden_score_upper_p",
    ):
        tested[col] = np.nan

    if len(burden_idx) >= 2:
        burden_sub = tested.loc[burden_idx].copy()
        # Re-rank within the burden-eligible subpopulation (not the whole graph),
        # so ranks reflect the population the axis is defined on.
        for c in burden_present:
            burden_sub[f"{c}_burden_pct"] = burden_sub.groupby(
                "window_idx", dropna=False
            )[c].rank(pct=True, method="average")
        burden_sub_pct = [f"{c}_burden_pct" for c in burden_present]
        burden_sub["burden_score"], burden_sub["burden_score_n"] = mean_with_count(
            burden_sub, burden_sub_pct
        )
        # Sparse subpopulation: calibrate with coarse strata (clade else pooled),
        # not per-window.
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
        )
        for col in (
            "burden_score",
            "burden_score_n",
            "burden_score_null_mean",
            "burden_score_null_sd",
            "burden_score_null_z",
            "burden_score_upper_p",
        ):
            tested.loc[burden_idx, col] = burden_sub[col]

    # ---- onward-spread SHAPE: characterisation only, NOT a detection axis -----
    # Dispersal (the shape of onward spread) was evaluated as a third detection
    # axis but demoted to characterisation. The signal is dominated by
    # out_degree == 2 nodes, where onward_entropy_norm is near-binary (balanced
    # vs lopsided pair); its rank-based p-value is then lumpy and poorly
    # calibrated, and ~91% of the candidates it would uniquely contribute are
    # exactly these tie-prone two-successor nodes. Rather than detect on a signal
    # we cannot calibrate, we retain the raw onward-shape descriptors —
    # onward_entropy_norm (evenness of spread) and dominant_successor_frac
    # (concentration in the top successor) — as characterisation features for the
    # archetype classifier (multiplicative fan-out vs concentrated hand-off). They
    # are already present on the table and play no part in scoring or tiering.

    # ------------------------------------------------------------------
    # 6. Candidate ranking and review tiers (two detection axes).
    # ------------------------------------------------------------------
    # Two near-orthogonal calibrated detection axes: burst (local size/excess,
    # universal) and burden (per-capita onward propagation, propagating subset).
    # Onward-spread shape (dispersal) was demoted to characterisation (see above),
    # so it does not contribute to scoring, ranking or tiering.
    #   - rank on the row-wise max calibrated z across the two axes;
    #   - assign a primary tier by which axis (if any) is significant;
    #   - record the set of significant axes as an "axes_fired" signature string
    #     for the archetype classifier.
    axis_z_cols = ["burst_score_null_z", "burden_score_null_z"]
    tested["max_axis_null_z"] = tested[axis_z_cols].max(axis=1)
    tested["candidate_rank"] = tested["max_axis_null_z"].rank(
        ascending=False,
        method="min",
        na_option="bottom",
    )

    enough_sequences = tested["cluster_size"].ge(min_cluster_size)
    burden_eligible = tested["burden_eligible"].fillna(False)

    # Per-axis significance (size-gated). burden also requires its subpopulation
    # eligibility — a node with no onward burden cannot be significant on an axis
    # that is undefined for it.
    burst_sig = tested["burst_score_upper_p"].le(0.05) & enough_sequences
    burden_sig = (
        tested["burden_score_upper_p"].le(0.05)
        & enough_sequences
        & burden_eligible
    )
    any_sig = burst_sig | burden_sig

    # axes_fired signature: which axes are significant, joined (e.g. "burst+burden").
    def _fired(row_burst, row_burden):
        fired = []
        if row_burst:
            fired.append("burst")
        if row_burden:
            fired.append("burden")
        return "+".join(fired) if fired else "none"

    tested["axes_fired"] = [
        _fired(b, u) for b, u in zip(burst_sig, burden_sig)
    ]

    # Weaker (possible-review) band on either axis, not already high-priority.
    burst_possible = tested["burst_score_upper_p"].le(0.10)
    burden_possible = tested["burden_score_upper_p"].le(0.10) & burden_eligible
    any_possible = (
        (burst_possible | burden_possible) & enough_sequences & ~any_sig
    )

    # Uncalibrated fallback: high observed burst or burden but no null p (stratum
    # too small to calibrate).
    uncalibrated = (
        (
            (tested["burst_score"].ge(0.90) & tested["burst_score_upper_p"].isna())
            | (tested["burden_score"].ge(0.90) & tested["burden_score_upper_p"].isna())
        )
        & enough_sequences
        & ~any_sig
    )

    # Primary tier: both-axis first, then single-axis, then weaker / fallback.
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

    # Record burden applicability explicitly for downstream interpretation.
    # Distinguishes "burden tested but not significant" from "no onward burden".
    tested["burden_status"] = np.where(
        burden_eligible, "burden_tested", "burden_not_applicable"
    )

    # ------------------------------------------------------------------
    # 7. Characterisation features (NOT used in detection).
    # ------------------------------------------------------------------
    # These describe detected candidates and feed the archetype classifier; they
    # play no part in the detection axes or tiering.
    #   - burst_vs_burden_contrast: signed contained-burst <-> onward-spreader
    #     axis on the propagating subset. Among clusters with onward burden, burst
    #     and burden oppose (r ~ -0.21): large-for-context clusters tend to be the
    #     more CONTAINED ones, while the real onward drivers are modest in size.
    #     Positive = burst-dominant (contained outbreak); negative = burden-dominant
    #     (spreader). NaN where no onward burden is defined.
    #   - onward_entropy_norm, dominant_successor_frac: onward-spread SHAPE
    #     (evenness; concentration in top successor). Demoted from a detection
    #     axis (out_degree==2 dominated, near-binary, poorly calibrated) and kept
    #     raw to distinguish multiplicative fan-out from concentrated hand-off.
    #   - novelty_fraction: fresh-emergence vs inherited continuation.
    #   - onward_retention_ratio: per-capita onward weight (descriptive).
    #   - coherence_score + the entropy-z columns: geographic/socio-demographic
    #     concentration.

    # Signed contained <-> spreader contrast on the burden-eligible subpopulation.
    # NaN where no onward burden is defined (a non-propagating cluster has no
    # contrast).
    tested["burst_vs_burden_contrast"] = np.where(
        burden_eligible,
        tested["burst_score_null_z"] - tested["burden_score_null_z"],
        np.nan,
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
    cluster_table = build_cluster_table(
        cluster_att,
        cluster_stats,
        edge_table,
        sequence_df=df,
    )

    print("...........detecting SSEs")
    sse_df = add_sse_node_metrics(cluster_table)

    out_path = PROJECT_ROOT / "sse_detection/results/sse_outputs"
    out_path.mkdir(exist_ok=True)
    sse_df.to_parquet(out_path / "cluster_table.parquet")
    edge_table.to_parquet(out_path / "edge_table.parquet")

    print("Done.")

if __name__ == "__main__":
    main()
