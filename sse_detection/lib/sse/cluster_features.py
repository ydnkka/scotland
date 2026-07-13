"""Cluster composition and node-feature assembly for SSE detection."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import DETECTION_RANDOM_SEED, N_ENTROPY_DRAWS
from .entropy import (
    add_mixing_tertiles,
    cluster_socio_demo_entropy,
)
from .transition_graph import (
    add_downstream_burden,
    add_edge_flow_metrics,
    add_new_downstream_metrics,
)


def build_cluster_stats(
    df: pd.DataFrame,
    *,
    n_entropy_draws: int = N_ENTROPY_DRAWS,
    random_state: int = DETECTION_RANDOM_SEED,
) -> pd.DataFrame:
    """Compute cluster-level composition and mixing-entropy stats."""
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
        "age_group",
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
    )

    return add_mixing_tertiles(cluster_stats)


def safe_mode(values: pd.Series):
    """Return the first modal non-missing value, or NaN when none exists."""
    values = values.dropna()
    if values.empty:
        return np.nan
    return values.mode().iloc[0]


def build_cluster_attributes(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse sequence-window rows to one row per cluster node."""
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
    return (
        df.groupby(["cluster_id", "window_id", "window_idx"], as_index=False)
        .agg(**agg)
        .sort_values(["window_idx", "cluster_size"], ascending=[True, False])
        .reset_index(drop=True)
    )


def build_cluster_table(
    cluster_att: pd.DataFrame,
    cluster_stats: pd.DataFrame,
    edge_table: pd.DataFrame,
    sequence_df: pd.DataFrame,
    *,
    min_supported_shared_sequences: int = 2,
) -> pd.DataFrame:
    """Join node attributes, composition stats, and transition-derived metrics."""
    cluster_table = cluster_att.merge(cluster_stats, on="cluster_id", how="left")
    cluster_table = add_edge_flow_metrics(cluster_table, edge_table)
    cluster_table = add_downstream_burden(cluster_table, edge_table)
    return add_new_downstream_metrics(
        cluster_table,
        sequence_df,
        edge_table,
        min_shared_sequences=min_supported_shared_sequences,
    )
