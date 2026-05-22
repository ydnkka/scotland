import pandas as pd

from sse_detection.lib.stats import (
    add_sse_node_metrics,
    categorise_sse_nodes,
    flag_sse,
)
from sse_detection.lib.palettes import (
    SSE_CATEGORY_ORDER,
    sse_category_palette_from,
)


def _candidate_frame(rows):
    base = {
        "window_idx": 1,
        "local_amplification_score_pct_window": 1.0,
        "cluster_size": 6,
        "in_degree": 0,
        "out_degree": 0,
        "in_strength": 0,
        "out_strength": 0,
        "novelty_fraction_pct_window": 0,
        "log_excess_over_upstream_pct_window": 0,
        "downstream_expansion_proxy_pct_window": 0,
        "downstream_expansion_proxy_pct_onward_window": 0.5,
        "out_strength_pct_window": 0,
        "downstream_entropy_norm": 0,
        "dominant_successor_frac": 0,
        "mixing_score": 0,
        "death": False,
        "isolated": False,
        "birth": False,
        "continuation": False,
        "merging": False,
    }
    return pd.DataFrame([{**base, **row} for row in rows])


def test_core_amplification_is_ranked_within_window():
    df = pd.DataFrame(
        {
            "window_idx": [1, 1, 2, 2],
            "local_amplification_score": [0.10, 0.20, 0.80, 0.90],
            "local_amplification_score_pct_window": [0.5, 1.0, 0.5, 1.0],
            "cluster_size": [6, 6, 6, 6],
            "in_degree": [0, 0, 0, 0],
            "out_degree": [0, 0, 0, 0],
            "in_strength": [0, 0, 0, 0],
            "out_strength": [0, 0, 0, 0],
            "novelty_fraction": [0, 0, 0, 0],
        }
    )

    out = categorise_sse_nodes(df, high_q=0.75, min_cluster_size=6)

    assert out.loc[1, "sse_candidate"]
    assert out.loc[3, "sse_candidate"]
    assert not out.loc[0, "sse_candidate"]
    assert not out.loc[2, "sse_candidate"]
    assert out.loc[0, "sse_category"] == "not_sse_like"
    assert out.loc[0, "sse_graph_category"] == "not_sse_like"


def test_onward_expansion_percentile_uses_onward_nodes_only():
    df = pd.DataFrame(
        {
            "window_idx": [1, 1, 1],
            "cluster_size": [10, 10, 10],
            "in_degree": [1, 1, 1],
            "out_degree": [0, 1, 1],
            "in_strength": [5, 5, 5],
            "out_strength": [0, 2, 5],
        }
    )

    out = add_sse_node_metrics(df)

    assert pd.isna(out.loc[0, "downstream_expansion_proxy_pct_onward_window"])
    assert out.loc[1, "downstream_expansion_proxy_pct_onward_window"] == 0.5
    assert out.loc[2, "downstream_expansion_proxy_pct_onward_window"] == 1.0


def test_contained_burst_is_reachable_for_terminal_candidate_with_burden():
    df = pd.DataFrame(
        {
            "window_idx": [1, 1],
            "local_amplification_score": [0.90, 0.95],
            "local_amplification_score_pct_window": [0.5, 1.0],
            "cluster_size": [6, 6],
            "in_degree": [1, 0],
            "out_degree": [0, 0],
            "in_strength": [3, 0],
            "out_strength": [0, 0],
            "novelty_fraction": [0, 0],
            "death": [True, True],
            "isolated": [False, True],
            "birth": [False, True],
            "continuation": [False, False],
            "merging": [False, False],
        }
    )

    out = categorise_sse_nodes(df, high_q=0.5, min_cluster_size=6)

    assert out.loc[0, "sse_candidate"]
    assert out.loc[0, "sse_onward_dynamic"] == "contained_burst"
    assert out.loc[0, "sse_graph_category"] == "terminal_sink__contained_burst"
    assert out.loc[0, "sse_category"] == "contained_local_burst"
    assert out.loc[1, "sse_candidate"]
    assert out.loc[1, "sse_onward_dynamic"] == "no_observed_onward_spread"
    assert out.loc[1, "sse_graph_category"] == "isolated_burst__no_observed_onward_spread"
    assert out.loc[1, "sse_category"] == "contained_local_burst"


def test_epidemiological_category_mapping_preserves_graph_category():
    rows = [
        {
            "expected_dynamic": "no_observed_onward_spread",
            "expected_category": "contained_local_burst",
        },
        {
            "in_strength": 2,
            "death": True,
            "expected_dynamic": "contained_burst",
            "expected_category": "contained_local_burst",
        },
        {
            "out_degree": 1,
            "out_strength": 2,
            "expected_dynamic": "single_successor_chain",
            "expected_category": "sustained_single_chain",
        },
        {
            "out_degree": 2,
            "out_strength": 2,
            "downstream_entropy_norm": 0.2,
            "expected_dynamic": "dominant_branch",
            "expected_category": "focused_branching_transmission",
        },
        {
            "out_degree": 2,
            "out_strength": 2,
            "downstream_entropy_norm": 0.8,
            "downstream_expansion_proxy_pct_window": 1.0,
            "expected_dynamic": "multi_branch_expander",
            "expected_category": "diffuse_branching_transmission",
        },
        {
            "out_degree": 2,
            "out_strength": 2,
            "downstream_entropy_norm": 0.8,
            "expected_dynamic": "multi_branch_seeder",
            "expected_category": "diffuse_branching_transmission",
        },
        {
            "out_degree": 2,
            "out_strength": 2,
            "downstream_entropy_norm": 0.8,
            "downstream_expansion_proxy_pct_window": 1.0,
            "mixing_score": 0.8,
            "expected_dynamic": "diverse_population_broadcaster",
            "expected_category": "mixed_population_dissemination",
        },
        {
            "out_degree": 2,
            "out_strength": 2,
            "downstream_entropy_norm": 0.5,
            "dominant_successor_frac": 0.5,
            "out_strength_pct_window": 1.0,
            "expected_dynamic": "high_volume_onward_spread",
            "expected_category": "high_volume_onward_transmission",
        },
        {
            "out_degree": 2,
            "out_strength": 2,
            "downstream_entropy_norm": 0.5,
            "dominant_successor_frac": 0.5,
            "expected_dynamic": "weak_or_ambiguous_onward_spread",
            "expected_category": "ambiguous_amplification_signal",
        },
    ]
    expected_dynamic = [row.pop("expected_dynamic") for row in rows]
    expected_category = [row.pop("expected_category") for row in rows]

    out = categorise_sse_nodes(_candidate_frame(rows), min_cluster_size=6)

    assert out["sse_onward_dynamic"].tolist() == expected_dynamic
    assert out["sse_category"].tolist() == expected_category
    assert out["sse_graph_category"].tolist() == [
        f"{role}__{dynamic}"
        for role, dynamic in zip(out["sse_role"], expected_dynamic)
    ]


def test_mixed_population_category_outranks_relay_status():
    df = _candidate_frame(
        [
            {
                "in_strength": 2,
                "out_strength": 4,
                "out_degree": 2,
                "continuation": True,
                "log_excess_over_upstream_pct_window": 1.0,
                "downstream_entropy_norm": 0.8,
                "downstream_expansion_proxy_pct_window": 1.0,
                "mixing_score": 0.8,
            }
        ]
    )

    out = categorise_sse_nodes(df, min_cluster_size=6)

    assert out.loc[0, "sse_role"] == "relay_amplifier"
    assert out.loc[0, "sse_onward_dynamic"] == "diverse_population_broadcaster"
    assert out.loc[0, "sse_category"] == "mixed_population_dissemination"
    assert out.loc[0, "sse_graph_category"] == (
        "relay_amplifier__diverse_population_broadcaster"
    )


def test_strict_birth_maps_to_putative_introduction_category():
    df = _candidate_frame(
        [
            {
                "birth": True,
                "out_degree": 1,
                "out_strength": 2,
                "novelty_fraction_pct_window": 1.0,
            }
        ]
    )

    out = categorise_sse_nodes(df, min_cluster_size=6)

    assert out.loc[0, "sse_role"] == "putative_birth"
    assert out.loc[0, "sse_onward_dynamic"] == "single_successor_chain"
    assert out.loc[0, "sse_category"] == "putative_introduction_burst"
    assert out.loc[0, "sse_graph_category"] == "putative_birth__single_successor_chain"


def test_relay_roles_map_to_secondary_relay_category():
    df = _candidate_frame(
        [
            {
                "in_strength": 2,
                "out_strength": 4,
                "out_degree": 1,
                "continuation": True,
                "log_excess_over_upstream_pct_window": 1.0,
            },
            {
                "in_degree": 2,
                "in_strength": 3,
                "out_strength": 4,
                "out_degree": 1,
                "merging": True,
            },
        ]
    )

    out = categorise_sse_nodes(df, min_cluster_size=6)

    assert out["sse_role"].tolist() == ["relay_amplifier", "merged_relay"]
    assert out["sse_category"].tolist() == [
        "secondary_relay_amplification",
        "secondary_relay_amplification",
    ]


def test_compact_sse_category_palette_covers_ordered_categories():
    palette = sse_category_palette_from(SSE_CATEGORY_ORDER)

    assert list(palette) == SSE_CATEGORY_ORDER
    assert set(palette) == set(SSE_CATEGORY_ORDER)
    assert palette["not_sse_like"] != palette["mixed_population_dissemination"]


def test_flag_sse_splits_first_observed_from_growth_and_skips_pre_observation_rows():
    df = pd.DataFrame(
        {
            "sequence_id": [f"s{i}" for i in range(15)],
            "collection_date": pd.to_datetime(
                [
                    "2024-01-08",
                    "2024-01-08",
                    "2024-01-08",
                    "2024-01-15",
                    "2024-01-15",
                    "2024-01-15",
                    "2024-01-15",
                    "2024-01-15",
                    "2024-01-15",
                    "2024-01-15",
                    "2024-01-15",
                    "2024-01-15",
                    "2024-01-15",
                    "2024-01-01",
                    "2024-01-15",
                ]
            ),
            "meta_cluster_id": ["A"] * 13 + ["B", "B"],
            "clade": ["X"] * 15,
        }
    )

    out = flag_sse(df, threshold=2, drop_incomplete_last_week=False)
    a_rows = out.loc[out["meta_cluster_id"].eq("A")].sort_values("week")

    assert a_rows["week"].min() == pd.Timestamp("2024-01-08")
    assert a_rows.iloc[0]["first_observed_burst"]
    assert not a_rows.iloc[0]["growth_burst"]
    assert a_rows.iloc[1]["growth_burst"]
    assert not a_rows.iloc[1]["first_observed_burst"]
    assert set(out.loc[out["is_sse"], "sse_flag_type"]) == {
        "first_observed_burst",
        "growth_burst",
    }
