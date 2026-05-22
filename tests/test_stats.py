import pandas as pd

from sse_detection.lib.stats import (
    add_sse_node_metrics,
    categorise_sse_nodes,
    flag_sse,
)


def test_core_amplification_is_ranked_within_window():
    df = pd.DataFrame(
        {
            "window_idx": [1, 1, 2, 2],
            "core_amplification_score": [0.10, 0.20, 0.80, 0.90],
            "core_amplification_score_pct_window": [0.5, 1.0, 0.5, 1.0],
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
            "core_amplification_score": [0.90, 0.95],
            "core_amplification_score_pct_window": [0.5, 1.0],
            "cluster_size": [6, 6],
            "in_degree": [1, 0],
            "out_degree": [0, 0],
            "in_strength": [3, 0],
            "out_strength": [0, 0],
            "novelty_fraction": [0, 0],
            "death": [True, True],
            "isolated": [False, True],
            "birth_like": [False, True],
            "continuation": [False, False],
            "merging": [False, False],
        }
    )

    out = categorise_sse_nodes(df, high_q=0.5, min_cluster_size=6)

    assert out.loc[0, "sse_candidate"]
    assert out.loc[0, "sse_onward_dynamic"] == "contained_burst"
    assert out.loc[1, "sse_candidate"]
    assert out.loc[1, "sse_onward_dynamic"] == "no_observed_onward_spread"


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
