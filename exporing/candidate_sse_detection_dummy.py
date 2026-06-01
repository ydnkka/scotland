#!/usr/bin/env python3
"""
Candidate superspreading event (SSE) detection on a sparse temporal cluster-transition graph.

This script uses dummy data to demonstrate a robust workflow for ranking candidate
SSE-like cluster nodes in overlapping epidemic windows.

The workflow:
1. Simulate temporal genomic clusters with socio-geodemographic metadata.
2. Build a sparse transition graph where nodes are clusters and directed edges link
   adjacent windows with weights equal to shared sequence counts.
3. Compute node-level features:
   - cluster size
   - incoming/outgoing transition weight
   - weighted out-degree
   - outgoing entropy
   - downstream burden
   - temporal compactness
   - socio-geodemographic coherence
   - sampling-adjusted excess size
4. Build a stratified permutation null within window/lineage/region strata.
5. Compute empirical p-values and anomaly z-scores.
6. Rank nodes as candidate SSEs.

Requirements:
    pip install pandas numpy networkx matplotlib

Optional:
    pip install scikit-learn

Run:
    python candidate_sse_detection_dummy.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# -----------------------------
# Configuration
# -----------------------------

RANDOM_SEED = 42
N_WINDOWS = 8
N_CLUSTERS_PER_WINDOW = (10, 18)
N_PERMUTATIONS = 250
OUTPUT_DIR = Path("sse_dummy_outputs")


@dataclass
class SimulationConfig:
    seed: int = RANDOM_SEED
    n_windows: int = N_WINDOWS
    n_clusters_min: int = N_CLUSTERS_PER_WINDOW[0]
    n_clusters_max: int = N_CLUSTERS_PER_WINDOW[1]
    sse_nodes: Tuple[Tuple[int, int], ...] = ((2, 3), (4, 7))  # (window, within-window cluster index)


# -----------------------------
# Utilities
# -----------------------------

def entropy_from_counts(counts: Iterable[int | float]) -> float:
    """Shannon entropy for positive counts."""
    arr = np.asarray(list(counts), dtype=float)
    arr = arr[arr > 0]
    if arr.size == 0:
        return 0.0
    p = arr / arr.sum()
    return float(-(p * np.log(p)).sum())


def normalise_series(s: pd.Series) -> pd.Series:
    """Z-score with protection against zero variance."""
    sd = s.std(ddof=0)
    if sd == 0 or np.isnan(sd):
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.mean()) / sd


def empirical_upper_tail_p(obs: float, null_values: np.ndarray) -> float:
    """
    Upper-tail empirical p-value with +1 smoothing.
    Smaller values indicate observed score is unusually high.
    """
    return float((1 + np.sum(null_values >= obs)) / (len(null_values) + 1))


# -----------------------------
# Dummy data simulation
# -----------------------------

def simulate_dummy_clusters(config: SimulationConfig) -> pd.DataFrame:
    """
    Simulate cluster nodes in overlapping epidemic windows.

    Each row represents one cluster node. The columns imitate node attributes you
    might have in a real analysis: size, lineage, region, sampling intensity,
    age composition, deprivation profile, setting enrichment and date range.
    """
    rng = np.random.default_rng(config.seed)

    regions = ["Greater Glasgow", "Lothian", "Tayside", "Grampian", "Lanarkshire"]
    lineages = ["B.1.1.7", "AY.4", "BA.1", "BA.2"]
    settings = ["community", "workplace", "care_home", "university", "hospital"]

    rows = []
    node_id = 0

    for window in range(config.n_windows):
        n_clusters = int(rng.integers(config.n_clusters_min, config.n_clusters_max + 1))

        for idx in range(n_clusters):
            node = f"w{window}_c{idx}"

            region = rng.choice(regions, p=[0.28, 0.24, 0.16, 0.16, 0.16])
            lineage = rng.choice(lineages, p=[0.25, 0.30, 0.25, 0.20])
            setting = rng.choice(settings, p=[0.58, 0.14, 0.08, 0.12, 0.08])

            # Background sequencing intensity varies by region and time.
            region_multiplier = {
                "Greater Glasgow": 1.40,
                "Lothian": 1.25,
                "Tayside": 0.90,
                "Grampian": 0.80,
                "Lanarkshire": 0.95,
            }[region]
            time_multiplier = 0.75 + 0.12 * window + rng.normal(0, 0.04)
            sampling_intensity = max(0.05, region_multiplier * time_multiplier)

            # Most clusters are small; sparse genomic cluster graphs often have many singletons.
            size = int(1 + rng.negative_binomial(n=1.2, p=0.55))

            # Inject a few SSE-like nodes: compact, larger, coherent and onward-seeding.
            is_injected_sse = (window, idx) in config.sse_nodes
            if is_injected_sse:
                size = int(rng.integers(25, 45))
                setting = rng.choice(["workplace", "university", "care_home"])
                sampling_intensity *= rng.uniform(0.9, 1.2)

            temporal_span_days = int(max(1, rng.gamma(shape=2.0, scale=2.0)))
            if is_injected_sse:
                temporal_span_days = int(rng.integers(1, 4))

            # Coherence-like attributes. Larger means more concentrated within one group/setting.
            location_homogeneity = np.clip(rng.beta(5, 2), 0, 1)
            age_homogeneity = np.clip(rng.beta(3, 3), 0, 1)
            deprivation_homogeneity = np.clip(rng.beta(2.5, 3.5), 0, 1)
            setting_enrichment = np.clip(rng.beta(2, 5), 0, 1)

            if is_injected_sse:
                location_homogeneity = np.clip(rng.normal(0.92, 0.04), 0, 1)
                age_homogeneity = np.clip(rng.normal(0.78, 0.08), 0, 1)
                deprivation_homogeneity = np.clip(rng.normal(0.65, 0.10), 0, 1)
                setting_enrichment = np.clip(rng.normal(0.90, 0.05), 0, 1)

            rows.append(
                {
                    "node": node,
                    "window": window,
                    "cluster_index": idx,
                    "node_number": node_id,
                    "cluster_size": size,
                    "region": region,
                    "lineage": lineage,
                    "dominant_setting": setting,
                    "sampling_intensity": sampling_intensity,
                    "temporal_span_days": temporal_span_days,
                    "location_homogeneity": location_homogeneity,
                    "age_homogeneity": age_homogeneity,
                    "deprivation_homogeneity": deprivation_homogeneity,
                    "setting_enrichment": setting_enrichment,
                    "injected_sse": is_injected_sse,
                }
            )
            node_id += 1

    nodes = pd.DataFrame(rows)

    # Sampling-adjusted expected size. In real work this could come from a GLM,
    # region-time sequencing denominators, case counts or a Bayesian observation model.
    nodes["expected_size_sampling"] = (
        nodes.groupby(["window", "lineage", "region"])["cluster_size"]
        .transform("mean")
        .fillna(nodes["cluster_size"].mean())
    )

    # Smooth the expected value towards global sampling intensity.
    global_mean = nodes["cluster_size"].mean()
    nodes["expected_size_sampling"] = (
        0.5 * nodes["expected_size_sampling"]
        + 0.5 * global_mean * nodes["sampling_intensity"] / nodes["sampling_intensity"].mean()
    ).clip(lower=0.5)

    return nodes


def simulate_transition_edges(nodes: pd.DataFrame, config: SimulationConfig) -> pd.DataFrame:
    """
    Simulate directed weighted edges between adjacent epidemic windows.

    Edge weights represent the number of shared sequences between cluster nodes in
    adjacent overlapping windows. The graph is intentionally sparse.
    """
    rng = np.random.default_rng(config.seed + 1)

    edges = []

    for window in range(config.n_windows - 1):
        current = nodes[nodes["window"] == window].copy()
        nxt = nodes[nodes["window"] == window + 1].copy()

        for _, src in current.iterrows():
            # Many nodes are isolated or terminal.
            base_prob = 0.08
            size_effect = min(0.25, src["cluster_size"] / 150)
            p_has_edge = base_prob + size_effect

            if src["injected_sse"]:
                p_has_edge = 0.95

            if rng.random() > p_has_edge:
                continue

            if src["injected_sse"]:
                n_targets = int(rng.integers(3, 7))
            else:
                n_targets = int(rng.choice([1, 1, 1, 2, 3], p=[0.55, 0.20, 0.10, 0.10, 0.05]))

            # Prefer targets sharing lineage or region, but allow spillover.
            target_scores = np.ones(len(nxt), dtype=float)
            target_scores += (nxt["lineage"].values == src["lineage"]) * 3.0
            target_scores += (nxt["region"].values == src["region"]) * 2.0
            target_scores += nxt["cluster_size"].values / max(1, nxt["cluster_size"].max())
            target_scores = target_scores / target_scores.sum()

            target_indices = rng.choice(
                np.arange(len(nxt)),
                size=min(n_targets, len(nxt)),
                replace=False,
                p=target_scores,
            )

            for target_idx in target_indices:
                tgt = nxt.iloc[target_idx]

                if src["injected_sse"]:
                    weight = int(rng.integers(3, min(15, src["cluster_size"]) + 1))
                else:
                    weight = int(rng.choice([1, 1, 1, 2, 3, 4], p=[0.45, 0.25, 0.12, 0.10, 0.05, 0.03]))

                # Shared sequences cannot exceed source or target cluster sizes.
                weight = max(1, min(weight, int(src["cluster_size"]), int(tgt["cluster_size"])))

                edges.append(
                    {
                        "source": src["node"],
                        "target": tgt["node"],
                        "source_window": int(src["window"]),
                        "target_window": int(tgt["window"]),
                        "weight": weight,
                    }
                )

    return pd.DataFrame(edges)


# -----------------------------
# Graph and feature engineering
# -----------------------------

def build_transition_graph(nodes: pd.DataFrame, edges: pd.DataFrame) -> nx.DiGraph:
    """Build directed weighted transition graph."""
    G = nx.DiGraph()

    for _, row in nodes.iterrows():
        G.add_node(row["node"], **row.to_dict())

    for _, row in edges.iterrows():
        G.add_edge(row["source"], row["target"], weight=float(row["weight"]))

    return G


def compute_node_features(nodes: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    """
    Compute node-level features for candidate SSE detection.
    """
    features = nodes.copy()

    if edges.empty:
        features["incoming_weight"] = 0
        features["outgoing_weight"] = 0
        features["weighted_outdegree"] = 0
        features["weighted_indegree"] = 0
        features["outgoing_entropy"] = 0.0
        features["downstream_burden"] = 0
    else:
        outgoing = edges.groupby("source").agg(
            outgoing_weight=("weight", "sum"),
            weighted_outdegree=("target", "nunique"),
            outgoing_entropy=("weight", entropy_from_counts),
        )
        incoming = edges.groupby("target").agg(
            incoming_weight=("weight", "sum"),
            weighted_indegree=("source", "nunique"),
        )

        features = features.merge(outgoing, left_on="node", right_index=True, how="left")
        features = features.merge(incoming, left_on="node", right_index=True, how="left")

        target_sizes = features[["node", "cluster_size"]].rename(
            columns={"node": "target", "cluster_size": "target_cluster_size"}
        )
        downstream = (
            edges.merge(target_sizes, on="target", how="left")
            .groupby("source")
            .agg(downstream_burden=("target_cluster_size", "sum"))
        )
        features = features.merge(downstream, left_on="node", right_index=True, how="left")

    for col in [
        "incoming_weight",
        "outgoing_weight",
        "weighted_outdegree",
        "weighted_indegree",
        "outgoing_entropy",
        "downstream_burden",
    ]:
        features[col] = features[col].fillna(0)

    features["sampling_adjusted_excess_size"] = (
        features["cluster_size"] - features["expected_size_sampling"]
    ) / np.sqrt(features["expected_size_sampling"].clip(lower=0.5))

    # Compactness: shorter span is more SSE-like, but use a bounded transform.
    features["temporal_compactness"] = 1 / (1 + features["temporal_span_days"])

    features["epi_coherence"] = (
        0.30 * features["location_homogeneity"]
        + 0.20 * features["age_homogeneity"]
        + 0.15 * features["deprivation_homogeneity"]
        + 0.25 * features["setting_enrichment"]
        + 0.10 * features["temporal_compactness"]
    )

    # Raw composite is interpretable but not yet null-calibrated.
    score_components = [
        "sampling_adjusted_excess_size",
        "outgoing_weight",
        "weighted_outdegree",
        "outgoing_entropy",
        "downstream_burden",
        "epi_coherence",
    ]

    for col in score_components:
        features[f"z_{col}"] = normalise_series(features[col])

    features["raw_sse_score"] = (
        1.25 * features["z_sampling_adjusted_excess_size"]
        + 1.10 * features["z_outgoing_weight"]
        + 0.85 * features["z_weighted_outdegree"]
        + 0.50 * features["z_outgoing_entropy"]
        + 1.00 * features["z_downstream_burden"]
        + 0.80 * features["z_epi_coherence"]
    )

    return features


# -----------------------------
# Null model
# -----------------------------

def permute_features_within_strata(
    features: pd.DataFrame,
    rng: np.random.Generator,
    strata_cols: Tuple[str, ...] = ("window", "lineage", "region"),
) -> pd.DataFrame:
    """
    Permutation null.

    This keeps the graph-derived node positions fixed, but permutes cluster size
    and epidemiological coherence-related attributes within window/lineage/region
    strata where possible. For small strata, falls back to window/lineage, then window.

    This answers:
        Is this node unusually SSE-like compared with other nodes sampled in a similar
        temporal, lineage and regional context?
    """
    perm = features.copy()

    cols_to_permute = [
        "cluster_size",
        "expected_size_sampling",
        "sampling_adjusted_excess_size",
        "temporal_span_days",
        "temporal_compactness",
        "location_homogeneity",
        "age_homogeneity",
        "deprivation_homogeneity",
        "setting_enrichment",
        "epi_coherence",
    ]

    # Use adaptive strata so we do not create many singleton strata.
    def choose_strata(df: pd.DataFrame) -> pd.Series:
        key3 = df[list(strata_cols)].astype(str).agg("|".join, axis=1)
        counts3 = key3.map(key3.value_counts())
        key2 = df[["window", "lineage"]].astype(str).agg("|".join, axis=1)
        counts2 = key2.map(key2.value_counts())
        key1 = df[["window"]].astype(str).agg("|".join, axis=1)
        return np.where(counts3 >= 4, key3, np.where(counts2 >= 4, key2, key1))

    perm["_perm_stratum"] = choose_strata(perm)

    for _, idx in perm.groupby("_perm_stratum").groups.items():
        idx = list(idx)
        if len(idx) <= 1:
            continue
        for col in cols_to_permute:
            perm.loc[idx, col] = rng.permutation(perm.loc[idx, col].values)

    # Recompute z columns and score using same formula.
    score_components = [
        "sampling_adjusted_excess_size",
        "outgoing_weight",
        "weighted_outdegree",
        "outgoing_entropy",
        "downstream_burden",
        "epi_coherence",
    ]
    for col in score_components:
        perm[f"z_{col}"] = normalise_series(perm[col])

    perm["raw_sse_score"] = (
        1.25 * perm["z_sampling_adjusted_excess_size"]
        + 1.10 * perm["z_outgoing_weight"]
        + 0.85 * perm["z_weighted_outdegree"]
        + 0.50 * perm["z_outgoing_entropy"]
        + 1.00 * perm["z_downstream_burden"]
        + 0.80 * perm["z_epi_coherence"]
    )

    return perm.drop(columns=["_perm_stratum"])


def add_permutation_p_values(
    features: pd.DataFrame,
    n_permutations: int = N_PERMUTATIONS,
    seed: int = RANDOM_SEED + 2,
) -> pd.DataFrame:
    """
    Compute empirical p-values for each node's composite SSE score.
    """
    rng = np.random.default_rng(seed)

    null_scores = np.zeros((len(features), n_permutations))
    for b in range(n_permutations):
        perm = permute_features_within_strata(features, rng)
        null_scores[:, b] = perm["raw_sse_score"].values

    obs = features["raw_sse_score"].values
    pvals = np.array([
        empirical_upper_tail_p(obs[i], null_scores[i, :])
        for i in range(len(features))
    ])

    out = features.copy()
    out["empirical_p"] = pvals
    out["null_mean_score"] = null_scores.mean(axis=1)
    out["null_sd_score"] = null_scores.std(axis=1, ddof=0)
    out["null_adjusted_score"] = (
        out["raw_sse_score"] - out["null_mean_score"]
    ) / out["null_sd_score"].replace(0, np.nan)
    out["null_adjusted_score"] = out["null_adjusted_score"].fillna(0)

    # A pragmatic classification, not a definitive causal claim.
    out["candidate_class"] = "unlikely_or_uninformative"

    large_or_connected = (
        (out["cluster_size"] >= 3)
        | (out["outgoing_weight"] >= 2)
        | (out["downstream_burden"] >= 5)
    )

    high_conf = (
        large_or_connected
        & (out["empirical_p"] <= 0.05)
        & (out["sampling_adjusted_excess_size"] > 1.0)
        & ((out["outgoing_weight"] > 0) | (out["downstream_burden"] > 0))
        & (out["epi_coherence"] >= out["epi_coherence"].quantile(0.60))
    )

    possible = (
        large_or_connected
        & ~high_conf
        & (
            (out["empirical_p"] <= 0.10)
            | (
                (out["sampling_adjusted_excess_size"] > 1.5)
                & (out["epi_coherence"] >= out["epi_coherence"].quantile(0.50))
            )
            | (
                (out["outgoing_weight"] >= out["outgoing_weight"].quantile(0.90))
                & (out["downstream_burden"] >= out["downstream_burden"].quantile(0.80))
            )
        )
    )

    terminal_large = (
        (out["cluster_size"] >= out["cluster_size"].quantile(0.90))
        & (out["outgoing_weight"] == 0)
    )

    out.loc[possible, "candidate_class"] = "possible_candidate_sse"
    out.loc[terminal_large, "candidate_class"] = "large_terminal_cluster"
    out.loc[high_conf, "candidate_class"] = "high_confidence_candidate_sse"

    return out.sort_values(
        ["candidate_class", "null_adjusted_score", "raw_sse_score"],
        ascending=[True, False, False],
    )


# -----------------------------
# Plotting
# -----------------------------

def plot_graph(G: nx.DiGraph, ranked: pd.DataFrame, output_dir: Path) -> None:
    """
    Plot sparse temporal transition graph.
    """
    output_dir.mkdir(exist_ok=True, parents=True)

    # Layout by epidemic window on x-axis; within-window index on y-axis.
    pos = {}
    for node, data in G.nodes(data=True):
        pos[node] = (data["window"], -data["cluster_index"])

    classes = ranked.set_index("node")["candidate_class"].to_dict()
    sizes = ranked.set_index("node")["cluster_size"].to_dict()

    node_sizes = [
        40 + 18 * math.sqrt(sizes.get(n, 1))
        for n in G.nodes()
    ]

    # Avoid specifying colours manually; use default cycle through grouped draws.
    plt.figure(figsize=(13, 8))
    for cls in ranked["candidate_class"].unique():
        nodelist = [n for n in G.nodes() if classes.get(n) == cls]
        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=nodelist,
            node_size=[node_sizes[list(G.nodes()).index(n)] for n in nodelist],
            label=cls,
            alpha=0.85,
        )

    edge_widths = [0.4 + 0.5 * G[u][v]["weight"] for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.35, arrows=True, arrowsize=8)
    nx.draw_networkx_labels(
        G,
        pos,
        labels={n: n if classes.get(n) == "high_confidence_candidate_sse" else "" for n in G.nodes()},
        font_size=8,
    )

    plt.title("Dummy sparse temporal cluster-transition graph")
    plt.xlabel("Epidemic window")
    plt.ylabel("Cluster index within window")
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(output_dir / "transition_graph.png", dpi=180)
    plt.close()


def plot_ranked_scores(ranked: pd.DataFrame, output_dir: Path, top_n: int = 25) -> None:
    """
    Plot top-ranked node scores.
    """
    output_dir.mkdir(exist_ok=True, parents=True)

    top = ranked.sort_values("null_adjusted_score", ascending=False).head(top_n).copy()
    top = top.sort_values("null_adjusted_score", ascending=True)

    plt.figure(figsize=(10, 7))
    plt.barh(top["node"], top["null_adjusted_score"])
    plt.xlabel("Null-adjusted SSE score")
    plt.ylabel("Cluster node")
    plt.title(f"Top {top_n} candidate SSE-like nodes")
    plt.tight_layout()
    plt.savefig(output_dir / "ranked_candidate_scores.png", dpi=180)
    plt.close()


# -----------------------------
# Main
# -----------------------------

def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

    config = SimulationConfig()
    nodes = simulate_dummy_clusters(config)
    edges = simulate_transition_edges(nodes, config)
    G = build_transition_graph(nodes, edges)

    features = compute_node_features(nodes, edges)
    ranked = add_permutation_p_values(features, n_permutations=N_PERMUTATIONS)

    nodes.to_csv(OUTPUT_DIR / "dummy_cluster_nodes.csv", index=False)
    edges.to_csv(OUTPUT_DIR / "dummy_transition_edges.csv", index=False)
    features.to_csv(OUTPUT_DIR / "node_features.csv", index=False)
    ranked.to_csv(OUTPUT_DIR / "ranked_candidate_sse_nodes.csv", index=False)

    plot_graph(G, ranked, OUTPUT_DIR)
    plot_ranked_scores(ranked, OUTPUT_DIR)

    cols = [
        "node",
        "window",
        "cluster_size",
        "region",
        "lineage",
        "dominant_setting",
        "outgoing_weight",
        "weighted_outdegree",
        "downstream_burden",
        "sampling_adjusted_excess_size",
        "epi_coherence",
        "raw_sse_score",
        "null_adjusted_score",
        "empirical_p",
        "candidate_class",
        "injected_sse",
    ]

    print("\nTop candidate SSE-like nodes\n" + "-" * 34)
    print(
        ranked.sort_values("null_adjusted_score", ascending=False)[cols]
        .head(15)
        .to_string(index=False)
    )

    print(f"\nOutputs written to: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
