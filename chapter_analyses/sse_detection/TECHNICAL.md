# SSE Detector: Technical Reference

## Entry point and configuration

Run:

```bash
python -m chapter_analyses.sse_detection.lib.sse.detection
```

The command loads sequence-window rows, builds the transition graph, assembles node features, scores/calibrates two axes, assigns tiers, validates regression columns, and writes detector plus graph-summary outputs.

It has no arguments. Settings in `lib/sse/config.py` are:

| Constant | Value | Use |
| --- | ---: | --- |
| `TRANSITION_WINDOW_STRIDE` | 2 | Retain alternate sorted source windows |
| `MIN_CLUSTER_SIZE` | 6 | Candidate-testing floor |
| `N_ENTROPY_DRAWS` | 1,000 | Composition-entropy null draws |
| `N_PERMUTATIONS` | 1,000 | Axis-profile permutations |
| `DETECTION_RANDOM_SEED` | 42 | Entropy and calibration RNG |

`load_sequence_data()` calls `utils.load_analysis_columns(...)` with policy variables and this stride. The shared loader retains every second sorted source window, then renumbers retained `window_idx` and `window_id` values consecutively; `cluster_id` remains the source cluster identifier.

## Input and transition graph

Input rows are sequence-window memberships. Transition construction requires `sequence_id`, `cluster_id`, `window_id`, and `window_idx`; configured context fields provide dates, lineage, cluster size, policy, surveillance, and socio-geodemographic attributes.

A node is one `cluster_id`. An edge `A -> B` exists when different clusters share at least one unique `sequence_id` and `B.window_idx = A.window_idx + 1` after stride filtering/renumbering. Its weight is `n_shared_sequences`.

The graph should be acyclic because edges move forward one retained window. Cumulative-burden construction raises if a cycle is present.

Per-node flow/shape fields include:

- `in_degree`, `out_degree`, `in_strength`, and `out_strength`;
- outgoing Shannon entropy, normalised entropy, effective successors, and dominant-successor fraction;
- weak-component ID;
- source/sink, branching/merging, primary-role, and observation-boundary flags.

Shape metrics are descriptive and do not enter candidate scores.

## Cluster features

`cluster_features.py` reduces the input to one row per cluster and joins:

- window, policy, lineage, size, duration, dates, and surveillance context;
- observed and null-standardised entropy for sex, age, SIMD, Data Zone, local authority, Health Board, and urban/rural composition;
- transition flow/shape, upstream novelty, and downstream burdens.

Socio-demographic and geographic fields are excluded from scoring.

### Upstream novelty

For target sequence set (S_i), direct-parent memberships are unioned as (P_i):

```text
unique_upstream_sequences = |S_i intersect P_i|
unique_local_new_sequences = |S_i minus P_i|
unique_local_new_sequences_ratio = |S_i minus P_i| / |S_i|
```

These values exist only when at least one parent is observed. Parentless nodes have `upstream_novelty_eligible = False` and missing novelty values.

### Downstream burden

Descriptive immediate fields count successor cluster size, successor sequences absent from the source, and a supported variant restricted to edges with at least 2 shared sequences.

For each source-target edge:

```text
attribution weight = edge shared sequences / target total incoming strength
source contribution = target sequences absent from source * attribution weight
```

Summing contributions gives fractional `source_attributable_new_downstream_burden`. It is graph credit, not causal attribution.

`cumulative_unique_new_sequences` unions all graph-reachable descendant memberships, removes source members, and counts each remaining identifier once.

## Detection axes

Only clusters with `cluster_size >= 6` are tested.

### Burst

The detector transforms size with `log1p`. It starts with the global tested mean/SD, then replaces them where groups have at least 20 observations and positive SD, first by `window_idx`, then by `window_idx + clade`:

```text
sampling_adjusted_excess_size =
  (log1p(cluster_size) - contextual mean) / contextual SD
```

Within each window it percentile-ranks:

```text
sampling_adjusted_excess_size
unique_local_new_sequences_ratio
```

`burst_score` is their available-component mean. Parentless nodes therefore have a one-component score. Raw log size and `log_excess_over_upstream` remain descriptive.

### Burden

The two components are:

```text
log1p(source_attributable_new_downstream_burden) - log1p(cluster_size)
log1p(cumulative_unique_new_sequences) - log1p(cluster_size)
```

A node is burden-eligible if either underlying burden is positive. Components are percentile-ranked within window among burden-eligible nodes, and `burden_score` is their mean. Other nodes retain missing burden scores.

## Permutation calibration

`add_composite_null_scores(...)` permutes raw component profiles within strata and recomputes within-window percentile composites. The production `null_mode="profile"` preserves dependence and component missingness; `"independent"` is available only through direct function use.

`choose_permutation_strata(...)` always uses its first available stratum field. It adds each subsequent nested field only for rows whose resulting group has at least 20 observations. Thus burst starts with `upstream_novelty_eligible`, then may refine by window and clade; burden uses clade as its sole stratum field. The 20-row floor does not collapse a sparse first field to an overall stratum.

For observed score (s_i) and (B) valid null scores:

```text
conservative p = (1 + count(null >= s_i)) / (1 + B)
randomized p   = (1 + count(null > s_i) + U_i * count(null == s_i)) / (1 + B)
```

The seeded randomised value is stored in `*_upper_p` and used operationally; both explicit versions are retained. Null mean, SD, and z-score are also written.

## Candidate tiers

| Tier | Rule |
| --- | --- |
| `size_ineligible` | cluster size below 6 |
| `high_priority_both_axes` | both applicable p-values `<= 0.05` |
| `high_priority_burst` | burst only `<= 0.05` |
| `high_priority_burden` | burden only `<= 0.05` |
| `possible_review` | no high-priority axis, but either applicable p-value `<= 0.10` |
| `high_score_uncalibrated` | score `>= 0.90` with its p-value missing |
| `background_or_low_information` | otherwise |

`axes_fired` records only axes significant at 0.05. `candidate_rank` ranks the maximum available axis null z-score. High-priority tiers, not `possible_review`, define the binary candidate outcome used downstream.

## Outputs

The detector writes to `results/sse_outputs/`:

| File | Unit |
| --- | --- |
| `cluster_table.parquet` | One scored row per cluster; authoritative detector/regression input |
| `edge_table.parquet` | One directed adjacent-retained-window edge |
| `transition_node_table.parquet` | Descriptive graph-role table |
| `transition_graph_summary.{csv,parquet}` | Whole graph |
| `transition_window_summary.{csv,parquet}` | Retained window |
| `transition_component_summary.{csv,parquet}` | Weak component |

Before writing, `validate_regression_cluster_columns(...)` checks the cluster table against the saved Bayesian specifications.

`build_composition_tables.py` separately writes `cluster_composition_<variable>.{csv,parquet}` for clusters meeting its size threshold. Level columns are proportions among non-missing sequence values; high-priority tiers map to `sse_status = candidate`.

## Verification and interpretation

Check adjacency, acyclicity, score missingness, component counts, calibration by size/window, tier sizes, burst-burden correlation on their joint domain, and right-censoring by source window. Regenerate figures, composition tables, and Bayesian results after detector changes.

Edges are shared membership, clusters depend on sampling and method parameters, source attribution is not causal, and absence of parents/descendants may reflect observation boundaries. Candidate status prioritises review; it does not establish a transmission chain or verified SSE.
