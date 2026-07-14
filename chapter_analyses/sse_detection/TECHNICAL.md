# SSE Detector: Technical Reference <!-- omit in toc -->

This document is the implementation contract for the Chapter 5 SSE detector. It describes how sequence-window records become a directed transition graph, how burst and burden scores are constructed, how null calibration assigns candidate tiers, and which files the pipeline writes.

For methodological interpretation and limitations, see `DETECTION_RATIONALE.md`. For Bayesian characterisation after detection, see `BAYESIAN_MODELS.md`. For commands and the directory map, see `README.md`.

## 1. Content

<!-- TOC tocDepth:2..3 chapterDepth:2..6 -->

- [1. Content](#1-content)
- [2. Pipeline entry point](#2-pipeline-entry-point)
- [3. Configuration](#3-configuration)
- [4. Input contract](#4-input-contract)
- [5. Transition graph](#5-transition-graph)
  - [5.1. 1 Nodes and edges](#51-1-nodes-and-edges)
  - [5.2. 2 Basic flow and shape metrics](#52-2-basic-flow-and-shape-metrics)
- [6. Cluster feature assembly](#6-cluster-feature-assembly)
- [7. Upstream novelty](#7-upstream-novelty)
- [8. Burst axis](#8-burst-axis)
  - [8.1. 1 Context-adjusted excess size](#81-1-context-adjusted-excess-size)
  - [8.2. 2 Components](#82-2-components)
- [9. Downstream burden](#9-downstream-burden)
  - [9.1. 1 Descriptive immediate measures](#91-1-descriptive-immediate-measures)
  - [9.2. 2 Source-attributable immediate burden](#92-2-source-attributable-immediate-burden)
  - [9.3. 3 Cumulative unique future burden](#93-3-cumulative-unique-future-burden)
  - [9.4. 4 Components and applicability](#94-4-components-and-applicability)
- [10. Null calibration](#10-null-calibration)
- [11. Candidate tiers](#11-candidate-tiers)
- [12. Boundary and censoring fields](#12-boundary-and-censoring-fields)
- [13. Output contract](#13-output-contract)
- [14. Downstream workflow](#14-downstream-workflow)
- [15. Verification and sensitivity checks](#15-verification-and-sensitivity-checks)
- [16. Interpretation limits](#16-interpretation-limits)

<!-- /TOC -->

## 2. Pipeline entry point

Run from the Scotland repository root in the project environment:

```bash
python -m chapter_analyses.sse_detection.lib.sse.detection
```

The entry point in `lib/sse/detection.py` executes these stages in order:

1. Load the required sequence-window columns and policy variables.
2. Retain the configured transition-window stride.
3. Build directed adjacent-window transition edges.
4. Assemble cluster attributes, composition features, graph metrics, upstream
   novelty, and downstream burden metrics.
5. Calculate and calibrate the burst and burden axes.
6. Assign candidate tiers.
7. Validate the regression-facing cluster-table schema.
8. Write detector and descriptive transition-graph outputs.

The detector has no command-line parameter layer. Its reproducibility settings are constants in `lib/sse/config.py`; changing them requires regenerating all detector outputs and downstream results.

## 3. Configuration

Current detector constants are:

| Constant                   |  Value | Purpose                                                                  |
| -------------------------- | -----: | ------------------------------------------------------------------------ |
| `TRANSITION_WINDOW_STRIDE` |    `2` | Retain alternate source windows before constructing adjacent transitions |
| `MIN_CLUSTER_SIZE`         |    `6` | Minimum cluster size eligible for candidate testing                      |
| `N_ENTROPY_DRAWS`          | `1000` | Monte Carlo draws for composition-entropy nulls                          |
| `N_PERMUTATIONS`           | `1000` | Profile permutations for each detection axis                             |
| `DETECTION_RANDOM_SEED`    |   `42` | Reproducible entropy and detector calibration                            |

`add_sse_node_metrics(...)` additionally defaults to an adaptive permutation stratum floor of 20 observations.

## 4. Input contract

`load_sequence_data()` calls `utils.load_analysis_columns(...)` with the columns declared in `ANALYSIS_COLUMNS`, policy variables enabled, and the configured window stride. The data are sequence-window membership records, not one row per person or infection.

Fields required by transition construction are:

```text
sequence_id
cluster_id
window_id
window_idx
```

The remaining configured columns supply cluster size, dates, clade/lineage, epidemic context, and socio-geodemographic characterisation. `cluster_id` must identify one window-specific cluster consistently throughout the input.

## 5. Transition graph

### 5.1. 1 Nodes and edges

A node is a window-specific sequence cluster. A directed edge `A -> B` exists when:

- `B` occurs exactly one retained window after `A`;
- `A` and `B` are different clusters; and
- at least one sequence belongs to both clusters.

`n_shared_sequences` is the number of unique shared sequence identifiers and becomes the edge weight. Edges represent continuity between clustering solutions; they do not represent observed transmission.

Because every edge moves forward by one retained `window_idx`, the transition graph is expected to be a directed acyclic graph. Cumulative burden calculation fails explicitly if a cycle is present.

### 5.2. 2 Basic flow and shape metrics

For each node, the pipeline records:

| Metric                         | Definition                                                      |
| ------------------------------ | --------------------------------------------------------------- |
| `in_degree` / `out_degree`     | Number of direct parent/successor clusters                      |
| `in_strength` / `out_strength` | Sum of incoming/outgoing shared-sequence edge weights           |
| `onward_entropy`               | Shannon entropy of outgoing edge weights                        |
| `onward_entropy_norm`          | Entropy divided by its maximum for the observed successor count |
| `effective_successors`         | Exponentiated entropy (Hill number of order 1)                  |
| `dominant_successor_frac`      | Largest outgoing edge weight divided by total outgoing weight   |

Degree and entropy describe branching shape. They are retained for characterisation and are not detection-axis components.

## 6. Cluster feature assembly

`lib/sse/cluster_features.py` produces one row per cluster by joining:

- temporal, lineage, policy, cluster-size, and epidemic-context attributes;
- socio-demographic composition and entropy measures;
- incoming/outgoing graph flow;
- direct successor size summaries;
- deduplicated upstream novelty;
- immediate and cumulative downstream burden.

Socio-demographic variables are never used to assign candidates. They are carried forward for descriptive and Bayesian characterisation.

## 7. Upstream novelty

Upstream novelty is calculated only for nodes with at least one observed direct parent.

For target cluster `i`, let `S_i` be its sequence set and let `P_i` be the union of the sequence sets of all direct parents. Parent sets are unioned before comparison, preventing the same inherited sequence from being counted more than once.

```text
unique_upstream_sequences = |S_i intersect P_i|
unique_local_new_sequences = |S_i minus P_i|
unique_local_new_sequences_ratio = |S_i minus P_i| / |S_i|
```

`upstream_novelty_eligible` is true when a direct parent is observed. For parentless nodes, novelty counts and the ratio remain missing: absence of an observed parent is not treated as evidence that every sequence is locally new.

## 8. Burst axis

Only clusters with `cluster_size >= 6` are tested.

### 8.1. 1 Context-adjusted excess size

Cluster size is transformed as `log1p(cluster_size)`. Expected log size and its standard deviation are estimated adaptively, beginning with the global tested population and using sufficiently populated `window_idx` and then `window_idx + clade` groups. A group is used only when it has at least 20 rows and a positive finite standard deviation.

```text
sampling_adjusted_excess_size =
    (log1p(cluster_size) - contextual mean) / contextual SD
```

### 8.2. 2 Components

The burst composite contains:

```text
sampling_adjusted_excess_size
unique_local_new_sequences_ratio
```

Each component is percentile-ranked within `window_idx`; `burst_score` is the row-wise mean of the available percentile ranks. Parentless nodes therefore have a one-component burst score based on contextual magnitude. Nodes with observed parents have both components.

`log_cluster_size` and `log_excess_over_upstream` remain descriptive columns but are excluded from the composite. Raw size strongly duplicates contextual excess size, while summed incoming strength can count support from multiple parents and does not provide a deduplicated novelty measure.

## 9. Downstream burden

### 9.1. 1 Descriptive immediate measures

For a source sequence set `S` and its direct successors, the pipeline retains:

- `new_downstream_burden`: unique successor sequences absent from `S`;
- `supported_new_downstream_burden`: the same count after excluding edges with
  fewer than two shared sequences;
- `new_downstream_children` and `supported_new_downstream_children`;
- `mean_successor_new_sequences`;
- `downstream_cluster_burden` and `mean_successor_cluster_size`.

These measures are descriptive or sensitivity outputs and do not receive additional weight in the current burden composite.

### 9.2. 2 Source-attributable immediate burden

For each direct source-to-target edge, the target's new sequence count is multiplied by that edge's share of the target's total incoming strength:

```text
edge attribution weight = edge n_shared_sequences / target in_strength
source contribution = target new sequences * edge attribution weight
source_attributable_new_downstream_burden = sum(source contributions)
```

The measure can be fractional. It allocates direct downstream accumulation across competing parents according to graph support; it is not a causal transmission attribution.

### 9.3. 3 Cumulative unique future burden

For each source, all graph-reachable descendants are found in topological order. Their sequence sets are unioned, source sequences are removed, and each remaining identifier is counted once:

```text
cumulative_unique_new_sequences =
    |union(all descendant sequences) minus source sequences|
```

This measures observed longer-term reach. It is right-censored: nodes near the end of the observation period have fewer possible future windows.

### 9.4. 4 Components and applicability

The burden composite contains source-size-normalised log ratios:

```text
log_source_attributable_new_downstream_burden_ratio =
    log1p(source_attributable_new_downstream_burden) - log1p(cluster_size)

log_cumulative_unique_new_sequences_ratio =
    log1p(cumulative_unique_new_sequences) - log1p(cluster_size)
```

A node is burden-eligible when either underlying burden count is positive. Eligible values are percentile-ranked within window, and `burden_score` is the mean of the component percentiles. Non-propagating nodes retain missing burden scores rather than structural zero scores.

## 10. Null calibration

For each axis, the detector permutes intact multivariate component profiles. Profile permutation preserves empirical dependence and component-missingness patterns. An independent-component alternative, `null_mode="independent"`, exists for sensitivity work but is not used by the pipeline entry point.

Permutation labels are adaptive: more detailed contextual strata are retained only where they contain at least 20 observations. Burst strata begin with `upstream_novelty_eligible`, ensuring parentless one-component profiles are not permuted against two-component profiles, followed by window and clade context. Burden is calibrated among burden-eligible nodes, with clade used where the adaptive size criterion permits.

For observed score `s_i` and `B` valid null scores, the smoothed one-sided upper-tail probability is:

```text
p_i = (1 + number(null score >= s_i)) / (1 + B)
```

The associated null mean, null standard deviation, and null z-score are also written. The test asks whether a score is unusually high; low scores naturally produce p-values near one.

## 11. Candidate tiers

Candidate assignment is demographic-blind and uses axis-specific p-values:

| Tier                            | Rule                                                          |
| ------------------------------- | ------------------------------------------------------------- |
| `size_ineligible`               | `cluster_size < 6`                                            |
| `high_priority_both_axes`       | Burst and burden `p <= 0.05`                                  |
| `high_priority_burst`           | Burst only `p <= 0.05`                                        |
| `high_priority_burden`          | Burden only `p <= 0.05`                                       |
| `possible_review`               | Neither high priority, but either applicable axis `p <= 0.10` |
| `high_score_uncalibrated`       | Score at least 0.90 but its p-value is unavailable            |
| `background_or_low_information` | None of the above                                             |

`axes_fired` records `burst`, `burden`, `burst+burden`, or `none`. `candidate_rank` orders nodes by the maximum available axis null z-score.

## 12. Boundary and censoring fields

Graph-role fields distinguish isolated, source, continuing, branching, merging, and sink nodes. The node table also records observation-boundary censoring based on the first and last retained windows.

Interpretation must distinguish:

- parentless nodes, for which upstream novelty is not observable;
- last-window and sink nodes, for which future burden may be unobservable;
- genuine structural termination from termination caused by sequencing or
  clustering resolution.

The current candidate rule handles upstream missingness explicitly. Cumulative burden remains sensitive to unequal follow-up and must be checked by window in robustness analyses.

## 13. Output contract

The detector writes to `chapter_analyses/sse_detection/results/sse_outputs/`:

| File                                         | Granularity and use                                                   |
| -------------------------------------------- | --------------------------------------------------------------------- |
| `cluster_table.parquet`                      | One row per cluster; detector, characterisation, and regression input |
| `edge_table.parquet`                         | One row per directed adjacent-window transition                       |
| `transition_node_table.parquet`              | Descriptive node and graph-role table                                 |
| `transition_graph_summary.{csv,parquet}`     | Whole-graph summary                                                   |
| `transition_window_summary.{csv,parquet}`    | Window-level transition summary                                       |
| `transition_component_summary.{csv,parquet}` | Weak-component summary                                                |

The cluster table is the authoritative detector output. Important column families include:

- identifiers, dates, window, policy, lineage, and cluster size;
- socio-demographic composition and mixing entropy;
- degree, strength, branching shape, and graph role;
- upstream novelty and immediate/cumulative downstream burden;
- component percentiles, scores, null summaries, p-values, applicability;
- candidate tier, fired axes, contrast, and rank.

`validate_regression_cluster_columns(...)` runs before writing, so a detector run fails rather than silently producing a cluster table incompatible with the saved Bayesian specifications.

## 14. Downstream workflow

After regenerating detector outputs:

```bash
python -m chapter_analyses.sse_detection.make_figures --skip-missing
python -m chapter_analyses.sse_detection.build_composition_tables
```

Bayesian models must be refitted when candidate labels, burst scores, burden scores, eligibility, or regression-facing feature definitions change. See `BAYESIAN_MODELS.md` for dry runs, fitting commands, priors, diagnostics, and saved-model outputs.

## 15. Verification and sensitivity checks

Before reporting results:

1. Confirm edge window deltas are exactly one retained step.
2. Confirm the transition graph is acyclic.
3. Inspect missingness and component counts separately for both axes.
4. Inspect background upper-tail p-values overall, by size, and by component
   availability using `lib/sse/diagnostics.py`.
5. Check candidate counts and size distributions by tier.
6. Compare burst and burden correlations only where both are applicable.
7. Assess cumulative burden and candidate rates by source window to identify
   right-censoring.
8. Compare descriptive supported/unfiltered direct burden variants.
9. Rebuild figures, tables, and Bayesian outputs from the same detector run.

## 16. Interpretation limits

- A transition edge is shared membership, not a transmission link.
- Source attribution is proportional graph credit, not causal attribution.
- Sequence counts reflect surveillance and sequencing coverage.
- Community detection can split or merge lineage structure across windows.
- Missing parents do not establish local novelty.
- Missing descendants do not establish containment.
- Cumulative measures depend on available follow-up.
- Candidate status prioritises clusters for review; it does not verify an SSE.
