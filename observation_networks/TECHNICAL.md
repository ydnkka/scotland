# Chapter 4 Technical Contract

## Analysis Boundary

This project describes the observed genomic surveillance record and the network
objects produced by EpiLink. It does not detect superspreading-compatible
signals. Terms such as `candidate`, `candidate_tier`, burst score, and burden
score belong to `sse_detection`.

The Chapter 4 analysis uses:

- good-QC rows by default;
- Leiden resolution `0.3`;
- population-weighted SIMD groupings from `utils.load_analysis_columns`;
- all rolling windows for observation, coverage, cluster, and compatibility
  network summaries, with compatibility summaries partitioned by Pango lineage
  within each window;
- alternate retained windows, stride `2`, for the temporal transition graph so
  that its input matches Chapter 5.

## Outputs

`build_tables.py` writes:

- `cohort_summary`: high-level observed-record counts.
- `window_coverage`: one row per rolling window with sequenced counts, positive
  test denominators, and sequencing proportion.
- `window_denominator_contrasts`: policy-period summaries of rolling-window
  denominator contrasts.
- `clade_window_counts`: sequence counts and proportions by clade and window.
- `sequence_composition_by_policy`: sequence-level composition by policy period.
- `cluster_table`: one row per window-level EpiLink cluster.
- `cluster_window_summary`: cluster count, size, duration, and geographic-spread
  summaries by window.
- `cluster_period_summary`: cluster summaries by policy period.
- `cluster_attribute_composition`: modal cluster-attribute composition.
- `transition_edge_table`: directed adjacent-window cluster-transition edges.
- `transition_node_table`: cluster nodes with in/out degree, strength,
  component, role, and downstream burden summaries.
- `transition_graph_summary`: scalar baseline graph summaries.
- `transition_window_summary`: transition graph node and outgoing-edge summaries
  by retained window.
- `transition_component_summary`: weak-component summaries.
- `transition_mixing_matrix`: directed cluster-level mixing matrices for the
  transition graph.
- `transition_assortativity`: scalar assortativity summaries derived from the
  transition mixing matrices.

`build_mixing.py` writes:

- `compatibility_mixing_matrix`: weighted categorical mixing matrices for the
  within-window, within-lineage compatibility network.
- `compatibility_assortativity`: scalar assortativity summaries derived from
  those matrices, with deterministic node-block jackknife uncertainty columns
  by default.
- `compatibility_degree_assortativity`: topology diagnostics for each
  compatibility window-lineage network, including degree assortativity,
  edge-weighted degree assortativity, and weighted strength assortativity.

The compatibility build is chunked by the per-group pairwise parquet files in
`data/processed/pairwise_distances_dataset`. Each worker reads one physical
pairwise file with `utils.load_pairwise_edges`, computes the three compatibility
tables for that window-lineage group, and writes same-stem intermediate parquet
chunks under:

- `observation_networks/results/intermediate/mixing_matrix`;
- `observation_networks/results/intermediate/comp_assortativity`;
- `observation_networks/results/intermediate/deg_assortativity`.

Existing same-stem intermediate chunks are skipped unless `--force` is passed.
After worker completion, the selected intermediate chunks are concatenated into
the final `compatibility_*` tables under `observation_networks/results/tables`.
Pairwise files at or above `--giant-threshold` (`50,000,000` sparse edges by
default, with unknown costs treated as giant) are skipped by default. Pass
`--include-giants` to run the giant-file phase and include those chunks in the
final compatibility tables.

`build_simd_validation.py` writes:

- `simd_population_weighting_datazone_assignments`: one row per Data Zone with
  stored/equal-Data-Zone and population-weighted SIMD group assignments.
- `simd_population_weighting_group_summary`: group-level population share, Data
  Zone count, and SIMD rank boundaries for stored, equal-Data-Zone, and
  population-weighted groupings.
- `simd_population_weighting_movement`: cross-tabulation of Data Zone movement
  from stored/equal-Data-Zone groupings to population-weighted groups.
- `simd_population_weighting_change_summary`: counts of Data Zones moving by
  group difference.
- `simd_population_weighting_diagnostics`: reproducibility checks and maximum
  population-share deviations.
The compact `simd_population_weighting_appendix_table.tex` LaTeX fragment is
rendered from `simd_population_weighting_group_summary` by
`make_figures.py`/`lib/figs/tables.py`.

## Mixing And Assortativity

Mixing is the full categorical edge-pairing pattern. For an attribute such as
age band, health board, or SIMD quintile, the mixing matrix records weighted
edge mass from each source category to each target category.

Assortativity is a one-number summary of that matrix:

```text
r = (observed same-category edge mass - expected same-category edge mass)
    / (1 - expected same-category edge mass)
```

For undirected compatibility networks, each edge is represented in both
directions before computing the matrix. For the temporal transition graph, edge
orientation is retained.

Compatibility assortativity uncertainty is estimated with a deterministic
node jackknife. For each attribute within each physical `(window_id,
pango_lineage)` pairwise file, missing labels are handled first and the
contributing labelled vertices are counted. Attributes with up to 1,000
contributing vertices use the standard leave-one-node jackknife. Larger
attributes use balanced hash blocks with an adaptive block count:

```text
K = min(--jackknife-blocks, max(50, ceil(sqrt(n_vertices)), ceil(n_vertices / 1000)))
```

with `--jackknife-blocks` defaulting to `1,000`. The observed edge-weighted
mixing matrix is kept fixed as the point estimate, then each node or block is
left out by subtracting all edge mass touching vertices in that unit.
Assortativity is recomputed from each leave-one-unit-out matrix.

The standard error is:

```text
SE = sqrt((K - 1) / K * sum_k (r_k - mean(r_k))^2)
```

where `K` is the number of finite leave-unit-out estimates. At least five finite
replicates are required before an SE or interval is reported; otherwise the
uncertainty columns are retained but set to `NaN`. Approximate interval columns
use `r_observed +/- 1.96 * SE`. This estimates uncertainty in the observed
sequence-level compatibility-network assortativity; it is not a random-label
significance test and does not test against a topology-preserving null. This is
intentional: pairwise edges share sequences, so edge-level resampling is not a
good default for the largest files.

The compatibility assortativity table includes `uncertainty_method`,
`jackknife_vertices_used`, `jackknife_blocks_used`, `jackknife_replicates`,
`jackknife_assortativity_mean`, `assortativity_se`,
`assortativity_ci_low`, and `assortativity_ci_high`. Missing attribute labels
are dropped from the edge-level calculation by default; passing
`--missing-label` keeps them as an explicit category. Pass
`--jackknife-blocks 0` to skip uncertainty columns for development runs.

Degree and strength assortativity are separate topology diagnostics. They
summarise whether highly connected or high-strength sequences connect to
similarly connected/high-strength sequences, rather than whether metadata
categories mix assortatively.

## Disclosure

Small-cell flags are added to publication-facing composition tables using the
default threshold in `lib/config.py`. Raw analysis tables are retained for
internal reproducibility. Before thesis or manuscript use, tables that expose
cohort, geography, or fine categorical combinations must be checked against the
active PHS disclosure rules.

## Figure Inputs

Figures are built from saved tables only. This keeps manuscript figure work
decoupled from expensive parquet and pairwise-edge scans.
