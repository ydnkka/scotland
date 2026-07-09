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
  network summaries;
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
  within-window compatibility network.
- `compatibility_assortativity`: scalar assortativity summaries derived from
  those matrices. When `build_mixing.py` is run with `--n-permutations`,
  this table also includes empirical permutation p-values and null-distribution
  summary columns.
- `compatibility_degree_assortativity`: topology diagnostics for each
  compatibility window, including degree assortativity, edge-weighted degree
  assortativity, and weighted strength assortativity.

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
- `simd_population_weighting_appendix_table.tex`: compact LaTeX table for the
  observation/network appendix.

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

When `build_mixing.py` is run with `--n-permutations B`, the compatibility
assortativity table also includes an empirical permutation test for each
window-attribute pair. The test holds the observed compatibility edge list,
edge weights, and category counts fixed, then randomly permutes vertex labels
across the vertices used by that window before recomputing assortativity. This
tests whether the observed edge-weighted same-category mixing is stronger than
would be expected from the same weighted network topology and the same marginal
label distribution. The reported p-value is two-sided:

```text
p = (count(|r_permuted| >= |r_observed|) + 1) / (B + 1)
```

Permutation output columns are added only when `B > 0`: `n_permutations`,
`permutation_p_value`, `null_assortativity_mean`, and
`null_assortativity_std`. Missing attribute labels are dropped from the
edge-level calculation by default; passing `--missing-label` keeps them as an
explicit category. The base `--permutation-seed` is made stable per window and
attribute, so reruns with the same inputs, seed, and permutation count are
deterministic while parallel window execution can complete in any order.

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
