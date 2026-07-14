# Chapter 4 Genomic Surveillance and Compatibility Networks: Technical Reference

This document is the implementation and output contract for the Chapter 4 observation and compatibility-network analysis. For commands and the directory map, see `README.md`.

## 1. Analysis boundary

This package describes the observed Scottish genomic-surveillance record and the window-specific network objects produced by EpiLink. It does not construct temporal cluster transitions or detect superspreading-compatible signals. Candidate labels, upstream novelty, burst scores, burden scores, null calibration, and candidate tiers belong to `chapter_analyses.sse_detection`.

The default Chapter 4 analysis uses:

- good-QC sequence rows (`QC_FILTER = "good"`);
- Leiden resolution `0.3` (`ANALYSIS_RESOLUTION = 0.3`);
- compatibility threshold `0.001` (`SPARSIFICATION_THRESHOLD = 0.001`);
- population-weighted SIMD categories attached by `utils.load_analysis_columns`;
- all rolling windows for cohort, coverage, vaccination, and cluster summaries;
- compatibility networks partitioned by physical `(window_id, pango_lineage)` pairwise files;
- a publication-facing small-cell flag below 5 observations (`DISCLOSURE_MIN_CELL = 5`).

`TRANSITION_WINDOW_STRIDE = 2` remains in the Chapter 4 configuration for shared context, but transition construction itself is owned by Chapter 5.

## 2. Data sources and units of analysis

### 2.1 Sequence-window analysis data

`load_chapter4_sequence_data()` loads the columns declared in `CHAPTER4_COLUMNS` through repository utilities and attaches policy variables. The input contains repeated sequence-window memberships, so functions explicitly deduplicate when their target unit is a unique sequence.

The main units are:

- sequence-window row: one sequence membership in one rolling-window cluster;
- unique sequence: one earliest/ordered record selected by `sequence_level_frame(...)`;
- cluster: one window-specific `cluster_id` within `window_id` and `window_idx`;
- window-lineage compatibility network: one physical pairwise parquet file for a `(window_id, pango_lineage)` group;
- Data Zone: the geographic unit used by SIMD population-weighting validation.

### 2.2 Pairwise compatibility data

Compatibility scans read `data/processed/pairwise_distances_dataset/*.parquet` with `utils.load_pairwise_edges`. Required fields for mixing are `id1`, `id2`, and `epilink_compatibility`; pairwise-distance summaries additionally use `snp_distance` and `temporal_distance`.

`data/processed/sparsified_edge_counts_by_window_lineage.parquet` supplies `pairwise_stem` and `sparse_edges` scheduling metadata. Manifest misses fall back to file size where possible; unknown costs are classified as giant.

### 2.3 Multi-resolution clustering data

Leiden sensitivity uses `data/processed/scotland_clustering_analysis_dataset.parquet`, retaining resolution, QC, window, sequence, cluster-size, geographic-spread, and duration fields.

## 3. Core cohort and cluster build

Run:

```bash
python -m chapter_analyses.genomic_networks.build_cluster_tables
```

The command loads the Chapter 4 sequence-window data once, optionally restricts the first `--max-windows N` windows, and writes all core tables. `--skip-transition`, `--max-transition-windows`, and `--transition-window-stride` are deprecated no-ops retained for compatibility.

### 3.1 Unique-sequence cohort summaries

`sequence_level_frame(...)` sorts by available collection date and window index, then retains one record per `sequence_id`. `cohort_summary` reports sequence-window rows, unique sequences, unique patients, date range, windows, clusters, clades, and Pango lineages where available.

### 3.2 Window coverage and denominators

`window_coverage` contains one row per rolling window, including dates, sequence count, positive-test denominator, sequencing proportion, policy assignment, and:

```text
sequences_per_positive_test = wn_no_sequences / wn_positive_tests
```

`window_denominator_contrasts` aggregates window counts and sequencing proportions by policy period using medians and ranges.

### 3.3 Sequence composition

`sequence_composition_by_policy` is a long table over the configured sex, age-band, age-group, SIMD, urban/rural, local-authority, and health-board attributes. Missing values are retained as an explicit `Missing` category for composition counts. Each row records the category count, group total, proportion, and `small_cell` flag.

### 3.4 Vaccination context

Vaccination summaries classify unique sequences as unvaccinated, one dose, two doses, booster/three-plus doses, or vaccinated with dose unknown. `vaccination_context_by_policy` summarises unique sequences by policy period; `vaccination_window_context` deduplicates within each window and summarises the same fields by rolling window. These are descriptive properties of sequenced records, not vaccine-effectiveness estimates.

### 3.5 Window-specific clusters

`cluster_table` contains one row per `(cluster_id, window_id, window_idx)`. It records cluster size, observed unique sequences and Data Zones, duration, first/last collection dates, lineage/clade/VOC, policy context, and modal categorical attributes. For each configured categorical attribute it also records the modal fraction and number of observed levels.

`cluster_window_summary` reports cluster counts, sequence memberships, median/p90/maximum size, duration, Data Zone spread, and clusters per 1,000 sequence memberships by window. `cluster_period_summary` reports analogous summaries by policy period. `cluster_attribute_composition` counts modal cluster categories by policy period and adds disclosure flags.

## 4. Compatibility-network construction

Run a full scan with:

```bash
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
```

For each physical window-lineage pairwise file, compatibility edges above `--compatibility-threshold` are loaded. Sequence identifiers are network vertices, `epilink_compatibility` is the edge weight, and node attributes come from the Chapter 4 sequence data.

The default categorical attributes are sex, age band, age group, SIMD quintile, urban/rural class, local authority, and health board. Missing endpoint labels are dropped by default; `--missing-label LABEL` retains them as an explicit category.

## 5. Mixing matrices and categorical assortativity

For an attribute, the weighted mixing matrix records compatibility-edge weight between source and target categories. Undirected edges are represented symmetrically so both directions contribute to the matrix.

Let `e` be the mixing matrix normalised to sum to one, `a` its row marginals, and `b` its column marginals. Nominal assortativity is:

```text
r = (trace(e) - sum(a * b)) / (1 - sum(a * b))
```

`compatibility_mixing_matrix` stores category-pair edge weight, contribution count, weight proportion, window, lineage, and physical pairwise stem. `compatibility_assortativity` stores the scalar point estimate, observed and expected same-category weights, contributing edges/categories, uncertainty fields, and pairwise stem.

Assortativity describes categorical pairing conditional on the observed compatibility network. It is not evidence of transmission, causality, or preferential contact without additional assumptions.

## 6. Jackknife uncertainty

Uncertainty is estimated over contributing labelled vertices rather than independent edges because many pairwise edges share endpoints.

- Attributes with at most 1,000 contributing vertices use leave-one-node jackknife.
- Larger attributes use deterministic balanced hash blocks.
- The adaptive block count is capped by `--jackknife-blocks`, which defaults to 1,000.
- `--jackknife-seed` defaults to 42 and participates in deterministic block assignment.
- At least five finite leave-unit-out estimates are required.
- `--jackknife-blocks 0` skips uncertainty calculation.

For larger networks, the requested block count is:

```text
K = min(--jackknife-blocks, max(50, ceil(sqrt(n_vertices)), ceil(n_vertices / 1000)))
```

Each replicate removes all mixing-matrix mass touching the omitted node or block and recomputes assortativity. With finite replicate estimates `r_k`, the standard error is:

```text
SE = sqrt((K - 1) / K * sum((r_k - mean(r_k))^2))
```

Approximate intervals are `r_observed ± 1.96 × SE`. These intervals quantify vertex-jackknife uncertainty in the observed compatibility network; they are not random-label or topology-preserving null tests.

The uncertainty contract includes `uncertainty_method`, `jackknife_vertices_used`, `jackknife_blocks_used`, `jackknife_replicates`, `jackknife_assortativity_mean`, `assortativity_se`, `assortativity_ci_low`, and `assortativity_ci_high`.

## 7. Degree and strength assortativity

`compatibility_degree_assortativity` provides topology diagnostics distinct from categorical assortativity:

- `degree_assortativity`: correlation of unweighted endpoint degrees;
- `weighted_degree_assortativity`: endpoint-degree correlation weighted by compatibility strength;
- `strength_assortativity`: correlation of endpoint weighted strengths;
- network size, edge count, total weight, mean/maximum degree, and mean/maximum strength.

These statistics ask whether highly connected or high-strength sequences connect to similarly connected or high-strength sequences. They do not measure categorical mixing.

## 8. Chunking, scheduling, and restart behaviour

The mixing build creates one task per physical pairwise parquet file. Small and giant tasks are scheduled separately using an edge-count threshold of 50,000,000 by default.

- Small tasks use `--workers`.
- Giant tasks are skipped unless `--include-giants` is supplied.
- Included giant tasks use `--giant-workers`, defaulting to one.
- `--dry-run` prints the schedule without scanning pairwise files.
- `--progress-every N` controls aggregate INFO logging; per-file events use DEBUG.

Each task writes same-stem parquet chunks under `results/intermediate/mixing_matrix`, `results/intermediate/comp_assortativity`, and `results/intermediate/deg_assortativity`. Existing complete chunks are reused unless `--force` is passed. After processing, chunks for the selected tasks are concatenated into final tables under `results/tables`.

Intermediate outputs are configuration-dependent. When changing the threshold, attributes, missing-label policy, or uncertainty settings, use a clean compatible intermediate set or `--force`; do not assume same-stem chunks from a different configuration are interchangeable.

## 9. Cluster pairwise-distance summary

`build_cluster_pairwise_distance_summary.py` selects non-singleton clusters at the requested resolution and QC filter, combines selected members within each window-lineage group, loads the corresponding physical pairwise file, and summarises within-selection SNP and temporal distances.

The output contains selected cluster/sequence counts, possible and observed pair counts, status, and median, quartile, and IQR summaries for SNP and temporal distance. Defaults are resolution 0.3, good QC, and minimum cluster size 2.

When `--max-clusters-per-window-lineage` is set, `--cluster-selection` chooses the largest, smallest, or cluster-ID-ordered subset. Development caps change the estimand and must not be presented as full-analysis results.

## 10. Leiden-resolution sensitivity

`build_sensitivity_tables.py --only leiden` compares window-level cluster structure across available Leiden resolutions against baseline resolution 0.3.

`leiden_resolution_window_sensitivity` records cluster count and density, size distribution, singleton cluster/sequence shares, duration, geographic spread, adjusted Rand index against baseline, and absolute/relative differences from baseline. Exact adjusted mutual information is optional through `--include-ami` because it is substantially slower on large high-resolution windows.

`leiden_resolution_sensitivity_summary` aggregates the window-level results by resolution using medians, quartiles, and extrema.

## 11. Sparsification sensitivity

`build_sensitivity_tables.py --only sparsification` evaluates compatibility thresholds against baseline `0.001`. The default grid is:

```text
0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1
```

For each physical pairwise group and threshold, `sparsification_threshold_sensitivity` records retained edge/weight totals and fractions, retained mean degree, scan coverage, partial-scan status, and ratios/differences from baseline. `sparsification_threshold_sensitivity_summary` aggregates group and pooled retention summaries by threshold.

`--max-row-groups-per-file` and `--max-rows-per-file` create approximate scans. Any output with `estimated_from_partial_scan = True` must be identified as approximate and should not replace the full sensitivity analysis.

## 12. SIMD population-weighting validation

`build_simd_validation.py` reconstructs equal-Data-Zone and population-weighted SIMD groups and compares them with stored assignments. The default is five groups; `--n-groups` also supports 10 and 20.

Outputs are:

- `simd_population_weighting_datazone_assignments`: one row per Data Zone with stored, equal-Data-Zone, and population-weighted assignments;
- `simd_population_weighting_group_summary`: population share, Data Zone count, and SIMD-rank boundaries by method and group;
- `simd_population_weighting_movement`: cross-tabulation of assignment movement;
- `simd_population_weighting_change_summary`: counts by magnitude/direction of group change;
- `simd_population_weighting_diagnostics`: reproducibility and population-balance checks.

## 13. Output contract

### 13.1 Core tables

`build_cluster_tables.py` writes:

- `cohort_summary.{csv,parquet}`;
- `window_coverage.{csv,parquet}`;
- `window_denominator_contrasts.{csv,parquet}`;
- `clade_window_counts.parquet`;
- `sequence_composition_by_policy.{csv,parquet}`;
- `vaccination_context_by_policy.{csv,parquet}`;
- `vaccination_window_context.{csv,parquet}`;
- `cluster_table.parquet`;
- `cluster_window_summary.{csv,parquet}`;
- `cluster_period_summary.{csv,parquet}`;
- `cluster_attribute_composition.{csv,parquet}`.

### 13.2 Compatibility tables

`build_mixing.py` writes:

- `compatibility_mixing_matrix.parquet`;
- `compatibility_assortativity.{csv,parquet}`;
- `compatibility_degree_assortativity.{csv,parquet}`.

### 13.3 Supplementary analysis tables

- `cluster_pairwise_distance_summary.{csv,parquet}`;
- `leiden_resolution_window_sensitivity.{csv,parquet}`;
- `leiden_resolution_sensitivity_summary.{csv,parquet}`;
- `sparsification_threshold_sensitivity.{csv,parquet}`;
- `sparsification_threshold_sensitivity_summary.{csv,parquet}`;
- the five SIMD validation tables described above.

Table writers use parquet and/or CSV explicitly by artifact. Parquet is the authoritative typed format; CSV is provided for inspection and reporting where configured.

## 14. Figures and LaTeX tables

`make_figures.py` reads saved result tables only; it does not rescan raw analysis or pairwise data. `--skip-missing` skips artifacts whose required tables are absent.

The current figure set covers sequence composition, cluster landscape, baseline assortativity, selected mixing matrices, vaccination context, parameter sensitivity, compatibility topology, SIMD population weighting, assortativity confidence intervals, and within-cluster pairwise distances.

LaTeX builders write cohort objects, policy denominators, vaccination context, cluster-period summary, compatibility-assortativity summary, and SIMD population-weighting fragments under `results/figures/`.

## 15. Disclosure and missingness

Publication-facing composition tables flag counts from 1 through 4 using the configured threshold of 5. Raw internal tables are not automatically disclosure-controlled outputs. Any table exposing cohort, geography, fine categories, or combinations must be reviewed against the active PHS disclosure rules before dissemination.

Composition tables retain missing values as an explicit category, whereas compatibility mixing drops missing endpoint attributes unless `--missing-label` is supplied. This distinction must be documented when comparing composition and mixing results.

## 16. Reproducibility and quality control

Before reporting Chapter 4 results:

1. Record the input dataset versions, QC filter, Leiden resolution, compatibility threshold, window selection, and attribute set.
2. Run `build_mixing --dry-run` and record which giant files are included or excluded.
3. Confirm final compatibility tables contain the expected physical pairwise stems.
4. Confirm uncertainty was computed with the intended jackknife configuration and inspect missing/finite interval counts.
5. Confirm pairwise-distance and sensitivity outputs were not produced from development caps or partial scans.
6. Review SIMD diagnostics and disclosure flags.
7. Regenerate all figures and LaTeX fragments from the same saved-table set.
8. Keep Chapter 4 outputs free of Chapter 5 candidate fields.

## 17. Interpretation limits

- The observed sequenced cohort is not the complete infected population.
- Window coverage and geographic sequencing intensity vary over time.
- Compatibility edges indicate EpiLink compatibility, not transmission.
- Rolling windows repeat sequences and must not be treated as independent cohorts.
- Community-detection output depends on resolution and graph sparsification.
- Categorical assortativity depends on category prevalence, missingness, and the observed network.
- Jackknife intervals quantify observed-network sensitivity to vertices or vertex blocks; they are not null-hypothesis significance tests.
- Vaccination tables describe the sequenced record and do not estimate vaccine effectiveness.
- Chapter 4 cluster summaries are window-specific and should not be interpreted as persistent outbreak entities without the Chapter 5 transition framework.
