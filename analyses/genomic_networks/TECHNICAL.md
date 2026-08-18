# Genomic Surveillance and Compatibility Networks: Technical Reference

## Analysis contract

The production defaults in `lib/config.py` are good Nextclade QC, Leiden resolution 0.3, compatibility threshold 0.001, population-weighted SIMD groups, and small-cell flagging for counts 1–4. These genomic-network constants are separate from `config.yaml` even where the values match.

Main units:

- sequence-window row: one sequence in one rolling-window cluster;
- unique sequence: one row selected after collection-date/window ordering;
- cluster: one window-specific `cluster_id`;
- compatibility network: one physical `(window_id, pango_lineage)` pairwise parquet;
- Data Zone: unit for SIMD validation.

`load_sequence_data()` loads the configured analysis columns, filters resolution/QC, attaches policy variables, and recomputes requested SIMD groups using population weights. Pairwise analysis reads `data/processed/pairwise_distances_dataset/*.parquet`.

## Cluster summaries

Run:

```bash
python -m analyses.genomic_networks.build_cluster_summaries
```

The command centralises cluster-derived outputs:

| Table                                        | Unit and content                                                                                                                                                                                                            |
| -------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cohort_summary`                             | Overall sequence-window, sequence, patient, date, window, total/singleton/non-singleton cluster, clade, and lineage counts                                                                                                  |
| `window_coverage`                            | Window dates, sequences, positive-test denominator, sequencing proportion, policy, and sequences per positive test                                                                                                          |
| `window_denominator_contrasts`               | Policy-period medians/ranges of window denominators and coverage                                                                                                                                                            |
| `sequence_composition_by_policy`             | Counts/proportions for configured categorical attributes, including `Missing` and small-cell flags                                                                                                                          |
| `test_reason_by_policy_era`                  | Unique-sequence test-reason counts by epidemic era                                                                                                                                                                          |
| `vaccination_window_context`                 | Within-window deduplicated vaccination categories                                                                                                                                                                           |
| `cluster_table`                              | One row per window-specific cluster with size, Data Zone spread, duration, residential centroid distance, lineage/policy, and modal attributes                                                                              |
| `cluster_window_summary`                     | Total, singleton, and non-singleton cluster counts plus size, duration, Data Zone spread, spatial distance, and rate summaries by window; duration and explicitly labelled non-singleton characteristics exclude singletons |
| `cluster_period_summary`                     | Total, singleton, and non-singleton counts by canonical policy period, with size, duration, Data Zone spread, and spatial distance restricted to non-singleton clusters                                                     |
| `cluster_pairwise_genetic_temporal_distance` | Selected/eligible clusters and sequences, possible/observed pair counts, status, and quartile/IQR summaries for within-cluster SNP and temporal distance by window-lineage                                                  |
| `cluster_pairwise_distance_summary`          | Overall unweighted and pair-count-weighted quartiles of window-lineage median SNP and temporal distance                                                                                                                     |

### Pairwise summarises

For each selected window-lineage group, the builder loads the physical pairwise file and summarises only pairs whose endpoints share a selected cluster. Genetic and temporal rollups keep rows with `status == "ok"` and at least 10 observed within-cluster pairwise rows by default. Spatial distance is based on pairwise Euclidean distances, in kilometres, between member residential Data Zone centroids (`dz_xcoord`, `dz_ycoord`) within a cluster. Cluster, window, or group caps are development options and alter the target population.

## Compatibility mixing

Run a full build with:

```bash
python -m analyses.genomic_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
```

Each task reads one pairwise parquet, retains `epilink_compatibility > threshold`, and joins endpoint attributes from the genomic-network sequence data. Default attributes are sex, age band, age group, SIMD quintile, urban/rural class, local authority, and Health Board. Pairs with a missing endpoint label are dropped attribute-by-attribute unless `--missing-label` is supplied. `n_edges_observed` is the retained compatibility-edge count before attribute missingness filtering; `n_edges_used` is the count used for that attribute after endpoint labels are checked. `--min-edges` keeps rows below the threshold but reports `NaN` estimates with `skipped_reason`.

For each attribute, an undirected edge contributes its weight symmetrically to a category mixing matrix. With the matrix normalised to `e` and row/column marginals `(a, b)`, nominal assortativity is:

```text
r = (trace(e) - sum(a * b)) / (1 - sum(a * b))
```

Outputs:

- `compatibility_mixing_matrix_bootstrap.parquet`: non-zero source/target category cells by attribute, window, lineage, and pairwise stem;
- `compatibility_assortativity_bootstrap`: categorical assortativity point estimates, same-category observed/expected weights, edge counts, categories, and bootstrap uncertainty;
- `compatibility_degree_assortativity_bootstrap`: unweighted degree, compatibility-weighted degree, and strength assortativity plus topology summaries; bootstrap uncertainty is reported for strength assortativity.

### Multiplier-bootstrap uncertainty

Assortativity uncertainty is estimated with an edge-weight multiplier bootstrap, not a vertex jackknife.

For each attribute within each window-lineage graph, the observed point estimate is computed once from the retained EpiLink weights. For each bootstrap replicate (b), the original edge weights (w_e) are perturbed by independent exponential multipliers:

```text
g_be ~ Exponential(1)
w_be* = w_e * g_be
r_b* = assortativity(w_be*)
```

The default is 500 replicates (`--bootstrap-replicates 500`), a 95% percentile interval (`--bootstrap-alpha 0.05`), and base seed 123 (`--bootstrap-seed 123`). Use `--bootstrap-replicates 0` to skip uncertainty while retaining point estimates. Finite bootstrap replicates are counted in `bootstrap_finite_replicates`; intervals and standard errors are `NaN` if no finite replicate estimates are available.

For (K) finite bootstrap estimates (r_b*):

```text
SE = sample_sd(r_b*)
CI = quantile(r_b*, alpha / 2), quantile(r_b*, 1 - alpha / 2)
```

For categorical assortativity, bootstrap replicates recompute the full weighted category matrix after the attribute-specific missing-label filter. For strength assortativity, each replicate recomputes node strengths from the perturbed weights before recomputing numeric assortativity. The unweighted degree and compatibility-weighted degree assortativity diagnostics are point estimates only.

These intervals describe sensitivity of the edge-weighted statistic to perturbing the retained edge weights. They are not a random-mixing null, do not model unsequenced infections, and do not remove whole vertices or blocks; dependence from shared endpoints, repeated rolling-window sequences, and lineage/window reuse therefore remains an interpretation limit.

### Scheduling and restart behaviour

The sparse-edge manifest supplies `pairwise_stem` and `sparse_edges`; misses fall back to file size and unknown costs are treated as giant. Small and giant files use separate worker pools.

Each task writes same-stem parquets below:

```text
results/intermediate/mixing_matrix_bootstrap/
results/intermediate/comp_assortativity_bootstrap/
results/intermediate/deg_assortativity_bootstrap/
```

Existing complete chunks are reused unless `--force` is passed. Chunk filenames do not encode configuration, so settings must not be mixed in one intermediate set. This includes attributes, threshold, missing-label handling, bootstrap settings, minimum-edge filtering, and giant-file inclusion.

## Sensitivity analyses

`build_sensitivity_tables.py --only leiden` compares each available resolution with baseline 0.3. Window-level outputs include counts, size/singleton summaries, duration/spread, adjusted Rand index, and baseline differences. `--include-ami` adds exact adjusted mutual information.

`--only sparsification` scans this default threshold grid:

```text
0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1
```

It reports edge/weight retention, retained mean degree, scan coverage, partial-scan status, and baseline ratios/differences per group and in aggregate. Outputs flagged `estimated_from_partial_scan` are approximate.

Tables:

```text
leiden_resolution_window_sensitivity
leiden_resolution_sensitivity_summary
sparsification_threshold_sensitivity
sparsification_threshold_sensitivity_summary
```

## SIMD validation

`build_simd_validation.py` compares stored, equal-Data-Zone, and population-weighted rank groups. It writes:

- `simd_population_weighting_datazone_assignments`;
- `simd_population_weighting_group_summary`;
- `simd_population_weighting_movement`;
- `simd_population_weighting_change_summary`;
- `simd_population_weighting_diagnostics`.

## Figures, disclosure, and checks

`python -m results.make_figures --domain genomic_networks` and
`python -m results.make_tables --domain genomic_networks` read saved
analysis-local tables and write project-level figures/LaTeX fragments under
`results/figures/` and `results/tables/`; they do not rescan source data.

Counts below 5 are flagged only in publication-facing composition tables. Raw/internal outputs remain potentially disclosive. Mixing drops missing endpoint labels by default, whereas composition retains a `Missing` level.

Before reporting, record input versions, filters, thresholds, resolutions, windows, attributes, bootstrap settings, minimum-edge filtering, giant-file inclusion, and any caps/partial scans. Confirm pairwise stems, uncertainty completeness, SIMD diagnostics, disclosure flags, and figure/table provenance.

## Interpretation limits

- The sequenced cohort is selected and coverage varies over time and place.
- Rolling windows repeat sequences.
- Compatibility edges are not transmission links.
- Clustering depends on sparsification and Leiden resolution.
- Assortativity depends on network/category prevalence and missingness.
- Window-specific clusters are not persistent outbreaks; temporal continuity is analysed separately by SSE detection.
