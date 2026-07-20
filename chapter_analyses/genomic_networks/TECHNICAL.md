# Chapter 4 Genomic Surveillance and Compatibility Networks: Technical Reference

## Analysis contract

The production defaults in `lib/config.py` are good Nextclade QC, Leiden resolution 0.3, compatibility threshold 0.001, population-weighted SIMD groups, and small-cell flagging for counts 1–4. These Chapter 4 constants are separate from `config.yaml` even where the values match.

Main units:

- sequence-window row: one sequence in one rolling-window cluster;
- unique sequence: one row selected after collection-date/window ordering;
- cluster: one window-specific `cluster_id`;
- compatibility network: one physical `(window_id, pango_lineage)` pairwise parquet;
- Data Zone: unit for SIMD validation.

`load_chapter4_sequence_data()` loads the configured analysis columns, filters resolution/QC, attaches policy variables, and recomputes requested SIMD groups using population weights. Pairwise analysis reads `data/processed/pairwise_distances_dataset/*.parquet`.

## Core tables

Run:

```bash
python -m chapter_analyses.genomic_networks.build_cluster_tables
```

The command loads data once and writes:

| Table | Unit and content |
| --- | --- |
| `cohort_summary` | Overall sequence-window, sequence, patient, date, window, total/singleton/non-singleton cluster, clade, and lineage counts |
| `window_coverage` | Window dates, sequences, positive-test denominator, sequencing proportion, policy, and sequences per positive test |
| `window_denominator_contrasts` | Policy-period medians/ranges of window denominators and coverage |
| `clade_window_counts` | Clade counts by window |
| `sequence_composition_by_policy` | Counts/proportions for configured categorical attributes, including `Missing` and small-cell flags |
| `test_reason_by_policy_era` | Unique-sequence test-reason counts by epidemic era |
| `vaccination_context_by_policy` | Unique-sequence vaccination categories by policy period |
| `vaccination_window_context` | Within-window deduplicated vaccination categories |
| `cluster_table` | One row per window-specific cluster with size, spread, duration, lineage/policy, and modal attributes |
| `cluster_window_summary` | Total, singleton, and non-singleton cluster counts plus size, duration, spread, and rate summaries by window; duration and explicitly labelled non-singleton characteristics exclude singletons |
| `cluster_period_summary` | Total, singleton, and non-singleton counts by policy period, with size, duration, and spread restricted to non-singleton clusters |
| `cluster_attribute_composition` | Modal cluster attributes by policy period, split into total, singleton, and non-singleton counts and proportions with corresponding small-cell flags |

All are CSV and parquet except `clade_window_counts` and `cluster_table`, which are parquet only.

Vaccination results describe sequenced cases; they do not estimate vaccine effectiveness.

## Compatibility mixing

Run a full build with:

```bash
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
```

Each task reads one pairwise parquet, retains `epilink_compatibility > threshold`, and joins endpoint attributes from the Chapter 4 sequence data. Default attributes are sex, age band, age group, SIMD quintile, urban/rural class, local authority, and Health Board. Pairs with a missing endpoint label are dropped unless `--missing-label` is supplied.

For each attribute, an undirected edge contributes its weight symmetrically to a category mixing matrix. With the matrix normalised to (e), row/column marginals (a,b), nominal assortativity is:

```text
r = (trace(e) - sum(a * b)) / (1 - sum(a * b))
```

Outputs:

- `compatibility_mixing_matrix.parquet`: category-pair weights, counts, proportions, window, lineage, and pairwise stem;
- `compatibility_assortativity.{csv,parquet}`: point estimates, same-category weights, contributing edges/categories, and uncertainty;
- `compatibility_degree_assortativity.{csv,parquet}`: unweighted degree, compatibility-weighted degree, and strength assortativity plus topology summaries.

### Vertex jackknife

Uncertainty removes labelled vertices rather than treating incident edges as independent.

- Up to 1,000 vertices: leave one vertex out.
- Larger networks: deterministic balanced blocks, capped by `--jackknife-blocks` (default 1,000).
- Requested block count: `min(cap, max(50, ceil(sqrt(n)), ceil(n/1000)))`.
- At least five finite replicates are required; `--jackknife-blocks 0` disables uncertainty.

For (K) finite replicate estimates (r_k):

```text
SE = sqrt((K - 1) / K * sum((r_k - mean(r_k))^2))
CI = observed r ± 1.96 * SE
```

These intervals measure sensitivity to observed vertices/blocks, not a random-mixing null.

### Scheduling and restart behaviour

The sparse-edge manifest supplies `pairwise_stem` and `sparse_edges`; misses fall back to file size and unknown costs are treated as giant. Small and giant files use separate worker pools.

Each task writes same-stem parquets below:

```text
results/intermediate/mixing_matrix/
results/intermediate/comp_assortativity/
results/intermediate/deg_assortativity/
```

Existing complete chunks are reused unless `--force` is passed. Chunk filenames do not encode configuration, so settings must not be mixed in one intermediate set.

## Pairwise-distance summary

`build_cluster_pairwise_distance_summary.py` selects clusters at the requested resolution/QC and minimum observed cluster size (defaults: 0.3, good, 2). For each selected window-lineage group it loads the physical pairwise file and summarises only pairs whose endpoints share a selected cluster.

`cluster_pairwise_distance_summary.{csv,parquet}` records selected/eligible clusters and sequences, possible/observed pair counts, status, and quartile/IQR summaries for SNP and temporal distance. Cluster or group caps alter the target population and are development options.

## Sensitivity analyses

`build_sensitivity_tables.py --only leiden` compares each available resolution with baseline 0.3. Window-level outputs include counts, size/singleton summaries, duration/spread, adjusted Rand index, and baseline differences. `--include-ami` adds exact adjusted mutual information.

`--only sparsification` scans this default threshold grid:

```text
0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1
```

It reports edge/weight retention, retained mean degree, scan coverage, partial-scan status, and baseline ratios/differences per group and in aggregate. Outputs flagged `estimated_from_partial_scan` are approximate.

Tables:

```text
leiden_resolution_window_sensitivity.{csv,parquet}
leiden_resolution_sensitivity_summary.{csv,parquet}
sparsification_threshold_sensitivity.{csv,parquet}
sparsification_threshold_sensitivity_summary.{csv,parquet}
```

## SIMD validation

`build_simd_validation.py` compares stored, equal-Data-Zone, and population-weighted rank groups. It writes:

- `simd_population_weighting_datazone_assignments.parquet`;
- `simd_population_weighting_group_summary.{csv,parquet}`;
- `simd_population_weighting_movement.{csv,parquet}`;
- `simd_population_weighting_change_summary.{csv,parquet}`;
- `simd_population_weighting_diagnostics.{csv,parquet}`.

## Figures, disclosure, and checks

`make_figures.py` reads saved result tables and writes current Chapter 4 figures plus LaTeX fragments under `results/figures/`; it does not rescan source data.

Counts below 5 are flagged only in publication-facing composition tables. Raw/internal outputs remain potentially disclosive. Mixing drops missing endpoint labels by default, whereas composition retains a `Missing` level.

Before reporting, record input versions, filters, thresholds, resolutions, windows, attributes, jackknife settings, giant-file inclusion, and any caps/partial scans. Confirm pairwise stems, uncertainty completeness, SIMD diagnostics, disclosure flags, and figure/table provenance.

## Interpretation limits

- The sequenced cohort is selected and coverage varies over time and place.
- Rolling windows repeat sequences.
- Compatibility edges are not transmission links.
- Clustering depends on sparsification and Leiden resolution.
- Assortativity depends on network/category prevalence and missingness.
- Window-specific clusters are not persistent outbreaks; temporal continuity is analysed separately in Chapter 5.
