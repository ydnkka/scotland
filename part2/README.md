# Part 2 Cluster Categorisation And Vaccination

This directory contains the Part 2 cluster-categorisation and vaccination
descriptive passes.

## Cluster Categories

`categorise_clusters.py` owns the shared Part 2 category rules and can also be
run as a standalone cached-table pass. It reads the cached Part 1 cluster table
through `utils.data.load_main_cluster_table` and writes cluster-level categories
for non-singleton clusters only. The singleton filter is applied before category
thresholds are estimated, so the size and geography categories describe
variation among reconstructed multi-case clusters rather than separating
singletons from everything else.

It writes categories for:

- cluster size: `small/moderate`, `large`, `very large`
- geographic dispersion, based on `cluster_n_datazones`
- mean-cluster SIMD quintile, where quintile 1 is most deprived
- SIMD, age, sex, and joint SIMD-age-sex excess mixing:
  `less mix`, `baseline`, `more mix`, or `not available`

Run from the repository root:

```bash
conda run -n PhD python part2/categorise_clusters.py
```

Default thresholds, estimated after filtering to non-singleton clusters:

- large size: `cluster_size >= ceil(90th percentile)`
- very large size: `cluster_size >= ceil(99th percentile)`
- large geographic dispersion: `cluster_n_datazones >= ceil(90th percentile)`
- very large geographic dispersion: `cluster_n_datazones >= ceil(99th percentile)`
- baseline mixing: excess discordance within `+/-0.01`

Outputs:

- `cache/cluster_categories.parquet`: full cached cluster table plus categories
- `tables/cluster_categories.csv`: compact row-level category table
- `tables/cluster_category_summary.csv`: one-way category counts
- `tables/cluster_category_thresholds.csv`: threshold/rule audit table
- `tables/cluster_category_combinations.csv`: cross-category cluster counts

Thresholds can be overridden at the command line, for example:

```bash
conda run -n PhD python part2/categorise_clusters.py \
  --large-size-min 10 \
  --very-large-size-min 50 \
  --mixing-baseline-band 0.02
```

## Cluster Characterisation

`cluster_characterisation.py` rebuilds primary-resolution cluster aggregates from
`data/processed/scotland_clustering_analysis_dataset.parquet` through
`utils.data.load_analysis_columns_pandas` instead of using the cached Part 1
table. It imports the shared size, geography, and SIMD category helpers from
`categorise_clusters.py`, then adds vaccination, demographic, and SIMD-domain
summaries including:

- cluster proportion vaccinated
- cluster vaccination profile: none, mixed-membership, or all vaccinated
- vaccination-status mixing category: homogeneous, baseline, mixed, or not available
- mean/median dose among vaccinated cluster members
- booster proportion
- mean/median days since latest prior vaccination
- index-case vaccination status, dose, and days since vaccination
- cluster size/geographic/SIMD categories for wave summaries
- SIMD domain quintiles for income, employment, education, health, geographic
  access, crime, and housing

The full vaccination cluster cache keeps all clusters, including singletons, but
category thresholds and category summary tables are based on non-singleton
clusters.

Age-stratified vaccination summaries use rollout-informed groups that approximate
the Scottish UK's Joint Committee on Vaccination and Immunisation (JCVI) eligibility sequence with the available five-year age bands:
`00-14`, `15-19`, `20-29`, `30-39`, `40-49`, `50-54`, `55-59`, `60-64`,
`65-69`, `70-74`, and `75+`. These cannot exactly distinguish `12-15`,
`16-17`, `18-19`, or `80+` because those groups are not separately identifiable
in the source `age_band` field.

Run:

```bash
conda run -n PhD python part2/cluster_characterisation.py
```

Main outputs:

- `cache/vaccination_cluster_table.parquet`
- `tables/vaccination_cluster_table.csv`
- `tables/vaccination_case_weekly_summary.csv`
- `tables/vaccination_case_weekly_simd_domain_summary.csv`
- `tables/vaccination_cluster_wave_category_summary.csv`
- `tables/vaccination_cluster_wave_simd_domain_summary.csv`
- `tables/vaccination_cluster_weekly_category_summary.csv`
- `tables/vaccination_descriptives.csv`
- `tables/vaccination_key_questions.csv`

The Part 2 framing is descriptive cluster characterisation rather than causal
vaccine-effectiveness estimation. Vaccination status among sequenced cases is
strongly confounded by rollout time, age eligibility, variant wave, prior
immunity, testing behaviour, and sequencing selection.

A useful primary descriptive question is:

> How did vaccination history, dose profile, and vaccination-status mixing
> characterise Scottish genomic clusters across waves, size categories,
> geographic-dispersion categories, and social groups?

The vaccination-status mixing category uses observed-minus-expected pairwise
discordance for binary vaccination status. Expected discordance is computed
within the same `window_id + pango_lineage` strata. Categories are:

- `homogeneous`: less vaccinated/unvaccinated mixing than expected
- `baseline`: within `+/-0.01` of expected mixing
- `mixed`: more vaccinated/unvaccinated mixing than expected
- `not available`: usually singleton clusters or missing expectation
