# Part 1 Main Analysis

This directory contains the Part 1 modelling pass using outcome-specific
models for the primary cluster outcomes.

## Question

After accounting for lineage, calendar time, local incidence, sequencing
intensity, and test positivity, are socioeconomic deprivation and surveillance
conditions associated with larger, more geographically dispersed, or more
socially mixed SARS-CoV-2 genomic clusters in Scotland?

## Analysis Population

- Input: `data/processed/scotland_clustering_analysis_dataset.parquet`
- QC filter: `nextclade_qc == "good"`
- Primary Leiden resolution: `0.3`
- Unit of analysis: one inferred cluster at the primary resolution
- Cluster rows after model-field filtering: 193,112
- Non-singleton clusters for mixing models: 84,067

The main analysis uses a single primary Leiden resolution so that the same
underlying sampled genomes are not repeatedly analysed across resolution values.

## Covariate Adjustment

The models adjust for:

- Mean SIMD deprivation, with higher values meaning more deprived clusters
- Local cumulative incidence
- Local sequencing fraction
- Window-level sequencing proportion
- Local test positivity
- Pango lineage, with rare lineages pooled into `Other rare lineages`
- Calendar time using an 8 df B-spline over `window_idx`

Calendar time is modelled with a spline rather than `window_id` fixed effects
because `window_id` fixed effects would absorb the window-level sequencing
proportion, making that surveillance coefficient non-identifiable.

## Count Outcomes

Cluster size and number of datazones have large structural masses at their
minimum values. They are therefore modelled with two-part hurdle models:

- Binary hurdle component: binomial GLM with logit link
- Positive count component: zero-truncated negative binomial model

Outcome definitions:

- `cluster_size`: hurdle is `cluster_size > 1`; positive count is
  `cluster_size - 1`
- `geographic_dispersion`: hurdle is `cluster_n_datazones > 1`; positive count
  is `cluster_n_datazones - 1`

The size-adjusted geographic-dispersion sensitivity is fitted only for the
positive zero-truncated count component. The binary hurdle component is not used
for the size-adjusted sensitivity model because cluster size partly defines
whether a cluster can exceed the structural minimum.

Duration remains available only in descriptive summaries and the supplementary
outcome-distribution figure because the fixed three-week clustering windows
mechanically constrain the observed span.

## Mixing Outcomes

Within-cluster mixing is measured as observed pairwise discordance minus the
expected pairwise discordance among sampled cases from the same lineage and
calendar window. Positive values mean clusters are more mixed than expected;
negative values mean clusters are more homogeneous than expected.

The main mixing outcomes are:

- SIMD quintile excess mixing
- Age-band excess mixing
- Sex excess mixing
- Joint SIMD-age-sex profile excess mixing

Mixing outcomes are fitted with linear models with the same adjustment set plus
cluster size.

## Mixing Predictors

The same observed-minus-expected excess mixing metrics are also fitted as
predictors in sensitivity models for the cluster count outcomes. These models
are reported separately from the primary hurdle/ZTNB tables because mixing is
only defined for clusters with at least two valid cases. The cluster-size hurdle
component is therefore skipped in the mixing-predictor tables; positive count
components and the geographic-spread hurdle component are estimated on the
complete-case population with non-missing mixing metrics.

## Outputs

Tables:

- `tables/dataset_descriptives.csv`
- `tables/covariate_scaling.csv`
- `tables/hurdle_count_model_results.csv`
- `tables/hurdle_count_model_diagnostics.csv`
- `tables/mixing_predictor_hurdle_count_model_results.csv`
- `tables/mixing_predictor_hurdle_count_model_diagnostics.csv`
- `tables/mixing_model_results.csv`
- `tables/mixing_model_diagnostics.csv`
- `tables/loglinear_count_model_results.csv`
- `tables/mixing_predictor_loglinear_count_model_results.csv`
- `tables/simd_domain_hurdle_count_model_results.csv`
- `tables/simd_domain_mixing_predictor_hurdle_count_model_results.csv`
- `tables/simd_domain_mixing_predictor_hurdle_count_model_diagnostics.csv`
- `tables/simd_domain_quintile_mixing_model_results.csv`
- `tables/simd_domain_demographic_mixing_model_results.csv`
- `tables/wave_specific_domain_demographic_mixing_model_results.csv`
- `tables/wave_specific_hurdle_count_model_results.csv`
- `tables/wave_specific_hurdle_count_model_diagnostics.csv`
- `tables/wave_specific_mixing_predictor_hurdle_count_model_results.csv`
- `tables/wave_specific_mixing_predictor_hurdle_count_model_diagnostics.csv`
- `tables/wave_cluster_outcome_descriptives.csv`
- `tables/observed_expected_mixing_matrices.csv`

Figures:

- `figures/hurdle_count_effects.png`
- `figures/hurdle_count_effects.pdf`
- `figures/mixing_predictor_hurdle_count_effects.png`
- `figures/mixing_predictor_hurdle_count_effects.pdf`
- `figures/mixing_effects.png`
- `figures/mixing_effects.pdf`

Cache:

- `cache/cluster_table.parquet`
- `cache/domain_wave_cluster_table.parquet`

## Domain And Wave Extensions

`domain_wave_analysis.py` refits the SIMD-domain and wave-specific analyses
under the same main formulation: primary resolution 0.3, good-QC genomes, rare
lineage pooling at 50 clusters, 8 df calendar spline, and window-clustered
standard errors.

Domain count outcomes use the same hurdle/ZTNB structure as the main count
models. Domain-quintile mixing, domain-demographic mixing, and wave-specific
domain-demographic mixing use observed-minus-expected pairwise discordance
linear models with cluster size included. BA.4, BQ.1, and XBB are retained in
descriptives and matrices but skipped in the wave-specific regression table when
they fall below the minimum non-singleton cluster/window threshold.

`wave_outcome_analysis.py` also refits the main hurdle/ZTNB cluster-outcome
models separately by wave for cluster size and geographic spread. These
wave-stratified outcome models retain the main calendar, incidence,
sequencing, test-positivity, and SIMD-deprivation covariates. Within-wave
lineage dummies are included where estimable, with rank-redundant columns
dropped after the substantive covariates and calendar spline terms are kept.
XBB is retained in descriptives but skipped in the regression table because it
falls below the minimum cluster count.

## Diagnostic Note

The positive cluster-size ZTNB component converges, but the dispersion parameter
hits the upper bound. This is a useful warning about the extreme right tail of
cluster size. Interpret positive cluster-size ratios as a heavy-tail sensitivity,
not as a perfectly well-behaved parametric count model.
