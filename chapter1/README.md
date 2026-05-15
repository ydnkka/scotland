# Chapter 1

Are larger or more geographically dispersed SARS-CoV-2 genomic clusters
characterised by sociodemographic mixing beyond what would be expected
given their size and the population composition of the lineage and
analysis window from which they were drawn?

## Question

Cluster *scale* — `cluster_size` and `cluster_n_datazones` — is the
outcome.  Per-cluster *excess sociodemographic mixing* (observed minus
expected pairwise discordance, computed within lineage × analysis-window
strata) is the exposure of interest.  The hypothesis is that larger or
more dispersed clusters, plausibly reflecting superspreading or rapid
transmission expansion, are more mixed across age, sex, and area
deprivation than a random-assembly null would predict.

## Analysis population

- Input: `data/processed/scotland_clustering_analysis_dataset.parquet`
- QC filter: `nextclade_qc == "good"`
- Primary Leiden resolution: `0.3`
- Unit of analysis: one inferred cluster at the primary resolution
- Excess-mixing predictors are only defined for clusters with at least
  two valid cases for the variable in question, so every fit is on the
  non-singleton sub-population.  Both outcomes — `cluster_size` and
  `cluster_n_datazones` — are modelled as **zero-truncated negative
  binomial only** and no hurdle component is fit.  `cluster_size` uses
  `cluster_size − 1`; geographic spread uses the raw unique datazone
  count because some non-singleton clusters still occupy one unique
  datazone.
- Pango lineages are pooled if they have fewer than 30 non-singleton
  clusters.

## Models

The same numerical pipeline is reused across all fits.  Each fit is a
ZTNB on the non-singleton positive count, with cluster-robust standard
errors clustered by analysis window.

| Fit | Outcome | Excess-mixing predictors | Notes |
|---|---|---|---|
| Main effects | size · spread | age, sex, SIMD | lineage-adjusted |
| Size-adjusted spread | spread | age, sex, SIMD | + `log_cluster_size_z` |
| Wave interactions | size · spread | age, sex, SIMD, all × wave | wave dummies replace lineage |
| Size-spline sensitivity | spread | age, sex, SIMD | 4-df B-spline on `log(size)` |
| SIMD-decile sensitivity | size · spread | age, sex, SIMD decile | sensitivity to SIMD bucketing |
| Finite-sample mixing sensitivity | size · spread | finite-sample standardised age, sex, SIMD | sensitivity to pair-count precision |
| Joint-profile adjusted sensitivity | size · spread | age, sex, SIMD, joint age × sex × SIMD profile | predictor-set sensitivity |
| Null-residual sensitivity | size · spread | residuals from `obs ∼ size + entropy + lineage + window` | sensitivity to the random-assembly null |
| Demographic profile | size · spread | joint age × sex profile | single-predictor supplement |
| Sociodemographic profile | size · spread | joint age × sex × SIMD profile | single-predictor supplement |
| Domain stratification | size · spread | age, sex, domain-specific quintile mixing | one model per SIMD domain |
| Wave stratification | size · spread | age, sex, SIMD | one model per wave |

All adjustment covariates (`deprivation_z`, `local_incidence_z`,
`local_seq_fraction_z`, `window_seq_fraction_z`, `test_positivity_z`)
are z-scored and an 8-df B-spline on `window_idx` captures calendar
time.

## Running

```
conda run -n PhD python chapter1/overall_analysis.py
conda run -n PhD python chapter1/domain_analysis.py
conda run -n PhD python chapter1/wave_analysis.py
```

Add `--sample-clusters N` to `overall_analysis.py` to fit on a
subsample (useful while iterating).

The main sensitivity runs can be generated with:

```
bash chapter1/run_chapter1.sh
```

The overall-analysis script also accepts `--window-stride 3`,
`--winsorise-quantile 0.99`, and `--exclude-tail-quantile 0.995` for
targeted runs into separate output directories.

## Outputs

Tables (`chapter1/tables/`):

- `dataset_descriptives.csv`
- `covariate_scaling.csv`
- `main_effects_results.csv` / `main_effects_diagnostics.csv`
- `size_adjusted_spread_results.csv` / `size_adjusted_spread_diagnostics.csv`
- `wave_interaction_results.csv` / `wave_interaction_diagnostics.csv`
- `size_spline_sensitivity_results.csv` / `size_spline_sensitivity_diagnostics.csv`
- `simd_decile_sensitivity_results.csv` / `simd_decile_sensitivity_diagnostics.csv`
- `finite_sample_mixing_sensitivity_results.csv` / `finite_sample_mixing_sensitivity_diagnostics.csv`
- `joint_profile_adjusted_results.csv` / `joint_profile_adjusted_diagnostics.csv`
- `null_residual_sensitivity_results.csv` / `null_residual_sensitivity_diagnostics.csv`
- `profile_predictor_results.csv` / `profile_predictor_diagnostics.csv`
- `domain_main_effects_results.csv` / `domain_main_effects_diagnostics.csv`
- `wave_stratified_results.csv` / `wave_stratified_diagnostics.csv`

Figures (`chapter1/figures/`):

- `main_effects_forest.{png,pdf}`
- `wave_interaction_slopes.{png,pdf}`
- `domain_forest.{png,pdf}`
- `wave_stratified_forest.{png,pdf}`

Cache (`chapter1/cache/`):

- `cluster_table.parquet`
- `domain_cluster_table.parquet`
- `wave_cluster_table.parquet`

## Differences from `part1/`

`part1/` interleaves "Line 1" (deprivation → outcomes) and "Line 2"
(mixing → outcomes) in one orchestrator, with several outcome variants
and sensitivities woven together.  Chapter 1 keeps only the Line-2
framing — *mixing as exposure, scale as outcome* — and the question
this chapter answers is whether boundary-crossing transmission scales
super-linearly relative to a within-lineage-window null.

See `analysis_documentation.md` for the full model specification and
interpretation guide.
