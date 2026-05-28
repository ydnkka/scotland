# Candidate SSE Socio-Geodemographic Association Notes

## Purpose

This analysis tests whether socio-geodemographic variables are associated with candidate superspreading-signature node membership. Candidate status is defined upstream in `sse_detection.ipynb` from temporal cluster-transition graph features. The association analysis asks a downstream descriptive question: after restricting the background to comparable cluster sizes, do candidate nodes differ in sequence composition or internal population mixing?

The main notebook is intentionally light. It configures and runs `sse_detection.lib.association_pipeline`, then displays the exported summary tables. The modelling logic lives in the reusable pipeline so the main analysis and sensitivity notebooks use the same preparation, fitting, adjustment, and export code.

## Inputs

- `sse_detection/results/sse_outputs/node_stats.parquet`: node-level candidate labels, cluster summaries, entropy metrics, window metadata, and variant labels.
- `data/processed/scotland_clustering_analysis_dataset.parquet`: sequence-window rows used for the composition models.
- `utils.data.CLADES`: curated top-clade mapping used by clade sensitivity analyses.
- `sse_detection/lib/regression.py`: Firth logistic, GLM, conditional-logit, Wald, odds-ratio, and fit-stat helpers.
- `sse_detection/lib/association_pipeline.py`: shared data-preparation, model-running, and CSV-export orchestration.

## Outputs

The main analysis saves CSV tables to `sse_detection/results/association_outputs`:

- `composition_wald.csv`
- `composition_odds_ratios.csv`
- `composition_fit_stats.csv`
- `mixing_wald.csv`
- `mixing_odds_ratios.csv`
- `mixing_fit_stats.csv`
- `cluster_diagnostics.csv`
- `model_failures.csv`, only when one or more fits fail

Sensitivity notebooks write equivalent tables to their own output directories, such as `sse_detection/results/sensitivity_clade` and `sse_detection/results/sensitivity_window`.

The policy-era context analysis saves CSV tables to `sse_detection/results/policy_outputs`:

- `policy_era_candidate_summary.csv`
- `policy_era_category_summary.csv`
- `policy_wald.csv`
- `policy_odds_ratios.csv`
- `policy_fit_stats.csv`

The vaccination-context analysis saves CSV tables to `sse_detection/results/vaccination_outputs`:

- `vaccination_candidate_summary.csv`
- `vaccination_composition_wald.csv`
- `vaccination_composition_odds_ratios.csv`
- `vaccination_composition_fit_stats.csv`
- `vaccination_node_wald.csv`
- `vaccination_node_odds_ratios.csv`
- `vaccination_node_fit_stats.csv`
- `vaccination_mixing_age_conditional_node_features.csv`
- `vaccination_mixing_age_conditional_summary.csv`
- `vaccination_mixing_age_conditional_category_summary.csv`
- `vaccination_mixing_age_conditional_wald.csv`
- `vaccination_mixing_age_conditional_odds_ratios.csv`
- `vaccination_mixing_age_conditional_fit_stats.csv`

## Eligible Analysis Set

The background set is restricted to nodes with `cluster_size >= min(candidate cluster size)`. This prevents candidate nodes from being compared against all small singleton-like background nodes. The candidate label is stored as `candidate`, a binary version of `sse_candidate`.

For composition models, sequence-window rows are loaded with `utils.data.load_analysis_columns(..., window_stride=2)`, matching the odd-window retention used by the detection workflow. These rows are joined to eligible node status by `cluster_id`.

For mixing models, the node-level model frame is built directly from eligible `node_stats` rows.

## Model Families

### Composition Models

Composition models use sequence-window rows joined to node-level candidate status. They test whether sequences in candidate nodes differ from eligible background nodes by:

- sex
- age band
- SIMD quintile
- urban/rural class
- health board

Each composition predictor is fitted in a single-predictor model and then all composition predictors are fitted jointly. Odds ratios are sequence-level associations with candidate-node membership.

### Node-Level Mixing Models

Mixing models use node-level entropy features from `node_stats`. The default main analysis uses entropy null-model z-scores:

- `sex_entropy_z`
- `age_entropy_z`
- `simd_entropy_z`
- `urban_rural_entropy_z`
- `health_board_entropy_z`

These coefficients describe differences in candidate odds per one null-model SD increase in entropy. In practice, the fitted feature is a z-score:

`entropy_z = (observed entropy - null mean entropy) / null SD entropy`

where the null mean and SD are estimated for nodes of the same size within the same window. A higher z-score means the node is more internally mixed than expected for its size and window under the null model. The observed-entropy sensitivity notebook instead uses `*_entropy_obs_x10` features, so odds ratios are per 0.1 increase in observed normalised entropy.

Do not adjust the main mixing models for candidate-defining graph quantities such as `cluster_size`, `local_amplification_score`, `out_strength`, or `onward_dissemination_score`, because that conditions on variables used to define the outcome.

## Adjustment Sets

The main primary adjustment set is:

- `C(window_idx)`
- `C(clade)` by default

The expanded adjustment set adds standardised local surveillance and epidemic-burden context:

- `z_dz_cum_prop_sequenced`
- `z_dz_cum_incidence_per_capita`
- `z_dz_7d_test_positivity`
- `z_log1p_dz_cum_positive_tests`

The window-sensitivity analysis does not include `C(window_idx)`. It uses alternative window-level surveillance adjusters instead:

- `z_wn_prop_sequenced`
- `z_log1p_wn_positive_tests`

Do not combine the window-level surveillance adjusters with `C(window_idx)` in the same model. They are constants within window strata, so this creates rank collinearity and can destabilise fitting.

The clade-sensitivity notebooks fit models within `clade_group`, using `utils.data.CLADES` plus `Other`, and deliberately remove clade from the adjustment set.

## Fitting Method

The default method is Firth-penalised logistic regression. Exact conditional logistic regression stratified by `window_idx` was considered, but the sequence-level window strata are too large for the recursive exact-likelihood implementation in `statsmodels`. Firth penalisation reduces sparse-data and separation bias while retaining explicit adjustment terms.

The pipeline also supports:

- `glm_clustered`: binomial GLM with cluster-robust standard errors.
- `conditional_logit_by_window`: conditional logistic model stratified by `window_idx`, mainly useful for smaller model frames.

## Filtering and Dropped Strata

Before fitting, the pipeline performs complete-case filtering for each model's outcome, predictors, and adjustment terms. For fixed-effect window models, windows with no candidate/background variation are dropped. Dropped row and stratum counts are carried into the exported tables through:

- `dropped_nonvarying_rows`
- `dropped_nonvarying_strata`
- `dropped_nonvarying_detail`

This makes sparse-window filtering visible in every result table.

## Multiple Testing

Benjamini-Hochberg correction is applied within each `domain` by `model_set` by `predictor_set` family. For grouped sensitivity analyses, adjustment is additionally performed within the grouping column. For vaccination joint blocks, correction is also kept within `joint_model` so each pre-specified block remains a separate robustness family.

The main confirmatory family is the single-predictor primary model set. Expanded and joint models are robustness or sensitivity families and should be interpreted separately rather than pooled with the primary confirmatory family.

## Interpretation Order

Use the tables in this order:

1. Single-predictor primary models: main screening family for each socio-geodemographic variable, adjusted for window and variant context.
2. Single-predictor expanded models: sensitivity to local surveillance and epidemic-burden adjustment.
3. Joint primary models: whether each predictor retains signal after mutual adjustment for the other socio-geodemographic predictors.
4. Joint expanded models: the most conditional specification, useful as a robustness check rather than the simplest effect summary.
5. McFadden pseudo-R2 and likelihood summaries: diagnostics for comparing models within the same outcome family and analysis frame.

Composition models describe who and where the sequences in candidate nodes come from. Node-level mixing models describe whether candidate nodes are unusually internally diverse or concentrated for their size and window.

## Policy Context Analysis

Policy era is treated as a contextual association rather than a formal causal effect. The analysis keeps the existing candidate definition and eligible-node restriction, then adds grouped policy eras derived from `policy_period`.

The policy model does not include `C(window_idx)` because policy periods are calendar-time variables and would be absorbed by window fixed effects. Its primary model adjusts for `C(clade_group)` plus window-level surveillance intensity; its expanded model adds local datazone surveillance/incidence context. `policy_intensity` is fitted as a single-exposure sensitivity in both policy adjustment sets rather than jointly with `policy_era`.

## Vaccination Context Analysis

Vaccination is treated as a contextual association rather than a formal causal effect. The analysis adds sequence-level vaccination predictors (`is_vaccinated`, dose category, booster status, days-since-vaccination category, and datazone vaccination-event coverage), node-level vaccination summaries, and age-conditional vaccination-mixing features aggregated across sequences within each eligible node.

Vaccination composition, node, and age-conditional mixing models retain the existing `C(window_idx) + C(clade)` primary adjustment because vaccination varies within retained windows. Their expanded models add non-vaccination context:

- composition expanded models add `C(age_band)`, `C(sex)`, `C(dz_simd_quintile)`, `C(dz_urban_rural_class)`, and `C(dz_health_board)`;
- node expanded models add standardised datazone surveillance/incidence context.

Vaccination outputs include single-exposure models plus selected joint blocks. Composition joint blocks pair datazone vaccination-event coverage with vaccination status, dose category, booster status, or days-since-vaccination category. Node joint models are deliberately conservative and pair node proportion vaccinated with node mean datazone vaccination-event coverage. Joint models are robustness checks, not replacements for the primary single-exposure summaries.
