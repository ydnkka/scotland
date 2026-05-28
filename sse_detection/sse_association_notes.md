# Candidate SSE Association Notes

## Purpose

This document covers the association analyses downstream of candidate
superspreading-signature detection. Candidate status is defined upstream in
`sse_detection.ipynb` from temporal cluster-transition graph features. The
association analyses ask a descriptive question: after restricting the
background to comparable cluster sizes, are candidate nodes associated with
differences in sequence composition, geographic context, internal population
mixing, policy era, or vaccination context?

These models are not causal effect estimates. They describe associations with
candidate-node membership under pre-specified adjustment sets and sensitivity
analyses.

The main notebooks stay intentionally light. They configure and run
`sse_detection.lib.association_pipeline`, then display the exported summary
tables. The preparation, fitting, multiple-testing adjustment, diagnostics, and
CSV export logic live in the reusable pipeline so the main analysis and
sensitivity/context analyses use the same implementation.

## Association Domains

| Domain | Unit | Main exposures | Main question |
| --- | --- | --- | --- |
| Sociodemographic composition | sequence-window rows | sex, age band, SIMD quintile | Do sequences inside candidate nodes come from different individual or area deprivation groups than eligible background nodes? |
| Geographic composition | sequence-window rows | urban/rural class, health board | Are candidate-node sequences geographically concentrated in different settlement or health-board contexts? |
| Local surveillance and epidemic context | sequence or node rows | sequencing proportion, incidence, positivity, positive tests | Are associations robust after accounting for local testing, sequencing, and epidemic burden? |
| Node-level internal mixing | eligible nodes | sex, age, SIMD, urban/rural, and health-board entropy | Are candidate nodes more internally mixed or more homogeneous than expected for their size and window? |
| Policy context | eligible nodes | grouped policy era, policy intensity | Are candidate nodes more frequent in particular policy periods after variant and surveillance adjustment? |
| Vaccination composition | sequence-window rows | vaccination status, dose category, booster status, vaccination recency, datazone vaccination-event coverage | Do candidate-node sequences have different vaccination context than eligible background sequences? |
| Vaccination node context | eligible nodes | node vaccination proportions, dose summaries, recency summaries, area vaccination coverage | Do candidate nodes differ in aggregate vaccination context? |
| Age-conditional vaccination mixing | eligible nodes | vaccination-mixing entropy within age strata | Are candidate nodes unusually mixed by vaccination status after conditioning on age structure? |

## Inputs

- `sse_detection/results/sse_outputs/node_stats.parquet`: node-level candidate
  labels, cluster summaries, entropy metrics, window metadata, variant labels,
  policy-period metadata, and local context summaries.
- `data/processed/scotland_clustering_analysis_dataset.parquet`:
  sequence-window rows used for composition, geographic, policy-adjacent, and
  vaccination-context models.
- `utils.data.CLADES`: curated top-clade mapping used by clade sensitivity
  analyses.
- `utils.data.PERIOD_INTENSITY`: numeric policy-intensity scores used by the
  policy-era sensitivity model.
- `sse_detection/lib/regression.py`: Firth logistic, GLM,
  conditional-logit, Wald, odds-ratio, fit-statistic, and multiple-testing
  helpers.
- `sse_detection/lib/association_pipeline.py`: shared data-preparation,
  model-running, diagnostic, and CSV-export orchestration.

## Outputs

The primary demographic, geographic, and node-mixing analysis saves CSV tables
to `sse_detection/results/association_outputs`:

- `composition_wald.csv`
- `composition_odds_ratios.csv`
- `composition_fit_stats.csv`
- `mixing_wald.csv`
- `mixing_odds_ratios.csv`
- `mixing_fit_stats.csv`
- `cluster_diagnostics.csv`
- `model_failures.csv`, only when one or more fits fail

Sensitivity analyses write equivalent association tables to their own output
directories:

- `sse_detection/results/sensitivity_clade`
- `sse_detection/results/sensitivity_window`
- `sse_detection/results/sensitivity_observed_entropy`

The policy-context analysis saves CSV tables to
`sse_detection/results/policy_outputs`:

- `policy_era_candidate_summary.csv`
- `policy_era_category_summary.csv`
- `policy_wald.csv`
- `policy_odds_ratios.csv`
- `policy_fit_stats.csv`
- `cluster_diagnostics.csv`
- `model_failures.csv`, only when one or more fits fail

The vaccination-context analysis saves CSV tables to
`sse_detection/results/vaccination_outputs`:

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
- `cluster_diagnostics.csv`
- `model_failures.csv`, only when one or more fits fail

## Eligible Analysis Set

All association analyses use the same candidate outcome and comparable
background rule. The background set is restricted to nodes with
`cluster_size >= min(candidate cluster size)`. This prevents candidate nodes
from being compared against all small singleton-like background nodes. The
binary outcome is stored as `candidate`, derived from `sse_candidate`.

For sequence composition models, sequence-window rows are loaded with
`utils.data.load_analysis_columns(..., window_stride=2)`, matching the
odd-window retention used by the detection workflow. These rows are joined to
eligible node status by `cluster_id`.

For node-level models, the model frame is built directly from eligible
`node_stats` rows. Node-level models include the internal mixing analysis,
policy-context analysis, vaccination node-context analysis, and vaccination
mixing analysis.

## Model Families

### Sociodemographic And Geographic Composition

Composition models use sequence-window rows joined to node-level candidate
status. They test whether sequences in candidate nodes differ from eligible
background nodes by:

- sex
- age band
- SIMD quintile
- urban/rural class
- health board

Sex and age band are individual-level sociodemographic variables. SIMD quintile
is an area deprivation variable. Urban/rural class and health board are the
primary geographic composition variables. The resulting odds ratios are
sequence-level associations with candidate-node membership, not individual
risk estimates.

Each composition predictor is fitted in a single-predictor model and then all
composition predictors are fitted jointly. Single-predictor models are the main
screening family; joint models ask whether a predictor retains signal after
mutual adjustment for the other sociodemographic and geographic predictors.

### Local Surveillance And Epidemic Context

The expanded adjustment sets add datazone-level surveillance and epidemic
burden context:

- `z_dz_cum_prop_sequenced`
- `z_dz_cum_incidence_per_capita`
- `z_dz_7d_test_positivity`
- `z_log1p_dz_cum_positive_tests`

These are used as context adjusters rather than primary exposures in the main
composition and mixing models. They help assess whether associations with
candidate status are robust to local sequencing intensity, infection burden,
test positivity, and positive-test volume.

The window-sensitivity analysis replaces fixed window effects with
window-level surveillance adjusters:

- `z_wn_prop_sequenced`
- `z_log1p_wn_positive_tests`

Do not combine the window-level surveillance adjusters with `C(window_idx)` in
the same model. They are constants within window strata, so this creates rank
collinearity and can destabilise fitting.

### Node-Level Internal Mixing

Node-level mixing models use entropy features from `node_stats`. The default
main analysis uses null-model z-scores:

- `sex_entropy_z`
- `age_entropy_z`
- `simd_entropy_z`
- `urban_rural_entropy_z`
- `health_board_entropy_z`

These coefficients describe differences in candidate odds per one null-model
standard deviation increase in entropy. The fitted feature is:

`entropy_z = (observed entropy - null mean entropy) / null SD entropy`

where the null mean and SD are estimated for nodes of the same size within the
same window. A higher z-score means the node is more internally mixed than
expected for its size and window under the null model. A lower z-score means
the node is more homogeneous than expected.

The observed-entropy sensitivity analysis instead uses
`*_entropy_obs_x10` features, so odds ratios are per 0.1 increase in observed
normalised entropy.

Do not adjust the main mixing models for candidate-defining graph quantities
such as `cluster_size`, `local_amplification_score`, `out_strength`, or
`onward_dissemination_score`, because that conditions on variables used to
define the outcome.

### Policy Context

Policy era is treated as contextual association, not a formal causal effect of
policy. The policy analysis keeps the existing candidate definition and
eligible-node restriction, then adds grouped policy eras derived from
`policy_period`.

Policy eras are ordered as:

- `early_restriction_easing`
- `autumn_winter_restrictions`
- `spring_summer_2021_easing`
- `near_normal_delta`
- `omicron_response`
- `post_restriction`

The primary policy model uses `post_restriction` as the reference era and
adjusts for `C(clade_group)` plus window-level surveillance intensity. It does
not include `C(window_idx)` because policy periods are calendar-time variables
and would be absorbed by window fixed effects. The expanded policy model adds
the datazone surveillance and epidemic-burden adjusters.

`policy_intensity` is fitted as a single-exposure sensitivity in the primary
and expanded policy adjustment sets. It is not fitted jointly with
`policy_era`, because both describe the same underlying policy-period context
at different resolutions.

### Vaccination Composition

Vaccination composition models use sequence-window rows joined to candidate
status. They test:

- individual vaccination status, `is_vaccinated`
- prior dose category, `vacc_dose_cat`
- prior booster status, `vacc_booster_status`
- days since prior vaccination, `days_since_vaccination_cat`
- datazone cumulative vaccination-event coverage, `z_dz_cum_prop_vaccinated`

Primary vaccination composition models adjust for `C(window_idx) + C(clade)`.
Expanded vaccination composition models add non-vaccination composition
context:

- `C(age_band)`
- `C(sex)`
- `C(dz_simd_quintile)`
- `C(dz_urban_rural_class)`
- `C(dz_health_board)`

Vaccination composition joint blocks pair datazone vaccination-event coverage
with vaccination status, dose category, booster status, or vaccination recency.
These joint blocks separate individual vaccination context from area-level
vaccination-event coverage as robustness checks.

### Vaccination Node Context

Vaccination node models aggregate sequence-level vaccination context onto each
eligible node and fit node-level candidate associations. The node exposures
are:

- node proportion vaccinated
- node proportion boosted
- node mean prior dose count
- node median days since prior vaccination
- node mean datazone vaccination-event coverage

Continuous node vaccination features are standardised, so odds ratios are
reported per one standard deviation increase unless a feature-specific
reference says otherwise. Primary models adjust for `C(window_idx) + C(clade)`.
Expanded node models add the datazone surveillance and epidemic-burden
adjusters. The pre-specified joint node block pairs node proportion vaccinated
with node mean datazone vaccination-event coverage.

### Age-Conditional Vaccination Mixing

Age-conditional vaccination-mixing models test whether candidate nodes are
unusually mixed by vaccination status after conditioning on age structure.
The pipeline computes a binary vaccination entropy feature within age-window
null structure and exports both the continuous entropy z-score and an ordered
tertile:

- `vaccination_mix_entropy_z`
- `vaccination_mix_tertile`

The ordered tertile levels are:

- `more_homogeneous`
- `as_expected`
- `more_mixed`

The primary vaccination-mixing model adjusts for `C(window_idx) + C(clade)`.
The `age_mixing` model additionally adjusts for `age_entropy_z`. The expanded
model adds `age_entropy_z` plus datazone surveillance and epidemic-burden
context. This separates vaccination mixing from general age mixing and local
testing or incidence context.

## Adjustment Sets

The main primary adjustment set is:

- `C(window_idx)`
- `C(clade)`

The main expanded adjustment set adds:

- `z_dz_cum_prop_sequenced`
- `z_dz_cum_incidence_per_capita`
- `z_dz_7d_test_positivity`
- `z_log1p_dz_cum_positive_tests`

The clade-sensitivity analyses fit models within `clade_group`, using
`utils.data.CLADES` plus `Other`, and deliberately remove clade from the
adjustment set.

The window-sensitivity analysis does not include `C(window_idx)`. It uses
window-level surveillance adjustment instead, which is useful when the window
fixed effects are too restrictive for a given interpretation.

## Fitting Method

The default method is Firth-penalised logistic regression. Exact conditional
logistic regression stratified by `window_idx` was considered, but the
sequence-level window strata are too large for the recursive exact-likelihood
implementation in `statsmodels`. Firth penalisation reduces sparse-data and
separation bias while retaining explicit adjustment terms.

The pipeline also supports:

- `glm_clustered`: binomial GLM with cluster-robust standard errors.
- `conditional_logit_by_window`: conditional logistic model stratified by
  `window_idx`, mainly useful for smaller model frames.

## Filtering And Dropped Strata

Before fitting, the pipeline performs complete-case filtering for each model's
outcome, predictors, and adjustment terms. For fixed-effect window models,
windows with no candidate/background variation are dropped. Dropped row and
stratum counts are carried into the exported tables through:

- `dropped_nonvarying_rows`
- `dropped_nonvarying_strata`
- `dropped_nonvarying_detail`

This makes sparse-window filtering visible in every result table.

## Multiple Testing

Benjamini-Hochberg correction is applied within each `domain` by `model_set`
by `predictor_set` family. For grouped sensitivity analyses, adjustment is
additionally performed within the grouping column. For vaccination joint
blocks, correction is also kept within `joint_model` so each pre-specified
block remains a separate robustness family.

The main confirmatory family is the single-predictor primary model set.
Expanded, joint, policy, vaccination, and sensitivity models are robustness or
context families and should be interpreted separately rather than pooled with
the primary screening family.

## Interpretation Order

Use the tables in this order:

1. Main single-predictor primary models: first-pass demographic, geographic,
   and mixing screening family, adjusted for window and clade context.
2. Main single-predictor expanded models: robustness to local surveillance and
   epidemic-burden adjustment.
3. Main joint primary models: whether each demographic, geographic, or mixing
   predictor retains signal after mutual adjustment within its model family.
4. Main joint expanded models: the most conditional main specification, useful
   as a robustness check rather than the simplest effect summary.
5. Clade, window, and observed-entropy sensitivity outputs: whether the main
   results depend on variant grouping, window adjustment strategy, or entropy
   scale.
6. Policy outputs: contextual timing patterns by grouped policy era and policy
   intensity.
7. Vaccination outputs: sequence, node, and age-conditional vaccination-mixing
   context.
8. McFadden pseudo-R2 and likelihood summaries: diagnostics for comparing
   models within the same outcome family and analysis frame.

Composition models describe who and where the sequences in candidate nodes
come from. Node-level mixing models describe whether candidate nodes are
unusually internally diverse or concentrated for their size and window. Policy
and vaccination models describe broader context around candidate-node
membership and should be reported as contextual associations.

## Reporting Guardrails

- Report odds ratios as associations with candidate-node membership, not as
  individual infection risk or causal policy/vaccination effects.
- Keep sequence-level and node-level interpretations separate.
- Mention the eligible-node restriction whenever comparing candidate and
  background nodes.
- Use primary single-predictor results for main screening statements and use
  expanded, joint, and sensitivity models to qualify robustness.
- Treat geographic variables carefully: health board and urban/rural class are
  composition/context variables, while datazone surveillance and epidemic
  variables are adjustment context.
- Do not infer that vaccination or policy caused candidate status without a
  causal design beyond these descriptive models.
