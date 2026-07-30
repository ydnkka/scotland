# Bayesian SSE Regression Models

This workflow characterises detector outcomes with sociogeodemographic mixing or sequence composition. It estimates associations, not causal effects.

## Code and inputs

- `lib/model/prep.py`: eligibility, predictors, complete-case/sampled frames, formulas, and paths.
- `lib/model/bayesian.py`: Bambi/PyMC models, priors, summaries, diagnostics, and output writing.
- `lib/model/runner.py`: CLI selection, fitting, logs, locks, and manifests.
- `lib/model/mixing.py` and `composition.py`: domain entry points.

The loader requires both detector files:

```text
results/sse_outputs/cluster_table.parquet
results/sse_outputs/edge_table.parquet
```

Model preparation uses `cluster_table`; `edge_table` is loaded as part of the `SseOutputs` contract but is not a formula input. Sequence-level composition rows are reloaded with the detector's window stride and joined to cluster outcomes.

High-priority burst, burden, or both-axis tiers define `candidate = 1`. Eligible nodes have `cluster_size` at least the smallest observed high-priority candidate size; this can be greater than the detector's six-sequence testing floor. The workflow fails if no high-priority candidates exist.

Two frames are produced:

- `eligible_nodes`: one cluster per row, used for mixing models;
- `eligible_sequence_data`: sequence-window rows inheriting cluster outcomes, used for composition models.

## Models

All formulas include varying intercepts for `policy_period` and `clade`:

```text
outcome ~ focal terms + optional context terms + (1|policy_period) + (1|clade)
```

| Family   | Likelihood | Outcome                         |
| -------- | ---------- | ------------------------------- |
| Logistic | Bernoulli  | `candidate`                     |
| Linear   | Gaussian   | `burst_score` or `burden_score` |

Rows missing the outcome, any selected predictor, or either grouping variable are removed separately for each specification. Burden models therefore retain only rows with a non-missing burden score.

### Mixing domain

Node-level model sets:

| Set                 | Fixed effects                                                                     |
| ------------------- | --------------------------------------------------------------------------------- |
| `null_primary`      | Null-standardised entropy for sex, age, SIMD, urban/rural class, and Health Board |
| `null_expanded`     | `null_primary` plus context adjusters                                             |
| `observed_primary`  | Observed entropy for the same five attributes, scaled ×10                         |
| `observed_expanded` | `observed_primary` plus context adjusters                                         |

Data Zone and local-authority entropy are computed by the detector but explicitly excluded from the current saved mixing specifications.

### Composition domain

Sequence-level model sets:

| Set        | Fixed effects                                                      |
| ---------- | ------------------------------------------------------------------ |
| `primary`  | Sex, age group, SIMD quintile, urban/rural class, and Health Board |
| `expanded` | `primary` plus context adjusters                                   |

Treatment-code references are defined in `COMPOSITION_SPECS`: Male, age 25–64, SIMD quintile 1, Large Urban Areas, and Greater Glasgow and Clyde.

Expanded models add standardised:

```text
wn_prop_sequenced_z
dz_cum_incidence_per_capita_z
dz_cum_prop_sequenced_z
```

## Commands

List selectors:

```bash
python -m chapter_analyses.sse_detection.lib.model.mixing --list-models
python -m chapter_analyses.sse_detection.lib.model.composition --list-models
```

Inspect formulas and fit-frame sizes without fitting:

```bash
python -m chapter_analyses.sse_detection.lib.model.mixing --dry-run --family logistic --model-set null_primary
python -m chapter_analyses.sse_detection.lib.model.composition --dry-run --family linear --outcome burden_score --model-set expanded
```

Fit selected models:

```bash
python -m chapter_analyses.sse_detection.lib.model.mixing --family logistic --model-set observed_expanded
python -m chapter_analyses.sse_detection.lib.model.composition --family linear --outcome burst_score --model-set primary
```

With no selectors, either domain command runs both families and every applicable outcome/model set.

Important options:

- selection: `--family`, `--outcome`, `--model-set`;
- sampling: `--sample-fraction` (default 1.0), `--sample-rows`, `--positive-fraction`, `--no-sample`;
- fitting: `--draws`, `--tune`, `--chains`, `--cores`, `--target-accept`, `--random-seed`, `--inference-method`;
- priors: `--fixed-prior-sigma`, `--intercept-prior-sigma`, `--random-effect-sigma`, `--residual-sigma`;
- execution: `--skip-existing`, `--continue-on-error`, `--save-idata`, `--live-progress`, `--centered`, and `--no-log-likelihood`.

Use `--help` for path, table-writing, diagnostic-display, log, and JAX options.

Logistic development samples preserve the requested positive fraction; when unspecified, the observed candidate rate is used. Composition frames can be much larger than node frames. Compare models only after checking their complete-case and sampling populations.

## Priors and sampling defaults

| Term               | Logistic                            | Linear                       |
| ------------------ | ----------------------------------- | ---------------------------- |
| Intercept          | `Normal(logit(observed rate), 1.5)` | `Normal(observed mean, 1.0)` |
| Common effects     | `Normal(0, 1.0)`                    | `Normal(0, 0.5)`             |
| Varying intercepts | `Normal(0, HalfNormal(1.0))`        | `Normal(0, HalfNormal(0.5))` |
| Residual SD        | —                                   | `HalfNormal(0.5)`            |

Defaults are 2,000 draws, 2,000 tuning iterations, 4 chains, 4 cores, target acceptance 0.99, PyMC inference, random seed 123, non-centred group effects, and log-likelihood calculation. Flags can override all listed scales/settings.

If a multi-core PyMC fit fails with the handled EOF/end-of-file condition, the fitter retries with one core. Formulas containing categorical terms use PyMC instead of `nutpie` to avoid coefficient-dimension conversion problems.

## Outputs

Default root:

```text
chapter_analyses/sse_detection/results/bayesian_outputs/
```

Root files include `eligibility_summary.csv`, `<domain>_selected_model_grid.csv`, `saved_model_manifest.csv`, `last_saved_model_manifest.csv`, and family-specific manifests.

Each `logistic/` or `linear/` directory contains `run_config.csv`, `model_grid.csv`, `fit_frame_summary.csv`, domain-specific grid/summary files, and a manifest. Linear paths include the outcome; logistic paths do not:

```text
logistic/<domain>/<model_set>/
linear/<domain>/<outcome>/<model_set>/
```

Consolidate saved summaries, metadata, fit-frame counts, and headline diagnostics into family/domain result tables with:

```bash
python -m chapter_analyses.sse_detection.combine_bayesian_results
```

The default output directory is `results/bayesian_outputs/consolidated_tables/`. It writes four CSV/parquet table pairs:

```text
mixing_logistic_consolidated_results
mixing_linear_consolidated_results
composition_logistic_consolidated_results
composition_linear_consolidated_results
```

Consolidated tables include `plot_*` columns for downstream figures. The single display field is `plot_label`, and it is populated for fixed effects, intercepts, policy-period/clade random intercepts, group SDs, and linear residual SDs. Composition model rows parse Patsy categorical terms with local label logic and include contrast references and panels. Mixing tables omit composition-only contrast columns such as `plot_level` and `plot_reference`; their `plot_scale` uses `null_standardised` instead of `null` so CSV readers do not treat it as missing. Sample columns are family-specific: logistic tables include candidate counts/rate, while linear tables include outcome mean/SD fields. The script also writes three thesis-facing summary tables: diagnostics, directional estimates, and directionless random-effect/residual SD components.

Each fitted model writes:

| File              | Content                                                              |
| ----------------- | -------------------------------------------------------------------- |
| `summary.csv`     | Posterior parameter, HDI, R-hat, ESS, scale, and direction summaries |
| `diagnostics.csv` | Divergence, BFMI, R-hat, ESS, and tree-depth checks                  |
| `metadata.csv`    | Formula, rows, outcome mean/rate, and sampling flag                  |
| `fit.log`         | Structured header, diagnostic report, saved paths, or traceback      |
| `idata.nc`        | Optional inference data with `--save-idata`                          |

Backend output is suppressed by default and not copied to `fit.log`; `--live-progress` shows it in the terminal. CSV, parquet, and NetCDF model artifacts use same-directory temporary files and atomic replacement.

Advisory locks protect shared prepared tables/manifests, and `<model_dir>/.fit.lock` prevents two processes fitting the same specification. A stale fit lock should be removed only after confirming no matching process is running.

## Interpretation

Logistic coefficients are log odds; the summary also provides odds ratios and posterior probabilities `P(OR > 1 | data)` and `P(OR < 1 | data)`. Linear coefficients are conditional expected differences in the relevant score with `P(beta > 0 | data)` and `P(beta < 0 | data)`.

Categorical coefficients compare a level with its configured treatment reference. Standardised continuous effects correspond to a one-SD predictor increase. These are conditional coefficients, not marginal contrasts.

Diagnostic defaults flag any divergence, minimum BFMI below 0.3, maximum R-hat above 1.01, or minimum bulk/tail ESS below 400. Tree depth is reported informationally. Do not interpret coefficients until fit-frame, convergence, and sampling differences have been reviewed.

## Reporting checklist

1. Confirm formula and output path in the model grid.
2. Check `fit_frame_summary.csv` and `metadata.csv` for complete-case and sampled rows.
3. Review `fit.log` and `diagnostics.csv`.
4. Interpret logistic effects as ORs and linear effects as score differences.
5. Compare primary/expanded models only on compatible fitted populations.
6. Report effect sizes, credible intervals, diagnostics, and the repeated cluster-outcome structure of sequence-level composition models.

To change specifications, edit constants and `default_model_specs(...)` in `model/prep.py`, dry-run the affected selectors, refit, then regenerate result tables and figures.
