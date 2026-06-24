# Bayesian SSE Regression Models: Technical Documentation <!-- omit in toc -->

This document describes the socio-geodemographic Bayesian regression workflow used for SSE candidate, burst-score, and burden-score analyses.

It covers:

- How the regression code is structured.
- Which data frames, outcomes, predictors, and model families are used.
- How to run one model version or all versions.
- How outputs are organized.
- How to interpret posterior summaries and diagnostics.
- Maintenance notes for changing or extending the models.

## Content <!-- omit in toc -->

<!-- TOC tocDepth:2..3 chapterDepth:2..6 -->

- [1. Workflow Structure](#1-workflow-structure)
- [2. Data Inputs](#2-data-inputs)
- [3. Model Domains](#3-model-domains)
  - [3.1. Mixing Models](#31-mixing-models)
  - [3.2. Composition Models](#32-composition-models)
- [4. Outcomes and Families](#4-outcomes-and-families)
- [5. Predictor Groups](#5-predictor-groups)
- [6. Model Specification Inspection](#6-model-specification-inspection)
- [7. Running Models](#7-running-models)
- [8. Sampling and Fit Frames](#8-sampling-and-fit-frames)
- [9. Output Structure](#9-output-structure)
- [10. Priors and Fitting Defaults](#10-priors-and-fitting-defaults)
  - [10.1. Logistic Models](#101-logistic-models)
  - [10.2. Linear Models](#102-linear-models)
- [11. Interpreting `summary.csv`](#11-interpreting-summarycsv)
  - [11.1. Logistic Models](#111-logistic-models)
  - [11.2. Linear Models](#112-linear-models)
- [12. Interpreting Categorical Composition Terms](#12-interpreting-categorical-composition-terms)
- [13. Interpreting Diagnostics](#13-interpreting-diagnostics)
- [14. Comparing Primary and Expanded Models](#14-comparing-primary-and-expanded-models)
- [15. Interpretation Caveats](#15-interpretation-caveats)
- [16. Quality-Control Checklist](#16-quality-control-checklist)
- [17. Maintenance Guide](#17-maintenance-guide)

<!-- /TOC -->

## 1. Workflow Structure

The regression workflow is split into three layers:

1. `sse_detection/lib/regression_prep.py`

   Prepares eligible node and sequence frames, defines predictor sets, builds formulas, and creates deterministic output paths.

2. `sse_detection/lib/bayesian_models.py`

   Fits Bambi/PyMC models, builds posterior summaries, computes diagnostics, and writes standard result files.

3. `sse_detection/lib/regression_runner.py`, `mixing_models.py`, and `composition_models.py`

   Provide repeatable command-line orchestration for running one model version or all model versions, with per-model logs.

## 2. Data Inputs

The regression workflow starts from:

```text
sse_detection/results/sse_outputs/cluster_table.parquet
sse_detection/results/sse_outputs/edge_table.parquet
```

These are loaded with:

```python
from sse_detection.lib import load_sse_outputs, prepare_regression_data

sse_outputs = load_sse_outputs("sse_detection/results/sse_outputs")
regression_data = prepare_regression_data(sse_outputs)
```

`prepare_regression_data(...)` creates two aligned modelling frames:

- `eligible_nodes`: one row per eligible cluster/node. Used by mixing models.
- `eligible_sequence_data`: one row per eligible sequence-window record. Used by composition models.

Eligibility is based on cluster size. Nodes are retained if `cluster_size` is at least the smallest high-priority candidate cluster size. Candidate status is defined by `candidate_tier` membership in the high-priority burst/burden tiers.

The preparation step also adds:

- observed mixing entropy scales where needed;
- standardized epidemic-context adjusters;
- candidate labels;
- cluster-level burst and burden score outcomes merged onto sequence-level rows.

## 3. Model Domains

### 3.1. Mixing Models

Mixing models are node-level regressions. Each row is an eligible cluster/node.

The focal predictors describe socio-geodemographic mixing entropy.

Available model sets:

| Model set           | Predictors                                                                |
| ------------------- | ------------------------------------------------------------------------- |
| `null_primary`      | Null-standardized mixing entropy features                                 |
| `null_expanded`     | Null-standardized mixing entropy features plus epidemic-context adjusters |
| `observed_primary`  | Observed mixing entropy scales                                            |
| `observed_expanded` | Observed mixing entropy scales plus epidemic-context adjusters            |

### 3.2. Composition Models

Composition models are sequence-level regressions. Each row is a sequence-window record that inherits candidate and score outcomes from its cluster.

The focal predictors describe the sequence's socio-geodemographic attributes.

Available model sets:

| Model set  | Predictors                                             |
| ---------- | ------------------------------------------------------ |
| `primary`  | Composition predictors only                            |
| `expanded` | Composition predictors plus epidemic-context adjusters |

Composition predictors are categorical and use treatment coding. Each coefficient is interpreted relative to the configured reference category.

## 4. Outcomes and Families

There are two regression families:

| Family   | Likelihood | Outcome(s)                    | Interpretation                         |
| -------- | ---------- | ----------------------------- | -------------------------------------- |
| Logistic | Bernoulli  | `candidate`                   | Candidate vs non-candidate association |
| Linear   | Gaussian   | `burst_score`, `burden_score` | Expected score difference              |

The model formula pattern is:

```text
outcome ~ focal_predictors + optional_context_adjusters + (1|policy_period) + (1|clade)
```

All models include varying intercepts for:

- `policy_period`
- `clade`

These group intercepts adjust for broad policy-period and variant/clade structure.

## 5. Predictor Groups

The definitive predictor definitions live in `sse_detection/lib/regression_prep.py`.

Key groups:

- `NULL_MIXING_FEATURES`: null-standardized mixing entropy features.
- `OBSERVED_MIXING_FEATURES`: observed entropy/mixing scales.
- `COMPOSITION_FEATURES`: sequence-level categorical socio-geodemographic predictors and reference categories.
- `EPIDEMIC_CONTEXT_ADJUSTERS`: surveillance and epidemic-context covariates.
- `GROUP_VARS`: varying-intercept grouping variables.

Current epidemic-context adjusters:

```text
wn_prop_sequenced_z
dz_cum_incidence_per_capita_z
dz_cum_prop_sequenced_z
```

Composition features currently include:

```text
sex
age_band
dz_simd_quintile
dz_urban_rural_class
dz_health_board
```

The references for categorical composition contrasts are defined in `COMPOSITION_SPECS` and exposed through `COMPOSITION_FEATURES`.

## 6. Model Specification Inspection

To inspect the live model specifications from Python:

```python
import pandas as pd
from sse_detection.lib import default_model_specs

def model_spec_table(domain: str) -> pd.DataFrame:
    rows = []
    for spec in default_model_specs(domain):
        rows.append(
            {
                "domain": spec.domain,
                "model_set": spec.model_set,
                "predictor_group": spec.predictor,
                "n_predictors": len(spec.predictors),
                "predictors": ", ".join(spec.predictors),
                "formula_terms": " + ".join(spec.terms),
            }
        )
    return pd.DataFrame(rows)

print(model_spec_table("mixing"))
print(model_spec_table("composition"))
```

## 7. Running Models

List available model selectors:

```bash
python -m sse_detection.lib.mixing_models --list-models
python -m sse_detection.lib.composition_models --list-models
```

Dry-run a selection without fitting Bambi/PyMC:

```bash
python -m sse_detection.lib.mixing_models --dry-run --family logistic --model-set null_primary
python -m sse_detection.lib.composition_models --dry-run --family linear --outcome burden_score --model-set expanded
```

Fit one model version:

```bash
python -m sse_detection.lib.mixing_models --family logistic --model-set observed_expanded
python -m sse_detection.lib.composition_models --family linear --outcome burst_score --model-set primary
```

Fit all models for one domain:

```bash
python -m sse_detection.lib.mixing_models
python -m sse_detection.lib.composition_models
```

Common options:

| Option                                     | Purpose                                   |
| ------------------------------------------ | ----------------------------------------- |
| `--family logistic`                        | Select regression family                  |
| `--model-set ...`                          | Select one or more model versions         |
| `--outcome candidate`                      | Select outcome                            |
| `--sample-fraction ...`                    | Fit a fraction of complete-case rows      |
| `--sample-rows ...`                        | Fit a fixed maximum number of rows        |
| `--skip-existing`                          | Avoid refitting models with saved outputs |
| `--continue-on-error`                      | Continue after a failed fit               |
| `--save-idata`                             | Save ArviZ inference data as `idata.nc`   |
| `--draws`, `--tune`, `--chains`, `--cores` | Sampler size and parallelism              |
| `--target-accept`                          | PyMC NUTS target acceptance probability   |

## 8. Sampling and Fit Frames

The scripts default to `--sample-fraction 1.0`, so the complete-case data are used unless another sampling option is specified.

For logistic models, sampling can preserve the observed positive fraction using `SampleSpec.positive_fraction`. If no positive fraction is supplied, the runner uses the observed candidate rate for the selected domain.

For composition models, the sequence-level frame can be much larger than the node-level frame. Use `--sample-rows` or `--sample-fraction` for development fits, then rerun final models on the intended full or sampled frame.

Always inspect:

```text
fit_frame_summary.csv
metadata.csv
```

before comparing estimates across models. Complete-case filtering and sampling can change the fitted population.

## 9. Output Structure

The default output root is:

```text
sse_detection/results/bayesian_socio_geo_demo_regression/
```

Run-level files:

| File                                  | Meaning                                                              |
| ------------------------------------- | -------------------------------------------------------------------- |
| `eligibility_summary.csv`             | Row counts and candidate rates for eligible node and sequence frames |
| `mixing_selected_model_grid.csv`      | Latest selected model grid for the mixing script                     |
| `composition_selected_model_grid.csv` | Latest selected model grid for the composition script                |
| `saved_model_manifest.csv`            | Accumulated manifest of saved, skipped, or failed model attempts     |
| `last_saved_model_manifest.csv`       | Manifest for the latest invocation only                              |
| `logistic/`                           | Logistic prepared-run tables and model outputs                       |
| `linear/`                             | Linear prepared-run tables and model outputs                         |

Family-level files under `logistic/` and `linear/`:

| File                             | Meaning                                                                   |
| -------------------------------- | ------------------------------------------------------------------------- |
| `run_config.csv`                 | Preparation settings                                                      |
| `model_grid.csv`                 | Prepared formulas and output directories                                  |
| `fit_frame_summary.csv`          | Complete-case rows, sampled rows, outcome rates/means, and sampling flags |
| `{domain}_model_grid.csv`        | Domain-specific prepared formulas                                         |
| `{domain}_fit_frame_summary.csv` | Domain-specific fit-frame summary                                         |
| `saved_model_manifest.csv`       | Family-specific saved model manifest                                      |

Per-model files:

| File              | Meaning                                                                           |
| ----------------- | --------------------------------------------------------------------------------- |
| `summary.csv`     | Posterior summaries for model parameters                                          |
| `diagnostics.csv` | Divergences, BFMI, R-hat, ESS, and tree-depth checks                              |
| `metadata.csv`    | Family, domain, outcome, formula, row count, outcome mean/rate, and sampling flag |
| `fit.log`         | Captured stdout/stderr for the model fit                                          |
| `idata.nc`        | Optional ArviZ inference data when `--save-idata` is used                         |

Typical path layout:

```text
sse_detection/results/bayesian_socio_geo_demo_regression/
  logistic/
    mixing/
      null_primary/
        summary.csv
        diagnostics.csv
        metadata.csv
        fit.log
      observed_expanded/
        ...
    composition/
      primary/
        ...
  linear/
    mixing/
      burst_score/
        null_primary/
          ...
      burden_score/
        observed_expanded/
          ...
    composition/
      burst_score/
        primary/
          ...
      burden_score/
        expanded/
          ...
```

## 10. Priors and Fitting Defaults

Priors are configured in `BayesianFitConfig` and applied by `fit_bayesian_model(...)`.

### 10.1. Logistic Models

- Bernoulli likelihood.
- Intercept prior centered at the logit of the observed outcome rate.
- Fixed effects use Normal priors centered at 0.
- Varying intercepts use Normal priors with HalfNormal SD priors.

### 10.2. Linear Models

- Gaussian likelihood.
- Intercept prior centered at the observed outcome mean.
- Fixed effects use Normal priors centered at 0.
- Varying intercepts use Normal priors with HalfNormal SD priors.
- Residual `sigma` uses a HalfNormal prior.

Default sampler settings:

```text
draws = 2000
tune = 2000
chains = 4
cores = 4
target_accept = 0.99
random_seed = 123
inference_method = pymc
```

Override these with script flags when needed.

## 11. Interpreting `summary.csv`

`summary.csv` is generated from ArviZ posterior summaries, with additional scale and direction summaries.

Common columns:

| Column                 | Meaning                                  |
| ---------------------- | ---------------------------------------- |
| `parameter`            | Posterior variable name                  |
| `mean`                 | Posterior mean                           |
| `sd`                   | Posterior standard deviation             |
| HDI columns            | Highest-density credible interval bounds |
| `r_hat`                | Chain convergence diagnostic             |
| `ess_bulk`, `ess_tail` | Effective sample-size diagnostics        |

### 11.1. Logistic Models

Logistic coefficients are on the log-odds scale.

The summary also includes odds-ratio columns, created by exponentiating posterior coefficients.

Interpretation:

- `OR_mean > 1`: higher candidate odds for larger predictor values, conditional on the other model terms.
- `OR_mean < 1`: lower candidate odds for larger predictor values.
- `P(OR > 1 | data)`: posterior probability that the association is positive on the odds scale.
- `P(OR < 1 | data)`: posterior probability that the association is negative on the odds scale.

For standardized continuous predictors ending in `_z`, the odds ratio corresponds to a one-standard-deviation predictor increase.

For categorical predictors, odds ratios compare each level to the configured reference level.

### 11.2. Linear Models

Linear coefficients are expected changes in `burst_score` or `burden_score`, conditional on the other model terms.

Interpretation:

- `mean > 0`: higher expected score.
- `mean < 0`: lower expected score.
- `P(beta > 0 | data)`: posterior probability of a positive association.
- `P(beta < 0 | data)`: posterior probability of a negative association.

For standardized predictors ending in `_z`, the coefficient is the expected outcome change for a one-standard-deviation predictor increase.

For categorical predictors, coefficients compare each level to the configured reference level.

## 12. Interpreting Categorical Composition Terms

Composition terms use treatment coding, for example:

```text
C(sex, Treatment(reference='Male'))
```

Each coefficient is a conditional contrast against the reference category.

Example interpretation:

```text
C(sex, Treatment(reference='Male'))[T.Female]
```

means the expected difference for `Female` relative to `Male`, conditional on the other predictors and varying intercepts.

These are not marginal population differences unless the model has been post-processed to estimate marginal contrasts.

## 13. Interpreting Diagnostics

The diagnostics table screens whether posterior summaries are trustworthy enough to interpret.

| Diagnostic     | Desired pattern                | Concern                                      |
| -------------- | ------------------------------ | -------------------------------------------- |
| Divergences    | Usually zero                   | Posterior geometry problems                  |
| BFMI           | Usually above about 0.3        | Poor energy exploration                      |
| Max R-hat      | Usually `<= 1.01`              | Chains have not mixed well                   |
| Bulk ESS       | Above threshold, currently 400 | Too few effective central posterior draws    |
| Tail ESS       | Above threshold, currently 400 | Too few effective tail posterior draws       |
| Max tree depth | Informational                  | Concerning if many transitions hit max depth |

If diagnostics are weak:

1. Inspect `fit.log` for sampler warnings.
2. Increase `--tune`.
3. Raise `--target-accept`.
4. Increase `--draws` or `--chains`.
5. Check sparse categories, quasi-separation, or highly collinear predictors.
6. Consider simplifying the model.

## 14. Comparing Primary and Expanded Models

Primary models estimate associations using only focal mixing or composition predictors plus group intercepts.

Expanded models add epidemic-context adjusters:

```text
wn_prop_sequenced_z
dz_cum_incidence_per_capita_z
dz_cum_prop_sequenced_z
```

Use the expanded model to assess whether focal associations are robust to these surveillance and epidemic-context covariates.

Changes between primary and expanded estimates can indicate confounding or mediation-like structure. They should not be treated as causal decomposition without a separate causal design.

## 15. Interpretation Caveats

- These are observational association models, not causal effect estimates.
- Sequence-level composition models repeat cluster-level outcomes across sequence-window rows, so rows are not independent biological outcomes in the same way as unique clusters.
- Varying intercepts for `policy_period` and `clade` adjust for broad temporal/variant structure, but they do not remove all possible confounding.
- Complete-case filtering means the fitted population can differ across model sets when predictors have different missingness.
- Sampling settings affect the fitted frame. Always inspect `fit_frame_summary.csv` and `metadata.csv`.
- Large posterior directional probabilities are evidence of association under the model, not proof of practical importance.
- Effect sizes, credible intervals, diagnostics, and fitted-row definitions should be reviewed together.

## 16. Quality-Control Checklist

Before reporting a model:

1. Confirm the intended model appears in `model_grid.csv`.
2. Confirm `fit_frame_summary.csv` has the expected row count and sampling status.
3. Check `metadata.csv` for formula, row count, and outcome rate/mean.
4. Read `fit.log` for warnings or fit failures.
5. Confirm `diagnostics.csv` has acceptable divergences, BFMI, R-hat, and ESS.
6. Interpret `summary.csv` on the right scale: ORs for logistic, coefficients for linear.
7. Compare primary and expanded models only after checking that their complete-case/sampled frames are comparable.

## 17. Maintenance Guide

To add or change models:

1. Edit predictor constants and `default_model_specs(...)` in `sse_detection/lib/regression_prep.py`.
2. Run a dry run with `--dry-run` and inspect the selected model grid.
3. Rebuild prepared model grids and inspect complete-case rows.
4. Fit the selected model script with `--skip-existing` when adding new models to an existing output directory.
5. Review `fit.log`, `diagnostics.csv`, and `summary.csv` before using estimates in figures or text.

Keep responsibilities separated:

- `regression_prep.py`: data frames, predictor definitions, formulas, output paths.
- `bayesian_models.py`: fitting, posterior summaries, diagnostics, result writing.
- `regression_runner.py`: command-line orchestration, selection, logging, manifests.
- `mixing_models.py` and `composition_models.py`: thin domain-specific entrypoints.
