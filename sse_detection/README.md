# SSE Detection

Code and technical outputs for thesis Chapter 5:

> Superspreading-Compatible Signals in EpiLink-Derived Temporal Cluster-Transition Graphs

This directory owns the superspreading-compatible signal detector and the
Bayesian characterisation workflow. Candidate status is assigned by a
graph-first detector before demographic, geographic, deprivation, policy,
or mixing variables are used for characterisation.

## Structure

- `DETECTION_RATIONALE.md`: methodological rationale for the two-axis
  detector, including the detection/characterisation boundary.
- `TECHNICAL.md`: Bayesian regression workflow and output contract.
- `lib/sse/`: detector pipeline package: configuration, output loading,
  entropy utilities, cluster features, scoring, transition graph construction,
  detector orchestration, and diagnostics.
- `lib/model/`: Bayesian model package: data preparation, Bambi/PyMC fitting,
  runner orchestration, and domain CLIs.
- `lib/forest.py`: forest-plot helpers for saved Bayesian outputs.
- `lib/figures.py` and `lib/figs/`: detector figure orchestration and figure
  builders.
- `sse_detection.ipynb`: exploratory detector development notebook.
- `composition_mixing_variable_correlations.ipynb`: exploratory correlation
  checks between composition and mixing variables.
- `bayesian_regression_forest_plots.ipynb`: figure-generation notebook for
  selected Bayesian forest plots.
- `make_figures.py`: rebuilds detector figures from saved SSE output tables.

Generated outputs are written under `sse_detection/results/`.

## Documentation Map

- Use this README for the directory map, main detector command, and top-level
  output locations.
- Use `DETECTION_RATIONALE.md` for detector design, null calibration,
  candidate interpretation, and method limitations.
- Use `TECHNICAL.md` for Bayesian model data frames, formulas, priors, run
  options, diagnostics, output files, and maintenance notes.

## Main Command

Run from the Scotland repository root.

Build the detector outputs:

```bash
python -m sse_detection.lib.sse.detection
```

Regenerate detector figures and Bayesian result tables from saved outputs:

```bash
python -m sse_detection.make_figures --skip-missing
```

The Bayesian mixing-model tables can also be regenerated independently:

```bash
python -m sse_detection.make_tables
```

Build cluster-level socio-demographic composition tables for downstream SSE
characterisation:

```bash
python -m sse_detection.build_composition_tables
```

See `TECHNICAL.md` for Bayesian dry-run examples, model selectors, fitting
commands, options, diagnostics, and interpretation guidance.

## Output Layout

- `results/sse_outputs/`: detector outputs, including the regression-facing
  cluster table, directed adjacent-window edge table, and transition
  node/window/component summaries.
- `results/figures/`: detector figures built from `results/sse_outputs/`.
- `results/bayesian_outputs/`: Bayesian candidate, burst-score, and
  burden-score characterisation outputs.
- `results/bayesian_outputs/figures/`: forest plots generated from selected
  saved Bayesian outputs.

`DETECTION_RATIONALE.md` defines which graph-derived quantities enter candidate
detection. `TECHNICAL.md` defines the fitted-model output contract.

## Relationship To Chapter 4

Chapter 4, implemented in `observation_networks`, describes the observed
sequenced record, compatibility-network mixing/assortativity, and window-level
clusters. This directory owns the alternate-window transition graph, candidate
prioritisation, and Bayesian characterisation. Results should not be
interpreted as confirmed superspreading events without external epidemiological
evidence.
