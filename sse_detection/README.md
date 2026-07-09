# SSE Detection

Code and technical outputs for thesis Chapter 5:

> Superspreading-Compatible Signals in EpiLink-Derived Temporal Cluster-Transition Graphs

This directory owns the superspreading-compatible signal detector and the
Bayesian characterisation workflow. Candidate status is assigned by a
graph-first detector before demographic, geographic, deprivation, policy,
vaccination, or mixing variables are used for characterisation.

## Structure

- `DETECTION_RATIONALE.md`: methodological rationale for the two-axis
  detector, including the detection/characterisation boundary.
- `TECHNICAL.md`: Bayesian regression workflow and output contract.
- `lib/sse_detection.py`: builds the alternate-window transition graph,
  cluster/node table, edge table, local-burst and onward-burden scores, null
  calibrations, and candidate tiers.
- `lib/io.py`: loaders for static SSE outputs.
- `lib/entropy.py`: cluster composition, entropy, and onward-edge entropy
  utilities.
- `lib/regression_prep.py`: model frames, formulas, predictors, outcomes, and
  output paths for Bayesian models.
- `lib/bayesian_models.py`: Bambi/PyMC fitting, posterior summaries, and
  diagnostics.
- `lib/regression_runner.py`: shared CLI orchestration for model fitting.
- `lib/mixing_models.py`: node-level mixing-regression CLI.
- `lib/composition_models.py`: sequence-level composition-regression CLI.
- `lib/forest.py`: forest-plot helpers for saved Bayesian outputs.
- `sse_detection.ipynb`: exploratory detector development notebook.
- `composition_mixing_variable_correlations.ipynb`: exploratory correlation
  checks between composition and mixing variables.
- `bayesian_regression_forest_plots.ipynb`: figure-generation notebook for
  selected Bayesian forest plots.

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
python -m sse_detection.lib.sse_detection
```

See `TECHNICAL.md` for Bayesian dry-run examples, model selectors, fitting
commands, options, diagnostics, and interpretation guidance.

## Output Layout

- `results/sse_outputs/`: detector outputs, including the cluster/node table
  and directed adjacent-window edge table.
- `results/bayesian_outputs/`: Bayesian candidate, burst-score, and
  burden-score characterisation outputs.
- `results/bayesian_outputs/figures/`: forest plots generated from selected
  saved Bayesian outputs.

`DETECTION_RATIONALE.md` defines which graph-derived quantities enter candidate
detection. `TECHNICAL.md` defines the fitted-model output contract.

## Relationship To Chapter 4

Chapter 4, implemented in `observation_networks`, describes the observed
sequenced record, compatibility-network mixing, window-level clusters, and the
transition-graph baseline. This directory starts from that graph-ready view and
adds candidate prioritisation plus Bayesian characterisation. Results should not
be interpreted as confirmed superspreading events without external
epidemiological evidence.
