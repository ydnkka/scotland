# SSE Detection

Chapter 5 code for an alternate-window cluster-transition graph, a demographic-blind superspreading-compatible signal detector, composition summaries, and Bayesian characterisation.

Candidates are graph-derived priorities for review, not confirmed superspreading events or transmission chains.

## Documentation

- [DETECTION_RATIONALE.md](DETECTION_RATIONALE.md): why burst and burden are separate detection axes.
- [TECHNICAL.md](TECHNICAL.md): implemented graph, scores, calibration, tiers, and outputs.
- [BAYESIAN_MODELS.md](BAYESIAN_MODELS.md): model frames, formulas, commands, priors, and diagnostics.

## Commands

Run from the repository root:

```bash
# Rebuild detector outputs
python -m chapter_analyses.sse_detection.lib.sse.detection

# Build wide cluster-composition tables
python -m chapter_analyses.sse_detection.build_composition_tables

# Rebuild available detector figures and Bayesian result tables
python -m chapter_analyses.sse_detection.make_figures --skip-missing

# Rebuild Bayesian result tables only
python -m chapter_analyses.sse_detection.make_tables
```

The detector has no CLI options; its settings are constants in `lib/sse/config.py`. Bayesian fitting commands are listed in [BAYESIAN_MODELS.md](BAYESIAN_MODELS.md).

## Layout

- `lib/sse/`: data loading, transition graph, composition entropy, feature assembly, scoring, and diagnostics.
- `lib/model/`: Bayesian preparation, fitting, and CLI orchestration.
- `lib/figs/` and `lib/figures.py`: figure and result-table builders.
- `results/sse_outputs/`: authoritative detector cluster/edge tables and graph summaries.
- `results/tables/`: composition and publication tables.
- `results/figures/`: detector and Bayesian figures.
- `results/bayesian_outputs/`: fitted model summaries, diagnostics, metadata, logs, and optional inference data.

Chapter 4 owns the underlying window-specific clusters and compatibility networks; this package owns their temporal continuity and downstream characterisation.
