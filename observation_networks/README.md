# Observation Networks

Code and technical outputs for thesis Chapter 4:

> Observation and Network Structure in Scottish SARS-CoV-2 Genomic Surveillance

This directory owns the descriptive observation/network analysis that sits
between the EpiLink pipeline and the Chapter 5 superspreading-compatible signal
detector. It should not create or consume candidate labels.

## Structure

- `lib/config.py`: shared Chapter 4 constants, output paths, and categorical
  attribute definitions.
- `lib/io.py`: loaders and writers built around `scotland/utils`.
- `lib/cohort.py`: observed cohort, rolling-window coverage, denominator, and
  sequence-composition summaries.
- `lib/clusters.py`: window-level cluster summaries and modal cluster
  composition.
- `lib/mixing.py`: weighted categorical mixing matrices and assortativity
  coefficients.
- `lib/simd.py`: SIMD population-weighting validation tables for the appendix.
- `lib/transition_graph.py`: alternate-window temporal cluster-transition graph
  construction and baseline summaries.
- `lib/figures.py`: reusable figure builders using `utils.style`.
- `build_tables.py`: builds the standard Chapter 4 tables.
- `build_mixing.py`: builds compatibility-network mixing matrices from sparse
  pairwise EpiLink edges.
- `build_simd_validation.py`: builds compact population-weighted SIMD validation
  tables for the appendix.
- `make_figures.py`: rebuilds figures from saved tables.
- `TECHNICAL.md`: analysis contract and output definitions.

Generated files are written under `observation_networks/results/`, which is
ignored by git and recreated by the build commands when needed.

## Main Commands

Run from the Scotland repository root.

Run the whole Chapter 4 build:

```bash
bash observation_networks/run_all.sh --workers 5
```

Using the project conda environment:

```bash
bash observation_networks/run_all.sh --conda-env PhD --workers 5
```

```bash
python -m observation_networks.build_tables
```

Build compatibility-network mixing for a small development window set:

```bash
python -m observation_networks.build_mixing --windows W080 W081
```

Add empirical permutation p-values for assortativity:

```bash
python -m observation_networks.build_mixing --windows W080 --n-permutations 1000 --missing-label Unknown
```

Build compatibility-network mixing for all retained windows:

```bash
python -m observation_networks.build_mixing --all-windows --workers 5
```

This also writes the compatibility-network degree/strength assortativity
diagnostic table used for the topology figure.

Build the SIMD population-weighting validation tables:

```bash
python -m observation_networks.build_simd_validation
```

Regenerate figures from available tables:

```bash
python -m observation_networks.make_figures --skip-missing
```

## Relationship To Chapter 5

Chapter 4 describes the observed sequenced record, window-level clusters,
compatibility-network mixing, and the transition-graph baseline. Chapter 5 may
reuse the graph-ready transition inputs, but candidate assignment and
superspreading-compatible scoring remain outside this directory.
