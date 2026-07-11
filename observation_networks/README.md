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
- `lib/figures.py`: orchestration layer for the thesis figure scripts in
  `lib/figs/`.
- `lib/figs/`: thesis figure and LaTeX table builders using `utils.style`.
- `build_cluster_tables.py`: builds the standard Chapter 4 tables.
- `build_mixing.py`: builds compatibility-network mixing matrices from sparse
  pairwise EpiLink edges, processing one window-lineage pairwise parquet file
  per worker and concatenating the intermediate chunks.
- `build_cluster_pairwise_distance_summary.py`: summarises pairwise SNP and
  temporal distances among sequences from selected non-singleton clusters
  within each window-lineage pairwise file.
- `build_simd_validation.py`: builds compact population-weighted SIMD validation
  tables for the appendix.
- `make_figures.py`: rebuilds figures and LaTeX table fragments from saved
  tables.
- `TECHNICAL.md`: analysis contract and output definitions.

Generated files are written under `observation_networks/results/`, which is
ignored by git and recreated by the build commands when needed.

## Main Commands

Run from the Scotland repository root.

Run the core Chapter 4 tables:

```bash
python -m observation_networks.build_cluster_tables
```

Run the full observation-network build by executing the table, mixing, SIMD
validation, and figure builders in sequence:

```bash
python -m observation_networks.build_cluster_tables
python -m observation_networks.build_mixing --all-windows --workers 4
python -m observation_networks.build_simd_validation
python -m observation_networks.make_figures --skip-missing
```

Build compatibility-network mixing for a small development window set:

```bash
python -m observation_networks.build_mixing --windows W080 W081
```

This processes every pairwise lineage file inside the requested windows. Use
`--max-windows N` to keep all lineages from only the first `N` selected windows
for development runs.

Summarise pairwise SNP and temporal distances for non-singleton cluster
sequences in a small development window set:

```bash
python -m observation_networks.build_cluster_pairwise_distance_summary \
  --windows W080 W081 --max-clusters-per-window-lineage 25
```

Cap adaptive large-file jackknife blocks for quick uncertainty checks:

```bash
python -m observation_networks.build_mixing --windows W080 --jackknife-blocks 100
```

Build compatibility-network mixing for all retained windows:

```bash
python -m observation_networks.build_mixing --all-windows --workers 4
```

By default this skips giant pairwise files, defined by `--giant-threshold`
(`50,000,000` sparse edges by default). For the full run including those files:

```bash
python -m observation_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
```

Each worker handles one `data/processed/pairwise_distances_dataset/*.parquet`
file. Per-file intermediate parquet outputs are written under
`observation_networks/results/intermediate/`, and the final concatenated
compatibility-network mixing, assortativity, and degree/strength assortativity
tables are written under `observation_networks/results/tables/`. The
assortativity table includes deterministic jackknife uncertainty by default:
attributes with up to 1,000 contributing vertices use leave-one-node jackknife,
and larger attributes use adaptive node blocks capped by `--jackknife-blocks`.
Pass `--jackknife-blocks 0` to skip uncertainty columns. At `INFO`, progress is
logged every 100 completed pairwise files by default; pass `--progress-every N`
to `build_mixing.py` to tune this. Per-file progress is available with
`--log-level DEBUG`.

Build the SIMD population-weighting validation tables:

```bash
python -m observation_networks.build_simd_validation
```

Regenerate figures from available tables:

```bash
python -m observation_networks.make_figures --skip-missing
```

## Relationship To Chapter 5

Chapter 4 describes the observed sequenced record, window-level clusters, and
compatibility-network mixing/assortativity. The alternate-window temporal
cluster-transition graph, candidate assignment, and superspreading-compatible
scoring are owned by `sse_detection`.
