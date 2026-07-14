# Genomic Surveillance and Compatibility Networks

Code and technical outputs for thesis Chapter 4:

> Observation and Network Structure in Scottish SARS-CoV-2 Genomic Surveillance

This directory owns the descriptive observation and compatibility-network analysis between the EpiLink pipeline and the Chapter 5 superspreading-compatible signal detector. It must not create or consume candidate labels.

## Documentation map

- Use this README for the directory map, execution order, common commands, and output locations.
- Use `TECHNICAL.md` for the analysis boundary, data contract, algorithms, uncertainty methods, table schemas, sensitivity analyses, and interpretation limits.
- Use `../sse_detection/README.md` for the temporal cluster-transition graph and Chapter 5 candidate detector.

## Structure

- `lib/config.py`: Chapter 4 paths, analysis constants, disclosure threshold, and categorical attribute specifications.
- `lib/io.py`: standard loaders, result-directory creation, and CSV/parquet table I/O.
- `lib/cohort.py`: cohort, rolling-window coverage, denominator, composition, and vaccination-context summaries.
- `lib/clusters.py`: window-specific EpiLink cluster tables and cluster composition summaries.
- `lib/mixing.py`: weighted categorical mixing, nominal assortativity, node-jackknife uncertainty, and degree/strength assortativity.
- `lib/simd.py`: population-weighted SIMD grouping validation.
- `lib/figures.py` and `lib/figs/`: figure and LaTeX-table orchestration from saved result tables.
- `build_cluster_tables.py`: builds the core cohort, coverage, vaccination, cluster, and composition tables.
- `build_mixing.py`: scans per-window/per-lineage sparse pairwise files and builds compatibility mixing and assortativity outputs.
- `build_cluster_pairwise_distance_summary.py`: summarises SNP and temporal distances among sequences in selected non-singleton clusters.
- `build_sensitivity_tables.py`: builds Leiden-resolution and compatibility-sparsification sensitivity tables.
- `build_simd_validation.py`: builds population-weighted SIMD validation tables.
- `make_figures.py`: regenerates all available Chapter 4 figures and LaTeX table fragments from saved tables.

Generated outputs are written under `chapter_analyses/genomic_networks/results/`.

## Environment and working directory

Run commands from the Scotland repository root using the project environment that provides pandas, PyArrow, NumPy, scikit-learn, plotting libraries, and the repository utilities.

## Recommended execution order

Build the complete Chapter 4 analysis in this order:

```bash
python -m chapter_analyses.genomic_networks.build_cluster_tables
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
python -m chapter_analyses.genomic_networks.build_cluster_pairwise_distance_summary --all-windows
python -m chapter_analyses.genomic_networks.build_sensitivity_tables
python -m chapter_analyses.genomic_networks.build_simd_validation
python -m chapter_analyses.genomic_networks.make_figures --skip-missing
```

The mixing command can be run without `--include-giants` when memory-heavy pairwise files should be skipped. Figures should be regenerated only after their required tables have been refreshed.

## Core tables

Build cohort, coverage, denominator, vaccination, cluster, and composition tables:

```bash
python -m chapter_analyses.genomic_networks.build_cluster_tables
```

Use `--max-windows N` for a development subset. The transition-related flags retained by this command are deprecated no-ops because Chapter 5 owns transition construction.

## Compatibility mixing and assortativity

Inspect the edge-count-aware schedule without processing files:

```bash
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --dry-run
```

Run a small development selection:

```bash
python -m chapter_analyses.genomic_networks.build_mixing --windows W080 W081 --workers 2
```

Run all selected non-giant files:

```bash
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4
```

Run the full build, including giant files in a conservative separate worker pool:

```bash
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
```

Each worker processes one physical `data/processed/pairwise_distances_dataset/*.parquet` file and writes same-stem intermediate chunks under `results/intermediate/`. Existing chunks are reused unless `--force` is supplied. Files with at least 50,000,000 sparse edges are classified as giant by default; unknown costs are also treated as giant.

Useful options include `--attributes`, `--compatibility-threshold`, `--max-windows`, `--missing-label`, `--jackknife-blocks`, `--jackknife-seed`, `--giant-threshold`, `--edge-manifest`, `--progress-every`, and `--log-level`. Use `--jackknife-blocks 0` only for development runs where uncertainty columns are not required.

## Pairwise distance summary

Build the complete non-singleton-cluster SNP and temporal-distance summary:

```bash
python -m chapter_analyses.genomic_networks.build_cluster_pairwise_distance_summary --all-windows
```

Run a capped development selection:

```bash
python -m chapter_analyses.genomic_networks.build_cluster_pairwise_distance_summary --windows W080 W081 --max-clusters-per-window-lineage 25
```

The command supports window and lineage filters, cluster-size and selection controls, resolution/QC overrides, custom output names and directories, and a dry run.

## Sensitivity tables

Build both Leiden-resolution and sparsification-threshold sensitivity families:

```bash
python -m chapter_analyses.genomic_networks.build_sensitivity_tables
```

Build one family only:

```bash
python -m chapter_analyses.genomic_networks.build_sensitivity_tables --only leiden
python -m chapter_analyses.genomic_networks.build_sensitivity_tables --only sparsification
```

For a quick approximate scan:

```bash
python -m chapter_analyses.genomic_networks.build_sensitivity_tables --max-windows 3 --max-files 10 --max-row-groups-per-file 2
```

## SIMD validation

Build population-weighted SIMD validation tables using quintiles:

```bash
python -m chapter_analyses.genomic_networks.build_simd_validation
```

The `--n-groups` option also supports 10 and 20 groups.

## Figures and LaTeX tables

Regenerate figures and LaTeX fragments from available saved tables:

```bash
python -m chapter_analyses.genomic_networks.make_figures --skip-missing
```

Without `--skip-missing`, the command fails when a required input table is absent.

## Output layout

- `results/tables/`: final CSV/parquet analysis tables.
- `results/intermediate/mixing_matrix/`: per-file compatibility mixing chunks.
- `results/intermediate/comp_assortativity/`: per-file categorical assortativity chunks.
- `results/intermediate/deg_assortativity/`: per-file topology-assortativity chunks.
- `results/figures/`: PNG/PDF figures and `.tex` table fragments.

## Relationship to Chapter 5

Chapter 4 describes the observed sequenced record, coverage, window-specific clusters, compatibility-network mixing, assortativity, topology, vaccination context, and parameter sensitivity. `chapter_analyses.sse_detection` owns alternate-window cluster transitions, upstream novelty, burst and onward-burden scoring, null calibration, candidate tiers, and Bayesian candidate characterisation.
