# Genomic Surveillance and Compatibility Networks

Chapter 4 code for cohort description, window-specific clusters, EpiLink compatibility networks, mixing, topology, and sensitivity analysis. Candidate labels and temporal transitions belong to [`../sse_detection`](../sse_detection/README.md).

See [TECHNICAL.md](TECHNICAL.md) for algorithms and table definitions.

## Build order

Run from the repository root:

```bash
python -m chapter_analyses.genomic_networks.build_cluster_tables
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
python -m chapter_analyses.genomic_networks.build_cluster_pairwise_distance_summary --all-windows
python -m chapter_analyses.genomic_networks.build_sensitivity_tables
python -m chapter_analyses.genomic_networks.build_simd_validation
python -m chapter_analyses.genomic_networks.make_figures --skip-missing
```

Build `data/processed/sparsified_edge_counts_by_window_lineage.parquet` first with `python3 method/build_sparsified_edge_manifest.py` for edge-aware mixing scheduling.

## Commands

`build_cluster_tables` writes the core cohort, coverage, composition, test-reason, vaccination, and cluster tables. `--max-windows N` is a development cap. Its transition flags are deprecated no-ops.

`build_mixing` requires either `--all-windows` or `--windows ...`. Useful invocations:

```bash
# Inspect scheduling only
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --dry-run

# Small selection
python -m chapter_analyses.genomic_networks.build_mixing --windows W080 W081 --workers 2

# All non-giant files; giant files are skipped by default
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4
```

Files with at least 50,000,000 sparse edges, and files with unknown cost, are giant by default. Use `--include-giants` to process them in the separate `--giant-workers` pool. Same-stem intermediate chunks are reused unless `--force` is supplied; do not reuse them after changing attributes, threshold, missing-value handling, or jackknife settings.

`build_cluster_pairwise_distance_summary` also requires `--all-windows` or `--windows`. It supports lineage, cluster-size, selection, resolution, QC, output, and development filters.

`build_sensitivity_tables` builds Leiden and sparsification sensitivity families. Select one with `--only leiden` or `--only sparsification`. Row/file caps and partial row-group scans are for development and can change or approximate the estimand.

`build_simd_validation` defaults to quintiles; `--n-groups` accepts 5, 10, or 20.

`make_figures` reads saved tables. Without `--skip-missing`, a missing required table fails the run.

## Layout

- `lib/`: configuration, I/O, cohort, cluster, mixing, SIMD, figure, and table logic.
- `build_*.py`: table-building entry points.
- `make_figures.py`: saved-table figure and LaTeX orchestration.
- `results/tables/`: final tables.
- `results/intermediate/`: per-pairwise-file mixing, assortativity, and topology chunks.
- `results/figures/`: figures and `.tex` fragments.

Compatibility edges are EpiLink-weighted plausible links, not observed transmission. Rolling-window clusters are not persistent outbreak entities.
