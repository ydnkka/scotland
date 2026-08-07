# Genomic Surveillance and Compatibility Networks

Chapter 4 code for cohort description, window-specific clusters, EpiLink compatibility networks, mixing, topology, and sensitivity analysis. Candidate labels and temporal transitions belong to [`../sse_detection`](../sse_detection/README.md).

See [TECHNICAL.md](TECHNICAL.md) for algorithms and table definitions.

## Build order

Run from the repository root:

```bash
python -m chapter_analyses.genomic_networks.build_cluster_summaries
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
python -m chapter_analyses.genomic_networks.build_sensitivity_tables
python -m chapter_analyses.genomic_networks.build_simd_validation
python -m chapter_analyses.genomic_networks.make_figures --skip-missing
```

Build `data/processed/sparsified_edge_counts_by_window_lineage.parquet` first with `python3 method/build_sparsified_edge_manifest.py` for edge-aware mixing scheduling.

## Commands

`build_cluster_summaries` writes the cohort, context, and cluster summary tables including within-cluster pairwise distance summaries, the overall pairwise-distance summary, and the combined period-level typical-summary table with duration plus genetic, temporal, and residential spatial distance. It processes all windows by default; use `--windows ...` or `--max-windows N` only for a pairwise-distance subset, `--reuse-input-tables` to refresh from existing cluster/pairwise inputs, or `--skip-pairwise-distances` for cluster-only outputs.

`build_mixing` requires either `--all-windows` or `--windows ...`. Useful invocations:

```bash
# Inspect scheduling only
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --dry-run

# Small selection
python -m chapter_analyses.genomic_networks.build_mixing --windows W080 W081 --workers 2

# All non-giant files; giant files are skipped by default
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4
```

Files with at least 50,000,000 sparse edges, and files with unknown cost, are giant by default. Use `--include-giants` to process them in the separate `--giant-workers` pool. Same-stem intermediate chunks are reused unless `--force` is supplied; do not reuse them after changing attributes, threshold, missing-value handling, bootstrap settings, or minimum-edge filtering.

Assortativity uncertainty uses an edge-weight multiplier bootstrap by default (`--bootstrap-replicates 500`, `--bootstrap-alpha 0.05`, `--bootstrap-seed 123`). Use `--bootstrap-replicates 0` to skip uncertainty. `--min-edges` keeps below-threshold analyses in the output with `NaN` estimates and a `skipped_reason`.

`build_sensitivity_tables` builds Leiden and sparsification sensitivity families. Select one with `--only leiden` or `--only sparsification`. Row/file caps and partial row-group scans are for development and can change or approximate the estimand.

`build_simd_validation` defaults to quintiles; `--n-groups` accepts 5, 10, or 20.

`make_figures` reads saved tables. Without `--skip-missing`, a missing required table fails the run.

## Layout

- `lib/`: configuration, I/O, cohort, cluster-table, cluster-rollup, pairwise-distance, mixing, SIMD, figure, and table logic.
- `build_*.py`: table-building entry points.
- `make_figures.py`: saved-table figure and LaTeX orchestration.
- `results/tables/`: final tables.
- `results/intermediate/`: per-pairwise-file bootstrap mixing, assortativity, and topology chunks.
- `results/figures/`: figures and `.tex` fragments.

Compatibility edges are EpiLink-weighted plausible links, not observed transmission. Rolling-window clusters are not persistent outbreak entities.
