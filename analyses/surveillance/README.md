# Genomic Surveillance and Policy Timeline

Descriptive analysis of sequence volume, Scottish restriction periods, lineage replacement, and sequencing coverage. It does not estimate policy effects, transmission, or vaccine effectiveness.

See [TECHNICAL.md](TECHNICAL.md) for definitions and output contracts.

## Commands

Run from the repository root:

```bash
# Build PNG/PDF figures under results/figures
python -m results.make_figures --domain surveillance

# Refresh companion CSV/parquet tables under analyses/surveillance/results/tables
python -m results.make_tables --domain surveillance
```

The individual figure builders live in `lib/figs/fig01.py` and
`lib/figs/fig02.py`. The central `results` commands provide the build entry
points; this analysis does not currently define root-level LaTeX table
fragments.

## Outputs

`policy_sequences_over_time` writes analysis-local CSV/parquet tables for raw
and smoothed clade frequencies, clade counts, dominance, overtake events, and
sequencing proportion, plus project-level `policy_sequences_over_time.{png,pdf}`.

`policy_index_comparison` writes `policy_indices_daily` and
`policy_index_correlation` as analysis-local CSV/parquet tables, plus
project-level `policy_index_comparison.{png,pdf}`.

Default output locations are:

```text
analyses/surveillance/results/tables/
results/figures/
```

Policy-period, window-policy, and lineage/clade lookup helpers live in the shared `utils` package and work from the processed daily policy calendar produced by the pipeline.
