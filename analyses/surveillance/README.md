# Genomic Surveillance and Policy Timeline

Descriptive analysis of sequence volume, Scottish restriction periods, lineage replacement, and sequencing coverage. It does not estimate policy effects, transmission, or vaccine effectiveness.

See [TECHNICAL.md](TECHNICAL.md) for definitions and output contracts.

## Commands

Run from the repository root:

```bash
python -m results.make_figures --domain surveillance
python -m results.make_tables --domain surveillance
```

The individual figure builders live in `lib/figs/fig01.py` and
`lib/figs/fig02.py`. Project-level publication outputs are built through the
root `results` commands.

## Outputs

`policy_sequences_over_time` writes CSV/parquet tables for raw and smoothed clade frequencies, clade counts, dominance, overtake events, and sequencing proportion, plus `policy_sequences_over_time.{png,pdf}`.

`policy_index_comparison` writes `policy_indices_daily` and `policy_index_correlation` as CSV/parquet, plus `policy_index_comparison.{png,pdf}`.

Default output locations are:

```text
analyses/surveillance/results/tables/
results/figures/
```

Policy-period, window-policy, and lineage/clade lookup helpers live in the shared `utils` package and work from the processed daily policy calendar produced by the pipeline.
