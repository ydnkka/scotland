# Genomic Surveillance and Policy Timeline

Descriptive analysis of sequence volume, Scottish restriction periods, lineage replacement, and sequencing coverage. It does not estimate policy effects, transmission, or vaccine effectiveness.

See [TECHNICAL.md](TECHNICAL.md) for definitions and output contracts.

## Commands

Run from the repository root:

```bash
python -m analyses.surveillance.policy_sequences_over_time
python -m analyses.surveillance.policy_index_comparison
```

The timeline command accepts `--smooth-window`, `--window-stride`, `--figure-dir`, `--table-dir`, and `--log-level`. The policy comparison accepts date, region, and output-directory overrides; use `--help` for the live interface.

## Outputs

`policy_sequences_over_time` writes CSV/parquet tables for raw and smoothed clade frequencies, clade counts, dominance, overtake events, and sequencing proportion, plus `policy_sequences_over_time.{png,pdf}`.

`policy_index_comparison` writes `policy_indices_daily` and `policy_index_correlation` as CSV/parquet, plus `policy_index_comparison.{png,pdf}`.

Default locations are:

```text
analyses/surveillance/results/tables/
analyses/surveillance/results/figures/
```

Policy-period, window-policy, and lineage/clade lookup helpers live in the shared `utils` package and work from the processed daily policy calendar produced by the pipeline.
