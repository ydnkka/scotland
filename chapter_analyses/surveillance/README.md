# Genomic Surveillance and Policy Timeline

Descriptive surveillance analysis for the Scottish SARS-CoV-2 sequenced record, restriction periods, sequencing intensity, and lineage replacement over time.

## Documentation map

- Use this README for the package map, main command, and output locations.
- Use `TECHNICAL.md` for data units, policy/stringency provenance, calculations, output definitions, and interpretation limits.
- Use `../genomic_networks/README.md` for Chapter 4 clusters and EpiLink compatibility networks.
- Use `../sse_detection/README.md` for Chapter 5 transition graphs and superspreading-compatible candidate detection.

## Structure

- `policy_sequences_over_time.py`: builds surveillance tables and the combined policy, sequence-count, lineage-frequency, and sequencing-proportion figure.
- `policy_index_comparison.py`: builds the two-panel Stringency versus Containment and Health Index validation figure and its supporting tables.
- `lib/config.py`: stable project/result paths and default smoothing/window-stride constants.
- `results/tables/`: index-free CSV and parquet surveillance tables.
- `results/figures/`: PNG and PDF surveillance figures.

## Main command

Run from the Scotland repository root in the project environment:

```bash
python -m chapter_analyses.surveillance.policy_sequences_over_time
```

Development and output options include:

```bash
python -m chapter_analyses.surveillance.policy_sequences_over_time --smooth-window 7 --window-stride 3 --log-level INFO
```

Use `--figure-dir` and `--table-dir` to override output locations.

Build the policy-index comparison figure and supporting correlation tables:

```bash
python -m chapter_analyses.surveillance.policy_index_comparison
```

## Outputs

The command writes:

- `results/tables/clade_frequency_by_period.{csv,parquet}`;
- `results/tables/clade_frequency_by_period_smoothed.{csv,parquet}`;
- `results/tables/clade_counts_by_period.{csv,parquet}`;
- `results/tables/clade_dominance_by_period.{csv,parquet}`;
- `results/tables/clade_overtake_events.{csv,parquet}`;
- `results/tables/sequencing_proportion_by_period.{csv,parquet}`;
- `results/figures/policy_sequences_over_time.{png,pdf}`.
- `results/tables/policy_indices_daily.{csv,parquet}`;
- `results/tables/policy_index_correlation.{csv,parquet}`;
- `results/figures/policy_index_comparison.{png,pdf}`.

## Policy stringency source

`utils/policy.py` reads the Scotland rows directly from the Stringency and Containment and Health Index tables under `data/raw/oxcgrt/`, reshapes the dated columns into daily tables, and calculates inclusive means for each predefined Scottish policy period. No hard-coded period index values remain.

## Scope

This package describes the sequenced surveillance record and policy/lineage timing. It does not estimate policy effects, transmission, vaccine effectiveness, compatibility-network mixing, or superspreading events.
