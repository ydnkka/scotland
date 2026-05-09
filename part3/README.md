# Part 3 Policy Periods And Variant-Phase Context

This directory rebuilds Part 3 as a descriptive policy chapter for the Scotland
SARS-CoV-2 genomic clustering project.  The chapter does not estimate causal
effects of Scottish COVID-19 policy from genomic cluster data alone.  Policy
periods are interpreted through variant phase, immunity, testing, sequencing,
and calendar time.

The main argument is:

> Genomic clusters can show policy-relevant structural changes, but policy
> signals are only interpretable through variant phase, surveillance, and
> timing.

## Scope

The analysis focuses on four selected policy phases:

- Autumn 2020, `P3 -> T1`: a relatively clean tightening during the B.1.177
  period.
- Winter 2020/21, `F5 -> L2`: the Alpha emergence centrepiece.
- Spring 2021, `L2 -> SL`: a warning example where easing overlaps with the
  Alpha tail.
- Summer 2021, `L0 -> NN`: the clearest policy-consistent increase in cluster
  geographic structure during Delta.

Other policy periods, including `OM`, `FE`, and `PR`, are retained in the full
context tables but are not developed as main narrative sections.

## Run

From the repository root:

```bash
conda run -n PhD python part3/part3_analysis.py
conda run -n PhD python part3/manuscript/make_figures.py
```

## Inputs

- `part1/main/cache/main_cluster_table.parquet`
- `data/processed/scotland_clustering_analysis_dataset.parquet`
- raw Nextclade TSV from `config.yaml` when present, with fallback to
  `data/raw/cog_all_scotland_nextclade.tsv`
- shared helpers in `utils/data.py`, `utils/policy.py`, and `utils/style.py`

The primary cluster filter follows Part 1: Leiden resolution `0.3` and
Nextclade QC `good`.

## Outputs

Main tables are written to `part3/tables/`, including:

- `period_descriptives.csv`
- `weekly_summaries.csv`
- `intensity_correlations.csv`
- `its_coefficients.csv`
- `its_weekly_t1_onset.csv`
- `its_weekly_l2_to_sl.csv`
- `its_weekly_nn_onset.csv`
- `alpha_phase_summary.csv`
- `alpha_mutation_trajectories.csv`
- `alpha_growth_params.csv`
- `alpha_growth_model_sensitivity.csv`
- `alpha_counterfactual_projections.csv`

Additional Alpha support tables are also written for health-board, local
authority, cluster-emergence, and counterfactual-trajectory plotting.

Manuscript figures are written to `part3/manuscript/figures/` as PDF, PNG, and
TIFF files.

## Interpretation

All policy-period analyses are associational.  The ITS models use policy dates
as descriptive anchors, not as exogenous shocks.  The Alpha counterfactuals are
fitted timing scenarios based on observed mutation frequencies and should not
be read as causal estimates of restriction effects.
