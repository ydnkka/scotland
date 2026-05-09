# Part 3 Analysis Documentation

## Overview

Part 3 asks how selected Scottish COVID-19 policy phases coincided with genomic
cluster structure.  It is deliberately not an exhaustive period catalogue.  The
analysis keeps the full 16-period context table, then develops four policy
moments where policy, variant advantage, and cluster structure are especially
informative.

The implementation is in `part3/part3_analysis.py`.

## Data Sources And Filters

Cluster-level policy summaries use the Part 1 cache:

- `part1/main/cache/main_cluster_table.parquet`

This cache is treated as the primary cluster-level input because it already
contains the Part 1 cluster outcomes and excess-mixing metrics.  The unit is
one inferred cluster/window row.

Sequence-level Alpha and mutation analyses use:

- `data/processed/scotland_clustering_analysis_dataset.parquet`
- raw Nextclade TSV from `config.yaml`, with fallback to
  `data/raw/cog_all_scotland_nextclade.tsv`

The sequence-level analysis filters to Leiden resolution `0.3` and
`nextclade_qc == "good"`, matching Part 1.  Alpha is identified from Pango
lineages beginning `B.1.1.7`.  Mutation trajectories are rebuilt from
Nextclade amino-acid substitutions for `S:N501Y`, `S:A222V`, and supporting
Alpha markers.

## Policy Attachment

Policy periods are attached with `utils.policy.attach_period_pandas` using
`wn_mid_date`.  Period definitions are not duplicated in the Part 3 scripts.

The selected periods are:

- `P3`, `T1`, `F5`, `L2`, `SL`, `L0`, and `NN`

The full context table also retains all other periods, including `OM`, `FE`,
and `PR`.

## Whole-Epidemic Policy Context

The script writes:

- `period_descriptives.csv`
- `weekly_summaries.csv`
- `intensity_correlations.csv`

Weekly summaries include total clusters, non-singleton clusters, singleton
fraction, median log cluster size, median log datazones, SIMD excess
discordance, age excess discordance, sex excess discordance, joint-profile
excess discordance, dominant policy period, and policy intensity.

Spearman correlations between policy intensity and weekly outcomes are reported
only as descriptive, confounded context.

## Interrupted Time-Series Summaries

The selected ITS transitions are:

- `t1_onset`: 2020-10-02, `P3 -> T1`
- `l2_to_sl`: 2021-04-02, `L2 -> SL`
- `nn_onset`: 2021-08-09, `L0 -> NN`

The model is:

```text
y_t = beta0 + beta1*t + beta2*post + beta3*(post*t) + error
```

The primary window is plus/minus 8 weeks around the transition.  Sensitivity
fits are run at plus/minus 6, 10, and 12 weeks.  Coefficients are written to
`its_coefficients.csv`; primary fitted weekly data are written to
`its_weekly_<transition>.csv`.

The main outcomes are:

- median log cluster size
- median log datazones
- mean SIMD excess discordance
- mean age excess discordance

The main figure uses the two structural outcomes.  Mixing outcomes are included
in the supplementary figure.

## Alpha/F5-L2 Case Study

Alpha emergence is rebuilt from the primary sequence-level table.  The three
reported phases are selected by observed window mid-dates:

- `cryptic_early`: W016-W021
- `multi_region_expansion`: W022-W024
- `f5_l2_bridge`: W025

The script writes:

- `alpha_phase_summary.csv`
- `alpha_cluster_emergence.csv`
- `alpha_health_board_weekly.csv`
- `alpha_local_authority_weekly.csv`
- `alpha_mutation_trajectories.csv`

Mutation trajectories are estimated from Nextclade markers joined to processed
sequence IDs.  `S:N501Y` is used as the Alpha marker and `S:A222V` as the
B.1.177 comparator.

## Growth Models And Counterfactuals

Binomial GLMs are fitted to weekly marker frequencies:

- Alpha marker `S:N501Y` during `F5`
- Alpha marker `S:N501Y` during `L2`
- B.1.177 marker `S:A222V` during `L2`

The primary model uses weekly positive tests as frequency weights.  Sensitivity
models use sequence-count weights, equal weekly weights, and coverage-adjusted
weights.  Results are written to:

- `alpha_growth_params.csv`
- `alpha_growth_model_sensitivity.csv`

The counterfactual timing analysis projects `S:N501Y` frequency under earlier
switches to the fitted L2 growth rate:

- actual L2 start, 2021-01-05
- 2020-12-08
- 2020-12-02, using the nearest observed expansion window
- 2020-11-02, the F5 start

Outputs are:

- `alpha_counterfactual_projections.csv`
- `alpha_counterfactual_trajectories.csv`

These are descriptive fitted scenarios, not causal policy estimates.

## Caveats

Policy periods are strongly entangled with variant replacement, changing
population immunity, testing eligibility and behaviour, sequencing coverage,
hospital burden, and seasonality.  The chapter should avoid claims that policy
caused a cluster outcome, that sequence data identify a specific
superspreading event, or that an earlier restriction scenario would have
prevented Alpha establishment.
