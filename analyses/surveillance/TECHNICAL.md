# Genomic Surveillance and Policy Timeline: Technical Reference

## Scope and inputs

`lib/figs/fig01.py` builds the sequence-over-time figure. It loads
`sequence_id`, `collection_date`, `clade`, and `wn_prop_sequenced` through
`utils.load_analysis_columns(...)`. It applies no Nextclade QC filter and, by
default, retains every third sorted analysis window. Retained windows are
renumbered by the shared loader. Unmatched clades are labelled `Other`.

Results describe observed sequenced records; the package does not fit causal or transmission models.

## Calculations

- Daily sequence count: unique `sequence_id` values by collection date, with missing calendar dates filled as zero.
- Daily trend: centred rolling mean, default 7 days (`--smooth-window`).
- Policy period: ordered, inclusive date intervals derived from the processed daily policy calendar via shared `utils` helpers.
- Lineage composition: weekly unique-sequence counts and within-week frequencies. The plotted frequency is a trailing 3-week mean fixed inside `plot_lineage_frequency_and_overtakes(...)`.
- Dominance: clade with the largest plotted frequency; an overtake is a change from the previous retained week.
- Sequencing proportion: weekly mean of row-level `wn_prop_sequenced`, aligned to retained lineage-frequency weeks; the plotted value uses the same trailing 3-week mean.

Because rolling-window rows repeat sequences, the loader stride and deduplication rules are part of the estimand and must be recorded.

## Policy indices

The shared policy helpers read the processed daily calendar at `data/processed/scotland_policy.parquet`, normalise the policy-period columns to `policy_period`, `policy_period_label`, `policy_period_order`, and `policy_era`, and expose the same date-indexed information used by the figure code. The parquet is produced upstream by `method/01_prep_metadata.py`, so a stale or missing policy build still fails when the pipeline is rebuilt.

`lib/figs/fig02.py` inner-joins the two daily policy-index series, restricts
them to the configured dates, and reports complete-day Pearson correlation,
Spearman correlation, an ordinary linear slope/intercept, and Pearson (r^2). It
intentionally omits independence-based p-values because adjacent daily policy
values are serially dependent.

## Output contract

When built through the project-level `results` commands, tables are index-free
CSV and parquet files under `analyses/surveillance/results/tables/`:

- `clade_frequency_by_period`: unsmoothed weekly frequencies;
- `clade_frequency_by_period_smoothed`: plotted trailing means;
- `clade_counts_by_period`: total and clade counts;
- `clade_dominance_by_period` and `clade_overtake_events`;
- `sequencing_proportion_by_period`: observed and plotted coverage;
- `policy_indices_daily` and `policy_index_correlation`.

Figures are PNG and PDF only under project-level `results/figures/`.

## Checks and limits

Record input versions, stride, smoothing, and policy dates. Check unmatched clades, missing coverage, complete OxCGRT days, and CSV/parquet agreement.

Sequence coverage, testing, variants, vaccination, behaviour, and policy all vary over time. Policy indices measure formal restrictions rather than adherence, and observed lineage dominance depends on sampling and aggregation.
