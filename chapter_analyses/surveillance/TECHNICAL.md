# Genomic Surveillance and Policy Timeline: Technical Reference

## Scope and inputs

`policy_sequences_over_time.py` loads `sequence_id`, `collection_date`, `clade`, and `wn_prop_sequenced` through `utils.load_analysis_columns(...)`. It applies no Nextclade QC filter and, by default, retains every third sorted analysis window (`--window-stride 3`). Retained windows are renumbered by the shared loader. Unmatched clades are labelled `Other`.

Results describe observed sequenced records; the package does not fit causal or transmission models.

## Calculations

- Daily sequence count: unique `sequence_id` values by collection date, with missing calendar dates filled as zero.
- Daily trend: centred rolling mean, default 7 days (`--smooth-window`).
- Policy period: ordered, inclusive date intervals in `utils/policy.py`.
- Lineage composition: weekly unique-sequence counts and within-week frequencies. The plotted frequency is a trailing 3-week mean fixed inside `plot_lineage_frequency_and_overtakes(...)`.
- Dominance: clade with the largest plotted frequency; an overtake is a change from the previous retained week.
- Sequencing proportion: weekly mean of row-level `wn_prop_sequenced`, aligned to retained lineage-frequency weeks; the plotted value uses the same trailing 3-week mean.

Because rolling-window rows repeat sequences, the loader stride and deduplication rules are part of the estimand and must be recorded.

## Policy indices

`utils/policy.py` reads exactly one `RegionName == "Scotland"` row from:

```text
data/raw/oxcgrt/OxCGRT_timeseries_StringencyIndex_v1.csv
data/raw/oxcgrt/OxCGRT_timeseries_ContainmentHealthIndex_v1.csv
```

Columns matching `%d%b%Y` are reshaped to daily observations. Each policy period receives the arithmetic mean over its inclusive start/end dates. Missing files, schema errors, or a non-unique Scotland row fail during import.

`policy_index_comparison.py` inner-joins the two daily series, restricts them to the configured dates, and reports complete-day Pearson correlation, Spearman correlation, an ordinary linear slope/intercept, and Pearson (r^2). It intentionally omits independence-based p-values because adjacent daily policy values are serially dependent.

## Output contract

Tables are index-free CSV and parquet files under `results/tables/`:

- `clade_frequency_by_period`: unsmoothed weekly frequencies;
- `clade_frequency_by_period_smoothed`: plotted trailing means;
- `clade_counts_by_period`: total and clade counts;
- `clade_dominance_by_period` and `clade_overtake_events`;
- `sequencing_proportion_by_period`: observed and plotted coverage;
- `policy_indices_daily` and `policy_index_correlation`.

Figures are PNG and PDF only under `results/figures/`.

## Checks and limits

Record input versions, stride, smoothing, and policy dates. Check unmatched clades, missing coverage, complete OxCGRT days, and CSV/parquet agreement.

Sequence coverage, testing, variants, vaccination, behaviour, and policy all vary over time. Policy indices measure formal restrictions rather than adherence, and observed lineage dominance depends on sampling and aggregation.
