# Genomic Surveillance and Policy Timeline: Technical Reference

This document defines the production data flow and outputs for the `chapter_analyses.surveillance` package. Each analysis result is descriptive of observed sequenced records.

## 1. Analysis boundary

The package summarises daily unique sequences, policy-period stringency, weekly lineage composition and dominance, lineage-overtake events, and the proportion of cases sequenced. `chapter_analyses.genomic_networks` owns cluster and compatibility-network analysis; `chapter_analyses.sse_detection` owns temporal transitions and candidate detection.

The production entry point is:

```bash
python -m chapter_analyses.surveillance.policy_sequences_over_time
```

## 2. Input sequence data

The entry point calls `utils.load_analysis_columns(...)` for `sequence_id`, `collection_date`, `clade`, and `wn_prop_sequenced`, with no QC restriction and a default analysis-window stride of 3. The stride is configurable through `--window-stride` and must be recorded with reported outputs.

Clades are mapped through the shared `utils.CLADES` labels, with unmatched values assigned to `Other`.

## 3. Daily sequence counts

Daily counts are the number of unique `sequence_id` values by `collection_date`. Missing calendar dates between the first and last observations are inserted with zero counts.

The displayed trend is a centred rolling mean:

```text
smoothed_count(date) = mean(daily count over the centred smoothing window)
```

The default smoothing window is 7 days and can be changed with `--smooth-window`.

## 4. Policy periods

Policy periods are defined centrally in `utils/policy.py` by ordered codes, labels, and inclusive start/end dates. `assign_period(...)` assigns dates to these intervals, while `attach_period(...)` adds code, label, mean period stringency, and mean period containment.

## 5. OxCGRT policy-index derivation

The production sources are:

```text
data/raw/oxcgrt/OxCGRT_timeseries_StringencyIndex_v1.csv
data/raw/oxcgrt/OxCGRT_timeseries_ContainmentHealthIndex_v1.csv
```

`load_oxcgrt_stringency(...)` and `load_oxcgrt_containment(...)` each require exactly one row with `RegionName == "Scotland"`, identify date columns using the `%d%b%Y` format, and return a long daily index table.

For policy period `p` with inclusive dates `[start_p, end_p]`:

```text
policy_stringency_p = mean(daily OxCGRT stringency_index where start_p <= date <= end_p)
policy_containment_p = mean(daily OxCGRT containment_index where start_p <= date <= end_p)
```

`derive_period_stringency(...)` and `derive_period_containment(...)` calculate these means when `utils.policy` is imported. A missing table, missing schema, or non-unique Scotland row fails explicitly rather than silently using stale constants.

## 6. Policy-index concordance analysis

`policy_index_comparison.py` aligns complete daily Stringency and Containment and Health Index observations over the defined study-policy period. Its two-column figure shows both time series and their daily scatter with a fitted linear relationship and identity reference.

The supporting `policy_index_correlation` table records complete-day count, date range, Pearson correlation, Spearman rank correlation, linear slope/intercept, and squared Pearson correlation. No conventional correlation p-value is reported because adjacent daily policy-index observations are serially dependent and should not be treated as independent replicates.

The observed strong correlation supports using Stringency Index as the selected summary policy measure in downstream descriptive work, while the differing levels and definitions mean the two indices should be described as strongly concordant rather than interchangeable.

## 7. Policy display

The figure uses the reversed `RdYlGn` colour map over stringency values from 1 to 100. Each policy interval is shown as a labelled strip and as low-alpha background shading behind daily sequence counts. Stringency is contextual annotation and is not used as a causal predictor in this package.

## 8. Lineage-frequency summaries

The analysis groups sequences into regular time periods, counts clades, converts counts to within-period frequencies, and smooths frequencies for plotting. Output tables retain both raw and plotted frequency series.

The dominant clade is the clade with the largest frequency in a period. A lineage-overtake event is recorded when the dominant clade differs from the preceding period's dominant clade. These are descriptive changes in the sequenced record and may be affected by sampling intensity.

## 9. Sequencing proportion

`sequencing_proportion_by_period` summarises the available `wn_prop_sequenced` context over the plotting periods and retains the plotted series. It describes surveillance coverage, not infection incidence.

## 10. Output contract

Every table is written as index-free CSV and parquet under `chapter_analyses/surveillance/results/tables/`:

- `clade_frequency_by_period`: unsmoothed within-period clade frequencies;
- `clade_frequency_by_period_smoothed`: frequency series used for plotting;
- `clade_counts_by_period`: total and clade-specific sequence counts;
- `clade_dominance_by_period`: dominant clade, its frequency, and previous dominant clade;
- `clade_overtake_events`: periods where the dominant clade changes;
- `sequencing_proportion_by_period`: observed and plotted sequencing proportions.
- `policy_indices_daily`: aligned daily Stringency and Containment and Health Index values;
- `policy_index_correlation`: correlation and fitted-line summary supporting index selection.

The combined surveillance figure and policy-index validation figure are written under `chapter_analyses/surveillance/results/figures/` in PNG, PDF, and TIFF formats.

## 11. Reproducibility checks

Before reporting outputs:

1. Record the processed analysis dataset and raw OxCGRT table versions.
2. Confirm exactly one Scotland row is selected from each OxCGRT index table.
3. Confirm policy-period date boundaries, daily alignment, and inclusive averaging.
4. Record the smoothing window and sequence-window stride.
5. Inspect missing dates, unmatched clades, and sequencing-proportion missingness.
6. Confirm CSV and parquet tables have identical row counts and fields.
7. Rebuild the figure from the same table-producing run.

## 12. Interpretation limits

- The sequenced record is a selected subset of infections.
- Policy periods are descriptive temporal strata and overlap epidemic, variant, behavioural, testing, and vaccination changes.
- OxCGRT stringency is an aggregate policy index, not measured adherence or individual exposure.
- Lineage dominance and overtakes depend on sequencing coverage and the chosen temporal aggregation/smoothing.
- The analysis does not estimate causal effects of policy on cases, sequences, lineage replacement, or transmission.
