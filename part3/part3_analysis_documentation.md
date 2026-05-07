# Part 3 Analysis Documentation

## 1. Overview

Part 3 characterises the associations between Scottish government COVID-19
policy restriction periods and SARS-CoV-2 genomic cluster structure. The
analysis uses the same Part 1 cluster-level dataset (Leiden resolution 0.3,
`nextclade_qc == "good"`) and the policy period definitions in
`utils/policy.py`, which encodes 16 chronological periods from emergence
(2020-03-01) through post-restriction (end 2023-05-05).

The analysis is explicitly descriptive and associational. Policy periods are
strongly confounded with variant waves, testing behaviour, and secular trends
in the epidemic; causal inference about policy effects is not possible without
a dedicated quasi-experimental design that this observational data cannot
support. The value of Part 3 is to contextualise the Part 1 and Part 2
findings within the policy timeline and to identify whether period-level
transitions carry detectable association signals in the cluster structure
data even after acknowledging the wave confounding.

## 2. Data Sources and Preparation

### 2.1 Input data

- `part1/main/cache/main_cluster_table.parquet` — the primary cluster-level
  analysis table produced by Part 1. One row per genomic cluster at the
  primary Leiden resolution (0.3). Contains 193,112 observations covering
  2020-07-14 to 2023-01-31 (window midpoint dates).
- `utils/policy.py` — 16 policy period definitions with `period_code`,
  `period_label`, `start_date`, `end_date`, and `intensity` (0–100 scale
  of restriction severity).

### 2.2 Study period

The sequence data begins on 2020-07-14, which falls within P3 (Route-map
phase 3, 2020-07-10 to 2020-10-01). Periods E0, L1, P1, and P2 pre-date the
available genomic data and are excluded from all analyses. The study window
therefore covers 12 of the 16 defined policy periods: P3 through PR.

### 2.3 Policy period assignment

Each cluster is assigned the policy period whose `[start_date, end_date]`
interval contains the cluster's window midpoint date (`wn_mid_date`). This
is performed by the `attach_period_pandas()` function in `utils/policy.py`,
a pandas-compatible wrapper around the polars `attach_period()` function.
Clusters whose window midpoint falls in an inter-period gap (none exists in
the current definitions — all 16 periods tile continuously) receive a missing
period code.

### 2.4 Derived columns

- `log_datazones`: log(clip(cluster_n_datazones, 1, ∞)) — log geographic
  dispersion (singletons get log(1) = 0).
- `is_non_singleton`: Boolean indicator, `cluster_size > 1`.
- `week_start`: ISO week start date derived from `wn_mid_date`.

## 3. Analysis Components

### 3.1 Period-level descriptive table

Computed separately for each of the 12 observed policy periods. Variables
reported:

| Variable | Definition |
|---|---|
| `n_clusters_total` | All clusters (singletons + non-singletons) |
| `n_clusters_nonsingleton` | Clusters with size > 1 |
| `pct_singleton` | 100 × (1 − n_ns / n_total) |
| `median_cluster_size` | Median cluster_size (non-singletons) |
| `iqr_cluster_size_lo/hi` | 25th/75th percentiles of cluster_size (non-singletons) |
| `median_datazones` | Median cluster_n_datazones (non-singletons) |
| `iqr_datazones_lo/hi` | 25th/75th percentiles of cluster_n_datazones (non-singletons) |
| `mean_simd_excess_discordance` | Mean simd_excess_discordance (non-singletons, valid) |
| `mean_age_excess_discordance` | Mean age_excess_discordance (non-singletons, valid) |

Saved to `part3/tables/period_descriptives.csv`.

### 3.2 Weekly aggregate series

Cluster outcomes are aggregated by ISO week (week_start). The weekly table
contains:

- **All-cluster aggregates**: total cluster count, singleton fraction.
- **Non-singleton aggregates**: median cluster_size, median log cluster_size,
  median datazones, median log datazones, mean SIMD excess discordance (valid
  non-singletons), mean age excess discordance (valid non-singletons).
- **Policy annotation**: dominant policy period code (modal value within the
  week), its label, and its intensity.

Saved to `part3/tables/weekly_summaries.csv`. This table is the primary input
for all figures.

### 3.3 Policy intensity correlation

Spearman rank correlations between the weekly dominant policy intensity and
each of the four weekly outcome variables are computed pooled across all
134 ISO weeks. These are descriptive summaries of the raw
intensity–outcome relationship before controlling for wave or calendar time.

Saved to `part3/tables/intensity_correlations.csv`.

**Observed correlations (Spearman ρ, n = 134 weeks):**

| Outcome | ρ |
|---|---|
| Median log cluster size | 0.74 |
| Median log datazones | 0.58 |
| Mean SIMD excess discordance | 0.02 |
| Mean age excess discordance | 0.59 |

The strongly positive correlations for cluster size and datazones are
predominantly driven by wave confounding: the same calendar periods that
carry high intensity (first and second lockdowns) also correspond to
B.1.177 and Alpha waves, which produced large, geographically constrained
clusters. The near-zero SIMD mixing correlation is consistent with the
Part 1 finding that SIMD deprivation mixing is not strongly associated with
the overall epidemic trajectory.

### 3.4 Interrupted time-series (ITS) analyses

#### 3.4.1 Transition selection

Three transition dates are chosen for ITS analyses on the basis that each
represents a meaningful policy shift occurring within a period of relative
variant-wave stability, minimising the most acute wave confounding:

| Transition | Date | Pre-period | Post-period | Context |
|---|---|---|---|---|
| T1-onset | 2020-10-02 | P3 (intensity 30) | T1 (intensity 55) | Restriction tightening; B.1.177 wave, pre-vaccine |
| L2→SL | 2021-04-02 | L2 (intensity 95) | SL (intensity 65) | Second lockdown lifted; Alpha-dominant, primary vaccine rollout in progress |
| NN-onset | 2021-08-09 | L0 (intensity 20) | NN (intensity 10) | Full removal of legal distancing; Delta-dominant, primary rollout largely complete |

T1-onset is the cleanest pre-vaccine within-variant test, occurring entirely
within the B.1.177 era before Alpha emerged. L2→SL occurs within the Alpha
wave but coincides with the late phase of Alpha when case numbers were
declining. NN-onset occurs within the Delta wave and captures the full legal
easing signal.

The second lockdown onset (L2, Jan 5 2021) was not chosen because it coincides
precisely with the emergence of Alpha, making variant-wave and policy confounding
inseparable. The Omicron restrictions (OM, Nov 2021) were similarly avoided
because BA.1 emergence and the policy tightening are nearly simultaneous.

#### 3.4.2 ITS window and aggregation

For each transition, a ±8 ISO-week window around the transition date is
extracted from the full cluster table. Non-singleton clusters within the window
are aggregated by ISO week. The four outcome variables are:

- `log_cluster_size`: weekly median log(cluster_size) (non-singletons)
- `log_datazones`: weekly median log(cluster_n_datazones) (non-singletons)
- `simd_excess_discordance`: weekly mean (non-singletons with valid scores)
- `age_excess_discordance`: weekly mean (non-singletons with valid scores)

#### 3.4.3 ITS model

A standard segmented OLS regression is fit for each outcome × transition
combination:

```
y_t = β0 + β1·t + β2·D_t + β3·(D_t · t) + ε_t
```

where:

- `t` = signed week offset from transition (negative = pre, 0 = first
  post-transition week)
- `D_t` = 0/1 post-transition indicator
- `D_t · t` = interaction capturing the slope change after transition

Parameter interpretations:
- β0: baseline level at the transition week under the pre-transition trend
- β1: pre-transition weekly trend
- β2: level change at transition (primary quantity of interest)
- β3: slope change after transition

Coefficients, 95 % CIs (OLS normal approximation), and p-values are saved
to `part3/tables/its_coefficients.csv`. Per-transition weekly data are saved
to `part3/tables/its_weekly_{label}.csv`.

#### 3.4.4 ITS results

**Level-change estimates (β2) with 95 % CI:**

| Transition | Outcome | β2 | 95 % CI | p |
|---|---|---:|---|---|
| T1-onset | log cluster size | −0.13 | [−0.34, +0.09] | 0.22 |
| T1-onset | log datazones | −0.07 | [−0.35, +0.21] | 0.59 |
| T1-onset | SIMD excess mixing | +0.02 | [−0.03, +0.07] | 0.49 |
| T1-onset | Age excess mixing | +0.00 | [−0.04, +0.05] | 0.84 |
| L2→SL | log cluster size | **−0.21** | [−0.41, −0.02] | **0.034** |
| L2→SL | log datazones | **−0.36** | [−0.60, −0.12] | **0.006** |
| L2→SL | SIMD excess mixing | +0.01 | [−0.03, +0.04] | 0.72 |
| L2→SL | Age excess mixing | +0.01 | [−0.01, +0.03] | 0.27 |
| NN-onset | log cluster size | +0.12 | [−0.01, +0.26] | 0.068 |
| NN-onset | log datazones | **+0.32** | [+0.07, +0.56] | **0.015** |
| NN-onset | SIMD excess mixing | −0.01 | [−0.06, +0.04] | 0.75 |
| NN-onset | Age excess mixing | −0.01 | [−0.04, +0.03] | 0.68 |

*Starred (bold) entries: p < 0.05.*

## 4. Interpretation Notes

### 4.1 T1-onset (Oct 2020): null signal

The introduction of a tiered restriction framework in October 2020 shows no
significant change in cluster size, geographic spread, or mixing metrics. This
may reflect that the actual restriction *experience* changed little in practice
between P3 and T1 (both had community gathering restrictions), or that any
effect on transmission takes several weeks to appear in genomic clusters
(reflecting earlier transmission events). The B.1.177 wave was mature by
October 2020 and cluster structure may have been more driven by the
heterogeneous national spread of the lineage than by the policy change per se.

### 4.2 L2→SL (Apr 2021): counter-intuitive reduction in cluster scale

The significant negative level changes for cluster size and datazones at the
L2→SL transition (lockdown lift) are initially counter-intuitive: why would
clusters get smaller as restrictions ease? Several non-exclusive explanations
apply:

1. **Alpha wave tail confounding**: The L2→SL transition occurs at the end of
   the Alpha wave peak. Case numbers and sequencing volume were declining
   independently of the policy change. Smaller sequential clusters late in a
   wave are expected under natural epidemic dynamics.
2. **Behavioural lag**: The Stay-local instruction (SL, Apr 2–25) is a very
   short period (24 days). Population mixing patterns may not have changed
   substantially over this interval.
3. **Vaccination context**: Primary vaccination was actively rolling out for
   50–65 year olds during this period. Reduced household secondary attack
   rates in vaccinated contacts would compress cluster sizes.

No causal interpretation of this finding is warranted without a counterfactual
study design.

### 4.3 NN-onset (Aug 2021): geographic dispersal increases

The significant positive level change in log datazones at the NN-onset
transition (L0 → Near-normal, full legal restrictions lifted) is consistent
with a priori expectations: removal of gathering limits, physical distancing
requirements, and capacity limits would allow transmission events spanning
more geographic areas. This is the most interpretable of the three ITS results,
though it coincides with the Delta wave at near-peak incidence (Delta peaked
in Scotland in late August–early September 2021), again limiting causal
attribution to the policy change alone.

### 4.4 Mixing metrics

None of the three transitions produced significant changes in SIMD or age
excess discordance at conventional levels. This is consistent with the Part 1
finding that demographic mixing patterns were more strongly associated with
structural and surveillance factors than with the acute policy context.

## 5. Limitations

1. **Wave confounding**: All three transitions occur within variant waves that
   independently drive cluster structure. ITS estimates conflate policy effects
   with wave dynamics and cannot be cleanly separated.
2. **Short observation windows**: ±8 weeks yields 16 weekly observations per
   transition, providing limited power for the ITS models.
3. **Ecological (cluster-level) analysis**: ITS outcomes are weekly medians of
   cluster-level statistics, not individual-level outcomes. Regression to the
   mean and aggregation artefacts apply.
4. **No comparison group**: A synthetic control or difference-in-differences
   design (e.g., using a comparator region with different policy timing) would
   be required for stronger causal inference. No such comparator is available
   within Scotland.
5. **Testing and sequencing behaviour**: Policy changes alter testing behaviour
   (people may test more or less after a policy announcement), which can
   change the sequencing denominator independently of transmission.

## 6. File Manifest

| Path | Description |
|---|---|
| `part3/part3_analysis.py` | Main analysis script |
| `part3/manuscript/make_figures.py` | Figure generation script |
| `part3/tables/period_descriptives.csv` | Per-period cluster outcome summary |
| `part3/tables/weekly_summaries.csv` | ISO-week cluster outcomes + policy annotation |
| `part3/tables/intensity_correlations.csv` | Spearman ρ: intensity vs outcomes |
| `part3/tables/its_coefficients.csv` | ITS OLS coefficients, CIs, p-values |
| `part3/tables/its_weekly_{label}.csv` | Per-transition weekly ITS data (3 files) |
| `part3/manuscript/figures/` | PDF, PNG, TIFF figure outputs |
| `utils/policy.py` | Policy period definitions + pandas helpers |
