# Paper 2 — Demographic and vaccination correlates of SARS-CoV-2 transmission clusters in Scotland

## One-line scope

How the age distribution, sex composition, and vaccination status of sequenced cases relate to the size and within-cluster homogeneity of SARS-CoV-2 transmission clusters in Scotland (Jul 2020 – Feb 2023), and how those relationships shifted across the Alpha → Delta → Omicron waves.

## Primary hypotheses

1. Clusters with a younger median age are larger than those with an older median age, after adjusting for VOC and surveillance intensity, reflecting higher-contact age groups driving onward transmission.
2. A higher proportion of vaccinated members within a cluster is associated with smaller cluster size, with the effect decomposing into pre-Omicron (strong) and Omicron (attenuated) eras.
3. Within-cluster age diversity (std of age midpoints) increases with cluster size, consistent with household-to-household spread broadening the age signature of the chain.

## Unit of analysis

One row per `(window_id, cluster_id)` at Leiden resolution 0.5. Uses `cluster_demographic_features.parquet` pre-computed by `analysis/demographic.py`, joined to the window-level `wn_prop_sequenced` from the master dataset.

## Figure list

| # | File                                        | Description |
|---|---------------------------------------------|-------------|
| 1 | `figures/fig1_age_over_time_by_epoch.py`    | Ridgeline/density of median cluster age by VOC epoch |
| 2 | `figures/fig2_vaccination_vs_cluster_size.py` | Vaccination prevalence vs. cluster size (binned) with bootstrap CI bands, per epoch |
| 3 | `figures/fig3_age_homogeneity.py`           | Within-cluster age std vs. cluster size; scatter + hex-binned density |
| 4 | `figures/fig4_sex_composition.py`           | Sex ratio (`frac_female`) vs. cluster size, stratified by VOC |
| 5 | `figures/fig5_demographic_forest.py`        | Forest plot of demographic IRRs / ORs on cluster size and on singleton status |
| 6 | `figures/fig6_voc_stratified_shifts.py`     | Paired panels showing the shift in demographic-size association across Alpha, Delta, Omicron |

## Statistical conventions (read once, then skim the figure notes)

- **IRR (incidence rate ratio)** — exponentiated coefficient from a negative-binomial GLM of `n_sequences` on the standardised predictor. Because predictors are standardised, an IRR = 1.15 means **"a 1-SD increase in the predictor is associated with 15% larger expected cluster size, holding other predictors constant"**. The `log(wn_prop_sequenced)` offset makes effects "per sequenced case" rather than "per observed case".
- **OR (odds ratio)** — exponentiated coefficient from a logistic GLM of `is_singleton` on the same standardised predictor set. An OR = 0.80 on `frac_vaccinated` means **"a 1-SD increase in the cluster's vaccinated fraction is associated with 20% lower odds of being a singleton"** (i.e. that cluster is *more* likely to be part of an onward-transmission chain).
- **Per 1-SD.** Every continuous demographic predictor (`median_age`, `frac_female`, `frac_vaccinated`, `mean_vacc_dose`) is Z-scored within the model frame before fitting. IRRs and ORs are therefore directly comparable across predictors — a bigger |log(estimate)| means a stronger effect in standardised units. Raw per-year or per-percentage-point effects can be recovered by dividing the log estimate by the column's SD.
- **95% CI** (`conf_low`, `conf_high`) are Wald intervals back-transformed; a CI that does not cross 1 is the usual α = 0.05 threshold.
- **Reference level** for `who_voc` is **Omicron** (combined BA.1 / BA.2+); VOC dummies appear in tables as `voc_Alpha`, `voc_Delta`, `voc_Pre-VOC`, etc.
- **Binomial-style SE** on mean `frac_female` is computed as `SD / √N` per bin; bands show ±1.96·SE — not bootstrap.
- **Bootstrap CIs** (Fig. 2) are percentile intervals from 800 resamples of the within-bin cluster-level `frac_vaccinated` values.
- **Age-homogeneity null** (Fig. 3) draws 500 random same-size subsets of `age_midpoint` from the entire epoch's sequence pool and takes the mean of their SDs — "what size-vs-SD curve would we see if cluster membership were independent of age?".
- **Log-scale axes** are used on all forest plots; the dashed reference line is at IRR = OR = 1.

## How to read each figure

### Fig. 1 — Median cluster age shifts across VOC epochs (ridgeline)

One ridge per VOC epoch, stacked vertically, all sharing the same x-axis (median age of cluster members, 0–100 years). Each ridge is a Gaussian KDE of the per-cluster **median ages**, normalised so every ridge peaks at the same height (the plot shows distribution *shape*, not *density*, across epochs). A short black tick inside each ridge marks the median of medians; the inline label reports `n` clusters and the epoch median. **Read it for shifts in modal age, not area.** Lower median under Omicron vs. Alpha is consistent with the school- and young-adult-driven Omicron wave; a bimodal ridge is the tell for simultaneous household and workplace transmission. Singletons are excluded (the median of one value isn't informative).

### Fig. 2 — Vaccination prevalence vs. cluster size, by epoch

X-axis: cluster size on a log scale, discretised into log-spaced bins per epoch. Y-axis: mean **`frac_vaccinated`** within the bin (values between 0 and 1). The band is a percentile bootstrap 95% CI over the bin's per-cluster `frac_vaccinated`. A **downward-sloping curve within an epoch ⇒ larger clusters are composed of less-vaccinated members**, the Paper 2 hypothesis. Compare epochs by vertical offset: Pre-VOC and Alpha curves sit near zero (no vaccine programme yet); Delta is the sweet spot for the "vaccine reduces onward transmission" signal; Omicron curves flatten towards the population mean (saturation). Companion table `tables/fig2_vacc_vs_size.csv` has one row per (epoch × size-bin) with `size_mid`, `mean`, `lo`/`hi` (bootstrap 95%), `n`.

### Fig. 3 — Within-cluster age diversity vs. cluster size

A hexbin density of clusters on (log1p(size), age-std), overlaid with (i) a running observed median with a 10th–90th percentile band and (ii) a dashed **random-draw null**: the expected age SD if a cluster of the same size were a random sample from that epoch's age pool. Clusters whose age-diversity sits **below** the null curve are *more age-homogeneous* than chance — the signature of household / nursing-home / school-cohort transmission. Clusters above the null are *more age-diverse* than chance, suggestive of mixed-contact settings. The hex colour scale (viridis) is cluster count per hex — use it to weight visual impressions toward the denser (smaller-cluster) regions.

### Fig. 4 — Sex composition vs. cluster size, by epoch

Each cluster contributes one dot (grey cloud, faint) at (cluster size, fraction female). Per-epoch lines with ±1.96·SE bands are rolling means within log-spaced size bins — the SE assumes within-bin binomial-style variation and will be slightly optimistic when clusters are non-independent. The horizontal dashed line at **0.51** is the Scotland population female share. **Systematic deviation from 0.51 in larger clusters ⇒ a setting-specific bias** (e.g. female-skewed care-home outbreaks, male-skewed workplace outbreaks); a line that hugs 0.51 means no aggregate sex signal. Expect Pre-VOC / Alpha in care-home settings to run above 0.51 at size > 50; Delta's return-to-work outbreaks to run below.

### Fig. 5 — Forest plot of demographic predictors

Two side-by-side forest panels, both with a log-scaled x-axis and a dashed null line at 1.00. **Left panel (IRR per 1-SD, cluster size)** — the NB GLM `n_sequences ~ demographics + VOC + offset(log wn_prop_sequenced)`. **Right panel (OR per 1-SD, singleton)** — the logistic GLM on `is_singleton`. Each row is one standardised demographic predictor (median age, frac. female, frac. vaccinated, mean vacc. dose). **Left panel rule of thumb:** IRR > 1 ⇒ higher value of the predictor is associated with *larger* clusters; IRR < 1 ⇒ *smaller*. **Right panel rule of thumb:** OR < 1 ⇒ higher value of the predictor is associated with being *less* likely to be a singleton (= more likely to be in a transmission chain). Companion tables: `tables/fig5_irr_size.csv` (NB) and `tables/fig5_or_singleton.csv` (logit). Columns: `term`, `estimate` (IRR or OR), `std_error` (log-scale SE), `conf_low`/`conf_high` (95% Wald CI), `z`, `p_value`.

### Fig. 6 — Within-epoch shifts of demographic IRRs on cluster size

Four side-by-side panels, one per predictor. In each panel, the y-axis lists VOC epochs; the x-axis shows the IRR (log scale) from **a separate NB GLM fit on clusters from that epoch only**, with the predictor's SE bars. Dashed reference line at IRR = 1. **Read it for trajectories, not levels.** A predictor whose IRR drifts from >1 (Alpha) towards 1 (Omicron) is consistent with effect-attenuation as population-level immunity saturates (Paper 2 hypothesis 2). Overlapping CIs across two epochs do not mean "no difference" — a formal interaction test would, but the figure is intended as an exploratory visual. Companion table `tables/fig6_voc_stratified_irrs.csv` columns: `term`, `epoch`, `estimate`, `std_error`, `conf_low`/`conf_high`, `z`, `p_value`, `n_obs`.

## Inputs

- `data/processed/scotland_clustering_analysis_dataset.parquet`
- `data/processed/cluster_demographic_features.parquet`
- `data/processed/cluster_summary.parquet`

## Running

```bash
python -m manuscripts.paper2_demographic.make_figures --figures manuscripts/paper2_demographic/figures
```

## Statistical models

`models/cluster_demographics.py`:

- **Cluster-size NB GLM** — `n_sequences ~ median_age + frac_female + frac_vaccinated + mean_vacc_dose + voc`, offset `log(wn_prop_sequenced)`.
- **Singleton logistic GLM** — same predictors, outcome `is_singleton`.
- **Age diversity linear GLM** — `age_diversity ~ log(n_sequences) + median_age + voc`.

## Target journals

*Virus Evolution*, *Eurosurveillance*, *Lancet Microbe*.

## Known caveats

1. Cluster-level aggregates of age and vaccination status are ecological; individual-level inference is not claimed.
2. The vaccination programme's age-staggered roll-out means `mean_vacc_dose` and `median_age` are strongly confounded — their joint interpretation requires an interaction term with calendar time.
3. Omicron cohort exposure is saturated; many Omicron clusters have `frac_vaccinated` near the population mean, reducing variance for inference.
