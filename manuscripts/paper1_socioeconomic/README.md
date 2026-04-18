# Paper 1 — Socioeconomic deprivation and SARS-CoV-2 transmission clustering in Scotland

## One-line scope

Whether area-level deprivation, measured by the 2020 Scottish Index of Multiple Deprivation (SIMD), predicts the size, persistence, and composition of SARS-CoV-2 transmission clusters inferred from COG-UK sequences in Scotland (Jul 2020 – Feb 2023), and whether that relationship changed across the Alpha → Delta → Omicron waves.

## Primary hypothesis

Clusters whose members predominantly reside in more deprived data zones are larger, less likely to be singletons, and have a broader geographic footprint than clusters drawn from less deprived areas, after adjusting for surveillance intensity (`wn_prop_sequenced`) and lineage. We expect the magnitude of the deprivation gradient to be strongest during Alpha and Delta and to attenuate during Omicron, when widespread population infection saturated the deprivation signal.

## Unit of analysis

One row per `(window_id, cluster_id)` at a single primary Leiden resolution (0.3). Sensitivity analyses sweep the full resolution grid `[0.1 … 0.8]`.

## Figure list

| # | File                                       | Description |
|---|--------------------------------------------|-------------|
| 1 | `figures/fig1_sequences_by_simd_over_time.py` | Weekly sequence counts stratified by SIMD quintile, with epoch shading and a `wn_prop_sequenced` rug |
| 2 | `figures/fig2_cluster_size_by_simd.py`     | Violin/box of cluster size distribution by SIMD quintile, faceted by VOC epoch |
| 3 | `figures/fig3_simd_domain_forest.py`       | Forest plot of IRRs for cluster size by SIMD domain (overall, income, employment, health, …) |
| 4 | `figures/fig4_singleton_odds_by_epoch.py`  | Odds ratio of singleton vs. non-singleton by SIMD quintile, one panel per VOC epoch |
| 5 | `figures/fig5_deprivation_lineage_heatmap.py` | Heat-map of median SIMD rank over time × top lineages |
| 6 | `figures/fig6_domain_decomposition.py`     | Proportion of the overall SIMD–cluster-size effect attributable to each SIMD domain via a variance-decomposition bar chart |

## Statistical conventions (read once, then skim the figure notes)

- **IRR (incidence rate ratio)** — the ratio of expected cluster sizes under a 1-unit change in the predictor (with all other predictors held constant) from a negative-binomial GLM whose log link makes effects multiplicative. We exponentiate the fitted coefficient so "estimate" is already an IRR. **IRR = 1.00** ⇒ no association; **IRR > 1** ⇒ larger clusters; **IRR < 1** ⇒ smaller clusters. A 95% CI that straddles 1 is not statistically distinguishable from "no effect" at α = 0.05.
- **OR (odds ratio)** — the multiplicative change in the odds of a binary outcome (here `is_singleton` = 1) from a logistic GLM. **OR < 1** on `is_singleton` means the predictor group is *less* likely to be a singleton (i.e. *more* likely to be part of a genetically linked onward chain).
- **Per 1-SD** — when a predictor is standardised before fitting, the IRR/OR is the effect of a one-standard-deviation increase in that predictor. Signed such that **positive = more deprived** for SIMD predictors.
- **Offset** — cluster-size NB GLMs include `log(wn_prop_sequenced)` as an offset, so an IRR of 1.30 means *30% more sequenced cases per sequenced case*, i.e. the effect is per sampling effort rather than per case observed.
- **95% CI** (`conf_low`, `conf_high`) are Wald normal-approximation intervals on the log scale, then back-transformed. `p_value` is the two-sided Z-test on the coefficient's log-scale estimate.
- **Log-scaled axes** are used on every forest plot; the reference vertical dashed line sits at 1.00 (the null IRR/OR).
- **Kruskal-Wallis p-values** (Fig. 2) test "are cluster-size distributions identical across SIMD quintiles within this epoch?" without assuming normality.
- **Reference levels.** `voc` uses Omicron as the reference; SIMD quintile forests use Q5 (least deprived) as the reference.

## How to read each figure

### Fig. 1 — Sequences per week by SIMD quintile, plus surveillance intensity

Two stacked panels sharing an x-axis. **Panel A** — weekly sequenced case counts, one line per SIMD quintile (Q1 = most deprived, red-hued; Q5 = least deprived, green-hued). Shaded vertical bands mark VOC epochs. Read this to check (i) whether the deprivation-quintile ordering is visually stable or reshuffles between waves, and (ii) whether any quintile is under-represented in a given window — that bias propagates into every downstream regression. **Panel B** — proportion of positive SARS-CoV-2 tests sequenced (`wn_prop_sequenced`) for that week. This is a surveillance-intensity context panel, not a regression output; low values (e.g. <1%) mean windows whose cluster-size inference is more fragile, and are why every NB GLM carries a `log(wn_prop_sequenced)` offset.

### Fig. 2 — Cluster size distribution by SIMD quintile, faceted by VOC epoch

Five violins per panel (one per SIMD quintile), one panel per VOC epoch. The y-axis is `log1p(n_sequences)`, not raw size — cluster-size distributions are heavily skewed, with singletons (y = 0) dominating and a long right tail. The notched boxplot inside each violin shows median and IQR; non-overlapping notches across quintiles hint at a credible median difference. The subtitle on each panel reports the Kruskal-Wallis p-value testing **"are the five quintile distributions all drawn from the same distribution within this epoch?"** A small p (p < 0.001) confirms some quintile differs, but doesn't by itself tell you it's Q1; look for a monotonic ordering of the violins for a deprivation-gradient story.

### Fig. 3 — Forest plot of IRRs for cluster size by SIMD domain

One row per SIMD domain (overall + 7 component domains). Each row shows the **IRR per 1-SD increase in deprivation** (sign-flipped, so more deprivation is positive) from a separate NB GLM that also adjusts for VOC and the log-sequencing-proportion offset. **IRR > 1 ⇒ more deprivation is associated with larger clusters.** The dashed vertical line is the null (IRR = 1). Rows are sorted by effect size for legibility. The right-margin annotation reports the `n` complete cases and the p-value. Companion table `tables/fig3_domain_irrs.csv` columns: `domain` (index), `estimate` (IRR), `conf_low`/`conf_high` (95% Wald CI), `p_value` (two-sided Z), `n` (complete cases). Rule of thumb: the "overall" row should sit between the income / employment rows (dominant signals) and the housing / crime rows (noisier signals); if it does not, the aggregate SIMD rank is being pulled by one domain worth naming in the Discussion.

### Fig. 4 — Singleton-cluster odds by SIMD quintile, per VOC epoch

One panel per VOC epoch. Each panel shows, from a logistic GLM that also adjusts for the standardised log-sequencing-proportion, the **OR of being a singleton** (vs. a multi-member cluster) for Q1–Q4 relative to **Q5** (least-deprived, not plotted because it is the reference). X-axis has a vertical dashed line at OR = 1. **OR < 1 for Q1 ⇒ the most deprived quintile is *less* likely to be a singleton, i.e. their sequences are more likely to sit inside a genetically linked onward chain** — the epidemiologically interesting direction for Paper 1. Epochs where Q1's point lies close to 1 with a wide CI are consistent with the hypothesis that widespread population infection (Omicron era) attenuates the deprivation gradient. Companion table `tables/fig4_singleton_ors.csv`: `epoch`, `quintile`, `estimate` (OR vs Q5), `conf_low`/`conf_high`, `p_value`, `z`, `std_error`.

### Fig. 5 — Median SIMD rank over time × top lineages (heat-map)

Rows are the top-20 PANGO lineages by count; columns are monthly bins. Cell colour = median **SIMD rank percentile** of sequences assigned to that (lineage, month), where **0 = most deprived, 1 = least**; the diverging palette is centred at the Scotland median (0.5). Empty cells mean the lineage was not sampled that month. The right-margin bar shows each lineage's total sequence count on a log scale — use it to temper claims about rows with little mass. A lineage-row that starts red (most-deprived communities) and migrates to blue (least-deprived) over its lifetime is the typical signature of a variant seeded in deprived areas and diffusing upward through the SIMD gradient.

### Fig. 6 — SIMD-domain decomposition in a mutually-adjusted model

A single NB GLM fitted with **all seven SIMD domain ranks simultaneously** plus VOC and offset. We plot each domain's share of `Σ |standardised coefficient|` (as a percent) — i.e. "when every domain is allowed to compete, what fraction of the signal does each carry?". The bar label shows the signed standardised coefficient (not exponentiated — this is on the log-IRR scale) and the p-value. **A large share with a p-value < 0.05 ⇒ that domain is doing independent work**; domains with large shares but non-significant p-values are collinear with another large-share domain. Shares are normalised so the seven bars sum to 100%. Companion table `tables/fig6_domain_decomposition.csv`: `domain`, `estimate` (standardised log-IRR coefficient), `std_error`, `conf_low`/`conf_high` (CI on the log-IRR scale), `p_value`, `abs_std_coef`, `share` (sums to 1).

## Inputs

- `data/processed/scotland_clustering_analysis_dataset.parquet`
- `data/processed/cluster_summary.parquet`
- `data/processed/cluster_simd_features.parquet` (auto-derived on first use if missing)

## Running

```bash
# From the repository root, with analysis scripts already run:
python -m manuscripts.paper1_socioeconomic.make_figures --output manuscripts/paper1_socioeconomic/output

# Or one at a time:
python -m manuscripts.paper1_socioeconomic.figures.fig1_sequences_by_simd_over_time
```

Outputs are written to `output/` in both PDF (for submission) and PNG (for quick review) formats.

## Statistical models

All regression specifications live in `models/simd_models.py`; every figure that needs a fit pulls it from there so the three NB panels (Figs 3, 5-companion, 6) and the one logistic panel (Fig 4) share design-matrix conventions.

- **Cluster size (headline)** — `cluster_size_model`. Negative-binomial GLM on `n_sequences` with predictors **SIMD quintile dummies** (Q1–Q4 vs. Q5 reference) *or* a standardised continuous rank when `deprivation_measure != "simd_quintile_mode"`, **VOC dummies** (Omicron reference), and a **natural cubic regression spline** (`patsy.cr`, 5 df) on days-since-study-start built from `wn_mid_date`. Offset is `log(wn_prop_sequenced)` (applied inside `stats.negbin_cluster_size`).
- **Per-domain forest** — `build_domain_forest_table`. One NB GLM per SIMD domain: `n_sequences ~ deprivation_sd (1-SD standardised, sign-flipped so +1 SD = more deprived) + VOC + cr(wn_mid_date) + offset(log wn_prop_sequenced)`.
- **Domain decomposition (joint)** — `build_domain_decomposition_table`. A single NB GLM with **all seven** domain ranks as standardised predictors simultaneously, plus VOC and spline. Returns standardised log-IRR coefficients and their share of `Σ|coef|`.
- **Singleton status** — `build_singleton_epoch_table`. Logistic GLM fit **separately within each VOC epoch** (VOC dummies would collapse) on `is_singleton ~ Q1..Q4 (vs Q5) + log_prop_seq_z + cr(wn_mid_date, df=3)`. The smaller df is because an epoch's temporal span is shorter than the full study window.

All NB / logistic fits are produced by `stats.negbin_cluster_size` and `stats.logit_singleton`, and all tidy tables are produced by `stats.tidy_glm` (Wald 95% CIs, two-sided Z p-values, IRR/OR exponentiation unless explicitly suppressed for the decomposition shares).

## Target journals

*Lancet Regional Health Europe*, *International Journal of Epidemiology*, *SSM – Population Health*. STROBE and a tailored genomic-surveillance methods supplement are expected.

## Known limitations to state up front

1. Sequencing coverage (`wn_prop_sequenced`) is not uniform across SIMD quintiles; any result must be shown to be robust to this as a covariate or inverse-probability weight.
2. Cluster assignments are genomic proxies for transmission, not observed transmission events.
3. Sliding windows overlap by two weeks, so the primary analysis uses one canonical window per sequence (the one whose mid-date is closest to `collection_date`) to avoid pseudo-replication.
