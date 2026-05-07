# Part 1 Main Results And Figure Guide

Generated for the main Part 1 analysis on 5 May 2026.

## Analysis Question

Part 1 asks whether socioeconomic deprivation and local surveillance conditions
are associated with SARS-CoV-2 genomic cluster characteristics in Scotland after
accounting for lineage, calendar time, local incidence, sequencing intensity,
and test positivity, using one primary Leiden clustering resolution.

The specific outcomes are:

- cluster size
- geographic spread, measured by the number of distinct datazones in a cluster
- observed-minus-expected socioeconomic and demographic mixing within clusters

The main manuscript figures are generated from `part1/main/tables` and
`part1/main/cache`, using the publication styling in `utils/style.py`.

## Analysis Population And Modelling Frame

The main analysis uses good-QC sequences at Leiden resolution 0.3. This produces
789,347 sequence rows and 193,112 inferred genomic clusters across 134 time
windows. Rare Pango lineages are pooled at a threshold of 50 clusters, giving
183 lineage model levels in the main pooled analysis.

The modelled count outcomes have large structural masses at their minimum
values, so the main count analysis uses two-part hurdle models:

- a binomial hurdle component for whether a cluster exceeds the structural
  minimum
- a zero-truncated negative-binomial positive-count component for the outcome
  among clusters exceeding the structural minimum

For cluster size, the hurdle is `cluster_size > 1` and the positive count is
`cluster_size - 1`. For geographic spread, the hurdle is
`cluster_n_datazones > 1` and the positive count is `cluster_n_datazones - 1`.

Cluster duration remains in the descriptive cluster summaries and the
supplementary outcome-distribution figure, but it is not modelled in the
current main count analysis because the fixed three-week clustering windows
mechanically constrain the observed span.

Mixing outcomes are fitted among non-singleton clusters. They are defined as
observed pairwise discordance within a cluster minus expected discordance among
sampled cases from the same lineage and calendar window. Positive values mean
more mixing than expected; negative values mean more homogeneity than expected.

Throughout the figures, higher `SIMD deprivation` means a more deprived cluster
on the transformed and standardised SIMD rank scale.

## Overall Results Summary

The main analysis does not support the simple hypothesis that more deprived
clusters are generally larger and more geographically dispersed.
In the pooled main models, higher SIMD deprivation is associated with slightly
lower odds of being non-singleton, smaller positive cluster sizes, and
substantially lower positive geographic spread.

Surveillance and epidemic-intensity variables show much stronger and more
consistent associations with cluster outcomes. Higher local incidence, higher
window-level sequencing proportion, and higher test positivity are associated
with larger and more geographically dispersed clusters. Local sequencing
fraction is especially strongly associated with larger positive cluster size and
positive geographic spread, consistent with the idea that local sampling
intensity affects the apparent scale of reconstructed genomic clusters.

Mixing results tell a different story. Overall SIMD deprivation is not clearly
associated with SIMD-quintile excess mixing, but it is associated with greater
age mixing and greater joint socioeconomic-demographic profile mixing, and with
lower sex mixing. Different SIMD domains show different patterns: education and
crime deprivation tend to increase domain-quintile mixing, whereas access and
housing deprivation are associated with lower domain-quintile mixing.

Wave-specific results are heterogeneous. Delta shows the clearest negative
association between deprivation and cluster size/geographic spread. BA.2 and
BA.4 show positive associations for positive cluster size and geographic spread,
but BA.4 has a much smaller sample and should be treated cautiously. XBB is
included in descriptives but not in wave-specific regression models because it
falls below the minimum cluster count.

---

## Main Figure 1: Main Cluster Outcome Models

**File:** `fig1_main_cluster_outcomes`

**What it shows.** A coefficient plot of adjusted ratios from the pooled hurdle
and zero-truncated negative-binomial (ZTNB) models for cluster size and
geographic spread. The four panels correspond to the cluster size hurdle (odds
ratio), positive cluster size (ZTNB count ratio), geographic spread hurdle (odds
ratio), and positive geographic spread (ZTNB count ratio). Each point is an
adjusted ratio per 1 SD higher covariate, with 95% confidence intervals. The
main exposure (SIMD deprivation) and surveillance and incidence covariates
(local incidence, local sequencing fraction, window sequencing proportion, test
positivity) are all shown alongside lineage and calendar time adjustments.

**Key visual patterns.** The SIMD deprivation point falls clearly below 1.0 for
positive cluster size (count ratio 0.926) and positive geographic spread (count
ratio 0.851), and modestly below 1.0 for the cluster size hurdle (OR 0.971),
while the geographic spread hurdle estimate sits close to 1.0. Surveillance
covariates far dominate the SIMD signal: test positivity shows count ratios of
2.65 for positive cluster size and 3.00 for positive geographic spread; local
sequencing fraction shows count ratios of 3.24 and 2.27 respectively. Higher
local incidence and window sequencing proportion are also consistently and
visibly associated with larger and more dispersed clusters.

**Suggested results paragraph.**
> In the pooled hurdle/ZTNB models, higher cluster-level SIMD deprivation was
> associated with slightly lower odds of being non-singleton (OR 0.971, 95% CI
> 0.960–0.983, p < 0.001) and smaller positive cluster size (count ratio 0.926,
> 95% CI 0.869–0.987, p = 0.018), but was not associated with geographic spread
> at the hurdle stage (OR 1.004, 95% CI 0.992–1.016, p = 0.522). Among clusters
> exceeding the structural minimum, positive geographic spread was substantially
> lower with higher deprivation (count ratio 0.851, 95% CI 0.792–0.915,
> p < 0.001). By contrast, local sequencing fraction (count ratio 3.24 for
> positive cluster size; 2.27 for positive geographic spread), test positivity
> (2.65 and 3.00 respectively), and window-level sequencing proportion showed far
> larger positive associations with cluster scale, consistent with surveillance
> and epidemic intensity being the dominant structural determinants of observed
> cluster size and geographic spread.

---

## Main Figure 2: Main Cluster Mixing Models

**File:** `fig2_main_cluster_mixing`

**What it shows.** A coefficient plot of adjusted changes in excess mixing, in
percentage points per 1 SD higher covariate. The four panels show results for
SIMD-quintile excess mixing, age-band excess mixing, sex excess mixing, and
joint SIMD-age-sex profile excess mixing, all modelled among non-singleton
clusters. Each point is an adjusted change in excess discordance (observed minus
expected within lineage-calendar stratum) per 1 SD covariate increase, with
95% confidence intervals.

**Key visual patterns.** The SIMD deprivation estimate for age excess mixing
(+1.66 pp) is the largest and most precisely estimated positive estimate across
the figure. Sex excess mixing shows a clear negative SIMD estimate (−0.78 pp).
The SIMD-quintile mixing estimate straddles zero (+0.31 pp), while joint profile
mixing is positive but smaller (+0.48 pp). Cluster size is a visible driver of
mixing: larger clusters show consistently higher SIMD, age, and joint-profile
mixing but lower sex mixing. Local incidence is positively associated with SIMD,
age, and joint-profile mixing. Sequencing and test-positivity covariates
generally push mixing estimates downward.

**Suggested results paragraph.**
> Higher cluster-level SIMD deprivation was not significantly associated with
> SIMD-quintile excess mixing (+0.31 pp, 95% CI −0.18 to +0.80, p = 0.208), but
> was positively associated with age excess mixing (+1.66 pp, 95% CI +1.29 to
> +2.03, p < 0.001) and joint SIMD-age-sex profile mixing (+0.48 pp, 95% CI
> +0.29 to +0.67, p < 0.001), and negatively associated with sex mixing (−0.78
> pp, 95% CI −1.16 to −0.39, p < 0.001). Cluster size was the strongest
> cross-outcome driver of excess mixing: larger clusters showed higher SIMD, age,
> and joint-profile mixing but lower sex mixing. Surveillance covariates
> (sequencing proportion, test positivity) generally reduced mixing estimates,
> consistent with sampling conditions influencing the observed within-cluster
> demographic composition.

---

## Main Figure 3: SIMD-Domain Mixing

**File:** `fig3_simd_domain_mixing`

**What it shows.** A four-panel coefficient plot comparing the overall SIMD
deprivation measure with each individual SIMD domain (income, employment,
education, health, crime, access, housing) for domain-quintile mixing (A), age
mixing (B), sex mixing (C), and joint age-sex profile mixing (D). Each point is
an adjusted change in excess discordance per 1 SD domain-deprivation increase
among non-singleton clusters, with 95% confidence intervals.

**Key visual patterns.** In panel A, two domains (access, −2.08 pp; housing,
−1.19 pp) are clearly left of zero, while education (+1.17 pp) and crime
(+1.09 pp) are clearly right of zero — a contrast that spans roughly 3 pp and
is the dominant visual feature of the panel. The overall SIMD estimate sits
between them, close to zero. In panel B (age mixing), most domain estimates
including overall SIMD are positive, with access deprivation as an outlier in
the negative direction. Panel C (sex mixing) shows predominantly negative
estimates across domains, broadly in line with the main SIMD result. Panel D
(joint age-sex mixing) is mostly positive with more variable magnitudes.

**Suggested results paragraph.**
> SIMD-domain models revealed heterogeneity in domain-quintile mixing that the
> overall SIMD measure conceals. Education deprivation (+1.17 pp, 95% CI +0.64
> to +1.69) and crime deprivation (+1.09 pp, 95% CI +0.55 to +1.64) were
> associated with greater within-cluster mixing across domain quintiles, whereas
> access deprivation (−2.08 pp, 95% CI −2.50 to −1.66) and housing deprivation
> (−1.19 pp, 95% CI −1.77 to −0.62) were associated with less mixing — a
> direction opposite to education and crime. Demographic mixing patterns (panels
> B–D) were broadly consistent with the main SIMD estimates: most domains showed
> positive age mixing and negative sex mixing, with access deprivation again
> behaving as an outlier, showing lower age mixing in contrast to all other
> domains.

---

## Main Figure 4: Wave-Specific SIMD Effects On Cluster Outcomes

**File:** `fig4_wave_specific_cluster_outcomes`

**What it shows.** A multi-wave coefficient plot showing the per-wave overall
SIMD deprivation effect on cluster size and geographic spread. Each wave
contributes four estimates — cluster size hurdle (OR), positive cluster size
(count ratio), geographic spread hurdle (OR), and positive geographic spread
(count ratio) — fitted using the same hurdle/ZTNB framing as the pooled main
analysis. Wave-specific models adjust for the same covariates as the main model;
within-wave lineage dummies are included where estimable.

**Key visual patterns.** Delta stands out visually with all four estimates
clearly below their null values: cluster size hurdle OR ~0.93, positive cluster
size count ratio ~0.80, geographic spread hurdle OR ~0.96, positive geographic
spread count ratio ~0.78. BA.2 and BA.4 show positive positive-count estimates
(BA.2 positive cluster size count ratio 1.19; BA.4 1.68), in contrast to Delta.
Earlier waves (B.1.177, Alpha) show estimates closer to the null with wider
confidence intervals. BA.4 has the widest intervals, reflecting its small sample
(2,669 clusters). B.1.177, Alpha, BA.5, and BQ.1 show more mixed or
outcome-specific patterns.

**Suggested results paragraph.**
> Per-wave SIMD models showed that the pooled negative association with cluster
> size and geographic spread was not stable across the epidemic. During the Delta
> wave, deprivation was consistently associated with lower odds of non-singleton
> clusters (OR 0.934, 95% CI 0.921–0.947), smaller positive cluster size (count
> ratio 0.797, 95% CI 0.725–0.876), lower odds of multi-datazone clusters (OR
> 0.958, 95% CI 0.945–0.971), and lower positive geographic spread (count ratio
> 0.781, 95% CI 0.703–0.867). BA.2 showed a contrasting pattern, with higher
> positive cluster size (count ratio 1.19, 95% CI 1.07–1.32) and weakly higher
> positive geographic spread (count ratio 1.10, 95% CI 0.97–1.24). BA.4 also
> showed positive positive-count associations, but with wide confidence intervals
> reflecting the small wave sample (n = 2,669). B.1.177 and Alpha showed weaker
> or outcome-specific associations. These wave-specific patterns argue against a
> single stable deprivation effect and suggest that epidemic and lineage context
> substantially modifies how deprivation relates to genomic cluster structure.

---

## Supplementary Figure 1: Outcome Distributions

**File:** `supp_fig1_outcome_distributions`

**What it shows.** A two-row distributional summary for the non-singleton cluster
population. The first row shows the distributions of cluster size, duration, and
distinct datazones among the 84,067 non-singleton clusters. The second row shows
the distributions of observed-minus-expected excess mixing for age, sex, and
SIMD deprivation quintile composition. All panels use the same non-singleton
population.

**Key visual patterns.** All three count outcomes (size, duration, datazones)
show long right tails even after restricting to non-singletons, with modal peaks
at the minimum (size 2, duration 0 days, datazones 1) and medians at 3, 4 days,
and 3 respectively. Maximum values (2,792 sequences, 19 days, 2,100 datazones)
are far above the median. In the mixing row, age and sex excess mixing are
centred slightly above zero, while SIMD deprivation mixing is centred slightly
below zero — visually establishing that same-quintile SIMD pairs are
over-represented within clusters before any regression adjustment.

**Suggested results paragraph.**
> Among 84,067 non-singleton clusters, the distributions of cluster size,
> duration, and geographic spread were highly right-skewed. Median cluster size
> was 3 sequences, with 38.8% of non-singleton clusters having size 2; median
> duration was 4 days, with 15.3% still at duration zero; and median distinct
> datazones was 3, with 12.0% confined to a single datazone. Maximum observed
> values reached 2,792 sequences, 19 days, and 2,100 datazones. Observed-minus-
> expected excess mixing was centred slightly above zero for age and sex and
> slightly below zero for SIMD deprivation quintile, indicating that same-quintile
> SIMD pairs are modestly over-represented within clusters relative to the
> lineage-window expectation, motivating the regression mixing analyses.

---

## Supplementary Figure 2: Size-Adjusted Positive Counts

**File:** `supp_fig2_size_adjusted_positive_counts`

**What it shows.** A coefficient plot comparing the primary positive-count model
for geographic spread with a size-adjusted version that additionally includes
log(cluster size) as a covariate. The SIMD deprivation estimate is highlighted
in both versions to show how the association changes after conditioning on
cluster size.

**Key visual patterns.** The SIMD deprivation estimate for positive geographic
spread visibly flips sign between the primary model (count ratio ~0.85, clearly
left of 1.0) and the size-adjusted model (count ratio ~1.03, right of 1.0). This
sign reversal is the central visual message of the figure: the unadjusted
negative geographic-spread result is explained by deprivation's association with
cluster size, not by reduced geographic diffusion at comparable cluster sizes.

**Suggested results paragraph.**
> After additionally adjusting the positive-count geographic spread model for
> cluster size, the overall SIMD deprivation estimate changed direction from the
> primary model (count ratio 0.851, 95% CI 0.792–0.915) to a weakly positive
> association (count ratio 1.027, 95% CI 1.010–1.044, p < 0.001). This sign
> reversal indicates that more deprived clusters are generally smaller, and that
> the unadjusted negative geographic-spread result is driven primarily by
> cluster-size differences. Among clusters of comparable size, there is a
> marginally wider geographic spread with higher deprivation, suggesting that
> deprivation influences cluster scale rather than geographic diffusion per se.

---

## Supplementary Figure 3: Log-Linear Versus Hurdle/ZTNB

**File:** `supp_fig3_loglinear_vs_hurdle_ztnb`

**What it shows.** A side-by-side coefficient comparison for the SIMD deprivation
association with cluster size and geographic spread, contrasting the log-linear
geometric mean ratio model with the hurdle (odds ratio) and ZTNB (count ratio)
components of the two-part main model. All estimates are on a log-ratio scale
and adjusted for the same set of covariates.

**Key visual patterns.** The log-linear estimates for SIMD deprivation are both
close to 1.0 (cluster size geometric mean ratio 0.992; geographic spread
geometric mean ratio 1.001), with narrow intervals straddling the null. In
contrast, the hurdle and positive-count components show clearly differentiated
estimates in both directions. The figure makes visually explicit that the
log-linear model averages over the structural mass at the minimum and thereby
masks the within-component patterns that the two-part model separates.

**Suggested results paragraph.**
> Log-linear models fitted for comparison produced SIMD deprivation estimates
> substantially attenuated relative to the hurdle/ZTNB main models: geometric
> mean ratios of 0.992 (95% CI 0.987–0.997) for cluster size and 1.001 (95% CI
> 0.996–1.006) for geographic spread. These near-null log-linear estimates
> reflect the averaging of hurdle and positive-count components across the
> structural mass at the count minimum, masking the negative positive-count
> associations visible in the two-part models. This comparison supports the
> decision to use hurdle/ZTNB models as the primary count analytical framework.

---

## Supplementary Figure 4: SIMD-Domain Cluster Outcomes

**File:** `supp_fig4_simd_domain_cluster_outcomes`

**What it shows.** A multi-domain extension of Main Figure 1, showing each SIMD
domain's (income, employment, education, health, crime, access, housing) adjusted
effect on all four count model components: cluster size hurdle, positive cluster
size, geographic spread hurdle, and positive geographic spread. The two-part
count structure is the same as in the main analysis.

**Key visual patterns.** Housing and crime deprivation stand out with the
clearest negative count-component associations: housing is associated with lower
odds of non-singleton clusters, lower odds of multi-datazone clusters, and lower
positive geographic spread; crime with smaller positive cluster size and lower
positive geographic spread. Access deprivation behaves distinctly from all other
domains — it is positively associated with the geographic spread hurdle and with
positive geographic spread — making it visually separable from the leftward
cluster of other domains. Overall SIMD and income/employment estimates are
intermediate in magnitude.

**Suggested results paragraph.**
> SIMD-domain count models showed that housing and crime deprivation had the
> strongest negative associations with positive cluster size and positive
> geographic spread, while access deprivation was associated with higher odds of
> multi-datazone clusters and higher positive geographic spread, in the opposite
> direction to most other domains. These domain-specific differences reinforce the
> conclusion that overall SIMD deprivation conflates distinct social mechanisms
> with different implications for genomic cluster structure, and that the access
> domain captures a different geography of deprivation from housing or crime.

---

## Supplementary Figure 5: SIMD-Domain Quintile Mixing

**File:** `supp_fig5_simd_domain_quintile_mixing`

**What it shows.** A single-panel coefficient plot isolating the domain-quintile
mixing estimates for each SIMD domain, analogous to panel A of Main Figure 3
but without the demographic mixing panels. The SIMD domains (income, employment,
education, health, crime, access, housing) are displayed alongside the overall
SIMD estimate for reference.

**Key visual patterns.** The education/crime versus access/housing contrast is
the dominant visual feature: education (+1.17 pp) and crime (+1.09 pp) are
clearly right of zero, access (−2.08 pp) and housing (−1.19 pp) clearly left,
spanning the widest range of any figure in the mixing series. Overall SIMD sits
near zero between these extremes. Income, employment, and health estimates are
closer to the null.

**Suggested results paragraph.**
> The domain-quintile mixing figure (Supplementary Figure 5) isolates the
> socioeconomic mixing contrast found in Main Figure 3. Education and crime
> deprivation were associated with greater within-cluster mixing across domain
> quintiles (education +1.17 pp, 95% CI +0.64 to +1.69; crime +1.09 pp, 95% CI
> +0.55 to +1.64), while access deprivation (−2.08 pp, 95% CI −2.50 to −1.66)
> and housing deprivation (−1.19 pp, 95% CI −1.77 to −0.62) showed the opposite
> pattern. This contrast supports interpreting domain-specific deprivation as
> capturing qualitatively different social geographies rather than a single
> deprivation mechanism.

---

## Supplementary Figure 6: SIMD-Domain Demographic Mixing

**File:** `supp_fig6_simd_domain_demographic_mixing`

**What it shows.** A three-panel coefficient plot showing age excess mixing (A),
sex excess mixing (B), and joint age-sex profile excess mixing (C) for each SIMD
domain and for overall SIMD. Each estimate is the adjusted change in excess
discordance per 1 SD domain-deprivation increase among non-singleton clusters,
with 95% confidence intervals.

**Key visual patterns.** Age mixing (panel A) is positive for most domains,
forming a cluster of estimates above zero; access deprivation is a clear visual
outlier below zero. Sex mixing (panel B) is predominantly negative across
domains, mirroring the main SIMD result. Joint age-sex profile mixing (panel C)
is broadly positive but more variable in magnitude than age mixing alone. Access
deprivation consistently behaves differently from other domains across all three
panels.

**Suggested results paragraph.**
> Across SIMD domains, age excess mixing was positively associated with
> deprivation for most domains (income, employment, education, health, crime,
> housing), with access deprivation as the exception showing lower age mixing.
> Sex mixing showed the reverse pattern, with most domains showing negative
> associations. Joint age-sex profile mixing was broadly positive for most
> domains but more variable across them. The consistent behavioural exception of
> the access domain across all three demographic mixing outcomes is notable and
> may reflect the distinct geographic and social processes captured by the
> access-to-services dimension of deprivation.

---

## Supplementary Figure 7: Wave-Specific Domain-Demographic Mixing

**File:** `supp_fig7_wave_specific_domain_demographic_mixing`

**What it shows.** A heatmap showing domain-specific adjusted changes in age
mixing, sex mixing, and joint age-sex profile mixing for each SIMD domain across
epidemic waves. Each cell represents the adjusted excess-mixing estimate for one
domain × wave combination, coloured by direction and magnitude. This is an
exploratory heterogeneity analysis.

**Key visual patterns.** Age mixing shows the most consistent positive domain-
deprivation pattern across waves, with positive values visible in most domain ×
wave cells for B.1.177, Alpha, and Delta. Access deprivation is predominantly
negative for age mixing across all waves, making it a consistent exception.
Sex mixing is more variable: BA.5 shows several negative estimates across
domains, whereas access in BA.5 is positive. Joint age-sex mixing is generally
positive in earlier waves and Delta for most domains, but varies more in Omicron
subwaves. Cells for smaller waves (BA.4, BQ.1, XBB) are visibly noisier.

**Suggested results paragraph.**
> Wave-specific domain-demographic mixing patterns were broadly consistent with
> the pooled estimates but showed heterogeneity across waves. Positive age mixing
> associations were visible in most domain × wave cells, particularly for
> B.1.177, Alpha, and Delta, with access deprivation as a consistent exception
> across all waves. Sex mixing was more variable: BA.5 showed negative estimates
> for several domains. The results for smaller waves (BA.4, BQ.1, XBB) should be
> interpreted cautiously given reduced cluster counts, and the wave-specific
> heatmap is best read as an exploratory check on pooled-estimate stability
> rather than a definitive wave-stratified analysis.

---

## Supplementary Figure 8: Observed-Versus-Expected Mixing Matrices

**File:** `supp_fig8_observed_expected_mixing_matrices`

**What it shows.** Two heatmaps of observed-minus-expected pairwise probability
matrices: one for SIMD quintile pairs (left) and one for age-band pairs (right).
Each cell shows the mean excess probability for a given pair type across all
non-singleton clusters, expressed in percentage points above or below expectation
within the same lineage and calendar window. The matrices give a descriptive view
of which pairings drive the aggregate discordance metrics.

**Key visual patterns.** For SIMD quintiles, the diagonal (same-quintile pairs)
is positive throughout, with the strongest excess at quintile 1 × quintile 1
(most deprived, +0.28 pp). Off-diagonal cells are close to zero or slightly
negative, particularly for quintile 4 × 5 and quintile 1 × 5 pairings. For age,
the within-band diagonal is positive across bands, peaking at 20–24 × 20–24
(+0.21 pp); adjacent young-adult pairs (20–24 × 25–29) are also positive. Some
child-to-adult pairings (0–4, 5–12 with adult bands) fall slightly below zero.

**Suggested results paragraph.**
> The observed-minus-expected mixing matrices (Supplementary Figure 8) illustrate
> the pairwise structure underlying the aggregate excess-discordance metrics. For
> SIMD quintiles, same-quintile pairs showed positive excess in all quintiles,
> with the most-deprived quintile 1 pairs showing the largest excess (+0.28 pp),
> consistent with spatially concentrated transmission within deprived areas. For
> age, within-band pairs showed the strongest positive excess in young adults
> (20–24: +0.21 pp), with adjacent young-adult pairs also above expectation.
> Some child-to-adult pairings were slightly below expectation. These patterns
> are consistent with genomic clusters preferentially capturing transmission
> chains within socioeconomically and age-similar groups, even after accounting
> for lineage and calendar-window composition.

---

## Sensitivity Analysis Results

Five sensitivity runs were generated with `part1_sensitivities.sh`. Each run
has its own result tables and manuscript figures so that the primary outputs are
not overwritten:

| Sensitivity | Tables | Figures |
|---|---|---|
| Health-board clustered standard errors | `part1/main/tables_health_board` | `part1/main/manuscript/figures_health_board` |
| Cluster-size positive-count offset | `part1/main/tables_size_offset` | `part1/main/manuscript/figures_size_offset` |
| Index-case SIMD exposure | `part1/main/tables_index_simd` | `part1/main/manuscript/figures_index_simd` |
| 99th-percentile positive-count winsorisation | `part1/main/tables_winsorise99` | `part1/main/manuscript/figures_winsorise99` |
| Approximately non-overlapping windows | `part1/main/tables_stride3` | `part1/main/manuscript/figures_stride3` |

### Health-Board Clustered Standard Errors

This sensitivity keeps the fitted coefficients unchanged but clusters standard
errors by health board rather than by sliding window. It is therefore a
standard-error sensitivity, not a point-estimate sensitivity.

The main count point estimates are unchanged, but uncertainty is much wider for
the primary unadjusted count components:

- Cluster size hurdle: OR 0.971, 95% CI 0.914 to 1.032.
- Positive cluster size: count ratio 0.926, 95% CI 0.664 to 1.293.
- Positive geographic spread: count ratio 0.851, 95% CI 0.577 to 1.255.

The size-adjusted positive-count result remains positive:

- Positive geographic spread, size-adjusted: count ratio 1.027, 95% CI 1.012
  to 1.041.

For mixing, age and sex associations remain clearly different from zero under
health-board clustering, while the joint-profile estimate weakens:

- SIMD mixing: +0.31 pp, 95% CI -3.36 to +3.98.
- Age mixing: +1.66 pp, 95% CI +0.98 to +2.34.
- Sex mixing: -0.78 pp, 95% CI -1.28 to -0.27.
- Joint profile mixing: +0.48 pp, 95% CI -0.04 to +0.99.

Interpretation: conclusions based on exact statistical significance are
sensitive to clustering level, but the direction of the main deprivation
estimates is unchanged. The most robust mixing signals are age and sex.

### Cluster-Size Positive-Count Offset

This sensitivity adds `log(wn_no_sequences)` as an offset in the cluster-size
positive-count model, changing the estimand from raw positive cluster size to
positive cluster size relative to the analysis-window sequencing pool.

Only the cluster-size positive model changes. The SIMD deprivation estimate is
almost identical to the primary result:

- Primary positive cluster size: count ratio 0.926, 95% CI 0.869 to 0.987.
- Offset positive cluster size: count ratio 0.925, 95% CI 0.868 to 0.985.

Interpretation: the negative deprivation association for positive cluster size
is not explained by differences in the number of sequences available in the
analysis window.

### Index-Case SIMD Exposure

This sensitivity replaces mean cluster SIMD deprivation with the SIMD
deprivation of the earliest collected sequence in the cluster. The corrected
index-SIMD sensitivity output applies this exposure consistently to the count
models, mixing models, and sensitivity figures.

For the count models, index-case SIMD gives a different pattern from mean
cluster SIMD:

- Cluster size hurdle: OR 0.967, 95% CI 0.956 to 0.978.
- Positive cluster size: count ratio 0.996, 95% CI 0.957 to 1.036.
- Geographic spread hurdle: OR 0.994, 95% CI 0.983 to 1.006.
- Positive geographic spread: count ratio 0.989, 95% CI 0.951 to 1.028.
- Positive geographic spread, size-adjusted: count ratio 1.002, 95% CI 0.998
  to 1.006.

The index-case SIMD mixing estimates are also attenuated relative to the main
mean-cluster SIMD estimates:

- SIMD mixing: +0.13 pp, 95% CI -0.16 to +0.42.
- Age mixing: +0.83 pp, 95% CI +0.61 to +1.06.
- Sex mixing: -0.58 pp, 95% CI -0.84 to -0.31.
- Joint profile mixing: +0.26 pp, 95% CI +0.15 to +0.37.

Interpretation: the negative positive-count associations and demographic
mixing associations for mean cluster SIMD are much weaker when the exposure is
the index-case SIMD. This suggests that mean cluster deprivation is capturing
the composition of the whole cluster, not only the deprivation context of the
earliest observed case.

### Winsorising Positive Counts At The 99th Percentile

This sensitivity caps positive count outcomes at the 99th percentile before
fitting the ZTNB models. Binary hurdle components are unchanged.

The main effect is attenuation of the positive cluster-size association:

- Primary positive cluster size: count ratio 0.926, 95% CI 0.869 to 0.987.
- Winsorised positive cluster size: count ratio 0.952, 95% CI 0.900 to 1.006.

The negative positive geographic-spread association remains:

- Primary positive geographic spread: count ratio 0.851, 95% CI 0.792 to 0.915.
- Winsorised positive geographic spread: count ratio 0.889, 95% CI 0.837 to
  0.944.

The size-adjusted positive-count result remains positive:

- Positive geographic spread, size-adjusted: count ratio 1.027, 95% CI 1.019
  to 1.036.

Interpretation: the positive cluster-size deprivation estimate is sensitive to
the extreme right tail, while the negative positive geographic-spread result and
the size-adjusted positive results are less sensitive to winsorisation.

### Approximately Non-Overlapping Windows

This sensitivity keeps only windows where `window_idx % 3 == 0`, reducing the
cluster table from 193,112 clusters to 63,991 clusters and the mixing-model
population from 84,067 to 27,897 non-singleton clusters.

Count-model results are directionally similar but less precise:

- Cluster size hurdle: OR 0.972, 95% CI 0.954 to 0.991.
- Positive cluster size: count ratio 0.920, 95% CI 0.811 to 1.044.
- Geographic spread hurdle: OR 1.003, 95% CI 0.983 to 1.024.
- Positive geographic spread: count ratio 0.861, 95% CI 0.745 to 0.996.
- Positive geographic spread, size-adjusted: count ratio 1.027, 95% CI 1.013
  to 1.041.

Mixing-model results remain close to the primary estimates for age, sex, and
joint profile:

- SIMD mixing: +0.24 pp, 95% CI -0.59 to +1.07.
- Age mixing: +1.63 pp, 95% CI +0.92 to +2.34.
- Sex mixing: -1.02 pp, 95% CI -1.63 to -0.40.
- Joint profile mixing: +0.44 pp, 95% CI +0.16 to +0.72.

Interpretation: reducing repeated-window dependence weakens precision, but the
main qualitative story is stable: no clear SIMD-mixing effect, positive age and
joint-profile mixing effects, negative sex-mixing effect, and only small or
negative deprivation associations for unadjusted count magnitude.

### Overall Sensitivity Interpretation

The sensitivity analyses support a cautious version of the main findings:

- The direction of the main mean-SIMD count estimates is generally stable, but
  their precision depends on the standard-error clustering level and on how the
  positive-count tail is handled.
- The strongest unadjusted count result that persists across several
  sensitivities is lower positive geographic spread with higher mean cluster
  deprivation.
- Size-adjusted positive geographic spread remains weakly positive under
  health-board clustering, winsorisation, and approximately non-overlapping
  windows.
- Index-case SIMD does not reproduce the positive-count associations seen with
  mean cluster SIMD, implying that the cluster-level composition exposure is
  scientifically different from the earliest-observed-case exposure.
- Mixing results are comparatively stable for age and sex. SIMD-quintile mixing
  remains near null, and joint-profile mixing is directionally positive but
  less robust to health-board clustering.

## Overall Suggested Results Paragraph

In the main hurdle/ZTNB models, higher cluster-level SIMD deprivation was not
associated with larger or more geographically dispersed genomic clusters after
adjustment for lineage, calendar time, local incidence, sequencing intensity,
and test positivity. Instead, deprivation was associated with slightly lower
odds of being non-singleton and smaller positive cluster size, and with
substantially lower positive geographic spread. Local epidemic and surveillance
conditions were much more strongly associated with cluster scale: higher
incidence, window-level sequencing proportion, and test positivity were
consistently associated with larger and more geographically dispersed clusters.
Mixing analyses showed a more nuanced socioeconomic pattern. Overall
deprivation was not clearly associated with SIMD-quintile mixing, but was
associated with greater age and joint socio-demographic mixing and lower sex
mixing. SIMD-domain models showed that education and crime deprivation were
associated with greater domain-quintile mixing, whereas access and housing
deprivation were associated with lower domain-quintile mixing. Per-wave outcome
models indicated that deprivation effects varied over time, with the strongest
negative associations seen during Delta and more heterogeneous patterns in
Omicron subwaves. Sensitivity analyses supported the qualitative pattern but
showed that statistical precision depends on the clustering level and on the
positive-count tail: health-board clustered standard errors widened several
count-outcome intervals, 99th-percentile winsorisation attenuated the positive
cluster-size association, and the index-case SIMD exposure did not reproduce
the mean-cluster-SIMD positive-count associations.

## Interpretation For The Part 1 Question

The answer to the original question is therefore mixed:

- Socioeconomic deprivation is associated with some cluster outcomes, but not in
  the simple direction of larger or more dispersed clusters.
- Surveillance and epidemic-intensity variables are strongly associated with
  apparent cluster scale, reinforcing the need to adjust for sequencing and
  testing context.
- SIMD domain results show that deprivation is not a single mechanism. Access,
  housing, education, crime, income, employment, and health domains capture
  different kinds of local social structure.
- Demographic mixing is where deprivation signals are clearest, especially for
  age mixing and joint profile mixing.
- Wave-specific analyses suggest that lineage and epidemic context modify these
  associations.
- Sensitivity analyses reinforce the need to separate point-estimate stability
  from statistical precision: directions are mostly stable, but inference is
  weaker under health-board clustered standard errors and when the heaviest
  positive-count tail is capped.

These results support a cautious conclusion: deprivation and local surveillance
conditions are associated with observed genomic cluster structure, but the
associations are outcome-specific, domain-specific, and wave-specific rather
than a single monotonic deprivation effect.
