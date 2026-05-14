# Part 1 Main Results And Figure Guide

Updated for the main Part 1 analysis on 13 May 2026.

## Analysis Question

Part 1 asks whether socioeconomic deprivation and local surveillance conditions
are associated with SARS-CoV-2 genomic cluster characteristics in Scotland after
accounting for lineage, calendar time, local incidence, sequencing intensity,
and test positivity, using one primary Leiden clustering resolution.

The specific outcomes are:

- cluster size
- geographic spread, measured by the number of distinct datazones in a cluster
- observed-minus-expected socioeconomic and demographic mixing within clusters

The main manuscript figures are generated from `part1/tables` and
`part1/cache`, with exported figure files in `part1/manuscript/figures`.

## Analysis Population And Modelling Frame

The main analysis uses sequences passing the overall Nextclade quality filter at
Leiden resolution $\gamma = 0.3$. This produces
789,347 sequence rows and 193,112 inferred genomic clusters across 134 time
windows. Rare Pango lineages are pooled at a threshold of 50 clusters, giving
183 lineage model levels in the main pooled analysis.

The modelled count outcomes have large structural masses at their minimum
values, so the main count analysis uses two-part hurdle models:

- a binomial hurdle component for whether a cluster exceeds the structural
  minimum
- a zero-truncated negative-binomial positive-count component for the outcome
  among clusters exceeding the structural minimum

For cluster size, the hurdle distinguishes singleton from non-singleton
clusters and the positive count is the number of additional sequences beyond a
singleton. For geographic spread, the hurdle distinguishes single-datazone from
multi-datazone clusters and the positive count is the number of additional
datazones beyond one.

Cluster duration remains in the descriptive cluster summaries and the
supplementary outcome-distribution figure, but it is not modelled in the
current main count analysis because the fixed three-week clustering windows
mechanically constrain the observed span.

Mixing outcomes are fitted among non-singleton clusters. They are defined as
observed pairwise discordance within a cluster minus expected discordance among
sampled cases from the same lineage and calendar window. Positive values mean
more mixing than expected; negative values mean more homogeneity than expected.

Throughout the figures, higher SIMD deprivation means a more deprived cluster
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

## Main Figure 1: Deprivation As Exposure — Overall Effects On Cluster Outcomes And Mixing

**File:** `fig1_deprivation_overall`

**What it shows.** A combined eight-panel figure of the pooled deprivation-as-
exposure models. The top row shows the four count-model components for cluster
size and geographic spread — cluster-size hurdle (odds ratio), positive cluster
size (ZTNB count ratio), geographic-spread hurdle (OR), and positive geographic
spread (ZTNB count ratio) — with adjusted ratios per 1 SD higher covariate and
95% CIs. The bottom row shows adjusted percentage-point differences in
observed-minus-expected pairwise discordance for the four mixing outcomes
(SIMD-quintile, age-band, sex, joint SIMD-age-sex profile) among non-singleton
clusters. Each panel displays SIMD deprivation alongside the four surveillance
covariates; the mixing panels additionally include log cluster size.

**Key visual patterns.** In the top row, the SIMD deprivation point falls
clearly below the null for positive cluster size (count ratio 0.926) and
positive geographic spread (count ratio 0.851), modestly below the null for the
cluster-size hurdle (OR 0.971), and is close to the null for the
geographic-spread hurdle.
Surveillance covariates dominate the panel visually: test positivity gives
positive-count ratios of 2.649 (cluster size) and 2.999 (geographic spread);
local sequencing fraction 3.240 and 2.269. In the bottom row, SIMD deprivation produces
its largest mixing estimate for age mixing (+1.66 pp) and a clearly negative
estimate for sex mixing (−0.78 pp); SIMD-quintile mixing straddles zero (+0.31
pp), and joint-profile mixing is positive but smaller (+0.48 pp). Log cluster
size is a strong cross-outcome driver of all four mixing components.

**Suggested results paragraph.**
> In the pooled hurdle/ZTNB count models (top row), higher cluster-level SIMD
> deprivation was associated with slightly lower odds of being non-singleton
> (OR 0.971, 95% CI 0.960–0.983, p < 0.001), smaller positive cluster size
> (count ratio 0.926, 95% CI 0.869–0.987, p = 0.018), and substantially lower
> positive geographic spread (count ratio 0.851, 95% CI 0.792–0.915, p < 0.001).
> Local sequencing fraction (count ratio 3.240 for positive cluster size; 2.269
> for positive geographic spread), test positivity (2.649 and 2.999), and
> window-level sequencing proportion showed far larger positive associations
> with cluster scale, consistent with surveillance and epidemic intensity being
> the dominant structural determinants of observed cluster size and geographic
> spread. In the linear mixing models (bottom row), deprivation was not clearly
> associated with SIMD-quintile mixing (+0.31 pp, 95% CI −0.18 to +0.80), but
> was positively associated with age mixing (+1.66 pp, 95% CI +1.29 to +2.03),
> negatively associated with sex mixing (−0.78 pp, 95% CI −1.16 to −0.39), and
> weakly positively associated with joint sociodemographic profile mixing
> (+0.48 pp, 95% CI +0.29 to +0.67).

---

## Main Figure 2: Deprivation As Exposure — Wave-Specific Effects On Cluster Outcomes

**File:** `fig2_deprivation_wave_specific`

**What it shows.** A multi-wave coefficient plot showing the per-wave overall
SIMD deprivation effect on the four count-model components: cluster-size hurdle
(OR), positive cluster size (ZTNB count ratio), geographic-spread hurdle (OR),
and positive geographic spread (ZTNB count ratio). Wave-specific models are
fitted separately for the eight dominant Pango-lineage waves (B.1.177, Alpha,
Delta, BA.1, BA.2, BA.4, BA.5, BQ.1) and adjust for the same surveillance and
calendar-time covariates as the pooled main model; within-wave lineage dummies
are included where estimable.

**Key visual patterns.** Delta stands out with all four estimates clearly below
their null values: cluster-size hurdle OR 0.934, positive cluster size count
ratio 0.797, geographic-spread hurdle OR 0.958, positive geographic-spread
count ratio 0.781. BA.2 and BA.4 show positive positive-count estimates (BA.2
positive cluster size count ratio 1.176; BA.4 1.674), in contrast to Delta.
Earlier waves (B.1.177, Alpha) show estimates closer to the null with wider
confidence intervals. BA.4 has the widest intervals, reflecting its small
sample (2,669 clusters).

**Suggested results paragraph.**
> Per-wave SIMD models showed that the pooled negative association with cluster
> size and geographic spread was not stable across the epidemic. During the
> Delta wave, deprivation was consistently associated with lower odds of
> non-singleton clusters (OR 0.934, 95% CI 0.921–0.947), smaller positive
> cluster size (count ratio 0.797, 95% CI 0.725–0.876), lower odds of
> multi-datazone clusters (OR 0.958, 95% CI 0.945–0.971), and lower positive
> geographic spread (count ratio 0.781, 95% CI 0.703–0.867). BA.2 showed a
> contrasting pattern, with higher positive cluster size (count ratio 1.176,
> 95% CI 1.061–1.303) and weakly higher positive geographic spread (count ratio
> 1.106, 95% CI 0.978–1.252). BA.4 also showed positive positive-count
> associations, but with wide confidence intervals reflecting the small wave
> sample (n = 2,669). B.1.177 and Alpha showed weaker or outcome-specific
> associations. These wave-specific patterns argue against a single stable
> deprivation effect and suggest that epidemic and lineage context
> substantially modifies how deprivation relates to genomic cluster structure.

---

## Main Figure 3: Excess Mixing As Predictor — Overall Effects On Cluster Scale

**File:** `fig3_mixing_overall`

**What it shows.** A three-panel coefficient plot of the four cluster-level
excess-mixing scores (SIMD-quintile, age, sex, joint SIMD-age-sex profile)
entered jointly as predictors of cluster scale alongside SIMD deprivation, the
four surveillance covariates, lineage fixed effects, and the calendar B-spline.
Panels: (A) positive cluster size (ZTNB count ratio; n = 84,067), (B)
geographic-spread hurdle (OR; multi- vs single-datazone among non-singletons;
n = 84,067), (C) positive geographic spread (ZTNB count ratio; n = 74,010). The
cluster-size hurdle is omitted because mixing scores are undefined for
singletons (see Methods). Each point is an adjusted ratio per 1 SD higher
mixing score with 95% CIs.

**Key visual patterns.** SIMD-quintile excess discordance is the dominant
positive predictor in all three panels: positive cluster size CR 3.478,
positive geographic spread CR 3.029, and geographic-spread hurdle OR 22.107. Age
excess discordance is also strongly positive on the ZTNB components (CR 1.668
positive cluster size; 1.966 positive geographic spread) and weakly positive on
the hurdle (OR 1.283). Sex excess discordance is small but negative on positive
cluster size and on the hurdle; joint SIMD-age-sex profile excess discordance
is negative on both ZTNB components (CR 0.812 positive cluster size; 0.718
positive geographic spread).

**Suggested results paragraph.**
> Among non-singleton clusters, SIMD-quintile excess discordance was the
> strongest mixing-side predictor of cluster scale: per 1 SD higher score,
> positive cluster size was 3.478-fold higher (95% CI 3.245–3.727), positive
> geographic spread 3.029-fold higher (95% CI 2.793–3.285), and the odds of
> multi-datazone spread 22.107-fold higher (95% CI 18.977–25.753). Age excess
> discordance was also substantially positively associated with both ZTNB
> components (positive cluster size CR 1.668, 95% CI 1.560–1.783; positive
> geographic spread CR 1.966, 95% CI 1.800–2.148) and weakly positively associated
> with the geographic-spread hurdle (OR 1.283, 95% CI 1.246–1.322). Joint
> SIMD-age-sex profile excess discordance was negatively associated with both
> ZTNB components (cluster size CR 0.812, 95% CI 0.763–0.865; geographic spread
> CR 0.718, 95% CI 0.665–0.775). Clusters bridging more across SIMD quintiles and
> age bands than the lineage-window baseline expectation are therefore detected
> as substantially larger and more geographically dispersed, while clusters
> whose excess discordance is concentrated in fine-grained joint profile
> coordinates are detected as smaller and less dispersed.

---

## Main Figure 4: Excess Mixing As Predictor — Wave-Specific Effects On Cluster Scale (ZTNB)

**File:** `fig4_mixing_wave_specific`

**What it shows.** A two-panel heatmap of per-wave adjusted ZTNB count ratios
per 1 SD higher excess-mixing score for the four mixing predictors
(SIMD-quintile, age, sex, joint SIMD-age-sex profile excess discordance),
entered jointly with SIMD deprivation, the four surveillance covariates,
lineage fixed effects, and the calendar B-spline. Rows are the eight dominant
Pango-lineage waves; columns are the four mixing predictors. Panels: (A)
positive cluster size (n = 84,067) and (B) positive geographic spread
(n = 74,010). Cells are coloured on a shared symmetric ratio scale (centred at
1, capped at ratio 5) and annotated with the count ratio. The corresponding
geographic-spread hurdle component is reported in Supplementary Table 1 rather
than as a third heatmap panel, because the multiple-datazone outcome is
positive for 88.0% of non-singleton clusters and combines with the strong
SIMD-quintile excess-mixing predictor to drive the binomial logistic component
toward quasi-separation in some waves, producing adjusted odds ratios
(approximately 29,000 in Alpha) that are uninformative on a heatmap.

**Key visual patterns.** Positive SIMD-quintile and age-band mixing-predictor
effects on the ZTNB components are directionally stable across waves but
largest in Alpha, Delta, BA.1, and BA.2 (SIMD column ratios 3.2–4.4 in panel A;
2.8–4.1 in panel B). Sex excess discordance is consistently negative or near
null, with the most extreme cell at Alpha geographic-spread ZTNB ratio 0.075.
Joint-profile excess discordance is mostly near 1 in panel A and ranges from
0.54 to 1.6 in panel B. Late Omicron subwaves (BA.4, BA.5, BQ.1) show
attenuated and noisier estimates.

**Suggested results paragraph.**
> Wave-stratified mixing-predictor count models (Figure 4 ZTNB components;
> Supplementary Table 1 for the companion geographic-spread hurdle) showed
> that the pooled SIMD-quintile and age-mixing associations with larger and
> more geographically dispersed clusters were directionally stable across
> waves but varied in magnitude. The largest positive SIMD-quintile mixing
> effects on positive cluster size (count ratios 3.4–4.4) and positive
> geographic spread (count ratios 2.8–4.1) were observed in waves with the
> deepest cluster recruitment (Alpha, Delta, BA.1, BA.2), mirroring the
> per-wave deprivation-as-exposure pattern in Figure 2. Late Omicron subwaves
> (BA.4, BA.5, BQ.1) produced smaller and noisier mixing-predictor effects.
> The geographic-spread hurdle component is reported as Supplementary
> Table 1 because the heavily imbalanced multiple-datazone outcome, positive
> in 88.0% of non-singleton clusters, combines with the strong SIMD-quintile
> excess-mixing predictor to produce implausibly large adjusted odds ratios
> (e.g. approximately 29,000 in Alpha) that should be read as evidence of strong direction rather than as
> interpretable effect sizes.

---

## Supplementary Figure 1: Outcome Distributions

**File:** `supp_fig1_outcome_distributions`

**What it shows.** A two-row distributional summary for the 84,067 non-singleton
clusters. Top row: histograms of cluster size, cluster duration, and number of
distinct datazones. Bottom row: histograms of observed-minus-expected excess
mixing for age, sex, and SIMD-deprivation quintile composition. All panels use
the same non-singleton population.

**Key visual patterns.** All three count outcomes show long right tails even
after restricting to non-singletons, with modal peaks at the structural minima
(size 2, duration 0 days, datazones 1) and medians at 3, 4 days, and 3
respectively. Maximum observed values reach 2,792 sequences, 19 days, and 2,100
datazones. Age and sex excess mixing centre slightly above zero; SIMD
deprivation excess mixing centres slightly below zero, visually establishing
that same-quintile SIMD pairs are over-represented within clusters before any
regression adjustment.

**Suggested results paragraph.**
> Among 84,067 non-singleton clusters, the distributions of cluster size,
> duration, and geographic spread were highly right-skewed. Median cluster size
> was 3 sequences, with 38.8% of non-singleton clusters having size 2; median
> duration was 4 days, with 15.3% still at duration zero; and median distinct
> datazones was 3, with 12.0% confined to a single datazone. Maximum observed
> values reached 2,792 sequences, 19 days, and 2,100 datazones. Observed-minus-
> expected excess mixing was centred slightly above zero for age and sex and
> slightly below zero for SIMD deprivation quintile, indicating that
> same-quintile SIMD pairs are modestly over-represented within clusters
> relative to the lineage-window expectation, motivating the regression mixing
> analyses.

---

## Supplementary Figure 2: Excess Mixing Predictor Distributions

**File:** `supp_fig2_mixing_distributions`

**What it shows.** Histograms of the four cluster-level mixing predictors used
as explanatory variables in the mixing-predictor count models (Main Figures 3
and 4): SIMD-quintile, age, sex, and joint age-sex profile excess-mixing
scores. Each score is the observed-minus-expected pair-discordance value for
that attribute within non-singleton clusters, on the percentage-point scale.

**Key visual patterns.** SIMD-quintile and age scores are unimodal, centred
near zero with longer positive than negative tails, indicating that clusters
more often bridge across quintiles or age bands than the lineage-window
expectation predicts. Sex excess discordance is more symmetric and narrower.
Joint SIMD-age-sex profile excess discordance has the heaviest right tail,
reflecting that fine-grained recombinations are over-represented in a small
number of large clusters.

**Suggested results paragraph.**
> The distributions of the four cluster-level excess-mixing scores
> (Supplementary Figure 2) show unimodal centred distributions with positive
> skew for SIMD-quintile, age, and joint sociodemographic profile scores, and
> a more symmetric distribution for sex. These distributions motivate the
> linear entry of each mixing score (on its native percentage-point scale) as
> a predictor in the mixing-predictor count models.

---

## Supplementary Figure 3: Observed-Versus-Expected Pair-Probability Matrices

**File:** `supp_fig3_observed_expected_matrices`

**What it shows.** Two heatmaps of observed-minus-expected pairwise probability
matrices: one for SIMD-quintile pairs (left) and one for age-band pairs
(right). Cells show the mean excess pair probability in percentage points
relative to the lineage- and calendar-window-matched expectation, with each
cell annotated with the numeric value.

**Key visual patterns.** For SIMD quintiles, the diagonal (same-quintile
pairs) is positive throughout, with the strongest excess at pairs of two
most-deprived-quintile cases (+0.28 pp). Off-diagonal cells are close to zero or
slightly negative. For age, the within-band diagonal is positive across bands,
peaking at pairs of two 20-24-year-old cases (+0.21 pp); adjacent young-adult pairs (20–24 ×
25–29) are also positive. Some child-to-adult pairings fall slightly below
zero.

**Suggested results paragraph.**
> The observed-minus-expected mixing matrices (Supplementary Figure 3)
> illustrate the pairwise structure underlying the aggregate excess-discordance
> metrics. Same-quintile SIMD pairs showed positive excess in all quintiles,
> with the most-deprived quintile pairs showing the largest excess
> (+0.28 pp), consistent with spatially concentrated transmission within
> deprived areas. Within-age-band pairs showed the strongest positive excess
> in young adults (20-24 years: +0.21 pp), with adjacent young-adult pairs also above
> expectation. These patterns are consistent with genomic clusters
> preferentially capturing transmission chains within socioeconomically and
> age-similar groups, even after accounting for lineage and calendar-window
> composition.

---

## Supplementary Figure 4: SIMD-Domain Deprivation Effects On Cluster Outcomes

**File:** `supp_fig4_deprivation_domain_outcomes`

**What it shows.** A multi-domain extension of Main Figure 1 (top row), showing
each SIMD subdomain's (income, employment, education, health, access, crime,
housing) adjusted effect on all four count-model components: cluster-size
hurdle, positive cluster size, geographic-spread hurdle, and positive
geographic spread. The two-part count structure is the same as in the main
analysis.

**Key visual patterns.** Housing and crime deprivation stand out with the
clearest negative count-component associations: housing is associated with
lower odds of non-singleton clusters, lower odds of multi-datazone clusters,
and lower positive geographic spread; crime with smaller positive cluster size
and lower positive geographic spread. Access deprivation behaves distinctly
from all other domains — it is positively associated with the geographic-spread
hurdle and with positive geographic spread — making it visually separable from
the leftward cluster of other domains. Overall SIMD and income/employment
estimates are intermediate in magnitude.

**Suggested results paragraph.**
> SIMD-domain count models showed that housing and crime deprivation had the
> strongest negative associations with positive cluster size and positive
> geographic spread, while access deprivation was associated with higher odds
> of multi-datazone clusters and higher positive geographic spread, in the
> opposite direction to most other domains. These domain-specific differences
> reinforce the conclusion that overall SIMD deprivation conflates distinct
> social mechanisms with different implications for genomic cluster structure.

---

## Supplementary Figure 5: SIMD-Domain Deprivation Effects On Mixing Outcomes

**File:** `supp_fig5_deprivation_domain_mixing`

**What it shows.** A four-panel coefficient plot showing the per-1-SD effect of
each SIMD subdomain on (A) domain-quintile excess mixing, (B) age excess
mixing, (C) sex excess mixing, and (D) joint age-sex profile excess mixing.
Each estimate is the adjusted percentage-point change in excess discordance
among non-singleton clusters, with 95% CIs.

**Key visual patterns.** In panel A, two domains (access, −2.08 pp; housing,
−1.19 pp) are clearly left of zero, while education (+1.17 pp) and crime
(+1.09 pp) are clearly right of zero — a contrast spanning roughly 3.0 pp. The
overall SIMD estimate sits between them, close to zero. In panel B (age
mixing), most domain estimates including overall SIMD are positive, with
access deprivation as an outlier in the negative direction. Panel C (sex
mixing) shows predominantly negative estimates across domains. Panel D (joint
age-sex mixing) is mostly positive with more variable magnitudes.

**Suggested results paragraph.**
> SIMD-domain mixing models revealed heterogeneity that the overall SIMD
> measure conceals. Education deprivation (+1.17 pp, 95% CI +0.64 to +1.69)
> and crime deprivation (+1.09 pp, 95% CI +0.55 to +1.64) were associated with
> greater within-cluster mixing across domain quintiles, whereas access
> deprivation (−2.08 pp, 95% CI −2.50 to −1.66) and housing deprivation
> (−1.19 pp, 95% CI −1.77 to −0.62) were associated with less mixing — opposite
> in direction to education and crime. Demographic mixing patterns were
> broadly consistent with the main SIMD estimates: most domains showed
> positive age mixing and negative sex mixing, with access deprivation
> behaving as an outlier across all three demographic mixing outcomes.

---

## Supplementary Figure 6: Wave-Specific SIMD-Domain Deprivation Effects On Demographic Mixing

**File:** `supp_fig6_deprivation_domain_wave_mixing`

**What it shows.** Three stacked heatmaps showing the per-1-SD effect of each
SIMD subdomain on (A) age, (B) sex, and (C) joint age-sex excess mixing for
each epidemic wave. Rows are SIMD subdomains (Overall, Income, Employment,
Education, Health, Access, Crime, Housing); columns are wave groups (B.1.177,
Alpha, Delta, BA.1, BA.2, BA.5). Cells are coloured by adjusted percentage-
point change on a shared symmetric scale (±5 pp) and annotated with the
numeric value.

**Key visual patterns.** Age mixing shows the most consistent positive
deprivation pattern across waves, with positive values in most domain × wave
cells for B.1.177, Alpha, and Delta. Access deprivation is predominantly
negative for age mixing across all waves, making it a consistent exception.
Sex mixing turns negative in BA.2 and BA.5 for most domains. Joint age-sex
mixing is generally positive in earlier waves and Delta for most domains, but
varies more in Omicron subwaves.

**Suggested results paragraph.**
> Wave-specific domain-demographic mixing patterns were broadly consistent
> with the pooled estimates but showed heterogeneity across waves. Positive
> age mixing associations were visible in most domain × wave cells,
> particularly for B.1.177, Alpha, and Delta, with access deprivation as a
> consistent exception across all waves. Sex mixing was more variable: BA.5
> showed negative estimates for several domains. The wave-specific heatmap is
> best read as an exploratory check on pooled-estimate stability rather than
> a definitive wave-stratified analysis.

---

## Supplementary Figure 7: Domain-Specific Mixing-Predictor Effects On ZTNB Cluster Outcomes

**File:** `supp_fig7_mixing_domain_outcomes`

**What it shows.** A two-panel heatmap of the per-1-SD effect of four
cluster-level mixing predictors (domain-quintile, age, sex, age-sex profile
excess mixing) on (A) the positive cluster-size ZTNB count ratio and (B) the
positive geographic-spread ZTNB count ratio, fit separately in each of eight
SIMD subdomains (rows: Overall, Income, Employment, Education, Health, Access,
Crime, Housing). A single ratio-scale colour bar is shared across panels, with
ticks centred at 1. Each cell is annotated with the ratio. The hurdle component of geographic spread is
omitted from this figure and reported in Supplementary Table 2 because the
heavily imbalanced multiple-datazone outcome combined with the strong
domain-quintile excess-mixing predictor drives the binomial logistic
component toward quasi-separation, producing odds ratios (~22) of comparable
size across all domains but with cluster-robust SE failures for two of them.

**Key visual patterns.** The domain-quintile column saturates the upper
triangle of the colour scale in both panels (ratios ~3.4 in panel A and ~2.8
in panel B), reflecting that domain-quintile excess discordance is the
dominant mixing-side predictor of cluster scale in every per-domain model.
Age mixing is positive but smaller (ratios ~1.5 panel A; ~2.1 panel B); sex
mixing is mildly negative; age-sex joint mixing is near 1. Between-row
variation is very small because each row represents a separate per-subdomain
model in which the age/sex/age-sex mixing predictors are identical observed
variables; only the domain-specific deprivation and domain-excess-mixing pair
differs across rows, and Scotland's SIMD subdomains are highly correlated so
per-domain coefficients land in essentially the same place.

**Suggested results paragraph.**
> Substituting each of the seven SIMD-domain quintile mixing scores for the
> overall SIMD-quintile mixing score, in turn, while retaining age, sex, and
> joint age-sex excess discordance, revealed substantial heterogeneity in the
> strength of domain-quintile mixing as a predictor of cluster scale
> (Supplementary Figure 7, ZTNB components; Supplementary Table 2 for the
> companion geographic-spread hurdle component). Domain-quintile mixing
> associations with cluster size and geographic spread were consistently
> positive across all domains, mirroring the pooled SIMD-quintile result.
> Education- and crime-quintile mixing were among the strongest domain-side
> predictors of cluster scale, whereas access- and housing-quintile mixing
> produced more attenuated estimates — paralleling the deprivation-as-exposure
> domain heterogeneity in Supplementary Figure 5. Crime and Education domains
> in Supplementary Table 2 lack cluster-robust standard errors because the
> window-clustered sandwich variance estimator failed numerically for those
> two hurdle fits (singular Hessian); point estimates remain valid.

---

## Supplementary Figure 8: Size-Adjusted Positive Geographic Spread

**File:** `supp_fig8_deprivation_size_adjusted`

**What it shows.** A coefficient plot comparing the primary positive-count
model for geographic spread with a size-adjusted version that additionally
includes log(cluster size) as a covariate. The SIMD deprivation estimate is
highlighted in both versions to show how the association changes after
conditioning on cluster size.

**Key visual patterns.** The SIMD deprivation estimate for positive
geographic spread visibly flips sign between the primary model (count ratio
0.851, clearly below the null) and the size-adjusted model (count ratio 1.027,
above the null). This sign reversal is the central visual message of the figure.

**Suggested results paragraph.**
> After additionally adjusting the positive-count geographic-spread model for
> cluster size, the overall SIMD deprivation estimate changed direction from
> the primary model (count ratio 0.851, 95% CI 0.792–0.915) to a weakly
> positive association (count ratio 1.027, 95% CI 1.019–1.035, p < 0.001).
> This sign reversal indicates that more deprived clusters are generally
> smaller, and that the unadjusted negative geographic-spread result is
> driven primarily by cluster-size differences. Among clusters of comparable
> size, there is a marginally wider geographic spread with higher deprivation,
> suggesting that deprivation influences cluster scale rather than geographic
> diffusion per se.

---

## Supplementary Figure 9: Deprivation Log-Linear Versus Hurdle/ZTNB Count Models

**File:** `supp_fig9_deprivation_loglinear`

**What it shows.** A side-by-side coefficient comparison for the SIMD
deprivation association with cluster size and geographic spread, contrasting
the log-linear geometric mean ratio model with the hurdle (odds ratio) and
ZTNB (count ratio) components of the two-part main model.

**Key visual patterns.** The log-linear SIMD estimates are both close to 1.000
(cluster size geometric mean ratio 0.992; geographic spread 1.001), with
narrow intervals straddling the null. The hurdle and positive-count
components show clearly differentiated estimates in both directions, making
visually explicit that the log-linear model averages over the structural mass
at the minimum and thereby masks the within-component patterns the two-part
model separates.

**Suggested results paragraph.**
> Log-linear models fitted for comparison produced SIMD deprivation estimates
> substantially attenuated relative to the hurdle/ZTNB main models: geometric
> mean ratios of 0.992 (95% CI 0.987–0.997) for cluster size and 1.001 (95%
> CI 0.996–1.006) for geographic spread. These near-null log-linear estimates
> reflect the averaging of hurdle and positive-count components across the
> structural mass at the count minimum, masking the negative positive-count
> associations visible in the two-part models. This comparison supports the
> decision to use hurdle/ZTNB models as the primary count analytical
> framework.

---

## Supplementary Figure 10: Mixing-Predictor Log-Linear Versus Hurdle/ZTNB Count Models

**File:** `supp_fig10_mixing_loglinear`

**What it shows.** The same comparison as Supplementary Figure 9 but for the
four cluster-level mixing predictors (SIMD-quintile, age, sex, joint
SIMD-age-sex profile excess mixing). Side-by-side log-linear (geometric mean
ratio) and hurdle/ZTNB estimates per 1 SD higher excess-mixing score.

**Key visual patterns.** Log-linear mixing-predictor estimates attenuate
toward 1 relative to the hurdle/ZTNB main models, in parallel with the
deprivation pattern in Supplementary Figure 9. The hurdle/ZTNB components
preserve the strong positive SIMD-quintile and age mixing-predictor effects on
cluster size and geographic spread, while the log-linear summary masks them.

**Suggested results paragraph.**
> Log-linear mixing-predictor count models attenuated the SIMD-quintile and
> age mixing-predictor effects on cluster size and geographic spread toward
> the null, again because the log-linear summary averages over the hurdle and
> positive-count components of the cluster-size and geographic-spread
> distributions. The hurdle/ZTNB main models preserve these positive
> associations, supporting the two-part formulation for the mixing-predictor
> count analyses as well as for the deprivation-as-exposure analyses.

---

## Supplementary Table 1: Wave-Specific Mixing-Predictor Effects On The Geographic-Spread Hurdle (Companion To Figure 4)

**File:** `supp_table_fig4_wave_mixing_hurdle_geographic_spread.csv`

**What it contains.** Window-clustered hurdle (binomial GLM with logit link)
odds ratios with 95% CIs and p-values for the four mixing predictors
(SIMD-quintile, age, sex, joint SIMD-age-sex profile excess mixing) by epidemic
wave. One row per wave × predictor combination. Columns: `wave_group`,
`term`, `coefficient`, `std_error_clustered_by_window`, `ratio`,
`ratio_ci_low`, `ratio_ci_high`, `p_value`, `n_observations`, `n_events`,
`notes`.

**Why a table instead of a heatmap.** The multiple-datazone hurdle outcome is
positive for 88.0% of non-singleton clusters and is therefore heavily
imbalanced. Combined with the strong SIMD-quintile excess-mixing predictor, the
binomial logistic component approaches quasi-separation in some waves,
producing implausibly large adjusted odds ratios (e.g. OR approximately 29,000 in the
Alpha wave, 95% CI 6,450–129,898). Such cells dominate a heatmap and make the
other predictor cells visually unreadable, so the wave-specific hurdle
estimates are reported as a table rather than as a third heatmap panel in
Figure 4.

**Suggested results paragraph.**
> Wave-specific geographic-spread hurdle odds ratios for the four mixing
> predictors are reported in Supplementary Table 1; the SIMD-quintile
> mixing-predictor hurdle OR was clearly positive in every wave but reached
> implausibly large values in Alpha (approximately 29,000), consistent with quasi-separation
> in the heavily imbalanced binary outcome. The age, sex, and joint-profile
> hurdle estimates were of more interpretable magnitude across waves and were
> broadly consistent in direction with the pooled Figure 3 estimates.

---

## Supplementary Table 2: SIMD-Domain Mixing-Predictor Effects On The Geographic-Spread Hurdle (Companion To Supplementary Figure 7)

**File:** `supp_table_fig7_domain_mixing_hurdle_geographic_spread.csv`

**What it contains.** Window-clustered hurdle odds ratios for the four
mixing predictors (domain-quintile, age, sex, age-sex profile excess mixing)
in each per-domain hurdle model. One row per domain × predictor combination.
Same column structure as Supplementary Table 1, with the addition of a
`predictor` key column mapping each domain-specific term to its conceptual
category. The crime and education rows carry the note "cluster-robust SE
unavailable (Hessian singular); point estimate only" because the
window-clustered sandwich variance estimator failed numerically for those two
hurdle fits — the maximum-likelihood point estimates remain valid.

**Why a table instead of a heatmap.** Same reasoning as Supplementary
Table 1: the geographic-spread hurdle outcome is positive in 88.0%
of non-singleton clusters and, when combined with the strong domain-quintile
excess-mixing predictor, produces domain-quintile odds ratios in the range
21-23. These comparable large effects saturate any practical heatmap colour
scale and obscure the smaller age, sex, and age-sex predictor cells.

**Suggested results paragraph.**
> Per-domain geographic-spread hurdle odds ratios for the four mixing
> predictors are reported in Supplementary Table 2. Domain-quintile excess
> mixing produced adjusted hurdle ORs of approximately 22 across all eight
> domain models, with very narrow between-domain variation — a consistency
> that follows directly from the highly correlated SIMD subdomains. Age
> excess discordance was positively associated with the hurdle outcome, sex was
> mildly negative, and joint age-sex was close to the null. Crime and Education rows
> are reported as point estimates only, because the cluster-robust sandwich
> variance estimator failed numerically for those two hurdle fits (singular
> Hessian under the heavy outcome imbalance and near-collinear domain-mixing
> predictor configuration); point estimates remain valid.

---

## Sensitivity Analysis Results

Five sensitivity runs were generated with `part1_sensitivities.sh`. Each run
has its own result tables and manuscript figures so that the primary outputs are
not overwritten:

| Sensitivity | Tables | Figures |
|---|---|---|
| Health-board clustered standard errors | `part1/sensitivity/tables_health_board` | `part1/sensitivity/figures_health_board` |
| Cluster-size positive-count offset | `part1/sensitivity/tables_size_offset` | `part1/sensitivity/figures_size_offset` |
| Index-case SIMD exposure | `part1/sensitivity/tables_index_simd` | `part1/sensitivity/figures_index_simd` |
| 99th-percentile positive-count winsorisation | `part1/sensitivity/tables_winsorise99` | `part1/sensitivity/figures_winsorise99` |
| Approximately non-overlapping windows | `part1/sensitivity/tables_stride3` | `part1/sensitivity/figures_stride3` |

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

This sensitivity adds the log number of sequences in the analysis window as an
offset in the cluster-size positive-count model, changing the estimand from raw
positive cluster size to positive cluster size relative to the analysis-window
sequencing pool.

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
  to 1.035.

Interpretation: the positive cluster-size deprivation estimate is sensitive to
the extreme right tail, while the negative positive geographic-spread result and
the size-adjusted positive results are less sensitive to winsorisation.

### Approximately Non-Overlapping Windows

This sensitivity keeps approximately every third analysis window, reducing the
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

In the main hurdle/ZTNB models (Figure 1), higher cluster-level SIMD
deprivation was not associated with larger or more geographically dispersed
genomic clusters after adjustment for lineage, calendar time, local
incidence, sequencing intensity, and test positivity. Instead, deprivation
was associated with slightly lower odds of being non-singleton and smaller
positive cluster size, and with substantially lower positive geographic
spread. Local epidemic and surveillance conditions were much more strongly
associated with cluster scale: higher incidence, window-level sequencing
proportion, and test positivity were consistently associated with larger and
more geographically dispersed clusters. Mixing analyses (Figure 1 bottom row)
showed a more nuanced socioeconomic pattern. Overall deprivation was not
clearly associated with SIMD-quintile mixing, but was associated with greater
age and joint socio-demographic mixing and lower sex mixing. SIMD-domain
models (Supplementary Figures 4–5) showed that education and crime
deprivation were associated with greater domain-quintile mixing, whereas
access and housing deprivation were associated with lower domain-quintile
mixing. Per-wave outcome models (Figure 2) indicated that deprivation
effects varied over time, with the strongest negative associations seen
during Delta and more heterogeneous patterns in Omicron subwaves.
Reversing the modelling direction (Figure 3; wave-specific in Figure 4
ZTNB components, with the companion geographic-spread hurdle in Supplementary
Table 1), within-cluster excess mixing was itself a substantial predictor of
cluster scale: per 1 SD higher SIMD-quintile excess discordance, positive
cluster size was 3.478-fold higher, positive geographic spread was 3.029-fold
higher, and the odds of multi-datazone spread were 22.107-fold higher. Age
excess discordance was also strongly positively associated with both ZTNB
components. The mixing-predictor effects were directionally stable across
SIMD subdomains (Supplementary Figure 7; Supplementary Table 2 for the
companion hurdle) and across epidemic waves, with the largest magnitudes in
Alpha, Delta, BA.1, and BA.2. Sensitivity analyses supported the qualitative
patterns: health-board clustered standard errors widened several count-
outcome intervals, 99th-percentile winsorisation attenuated the positive
cluster-size association, and the index-case SIMD exposure did not reproduce
the mean-cluster-SIMD positive-count associations.

## Interpretation For The Part 1 Question

The answer to the original question is therefore mixed:

- Socioeconomic deprivation is associated with some cluster outcomes, but not
  in the simple direction of larger or more dispersed clusters.
- Surveillance and epidemic-intensity variables are strongly associated with
  apparent cluster scale, reinforcing the need to adjust for sequencing and
  testing context.
- SIMD domain results show that deprivation is not a single mechanism.
  Access, housing, education, crime, income, employment, and health domains
  capture different kinds of local social structure.
- Demographic mixing is where deprivation-as-exposure signals are clearest,
  especially for age mixing and joint profile mixing.
- Reversing the modelling direction, within-cluster SIMD-quintile and age
  excess discordance are themselves substantial positive predictors of
  cluster size and geographic spread, consistent with a bridging-transmission
  account of how cluster scale arises. The hurdle component of geographic
  spread is reported as Supplementary Tables 1 and 2 rather than as a heatmap
  because the heavily imbalanced binary outcome combines with the strong
  mixing predictor to drive the binomial logistic component toward quasi-
  separation, producing odds ratios that should be read as evidence of strong
  direction rather than as interpretable effect sizes.
- Wave-specific analyses (both the deprivation-as-exposure and the
  mixing-as-predictor lines) suggest that lineage and epidemic context modify
  these associations.
- Sensitivity analyses reinforce the need to separate point-estimate
  stability from statistical precision: directions are mostly stable, but
  inference is weaker under health-board clustered standard errors and when
  the heaviest positive-count tail is capped.

These results support a cautious conclusion: deprivation and local
surveillance conditions are associated with observed genomic cluster
structure, but the associations are outcome-specific, domain-specific, and
wave-specific rather than a single monotonic deprivation effect. Cluster
scale in genomic surveillance reflects both the deprivation profile of the
cases involved and the bridging structure of contact and transmission that
drew them into the same cluster.
