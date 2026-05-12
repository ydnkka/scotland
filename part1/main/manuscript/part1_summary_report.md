# Scottish SARS-CoV-2 Genomic Clustering: Part 1 Summary Report

Prepared: 6 May 2026

## 1. Introduction

This Part 1 analysis asks whether socioeconomic deprivation and local
surveillance conditions are associated with the structure of SARS-CoV-2 genomic
clusters in Scotland. The main question is whether clusters linked to more
deprived areas are larger, more geographically dispersed, or
more socially/demographically mixed after accounting for lineage, calendar time,
local incidence, sequencing intensity, and test positivity.

The core conclusion is cautious: deprivation is associated with some aspects of
cluster structure, but not in the simple direction of larger or
more geographically dispersed clusters. Surveillance and epidemic context are
strong drivers of apparent cluster scale, while deprivation signals are clearest
in demographic mixing and in selected SIMD domains.

## 2. Analysis Population

The analysis used a primary Leiden clustering resolution of `0.3`. Below is some
summary information about the dataset:

- QC filter: `nextclade_qc == "good"`
- Unit of analysis: one inferred genomic cluster at the primary resolution within calendar window
- Sequence rows used: 789,347 (281,320 unique sequences)
- Inferred clusters before model-field filtering: 193,160
- Cluster rows used in the main regression models: 193,112
- Sliding time windows: 134
- Raw Pango lineages: 788
- Modelled lineage levels after rare-lineage pooling: 183
- Non-singleton clusters used for mixing models: 84,067

A single primary Leiden resolution was used because from the simulation study, there is minimal variation
between clusters at the resolutions 0.1–0.6.

## 3. Key Variables

### 3.1. Main Exposure

The main exposure is mean cluster-level SIMD deprivation. SIMD rank was
transformed so that higher values represent greater deprivation, then
standardised. The main pooled models use overall SIMD deprivation; domain
analyses separately examine:

- Overall SIMD
- Income deprivation
- Employment deprivation
- Education deprivation
- Health deprivation
- Geographic access deprivation
- Crime deprivation
- Housing deprivation

### 3.2. Count Outcomes

Three cluster-level count descriptors were constructed:

- **Cluster size:** number of sequences in the inferred genomic cluster.
- **Duration:** days between the earliest and latest sampled case in the
  cluster. This is retained descriptively, but not modelled in the current main count analysis because the fixed three-week clustering windows mechanically
  constrain the observed span.
- **Geographic spread:** number of distinct datazones represented in the
  cluster.

These outcomes have large structural masses at their minimum values:

- 56.5% of clusters are singletons.
- 63.1% have duration zero days.
- 61.7% are observed in a single datazone.
- Maximum observed values are 2,792 sequences, 19 days, and 2,100 datazones.

### 3.3. Mixing Outcomes

Within-cluster mixing was measured as observed pairwise discordance minus
expected pairwise discordance among sampled cases from the same lineage and
calendar window. Positive values mean the cluster is more mixed than expected;
negative values mean it is more homogeneous than expected.

The main mixing outcomes were:

- SIMD quintile excess mixing
- Age-band excess mixing
- Sex excess mixing
- Joint SIMD-age-sex profile excess mixing

In the domain analysis, domain-quintile excess mixing was also estimated for
each SIMD domain.

### 3.4. Adjustment Variables

The models adjusted for:

- Local cumulative incidence
- Local sequencing fraction
- Window-level sequencing proportion
- Local test positivity
- Pango lineage, with rare lineages pooled into `Other` rare lineages
- Calendar time using an 8 df B-spline over `window_idx`

Mixing models additionally adjusted for cluster size.

## 4. Modelling Overview

### 4.1. Count Models

Cluster size and geographic spread were modelled using two-part
hurdle models:

- A binary hurdle component: binomial GLM with logit link.
- A positive-count component: zero-truncated negative binomial model.

Outcome-specific definitions:

- **Cluster size:** hurdle is `cluster_size > 1`; positive count is
  `cluster_size - 1`.
- **Geographic spread:** hurdle is `cluster_n_datazones > 1`; positive count is
  `cluster_n_datazones - 1`.

This formulation separates the probability of exceeding the structural minimum
from the magnitude of the outcome among clusters that exceed it.

Duration remains in the descriptive summaries and the supplementary outcome
distribution figure only.

### 4.2. Mixing Models

Mixing outcomes were analysed among non-singleton clusters using linear models
with the same adjustment set as the count models, plus cluster size. The
coefficient is interpreted as the adjusted change in excess mixing, in
percentage points, per 1 SD higher covariate.

### 4.3. Mixing-Predictor Count Models

A second line of inquiry reverses the modelling direction. The four
cluster-level excess-mixing scores (SIMD-quintile, age, sex, joint
SIMD-age-sex profile) are entered as predictors of cluster size and
geographic spread, alongside SIMD deprivation, the four surveillance
covariates, lineage fixed effects, and the calendar B-spline. The model
form is otherwise identical to the deprivation-as-exposure count models,
retaining the hurdle/ZTNB specification. Because excess discordance is
undefined for singletons, these models are fitted among non-singleton
clusters only, so the cluster-size hurdle is not estimable.

### 4.4. Extensions And Sensitivities

Additional analyses included:

- SIMD-domain versions of the count and mixing models (deprivation as
  exposure).
- Wave-specific versions of the main count models (deprivation as exposure).
- Wave-specific domain-demographic mixing models.
- Wave-specific and SIMD-domain versions of the mixing-predictor count
  models.
- Size-adjusted positive-count sensitivity for geographic spread.
- Log-linear versus hurdle/ZTNB comparisons (for both lines).
- Health-board clustered standard errors.
- Index-case SIMD exposure instead of mean cluster SIMD.
- 99th-percentile positive-count winsorisation.
- Approximately non-overlapping window sensitivity.

## 5. Main Figures

The Part 1 manuscript organises results around two complementary lines of
inquiry. Lines 1–2 of the figures cover the deprivation-as-exposure analyses;
lines 3–4 cover the reverse direction, with the four cluster-level
excess-mixing scores entered as predictors of cluster size and geographic
spread. Supplementary figures and tables extend each main figure to SIMD
subdomain breakdowns, sensitivities, and model-form comparators.

### 5.1. Figure 1: Deprivation As Exposure — Overall Effects On Cluster Outcomes And Mixing

This is the core figure for the deprivation-as-exposure analyses. The top row
shows adjusted odds ratios (hurdle components) and ZTNB count ratios (positive
components) for cluster size and geographic spread; the bottom row shows
adjusted percentage-point changes in observed-minus-expected excess mixing for
SIMD-quintile, age, sex, and joint sociodemographic profile.

![Figure 1](figures/fig1_deprivation_overall.png)

### 5.2. Figure 2: Deprivation As Exposure — Wave-Specific Effects On Cluster Outcomes

This figure shows that deprivation effects on the four count-model components
vary by epidemic wave. It helps avoid overinterpreting the pooled deprivation
effect as a single stable process. Delta shows the clearest negative
associations; BA.2 and BA.4 contrast with positive positive-count estimates.

![Figure 2](figures/fig2_deprivation_wave_specific.png)

### 5.3. Figure 3: Excess Mixing As Predictor — Overall Effects On Cluster Scale

Three-panel coefficient plot reversing the modelling direction: the four
cluster-level excess-mixing scores (SIMD-quintile, age, sex, joint profile)
enter as predictors alongside SIMD deprivation and the surveillance
covariates. SIMD-quintile and age excess discordance are strongly positively
associated with positive cluster size and positive geographic spread; joint
profile mixing is negatively associated. The cluster-size hurdle is omitted
because mixing scores are undefined for singletons.

![Figure 3](figures/fig3_mixing_overall.png)

### 5.4. Figure 4: Excess Mixing As Predictor — Wave-Specific Effects On Cluster Scale (ZTNB)

Two-panel heatmap of per-wave mixing-predictor effects on the ZTNB cluster-size
and ZTNB geographic-spread components. A single shared ratio-scale colour bar
spans 0.2–5 with `extend="both"` triangles; cells are annotated with the raw
count ratio. The geographic-spread hurdle component is reported separately in
Supplementary Table 1 because the heavily imbalanced binary outcome combines
with the strong SIMD-excess-mixing predictor to produce implausibly large
adjusted odds ratios (~29,000 in Alpha) that obscure heatmap interpretation.

![Figure 4](figures/fig4_mixing_wave_specific.png)

### 5.5. Supplementary Highlights

Two supplementary figures are particularly useful for reading the main
results in context:

- **Supplementary Figure 1** (outcome distributions, `supp_fig1_outcome_distributions.png`)
  motivates the hurdle/ZTNB formulation by showing the structural mass at the
  count minima and the long right tails of the non-singleton distributions.
- **Supplementary Figure 8** (size-adjusted positive counts, `supp_fig8_deprivation_size_adjusted.png`)
  shows that the unadjusted negative SIMD-deprivation association with
  positive geographic spread flips to a weakly positive association after
  conditioning on cluster size — deprivation affects cluster scale, not
  geographic diffusion at comparable cluster sizes.

Full panel-by-panel descriptions of all 10 supplementary figures and the two
supplementary tables (Supplementary Table 1, companion to Figure 4;
Supplementary Table 2, companion to Supplementary Figure 7) are in
`part1_results_and_figures_description.md`; concise captions are also
available in `figures/part1_supplementary_files.md`.

## 6. Key Findings

### 6.1. More Deprived Clusters Were Not Generally Larger Or More Dispersed

The main pooled count models do not support the simple hypothesis that higher
SIMD deprivation is associated with larger or more geographically dispersed
genomic clusters.

For SIMD deprivation in the main models:

- Cluster size hurdle: OR 0.971, 95% CI 0.960 to 0.983.
- Positive cluster size: count ratio 0.926, 95% CI 0.869 to 0.987.
- Geographic spread hurdle: OR 1.004, 95% CI 0.992 to 1.016.
- Positive geographic spread: count ratio 0.851, 95% CI 0.792 to 0.915.

The strongest deprivation-associated count result is therefore lower positive
geographic spread, not greater spread.

### 6.2. Surveillance And Epidemic Context Strongly Shape Apparent Cluster Scale

Local incidence, sequencing intensity, window-level sequencing proportion, and
test positivity were more consistently associated with larger apparent clusters
than SIMD deprivation.

Examples from the positive-count components:

- Local sequencing fraction was strongly associated with positive cluster size
  and positive geographic spread.
- Window-level sequencing proportion was positively associated with all count
  components.
- Test positivity showed large positive associations with positive cluster size
  and geographic spread.

Interpretation: reconstructed genomic clusters are partly a function of
transmission, but also of how many infections are detected, sequenced, and
linked within a time window and locality.

### 6.3. Deprivation Signals Are Clearer In Demographic Mixing Than In Cluster Size

Overall SIMD deprivation was not clearly associated with SIMD-quintile excess
mixing itself:

- SIMD-quintile excess mixing: +0.31 percentage points, 95% CI -0.18 to 0.80.

However, deprivation was associated with several demographic mixing outcomes:

- Age excess mixing: +1.66 percentage points, 95% CI 1.29 to 2.03.
- Sex excess mixing: -0.78 percentage points, 95% CI -1.16 to -0.39.
- Joint SIMD-age-sex profile excess mixing: +0.48 percentage points, 95% CI
  0.29 to 0.67.

This suggests that deprivation may be more strongly related to the demographic
composition of clusters than to simple mixing across SIMD quintiles.

### 6.4. SIMD Domains Behave Differently

The SIMD-domain analyses show that the overall SIMD result masks domain-specific
patterns.

For domain-quintile mixing:

- Education deprivation was associated with higher domain-quintile mixing.
- Crime deprivation was associated with higher domain-quintile mixing.
- Access deprivation was associated with lower domain-quintile mixing.
- Housing deprivation was associated with lower domain-quintile mixing.

For demographic mixing:

- Age mixing was higher for most deprivation domains.
- Sex mixing tended to be lower with higher deprivation across several domains.
- Joint age-sex mixing was generally higher for overall, income, employment,
  education, health, and housing deprivation.
- Access deprivation was an exception in several analyses.

Interpretation: deprivation should not be treated as a single homogeneous
mechanism. Different SIMD domains likely capture different social and spatial
processes relevant to transmission, observation, and cluster composition.

### 6.5. Deprivation Effects Vary By Epidemic Wave

Wave-specific models show that deprivation associations are not stable across
the pandemic.

Delta showed the clearest negative deprivation associations:

- Lower odds of non-singleton clusters.
- Smaller positive cluster sizes.
- Lower odds of multi-datazone spread.
- Lower positive geographic spread.

Later Omicron subwaves were more heterogeneous. BA.2 and BA.4 showed some
positive positive-count associations, but BA.4 had a much smaller sample and
should be interpreted cautiously.

Interpretation: the pooled deprivation effect should be read as an average over
changing epidemic contexts rather than as a stable biological or social effect.

### 6.6. Within-Cluster Mixing Strongly Predicted Cluster Scale

Reversing the modelling direction (Figure 3; wave-specific in Figure 4),
within-cluster excess mixing was itself a substantial predictor of cluster
size and geographic spread among non-singleton clusters. Per 1 SD higher
SIMD-quintile excess discordance:

- Positive cluster size was 3.48-fold higher (95% CI 3.25–3.73).
- Positive geographic spread was 3.03-fold higher (95% CI 2.79–3.29).
- Odds of multi-datazone spread were 22.1-fold higher (95% CI 19.0–25.8).

Age excess discordance was also strongly positively associated with both ZTNB
components (positive cluster size CR 1.67; positive geographic spread CR
1.97). Joint SIMD-age-sex profile excess discordance was negatively associated
with both components.

Wave-specific extension (Figure 4, ZTNB only) showed the same directional
pattern, with the largest mixing-predictor effects in waves with the deepest
cluster recruitment (Alpha, Delta, BA.1, BA.2). The geographic-spread hurdle
component is reported in Supplementary Table 1 rather than as a heatmap
because the heavily imbalanced binary outcome (`datazones_gt1` positive in
88% of clusters) combined with the strong SIMD-excess-mixing predictor
produced implausibly large adjusted odds ratios (~29,000 in Alpha).

Interpretation: clusters that bridge across SIMD quintiles and age bands more
than the lineage-window baseline expectation are detected as substantially
larger and more geographically dispersed, consistent with a bridging-
transmission account of cluster scale. The pattern is associational rather
than causal.

### 6.7. SIMD-Domain Mixing-Predictor Heterogeneity

Substituting each SIMD subdomain's quintile excess-mixing score for the
overall SIMD-quintile score (Supplementary Figure 7) confirmed that all
seven SIMD subdomains gave positive domain-quintile mixing-predictor effects
on cluster size and geographic spread (Supplementary Figure 7) with very
similar magnitudes — the per-domain models share the same age/sex/joint
demographic mixing variables and only differ in which deprivation/excess-
mixing pair is included, so the results reflect the consistency of the
mixing-predictor finding across SIMD subdomains. The companion
Supplementary Table 2 reports the corresponding geographic-spread hurdle
results; crime and education rows are reported as point estimates only
because the cluster-robust sandwich variance estimator failed numerically
for those two hurdle fits.

### 6.8. Sensitivity Analyses Support A Cautious Interpretation

Sensitivity analyses generally supported the qualitative pattern but showed that
some inferential strength depends on modelling choices.

Key points:

- Health-board clustered standard errors widened several confidence intervals,
  especially for count outcomes.
- Winsorising positive counts attenuated the positive cluster-size association
  but did not remove the lower positive geographic-spread result.
- Index-case SIMD did not reproduce the mean-cluster-SIMD positive-count
  associations, suggesting that mean cluster deprivation captures whole-cluster
  composition rather than simply the deprivation context of the earliest
  observed case.
- Approximately non-overlapping windows reduced precision but retained the
  broad qualitative pattern.
- Age and sex mixing findings were among the more stable results.

## 7. Interpretation

The analysis provides evidence that observed genomic cluster structure is
shaped jointly by social context, surveillance context, and the bridging
structure of within-cluster mixing. Three threads run through the results.

First, higher overall SIMD deprivation does not produce a simple pattern of
larger or more geographically dispersed clusters. Deprivation is linked to
smaller positive cluster size and lower positive geographic spread in the
main pooled models. Among clusters of comparable size, the size-adjusted
geographic-spread sensitivity suggests only a weak positive association with
spread. The strongest and most consistent drivers of apparent cluster size
and geographic spread are surveillance and epidemic-intensity variables.

Second, the most interpretable deprivation signal in the
deprivation-as-exposure analyses is in cluster composition: more deprived
clusters show higher age mixing and joint socio-demographic profile mixing,
lower sex mixing, and little evidence of increased SIMD-quintile mixing.
Domain analyses show that education, crime, access, and housing deprivation
behave differently, reinforcing that overall SIMD is not a single mechanism.

Third, reversing the modelling direction (mixing as predictor of cluster
scale) revealed that clusters which bridge across SIMD quintiles and age
bands more than the lineage-window baseline expectation are detected as
substantially larger and more geographically dispersed. SIMD-quintile excess
discordance was the strongest mixing-side predictor in every panel of
Figure 3, with effects 3- to 22-fold larger in magnitude per 1 SD than the
deprivation-as-exposure estimates. The hurdle component of geographic spread
is reported as a supplementary table (Supplementary Tables 1 and 2) rather
than as a heatmap because the heavily imbalanced binary outcome combined with
the strong mixing predictor drives the binomial logistic component toward
quasi-separation, producing implausibly large adjusted odds ratios that
should be read as evidence of strong direction rather than as interpretable
effect sizes.

The results should be framed as outcome-specific and context-dependent rather
than as evidence for one monotonic deprivation effect. The findings also
highlight the importance of adjusting for surveillance variables when using
genomic clusters as epidemiological outcomes, and the value of treating
within-cluster mixing structure as both an outcome (in the
deprivation-as-exposure line) and a predictor (in the reverse line).

## 8. Takeaway

In Scottish SARS-CoV-2 genomic clusters, socioeconomic deprivation was not
associated with generally larger or more geographically dispersed clusters
after adjustment for lineage, calendar time, local incidence, sequencing
intensity, and test positivity. Apparent cluster scale was more strongly
related to local epidemic and surveillance conditions. Deprivation-related
signals were clearest in cluster composition: higher deprivation was
associated with greater age mixing, lower sex mixing, and greater joint
socio-demographic profile mixing, with substantial heterogeneity across SIMD
domains and epidemic waves. Reversing the modelling direction, within-cluster
SIMD-quintile and age excess discordance were themselves substantial
positive predictors of cluster size and geographic spread, consistent with a
bridging-transmission account of how cluster scale arises. Cluster scale in
genomic surveillance therefore reflects both the deprivation profile of the
cases involved and the bridging structure of contact and transmission that
drew them into the same cluster.
