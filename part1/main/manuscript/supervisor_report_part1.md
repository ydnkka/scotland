# Scottish SARS-CoV-2 Genomic Clustering: Part 1 Summary Report

Prepared: 6 May 2026

## Purpose

This Part 1 analysis asks whether socioeconomic deprivation and local
surveillance conditions are associated with the structure of SARS-CoV-2 genomic
clusters in Scotland. The main question is whether clusters linked to more
deprived areas are larger, longer-lasting, more geographically dispersed, or
more socially/demographically mixed after accounting for lineage, calendar time,
local incidence, sequencing intensity, and test positivity.

The core conclusion is cautious: deprivation is associated with some aspects of
cluster structure, but not in the simple direction of larger, longer-lasting, or
more geographically dispersed clusters. Surveillance and epidemic context are
strong drivers of apparent cluster scale, while deprivation signals are clearest
in demographic mixing and in selected SIMD domains.

## Analysis Population

The analysis uses the processed Scottish SARS-CoV-2 clustering dataset:

- Input dataset: `data/processed/scotland_clustering_analysis_dataset.parquet`
- QC filter: `nextclade_qc == "good"`
- Primary Leiden clustering resolution: `0.3`
- Unit of analysis: one inferred genomic cluster at the primary resolution
- Sequence rows used: 789,347 (281,320 unique sequences)
- Inferred clusters: 193,160
- Sliding time windows: 134
- Raw Pango lineages: 788
- Modelled lineage levels after rare-lineage pooling: 183
- Non-singleton clusters used for mixing models: 84,067

A single primary Leiden resolution was used so that the same sampled genomes
were not repeatedly analysed across multiple clustering resolutions.

## Key Variables

### Main Exposure

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

### Count Outcomes

Three cluster-level count outcomes were analysed:

- **Cluster size:** number of sequences in the inferred genomic cluster.
- **Duration:** days between the earliest and latest sampled case in the
  cluster.
- **Geographic spread:** number of distinct datazones represented in the
  cluster.

These outcomes have large structural masses at their minimum values. In the
main cluster table:

- 56.5% of clusters are singletons.
- 63.1% have duration zero days.
- 61.7% are observed in a single datazone.
- Maximum observed values are 2,792 sequences, 19 days, and 2,100 datazones.

### Mixing Outcomes

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

### Adjustment Variables

The models adjusted for:

- Local cumulative incidence
- Local sequencing fraction
- Window-level sequencing proportion
- Local test positivity
- Pango lineage, with rare lineages pooled into `Other rare lineages`
- Calendar time using an 8 df B-spline over `window_idx`

Mixing models additionally adjusted for cluster size.

## Modelling Overview

### Count Models

Cluster size, duration, and geographic spread were modelled using two-part
hurdle models:

- A binary hurdle component: binomial GLM with logit link.
- A positive-count component: zero-truncated negative binomial model.

Outcome-specific definitions:

- **Cluster size:** hurdle is `cluster_size > 1`; positive count is
  `cluster_size - 1`.
- **Duration:** hurdle is `duration_days > 0`; positive count is
  `duration_days`.
- **Geographic spread:** hurdle is `cluster_n_datazones > 1`; positive count is
  `cluster_n_datazones - 1`.

This formulation separates the probability of exceeding the structural minimum
from the magnitude of the outcome among clusters that exceed it.

### Mixing Models

Mixing outcomes were analysed among non-singleton clusters using linear models
with the same adjustment set as the count models, plus cluster size. The
coefficient is interpreted as the adjusted change in excess mixing, in
percentage points, per 1 SD higher covariate.

### Extensions And Sensitivities

Additional analyses included:

- SIMD-domain versions of the count and mixing models.
- Wave-specific versions of the main count models.
- Wave-specific domain-demographic mixing models.
- Size-adjusted positive-count sensitivities for duration and geographic spread.
- Log-linear versus hurdle/ZTNB comparisons.
- Health-board clustered standard errors.
- Index-case SIMD exposure instead of mean cluster SIMD.
- 99th-percentile positive-count winsorisation.
- Approximately non-overlapping window sensitivity.

## Recommended Figure Set For Supervisors

The following figures are the clearest set for a supervisor-facing report or
meeting.

### 1. Outcome And Mixing Distributions

Use: `supp_fig1_outcome_distributions`

This figure shows why the modelling uses hurdle/ZTNB count models and why
mixing is analysed among non-singleton clusters. The first row shows cluster
size, duration, and distinct datazones among clusters with size greater than 1.
The second row shows age, sex, and deprivation excess mixing distributions.

![Supplementary Figure 1](figures/supp_fig1_outcome_distributions.png)

### 2. Main Cluster Outcome Models

Use: `fig1_main_cluster_outcomes`

This is the core figure for the count outcomes. It shows adjusted odds ratios
and positive-count ratios for cluster size, duration, and geographic spread.

![Figure 1](figures/fig1_main_cluster_outcomes.png)

### 3. Main Cluster Mixing Models

Use: `fig2_main_cluster_mixing`

This is the core figure for social and demographic composition. It shows how
SIMD deprivation and surveillance covariates relate to SIMD, age, sex, and joint
profile excess mixing.

![Figure 2](figures/fig2_main_cluster_mixing.png)

### 4. SIMD-Domain Mixing

Use: `fig3_simd_domain_mixing`

This figure shows that SIMD domains do not behave as interchangeable measures
of deprivation. It is important for interpretation because education, crime,
access, and housing domains show different mixing patterns.

![Figure 3](figures/fig3_simd_domain_mixing.png)

### 5. Wave-Specific Cluster Outcomes

Use: `fig4_wave_specific_cluster_outcomes`

This figure shows that deprivation effects vary by epidemic wave. It helps
avoid overinterpreting the pooled deprivation effect as a single stable process.

![Figure 4](figures/fig4_wave_specific_cluster_outcomes.png)

### 6. Size-Adjusted Positive Counts

Use: `supp_fig2_size_adjusted_positive_counts`

This figure is useful as a sensitivity result. It shows that after additionally
adjusting positive duration and geographic spread models for cluster size,
deprivation has weak positive associations with these outcomes.

![Supplementary Figure 2](figures/supp_fig2_size_adjusted_positive_counts.png)

## Key Findings

### 1. More Deprived Clusters Were Not Generally Larger Or More Dispersed

The main pooled count models do not support the simple hypothesis that higher
SIMD deprivation is associated with larger, longer-lasting, or more
geographically dispersed genomic clusters.

For SIMD deprivation in the main models:

- Cluster size hurdle: OR 0.971, 95% CI 0.960 to 0.983.
- Positive cluster size: count ratio 0.926, 95% CI 0.869 to 0.987.
- Duration hurdle: OR 0.992, 95% CI 0.979 to 1.004.
- Positive duration: count ratio 1.003, 95% CI 0.995 to 1.011.
- Geographic spread hurdle: OR 1.004, 95% CI 0.992 to 1.016.
- Positive geographic spread: count ratio 0.851, 95% CI 0.792 to 0.915.

The strongest deprivation-associated count result is therefore lower positive
geographic spread, not greater spread.

### 2. Surveillance And Epidemic Context Strongly Shape Apparent Cluster Scale

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

### 3. Deprivation Signals Are Clearer In Demographic Mixing Than In Cluster Size

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

### 4. SIMD Domains Behave Differently

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

### 5. Deprivation Effects Vary By Epidemic Wave

Wave-specific models show that deprivation associations are not stable across
the pandemic.

Delta showed the clearest negative deprivation associations:

- Lower odds of non-singleton clusters.
- Smaller positive cluster sizes.
- Lower odds of positive duration.
- Lower odds of multi-datazone spread.
- Lower positive geographic spread.

Later Omicron subwaves were more heterogeneous. BA.2 and BA.4 showed some
positive positive-count associations, but BA.4 had a much smaller sample and
should be interpreted cautiously.

Interpretation: the pooled deprivation effect should be read as an average over
changing epidemic contexts rather than as a stable biological or social effect.

### 6. Sensitivity Analyses Support A Cautious Interpretation

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

## Interpretation For Supervisors

The analysis provides evidence that observed genomic cluster structure is shaped
by both social context and surveillance context, but the strongest and most
consistent drivers of apparent cluster size and geographic spread are epidemic
and sampling variables.

Higher overall SIMD deprivation does not produce a simple pattern of larger,
longer, or more geographically dispersed clusters. Instead, deprivation is linked
to smaller positive cluster size and lower positive geographic spread in the main
pooled models. However, among clusters of comparable size, size-adjusted
sensitivities suggest weak positive associations with duration and geographic
spread.

The most interpretable deprivation signal is in cluster composition: more
deprived clusters show higher age mixing and joint socio-demographic profile
mixing, lower sex mixing, and little evidence of increased SIMD-quintile mixing.
Domain analyses show that education, crime, access, and housing deprivation
behave differently, reinforcing that overall SIMD is not a single mechanism.

The results should be framed as outcome-specific and context-dependent rather
than as evidence for one monotonic deprivation effect. The findings also
highlight the importance of adjusting for surveillance variables when using
genomic clusters as epidemiological outcomes.

## Suggested High-Level Takeaway

In Scottish SARS-CoV-2 genomic clusters, socioeconomic deprivation was not
associated with generally larger or more geographically dispersed clusters after
adjustment for lineage, calendar time, local incidence, sequencing intensity,
and test positivity. Apparent cluster scale was more strongly related to local
epidemic and surveillance conditions. Deprivation-related signals were clearest
in demographic composition: higher deprivation was associated with greater age
mixing, lower sex mixing, and greater joint socio-demographic profile mixing,
with substantial heterogeneity across SIMD domains and epidemic waves.

