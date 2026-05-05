# Part 1 Main Results And Figure Guide

Generated for the main Part 1 analysis on 5 May 2026.

## Analysis Question

Part 1 asks whether socioeconomic deprivation and local surveillance conditions
are associated with SARS-CoV-2 genomic cluster characteristics in Scotland after
accounting for lineage, calendar time, local incidence, sequencing intensity,
test positivity, and clustering resolution.

The specific outcomes are:

- cluster size
- cluster duration
- geographic spread, measured by the number of distinct datazones in a cluster
- observed-minus-expected socioeconomic and demographic mixing within clusters

The main manuscript figures are generated from `part1/main/tables` and
`part1/main/cache`, using the publication styling in `utils/style.py`.

## Analysis Population And Modelling Frame

The main analysis uses good-QC sequences at Leiden resolution 0.3. This produces
789,347 sequence rows and 193,112 inferred genomic clusters across 134 time
windows. Rare Pango lineages are pooled at a threshold of 50 clusters, giving
183 lineage model levels in the main pooled analysis.

The count outcomes have large structural masses at their minimum values, so the
main count analysis uses two-part hurdle models:

- a binomial hurdle component for whether a cluster exceeds the structural
  minimum
- a zero-truncated negative-binomial positive-count component for the outcome
  among clusters exceeding the structural minimum

For cluster size, the hurdle is `cluster_size > 1` and the positive count is
`cluster_size - 1`. For duration, the hurdle is `duration_days > 0` and the
positive count is `duration_days`. For geographic spread, the hurdle is
`cluster_n_datazones > 1` and the positive count is `cluster_n_datazones - 1`.

Mixing outcomes are fitted among non-singleton clusters. They are defined as
observed pairwise discordance within a cluster minus expected discordance among
sampled cases from the same lineage and calendar window. Positive values mean
more mixing than expected; negative values mean more homogeneity than expected.

Throughout the figures, higher `SIMD deprivation` means a more deprived cluster
on the transformed and standardised SIMD rank scale.

## Overall Results Summary

The main analysis does not support the simple hypothesis that more deprived
clusters are generally larger, longer lasting, and more geographically dispersed.
In the pooled main models, higher SIMD deprivation is associated with slightly
lower odds of being non-singleton, smaller positive cluster sizes, and
substantially lower positive geographic spread. Duration is not materially
associated with overall SIMD deprivation in the primary count models.

Surveillance and epidemic-intensity variables show much stronger and more
consistent associations with cluster outcomes. Higher local incidence, higher
window-level sequencing proportion, and higher test positivity are associated
with larger, longer-lasting, and more geographically dispersed clusters. Local
sequencing fraction is especially strongly associated with larger positive
cluster size and positive geographic spread, consistent with the idea that
local sampling intensity affects the apparent scale of reconstructed genomic
clusters.

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

## Main Figure 1: Main Cluster Outcome Models

File stem: `fig1_main_cluster_outcomes`

This figure shows the main hurdle and zero-truncated negative-binomial results
for cluster size, duration, and geographic spread. Each point is an adjusted
ratio per 1 SD higher covariate. For hurdle panels, the ratio is an odds ratio.
For positive-count panels, the ratio is a ZTNB count ratio.

### SIMD Deprivation

In the pooled main model:

- Cluster size hurdle: OR 0.971, 95% CI 0.960 to 0.983, p = 1.43e-06.
- Positive cluster size: ZTNB count ratio 0.926, 95% CI 0.891 to 0.963,
  p = 0.000126.
- Duration hurdle: OR 0.992, 95% CI 0.979 to 1.004, p = 0.192.
- Positive duration: ZTNB count ratio 1.003, 95% CI 0.994 to 1.012,
  p = 0.520.
- Geographic spread hurdle: OR 1.004, 95% CI 0.992 to 1.016, p = 0.522.
- Positive geographic spread: ZTNB count ratio 0.851, 95% CI 0.813 to 0.891,
  p = 3.36e-12.

The strongest SIMD deprivation result is therefore not larger clusters, but
smaller positive cluster size and lower positive geographic spread after
adjustment.

### Incidence And Surveillance Conditions

Local incidence, local sequencing fraction, window sequencing proportion, and
test positivity are much more consistently associated with larger apparent
clusters:

- Higher local incidence is associated with higher odds of exceeding the
  structural minimum for all three count outcomes and with larger positive
  cluster size and geographic spread.
- Higher local sequencing fraction is strongly associated with positive cluster
  size, count ratio 3.24, and positive geographic spread, count ratio 2.27.
- Higher window sequencing proportion is positively associated with every count
  component.
- Higher test positivity has the largest and most consistent positive
  associations, including positive cluster size count ratio 2.65 and positive
  geographic spread count ratio 3.00.

Interpretation: local epidemic intensity and sampling conditions have a major
relationship with reconstructed cluster scale. These covariates are not merely
technical nuisances; they shape what a genomic cluster looks like in the data.

## Main Figure 2: Main Cluster Mixing Models

File stem: `fig2_main_cluster_mixing`

This figure shows adjusted changes in excess mixing, in percentage points per
1 SD higher covariate. The outcomes are SIMD-quintile mixing, age-band mixing,
sex mixing, and joint SIMD-age-sex profile mixing.

### SIMD Deprivation

Overall SIMD deprivation has outcome-specific mixing associations:

- SIMD-quintile excess mixing: +0.31 pp, 95% CI -0.18 to 0.80, p = 0.208.
- Age excess mixing: +1.66 pp, 95% CI 1.29 to 2.03, p = 2.57e-18.
- Sex excess mixing: -0.78 pp, 95% CI -1.16 to -0.39, p = 8.57e-05.
- Joint SIMD-age-sex profile excess mixing: +0.48 pp, 95% CI 0.29 to 0.67,
  p = 8.27e-07.

Thus, deprivation is not clearly associated with more mixing across SIMD
quintiles themselves, but it is associated with more age mixing and more joint
socio-demographic profile mixing.

### Other Covariates

Cluster size is a strong driver of mixing estimates. Larger clusters have more
SIMD, age, and joint-profile excess mixing, but lower sex excess mixing. Local
incidence is positively associated with SIMD, age, and joint-profile mixing.
Sequencing and test-positivity covariates generally move SIMD and joint-profile
mixing downward, which is consistent with surveillance conditions influencing
the observed composition of genomic clusters.

## Main Figure 3: SIMD-Domain Mixing

File stem: `fig3_simd_domain_mixing`

This figure compares the overall SIMD deprivation measure with individual SIMD
domains. It has four panels: domain-quintile mixing, age mixing, sex mixing, and
joint age-sex mixing.

### Domain-Quintile Mixing

The domain-specific effects are not interchangeable:

- Education deprivation is associated with higher domain-quintile excess mixing:
  +1.17 pp, 95% CI 0.64 to 1.69.
- Crime deprivation is also associated with higher domain-quintile excess
  mixing: +1.09 pp, 95% CI 0.55 to 1.64.
- Access deprivation is associated with lower domain-quintile excess mixing:
  -2.08 pp, 95% CI -2.50 to -1.66.
- Housing deprivation is also associated with lower domain-quintile excess
  mixing: -1.19 pp, 95% CI -1.77 to -0.62.
- Overall, income, employment, and health deprivation have smaller or less
  certain domain-quintile mixing estimates.

This suggests that the "SIMD effect" depends on the social process represented
by the domain. Education and crime domains behave differently from access and
housing domains.

### Domain-Demographic Mixing

Age mixing is higher for most deprivation domains, including overall SIMD,
income, employment, education, health, crime, and housing. Access deprivation is
the exception and is associated with lower age mixing.

Sex mixing tends to be lower with higher deprivation for overall SIMD and most
domains. Joint age-sex mixing is generally higher for overall, income,
employment, education, health, and housing deprivation, but weaker for access
and crime.

## Main Figure 4: Wave-Specific SIMD Effects On Cluster Outcomes

File stem: `fig4_wave_specific_cluster_outcomes`

This figure shows the per-wave overall SIMD deprivation effect on cluster size,
duration, and geographic spread. It uses the same hurdle/ZTNB framing as the
main count analysis, but models are fitted separately by wave.

The wave-specific outcome models adjust for the same main covariates: calendar
time spline, local incidence, local sequencing fraction, window sequencing
proportion, test positivity, and SIMD deprivation. Within-wave lineage dummies
are omitted because the wave strata are already lineage-defined and nested
lineage terms can cause separation or singular covariance estimates.

### Main Patterns

Delta shows the clearest negative deprivation association:

- Cluster size hurdle: OR 0.938, 95% CI 0.925 to 0.951.
- Positive cluster size: count ratio 0.818, 95% CI 0.761 to 0.879.
- Duration hurdle: OR 0.949, 95% CI 0.935 to 0.963.
- Geographic spread hurdle: OR 0.962, 95% CI 0.949 to 0.976.
- Positive geographic spread: count ratio 0.789, 95% CI 0.725 to 0.860.

BA.2 shows a different pattern among clusters exceeding the structural minimum:

- Positive cluster size: count ratio 1.26, 95% CI 1.18 to 1.34.
- Positive geographic spread: count ratio 1.17, 95% CI 1.09 to 1.26.

BA.4 also shows positive positive-count effects:

- Positive cluster size: count ratio 1.75, 95% CI 1.17 to 2.64.
- Positive geographic spread: count ratio 1.67, 95% CI 1.11 to 2.51.

However, BA.4 has only 2,669 clusters, so these estimates should be interpreted
as exploratory and more sensitive to sparse-wave structure.

BA.5 has lower odds of non-singleton clusters and lower positive geographic
spread with higher deprivation. BQ.1 has lower odds of non-singleton clusters.
B.1.177 and Alpha show weaker or outcome-specific associations.

### Interpretation

The per-wave figure argues against a single stable deprivation effect across
the pandemic. Instead, deprivation effects depend on the lineage and epidemic
context. The pooled negative association with cluster size and geographic
spread is heavily compatible with the Delta pattern, whereas later Omicron
subwaves show more heterogeneous associations.

## Supplementary Figure 1: Outcome Distributions

File stem: `supp_fig1_outcome_distributions`

This figure documents why hurdle/ZTNB models are used. The three count outcomes
are heavily concentrated at their structural minima and have long right tails.

In the main cluster table:

- 56.5% of clusters are singletons.
- 63.1% have duration zero days.
- 61.7% are observed in a single datazone.
- The median is 1 for cluster size, 0 for duration, and 1 for distinct
  datazones.
- The maximum observed values are 2,792 sequences, 19 days, and 2,100
  datazones.

The distributional shape makes ordinary log-linear modelling too compressed for
the main question, because it blends the structural-minimum process with the
positive-count tail.

## Supplementary Figure 2: Size-Adjusted Positive Counts

File stem: `supp_fig2_size_adjusted_positive_counts`

This sensitivity analysis asks whether duration and geographic spread effects
remain after additionally adjusting the positive-count models for cluster size.

After size adjustment, overall SIMD deprivation is weakly positively associated
with:

- positive duration: count ratio 1.01, 95% CI 1.01 to 1.02, p = 4.26e-05
- positive geographic spread: count ratio 1.03, 95% CI 1.01 to 1.04,
  p = 8.97e-05

This is important because the unadjusted positive geographic-spread model shows
a negative deprivation association. The size-adjusted sensitivity suggests that
more deprived clusters are generally smaller, but among clusters of comparable
size there may be slightly longer duration or wider geographic spread.

## Supplementary Figure 3: Log-Linear Versus Hurdle/ZTNB

File stem: `supp_fig3_loglinear_vs_hurdle_ztnb`

This figure compares the previous log-linear framing with the main hurdle/ZTNB
formulation.

For SIMD deprivation, the log-linear estimates are:

- cluster size geometric mean ratio 0.992, 95% CI 0.987 to 0.997
- duration geometric mean ratio 0.999, 95% CI 0.994 to 1.004
- geographic spread geometric mean ratio 1.001, 95% CI 0.996 to 1.006

The log-linear model broadly agrees that duration is not meaningfully associated
with overall deprivation, but it substantially attenuates the positive-count
patterns seen in the ZTNB models. This supports the decision to use hurdle/ZTNB
models for the main count outcomes.

## Supplementary Figure 4: SIMD-Domain Cluster Outcomes

File stem: `supp_fig4_simd_domain_cluster_outcomes`

This figure extends the main count models to each SIMD domain. It uses the same
two-part count structure as Figure 1.

The strongest negative domain effects are for housing and crime deprivation.
Housing deprivation is associated with lower odds of non-singleton clusters,
lower odds of positive duration, lower odds of multi-datazone clusters, and
lower positive geographic spread. Crime deprivation is associated with smaller
positive cluster size and lower positive geographic spread.

Access deprivation behaves differently. It is associated with slightly higher
odds of non-singleton clusters and higher positive geographic spread, while
several other domains show negative positive-count associations. This supports
reporting SIMD domains rather than treating overall SIMD as a single mechanism.

## Supplementary Figure 5: SIMD-Domain Quintile Mixing

File stem: `supp_fig5_simd_domain_quintile_mixing`

This one-panel figure isolates the domain-quintile mixing result from Figure 3.
It is useful as a cleaner visual if the manuscript needs to focus specifically
on socioeconomic mixing rather than demographic mixing.

The key contrast remains education/crime versus access/housing:

- education and crime deprivation are associated with more domain-quintile
  excess mixing
- access and housing deprivation are associated with less domain-quintile
  excess mixing

## Supplementary Figure 6: SIMD-Domain Demographic Mixing

File stem: `supp_fig6_simd_domain_demographic_mixing`

This figure separates age, sex, and joint age-sex mixing for each SIMD domain.
It shows that most domains have positive age-mixing associations, most have
negative sex-mixing associations, and several have positive joint age-sex
associations.

Access deprivation again behaves differently, with negative age mixing and no
clear sex or joint age-sex association in the pooled domain model.

## Supplementary Figure 7: Wave-Specific Domain-Demographic Mixing

File stem: `supp_fig7_wave_specific_domain_demographic_mixing`

This heatmap shows whether the domain-demographic mixing associations are stable
across waves. They are not fully stable.

Age mixing shows the most consistent positive domain-deprivation pattern across
B.1.177, Alpha, Delta, BA.1, and BA.5 for many domains. Access is often an
exception and can be negative. Sex mixing is more variable across waves; BA.5
shows several negative sex-mixing estimates, whereas access in BA.5 is positive.
Joint age-sex mixing is generally positive in earlier waves and Delta for
several domains, but again varies by domain and wave.

This figure should be interpreted as an exploratory heterogeneity analysis,
especially for smaller waves.

## Supplementary Figure 8: Observed-Versus-Expected Mixing Matrices

File stem: `supp_fig8_observed_expected_mixing_matrices`

This figure shows the observed-minus-expected pair probability matrices for
SIMD quintile and age band. It gives a descriptive view of which pairs drive the
excess-mixing summaries.

For SIMD quintiles, the largest overall excess is same-quintile 1 pairs, about
+0.28 percentage points above expected. Some cross-quintile pairs, especially
quintile 4 with 5 and quintile 1 with 5, are slightly below expected. The
absolute pair-probability differences are small because they are measured over
the full pair matrix.

For age, the largest overall excess is within age band 20-24, about +0.21
percentage points above expected. Adjacent young-adult pairs, such as 20-24
with 25-29, also show positive excess. Several child-to-adult pairings are
slightly below expected.

The matrices are descriptive rather than inferential. They help make the
regression mixing outcomes more interpretable by showing the pairwise structure
behind the aggregate discordance metrics.

## Suggested Results Paragraph

In the main hurdle/ZTNB models, higher cluster-level SIMD deprivation was not
associated with larger or longer-lasting genomic clusters after adjustment for
lineage, calendar time, local incidence, sequencing intensity, and test
positivity. Instead, deprivation was associated with slightly lower odds of
being non-singleton and smaller positive cluster size, and with substantially
lower positive geographic spread. Duration showed little evidence of an overall
deprivation association. Local epidemic and surveillance conditions were much
more strongly associated with cluster scale: higher incidence, window-level
sequencing proportion, and test positivity were consistently associated with
larger, longer-lasting, and more geographically dispersed clusters. Mixing
analyses showed a more nuanced socioeconomic pattern. Overall deprivation was
not clearly associated with SIMD-quintile mixing, but was associated with
greater age and joint socio-demographic mixing and lower sex mixing. SIMD-domain
models showed that education and crime deprivation were associated with greater
domain-quintile mixing, whereas access and housing deprivation were associated
with lower domain-quintile mixing. Per-wave outcome models indicated that
deprivation effects varied over time, with the strongest negative associations
seen during Delta and more heterogeneous patterns in Omicron subwaves.

## Interpretation For The Part 1 Question

The answer to the original question is therefore mixed:

- Socioeconomic deprivation is associated with some cluster outcomes, but not in
  the simple direction of larger, longer-lasting, or more dispersed clusters.
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

These results support a cautious conclusion: deprivation and local surveillance
conditions are associated with observed genomic cluster structure, but the
associations are outcome-specific, domain-specific, and wave-specific rather
than a single monotonic deprivation effect.
