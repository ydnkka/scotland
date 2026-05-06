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

## Main Figure 1: Main Cluster Outcome Models

File stem: `fig1_main_cluster_outcomes`

This figure shows the main hurdle and zero-truncated negative-binomial results
for cluster size and geographic spread. Each point is an adjusted ratio per
1 SD higher covariate. For hurdle panels, the ratio is an odds ratio. For
positive-count panels, the ratio is a ZTNB count ratio.

### SIMD Deprivation

In the pooled main model:

- Cluster size hurdle: OR 0.971, 95% CI 0.960 to 0.983, p = 1.43e-06.
- Positive cluster size: ZTNB count ratio 0.926, 95% CI 0.869 to 0.987,
  p = 0.018.
- Geographic spread hurdle: OR 1.004, 95% CI 0.992 to 1.016, p = 0.522.
- Positive geographic spread: ZTNB count ratio 0.851, 95% CI 0.792 to 0.915,
  p = 1.24e-05.

The strongest SIMD deprivation result is therefore not larger clusters, but
smaller positive cluster size and lower positive geographic spread after
adjustment.

### Incidence And Surveillance Conditions

Local incidence, local sequencing fraction, window sequencing proportion, and
test positivity are much more consistently associated with larger apparent
clusters:

- Higher local incidence is associated with higher odds of exceeding the
  structural minimum for both modelled count outcomes and with larger positive
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

This figure shows the per-wave overall SIMD deprivation effect on cluster size
and geographic spread. It uses the same hurdle/ZTNB framing as the main count
analysis, but models are fitted separately by wave.

The wave-specific outcome models adjust for the same main covariates: calendar
time spline, local incidence, local sequencing fraction, window sequencing
proportion, test positivity, and SIMD deprivation. Within-wave lineage dummies
are included where estimable, with rank-redundant columns dropped after the
substantive covariates and calendar spline terms are kept.

### Main Patterns

Delta shows the clearest negative deprivation association:

- Cluster size hurdle: OR 0.934, 95% CI 0.921 to 0.947.
- Positive cluster size: count ratio 0.797, 95% CI 0.725 to 0.876.
- Geographic spread hurdle: OR 0.958, 95% CI 0.945 to 0.971.
- Positive geographic spread: count ratio 0.781, 95% CI 0.703 to 0.867.

BA.2 shows a different pattern among clusters exceeding the structural minimum:

- Positive cluster size: count ratio 1.19, 95% CI 1.07 to 1.32.
- Positive geographic spread: count ratio 1.10, 95% CI 0.97 to 1.24.

BA.4 also shows positive positive-count effects:

- Positive cluster size: count ratio 1.68, 95% CI 0.27 to 10.28.
- Positive geographic spread: count ratio 1.80, 95% CI 1.26 to 2.58.

However, BA.4 has only 2,669 clusters, so these estimates should be interpreted
as exploratory and more sensitive to sparse-wave structure.

BA.5 has lower odds of non-singleton clusters and lower positive geographic
spread with higher deprivation. BQ.1 has lower odds of non-singleton clusters
but higher positive cluster size and higher positive geographic spread among
clusters exceeding the structural minimum.
B.1.177 and Alpha show weaker or outcome-specific associations.

### Interpretation

The per-wave figure argues against a single stable deprivation effect across
the pandemic. Instead, deprivation effects depend on the lineage and epidemic
context. The pooled negative association with cluster size and geographic
spread is heavily compatible with the Delta pattern, whereas later Omicron
subwaves show more heterogeneous associations.

## Supplementary Figure 1: Outcome Distributions

File stem: `supp_fig1_outcome_distributions`

This figure now summarises the model outcomes among non-singleton clusters
(`cluster_size > 1`). The first row shows the distributions of cluster size,
duration, and distinct datazones in this non-singleton population. The second
row shows observed-minus-expected excess mixing for age, sex, and SIMD
deprivation quintile composition.

In the main cluster table, 84,067 clusters are non-singletons (43.5% of all
clusters). Among these non-singleton clusters:

- The median cluster size is 3 sequences, and 38.8% have size 2.
- The median duration is 4 days; 15.3% still have duration zero days.
- The median geographic spread is 3 distinct datazones; 12.0% remain within a
  single datazone.
- The maximum observed values are 2,792 sequences, 19 days, and 2,100
  datazones.

The count outcomes retain long right tails even after removing singletons. The
mixing panels show that the excess-mixing outcomes are also heterogeneous, with
age and sex mixing centred slightly above zero and SIMD deprivation mixing
centred slightly below zero.

## Supplementary Figure 2: Size-Adjusted Positive Counts

File stem: `supp_fig2_size_adjusted_positive_counts`

This sensitivity analysis asks whether the geographic spread effect remains
after additionally adjusting the positive-count model for cluster size.

After size adjustment, overall SIMD deprivation is weakly positively associated
with:

- positive geographic spread: count ratio 1.03, 95% CI 1.01 to 1.04,
  p = 4.67e-11

This is important because the unadjusted positive geographic-spread model shows
a negative deprivation association. The size-adjusted sensitivity suggests that
more deprived clusters are generally smaller, but among clusters of comparable
size there may be slightly wider geographic spread.

## Supplementary Figure 3: Log-Linear Versus Hurdle/ZTNB

File stem: `supp_fig3_loglinear_vs_hurdle_ztnb`

This figure compares the previous log-linear framing with the main hurdle/ZTNB
formulation.

For SIMD deprivation, the log-linear estimates are:

- cluster size geometric mean ratio 0.992, 95% CI 0.987 to 0.997
- geographic spread geometric mean ratio 1.001, 95% CI 0.996 to 1.006

The log-linear model substantially attenuates the positive-count patterns seen
in the ZTNB models, especially for cluster size and geographic spread. This
supports the decision to use hurdle/ZTNB models for the main count outcomes.

## Supplementary Figure 4: SIMD-Domain Cluster Outcomes

File stem: `supp_fig4_simd_domain_cluster_outcomes`

This figure extends the main count models to each SIMD domain. It uses the same
two-part count structure as Figure 1.

The strongest negative domain effects are for housing and crime deprivation.
Housing deprivation is associated with lower odds of non-singleton clusters,
lower odds of multi-datazone clusters, and lower positive geographic spread.
Crime deprivation is associated with smaller positive cluster size and lower
positive geographic spread.

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

## Suggested Results Paragraph

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
