# Part 2 Results And Figure Guide

Generated for the Part 2 vaccination characterisation analysis on 7 May 2026.

## Analysis Question

Part 2 asks how vaccination status, dose profile, and vaccination-status mixing
characterised Scottish SARS-CoV-2 genomic clusters across epidemic waves, cluster
size and geographic-dispersion categories, and SIMD deprivation groups.

The specific descriptive outcomes are:

- weekly vaccinated-case proportions by JCVI rollout age group and SIMD
  deprivation quintile
- cluster vaccination profile (none / mixed / all vaccinated) by wave
- vaccination-status mixing (homogeneous / baseline / mixed) by wave
- dose history (dose number, booster coverage, days since last dose) by wave
  and SIMD quintile
- demographic and geographic mixing categories by wave and cluster size

The manuscript figures are generated from `part2/tables` CSV files and
`part2/cache` parquet files by `part2/manuscript/make_figures.py`.

## Analysis Population And Descriptive Frame

The analysis uses good-QC sequences at Leiden resolution 0.3. This produces
193,160 inferred genomic clusters, of which 84,080 are non-singleton (cluster
size > 1). All category thresholds are estimated within the non-singleton
population.

Cluster size categories (90th/99th non-singleton percentile thresholds): small/
moderate (< 13 sequences), large (13–73), very large (≥ 74). Geographic
dispersion categories: low/moderate (< 11 datazones), large (11–60), very large
(≥ 61). Overall SIMD quintile 1 is most deprived.

Vaccination-status mixing categories use observed-minus-expected pairwise
vaccination-status discordance within the same `window_id × pango_lineage`
stratum: homogeneous (< −0.01), baseline (±0.01), mixed (> +0.01). Demographic
mixing categories for SIMD, age, sex, and joint profile use the same ±0.01
threshold applied to the corresponding Part 1 excess-discordance scores.

## Overall Descriptive Summary

Across 193,160 clusters, mean cluster proportion vaccinated was 68.9%.
Vaccination-status mixing was homogeneous in 53.3% of non-singleton clusters,
at baseline in 7.1%, and mixed in 39.6% (singletons excluded). Median days
since last vaccination among vaccinated cluster members was 120 days.

Cluster vaccination profiles followed the national rollout trajectory
faithfully. None-vaccinated clusters fell from 89.5% in B.1.177 to 0% in XBB;
all-vaccinated clusters rose from 0.5% to 80.6% over the same period.
Mixed-vaccination clusters peaked at 65.3% during Delta, when the rollout was
spreading heterogeneously across the adult population.

Demographic mixing among non-singleton clusters was predominantly "more mix"
by age (54.2%) and sex (56.6%), but predominantly "less mix" by SIMD
deprivation quintile (54.4%). Joint profile mixing was mostly at baseline
(75.7%), reflecting the high specificity of the combined SIMD-age-sex category.

---

## Figure 1: Vaccinated-Case Proportions Over Time

File stem: `fig1_vaccinated_cases_over_time`

This figure tracks the weekly proportion of sequenced cases that were vaccinated
throughout the study period. It is the entry point into the Part 2 narrative
because it shows how the sequenced infected population changed character as the
rollout progressed.

### Panel A — by JCVI rollout age group

Older age groups (65–74, 75+) reached high breakthrough-case proportions first,
consistent with the JCVI priority sequencing from oldest to youngest. Lines are
masked where the weekly stratum count falls below 20. The younger groups (00–14,
15–19) started later and rose more gradually, with the 00–14 group never
reaching the same proportions as adults because most of that group was not
eligible for primary vaccination within the study window.

### Panel B — by overall SIMD deprivation quintile

There is a visible deprivation gradient in the pace at which breakthrough cases
accumulated during the primary-course rollout phase: the least-deprived quintile
(Q5) shows elevated proportions earlier than the most-deprived (Q1). This is
consistent with differential uptake and access during the early rollout. Lines
converge by the booster phase and remain broadly similar through the Omicron
subwaves.

### Key interpretive note

The y-axis shows vaccinated cases as a proportion of sequenced cases, not as a
proportion of all infected cases or of the vaccinated population. Rising
proportions reflect a combination of increasing vaccine coverage, changing
testing behaviour, and sequencing patterns — not vaccine effectiveness.

---

## Figure 2: Cluster Vaccination Profile And Mean Vaccination By Wave And Size

File stem: `fig2_cluster_vaccination_by_wave_and_category`

### Panel A — cluster vaccination profile by wave

Stacked proportional bars for non-singleton clusters. The wave progression
shows the rollout trajectory in the cluster population:

- B.1.177: 89.5% none vaccinated — effectively a pre-rollout baseline.
- Alpha: 41.5% none / 56.2% mixed — rollout reaching younger adults.
- Delta: 10.3% none / 65.3% mixed / 24.4% all vaccinated — peak mixed-coverage
  period.
- BA.1: 4.4% none / 49.9% mixed / 45.7% all vaccinated — transition wave.
- BA.2 onwards: all-vaccinated clusters are the majority. By BQ.1, 78.9% of
  non-singleton clusters were entirely vaccinated.
- XBB: 80.6% all vaccinated, 0% none vaccinated.

### Panel B — mean cluster proportion vaccinated by wave and cluster size

Among all waves, mean vaccination within a cluster increases monotonically from
B.1.177 to XBB. By wave, very large clusters tend to show slightly lower mean
vaccination proportions than small/moderate clusters in some waves (particularly
BA.1–BA.2), consistent with larger clusters having greater within-cluster
diversity of vaccination status. The dot-plot format allows the three size
categories to be compared at each wave without stacking.

---

## Figure 3: Vaccination-Status Mixing By Wave

File stem: `fig3_vaccination_mixing_by_wave`

This figure is the vaccination analogue of the excess-mixing figures in Part 1.
It shows, for each wave, what fraction of non-singleton clusters were more
homogeneous, at baseline, or more mixed by vaccination status than expected from
cases in the same lineage and calendar window.

### Key wave patterns

- **B.1.177**: 53.9% at baseline — before rollout, vaccination status is nearly
  uniform (almost no one vaccinated), so there is little scope for excess mixing.
- **Alpha**: 57.5% homogeneous, 39.9% mixed — the first wave with substantial
  heterogeneity in vaccination status. Nearly half of non-singleton clusters
  showed greater vaccination-status mixing than expected, consistent with rapid
  rollout creating diverse exposure groups during Alpha transmission.
- **Delta**: 46.9% homogeneous, 49.1% mixed — the peak mixed-mixing wave.
  Mixed vaccination-status clusters were nearly as common as homogeneous ones.
- **BA.1–BA.5**: homogeneous clusters increasingly dominate (57–66%) as
  vaccination becomes near-universal in the sequenced infected population.
- **BQ.1–XBB**: shift back toward baseline as vaccine status loses discriminating
  power in a nearly fully vaccinated population.

The n labels above each bar show the total non-singleton cluster count per wave.

---

## Figure 4: Demographic Mixing Categories By Wave

File stem: `fig4_demographic_mixing_by_wave`

This 2×2 figure shows stacked proportional bars for SIMD deprivation mixing (A),
age mixing (B), sex mixing (C), and joint SIMD-age-sex profile mixing (D),
by epidemic wave. Each panel uses the three-category scheme (less mix / baseline
/ more mix) with the ±0.01 threshold applied to the Part 1 excess-discordance
scores.

### SIMD deprivation mixing (Panel A)

SIMD mixing shows the highest fraction of "less mix" clusters across all waves:
54.4% overall. Alpha has the highest "less mix" fraction at 67.9%, consistent
with early-rollout transmission being geographically concentrated in areas of
similar deprivation. The "more mix" fraction increases in later Omicron subwaves
(BA.4–XBB: 43–47%), perhaps reflecting broader geographic spread of XBB.

### Age mixing (Panel B)

Age mixing shows a slight majority of "more mix" clusters (54.2% overall) in
most waves. B.1.177 (55.8%) and Delta (53.5%) show the highest "more mix"
fractions. XBB has the lowest (50.4%). This pattern is consistent with Part 1
finding that higher deprivation was associated with greater age mixing.

### Sex mixing (Panel C)

Sex mixing also shows a "more mix" majority (56.6% overall). The pattern is
broadly stable across waves, with the exception of XBB and Other, which have
higher baseline fractions due to their smaller cluster counts. The finding
that most clusters are more sex-mixed than expected is consistent with
community transmission involving both sexes rather than sex-segregated
transmission.

### Joint profile mixing (Panel D)

Profile mixing is dominated by baseline (75.7% overall), with "more mix"
accounting for only 9.4% and "less mix" 14.9%. The BQ.1 wave shows the highest
"more mix" fraction (45.8%), and BA.4 is also high at 40.3%. The large baseline
fraction reflects that the joint SIMD-age-sex combination is sufficiently
specific that most clusters do not deviate substantially from within-stratum
expectations.

---

## Figure 5: Geographic Dispersion Category By Wave

File stem: `fig5_geographic_dispersion_by_wave`

### Panel A — geographic dispersion category proportions by wave

Stacked proportional bars showing the fraction of clusters in each dispersion
category by wave. Low/moderate dispersion clusters are the overwhelming majority
in all waves (89.8% overall). Large and very large dispersion clusters are most
prevalent in Delta and B.1.177, consistent with the scale of those outbreaks.
Later Omicron subwaves show a higher fraction of low/moderate dispersion
clusters, possibly reflecting more localised transmission chains.

### Panel B — mean cluster proportion vaccinated by wave and dispersion category

Very large dispersion clusters show higher mean vaccination in later waves
compared with low/moderate dispersion clusters. This likely reflects geographic
spread being driven primarily by vaccinated-only clusters in the booster era,
when coverage was near-universal.

---

## Figure 6: Booster Coverage And Dose Recency By Wave And SIMD

File stem: `fig6_dose_recency_by_simd`

This double-panel heatmap summarises the wave × SIMD quintile structure of two
dose-history variables among non-singleton clusters.

### Panel A — mean booster coverage among vaccinated members (%)

Booster coverage was negligible before BA.1 (B.1.177 and Alpha: 0%). From BA.1
onwards it rose rapidly across all quintiles. In the BA.4–XBB period, booster
coverage among vaccinated cluster members approached 72–85% across most quintiles.
There are modest SIMD gradients within waves: less-deprived quintiles (Q4–Q5)
tend to show slightly higher booster coverage than Q1–Q2, but the absolute
differences are small.

Cells are masked (shown grey) where mean cluster vaccination proportion was
below 5%, representing waves where the vaccinated cluster population was too
small to compute meaningful booster estimates.

### Panel B — mean days since last prior vaccination dose

Days since last dose increase across waves from roughly 10–30 days in B.1.177/
Alpha to over 200 days in XBB. The gradients across SIMD quintiles are modest
within most waves. The BA.4–XBB period shows high and broadly similar values
across all quintiles (200–230 days), consistent with booster campaigns having
been completed many months prior and no major new rollout occurring during these
waves.

---

## Supplementary Figure 1: Weekly Vaccination-Status Mixing Evolution

File stem: `supp_fig1_weekly_mixing_evolution`

This figure shows the 4-week rolling-mean fractions of homogeneous, baseline,
and mixed vaccination-status clusters as a stacked area chart over the full
calendar period. It provides a continuous temporal view of the wave-level
patterns summarised in Figure 3.

The transition from baseline dominance (pre-rollout) through the mixed peak
(Delta, 2021 Q3–Q4) to homogeneous dominance (BA.2 onwards) is clearly visible
as a single continuous narrative. The chart ends with a brief return toward
baseline fractions in the XBB period, consistent with vaccination status
becoming uninformative at near-universal coverage.

---

## Supplementary Figure 2: Deprivation Gradient In Dose Metrics By SIMD Domain

File stem: `supp_fig2_domain_dose_gradient`

This double-panel heatmap shows, for each SIMD domain and wave, the Q5 − Q1
gradient in booster coverage (panel A) and mean days since last dose (panel B).
Positive values in panel A mean the least-deprived quintile had higher booster
coverage than the most deprived; positive values in panel B mean the
least-deprived quintile had longer time since last dose.

### Booster coverage gradient (Panel A)

Most domain × wave combinations show a modest positive gradient (less-deprived
Q5 had higher booster coverage than most-deprived Q1), with the largest
gradients during BA.1 and BA.2. Access deprivation stands out with negative
gradients in several waves, meaning most-access-deprived clusters had higher
booster coverage than least-access-deprived clusters in those waves. This is
consistent with the access domain capturing a different social geography than
other deprivation dimensions.

### Days since last dose gradient (Panel B)

The days-since-dose gradient is more variable in sign and larger in magnitude
than the booster gradient. For most domains in BA.1, less-deprived clusters
had fewer days since their last dose (negative gradient), consistent with the
booster campaign reaching less-deprived areas first. By BA.4–XBB the gradient
largely flattens. Access deprivation again shows a different pattern from the
other domains.

---

## Supplementary Figure 3: Cross-Category Heatmap — "More Mix"

File stem: `supp_fig3_cross_category_heatmap`

This 2×2 heatmap shows, for each demographic mixing dimension (SIMD, age, sex,
joint profile), the fraction of non-singleton clusters classified as "more mix"
in each SIMD quintile × cluster size stratum.

### Key patterns

**SIMD mixing (Panel A):** Small/moderate clusters show moderate "more mix"
fractions (15–55% across quintiles), peaking at Q3. Large and very large
clusters are almost entirely "less mix" (see supp_fig4), so "more mix" fractions
are near zero for large and very large clusters. The grey "—" cells indicate
strata where few or no clusters were classified due to sparse cell counts.

**Age mixing (Panel B):** Small/moderate clusters show consistently high "more
mix" fractions (57–65% across quintiles), while large clusters show 7–22%.
This supports the interpretation that large clusters are more geographically
concentrated and therefore sample a narrower age range than chance would predict.

**Sex mixing (Panel C):** The pattern is similar across size categories and
quintiles: most small/moderate clusters show "more mix" (57–59%). Large clusters
show moderate "more mix" fractions (21–53%). Sex mixing is less dependent on
cluster size than SIMD or age mixing.

**Joint profile mixing (Panel D):** "More mix" fractions are low across all
strata (1–11% for small/moderate, slightly higher for some large-cluster
quintiles). This reflects that joint profile mixing is predominantly at
baseline, as shown in Figure 4D.

---

## Supplementary Figure 4: Cross-Category Heatmap — "Less Mix"

File stem: `supp_fig4_cross_category_less_mix_heatmap`

The complement to supp_fig3, showing the fraction of non-singleton clusters
classified as "less mix" (more homogeneous than expected).

### Key patterns

**SIMD mixing (Panel A):** Large clusters are overwhelmingly "less mix" across
all quintiles (82–100%), confirming that larger genomic clusters are highly
concentrated by SIMD deprivation quintile. Even among small/moderate clusters,
16–83% are "less mix" by SIMD. This is the strongest size-gradient result across
all mixing dimensions.

**Age mixing (Panel B):** Less-mix fractions for small/moderate clusters are
32–39% across quintiles. For large clusters, 49–89% are "less mix", particularly
in Q5 (least deprived, 89%). Age homogeneity in large, least-deprived clusters
may reflect spatial concentration of outbreaks in workplace or community settings
with narrow age ranges.

**Sex mixing (Panel C):** Less-mix fractions are relatively uniform across size
and quintile (21–59% for large clusters). Sex mixing is less differentiated by
cluster size than SIMD or age mixing.

**Joint profile mixing (Panel D):** Large clusters show striking "less mix"
fractions in some quintiles: Q5 (89%), Q1 (79%). This means that in large
clusters — which are both geographically concentrated (from supp_fig2) and
SIMD-homogeneous — the full joint SIMD-age-sex profile is also homogeneous
relative to expectation. This is the strongest signal for compositional
clustering in the largest transmission events.

---

## Suggested Results Paragraph

The proportion of sequenced cases who were vaccinated rose from near zero in
mid-2021 to above 90% in the XBB period, with older JCVI priority groups and
less-deprived areas accumulating breakthrough cases earlier than younger and more
deprived groups. Cluster vaccination profiles tracked the national rollout:
none-vaccinated clusters fell from 89.5% of non-singleton clusters in B.1.177
to 0% in XBB, while all-vaccinated clusters rose from 0.5% to 80.6% over the
same period. Mixed vaccination-status clusters were most prevalent during Delta
(65.3%) and Alpha (56.2%), when the rollout was spreading heterogeneously across
the adult population. Vaccination-status mixing within clusters followed the
same trajectory: mixed-mixing clusters were most common during Delta and BA.1,
shifting to homogeneous dominance from BA.2 onwards as vaccination became the
near-universal background condition. Booster coverage among vaccinated cluster
members reached 72–85% in the BA.4–XBB period, with modest but consistent SIMD
gradients. Median days since last dose plateaued at roughly 210–220 days in the
late Omicron subwaves. Cross-category analyses showed that large clusters were
strongly homogeneous by SIMD deprivation quintile (82–100% "less mix") regardless
of the deprivation stratum, and that this size-dependent SIMD homogeneity
extended to the joint SIMD-age-sex profile in the largest transmission events.

## Interpretation In Context Of Part 1

The vaccination characterisation in Part 2 reinforces the Part 1 conclusion
without overturning it. Part 1 showed that deprivation was not associated with
larger or more dispersed clusters in a simple positive direction, and that
epidemic and surveillance context dominated apparent cluster scale. Part 2 shows
that the vaccination profiles of those clusters were primarily determined by
wave timing and rollout calendar, not by cluster-structural differences related
to deprivation.

The key addition from Part 2 is that large clusters were markedly more SIMD-
homogeneous than small/moderate clusters. This is consistent with large clusters
representing transmission events concentrated in specific geographic localities
with narrow socioeconomic ranges, rather than diffuse community transmission
spanning deprivation quintiles. It also explains why the Part 1 mixing models
found cluster size to be a strong predictor of SIMD excess mixing: large clusters
generate more "less mix" outcomes by SIMD almost regardless of their deprivation
quintile.

The deprivation gradients in vaccination coverage and dose recency were
measurable but modest. They did not translate into the large deprivation
differentials in cluster structure that a simple vaccine-protection hypothesis
would predict, which is consistent with the Part 1 finding that surveillance and
epidemic context — not deprivation per se — are the primary structural
determinants of observed genomic cluster scale.
