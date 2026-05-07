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

**File:** `fig1_vaccinated_cases_over_time`

**What it shows.** Two-panel weekly time series covering the full study period.
Panel A shows the proportion of sequenced cases who were vaccinated, stratified
by JCVI rollout age group (00–14, 15–19, 20–29, 30–39, 40–49, 50–64, 65–74,
75+). Panel B shows the same proportion stratified by overall SIMD deprivation
quintile (Q1 most deprived to Q5 least deprived). Lines are masked where the
weekly stratum count falls below 20. Both panels share the same calendar axis.

**Key visual patterns.** Panel A shows a clear age-gradient in timing: the
oldest groups (65–74, 75+) reach high vaccinated-case proportions first,
consistent with the JCVI priority rollout sequence. Younger groups rise later
and more gradually; the 00–14 group never reaches adult proportions within the
study window. Panel B shows a deprivation gradient during the primary-course
rollout phase: the least-deprived quintile (Q5) shows elevated proportions
earlier than the most-deprived quintile (Q1). Lines converge across quintiles
by the booster phase and remain broadly similar through the Omicron subwaves.

**Suggested results paragraph.**
> Figure 1 shows the weekly proportion of sequenced cases who were vaccinated
> throughout the study period. Older JCVI priority groups accumulated
> breakthrough cases first: the 65–74 and 75+ groups reached high proportions
> well before younger groups, and the 00–14 group never reached adult-level
> proportions within the study window (Panel A). A deprivation gradient was
> visible during the primary rollout: the least-deprived quintile (Q5) showed
> elevated vaccinated-case proportions earlier than the most-deprived quintile
> (Q1), consistent with differential uptake and earlier rollout access in less-
> deprived areas (Panel B). By the booster phase, proportions across quintiles
> had largely converged.

---

## Figure 2: Cluster Vaccination Profile And Mean Vaccination By Wave And Size

**File:** `fig2_cluster_vaccination_by_wave_and_category`

**What it shows.** Panel A shows stacked proportional bars for non-singleton
clusters, with each bar representing one epidemic wave and each stack segment
showing the fraction with none vaccinated, mixed, or all vaccinated. Panel B
shows a dot plot of mean cluster proportion vaccinated by wave and cluster size
category (small/moderate, large, very large), allowing size-category comparisons
within and across waves.

**Key visual patterns.** Panel A tracks the rollout trajectory across waves:
B.1.177 is almost entirely none-vaccinated (89.5%); Delta is the peak
mixed-coverage wave (65.3% mixed); from BA.2 onwards all-vaccinated clusters
dominate, reaching 80.6% in XBB. Panel B shows that within waves, very large
clusters tend to have slightly lower mean vaccination than small/moderate
clusters during BA.1–BA.2, consistent with larger clusters having greater
within-cluster heterogeneity in vaccination status. By later subwaves
differences across size categories shrink as coverage approaches saturation.

**Suggested results paragraph.**
> Cluster vaccination profiles tracked the national rollout trajectory across
> waves (Panel A). None-vaccinated clusters fell from 89.5% of non-singleton
> clusters in B.1.177 to 0% in XBB, while all-vaccinated clusters rose from
> 0.5% to 80.6% over the same period. Mixed-vaccination clusters peaked during
> Delta (65.3%) and Alpha (56.2%), when the rollout was spreading
> heterogeneously across the adult population. Within waves, very large clusters
> showed slightly lower mean vaccination proportions than small/moderate clusters
> during BA.1–BA.2, consistent with within-cluster heterogeneity in vaccination
> status being higher in larger clusters (Panel B). By late Omicron subwaves,
> differences across size categories were negligible.

---

## Figure 3: Vaccination-Status Mixing By Wave

**File:** `fig3_vaccination_mixing_by_wave`

**What it shows.** Stacked proportional bars for non-singleton clusters, one bar
per epidemic wave, showing the fraction of clusters classified as homogeneous
(more same-vaccination-status pairs than expected), at baseline, or mixed (more
cross-vaccination-status pairs than expected), using observed-minus-expected
pairwise vaccination-status discordance within the same lineage and calendar
window. Cluster counts per wave are annotated above each bar.

**Key visual patterns.** B.1.177 is predominantly at baseline (53.9%), as near-
universal unvaccinated status leaves little scope for excess mixing. Alpha shows
the first substantial mixed-cluster fraction (39.9%). Delta is the peak mixed-
mixing wave (49.1% mixed, 46.9% homogeneous) — the only wave where mixed and
homogeneous categories are nearly equal. From BA.2 onwards, homogeneous clusters
dominate with increasing fractions (57–66%), reflecting near-universal
vaccination. BQ.1 and XBB shift back toward baseline as vaccination status
loses discriminating power.

**Suggested results paragraph.**
> Vaccination-status mixing within clusters followed the rollout trajectory.
> Before vaccination, clusters were predominantly at baseline (53.9% in
> B.1.177). During Alpha, 39.9% of non-singleton clusters showed mixed
> vaccination-status composition. Mixed-mixing clusters were most prevalent
> during Delta (49.1%), coinciding with the period of greatest heterogeneity in
> population vaccination coverage. From BA.2 onwards, homogeneous clusters
> dominated (57–66%), reflecting near-universal vaccination among the sequenced
> infected population. BQ.1 and XBB showed a partial return toward baseline
> fractions as vaccination status became uninformative at near-universal coverage.

---

## Figure 4: Demographic Mixing Categories By Wave

**File:** `fig4_demographic_mixing_by_wave`

**What it shows.** A 2×2 figure of stacked proportional bars for non-singleton
clusters, one panel per demographic mixing dimension: SIMD deprivation quintile
mixing (A), age-band mixing (B), sex mixing (C), and joint SIMD-age-sex profile
mixing (D). Each bar represents one epidemic wave and each stack shows the
fraction of clusters classified as less mix, baseline, or more mix using the
±0.01 threshold applied to the Part 1 excess-discordance scores.

**Key visual patterns.** Panel A (SIMD mixing) is the only dimension with a
"less mix" majority (54.4% overall): Alpha stands out with 67.9% less mix,
consistent with early-rollout transmission being geographically concentrated in
areas of similar deprivation. The "more mix" fraction increases in later Omicron
subwaves. Panel B (age mixing) shows a "more mix" majority (54.2% overall) that
is broadly stable across waves. Panel C (sex mixing) similarly shows a "more
mix" majority (56.6% overall) across all waves. Panel D (joint profile mixing)
is dominated by the baseline category (75.7% overall), with deviations most
prominent in BQ.1 (45.8% more mix) and BA.4 (40.3% more mix).

**Suggested results paragraph.**
> SIMD deprivation mixing showed the highest "less mix" fraction of any
> demographic dimension (54.4% overall), with the Alpha wave having the most
> spatially homogeneous SIMD composition (67.9% less mix), consistent with
> early rollout transmission being geographically concentrated (Panel A). The
> "more mix" fraction increased in later Omicron subwaves (BA.4–XBB: 43–47%).
> Age mixing was predominantly "more mix" (54.2% overall), broadly stable across
> waves (Panel B). Sex mixing similarly showed a "more mix" majority (56.6%
> overall), stable across waves (Panel C). Joint SIMD-age-sex profile mixing was
> dominated by the baseline category (75.7% overall), with deviations most
> notable in BQ.1 and BA.4, which showed the highest "more mix" fractions
> (45.8% and 40.3% respectively; Panel D).

---

## Figure 5: Geographic Dispersion Category By Wave

**File:** `fig5_geographic_dispersion_by_wave`

**What it shows.** Panel A shows stacked proportional bars of geographic
dispersion category (low/moderate, large, very large) for non-singleton clusters
by epidemic wave. Panel B shows a dot plot of mean cluster proportion vaccinated
by wave and dispersion category, allowing comparison of vaccination levels across
dispersion categories within each wave.

**Key visual patterns.** Panel A shows that low/moderate dispersion clusters
dominate in every wave (89.8% overall), with large and very large dispersion
clusters most prevalent during Delta and B.1.177. Omicron subwaves show a higher
fraction of low/moderate dispersion clusters. In Panel B, very large dispersion
clusters show higher mean vaccination proportions than low/moderate clusters
in later waves (from BA.2 onwards), while the pattern reverses or is flat in
earlier waves.

**Suggested results paragraph.**
> Low/moderate geographic dispersion clusters were the overwhelming majority in
> all waves (89.8% overall; Panel A). Large and very large dispersion clusters
> were most prevalent during Delta and B.1.177, consistent with the scale of
> those outbreaks, while later Omicron subwaves showed a higher fraction of
> low/moderate dispersion clusters, possibly reflecting more localised
> transmission chains. In the booster era (BA.2 onwards), very large dispersion
> clusters showed higher mean vaccination proportions than low/moderate dispersion
> clusters (Panel B), consistent with geographically dispersed late-pandemic
> transmission events being driven primarily by fully vaccinated cluster members.

---

## Figure 6: Booster Coverage And Dose Recency By Wave And SIMD

**File:** `fig6_dose_recency_by_simd`

**What it shows.** A double-panel heatmap. Panel A shows mean booster coverage
(%) among vaccinated cluster members for each wave × SIMD quintile cell. Panel B
shows mean days since last prior vaccination dose for each wave × SIMD quintile
cell. Cells where mean cluster vaccination proportion was below 5% are masked
(shown grey), representing waves where the vaccinated cluster population was too
small for meaningful booster estimates.

**Key visual patterns.** Panel A shows booster coverage rising sharply from
negligible levels before BA.1 to 72–85% across all quintiles in the BA.4–XBB
period. Modest SIMD gradients are visible within waves: less-deprived quintiles
(Q4–Q5) tend to show slightly higher booster coverage than Q1–Q2, but absolute
differences are small. Panel B shows days since last dose increasing steadily
from roughly 10–30 days in B.1.177/Alpha to over 200 days in XBB. SIMD gradients
within late waves are largely flat, and the BA.4–XBB period shows broadly similar
values across quintiles (200–230 days).

**Suggested results paragraph.**
> Booster coverage among vaccinated cluster members was negligible before BA.1
> and rose to 72–85% during the BA.4–XBB period (Panel A). Modest SIMD gradients
> were present within waves, with less-deprived quintiles showing slightly higher
> booster coverage than most-deprived quintiles, but absolute differences were
> small across all waves. Mean days since last vaccination dose increased from
> approximately 10–30 days in B.1.177/Alpha to over 200 days in XBB (Panel B),
> with SIMD gradients within late waves largely flat, indicating that the booster
> campaign timing was broadly similar across deprivation quintiles by the late
> Omicron subwave period. Median days since last dose plateaued at roughly
> 210–220 days in the late Omicron subwaves.

---

## Supplementary Figure 1: Weekly Vaccination-Status Mixing Evolution

**File:** `supp_fig1_weekly_mixing_evolution`

**What it shows.** A stacked area chart of 4-week rolling-mean fractions of
homogeneous, baseline, and mixed vaccination-status clusters over the full
calendar period (mid-2020 to early 2023). Each area represents one mixing
category; the three sum to 1.0 at each week. The chart provides a continuous
temporal view of the wave-level patterns summarised in Figure 3.

**Key visual patterns.** The figure shows a clear three-phase narrative: a
baseline-dominated period before the rollout; a mixed-cluster peak coinciding
with Delta (2021 Q3–Q4); and a homogeneous-dominant period from BA.2 onwards.
A brief return toward baseline fractions in the XBB period is visible at the
rightmost end of the chart. The continuity of the area chart makes the gradual
transitions between phases clearer than the wave-by-wave bar chart in Figure 3.

**Suggested results paragraph.**
> Supplementary Figure 1 shows the continuous temporal evolution of vaccination-
> status mixing categories over the full study period. The baseline-dominated
> pre-rollout phase, the rise of mixed-cluster fractions during the primary
> rollout and Alpha wave, the Delta peak of mixed-mixing clusters, and the
> subsequent homogeneous-dominant phase from BA.2 onwards are all clearly visible
> as a continuous narrative. A brief return toward baseline fractions appears at
> the end of the study period (XBB), consistent with vaccination status losing
> discriminating power at near-universal coverage.

---

## Supplementary Figure 2: Deprivation Gradient In Dose Metrics By SIMD Domain

**File:** `supp_fig2_domain_dose_gradient`

**What it shows.** A double-panel heatmap showing, for each SIMD domain and
wave, the Q5 − Q1 gradient in booster coverage (Panel A) and mean days since
last dose (Panel B) among vaccinated cluster members. Positive values in Panel A
mean the least-deprived quintile had higher booster coverage than the most
deprived; positive values in Panel B mean the least-deprived quintile had longer
time since last dose.

**Key visual patterns.** Panel A shows mostly positive gradients (less-deprived
Q5 had higher booster coverage than most-deprived Q1), with the largest
gradients during BA.1 and BA.2. Access deprivation stands out with negative
gradients in several waves, indicating that the most-access-deprived clusters
had higher booster coverage than the least-access-deprived in those waves. Panel
B shows more variable and larger-magnitude gradients by sign. For most domains
in BA.1, less-deprived clusters had fewer days since last dose (negative
gradient), consistent with the booster campaign reaching less-deprived areas
first. By BA.4–XBB the gradient largely flattens across all domains.

**Suggested results paragraph.**
> Supplementary Figure 2 shows the Q5 − Q1 deprivation gradient in booster
> coverage and dose recency, stratified by SIMD domain and wave. Most domain ×
> wave combinations showed a modest positive booster-coverage gradient, with
> less-deprived clusters having higher booster coverage, and the largest
> gradients in BA.1 and BA.2 (Panel A). Access deprivation was a consistent
> exception, showing negative gradients in several waves — meaning most-access-
> deprived clusters had higher booster coverage than least-access-deprived
> clusters. The days-since-dose gradient (Panel B) was more variable in sign and
> magnitude; in BA.1, less-deprived clusters tended to have fewer days since
> their last dose, consistent with the booster campaign reaching those areas
> first. By BA.4–XBB both gradients were largely flat across all domains.

---

## Supplementary Figure 3: Cross-Category Heatmap — "More Mix"

**File:** `supp_fig3_cross_category_heatmap`

**What it shows.** A 2×2 heatmap showing, for each demographic mixing dimension
(SIMD quintile mixing, age mixing, sex mixing, joint SIMD-age-sex profile
mixing), the fraction of non-singleton clusters classified as "more mix" in each
SIMD quintile × cluster size category stratum. Grey "—" cells indicate strata
where cluster counts were too sparse for reliable estimates.

**Key visual patterns.** SIMD mixing (Panel A): small/moderate clusters show
moderate "more mix" fractions (15–55% across quintiles), but large and very
large clusters are near zero, as these are overwhelmingly "less mix" (see
Supplementary Figure 4). Age mixing (Panel B): small/moderate clusters show
consistently high "more mix" fractions (57–65% across quintiles), while large
clusters drop to 7–22%, consistent with larger clusters sampling a narrower age
range. Sex mixing (Panel C): more uniform across size and quintile, with
small/moderate clusters at 57–59% and large clusters at 21–53%. Joint profile
mixing (Panel D): "more mix" fractions are low across all strata (1–11% for
small/moderate), reflecting that joint profile mixing is predominantly at
baseline.

**Suggested results paragraph.**
> Supplementary Figure 3 shows the fraction of clusters classified as "more mix"
> by demographic dimension across SIMD quintile and cluster size categories.
> SIMD-quintile "more mix" fractions were moderate in small/moderate clusters
> (15–55% across quintiles) but near zero in large and very large clusters,
> which are overwhelmingly SIMD-homogeneous (Panel A). Age "more mix" fractions
> were consistently high in small/moderate clusters (57–65%) and fell markedly
> in large clusters (7–22%), consistent with larger clusters sampling narrower
> age ranges from geographically concentrated outbreaks (Panel B). Sex mixing
> "more mix" fractions were broadly uniform across size and quintile (Panel C).
> Joint profile "more mix" fractions were low across all strata, reflecting the
> near-universal baseline classification for this high-dimensional category
> (Panel D).

---

## Supplementary Figure 4: Cross-Category Heatmap — "Less Mix"

**File:** `supp_fig4_cross_category_less_mix_heatmap`

**What it shows.** The complement to Supplementary Figure 3, showing the
fraction of non-singleton clusters classified as "less mix" (more homogeneous
than expected within lineage-calendar stratum) in each SIMD quintile × cluster
size category stratum for each demographic mixing dimension.

**Key visual patterns.** SIMD mixing (Panel A): large clusters are
overwhelmingly "less mix" across all quintiles (82–100%), the strongest
size-gradient result across all mixing dimensions. Even small/moderate clusters
show 16–83% "less mix" by SIMD. Age mixing (Panel B): large clusters show
49–89% "less mix", with particularly high age homogeneity in large, least-
deprived clusters (Q5: 89%). Sex mixing (Panel C): "less mix" fractions are more
uniform across size and quintile (21–59% for large clusters). Joint profile
mixing (Panel D): large clusters show striking "less mix" fractions in some
quintiles (Q5: 89%, Q1: 79%), the strongest compositional clustering signal
across all dimensions.

**Suggested results paragraph.**
> Cross-category analyses showed that large clusters were strongly homogeneous by
> SIMD deprivation quintile (82–100% "less mix" regardless of the deprivation
> stratum), and that this size-dependent SIMD homogeneity extended to the joint
> SIMD-age-sex profile in the largest transmission events (Q5 large clusters:
> 89% less mix; Q1: 79%; Panel D). Age homogeneity was also high in large,
> least-deprived clusters (89% less mix; Panel B), potentially reflecting
> workplace or community outbreak settings with narrow age ranges. Sex mixing
> showed a more uniform "less mix" distribution across size and quintile
> categories (Panel C). These patterns are consistent with large clusters
> representing spatially concentrated transmission events rather than diffuse
> community spread spanning multiple socioeconomic and demographic groups.

---

## Overall Suggested Results Paragraph

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
