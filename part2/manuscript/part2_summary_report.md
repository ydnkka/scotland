# Scottish SARS-CoV-2 Genomic Clustering: Part 2 Summary Report

Prepared: 7 May 2026

## 1. Introduction

This Part 2 analysis characterises the vaccination profiles of SARS-CoV-2
genomic clusters in Scotland over the full study period (B.1.177 through XBB).
The central question is how vaccination status, dose history, and
vaccination-status mixing characterised clusters across epidemic waves, cluster
size categories, geographic-dispersion categories, and social groups defined by
SIMD deprivation quintile.

Part 2 is explicitly descriptive. The goal is not to estimate vaccine
effectiveness against infection or transmission, but to understand the
vaccination composition of the sequenced infected population as it changed
during the rollout and across successive variant waves. Vaccination status among
sequenced cases is strongly confounded by rollout timing, age eligibility,
testing behaviour, prior immunity, and sequencing selection. These confounders
make causal inference inappropriate without dedicated test-negative or cohort
study designs.

The analysis connects directly to Part 1. Where Part 1 showed that deprivation
was not associated with larger or more geographically dispersed clusters in a
simple positive direction, Part 2 examines whether the vaccination profiles of
those clusters — and their social and demographic mixing — changed as vaccine
coverage expanded through the population.

## 2. Analysis Population

Part 2 uses the same primary analysis dataset as Part 1:

- QC filter: `nextclade_qc == "good"`
- Primary Leiden resolution: `0.3`
- Unique sequences: 281,320
- Total inferred genomic clusters: 193,160
- Singleton clusters: 109,080 (56.5%)
- Non-singleton clusters: 84,080 (43.5%)
- Epidemic waves: B.1.177, Alpha, Delta, BA.1, BA.2, BA.4, BA.5, BQ.1, XBB, Other

The non-singleton population of 84,080 clusters is the primary analysis
population for mixing and category analyses. Category thresholds are estimated
within this population.

## 3. Key Variables

### 3.1 Vaccination characterisation variables

For each cluster:

- `cluster_prop_vaccinated`: proportion of cluster members with known vaccinated
  status
- `cluster_vaccination_profile`: none vaccinated / mixed vaccination / all
  vaccinated
- `cluster_prop_boosted_all_members`: proportion of all members with a booster
  dose recorded
- `mean_dose_vaccinated_members`: mean dose number among vaccinated members
- `mean_days_since_vaccination`: mean days since last prior vaccination among
  vaccinated members (dose recency)
- `vaccination_mixing_category`: homogeneous / baseline / mixed, based on
  observed-minus-expected pairwise discordance for vaccination status

### 3.2 Cluster categories

Cluster size and geographic-dispersion categories are based on non-singleton
percentile thresholds:

| Category           | Size rule         | Geography rule              |
|--------------------|-------------------|-----------------------------|
| small/moderate     | cluster_size < 13 | datazones < 11              |
| large              | 13 ≤ size < 74    | 11 ≤ datazones < 61         |
| very large         | size ≥ 74         | datazones ≥ 61              |

Overall SIMD quintile: quintile 1 is most deprived, quintile 5 is least
deprived. Non-singleton clusters are concentrated in Q3 (37.6%) with a
right-skewed distribution towards less-deprived quintiles (Q4: 20.9%, Q5: 7.5%)
and a smaller tail in Q1 (8.9%) and Q2 (25.1%).

### 3.3 Demographic mixing categories

The four demographic mixing dimensions (SIMD deprivation, age, sex, and joint
SIMD-age-sex profile) are categorised from the Part 1 excess-discordance scores
using a ±0.01 baseline band. Across all non-singleton clusters:

| Mixing dimension | Less mix | Baseline | More mix |
|------------------|:--------:|:--------:|:--------:|
| SIMD             |   54.4%  |    6.8%  |   38.8%  |
| Age              |   38.2%  |    7.7%  |   54.2%  |
| Sex              |   34.5%  |    8.9%  |   56.6%  |
| Joint profile    |   14.9%  |   75.7%  |    9.4%  |

The majority of non-singleton clusters are more mixed than expected by age and
sex, while the majority are less mixed than expected by SIMD deprivation
quintile. Joint profile mixing is predominantly at baseline, reflecting that
the full SIMD-age-sex combination is a highly specific category where random
within-stratum variation dominates.

## 4. Key Findings

### 4.1 Breakthrough-case frequency rose sharply and earlier in older and less deprived groups

The weekly proportion of sequenced cases who were vaccinated increased rapidly
from near zero in 2021 Q1 to above 80–90% by 2022 across most age and SIMD
groups. Consistent with the JCVI rollout sequence, older age groups (65–74 and
75+) reached high breakthrough proportions first. There were clear deprivation
gradients in the pace of breakthrough-case accumulation during the initial
primary-course rollout: less-deprived quintiles (Q4–Q5) showed elevated
proportions earlier, while the most deprived quintile (Q1) lagged behind during
primary vaccination. Gradients converged by the booster phase but differences
by deprivation were visible throughout.

### 4.2 Cluster vaccination profiles shifted monotonically across waves

In B.1.177, 89.5% of non-singleton clusters contained no vaccinated members
and only 0.5% were all-vaccinated. By XBB, these proportions had entirely
reversed: 80.6% of non-singleton clusters were all-vaccinated and none contained
exclusively unvaccinated members.

Mixed-vaccination clusters were most prevalent during Delta (65.3%) and Alpha
(56.2%), when the rollout was spreading rapidly across the adult population and
within-cluster heterogeneity in vaccination status was maximal. From BA.2
onwards, all-vaccinated clusters were the majority and continued to grow.

| Wave    | None vaccinated | Mixed vaccination | All vaccinated |
|---------|----------------:|------------------:|---------------:|
| B.1.177 |           89.5% |             10.0% |           0.5% |
| Alpha   |           41.5% |             56.2% |           2.3% |
| Delta   |           10.3% |             65.3% |          24.4% |
| BA.1    |            4.4% |             49.9% |          45.7% |
| BA.2    |            2.6% |             41.0% |          56.4% |
| BA.5    |            1.8% |             27.1% |          71.0% |
| BQ.1    |            0.6% |             20.5% |          78.9% |
| XBB     |            0.0% |             19.4% |          80.6% |

### 4.3 Vaccination-status mixing within clusters peaked during mixed-coverage waves

Vaccination-status mixing categories show a clear wave pattern that tracks the
trajectory of population vaccine coverage.

In the pre-rollout B.1.177 era, 53.9% of clusters were at baseline — expected,
since near-zero coverage means expected discordance was already near zero and
there was little scope for excess mixing. In Alpha and Delta, when the rollout
was actively distributing vaccine heterogeneously across the population, 39.9%
and 49.1% of non-singleton clusters respectively showed mixed vaccination-status
composition (more mixed than expected). Homogeneous clusters (less mixed than
expected) became the dominant category from BA.1 onwards, reaching 66.3% in
BA.5, as the majority of sequenced cases were vaccinated and vaccination became
the near-universal background condition. XBB shifted back toward baseline
(56.6%), consistent with a highly vaccinated sequenced population where
vaccination status has little residual discriminatory power.

| Wave    | Homogeneous | Baseline | Mixed |
|---------|------------:|---------:|------:|
| B.1.177 |       36.7% |    53.9% |  9.5% |
| Alpha   |       57.5% |     2.6% | 39.9% |
| Delta   |       46.9% |     4.0% | 49.1% |
| BA.1    |       57.4% |     2.0% | 40.6% |
| BA.2    |       62.2% |     3.8% | 34.0% |
| BA.5    |       66.3% |    10.3% | 23.5% |
| BQ.1    |       51.7% |    29.7% | 18.6% |
| XBB     |       29.5% |    56.6% | 14.0% |

### 4.4 Dose recency increased monotonically, plateauing in late Omicron

Median days since last vaccination among vaccinated cluster members increased
from 13 days in B.1.177 to 218 days in XBB. The largest jump occurred between
Alpha (28 days) and Delta (108 days), reflecting the rollout progressing from
initial doses to the period when primary-course doses were several months in
the past. In BA.4 through XBB, median dose recency stabilised at roughly 208–218
days — consistent with a period when the most recent booster campaigns had been
completed months earlier and no new major rollout had occurred.

### 4.5 Deprivation gradients in booster coverage and dose recency are modest but present

Among all-wave non-singleton clusters, there is modest variation in booster
coverage and days since vaccination across SIMD quintiles within waves. The
Q5 − Q1 gradient in booster coverage (proportion of vaccinated members with a
booster dose) is positive for most post-Delta waves, meaning less-deprived
clusters had higher booster coverage. The gradient in mean days since last dose
is largely driven by rollout timing, but access deprivation in particular shows
a different pattern from overall SIMD in some domains.

The SIMD-domain analysis shows that these dose-recency and booster-coverage
gradients are not uniform across deprivation domains. Income, employment, and
health domains broadly follow the overall SIMD gradient, while geographic access
deprivation is systematically different in direction and magnitude from the
other domains in several waves.

### 4.6 Large clusters show stronger homogeneous mixing by SIMD and size

Cross-category heatmaps show that SIMD-deprivation "less mix" fractions are
substantially higher for large and very large clusters than for small/moderate
clusters. Among large clusters, 85–100% of non-singleton clusters are classified
as less mixed than expected by SIMD, compared with 16–55% of small/moderate
clusters. This is consistent with larger genomic clusters being more geographically
concentrated and therefore more socioeconomically homogeneous by SIMD quintile.

Age mixing shows an opposite pattern in small/moderate clusters: the majority
(57–65% across quintiles) are classified as more mixed than expected by age,
consistent with broader age exposure in smaller clusters. Sex mixing shows
little systematic gradient across size or SIMD.

Joint profile mixing is predominantly at baseline or less-mixed for large
clusters (79–89% less mix in large clusters across some quintiles), reflecting
that large clusters tend to be internally homogeneous on the combined
SIMD-age-sex profile.

## 5. Interpretation

The Part 2 results should be read alongside Part 1 rather than in isolation.
Part 1 showed that deprivation was not associated with larger or more
geographically dispersed clusters after adjustment for surveillance conditions,
and that deprivation signals were clearest in demographic mixing rather than
cluster scale. Part 2 adds context to these findings:

The vaccination profile of clusters tracked the population rollout faithfully.
Clusters were not exceptional relative to background vaccination coverage in
their wave and lineage stratum; instead, the proportion of vaccinated members
reflected where each wave fell in the rollout calendar. Mixed vaccination-status
clusters peaked when vaccination was most heterogeneously distributed across the
population, not because clusters were intrinsically different.

The homogeneous mixing trend from BA.1 onwards — where most clusters were less
vaccinated-status-mixed than expected — does not imply transmission was
occurring preferentially within vaccination groups. Rather, it reflects that in
a predominantly vaccinated population, the expected discordance within a sampled
stratum is already high, and clusters are more likely to fall below that
expectation simply because vaccination status is no longer a strong discriminant.

The deprivation gradients in breakthrough-case rates, booster coverage, and dose
recency are consistent with the Part 1 findings of outcome-specific and
domain-specific deprivation signals. They suggest that unequal vaccine access
or uptake among deprived groups was measurable in the sequenced infected
population, but that these gradients were not large enough to produce the
simple pattern of more and larger clusters in deprived areas that might be
predicted by a naive vaccine-protection hypothesis.

## 6. Takeaway

In Scottish SARS-CoV-2 genomic clusters, vaccination profiles reflected the
national rollout trajectory rather than any exceptional clustering of vaccinated
or unvaccinated cases. Mixed vaccination-status clusters were most common during
Delta and BA.1, when the rollout was most heterogeneous. From BA.2 onwards,
all-vaccinated clusters dominated and vaccination-status homogeneity increased.
Modest SIMD gradients in booster coverage and dose recency were present but
varied by SIMD domain, with geographic access deprivation behaving differently
from income or employment deprivation. Large clusters showed stronger
SIMD-deprivation homogeneity, consistent with geographic concentration of large
transmission events. These findings complement the Part 1 conclusion that
deprivation and surveillance conditions, rather than vaccination status per se,
are the primary structural determinants of observed genomic cluster scale.
