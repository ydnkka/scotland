# Scottish SARS-CoV-2 Genomic Clustering: Part 3 Summary Report

Prepared: 7 May 2026

## 1. Introduction

Part 3 examines the associations between Scottish government COVID-19 policy
restriction periods and the structure of SARS-CoV-2 genomic clusters. The
central question is whether policy restriction intensity — the severity of the
restrictions in force during each epidemic period — is associated with cluster
size, geographic dispersion, or demographic mixing, after acknowledging that
policy periods are strongly confounded with variant waves and calendar time.

This analysis is explicitly descriptive and associational. Scotland implemented
16 distinct restriction periods between March 2020 and April 2022, ranging from
the initial emergence period through first and second national lockdowns, a
tiered framework, the Omicron-wave restrictions, and final easing. The
genomic data begin in July 2020, covering 12 of these 16 periods. Policy
periods and variant waves are nearly perfectly collinear in calendar time:
the second lockdown coincides with Alpha's emergence, near-normal conditions
coincide with Delta's peak, and the Omicron-wave restrictions coincide with
BA.1's emergence. Causal inference about policy effects is therefore not
possible, and Part 3 makes no such claims.

The value of Part 3 is threefold: (1) it completes the descriptive
characterisation of how cluster structure varied across the full arc of the
epidemic, anchored to the policy timeline rather than to variant waves; (2) it
provides three targeted interrupted-time-series (ITS) analyses at transitions
selected for minimal wave confounding; (3) it contextualises the Part 1 and
Part 2 findings by showing where in the policy calendar the key structural
changes in clusters occurred.

## 2. Analysis Population

Part 3 uses the Part 1 primary analysis dataset:

- QC filter: `nextclade_qc == "good"`
- Primary Leiden resolution: 0.3
- Policy period assignment: `wn_mid_date` matched to `utils/policy.py`
- Study start: 2020-07-10 (P3, Route-map phase 3)
- Total clusters: 193,112
- Non-singleton clusters: 84,067 (43.5%)
- Policy periods observed: P3, T1, F5, L2, SL, L3, L21, L0, NN, OM, FE, PR

## 3. Policy Period Context

| Code | Label | Intensity | n clusters | Median size (NS) | Median DZ (NS) |
|---|---|---:|---:|---:|---:|
| P3 | Route-map phase 3 | 30 | 1,570 | 4.0 | 3.0 |
| T1 | Pre-tier tightening | 55 | 1,379 | 4.0 | 3.0 |
| F5 | Five-tier framework | 65 | 2,428 | 3.0 | 3.0 |
| L2 | Second lockdown | 95 | 9,247 | 4.0 | 3.0 |
| SL | Stay local — Level 3 | 65 | 1,689 | 5.0 | 3.0 |
| L3 | Level 3 | 55 | 1,169 | 5.0 | 3.0 |
| L21 | Level 2 / Level 1 | 38 | 10,184 | 4.0 | 3.0 |
| L0 | Level 0 | 20 | 5,245 | 3.0 | 3.0 |
| NN | Near-normal | 10 | 48,407 | 3.0 | 3.0 |
| OM | Omicron wave | 42 | 21,394 | 3.0 | 3.0 |
| FE | Final easing | 15 | 38,699 | 3.0 | 3.0 |
| PR | Post-restriction | 3 | 51,701 | 3.0 | 2.0 |

*NS = non-singleton clusters; DZ = datazones.*

The Near-normal (NN, Delta) and Post-restriction (PR) periods account for the
majority of sequenced clusters by volume (48,407 and 51,701 respectively),
reflecting both the high incidence of these waves and the substantial sequencing
effort sustained through 2021–2022.

## 4. Key Findings

### 4.1 Policy intensity shows strong positive correlation with cluster size but not SIMD mixing

Spearman correlations between weekly median log cluster size and weekly policy
intensity are strongly positive (ρ = 0.74, n = 134 weeks). Weekly log
datazones and intensity are moderately positively correlated (ρ = 0.58). Mean
age excess discordance correlates positively with intensity (ρ = 0.59). Mean
SIMD excess discordance shows effectively no correlation with intensity
(ρ = 0.02). These pooled correlations are largely driven by wave confounding —
the same calendar windows with high restriction intensity also correspond to
waves (B.1.177, Alpha) that produced the largest, most geographically
concentrated clusters — and do not represent causal policy effects.

### 4.2 No detectable change in cluster structure at T1-onset (Oct 2020)

The introduction of the five-level tiered framework (Pre-tier tightening, T1)
in October 2020 was not associated with significant changes in cluster size
(β = −0.13, 95 % CI −0.34 to +0.09, p = 0.22), geographic dispersion (β =
−0.07, p = 0.59), or mixing metrics (all p > 0.48). This is the cleanest
within-variant pre-vaccine ITS test in the dataset, occurring entirely within
the B.1.177 era. The null result may reflect a lack of substantive change in
actual contact behaviour between P3 and T1, a lag between policy change and
genomic cluster formation, or insufficient statistical power in the 16-week
ITS window.

### 4.3 Cluster size and geographic spread declined at the L2→SL transition (Apr 2021)

When the second lockdown eased to Stay-local Level 3 on 2 April 2021, there
were statistically significant reductions in both median log cluster size
(β = −0.21, 95 % CI −0.41 to −0.02, p = 0.034) and log datazones
(β = −0.36, 95 % CI −0.60 to −0.12, p = 0.006). This is counter-intuitive
in the direction of the association (restrictions eased but clusters got
smaller), and the most likely explanation is the Alpha wave tail confounding:
the April 2021 transition coincides with the declining phase of the Alpha
wave, when cluster sizes were already contracting independently of the policy
change. The short duration of the SL period (24 days) limits the interpretive
value of this transition.

### 4.4 Geographic dispersal increased significantly at NN-onset (Aug 2021)

The full removal of legal distancing requirements at the Near-normal (NN)
onset on 9 August 2021 was associated with a significant increase in cluster
geographic spread: log datazones increased by β = +0.32 (95 % CI +0.07 to
+0.56, p = 0.015). The corresponding change in log cluster size was positive
but fell short of significance (β = +0.12, p = 0.068). SIMD and age mixing
metrics were not significantly affected (both p > 0.67). The geographic
dispersal signal is consistent with the hypothesis that removal of gathering
restrictions allowed transmission to span a wider geographic footprint within
the Delta wave, though simultaneous high Delta incidence means the policy and
wave effects cannot be disentangled.

### 4.5 Demographic mixing metrics were unaffected by any of the three transitions

Across all three ITS transitions, neither SIMD excess discordance nor age
excess discordance showed statistically significant level changes (all p > 0.27
for level-change terms). This is consistent with the Part 1 finding that
demographic mixing within clusters is not strongly driven by the overall
restriction context but instead reflects the local social structure of
transmission events. Policy changes that affect contact rates globally do not
appear to substantially alter the socioeconomic composition of who transmits
to whom within a cluster.

## 5. Interpretation

The Part 3 results complement Parts 1 and 2 rather than superseding them.

The strong pooled correlation between policy intensity and cluster size (ρ =
0.74) is informative precisely because it is largely confounded: it reflects
the fact that the periods of highest restriction produced the largest, most
genomically coherent clusters (household transmission during lockdowns), while
the post-restriction period produced many small, geographically diffuse
clusters (community transmission in a largely immune population). This is the
expected signature of epidemic suppression and should not be interpreted as
evidence that restrictions per se generated larger clusters.

The ITS analyses are more informative at the transition level. The null result
at T1-onset suggests that the formal introduction of a tiered policy framework
in late 2020 did not produce an acute change in genomic cluster structure
detectable at the weekly scale. The NN-onset result — increased geographic
dispersal at full legal easing — is the most a priori consistent finding,
and its statistical significance survives the short observation window.

The absence of any mixing-metric signal across all three transitions reinforces
the Part 1 conclusion that demographic mixing within clusters is determined by
the local social structure of transmission (which household, neighbourhood, and
community you belong to) rather than by the overall policy environment. Policies
that restrict contacts broadly do not appear to change the socioeconomic
composition of the contacts that still occur.

## 6. Takeaway

In Scottish SARS-CoV-2 genomic clusters, policy restriction intensity is
positively correlated with cluster size and geographic spread across the
epidemic as a whole, but this association is largely a consequence of wave
confounding rather than a direct policy effect. Targeted ITS analyses at three
within-variant policy transitions find: no significant structural change at the
October 2020 tier introduction; a paradoxical reduction in cluster scale at
the April 2021 lockdown lift (interpreted as Alpha wave tail dynamics); and a
significant increase in cluster geographic dispersal at the August 2021 full
legal easing (consistent with removal of geographic mixing constraints at the
start of Delta's community spread). Demographic mixing metrics were unaffected
by all three transitions, supporting the Part 1 conclusion that socioeconomic
composition of clusters is driven by local social structure rather than the
broad policy environment.
