# Scottish SARS-CoV-2 Genomic Clustering: Part 4 Summary Report

Prepared: 8 May 2026

## 1. Introduction

Part 4 focuses on the emergence of Alpha in Scotland as a case study in how a
new variant can arise within, and then overtake, an existing epidemic and
policy context. The analysis follows Alpha from its cryptic early signal after
the easing phase of the first lockdown, through the decline of B.1.177, the
introduction of the five-tier framework, and the subsequent introduction of the
second national lockdown.

The central question is whether the timing of restriction measures could have
changed the Alpha trajectory. This is addressed descriptively using mutation
frequency trajectories, overlapping-window cluster evidence, demographic
profiles, regional footprint summaries, and simple fitted counterfactuals. The analysis does
not identify a specific transmission event and does not make causal claims
about policy effects.

## 2. Analysis inputs and definitions

Part 4 uses the same primary genomic filtering conventions as the earlier
parts:

- QC filter: `nextclade_qc == "good"`
- Primary Leiden resolution: 0.3
- Alpha definition for Figure 1 and phase summaries: Pango lineage beginning
  `B.1.1.7`
- Alpha marker for time-series modelling: S:N501Y
- B.1.177 marker: S:A222V
The key restriction periods for this analysis are:

| Code | Period | Start | End | Intensity |
|---|---|---:|---:|---:|
| P3 | Route-map phase 3 | 2020-07-10 | 2020-10-01 | 30 |
| T1 | Pre-tier tightening | 2020-10-02 | 2020-11-01 | 55 |
| F5 | Five-tier framework | 2020-11-02 | 2021-01-04 | 65 |
| L2 | Second lockdown | 2021-01-05 | 2021-04-01 | 95 |

## 3. Early Alpha phase summaries

The phase summaries below are based on unique sequences for demographic
interpretation. Cluster-weighted counts are still retained in the source table
because they describe overlapping-window cluster persistence, but they can
duplicate the same sequence across adjacent windows.

| Phase | Windows | Clusters | Summed cluster size | Unique sequences | Duplicate memberships | Unique demographic profile |
|---|---|---:|---:|---:|---:|---|
| Cryptic GGC chain | W016-W021 | 15 | 96 | 51 | 45 | Mean age 51.1; 54.9% female; 47.1% SIMD 1; 100.0% unvaccinated |
| Multi-region expansion | W022-W024 | 113 | 577 | 291 | 286 | Mean age 49.5; 54.0% female; 28.5% SIMD 1; 95.5% unvaccinated |
| F5/L2 bridge | W025 | 102 | 458 | 458 | 0 | Mean age 45.0; 53.5% female; 26.2% SIMD 1; 95.9% unvaccinated |

The cryptic phase is the most sensitive to duplicate counting. In the
cluster-weighted view it appears older, more female, and more deprived
(weighted mean age 52.7; 74.0% female; 93.8% SIMD 1). In the unique-sequence
view, the same phase remains enriched for more deprived areas but is much less
demographically extreme (mean age 51.1; 54.9% female; 47.1% SIMD 1). The
unique-sequence profile is therefore the preferred demographic interpretation.

## 4. Main figures

### 4.1 Regional expansion of Alpha

This figure summarises the regional footprint of unique Alpha sequences across
three early phases: W016-W021, W022-W024, and W025. Panel A shows the expansion
in total unique sequences and number of affected health boards; Panel B shows
the health-board-by-phase heatmap.

![Figure 1](figures/fig1_alpha_geographic_emergence.png)

### 4.2 Alpha/B.1.177 transition and hospital pressure

This figure places Alpha growth alongside the decline of B.1.177, policy
periods, and Scotland hospital occupancy.

![Figure 2](figures/fig2_alpha_b1177_transition.png)

### 4.3 Counterfactual L2 timing and fitted growth rates

This figure compares observed Alpha growth with counterfactual earlier L2
timings, and summarises the fitted growth/decline rates used in those
scenarios.

![Figure 3](figures/fig3_counterfactual_growth.png)

## 5. Key findings

### 5.1 Alpha was initially detected as a small, demographically specific signal

The earliest Alpha phase contained 51 unique sequences across 47 datazones
and 7 health boards, with a unique-sequence profile of mean age 51.1 years,
54.9% female, 47.1% SIMD 1, and 100% unvaccinated. This supports wording that
the early signal was associated with more deprived areas and a pre-vaccine
population, but it should not be described as exclusively older or female once
duplicate sequence memberships are removed.

### 5.2 The overlap pattern supports a possible superspreading or rapid-expansion window

The strongest overlap evidence occurs across the W022-W024 expansion period:
`W022|B.1.1.7|R0.3|C6` shares 60 sequences with
`W023|B.1.1.7|R0.3|C4` (95.2% of W022/C6), and W023/C4 shares 86 sequences
with `W024|B.1.1.7|R0.3|C1` (93.5% of W023/C4; 98.9% of W024/C1). This is
consistent with a rapid expansion or superspreading-like transmission window
around the 8 December rise, although the sequence data cannot identify the
specific setting or event.

### 5.3 Alpha replaced B.1.177 during F5 and was already established by L2

The Alpha marker S:N501Y increased from 3.2% in W021 to 17.7% in W022, while
the B.1.177 marker S:A222V declined through December. By the start of L2 on
5 January 2021, Alpha was already near or above majority frequency in the
observed and fitted trajectories. This timing is important for interpretation:
the second lockdown began after the key Alpha expansion was underway.

### 5.4 L2 was associated with slower Alpha growth, but not reversal

The primary model fitted a binomial GLM to the sequenced S:N501Y counts, with
additional window-level weighting by confirmed positive-test volume. The fitted
Alpha growth rate was 0.0812/day under F5 (95% CI 0.0711-0.0913; doubling time
8.5 days) and 0.0661/day under L2 (95% CI 0.0619-0.0702; doubling time
10.5 days). The L2/F5 rate ratio was 81.4%, so the fitted Alpha growth rate
under L2 was 19% slower than under F5. B.1.177 declined under L2 with an
estimated halving time of 10.6 days.

Sensitivity analyses compared this primary model with two alternatives:

| Model | F5 r/day | L2 r/day | L2/F5 | L2 slower |
|---|---:|---:|---:|---:|
| Unweighted binomial GLM | 0.0789 | 0.0694 | 87.9% | 12.1% |
| Positive-test weighted GLM | 0.0812 | 0.0661 | 81.4% | 18.6% |
| Coverage-adjusted GLM | 0.0865 | 0.0433 | 50.0% | 50.0% |

The direction of the result is stable across models: Alpha grew more slowly in
L2 than in F5. The estimated magnitude is sensitive to how sequencing coverage
is handled.

### 5.5 Earlier L2 timing would delay, not prevent, Alpha dominance in the fitted counterfactuals

In the counterfactual projections, Alpha reaches 50% frequency on:

| Scenario | 50% Alpha date |
|---|---:|
| Actual F5 -> L2 on 5 Jan | 2021-01-05 |
| L2 from 8 Dec | 2021-01-12 |
| L2 from 2 Dec | 2021-01-19 |
| L2 from 2 Nov | 2021-01-26 |

These projections suggest that timely restriction measures may have changed
the trajectory by delaying Alpha dominance and reducing the overlap between
the Alpha rise and the B.1.177 winter burden. Under the fitted assumptions,
however, even very early L2-level restrictions delay dominance by weeks rather
than preventing it.

## 6. Manuscript interpretation

The strongest Part 4 framing is that Alpha did not emerge into a neutral
background. It emerged after Phase 3 easing, during a B.1.177-dominated autumn
wave, through the tiered and five-tier frameworks, and before the second
lockdown could fully affect transmission. The demographic analysis helps
separate real early enrichment from duplicated overlapping-window counts, while
the counterfactual model supports a cautious conclusion: earlier restrictions
could plausibly have shifted the timing and burden of the Alpha transition,
but the available genomic data do not support a claim that policy timing alone
would have prevented Alpha establishment.
