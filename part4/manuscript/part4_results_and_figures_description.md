# Part 4 Results and Figures Description

This document provides figure-by-figure narrative guidance for the Part 4
manuscript. Each section describes what the figure shows, what the reader
should take from it, and a suggested results paragraph.

---

## Figure 1: Regional expansion of Alpha
**File:** `fig1_alpha_geographic_emergence`

**What it shows.** A two-panel regional expansion summary for Alpha across the
early seeding, expansion, and F5/L2 bridge phases. Panel A shows the total
number of unique Alpha sequences in each phase, with the number of health
boards overlaid as a line. Panel B is a health-board-by-phase heatmap of
unique Alpha sequence counts. Counts are unique sequences, not summed cluster
memberships, so sequences duplicated across overlapping cluster windows are
counted only once within each phase.

**Key visual patterns.** The cryptic phase (W016-W021; 4 Nov-11 Dec 2020)
contains 51 unique Alpha sequences across 7 health boards, with the strongest
signal in Greater Glasgow and Clyde. The expansion phase (W022-W024;
30 Nov-1 Jan) increases to 291 unique sequences across 11 health boards, and
the heatmap shows broader contributions from Grampian, Lanarkshire, Lothian,
Dumfries and Galloway, Borders, and Tayside. By the F5/L2 bridge week (W025;
19 Dec-8 Jan), Alpha reaches 458 unique sequences across 11 health boards,
with Greater Glasgow and Clyde still dominant but a visibly wider regional
footprint.

**Suggested results paragraph.**
> Figure 1 summarises the regional emergence of Alpha across the three early
> phases defined from the overlapping-window cluster analysis. During the
> cryptic GGC chain (W016-W021), Alpha was visible as 51 unique sequences
> across 7 health boards, with the strongest signal in Greater Glasgow and
> Clyde. The subsequent expansion phase (W022-W024) broadened sharply, with
> 291 unique sequences detected across 11 health boards, and by W025 the F5/L2
> bridge phase included 458 unique sequences across 11 health boards. The
> heatmap supports a transition from localized or weakly detected seeding into
> multi-region expansion before the second national lockdown could fully alter
> transmission dynamics.

---

## Figure 2: Alpha/B.1.177 transition and hospital pressure
**File:** `fig2_alpha_b1177_transition`

**What it shows.** A two-panel time series from October 2020 to April 2021.
Panel A plots the frequency of S:N501Y, used here as an Alpha marker, and
S:A222V, used here as a B.1.177 marker. Panel B plots sequenced lineage
composition as a stacked area chart, with Scotland hospital occupancy overlaid
on a secondary axis. Policy-period shading and vertical markers indicate the
pre-tier tightening period (T1), the five-tier framework (F5), the second
lockdown (L2), and the stay-local period (SL).

**Key visual patterns.** S:N501Y remains low through late November, then rises
from 3.2% in W021 to 17.7% in W022, the 8 December expansion point annotated
in the figure. B.1.177 remains dominant during T1/F5 but is rapidly displaced
as Alpha rises through December and January. By the L2 transition on 5 January
2021, Alpha is already around the majority threshold in the fitted and observed
series, while hospital occupancy is rising toward its winter peak.

**Suggested results paragraph.**
> Figure 2 places the Alpha emergence in the context of the preceding B.1.177
> wave and the Scottish policy timeline. The Alpha marker S:N501Y rose from
> 3.2% in W021 to 17.7% in W022, marking the abrupt expansion around
> 8 December 2020, while the B.1.177-associated S:A222V marker declined over
> the same period. Sequenced lineage composition shows the same replacement
> dynamic: B.1.177 dominated through the tiered and five-tier periods, but
> Alpha crossed into dominance around the start of L2. Hospital occupancy was
> also rising through this transition, indicating that the second lockdown was
> introduced after Alpha had already established substantial frequency and
> while winter pressure was increasing.

---

## Figure 3: Counterfactual lockdown timing and growth-rate comparison
**File:** `fig3_counterfactual_growth`

**What it shows.** A two-panel counterfactual figure. Panel A compares the
observed S:N501Y trajectory with fitted logistic projections under the actual
F5 to L2 transition and three hypothetical earlier L2 timings: 2 November,
2 December, and 8 December 2020. Panel B compares the fitted Alpha growth
rate under F5, the fitted Alpha growth rate under L2, and the fitted B.1.177
decline rate under L2.

**Key visual patterns.** Alpha grows faster under F5 than under L2 in the
primary positive-test weighted model: the fitted F5 growth rate is 0.0812 per
day (doubling 8.5 days; 95% CI 7.6-9.7), whereas the fitted L2 growth rate is
0.0661 per day (doubling 10.5 days; 95% CI 9.9-11.2). The L2/F5 rate ratio is
81.4%, meaning L2 is estimated to be 19% slower than F5 for Alpha growth. The
model-sensitivity table gives a range from 12% slower in an unweighted
binomial GLM to 50% slower after adding proportion sequenced. B.1.177 declines
under L2 with an estimated halving time of 10.6 days. In the fitted
counterfactuals, 50% Alpha frequency occurs on 5 January under the actual
policy timing, 12 January if L2 begins from 8 December, 19 January if L2 begins
from 2 December, and 26 January if L2 begins from 2 November.

**Suggested results paragraph.**
> Figure 3 evaluates whether earlier restriction timing could plausibly have
> changed the Alpha trajectory. The primary model fitted a binomial GLM to the
> sequenced S:N501Y counts, with additional window-level weighting by confirmed
> positive-test volume. It estimated faster Alpha growth during the five-tier
> framework (r = 0.0812/day, doubling 8.5 days) than during L2 (r =
> 0.0661/day, doubling 10.5 days), corresponding to a 19% slower growth rate
> under L2. Counterfactual projections suggest
> that earlier imposition of L2-level restrictions would have delayed Alpha
> dominance, shifting the 50% date from 5 January under the observed timing to
> 12 January for L2 from 8 December, 19 January for L2 from 2 December, and
> 26 January for L2 from 2 November. These projections are descriptive rather
> than causal, but they suggest that timely restriction measures may have
> changed the short-term trajectory by delaying dominance and reducing
> concurrent B.1.177 burden, rather than preventing Alpha establishment.

---

## Supporting tables for manuscript text

- `part4_alpha_phase_demographic_summary.csv`: cluster-weighted and
  unique-sequence demographic summaries for the three early Alpha phases.
- `part4_alpha_chain_overlaps.csv`: sequence-overlap evidence linking Alpha
  clusters across adjacent windows.
- `part4_growth_params.csv`: fitted growth and decline rates, doubling or
  halving times, confidence intervals, and pseudo-R2 values.
- `part4_counterfactual_projections.csv`: daily projected Alpha frequencies
  under the actual and counterfactual L2-timing scenarios.
- `part4_mutation_trajectories.csv`: weekly mutation-marker frequencies used
  in Figures 2 and 3.
