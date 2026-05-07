# Part 3 Results and Figures Description

This document provides figure-by-figure narrative guidance for the Part 3
manuscript. Each section describes what the figure shows, what the reader
should take from it, and a suggested results paragraph.

---

## Figure 1: Weekly cluster outcomes and policy context
**File:** `fig1_weekly_time_series`

**What it shows.** Two-panel time series covering the full study period
(July 2020 – January 2023). Panel A shows the weekly median log cluster size
for non-singleton clusters; Panel B shows the policy intensity stepped series
with period code labels annotating each period. Both panels have background
shading coloured from blue (low intensity) to red (high intensity) for each
policy period. Three vertical dashed lines mark the ITS transition dates.

**Key visual patterns.** Panel A makes the wave structure of cluster sizes
immediately visible: log cluster size peaks during the second lockdown (L2,
Alpha wave, early 2021) and again during the Omicron BA.1 wave (OM, late 2021),
then declines steadily through 2022 as immunity accumulates. The shading shows
that the largest-cluster periods coincide with the highest-intensity periods —
the visual correlation that motivates the analysis while also illustrating
its confounded nature. Panel B shows the intensity staircase: the dramatic
rise to lockdown-level intensity in January 2021, the stepped relaxation
through 2021, and the low-intensity regime from August 2021 onwards.

**Suggested results paragraph.**
> Figure 1 illustrates the policy context of the Part 3 analysis. Median
> log cluster size among non-singleton clusters followed a clear wave
> structure, peaking during the second lockdown and early Omicron periods and
> declining through 2022 (Panel A). The policy intensity series (Panel B) shows
> the full arc from the initial restriction tightening to the final legal
> relaxation in April 2022. Visual inspection confirms the positive association
> between restriction intensity and cluster size noted in the Spearman
> correlation analysis (ρ = 0.74), but also illustrates that the periods of
> highest intensity correspond directly to the Alpha and Omicron wave peaks,
> making the association uninterpretable without wave adjustment. The three
> dashed lines mark the ITS transition dates analysed in Figure 2.

---

## Figure 2: Interrupted time-series at three policy transitions
**File:** `fig2_its_transitions`

**What it shows.** A 3×2 panel of ITS analyses. Rows correspond to the three
chosen transitions: T1-onset (Oct 2020), L2→SL (Apr 2021), and NN-onset
(Aug 2021). The left column shows weekly median log cluster size; the right
column shows weekly median log datazones. Each panel plots observed weekly
data points (blue = pre-transition, red = post-transition) with fitted OLS
segmented-regression lines. The level-change estimate (β_post) with 95 % CI
is annotated in each panel; an asterisk (*) denotes p < 0.05.

**Key visual patterns.** The T1-onset row (top) shows flat or slightly
declining trends in both outcomes, with the fitted pre- and post-lines nearly
continuous at the transition — consistent with the null ITS result. The
L2→SL row (middle) shows a downward level shift: post-transition points are
visibly lower than the pre-transition trend would predict, especially for
datazones. The NN-onset row (bottom) shows an upward level shift for datazones
in particular, with the post-transition points tracking above the pre-transition
line — the strongest and most interpretable signal in the figure.

**Suggested results paragraph.**
> Figure 2 presents the interrupted time-series analyses at the three selected
> policy transitions. At T1-onset (pre-tier tightening introduction in October
> 2020; top row), neither cluster size nor geographic spread showed a
> statistically significant level change (β_size = −0.13, 95 % CI −0.34 to
> +0.09, p = 0.22; β_DZ = −0.07, p = 0.59), consistent with the null
> hypothesis that this policy formalisation did not produce an acute change in
> genomic cluster structure. At the second lockdown lift (L2→SL, April 2021;
> middle row), significant downward level shifts were observed for both cluster
> size (β = −0.21, 95 % CI −0.41 to −0.02, p = 0.034) and datazones (β =
> −0.36, 95 % CI −0.60 to −0.12, p = 0.006); these counter-intuitive
> reductions are most plausibly attributed to the naturally declining Alpha
> wave rather than to the easing of restrictions per se. At NN-onset (full
> legal easing, August 2021; bottom row), log datazones showed a significant
> positive level shift (β = +0.32, 95 % CI +0.07 to +0.56, p = 0.015),
> suggesting that the removal of gathering and distancing rules was associated
> with more geographically dispersed transmission within the Delta wave.

---

## Figure 3: Policy-period cluster outcome dot chart
**File:** `fig3_period_outcomes`

**What it shows.** Horizontal dot chart with one point per observed policy
period. Left panel: median log cluster size (non-singletons). Right panel:
median log datazones (non-singletons). Points are coloured by policy intensity
using the RdYlBu_r scale (blue = low, red = high), with a shared colourbar
at the right. Each point is annotated with the non-singleton cluster count
for context on period volume. Periods are ordered chronologically from top
(P3) to bottom (PR).

**Key visual patterns.** Both panels show that the largest median values
occur in the middle of the chronological range — specifically during the
SL, L3, and L2 periods (Alpha-dominant, early 2021) — rather than at the
highest-intensity periods (L2, second lockdown). The PR period (post-
restriction) has the smallest median values on both outcomes, consistent with
the very large number of small, geographically compact post-Omicron clusters.
The colour gradient from blue (P3, L0, NN, PR) to red (L2, F5, T1) adds the
intensity context without implying a causal direction.

**Suggested results paragraph.**
> Figure 3 compares median log cluster size and log datazones across the 12
> observed policy periods. The largest cluster sizes were observed in the
> Stay-local (SL) and Level 3 (L3) periods during the Alpha wave (median
> cluster size 5 in both periods), not during the second lockdown peak (L2,
> median size 4), likely reflecting that Alpha clusters were larger due to
> the emerging variant's higher transmissibility rather than the policy context.
> The Post-restriction (PR) period produced the smallest median cluster sizes
> and datazones (median size 3, median datazones 2), consistent with endemic
> low-intensity transmission in a highly immune population. The intensity
> colourbar highlights that the high-volume periods (NN, PR) are low-intensity,
> while the periods with the largest clusters are broadly spread across the
> intensity range, further illustrating the confounded nature of the
> intensity–cluster-size association.

---

## Supplementary Figure 1: Weekly mixing metrics with policy overlay
**File:** `supp_fig1_weekly_mixing`

**What it shows.** Four-panel weekly time series (2×2). Panel A: mean SIMD
excess discordance. Panel B: mean age excess discordance. Panel C: median log
datazones. Panel D: policy intensity stepped series. All panels share the same
calendar axis (July 2020 – January 2023) with background shading and ITS
transition marker lines.

**Key visual patterns.** The SIMD excess discordance series (Panel A) is
notably flat and close to zero throughout the study period, with no visible
wave structure or policy-period trend — consistent with the near-zero
Spearman ρ of 0.02 and with the Part 1 finding that SIMD mixing is not
strongly associated with the epidemic trajectory. The age excess discordance
series (Panel B) shows more temporal variation, with elevated values during
the pre-vaccine and Alpha periods, declining through 2022. Panel C mirrors
Figure 1A. Panel D provides the intensity reference.

**Suggested results paragraph.**
> Supplementary Figure 1 shows weekly mixing metric time series alongside
> the policy intensity series. SIMD excess discordance (Panel A) showed no
> systematic temporal trend across the study period (Spearman ρ with intensity:
> 0.02), consistent with the absence of significant ITS-level changes at any
> of the three policy transitions. Age excess discordance (Panel B) showed
> more temporal variation, correlating positively with intensity (ρ = 0.59)
> in a pattern driven primarily by the Alpha-period peak. No step changes
> coinciding with the ITS transition dates are visible for either mixing metric,
> reinforcing the conclusion that demographic mixing within clusters was not
> acutely responsive to policy-level changes in the restriction environment.
