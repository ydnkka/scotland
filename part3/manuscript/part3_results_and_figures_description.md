# Part 3 Results And Figure Guide

Generated for the Part 3 policy-phase and variant-emergence analysis on 9 May 2026.

## Analysis Question

Part 3 asks whether Scottish SARS-CoV-2 genomic cluster structure was
associated with restriction-intensity policy phases across the pre-Omicron
epidemic, and whether the timing of the second national lockdown (L2) relative
to the emergence of the Alpha (B.1.1.7) variant plausibly altered its
trajectory.

The specific descriptive and modelling questions are:

- Were weekly cluster-size and geographic-spread outcomes correlated with
  policy restriction intensity across the 16 observed policy periods?
- Did three selected policy transitions produce acute changes in genomic cluster
  structure, as estimated by interrupted time-series segmented regression?
- Did restriction context (five-tier framework F5 versus second lockdown L2)
  correspond to measurably different Alpha log-odds growth rates?
- What was the pre-L2 meta-cluster structure of the early Alpha expansion, and
  what does it imply for the likely effect of earlier L2 imposition?

The main manuscript figures are generated from `part3/tables`, `part3/cache`,
and `part3/notebooks/tables`, using the publication styling in `utils/style.py`.

## Analysis Population And Modelling Frame

The analysis uses good-QC sequences at Leiden resolution 0.3 across 16 observed
policy periods from July 2020 to April 2022. The full analytical dataset
contains 48,407 clusters and 183,190 sequences in the NN period alone; the
complete cross-period population is described in `part3/tables/period_descriptives.csv`.

Policy periods are assigned a scalar restriction-intensity score (P3 = 30,
T1 = 55, F5 = 65, L2 = 95, SL = 65, NN = 10, etc.) reflecting the legal
restriction level in force. Spearman rank correlations between period-level
intensity scores and period-median cluster outcomes are used as a descriptive
summary of the intensity–cluster-structure association.

The interrupted time-series (ITS) analysis uses segmented OLS regression of the
form:

> y_t = β₀ + β₁·t + β₂·post + β₃·(post·t) + ε

where y_t is the weekly median log cluster size or median log datazones, t is a
continuous week counter, post is an indicator for the post-transition period, and
β₂ (the level-change estimate) is the primary ITS parameter of interest. Three
transitions are analysed at a primary ±8-week window: T1-onset (October 2020),
L2→SL (April 2021), and NN-onset (August 2021). ITS windows of ±6, ±10, and
±12 weeks are used for sensitivity analyses.

Alpha and B.1.177 growth models are binomial GLMs fitted to sequenced variant-
marker counts (S:N501Y for Alpha; S:A222V for B.1.177) with positive-test-volume
weighting as the primary specification. The fitted log-odds coefficient is
expressed as a growth rate r (per day), a doubling time (days), and an
odds-ratio per week. Phase-specific models are fitted separately for the F5 and
L2 policy periods.

The pre-L2 Alpha meta-cluster network is constructed by (i) identifying all
Alpha-containing rolling-window cluster assignments, (ii) linking adjacent
windows when they share at least one sequence, and (iii) extracting connected
components. Each connected component is a meta-cluster, labelled AM001 (largest)
through AM078. All unique-sequence counts are deduplicated within each
meta-cluster.

## Overall Results Summary

Because cluster size is highly right-skewed across all policy periods, two
complementary per-period metrics — clustering rate and the negative-binomial
dispersion parameter k̂ — are reported alongside the Spearman correlations as
the primary distributional summaries:

- **Clustering rate** = (total sequences − total clusters) / total sequences.
  This is the fraction of sequenced isolates that are secondary cluster members
  (as opposed to singletons). It ranges from 0 (all singletons) toward 1 (one
  large chain).
- **Dispersion parameter k̂** (MME) = x̄² / (s² − x̄), where x̄ is the
  period mean cluster size and s² is the period variance. Lower k indicates
  greater dominance by a few large "superspreading-like" clusters; k → ∞
  is Poisson.

Cross-period Spearman correlations show a positive association between
restriction intensity and median log cluster size (ρ = 0.741) and a negative
association with singleton fraction (ρ = −0.621). The association with mean SIMD
excess discordance is negligible (ρ = 0.019), consistent with Part 1's finding
that socioeconomic mixing within clusters is not strongly linked to the epidemic
trajectory.

ITS analyses at three selected transitions produce a heterogeneous picture.
The T1-onset (October 2020 tier tightening) shows no significant level change in
cluster size or geographic spread, consistent with a null acute effect. The
L2→SL transition (April 2021 lockdown lift) shows a significant decrease in
geographic spread (β = −0.355, p < 0.001), which is interpreted as a
confounded signal: easing coincided with the natural decline of the Alpha wave.
The NN-onset (August 2021 near-normal reopening) shows positive level changes in
both outcomes (cluster size β = +0.116, p = 0.016; datazones β = +0.312,
p = 0.036), the most policy-consistent genomic signal in the chapter, occurring
within the stable Delta-dominant phase.

Alpha log-odds growth was faster under the five-tier framework F5 (r = 0.0851/
day, doubling time 8.1 days, OR/week = 1.815) than under the second lockdown L2
(r = 0.0618/day, doubling time 11.2 days, OR/week = 1.542). Counterfactual
projections suggest that L2-level restrictions imposed 4, 34, or 64 days earlier
would have delayed 50% Alpha dominance by 11, 13, or 25 days respectively, but
would not have prevented establishment.

Pre-L2 Alpha meta-cluster analysis identifies 78 connected components covering
442 unique pre-L2 Alpha sequences. The largest component, AM001, contained 234
sequences (52.9% of the pre-L2 Alpha total), was dominated by Greater Glasgow
and Clyde, and was strongly enriched for the `ORF1a:L730F` mutation (85% of
AM001 vs 18.3% of non-AM001 Alpha). This concentration implies that most of the
pre-L2 Alpha burden was driven by a small number of high-amplification
components that had established across multiple regions before L2 could take
effect.

---

## Figure 1: Policy Timeline And Weekly Cluster Structure

**File:** `fig1_policy_timeline_cluster_structure`

**What it shows.** Four-panel figure covering the full study period
(July 2020 – April 2022). A thin policy-intensity colour strip at the top
labels each period by code. Panel A shows weekly median cluster size with an
IQR shaded band (Q25–Q75 across clusters within each window). Panel B shows
weekly clustering rate — (total sequences − total clusters) / total sequences
— the fraction of sequenced isolates that are secondary chain members rather
than isolated detections. Panel C shows the weekly dispersion parameter k̂ (MME
fitted to all clusters including singletons) on a log y-axis. Vertical dashed
lines in all three panels mark the four key policy transition dates (T1, L2, SL,
NN).

**Key visual patterns.** Panel A shows the wave structure of cluster size but
also its right-skewed character: the IQR band is wide during the Alpha wave
(L2/SL) and compressed toward size 1 in the post-restriction period, making the
contrast between periods much more legible than the log-median alone. The Q75
reaches 10–15 sequences during the Alpha wave peak, while Q25 stays at 1
throughout. Panel B (clustering rate) peaks near 0.88 during the Alpha-wave
periods (L3/SL), meaning 88% of sequences were secondary chain members at that
point, and falls to its lowest values in the post-restriction phase (~0.59) and
the Level 0 transition (~0.63). Panel C (k̂) is below 1.0 in every window,
confirming that no period shows Poisson-like spread; the log scale makes the
large contrast visible: k̂ drops to near 0.01 in the Omicron and long
multi-variant transition periods (OM, L21), while the Alpha-wave periods (F5,
SL) show the relatively highest k̂ values (0.27–0.36) of the selected phases.

**Suggested results paragraph.**
> Figure 1 presents three complementary weekly summaries of cluster-size
> structure across the Scottish policy timeline. Because cluster size was highly
> right-skewed in all periods, we report median cluster size with the
> interquartile range (Panel A), clustering rate (Panel B), and the
> negative-binomial dispersion parameter k̂ (Panel C) alongside the policy
> intensity colour strip. The IQR band in Panel A shows that the Q75 reached
> 10–15 sequences during the Alpha wave (L2/SL) while the Q25 remained at 1
> throughout, illustrating the persistent dominance of singletons alongside
> sporadic large clusters. Clustering rate was highest during the Alpha-wave
> periods (L3/SL: ~0.88) and lowest in the post-restriction phase (~0.59),
> tracking the epidemic wave structure closely. The dispersion parameter k̂ was
> below 1.0 in every week (persistent strong overdispersion), with the most
> extreme values (k̂ near 0.01) occurring during the long Omicron and
> post-Alpha multi-variant periods, and relatively higher values (k̂ 0.27–0.36)
> during the Alpha wave — consistent with a moderately concentrated rather than
> extreme size distribution during that period.

---

## Figure 2: Interrupted Time-Series At Three Policy Transitions (Main)

**File:** `fig2_selected_policy_transitions`

**What it shows.** A 3×2 panel of ITS analyses at the primary ±8-week window.
Rows correspond to the three selected transitions: T1-onset (October 2020),
L2→SL (April 2021), and NN-onset (August 2021). The left column shows weekly
median log cluster size; the right column shows weekly median log datazones.
Each panel plots weekly observed data points with IQR error bars (Q25–Q75 of
the within-window cluster-size or datazone distribution, log-transformed) and a
fitted OLS segmented-regression line. The β_post level-change estimate and
transition date dashed line are shown in each panel.

**Key visual patterns.** The IQR error bars on each point make the distributional
width visible within the ITS frame: in the T1-onset row (top), bars are
moderate and both outcomes show flat or slightly declining trends, with no
visible discontinuity at the transition. In the L2→SL row (middle), the
downward level shift in datazones is clear even against the wide IQR bars —
consistent with a sharp geographic contraction as the Alpha wave wound down.
The NN-onset row (bottom) shows both outcomes rising after the transition,
with the post-transition IQR bars also widening slightly, consistent with
larger and more dispersed Delta-wave clusters.

**Suggested results paragraph.**
> Figure 2 presents the ITS analyses at three selected policy transitions. At
> T1-onset (October 2020; top row), neither median log cluster size nor log
> datazones showed a statistically significant immediate level change
> (cluster size β_post = −0.080, p = 0.197), consistent with a null acute
> effect of this policy formalisation. At the second lockdown lift (L2→SL,
> April 2021; middle row), a significant downward level shift was observed for
> geographic spread (datazones β_post = −0.355, p < 0.001); this
> counter-intuitive reduction during easing is most plausibly attributed to the
> concurrent natural decline of the Alpha wave. At NN-onset (August 2021;
> bottom row), positive level changes were observed for both cluster size
> (β_post = +0.116, p = 0.016) and geographic spread (β_post = +0.312,
> p = 0.036), the most policy-consistent genomic signal in the chapter,
> occurring within the stable Delta-dominant phase. Supplementary Figures 2a
> and 2b show the same ITS framework applied to clustering rate and dispersion
> k̂, respectively, as distributional complements to the log-median outcomes.

---

## Figure 3: Alpha Emergence And Regional Expansion During F5/L2

**File:** `fig3_alpha_emergence`

**What it shows.** A two-panel summary of Alpha expansion across the pre-L2
phases. Panel A shows a weekly time series of S:N501Y (Alpha marker) and
S:A222V (B.1.177 marker) frequencies, with policy-period shading marking the
T1, F5, L2, and SL periods and vertical markers at the three Alpha phase
boundaries. Panel B is a health-board-by-phase heatmap of unique Alpha sequence
counts across the three phases: cryptic GGC-dominated chain (W016–W021),
expansion phase (W022–W024), and F5/L2 bridge (W025). Counts are unique
sequences, deduplicated across overlapping cluster windows.

**Key visual patterns.** Panel A shows S:N501Y rising from near zero in early
November 2020 to a visible fraction by late November (W022 expansion point),
then accelerating through December as B.1.177 (S:A222V) declines. By the L2
transition date (5 January 2021), Alpha has already achieved majority frequency
in the fitted and observed series. The three annotated phase boundaries
correspond to the genomic phase definitions used in the heatmap. Panel B shows
that the cryptic phase (W016–W021) captured 51 unique Alpha sequences across 7
health boards, dominated by Greater Glasgow and Clyde (GGC). The expansion phase
(W022–W024) broadened sharply to 291 sequences across 11 health boards, with
Grampian, Lanarkshire, Lothian, and Tayside contributing visibly. By the F5/L2
bridge week (W025), Alpha reached 458 sequences across 11 health boards, with
GGC still dominant but the heatmap showing a wider regional footprint.

**Suggested results paragraph.**
> Figure 3 summarises the emergence of Alpha and its B.1.177 predecessor during
> the F5–L2 policy transition. The Alpha marker S:N501Y rose from 2.1% of
> sequenced cases during the cryptic phase (W016–W021; 51 unique sequences, 7
> health boards) to 21.0% during the expansion phase (W022–W024; 291 sequences,
> 11 health boards), and reached 43.0% by the F5/L2 bridge week (W025; 458
> sequences, 11 health boards; Panel A). This trajectory indicates that Alpha had
> achieved majority frequency and multi-regional penetration before the second
> national lockdown was imposed on 5 January 2021. The health-board heatmap
> (Panel B) shows that Greater Glasgow and Clyde led the cryptic expansion before
> other boards, but that ten additional boards were represented during the F5/L2
> bridge week, consistent with a transition from localised seeding into national
> multi-region establishment ahead of L2.

---

## Figure 4: Counterfactual L2 Timing And Growth-Rate Comparison

**File:** `fig4_counterfactual_growth`

**What it shows.** A two-panel counterfactual figure. Panel A plots observed
S:N501Y frequencies alongside fitted logistic projections under the actual
F5-to-L2 transition date and three hypothetical earlier L2 timings: 4 weeks
earlier, 34 days earlier, and 64 days earlier. Panel B compares the fitted
Alpha growth rate during F5, the fitted Alpha growth rate during L2, and the
fitted B.1.177 decline rate during L2, displayed as rate estimates (per day)
with 95% confidence intervals and accompanying doubling or halving times.

**Key visual patterns.** Panel B shows that Alpha grew faster during F5 than
during L2 in the primary positive-test-weighted model: the F5 growth rate
(r = 0.0851/day, doubling time 8.1 days) is clearly to the right of the L2
estimate (r = 0.0618/day, doubling time 11.2 days), and B.1.177's decline rate
under L2 (r = −0.0655/day, halving time 10.6 days) is at the negative end of
the panel. In Panel A, the three counterfactual curves shift the projected 50%
dominance date rightward: the actual timing produces 50% Alpha frequency on the
observed date, while progressively earlier L2 imposition shifts this by 11, 13,
and 25 days respectively. By the latest scenario (L2 beginning 64 days earlier),
Alpha dominance is still achieved within the panel's time frame — the
counterfactuals delay but do not prevent establishment.

**Suggested results paragraph.**
> Figure 4 evaluates whether earlier restriction timing could plausibly have
> altered the Alpha trajectory. The primary positive-test-weighted binomial GLM
> estimated faster Alpha log-odds growth during the five-tier framework (F5:
> r = 0.0851/day, doubling time 8.1 days, OR/week = 1.815) than during the
> second lockdown (L2: r = 0.0618/day, doubling time 11.2 days, OR/week = 1.542;
> Panel B). B.1.177 declined under L2 at r = −0.0655/day (halving time 10.6 days).
> Counterfactual projections (Panel A) suggest that imposing L2-level restrictions
> 4 weeks earlier would have delayed Alpha achieving 50% frequency by 11 days;
> 34 days earlier by 13 days; and 64 days earlier by 25 days. These projections
> are descriptive rather than causal, and all three scenarios result in Alpha
> dominance within the study window, suggesting that the second lockdown may have
> slowed the variant's initial expansion but could not have prevented its
> establishment given the pre-L2 genomic evidence reviewed in Figure 3 and
> Supplementary Figures 2–3.

---

## Supplementary Figure 1: Mixing Outcomes At ITS Transitions

**File:** `supp_fig1_weekly_mixing`

**What it shows.** A 3×2 panel showing the same ITS framework as Figure 2, but
for mixing outcomes rather than cluster-structural outcomes. Rows correspond to
the three selected transitions (T1-onset, L2→SL, NN-onset); the left column
shows mean SIMD excess discordance, and the right column shows mean age excess
discordance. Each panel plots the observed weekly mixing means with fitted
segmented-regression lines and annotated β_post estimates.

**Key visual patterns.** The SIMD excess-discordance series is notably flat and
close to zero throughout all three transition windows, with narrow confidence
intervals on the β_post estimates that all straddle zero. This is consistent
with the near-zero cross-period Spearman correlation (ρ = 0.019) and with Part
1's finding that SIMD mixing within clusters is not associated with epidemic
trajectory. The age excess-discordance series shows more temporal variation —
elevated values during the pre-vaccine and Alpha periods — but also produces
β_post estimates that are weaker and less consistent than the cluster-size and
datazones outcomes in Figure 2. The supplementary location of these panels
reflects their interpretive limitations: mixing metrics are not acutely
responsive to policy-level changes.

**Suggested results paragraph.**
> Supplementary Figure 1 shows the ITS analysis repeated for demographic mixing
> outcomes. SIMD excess discordance showed no significant level change at any of
> the three transitions (all β_post 95% CIs included zero), consistent with the
> near-zero cross-period Spearman correlation (ρ = 0.019) and with Part 1's
> conclusion that socioeconomic mixing within clusters is not associated with
> the epidemic trajectory. Age excess discordance showed greater temporal
> variation — elevated during the pre-vaccine and Alpha phases — but no
> statistically significant or directionally consistent acute level change
> coinciding with any of the three policy transitions. These mixing outcomes are
> kept in supplementary materials because they do not add clear evidence of
> policy-linked demographic-mixing effects beyond the cluster-structural signals
> in Figure 2.

---

## Supplementary Figure 2a: ITS At Three Transitions — Clustering Rate

**File:** `supp_fig2a_its_clustering_rate`

**What it shows.** A 3×2 panel with the same three-transition ITS frame as
Figure 2, but with clustering rate as both columns:
- Left: **cluster-size clustering rate** = (total sequences − clusters) /
  total sequences — the fraction of sequences that are secondary chain members.
- Right: **geographic clustering rate** = (sum of datazones − clusters) /
  sum of datazones — the fraction of "datazone-slots" that are secondary (i.e.,
  the cluster spans more than one datazone).

Both outcomes are on a 0–1 proportion scale. Observed weekly values are plotted
as scatter points with the OLS segmented-regression fitted line in each panel.

**Key visual patterns.** At T1-onset (top row), both clustering rates show
no clear immediate discontinuity, consistent with the null ITS result in the
main figure. At L2→SL (middle row), the geographic clustering rate shows a
visible downward shift (β_post ≈ −0.115), sharper than the cluster-size rate
shift (β_post ≈ −0.095), reflecting that geographic spread contracted more
acutely during the Alpha tail than raw chain membership. At NN-onset (bottom
row), both rates increase, with the geographic clustering rate showing the
larger level change (β_post ≈ +0.210 vs +0.191 for cluster size), consistent
with the Delta wave increasing both the fraction of sequences in chains and
the fraction of datazones belonging to multi-datazone clusters.

**Suggested results paragraph.**
> Supplementary Figure 2a presents the ITS analysis in terms of clustering rate
> for both cluster size (left) and geographic spread (right). At T1-onset, both
> rates showed no significant immediate level change. At L2→SL, both rates
> declined after easing, with the geographic rate showing the larger shift
> (β_post ≈ −0.115 vs −0.095 for cluster size), reinforcing the interpretation
> that geographic contraction during the Alpha decline was the more acute signal.
> At NN-onset, both rates increased, confirming that near-normal reopening in the
> Delta-dominant phase was associated with a higher fraction of sequences in
> transmission chains and a higher fraction of datazones belonging to
> multi-datazone clusters.

---

## Supplementary Figure 2b: ITS At Three Transitions — Dispersion k̂

**File:** `supp_fig2b_its_dispersion`

**What it shows.** A 3×2 panel with the same three-transition ITS frame as
Figure 2, but with the MME dispersion parameter k̂ as both columns:
- Left: **k̂ for cluster size** — fitted to the per-window distribution of
  cluster sizes (all clusters including singletons).
- Right: **k̂ for geographic spread** — fitted to the per-window distribution
  of `cluster_n_datazones`.

Both axes use a log scale because k̂ varies over more than an order of magnitude
across windows. Lower values indicate greater dominance by a few large events;
higher values indicate a more Poisson-like spread across clusters.

**Key visual patterns.** k̂ is volatile at the window level (each point reflects
a single week's cluster distribution), so scatter around the fitted line is
wider than in Figure 2 or Supplementary Figure 2a. Nonetheless, directional
patterns are visible: at T1-onset (top), no clear break in either k̂ series.
At L2→SL (middle), both k̂ series show an upward shift after easing
(k̂_size β_post ≈ +0.23; k̂_geo β_post ≈ +0.23 on the log scale), meaning
the cluster-size and datazone distributions became somewhat less concentrated
after the Alpha wave declined — consistent with fewer dominant large clusters
in the SL period. At NN-onset (bottom), k̂_geo shows a downward level shift
(β_post ≈ −0.15) while k̂_size shows less consistent movement, suggesting
that near-normal reopening increased concentration in the geographic spread
distribution (a few large multi-datazone Delta clusters becoming more dominant)
even as cluster sizes were rising overall.

**Suggested results paragraph.**
> Supplementary Figure 2b shows how the per-window dispersion parameter k̂
> changed at the three ITS transitions. At T1-onset, neither k̂ for cluster
> size nor for geographic spread showed a systematic level change. At L2→SL,
> both k̂ estimates increased after easing (log-scale β_post ≈ +0.23 in each),
> indicating that the Alpha-decline period produced a more evenly distributed
> size structure — fewer dominant large clusters — relative to the pre-easing
> window. At NN-onset, k̂ for geographic spread declined (β_post ≈ −0.15),
> consistent with a small number of large geographically dispersed Delta clusters
> becoming more dominant as restrictions were removed, even as overall cluster
> volume increased.

---

## Supplementary Figure 2: Pre-L2 Alpha Meta-Cluster Amplification

**File:** `supp_fig2_meta_cluster_amplification`

**What it shows.** A four-panel summary of the rolling-window meta-cluster
network for pre-L2 Alpha sequences. Panel A shows the rank-size distribution of
all 78 inferred pre-L2 Alpha meta-clusters (unique-sequence count by rank).
Panel B shows weekly unique-sequence counts by meta-cluster, stacked,
highlighting AM001 against smaller components. Panel C shows the cumulative
pre-L2 unique-sequence contribution by meta-cluster rank. Panel D tracks
candidate meta-cluster signature mutations — specifically `ORF1a:L730F` — among
Alpha sequences after L2, with AM001 and non-AM001 series shown separately.

**Key visual patterns.** Panel A is highly right-skewed: 49 of 78 meta-clusters
contain a single sequence, while 6 contain ≥ 10 sequences. AM001 dominates with
234 sequences out of 442 total (52.9%). Panel B shows AM001 growing from early
November and accounting for the majority of the visible pre-L2 Alpha burden;
smaller components appear mainly in December. Panel C shows that the six
largest components collectively contain 312/442 (70.6%) of pre-L2 Alpha
sequences, while AM001 alone accounts for 53%. Panel D shows that `ORF1a:L730F`
is present in 85% of AM001 sequences versus 18.3% of non-AM001 Alpha sequences,
but the marker persists after L2 without being fully private to AM001 — it
declines gradually through the Alpha wave.

**Suggested results paragraph.**
> The pre-L2 Alpha expansion was not evenly distributed across inferred
> introduction events. Of 78 meta-clusters spanning 442 unique pre-L2 Alpha
> sequences, 49 were singletons and 6 contained at least 10 sequences
> (Supplementary Figure 2, Panel A). The largest component — AM001 — contained
> 234 sequences (52.9% of the pre-L2 Alpha total) and grew continuously from
> early November 2020, accounting for the majority of the visible pre-L2 Alpha
> burden in Panel B. The six largest meta-clusters together contained 312/442
> sequences (70.6%; Panel C). AM001 carried a strongly enriched `ORF1a:L730F`
> signature (85% of AM001 vs 18.3% of non-AM001 Alpha; Panel D), which persisted
> through the Alpha wave without being fully private to AM001. These findings are
> consistent with superspreading-like amplification before L2, but cannot
> identify a specific event from sequence data alone.

---

## Supplementary Figure 3: Context Of The Six Largest Alpha Meta-Clusters

**File:** `supp_fig3_meta_cluster_context`

**What it shows.** A four-panel contextual summary of the six largest pre-L2
Alpha meta-clusters, using one row per unique sequence. Panel A shows the health-
board composition of each meta-cluster. Panel B shows age-group composition.
Panel C shows SIMD quintile composition. Panel D shows grouped test-reason
composition. Each panel displays stacked proportional bars, one bar per
meta-cluster, ordered by size (AM001 leftmost).

**Key visual patterns.** Panel A shows that AM001 is dominated by Greater Glasgow
and Clyde sequences, AM003 by Grampian, AM035 by Highland, and AM024 by Borders.
This multi-board geographic spread at the meta-cluster level confirms that the
largest pre-L2 Alpha burden was both regionally structured and multi-regional by
the L2 boundary. Panel B shows that AM001 includes a substantial older component
(65+), while several smaller late-December components are more weighted toward
working-age groups. Panel C shows that AM001 has the strongest contribution from
SIMD Q1 (most deprived) among the six largest components, whereas AM003 is more
weighted toward less-deprived quintiles — patterns that are geographically
entangled and should not be read as independent social-risk effects. Panel D
shows that recorded test reasons are dominated by symptomatic testing where
available but that missing test-reason data are common, particularly for AM001
and AM024.

**Suggested results paragraph.**
> The six largest pre-L2 Alpha meta-clusters differed systematically in
> geography, age structure, deprivation composition, and test-reason
> completeness (Supplementary Figure 3). AM001 was Greater Glasgow and
> Clyde-dominated, had the largest older-age component, and showed the highest
> SIMD Q1 fraction of the six components. AM003 was Grampian-dominated and
> skewed toward less-deprived quintiles; AM035 and AM024 reflected Highland and
> Borders respectively. Age composition of smaller late-December components
> skewed toward working-age or younger groups, though counts are small and these
> patterns should be treated descriptively. The geographic and demographic
> heterogeneity across these six components reinforces the central Part 3
> caution: genomic meta-cluster structure can show policy-relevant context, but
> interpretation depends on regional timing, surveillance, and the particular
> variant phase.

---

## Supplementary Tables

### Supplementary Table 1: Lagged Intensity Correlations

**File:** `part3/tables/supp_lagged_intensity_correlations.csv`

**What it shows.** Spearman ρ between weekly policy intensity (lagged 0–4 weeks
relative to each outcome) and five cluster outcomes: median log cluster size,
median log datazones, singleton fraction, mean SIMD excess discordance, and mean
age excess discordance. This table quantifies whether the intensity–outcome
correlations strengthen when intensity is lagged to account for the delay
between policy changes and detectable changes in genomic cluster structure.

**Key findings.** Correlations with cluster size, datazones, and singleton
fraction strengthen moderately as intensity is lagged (e.g. cluster size ρ
rises from 0.741 at lag 0 to approximately 0.78 at lag 4), consistent with a
2–4 week delay between restriction changes and observable changes in genomic
cluster structure. SIMD excess discordance remains near zero at all lags
(0.019 at lag 0, declining toward −0.12 at lag 4), confirming that
socioeconomic mixing within clusters is not captured by the policy intensity
signal at any lag tested. Age excess discordance peaks at shorter lags and
weakens slightly at longer lags.

**Suggested results paragraph.**
> We tested whether intensity–outcome correlations strengthened when intensity
> was lagged by 1–4 weeks to account for the epidemiological delay between
> policy changes and observable changes in genomic cluster structure
> (Supplementary Table 1). Spearman ρ for median log cluster size increased from
> 0.741 at lag 0, and the singleton fraction correlation strengthened from −0.621
> at lag 0, consistent with a short but detectable delay between restriction
> intensity and cluster-structural changes. SIMD excess discordance remained near
> zero at all lags (range: 0.019 to approximately −0.12), providing no evidence
> that socioeconomic mixing within clusters is associated with policy intensity
> even with lagged timing.

---

### Supplementary Table 2: ITS Window Sensitivity

**File:** `part3/tables/supp_its_window_sensitivity.csv`

**What it shows.** ITS level-change estimates (β_post) and p-values for all
three transitions and four outcomes (log cluster size, log datazones, SIMD
excess discordance, age excess discordance), repeated at ±6, ±8, ±10, and
±12-week windows around each transition date.

**Key findings.** The null result at T1-onset is robust: no window produces a
significant level change for cluster size or datazones. The datazones reduction
at L2→SL is significant across all four windows, whereas the cluster-size
reduction is significant only at narrower windows, suggesting sensitivity to
window choice for that specific estimate. The datazones increase at NN-onset is
significant across all windows tested and strengthens at wider windows, supporting
a robust and growing geographic-dispersal effect in the Delta-wave context.
Mixing outcomes (SIMD, age excess discordance) produce no significant level
changes at any transition or window combination.

**Suggested results paragraph.**
> To assess sensitivity of the primary ITS results to the ±8-week window choice,
> all fits were repeated at windows of ±6, ±10, and ±12 weeks (Supplementary
> Table 2). The null result at T1-onset was consistent across all window widths
> for both structural outcomes. The geographic-dispersal increase at NN-onset
> (log datazones β_post = +0.312 at ±8 weeks) was significant at all windows
> tested and strengthened at wider windows (β range: approximately +0.27 to
> +0.38), supporting a robust effect. The datazones reduction at L2→SL was also
> robust across windows, while the cluster-size reduction at that transition was
> significant only at the primary ±8-week window, indicating greater sensitivity
> to window choice for that specific estimate. Mixing outcomes showed no
> significant level changes at any window width.

---

### Supplementary Table 3: Policy-Period Lineage Context

**File:** `part3/tables/supp_policy_lineage_context.csv`

**What it shows.** For each observed policy period, the dominant lineage group
at the period start date, its frequency in the surveillance data, and the timing
of the nearest previous and next variant overtake event. This table documents
the degree of variant-change confounding present at each policy moment, and
justifies the selection of the three ITS transition dates.

**Key findings.** The three ITS transitions were chosen to minimise variant
confounding, and this table confirms that selection: T1-onset falls 101 days
before Alpha's dominance; L2→SL falls within the stable Alpha-dominant phase
(52 days before the Alpha-to-Delta overtake); NN-onset falls 133 days before
BA.1's emergence. The most confounded policy moments in the dataset are the L2
onset (Alpha overtake just 6 days later) and the Omicron-wave restriction onset
(BA.1 overtake 21 days later) — neither of which was selected as an ITS
transition point.

**Suggested results paragraph.**
> To contextualise the ITS transition selection, each policy period start date
> was annotated with the dominant lineage at that moment and the timing of the
> nearest variant overtake event (Supplementary Table 3). The three selected
> transitions were among the least confounded by simultaneous variant change:
> T1-onset occurred 101 days before Alpha's establishment, L2→SL occurred within
> the stable Alpha-dominant phase (52 days before the Alpha-to-Delta overtake),
> and NN-onset occurred within the stable Delta-dominant phase (133 days before
> BA.1's emergence). By contrast, the start of the second lockdown (L2) coincides
> with the Alpha wave onset (overtake within 6 days), and the Omicron-wave
> restriction period begins 21 days before BA.1 dominated — both cases where ITS
> analysis would be maximally confounded by concurrent variant change, and which
> were therefore excluded from the primary ITS analysis.

---

### Supplementary Table 4: Alpha Growth Model Sensitivity

**File:** `part3/tables/alpha_growth_model_sensitivity.csv`

**What it shows.** Fitted Alpha (S:N501Y) log-odds growth rate (r, per day),
doubling time (days), and OR/week under four model specifications: (1) the
primary positive-test-weighted binomial GLM, (2) an unweighted binomial GLM,
(3) a model additionally including window-level proportion sequenced as a
covariate, and (4) a model including both weights and the proportion-sequenced
covariate. Estimates are given separately for the F5 and L2 phases.

**Key findings.** The primary positive-test-weighted estimates (F5: doubling
8.1 days; L2: doubling 11.2 days) are broadly consistent across sensitivity
specifications, with the F5/L2 rate ratio suggesting Alpha growth was 19–50%
slower under L2 depending on specification. The direction of the difference —
faster growth under F5 than L2 — is robust across all four specifications. The
magnitude varies: the unweighted binomial GLM produces the smallest F5/L2
difference, while adding proportion-sequenced as a covariate increases the
apparent slowing effect of L2.

**Suggested results paragraph.**
> Supplementary Table 4 examines the sensitivity of the estimated F5 versus L2
> Alpha growth-rate difference to model specification. The primary
> positive-test-weighted binomial GLM estimated faster growth under F5
> (r = 0.0851/day, doubling 8.1 days) than under L2 (r = 0.0618/day, doubling
> 11.2 days). This direction — faster growth during F5 than L2 — was consistent
> across all four specifications tested. The magnitude of the difference varied:
> the unweighted binomial GLM produced the smallest L2/F5 rate ratio, suggesting
> that weighting by positive-test volume matters for the magnitude but not the
> direction of the estimate. These results support treating the growth-rate
> comparison as directionally robust while acknowledging model-specification
> uncertainty in the exact magnitude.

---

### Supplementary Table 5: Per-Period Cluster Outcome Descriptives And Dispersion

**Files:** `part3/tables/period_descriptives.csv` (counts, medians, mixing means);
`part3/tables/period_clustering_dispersion.csv` (clustering rate, k̂)

**What it shows.** For each of the 16 observed policy periods: the period code,
start date, restriction intensity score, total cluster count, singleton fraction,
clustering rate, mean cluster size, variance of cluster sizes, negative-binomial
dispersion parameter k̂ fitted to all clusters (k_all) and to non-singleton
clusters only (k_non_singleton), mean SIMD excess discordance, and mean age
excess discordance. Clustering rate = (total sequences − total clusters) / total
sequences. k̂_all = x̄² / (s² − x̄) estimated from the per-period cluster-size
distribution.

**Key findings.** Clustering rate peaked at 0.88 during the Alpha-wave periods
(L3, SL), meaning 88% of sequences in those periods were secondary chain members.
It was lowest in the post-restriction phase (PR: 0.59) and the Level 0 transition
(L0: 0.63). k̂ was below 1.0 in every observed period, confirming strong
overdispersion throughout. The most extreme overdispersion occurred in the long
multi-variant Omicron and transition periods (OM k = 0.010; L21 k = 0.017),
where variance in cluster size was very large (OM variance 2,693; L21 variance
3,119), driven by mixing of many small endemic clusters with a few very large
outbreak clusters. Relatively higher k values during F5 (k = 0.265) and SL
(k = 0.362) indicate that the Alpha-wave size distribution was less dominated by
extreme large clusters than in the subsequent multi-variant and Omicron periods.

**Suggested results paragraph.**
> Supplementary Table 5 provides the per-period cluster-size distributional
> summary underlying the cross-period intensity analysis. Because cluster size was
> strongly right-skewed, we report clustering rate and the MME dispersion
> parameter k̂ alongside the conventional singleton fraction and mean cluster size.
> Clustering rate ranged from 0.59 (PR, post-restriction) to 0.88 (L3/SL, Alpha
> wave), reflecting the epidemic wave structure. k̂ was below 1.0 in every
> period, indicating persistent overdispersion; the most extreme values occurred
> in the Omicron wave (k = 0.010) and the Level 2/1 transition (k = 0.017), where
> variance exceeded 2,600 sequences². By contrast, F5 (k = 0.265) and SL
> (k = 0.362) had relatively moderate overdispersion during the Alpha wave. Mean
> SIMD excess discordance was small and variable across all periods (range
> approximately −0.26 to +0.06), consistent with the near-zero Spearman
> correlation (ρ = 0.019) between intensity and SIMD mixing.

---

## Overall Suggested Results Paragraph

Because cluster size and geographic spread were both strongly right-skewed
across all policy periods, we characterised each per-window distribution using
three complementary metrics — (i) log-median with IQR, (ii) clustering rate
= (total minus clusters) / total, and (iii) the negative-binomial dispersion
parameter k̂ (MME) — applied to both cluster size and geographic spread
(datazones). Cluster-size clustering rate ranged from 0.59 (PR) to 0.88
(L3/SL); geographic clustering rate showed a similar pattern, reflecting that
geographic spread tracked the epidemic wave structure in parallel with chain
membership. k̂ was below 1.0 in every window for both outcomes (persistent
strong overdispersion); the most extreme values occurred in the long Omicron
and multi-variant transition periods (OM k̂_size = 0.010; L21 k̂_size = 0.017),
while the Alpha-wave periods showed relatively higher k̂ (F5: 0.265; SL: 0.362),
consistent with a moderately rather than extremely concentrated size distribution
during that phase. Cross-period Spearman correlations confirmed a positive
association between restriction intensity and median log cluster size (ρ = 0.741)
and a negative association with singleton fraction (ρ = −0.621), while SIMD
excess discordance within clusters showed no meaningful correlation with
intensity (ρ = 0.019). ITS analyses at three transitions produced directionally
consistent results across all three outcome classes: null effects at T1-onset,
a downward shift in geographic spread at L2→SL driven by the Alpha decline,
and the most policy-consistent signal at NN-onset — positive level changes in
cluster-size clustering rate (+0.191), geographic clustering rate (+0.210), log
median datazones (+0.312), and a corresponding decrease in geographic k̂ (−0.15
log-units), all consistent with near-normal reopening increasing both the
fraction of sequences in chains and the dominance of a few large geographically
dispersed Delta clusters.
Lagged correlations strengthened slightly for cluster-structural outcomes,
consistent with a short epidemiological delay between policy changes and
detectable genomic cluster changes. Interrupted time-series analyses at three
transitions with minimal variant-change confounding produced a heterogeneous
picture: the T1-onset produced no acute change in cluster structure; the
L2→SL transition showed a significant reduction in geographic spread (β_post =
−0.355, p < 0.001) that is most plausibly attributed to the concurrent Alpha
wave decline; and the NN-onset showed positive level changes in both cluster
size (β_post = +0.116, p = 0.016) and geographic spread (β_post = +0.312,
p = 0.036) within the stable Delta-dominant phase, the most policy-consistent
genomic signal in the chapter. Alpha log-odds growth was faster during F5 than
L2 (doubling 8.1 days vs 11.2 days; OR/week 1.815 vs 1.542), but Alpha had
already achieved multi-regional penetration (11 health boards, 458 unique
sequences) before L2 was imposed. Counterfactual projections suggest earlier L2
imposition would have delayed Alpha dominance by 11–25 days but not prevented
it. Pre-L2 meta-cluster analysis identified AM001 — a Greater Glasgow and
Clyde-dominated high-amplification component containing 52.9% of pre-L2 Alpha
sequences and strongly enriched for `ORF1a:L730F` (85% vs 18.3%) — as the
dominant driver of the pre-L2 Alpha expansion.

## Interpretation In Context Of Parts 1 And 2

The Part 3 findings add a policy-phase and variant-emergence layer to the
socioeconomic and vaccination characterisation from Parts 1 and 2.

Part 1 established that deprivation was not associated with larger or more
dispersed clusters in a simple positive direction, and that epidemic and
surveillance conditions — not deprivation per se — were the primary structural
determinants of observed cluster scale. Part 3's cross-period correlation
analysis is consistent with this: the positive association between restriction
intensity and cluster size is substantially explained by the coincidence of
high-restriction periods with the largest epidemic wave (Alpha), not by a direct
policy effect on transmission chain size. The negligible SIMD mixing correlation
with intensity further confirms Part 1's finding that socioeconomic mixing
within clusters is not tied to the policy or epidemic trajectory.

Part 2 showed that cluster vaccination profiles tracked the national rollout
trajectory and that mixed-vaccination mixing peaked at Delta. Part 3's ITS
analysis identifies the NN-onset (Delta-wave context) as the one policy
transition with a positive and robust cluster-dispersal signal, consistent with
the Delta phase being a period of both high transmission and heterogeneous
population immunity — the context in which relaxing gathering restrictions
would plausibly have the greatest impact on observable genomic cluster spread.

The Alpha counterfactual analysis reinforces a recurring theme: the key driver
of the pre-L2 Alpha burden was a small number of high-amplification meta-clusters
that had already seeded multiple health boards before any restriction change
could take effect. This is structurally parallel to Part 1's finding that
surveillance and epidemic context dominate apparent cluster structure, and to
Part 2's finding that large clusters are strongly socioeconomically homogeneous —
large, geographically concentrated transmission events shape the measurable
genomic signal more than policy-phase or socioeconomic context alone.

Together, the three parts converge on the same methodological caution: genomic
cluster structure is shaped by overlapping epidemic, surveillance, policy, and
social processes that are difficult to disentangle with observational genomic
data alone. Results should be treated as descriptive of those overlapping
processes rather than as causal estimates of policy or deprivation effects.
