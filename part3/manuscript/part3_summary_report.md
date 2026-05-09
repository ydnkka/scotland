# Part 3 Summary Report

## Core Interpretation

Scottish COVID-19 policy periods were strongly entangled with variant waves,
immunity, testing, sequencing, and calendar time.  Part 3 therefore treats
policy as descriptive epidemic context, not as a source of causal policy
effects.  Selected periods are still useful because they show how restrictions,
variant advantage, and genomic cluster structure interacted at key moments.

## Main Findings

### Autumn 2020, P3 -> T1

The `P3 -> T1` transition is the cleanest tightening period without a major
variant replacement.  In the primary plus/minus 8-week ITS, the immediate
post-T1 terms were small and not statistically clear for median log cluster
size (`post = -0.080`, p = 0.20) or median log datazones (`post = -0.027`,
p = 0.87).  The post-transition slopes drift downward, but this is best read
as a weak/null tightening signal rather than causal evidence.

### Winter 2020/21, F5 -> L2 and Alpha

The Alpha case study is the centrepiece.  Pango-defined Alpha increased from
51 sequences in W016-W021 to 291 in W022-W024, then 458 in W025.  By W025,
Alpha represented 43.0% of primary good-QC sequences in that phase summary and
was present across 11 health boards and 29 local authorities.

The mutation trajectory shows `S:N501Y` rising from 17.7% in W022
(2020-12-08) to 47.5% in W025 (2020-12-29), then 60.4% in W026
(2021-01-05), the L2 start window.  Positive-test-weighted growth models
estimate Alpha `S:N501Y` growth of 0.596 log-odds per week during F5 and 0.433
during L2.  The B.1.177 comparator `S:A222V` declines during L2 at -0.459
log-odds per week.

The counterfactual timing scenarios suggest delays rather than prevention:
switching to the fitted L2 growth rate on 2020-12-08 projects 50% `S:N501Y`
around 2021-01-09; switching on 2020-12-02 projects around 2021-01-11; switching
from the F5 start projects around 2021-01-23.  These are descriptive fitted
scenarios, not causal policy estimates.

The supplementary pre-L2 meta-cluster analysis decomposes this rise into
connected components of Alpha-containing rolling-window cluster calls.  Among
442 unique pre-L2 Alpha sequences, the graph contains 78 meta-clusters, but the
distribution is highly skewed.  Six components contain at least 10 sequences
and together account for 312/442 pre-L2 Alpha sequences (70.6%).  AM001 alone
contains 234/442 sequences (52.9%).

AM001 is a high-amplification, Greater Glasgow and Clyde-dominated component
that is strongly enriched for `ORF1a:L730F` (199/234 AM001 sequences, compared
with 38/208 non-AM001 pre-L2 Alpha sequences).  The marker is not private to
AM001, so it should be treated as an enriched signature rather than a unique
lineage definition.  The AM001 signature remains common through the later Alpha
wave, supporting the interpretation that one large pre-L2 component made a
substantial contribution to Alpha establishment.

The top-six contextual analysis uses one row per unique sequence, avoiding
overlapping-window inflation.  AM001 is older and more deprived-skewed than
several smaller components: 85/234 sequences are aged 40-64, 79/234 are aged
65+, and 86/234 are in SIMD Q1.  AM003 is mainly Grampian and weighted toward
less-deprived quintiles; AM035 is mainly Highland; AM024 is Borders.  Recorded
test reasons are dominated by symptomatic testing where present, but missing
test-reason data are common, especially in AM001 and AM024.

These findings refine the F5/L2 interpretation.  Alpha establishment before L2
was not just a smooth population-level curve; it was disproportionately built
from a small number of high-amplification genomic components.  This is
consistent with superspreading-like amplification and regional seeding before
L2, but it does not identify a specific event from sequence data alone.  It
also supports the counterfactual conclusion that earlier L2-level restrictions
may plausibly have delayed dominance and reduced burden, while complete
prevention of Alpha establishment was unlikely once these components were
already expanding.

### Spring 2021, L2 -> SL

After the move from L2 to stay-local/Level 3, median log cluster size continues
to decline in the primary ITS (`post_t = -0.072` per week, p = 0.001).  This is
a warning example: apparent policy associations can run opposite to expectation
when wave phase dominates.  Easing coincided with the Alpha tail.

### Summer 2021, L0 -> NN

The `L0 -> NN` transition gives the clearest policy-consistent genomic signal.
Median log cluster size increased immediately (`post = 0.116`, p = 0.016) and
continued rising after the transition (`post_t = 0.051` per week, p < 0.001).
Median log datazones also increased (`post = 0.312`, p = 0.036; `post_t =
0.048` per week, p = 0.038).  Mixing metrics were weak or null, so the most
interpretable signal is geographic and cluster-size structure during Delta.

## Whole-Epidemic Context

Weekly policy intensity is correlated with several cluster outcomes, including
median log cluster size (Spearman rho = 0.74) and median log datazones
(rho = 0.58).  These correlations are explicitly confounded by variant phase,
testing, sequencing, immunity, and calendar time.  SIMD excess-discordance is
near-null in the whole-epidemic correlation table (rho = 0.02).

## Caveats

This analysis does not identify causal policy effects.  It does not identify
specific superspreading events from sequence data alone.  The Alpha
counterfactuals should be described as timing projections under fitted marker
growth rates, not estimates of what restrictions would have caused.  The
pre-L2 Alpha meta-clusters should similarly be described as inferred
high-amplification components, not confirmed introductions or events.
