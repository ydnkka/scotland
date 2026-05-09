# Part 3 Results And Figures Description

## Figure 1. Policy timeline and weekly cluster structure

Figure 1 places weekly median log cluster size against the full Scottish policy
timeline.  The figure is intended as context rather than evidence of causal
policy effects.  Cluster structure varies across periods in ways that are
visibly entangled with variant waves and surveillance intensity.

## Figure 2. Selected policy transitions

Figure 2 shows the primary plus/minus 8-week segmented OLS summaries for three
selected policy transitions.

For `P3 -> T1`, there is no clear immediate tightening-associated break in
median log cluster size or log datazones.  The fitted slopes drift downward
after T1, but this is not interpreted as a causal restriction effect.

For `L2 -> SL`, cluster structure declines after easing.  This is the warning
example: apparent associations run opposite to the simple policy expectation
because easing coincides with the Alpha tail.

For `L0 -> NN`, both median log cluster size and median log datazones increase
around near-normal reopening.  This is the most policy-consistent genomic
signal in the chapter, although still descriptive and Delta-phase dependent.

## Figure 3. Alpha emergence during F5/L2

Figure 3 shows Alpha marker expansion (`S:N501Y`) against the B.1.177 marker
(`S:A222V`) and displays Alpha spread across health boards.  Alpha expanded
during `F5`, was a large fraction of marker-positive sequences by W025, and was
already established at the start of `L2`.

The figure supports a timing interpretation: L2 may have slowed Alpha growth,
but the variant had already bridged multiple regions before L2 could fully
take effect.

## Figure 4. Counterfactual L2 timing and growth rates

Figure 4 combines observed `S:N501Y` frequencies, fitted trajectories, and
earlier-L2 timing scenarios.  The primary positive-test-weighted model estimates
faster Alpha log-odds growth during `F5` than during `L2`, while the `S:A222V`
comparator declines under `L2`.

The counterfactual curves are descriptive fitted scenarios.  They suggest that
earlier L2-level restrictions could have delayed Alpha dominance and reduced
B.1.177 burden, but they do not support a claim that Alpha establishment would
have been prevented.

## Supplementary Figure 1. Mixing outcomes

The supplementary ITS figure shows SIMD and age excess-discordance outcomes for
the same selected transitions.  These mixing outcomes are weaker and less
consistent than cluster size or geographic spread, so they are kept out of the
main transition figure.

## Supplementary Figure 2. Pre-L2 Alpha meta-cluster amplification

Supplementary Figure 2 extends Figure 3 by decomposing the early Alpha rise into
connected components of rolling-window genomic clusters.  Nodes in the
underlying analysis are Alpha-containing cluster/window assignments; adjacent
windows are linked when they share a sequence, and connected components are
treated as Alpha meta-clusters.

Panel A shows the rank-size distribution of all 78 inferred pre-L2 Alpha
meta-clusters.  The distribution is highly skewed: 49 components contain a
single sequence, whereas six components contain at least 10 sequences.  AM001
dominates the pre-L2 Alpha population with 234 of 442 unique sequences.

Panel B shows weekly unique-sequence counts before L2.  AM001 grows from early
November and accounts for most of the visible pre-L2 Alpha burden.  Smaller
components appear mainly in December, consistent with multiple introductions or
secondary high-amplification expansions emerging shortly before L2.

Panel C shows cumulative pre-L2 expansion.  The six largest components contain
312/442 pre-L2 Alpha sequences, while AM001 alone contains 234/442.  This
supports the interpretation that a small number of high-amplification
components had already set up much of the Alpha wave before L2 could take
effect.

Panel D tracks candidate meta-cluster signature mutations after L2.  AM001 is
strongly enriched for `ORF1a:L730F`, which remains common among Alpha sequences
through the subsequent wave, but the marker is not private to AM001.  The
trajectory therefore supports AM001's large contribution without reducing the
wave to a single genetically unique event.

Suggested interpretation:

> The pre-L2 Alpha expansion was not evenly distributed across inferred
> introductions.  Instead, most early Alpha burden was concentrated in a small
> set of connected genomic meta-clusters.  AM001, a Greater Glasgow and
> Clyde-dominated high-amplification component, contained 52.9% of unique
> pre-L2 Alpha sequences and carried an enriched `ORF1a:L730F` signature that
> persisted through the Alpha wave.  These findings are consistent with
> superspreading-like amplification before L2, but do not identify a specific
> event from sequence data alone.

## Supplementary Figure 3. Context of the six largest Alpha meta-clusters

Supplementary Figure 3 describes the six largest pre-L2 Alpha meta-clusters
using one row per unique sequence, not overlapping rolling-window rows.

Panel A shows geography.  AM001 is dominated by Greater Glasgow and Clyde
sequences, AM003 by Grampian, AM035 by Highland, and AM024 by Borders.  This
supports the interpretation that the largest pre-L2 Alpha burden was both
regionally structured and multi-regional by the L2 boundary.

Panel B shows age composition.  AM001 includes a substantial older component
(65+), while several smaller late-December components skew toward working-age
or younger groups.  Counts are small outside AM001, so these patterns should be
treated descriptively.

Panel C shows SIMD quintile composition.  AM001 has the strongest contribution
from the most-deprived quintile among the six largest components, whereas
AM003 is more weighted toward less-deprived quintiles.  These differences are
geographically entangled and should not be interpreted as independent social
risk effects.

Panel D shows grouped test reasons.  Recorded reasons are dominated by
symptomatic testing where available, but missing test-reason data are common,
especially for AM001 and AM024.  Testing composition therefore provides
surveillance context rather than a complete ascertainment model.

Suggested interpretation:

> The six largest pre-L2 Alpha meta-clusters differed in geography,
> deprivation, age structure, and test-reason completeness.  AM001 was large,
> Greater Glasgow and Clyde-dominated, older, and more weighted toward SIMD Q1
> than several smaller components.  These contextual differences reinforce the
> central Part 3 caution: genomic clusters can show policy-relevant structure,
> but interpretation depends on variant phase, surveillance, and regional
> timing.
