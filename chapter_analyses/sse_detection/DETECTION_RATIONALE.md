# SSE Detection Rationale

## Overview

We detect candidate superspreading events (SSEs) from a cluster-transition network in which nodes are clusters of genetically similar, temporally proximate SARS-CoV-2 sequences, and directed edges link clusters in adjacent overlapping epidemic windows that share sequences. An edge therefore encodes sequence co-membership across window-specific clustering solutions; it is a marker of lineage continuity, not a directly observed transmission event. Candidate detection is built on network topology and cluster-size behaviour, and is deliberately blind to socio-demographic information, which is reserved for the separate task of characterising detected candidates.

Detection rests on two calibrated graph-derived axes: **local burst** and **onward burden**. These axes are calibrated against permutation null expectations and are scored only where the corresponding quantity is meaningfully defined. Onward-spread shape, including dispersal across successors, is retained as a characterisation feature rather than a candidate-defining axis.

## Why Detection Is Separated From Characterisation

A recurring difficulty in network-based outbreak detection is the conflation of two questions: whether a cluster is anomalous at all, and what kind of event it represents. The second question is naturally answered with socio-demographic, geographic, and onward-shape information: a geographically concentrated, age-homogeneous cluster may suggest an institutional outbreak, while a cluster with broad onward branching may suggest wider dissemination. Folding those features into the detector, however, ties the candidate set to sampling structure and makes the detection step harder to interpret in isolation.

We therefore restrict detection to demographic-blind signals derived from cluster size and onward sequence burden, and reserve socio-demographic concentration measures, geography, age, and onward-spread shape as descriptive features attached to detected candidates. This keeps detector sensitivity dependent on sequence clustering and network structure, and ensures that which clusters are flagged is not determined by who happened to be sequenced. The descriptive features then do the work they are best suited to: sorting detected candidates into interpretable archetypes at the characterisation stage.

## Why The Detection Axes Are Distinct

The intuitive single notion of "magnitude" conflates quantities that behave differently in this data. A large local accumulation and a high onward sequence burden are related, but they are not the same signal.

The first detection axis is **local burst**: how large or novel a cluster is relative to its epidemiological context and upstream support. It combines context-adjusted excess size with the fraction of cluster sequences not observed in the union of its direct parents. Parent memberships are unioned before comparison so sequences shared with multiple parents are counted once. Parentless clusters retain a context-adjusted size score, but upstream novelty is marked not applicable rather than treating every sequence as novel. Because cluster sizes are heavy-tailed, the context adjustment is computed in log space and standardised within epidemic-window and clade strata rather than with a Poisson-style assumption that would understate the tail. Raw log cluster size and the former inflow ratio remain descriptive measures but are excluded from the burst composite to avoid repeatedly weighting size.

The second detection axis is **onward burden**: how much sequence mass follows a cluster through the transition graph. It combines two source-size-normalised quantities with distinct interpretations. Source-attributable direct burden apportions each immediate successor's new sequences among its parents according to their shares of the successor's incoming shared-sequence support. Cumulative unique future burden counts sequences absent from the source but present in any graph-reachable descendant, deduplicating sequences that recur across later windows or branches. The older unfiltered and edge-thresholded direct-burden ratios remain descriptive sensitivity measures; they are not additional components of the burden score because they largely repeat the same direct accumulation signal.

Onward-spread shape remains useful, but it answers a different question. It asks whether outgoing support is concentrated in one successor or distributed across multiple successors. That shape helps describe the phenotype of a candidate, but it is too structurally dependent on branching availability to serve as a stable primary detection axis in this analysis.

Measured on the subpopulations where both are jointly defined, local burst and onward burden are only partly correlated. That justifies treating them as separate dimensions rather than collapsing them into a single averaged score. The separation is also substantively useful: a single magnitude score would tend to rank clusters by size and could bury modest clusters that carry disproportionate onward burden.

## Why Some Signals Are Scoped To Subpopulations

The network is overwhelmingly non-propagating. Many clusters either terminate or continue linearly to a single successor, and only a minority carry appreciable onward burden. This is not a defect to be corrected but a structural property of lineage continuity under overlapping-window clustering, and it has direct consequences for scoring.

Onward burden is undefined, or trivially zero, for clusters that do not propagate. A cluster with no graph-reachable new downstream sequences has no onward burden to be anomalous about. Scoring burden across all clusters would rank a large mass of structural zeros and make the ranks of genuinely propagating clusters less meaningful.

The framework therefore records applicability explicitly. Local burst is defined for every size-eligible cluster and is scored across the full tested population. Onward burden is interpreted for clusters with positive source-attributable direct burden or positive cumulative unique future burden. Clusters outside an axis's eligible set are marked as not applicable for that axis rather than assigned a degenerate score.

## How Significance Is Assessed

For each detection axis, an observed score is formed from within-window percentile ranks of its components. A null distribution is then generated by profile permutation within adaptive contextual strata. Permuting intact multivariate profiles preserves the correlation structure among components and tests whether an assembled SSE-like profile is unusually placed in its context. This is more conservative than permuting components independently. The reported operational p-value is a seeded randomized, smoothed empirical upper-tail probability. Random allocation within exact null-score ties prevents the discrete percentile composites from producing artificial point masses at one; the corresponding conservative p-value that counts all ties in the upper tail is retained for audit. The test is one-sided because the alternative of interest is that a cluster scores higher than its context, not lower.

Single-component descriptive quantities, including onward-spread shape descriptors, are not promoted to detection axes merely because they can be ranked. In this version they are carried forward for interpretation rather than used to define candidate status.

## Interpreting One-Sided Calibration

Because detection is one-sided, the distribution of upper-tail p-values should be interpreted carefully. Under the conservative audit definition, clusters that are small or low-burden for their context can accumulate near `p = 1` because most of the null distribution, including exact ties, meets or exceeds their observed score. The operational randomized definition spreads those ties across their attainable interval. False-positive inflation would instead appear as an excess near zero among all tested nodes.

Calibration should therefore be assessed among all tested nodes and within meaningful strata, especially by cluster size and component availability. Candidate status may be used to colour the histogram but must not remove tested nodes from the calibration distribution. The randomized operational p-values should be approximately uniform under the permutation null; the conservative audit values may retain upper-tail point masses caused by discreteness.

## From Axes To Candidates And Archetypes

A cluster is treated as a candidate if it is significant on either applicable detection axis, and the set of significant axes is recorded as its signature. The framework does not require every candidate to be significant on both axes, because the axes identify different kinds of events. Both-axis candidates are retained where they occur, but single-axis candidates remain meaningful and should not be discarded merely because they express only one part of the SSE phenotype.

The relationship between local burst and onward burden provides the interpretive spine for archetype classification. Some candidates are large local bursts with limited detected onward propagation; others are modest clusters that carry disproportionate onward burden. Onward-spread shape, together with reserved socio-demographic and geographic features, supports downstream classification into outbreak-like, relay-like, and broadly disseminating archetypes without changing which clusters were detected.

## Scope And Limitations

Two limitations follow directly from the construction. First, edges record genomic co-membership across windows, not transmission. The axes measure properties of the clustering and transition network that may correlate with transmission intensity but remain affected by sampling, sequencing density, and community-detection resolution. Second, a cluster that registers as a local burst with no detected onward signal should be read as showing no detected onward propagation, not as a proven contained outbreak. Onward spread that left no cross-window shared-sequence trace, whether because of sampling gaps or clustering resolution, would be invisible to the network.

The framework therefore yields a calibrated, demographic-blind prioritisation of candidate clusters for review. Its relationship to verified superspreading should be evaluated against external epidemiological evidence rather than assumed from the graph alone.

## Statistical References

- Habiger, J. D. and Peña, E. A. (2011). Randomised p-values and
  nonparametric procedures in multiple testing. *Journal of Nonparametric
  Statistics*, 23(3), 583--604.
  <https://doi.org/10.1080/10485252.2010.482154>
- Hemerik, J. and Goeman, J. (2018). Exact testing with random permutations.
  *TEST*, 27, 811--825. <https://doi.org/10.1007/s11749-017-0571-1>
- Phipson, B. and Smyth, G. K. (2010). Permutation p-values should never be
  zero: calculating exact p-values when permutations are randomly drawn.
  *Statistical Applications in Genetics and Molecular Biology*, 9(1), Article
  39. <https://doi.org/10.2202/1544-6115.1585>
