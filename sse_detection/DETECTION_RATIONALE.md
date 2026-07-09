# SSE Detection Rationale

## Overview

We detect candidate superspreading events (SSEs) from a cluster-transition network in which nodes are clusters of genetically similar, temporally proximate SARS-CoV-2 sequences, and directed edges link clusters in adjacent overlapping epidemic windows that share sequences. An edge therefore encodes sequence co-membership across window-specific clustering solutions; it is a marker of lineage continuity, not a directly observed transmission event. Candidate detection is built on network topology and cluster-size behaviour, and is deliberately blind to socio-demographic information, which is reserved for the separate task of characterising detected candidates.

Detection now rests on two calibrated graph-derived axes: **local burst** and **onward burden**. These axes are calibrated against permutation null expectations and are scored only where the corresponding quantity is meaningfully defined. Onward-spread shape, including dispersal across successors, is retained as a characterisation feature rather than a candidate-defining axis.

## Why Detection Is Separated From Characterisation

A recurring difficulty in network-based outbreak detection is the conflation of two questions: whether a cluster is anomalous at all, and what kind of event it represents. The second question is naturally answered with socio-demographic, geographic, and onward-shape information: a geographically concentrated, age-homogeneous cluster may suggest an institutional outbreak, while a cluster with broad onward branching may suggest wider dissemination. Folding those features into the detector, however, ties the candidate set to sampling structure and makes the detection step harder to interpret in isolation.

We therefore restrict detection to demographic-blind signals derived from cluster size and onward sequence burden, and reserve socio-demographic concentration measures, geography, age, and onward-spread shape as descriptive features attached to detected candidates. This keeps detector sensitivity dependent on sequence clustering and network structure, and ensures that which clusters are flagged is not determined by who happened to be sequenced. The descriptive features then do the work they are best suited to: sorting detected candidates into interpretable archetypes at the characterisation stage.

## Why The Detection Axes Are Distinct

The intuitive single notion of "magnitude" conflates quantities that behave differently in this data. A large local accumulation and a high onward sequence burden are related, but they are not the same signal.

The first detection axis is **local burst**: how large or novel a cluster is relative to its epidemiological context and upstream support. This is captured by context-adjusted excess size, excess over upstream inflow, and the fraction of the cluster not explained by incoming overlap. Because cluster sizes are heavy-tailed, the context adjustment is computed in log space and standardised within epidemic-window and clade strata rather than with a Poisson-style assumption that would understate the tail. Local burst measures the intensity of local accumulation, which is often the first question in public-health review.

The second detection axis is **onward burden**: how much sequence mass a cluster carries into adjacent-window successors. The cleanest burden measures are downstream quantities scaled against source size or adjusted for overlap with the source. Raw downstream counts can mechanically track cluster size and should therefore be interpreted cautiously unless they are ratio-normalised or calibrated alongside local-burst terms.

Onward-spread shape remains useful, but it answers a different question. It asks whether outgoing support is concentrated in one successor or distributed across multiple successors. That shape helps describe the phenotype of a candidate, but it is too structurally dependent on branching availability to serve as a stable primary detection axis in this analysis.

Measured on the subpopulations where both are jointly defined, local burst and onward burden are only partly correlated. That justifies treating them as separate dimensions rather than collapsing them into a single averaged score. The separation is also substantively useful: a single magnitude score would tend to rank clusters by size and could bury modest clusters that carry disproportionate onward burden.

## Why Some Signals Are Scoped To Subpopulations

The network is overwhelmingly non-propagating. Many clusters either terminate or continue linearly to a single successor, and only a minority carry appreciable onward burden. This is not a defect to be corrected but a structural property of lineage continuity under overlapping-window clustering, and it has direct consequences for scoring.

Onward burden is undefined, or trivially zero, for clusters that do not propagate. A cluster with no onward edges has no onward burden to be anomalous about. Scoring burden across all clusters would rank a large mass of structural zeros and make the ranks of genuinely propagating clusters less meaningful.

The framework therefore records applicability explicitly. Local burst is defined for every size-eligible cluster and is scored across the full tested population. Onward burden is interpreted for clusters with observed onward support. Clusters outside an axis's eligible set are marked as not applicable for that axis rather than assigned a degenerate score.

## How Significance Is Assessed

For each detection axis, an observed score is formed from within-window percentile ranks of its components. A null distribution is then generated by profile permutation within adaptive contextual strata. Permuting intact multivariate profiles preserves the correlation structure among components and tests whether an assembled SSE-like profile is unusually placed in its context. This is more conservative than permuting components independently. The reported p-value is a smoothed empirical upper-tail probability, and the test is one-sided because the alternative of interest is that a cluster scores higher than its context, not lower.

Single-component descriptive quantities, including onward-spread shape descriptors, are not promoted to detection axes merely because they can be ranked. In this version they are carried forward for interpretation rather than used to define candidate status.

## Interpreting One-Sided Calibration

Because detection is one-sided, the background distribution of upper-tail p-values should be interpreted carefully. Clusters that are small or low-burden for their context naturally fall near p = 1, because most of the null distribution exceeds their observed score. This is conservative behaviour, not false-positive inflation. False-positive inflation would instead appear as an excess of background p-values near zero.

Calibration should therefore be assessed within meaningful strata, especially by cluster size and component availability. Among clusters large enough to express the relevant signal, background upper-tail p-values should be approximately uniform after excluding flagged candidates. At the size floor, a concentration near p = 1 is expected and should be reported as conservative under-calling of the least anomalous clusters rather than as evidence of lower-tail inflation.

## From Axes To Candidates And Archetypes

A cluster is treated as a candidate if it is significant on either applicable detection axis, and the set of significant axes is recorded as its signature. The framework does not require every candidate to be significant on both axes, because the axes identify different kinds of events. Both-axis candidates are retained where they occur, but single-axis candidates remain meaningful and should not be discarded merely because they express only one part of the SSE phenotype.

The relationship between local burst and onward burden provides the interpretive spine for archetype classification. Some candidates are large local bursts with limited detected onward propagation; others are modest clusters that carry disproportionate onward burden. Onward-spread shape, together with reserved socio-demographic and geographic features, supports downstream classification into outbreak-like, relay-like, and broadly disseminating archetypes without changing which clusters were detected.

## Scope And Limitations

Two limitations follow directly from the construction. First, edges record genomic co-membership across windows, not transmission. The axes measure properties of the clustering and transition network that may correlate with transmission intensity but remain affected by sampling, sequencing density, and community-detection resolution. Second, a cluster that registers as a local burst with no detected onward signal should be read as showing no detected onward propagation, not as a proven contained outbreak. Onward spread that left no cross-window shared-sequence trace, whether because of sampling gaps or clustering resolution, would be invisible to the network.

The framework therefore yields a calibrated, demographic-blind prioritisation of candidate clusters for review. Its relationship to verified superspreading should be evaluated against external epidemiological evidence rather than assumed from the graph alone.
