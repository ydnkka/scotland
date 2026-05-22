# Superspreading Signature Categories

These labels describe node-level signatures in the temporal cluster-transition
network. They are heuristic categories for review, not confirmed exposure
events.

`sse_category` is now the compact epidemiological category used for headline
interpretation and plotting. The previous graph-diagnostic composite label is
preserved in `sse_graph_category` with the form:

`<role>__<onward_dynamic>`

Rows labelled `not_sse_like` did not pass the amplification screen.

## Windowing and Burden Screen

Node-clusters are built in overlapping 3-week windows advanced in 1-week
increments in the upstream clustering dataset. This detection notebook keeps
every other original window after construction, so retained node windows are
2 weeks apart but still arise from overlapping 3-week sequence windows.

Edges connect adjacent retained windows when sequences are shared across
node-clusters. Because the underlying windows overlap, `in_strength`,
`out_strength`, and adjacency edges partly measure window carry-over and
sampling continuity rather than pure onward transmission. They are therefore
used as SSE-like signatures for review, not as confirmed epidemiological
transmission events.

The candidate screen requires at least 6 sampled sequences in a node-cluster.
This modest absolute burden floor is a deliberate scientific choice: clusters
of only 3-5 sampled sequences can be high-novelty introductions in molecular
data, but they are too small to label as candidate SSE-like nodes by default.
Their novelty and amplification metrics remain available for sensitivity
review.

## Epidemiological Categories

Categories are mutually exclusive and assigned in priority order after the
candidate screen:

| `sse_category` | Interpretation | Typical graph evidence |
| --- | --- | --- |
| `not_sse_like` | Node did not pass the amplification screen. | Candidate screen is false. |
| `mixed_population_dissemination` | Candidate associated with onward spread into a socio-geodemographically diverse population. | `sse_onward_dynamic == diverse_population_broadcaster`. |
| `putative_introduction_burst` | Newly observed high-novelty candidate with observed onward continuity. | `sse_role == putative_birth`; strict zero incoming degree, high novelty, and at least one outgoing edge. |
| `secondary_relay_amplification` | Candidate is not necessarily the original introduction, but appears to be an important secondary amplifier. | `sse_role` is `relay_amplifier` or `merged_relay`, unless mixed-population dissemination applies. |
| `diffuse_branching_transmission` | Candidate seeds multiple descendant clusters. | `sse_onward_dynamic` is `multi_branch_seeder` or `multi_branch_expander`. |
| `focused_branching_transmission` | Candidate seeds multiple successors, but observed onward continuity is concentrated in one branch. | `sse_onward_dynamic == dominant_branch`. |
| `sustained_single_chain` | Candidate feeds a persistent but narrow observed continuity chain. | `sse_onward_dynamic == single_successor_chain`. |
| `contained_local_burst` | Large short-lived burst with limited detected onward transmission. | `sse_onward_dynamic` is `contained_burst` or `no_observed_onward_spread`. |
| `high_volume_onward_transmission` | Candidate has high outgoing burden but not enough entropy evidence for a branching label. | `sse_onward_dynamic == high_volume_onward_spread`. |
| `ambiguous_amplification_signal` | Candidate merits review, but onward evidence is weak or conflicting. | Remaining candidate signatures, usually `weak_or_ambiguous_onward_spread`. |

## Graph-Diagnostic Role Labels

`sse_role` remains available for diagnostic interpretation.

| Role | Interpretation | Typical evidence |
| --- | --- | --- |
| `putative_birth` | A newly observed node that may represent the start of a visible burst. | No incoming graph edge, high novelty, high amplification, and observed onward continuity. |
| `relay_amplifier` | A node with inherited upstream overlap that appears to amplify before continuing onward. | Incoming and outgoing overlap present, high net amplification. |
| `merged_relay` | A node where multiple upstream sources converge and then continue onward. | Multiple incoming edges, outgoing overlap present. |
| `terminal_sink` | A node with upstream overlap but no observed onward continuity. | Incoming overlap present, zero outgoing strength. |
| `isolated_burst` | A large or unusual node with no observed network continuity on either side. | Zero incoming and outgoing strength; possible sampling gap or event outside observed windows. |
| `unclear_origin` | Candidate passes the amplification screen, but the network role is not cleanly identifiable. | Mixed or borderline lifecycle evidence. |

## Graph-Diagnostic Onward Dynamic Labels

`sse_onward_dynamic` remains available for diagnostic interpretation.

| Onward dynamic | Interpretation | Typical evidence |
| --- | --- | --- |
| `no_observed_onward_spread` | No successor node is observed in the adjacent-window graph. | Zero `out_strength`; interpret carefully if right-censored. |
| `contained_burst` | Candidate is large or unusual but has weak onward continuity. | Terminal node with incoming burden, or low downstream expansion ranked among nodes with observed onward spread. |
| `single_successor_chain` | Candidate continues through one observed successor. | Exactly one outgoing successor; no branching choice is observed. |
| `dominant_branch` | Candidate has multiple successors, but one branch dominates. | At least two outgoing successors plus low downstream entropy or high `dominant_successor_frac`. |
| `multi_branch_seeder` | Candidate seeds multiple successor clusters without necessarily producing high total onward burden. | Multiple outgoing edges and high downstream entropy. |
| `multi_branch_expander` | Candidate seeds multiple substantial successor branches. | High downstream entropy plus high downstream expansion. |
| `diverse_population_broadcaster` | Candidate branches onward into a socio-geodemographically diverse population. | High downstream entropy, high downstream expansion, and high observed-normalised `mixing_score`. |
| `high_volume_onward_spread` | Candidate has high outgoing burden but not enough entropy evidence for a branching label. | High `out_strength` with multiple successors, but downstream evenness is ambiguous. |
| `weak_or_ambiguous_onward_spread` | Candidate passes the amplification screen, but onward evidence is weak or conflicting. | Some onward continuity, but no clear concentration, branching, or broadcaster signature. |

## Censoring Notes

| Note | Meaning |
| --- | --- |
| `not_censored` | Node is not at the full-dataset or VOC-epoch boundary. |
| `left_censored_origin_uncertain` | Origin may have occurred before the observed dataset or VOC epoch. |
| `right_censored_onward_uncertain` | Onward spread may continue beyond the observed dataset or VOC epoch. |
| `both_left_and_right_censored` | Both origin and onward interpretation are boundary-limited. |
