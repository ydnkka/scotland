# Superspreading Signature Categories

These labels describe node-level signatures in the temporal cluster-transition
network. They are heuristic categories for review, not confirmed exposure
events. Each candidate label has the form:

`<role>__<onward_dynamic>`

Rows labelled `not_sse_like` did not pass the amplification screen.

## Role Labels

| Role | Interpretation | Typical evidence |
|---|---|---|
| `putative_birth` | A newly observed or near-new node that may represent the start of a visible burst. | Little or no incoming overlap, high novelty, high amplification. |
| `relay_amplifier` | A node with inherited upstream overlap that appears to amplify before continuing onward. | Incoming and outgoing overlap present, high net amplification. |
| `merged_relay` | A node where multiple upstream sources converge and then continue onward. | Multiple incoming edges, outgoing overlap present. |
| `terminal_sink` | A node with upstream overlap but no observed onward continuity. | Incoming overlap present, zero outgoing strength. |
| `isolated_burst` | A large or unusual node with no observed network continuity on either side. | Zero incoming and outgoing strength; possible sampling gap or event outside observed windows. |
| `unclear_origin` | Candidate passes the amplification screen, but the network role is not cleanly identifiable. | Mixed or borderline lifecycle evidence. |

## Onward Dynamic Labels

| Onward dynamic | Interpretation | Typical evidence |
|---|---|---|
| `no_observed_onward_spread` | No successor node is observed in the adjacent-window graph. | Zero `out_strength`; interpret carefully if right-censored. |
| `contained_burst` | Candidate is large or unusual but has weak onward continuity. | Low `out_degree`, low downstream expansion, or terminal lifecycle. |
| `single_dominant_chain` | Onward spread is mostly concentrated into one successor. | Low downstream entropy or high `dominant_successor_frac`. |
| `multi_branch_seeder` | Candidate seeds multiple successor clusters without necessarily producing high total onward burden. | Multiple outgoing edges and high downstream entropy. |
| `multi_branch_expander` | Candidate seeds multiple substantial successor branches. | High downstream entropy plus high downstream expansion. |
| `diffuse_spatial_broadcaster` | Candidate branches onward and is spatially or demographically broad. | High downstream entropy, high downstream expansion, and broad spatial or mixing evidence. |
| `high_volume_onward_spread` | Candidate has high outgoing burden but not enough entropy evidence for a branching label. | High `out_strength` with multiple successors, but downstream evenness is ambiguous. |
| `weak_or_ambiguous_onward_spread` | Candidate passes the amplification screen, but onward evidence is weak or conflicting. | Some onward continuity, but no clear concentration, branching, or broadcaster signature. |

## Common Composite Categories

| Category | Short description | Typical signature |
|---|---|---|
| `putative_birth__contained_burst` | New large cluster with limited observed onward spread. | High amplification, high novelty, weak downstream expansion. |
| `putative_birth__single_dominant_chain` | New burst that mainly continues through one successor lineage. | High novelty, one dominant outgoing successor. |
| `putative_birth__multi_branch_seeder` | New burst that seeds several successor clusters. | High novelty, multiple outgoing edges, high downstream entropy. |
| `putative_birth__multi_branch_expander` | New burst that seeds multiple substantial onward branches. | High novelty, high `out_strength`, high downstream entropy and expansion. |
| `putative_birth__diffuse_spatial_broadcaster` | New burst followed by broad onward dissemination. | High novelty, high downstream entropy, broad spatial or mixing evidence. |
| `relay_amplifier__contained_burst` | Secondary amplification with limited onward spread. | Incoming overlap, high net amplification, weak downstream expansion. |
| `relay_amplifier__single_dominant_chain` | Secondary amplification that mainly feeds one continuing chain. | Incoming overlap, high net amplification, one dominant successor. |
| `relay_amplifier__multi_branch_seeder` | Secondary amplification that splits into several successors. | Incoming overlap, high amplification, high downstream entropy. |
| `relay_amplifier__multi_branch_expander` | High-impact secondary amplifier with substantial branching. | Incoming overlap, high `out_strength`, high downstream expansion. |
| `relay_amplifier__diffuse_spatial_broadcaster` | Secondary amplifier with broad onward dissemination. | High net amplification, high downstream entropy, broad spatial or mixing evidence. |
| `merged_relay__multi_branch_expander` | Multiple incoming sources merge and then expand onward. | High `in_degree`, high `out_degree`, high `out_strength`. |
| `merged_relay__single_dominant_chain` | Merged node continues mainly through one dominant successor. | High `in_degree`, low downstream entropy, high `dominant_successor_frac`. |
| `isolated_burst__no_observed_onward_spread` | Large or unusual node with no observed graph continuity. | No incoming or outgoing overlap; check sampling and temporal boundaries. |
| `terminal_sink__no_observed_onward_spread` | Upstream cluster continuity ends at this node. | Nonzero `in_strength`, zero `out_strength`; often right-censored or dying out. |
| `unclear_origin__weak_or_ambiguous_onward_spread` | Amplification signal exists but origin and onward dynamics are unclear. | Candidate score is high, but lifecycle and downstream features conflict. |

## Family-Level Grouping

| Family | Includes | Interpretation |
|---|---|---|
| Contained SSE | `*__contained_burst`, `*__no_observed_onward_spread` | Large short-lived burst with limited detected onward transmission. |
| Chain SSE | `*__single_dominant_chain` | Event seeded a persistent but narrow transmission chain. |
| Branching SSE | `*__multi_branch_seeder`, `*__multi_branch_expander` | Event seeded multiple descendant clusters. |
| Broadcaster SSE | `*__diffuse_spatial_broadcaster` | Event associated with broad spatial or demographic dissemination. |
| Relay SSE | `relay_amplifier__*`, `merged_relay__*` | Not necessarily the original introduction, but important secondary amplification. |
| Ambiguous SSE | `unclear_origin__*`, `*__weak_or_ambiguous_onward_spread`, `*__high_volume_onward_spread` | Candidate merits review, but evidence does not support a narrower label. |

## Censoring Notes

| Note | Meaning |
|---|---|
| `not_censored` | Node is not at the full-dataset or VOC-epoch boundary. |
| `left_censored_origin_uncertain` | Origin may have occurred before the observed dataset or VOC epoch. |
| `right_censored_onward_uncertain` | Onward spread may continue beyond the observed dataset or VOC epoch. |
| `both_left_and_right_censored` | Both origin and onward interpretation are boundary-limited. |
