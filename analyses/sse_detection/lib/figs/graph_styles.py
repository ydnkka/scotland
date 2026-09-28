"""Shared transition-graph role names and colours."""

CLUSTER_ROLE_GROUPS = {
    "isolated": "Isolated",
    "single_outgoing_source": "Source",
    "source_branching": "Branching source",
    "single_incoming_sink": "Sink",
    "merging_sink": "Merging sink",
    "linear_continuation": "Linear continuation",
    "internal_branching": "Internal branching",
    "internal_merging": "Internal merging",
    "merge_and_branch": "Merge and branch",
    "other": "Other",
}

ROLE_COLORS = {
    "Isolated": "#4c78a8",
    "Source": "#f58518",
    "Branching source": "#e45756",
    "Sink": "#54a24b",
    "Merging sink": "#e15ebc",
    "Linear continuation": "#8e6c8a",
    "Internal branching": "#ff9da6",
    "Internal merging": "#9d755d",
    "Merge and branch": "#72b7b2",
    "Other": "#bab0ac",
}
