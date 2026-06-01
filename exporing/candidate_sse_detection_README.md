# Candidate SSE detection dummy example

Files included:

- `candidate_sse_detection_dummy.py`: standalone Python script.
- `candidate_sse_detection_dummy.ipynb`: Jupyter notebook version.

The example simulates a sparse temporal cluster-transition graph and ranks candidate SSE-like nodes using:

- sampling-adjusted excess cluster size;
- outgoing transition weight;
- weighted out-degree;
- outgoing edge entropy;
- downstream burden;
- temporal compactness;
- socio-geodemographic coherence;
- stratified permutation null scores.

## Run the script

```bash
pip install pandas numpy networkx matplotlib
python candidate_sse_detection_dummy.py
```

The script writes outputs to `sse_dummy_outputs/`.

## Main outputs

- `dummy_cluster_nodes.csv`
- `dummy_transition_edges.csv`
- `node_features.csv`
- `ranked_candidate_sse_nodes.csv`
- `transition_graph.png`
- `ranked_candidate_scores.png`

## How to adapt to real data

Replace the dummy `nodes` and `edges` tables.

Your node table should include at least:

```text
node
window
cluster_size
region / health board
lineage
sampling intensity or sequencing denominator
temporal span
socio-geodemographic coherence variables
```

Your edge table should include:

```text
source
target
source_window
target_window
weight
```

The most important real-data improvement is to replace the toy `expected_size_sampling`
calculation with a proper sampling model based on case counts, sequencing proportions,
region, lineage and epidemic window.
