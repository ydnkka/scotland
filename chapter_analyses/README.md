# Thesis Chapter Analyses

This namespace contains the sequential surveillance, genomic-network, and superspreading-compatible signal analyses used by the thesis.

## Packages

- `surveillance/`: descriptive sequence surveillance, Scottish policy periods, OxCGRT index validation, lineage replacement, and sequencing coverage.
- `genomic_networks/`: observed cohort, window-specific EpiLink clusters, compatibility-network mixing and assortativity, topology, and sensitivity analyses.
- `sse_detection/`: temporal cluster-transition graph, burst and onward-burden candidate detection, figures, composition tables, and Bayesian characterisation.

Each package contains its own README, technical reference, results directory, and module entry points.

## Main entry points

```bash
python -m chapter_analyses.surveillance.policy_sequences_over_time
python -m chapter_analyses.surveillance.policy_index_comparison
python -m chapter_analyses.genomic_networks.build_cluster_tables
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4
python -m chapter_analyses.sse_detection.lib.sse.detection
```

Run all commands from the Scotland repository root in the project environment.
