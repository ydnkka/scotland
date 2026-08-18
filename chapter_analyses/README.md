# Analysis Packages

Importable analysis packages built from `data/processed/scotland_clustering_analysis_dataset.parquet` and, where required, the pairwise compatibility parquets.

| Package | Scope | Main entry point |
| --- | --- | --- |
| [`surveillance`](surveillance/README.md) | Policy timeline, sequence counts, lineage replacement, and sequencing coverage | `python -m chapter_analyses.surveillance.policy_sequences_over_time` |
| [`genomic_networks`](genomic_networks/README.md) | Cohort, clusters, compatibility mixing/topology, and sensitivity analyses | `python -m chapter_analyses.genomic_networks.build_cluster_tables` |
| [`sse_detection`](sse_detection/README.md) | Transition graph, candidate detector, composition, and Bayesian characterisation | `python -m chapter_analyses.sse_detection.lib.sse.detection` |

Run modules from the repository root. Each package writes generated artifacts to its own `results/` directory.
