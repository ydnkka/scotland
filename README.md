# Scotland SARS-CoV-2 Genomic Surveillance and Network Analysis

This repository links Scottish SARS-CoV-2 sequence metadata to rolling-window genetic compatibility, Leiden clusters, descriptive surveillance, compatibility-network analysis, and superspreading-compatible (SSE) candidate characterisation.

EpiLink edges and cluster transitions are compatibility and continuity signals, not observed transmission events.

## Workflow

```text
raw sequence and linked metadata
  -> metadata preparation
  -> 3-week window × Pango-lineage groups
  -> TN93 pairwise distances
  -> EpiLink compatibility and sparsification
  -> Leiden clusters
  -> consolidated analysis dataset
  -> surveillance | genomic networks | SSE detection and characterisation
```

## Repository map

- `config.yaml`: method paths and parameters.
- `method/`: preprocessing, pairwise scoring, clustering, consolidation, and batch helpers.
- `utils/`: shared data, policy, and plotting helpers.
- `chapter_analyses/surveillance/`: policy timeline and lineage surveillance.
- `chapter_analyses/genomic_networks/`: Chapter 4 cohort, cluster, mixing, topology, and sensitivity analyses.
- `chapter_analyses/sse_detection/`: Chapter 5 transition graph, detector, figures, and Bayesian models.
- `assets/`: method schematics and research notes.
- `data/`: local raw and processed data; excluded from version control.

Detailed documentation:

- [Method pipeline](method/PIPELINE.md)
- [Chapter analyses](chapter_analyses/README.md)
- [Surveillance](chapter_analyses/surveillance/README.md) ([technical reference](chapter_analyses/surveillance/TECHNICAL.md))
- [Genomic networks](chapter_analyses/genomic_networks/README.md) ([technical reference](chapter_analyses/genomic_networks/TECHNICAL.md))
- [SSE detection](chapter_analyses/sse_detection/README.md) ([technical reference](chapter_analyses/sse_detection/TECHNICAL.md), [rationale](chapter_analyses/sse_detection/DETECTION_RATIONALE.md), [Bayesian models](chapter_analyses/sse_detection/BAYESIAN_MODELS.md))

## Environment

Run commands from the repository root. Install the Python dependencies, plus the external method tools `samtools`, `tn93`, and GNU parallel:

```bash
python -m pip install -r requirements.txt
```

`requirements.txt` includes the Bayesian fitting stack and CairoSVG used by the method-schematic builders.

## Method pipeline

The configured defaults are 3-week windows stepped weekly, groups of at least 2 sequences, a 29,903-base alignment, Leiden resolutions 0.1-0.8, seed 42, and compatibility sparsification at 0.001. Chapter 4 separately fixes its primary Leiden resolution at 0.3.

Run the method stages in order:

```bash
python3 method/01_prep_metadata.py
python3 method/02_gen_tn93_commands.py
./method/parallel_run.sh -c data/processed/group_fastas/tn93_commands.txt
python3 method/03_build_pairwise_dataset.py
python3 method/04_gen_cluster_commands.py
./method/parallel_run.sh -c data/processed/cluster_commands.txt
python3 method/05_consolidate.py
```

The shared downstream input is `data/processed/scotland_clustering_analysis_dataset.parquet`. Build the optional sparse-edge scheduling manifest before a full Chapter 4 mixing run:

```bash
python3 method/build_sparsified_edge_manifest.py
```

See [method/PIPELINE.md](method/PIPELINE.md) for inputs, outputs, restart behaviour, and CLI options.

## Analysis entry points

```bash
# Surveillance
python -m chapter_analyses.surveillance.policy_sequences_over_time
python -m chapter_analyses.surveillance.policy_index_comparison

# Chapter 4
python -m chapter_analyses.genomic_networks.build_cluster_tables
python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
python -m chapter_analyses.genomic_networks.build_cluster_pairwise_distance_summary --all-windows
python -m chapter_analyses.genomic_networks.build_sensitivity_tables
python -m chapter_analyses.genomic_networks.build_simd_validation
python -m chapter_analyses.genomic_networks.make_figures --skip-missing

# Chapter 5
python -m chapter_analyses.sse_detection.lib.sse.detection
python -m chapter_analyses.sse_detection.build_composition_tables
python -m chapter_analyses.sse_detection.lib.model.mixing
python -m chapter_analyses.sse_detection.lib.model.composition
python -m chapter_analyses.sse_detection.make_figures --skip-missing
```

Generated analysis outputs are written below each package's `results/` directory and are excluded from version control.

## Governance and reproducibility

PHS-linked raw inputs are confidential. Do not commit person-level data, credentials, machine-specific paths, or disclosive outputs. Small-cell flags are aids, not substitutes for formal disclosure review.

For a reproducible run, record the Git revision, environment and external-tool versions, `config.yaml`, raw-data versions, filters, seeds, thresholds, resolutions, window selection, development caps, and batch logs. Rebuild downstream outputs whenever cluster, transition, detector, or policy definitions change.
