# Scotland SARS-CoV-2 Genomic Surveillance and Network Analysis

This repository links Scottish SARS-CoV-2 sequence metadata to rolling-window genetic compatibility, Leiden clusters, descriptive surveillance, compatibility-network analysis, and superspreading-compatible candidate characterisation.

EpiLink edges and cluster transitions are compatibility and continuity signals, not observed transmission events. The EpiLink project is maintained at [ydnkka/epilink](https://github.com/ydnkka/epilink); use its [`evaluation` branch](https://github.com/ydnkka/epilink/tree/evaluation) for the version associated with this analysis.

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
- `utils/`: shared data, policy, mapping, and plotting helpers.
- `analyses/surveillance/`: policy timeline and lineage surveillance.
- `analyses/genomic_networks/`: cohort, cluster, mixing, topology, and sensitivity analyses.
- `analyses/sse_detection/`: transition graph, detector, figures, and Bayesian models.
- `results/`: central publication figure and LaTeX-table builders.
- `data/`: local raw and processed data; excluded from version control.

Detailed documentation:

- [Clustering pipeline](method/PIPELINE.md)
- [Pipeline output](data/processed/DATASET_DESCRIPTION.md)
- [Output analyses](analyses/README.md)
- [Surveillance](analyses/surveillance/README.md) ([technical reference](analyses/surveillance/TECHNICAL.md))
- [Genomic networks](analyses/genomic_networks/README.md) ([technical reference](analyses/genomic_networks/TECHNICAL.md))
- [SSE detection](analyses/sse_detection/README.md) ([technical reference](analyses/sse_detection/TECHNICAL.md), [rationale](analyses/sse_detection/README.md#detection-rationale), [Bayesian models](analyses/sse_detection/BAYESIAN_MODELS.md))
- [Publication results](results/README.md)

## Environment

Run commands from the repository root. Install the Python dependencies, plus the external method tools `samtools`, `tn93`, and GNU `parallel`:

```bash
python -m pip install -r requirements.txt
```

`requirements.txt` includes the Bayesian fitting stack and CairoSVG used by the method-schematic builders.

## Clustering pipeline

The configured defaults are 3-week windows stepped weekly, groups of at least 2 sequences, a 29,903-base alignment, Leiden resolutions 0.1-0.8, seed 42, and compatibility sparsification at 0.001. The genomic-network analysis separately fixes its primary Leiden resolution at 0.3.

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

The shared downstream input is `data/processed/scotland_clustering_analysis_dataset.parquet`. Build the optional sparse-edge scheduling manifest before a full genomic-network mixing run:

```bash
python3 method/build_sparsified_edge_manifest.py
```

See [method/PIPELINE.md](method/PIPELINE.md) for inputs, outputs, restart behaviour, and CLI options.

## Analysis entry points

```bash
# Analysis-local table and intermediate builders
python -m analyses.genomic_networks.build_cluster_summaries
python -m analyses.genomic_networks.build_mixing --all-windows --workers 4 --include-giants --giant-workers 1
python -m analyses.genomic_networks.build_sensitivity_tables
python -m analyses.genomic_networks.build_simd_validation
python -m analyses.sse_detection.lib.sse.detection
python -m analyses.sse_detection.build_composition_tables
python -m analyses.sse_detection.lib.model.mixing
python -m analyses.sse_detection.lib.model.composition

# Publication figures and LaTeX table fragments
python -m results.make_figures --skip-missing
python -m results.make_tables --skip-missing
```

Use `python -m results.make_figures --list` or `python -m results.make_tables --list`
to inspect available builders. Individual builders can be selected by fully
qualified name, for example:

```bash
python -m results.make_figures surveillance:policy_sequences_over_time
python -m results.make_tables genomic_networks:tab_policy_denominators
```

Analysis-local CSV/parquet outputs are written below each package's `results/`
directory. Project-level publication figures are written to `results/figures/`,
and project-level LaTeX table fragments are written to `results/tables/`.

## Governance and reproducibility

PHS-linked raw inputs are confidential. Do not commit person-level data, credentials, machine-specific paths, or disclosive outputs. Small-cell flags are aids, not substitutes for formal disclosure review.

For a reproducible run, record the Git revision, environment and external-tool versions, `config.yaml`, raw-data versions, filters, seeds, thresholds, resolutions, window selection, development caps, and batch logs. Rebuild downstream outputs whenever cluster, transition, detector, or policy definitions change.
