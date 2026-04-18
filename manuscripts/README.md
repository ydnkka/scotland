# Manuscripts — Scotland SARS-CoV-2 sequence-level clustering

This directory holds three self-contained manuscript projects, each with its own scope, figure code, and regression models. Each paper draws on the sequence-level analysis dataset and the cluster-level feature tables produced by the scripts in `analysis/`.

## Papers

| Dir | Title | Primary question |
|-----|-------|------------------|
| `paper1_socioeconomic/` | *Socioeconomic deprivation and SARS-CoV-2 transmission clustering in Scotland* | Does area-level deprivation (SIMD) predict cluster size and the probability of cluster membership, across VOC epochs? |
| `paper2_demographic/`   | *Demographic and vaccination correlates of SARS-CoV-2 transmission clusters in Scotland* | How do age, sex, and vaccination status shape cluster size, age diversity, and singleton status? |
| `paper3_spatial/`       | *Spatial structure of SARS-CoV-2 transmission clusters in Scotland* | How does geographic footprint vary with deprivation, urban-rural status, and VOC? |

Each paper folder has its own `README.md` with the full scope, hypotheses, figure list, model specifications, and target journals.

## Directory layout

```
manuscripts/
├── _common/                       Shared helpers
│   ├── data.py                    Config-aware parquet loaders, VOC epochs
│   ├── stats.py                   NB/logistic GLM helpers, bootstrap CIs
│   └── style.py                   Publication theme, palettes, figure sizes
├── paper1_socioeconomic/
│   ├── README.md
│   ├── Makefile
│   ├── make_figures.py            Orchestrator for all figures
│   ├── figures/                   One script per figure
│   ├── models/                    Regression-frame builders
│   ├── tables/                    CSV outputs (tidy regression results)
│   └── output/                    PDF + PNG figure outputs
├── paper2_demographic/            (same layout)
└── paper3_spatial/                (same layout)
```

## Running everything

From the repository root (where `config.yaml` lives):

```bash
# Install deps (only needed once)
pip install -r requirements.txt statsmodels scipy pyarrow matplotlib

# Build all figures for all three papers
make -C manuscripts all

# Or paper-by-paper
make -C manuscripts/paper1_socioeconomic figs
make -C manuscripts/paper2_demographic figs
make -C manuscripts/paper3_spatial figs

# Or a single figure
python -m manuscripts.paper1_socioeconomic.figures.fig3_simd_domain_forest
```

## Data contract

The figure scripts expect these parquets to exist (the first three are produced by scripts in `analysis/`):

- `data/processed/scotland_clustering_analysis_dataset.parquet`
- `data/processed/cluster_summary.parquet`
- `data/processed/cluster_demographic_features.parquet`
- `data/processed/cluster_simd_features.parquet` *(auto-derived on first use if absent)*

See `data/processed/analysis_dataset_description.md` for the full column dictionary.

## Conventions

- **Primary resolution.** Headline results use Leiden resolution `0.5` (`manuscripts._common.data.PRIMARY_RESOLUTION`). Resolution sensitivity appears as supplementary figures.
- **QC filter.** All analyses exclude sequences with `nextclade_qc` in `{mediocre, bad}`.
- **VOC epochs.** Derived from the data by `manuscripts.common.data.derive_voc_epochs_from_data()` — contiguous weeks where one WHO VOC holds ≥50% of sequenced cases, with Omicron split into BA.1 vs BA.2+ via `pango_lineage`. Exposed as `manuscripts.common.data.VOC_EPOCHS` for consistent stratification across all papers; falls back to the hardcoded `VOC_EPOCHS_DEFAULT` if the master parquet is unavailable.
- **Offset.** Cluster-size models use `log(wn_prop_sequenced)` as an offset so effect sizes are interpreted per sequenced case rather than per observed case.
- **Output formats.** Every figure is written to both `.pdf` (for submission) and `.png` (for quick review) at 400 dpi. TrueType fonts are embedded (`pdf.fonttype = 42`) so that panels remain editable in Illustrator.

## Reproducibility checklist before submission

- [ ] Pin versions in `requirements.txt` and capture a `pip freeze > requirements.lock`.
- [ ] Record the exact commit SHA used to build each final figure inside the figure-output directory (e.g. `git rev-parse HEAD > output/GIT_SHA`).
- [ ] Re-derive at least one headline number for each paper from the raw analysis dataset in a separate notebook.
- [ ] Confirm cell-size suppression for any small-count maps or tables (genomic surveillance data governance).
- [ ] Ensure resolution-sensitivity supplementary tables are generated with `--only` in the orchestrator loop.
