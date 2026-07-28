# Scotland SARS-CoV-2 Data Integration and Clustering Pipeline

The pipeline prepares linked metadata, groups sequences by overlapping time window and Pango lineage, computes TN93 distances, scores EpiLink compatibility, applies Leiden clustering, and consolidates the results.

All configured paths are relative to the repository root unless absolute. Run commands from that root.

## Execution

```bash
python3 method/01_prep_metadata.py
python3 method/02_gen_tn93_commands.py
./method/parallel_run.sh -c data/processed/group_fastas/tn93_commands.txt
python3 method/03_build_pairwise_dataset.py
python3 method/04_gen_cluster_commands.py
./method/parallel_run.sh -c data/processed/cluster_commands.txt
python3 method/05_consolidate.py
```

TN93 jobs must finish before stage 03; clustering jobs must finish before stage 05. Keep `data.processed.group_fasta_dir` through clustering because its `.ids` files retain sequences with no edge above the sparsification threshold.

Most Python stages accept `--config`, `--root`, and `--log-level`; inspect `--help` for stage-specific filters and overrides.

## Configuration

`config.yaml` has four top-level sections:

- `pipeline`: window width/step, minimum group size, alignment length, Leiden resolutions, seed, and sparsification threshold;
- `tn93`: TN93 distance threshold, ambiguous-site mode and fraction, minimum overlap, quiet flag, and the required Nextclade overall QC status;
- `data`: raw-input and processed-output paths;
- `policy`: the Scotland region label and the ordered, contiguous policy-period calendar used to annotate samples and downstream analyses.

The current defaults are 3-week half-open windows stepped weekly, minimum group size 2, alignment length 29,903, Leiden resolutions 0.1–0.8, seed 42, and sparsification `epilink_compatibility > 0.001`. TN93 reports distances below 0.0005, uses `resolve` with a 0.05 ambiguity fraction and 500-nt minimum overlap, and receives only Nextclade `good` sequences.

## Stage contracts

### 01: metadata preparation

`01_prep_metadata.py` builds all tabular inputs used by the later method stages. It links confidential PHS records to COG-UK/Nextclade annotations, prepares area-level surveillance and sociodemographic covariates, and constructs the daily Scotland policy calendar.

The stage reads:

| Config key                                        | Source and use                                                                                             |
| ------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `metadata_csv`                                    | PHS sequenced-case records: specimen/patient identifiers, collection date, demographics, and Data Zone     |
| `nextclade_tsv`                                   | COG-UK sequence name, clade, WHO variant label, Pango lineage, and `qc.overallStatus`                      |
| `testing_csv`                                     | PHS test records used for specimen attributes, reinfection history, and daily Data Zone testing totals     |
| `vaccination_csv`                                 | PHS vaccination records used for each sample's latest prior dose and daily Data Zone vaccination summaries |
| `simd_csv`                                        | SIMD 2020v2 population, deprivation, urban/rural, Local Authority, and Health Board attributes             |
| `geography_shp`                                   | 2011 Data Zone boundaries used for geometry, centroid coordinates, and area                                |
| `daily_hb_trends_csv`                             | Daily Health Board cases, tests, hospital/ICU burden, and reinfections                                     |
| `oxcgrt_stringency_csv`, `oxcgrt_containment_csv` | Daily Scotland Oxford policy indices                                                                       |
| `policy.periods`                                  | Ordered period codes, labels, inclusive date bounds, and broader policy eras                               |

It writes seven processed datasets:

| Config key    | Content                                                                                                                                                                                                              |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `metadata`    | One row per retained sequenced specimen, with genomic annotations/QC, demographics, Data Zone coordinates, linked test attributes, policy period, reinfection status, and latest vaccination on or before collection |
| `testing`     | One row per collection date and Data Zone, with total, positive, negative, PCR-positive, LFD-positive, and care-home-linked test counts                                                                              |
| `vaccination` | One row per vaccination date and Data Zone, with unique patients vaccinated and age/dose-number summary statistics; these are vaccination-event summaries rather than deduplicated ever-vaccinated residents         |
| `simd`        | One row per retained Data Zone, with population, SIMD ranks/groups, urban/rural class, Local Authority, and Health Board fields available in the source release                                                      |
| `geography`   | GeoParquet with one row per Data Zone: boundary geometry, OSGB centroid coordinates, area in km², and joined SIMD attributes                                                                                         |
| `hb_trends`   | One row per date and Health Board, excluding the Scotland-wide aggregate, with case, testing, hospital, ICU, and reinfection measures                                                                                |
| `policy`      | One row per calendar date, with Scotland stringency and containment indices plus configured period code, label, inclusive bounds, order, and policy era                                                              |

#### Sequence metadata derivations

- Nextclade fields are joined by the PHS sequence identifier parsed from `seqName`. Immediately afterwards, records are restricted to `tn93.nextclade_qc` (currently `good`), so the processed metadata and every downstream stage use the same quality cohort.
- Duplicate specimens retain their earliest metadata record. Rows missing a required modelling field (`datazone`, collection date, patient/sequence identifiers, sex, age band, clade, lineage, or Nextclade QC) are removed.
- Source age bands are retained as `age_band`, converted to numeric `age_midpoint`, and collapsed to the analytical `age_group` categories used by downstream models.
- `test_type`, raw test reason, and S-gene status are linked by specimen. `test_reason_raw` preserves the source value; `test_reason` provides the stable analytical categories. Missing reasons remain `missing`, and previously unseen non-missing values are logged and assigned to `other`.
- `is_reinfection` identifies a positive specimen collected at least 90 days after the same patient's preceding positive test.
- Vaccination fields describe the most recent recorded dose on or before the sample date: date, dose number, product, booster flag, and derived vaccination indicator.
- Each sample date must map to a configured policy period. The metadata stores `policy_period`, `policy_period_label`, and `policy_era`; the daily stringency and containment values remain in the standalone `policy` dataset.

#### Policy calendar

The configured policy periods must be unique, contiguous, and non-overlapping. Their start and end dates are inclusive, and `period_order` follows chronological order. The two wide OxCGRT files are reshaped to daily Scotland series and outer-joined to this calendar. Missing index values are retained and logged, whereas an uncovered sample date is an error. Use `python3 method/01_prep_metadata.py --policy-only` to rebuild just the processed policy dataset.

SIMD rows containing a null in any retained field are dropped with a warning. Because the later SIMD merge is inner, this can also remove samples whose Data Zones are not retained in the processed SIMD dataset.

### 02: TN93 command generation

`02_gen_tn93_commands.py` reapplies `tn93.nextclade_qc` as a guard for metadata created by an earlier pipeline version, deduplicates by `sequence_id`, creates half-open windows (`start <= collection_date < end`), and partitions each window by Pango lineage. Groups smaller than `pipeline.min_group_size` are skipped.

For each retained group it writes `group_fastas/<stem>.ids` and adds a compound `samtools faidx && tn93` command to `tn93_commands_file`. The TN93 invocation uses the configured threshold (`-t`), ambiguous-site strategy (`-a`), ambiguity fraction (`-g`), minimum overlap (`-l`), and quiet mode (`-q`). Output stems have the form `W042_BA.2_317`. A sequence can occur in up to three windows with the default width and step.

### Batch runner

`parallel_run.sh` executes one command per line with GNU parallel:

```bash
./method/parallel_run.sh -c COMMANDS_FILE [-j N] [--retries N] [--timeout SECS]
  [--joblog FILE] [--resume-failed] [--tmpdir DIR] [--no-compress]
  [--progress] [--dry-run]
```

The default worker count is the detected CPU count. The default joblog is beside the command file. Job buffers use `--tmpdir`, then `$TMPDIR`, `/var/tmp`, `/tmp`, or a `.parallel-tmp` directory beside the joblog.

### 03: pairwise dataset

`03_build_pairwise_dataset.py` processes every matching TN93 CSV:

1. remove invalid/self pairs and deduplicate undirected endpoints;
2. calculate `snp_distance = round(tn93_distance * alignment_length)`;
3. calculate absolute sampling-date difference in days;
4. calculate EpiLink compatibility;
5. write one Zstd-compressed parquet per group.

The schema is:

```text
window_id, pango_lineage, nunique_sequences, id1, id2,
tn93_distance, snp_distance, temporal_distance, epilink_compatibility
```

Existing parquets are skipped unless `--force` is used. An EpiLink exception is logged and produces missing compatibility values for that group rather than aborting the stage.

### EpiLink wrapper

`epilink_wrapper.py` creates one EpiLink instance per process with default SARS-CoV-2 epidemic and inference target parameters.

The EpiLink seed is a code constant (`RNG_SEED`), not read from `config.yaml`.

### 04 and clustering: Leiden assignments

`04_gen_cluster_commands.py` creates one `cluster_pairwise.py` command per pairwise parquet and forwards the configured resolutions, seed, sparsification threshold, and matching `.ids` path. `--include` and `--exclude` restrict stems.

`cluster_pairwise.py`:

1. reads edges whose compatibility is strictly above the threshold;
2. deduplicates undirected pairs;
3. adds endpoint-free nodes from `.ids`;
4. runs weighted Leiden for every configured resolution with 10 iterations;
5. writes one row per sequence and resolution: `sequence_id`, `window_id`, `resolution`, `cluster_id`.

If no edge survives, every known sequence becomes a singleton. If `.ids` is absent, membership is inferred from all pairwise endpoints before thresholding; truly endpoint-free sequences cannot then be recovered.

### 05: consolidation

`05_consolidate.py` applies the configured Nextclade QC filter before constructing windows or singleton assignments, rejects stale cluster assignments outside that cohort, and joins the processed sequence metadata (including sample-date policy fields) to SIMD, testing, vaccination, and Health Board data. It verifies that every final row has the configured QC status and writes `data.processed.analysis_dataset`.

Important derived fields include theoretical window dates; cluster size, duration, and Data Zone count; same-day and 7-row testing positivity; cumulative Data Zone sequence/test counts; cumulative vaccination-event count divided by population; incidence; and the latest Health Board report on or before collection. The vaccination quantity counts recorded dose events, not unique vaccinated residents. Patient IDs are replaced by run-specific `P000001`-style labels.

The SIMD join is inner and can remove sequences whose Data Zones are absent after SIMD cleaning. Health Board trends are skipped with a warning if the SIMD input lacks `HBcode`.

## Optional sparse-edge manifest

Chapter 4 uses an edge-count manifest to schedule very large pairwise files:

```bash
python3 method/build_sparsified_edge_manifest.py
```

It scans `epilink_compatibility` in batches, counts rows strictly above the configured threshold, and writes parquet plus CSV at `data.processed.sparsified_edge_manifest`. This is a downstream scheduling artifact, not a prerequisite for method consolidation.

## Dependencies and reproducibility

Python dependencies are in `requirements.txt`; external executables are `samtools` (version 1.23), [`tn93`](https://github.com/veg/tn93), and GNU `parallel` (version 20210822).

Record `config.yaml`, raw inputs, software versions, and batch joblogs. Leiden uses `pipeline.seed`; EpiLink uses its separate fixed seed. Changing alignment length, grouping, QC status, thresholds, resolution, or either seed requires rebuilding the affected intermediate and downstream outputs. Remove stale group FASTA/ID, TN93, pairwise, and cluster outputs when changing the QC cohort; consolidation rejects cluster assignments containing sequences outside the configured cohort. Do not mix same-stem chunks created under different settings.
