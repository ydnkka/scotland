# Scotland SARS-CoV-2 Clustering Pipeline

## Overview

The pipeline groups a national genomic surveillance dataset into temporal windows, partitions each window by Pango lineage, computes pairwise genetic distances with TN93, converts those distances into epidemiological compatibility weights with EpiLink, and infers transmission clusters using the Leiden community-detection algorithm. Clinical, demographic, and area-level contextual data are then joined to produce a single rectangular analysis dataset.

All paths are declared in `config.yaml` (repo root).

---

## Execution order

```
Step 1  python3 method/01_prep_metadata.py
Step 2  python3 method/02_gen_tn93_commands.py
Step 3  ./method/parallel_run.sh -c <tn93_commands_file>
Step 4  python3 method/04_gen_cluster_commands.py
Step 5  ./method/parallel_run.sh -c <cluster_commands_file>
Step 6  python3 method/05_consolidate.py
```

Steps 3 and 5 run externally generated command files in parallel using GNU parallel. Step 3 must complete before Step 4, and Step 5 must complete before Step 6.

---

## File descriptions

### `config.yaml`  *(repo root)*

Central configuration for the entire pipeline. Contains three top-level sections:

- **`pipeline`** — window parameters (`window_size_weeks`, `step_weeks`), minimum group size for TN93/Leiden processing, SARS-CoV-2 genome length (`alignment_length`), Leiden resolution sweep, random seed, and edge-weight sparsification threshold.
- **`tn93`** — distance threshold, ambiguous-site handling, and verbosity flag passed to the `tn93` executable.
- **`data.raw`** — paths to raw input files (COG-UK FASTA, Nextclade TSV, PHS metadata/testing/vaccination CSVs, SIMD CSV, datazone shapefile).
- **`data.processed`** — paths for every intermediate and final output written by the pipeline.

---

### `01_prep_metadata.py`

**Role:** Reads all raw data sources and writes five cleaned parquet files used by all downstream steps.

**Inputs (from `data.raw`):**
- `metadata_csv` — PHS sequenced-case records (demographics, specimen IDs, datazones)
- `nextclade_tsv` — Nextclade QC calls, Pango lineage, WHO VOC label, clade
- `vaccination_csv` — individual vaccination events (patient ID, dose number, date, datazone)
- `testing_csv` — individual PCR/LFD test records (result, date, patient ID, datazone)
- `simd_csv` — SIMD 2020v2 domain ranks/deciles for all Scottish datazones
- `geography_shp` — datazone boundary shapefile (OSGB36/EPSG:27700)

**Outputs (to `data.processed`):**
| Key | Content |
|---|---|
| `metadata` | One row per sequenced case: sequence ID, collection date, demographics, Pango lineage, Nextclade QC, datazone centroid coordinates, and the patient's most-recent prior vaccination dose/date at time of sampling. |
| `testing` | Daily positive/negative/total test counts aggregated to datazone level. |
| `vaccination` | Daily vaccination counts and age/dose statistics aggregated to datazone level. |
| `simd` | One row per datazone: population, SIMD rank/quintile/decile/vigintile, and all seven domain ranks. |
| `geography` | Geoparquet with datazone centroids (easting/northing) and SIMD attributes. |

**Key logic:**
- Age-band strings (e.g. `30-34`, `85+`) are converted to numeric midpoints.
- Vaccination history: for each sequenced specimen the most-recent prior dose on or before the collection date is looked up by patient ID, producing `vacc_dose_number` and `vaccination_date`.
- SIMD null values are dropped with a logged warning rather than raising a hard error.

---

### `02_gen_tn93_commands.py`

**Role:** Partitions the sequence dataset into (window × lineage) groups and generates a command file that extracts per-group FASTAs and computes pairwise TN93 distances.

**Inputs:** `data.processed.metadata` parquet.

**Outputs:**
- One `.ids` text file per group inside `group_fasta_dir` — each line is a `sequence_id` (FASTA header) belonging to that group.
- `tn93_commands_file` — one shell command per group of the form:
  ```
  samtools faidx -o <group>.fasta <all_seqs>.fasta.gz -r <group>.ids && \
  tn93 -q -t <threshold> -a <ambig> -o <group>.csv <group>.fasta
  ```

**Windowing:** 3-week sliding windows with 1-week stride (configurable). Groups with fewer than `min_group_size` sequences are skipped.

**Naming convention:** `W{window_idx}_{lineage_slug}_{n_seqs}` — e.g. `W042_BA.2_317`.

---

### `parallel_run.sh`

**Role:** Executes a command file (one shell command per line) in parallel using GNU parallel.

**Usage:**
```bash
./method/parallel_run.sh -c <commands_file> [options]

Options:
  -j N            Number of parallel workers (default: all CPU cores)
  --retries N     Retry failed jobs N times
  --timeout SECS  Per-job timeout
  --resume-failed Re-run only previously failed jobs (uses joblog)
  --progress      Show live progress bar
  --dry-run       Print first 5 commands and exit
```

Used twice: once for TN93 distance computation (Step 3) and once for cluster inference (Step 5).

---

### `epilink_wrapper.py`

**Role:** Thin module-level wrapper around the `EpiLink` scorer. Called by `03_process_group.py` to convert (SNP distance, temporal distance) pairs into epidemiological compatibility weights in [0, 1].

**EpiLink model parameters** (all in `NATURAL_HISTORY`):
- Incubation-period gamma distribution (shape 5.807, scale 0.948)
- Latent-period shape 3.38; symptomatic rate 37%
- Transmission rate ratio 2.29 (pre-symptomatic vs symptomatic)
- SARS-CoV-2 substitution rate 1×10⁻³ substitutions/site/year
- Genome length 29903 nt — **must match `pipeline.alignment_length` in `config.yaml`**
- `MAXIMUM_DEPTH = 0`: only direct transmission links are scored
- `TARGET = ("ad(0)", "ca(0,0)")`: ancestor–descendant and common-ancestor at depth 0
- 10,000 Monte-Carlo samples

The EpiLink instance is constructed once per process (module-level singleton) so the MC lookup tables are built only on first call.

---

### `03_process_group.py`

**Role:** Processes a single (window × lineage) group end-to-end: TN93 distances → EpiLink weights → Leiden clusters. Designed to be called in parallel via `parallel_run.sh`.

**Inputs (all passed as CLI arguments generated by `04_gen_cluster_commands.py`):**
- `--tn93-csv` — CSV with columns `ID1`, `ID2`, `Distance` (proportional, not SNP count)
- `--seq-ids` — `.ids` file listing all sequence IDs in this group
- `--metadata` — processed metadata parquet (used only to look up collection dates)
- `--out-long-dir` — directory where the output parquet is written

**Processing steps:**
1. Load and clean TN93 distances; deduplicate undirected pairs.
2. Convert proportional distances to SNP counts: `round(Distance × alignment_length)`.
3. Look up collection dates from metadata; compute absolute temporal distances in days.
4. Score each pair with EpiLink → compatibility weight in [0, 1].
5. Drop edges with weight ≤ `sparsification` (default 1×10⁻⁴).
6. Build a weighted igraph graph; add any isolates (sequences with no edges above threshold).
7. Run Leiden community detection at each resolution in `leiden_resolutions`.
8. Write long-format parquet: one row per (sequence × resolution) with columns `sequence_id`, `window_id`, `resolution`, `cluster_id`.

**Fallback:** if TN93 produced no usable edges (empty group, all distances missing, or all weights below threshold), every sequence is assigned its own singleton cluster at every resolution.

**Output:** `<stem>.parquet` in `cluster_long_dir`.

---

### `04_gen_cluster_commands.py`

**Role:** Pairs each `.csv` in `tn93_results_dir` with its `.ids` file in `group_fasta_dir` and generates a command file that calls `03_process_group.py` once per group.

**Inputs:** `tn93_results_dir` and `group_fasta_dir` directories (populated by Steps 2–3).

**Output:** `cluster_commands_file` — one `python3 method/03_process_group.py ...` command per group, with all pipeline parameters forwarded from `config.yaml`.

Supports `--include` and `--exclude` regex flags for reprocessing subsets of groups without regenerating the full command file.

---

### `05_consolidate.py`

**Role:** Combines all per-group cluster parquets into a single wide analysis dataset by joining in all metadata, SIMD deprivation data, daily testing/vaccination aggregates, and health-board surveillance trends.

**Inputs:** all `data.processed` parquets produced by Steps 1–5.

**Processing steps:**
1. Concatenate all per-group long parquets from `cluster_long_dir`.
2. Add singleton assignments for lineage/window groups with exactly one sequence (skipped by TN93 step when `min_group_size ≥ 2`).
3. Inner-join with SIMD (logs sequences dropped from datazones absent in SIMD).
4. Rename `vaccination_date` → `vacc_date_prior` to preserve the individual-level vaccination lookback from Step 1.
5. Join daily testing and vaccination aggregates on `(collection_date, datazone)`.
6. Attach theoretical window date boundaries (`wn_start_date`, `wn_mid_date`, `wn_end_date`) from the configured window parameters — not from the empirical range of sample dates.
7. **Cumulative vaccination coverage** (`dz_cum_vaccinated`, `dz_cum_prop_vaccinated`): running total of vaccination events per datazone up to the sequence's collection date, attached via `pd.merge_asof(direction='backward')`.
8. **Cumulative sequencing fraction** (`dz_cum_sequences`, `dz_cum_positive_tests`, `dz_cum_prop_sequenced`): for each sequence, the proportion of all positive tests ever recorded in its datazone that have a linked genome, as of the collection date. A measure of local genomic surveillance intensity over time.
9. **7-day rolling test positivity** (`dz_7d_test_positivity`): sum of positive and total tests over the 7 most-recent testing-data rows in the datazone, joined by direct merge on `(collection_date, datazone)`.
10. **Cluster descriptors** (`cluster_size`, `cluster_n_datazones`, `cluster_start_date`, `cluster_end_date`, `cluster_duration_days`): per-cluster aggregates computed on the fully merged dataset and joined back on `cluster_id`.
11. **Derived variables:** `days_since_vaccination` (collection_date − vacc_date_prior, days; NaN if unvaccinated); `dz_test_positivity` (same-day positivity rate); `dz_cum_incidence_per_capita` (cumulative positive tests per head of datazone population).
12. **Health-board daily trends** (hospital admissions/occupancy, ICU, reinfections): joined via `pd.merge_asof(direction='backward', by='dz_health_board_code')` using the `dz_health_board_code` field from SIMD. Skipped with a warning if the SIMD release lacks `HBcode`.
13. Anonymise `patient_id` (replaced with zero-padded integers `P000001`, …).
14. Write final parquet to `analysis_dataset`.

**Output columns (~80):** window metadata, sequence/cluster identifiers, cluster descriptors, sample-level demographics, individual vaccination history (including product name, booster flag, and days since vaccination), specimen test attributes (test type, reason, S-gene status), datazone SIMD attributes (including local authority and health board codes), daily testing counts (PCR/LFD/care-home breakdown), test positivity (point and 7-day rolling), cumulative vaccination coverage, cumulative sequencing fraction, cumulative incidence per capita, and health-board daily surveillance trends.

---

## Data flow diagram

```
config.yaml
    │
    ▼
01_prep_metadata.py
    ├─► scotland_sequence_metadata.parquet
    ├─► scotland_testing.parquet
    ├─► scotland_datazone_vaccinations.parquet
    ├─► scotland_datazone_simd_data.parquet
    └─► scotland_geography.parquet
          │
          ▼
    02_gen_tn93_commands.py
          ├─► group_fastas/*.ids  (one per window×lineage group)
          └─► tn93_commands.txt
                    │
                    ▼
          parallel_run.sh  [GNU parallel]
                    └─► tn93_results/*.csv  (pairwise distances)
                              │
                              ▼
                    04_gen_cluster_commands.py
                              └─► cluster_commands.txt
                                        │
                                        ▼
                              parallel_run.sh  [GNU parallel]
                                        └─► cluster_assignments_long/*.parquet
                                                  │
                                                  ▼
                                        05_consolidate.py
                                                  └─► scotland_clustering_analysis_dataset.parquet
```

---

## Dependencies

| Package | Use |
|---|---|
| `pandas`, `numpy` | Data manipulation throughout |
| `geopandas` | Datazone shapefile reading and centroid computation (`01`) |
| `pyarrow` / `fastparquet` | Parquet I/O |
| `pyyaml` | Config parsing |
| `igraph` (`python-igraph`) | Graph construction and Leiden clustering (`03`) |
| `epilink` | Epidemiological compatibility scoring (`epilink_wrapper`, `03`) |
| `samtools` | Per-group FASTA extraction from compressed multi-sequence file (`02`→`03`) |
| `tn93` | Pairwise TN93 distance computation (`02`→`03`) |
| `GNU parallel` | Parallel execution of TN93 and clustering steps |

---

## Reproducibility notes

- All random operations (Leiden, EpiLink MC sampling) use `pipeline.seed` from `config.yaml` (default 42).
- Window boundaries are fully determined by `window_size_weeks`, `step_weeks`, and the date range of the metadata — rerunning with the same config and data produces identical windows.
- `genome_length` in `epilink_wrapper.py` and `pipeline.alignment_length` in `config.yaml` must always be equal (both default to 29903 nt for the SARS-CoV-2 reference genome). Changing one without the other will silently corrupt SNP-to-distance conversion.
