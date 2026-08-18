# Scotland SARS-CoV-2 Sequence-Level Analysis Dataset

`scotland_clustering_analysis_dataset.parquet`

---

## Overview

The analysis dataset contains sequence-window-resolution records for Scottish
SARS-CoV-2 genomes passing the configured Nextclade QC filter (`good` by
default). It links each retained sequence to genomic annotations, anonymised PHS
testing and vaccination fields, policy period, Data Zone attributes, and Health
Board context.

Rows are generated from three-week analysis windows stepped weekly. Within each
window, sequences are grouped by Pango lineage, scored with TN93 and EpiLink,
and clustered with Leiden across the configured resolutions. Clusters are
genomic-temporal compatibility groups, not observed transmission chains.

Dataset dimensions depend on the raw-data release and pipeline configuration;
calculate them from the generated Parquet rather than treating them as fixed.

---

## Row structure

One row = one sequence × one analysis window × one Leiden resolution.

A sequence collected in the overlap of consecutive windows appears in multiple windows. At eight resolutions per window, a sequence in a single window produces eight rows. The full combination for a typical sequence therefore ranges from 8 rows (one window only) to 24 rows (three windows). Filtering to a single resolution (e.g. `resolution == 0.3`) reduces the dataset to one row per sequence per window.

---

## Columns

### Window-level identifiers and summaries

| Column              | Type  | Description                                                                                                                              |
| ------------------- | ----- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `window_idx`        | int   | Sequential integer index of the analysis window (1, 2, 3, …)                                                                             |
| `window_id`         | str   | Zero-padded window label used in cluster IDs (e.g. `W001`, `W042`)                                                                       |
| `wn_start_date`     | date  | Theoretical inclusive start date of the window (from configured step/size, not empirical min of sample dates)                            |
| `wn_mid_date`       | date  | Midpoint date of the window                                                                                                              |
| `wn_end_date`       | date  | Theoretical exclusive end date of the window                                                                                             |
| `wn_no_sequences`   | int   | Number of unique QC-retained sequences with a collection date within this window                                                         |
| `wn_positive_tests` | int   | Total positive PCR/LFD tests recorded across all Data Zones within the window                                                            |
| `wn_prop_sequenced` | float | Fraction of window positive tests represented by QC-retained sequences (`wn_no_sequences / wn_positive_tests`); NaN if no positive tests |

---

### Sequence and cluster identifiers

| Column        | Type  | Description                                                                                                                                                                                                                                                                                                                                                                     |
| ------------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sequence_id` | str   | COG-UK sequence identifier (format `Scotland/LIVE-XXXXX/YYYY`); primary key to genomic data                                                                                                                                                                                                                                                                                     |
| `patient_id`  | str   | Anonymised patient identifier (format `P000001`); one patient may have multiple sequences                                                                                                                                                                                                                                                                                       |
| `resolution`  | float | Leiden resolution parameter at which `cluster_id` was assigned (one of 0.1, 0.2, …, 0.8)                                                                                                                                                                                                                                                                                        |
| `cluster_id`  | str   | Cluster assignment, formatted `{window_id}\|{lineage}\|R{resolution}\|C{n}` for multi-sequence clusters (e.g. `W042\|BA.2\|R0.3\|C017`) or `{window_id}\|{lineage}\|R{resolution}\|S0` for singletons (lineage groups of size 1, not processed by TN93/Leiden) or  `{window_id}\|{lineage}\|R{resolution}\|S{n}` for singletons (processed by TN93/Leiden). `n` starts from `1` |

---

### Cluster descriptors

These fields summarise the inferred transmission cluster to which this sequence belongs at the given resolution. They are computed from all sequences in the same `cluster_id` after the full metadata join.

| Column                  | Type | Description                                                                                                                                     |
| ----------------------- | ---- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `cluster_size`          | int  | Number of unique sequences assigned to this cluster                                                                                             |
| `cluster_n_datazones`   | int  | Number of distinct Data Zones represented in the cluster; a coarse measure of geographic spread                                                 |
| `cluster_start_date`    | date | Earliest collection date among sequences in the cluster                                                                                         |
| `cluster_end_date`      | date | Latest collection date among sequences in the cluster                                                                                           |
| `cluster_duration_days` | int  | `cluster_end_date − cluster_start_date` in whole days; 0 for clusters where all sequences were collected on the same day (including singletons) |

---

### Sample-level fields

#### Geography

| Column            | Type  | Description                                                                      |
| ----------------- | ----- | -------------------------------------------------------------------------------- |
| `collection_date` | date  | Date the sample was collected from the patient                                   |
| `datazone`        | str   | 2011 Scottish Data Zone (DZ) code of patient residence (e.g. `S01006506`)        |
| `dz_xcoord`       | float | Easting of the DZ polygon centroid, British National Grid (OSGB36 / EPSG:27700)  |
| `dz_ycoord`       | float | Northing of the DZ polygon centroid, British National Grid (OSGB36 / EPSG:27700) |

#### Demographics

| Column         | Type  | Description                                                                                                            |
| -------------- | ----- | ---------------------------------------------------------------------------------------------------------------------- |
| `sex`          | str   | Recorded patient sex (`Male`, `Female`, or other/unknown)                                                              |
| `is_female`    | float | Binary indicator: 1.0 if `sex == "Female"`, 0.0 otherwise                                                              |
| `age_band`     | str   | Age band at time of sampling (e.g. `30-34`, `75+`)                                                                     |
| `age_group`    | str   | Derived analytical age group: `00-04`, `05-14`, `15-24`, `25-64`, `65-74`, or `75+`                                    |
| `age_midpoint` | float | Midpoint of `age_band` in years (open-ended bands such as `75+` use `lower + 2.5`); intended as a continuous covariate |

The age-group mapping retains `00-04` and `75+`; combines `05-09` and `10-14` as `05-14`; combines `15-19` and `20-24` as `15-24`; combines bands from `25-29` through `60-64` as `25-64`; and combines `65-69` and `70-74` as `65-74`. The original `age_band` is retained alongside this derived variable.

#### Policy context

These sample-level fields are assigned by joining `collection_date` to `data/processed/scotland_policy.parquet`. They describe policy on the sample collection date, not the midpoint of the rolling analysis window. Analyses requiring window-level policy context join `wn_mid_date` to the policy calendar separately.

| Column                | Type | Description                                                                                                                                                                    |
| --------------------- | ---- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `policy_period`       | str  | Ordered compact policy-period code (e.g. `P3`, `L2`, `OM`, `PR`) assigned from `collection_date`                                                                               |
| `policy_period_label` | str  | Human-readable label corresponding to `policy_period` (e.g. `Route map phase 3`, `Second lockdown`, `Omicron wave`, `Post-restriction`)                                        |
| `policy_era`          | str  | Broader policy grouping: `early_restriction_easing`, `autumn_winter_restrictions`, `spring_summer_2021_easing`, `near_normal_delta`, `omicron_response`, or `post_restriction` |

#### Vaccination history

All vaccination fields reflect the patient's most-recent dose received on or before `collection_date`.

| Column                   | Type  | Description                                                                                                                                    |
| ------------------------ | ----- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `is_vaccinated`          | float | 1.0 if the patient had received ≥1 dose by `collection_date`, 0.0 otherwise                                                                    |
| `vacc_dose_number`       | float | Number of doses received by `collection_date` (0 if unvaccinated)                                                                              |
| `vacc_date_prior`        | date  | Date of the patient's most-recent vaccine dose on or before sample collection; NaT if unvaccinated                                             |
| `vacc_product_name`      | str   | Name of the vaccine product received at the most-recent prior dose (e.g. `COVID-19 mRNA Vaccine Moderna`); NaN if unvaccinated or not recorded |
| `vacc_booster`           | float | 1.0 if the most-recent prior dose was classified as a booster, 0.0 otherwise; NaN if unvaccinated or not recorded                              |
| `days_since_vaccination` | float | Days between `vacc_date_prior` and `collection_date`; NaN if unvaccinated                                                                      |

#### Test record

These fields are joined from the PHS testing record linked to the sequence by `specimen_id`.

| Column            | Type  | Description                                                                                                                                                                                                                                              |
| ----------------- | ----- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `test_type`       | str   | Type of test that yielded the sequenced positive (`PCR` or `ANTIGEN`)                                                                                                                                                                                    |
| `test_reason_raw` | str   | Original reason value from the linked PHS testing record; retained for provenance and NaN when no reason was recorded                                                                                                                                    |
| `test_reason`     | str   | Canonical analytical grouping of `test_reason_raw`: `symptomatic_citizen`, `symptomatic_essential_worker`, `contact_tracing`, `confirmatory`, `isolation_scheme`, `clinical`, `surveillance_research`, `local_outbreak`, `travel`, `other`, or `missing` |
| `s_gene_status`   | str   | S-gene target failure (SGTF) status from the linked PCR result; values typically `DETECTED`, `NOT_DETECTED`, or NaN where not applicable or not reported                                                                                                 |
| `is_reinfection`  | float | 1.0 if this positive test occurred ≥90 days after the same patient's most-recent prior positive test (standard PHS/ECDC reinfection threshold); 0.0 for first infections and non-positive tests                                                          |

#### Genomic annotation

| Column          | Type | Description                                                                                             |
| --------------- | ---- | ------------------------------------------------------------------------------------------------------- |
| `pango_lineage` | str  | PANGO lineage assigned by Nextclade (e.g. `BA.2`, `AY.4.2`)                                             |
| `clade`         | str  | Nextclade clade assignment (e.g. `21K`, `22B`)                                                          |
| `who_voc`       | str  | WHO Variant of Concern label where applicable (e.g. `Delta`, `Omicron`); NaN for non-VOC lineages       |
| `nextclade_qc`  | str  | Nextclade overall QC status; every final row matches configured `tn93.nextclade_qc` (`good` by default) |

---

### Data Zone sociodemographic attributes

These are static attributes of the patient's DZ of residence from SIMD 2020v2 and the SG 2011 datazone boundary file. They do not vary across rows for the same DZ.

#### Population and area

| Column                      | Type  | Description                                                                 |
| --------------------------- | ----- | --------------------------------------------------------------------------- |
| `dz_population`             | int   | Total resident population of the DZ (SIMD 2020v2)                           |
| `dz_working_age_population` | int   | Working-age (16–64) resident population of the DZ                           |
| `dz_area_km2`               | float | Area of the DZ in km² (from `StdAreaKm2` in the SG 2011 boundary shapefile) |
| `dz_population_density`     | float | Resident population per km² (`dz_population / dz_area_km2`)                 |

#### SIMD 2020v2 deprivation

All ranks run from 1 (most deprived) to 6,976 (least deprived) across the 6,976 Scottish DZs.

| Column                    | Type | Description                                            |
| ------------------------- | ---- | ------------------------------------------------------ |
| `dz_simd_rank`            | int  | Overall SIMD 2020v2 rank                               |
| `dz_simd_quintile`        | int  | Overall SIMD quintile (1 = most deprived, 5 = least)   |
| `dz_simd_decile`          | int  | Overall SIMD decile (1 = most deprived, 10 = least)    |
| `dz_simd_vigintile`       | int  | Overall SIMD vigintile (1 = most deprived, 20 = least) |
| `dz_simd_income_rank`     | int  | SIMD income domain rank                                |
| `dz_simd_employment_rank` | int  | SIMD employment domain rank                            |
| `dz_simd_education_rank`  | int  | SIMD education, skills and training domain rank        |
| `dz_simd_health_rank`     | int  | SIMD health domain rank                                |
| `dz_simd_access_rank`     | int  | SIMD geographic access domain rank                     |
| `dz_simd_crime_rank`      | int  | SIMD crime domain rank                                 |
| `dz_simd_housing_rank`    | int  | SIMD housing domain rank                               |

#### Geographic classification and administrative codes

| Column                    | Type | Description                                                                                |
| ------------------------- | ---- | ------------------------------------------------------------------------------------------ |
| `dz_urban_rural_class`    | str  | Scottish Government 6-fold urban-rural classification of the DZ (e.g. `Large Urban Areas`) |
| `dz_local_authority`      | str  | Local authority (council area) name containing the DZ                                      |
| `dz_local_authority_code` | str  | Local authority code (e.g. `S12000036`)                                                    |
| `dz_health_board`         | str  | NHS Health Board name containing the DZ                                                    |
| `dz_health_board_code`    | str  | NHS Health Board code (e.g. `S08000024`); used to join health board daily trends           |

---

### Data Zone testing counts (on `collection_date`)

These are daily aggregates for the patient's DZ on the specific `collection_date`. They represent testing activity on that single day, not cumulative totals.

| Column                  | Type  | Description                                                                                                                                                                                                                        |
| ----------------------- | ----- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dz_total_tests`        | int   | All SARS-CoV-2 tests (PCR + LFD) performed on DZ residents on `collection_date`                                                                                                                                                    |
| `dz_positive_tests`     | int   | Positive tests (PCR + LFD) on DZ residents on `collection_date`                                                                                                                                                                    |
| `dz_negative_tests`     | int   | Negative tests on DZ residents on `collection_date`                                                                                                                                                                                |
| `dz_pcr_positive_tests` | int   | PCR-confirmed positive tests on DZ residents on `collection_date`                                                                                                                                                                  |
| `dz_lfd_positive_tests` | int   | LFD (antigen)-confirmed positive tests on DZ residents on `collection_date`                                                                                                                                                        |
| `dz_care_home_tests`    | int   | Tests linked to a care home facility (care_home_id not null) on DZ residents on `collection_date`                                                                                                                                  |
| `dz_test_positivity`    | float | Same-day test positivity rate in the DZ (`dz_positive_tests / dz_total_tests`); NaN if no tests recorded                                                                                                                           |
| `dz_7d_test_positivity` | float | 7-day rolling test positivity rate in the DZ: sum of positive tests over the 7 most-recent reporting dates in the DZ divided by sum of total tests over the same period; computed on the daily testing aggregate (min_periods = 1) |

---

### Data Zone vaccination (on `collection_date`)

| Column                   | Type  | Description                                                                                                                                                                                                                                                             |
| ------------------------ | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dz_total_vaccinated`    | int   | Number of distinct patients in the DZ who received any vaccine dose on `collection_date` (new events on that day only, not cumulative)                                                                                                                                  |
| `dz_cum_vaccinated`      | float | Running total of vaccination events in the DZ on or before `collection_date` (via backward merge_asof on vaccination date). Note: counts vaccination events per dose event, not unique individuals ever vaccinated, so slightly over-counts distinct vaccinated persons |
| `dz_cum_prop_vaccinated` | float | `dz_cum_vaccinated / dz_population`; cumulative vaccination coverage proxy                                                                                                                                                                                              |

---

### Data Zone cumulative surveillance

These fields characterise the local genomic surveillance effort and epidemic burden in the DZ up to and including `collection_date`.

| Column                        | Type  | Description                                                                                                                                                                                                                                           |
| ----------------------------- | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dz_cum_sequences`            | float | Total unique QC-retained genomes collected from the DZ on or before `collection_date`                                                                                                                                                                 |
| `dz_cum_positive_tests`       | float | Total positive tests (PCR + LFD) recorded in the DZ on or before `collection_date`                                                                                                                                                                    |
| `dz_cum_prop_sequenced`       | float | `dz_cum_sequences / dz_cum_positive_tests`; the fraction of the DZ's cumulative positive tests represented by a QC-retained genome as of this date — a local measure of usable genomic surveillance intensity over time; NaN if no positive tests yet |
| `dz_cum_incidence_per_capita` | float | `dz_cum_positive_tests / dz_population`; cumulative positive tests per head of DZ population as of `collection_date`                                                                                                                                  |

---

### Health Board daily trends

These fields are joined from the PHS daily health board COVID surveillance report for the patient's NHS Health Board (`dz_health_board_code`) using a backward-looking merge: the most-recent report on or before `collection_date` is used. All HB fields are NaN if `dz_health_board_code` is absent from the SIMD release or if no HB report precedes the collection date.

| Column                   | Type  | Description                                                                                           |
| ------------------------ | ----- | ----------------------------------------------------------------------------------------------------- |
| `hb_daily_positive`      | float | New positive cases reported in the Health Board on the matched report date                            |
| `hb_cumulative_positive` | float | Cumulative positive cases in the Health Board as of the matched report date                           |
| `hb_hospital_admissions` | float | New hospital admissions with confirmed COVID-19 in the Health Board on the matched report date        |
| `hb_hospital_occupancy`  | float | COVID-19 inpatients in hospital in the Health Board on the matched report date                        |
| `hb_icu_admissions`      | float | New ICU admissions with confirmed COVID-19 in the Health Board on the matched report date             |
| `hb_icu_occupancy_lt28d` | float | COVID-19 ICU patients with length of stay <28 days in the Health Board on the matched report date     |
| `hb_icu_occupancy_ge28d` | float | COVID-19 ICU patients with length of stay ≥28 days in the Health Board on the matched report date     |
| `hb_daily_reinfections`  | float | New reinfection cases reported in the Health Board on the matched report date                         |
| `hb_reinfection_rate`    | float | Reinfections as a percentage of all new positive cases in the Health Board on the matched report date |

---

## Notes on missing values

- Fields sourced from daily DZ aggregates (`dz_*`) are NaN for sequences from DZs with no testing/vaccination activity on `collection_date`; this is expected for low-density rural DZs.
- `vacc_date_prior`, `vacc_product_name`, `vacc_booster`, and `days_since_vaccination` are NaN for unvaccinated patients.
- `s_gene_status` is NaN for LFD tests and PCR tests where S-gene result was not reported.
- `test_reason_raw` is NaN where no reason was recorded in the PHS testing record. The corresponding canonical `test_reason` is the explicit string `missing`, rather than NaN. Recognised residual reasons and previously unseen non-null values are grouped as `other`.
- `age_group`, `policy_period`, `policy_period_label`, and `policy_era` are complete for the collection dates and age bands represented in the current analysis dataset; preprocessing raises an error if an age band or collection date cannot be assigned.
- `nextclade_qc` is complete and restricted to the status configured by `tn93.nextclade_qc`; with the default configuration, every row is `good`.
- `who_voc` is NaN for lineages not designated as a WHO Variant of Concern or Interest.
- `dz_7d_test_positivity` uses all available prior rows (min_periods = 1), so values early in the pandemic are based on fewer than 7 days.
- `hb_*` fields are NaN if `dz_health_board_code` was absent in the SIMD release used, or if no Health Board report precedes the sequence's collection date.
- Sequences from Data Zones absent from the SIMD 2020v2 table are dropped before the final dataset is written (logged as a warning by the pipeline).

## Sources and pipeline

| Source                                                       | Content                                                                   |
| ------------------------------------------------------------ | ------------------------------------------------------------------------- |
| COG-UK aligned FASTA (`cog_all_scotland_aligned.fasta.gz`)   | Consensus genome sequences                                                |
| Nextclade TSV (`cog_all_scotland_nextclade.tsv`)             | QC calls, PANGO lineage, WHO VOC label, clade                             |
| PHS sequenced cases CSV (`AnnaSequencedCases.csv`)           | Specimen/patient IDs, collection date, demographics, DZ                   |
| PHS testing CSV (`PHS_CovidTesting_2023-02-22.csv`)          | Individual test records (type, reason, S-gene, result)                    |
| PHS vaccination CSV (`PHS_Vaccinations_2023-02-22.csv`)      | Individual dose records (date, product, dose number, booster)             |
| SIMD 2020v2 CSV (`2020v2_simd.csv`)                          | DZ-level deprivation ranks, population, HB/LA codes                       |
| SG datazone boundary shapefile (`sg_datazone_bdry_2011.shp`) | DZ polygon geometry, centroid coordinates, area                           |
| PHS daily HB trends CSV (`daily_trend_by_health_board.csv`)  | Daily hospital, ICU, reinfection, and case counts by Health Board         |
| OxCGRT Stringency and Containment/Health index CSVs          | Daily Scotland policy indices used to build the processed policy calendar |
| Policy-period specification (`config.yaml`)                  | Policy codes, labels, inclusive start/end dates, order, and broader eras  |

Pipeline steps: `01_prep_metadata.py` $\rightarrow$ `02_gen_tn93_commands.py` $\rightarrow$ TN93 (GNU parallel) $\rightarrow$ `03_build_pairwise_dataset.py` $\rightarrow$ `04_gen_cluster_commands.py` $\rightarrow$ Leiden clustering (GNU parallel) $\rightarrow$ `05_consolidate.py`. The Nextclade QC restriction is applied before window construction and retained through pairwise calculations, clustering, and final consolidation. See `method/PIPELINE.md` for full details.
