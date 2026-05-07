# Part 2 Analysis Documentation

## Vaccination characterisation of SARS-CoV-2 genomic clusters in Scotland

This document records the current Part 2 analysis implemented across two scripts:

```text
part2/cluster_characterisation.py
part2/categorise_clusters.py
```

It is written as an implementation record: what data were used, how cluster-level
vaccination and demographic variables were constructed, what category rules were
applied, what summary tables were written, and how the manuscript figures were
produced.

The Part 2 research questions are:

1. Among sequenced SARS-CoV-2 cases, how did breakthrough-case frequency change
   over calendar time, and did this differ by vaccination-rollout age group, sex,
   overall SIMD quintile, or SIMD domain quintile?
2. Were clusters more homogeneous or more mixed by vaccination status than
   expected for cases from the same lineage and calendar window?
3. Did vaccine dose profile and dose recency among sequenced infected cases vary
   across overall SIMD quintiles, SIMD domain quintiles, or epidemic waves?
4. Do all-vaccinated, mixed-vaccination, and unvaccinated clusters have different
   demographic and geographic mixing profiles?

Part 2 is explicitly descriptive. Vaccination status among sequenced cases is
strongly confounded by rollout time, age eligibility, prior immunity, testing
behaviour, and sequencing selection, so causal vaccine-effectiveness inference
is not attempted.

---

## 1. Source Data

### 1.1 Sequence-level data

`cluster_characterisation.py` reads directly from the repository configuration:

```text
config.yaml -> data.processed.analysis_dataset
```

which resolves to:

```text
data/processed/scotland_clustering_analysis_dataset.parquet
```

The same QC filter and primary resolution used in Part 1 are applied:

| Filter            | Value |
|-------------------|------:|
| Nextclade QC      | `good` |
| Leiden resolution | `0.3` |

Additional columns specific to Part 2 are read from the same source:

```text
is_vaccinated
vacc_dose_number
vacc_date_prior
vacc_product_name
vacc_booster
days_since_vaccination
dz_simd_income_rank
dz_simd_employment_rank
dz_simd_education_rank
dz_simd_health_rank
dz_simd_access_rank
dz_simd_crime_rank
dz_simd_housing_rank
```

### 1.2 Part 1 cluster cache

`categorise_clusters.py` reads the cached Part 1 cluster table through
`utils.data.load_main_cluster_table`, which points to:

```text
part1/main/cache/main_cluster_table.parquet
```

The Part 1 cluster table provides cluster-level excess-mixing scores (SIMD,
age, sex, joint profile) computed in `part1/main/main_analysis.py`. These
scores are the basis for the Part 2 mixing categories.

---

## 2. Why This Is Part 2 And Not Part 1

Part 1 modelled deprivation effects on cluster size, geographic spread, and
within-cluster socioeconomic and demographic mixing using regression models with
surveillance adjustments. Part 2 adds vaccination status and dose history as
cluster-characterisation variables without fitting causal models for two reasons:

1. The vaccination rollout is too strongly confounded with wave timing, age,
   immune history, and testing behaviour to support causal inference without
   dedicated test-negative or cohort designs.
2. The descriptive questions — how did vaccination profiles look across waves,
   sizes, and social groups? — are scientifically worthwhile in their own right
   as a characterisation of the sequenced infected population.

---

## 3. Cluster Table Construction

`cluster_characterisation.py` rebuilds primary-resolution cluster aggregates
from the sequence-level table, adding vaccination and SIMD-domain variables.
Because this script does not apply Part 1's model-field completeness filter, it
retains 193,160 clusters rather than the 193,112 clusters used in the main
Part 1 regression models.

### 3.1 Vaccination variables

For each cluster, vaccination fields are constructed from sequence-level records
where vaccination status was known (`n_vaccination_known`):

| Cluster field                              | Construction                                                               |
|--------------------------------------------|----------------------------------------------------------------------------|
| `n_vaccination_known`                      | Sequences with known vaccination status                                    |
| `n_vaccinated`                             | Sequences with known vaccinated status                                     |
| `cluster_prop_vaccinated`                  | `n_vaccinated / n_vaccination_known`                                       |
| `n_boosted`                                | Sequences with booster dose recorded                                       |
| `cluster_prop_boosted_all_members`         | `n_boosted / cluster_size`                                                 |
| `cluster_prop_boosted_vaccinated_members`  | `n_boosted / n_vaccinated` (among vaccinated only)                         |
| `mean_dose_all_members`                    | Mean dose number across all sequences                                      |
| `mean_dose_vaccinated_members`             | Mean dose number among vaccinated sequences                                |
| `median_dose_vaccinated_members`           | Median dose number among vaccinated sequences                              |
| `mean_days_since_vaccination`              | Mean days since last prior vaccination among vaccinated sequences          |
| `median_days_since_vaccination`            | Median days since last prior vaccination among vaccinated sequences        |
| `index_is_vaccinated`                      | Vaccination status of the earliest-collected sequence in the cluster       |
| `index_vacc_dose_number`                   | Dose number of the earliest-collected sequence                             |
| `index_vacc_booster`                       | Booster flag for the earliest-collected sequence                           |
| `index_days_since_vaccination`             | Days since vaccination for the earliest-collected sequence                 |

### 3.2 Vaccination profile

Each cluster is assigned a `cluster_vaccination_profile` based on the
distribution of vaccination status among its members:

| Profile             | Rule                                               |
|---------------------|----------------------------------------------------|
| `none vaccinated`   | `cluster_prop_vaccinated == 0`                     |
| `all vaccinated`    | `cluster_prop_vaccinated == 1`                     |
| `mixed vaccination` | `0 < cluster_prop_vaccinated < 1`                  |
| `vaccination unknown` | No sequences with known vaccination status       |

### 3.3 Vaccination-status mixing category

Vaccination-status mixing is measured by observed-minus-expected pairwise
discordance for binary vaccination status, using the same framework as the
Part 1 mixing outcomes. Expected discordance is computed within the same
`window_id × pango_lineage` stratum.

| Category        | Rule                                                   |
|-----------------|--------------------------------------------------------|
| `homogeneous`   | Excess discordance < −0.01 (less mixed than expected)  |
| `baseline`      | Excess discordance within ±0.01                        |
| `mixed`         | Excess discordance > +0.01 (more mixed than expected)  |
| `not available` | Singleton clusters or missing expected discordance     |

The ±0.01 baseline band matches the threshold used for the Part 2 demographic
mixing categories.

### 3.4 SIMD domain quintiles

For each of the seven SIMD domains, cluster-level mean domain ranks are computed
and divided into quintiles using the same rank-based approach as overall SIMD.
Quintile 1 is most deprived in each domain.

| Domain       | Rank column                    |
|--------------|-------------------------------|
| overall      | `dz_simd_rank`                |
| income       | `dz_simd_income_rank`         |
| employment   | `dz_simd_employment_rank`     |
| education    | `dz_simd_education_rank`      |
| health       | `dz_simd_health_rank`         |
| access       | `dz_simd_access_rank`         |
| crime        | `dz_simd_crime_rank`          |
| housing      | `dz_simd_housing_rank`        |

### 3.5 Age groups

Age-stratified summaries use JCVI rollout-informed groups that approximate
Scottish vaccination eligibility priority:

```text
00-14  15-19  20-29  30-39  40-49
50-54  55-59  60-64  65-69  70-74  75+
```

These groups cannot exactly distinguish `12–15`, `16–17`, `18–19`, or `80+`
because those boundaries are not separately identifiable in the source five-year
`age_band` field.

---

## 4. Cluster Category Rules

`categorise_clusters.py` assigns all category labels. All thresholds are
estimated among **non-singleton clusters only** (cluster_size > 1). The singleton
filter is applied before threshold estimation, so categories describe variation
within the multi-case cluster population.

### 4.1 Size categories

| Category         | Rule                        | Basis                         |
|------------------|-----------------------------|-------------------------------|
| `small/moderate` | `cluster_size < 13`         | Non-singleton 90th percentile |
| `large`          | `13 ≤ cluster_size < 74`    | Non-singleton 90th/99th       |
| `very large`     | `cluster_size ≥ 74`         | Non-singleton 99th percentile |

### 4.2 Geographic dispersion categories

| Category                  | Rule                              | Basis                         |
|---------------------------|-----------------------------------|-------------------------------|
| `low/moderate dispersion` | `cluster_n_datazones < 11`        | Non-singleton 90th percentile |
| `large dispersion`        | `11 ≤ cluster_n_datazones < 61`   | Non-singleton 90th/99th       |
| `very large dispersion`   | `cluster_n_datazones ≥ 61`        | Non-singleton 99th percentile |

### 4.3 SIMD quintile categories

Overall SIMD quintiles are assigned from mean cluster SIMD rank, dividing the
national rank range (1–6,976) into five equal-width bands. Quintile 1 is most
deprived.

### 4.4 Demographic mixing categories

The four demographic mixing dimensions use the Part 1 excess-discordance scores:

| Dimension     | Source column               | Category rules          |
|---------------|-----------------------------|-------------------------|
| SIMD          | `simd_excess_discordance`   | < −0.01 / ±0.01 / > +0.01 |
| Age           | `age_excess_discordance`    | < −0.01 / ±0.01 / > +0.01 |
| Sex           | `sex_excess_discordance`    | < −0.01 / ±0.01 / > +0.01 |
| Joint profile | `profile_excess_discordance`| < −0.01 / ±0.01 / > +0.01 |

Categories: `less mix`, `baseline`, `more mix`, `not available`.

---

## 5. Current Descriptive Results

### 5.1 Overall vaccination summary

| Measure                                              | Value      |
|------------------------------------------------------|------------|
| Total clusters                                       | 193,160    |
| Non-singleton clusters                               | 84,080     |
| Singleton clusters                                   | 109,080    |
| Mean cluster proportion vaccinated (all clusters)    | 68.9%      |
| Mean cluster proportion vaccinated (non-singleton)   | 65.7%      |
| Mean cluster proportion boosted (all members)        | 35.0%      |
| Median days since vaccination (vaccinated clusters)  | 120 days   |

### 5.2 Wave distribution

| Wave    | Total clusters |
|---------|---------------:|
| B.1.177 |          4,625 |
| Alpha   |         12,112 |
| Delta   |         74,272 |
| BA.1    |         32,928 |
| BA.2    |         38,899 |
| BA.4    |          2,669 |
| BA.5    |         16,440 |
| BQ.1    |          3,326 |
| XBB     |            509 |
| Other   |          7,380 |

### 5.3 Non-singleton category distribution

| Category variable                    | Category              | Clusters | Fraction |
|--------------------------------------|-----------------------|---------:|---------:|
| Cluster size                         | small/moderate        |  75,404  |   89.7%  |
| Cluster size                         | large                 |   7,810  |    9.3%  |
| Cluster size                         | very large            |     853  |    1.0%  |
| Geographic dispersion                | low/moderate          |  75,505  |   89.8%  |
| Geographic dispersion                | large                 |   7,706  |    9.2%  |
| Geographic dispersion                | very large            |     856  |    1.0%  |
| Overall SIMD quintile                | Q1 (most deprived)    |   7,517  |    8.9%  |
| Overall SIMD quintile                | Q2                    |  21,095  |   25.1%  |
| Overall SIMD quintile                | Q3                    |  31,602  |   37.6%  |
| Overall SIMD quintile                | Q4                    |  17,536  |   20.9%  |
| Overall SIMD quintile                | Q5 (least deprived)   |   6,317  |    7.5%  |
| Vaccination-status mixing            | homogeneous           |  44,837  |   53.3%  |
| Vaccination-status mixing            | baseline              |   5,980  |    7.1%  |
| Vaccination-status mixing            | mixed                 |  33,263  |   39.6%  |
| SIMD demographic mixing              | less mix              |  45,773  |   54.4%  |
| SIMD demographic mixing              | baseline              |   5,697  |    6.8%  |
| SIMD demographic mixing              | more mix              |  32,597  |   38.8%  |
| Age demographic mixing               | less mix              |  32,099  |   38.2%  |
| Age demographic mixing               | baseline              |   6,432  |    7.7%  |
| Age demographic mixing               | more mix              |  45,536  |   54.2%  |
| Sex demographic mixing               | less mix              |  29,033  |   34.5%  |
| Sex demographic mixing               | baseline              |   7,456  |    8.9%  |
| Sex demographic mixing               | more mix              |  47,578  |   56.6%  |
| Joint profile mixing                 | less mix              |  12,541  |   14.9%  |
| Joint profile mixing                 | baseline              |  63,636  |   75.7%  |
| Joint profile mixing                 | more mix              |   7,890  |    9.4%  |

### 5.4 Vaccination profile by wave (non-singleton clusters, %)

| Wave    | None vaccinated | Mixed vaccination | All vaccinated |
|---------|----------------:|------------------:|---------------:|
| B.1.177 |           89.5% |             10.0% |           0.5% |
| Alpha   |           41.5% |             56.2% |           2.3% |
| Delta   |           10.3% |             65.3% |          24.4% |
| BA.1    |            4.4% |             49.9% |          45.7% |
| BA.2    |            2.6% |             41.0% |          56.4% |
| BA.4    |            2.2% |             28.4% |          69.4% |
| BA.5    |            1.8% |             27.1% |          71.0% |
| BQ.1    |            0.6% |             20.5% |          78.9% |
| XBB     |            0.0% |             19.4% |          80.6% |

The proportion of all-vaccinated clusters rose monotonically across waves,
reaching 80.6% in XBB. Mixed-vaccination clusters peaked in Delta (65.3%) as
the rollout reached the adult population. None-vaccinated clusters fell from
89.5% in B.1.177 to essentially zero by BQ.1 and XBB.

### 5.5 Vaccination-status mixing by wave (non-singleton clusters, %)

| Wave    | Homogeneous | Baseline | Mixed |
|---------|------------:|---------:|------:|
| B.1.177 |       36.7% |    53.9% |  9.5% |
| Alpha   |       57.5% |     2.6% | 39.9% |
| Delta   |       46.9% |     4.0% | 49.1% |
| BA.1    |       57.4% |     2.0% | 40.6% |
| BA.2    |       62.2% |     3.8% | 34.0% |
| BA.4    |       63.2% |    13.0% | 23.8% |
| BA.5    |       66.3% |    10.3% | 23.5% |
| BQ.1    |       51.7% |    29.7% | 18.6% |
| XBB     |       29.5% |    56.6% | 14.0% |

The B.1.177 baseline reflects a pre-rollout era: most clusters fell at baseline
because almost no one was vaccinated, so expected discordance was near zero.
Homogeneous mixing dominated once vaccine coverage was high (BA.2 onwards),
while mixed vaccination-status clusters were proportionally largest during
Delta and BA.1, when the rollout was most heterogeneous across the adult
population.

### 5.6 Dose recency by wave

| Wave    | Non-singleton clusters | Median days since last dose |
|---------|----------------------:|----------------------------:|
| B.1.177 |                 2,098 |                       13    |
| Alpha   |                 6,898 |                       28    |
| Delta   |                33,974 |                      108    |
| BA.1    |                14,062 |                      100    |
| BA.2    |                16,986 |                      143    |
| BA.4    |                   894 |                      208    |
| BA.5    |                 5,449 |                      215    |
| BQ.1    |                   963 |                      214    |
| XBB     |                   129 |                      218    |

Median days since last dose among vaccinated cluster members increased from
roughly two weeks in B.1.177 to over seven months in the late Omicron
subwaves. The broadly similar values for BA.4–XBB suggest that dose recency
plateaued once most cases were occurring far beyond their most recent booster
dose.

---

## 6. Outputs

### 6.1 Cluster characterisation tables (part2/tables/)

| File                                                 | Contents                                                                      |
|------------------------------------------------------|-------------------------------------------------------------------------------|
| `vaccination_cluster_table.csv`                      | Full cluster-level vaccination and category variables                         |
| `vaccination_case_weekly_summary.csv`                | Weekly vaccinated-case proportions by stratum type and stratum                |
| `vaccination_case_weekly_simd_domain_summary.csv`    | Weekly vaccinated-case proportions by SIMD domain                             |
| `vaccination_cluster_wave_category_summary.csv`      | Wave × cluster category summaries (means, counts, vaccination metrics)        |
| `vaccination_cluster_wave_simd_domain_summary.csv`   | Wave × SIMD domain × SIMD quintile summaries                                  |
| `vaccination_cluster_weekly_category_summary.csv`    | Weekly cluster category fractions and vaccination means                       |
| `vaccination_descriptives.csv`                       | Overall descriptive counts and proportions                                    |
| `vaccination_key_questions.csv`                      | Priority research questions with suggested outputs and caveats                |

### 6.2 Cluster category tables (part2/tables/)

| File                                  | Contents                                                     |
|---------------------------------------|--------------------------------------------------------------|
| `cluster_categories.csv`             | Row-level category table (cluster_id + all category columns) |
| `cluster_category_summary.csv`       | One-way category counts and fractions                        |
| `cluster_category_thresholds.csv`    | Threshold and rule audit table                               |
| `cluster_category_combinations.csv`  | Cross-category cluster counts                                |

### 6.3 Supplementary question tables (part2/tables/)

These tables are generated by `part2/supplementary_questions.py` and answer
secondary interpretation questions rather than the core Part 2 figure set.

| File | Contents |
|------|----------|
| `supp_vaccination_missingness_summary.csv` | Vaccination-status completeness by wave, age group, sex, SIMD quintile, and sequencing-fraction quartile |
| `supp_vaccination_profile_cluster_mixing_summary.csv` | Cluster structure and Part 1 demographic-mixing summaries by wave and cluster vaccination profile |
| `supp_vaccination_mixing_demographic_summary.csv` | Cluster structure and demographic-mixing summaries by wave and vaccination-status mixing category |
| `supp_simd_domain_vaccination_by_quintile.csv` | Booster, dose-recency, and vaccination summaries by wave, SIMD domain, and quintile |
| `supp_simd_domain_vaccination_gradients.csv` | Wave-specific Q5 minus Q1 vaccination, booster, and dose-recency gradients for each SIMD domain |

### 6.4 Parquet caches (part2/cache/)

| File                                    | Contents                                                            |
|-----------------------------------------|---------------------------------------------------------------------|
| `vaccination_cluster_table.parquet`     | Full vaccination cluster table (all 193,160 clusters)               |
| `cluster_categories.parquet`           | Non-singleton clusters with Part 1 mixing scores plus Part 2 categories |

### 6.5 Manuscript figures (part2/manuscript/figures/)

Ten figures generated by `part2/manuscript/make_figures.py`. See `README.md` for
the full figure list. Figures are written as PDF, PNG, and LZW-compressed TIFF
at 600 dpi.

---

## 7. Reproducibility

Run the full Part 2 pipeline from the repository root:

```bash
# Step 1 — assign cluster categories (reads Part 1 cache)
conda run -n PhD python part2/categorise_clusters.py

# Step 2 — build vaccination cluster characterisation tables and caches
conda run -n PhD python part2/cluster_characterisation.py

# Step 3 — generate manuscript figures
conda run -n PhD python part2/manuscript/make_figures.py

# Optional — generate supplementary question tables
conda run -n PhD python part2/supplementary_questions.py
```

Important defaults:

| Script                       | Argument                   | Default                                   |
|------------------------------|----------------------------|-------------------------------------------|
| `cluster_characterisation.py`| `--qc`                     | `good`                                    |
| `cluster_characterisation.py`| `--resolution`             | `0.3`                                     |
| `categorise_clusters.py`     | `--large-size-min`         | ceil(non-singleton 90th percentile) = 13  |
| `categorise_clusters.py`     | `--very-large-size-min`    | ceil(non-singleton 99th percentile) = 74  |
| `categorise_clusters.py`     | `--mixing-baseline-band`   | `0.01`                                    |
| `make_figures.py`            | `--tables-dir`             | `part2/tables`                            |
| `make_figures.py`            | `--cache-dir`              | `part2/cache`                             |
| `make_figures.py`            | `--out-dir`                | `part2/manuscript/figures`                |

Category thresholds can be overridden at the command line:

```bash
conda run -n PhD python part2/categorise_clusters.py \
  --large-size-min 10 \
  --very-large-size-min 50 \
  --mixing-baseline-band 0.02
```

---

## 8. Caveats

1. Vaccination status among sequenced cases is **not** a random sample.
   Vaccinated people were tested and sequenced differently over the rollout
   period, and sequencing intensity varied by region, wave, and health-board
   capacity.
2. Vaccination-status mixing is measured relative to the sampled
   `window × lineage` stratum, not relative to the population contact matrix.
   It is a descriptive benchmark, not a causal counterfactual.
3. Dose recency is a proxy for rollout phase as much as for individual immune
   waning. Days since last dose in later waves largely reflects calendar time
   since booster campaigns rather than individual decisions.
4. SIMD is area-based and may not represent individual socioeconomic position.
5. JCVI age groups are approximate. The source five-year `age_band` field cannot
   separately identify `12–15`, `16–17`, or `18–19`.
6. BA.4 and XBB have small non-singleton cluster counts (894 and 129
   respectively) and should be interpreted with caution.
7. All Part 2 analyses share the same primary-resolution (0.3) cluster table
   as Part 1 to ensure internal consistency.
