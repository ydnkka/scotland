# Part 1 Analysis Documentation

## Socioeconomic deprivation, surveillance conditions, and SARS-CoV-2 genomic cluster structure in Scotland

This document describes the analyses implemented for Part 1 of the Scotland SARS-CoV-2 genomic clustering project. It is written as an implementation record: what data were used, how the cluster-level variables were constructed, what models were fitted, what outputs were produced, and what conclusions can be drawn.

The motivating research question was:

> After accounting for lineage, calendar time, Leiden resolution, local incidence, and sequencing intensity, are socioeconomic deprivation and local surveillance conditions associated with larger, longer-lasting, or more geographically dispersed SARS-CoV-2 genomic clusters in Scotland?

During the analysis, this question was expanded to include cluster composition:

> Do socioeconomic deprivation, SIMD domains, and local surveillance conditions shape who mixes with whom inside genomic clusters, in terms of socioeconomic and demographic composition?

The final Part 1 analysis therefore has two linked parts:

1. Cluster magnitude: whether clusters are larger, longer lasting, or more geographically dispersed.
2. Cluster composition: whether clusters are more or less mixed by SIMD, SIMD domains, age, sex, and joint demographic profiles than expected.

## 1. Source Data

The source dataset was:

```text
data/processed/scotland_clustering_analysis_dataset.parquet
```

This parquet is produced by the clustering pipeline described in:

```text
method/PIPELINE.md
method/05_consolidate.py
```

The source table is sequence-level but repeated across analysis windows and Leiden resolutions. Each row represents a sequenced case in a particular sliding time window and clustering resolution.

Key features of the source dataset:

| Feature | Value used in Part 1 |
|---|---:|
| Total rows in processed parquet | 7,287,320 |
| Rows retained after `nextclade_qc == "good"` | 6,314,776 |
| Cluster-level observations after aggregation | 2,482,667 |
| Non-singleton clusters used in mixing regressions | 767,494 |
| Sliding windows | 134 |
| Pango lineages represented in cluster outcome models | 788 |
| Leiden resolutions | 8 (`0.1` to `0.8`) |

Only sequences with `nextclade_qc == "good"` were used in the implemented analyses.

## 2. Clustering Context

Clusters were inferred before Part 1 analysis using the EpiLink/Leiden pipeline:

1. Sequences were grouped by sliding 3-week window and Pango lineage.
2. Pairwise TN93 genetic distances were computed within each group.
3. Distances were converted into EpiLink compatibility weights.
4. Weighted graphs were clustered using Leiden community detection.
5. The procedure was repeated across Leiden resolutions `0.1` to `0.8`.
6. Cluster assignments were merged with sequence metadata, SIMD data, testing data, vaccination data, and local surveillance variables.

The main cluster identifier was `cluster_id`, which encodes the window, lineage, Leiden resolution, and community assignment.

## 3. General Analytical Approach

The first important decision was to analyse clusters, not individual rows. Outcomes such as size, duration, geographic spread, and mixing are cluster-level properties. Treating every sequence row as an independent observation would over-weight larger clusters and misrepresent the unit of analysis.

For each unique `cluster_id`, sequence-level metadata were collapsed to one cluster-level row. Cluster-level averages were used for contextual variables such as mean deprivation, mean local incidence, and mean local sequencing coverage.

The analyses adjusted for:

| Adjustment | Implementation |
|---|---|
| Pango lineage | Fixed effects for `pango_lineage` |
| Calendar time | Fixed effects for `window_id` |
| Leiden resolution | Fixed effects for `resolution_label` |
| Local incidence | Standardised log cumulative local incidence |
| Local sequencing intensity | Standardised logit cumulative datazone sequencing fraction |
| Window-level sequencing intensity | Standardised logit window sequencing proportion |
| Local surveillance/testing pressure | Standardised logit local 7-day test positivity |
| Cluster size | Included in composition/mixing models and size-adjusted duration/spread models |

Continuous predictors were transformed and standardised so that coefficients correspond to a 1 standard deviation increase.

For SIMD and SIMD domains, ranks were negated before standardisation, so that higher values mean greater deprivation.

## 4. Cluster Outcome Analysis

Implemented in:

```text
part1/cluster_outcome_models.py
```

Outputs:

```text
part1/tables/cluster_outcome_model_results.csv
part1/tables/cluster_outcome_model_diagnostics.csv
part1/tables/cluster_outcome_descriptives.csv
part1/tables/cluster_outcome_covariate_scaling.csv
part1/figures/cluster_outcome_model_effects.png
part1/figures/cluster_outcome_model_effects.pdf
```

### 4.1 Outcomes

Three primary cluster outcomes were analysed:

| Outcome | Definition | Modelled variable |
|---|---|---|
| Cluster size | Number of unique sequences in the cluster | `log(cluster_size)` |
| Duration | Days between first and last collection date | `log(cluster_duration_days + 1)` |
| Geographic dispersion | Number of distinct datazones in the cluster | `log(cluster_n_datazones)` |

Two additional size-adjusted models were fitted:

| Model | Purpose |
|---|---|
| Duration, size-adjusted | Tests whether clusters are longer after accounting for their size |
| Geographic dispersion, size-adjusted | Tests whether clusters span more datazones after accounting for their size |

### 4.2 Model form

The implemented model was a log-linear fixed-effect model:

```text
log(outcome_i) =
    deprivation_i
  + local_incidence_i
  + local_sequencing_fraction_i
  + window_sequencing_fraction_i
  + test_positivity_i
  + lineage fixed effects
  + window fixed effects
  + Leiden resolution fixed effects
  + error_i
```

For the size-adjusted duration and geographic dispersion models, `log_cluster_size_z` was added.

Models were solved with least squares on the log-transformed outcome. Standard errors were clustered by `window_id`. Exponentiated coefficients are interpreted as adjusted geometric mean ratios.

### 4.3 Descriptive results

After QC filtering and aggregation:

| Quantity | Value |
|---|---:|
| Cluster rows | 2,482,667 |
| Singleton cluster fraction | 69.1 percent |
| Median cluster size | 1 |
| 90th percentile cluster size | 3 |
| 99th percentile cluster size | 22 |
| Maximum cluster size | 4,930 |
| Median duration | 0 days |
| 90th percentile duration | 4 days |
| Maximum duration | 20 days |
| Median number of datazones | 1 |
| 90th percentile number of datazones | 3 |
| Maximum number of datazones | 3,142 |

### 4.4 Main outcome results

For the three primary outcome models:

| Predictor | Cluster size GMR | Duration GMR | Geographic dispersion GMR |
|---|---:|---:|---:|
| Mean SIMD deprivation | 0.996 | 0.999 | 1.000 |
| Local cumulative incidence | 1.043 | 1.035 | 1.046 |
| Local sequencing fraction | 1.005 | 1.005 | 0.999 |
| Window sequencing proportion | 1.190 | 1.072 | 1.174 |
| Local test positivity | 1.076 | 1.053 | 1.054 |

Interpretation:

- Overall SIMD deprivation was not meaningfully associated with larger, longer-lasting, or more geographically dispersed clusters.
- The only primary deprivation association was a very small negative association with cluster size.
- Local incidence, test positivity, and window-level sequencing proportion were much stronger predictors of cluster magnitude.
- Local cumulative sequencing fraction was close to null in the primary cluster magnitude models.

In size-adjusted models, deprivation was associated with very small increases in duration and geographic dispersion, but these effects were tiny in magnitude and likely less important than cluster size itself.

## 5. Socioeconomic and Demographic Mixing Analysis

Implemented in:

```text
part1/cluster_mixing_analysis.py
```

Outputs:

```text
part1/tables/cluster_mixing_model_results.csv
part1/tables/cluster_mixing_model_diagnostics.csv
part1/tables/cluster_mixing_descriptives.csv
part1/tables/cluster_mixing_covariate_scaling.csv
part1/figures/cluster_mixing_model_effects.png
part1/figures/cluster_mixing_model_effects.pdf
```

### 5.1 Rationale

The outcome models ask whether deprived areas generate bigger clusters. That does not fully answer whether transmission structure is socially patterned. A cluster may be the same size across deprivation contexts but differ in who is connected to whom.

The mixing analysis therefore asks whether cases within the same genomic cluster are more or less mixed by socioeconomic or demographic category than expected.

### 5.2 Pairwise discordance

For a categorical variable such as SIMD quintile, within-cluster mixing was measured as pairwise discordance:

```text
discordance = probability that two different cases drawn from the same cluster
              belong to different categories
```

For a cluster with category counts `n_k` and total valid cases `n`, this was calculated as:

```text
discordance = 1 - sum_k n_k * (n_k - 1) / [n * (n - 1)]
```

This was computed for:

| Mixing outcome | Category variable |
|---|---|
| SIMD mixing | `dz_simd_quintile` |
| Age mixing | `age_band` |
| Sex mixing | `sex` |
| Joint socioeconomic-demographic profile mixing | `dz_simd_quintile + age_band + sex` |

Singleton clusters were excluded from mixing models because pairwise mixing is undefined for a cluster of size 1.

### 5.3 Expected mixing

Observed discordance was compared with the expected discordance among sampled cases from the same:

```text
window_id x pango_lineage x Leiden resolution
```

The main mixing outcome was:

```text
excess discordance = observed cluster discordance - expected stratum discordance
```

Interpretation:

| Value | Meaning |
|---|---|
| `0` | Cluster is as mixed as expected for that lineage, time window, and resolution |
| Positive | Cluster is more mixed than expected |
| Negative | Cluster is more homogeneous or assortative than expected |

### 5.4 Mixing model

The regression model was:

```text
excess_mixing_i =
    mean_SIMD_deprivation_i
  + local_incidence_i
  + local_sequencing_fraction_i
  + window_sequencing_fraction_i
  + test_positivity_i
  + log_cluster_size_i
  + lineage fixed effects
  + window fixed effects
  + Leiden resolution fixed effects
  + error_i
```

Coefficients are reported in percentage points.

### 5.5 Main mixing results

For overall SIMD deprivation:

| Mixing outcome | Effect per 1 SD higher deprivation |
|---|---:|
| SIMD quintile mixing | -0.69 percentage points |
| Age-band mixing | +1.55 percentage points |
| Sex mixing | -0.70 percentage points |
| Joint SIMD-age-sex profile mixing | +0.51 percentage points |

Interpretation:

- More deprived cluster composition was associated with slightly less SIMD mixing, meaning clusters were somewhat more SIMD-homogeneous.
- More deprived cluster composition was associated with more age mixing.
- More deprived cluster composition was associated with slightly less sex mixing.
- The effects were statistically precise because of the large dataset, but most were modest in magnitude.

Other covariates were also important:

- Higher local incidence was associated with more SIMD, age, and joint profile mixing.
- Higher window sequencing proportion was associated with more mixing.
- Higher test positivity was associated with less SIMD mixing but somewhat more sex mixing.
- Larger clusters were generally more mixed, especially for SIMD and age.

## 6. SIMD Domain Analysis

Implemented in:

```text
part1/simd_domain_analysis.py
```

Outputs:

```text
part1/tables/simd_domain_outcome_model_results.csv
part1/tables/simd_domain_outcome_model_diagnostics.csv
part1/tables/simd_domain_mixing_model_results.csv
part1/tables/simd_domain_mixing_model_diagnostics.csv
part1/tables/simd_domain_covariate_scaling.csv
part1/figures/simd_domain_outcome_effects.png
part1/figures/simd_domain_outcome_effects.pdf
part1/figures/simd_domain_mixing_effects.png
part1/figures/simd_domain_mixing_effects.pdf
```

### 6.1 Domains analysed

The following SIMD ranks were analysed:

| Domain | Source column |
|---|---|
| Overall | `dz_simd_rank` |
| Income | `dz_simd_income_rank` |
| Employment | `dz_simd_employment_rank` |
| Education | `dz_simd_education_rank` |
| Health | `dz_simd_health_rank` |
| Access | `dz_simd_access_rank` |
| Crime | `dz_simd_crime_rank` |
| Housing | `dz_simd_housing_rank` |

Each rank was negated and standardised so that higher values mean more domain-specific deprivation.

### 6.2 Domain effects on cluster magnitude

Each SIMD domain was entered one at a time into the cluster outcome models.

Main results:

- Domain effects on cluster size, duration, and geographic dispersion were generally small.
- Housing deprivation showed the clearest negative association with cluster size and duration.
- Access behaved differently from material domains: it showed a small positive association with cluster size and a small negative association with geographic dispersion.
- These domain effects were much smaller than the effects of incidence, test positivity, and sequencing effort.

### 6.3 Domain-quintile mixing

For each domain, ranks were converted into domain-specific quintiles. The analysis then asked whether clusters were more or less mixed by those domain quintiles than expected.

Effect of 1 SD higher domain deprivation on domain-quintile mixing:

| Domain | Effect on domain-quintile excess mixing |
|---|---:|
| Overall | -0.69 percentage points |
| Income | -1.09 percentage points |
| Employment | -0.65 percentage points |
| Education | -0.02 percentage points |
| Health | -0.85 percentage points |
| Access | -1.28 percentage points |
| Crime | -0.26 percentage points |
| Housing | -2.41 percentage points |

Interpretation:

- Housing, access, income, health, and employment deprivation were associated with more domain-specific homogeneity.
- Education and crime showed little clear association with domain-quintile mixing.
- Housing deprivation had the strongest association with reduced domain mixing.

## 7. SIMD Domains and Demographic Mixing

Implemented in:

```text
part1/simd_domain_demographic_mixing.py
```

Outputs:

```text
part1/tables/simd_domain_demographic_mixing_model_results.csv
part1/tables/simd_domain_demographic_mixing_model_diagnostics.csv
part1/tables/simd_domain_demographic_mixing_descriptives.csv
part1/tables/simd_domain_demographic_mixing_covariate_scaling.csv
part1/figures/simd_domain_demographic_mixing_effects.png
part1/figures/simd_domain_demographic_mixing_effects.pdf
```

### 7.1 Purpose

This analysis asked whether specific deprivation domains predicted demographic mixing within clusters.

For each SIMD domain, the exposure was:

```text
cluster mean domain deprivation, standardised
```

The outcomes were:

| Outcome | Meaning |
|---|---|
| Age-band excess mixing | More or less age mixing than expected |
| Sex excess mixing | More or less sex mixing than expected |
| Joint age-sex excess mixing | More or less joint age-sex profile mixing than expected |

Each SIMD domain was modelled one at a time.

### 7.2 Main results

Effect per 1 SD higher domain deprivation:

| Domain | Age mixing | Sex mixing | Age-sex profile mixing |
|---|---:|---:|---:|
| Overall | +1.55 pp | -0.70 pp | +0.71 pp |
| Income | +1.54 pp | -0.64 pp | +0.68 pp |
| Employment | +1.60 pp | -0.61 pp | +0.76 pp |
| Education | +1.58 pp | -0.65 pp | +0.77 pp |
| Health | +1.43 pp | -0.69 pp | +0.66 pp |
| Access | -0.43 pp | +0.14 pp | -0.13 pp |
| Crime | +0.57 pp | -0.62 pp | +0.20 pp |
| Housing | +0.94 pp | -0.51 pp | +0.29 pp |

Interpretation:

- Material deprivation domains were consistently associated with more age mixing.
- The same domains were generally associated with slightly less sex mixing.
- Joint age-sex profile mixing was generally higher in more deprived material-domain clusters.
- Access again behaved differently, suggesting it may capture rurality/remoteness rather than deprivation in the same sense as income, employment, education, health, or housing.

## 8. Wave-Specific Domain-Demographic Mixing

Implemented in:

```text
part1/wave_specific_domain_demographic_mixing.py
```

Outputs:

```text
part1/tables/wave_specific_domain_demographic_mixing_model_results.csv
part1/tables/wave_specific_domain_demographic_mixing_model_diagnostics.csv
part1/tables/wave_specific_domain_demographic_mixing_covariate_scaling.csv
part1/figures/wave_specific_domain_demographic_mixing_effects.png
part1/figures/wave_specific_domain_demographic_mixing_effects.pdf
```

### 8.1 Wave assignment

Pango lineages were mapped to broad wave groups:

| Wave group | Rule |
|---|---|
| B.1.177 | lineage starts with `B.1.177` |
| Alpha | `B.1.1.7` or sublineage |
| Delta | `AY.*` or `B.1.617.2` |
| BA.1 | lineage starts with `BA.1` |
| BA.2 | lineage starts with `BA.2` |
| BA.4 | lineage starts with `BA.4` |
| BA.5 | lineage starts with `BA.5` or `BE.*` |
| BQ.1 | lineage starts with `BQ.` |
| XBB | lineage starts with `XBB` |

Wave-specific models were fitted only for waves with enough non-singleton clusters. XBB had 819 non-singleton clusters and was skipped by the default threshold.

### 8.2 What was modelled

Within each wave, the domain-demographic mixing models from Section 7 were refitted separately. This tested whether the relationship between SIMD domains and demographic mixing was stable across epidemic phases.

### 8.3 Main findings

The domain-demographic mixing effects were strongly wave-dependent.

Examples of strongest domain effects by wave:

| Wave | Mixing outcome | Strongest domain effect |
|---|---|---:|
| B.1.177 | Age mixing | Income, +5.36 pp |
| B.1.177 | Age-sex mixing | Education, +4.05 pp |
| Alpha | Age mixing | Health, +2.35 pp |
| Delta | Age mixing | Education, +1.59 pp |
| BA.1 | Age mixing | Education, +0.95 pp |
| BA.2 | Sex mixing | Education, -1.67 pp |
| BA.5 | Age mixing | Overall SIMD, +3.08 pp |
| BQ.1 | Age mixing | Housing, +4.75 pp |

Interpretation:

- There is no single stable deprivation effect across the pandemic.
- Age-mixing effects were particularly strong in B.1.177 and re-emerged in some later Omicron periods.
- Sex-mixing effects were smaller and often negative, but varied by wave.
- Access behaved differently across waves, reinforcing the idea that it reflects spatial/rural structure as much as material deprivation.

## 9. Observed-vs-Expected Mixing Matrices

Implemented in:

```text
part1/observed_expected_mixing_matrices.py
```

Outputs:

```text
part1/tables/observed_expected_mixing_matrices.csv
part1/figures/observed_expected_simd_matrix_overall.png
part1/figures/observed_expected_simd_matrix_overall.pdf
part1/figures/observed_expected_age_matrix_overall.png
part1/figures/observed_expected_age_matrix_overall.pdf
part1/figures/observed_expected_simd_matrix_by_wave.png
part1/figures/observed_expected_simd_matrix_by_wave.pdf
```

### 9.1 Purpose

The scalar mixing scores indicate whether clusters are more or less mixed overall. The matrices show which groups are mixing more or less than expected.

Matrices were produced for:

| Matrix | Categories |
|---|---|
| SIMD | Quintiles 1 to 5 |
| Age | 5-year age bands from `00-04` to `75+` |

### 9.2 Observed pair probabilities

For each cluster, ordered pairs of cases were counted. For categories `i` and `j`, the number of within-cluster ordered pairs was:

```text
n_i * n_j              if i != j
n_i * (n_i - 1)        if i == j
```

These were summed across clusters.

### 9.3 Expected pair probabilities

Expected pairs were computed within the same:

```text
window_id x pango_lineage x Leiden resolution
```

Expected values were based on the marginal category composition of the sampled stratum, scaled by the number of within-cluster ordered pairs available in that stratum.

The output table reports:

| Column | Meaning |
|---|---|
| `observed_pairs` | Count of observed within-cluster ordered pairs |
| `expected_pairs` | Expected count under stratum-level random mixing |
| `observed_probability` | Observed pair probability |
| `expected_probability` | Expected pair probability |
| `excess_probability` | Observed minus expected probability |
| `excess_percentage_points` | Excess probability in percentage points |
| `observed_expected_ratio` | Observed probability divided by expected probability |

### 9.4 Matrix findings

Overall SIMD matrix:

- The largest excess cell was SIMD Q1-Q1: +0.22 percentage points overall.
- Q1-Q2 and Q2-Q1 were also slightly above expected.
- Pairs involving Q5, especially Q5 with Q1/Q3/Q4/Q5, tended to be below expected.
- Wave-specific matrices showed stronger patterns in some waves, for example Alpha Q1-Q1 was +2.03 percentage points.

Overall age matrix:

- The strongest excess was 20-24 with 20-24: +0.17 percentage points.
- Adjacent young adult combinations such as 20-24 with 25-29 were also above expected.
- This supports a young-adult assortativity signal within clusters.

## 10. Overall Conclusions from Part 1

The main conclusion is:

> Socioeconomic deprivation was not strongly associated with cluster magnitude, but it was associated with cluster composition.

More specifically:

1. Deprivation did not meaningfully predict larger, longer, or more geographically dispersed genomic clusters after adjustment.
2. Local epidemic and surveillance context mattered much more for cluster magnitude.
3. Clusters were not socially or demographically random.
4. SIMD mixing showed mild assortativity, especially around the most deprived quintiles.
5. Material deprivation domains were associated with more age mixing but slightly less sex mixing.
6. Housing, income, health, employment, and access domains were associated with more domain-specific homogeneity.
7. Access behaved differently from material domains, probably because it captures rurality/remoteness.
8. Domain-demographic mixing effects varied substantially by epidemic wave.
9. The strongest scientific story is about transmission composition and social mixing, not simply cluster size.

A concise manuscript-style interpretation would be:

> After adjustment for lineage, calendar time, Leiden resolution, local incidence, sequencing intensity, and test positivity, socioeconomic deprivation was not meaningfully associated with larger, longer-lasting, or more geographically dispersed genomic clusters. However, inferred clusters were socially and demographically structured. Clusters showed modest excess assortativity by SIMD and age, and specific SIMD domains were associated with altered age and sex mixing. These associations varied across epidemic waves, suggesting that the social structure of recent transmission was shaped by changing variant, behavioural, policy, and surveillance contexts rather than by a single fixed deprivation gradient.

## 11. Important Caveats

These analyses are descriptive and associational, not causal.

Important caveats:

- Clusters represent inferred recent genomic linkage, not confirmed direct transmission.
- All results depend on the EpiLink/Leiden clustering framework and the chosen window/resolution design.
- The same sequences appear across multiple windows and Leiden resolutions, so observations are not fully independent.
- Expected mixing is conditional on sampled cases within lineage, window, and resolution. It is not a population contact matrix.
- Sequencing intensity is both a covariate and a source of ascertainment bias.
- SIMD is area-based and may not represent individual socioeconomic position.
- SIMD domains are correlated, so domain-specific models were fitted one at a time.
- Large sample size makes small effects statistically precise; interpretation should focus on effect magnitude and consistency.

## 12. Reproducibility

The main scripts can be rerun from the repository root:

```bash
python3 part1/cluster_outcome_models.py
python3 part1/cluster_mixing_analysis.py
python3 part1/simd_domain_analysis.py
python3 part1/simd_domain_demographic_mixing.py
python3 part1/wave_specific_domain_demographic_mixing.py
python3 part1/observed_expected_mixing_matrices.py
```

The scripts write their outputs into:

```text
part1/tables/
part1/figures/
```

The key summary tables are:

```text
part1/tables/cluster_outcome_model_results.csv
part1/tables/cluster_mixing_model_results.csv
part1/tables/simd_domain_mixing_model_results.csv
part1/tables/simd_domain_demographic_mixing_model_results.csv
part1/tables/wave_specific_domain_demographic_mixing_model_results.csv
part1/tables/observed_expected_mixing_matrices.csv
```

The key figures are:

```text
part1/figures/cluster_outcome_model_effects.png
part1/figures/cluster_mixing_model_effects.png
part1/figures/simd_domain_demographic_mixing_effects.png
part1/figures/wave_specific_domain_demographic_mixing_effects.png
part1/figures/observed_expected_simd_matrix_overall.png
part1/figures/observed_expected_age_matrix_overall.png
part1/figures/observed_expected_simd_matrix_by_wave.png
```

