# Part 1 Main Analysis Documentation

## Socioeconomic deprivation, surveillance conditions, and SARS-CoV-2 genomic cluster structure in Scotland

This document records the current main Part 1 analysis implemented in:

```text
part1/main/main_analysis.py
```

It is written as an implementation record: what data were used, how sequence rows were collapsed to cluster rows, how outcomes and covariates were constructed, what models were fitted, what outputs were written, how the summary plots were produced, and what the current main results show.

The main research questions are:

1. Are socioeconomic deprivation and local surveillance conditions associated with the probability that a genomic cluster exceeds its structural minimum size or geographic spread?
2. Among clusters that exceed those structural minima, are the same covariates associated with the positive count magnitude?
3. Among non-singleton clusters, are socioeconomic deprivation and surveillance conditions associated with within-cluster socioeconomic and demographic mixing beyond what would be expected for the same lineage and analysis window?

The current main analysis deliberately uses one primary Leiden resolution, `0.3`, rather than treating repeated resolutions as independent observations.

## 1. Source Data

The source dataset is read from the repository configuration:

```text
config.yaml -> data.processed.analysis_dataset
```

In the current repository this resolves to:

```text
data/processed/scotland_clustering_analysis_dataset.parquet
```

The source table is sequence-level and includes repeated sequence appearances across sliding windows and Leiden resolutions. The main analysis filters it before aggregation:

| Filter            |  Value |
|-------------------|-------:|
| Nextclade QC      | `good` |
| Leiden resolution |  `0.3` |

Only the columns required for the main analysis are read from parquet. These include identifiers, window metadata, collection dates, geography, lineage, QC status, SIMD, age, sex, incidence, sequencing proportions, test positivity, window sequence counts, and health board.

Current retained analysis counts from `part1/main/tables/main_dataset_descriptives.csv`:

| Quantity                                        |   Value |
|-------------------------------------------------|--------:|
| Sequence rows used                              | 789,347 |
| Sequence rows dropped for missing model fields  |     105 |
| Cluster rows                                    | 193,112 |
| Non-singleton clusters used in mixing models    |  84,067 |
| Sliding windows                                 |     134 |
| Raw Pango lineages                              |     788 |
| Lineage model levels after rare-lineage pooling |     183 |
| Minimum clusters for separate lineage level     |      50 |

## 2. Why This Is the Main Analysis

Earlier exploratory Part 1 scripts fitted log-linear models across all Leiden resolutions. The current main analysis supersedes those scripts for the central manuscript-style estimates because it:

1. Uses a single pre-specified primary resolution (`0.3`) to avoid counting the same underlying cluster structure repeatedly across resolutions.
2. Models structural zeros explicitly with hurdle components instead of relying only on log-transformed outcomes.
3. Uses zero-truncated negative-binomial models for the positive count components, matching the heavy-tailed count nature of cluster size and geographic spread.
4. Keeps mixing models aligned with the same primary-resolution cluster table.
5. Writes all primary tables, diagnostics, cache files, and summary figures under `part1/main/`.

The older scripts remain useful as exploratory or supplementary analyses, but the outputs documented here are the main analysis outputs.

## 3. Cluster Table Construction

The function `build_cluster_table()` collapses sequence rows to one row per `cluster_id`.

### 3.1 Required sequence fields

Rows are dropped if they are missing fields required for model construction:

```text
cluster_id
sequence_id
window_id
window_idx
collection_date
datazone
pango_lineage
dz_simd_rank
dz_cum_incidence_per_capita
dz_cum_prop_sequenced
wn_prop_sequenced
dz_7d_test_positivity
```

### 3.2 Cluster-level aggregations

For each cluster, the script computes:

| Cluster field                     | Construction                                                |
|-----------------------------------|-------------------------------------------------------------|
| `cluster_size`                    | Number of unique sequence IDs                               |
| `cluster_n_datazones`             | Number of distinct datazones represented                    |
| `cluster_start_date`              | Minimum collection date                                     |
| `cluster_end_date`                | Maximum collection date                                     |
| `duration_days`                   | `cluster_end_date - cluster_start_date` in days             |
| `resolution`                      | First resolution value, after filtering to `0.3`            |
| `window_id`                       | First analysis-window ID                                    |
| `window_idx`                      | First numeric window index                                  |
| `wn_mid_date`                     | First window midpoint date                                  |
| `pango_lineage`                   | First lineage                                               |
| `mean_simd_rank`                  | Mean SIMD rank among cluster sequences                      |
| `mean_local_incidence_per_capita` | Mean cumulative local incidence                             |
| `mean_local_seq_fraction`         | Mean cumulative datazone sequencing fraction                |
| `mean_window_seq_fraction`        | Mean window sequencing proportion                           |
| `mean_test_positivity`            | Mean local 7-day test positivity                            |
| `wn_no_sequences`                 | Window-level sequence count                                 |
| `health_board`                    | Health board code from the first cluster row                |
| `index_simd_rank`                 | SIMD rank of the earliest collected sequence in the cluster |

The `index_simd_rank` field supports a sensitivity analysis in which index-case SIMD replaces mean cluster SIMD.

### 3.3 Hurdle and positive-count variables

The cluster descriptors have structural minima. The script therefore creates
binary indicators for exceeding those minima and positive-count variables for
clusters above the minima. In the current main analysis, the cluster-size and
geographic-dispersion variables feed the fitted hurdle models, while the
duration variables are retained for descriptive summaries and the supplementary
outcome-distribution figure because the fixed three-week clustering windows
mechanically constrain the observed span.

| Outcome               | Structural minimum | Binary component   | Positive component                           |
|-----------------------|-------------------:|--------------------|----------------------------------------------|
| Cluster size          |         1 sequence | `cluster_size_gt1` | `cluster_size_excess = cluster_size - 1`     |
| Duration              |             0 days | `duration_gt0`     | `duration_positive_days = duration_days`     |
| Geographic dispersion |         1 datazone | `datazones_gt1`    | `datazones_excess = cluster_n_datazones - 1` |

The positive component is fitted only among clusters where the positive
variable is greater than zero.

### 3.4 Descriptive cluster distributions

Current descriptive results:

| Measure            | Median | 75th percentile | 90th percentile | 99th percentile | Maximum |  Structural-minimum fraction |
|--------------------|-------:|----------------:|----------------:|----------------:|--------:|-----------------------------:|
| Cluster size       |      1 |               3 |               6 |              39 |   2,792 |       56.5 percent singleton |
| Duration           | 0 days |          3 days |          7 days |         12 days | 19 days |        63.1 percent zero-day |
| Distinct datazones |      1 |               2 |               5 |              32 |   2,100 | 61.7 percent single-datazone |

The distributions are strongly right-skewed, especially cluster size and geographic spread. This is the practical reason for using a hurdle model with zero-truncated negative-binomial positive components.

## 4. Covariates and Transformations

The primary covariates are:

| Model term              | Meaning                                          | Source transformation                                   |
|-------------------------|--------------------------------------------------|---------------------------------------------------------|
| `deprivation_z`         | Mean cluster SIMD deprivation                    | `-mean_simd_rank`, standardised                         |
| `local_incidence_z`     | Local cumulative incidence                       | `log1p(mean incidence per capita * 1000)`, standardised |
| `local_seq_fraction_z`  | Local cumulative datazone sequencing fraction    | Clipped logit, standardised                             |
| `window_seq_fraction_z` | Analysis-window sequencing proportion            | Clipped logit, standardised                             |
| `test_positivity_z`     | Local 7-day test positivity                      | Clipped logit, standardised                             |
| `log_cluster_size_z`    | Cluster size for size-adjusted and mixing models | `log(cluster_size)`, standardised                       |

SIMD ranks are negated before standardisation so that higher values mean greater deprivation.

Current standardisation values:

| Standardised column     | Source column               | Source mean | Source SD |
|-------------------------|-----------------------------|------------:|----------:|
| `deprivation_z`         | `deprivation_raw`           |   -3420.328 |  1793.221 |
| `index_deprivation_z`   | `index_deprivation_raw`     |   -3426.554 |  2031.553 |
| `local_incidence_z`     | `local_incidence_log`       |       5.269 |     0.965 |
| `local_seq_fraction_z`  | `local_seq_fraction_logit`  |      -1.714 |     0.859 |
| `window_seq_fraction_z` | `window_seq_fraction_logit` |      -1.727 |     0.652 |
| `test_positivity_z`     | `test_positivity_logit`     |      -2.270 |     0.840 |
| `log_cluster_size_z`    | `log_cluster_size`          |       0.604 |     0.888 |

All reported model coefficients therefore correspond to a 1 SD higher covariate.

## 5. Lineage and Calendar Adjustment

The main analysis adjusts for lineage and calendar time as follows:

| Adjustment                 | Implementation                                                             |
|----------------------------|----------------------------------------------------------------------------|
| Pango lineage              | Fixed-effect dummy variables for `lineage_model`                           |
| Rare lineages              | Lineages with fewer than 50 clusters are pooled into `Other rare lineages` |
| Calendar time              | Cubic B-spline basis over `window_idx`                                     |
| Calendar spline complexity | `calendar_spline_df = 8` by default                                        |

The current model matrix uses 183 lineage model levels and 8 calendar spline columns. One lineage level is the reference level after dummy coding.

The analysis does not include Leiden-resolution fixed effects because the data have already been filtered to one primary resolution.

## 6. Count Models

Count models are fitted by `fit_count_models()`. There are two primary outcomes
plus one size-adjusted positive-count model. Duration is not fitted in the
current main count analysis because the fixed three-week clustering windows
mechanically constrain the observed span.

### 6.1 Model specifications

| Model spec                            | Binary component   | Positive component       | Includes cluster size? |
|---------------------------------------|--------------------|--------------------------|------------------------|
| `cluster_size`                        | `cluster_size_gt1` | `cluster_size_excess`    | No                     |
| `geographic_dispersion`               | `datazones_gt1`    | `datazones_excess`       | No                     |
| `geographic_dispersion_size_adjusted` | Not fitted         | `datazones_excess`       | Yes                    |

For the size-adjusted model, `log_cluster_size_z` is included so that
geographic spread is assessed conditional on reconstructed cluster size.

### 6.2 Binary hurdle component

The binary component asks whether the cluster exceeds its structural minimum.

Model family:

```text
Binomial GLM with logit link
```

Model form:

```text
logit[P(Y_i > minimum)] =
    beta_0
  + beta_1 deprivation_z
  + beta_2 local_incidence_z
  + beta_3 local_seq_fraction_z
  + beta_4 window_seq_fraction_z
  + beta_5 test_positivity_z
  + lineage fixed effects
  + calendar spline terms
```

Binary component coefficients are exponentiated and reported as odds ratios.

### 6.3 Positive zero-truncated count component

The positive component asks how large the count is after the structural minimum has been exceeded.

Model family:

```text
Zero-truncated negative binomial
```

Model form:

```text
log[E(Y_i | Y_i > 0)] =
    beta_0
  + beta_1 deprivation_z
  + beta_2 local_incidence_z
  + beta_3 local_seq_fraction_z
  + beta_4 window_seq_fraction_z
  + beta_5 test_positivity_z
  + optional beta_6 log_cluster_size_z
  + lineage fixed effects
  + calendar spline terms
```

The zero-truncated negative-binomial likelihood is implemented directly in `main_analysis.py`. The optimiser is `scipy.optimize.minimize()` with L-BFGS-B. The dispersion parameter is represented as `log_alpha` and bounded between `-10` and `8`.

Cluster-robust standard errors for the custom ZTNB model are computed with a sandwich estimator:

1. Observation-level scores are computed analytically.
2. The bread is based on a numerical Hessian of the log likelihood at the maximum-likelihood estimate.
3. Scores are summed by the selected clustering variable.
4. The sandwich covariance uses the clustered score cross-product as meat.
5. A finite-sample correction is applied when the number of groups is greater than one.

Positive component coefficients are exponentiated and reported as count ratios.

### 6.4 Standard-error clustering

The default clustering variable is:

```text
window_id
```

This accounts for dependence among clusters from the same sliding analysis window. A command-line sensitivity option allows clustering by health board instead:

```bash
conda run -n PhD python part1/main/main_analysis.py --cluster-by health_board
```

The result table column remains named `std_error_clustered_by_window` for historical compatibility, while the diagnostics table records the actual `cluster_by` value used.

### 6.5 Current count-model diagnostics

Current diagnostics from `main_hurdle_count_model_diagnostics.csv`:

| Outcome                              | Component     | Observations |      Events / mean response | Converged |
|--------------------------------------|---------------|-------------:|----------------------------:|-----------|
| Cluster size                         | Binary hurdle |      193,112 |               84,067 events | Yes       |
| Cluster size                         | Positive ZTNB |       84,067 | Mean positive response 7.09 | Yes       |
| Geographic dispersion                | Binary hurdle |      193,112 |               74,010 events | Yes       |
| Geographic dispersion                | Positive ZTNB |       74,010 | Mean positive response 6.46 | Yes       |
| Geographic dispersion, size-adjusted | Positive ZTNB |       74,010 | Mean positive response 6.46 | Yes       |

The cluster-size positive component has a very large estimated `alpha` and is flagged at the upper bound. This reflects the extreme right tail of cluster-size excess counts and should be kept in mind when interpreting that positive-count component.

The geographic dispersion positive component also has a very large estimated `alpha` (approximately 2,840), close to but technically below the upper bound of `exp(8) ≈ 2,981`. It is not flagged at the upper bound, but given its proximity the same caution applies: the dispersion estimate is unreliable and the heavy right tail of datazone-excess counts should be kept in mind when interpreting that component.

## 7. Main Count Results

### 7.1 Deprivation effects

The main deprivation estimates are:

| Outcome                              | Component     | Ratio per 1 SD higher deprivation | 95 percent CI  | Interpretation                                           |
|--------------------------------------|---------------|----------------------------------:|----------------|----------------------------------------------------------|
| Cluster size                         | Binary hurdle |                             0.971 | 0.960 to 0.983 | Lower odds of being non-singleton                        |
| Cluster size                         | Positive ZTNB |                             0.926 | 0.869 to 0.987 | Fewer additional sequences among non-singletons          |
| Geographic dispersion                | Binary hurdle |                             1.004 | 0.992 to 1.016 | Near null                                                |
| Geographic dispersion                | Positive ZTNB |                             0.851 | 0.792 to 0.915 | Fewer additional datazones among multi-datazone clusters |
| Geographic dispersion, size-adjusted | Positive ZTNB |                             1.027 | 1.019 to 1.035 | Slightly more datazones conditional on size              |

The primary interpretation is that mean cluster deprivation is not associated
with a simple increase in cluster magnitude. It is slightly negatively
associated with the probability and positive magnitude of cluster size, and
with unadjusted positive geographic spread. After conditioning on cluster size,
deprivation is associated with a small increase in geographic dispersion,
suggesting that among clusters of similar reconstructed size, more deprived
clusters may spread across slightly more datazones, but the effect size is
small.

### 7.2 Surveillance and incidence effects

The surveillance and incidence covariates are much stronger predictors of count outcomes than deprivation.

Selected ratios per 1 SD higher covariate:

| Outcome/component              | Local incidence | Local sequencing fraction | Window sequencing proportion | Test positivity |
|--------------------------------|----------------:|--------------------------:|-----------------------------:|----------------:|
| Cluster size hurdle            |           1.173 |                     1.067 |                        1.252 |           1.448 |
| Cluster size positive          |           1.650 |                     3.240 |                        1.314 |           2.649 |
| Geographic dispersion hurdle   |           1.223 |                     1.047 |                        1.170 |           1.314 |
| Geographic dispersion positive |           1.699 |                     2.269 |                        1.274 |           2.999 |

These estimates support the interpretation that local epidemic intensity and
surveillance conditions strongly shape whether clusters are detected as larger
or more geographically dispersed.

## 8. Mixing Outcomes

The main analysis also fits mixing models among non-singleton clusters.

### 8.1 Pairwise discordance

Within each cluster, mixing is measured by pairwise discordance:

```text
discordance = probability that two different cases drawn from the same cluster
              belong to different categories
```

For category counts `n_k` and total valid cases `n`:

```text
discordance = 1 - sum_k n_k * (n_k - 1) / [n * (n - 1)]
```

This is computed for four categorical variables:

| Outcome              | Category variable                   |
|----------------------|-------------------------------------|
| SIMD mixing          | `dz_simd_quintile`                  |
| Age mixing           | `age_band`                          |
| Sex mixing           | `sex`                               |
| Joint profile mixing | `dz_simd_quintile + age_band + sex` |

Singleton clusters are excluded because pairwise discordance is undefined when `n = 1`.

### 8.2 Expected discordance

Observed cluster discordance is compared with expected discordance among sampled cases from the same:

```text
window_id x pango_lineage
```

Resolution is not included in the expected-mixing strata because the main analysis has already filtered to one primary resolution.

The modelled outcome is:

```text
excess discordance = observed cluster discordance - expected stratum discordance
```

Interpretation:

| Value    | Meaning                                                    |
|----------|------------------------------------------------------------|
| `0`      | Cluster is as mixed as expected for its lineage and window |
| Positive | Cluster is more mixed than expected                        |
| Negative | Cluster is more homogeneous or assortative than expected   |

### 8.3 Mixing model form

For each mixing outcome, the model is:

```text
excess_mixing_i =
    beta_0
  + beta_1 deprivation_z
  + beta_2 local_incidence_z
  + beta_3 local_seq_fraction_z
  + beta_4 window_seq_fraction_z
  + beta_5 test_positivity_z
  + beta_6 log_cluster_size_z
  + lineage fixed effects
  + calendar spline terms
  + error_i
```

Models are fitted as OLS regressions with clustered standard errors. Coefficients are reported as percentage-point changes in excess discordance per 1 SD higher covariate.

### 8.4 Current mixing descriptives

Among 84,067 non-singleton clusters:

| Outcome       | Mean excess discordance | Median | 10th percentile | 90th percentile |
|---------------|------------------------:|-------:|----------------:|----------------:|
| SIMD          |                  -0.161 | -0.047 |          -0.798 |           0.203 |
| Age           |                  -0.077 |  0.031 |          -0.333 |           0.083 |
| Sex           |                   0.010 |  0.035 |          -0.499 |           0.503 |
| Joint profile |                  -0.025 |  0.008 |          -0.039 |           0.010 |

The negative mean SIMD and age values indicate that, before covariate adjustment, clusters tend to be less mixed by those dimensions than expected from their lineage-window sampling strata.

## 9. Main Mixing Results

### 9.1 Deprivation effects

Effect per 1 SD higher mean cluster deprivation:

| Mixing outcome                    | Effect in percentage points | 95 percent CI  | Interpretation                  |
|-----------------------------------|----------------------------:|----------------|---------------------------------|
| SIMD quintile mixing              |                    +0.31 pp | -0.17 to +0.80 | Not clearly different from zero |
| Age-band mixing                   |                    +1.66 pp | +1.29 to +2.03 | More age mixing                 |
| Sex mixing                        |                    -0.78 pp | -1.16 to -0.39 | Less sex mixing                 |
| Joint SIMD-age-sex profile mixing |                    +0.48 pp | +0.29 to +0.66 | More joint profile mixing       |

The main mixing result is more nuanced than a simple deprivation-magnitude story. Deprivation does not clearly increase SIMD-quintile mixing in the current main model, but it is associated with more age mixing, less sex mixing, and slightly more joint profile mixing.

### 9.2 Other covariate effects

Selected effects in percentage points:

| Mixing outcome | Local incidence | Local sequencing fraction | Window sequencing proportion | Test positivity | Cluster size |
|----------------|----------------:|--------------------------:|-----------------------------:|----------------:|-------------:|
| SIMD           |           +4.86 |                     -1.22 |                        -3.79 |           -6.65 |        +7.49 |
| Age            |           +1.56 |                     -1.02 |                        -0.76 |           -0.73 |        +2.80 |
| Sex            |           -0.79 |                     -0.36 |                        +0.61 |           +1.15 |        -1.20 |
| Joint profile  |           +1.09 |                     -0.54 |                        -0.30 |           -0.98 |        +0.98 |

Cluster size is a strong predictor of SIMD and age mixing. Test positivity has a strong negative association with SIMD excess mixing and a positive association with sex excess mixing.

## 10. Wave-Specific Cluster Outcome Descriptives

Wave-level summaries of cluster outcomes are written to
`part1/main/tables/main_wave_cluster_outcome_descriptives.csv`. These describe
the distribution of cluster outcomes within each dominant variant wave and
provide context for the main regression results.

### 10.1 Cluster size by wave

| Wave    | Clusters | Singleton fraction | 75th pct | 90th pct | Among non-singletons: median | mean |
|---------|---------:|--------------------|----------:|----------:|-----------------------------:|-----:|
| B.1.177 |   4,621  | 54.6%              |       3.0 |       7.0 |                          2.0 | 6.95 |
| Alpha   |  12,112  | 43.0%              |       5.0 |      15.0 |                          3.0 | 10.04|
| Delta   |  74,272  | 54.3%              |       3.0 |       7.0 |                          2.0 | 6.62 |
| BA.1    |  32,928  | 57.3%              |       3.0 |       6.0 |                          2.0 | 9.41 |
| BA.2    |  38,893  | 56.3%              |       3.0 |       6.0 |                          2.0 | 6.81 |
| BA.4    |   2,669  | 66.5%              |       2.0 |       4.0 |                          1.0 | 4.30 |
| BA.5    |  16,423  | 66.9%              |       2.0 |       4.0 |                          1.0 | 3.77 |
| BQ.1    |   3,314  | 70.9%              |       2.0 |       3.0 |                          1.0 | 2.29 |
| XBB     |     509  | 74.7%              |       2.0 |       3.0 |                          1.0 | 2.86 |

Alpha has the lowest singleton fraction (43.0%) and the largest non-singleton
mean cluster size (10.04 additional sequences), reflecting the concentrated
pre-SL expansion documented in Part 3. Late Omicron subwaves (BA.4, BA.5,
BQ.1, XBB) show markedly higher singleton fractions (66–75%) and smaller
median non-singleton sizes, consistent with faster variant turnover and
lower sequencing depth.

### 10.2 Geographic dispersion by wave

| Wave    | Single-datazone fraction | Among non-singletons: median datazones | 90th pct |
|---------|--------------------------:|----------------------------------------:|---------:|
| B.1.177 |            58.4%          |                                     2.0 |     13.0 |
| Alpha   |            50.7%          |                                     3.0 |     17.0 |
| Delta   |            59.3%          |                                     2.0 |     11.0 |
| BA.1    |            61.8%          |                                     2.0 |      9.0 |
| BA.2    |            62.0%          |                                     2.0 |     10.0 |
| BA.4    |            70.8%          |                                     1.0 |      8.0 |
| BA.5    |            71.3%          |                                     2.0 |      7.9 |
| BQ.1    |            77.1%          |                                     1.0 |      5.0 |
| XBB     |            80.4%          |                                     2.0 |      7.1 |

Alpha has the lowest single-datazone fraction (50.7%) and the widest geographic
spread among non-singleton clusters (median 3.0 datazones, 90th percentile
17.0). The pattern mirrors cluster-size: late Omicron subwaves are markedly
more geographically concentrated.

## 11. Mixing-Predictor Count Models

The main script also fits mixing-predictor count models in which the four
excess-mixing scores (SIMD, age, sex, joint profile) are added alongside the
primary surveillance covariates. These models assess whether within-cluster
demographic mixing is associated with cluster size and geographic dispersion
beyond the primary covariate set. Results are written to:

```text
main_mixing_predictor_hurdle_count_model_results.csv
main_mixing_predictor_hurdle_count_model_diagnostics.csv
```

### 11.1 Key mixing-predictor count results

Selected count-ratio and odds-ratio estimates per 1 SD higher mixing score,
from the positive-count (ZTNB) components among non-singleton clusters
(n = 84,067 for cluster size; n = 74,010 for geographic dispersion):

| Mixing predictor         | Cluster size count ratio | Geographic dispersion count ratio | Geographic dispersion hurdle OR |
|--------------------------|-------------------------:|----------------------------------:|--------------------------------:|
| SIMD excess discordance  |       3.478 (3.245–3.727)|                    3.029 (2.793–3.285)|                  22.107 (18.977–25.753)|
| Age excess discordance   |       1.668 (1.560–1.783)|                    1.966 (1.800–2.148)|                   1.283 (1.246–1.322)|
| Sex excess discordance   |       0.852 (0.764–0.949)|                    1.060 (0.918–1.223)|                   0.770 (0.737–0.805)|
| Joint profile excess     |       0.812 (0.763–0.865)|                    0.718 (0.665–0.775)|                   1.029 (1.007–1.050)|

SIMD excess discordance is the strongest mixing predictor of both cluster size
and geographic dispersion. A 1 SD increase in SIMD excess discordance is
associated with a 3.48-fold increase in additional cluster sequences (among
non-singletons) and an odds ratio of 22.1 for exceeding a single datazone.
Age excess discordance is also a substantial positive predictor of both
outcomes. These results suggest that more socioeconomically and age-mixed
clusters are detected as larger and more geographically dispersed, consistent
with clusters that span more types of settings.

## 12. Sensitivity Options Built Into the Main Script

The command-line interface supports several pre-specified sensitivity runs.

### 12.1 Alternative standard-error clustering

```bash
conda run -n PhD python part1/main/main_analysis.py --cluster-by health_board
```

This clusters the sandwich standard errors by health board instead of window. The default is `window_id`.

### 12.2 Size-offset cluster-size positive model

```bash
conda run -n PhD python part1/main/main_analysis.py --use-size-offset
```

This includes `log(wn_no_sequences)` as an offset in the cluster-size positive-count model. The estimand changes from raw reconstructed cluster size to cluster size relative to the number of sequences available in the analysis window.

### 12.3 Tail winsorisation

```bash
conda run -n PhD python part1/main/main_analysis.py --winsorise-quantile 0.99
```

This caps positive count outcomes at the specified quantile before fitting the ZTNB models. It is intended to test sensitivity to very large right-tail clusters.

### 12.4 Index-case SIMD exposure

```bash
conda run -n PhD python part1/main/main_analysis.py --use-index-simd
```

This replaces mean cluster SIMD deprivation with deprivation for the earliest collected sequence in the cluster.

### 12.5 Approximately non-overlapping windows

```bash
conda run -n PhD python part1/main/main_analysis.py --window-stride 3
```

This keeps only clusters from windows where `window_idx % 3 == 0`. With 3-week windows advanced in 1-week steps, this approximates a non-overlapping window sensitivity.

### 12.6 Separate output directories for sensitivities

Sensitivity runs should use separate output directories to avoid overwriting primary results:

```bash
conda run -n PhD python part1/main/main_analysis.py \
  --cluster-by health_board \
  --tables-dir part1/main/tables_health_board \
  --figures-dir part1/main/figures_health_board \
  --cache-dir part1/main/cache_health_board
```

## 13. Outputs

The main run writes tables to:

```text
part1/main/tables/
```

Primary table outputs:

| File                                      | Contents                                          |
|-------------------------------------------|---------------------------------------------------|
| `main_covariate_scaling.csv`                                  | Means and SDs used for standardisation                           |
| `main_dataset_descriptives.csv`                               | Analysis counts and descriptive outcome summaries                |
| `main_wave_cluster_outcome_descriptives.csv`                  | Per-wave cluster size and geographic dispersion summaries        |
| `main_hurdle_count_model_results.csv`                         | Binary and positive count effect estimates (primary covariates)  |
| `main_hurdle_count_model_diagnostics.csv`                     | Count-model diagnostics                                          |
| `main_mixing_predictor_hurdle_count_model_results.csv`        | Count estimates with mixing scores added as predictors           |
| `main_mixing_predictor_hurdle_count_model_diagnostics.csv`    | Diagnostics for mixing-predictor count models                    |
| `main_mixing_predictor_loglinear_count_model_results.csv`     | Log-linear supplementary estimates with mixing scores            |
| `main_loglinear_count_model_results.csv`                      | Log-linear count estimates (supplementary)                       |
| `main_mixing_model_results.csv`                               | Mixing-model effect estimates                                    |
| `main_mixing_model_diagnostics.csv`                           | Mixing-model diagnostics                                         |
| `main_observed_expected_mixing_matrices.csv`                  | Observed and expected mixing matrices                            |
| `main_simd_domain_hurdle_count_model_results.csv`             | SIMD-domain hurdle count results                                 |
| `main_simd_domain_mixing_model_results.csv`                   | SIMD-domain mixing model results                                 |
| `main_wave_specific_hurdle_count_model_results.csv`           | Wave-specific hurdle count results                               |
| `main_wave_specific_mixing_predictor_hurdle_count_model_results.csv` | Wave-specific hurdle with mixing predictors             |

The main run also writes the cluster cache:

```text
part1/main/cache/main_cluster_table.parquet
```

Figures are written to:

```text
part1/main/figures/
```

Primary summary figures:

```text
main_hurdle_count_effects.png
main_hurdle_count_effects.pdf
main_mixing_effects.png
main_mixing_effects.pdf
```

## 14. Summary Plotting

The summary plots are generated by:

```text
plot_count_effects()
plot_mixing_effects()
```

Both functions now use the project style module:

```text
utils/style.py
```

Specifically, the plotting code:

1. Sets the non-interactive Matplotlib backend through `load_plot_style()`.
2. Imports `utils.style` from the repository root.
3. Uses `style.new_figure()` for publication-style figure sizes and rcParams.
4. Uses the shared `style.SIMD_DOMAIN_PALETTE` through `term_colours()`.
5. Uses `style.save_figure()` to write both PDF and PNG outputs.

The count summary plot:

| Feature          | Implementation                   |
|------------------|----------------------------------|
| Figure file stem | `main_hurdle_count_effects`      |
| Layout           | 2 outcomes x 2 components        |
| X-axis scale     | Log ratio scale                  |
| Reference line   | Ratio = 1                        |
| Major ticks      | `0.8`, `1`, `1.5`, `2`, `3`, `4` |
| Minor ticks      | Suppressed                       |
| Output formats   | PNG and PDF                      |

The mixing summary plot:

| Feature          | Implementation                                                       |
|------------------|----------------------------------------------------------------------|
| Figure file stem | `main_mixing_effects`                                                |
| Layout           | Single horizontal coefficient plot                                   |
| X-axis scale     | Percentage-point difference                                          |
| Reference line   | Difference = 0                                                       |
| Major ticks      | Symmetric ticks every 2 percentage points for the current axis range |
| Output formats   | PNG and PDF                                                          |

The x-axis ticks were explicitly controlled because Matplotlib's defaults on log-ratio and wide percentage-point axes produced less readable tick placement in the summary figures.

## 15. Reproducibility

The main analysis can be rerun from the repository root with:

```bash
conda run -n PhD python part1/main/main_analysis.py
```

Important default arguments:

| Argument                 | Default     |
|--------------------------|-------------|
| `--qc`                   | `good`      |
| `--primary-resolution`   | `0.3`       |
| `--lineage-min-clusters` | `50`        |
| `--calendar-spline-df`   | `8`         |
| `--maxiter`              | `1000`      |
| `--cluster-by`           | `window_id` |
| `--winsorise-quantile`   | `0.0`       |
| `--window-stride`        | `1`         |

To regenerate only the two summary figures from existing CSV results, use:

```bash
conda run -n PhD python -c "
from pathlib import Path
import pandas as pd
from part1.main.main_analysis import plot_count_effects, plot_mixing_effects

root = Path('.')
tables = root / 'part1/main/tables'
figures = root / 'part1/main/figures'
plot_count_effects(
    pd.read_csv(tables / 'main_hurdle_count_model_results.csv'),
    figures / 'main_hurdle_count_effects',
)
plot_mixing_effects(
    pd.read_csv(tables / 'main_mixing_model_results.csv'),
    figures / 'main_mixing_effects',
)
"
```

## 16. Current Main Interpretation

The current main results support three linked conclusions:

1. Deprivation is not associated with larger cluster magnitude in a simple positive direction. In the main count models, higher mean cluster deprivation is associated with slightly lower odds of a non-singleton cluster, fewer additional sequences among non-singletons, and fewer additional datazones among clusters spanning multiple datazones.
2. After conditioning on cluster size, deprivation is associated with a very small increase in positive geographic dispersion. This suggests that size is an important mediator or confounder for interpreting spatial spread.
3. Local epidemic and surveillance context is much more strongly associated with cluster magnitude than deprivation. Local incidence, test positivity, local sequencing fraction, and window sequencing proportion show larger and more consistent effects across count components.

For mixing, the current main results show:

1. Mean cluster deprivation is not clearly associated with SIMD-quintile excess mixing in the primary model.
2. Higher deprivation is associated with more age mixing, less sex mixing, and slightly more joint SIMD-age-sex profile mixing.
3. Cluster size and test positivity are important predictors of mixing outcomes, especially SIMD mixing.

Overall, the scientific story is not simply that deprived areas have bigger genomic clusters. The stronger story is that cluster detection and apparent magnitude are shaped by epidemic intensity and surveillance conditions, while cluster composition shows specific socioeconomic and demographic structure.

## 17. Caveats

Important caveats for interpreting the main analysis:

1. Genomic clusters represent inferred recent linkage, not confirmed direct transmission.
2. The analysis is descriptive and associational, not causal.
3. The primary analysis uses one Leiden resolution to avoid repeated-resolution dependence, but results still depend on the clustering framework and the chosen resolution.
4. Sliding analysis windows induce dependence; default standard errors are clustered by `window_id`, and non-overlapping-window sensitivity runs are available.
5. Expected mixing is conditional on sampled sequences within lineage and window; it is not a population contact matrix.
6. Sequencing intensity is both an adjustment variable and part of the ascertainment process.
7. SIMD is area-based and may not represent individual socioeconomic position.
8. Positive count components, especially cluster-size excess, are heavy-tailed; the cluster-size ZTNB dispersion estimate reaches its upper bound in the current main diagnostics.
9. Large sample size makes small effects statistically precise, so interpretation should emphasise effect size and consistency rather than p-values alone.
