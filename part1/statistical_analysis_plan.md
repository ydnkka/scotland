# Statistical Analysis Plan — Part 1

## Socioeconomic Deprivation, Local Surveillance, and SARS-CoV-2 Genomic Cluster Characteristics in Scotland

**Research question:** After accounting for lineage, calendar time, Leiden resolution, local incidence, and sequencing intensity, are socioeconomic deprivation and local surveillance conditions associated with larger, longer-lasting, or more geographically dispersed SARS-CoV-2 genomic clusters in Scotland?

---

## 1. Study aims

This analysis tests whether two structural features of a cluster's source population — the socioeconomic deprivation of the data zones (DZs) that contributed sequences to the cluster, and the local genomic surveillance intensity at the time of sampling — are independently associated with three cluster-level characteristics: size (number of sequences), duration (days from first to last sampled sequence), and geographic dispersal (number of distinct DZs represented). All models adjust for pango lineage, calendar time, Leiden resolution, and the local epidemiological context (positive-test burden and sequencing fraction) in which each cluster was detected.

---

## 2. Data and unit of analysis

### 2.1 Source dataset

The sequence-level parquet (`scotland_clustering_analysis_dataset.parquet`) contains approximately 7.3 million rows structured as one row per unique sequence × analysis window × Leiden resolution. Three-week analysis windows advance in one-week steps (overlap = two weeks). Eight Leiden resolutions (0.1–0.8) are applied within each window × lineage group.

### 2.2 Aggregation to the cluster level

All three outcomes are properties of whole clusters, not individual sequences. The modelling unit is therefore **one row per unique `cluster_id`**. Each cluster is identified by the composite key `{window_id}|{lineage}|R{resolution}|C{n}` (multi-sequence) or `|S{n}` (singleton). Aggregation proceeds as follows.

**Outcome columns** (`cluster_size`, `cluster_duration_days`, `cluster_n_datazones`) are already computed at cluster level and are identical for all sequences sharing a `cluster_id`; take any one value per cluster (e.g. first row after groupby).

**Exposure and covariate columns** require within-cluster aggregation:

| Variable | Aggregation rule | Rationale |
|---|---|---|
| `dz_simd_rank` (SIMD deprivation) | Mean across all sequences in the cluster | Captures the average socioeconomic context of the entire transmission event, not just the index case; use raw rank (1–6,976) on a continuous scale |
| `dz_cum_prop_sequenced` (DZ surveillance intensity) | Mean across all sequences | Reflects average local sequencing coverage at the time each case was sampled |
| `wn_prop_sequenced` (window-level sequencing fraction) | Constant per window; take the single value for the cluster's window | Measures Scotland-wide sequencing effort within the relevant time window |
| `dz_cum_incidence_per_capita` (cumulative local incidence) | Mean across sequences | Average epidemic burden in the DZs contributing to the cluster |
| `wn_positive_tests` (window-level case burden) | Constant per window | Epidemiological context at time of detection |
| `wn_no_sequences` (window size) | Constant per window | Used as a model offset for cluster size |
| `pango_lineage` | Constant per cluster (lineage is part of the cluster key) | Direct attribute |
| `window_idx` / `wn_mid_date` | Constant per cluster | Calendar time |
| `resolution` | Constant per cluster | Leiden granularity |
| `dz_health_board_code` | Modal value | For random-effect grouping |

### 2.3 Singleton classification

Sequences that appear as singletons in the `cluster_id` field come from two sources:

- **S0 singletons**: lineage × window groups of size 1 that were never processed by TN93 or Leiden (i.e. the whole lineage group had only one sequence in that window). These have `cluster_size = 1` but are not the product of the community-detection algorithm.
- **Sn singletons** (n ≥ 1): sequences that passed through TN93/Leiden but were not assigned to any multi-sequence community. These represent genuine epidemiological isolation within a resolved lineage group.

For the primary analysis, **all singletons are included**. This is necessary because their prevalence is itself informative — areas with high deprivation or poor surveillance may produce more singletons if transmission is either too fragmented or too under-sampled to form detectable clusters. Excluding singletons would induce selection bias towards detected, genetically-linked transmission events.

Singleton status (binary flag: `cluster_size == 1`) is included as a component of the model structure via the hurdle / zero-truncated framework described in Section 6, rather than being filtered out.

---

## 3. Outcomes

### 3.1 Cluster size (`cluster_size`)

**Definition:** Number of unique sequences assigned to the cluster. Minimum possible value is 1 (singletons). There are no structural zeros.

**Distribution:** Highly right-skewed positive integer, overdispersed relative to Poisson. The count excludes any notion of "true" transmission chain size; it reflects detected and sequenced cases within the cluster.

**Transformation / model family:** Zero-truncated negative binomial (ZTNB) regression. Because the minimum value is 1 (no cluster can have zero sequences), standard negative binomial with `log(wn_no_sequences)` as an offset is appropriate, but the zero-truncated parameterisation avoids estimating the probability mass at zero. In practice, given the very large proportion of size-1 clusters (singletons), a hurdle negative binomial may fit better: the binomial hurdle part models whether a cluster exceeds size 1, and the zero-truncated NB part models the size conditional on exceeding 1.

**Offset:** `log(wn_no_sequences)` — the number of sequences available in the window — is included as an offset so that the modelled quantity is the cluster's size *relative to the sampling pool*, not its raw count.

### 3.2 Cluster duration (`cluster_duration_days`)

**Definition:** Days between the earliest and latest sequence collection dates within the cluster. Singletons and clusters where all sequences share a collection date have `cluster_duration_days = 0`.

**Distribution:** Non-negative integer with an excess of zeros. The zero mass has two components: true singletons (structurally zero) and small multi-sequence clusters sampled on the same day (informative zero).

**Model family:** Hurdle negative binomial (also called a two-part model):
1. **Binary part (logistic regression):** Models `P(duration > 0 | covariates)`, i.e. whether any temporal spread is observed.
2. **Count part (truncated negative binomial):** Models `E(duration | duration > 0, covariates)`.

Both parts are fitted simultaneously. Alternatively, a zero-inflated negative binomial may be considered if the zero-generating process is thought to be distinct from the count-generating process; model selection between hurdle and ZINB will be guided by Vuong's test and AIC.

### 3.3 Geographic dispersal (`cluster_n_datazones`)

**Definition:** Number of distinct Scottish 2011 Data Zones (DZs) represented in the cluster. Minimum value is 1.

**Primary metric:** `cluster_n_datazones` (integer count, min 1, no structural zeros).

**Supplementary metric (optional computation):** The maximum pairwise Euclidean distance between DZ population-weighted centroids (`dz_xcoord`, `dz_ycoord`) among sequences in the cluster, in kilometres (OSGB36 / EPSG:27700 coordinates; divide by 1,000). This provides a continuous measure of the geographic footprint that is not dependent on DZ boundaries or density. This metric must be computed during aggregation from the sequence-level coordinate columns before collapsing to the cluster-level dataset.

**Model family:** Zero-truncated negative binomial for `cluster_n_datazones` (minimum 1, strong right skew). If using maximum pairwise distance (continuous, right-skewed, zero when cluster_size = 1), a two-part model analogous to cluster duration is more appropriate: logistic for `P(distance > 0)`, gamma or log-normal for `E(distance | distance > 0)`.

---

## 4. Primary exposures

### 4.1 Socioeconomic deprivation (SIMD)

**Variable:** Mean `dz_simd_rank` across sequences in the cluster (continuous, range 1–6,976; lower values indicate greater deprivation).

**Primary specification:** Continuous, linearly entered as `dz_simd_rank_mean` in all three models. A lower mean rank therefore corresponds to a cluster from more deprived areas. To aid interpretation of coefficients, the variable is rescaled: **`simd_rank_std = (mean_dz_simd_rank − μ) / σ`** (standardised to mean 0, SD 1 across clusters). A one-unit increase in `simd_rank_std` represents a one-standard-deviation shift towards less deprivation.

**Alternative specification (sensitivity):** SIMD quintile of the index case (the sequence with the earliest `collection_date` in the cluster), treated as an ordinal five-level factor or as a linear trend 1–5. This tests whether the result is driven by the full distribution of deprivation within clusters or by the characteristic deprivation of the case that likely seeded the cluster.

**Interpretation direction:** If more deprived clusters (lower SIMD rank; more negative `simd_rank_std`) are larger/longer/more dispersed, the coefficient on `simd_rank_std` will be positive (since higher rank = less deprived = smaller expected cluster).

### 4.2 Local surveillance conditions

**Variable:** Mean `dz_cum_prop_sequenced` across sequences in the cluster — the fraction of cumulative positive tests in each DZ that received a linked genome sequence, averaged across the cluster's contributing DZs at the time of each sequence's collection.

**Rationale:** A high `dz_cum_prop_sequenced` indicates that the DZ has been well-covered by genomic surveillance historically. Clusters arising in high-surveillance DZs may appear larger or longer simply because more of their constituent cases were sequenced, not because transmission was more extensive. Conversely, low-surveillance DZs may produce clusters that are truncated in size and spatial extent because many transmission links are invisible.

**Scaling:** Logit-transform or standardise (mean 0, SD 1). Because many DZs early in the pandemic have very low or very high proportions, logit-transformation is preferred if the distribution is strongly bimodal near 0 and 1.

**Secondary surveillance variable:** `wn_prop_sequenced` (window-level) captures Scotland-wide sequencing effort in the time window. This is collinear with calendar time (sequencing intensity varied over the pandemic) and should be used cautiously; include it alongside the DZ-level measure but check variance inflation.

---

## 5. Covariates

All models adjust for the following covariates to answer the question "after accounting for...":

### 5.1 Pango lineage

**Variable:** `pango_lineage` (categorical, hundreds of levels).

**Specification:** Fixed-effect dummy variables for the top *k* lineages by frequency (suggested k = 20, capturing > 90% of sequences in most waves), with all remaining lineages collapsed to an "Other" reference category. Alternatively, `who_voc` (6 levels: Alpha, Delta, Omicron, BA.2, etc., or NaN) can replace lineage to reduce degrees of freedom when the lineage-level fixed effects are not estimable due to sparse data.

**Rationale:** Different lineages have different intrinsic growth rates, within-host viral loads, and test sensitivities that directly affect cluster detection and apparent size independently of socioeconomic context.

### 5.2 Calendar time

**Variable:** `window_idx` (integer, 1 to N) or equivalently `wn_mid_date`.

**Specification:** Restricted cubic spline (RCS) with 4 knots placed at the 5th, 35th, 65th, and 95th percentiles of the window index distribution. This flexibly captures the non-linear trajectory of pandemic dynamics (policy changes, wave structure, vaccine rollout) without constraining the time effect to be linear. An alternative is a piecewise linear specification using the pre-defined wave boundaries in `method/wave_dates.parquet`.

**Interaction with lineage:** Because waves are largely defined by lineage turnover, the spline on time will partially colinear with lineage fixed effects. The model should be checked for variance inflation; if severe, replace the spline with wave fixed effects (wave × lineage interactions if necessary).

### 5.3 Leiden resolution

**Variable:** `resolution` (continuous: 0.1, 0.2, …, 0.8; or an 8-level factor).

**Specification:** 8-level factor (reference: 0.3, a commonly used intermediate resolution). Higher resolution values produce finer-grained communities (more clusters, smaller clusters on average), so resolution is a mechanical determinant of cluster size and must be adjusted for to avoid confounding.

**Note on primary vs. sensitivity analyses:** Because including all 8 resolutions in a single model requires correct handling of the within-window-within-lineage correlation structure (the same sequences contribute to all 8 resolution-specific clusters), the **primary analysis fixes resolution at 0.3** and adjusts for it trivially. The sensitivity analysis repeats the primary models at each of the 8 resolutions separately, and an omnibus sensitivity model includes all resolutions simultaneously with a random intercept for `window_id|lineage` to account for the dependency.

### 5.4 Local incidence

**Variable:** Mean `dz_cum_incidence_per_capita` across sequences in the cluster (cumulative positive tests per capita in each DZ at the time of sampling, averaged across cluster sequences).

**Rationale:** Areas with a higher epidemic burden have more circulating virus and more opportunities for sustained transmission; they will tend to produce larger and longer clusters regardless of deprivation. This variable decouples the effect of deprivation from the general level of virus circulation in the area.

**Scaling:** Right-skewed; log-transform (`log(x + ε)`, where ε is a small positive constant if zeros exist, or `log1p(x)`).

### 5.5 Window-level sequencing intensity

**Variable:** `wn_prop_sequenced` — the fraction of Scotland-wide positive tests that were sequenced in the window.

**Rationale:** High window-level sequencing effort inflates apparent cluster size and geographic coverage because more of the true transmission chain is observed. This is a sampling-effort covariate, not a substantive predictor.

**Scaling:** Logit-transform (values range 0–1).

---

## 6. Model specification

For each of the three outcomes, the primary model is fitted on the cluster-level dataset restricted to `resolution == 0.3`. All models are fitted using maximum likelihood; standard errors are clustered at the `dz_health_board_code` level using a sandwich estimator to account for residual spatial correlation between clusters from the same NHS Health Board region.

### 6.1 Model A — Cluster size

$$
\log(\mu_i) = \beta_0 + \beta_1 \cdot \text{simd\_rank\_std}_i + \beta_2 \cdot \text{logit}(\text{dz\_cum\_prop\_seq}_i) + \mathbf{z}_i^\top \boldsymbol{\gamma} + \log(\text{wn\_no\_sequences}_i)
$$

where $\mu_i = E(\text{cluster\_size}_i)$, $\mathbf{z}_i$ contains the lineage fixed effects, the time spline, and the log-incidence and logit-sequencing-intensity covariates, and the final term is the log-offset. The response distribution is **negative binomial** (NB2 parameterisation), allowing overdispersion parameter $\phi$ to be estimated from the data. If the proportion of size-1 clusters is very high (> 50%), a **hurdle negative binomial** is preferred, with the logistic part modelling $P(\text{size} > 1)$ and the truncated NB part modelling $E(\text{size} \mid \text{size} > 1)$.

### 6.2 Model B — Cluster duration

A **hurdle negative binomial** with two sub-models:

**Part 1 — Binary (logistic):** $\text{logit}[P(\text{duration} > 0)] = \beta_0 + \beta_1 \cdot \text{simd\_rank\_std} + \beta_2 \cdot \text{logit}(\text{dz\_cum\_prop\_seq}) + \mathbf{z}^\top \boldsymbol{\gamma}$

**Part 2 — Count (truncated NB):** $\log(\mu | \text{duration} > 0) = \alpha_0 + \alpha_1 \cdot \text{simd\_rank\_std} + \alpha_2 \cdot \text{logit}(\text{dz\_cum\_prop\_seq}) + \mathbf{z}^\top \boldsymbol{\delta}$

The same covariate vector $\mathbf{z}$ enters both parts (lineage, time spline, incidence, window sequencing fraction). Cluster size (`log(cluster_size)`) is additionally included as a covariate in the count part only, because longer clusters tend to be larger and size is a mechanical correlate of duration.

### 6.3 Model C — Geographic dispersal

**Primary:** Zero-truncated negative binomial for `cluster_n_datazones` (minimum 1):

$$\log(\mu_i) = \beta_0 + \beta_1 \cdot \text{simd\_rank\_std}_i + \beta_2 \cdot \text{logit}(\text{dz\_cum\_prop\_seq}_i) + \beta_3 \cdot \log(\text{cluster\_size}_i) + \mathbf{z}_i^\top \boldsymbol{\gamma}$$

Cluster size is included as a covariate because a larger cluster will span more DZs simply by chance; the coefficient on deprivation thus captures dispersal in excess of what is explained by size.

**Supplementary:** If maximum pairwise distance is computed, a hurdle model is fitted with logistic part for $P(\text{distance} > 0)$ and log-normal or gamma regression for $E(\text{distance} \mid \text{distance} > 0)$, with the same covariate structure.

---

## 7. Model validation and diagnostics

For each fitted model:

- **Rootogram or PIT histogram** to assess distributional fit of the count component.
- **Dispersion test** (ratio of Pearson residuals to degrees of freedom) to confirm the chosen NB family is appropriate.
- **Variance inflation factors (VIF)** for all covariates; flag any VIF > 5, particularly the time spline versus lineage and the two sequencing-intensity variables.
- **Scaled Schoenfeld-style residuals** (for duration models) to check whether the proportional effect of covariates is stable across the range of the response.
- **Cook's distance / influence diagnostics** to identify clusters with disproportionate leverage (e.g. very large multi-wave lineage clusters).
- **Cluster-robust standard errors** (sandwich estimator) clustered on `dz_health_board_code` as the primary approach; compare with naïve SEs to quantify the degree of spatial dependency.

---

## 8. Singleton-specific supplementary analysis

Since singletons dominate the dataset at lower resolutions and represent a distinct epidemiological category, a supplementary analysis compares the probability of being a singleton between deprived and less deprived areas using **logistic regression** with the binary outcome `is_singleton = (cluster_size == 1)`, adjusted for the same covariates. This answers the complementary question: are cases from deprived areas less likely to be genetically linked to other sequenced cases, which could reflect either fragmented transmission or surveillance gaps?

---

## 9. Sensitivity analyses

| Sensitivity | Description |
|---|---|
| **S1 — All 8 resolutions** | Repeat Models A–C at each resolution (0.1–0.8) separately; display coefficient estimates for `simd_rank_std` and `logit(dz_cum_prop_seq)` across resolutions to demonstrate robustness to the Leiden granularity choice. |
| **S2 — Index-case SIMD** | Replace mean cluster SIMD rank with the SIMD quintile of the index case (earliest collection date) as the deprivation exposure. |
| **S3 — SIMD domains** | Replace overall SIMD rank with individual SIMD domain ranks (income, employment, housing, access, health) entered separately to explore which dimension of deprivation drives any association. |
| **S4 — Exclude S0 singletons** | Restrict to clusters that passed through TN93/Leiden (i.e. exclude `|S0` clusters) to ensure findings are not driven by lineage-group singletons that were never subject to community detection. |
| **S5 — Lineage replaced by WHO VOC** | Replace pango lineage fixed effects with 6-level WHO VOC factor (Alpha, Delta, BA.1, BA.2, BA.4/5, non-VOC) to reduce dimensionality and check stability. |
| **S6 — Wave fixed effects** | Replace the time RCS with categorical wave fixed effects (derived from `method/wave_dates.parquet`) to assess sensitivity to the temporal smoothing assumption. |
| **S7 — Health Board random intercept** | Replace the cluster-robust sandwich estimator with a mixed-effects model including a random intercept for `dz_health_board_code` (14 NHS Health Boards); compare fixed-effect estimates. |
| **S8 — Maximum pairwise distance** | Use computed maximum within-cluster pairwise Euclidean distance (km) as the geographic dispersal outcome instead of `cluster_n_datazones`. |

---

## 10. Reporting

Results for each model will be reported as:

- **Incidence rate ratios (IRR)** or **odds ratios (OR)** with 95% confidence intervals (cluster-robust) for `simd_rank_std` and `logit(dz_cum_prop_seq)`, and for all other covariates.
- **Marginal effects at representative values:** predicted cluster size, duration, and dispersal at SIMD quintiles 1–5, holding all other covariates at their means, with 95% CIs from the delta method.
- **Model fit summary:** log-likelihood, AIC, BIC, and the test of the zero component for hurdle models.
- **Sensitivity S1 forest plot:** coefficient estimate and 95% CI for `simd_rank_std` across all 8 resolutions for each outcome.

---

## 11. Software

All analyses will be conducted in Python. Recommended libraries:

| Task | Library |
|---|---|
| Data manipulation and aggregation | `pandas`, `pyarrow` |
| Statistical modelling (GLM, NB, hurdle) | `statsmodels` (GLM with NB family; use `statsmodels.discrete.count_model.ZeroInflatedNegativeBinomialP` or `HurdleCountModel` when available) or `PyMC` for Bayesian variants |
| Spline basis construction | `patsy` (`cr()` or `bs()` splines) |
| Cluster-robust standard errors | `statsmodels` sandwich covariance (`HC3`, or cluster via `cov_type='cluster'`) |
| Geospatial distance computation | `numpy` (pairwise Euclidean from BNG coordinates) or `pyproj` / `geopandas` if reprojection is needed |
| Visualisation | `matplotlib`, `seaborn`, `plotnine` |

The analysis script will:
1. Load the parquet dataset using `pyarrow` and filter to `resolution == 0.3`.
2. Aggregate to cluster level following the rules in Section 2.2.
3. Compute derived variables (standardised SIMD, logit transforms, log-incidence, max pairwise distance if required).
4. Fit Models A, B, C and the singleton supplementary model.
5. Extract and format coefficient tables and marginal effects.
6. Loop over resolutions 0.1–0.8 for Sensitivity S1.

---

*Prepared: May 2026. This plan was written before any modelling was conducted; it constitutes the pre-specified analysis for Part 1.*
