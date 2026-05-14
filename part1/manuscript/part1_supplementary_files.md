# Part 1 Supplementary Files

Companion document to the Part 1 manuscript. Contains the figure captions for every supplementary figure produced by `make_figures.py`, and the supplementary tables for hurdle geographic-spread results that were dropped from Figure 4 and Supplementary Figure 7 (those are now ZTNB-only heatmaps). Figures are saved as PDF, PNG, and TIFF in this directory.

## Supplementary Figures

### Supplementary Figure 1: Cluster-outcome and excess-mixing distributions among non-singleton clusters

**File:** `supp_fig1_outcome_distributions` (PDF/PNG/TIFF)

Two-row distributional summary for the 84,067 non-singleton clusters. Top row: histograms of cluster size, cluster duration, and number of distinct datazones. Bottom row: histograms of observed-minus-expected excess mixing for age, sex, and SIMD-deprivation quintile composition. Modal peaks sit at the structural minima (size 2, duration 0 days, datazones 1) and medians at 3, 4 days, and 3 respectively. Age and sex excess mixing centre slightly above zero; SIMD-quintile excess mixing centres slightly below zero, motivating the regression mixing analyses.

### Supplementary Figure 2: Excess mixing predictor distributions

**File:** `supp_fig2_mixing_distributions` (PDF/PNG/TIFF)

Distributions of the four cluster-level mixing predictors used as explanatory variables in the mixing-predictor count models: SIMD-quintile, age, sex, and joint age-sex-profile excess-mixing scores. Each is the observed-minus-expected pair-discordance score for that attribute within non-singleton clusters, on the percentage-point scale.

### Supplementary Figure 3: Observed-minus-expected pair-probability matrices

**File:** `supp_fig3_observed_expected_matrices` (PDF/PNG/TIFF)

Heatmaps of observed-minus-expected pair probabilities for SIMD-quintile pairs (left) and age-band pairs (right), averaged across all non-singleton clusters. Cells show mean excess pair probability in percentage points relative to the lineage- and calendar-window-matched expectation. Same-quintile SIMD pairs are positive throughout, peaking at quintile 1 × quintile 1 (most deprived, +0.3 pp). Same-age-band pairs are positive on the diagonal, peaking among young adults (20-24, +0.2 pp). All cells are annotated with the numerical excess value (pp).

### Supplementary Figure 4: SIMD-domain deprivation effects on count outcomes

**File:** `supp_fig4_deprivation_domain_outcomes` (PDF/PNG/TIFF)

Coefficient plot of the per-1-SD effect of each SIMD subdomain (income, employment, education, health, access, crime, housing) on the four count-model components: cluster-size hurdle (odds ratio), positive cluster size (ZTNB count ratio), geographic-spread hurdle (odds ratio), and positive geographic spread (ZTNB count ratio). Coefficients are estimated in domain-specific models adjusting for the surveillance and lineage-window covariates. Housing and crime show the strongest negative associations with positive cluster size and positive geographic spread; access deprivation is the only domain positively associated with the geographic-spread hurdle and positive geographic spread.

### Supplementary Figure 5: SIMD-domain deprivation effects on mixing outcomes

**File:** `supp_fig5_deprivation_domain_mixing` (PDF/PNG/TIFF)

Four-panel coefficient plot of the per-1-SD effect of each SIMD subdomain on (A) domain-quintile excess mixing, (B) age excess mixing, (C) sex excess mixing, and (D) joint age-sex profile excess mixing. Mixing outcomes are observed-minus-expected pair-discordance scores in percentage points; the x-axis is the adjusted percentage-point change per 1 SD higher domain deprivation. Education and crime deprivation are associated with greater domain-quintile mixing; access and housing with the reverse pattern. Across age, sex, and joint age-sex mixing, access deprivation behaves opposite to the other six domains.

### Supplementary Figure 6: Wave-specific SIMD-domain deprivation effects on demographic mixing

**File:** `supp_fig6_deprivation_domain_wave_mixing` (PDF/PNG/TIFF)

Three stacked heatmaps showing the per-1-SD effect of each SIMD subdomain on (A) age, (B) sex, and (C) joint age-sex excess mixing for each epidemic wave. Rows are SIMD subdomains (Overall, Income, Employment, Education, Health, Access, Crime, Housing); columns are wave groups (B.1.177, Alpha, Delta, BA.1, BA.2, BA.5). Cells are coloured by adjusted percentage-point change with shared colour scale across all three panels (±5 pp), and annotated with the numeric value. Age mixing is positive for most domain × wave cells in earlier waves, with access as a consistent negative outlier; sex mixing turns negative in BA.2 and BA.5 for most domains.

### Supplementary Figure 7: Domain-specific mixing-predictor effects on ZTNB cluster outcomes

**File:** `supp_fig7_mixing_domain_outcomes` (PDF/PNG/TIFF)

Two-panel heatmap of the per-1-SD effect of four cluster-level mixing predictors (domain-quintile, age, sex, age-sex profile excess mixing) on (A) the positive cluster-size ZTNB count ratio and (B) the positive geographic-spread ZTNB count ratio, fit separately in each of eight SIMD subdomains (rows). A single ratio-scale colour bar is shared across panels and capped at ratio 5 (the largest observed cell ratio in either panel). The domain-quintile column saturates the upper triangle in both panels (cell values 2.8-3.4), and all cells are annotated with the raw ratio. The hurdle component of geographic spread is omitted from this figure and reported in Supplementary Table 2. Note: between-row variation is very small because each row represents a separate per-subdomain model in which the age/sex/age-sex mixing predictors are identical observed variables; only the domain-specific deprivation and domain-excess-mixing pair differs across rows, and Scotland's SIMD subdomains are highly correlated so per-domain coefficients land in essentially the same place.

### Supplementary Figure 8: Size-adjusted positive geographic spread

**File:** `supp_fig8_deprivation_size_adjusted` (PDF/PNG/TIFF)

Coefficient plot comparing SIMD-deprivation and surveillance-covariate effects on the positive ZTNB geographic-spread model with and without additional adjustment for log cluster size. The SIMD point estimate flips direction once cluster size is conditioned on (count ratio 0.851, 95% CI 0.792-0.915 unadjusted; 1.027, 95% CI 1.010-1.044 size-adjusted), showing that the unadjusted negative geographic-spread association is explained by deprivation's association with cluster size.

### Supplementary Figure 9: Deprivation log-linear vs hurdle/ZTNB count models

**File:** `supp_fig9_deprivation_loglinear` (PDF/PNG/TIFF)

Coefficient plot contrasting SIMD-deprivation effects on cluster size and geographic spread from a single-component log-linear (Poisson-style) model with the corresponding hurdle (odds ratio) and ZTNB (count ratio) components of the two-part main model. Log-linear estimates are substantially attenuated (cluster size geometric mean ratio 0.992, geographic spread 1.001) because they average over the structural mass at the count minimum, masking the within-component associations that the two-part model separates.

### Supplementary Figure 10: Mixing-predictor log-linear vs hurdle/ZTNB count models

**File:** `supp_fig10_mixing_loglinear` (PDF/PNG/TIFF)

Same contrast as Supplementary Figure 9 but for the four mixing predictors (SIMD-quintile, age, sex, age-sex profile excess mixing). Log-linear estimates attenuate the mixing-predictor effects on count outcomes, again because of averaging across the hurdle and positive-count components of the cluster-size / geographic-spread distributions.

## Supplementary Tables

### Supplementary Table 1: Wave-specific mixing-predictor effects on the geographic-spread hurdle (companion to Figure 4)

**File:** `supp_table_fig4_wave_mixing_hurdle_geographic_spread.csv`

Window-clustered hurdle (binomial GLM with logit link) odds ratios for the four mixing predictors (SIMD-quintile, age, sex, age-sex profile excess mixing) by epidemic wave. This component is omitted from Figure 4 because the SIMD coefficient reaches an odds ratio of ~29,000 in the Alpha wave (95% CI 6,450-129,898), making a heatmap uninformative.

> | Wave Group | Mixing predictor | Ratio (95% CI)      | p        | Notes |
> |------------|------------------|---------------------|----------|-------|
> | Alpha      | Age              | 1.41 (1.26–1.58)    | 1.5e-09  |       |
> | Alpha      | Joint profile    | 0.904 (0.833–0.981) | 0.0158   |       |
> | Alpha      | Sex              | 0.955 (0.895–1.02)  | 0.167    |       |
> | Alpha      | SIMD             | 28945 (6450–129898) | 5.2e-41  |       |
> | B.1.177    | Age              | 1.96 (1.56–2.47)    | 6.1e-09  |       |
> | B.1.177    | Joint profile    | 0.845 (0.697–1.02)  | 0.0865   |       |
> | B.1.177    | Sex              | 0.772 (0.619–0.963) | 0.0218   |       |
> | B.1.177    | SIMD             | 22.8 (13.5–38.5)    | 1.7e-31  |       |
> | BA.1       | Age              | 1.36 (1.29–1.43)    | 2.0e-32  |       |
> | BA.1       | Joint profile    | 1.04 (0.995–1.08)   | 0.088    |       |
> | BA.1       | Sex              | 0.791 (0.755–0.828) | 8.7e-24  |       |
> | BA.1       | SIMD             | 34.8 (22.2–54.6)    | 1.0e-53  |       |
> | BA.2       | Age              | 1.21 (1.14–1.3)     | 3.8e-09  |       |
> | BA.2       | Joint profile    | 1.05 (0.993–1.11)   | 0.0837   |       |
> | BA.2       | Sex              | 0.656 (0.626–0.688) | 3.4e-68  |       |
> | BA.2       | SIMD             | 27 (21.2–34.4)      | 5.5e-156 |       |
> | BA.4       | Age              | 1.67 (1.16–2.41)    | 0.00563  |       |
> | BA.4       | Joint profile    | 1.46 (1.16–1.82)    | 0.00104  |       |
> | BA.4       | Sex              | 0.478 (0.304–0.752) | 0.00142  |       |
> | BA.4       | SIMD             | 7.82 (5.34–11.4)    | 3.4e-26  |       |
> | BA.5       | Age              | 1.35 (1.19–1.54)    | 4.9e-06  |       |
> | BA.5       | Joint profile    | 1.28 (1.18–1.39)    | 2.3e-09  |       |
> | BA.5       | Sex              | 0.552 (0.479–0.637) | 3.5e-16  |       |
> | BA.5       | SIMD             | 11 (9.34–12.9)      | 3.9e-190 |       |
> | BQ.1       | Age              | 1.31 (0.996–1.72)   | 0.0536   |       |
> | BQ.1       | Joint profile    | 1.09 (0.819–1.46)   | 0.545    |       |
> | BQ.1       | Sex              | 0.641 (0.524–0.783) | 1.3e-05  |       |
> | BQ.1       | SIMD             | 7.05 (4.72–10.5)    | 1.2e-21  |       |
> | Delta      | Age              | 1.22 (1.15–1.3)     | 1.4e-10  |       |
> | Delta      | Joint profile    | 0.998 (0.973–1.02)  | 0.846    |       |
> | Delta      | Sex              | 0.841 (0.789–0.897) | 9.9e-08  |       |
> | Delta      | SIMD             | 50.5 (38.7–65.8)    | 7.1e-185 |       |

*Full coefficient, standard error, z, n_observations, and n_events columns are in the companion CSV.*

### Supplementary Table 2: SIMD-domain mixing-predictor effects on the geographic-spread hurdle (companion to Supplementary Figure 7)

**File:** `supp_table_fig7_domain_mixing_hurdle_geographic_spread.csv`

Window-clustered hurdle odds ratios for the four mixing predictors (domain-quintile, age, sex, age-sex profile excess mixing) in each per-domain model. The crime and education rows have point estimates but no cluster-robust SE, CI or p-value because the window-clustered sandwich variance estimator failed numerically for those two hurdle fits (Hessian effectively singular under the heavy outcome imbalance, 88% of clusters being multi-datazone). Point estimates remain valid.

> | Domain     | Mixing predictor | Ratio (95% CI)         | p       | Notes                                                                 |
> |------------|------------------|------------------------|---------|-----------------------------------------------------------------------|
> | Access     | Age              | 1.35 (1.31–1.4)        | 3.8e-69 |                                                                       |
> | Access     | Age-sex          | 0.995 (0.96–1.03)      | 0.789   |                                                                       |
> | Access     | Domain quintile  | 22.1 (19–25.8)         | <1e-300 |                                                                       |
> | Access     | Sex              | 0.802 (0.769–0.838)    | 9.1e-24 |                                                                       |
> | Crime      | Age              | 1.29 (CI unavailable)  | —       | cluster-robust SE unavailable (Hessian singular); point estimate only |
> | Crime      | Age-sex          | 1.04 (CI unavailable)  | —       | cluster-robust SE unavailable (Hessian singular); point estimate only |
> | Crime      | Domain quintile  | 23.1 (CI unavailable)  | —       | cluster-robust SE unavailable (Hessian singular); point estimate only |
> | Crime      | Sex              | 0.799 (CI unavailable) | —       | cluster-robust SE unavailable (Hessian singular); point estimate only |
> | Education  | Age              | 1.28 (CI unavailable)  | —       | cluster-robust SE unavailable (Hessian singular); point estimate only |
> | Education  | Age-sex          | 1.02 (CI unavailable)  | —       | cluster-robust SE unavailable (Hessian singular); point estimate only |
> | Education  | Domain quintile  | 22.1 (CI unavailable)  | —       | cluster-robust SE unavailable (Hessian singular); point estimate only |
> | Education  | Sex              | 0.794 (CI unavailable) | —       | cluster-robust SE unavailable (Hessian singular); point estimate only |
> | Employment | Age              | 1.29 (1.25–1.33)       | 2.6e-63 |                                                                       |
> | Employment | Age-sex          | 1.02 (0.99–1.05)       | 0.205   |                                                                       |
> | Employment | Domain quintile  | 21.8 (18.7–25.4)       | <1e-300 |                                                                       |
> | Employment | Sex              | 0.789 (0.759–0.821)    | 2.2e-32 |                                                                       |
> | Health     | Age              | 1.3 (1.26–1.34)        | 1.7e-62 |                                                                       |
> | Health     | Age-sex          | 1.02 (0.986–1.05)      | 0.289   |                                                                       |
> | Health     | Domain quintile  | 22.7 (19.4–26.5)       | <1e-300 |                                                                       |
> | Health     | Sex              | 0.805 (0.772–0.841)    | 4.8e-23 |                                                                       |
> | Housing    | Age              | 1.3 (1.25–1.34)        | 2.5e-54 |                                                                       |
> | Housing    | Age-sex          | 1.03 (1–1.06)          | 0.0265  |                                                                       |
> | Housing    | Domain quintile  | 23.4 (20.2–27.1)       | <1e-300 |                                                                       |
> | Housing    | Sex              | 0.778 (0.743–0.814)    | 1.6e-27 |                                                                       |
> | Income     | Age              | 1.29 (1.25–1.33)       | 6.6e-65 |                                                                       |
> | Income     | Age-sex          | 1.03 (1–1.06)          | 0.0425  |                                                                       |
> | Income     | Domain quintile  | 22.1 (19.1–25.5)       | <1e-300 |                                                                       |
> | Income     | Sex              | 0.774 (0.744–0.806)    | 3.1e-35 |                                                                       |
> | Overall    | Age              | 1.29 (1.25–1.33)       | 2.1e-61 |                                                                       |
> | Overall    | Age-sex          | 1.03 (1.01–1.06)       | 0.0151  |                                                                       |
> | Overall    | Domain quintile  | 22.2 (19.1–25.8)       | <1e-300 |                                                                       |
> | Overall    | Sex              | 0.771 (0.738–0.807)    | 5.2e-30 |                                                                       |

*Full coefficient, standard error, z, n_observations, and n_events columns are in the companion CSV.*

