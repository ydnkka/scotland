# Summary Report: Scottish SARS-CoV-2 Genomic Cluster Analyses

7 May 2026

## 1. Purpose

This report summarises the Scottish SARS-CoV-2 genomic cluster analyses intended
for three distinct thesis chapters. It outlines the focus of each chapter, the
main results, key figures and tables, and a short request for advice on how best
to shape the work into a manuscript.

## 2. Summary

The three analyses form a coherent thesis arc. Part 1 is the strongest
analytical chapter and examines deprivation, surveillance, and genomic cluster
structure. Part 2 is a descriptive vaccination chapter, showing how cluster
vaccination profiles changed through the rollout and across variant waves. Part
3 is a policy-context chapter, describing how cluster structure aligned with
Scottish restriction periods while carefully avoiding causal claims because of
strong wave and calendar-time confounding.

Across all three parts, the main insight is that genomic clusters reflect a
mixture of transmission, local social structure, variant wave, vaccination
rollout, surveillance intensity, and policy context. Deprivation matters, but
not as a simple pattern of larger or more geographically dispersed clusters in
more deprived areas. The clearest deprivation-related signals are in cluster
composition and demographic mixing.

## 3. Chapter Overview

| Thesis part | Main focus                                       | Core contribution                                                                                                                                                   | Manuscript potential                                                        |
|-------------|--------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------|
| Part 1      | Deprivation, surveillance, and cluster structure | Shows that apparent cluster size and spread are strongly shaped by surveillance and epidemic context, while deprivation signals are clearest in demographic mixing. | Strongest standalone manuscript candidate.                                  |
| Part 2      | Vaccination profiles within clusters             | Describes how clusters shifted from unvaccinated to mixed to predominantly vaccinated across rollout and variant waves.                                             | Useful as a descriptive chapter or supporting context.                      |
| Part 3      | Policy periods and cluster structure             | Places cluster dynamics on the Scottish policy timeline and uses targeted ITS analyses at selected transitions.                                                     | Best framed as contextual or secondary unless paired carefully with Part 1. |

## 4. Analysis Population

| Measure                                   |                                          Value |
|-------------------------------------------|-----------------------------------------------:|
| QC filter                                 |                       `nextclade_qc == "good"` |
| Primary clustering resolution             |                                     Leiden 0.3 |
| Unique sequences                          |                                        281,320 |
| Sequence rows used                        |                                        789,347 |
| Total inferred clusters                   |                                        193,160 |
| Cluster rows used in main Part 1/3 models |                                        193,112 |
| Non-singleton clusters                    | 84,067 in Part 1/3; 84,080 in Part 2 summaries |
| Sliding time windows                      |                                            134 |
| Raw Pango lineages                        |                                            788 |
| Modelled lineage levels after pooling     |                                            183 |

## 5. Part 1: Socioeconomic Deprivation, Surveillance, and Cluster Structure

### 5.1. Focus

Part 1 asks whether clusters linked to more deprived areas were larger, more
geographically dispersed, or more socially and demographically mixed. Outcomes
include cluster size, number of distinct datazones, and observed-minus-expected
mixing by SIMD, age, sex, and joint SIMD-age-sex profile.

### 5.2. Key Results

| Result                                          |                                       Estimate | Interpretation                                                                                                     |
|-------------------------------------------------|-----------------------------------------------:|--------------------------------------------------------------------------------------------------------------------|
| SIMD deprivation and non-singleton cluster odds |                OR 0.971, 95% CI 0.960 to 0.983 | More deprived clusters were slightly less likely to exceed singleton size.                                         |
| SIMD deprivation and positive cluster size      |       Count ratio 0.926, 95% CI 0.869 to 0.987 | Among non-singleton clusters, higher deprivation was associated with smaller clusters.                             |
| SIMD deprivation and multi-datazone spread      |                OR 1.004, 95% CI 0.992 to 1.016 | No clear association with crossing the one-datazone threshold.                                                     |
| SIMD deprivation and positive geographic spread |       Count ratio 0.851, 95% CI 0.792 to 0.915 | Higher deprivation was associated with lower geographic spread among clusters already spanning multiple datazones. |
| SIMD-quintile excess mixing                     |  +0.31 percentage points, 95% CI -0.18 to 0.80 | No clear evidence of increased SIMD-quintile mixing.                                                               |
| Age excess mixing                               |   +1.66 percentage points, 95% CI 1.29 to 2.03 | More deprived clusters showed greater age mixing.                                                                  |
| Sex excess mixing                               | -0.78 percentage points, 95% CI -1.16 to -0.39 | More deprived clusters showed lower sex mixing.                                                                    |
| Joint SIMD-age-sex excess mixing                |   +0.48 percentage points, 95% CI 0.29 to 0.67 | More deprived clusters showed greater joint demographic-profile mixing.                                            |

### 5.3. Main Insight

The main models do not support a simple hypothesis that deprivation produces
larger or more geographically dispersed genomic clusters. Surveillance and
epidemic context, including local incidence, sequencing intensity, window-level
sequencing proportion, and test positivity, were more consistent predictors of
apparent cluster scale. Deprivation-related signals were clearest in demographic
mixing and varied across SIMD domains.

### 5.4. Key Figures

| Figure          | What it shows                                                                                    |
|-----------------|--------------------------------------------------------------------------------------------------|
| Part 1 Figure 1 | Main hurdle and zero-truncated negative-binomial results for cluster size and geographic spread. |
| Part 1 Figure 2 | Adjusted associations with SIMD, age, sex, and joint demographic excess mixing.                  |
| Part 1 Figure 3 | How different SIMD domains relate to deprivation-domain and demographic mixing.                  |
| Part 1 Figure 4 | Wave-specific deprivation associations for cluster size and geographic spread.                   |

![Part 1 Figure 1: Main cluster outcomes](part1/main/manuscript/figures/fig1_main_cluster_outcomes.png)

![Part 1 Figure 2: Main cluster mixing](part1/main/manuscript/figures/fig2_main_cluster_mixing.png)

![Part 1 Figure 3: SIMD-domain mixing](part1/main/manuscript/figures/fig3_simd_domain_mixing.png)

![Part 1 Figure 4: Wave-specific cluster outcomes](part1/main/manuscript/figures/fig4_wave_specific_cluster_outcomes.png)

### 5.5. Key Tables

| Table                                  | File                                                                                                                             |
|----------------------------------------|----------------------------------------------------------------------------------------------------------------------------------|
| Dataset descriptives                   | [main_dataset_descriptives.csv](part1/main/tables/main_dataset_descriptives.csv)                                                 |
| Main hurdle count model results        | [main_hurdle_count_model_results.csv](part1/main/tables/main_hurdle_count_model_results.csv)                                     |
| Main mixing model results              | [main_mixing_model_results.csv](part1/main/tables/main_mixing_model_results.csv)                                                 |
| SIMD-domain demographic mixing results | [main_simd_domain_demographic_mixing_model_results.csv](part1/main/tables/main_simd_domain_demographic_mixing_model_results.csv) |
| Wave-specific hurdle count results     | [main_wave_specific_hurdle_count_model_results.csv](part1/main/tables/main_wave_specific_hurdle_count_model_results.csv)         |

## 6. Part 2: Vaccination Profiles Within Genomic Clusters

### 6.1. Focus

Part 2 characterises vaccination status within non-singleton genomic clusters
across epidemic waves, cluster size, geographic-dispersion categories, and SIMD
quintiles. This is descriptive rather than causal and should not be interpreted
as vaccine-effectiveness analysis.

### 6.2. Key Results

| Wave    | None vaccinated | Mixed vaccination | All vaccinated |
|---------|----------------:|------------------:|---------------:|
| B.1.177 |           89.5% |             10.0% |           0.5% |
| Alpha   |           41.5% |             56.2% |           2.3% |
| Delta   |           10.3% |             65.3% |          24.4% |
| BA.1    |            4.4% |             49.9% |          45.7% |
| BA.2    |            2.6% |             41.0% |          56.4% |
| BA.5    |            1.8% |             27.1% |          71.0% |
| BQ.1    |            0.6% |             20.5% |          78.9% |
| XBB     |            0.0% |             19.4% |          80.6% |

| Vaccination-status mixing result        | Pattern                                                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------------------------------|
| Mixed vaccination-status clusters       | Peaked in Delta, when 49.1% of non-singleton clusters were more mixed than expected.                                    |
| Homogeneous vaccination-status clusters | Became dominant from BA.1 onwards, reaching 66.3% in BA.5.                                                              |
| Booster coverage by SIMD                | Modest positive Q5-Q1 gradients in most post-Delta waves, indicating higher booster coverage in less deprived clusters. |
| Dose recency                            | Median days since vaccination increased from 13 days in B.1.177 to 218 days in XBB.                                     |
| Large clusters                          | More likely to be homogeneous by SIMD, consistent with geographic concentration of large transmission events.           |

### 6.3. Main Insight

Cluster vaccination profiles largely tracked the national rollout. Mixed
vaccination-status clusters were most common when vaccine coverage was changing
rapidly across the population. Later clusters were predominantly all-vaccinated
because vaccination had become the background condition among sequenced cases,
not because vaccinated cases formed exceptional clusters.

### 6.4. Key Figures

| Figure          | What it shows                                                                                     |
|-----------------|---------------------------------------------------------------------------------------------------|
| Part 2 Figure 1 | Weekly proportion of sequenced cases who were vaccinated, by rollout age group and SIMD quintile. |
| Part 2 Figure 2 | Shift from unvaccinated to mixed and all-vaccinated clusters across waves and size categories.    |
| Part 2 Figure 3 | Vaccination-status mixing categories by wave, relative to same-lineage/window expectations.       |
| Part 2 Figure 6 | Booster coverage and days since last vaccination by wave and SIMD quintile.                       |

![Part 2 Figure 1: Vaccinated cases over time](part2/manuscript/figures/fig1_vaccinated_cases_over_time.png)

![Part 2 Figure 2: Cluster vaccination by wave and category](part2/manuscript/figures/fig2_cluster_vaccination_by_wave_and_category.png)

![Part 2 Figure 3: Vaccination-status mixing by wave](part2/manuscript/figures/fig3_vaccination_mixing_by_wave.png)

![Part 2 Figure 6: Booster coverage and dose recency by SIMD](part2/manuscript/figures/fig6_dose_recency_by_simd.png)

### 6.5. Key Tables

| Table                                     | File                                                                                                        |
|-------------------------------------------|-------------------------------------------------------------------------------------------------------------|
| Vaccination descriptives                  | [vaccination_descriptives.csv](part2/tables/vaccination_descriptives.csv)                                   |
| Vaccination cluster table                 | [vaccination_cluster_table.csv](part2/tables/vaccination_cluster_table.csv)                                 |
| Cluster category summary                  | [cluster_category_summary.csv](part2/tables/cluster_category_summary.csv)                                   |
| Vaccination cluster wave/category summary | [vaccination_cluster_wave_category_summary.csv](part2/tables/vaccination_cluster_wave_category_summary.csv) |
| SIMD-domain vaccination gradients         | [supp_simd_domain_vaccination_gradients.csv](part2/tables/supp_simd_domain_vaccination_gradients.csv)       |

## 7. Part 3: Policy Restrictions and Cluster Structure

### 7.1. Focus

Part 3 examines how Scottish COVID-19 restriction periods align with genomic
cluster size, geographic dispersion, and demographic mixing. Because restriction
periods are strongly confounded with calendar time, variant waves, immunity, and
surveillance changes, the analysis is descriptive and associational.

### 7.2. Policy Period Codes

Policy intensity is an ordinal restriction score where higher values indicate
more restrictive conditions. The genomic cluster analysis begins in P3, so E0,
L1, P1, and P2 are included below for code completeness but are not observed in
the Part 3 cluster models.

| Code | Policy period label  | Dates                    | Intensity | In cluster analysis |
|------|----------------------|--------------------------|----------:|---------------------|
| E0   | Emergence            | 2020-03-01 to 2020-03-23 |        15 | No                  |
| L1   | First lockdown       | 2020-03-24 to 2020-05-28 |       100 | No                  |
| P1   | Route map phase 1    | 2020-05-29 to 2020-06-18 |        72 | No                  |
| P2   | Route map phase 2    | 2020-06-19 to 2020-07-09 |        52 | No                  |
| P3   | Route map phase 3    | 2020-07-10 to 2020-10-01 |        30 | Yes                 |
| T1   | Pre-tier tightening  | 2020-10-02 to 2020-11-01 |        55 | Yes                 |
| F5   | Five-tier framework  | 2020-11-02 to 2021-01-04 |        65 | Yes                 |
| L2   | Second lockdown      | 2021-01-05 to 2021-04-01 |        95 | Yes                 |
| SL   | Stay local - Level 3 | 2021-04-02 to 2021-04-25 |        65 | Yes                 |
| L3   | Level 3              | 2021-04-26 to 2021-05-16 |        55 | Yes                 |
| L21  | Level 2 / Level 1    | 2021-05-17 to 2021-07-18 |        38 | Yes                 |
| L0   | Level 0              | 2021-07-19 to 2021-08-08 |        20 | Yes                 |
| NN   | Near-normal          | 2021-08-09 to 2021-11-28 |        10 | Yes                 |
| OM   | Omicron wave         | 2021-11-29 to 2022-01-23 |        42 | Yes                 |
| FE   | Final easing         | 2022-01-24 to 2022-04-17 |        15 | Yes                 |
| PR   | Post-restriction     | 2022-04-18 to 2023-05-05 |         3 | Yes                 |

### 7.3. Key Results

| Analysis                                            |              Result | Interpretation                                                                   |
|-----------------------------------------------------|--------------------:|----------------------------------------------------------------------------------|
| Policy intensity and weekly median log cluster size | Spearman rho = 0.74 | Strong positive correlation, likely driven by wave confounding.                  |
| Policy intensity and weekly log datazones           | Spearman rho = 0.58 | Moderate positive correlation, also confounded by wave and surveillance context. |
| Policy intensity and mean SIMD excess discordance   | Spearman rho = 0.02 | No meaningful correlation with SIMD mixing.                                      |
| Policy intensity and mean age excess discordance    | Spearman rho = 0.59 | Positive correlation, but not interpretable causally.                            |

| ITS transition         |            Cluster size |       Geographic spread | Mixing metrics         | Interpretation                                                                                           |
|------------------------|------------------------:|------------------------:|------------------------|----------------------------------------------------------------------------------------------------------|
| T1 onset, October 2020 |  beta = -0.13, p = 0.22 |  beta = -0.07, p = 0.59 | No significant changes | No detectable structural change at the tier introduction.                                                |
| L2 to SL, April 2021   | beta = -0.21, p = 0.034 | beta = -0.36, p = 0.006 | No significant changes | Decline likely reflects Alpha-wave tail dynamics rather than easing itself.                              |
| NN onset, August 2021  | beta = +0.12, p = 0.068 | beta = +0.32, p = 0.015 | No significant changes | Geographic spread increased after legal distancing requirements ended, within strong Delta-wave context. |

### 7.4. Main Insight

Restriction intensity was correlated with cluster size and spread across the
whole epidemic, but these correlations are mainly contextual rather than causal.
The strongest transition-level signal was increased geographic spread after the
August 2021 move to near-normal restrictions. Demographic mixing did not change
significantly across the selected policy transitions, supporting the Part 1
interpretation that cluster composition reflects local social structure more
than broad policy intensity.

### 7.5. Key Figures

| Figure          | What it shows                                                                          |
|-----------------|----------------------------------------------------------------------------------------|
| Part 3 Figure 1 | Weekly cluster outcomes shown against policy-period shading and restriction intensity. |
| Part 3 Figure 2 | Interrupted time-series plots for the T1 onset, L2 to SL, and NN onset transitions.    |
| Part 3 Figure 3 | Median cluster size and geographic spread compared across observed policy periods.     |

![Part 3 Figure 1: Weekly time series and policy context](part3/manuscript/figures/fig1_weekly_time_series.png)

![Part 3 Figure 2: Interrupted time-series transitions](part3/manuscript/figures/fig2_its_transitions.png)

![Part 3 Figure 3: Policy-period outcomes](part3/manuscript/figures/fig3_period_outcomes.png)

### 7.6. Key Tables

| Table                         | File                                                                                          |
|-------------------------------|-----------------------------------------------------------------------------------------------|
| Policy-period descriptives    | [period_descriptives.csv](part3/tables/period_descriptives.csv)                               |
| Intensity correlations        | [intensity_correlations.csv](part3/tables/intensity_correlations.csv)                         |
| ITS coefficients              | [its_coefficients.csv](part3/tables/its_coefficients.csv)                                     |
| Lagged intensity correlations | [supp_lagged_intensity_correlations.csv](part3/tables/supp_lagged_intensity_correlations.csv) |
| ITS window sensitivity        | [supp_its_window_sensitivity.csv](part3/tables/supp_its_window_sensitivity.csv)               |

## 8. Overall Interpretation

The strongest recurring conclusion is that genomic clusters cannot be interpreted
as direct measures of transmission scale without accounting for surveillance and
epidemic context. Deprivation, vaccination, and policy context all matter, but
they operate through different parts of the system:

- Deprivation is most visible in cluster composition and demographic mixing.
- Vaccination status mainly reflects rollout timing and background population
  coverage among sequenced cases.
- Policy periods help contextualise the epidemic timeline, but policy effects
  cannot be separated cleanly from variant waves and calendar time in these
  data.

## 9. Manuscript Advice Request

I would appreciate advice on the best manuscript strategy. My current view is
that Part 1 is the strongest standalone paper because it has the clearest
inferential structure and central argument: socioeconomic deprivation does not
map simply onto larger or more geographically dispersed genomic clusters once
surveillance and epidemic context are accounted for, but it is associated with
demographic mixing and cluster composition.

One option is therefore to develop Part 1 as the main manuscript and use Parts 2
and 3 as thesis chapters or as contextual supporting material. Another option is
to draw selectively across all three parts to write a broader paper on the
social, vaccination, and policy context of SARS-CoV-2 genomic clusters in
Scotland. My concern is that combining all three parts may make the manuscript
too broad, while focusing only on Part 1 may underuse the richer vaccination and
policy analyses. I would be grateful for advice on which framing is most
publishable and whether Parts 2 and 3 should be kept separate, used as
supplementary context, or developed into future papers.