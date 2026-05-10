# Summary Report: Scottish SARS-CoV-2 Genomic Cluster Analyses

10 May 2026

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
Scottish restriction periods, characterising the emergence and growth of the Alpha
variant, and presenting a meta-cluster analysis of the major pre-lockdown Alpha
lineages.

Across all three parts, the main insight is that genomic clusters reflect a
mixture of transmission, local social structure, variant wave, vaccination
rollout, surveillance intensity, and policy context. Deprivation matters, but
not as a simple pattern of larger or more geographically dispersed clusters in
more deprived areas. The clearest deprivation-related signals are in cluster
composition and demographic mixing.

## 3. Chapter Overview

| Thesis part | Main focus | Core contribution | Manuscript potential |
|-------------|------------|-------------------|----------------------|
| Part 1 | Deprivation, surveillance, and cluster structure | Shows that apparent cluster size and spread are strongly shaped by surveillance and epidemic context, while deprivation signals are clearest in demographic mixing. | Strongest standalone manuscript candidate. |
| Part 2 | Vaccination profiles within clusters | Describes how clusters shifted from unvaccinated to mixed to predominantly vaccinated across rollout and variant waves. | Useful as a descriptive chapter or supporting context. |
| Part 3 | Policy periods, Alpha emergence, and cluster structure | Places cluster dynamics on the Scottish policy timeline, uses targeted ITS analyses at selected transitions, characterises Alpha emergence and growth advantage under different restriction levels, and maps pre-lockdown Alpha meta-cluster structure. | Best framed as contextual or secondary unless paired carefully with Part 1; Alpha emergence case study adds independent interest. |

## 4. Analysis Population

| Measure | Value |
|---------|------:|
| QC filter | `nextclade_qc == "good"` |
| Primary clustering resolution | Leiden 0.3 |
| Unique sequences | 281,320 |
| Sequence rows used | 789,347 |
| Total inferred clusters | 193,160 |
| Cluster rows used in main Part 1/3 models | 193,112 |
| Non-singleton clusters | 84,067 in Part 1/3; 84,080 in Part 2 summaries |
| Sliding time windows | 134 |
| Raw Pango lineages | 788 |
| Modelled lineage levels after pooling | 183 |

## 5. Part 1: Socioeconomic Deprivation, Surveillance, and Cluster Structure

### 5.1. Focus

Part 1 asks whether clusters linked to more deprived areas were larger, more
geographically dispersed, or more socially and demographically mixed. Outcomes
include cluster size, number of distinct datazones, and observed-minus-expected
mixing by SIMD, age, sex, and joint SIMD-age-sex profile.

### 5.2. Key Results

| Result | Estimate | Interpretation |
|--------|----------|----------------|
| SIMD deprivation and non-singleton cluster odds | OR 0.971, 95% CI 0.960 to 0.983 | More deprived clusters were slightly less likely to exceed singleton size. |
| SIMD deprivation and positive cluster size | Count ratio 0.926, 95% CI 0.869 to 0.987 | Among non-singleton clusters, higher deprivation was associated with smaller clusters. |
| SIMD deprivation and multi-datazone spread | OR 1.004, 95% CI 0.992 to 1.016 | No clear association with crossing the one-datazone threshold. |
| SIMD deprivation and positive geographic spread | Count ratio 0.851, 95% CI 0.792 to 0.915 | Higher deprivation was associated with lower geographic spread among clusters already spanning multiple datazones. |
| SIMD-quintile excess mixing | +0.31 percentage points, 95% CI -0.18 to 0.80 | No clear evidence of increased SIMD-quintile mixing. |
| Age excess mixing | +1.66 percentage points, 95% CI 1.29 to 2.03 | More deprived clusters showed greater age mixing. |
| Sex excess mixing | -0.78 percentage points, 95% CI -1.16 to -0.39 | More deprived clusters showed lower sex mixing. |
| Joint SIMD-age-sex excess mixing | +0.48 percentage points, 95% CI 0.29 to 0.67 | More deprived clusters showed greater joint demographic-profile mixing. |

### 5.3. Main Insight

The main models do not support a simple hypothesis that deprivation produces
larger or more geographically dispersed genomic clusters. Surveillance and
epidemic context, including local incidence, sequencing intensity, window-level
sequencing proportion, and test positivity, were more consistent predictors of
apparent cluster scale. Deprivation-related signals were clearest in demographic
mixing and varied across SIMD domains.

### 5.4. Key Figures

| Figure | What it shows |
|--------|---------------|
| Part 1 Figure 1 | Main hurdle and zero-truncated negative-binomial results for cluster size and geographic spread. |
| Part 1 Figure 2 | Adjusted associations with SIMD, age, sex, and joint demographic excess mixing. |
| Part 1 Figure 3 | How different SIMD domains relate to deprivation-domain and demographic mixing. |
| Part 1 Figure 4 | Wave-specific deprivation associations for cluster size and geographic spread. |

![Part 1 Figure 1: Main cluster outcomes](part1/main/manuscript/figures/fig1_main_cluster_outcomes.png)

![Part 1 Figure 2: Main cluster mixing](part1/main/manuscript/figures/fig2_main_cluster_mixing.png)

![Part 1 Figure 3: SIMD-domain mixing](part1/main/manuscript/figures/fig3_simd_domain_mixing.png)

![Part 1 Figure 4: Wave-specific cluster outcomes](part1/main/manuscript/figures/fig4_wave_specific_cluster_outcomes.png)

### 5.5. Key Tables

| Table | File |
|-------|------|
| Dataset descriptives | [main_dataset_descriptives.csv](part1/main/tables/main_dataset_descriptives.csv) |
| Main hurdle count model results | [main_hurdle_count_model_results.csv](part1/main/tables/main_hurdle_count_model_results.csv) |
| Main mixing model results | [main_mixing_model_results.csv](part1/main/tables/main_mixing_model_results.csv) |
| SIMD-domain demographic mixing results | [main_simd_domain_demographic_mixing_model_results.csv](part1/main/tables/main_simd_domain_demographic_mixing_model_results.csv) |
| Wave-specific hurdle count results | [main_wave_specific_hurdle_count_model_results.csv](part1/main/tables/main_wave_specific_hurdle_count_model_results.csv) |

## 6. Part 2: Vaccination Profiles Within Genomic Clusters

### 6.1. Focus

Part 2 characterises vaccination status within non-singleton genomic clusters
across epidemic waves, cluster size, geographic-dispersion categories, and SIMD
quintiles. This is descriptive rather than causal and should not be interpreted
as vaccine-effectiveness analysis.

### 6.2. Key Results

| Wave | None vaccinated | Mixed vaccination | All vaccinated |
|------|----------------:|------------------:|---------------:|
| B.1.177 | 89.5% | 10.0% | 0.5% |
| Alpha | 41.5% | 56.2% | 2.3% |
| Delta | 10.3% | 65.3% | 24.4% |
| BA.1 | 4.4% | 49.9% | 45.7% |
| BA.2 | 2.6% | 41.0% | 56.4% |
| BA.5 | 1.8% | 27.1% | 71.0% |
| BQ.1 | 0.6% | 20.5% | 78.9% |
| XBB | 0.0% | 19.4% | 80.6% |

| Vaccination-status mixing result | Pattern |
|----------------------------------|---------|
| Mixed vaccination-status clusters | Peaked in Delta, when 49.1% of non-singleton clusters were more mixed than expected. |
| Homogeneous vaccination-status clusters | Became dominant from BA.1 onwards, reaching 66.3% in BA.5. |
| Booster coverage by SIMD | Modest positive Q5-Q1 gradients in most post-Delta waves, indicating higher booster coverage in less deprived clusters. |
| Dose recency | Median days since vaccination increased from 13 days in B.1.177 to 218 days in XBB. |
| Large clusters | More likely to be homogeneous by SIMD, consistent with geographic concentration of large transmission events. |

### 6.3. Main Insight

Cluster vaccination profiles largely tracked the national rollout. Mixed
vaccination-status clusters were most common when vaccine coverage was changing
rapidly across the population. Later clusters were predominantly all-vaccinated
because vaccination had become the background condition among sequenced cases,
not because vaccinated cases formed exceptional clusters.

### 6.4. Key Figures

| Figure | What it shows |
|--------|---------------|
| Part 2 Figure 1 | Weekly proportion of sequenced cases who were vaccinated, by rollout age group and SIMD quintile. |
| Part 2 Figure 2 | Shift from unvaccinated to mixed and all-vaccinated clusters across waves and size categories. |
| Part 2 Figure 3 | Vaccination-status mixing categories by wave, relative to same-lineage/window expectations. |
| Part 2 Figure 6 | Booster coverage and days since last vaccination by wave and SIMD quintile. |

![Part 2 Figure 1: Vaccinated cases over time](part2/manuscript/figures/fig1_vaccinated_cases_over_time.png)

![Part 2 Figure 2: Cluster vaccination by wave and category](part2/manuscript/figures/fig2_cluster_vaccination_by_wave_and_category.png)

![Part 2 Figure 3: Vaccination-status mixing by wave](part2/manuscript/figures/fig3_vaccination_mixing_by_wave.png)

![Part 2 Figure 6: Booster coverage and dose recency by SIMD](part2/manuscript/figures/fig6_dose_recency_by_simd.png)

### 6.5. Key Tables

| Table | File |
|-------|------|
| Vaccination descriptives | [vaccination_descriptives.csv](part2/tables/vaccination_descriptives.csv) |
| Vaccination cluster table | [vaccination_cluster_table.csv](part2/tables/vaccination_cluster_table.csv) |
| Cluster category summary | [cluster_category_summary.csv](part2/tables/cluster_category_summary.csv) |
| Vaccination cluster wave/category summary | [vaccination_cluster_wave_category_summary.csv](part2/tables/vaccination_cluster_wave_category_summary.csv) |
| SIMD-domain vaccination gradients | [supp_simd_domain_vaccination_gradients.csv](part2/tables/supp_simd_domain_vaccination_gradients.csv) |

## 7. Part 3: Policy Restrictions, Alpha Emergence, and Cluster Structure

### 7.1. Focus

Part 3 examines how Scottish COVID-19 restriction periods align with genomic
cluster size, geographic dispersion, and demographic mixing; characterises the
emergence and comparative growth of Alpha (B.1.1.7) under different restriction
levels; and maps the major pre-lockdown Alpha lineages using a rolling-window
meta-cluster network. Because restriction periods are strongly confounded with
calendar time, variant waves, immunity, and surveillance changes, all
associations between policy and genomic outcomes are descriptive and
associational.

Two per-period summary metrics complement the standard log-median outcomes
throughout the chapter:

- **Clustering rate** = (total sequences − total clusters) / total sequences: the
  fraction of sequenced isolates that are secondary cluster members rather than
  isolated singletons. Higher values indicate more transmission-chain structure.
- **Dispersion parameter k̂** (method of moments) = x̄² / (s² − x̄), where x̄
  is the period mean cluster size and s² is the variance. Values below 1 indicate
  overdispersion (a few very large clusters dominate); k̂ → ∞ approaches Poisson
  spread. This parameter is reported on a log scale in figures because of the
  wide cross-period range.

### 7.2. Policy Period Codes

Policy intensity is an ordinal restriction score where higher values indicate
more restrictive conditions. The genomic cluster analysis begins in P3, so E0,
L1, P1, and P2 are included below for code completeness but are not observed in
the Part 3 cluster models.

| Code | Policy period label | Dates | Intensity | In cluster analysis |
|------|---------------------|-------|----------:|---------------------|
| E0 | Emergence | 2020-03-01 to 2020-03-23 | 15 | No |
| L1 | First lockdown | 2020-03-24 to 2020-05-28 | 100 | No |
| P1 | Route map phase 1 | 2020-05-29 to 2020-06-18 | 72 | No |
| P2 | Route map phase 2 | 2020-06-19 to 2020-07-09 | 52 | No |
| P3 | Route map phase 3 | 2020-07-10 to 2020-10-01 | 30 | Yes |
| T1 | Pre-tier tightening | 2020-10-02 to 2020-11-01 | 55 | Yes |
| F5 | Five-tier framework | 2020-11-02 to 2021-01-04 | 65 | Yes |
| L2 | Second lockdown | 2021-01-05 to 2021-04-01 | 95 | Yes |
| SL | Stay local — Level 3 | 2021-04-02 to 2021-04-25 | 65 | Yes |
| L3 | Level 3 | 2021-04-26 to 2021-05-16 | 55 | Yes |
| L21 | Level 2 / Level 1 | 2021-05-17 to 2021-07-18 | 38 | Yes |
| L0 | Level 0 | 2021-07-19 to 2021-08-08 | 20 | Yes |
| NN | Near-normal | 2021-08-09 to 2021-11-28 | 10 | Yes |
| OM | Omicron wave | 2021-11-29 to 2022-01-23 | 42 | Yes |
| FE | Final easing | 2022-01-24 to 2022-04-17 | 15 | Yes |
| PR | Post-restriction | 2022-04-18 to 2023-05-05 | 3 | Yes |

### 7.3. Descriptive Results: Policy Intensity Correlations

Cross-period Spearman correlations show that restriction intensity was
positively correlated with median log cluster size (ρ = 0.741) and negatively
correlated with singleton fraction (ρ = −0.621), but these associations are
mainly driven by the co-occurrence of high-restriction periods with large variant
waves. The association with mean SIMD excess discordance was negligible (ρ =
0.019), consistent with Part 1's finding that within-cluster socioeconomic mixing
is not strongly linked to the epidemic trajectory.

Per-period clustering rate and k̂ values for the selected phases are shown below.
k̂ was below 1.0 in all periods, indicating persistent overdispersion throughout
the epidemic. The Alpha-wave periods (F5, SL) showed the relatively highest k̂
values (0.26–0.36) among selected phases, while the long multi-variant periods
(OM k̂ = 0.010, L21 k̂ = 0.017) and the post-restriction phase (PR k̂ = 0.082)
showed the most extreme overdispersion.

| Period | Clustering rate | k̂ (all clusters) | k̂ (non-singletons) |
|--------|----------------:|------------------:|---------------------:|
| P3 | 0.852 | 0.108 | 0.206 |
| T1 | 0.780 | 0.172 | 0.359 |
| F5 | 0.675 | 0.265 | 0.472 |
| L2 | 0.831 | 0.144 | 0.247 |
| SL | 0.844 | 0.362 | 0.615 |
| L0 | 0.634 | 0.216 | 0.360 |
| NN | 0.736 | 0.143 | 0.247 |

### 7.4. Interrupted Time-Series Results

Three policy transitions were analysed using segmented OLS regression at a
primary ±8-week window. The level-change estimate β_post is the primary
parameter of interest.

#### Primary outcomes: log-median cluster size and geographic spread

| ITS transition | Cluster size β_post | Datazones β_post | Interpretation |
|----------------|--------------------:|-----------------:|----------------|
| T1 onset, October 2020 | −0.080, p = 0.197 | −0.082, p = 0.588 | No detectable structural change at the tier introduction. |
| L2 to SL, April 2021 | −0.171, p = 0.102 | −0.355, p < 0.001 | Geographic decline likely reflects Alpha-wave tail dynamics rather than easing itself; cluster size change no longer significant. |
| NN onset, August 2021 | +0.116, p = 0.016 | +0.312, p = 0.036 | Both outcomes increased after legal distancing requirements ended, within a stable Delta-wave context — the most policy-consistent signal in the chapter. |

#### Supplementary outcomes: clustering rate and dispersion k̂

The same ITS framework was applied to clustering rate and k̂ as distributional
complements. Level-change estimates at each transition were broadly consistent
with the log-median results:

| ITS transition | Clustering rate β_post (cluster) | Clustering rate β_post (geo) | k̂ β_post (geo, log scale) |
|----------------|----------------------------------:|-----------------------------:|---------------------------:|
| T1 onset | −0.017 | −0.015 | −0.290 |
| L2 to SL | −0.095 | −0.115 | +0.233 |
| NN onset | +0.191 | +0.210 | −0.149 |

The positive k̂ at L2→SL reflects a transient increase in geographic
concentration (fewer, relatively larger geo-spread clusters) as the Alpha wave
wound down. The negative k̂ at NN-onset, alongside the positive clustering rate,
indicates that Delta-wave clusters increased overall chain membership while also
becoming more geographically dispersed (more, moderately sized geo-spread
clusters rather than a single dominant lineage).

### 7.5. Alpha Emergence Case Study

Binomial GLM growth models fitted to S:N501Y (Alpha) and S:A222V (B.1.177)
marker frequencies show that Alpha expanded faster under the less restrictive
five-tier framework (F5) than under the subsequent second lockdown (L2):

| Phase | Growth rate r (per day) | Doubling time | OR per week |
|-------|------------------------:|-------------:|------------:|
| F5 (five-tier framework) | 0.0851 | 8.1 days | 1.815 |
| L2 (second lockdown) | 0.0618 | 11.2 days | 1.542 |

Counterfactual projections estimate the effect of imposing L2-level restrictions
earlier:

| Earlier L2 imposition | Estimated delay to 50% Alpha dominance |
|-----------------------|----------------------------------------:|
| 4 days earlier | 11 days |
| 34 days earlier | 13 days |
| 64 days earlier | 25 days |

These projections suggest that earlier restrictions would have modestly slowed
Alpha's ascent but could not have prevented establishment, given the growth
advantage already present.

### 7.6. Pre-L2 Alpha Meta-Cluster Analysis

The pre-L2 Alpha meta-cluster network was constructed by linking rolling-window
cluster assignments that share at least one sequence, then extracting connected
components. Key findings:

| Measure | Value |
|---------|------:|
| Connected components (meta-clusters) | 78 |
| Unique pre-L2 Alpha sequences covered | 442 |
| Largest component (AM001) — unique sequences | 234 (52.9% of pre-L2 Alpha total) |
| AM001 primary health board | Greater Glasgow and Clyde |
| ORF1a:L730F prevalence in AM001 | 85% |
| ORF1a:L730F prevalence in non-AM001 Alpha | 18.3% |

The strong enrichment of ORF1a:L730F in AM001 identifies this mutation as a
lineage marker for the dominant pre-L2 chain. The concentration of over half of
all pre-L2 Alpha sequences in a single meta-cluster implies that most of the
pre-lockdown Alpha burden was driven by a small number of high-amplification
components that had established across multiple regions before L2 could take
effect.

### 7.7. Main Insight

Restriction intensity was correlated with cluster size and spread across the
whole epidemic, but these correlations are primarily contextual rather than
causal because restriction periods co-occur with distinct variant waves and
surveillance regimes. The strongest transition-level signal was a coincident
increase in both cluster size and geographic spread after the August 2021 move
to near-normal restrictions, occurring within the stable Delta-dominant context
where wave confounding was least severe. Demographic mixing did not change
significantly across the selected policy transitions, supporting the Part 1
interpretation that cluster composition reflects local social structure more
than broad policy intensity.

The Alpha emergence analysis adds independent evidence that restriction level
modulated growth rate: Alpha expanded roughly 40% faster (by doubling time)
under the five-tier framework than under the second lockdown. However, the
meta-cluster results show that Alpha had already established a large, geographically
dispersed dominant chain (AM001) before L2 was imposed, limiting the potential
impact of even earlier intervention.

### 7.8. Key Figures

| Figure | What it shows |
|--------|---------------|
| Part 3 Figure 1 | Four-panel weekly time series: policy colour strip, median cluster size with IQR shading, clustering rate, and dispersion k̂ (log scale). |
| Part 3 Figure 2 | ITS analyses (primary) at T1 onset, L2→SL, and NN onset — log-median cluster size (left) and log datazones (right), with IQR error bars. |
| Part 3 Figure 3 | Alpha emergence: weekly S:N501Y and S:A222V frequencies, and health-board-by-phase heatmap of pre-L2 Alpha sequences. |
| Part 3 Figure 4 | Counterfactual projections for earlier L2 imposition on Alpha dominance timing. |
| Part 3 Supplementary Figure 2a | ITS clustering rate at the three transitions (cluster size left, geographic right). |
| Part 3 Supplementary Figure 2b | ITS dispersion k̂ at the three transitions (cluster size left, geographic right, log scale). |
| Part 3 Supplementary Figure (meta-cluster) | Pre-L2 Alpha meta-cluster top-component demographics and ORF1a:L730F enrichment. |

![Part 3 Figure 1: Weekly time series and policy context](part3/manuscript/figures/fig1_policy_timeline_cluster_structure.png)

![Part 3 Figure 2: Interrupted time-series transitions](part3/manuscript/figures/fig2_selected_policy_transitions.png)

![Part 3 Figure 3: Alpha emergence](part3/manuscript/figures/fig3_alpha_emergence_f5_l2.png)

![Part 3 Figure 4: Counterfactual projections](part3/manuscript/figures/fig4_alpha_counterfactual_timing.png)

![Part 3 Supplementary Figure 2a: ITS clustering rate](part3/manuscript/figures/supp_fig2a_its_clustering_rate.png)

![Part 3 Supplementary Figure 2b: ITS dispersion](part3/manuscript/figures/supp_fig2b_its_dispersion.png)

### 7.9. Key Tables

| Table | File |
|-------|------|
| Policy-period descriptives | [period_descriptives.csv](part3/tables/period_descriptives.csv) |
| Per-period clustering rate and k̂ | [period_clustering_dispersion.csv](part3/tables/period_clustering_dispersion.csv) |
| Intensity correlations | [intensity_correlations.csv](part3/tables/intensity_correlations.csv) |
| ITS coefficients | [its_coefficients.csv](part3/tables/its_coefficients.csv) |
| Lagged intensity correlations | [supp_lagged_intensity_correlations.csv](part3/tables/supp_lagged_intensity_correlations.csv) |
| ITS window sensitivity | [supp_its_window_sensitivity.csv](part3/tables/supp_its_window_sensitivity.csv) |

## 8. Overall Interpretation

The strongest recurring conclusion is that genomic clusters cannot be interpreted
as direct measures of transmission scale without accounting for surveillance and
epidemic context. Deprivation, vaccination, and policy context all matter, but
they operate through different parts of the system:

- Deprivation is most visible in cluster composition and demographic mixing,
  not in raw cluster size or geographic spread.
- Vaccination status mainly reflects rollout timing and background population
  coverage among sequenced cases.
- Policy periods help contextualise the epidemic timeline, but policy effects
  on cluster structure cannot be separated cleanly from variant waves and calendar
  time — with the partial exception of the NN-onset transition and the comparative
  Alpha growth rates under F5 and L2.

The Alpha emergence case study illustrates how the analysis moves beyond purely
descriptive correlation: the growth-model comparison provides a quantified
estimate of restriction-level effects on variant spread under explicitly stated
counterfactual assumptions, while the meta-cluster structure shows how established
transmission chains constrain the realistic impact of any single policy moment.

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
policy analyses.

For Part 3 specifically, the Alpha emergence case study — including the
counterfactual projections and the meta-cluster mutation enrichment — may have
enough independent interest to support a focused short paper, distinct from the
broader policy-context framing. I would be grateful for advice on which framing
is most publishable and whether Parts 2 and 3 should be kept separate, used as
supplementary context, or developed into future papers.
