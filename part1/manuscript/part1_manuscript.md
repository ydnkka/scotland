# Socioeconomic deprivation, local surveillance intensity, and the structure of SARS-CoV-2 genomic transmission clusters in Scotland, 2020–2023

Dominic Arthur$^1$, Christopher J. Banks$^1$, Rowland R. Kao$^2$\*

$^1$ The Roslin Institute, University of Edinburgh, Edinburgh, United Kingdom

$^2$ School of Physics and Astronomy, University of Edinburgh, Edinburgh, United Kingdom

\* Correspondence: rowland.kao@ed.ac.uk

---

## 1. Introduction

The COVID-19 pandemic accelerated the integration of pathogen genomic surveillance into routine public-health practice, generating unprecedented volumes of viral sequence data linked to demographic and epidemiological metadata. The COVID-19 Genomics UK (COG-UK) Consortium, in particular, established one of the world’s largest integrated genomic-surveillance systems, enabling near real-time reconstruction of viral introductions, lineage replacement, and outbreak structure across the United Kingdom [1,2]. More broadly, pathogen genomics is now recognised as a core component of infectious-disease surveillance, complementing conventional case-, admission-, and mortality-based epidemiology by providing information on genetic relatedness and transmission dynamics [3--5].

Within Scotland, genomic epidemiology contributed substantially to the national COVID-19 response. Early analyses demonstrated that the first epidemic wave was driven by multiple introductions, predominantly from continental Europe [6]. Subsequent studies investigated transmission in settings of particular epidemiological or policy relevance, including renal dialysis units [7], university accommodation [8], hospital-onset infections [9,10], and care homes [11]. Variant-focused analyses further showed that the Alpha variant was associated with increased clinical severity, while the Delta wave displayed marked age- and deprivation-patterned differences in hospitalisation risk [12,13]. Genomic surveillance data have also been used to evaluate travel policies and quantify importation pressure into Scotland [14]. Collectively, these studies established the value of genomic data for understanding transmission pathways in Scotland; however, they have largely focused on early introductions, specific institutional outbreaks, or discrete epidemic phases rather than on nationwide recent-transmission structure across the full epidemic.

A parallel literature in both the UK and internationally has demonstrated that the burden of COVID-19 was strongly socially patterned. Higher infection rates, more severe illness, and poorer outcomes were consistently observed in more deprived populations and among groups with greater occupational exposure, overcrowded housing, and higher reliance on essential work [15--18]. Scotland followed this broader pattern, with deprivation associated with inequalities in infection, hospital admission, and mortality [12,13,19,20,62]. Despite this, most inequalities research has relied on case counts, testing data, hospital admissions, or mortality records rather than on genome-informed representations of recent transmission [16,17,20]. Consequently, it remains unclear whether deprivation is associated solely with a higher burden of diagnosed disease, or whether it is also reflected in transmission structure itself---for example, through participation in larger recent-transmission clusters or in clusters characterised by particular sociodemographic mixing patterns.

Linking SARS-CoV-2 genomic data to recent-transmission structure at population scale presents substantial methodological challenges. Detailed phylodynamic and transmission-inference approaches, including outbreaker2 [21], TransPhylo [22], and SCOTTI [23], are well suited to reconstructing modest-scale outbreaks but become computationally demanding at national scale and often require assumptions that are difficult to justify in routine surveillance settings. By contrast, threshold- or distance-based clustering approaches, such as fixed single-nucleotide polymorphism (SNP) cut-offs or hierarchical genetic-distance clustering, are scalable and computationally efficient but have limited biological grounding. Small genetic distances do not map cleanly onto recent transmission, particularly when pathogen diversity is low, sampling is incomplete, or sampling time contains information ignored by static genetic thresholds [24--28]. These limitations are amplified in superspreading-driven epidemics, where many epidemiologically distinct infections may remain genetically similar because transmission occurs rapidly over short time intervals [29,30].

To generate interpretable nationwide cluster definitions from a large genomic-surveillance dataset, we applied EpiLink, a process-based pairwise compatibility model recently developed by our group [31]. EpiLink evaluates whether the observed temporal and genetic separation between sampled cases is compatible with simulated distributions generated under plausible recent-transmission scenarios, parameterised using published estimates of SARS-CoV-2 substitution rate, infectiousness profile, and disease progression dynamics [32–35]. Unlike static SNP thresholds, EpiLink produces graded pairwise compatibility scores that explicitly propagate uncertainty in incubation timing, testing delay, and mutation accumulation into downstream clustering. Combined with the Leiden community-detection algorithm [36], this framework provides an interpretable, threshold-free approach for reconstructing recent-transmission clusters across the full surveillance period.

A second, complementary question arises once both cluster composition and cluster scale are quantified. If deprivation partly shapes transmission through differences in contact structure across age groups, sex, and socioeconomic strata, then within-cluster sociodemographic mixing is not merely an outcome of deprivation but may itself influence how transmission propagates through the population. Social-contact studies consistently demonstrate strong age-assortative mixing, alongside additional structuring by setting and household composition [37,38]. Theoretical and empirical work on superspreading further suggests that transmission events bridging otherwise distinct social groups can seed disproportionately large outbreaks [30,31]. From this perspective, clusters whose sampled cases bridge demographic and socioeconomic strata more than expected are likely to appear larger and more geographically dispersed. Treating sociodemographic mixing as predictors of cluster size and geographic spread therefore addresses a distinct question: whether the bridging properties of clusters are associated with their reconstructed epidemiological scale.

In this study, we applied the EpiLink–Leiden framework to the Scottish national SARS-CoV-2 genomic-surveillance dataset to investigate two complementary sets of questions. First, treating deprivation as the exposure, we examined whether mean cluster-level area deprivation, measured using the Scottish Index of Multiple Deprivation (SIMD) [39], was associated with cluster size and geographic spread; whether deprivation was associated with within-cluster sociodemographic mixing, quantified as observed-minus-expected pairwise discordance across SIMD quintiles, age bands, sex, and joint sociodemographic profiles; and whether these associations differed by SIMD domain and epidemic wave. Second, treating excess mixing as the exposure of interest, we examined whether within-cluster bridging across sociodemographic strata predicted cluster size and geographic spread after adjustment for the same covariates, and whether these relationships varied across SIMD domains and epidemic waves.

Throughout, we treated local sequencing intensity and test positivity not only as potential confounders but also as substantive surveillance variables, because the genomic structure recovered by surveillance systems depends fundamentally on which infections are sampled, where, and when [40,41]. Our central premise is that the structure of inferred SARS-CoV-2 transmission clusters is jointly shaped by social, surveillance, and epidemic contexts. Consequently, interpretation of cluster-based inequalities in transmission should be considered context-dependent rather than monotonic or uniform across the epidemic.

---

## 2. Methods

### 2.1. Study population and source data

We analysed SARS-CoV-2 whole-genome sequences collected from Scottish residents between July 2020 and February 2023 through Public Health Scotland and the COVID-19 Genomics UK (COG-UK) surveillance infrastructure [1,6]. Consensus genomes were aligned to the SARS-CoV-2 reference genome (GenBank accession MN908947) and quality-controlled using Nextclade [42]; only records passing the Nextclade `good` overall quality filter were retained. Sequence records were linked to collection date, anonymised patient identifier, residential datazone, age band, sex, Pango lineage, and the Scottish Index of Multiple Deprivation 2020 version 2 (SIMD). SIMD linkage included the overall rank and seven domain ranks: income, employment, education, health, geographic access, crime, and housing [39]. Datazone-level (standard small-area geography of Scotland with $\approx$760 residents on average) testing aggregates and health-board-level surveillance trends were used to derive local surveillance covariates.

### 2.2. Cluster inference pipeline

Cluster inference followed the EpiLink-Leiden framework described elsewhere [31]. Briefly, sequences were grouped into 3-week sliding windows advanced in 1-week steps and then partitioned by Pango lineage. Within each window-lineage group containing at least two sequences, pairwise consensus genetic distances were computed under the TN93 substitution model [43] using the standalone `tn93` tool [44] and paired with specimen-collection time differences. Each pair was scored using EpiLink, a process-based compatibility model that evaluates whether the observed genetic and temporal separation is consistent with simulated recent-transmission scenarios under specified assumptions about mutation accumulation, disease progression, and infectiousness [31]. Pairwise compatibility scores were used as graph edge weights, and recent-transmission clusters were inferred using the Leiden community-detection algorithm with a modularity-based objective [36] at the primary resolution $\gamma = 0.3$ (partition are stability across $\gamma \in$ [0.1, 0.6] in the EpiLink benchmark). Sequences in window-lineage groups with fewer than two sequences were assigned to singleton clusters.

### 2.3. Cluster-level outcomes and covariates

The resulting long-format cluster table, with one row for each sequence, window, and Leiden resolution combination, was collapsed to one row per cluster at the primary resolution. For cluster $c$, the primary count outcomes were cluster size $N_c$ (number of unique sequences) and geographic dispersion $D_c$ (number of distinct residential datazones represented in the cluster). Cluster duration was calculated as the interval between the earliest and latest sampled case but retained only descriptively because the fixed 3-week window structure mechanically constrains observed duration.

Because both primary count outcomes are bounded below by 1, we modelled excess counts above the structural minimum:

$$
Y^{(N)}_c = N_c - 1,
\qquad
Y^{(D)}_c = D_c - 1.
$$

The associated hurdle indicators were:

$$
H^{(N)}_c=\mathbb{I}(N_c > 1),
\qquad
H^{(D)}_c=\mathbb{I}(D_c > 1).
$$

Mean cluster-level deprivation was defined as the mean SIMD rank across sampled cases in the cluster, multiplied by $-1$ so that higher values indicate greater deprivation. Local cumulative incidence was log-transformed after scaling per 1,000 population. Local cumulative sequencing fraction, window-level sequencing proportion, and local 7-day test positivity were logit-transformed after clipping probabilities away from 0 and 1. Continuous covariates were standardised to mean 0 and unit variance. Calendar time was represented using a cubic B-spline basis over the window index, and Pango lineages represented by fewer than 50 clusters were pooled into an "Other rare lineages" category.

### 2.4. Within-cluster sociodemographic mixing

Within-cluster sociodemographic mixing was measured using observed-minus-expected pairwise discordance. For categorical variable $V$, let $N_{ck}$ be the number of sampled cases in cluster $c$ belonging to category $k$, therefore $N_c=\sum_k N_{ck}$. Observed pairwise discordance was defined as:

$$
d^{\mathrm{obs}}_{c,V} = 1 - \sum_k\frac{N_{ck}(N_{ck}-1)}{N_c(N_c-1)}.
$$

Expected discordance was computed from all sampled cases in the same window-lineage stratum. If $m_{wlk}$ is the number of cases in window $w$, lineage $l$, and category $k$, with $m_{wl}=\sum_k m_{wlk}$, then:

$$
d^{\mathrm{exp}}_{wl,V} = 1 - \sum_k\frac{m_{wlk}(m_{wlk}-1)}{m_{wl}(m_{wl}-1)}.
$$

The cluster-level excess-mixing outcome was:

$$
e_{c,V} = d^{\mathrm{obs}}_{c,V} - d^{\mathrm{exp}}_{w(c)l(c),V}.
$$

Positive values indicate more cross-category mixing than expected for the cluster's window-lineage context, whereas negative values indicate more assortative composition. Excess discordance was computed for SIMD quintile, age band, sex, and a joint SIMD-age-sex profile.

### 2.5. Regression Analyses

The primary regression analysis examined deprivation as the exposure. For each count outcome, we used a two-part hurdle formulation to separate the probability of exceeding the structural minimum from the magnitude of the positive excess count [50,51; see Technical Implementation Appendix, Sections 9-10]. The first component modelled whether cluster $c$ exceeded the structural minimum:

$$
H_c \sim \operatorname{Bernoulli}(\pi_c),
\qquad
\operatorname{logit}(\pi_c) = \mathbf{x}_c^\top\gamma.
$$

The second component was fitted only among clusters with positive excess counts and used a zero-truncated negative-binomial model:

$$
Y_c \mid Y_c > 0 \sim \operatorname{ZTNB}(\mu_c,\alpha),
\qquad
\log(\mu_c) = \mathbf{x}_c^\top\beta.
$$

Here, $Y_c$ denotes the count above the structural minimum. For cluster size, this was additional sequences beyond a singleton cluster; for geographic dispersion, it was additional datazones beyond a single-datazone cluster. The covariate vector $\mathbf{x}_c$ included standardised deprivation, local cumulative incidence, local cumulative sequencing fraction, window-level sequencing proportion, local test positivity, calendar spline terms, and lineage fixed effects. For geographic dispersion, an additional positive-count specification adjusted for log cluster size to distinguish wider spatial spread from larger reconstructed clusters.

The ZTNB component used an NB2 negative-binomial parent distribution, with variance function [52,53]:

$$
\operatorname{Var}(Y_c) = \mu_c + \alpha\mu_c^2,
$$

where $\alpha>0$ was estimated on the log scale. If $f_{\mathrm{NB}}(y;\mu,\alpha)$ is the ordinary negative-binomial probability mass function and $p_0=f_{\mathrm{NB}}(0;\mu,\alpha)$, then the zero-truncated probability for $y=1,2,\ldots$ is:

$$
P(Y=y \mid Y>0) = \frac{f_{\mathrm{NB}}(y;\mu,\alpha)}{1-p_0}.
$$

Full details of the parameterisation, likelihood, analytical score, optimisation, and numerical checks are provided in the Technical Implementation Appendix.

For the excess-mixing outcomes, we fitted linear models among non-singleton clusters:

$$
e_{c,V} = \mathbf{x}_c^\top \delta + \varepsilon_c,
$$

where $\mathbf{x}_c$ included the same deprivation and surveillance covariates, calendar spline terms, lineage fixed effects, and log cluster size. Coefficients were reported as percentage-point differences in excess discordance per 1 standard deviation higher covariate value.

To address the complementary question of whether bridging across sociodemographic strata predicted cluster scale, we refitted the count models after adding four standardised excess-mixing predictors: SIMD quintile, age band, sex, and joint SIMD-age-sex profile. These models retained deprivation, surveillance covariates, calendar splines, and lineage fixed effects as adjustments. Because excess mixing is undefined for singleton clusters, the cluster-size hurdle component was not estimable in this analysis. The positive cluster-size ZTNB component, the geographic-dispersion hurdle and positive ZTNB components, and the size-adjusted positive geographic-dispersion component were fitted among non-singleton clusters.

### 2.6. Inference

All models used cluster-robust standard errors, clustered by analysis window in the primary analysis to account for dependence induced by overlapping windows [54,55]. A sensitivity analysis instead clustered standard errors by health board.

For the custom ZTNB fits, standard errors were obtained using a sandwich covariance estimator with an observed-information bread and score-based meat [54,55; see Technical Implementation Appendix, Section 11]. Let $\hat{\theta}$ denote the fitted parameter vector and $s_c(\hat{\theta})$ the observation-level score. For robust cluster $g$, the cluster-level score was:

$$
S_g(\hat{\theta}) = \sum_{c\in g}s_c(\hat{\theta}).
$$

The sandwich covariance estimator was:

$$
\widehat{\operatorname{Var}}(\hat{\theta}) =
\widehat{B}^{-1}
\left[\sum_g S_g(\hat{\theta})S_g(\hat{\theta})^\top\right]
\widehat{B}^{-1},
$$

where $\widehat{B}$ is the observed information matrix. A finite-sample correction,

$$
\frac{G}{G-1}\frac{n-1}{n-k},
$$

was applied when $G>1$ and $n>k$, where $G$ is the number of robust clusters, $n$ is the number of observations, and $k$ is the total number of estimated parameters, including the ZTNB dispersion parameter.

Hurdle coefficients were reported as adjusted odds ratios. ZTNB coefficients were exponentiated and reported as adjusted ZTNB ratios for the parent negative-binomial mean parameter; because zero truncation also changes the conditional positive mean through $1-p_0$, exact conditional mean contrasts require prediction from the fitted model. Linear mixing coefficients were reported as adjusted percentage-point differences.

### 2.7. Sensitivity and extension analyses

We conducted sensitivity analyses using health-board-clustered standard errors, a log offset for the number of sequences in the analysis window in the positive cluster-size model [NOTE TO SELF: `PRIMARY_TERMS` includes `window_seq_fraction_z`], size adjustment for positive geographic dispersion, index-case SIMD deprivation, 99th-percentile winsorisation of positive counts, an approximately non-overlapping window subset retaining every third window, and log-linear single-component models as comparators. We also fitted wave-stratified models for B.1.177, Alpha, Delta, BA.1, BA.2, BA.4, BA.5, and BQ.1, excluding XBB from regression models because of small sample size. Finally, domain-specific models were fitted for all seven SIMD domains and extended to the mixing-predictor count-model framework.

---

## 3. Results

### 3.1. Analysis dataset

After applying Nextclade QC and the primary Leiden resolution ($\gamma$ = 0.3), the analysis dataset comprised 789,347 sequence-rows representing 281,320 unique sequences from Scottish residents, distributed across 134 sliding three-week windows and 788 raw Pango lineages (183 levels after rare-lineage pooling). The cluster table contained 193,112 inferred genomic clusters, of which 84,067 (43.5%) were non-singletons and were retained for mixing analyses. Outcome distributions were strongly right-skewed: median cluster size was 1 with 56.5% singletons and a maximum of 2,792 sequences; median distinct datazones was 1 with 61.7% confined to a single datazone and a maximum of 2,100 datazones (Supplementary Figure 1). Among non-singleton clusters, median size was 3 and median distinct datazones was 3. Wave-level summaries showed a marked shift towards smaller, more spatially confined clusters in late Omicron subwaves (Table 1).

> **Table 1.** Analysis population and outcome distributions across the Scottish SARS-CoV-2 genomic surveillance dataset, July 2020 to February 2023, at primary Leiden resolution $\gamma$ = 0.3.
> 
>**Panel A. Data description**
> 
> | Quantity                                                                  |                  Value |
> |---------------------------------------------------------------------------|-----------------------:|
> | *Cohort*                                                                  |                        |
> | Sequence-rows used                                                        |                789,347 |
> | Unique sequences                                                          |                281,320 |
> | Sliding 3-week windows                                                    |                    134 |
> | Raw Pango lineages                                                        |                    788 |
> | Pango lineage levels modelled (after pooling lineages with < 50 clusters) |                    183 |
> | Inferred clusters ($\gamma$ = 0.3)                                        |                193,112 |
> | Non-singleton clusters used in mixing models                              |                 84,067 |
> | *Outcome distributions across all clusters*                               |                        |
> | Cluster size: structural minimum (singleton)                              |                  56.5% |
> | Cluster size: median / 75th / 90th / 99th percentile / maximum            | 1 / 3 / 6 / 39 / 2,792 |
> | Duration: structural minimum (0 days)                                     |                  63.1% |
> | Duration: median / 75th / 90th / 99th percentile / maximum (days)         |    0 / 3 / 7 / 12 / 19 |
> | Distinct datazones: structural minimum (1 datazone)                       |                  61.7% |
> | Distinct datazones: median / 75th / 90th / 99th percentile / maximum      | 1 / 2 / 5 / 32 / 2,100 |
> 
> **Panel B. Per-wave breakdown**  
> Dominant Pango-lineage waves used in wave-stratified analyses (additional clusters from inter-wave periods are not shown). "Mean cluster size (non-singletons)" is the mean of cluster sizes among clusters with size > 1, computed as 1 plus the mean positive-count value reported in the analysis tables. "Single-datazone %" is the fraction of non-singleton clusters confined to a single datazone, excluding singletons (which are single-datazone by construction).
>
> | Wave    | Clusters | Singleton % |                  Mean cluster size |                  Single-datazone % |
> |---------|---------:|------------:|-----------------------------------:|-----------------------------------:|
> | B.1.177 |    4,621 |        54.6 |                               7.95 |                                8.3 |
> | Alpha   |   12,112 |        43.0 |                              11.04 |                               13.4 |
> | Delta   |   74,272 |        54.3 |                               7.62 |                               11.1 |
> | BA.1    |   32,928 |        57.3 |                              10.41 |                               10.5 |
> | BA.2    |   38,893 |        56.3 |                               7.81 |                               13.1 |
> | BA.4    |    2,669 |        66.5 |                               5.30 |                               12.8 |
> | BA.5    |   16,423 |        66.9 |                               4.77 |                               13.4 |
> | BQ.1    |    3,314 |        70.9 |                               3.29 |                               21.2 |
> | XBB     |      509 |        74.7 |                               3.86 |                               22.5 |

### 3.2. Deprivation was not associated with larger or more geographically dispersed clusters

In the primary pooled hurdle/ZTNB models, there was little evidence that clusters with higher mean SIMD deprivation were larger or more geographically dispersed (Figure 1 A-D). Higher cluster-level deprivation was associated with slightly lower odds of being non-singleton (OR = 0.971, 95% CI 0.960–0.983, p < 0.001) and with smaller positive cluster sizes among non-singleton clusters (count ratio = 0.926, 95% CI 0.869–0.987, p = 0.018) [50,52]. For geographic spread, the hurdle component was close to the null (OR = 1.004, 95% CI 0.992–1.016, p = 0.522), whereas the positive-count component was significantly below 1.0 (count ratio = 0.851, 95% CI 0.792–0.915, p < 0.001). Thus, among clusters that spanned more than one data zone, clusters with higher mean deprivation tended to span fewer additional data zones.

![Figure 1](figures/fig1_deprivation_overall.png)

However, this negative association appeared to be explained by cluster size. When the positive geographic-spread model was additionally adjusted for log cluster size, the deprivation estimate reversed direction (count ratio = 1.027, 95% CI 1.019–1.035, p < 0.001; Supplementary Figure 8). This suggests that the unadjusted negative association between deprivation and geographic spread was driven primarily by the tendency for more deprived clusters to be smaller, rather than by evidence of reduced spatial diffusion among clusters of comparable size.

Surveillance and epidemic-intensity covariates showed substantially stronger associations with cluster scale than deprivation. Per 1 SD increase in each covariate, positive cluster size was most strongly associated with local sequencing fraction (count ratio = 3.24), test positivity (2.65), local incidence (1.65), and window-level sequencing proportion (1.31). A similar pattern was observed for positive geographic spread, with count ratios of 2.27 for local sequencing fraction, 3.00 for test positivity, 1.70 for local incidence, and 1.27 for window-level sequencing proportion. These effects are consistent with the expectation that observed genomic cluster size and spatial extent are shaped not only by transmission dynamics, but also by epidemic intensity and sequencing coverage [56, 57].

After adjustment for log cluster size, only local incidence retained a positive association with geographic spread, although the effect was substantially attenuated (count ratio = 1.06). In contrast, local sequencing fraction (0.94), test positivity (0.90), and window-level sequencing proportion (0.94) were negatively associated with size-adjusted geographic spread. This indicates that these covariates were primarily associated with the number of sampled cases in a cluster, rather than with greater geographic dispersion conditional on cluster size.

Overall, these results do not support the hypothesis that higher deprivation was associated with larger or more spatially dispersed SARS-CoV-2 genomic clusters. Instead, deprivation showed weak negative associations with cluster size and unadjusted geographic spread, while surveillance intensity and local epidemic burden were much stronger predictors of observed cluster scale. Because SIMD is an area-based measure of relative deprivation across Scottish data zones, these findings should be interpreted as associations with the deprivation profile of places, rather than as direct individual-level socioeconomic effects [39].

### 3.3. Deprivation was associated with cluster composition

Mixing analyses revealed a different pattern. In the primary models among 84,067 non-singleton clusters (Figure 1 E-H), mean cluster SIMD deprivation was not clearly associated with SIMD-quintile excess mixing (+0.31 percentage points, 95% CI -0.18 to +0.80, p = 0.208) but was clearly associated with greater age-band excess mixing (+1.66 pp, 95% CI +1.29 to +2.03, p < 0.001), with reduced sex excess mixing (−0.78 pp, 95% CI −1.16 to −0.39, p < 0.001), and with slightly greater joint SIMD-age-sex profile excess mixing (+0.48 pp, 95% CI +0.29 to +0.67, p < 0.001). Cluster size was the strongest cross-outcome predictor of mixing (per 1 SD log size: +7.49 pp for SIMD, +2.80 pp for age, −1.20 pp for sex, +0.98 pp for joint profile), and surveillance covariates generally reduced mixing point estimates. 

Observed-minus-expected pair-probability matrices showed non-random sociodemographic composition within genomic clusters (Supplementary Figure 3) . For SIMD quintile, same-quintile excess mixing was strongest in the most deprived quintile, with quintile 1 x 1 pairs observed 0.28 pp more often than expected. Quintile 2 x 2 and 3 x 3 pairs were also modestly above expectation, whereas pairings involving the least deprived quintile were generally below expectation, including quintile 5 x 5 (-0.09 pp). For age, excess pairing was concentrated among young adults: 20-24 x 20-24 showed the largest positive excess (+0.21 pp), with adjacent 15-19 x 20-24 and 20-24 x 25-29 pairings also above expectation. These patterns indicate that reconstructed clusters retain measurable socioeconomic and age structure beyond that expected from the background lineage-window composition [37,38], but they should be interpreted as sampled genomic-cluster composition rather than direct contact or transmission matrices.

Full numerical estimates for the main count and mixing models are given in Table 2.

> **Table 2.** Adjusted estimates from the main pooled hurdle / zero-truncated negative-binomial (ZTNB) count models and the main excess-mixing linear models. Cluster size and geographic spread (number of distinct datazones) are modelled in two parts: a binomial hurdle component (whether the cluster exceeds its structural minimum, reported as odds ratios) and a ZTNB positive-count component among clusters above the minimum (reported as count ratios). Excess-mixing models are linear regressions among 84,067 non-singleton clusters, with coefficients reported as adjusted percentage-point differences in observed-minus-expected pairwise discordance. All models adjust for Pango lineage (183 levels), calendar time (cubic B-spline), and the covariates shown. Mixing models additionally adjust for log cluster size. Standard errors are clustered by analysis window. Estimates per 1 SD higher covariate. 95% CIs and p-values are shown for the primary exposure (SIMD deprivation); for surveillance covariates, point estimates are shown and the corresponding figures (Figs 1–2) and supplementary tables include full inference.
> 
> ---
>
> **Panel A. Count models**  
> n = 193,112 clusters for hurdle components; n = 84,067 / 74,010 for positive cluster size and positive geographic spread respectively.
>
> | Covariate (per 1 SD)         |              Cluster size hurdle (OR) |          Cluster size positive (CR) | Geographic spread hurdle (OR) |       Geographic spread positive (CR) |
> |------------------------------|--------------------------------------:|------------------------------------:|------------------------------:|--------------------------------------:|
> | **SIMD deprivation**         | **0.971 (0.960–0.983)<sup>***</sup>** | **0.926 (0.869–0.987)<sup>*</sup>** |       **1.004 (0.992–1.016)** | **0.851 (0.792–0.915)<sup>***</sup>** |
> | Local cumulative incidence   |                                 1.173 |                               1.650 |                         1.223 |                                 1.699 |
> | Local sequencing fraction    |                                 1.067 |                               3.240 |                         1.047 |                                 2.269 |
> | Window sequencing proportion |                                 1.252 |                               1.314 |                         1.170 |                                 1.274 |
> | Local test positivity        |                                 1.448 |                               2.649 |                         1.314 |                                 2.999 |
>
> **Panel B. Mixing models**  
> n = 84,067 non-singleton clusters; outcomes in percentage-point excess discordance.
>
> | Covariate (per 1 SD)         | SIMD-quintile mixing (pp) |                   Age-band mixing (pp) |                        Sex mixing (pp) | Joint SIMD-age-sex profile mixing (pp) |
> |------------------------------|--------------------------:|---------------------------------------:|---------------------------------------:|---------------------------------------:|
> | **SIMD deprivation**         |  **+0.31 (−0.18, +0.80)** | **+1.66 (+1.29, +2.03)<sup>***</sup>** | **−0.78 (−1.16, −0.39)<sup>***</sup>** | **+0.48 (+0.29, +0.67)<sup>***</sup>** |
> | Local cumulative incidence   |                     +4.86 |                                  +1.56 |                                  −0.79 |                                  +1.09 |
> | Local sequencing fraction    |                     −1.22 |                                  −1.02 |                                  −0.36 |                                  −0.54 |
> | Window sequencing proportion |                     −3.79 |                                  −0.76 |                                  +0.61 |                                  −0.30 |
> | Local test positivity        |                     −6.65 |                                  −0.73 |                                  +1.15 |                                  −0.98 |
> | Log cluster size             |                     +7.49 |                                  +2.80 |                                  −1.20 |                                  +0.98 |
>
> **Panel C. Mixing-predictor count models**  
> n = 84,067 non-singleton clusters for cluster size and the geographic-spread hurdle; n = 74,010 for positive geographic spread. The four excess-mixing scores enter the same hurdle/ZTNB specification as Panel A, alongside SIMD deprivation, the four surveillance covariates, lineage fixed effects, and the calendar B-spline. Estimates per 1 SD higher mixing score.
>
> | Mixing predictor (per 1 SD)       |         Cluster size positive (CR) |         Geographic spread hurdle (OR) |    Geographic spread positive (CR) |
> |-----------------------------------|-----------------------------------:|--------------------------------------:|-----------------------------------:|
> | SIMD-quintile excess discordance  | **3.48 (3.25–3.73)<sup>***</sup>** | **22.11 (18.98–25.75)<sup>***</sup>** | **3.03 (2.79–3.29)<sup>***</sup>** |
> | Age-band excess discordance       | **1.67 (1.56–1.78)<sup>***</sup>** |    **1.28 (1.25–1.32)<sup>***</sup>** | **1.97 (1.80–2.15)<sup>***</sup>** |
> | Sex excess discordance            |  **0.85 (0.76–0.95)<sup>**</sup>** |    **0.77 (0.74–0.81)<sup>***</sup>** |                   1.06 (0.92–1.22) |
> | Joint SIMD-age-sex profile excess | **0.81 (0.76–0.87)<sup>***</sup>** |     **1.03 (1.01–1.05)<sup>**</sup>** | **0.72 (0.67–0.78)<sup>***</sup>** |
> ---
> OR = adjusted odds ratio; CR = adjusted count ratio; pp = percentage points. Significance markers: <sup>*</sup> p < 0.05; <sup>**</sup> p < 0.01; <sup>***</sup> p < 0.001. SIMD ranks were negated and standardised, so positive coefficients in mixing models, and ratios < 1 in count models, indicate associations in the direction of greater deprivation. Excess-mixing scores in Panel C are observed-minus-expected pairwise discordance within lineage × window strata, standardised to unit variance. Sensitivity analyses are summarised in the Discussion and reported in full in the supplementary tables.

### 3.4. SIMD domains behaved heterogeneously

Domain-specific models indicated that the overall SIMD association concealed marked heterogeneity across deprivation domains, consistent with the construction of SIMD as a multidimensional small-area index rather than a single-axis measure of socioeconomic position [39]. For domain-quintile excess mixing (Supplementary Figure 5), education deprivation (+1.17 pp, 95% CI +0.64 to +1.69) and crime deprivation (+1.09 pp, 95% CI +0.55 to +1.64) were associated with greater within-cluster mixing, whereas geographic access deprivation (−2.08 pp, 95% CI −2.50 to −1.66) and housing deprivation (−1.19 pp, 95% CI −1.77 to −0.62) were associated with lower mixing, giving an approximate 3 pp range across the seven SIMD domains. For age mixing, most domains showed positive associations with deprivation, with geographic access deprivation again emerging as a consistent negative outlier. Sex mixing showed the opposite pattern, with most domains negatively associated with deprivation. For joint demographic-profile mixing, associations were broadly positive and closely resembled those observed for age mixing. Domain-specific estimates for cluster size and geographic spread (Supplementary Figure 4) showed less heterogeneity: with the exception of geographic access deprivation, which was weakly associated with smaller positive cluster size and lower positive geographic spread, the remaining six domains mirrored the overall negative SIMD association reported above. These findings suggest that SIMD domains capture distinct social and spatial dimensions of deprivation relevant to transmission processes and surveillance ascertainment, rather than acting as interchangeable indicators of a single deprivation construct.

### 3.5. Deprivation effects were not stable across epidemic waves

Per-wave analyses (Figure 2) showed that the pooled negative association between deprivation and cluster size or geographic spread was not stable across the epidemic. The strongest and most consistent negative deprivation associations occurred during the Delta wave: cluster size hurdle OR 0.934 (95% CI 0.921–0.947), positive cluster size count ratio 0.797 (95% CI 0.725–0.876), geographic spread hurdle OR 0.958 (95% CI 0.945–0.971), and positive geographic spread count ratio 0.781 (95% CI 0.703–0.867). Omicron sublineage BA.2 showed a contrasting pattern, with positive cluster size count ratio 1.176 (95% CI 1.061–1.302) and weakly higher positive geographic spread. Omicron sublineage BA.4 also showed positive associations, but with very wide intervals reflecting the small per-wave sample (n = 2,669 clusters). Domain-wave heatmaps for demographic mixing (Supplementary Figure 6) showed broadly stable directional patterns for age and joint-profile mixing but with greater wave-to-wave variability.

![Figure 2](figures/fig2_deprivation_wave_specific.png)

### 3.6. Within-cluster mixing strongly predicted cluster scale

Reversing the modelling direction, we asked whether the four excess-mixing scores themselves predicted cluster size and geographic spread (Figure 3; Table 2 Panel C). Because mixing scores are undefined for singletons, the cluster-size hurdle is not estimable in these models; we therefore report three count components: positive cluster size (ZTNB), geographic-spread hurdle (multi- vs single-datazone among non-singletons), and positive geographic spread (ZTNB). Among non-singleton clusters, SIMD-quintile excess discordance was by far the strongest mixing-side predictor of all three components. Per 1 SD higher SIMD excess discordance, positive cluster size was 3.48-fold higher (95% CI 3.25–3.73), positive geographic spread was 3.03-fold higher (95% CI 2.79–3.29), and the odds of exceeding a single datazone were 22.1-fold higher (95% CI 18.98–25.75). Age excess discordance was also substantially positively associated with both count outcomes (positive cluster size CR 1.67, 95% CI 1.56–1.78; positive geographic spread CR 1.97, 95% CI 1.80–2.15), and weakly positively associated with the geographic-spread hurdle (OR 1.28, 95% CI 1.25–1.32). Sex excess discordance showed a small negative association with positive cluster size (CR 0.85, 95% CI 0.76–0.95) and with the geographic-spread hurdle (OR 0.77, 95% CI 0.74–0.81), and was near null for positive geographic spread. Joint SIMD-age-sex profile excess discordance was negatively associated with both positive cluster size (CR 0.81, 95% CI 0.76–0.87) and positive geographic spread (CR 0.72, 95% CI 0.67–0.78). Taken together, these estimates show that clusters bridging more across SIMD quintiles and age bands than the lineage-window baseline expectation are detected as substantially larger and more geographically dispersed, while clusters whose excess discordance is concentrated in joint profile coordinates — that is, in fine-grained sociodemographic recombinations beyond what marginal SIMD and age mixing already explain — are detected as smaller and less dispersed.

![Figure 3](figures/fig3_mixing_overall.png)

### 3.7. SIMD-domain mixing predictors showed heterogeneity

Substituting each of the seven SIMD-domain quintile mixing scores for the overall SIMD-quintile mixing score, in turn, while retaining age, sex, and joint age-sex excess discordance, revealed substantial heterogeneity in the strength of domain-quintile mixing as a predictor of cluster scale (Supplementary Figure 7, ZTNB components; Supplementary Table 2 for the corresponding geographic-spread hurdle component, which is omitted from the heatmap for the same imbalance/quasi-separation reasons described for Figure 3). The direction of the domain-quintile mixing → cluster size and → geographic spread associations was consistently positive across all domains, mirroring the pooled SIMD-quintile result, but the magnitude varied by domain in a way that paralleled the domain heterogeneity observed for deprivation as an exposure (Supplementary Figure 5). Education- and crime-quintile mixing, which were the domains most positively associated with within-cluster mixing in the deprivation-as-exposure analyses, were among the strongest domain-side predictors of cluster scale, whereas access- and housing-quintile mixing — the domains that behaved as outliers in the deprivation-as-exposure analyses — produced more attenuated mixing-predictor estimates. In Supplementary Table 2 the crime and education domains report point estimates without cluster-robust standard errors, confidence intervals, or p-values: the window-clustered sandwich variance estimator failed numerically for those two hurdle fits because the Hessian (the "bread" of the sandwich) was effectively singular under the heavy outcome imbalance and the near-collinear `<domain>_domain_excess_mixing_z` predictor configuration; the maximum-likelihood point estimates remain valid. Age and joint age-sex mixing predictors remained strongly positive across all domain specifications, and sex mixing remained close to null or weakly negative, indicating that the demographic-mixing predictors of cluster scale are not confounded with the choice of socioeconomic domain. This symmetry between the deprivation-as-exposure and the mixing-as-predictor analyses reinforces the conclusion that SIMD domains capture qualitatively different social geographies — not only of who appears in deprived clusters, but of how mixing across deprivation strata maps onto reconstructed cluster scale.

### 3.8. Wave-specific mixing-predictor effects

Stratifying the mixing-predictor count models by epidemic wave (Figure 4, ZTNB components; Supplementary Table 1 for the corresponding geographic-spread hurdle component) showed that the pooled SIMD-quintile and age-mixing associations with larger and more geographically dispersed clusters were directionally stable across waves but varied in magnitude. The largest positive SIMD-quintile mixing effects on positive cluster size and positive geographic spread were observed in waves with the deepest cluster recruitment and most active community transmission (Alpha, Delta, BA.1, BA.2), echoing the per-wave deprivation-as-exposure pattern (Figure 2). Later Omicron subwaves (BA.4, BA.5, BQ.1) produced smaller and noisier mixing-predictor effects, consistent with their smaller per-wave cluster counts and higher singleton fractions (Table 1) and with the more rapid variant turnover that shortened the effective cluster-resolution window. The negative association of joint SIMD-age-sex profile excess discordance with cluster scale was observed in most waves, indicating that the fine-grained recombination signature is not specific to any single epidemic phase. The geographic-spread hurdle component is shown for the overall SIMD analysis in Figure 3B (OR 22.1, 95% CI 19.0–25.8 per 1 SD SIMD-quintile excess discordance); the corresponding per-wave estimates are reported in Supplementary Table 1 rather than as a heatmap. The hurdle binary outcome is `datazones_gt1`, which is 1 for 88% of clusters — heavily imbalanced. Combined with the strong `simd_excess_mixing_z` predictor, the binomial logistic component approaches quasi-separation in some waves, producing implausibly large adjusted odds ratios (~29,000 in the Alpha wave, 95% CI 6,450–129,898) that obscure heatmap interpretation; the magnitudes should be read as evidence of strong direction rather than as interpretable effect sizes. Wave-specific mixing-predictor estimates should be interpreted alongside the per-wave deprivation-as-exposure estimates, because they describe complementary features of the same set of clusters observed under wave-specific surveillance and immunity conditions.

### 3.9. Sensitivity analyses

The hurdle/ZTNB formulation captured patterns invisible to a log-linear comparator: log-linear estimates for SIMD on cluster size and geographic spread were both close to 1.0 (geometric mean ratios 0.992 and 1.001), reflecting averaging over the hurdle and positive-count components (Supplementary Figure 3). The size-offset cluster-size positive model produced an SIMD count ratio (0.925, 95% CI 0.868–0.985) almost identical to the primary estimate (0.926), indicating that the negative association is not explained by differences in available sequence pool. Winsorisation at the 99th percentile attenuated the positive cluster-size SIMD estimate (count ratio 0.952, 95% CI 0.900–1.006) but did not abolish the negative positive geographic-spread estimate (count ratio 0.889, 95% CI 0.837–0.944), and left the size-adjusted positive geographic-spread estimate weakly positive. Index-case SIMD as exposure did not reproduce the mean-cluster-SIMD positive-count associations (positive cluster size count ratio 0.996, 95% CI 0.957–1.036; positive geographic spread count ratio 0.989, 95% CI 0.951–1.028), and gave attenuated mixing estimates, indicating that the mean-cluster exposure captures whole-cluster deprivation composition rather than only the deprivation of the earliest detected case. Health-board clustered standard errors widened many count-outcome confidence intervals (positive cluster size 95% CI 0.664–1.293; positive geographic spread 95% CI 0.577–1.255) but left the directions and the size-adjusted positive geographic-spread result unchanged, and left age and sex mixing point estimates clearly different from zero. The approximately non-overlapping window subsample reduced precision but reproduced the qualitative directional pattern across all outcomes.

---

## 4. Discussion

We used a process-based, threshold-free pairwise compatibility model [31] to construct national SARS-CoV-2 genomic clusters across Scotland between 2020 and 2023, and we asked whether area deprivation, after adjustment for lineage, calendar time, and local surveillance conditions, is associated with the size, geographic spread, and sociodemographic composition of those clusters. Three substantive findings emerge.

First, more deprived clusters were not generally larger or more geographically dispersed. The negative pooled association of deprivation with positive cluster size and positive geographic spread, and its sign reversal in the size-adjusted geographic-spread model, indicates that deprivation is associated with a particular cluster-size profile (smaller, but spread across somewhat more datazones when conditioning on size), rather than with the simple "deprivation amplifies transmission cluster scale" story that case-based inequalities literature might suggest [16,17,19]. The much larger associations of surveillance covariates — test positivity, sequencing intensity, and local incidence — with cluster scale support a long-standing concern in genomic epidemiology that the apparent structure of inferred clusters depends on which cases were sequenced where and when [40,41]. In a national surveillance setting with non-uniform sequencing intensity across areas, time, and lineages [1,40], conflating the size of reconstructed clusters with the true intensity of local transmission risks misattributing surveillance gradients to social ones.

Second, deprivation is more clearly visible in the composition of clusters than in their size. The age-mixing association (+1.66 pp per 1 SD deprivation), the negative sex-mixing association (−0.78 pp), and the slight positive joint sociodemographic profile-mixing association (+0.48 pp) are consistent with a picture in which more deprived clusters include a more age-heterogeneous and somewhat more sociodemographically mixed set of cases, while remaining strongly assortative on SIMD quintile (a pattern visible in the observed-versus-expected mixing matrices). Greater age mixing in more deprived clusters fits prior evidence on multi-generational household structures, higher household crowding, and broader intergenerational contact patterns in more deprived neighbourhoods [16,18,37], as well as observations of differential ability to isolate by socioeconomic position [38,63]. The reduced sex mixing — that is, more sex-homogeneous clusters — could plausibly arise from occupational settings with skewed sex composition (e.g. care and cleaning roles, transport, hospitality, warehousing) that have repeatedly been identified as high-exposure environments during the pandemic [15,18,64]. We emphasise, however, that we cannot infer specific transmission settings from these patterns alone, only that the composition of more deprived clusters differs measurably from what would be expected of size-matched clusters drawn from the same lineage and time window.

Third, the SIMD domains behaved heterogeneously across both mixing and count outcomes. Education and crime deprivation were associated with greater domain-quintile mixing, whereas access and housing deprivation were associated with less. Access deprivation, in particular, was a consistent outlier across the demographic mixing outcomes. This is biologically plausible: SIMD's geographic-access domain heavily weights distance to GP, primary school, post office, retail, and petrol stations [65], and is therefore largely a rural and small-town deprivation measure, with a fundamentally different geography from urban income, employment, or housing deprivation. Recent Scottish work has explicitly highlighted the limitations of using SIMD as a single composite construct for COVID-19 inequalities analyses [62], and similar concerns have been raised for other deprivation indices in pandemic contexts [63,64,66]. Our results provide a sharply genomic-epidemiology-relevant version of this argument: which SIMD domain one uses materially affects what one concludes about deprivation and recent transmission.

The wave-specific analyses extend this caution further. The strongest negative deprivation associations occurred during Delta — a period of intense local transmission, high test positivity, and substantial age-graded vaccination rollout [13,67] — while Omicron subwaves were more variable. Because surveillance intensity, sequencing strategy, and the population susceptibility landscape changed substantially across these waves [1,40,62], it is difficult to disentangle a stable "deprivation effect" from interactions with the changing surveillance and immunity environment. The wave-specific estimates should therefore be read as evidence that any single pooled deprivation estimate averages over substantively different epidemic regimes.

The mixing-predictor analyses point to a fourth, complementary finding: clusters with higher SIMD-quintile and age-band excess discordance were detected as substantially larger and more geographically dispersed than otherwise similar clusters from the same lineage and time window. The magnitudes are striking — a 3.5-fold positive cluster-size association and a 22-fold geographic-spread hurdle association per 1 SD higher SIMD excess discordance — and are consistent with a bridging-transmission account in which clusters that mix across socioeconomic strata reach larger numbers of cases and span wider geographies. This pattern fits the broader contact-pattern literature, in which population-level mixing is strongly age-assortative within settings but bridging events (workplaces, transport, mixed-age households, multi-setting superspreading) can connect otherwise weakly linked subpopulations [31,29,38], and is consistent with empirical superspreading work showing that mixed-venue events seed disproportionately large transmission clusters [31,32]. The age-mixing association reinforces this: clusters with greater age bridging than the lineage-window baseline plausibly capture inter-household and inter-generational transmission, which has been repeatedly linked to deprivation through household structure and occupational exposure [16,18,58,68]. The negative association of joint SIMD-age-sex profile excess discordance with cluster scale — opposite in sign to the marginal SIMD and age associations — indicates that excess fine-grained sociodemographic recombination beyond what marginal SIMD and age mixing already explain does not buy additional cluster scale and may instead mark more fragmented chains of transmission. The domain- and wave-specific extensions of the mixing-predictor models (Supplementary Figure 7 and Figure 4) parallel the deprivation-as-exposure heterogeneity in two ways. First, education- and crime-quintile mixing produced the strongest domain-side mixing-predictor effects on cluster scale, mirroring the same two domains' positive associations with within-cluster mixing in the deprivation-as-exposure analyses (Supplementary Figure 5), while access- and housing-quintile mixing produced attenuated mixing-predictor effects, again as outliers. Second, the SIMD-quintile and age-mixing positive mixing-predictor effects were largest in waves with the deepest cluster recruitment (Alpha, Delta, BA.1, BA.2) and attenuated in late Omicron subwaves, paralleling the wave-specific pattern in the deprivation-as-exposure count models. This symmetry strengthens the central message that domain and wave heterogeneity are features of how social structure intersects with transmission and surveillance rather than artefacts of one specific modelling direction. Two interpretational caveats apply. First, the mixing-predictor estimates are associational rather than causal: cluster scale and cluster composition are jointly shaped by transmission, sampling, and ascertainment processes, and the direction "more bridging → larger clusters" should not be read mechanistically without further work. Second, although excess discordance is constructed as observed-minus-expected within lineage × window strata, residual size dependence cannot be fully ruled out at the largest values of the count outcomes, so the SIMD-mixing estimates are best understood alongside the within-cluster composition analyses (Figure 1, bottom row) rather than in isolation. Taken with the rest of the paper, however, these results sharpen the central message: cluster scale in genomic surveillance reflects both the deprivation profile of the cases involved and the bridging structure of contact and transmission that drew them into the same cluster.

These findings sit at the intersection of three literatures that have not been well linked at scale in Scotland. The Scottish genomic literature has focused on introductions [6], institutional outbreaks [7,8,11], variant severity [12,13], and travel policy [14], with little attention to population-wide cluster structure and inequalities. The Scottish inequalities literature has documented gradients in infection and outcomes but has rarely incorporated genome-informed cluster definitions [16,19]. Preprint work has begun to combine Scottish sequences with deprivation and spatial covariates [71], but to our knowledge no published nation-wide analysis has used a scalable recent-linkage clustering framework to ask how deprivation and surveillance jointly shape cluster outcomes across multiple waves. Our use of EpiLink [31] addresses the principal methodological barrier — providing graded, threshold-free compatibility scores that scale to national volumes — and our explicit modelling of surveillance covariates addresses the principal interpretive barrier — the dependence of cluster scale on sequencing intensity.

This study has several limitations. First, inferred genomic clusters represent recent-transmission compatibility, not direct transmission [31]. We chose this framing precisely because direct transmission inference is intractable at national scale and because EpiLink's compatibility scores have a clear mechanistic interpretation, but our cluster-level estimates should be read as descriptions of observed cluster structure rather than as estimates of directly inferred transmission. Second, our analysis is descriptive and associational; we do not claim causal interpretation for the deprivation estimates. Third, SIMD is an area-level measure and may not capture individual-level socioeconomic position, particularly in heterogeneous urban neighbourhoods. Fourth, sliding analytic windows induce dependence across cluster rows; we addressed this with window-clustered standard errors and an approximately non-overlapping window sensitivity, but residual dependence may persist. Fifth, we treat sequencing intensity as both a confounder and a substantive surveillance variable; the strong dominance of surveillance covariates in count outcomes underscores that the line between adjustment and ascertainment is genuinely blurred in genomic surveillance, and our conclusions about deprivation are therefore conditional on the surveillance regime in which the data were generated. Sixth, the cluster-size positive-component dispersion estimate reaches its upper bound under the ZTNB, indicating that the heavy right tail is not fully captured by the model; winsorisation analyses suggest that the substantive conclusions for geographic spread are robust to this, but the positive cluster-size association is more tail-sensitive. Seventh, large sample size makes small effects statistically precise; we have emphasised effect-size and direction over p-values.

For practice, our results have three implications. First, claims that "more deprived areas have larger genomic clusters" should be interpreted cautiously in surveillance datasets where sequencing intensity is itself socially patterned, because the apparent magnitude of clusters is heavily shaped by surveillance and epidemic-intensity covariates. Second, inequalities-focused genomic-epidemiology analyses should report results separately by SIMD (or equivalent) domain rather than using only the composite index, because composite indices can mask domain-level contrasts as large as the overall effect. Third, where inequalities-related questions about cluster composition are of interest (e.g. age mixing, joint sociodemographic mixing), excess-discordance-style measures conditional on lineage and time provide a more interpretable signal than cluster size or geographic spread alone.

---

## 5. Conclusion

Across more than 280,000 SARS-CoV-2 sequences spanning three years and eight Pango-lineage-defined epidemic waves in Scotland, mean cluster-level area deprivation was not associated with larger or more geographically dispersed genomic clusters after adjustment for lineage, calendar time, and local surveillance conditions. Cluster scale was substantially more strongly shaped by surveillance and epidemic-intensity covariates than by deprivation. Deprivation-related signals were clearest in the within-cluster sociodemographic composition — particularly age mixing — and varied markedly across SIMD domains and epidemic waves. Reversing the modelling direction, within-cluster excess mixing across SIMD quintiles and age bands was itself a strong predictor of larger and more geographically dispersed clusters, consistent with a bridging-transmission interpretation in which mixed-strata clusters are detected as larger in surveillance. These findings argue against a single monotonic "deprivation amplifies transmission" interpretation drawn from cluster scale alone, support the use of domain-disaggregated and wave-stratified analyses in genomic-epidemiological inequalities work, and reinforce the importance of explicit surveillance adjustment in any cluster-based equity claim.

---

## Figures

The four main figures are organised around the two complementary lines of inquiry described in the Introduction. Figures 1 and 2 cover the deprivation-as-exposure line; Figures 3 and 4 cover the excess-mixing-as-predictor line.

**Figure 1.** Deprivation as exposure: overall effects on cluster size, geographic spread, and within-cluster sociodemographic mixing. Top row: adjusted odds ratios (binary hurdle components) and zero-truncated negative-binomial (ZTNB) count ratios (positive components) for cluster size and geographic spread (number of distinct datazones), per 1 SD higher covariate, across 193,112 clusters from 134 sliding three-week windows at primary Leiden resolution γ = 0.3. Bottom row: adjusted percentage-point differences in observed-minus-expected pairwise discordance for SIMD quintile, age band, sex, and joint SIMD-age-sex profile excess mixing, among 84,067 non-singleton clusters. SIMD deprivation effects are modest and largely negative for count components; among mixing outcomes, deprivation is positively associated with age mixing and joint-profile mixing, negatively associated with sex mixing, and near null for SIMD-quintile mixing. Surveillance and epidemic-intensity covariates show substantially larger positive associations with count components. `figures/fig1_deprivation_overall.pdf`

**Figure 2.** Deprivation as exposure: wave-specific effects on cluster outcomes. Per-wave adjusted estimates of mean cluster SIMD deprivation on the four count-model components (cluster size hurdle, positive cluster size, geographic spread hurdle, positive geographic spread), for the eight dominant Pango-lineage waves modelled (B.1.177, Alpha, Delta, BA.1, BA.2, BA.4, BA.5, BQ.1). Delta shows the clearest negative pooled association; BA.2 and BA.4 show contrasting positive positive-count associations; earlier waves and BA.5/BQ.1 are more heterogeneous. `figures/fig2_deprivation_wave_specific.pdf`

**Figure 3.** Excess mixing as predictor: overall effects on cluster scale. Adjusted ratios per 1 SD higher excess-mixing score for the four mixing predictors (SIMD-quintile, age-band, sex, joint SIMD-age-sex profile excess discordance), entered jointly with SIMD deprivation, the four surveillance covariates, lineage fixed effects, and the calendar B-spline. Three panels are shown: (A) cluster size positive ZTNB count ratio (n = 84,067), (B) geographic spread hurdle odds ratio (n = 84,067; multi- vs single-datazone among non-singletons), and (C) geographic spread positive ZTNB count ratio (n = 74,010). The cluster-size hurdle is omitted because mixing scores are undefined for singletons (see Methods). SIMD-quintile and age-band excess discordance are strongly positively associated with cluster scale; joint-profile excess discordance is negatively associated. `figures/fig3_mixing_overall.pdf`

**Figure 4.** Excess mixing as predictor: wave-specific effects on cluster scale (ZTNB components). Per-wave adjusted count ratios of the four excess-mixing scores on the two ZTNB count-model components, for the eight dominant Pango-lineage waves. Two panels: (A) positive cluster size (n = 84,067) and (B) positive geographic spread (n = 74,010). Cells are coloured on a shared symmetric ratio scale (capped at ratio 5) and annotated with the raw count ratio per 1 SD higher excess-mixing score. The corresponding geographic-spread hurdle component is reported in Supplementary Table 1 rather than as a heatmap: the hurdle outcome (`datazones_gt1`) is 1 for 88% of non-singleton clusters — heavily imbalanced — and combines with the strong `simd_excess_mixing_z` predictor to drive the binomial logistic component toward quasi-separation in some waves, producing implausibly large odds ratios (~29,000 in Alpha) that obscure heatmap interpretation; the overall SIMD value of OR 22.1 (95% CI 19.0–25.8) is shown in Figure 3B. Positive SIMD-quintile and age-band mixing-predictor effects on the ZTNB components are directionally stable across waves but largest in Alpha, Delta, BA.1, and BA.2; late Omicron subwaves (BA.4, BA.5, BQ.1) show attenuated and noisier estimates. `figures/fig4_mixing_wave_specific.pdf`

Supplementary Figures 1–10 cover outcome and mixing distributions, observed–expected matrices, domain-specific extensions of both lines of inquiry, the size-adjusted geographic-spread sensitivity, the wave × domain demographic-mixing heatmap, and the log-linear-vs-hurdle/ZTNB comparators. Filenames follow the `supp_figN_*` convention written by `make_figures.py`. Supplementary Table 1 (companion to Figure 4) and Supplementary Table 2 (companion to Supplementary Figure 7) report the geographic-spread hurdle results that were omitted from those heatmaps for the imbalance / quasi-separation reasons described above. Full panel-by-panel captions for every supplementary figure and the embedded supplementary tables are in `figures/part1_supplementary_files.md`; longer narrative descriptions are in the accompanying `part1_results_and_figures_description.md`.
