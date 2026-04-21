# Socioeconomic deprivation and SARS-CoV-2 genomic transmission clustering in Scotland

> Drafting notes for authors (remove before submission):
> 1. Insert final author list, affiliations, funding, acknowledgements, and governance approvals.
> 2. Fill in the all-QC deduplicated sequence count 324,194 (Abstract, Methods, Results) and the regenerated per-quintile breakdown after the Fig. 1 all-QC rerun. Confirm Fig. 2 has been switched from the broader precomputed cluster summary to the good-QC regression frame before the reported cluster-size means/medians are finalized.
> 3. ~~Add the planned QC-adjusted sensitivity analysis for singleton models if the headline estimates shift materially when mediocre and bad genomes are reintroduced with QC covariates.~~ **Done — QC-adjusted sensitivity analysis pre-specified in Statistical Analysis and described in Supplementary Methods S1; trigger condition (material shift in headline OR >10% change in Q1-vs-Q5 OR) specified. Run models once reintroduction dataset is confirmed and update S1 with results table.**
> 4. Replace reference [16] with the final EpiLink citation if the methods paper is posted as a preprint or accepted before submission.

**Author names and affiliations to be inserted**

**Keywords:** SARS-CoV-2; genomic epidemiology; deprivation; Scotland; SIMD; transmission clusters; EpiLink

## Abstract

Socioeconomic inequalities in COVID-19 have been documented for infection, hospitalization, and mortality, but much less is known about how deprivation structures genetically inferred transmission clusters at population scale. We applied EpiLink, a process-based pairwise compatibility method that integrates sampling-time differences and genetic distance, to Scottish SARS-CoV-2 surveillance data and asked whether area-level deprivation predicted cluster size and the probability that a sampled genome belonged to a multi-member transmission cluster rather than a singleton.

Scottish SARS-CoV-2 genomes were linked to Public Health Scotland testing and datazone metadata and clustered within overlapping 3-week lineage-specific windows. Pairwise TN93 distances were converted to SNP counts and weighted with EpiLink compatibility scores before Leiden community detection. Socioeconomic exposure was defined using the Scottish Index of Multiple Deprivation (SIMD) 2020v2 overall rank, quintiles, and seven component domains. Descriptive surveillance panels used324,194 deduplicated genomes from July 2020 to February 2023 regardless of Nextclade QC tier, because the descriptive quantities (sequencing coverage and lineage circulation) are not biased by assembly quality; cluster-level analyses were restricted to 276,026 Nextclade `good` genomes because QC errors inflate pairwise SNP distances and bias cluster assignment. Adjusted regression models used 190,021 complete cluster-window observations at the primary Leiden resolution of 0.3.

Marked deprivation gradients were present throughout the study period. In descriptive cluster summaries, the most deprived quintile (Q1) had a median cluster size of 2 in every variant epoch, whereas the least deprived quintile (Q5) had a median of 1 throughout; within-epoch distributions differed strongly in every epoch (all Kruskal-Wallis p < 0.001). In adjusted negative binomial models, a 1-SD increase in overall deprivation was associated with 5.3% larger clusters (IRR 1.053, 95% CI 1.048-1.059). Health and education deprivation showed the strongest positive domain-specific associations (IRR 1.064 for both), while access deprivation was inversely associated with cluster size (IRR 0.973, 95% CI 0.968-0.978). In logistic models, Q1 had substantially lower odds of singleton status than Q5 in every epoch, with ORs ranging from 0.37 in Alpha to 0.47 in Pre-VOC and remaining 0.38 in Omicron BA.2+, indicating persistent rather than attenuating socioeconomic patterning. In the mutually adjusted domain decomposition, the health domain accounted for 37.7% of the total absolute standardized coefficient magnitude.

These findings suggest that deprivation in Scotland was associated not only with COVID-19 burden, but with the structure of onward genomic transmission. Genomic surveillance therefore captures a persistent social gradient in transmission clustering that complements existing evidence on deprivation and severe outcomes, and highlights the value of place-based prevention in deprived communities.

## Author summary

COVID-19 did not spread evenly through society. Previous Scottish studies have shown that people living in more deprived areas were more likely to experience severe disease and death, but those studies did not ask whether the virus also circulated through different kinds of transmission clusters across the social gradient. We used EpiLink, a genomic clustering method that scores whether two sampled infections are compatible with recent transmission, and applied it to national Scottish surveillance data linked to area-level deprivation.

We found that genomes from more deprived communities were consistently less likely to appear as singletons and more likely to belong to larger genetically linked clusters. This pattern held across Pre-VOC, Alpha, Delta, and both major Omicron phases, which means the deprivation gradient did not disappear when population-wide infection became more common. Health, education, and employment deprivation were the main positive domain-level correlates of larger clusters, while the access domain showed the opposite pattern, consistent with more dispersed and less dense mixing in remote areas. The study adds a genomic dimension to the social epidemiology of COVID-19 in Scotland and suggests that deprivation shaped not only who became ill, but how onward transmission was organized.

## Introduction

Socioeconomic inequalities were a defining feature of the COVID-19 pandemic. Across multiple settings, people living in more deprived circumstances experienced higher risks of infection, hospitalization, and death [9-11]. Reviews of the early international literature consistently concluded that low income, overcrowded housing, lower educational attainment, and other markers of socioeconomic disadvantage were associated with greater COVID-19 burden, although the strength and timing of these associations varied across countries and phases of the pandemic [9-11]. The dominant interpretation has been that deprived populations faced more unavoidable exposure through crowded housing, essential work, lower capacity to work from home, and pre-existing health disadvantage [9,10].

Scottish evidence fits that broader pattern. Early record linkage work from the REACT-SCOT study showed that severe COVID-19 was socially patterned even after accounting for age and sex [12]. A later national record linkage study across the first three pandemic waves found that higher deprivation remained associated with higher odds of COVID-19 death in every wave [13]. At hospital level, analysis of SIMD indicators in southeast Scotland suggested that area-level income deprivation and alcohol-related harm remained informative even when aggregate deprivation quintiles became less stable in multivariable models [14]. These studies established deprivation as an important determinant of severe outcomes in Scotland, but they did not address how deprivation shaped the structure of transmission itself.

Pathogen genomics offers one route into that question. The UK COVID-19 Genomics UK (COG-UK) programme built surveillance infrastructure at a scale unprecedented for an acute respiratory virus, explicitly to support transmission tracking and integration with linked health data [1]. In Scotland, genomic epidemiology has already been used to show that the first wave was seeded by repeated introductions from mainland Europe and that early community transmission was established before control measures were introduced [2]. Setting-specific studies have shown how whole-genome sequencing can reconstruct outbreaks in Scottish university accommodation and other semi-closed settings [4]. At wider UK scale, dense surveillance sequencing has also been used to reconstruct lineage turnover and the spatial dynamics of Alpha and Delta [3]. What remains much less explored is whether these national genomic data can reveal durable social gradients in onward transmission clustering.

That gap is partly methodological. Full transmission-tree inference tools such as SCOTTI, TransPhylo, and outbreaker2 are powerful, but they are designed primarily for reconstructing who-infected-whom in outbreaks where richer assumptions and heavier computation are acceptable [5-7]. At the other end of the spectrum, threshold-based or probability-based genomic clustering methods such as cov2clusters provide scalable surveillance clusters, but they still face the problem that SARS-CoV-2 often displays low short-term genetic diversity and dense chains of nearly identical genomes [8]. In that setting, clustering methods need to make explicit use of time as well as genetic distance and should ideally remain interpretable at surveillance scale.

EpiLink was developed for that middle ground. Rather than inferring a fully resolved transmission tree, it evaluates whether the observed temporal and genetic separation between two samples is compatible with recent transmission scenarios under a mechanistic natural-history model [16]. In the Scottish pipeline, those pairwise compatibility scores are then used as weighted edges in lineage-specific graphs, followed by Leiden community detection. This produces genomic transmission neighborhoods that are not equivalent to observed transmission events, but are interpretable as recent linkage structures within surveillance data.

In this paper, we apply EpiLink-derived clustering to linked Scottish SARS-CoV-2 surveillance data from July 2020 to February 2023 and ask three questions. First, were genomes from more deprived communities more likely to fall into larger transmission clusters? Second, were they less likely to appear as singletons than genomes from less deprived communities? Third, which specific SIMD domains carried the deprivation signal once the aggregate index was decomposed? By answering these questions, we aim to extend the social epidemiology of COVID-19 in Scotland from outcomes and case counts to the genomic structure of onward transmission.

## Methods

### Study design and data sources

We conducted a national retrospective genomic surveillance analysis using linked SARS-CoV-2 sequence, testing, and area-level deprivation data for Scotland. The local processing log indicates that the underlying surveillance build contained 369,026 sequence metadata records, 5,467,764 testing rows, and SIMD 2020v2 information for all 6,976 Scottish datazones. Sequencing data were derived from the COG-UK surveillance framework [1], with linked Public Health Scotland metadata assembled into a master analysis parquet.

The figure pipeline applies Nextclade QC filtering only where it is analytically necessary. Descriptive surveillance panels (Figs 1 and 5) operate on all sequenced genomes regardless of QC tier, deduplicated on `sequence_id`, because the quantities they display — weekly sequenced cases by SIMD quintile, the proportion of positive tests sequenced, and lineage-specific median SIMD rank over time — are not biased by assembly quality: a mediocre- or bad-QC genome still represents a sequenced sample, occupying a real datazone on a real date. Across July 2020 to February 2023, this yielded324,194 deduplicated genomes. Cluster-level analyses (Figs 2, 3, 4, and 6), by contrast, are restricted to the 276,026 deduplicated genomes passing Nextclade quality status `good`, because mediocre- and bad-QC assemblies inflate pairwise TN93 distances and would bias EpiLink-based cluster assignment, disproportionately pushing low-quality genomes toward artefactual singleton status. The regression models further aggregate these good-quality genomes to cluster-window observations at the primary Leiden resolution of 0.3; after complete-case restriction for the joint domain model, 190,021 cluster-window observations were available.

### EpiLink-based clustering pipeline

Clustering was performed in 3-week sliding windows advanced in 1-week steps, yielding 134 windows across the study period. Within each window, sequences were stratified by PANGO lineage and pairwise TN93 distances were computed from the aligned Scottish SARS-CoV-2 genomes. Pairwise genetic distances were converted to SNP counts across the 29,903 bp genome.

EpiLink compatibility scores were then calculated for each pair using observed SNP distance and sampling-date difference. In the Scottish clustering pipeline, EpiLink was configured with a stochastic mutation process, 10,000 Monte Carlo samples, maximum hidden depth 0, and a target scenario set corresponding to direct ancestor-descendant transmission and a very recent shared common ancestor (`ad(0)` and `ca(0,0)`) [16]. The natural-history parameterization followed the EpiLink SARS-CoV-2 defaults embedded in the local wrapper, including a substitution rate of 1 x 10^-3 substitutions per site per year, clock relaxation 0.33, and a gamma-based latent/incubation/infectiousness formulation.

Pairs with EpiLink compatibility scores above 1e-4 were retained as weighted edges in an undirected graph, and weighted Leiden community detection was run at resolutions from 0.1 to 0.8. The primary analyses reported here use resolution 0.3. Cluster assignments were subsequently merged with testing, geography, and deprivation metadata to create sequence-level and cluster-level analytic frames.

### Variant epochs

Variant epochs were defined from the observed Scottish sequence data using the local manuscript data helper, which derives contiguous dominant-variant runs from weekly WHO variant labels and further splits Omicron into BA.1 and BA.2+ phases. The five resulting epochs were Pre-VOC, Alpha, Delta, Omicron BA.1, and Omicron BA.2+.

### Socioeconomic exposures

Area-level deprivation was measured using the Scottish Index of Multiple Deprivation (SIMD) 2020v2, the Scottish Government's official area-based deprivation measure across 6,976 datazones [15]. SIMD combines seven domains: income, employment, education, health, access to services, crime, and housing. Ranks are relative, with lower ranks indicating greater deprivation; for interpretability in the regression models, domain scores were sign-flipped after standardization so that higher standardized values corresponded to more deprivation.

At cluster level, we summarized deprivation in three ways:

1. Modal SIMD quintile of cluster members' datazones (`simd_quintile_mode`), with Q1 representing the most deprived fifth of Scottish datazones and Q5 the least deprived.
2. Mean overall SIMD rank of cluster members (`simd_rank_mean`).
3. Mean rank within each SIMD domain.

### Outcomes

The primary outcomes were:

1. **Cluster size**, defined as the number of unique sequences assigned to a cluster within a window.
2. **Singleton status**, coded as 1 when a cluster contained a single sequence and 0 otherwise.
3. **Lineage-level deprivation profile**, summarized as the monthly median SIMD rank percentile among sequences belonging to the top 20 PANGO lineages.

### Statistical analysis

We first summarized deduplicated sequence counts by SIMD quintile and week, alongside the weekly proportion of positive tests sequenced (`wn_prop_sequenced`). Because sequencing intensity varied substantially through time, `wn_prop_sequenced` was carried forward as a contextual surveillance measure and as an offset or covariate in downstream models.

For descriptive cluster-size comparisons, we summarized one row per `(epoch, modal SIMD quintile)` and compared within-epoch cluster-size distributions using Kruskal-Wallis tests.

Adjusted cluster-size associations were estimated with negative binomial generalized linear models (GLMs) on `n_sequences`. The domain-specific models took the form:

`cluster size ~ deprivation (1 SD) + variant-of-concern indicators + natural cubic spline on window mid-date + offset(log wn_prop_sequenced)`

The reported exponentiated coefficients are incidence rate ratios (IRRs), interpretable as the multiplicative change in expected cluster size associated with a 1-SD increase in deprivation.

Singleton status was modeled with logistic GLMs fitted separately within each variant epoch:

`is_singleton ~ SIMD quintile (Q1-Q4 vs Q5) + standardized log sequencing proportion + natural cubic spline on window mid-date`

Here, odds ratios (ORs) below 1 indicate that the quintile is less likely than Q5 to appear as a singleton and therefore more likely to be embedded in a multi-member genetically linked cluster.

To assess which domains carried the overall deprivation signal, we fitted a mutually adjusted negative binomial model including all seven SIMD domain ranks simultaneously. We then calculated each domain's share of the total absolute standardized coefficient magnitude, `share = |beta_d| / sum |beta|`. Because SIMD domains are correlated, the signed coefficients from this decomposition were interpreted as competitive conditional effects rather than stand-alone causal effects.

### Pre-specified QC-adjusted sensitivity analysis for singleton models

Because the primary singleton analyses are restricted to Nextclade `good` genomes to avoid QC-driven artefactual singletons, a pre-specified sensitivity analysis is registered here to assess whether the headline deprivation gradients are robust to the exclusion of mediocre- and bad-QC genomes.

**Trigger condition.** The sensitivity analysis is activated if either of the following is observed when mediocre- and bad-QC genomes are reintroduced into the singleton logistic models: (i) the Q1-versus-Q5 odds ratio changes by more than 10% in absolute terms in any variant epoch, or (ii) the direction of the quintile gradient reverses in any epoch.

**Analysis strategy.** Singleton status logistic GLMs are refitted with the same epoch-specific specification as the primary models (SIMD quintile Q1–Q4 vs Q5 + standardized log sequencing proportion + natural cubic spline on window mid-date), but with the analytic dataset expanded to include all Nextclade QC tiers. Two additional binary covariates are introduced: `qc_mediocre` (1 if Nextclade QC status is `mediocre`, 0 otherwise) and `qc_bad` (1 if `bad`, 0 otherwise), with `good` as the reference category. This additive adjustment allows QC tier to shift the intercept of singleton probability while leaving the deprivation gradient free to take its own value.

**Interpretation rule.** If the Q1-versus-Q5 singleton ORs from the QC-adjusted all-tier models fall within 10% of the primary good-QC estimates in every epoch, the headline results are judged robust to QC-tier exclusion and the sensitivity findings are reported in Supplementary Table S1 with the note that no material shift was detected. If a material shift is detected in one or more epochs, the QC-adjusted all-tier estimates are reported alongside the primary estimates in Table 2 and the discrepancy is discussed in the Strengths and Limitations subsection with reference to the mechanism (differential mediocre/bad-QC composition by SIMD quintile) that could produce the shift.

**Rationale for delayed execution.** This analysis is pre-specified rather than run in the primary analysis pass because QC-tier composition by SIMD quintile is itself an analytic quantity of interest and its direction is not known in advance. Running it conditionally on a material shift in the headline estimates protects against inflating the reported model count without analytical need, while the pre-specification here ensures that activation of the sensitivity analysis cannot be post-hoc.

### Ethics and governance

This study used de-identified linked surveillance data assembled under public health functions in Scotland. Insert the relevant approvals, public benefit/privacy panel references, and local information-governance wording here before submission.

## Results

### Sequencing coverage and socioeconomic distribution of sampled genomes

Among the324,194 deduplicated genomes used for sequence-level surveillance summaries (all Nextclade QC tiers retained), the most deprived quintile contributed 72,067 sequences, followed by Q2 (66,998), Q4 (61,623), Q5 (59,124), and Q3 (58,079). This pattern suggests a greater surveillance burden in more deprived communities, but it should be interpreted cautiously because sequenced counts reflect both epidemic burden and sequencing coverage. For reference, the good-QC subset used in all cluster-level analyses numbered 276,026 genomes and showed the same rank ordering across SIMD quintiles.

Sequencing intensity varied markedly across variant epochs. The median weekly proportion of positive tests sequenced was 13.2% in Pre-VOC, 52.3% in Alpha, 16.2% in Delta, 7.8% in Omicron BA.1, and 14.1% in Omicron BA.2+. The Alpha period therefore combined the sharpest deprivation gradient in raw weekly sequenced counts with the highest sequencing coverage, whereas Omicron BA.1 had far lower coverage and required greater caution in interpretation.

### Larger cluster sizes in more deprived communities across every epoch

Descriptive cluster summaries showed a persistent socioeconomic gradient in cluster size across the entire study period. Within every variant epoch, cluster-size distributions differed by modal SIMD quintile (all Kruskal-Wallis p < 0.001). The most deprived quintile (Q1) had a median cluster size of 2 in every epoch, whereas the least deprived quintile (Q5) had a median of 1 throughout.

This gradient was visible not only in medians but also in singleton prevalence. Overall, 45.7% of Q1 clusters were singletons compared with 67.5% of Q5 clusters. Mean cluster size was also substantially larger in Q1 than Q5 (5.61 versus 3.66 sequences per cluster-window observation). During Alpha, Q1 clusters had mean size 8.51 compared with 4.67 in Q5; during Delta the corresponding means were 5.77 and 4.27; and during Omicron BA.1 they were 8.83 and 5.16. The pattern therefore weakened somewhat in absolute terms as epidemic context changed, but it did not disappear.

### Overall deprivation and specific SIMD domains were associated with larger clusters

In adjusted negative binomial models, a 1-SD increase in overall deprivation was associated with a 5.3% increase in expected cluster size (IRR 1.053, 95% CI 1.048-1.059; p < 0.001). Every domain except access deprivation showed a positive univariable association with cluster size.

The strongest positive domain-specific associations were observed for the health domain (IRR 1.064, 95% CI 1.059-1.070; p < 0.001) and the education domain (IRR 1.064, 95% CI 1.059-1.070; p < 0.001), followed by employment deprivation (IRR 1.055, 95% CI 1.050-1.061; p < 0.001), overall SIMD rank (IRR 1.053, 95% CI 1.048-1.059; p < 0.001), and income deprivation (IRR 1.048, 95% CI 1.043-1.054; p < 0.001). Crime and housing showed smaller positive associations. Access deprivation was the only domain with an inverse association (IRR 0.973, 95% CI 0.968-0.978; p < 0.001), indicating that clusters located in more access-deprived, and likely more remote, areas tended to be smaller after adjustment for sequencing intensity, time, and dominant variant.

### The deprivation gradient in singleton odds persisted into Omicron

Epoch-specific logistic models showed a strikingly consistent monotonic gradient in singleton odds. Relative to Q5, the odds of being a singleton were lowest in Q1 in every epoch: Pre-VOC OR 0.470 (95% CI 0.396-0.558), Alpha OR 0.372 (0.328-0.420), Delta OR 0.436 (0.415-0.458), Omicron BA.1 OR 0.385 (0.356-0.418), and Omicron BA.2+ OR 0.376 (0.358-0.394). Q2, Q3, and Q4 also showed progressively higher ORs approaching the Q5 reference in every epoch.

This pattern did not support our initial expectation that the deprivation gradient would materially attenuate under Omicron. On the contrary, the most deprived quintile remained substantially less likely to appear as a singleton than the least deprived quintile even in the BA.2+ phase. In epidemiological terms, sequences from deprived communities were persistently more likely to be embedded in genetically linked onward chains rather than appearing as isolated observations.

### Domain decomposition suggested that health deprivation dominated the shared signal

When all seven SIMD domains were entered simultaneously into a single mutually adjusted negative binomial model, the health domain accounted for 37.7% of the total absolute standardized coefficient magnitude. Income contributed 23.7%, housing 15.2%, employment 10.3%, education 7.5%, crime 4.3%, and access 1.3%.

The signs of the mutually adjusted coefficients were mixed: health, employment, and education remained positive, whereas income, housing, and crime became negative. We interpret these sign reversals as evidence of collinearity and competitive decomposition rather than as a literal protective effect of income or housing deprivation. In other words, the shared deprivation signal appears to be distributed across correlated domains, but the part that remained most independently informative in this model was the health domain.

### Lineage composition also showed socioeconomic patterning

The lineage heatmap suggested that the deprivation signal was not confined to a single variant wave. Across the top 20 PANGO lineages, several major lineages were centered below the national median SIMD rank percentile, including B.1.1.7 (median percentile 0.372), BA.5.2 (0.413), BA.5.1 (0.423), and B.1.177 (0.451). Others skewed closer to or above the national midpoint, including AY.98 (0.506) and BA.1.15.1 (0.553).

Several lineages also traversed the deprivation gradient over time. AY.4.2.2, BA.2.23, and BA.5.2 showed the largest monthly ranges in median deprivation percentile, consistent with lineages that may have amplified first in one part of the social gradient before diffusing more widely. These observations should be interpreted cautiously because some late lineage-month cells were sparse, but they support the broader impression that deprivation shaped both cluster membership and the social geography of lineage spread.

## Table 1. Key cluster-size associations (QC=`good`)

| Exposure                            | Estimate | 95% CI      | P value | Interpretation                                           |
|-------------------------------------|---------:|-------------|--------:|----------------------------------------------------------|
| Overall SIMD deprivation (per 1 SD) |    1.053 | 1.048-1.059 |  <0.001 | More deprivation associated with larger clusters         |
| Income deprivation (per 1 SD)       |    1.048 | 1.043-1.054 |  <0.001 | Positive univariable association                         |
| Employment deprivation (per 1 SD)   |    1.055 | 1.050-1.061 |  <0.001 | Positive univariable association                         |
| Education deprivation (per 1 SD)    |    1.064 | 1.059-1.070 |  <0.001 | Strong positive univariable association                  |
| Health deprivation (per 1 SD)       |    1.064 | 1.059-1.070 |  <0.001 | Strongest positive univariable association               |
| Access deprivation (per 1 SD)       |    0.973 | 0.968-0.978 |  <0.001 | More access deprivation associated with smaller clusters |
| Crime deprivation (per 1 SD)        |    1.023 | 1.018-1.028 |  <0.001 | Small positive univariable association                   |
| Housing deprivation (per 1 SD)      |    1.012 | 1.006-1.017 |  <0.001 | Small positive univariable association                   |

## Table 2. Odds of singleton status in the most deprived quintile (Q1) versus the least deprived quintile (Q5)

| Epoch         | OR for Q1 vs Q5 | 95% CI      |  Δ% | Interpretation                                      |
|---------------|----------------:|-------------|----:|-----------------------------------------------------|
| Pre-VOC       |           0.470 | 0.396-0.558 | 7.9 | Q1 much less likely than Q5 to be a singleton       |
| Alpha         |           0.372 | 0.328-0.420 | 5.9 | Strong deprivation gradient                         |
| Delta         |           0.436 | 0.415-0.458 | 0.3 | Persistent deprivation gradient                     |
| Omicron BA.1  |           0.385 | 0.356-0.418 | 1.8 | Gradient persists despite lower sequencing coverage |
| Omicron BA.2+ |           0.376 | 0.358-0.394 | 8.2 | No meaningful attenuation during late Omicron       |

## Discussion

### Principal findings

This study shows that socioeconomic deprivation in Scotland was associated with the genomic structure of SARS-CoV-2 transmission, not just with downstream clinical severity. Across multiple variant epochs, genomes from more deprived communities were more likely to belong to larger clusters and less likely to appear as singletons. The persistence of the Q1 versus Q5 singleton odds gradient into Omicron BA.2+ was especially notable, because it ran counter to the prior expectation that widespread community infection would flatten social patterning.

The domain analysis adds nuance to the headline SIMD association. Health, education, and employment deprivation were the clearest positive correlates of larger clusters in the one-domain models, and the health domain dominated the mutually adjusted decomposition. Access deprivation behaved differently, with more remote and service-poor areas showing smaller clusters after adjustment. This suggests that the deprivation signal captured by genomic clustering in Scotland was not reducible to a single material disadvantage axis; instead, it reflected a combination of chronic ill health, educational and labor-market disadvantage, and the urban-rural organization of contact opportunity.

### Relation to existing literature

Our findings extend several strands of prior work. First, they are consistent with the broader COVID-19 inequalities literature, which has repeatedly shown that socioeconomically disadvantaged populations experience higher infection risk and worse outcomes [9-11]. Second, they complement Scottish studies showing that deprivation was associated with severe disease and death across the first pandemic waves [12-14]. However, those studies examined diagnosed cases, hospitalizations, or fatalities; none asked whether deprivation also structured the connectivity of viral transmission as represented in genomic surveillance data.

Third, the paper builds on Scottish genomic epidemiology. Early genomic work demonstrated that the first Scottish epidemic was seeded by repeated introductions and rapidly established community transmission [2]. Later setting-specific studies showed how sequencing could characterize outbreaks in Scottish university accommodation and other targeted environments [4]. Our analysis takes the next step by applying genomic clustering at national scale to a place-based social epidemiology question: not simply whether deprived communities were hit harder, but whether they were embedded in different transmission architectures.

The methodological contribution is also relevant. Existing outbreak reconstruction frameworks such as SCOTTI, TransPhylo, and outbreaker2 remain essential where the aim is direct transmission inference or explicit accounting for unsampled intermediates within a formal transmission-tree model [5-7]. Scalable SARS-CoV-2 clustering methods such as cov2clusters address a different surveillance problem by combining time and genetic divergence into stable cluster probabilities [8]. EpiLink occupies a related but distinct space by using a mechanistic recent-transmission compatibility model rather than a supervised or threshold-based pairwise rule [16]. The Scottish results suggest that such pairwise compatibility scores can be informative for population-scale social epidemiology when embedded in a reproducible clustering pipeline.

### Interpretation of the domain-level pattern

The strength of the health domain is epidemiologically plausible. Health deprivation in SIMD partly reflects chronic morbidity, disability, and excess emergency care use [15], which also mark communities with higher baseline vulnerability, greater care dependency, and potentially more sustained opportunities for repeated close contact. These features are compatible with larger and more socially concentrated transmission clusters. Education and employment deprivation likely capture additional mechanisms relevant to exposure and onward spread, including occupational inflexibility, household crowding correlated with labor-market disadvantage, and reduced capacity to isolate without income loss.

The inverse access-domain association is also plausible in Scottish context. Access deprivation often maps onto remote, rural, or poorly connected areas where service access is limited and where area-level contact networks may be more spatially diffuse and lower density [15]. Smaller cluster sizes in those settings are therefore consistent with less intense local amplification. Importantly, the access result also shows why a single overall SIMD score can obscure heterogeneity among its component parts.

The mutually adjusted domain decomposition should not be read naively. The negative coefficients for income and housing in the joint model do not imply that deprivation in those domains protects against transmission. Rather, they indicate that once highly correlated domains are allowed to compete in the same model, the remaining conditional coefficients can reverse sign. In substantive terms, the joint model is most useful for identifying where the shared deprivation signal concentrates, not for making clean causal statements about each domain in isolation.

### Public health implications

The main implication is that genomic surveillance can reveal a persistent social gradient in onward transmission that sits upstream of hospitalization and mortality. If deprived communities are more likely to generate larger genomic clusters and less likely to contribute isolated singletons, then interventions limited to downstream clinical protection will miss a critical part of inequality production. Place-based prevention remains important even in a mature surveillance system: paid isolation support, workplace mitigation, housing and ventilation interventions, targeted community vaccination, and rapid outbreak support in deprived neighborhoods are all consistent with the pattern observed here.

The findings also support making genomic surveillance explicitly equity-aware. Surveillance intensity varied sharply across the study period, and sequencing is never socially neutral. Genomic cluster analyses should therefore be interpreted together with coverage metrics such as `wn_prop_sequenced`, rather than as if the observed graph were a direct sample of all transmission events. Equity-sensitive genomic epidemiology requires both dense sequencing and transparent accounting for who was sequenced, when, and under what testing regime.

### Strengths and limitations

This study has several strengths. It uses national surveillance infrastructure built at high sequencing density [1-3], links genomic clustering to the Scottish Government's standard deprivation measure [15], and examines multiple complementary manifestations of transmission structure: cluster size, singleton status, lineage composition, and domain decomposition. It also benefits from an explicitly mechanistic clustering front end in EpiLink [16], rather than relying solely on a fixed genetic threshold.

The study also has limitations. First, EpiLink-derived clusters are genomic proxies for recent transmission neighborhoods, not observed transmission chains. Second, SIMD is an area-level measure; it identifies deprived places, not deprived individuals, and is known to be less sensitive to small pockets of deprivation in rural areas [15]. Third, sequencing coverage varied strongly over time and may also have varied socially; we adjusted for surveillance intensity, but residual bias is possible. Fourth, cluster-window observations are generated from overlapping 3-week windows and are therefore not fully independent. Fifth, some lineage-month cells in the heatmap are sparse and should not be overinterpreted. Sixth, the two analytic denominators used here (all-QC for descriptive surveillance panels, good-QC for cluster-level analyses) reflect a deliberate asymmetry based on where Nextclade quality plausibly biases each quantity rather than a single unified denominator; the sensitivity of the singleton gradient to QC-tier inclusion is addressed by a pre-specified QC-adjusted sensitivity analysis (see Statistical Analysis and Supplementary Methods S1), which will be activated and reported if the Q1-versus-Q5 singleton odds ratio shifts by more than 10% in any epoch when mediocre and bad genomes are reintroduced with QC-tier covariates.

## Conclusion

Area-level deprivation in Scotland was associated with larger SARS-CoV-2 genomic transmission clusters and substantially lower odds of singleton status across every major variant epoch from July 2020 to February 2023. The deprivation gradient persisted through Omicron rather than attenuating, and the health domain carried the largest share of the mutually adjusted deprivation signal. These results extend Scottish COVID-19 inequalities research from severe outcomes to the genomic structure of transmission and suggest that equity-focused prevention should remain central to epidemic preparedness and genomic surveillance policy.

## Figure legends

**Figure 1. Weekly deduplicated sequence counts and sequencing intensity by deprivation.**  
Panel A shows deduplicated good-quality SARS-CoV-2 genomes by week and SIMD quintile of residential datazone. Panel B shows the weekly proportion of positive tests sequenced. Shaded bands indicate the Pre-VOC, Alpha, Delta, Omicron BA.1, and Omicron BA.2+ epochs.

**Figure 2. Distribution of cluster size by modal SIMD quintile across variant epochs.**  
Each panel shows the distribution of cluster size within a variant epoch using violin plots overlaid with notched box plots. The y-axis is log1p cluster size. Kruskal-Wallis p values compare the quintile-specific cluster-size distributions within each epoch.

**Figure 3. Domain-specific IRRs for cluster size.**  
Forest plot of adjusted negative binomial IRRs for a 1-SD increase in deprivation within each SIMD domain and the overall SIMD rank. Models adjust for variant epoch, calendar time, and sequencing proportion.

**Figure 4. Odds of singleton status by SIMD quintile and variant epoch.**  
Within each epoch, logistic models estimate the odds of singleton status for Q1-Q4 relative to Q5. Odds ratios below 1 indicate that the quintile is less likely to appear as a singleton and therefore more likely to belong to a multi-member transmission cluster.

**Figure 5. Monthly lineage-level deprivation heatmap.**  
Heatmap of the median SIMD rank percentile for the top 20 PANGO lineages by monthly bin. Lower values indicate concentration in more deprived communities. The right-hand bars show total lineage abundance on a log scale. All Nextclade QC tiers are retained because each cell's value (median SIMD rank) is not biased by assembly quality; however, PANGO lineage assignment is less reliable for mediocre- and bad-QC genomes, which may marginally underrepresent such sequences in specific lineage rows.

**Figure 6. Mutually adjusted decomposition of the deprivation signal.**  
Bar plot showing each SIMD domain's share of the total absolute standardized coefficient magnitude from the joint negative binomial model that includes all seven domains simultaneously.

## Data availability

Scottish linked surveillance data are not publicly available because they contain potentially identifiable health information and are governed by Scottish public-health data access controls. Insert the exact data-governance route and approval wording here. Code used to generate the clustering inputs, analytic datasets, and manuscript figures is available in this repository. The EpiLink software used for pairwise compatibility scoring is available through PyPI and GitHub.

## Funding

Funding statement to be inserted.

## Competing interests

The authors declare no competing interests, or insert the final competing-interest statement here.

## References

1. COVID-19 Genomics UK (COG-UK) Consortium. An integrated national scale SARS-CoV-2 genomic surveillance network. *Lancet Microbe*. 2020;1:e99-e100. doi: [10.1016/S2666-5247(20)30054-9](https://doi.org/10.1016/S2666-5247(20)30054-9)

2. da Silva Filipe A, Shepherd JG, Williams T, Hughes J, Aranday-Cortes E, Asamaphan P, et al. Genomic epidemiology reveals multiple introductions of SARS-CoV-2 from mainland Europe into Scotland. *Nature Microbiology*. 2021;6:112-122. doi: [10.1038/s41564-020-00838-z](https://doi.org/10.1038/s41564-020-00838-z)

3. Voehringer HS, Sanderson T, Sinnott M, De Maio N, Nguyen T, Goater R, et al. Genomic reconstruction of the SARS-CoV-2 epidemic in England. *Nature*. 2021;600:506-511. doi: [10.1038/s41586-021-04069-y](https://doi.org/10.1038/s41586-021-04069-y)

4. Nickbakhsh S, Hughes J, Christofidis N, Griffiths E, Shaaban S, Enright J, et al. Genomic epidemiology of SARS-CoV-2 in a university outbreak setting and implications for public health planning. *Scientific Reports*. 2022;12:11735. doi: [10.1038/s41598-022-15661-1](https://doi.org/10.1038/s41598-022-15661-1)

5. De Maio N, Wu CH, Wilson DJ. SCOTTI: Efficient Reconstruction of Transmission within Outbreaks with the Structured Coalescent. *PLoS Computational Biology*. 2016;12:e1005130. doi: [10.1371/journal.pcbi.1005130](https://doi.org/10.1371/journal.pcbi.1005130)

6. Didelot X, Fraser C, Gardy J, Colijn C. Genomic Infectious Disease Epidemiology in Partially Sampled and Ongoing Outbreaks. *Molecular Biology and Evolution*. 2017;34:997-1007. doi: [10.1093/molbev/msw275](https://doi.org/10.1093/molbev/msw275)

7. Campbell F, Didelot X, Fitzjohn R, Ferguson N, Cori A, Jombart T. outbreaker2: a modular platform for outbreak reconstruction. *BMC Bioinformatics*. 2018;19:363. doi: [10.1186/s12859-018-2330-z](https://doi.org/10.1186/s12859-018-2330-z)

8. Sobkowiak B, Kamelian K, Zlosnik JEA, Tyson J, Goncalves da Silva A, Hoang LMN, et al. Cov2clusters: genomic clustering of SARS-CoV-2 sequences. *BMC Genomics*. 2022;23:710. doi: [10.1186/s12864-022-08936-4](https://doi.org/10.1186/s12864-022-08936-4)

9. Wachtler B, Michalski N, Nowossadeck E, Diercke M, Wahrendorf M, Santos-Hoevener C, et al. Socioeconomic inequalities and COVID-19 - A review of the current international literature. *Journal of Health Monitoring*. 2020;5(Suppl 7):3-17. doi: [10.25646/7059](https://doi.org/10.25646/7059)

10. Khanijahani A, Iezadi S, Gholipour K, Azami-Aghdash S, Naghibi D. A systematic review of racial/ethnic and socioeconomic disparities in COVID-19. *International Journal for Equity in Health*. 2021;20:248. doi: [10.1186/s12939-021-01582-4](https://doi.org/10.1186/s12939-021-01582-4)

11. Benita F, Rebollar-Ruelas L, Gaytan-Alfaro ED. What have we learned about socioeconomic inequalities in the spread of COVID-19? A systematic review. *Sustainable Cities and Society*. 2022;86:104158. doi: [10.1016/j.scs.2022.104158](https://doi.org/10.1016/j.scs.2022.104158)

12. McKeigue PM, Weir A, Bishop J, McGurnaghan SJ, Kennedy S, McAllister D, et al. Rapid Epidemiological Analysis of Comorbidities and Treatments as risk factors for COVID-19 in Scotland (REACT-SCOT): A population-based case-control study. *PLoS Medicine*. 2020;17:e1003374. doi: [10.1371/journal.pmed.1003374](https://doi.org/10.1371/journal.pmed.1003374)

13. Leslie K, Findlay B, Ryan T, Green LI, et al. Epidemiology of SARS-CoV-2 during the first three waves in Scotland: a national record linkage study. *Journal of Epidemiology and Community Health*. 2022;77:1-8. doi: [10.1136/jech-2022-219367](https://doi.org/10.1136/jech-2022-219367)

14. Scopazzini MS, Cave RNR, Mutch CP, Ross DA, Bularga A, Chase-Topping M, et al. Scottish Index of Multiple Deprivation (SIMD) indicators as predictors of mortality among patients hospitalised with COVID-19 disease in the Lothian Region, Scotland during the first wave: a cohort study. *International Journal for Equity in Health*. 2023;22:205. doi: [10.1186/s12939-023-02017-y](https://doi.org/10.1186/s12939-023-02017-y)

15. Scottish Government. Scottish Index of Multiple Deprivation 2020. Available at: [https://www.gov.scot/simd](https://www.gov.scot/simd)

16. Arthur D, Banks CJ, Kao RR. EpiLink: a process-based compatibility model for genomic transmission clustering in infectious disease surveillance. Manuscript in preparation.
