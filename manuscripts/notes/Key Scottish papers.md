## Key Papers on Genomic Transmission, Clustering, Phylogenetics, and Phylogeography in Scotland

### 1. Scotland-focused, peer-reviewed papers on genomic transmission, clustering, phylogenetics or phylogeography

1. da Silva Filipe et al. (2021)
    _Genomic epidemiology reveals multiple introductions of SARS-CoV-2 from mainland Europe into Scotland._
    Journal: Nature Microbiology
    Focus: National-scale introductions into Scotland; phylogenetics; phylogeography; early outbreak reconstruction. 
2. Li et al. (2021)
    _Genetic epidemiology of SARS-CoV-2 transmission in renal dialysis units – a high risk community-hospital interface._
    Journal: Journal of Infection
    Focus: Transmission clustering in Scottish renal dialysis units; hospital versus community acquisition.
3. Nickbakhsh et al. (2022)
    _Genomic epidemiology of SARS-CoV-2 in a university outbreak setting and implications for public health planning._
    Journal: Scientific Reports
    Focus: University of Glasgow outbreak; cluster structure; multiple introductions; outbreak reconstruction.
4. Cotton et al. (2023)
    _Investigation of hospital discharge cases and SARS-CoV-2 introduction into Lothian care homes._
    Journal: Journal of Hospital Infection
    Focus: Care-home introductions in Lothian; genomic investigation of possible hospital-seeded transmission.

### 2. Papers with Scottish genomic outbreak data, but not Scotland-only or not mainly national transmission mapping

1. Stirrup et al. (2021)
    _Rapid feedback on hospital onset SARS-CoV-2 infections combining epidemiological and sequencing data._
    Journal: eLife
    Focus: Hospital-onset infection analysis using sequencing plus epidemiology; includes Glasgow alongside Sheffield.
2. Stirrup et al. (2022)
    _Effectiveness of rapid SARS-CoV-2 genome sequencing in supporting infection control for hospital-onset COVID-19 infection: multicentre, prospective study._
    Journal: eLife
    Focus: UK multicentre hospital IPC genomics study with Scottish participation; not Scotland-only.

### 3. Scotland-focused genomic epidemiology papers adjacent to transmission mapping

1. McLachlan et al. (2024)
    _Evaluation of risk-based travel policy for the COVID-19 epidemic in Scotland: a population-based surveillance study._
    Journal: BMJ Open
    Focus: Importation risk, travel policy, and genomic evidence on variant importation into Scotland; closest to policy-oriented phylogeography.
2. Pascall et al. (2023)
    _The SARS-CoV-2 Alpha variant was associated with increased clinical severity of COVID-19 in Scotland: a genomics-based retrospective cohort analysis._
    Journal: PLOS ONE
    Focus: Scottish lineage-based genomics study; mainly severity rather than transmission clustering, but still uses Scottish genomic data and phylogenetic analysis.
3. Sheikh et al. (2021)
    _SARS-CoV-2 Delta VOC in Scotland: demographics, risk of hospital admission, and vaccine effectiveness._
    Journal: The Lancet
    Focus: Variant epidemiology in Scotland using viral genetic classification; includes demographic and deprivation-related analysis, though not primarily a clustering or phylogenetic paper.

### 4. Preprints / non-peer-reviewed work highly relevant to Scotland genomic transmission and socioeconomic factors

1. Lycett et al. (2021)
   _Epidemic waves of COVID-19 in Scotland: a genomic perspective on the impact of the introduction and relaxation of lockdown on SARS-CoV-2._
   Status: Preprint
   Focus: Scotland-wide genomic perspective on epidemic waves and lockdown effects; highly relevant phylodynamic/phylogeographic work.
2. Gamża et al. (preprint)
   _Spatial analysis of phylogenetic, population and deprivation data from Scottish SARS-CoV-2 outbreak reveals patterns of the community transmission._
   Status: Preprint / SSRN version
   Focus: Scottish SARS-CoV-2 sequences linked to deprivation and spatial data; directly relevant to socioeconomic patterning of genomic transmission.
3. Gamża et al. (preprint, later title/version)
   _Infector characteristics exposed by spatial analysis of SARS-CoV-2 sequence and demographic data analysed at fine geographical scales._
   Status: Preprint / arXiv version
   Focus: Fine-scale spatial, sequence, and demographic analysis in Scotland; likely the later version of the same project.

### 5. Best category for “socioeconomic factors + genetic data” specifically

The papers from the list above that fit this best are:

* Gamża et al. preprints — the most directly focused on deprivation/spatial patterning plus sequence data in Scotland.
* Sheikh et al. (Delta VOC in Scotland) — includes deprivation-related epidemiological analysis with variant classification.
* Li et al. — includes SIMD context, but deprivation is not the main analytical focus.

___ 

### Introduction

The COVID-19 pandemic prompted an unprecedented expansion of pathogen genomic surveillance and established SARS-CoV-2 as the clearest modern example of how viral sequence data can inform public health. At global scale, genomic analyses were used to reconstruct introductions, track variant emergence and spread, identify hospital and community outbreaks, and support phylodynamic inference on epidemic growth and replacement. In the United Kingdom, the COVID-19 Genomics UK (COG-UK) Consortium created one of the world’s largest integrated genomic surveillance systems, enabling near real-time linkage of viral genomes to epidemiological metadata across multiple epidemic phases.

In Scotland, genomic epidemiology played an important role from the earliest phase of the epidemic. Early national analyses showed that the first wave was seeded by numerous independent introductions, predominantly from mainland Europe, and demonstrated the value of combining viral genomes with travel and epidemiological data to distinguish importation from local spread. Subsequent Scottish studies used sequencing to investigate transmission in higher-risk or policy-relevant settings, including renal dialysis units, university outbreaks, hospital-onset infection, and care homes, while later work examined the role of travel policy in variant importation. Together, these studies established that genomic data could illuminate transmission pathways in Scotland, but they were largely focused on early introductions, specific outbreaks, or particular institutional settings rather than population-wide recent transmission structure across the full epidemic.

This pattern mirrors a broader international literature. Large-scale phylogenetic and phylodynamic analyses have been highly informative for reconstructing introductions, lineage turnover, and regional or international dissemination, but they are not always well matched to the problem of identifying recent transmission structure among contemporaneous cases within densely sampled epidemics. For SARS-CoV-2 in particular, low short-term sequence diversity, rapid epidemic growth, incomplete sampling, and the possibility of repeated importations of genetically similar viruses can all limit the resolution of phylogenies for inferring direct or near-direct transmission. As a result, many applied studies relied on simple SNP thresholds or ad hoc phylogenetic clustering rules, which are operationally convenient but can be insensitive to timing and may merge epidemiologically distinct events or split related ones depending on the threshold chosen.

To address these limitations, several groups developed probabilistic or semi-mechanistic frameworks that integrate genetic distances with sampling times, symptom onset, or serial interval information when assessing whether two cases are plausibly linked within a recent transmission chain. Examples include methods such as outbreaker2, transcluster, and A2B-COVID, which move beyond purely genetic thresholding by incorporating the temporal constraints imposed by infection and sampling processes. These approaches helped demonstrate that joint temporal-genetic modelling can improve recent-transmission inference, especially in settings such as hospitals or other defined outbreaks. However, most have been applied at relatively modest scale, within bounded settings, or as tools for outbreak investigation rather than as a basis for nation-wide clustering over hundreds of thousands of genomes spanning multiple waves and variants.

At the same time, a separate and extensive literature has shown that the burden of COVID-19 was socially patterned. Internationally, socioeconomic disadvantage has been associated with higher infection risk, greater exposure through occupation and household crowding, and worse clinical outcomes. Scotland followed this broader pattern: deprivation was associated with inequalities in infection, hospitalisation and mortality, and Scottish studies of Alpha and Delta highlighted the importance of demographic and contextual heterogeneity across epidemic phases. Yet this inequalities literature has generally been built from case counts, testing data, admissions, or deaths, rather than from genome-informed representations of recent transmission. Consequently, an important question remains unresolved: whether deprivation is associated merely with a higher burden of diagnosed disease, or whether it is also reflected in the structure of transmission, for example through participation in larger recent-transmission clusters or reduced probability of appearing as apparently isolated cases.

Within Scotland, very little published work has directly joined viral genetic data to deprivation at fine scale. The clearest example is the recent preprint literature linking Scottish SARS-CoV-2 sequences to deprivation and spatial covariates, which suggests that socio-spatial inequalities may shape inferred transmission patterns. However, this work remains limited relative to the scale of available Scottish genomic surveillance data, and there is still no published nation-wide analysis, to our knowledge, that uses a scalable recent-linkage framework to study how area deprivation relates to cluster size, singleton membership, and the persistence of these associations across successive variant-dominated epidemic epochs.

This gap matters for both epidemiological interpretation and public health practice. If deprivation is associated with larger recent-transmission clusters, that would suggest inequalities are expressed not only through differential susceptibility or access to testing, but through differences in contact structure, exposure contexts, or opportunities for onward spread. If, conversely, cases in less deprived areas are more often singletons, that may reflect distinct mobility, mixing, household, or ascertainment patterns. And if these relationships persist across periods dominated by different variants of concern, despite marked shifts in transmissibility, immunity, control measures and population behaviour, that would indicate a robust underlying social gradient in transmission opportunity rather than a transient feature of a single wave.

In this study, we address these questions using a national SARS-CoV-2 genomic dataset from Scotland comprising more than 350,000 sequences collected between 2020 and 2023 through COG-UK and linked to the Scottish Index of Multiple Deprivation (SIMD). We develop an interpretable pairwise linkage model that combines temporal and genetic distances to estimate compatibility with recent transmission, defined here as direct transmission or co-primary infection from a shared source. These compatibility scores define a weighted network over cases, which is sparsified by removing very weak edges before applying community detection to infer recent transmission clusters. Using this framework, we ask three questions: whether larger inferred clusters are associated with greater area deprivation; whether cases from less deprived areas are more likely to appear in singleton clusters; and whether any such gradient persists across distinct epidemic epochs characterised by successive dominance of major variants of concern. By integrating scalable genome-informed recent-transmission inference with deprivation at national scale, our study aims to connect Scotland’s genomic surveillance infrastructure to one of the central unresolved themes of the pandemic: the social patterning of transmission itself.