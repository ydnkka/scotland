# Part 4 Analysis Documentation

## Alpha Emergence, Super-Seeding, and the F5 vs L2 Counterfactual: Would Earlier Lockdown Have Changed Scotland's Pandemic Trajectory?

---

## 1. Overview and Research Questions

Part 4 is a genomic case study examining how Alpha emerged in Scotland at the junction between post-first-lockdown reopening, the autumn B.1.177 wave, the five-tier framework, and the second national lockdown. The analysis is organised around four linked questions:

1. **What epidemic context did Alpha enter?** Describe the transition from Route-map phase 3 easing, through B.1.177 dominance and decline, into the five-tier framework and second lockdown.
2. **Who and where was Alpha initially associated with in the genomic data?** Characterise the early Alpha clusters by geography, age, sex, area deprivation, vaccination status, and test reason while treating these as descriptive surveillance signals rather than causal risk factors.
3. **Was there evidence for a key super-seeding transition?** Identify whether the early GGC chain appears to have expanded through a focal high-amplification event or through multiple unrelated introductions.
4. **Could earlier or more timely restrictions have changed the trajectory?** Fit logistic growth models to S:N501Y frequency under F5 and L2 conditions and project counterfactual scenarios for Alpha dominance and B.1.177 burden.

The analysis uses the same Leiden 0.3 / nextclade_qc = "good" dataset as Parts 1–3, supplemented by the raw Nextclade TSV for mutation-level frequency tracking.

---

## 2. Policy Context

| Code | Period | Dates | Intensity |
|---|---|---|---|
| P3 | Route map phase 3 | 2020-07-10 → 2020-10-01 | 30 |
| T1 | Pre-tier tightening | 2020-10-02 → 2020-11-01 | 55 |
| **F5** | **Five-tier framework** | **2020-11-02 → 2021-01-04** | **65** |
| **L2** | **Second lockdown** | **2021-01-05 → 2021-04-01** | **95** |
| SL | Stay local — Level 3 | 2021-04-02 → 2021-04-25 | 65 |

Route-map phase 3 was the main post-first-lockdown easing period in the dataset: hospitality, tourism accommodation, retail and schools reopened during this phase, while face coverings and localised restrictions remained in place. B.1.177 became the dominant lineage during and after this reopening period and carried the autumn epidemic into pre-tier tightening and F5.

The five-tier framework was Scotland's primary policy response for the 9 weeks between the first confirmed Alpha importation (4 November 2020) and the second lockdown (5 January 2021). F5 therefore sits at a critical point in the epidemic: B.1.177 was beginning to decline, but Alpha was entering and then expanding within the same social and surveillance context. Scotland's L2 was announced on 4 January and took effect 5 January — coinciding with England's third national lockdown. At the point L2 was imposed, Alpha already constituted 60% of all sequenced cases in Scotland.

---

## 3. The Alpha Seeding Chain

### 3.1 Genomic index case

**Cluster:** `W016|B.1.1.7|R0.3|S0`  
**Sequence:** `Scotland/QEUH-B3184F/2020`  
**Date:** 4 November 2020  
**Location:** Glasgow City, Greater Glasgow and Clyde  
**Test reason:** Symptomatic essential worker  
**QC status:** Good  

This single sequence carries the complete canonical Alpha VOC haplotype at confirmed good QC: all 7 Alpha-defining Spike mutations (S:N501Y, S:P681H, S:A570D, S:T716I, S:S982A, S:D1118H, S:D614G), the N protein changes (N:D3L, N:R203K, N:G204R, N:S235F), and the ORF8 truncation (Q27\*, R52I, Y73C). Zero private labelled mutations — a clean importation, not a sequence that accumulated additional changes during local spread.

### 3.2 Cryptic GGC spread: 4 November – 7 December 2020

The 5 weeks following the index case show a chain of small clusters confined almost entirely to Greater Glasgow and Clyde:

| Cluster | n | Dates | HBs | Key local authorities | Evidence of linkage |
|---|---|---|---|---|---|
| W016/S0 | 1 | 4 Nov | 1 | Glasgow City | Index case |
| W017/C1 | 5 | 4–13 Nov | 2 | Glasgow, N. Lanarkshire, W. Dunbartonshire | 1 seq = same as W016/S0 |
| W018/C2 | 10 | 10–20 Nov | 2 | Adds E. Renfrewshire, Renfrewshire | 4/5 W017 seqs present (80%) |
| W019/C1 | 9 | 10–18 Nov | 2 | Same geography | 4/5 W017 seqs also here (80%) |
| W019/C2 | 8 | 20–27 Nov | 2 | Adds Argyll & Bute, E. Dunbartonshire | 8 seqs shared with W018 |
| W020/C2 | 19 | 20 Nov–4 Dec | 4 | Adds Falkirk | 8 seqs shared with W019 |
| W021/C1 | 33 | 25 Nov–11 Dec | 5 | 11 local authorities | 17/19 W020 seqs present |

The 80% sequence overlap between W017 and both W018 and W019 is direct evidence these are the same transmission chain captured across successive sliding windows — not independent importations. The chain doubles approximately every 8 days within GGC, consistent with the Scotland-wide F5 growth rate estimate.

During this entire period (4 Nov – 7 Dec), S:N501Y frequency in Scotland remained below 4%, and the cluster chain was geographically contained. The five-tier framework was in operation throughout, with GGC under Level 3–4 restrictions.

### 3.3 Initial demographic and area profile

The early Alpha signal was not demographically neutral. Because clusters are inferred in overlapping windows, the summed cluster-size values are cluster memberships across phase windows, not unique person counts. The table therefore reports both summed cluster memberships and unique sequences. Demographic interpretation should prioritise the unique-sequence columns, while the cluster-weighted columns describe which demographic profiles dominated the inferred cluster signal. The source table is `part4_alpha_phase_demographic_summary.csv`.

| Phase | Windows | Clusters | Summed memberships | Unique sequences | Cluster-weighted signal | Unique-sequence profile |
|---|---|---:|---:|---:|---|---|
| Cryptic GGC chain | W016–W021 | 15 | 96 | 51 | Mean age 52.7; 75+ dominant; female-dominant; SIMD 1 dominant | Mean age 51.1; 75+ 19.6%; female 54.9%; SIMD 1 47.1%; unvaccinated 100% |
| Multi-region expansion | W022–W024 | 113 | 577 | 291 | Mean age 49.3; 75+ dominant; female-dominant; SIMD 1 dominant | Mean age 49.5; 75+ 18.9%; female 54.0%; SIMD 1 28.5%; unvaccinated 95.5% |
| F5/L2 bridge | W025 | 102 | 458 | 458 | Mean age 45.0; 60–64 dominant; female-dominant; SIMD 1 dominant | Mean age 45.0; 75+ 12.7%; female 53.5%; SIMD 1 26.2%; unvaccinated 95.9% |

The index case was a symptomatic essential worker, but the subsequent early chain is better described as a more-deprived, still largely unvaccinated Central Belt signal with older-age representation, rather than as an occupationally specific outbreak. This distinction matters: the data support questions about demographic association and possible high-exposure settings, but they do not identify a workplace, care home, school, or household event without additional epidemiological linkage data.

### 3.4 Candidate super-seeding transition: the 8 December expansion

On 8 December 2020, S:N501Y frequency jumps from 3.2% (1 Dec) to **17.7%** — a 5.5× increase in a single window. This corresponds to the simultaneous appearance of Alpha clusters in multiple health boards for the first time:

- `W022|B.1.1.7|R0.3|C6`: n=67, 8 health boards, 16 local authorities
- `W022|B.1.1.7|R0.3|C1`: n=17, 4 health boards
- Additional W022/W023 Alpha clusters spreading beyond GGC

The Dec 8 window marks the transition from contained GGC spread to multi-regional community transmission. The overlap evidence points to a candidate high-amplification transition rather than a purely diffuse rise: `W022|B.1.1.7|R0.3|C6` shares 60 sequences with `W023|B.1.1.7|R0.3|C4` (95.2% of W022/C6), and `W023|B.1.1.7|R0.3|C4` shares 86 sequences with `W024|B.1.1.7|R0.3|C1` (93.5% of W023/C4 and 98.9% of W024/C1). In other words, the apparent explosion is not just a frequency artefact; it is a persistent transmission cluster observed across successive windows.

The precise exposure event cannot be named from sequence data alone. The defensible interpretation is that a previously localised Alpha chain crossed a threshold into multi-regional spread in early December, consistent with a super-seeding process: a small number of linked transmission chains generating a disproportionate share of onward dissemination. F5 did not prevent or detectably slow this transition, but the available data cannot prove whether a specific gathering, institution, occupational setting, or travel-linked introduction caused it.

### 3.5 National superspreading: focal F5/L2 bridge cluster

**Cluster:** `W025|B.1.1.7|R0.3|C25`  
**Size:** 135 sequences, 120 datazones, 10 health boards, 24 local authorities  
**Dates:** 28 December 2020 – 8 January 2021 (spans F5/L2 boundary)  

This cluster represents the peak of the pre-L2 Alpha wave and directly bridges the F5→L2 policy transition. By its dates, Alpha was already in every major Scottish health board. The second lockdown was imposed into an epidemic where Alpha was already nationally established.

---

## 4. Growth Rate Analysis

### 4.1 Model: positive-test weighted binomial GLM

S:N501Y frequency per sliding window was modelled with logistic growth, fitted separately for the F5 phase and the L2 phase using a binomial generalised linear model (GLM) via `statsmodels`. We fitted a binomial GLM to the sequenced S:N501Y counts, with additional window-level weighting by confirmed positive-test volume. The response is [n_with_S:N501Y, n_without_S:N501Y] per window; the predictor is days since the anchor date (3 November 2020).

Each window is weighted by its total positive-test count (`wn_positive_tests`), normalised to mean 1. This weights the fitted sequenced proportion by the confirmed epidemic volume represented by that window; it does not by itself adjust for sequencing coverage. Sequencing coverage is evaluated separately in a sensitivity GLM that includes the proportion of positives sequenced (`n_seqs / wn_positive_tests`) as an additional covariate. Confidence intervals are 95% Wald CI on the slope parameter.

B.1.177 decline under L2 is fitted as an exponential OLS on log(frequency). This model is retained as a descriptive comparator for the displacement phase.

| Phase | Growth rate r (/day) | 95% CI r | Doubling time | 95% CI | Pseudo-R² |
|---|---|---|---|---|---|
| Alpha under F5 | 0.0812 | 0.0711–0.0913 | **8.5 days** | 7.6–9.7d | 0.884 |
| Alpha under L2 | 0.0661 | 0.0619–0.0702 | **10.5 days** | 9.9–11.2d | 0.933 |
| B.1.177 under L2 | −0.0654 (decline) | — | halving **10.6 days** | — | 0.944 (OLS R²) |

### 4.2 Growth-model sensitivity

The source table is `part4_growth_model_sensitivity.csv`. Three S:N501Y GLM specifications were compared:

| Model | Formula | F5 r/day | L2 r/day | L2/F5 | L2 slower | F5 doubling | L2 doubling |
|---|---|---:|---:|---:|---:|---:|---:|
| Unweighted binomial GLM | logit(freq) ~ days | 0.0789 | 0.0694 | 87.9% | 12.1% | 8.8d | 10.0d |
| Positive-test weighted GLM | logit(freq) ~ days, weighted by positive tests | 0.0812 | 0.0661 | 81.4% | 18.6% | 8.5d | 10.5d |
| Coverage-adjusted GLM | logit(freq) ~ days + proportion sequenced | 0.0865 | 0.0433 | 50.0% | 50.0% | 8.0d | 16.0d |

The direction of inference is stable across specifications: Alpha grew more slowly in L2 than in F5. The magnitude is model-dependent. The unweighted binomial GLM estimates a modest 12% reduction, the positive-test weighted primary model estimates a 19% reduction, and the coverage-adjusted model estimates a larger 50% reduction after explicitly adjusting for changing sequencing intensity. The primary analysis therefore reports the middle specification and treats the coverage-adjusted model as a sensitivity that shows the conclusion is not weakened by accounting for sequencing coverage.

### 4.3 Interpretation

The L2 growth rate is 81% of the F5 rate in the primary positive-test weighted model — **L2 reduced Alpha's logistic growth rate by 19%**. This is the central quantitative finding used for the counterfactual analysis, while the sensitivity table gives the plausible range under alternative treatment of window weighting and sequencing coverage.

B.1.177, in contrast, declined rapidly under L2 with a halving time of 10.6 days. The divergence between Alpha's continued growth and B.1.177's sharp decline under the same policy regime (L2, intensity = 95) confirms that Alpha's sustained growth was driven by its transmission advantage (~50–75% above WT) rather than by permissive policy conditions.

---

## 5. Counterfactual Analysis

### 5.1 Scenarios

Three counterfactual scenarios are compared against the actual timeline, using the L2-phase growth rate from the counterfactual switch date onward. Projections use the primary positive-test weighted binomial GLM rates (F5: r = 0.0812/day; L2: r = 0.0661/day).

| Scenario | Switch date | Alpha % on 5 Jan | Alpha reaches 50% |
|---|---|---|---|
| **Actual (F5 → L2 on 5 Jan)** | 5 Jan 2021 | **51%** | **5 Jan** |
| L2 from 8 Dec (explosion) | 8 Dec 2020 | 40% | 12 Jan |
| L2 from 2 Dec | 2 Dec 2020 | 38% | 19 Jan |
| **L2 from 2 Nov (immediate)** | 2 Nov 2020 | **28%** | **26 Jan** |

### 5.2 Key finding: policy difference was modest for Alpha

Even under the most aggressive counterfactual — lockdown imposed immediately when the index case was detected on 2 November — Alpha would still have reached 28% by the date L2 was actually imposed (5 January), and would have reached 50% dominance only three weeks later (26 January vs 5 January actual). Even the mid-range scenario (L2 from 8 December, the date of the multi-region explosion) delays 50% dominance by only one week (12 January).

The fundamental reason is the 19% growth rate difference: the gap between F5 (intensity 65) and L2 (intensity 95) was insufficient to contain a variant with Alpha's transmission advantage. Once Alpha was seeded in Scotland, its eventual dominance was strongly constrained by virology and repeated opportunities for onward spread, not simply by whether Scotland used a tiered framework or a full lockdown.

However, "not preventing dominance" is not the same as "making no difference." The counterfactuals indicate that earlier restrictions could plausibly have changed the timing and burden of the transition even if they did not reverse the ordering of variants. The most realistic policy question is therefore not whether F5 versus L2 could have stopped Alpha altogether, but whether earlier L2-level controls could have reduced the size of the B.1.177 background epidemic, delayed the Alpha crossing point, or reduced the scale of the early December super-seeding transition.

### 5.3 Where earlier L2 would have made a difference: B.1.177

While earlier L2 would not have prevented the Alpha wave, it would have significantly reduced the **B.1.177 burden** during November–December 2020. Under L2, B.1.177 (S:A222V) declined with a halving time of 10.6 days. Hospital occupancy during F5 was driven primarily by B.1.177 cases (occupancy peaked at ~1,241 beds on 2 November as T1 ended) and only began its steep L2-driven decline after 5 January. 

Under counterfactual L2 from 2 November, B.1.177 would have been at approximately 1/16th of its actual frequency by late December — preventing the F5-period hospital burden and entering L2 from a much lower baseline. The Alpha-driven peak occupancy (2,049 beds on 21 January) might still have occurred, but the preceding weeks of B.1.177 strain would have been substantially reduced.

### 5.4 The S:A222V second rise

An unexpected finding is the reappearance of S:A222V (the B.1.177 signature Spike mutation) at up to 35% frequency during the NN/OM periods (August–November 2021). This is **not residual B.1.177** — it is convergent acquisition of S:A222V in Delta sub-lineages, principally AY.4.2.2 (n=7,711 in NN/OM) and AY.4.2 (n=4,572). AY.4.2 was the UK's autumn 2021 "Delta plus" variant of concern. The same Spike mutation that defined B.1.177's European summer 2020 wave reappears 13 months later in a phylogenetically unrelated lineage, illustrating convergent evolution at an immune/fitness-relevant site.

---

## 6. Answer to the Research Question

**Would earlier imposition of L2 instead of F5 have significantly changed the trajectory of the Alpha pandemic in Scotland?**

The answer is **probably not enough to prevent Alpha dominance, but enough to matter for timing, B.1.177 burden, and possibly the scale of early amplification.**

The five-tier framework and the second lockdown differed by 19% in their effectiveness at slowing Alpha's logistic frequency growth in the primary positive-test weighted binomial GLM. Sensitivity analyses give a range from 12% slower in an unweighted binomial GLM to 50% slower in a coverage-adjusted GLM; the qualitative conclusion is unchanged. Even under the most aggressive scenario — lockdown from 2 November — Alpha would still have reached 50% dominance by 26 January, compared with 5 January actual. Under more realistic scenarios (lockdown from the December explosion), the delay is only one week. These estimates argue against a strong claim that earlier L2 would have stopped Alpha.

The more meaningful counterfactual benefit of earlier L2 would have been in the suppression of B.1.177 transmission during November–December 2020, potentially reducing F5-period hospital occupancy and entering the Alpha wave from a lower baseline. Earlier restrictions may also have reduced the number of infectious opportunities available to the early Alpha chain before the 8 December expansion. That effect is not fully captured by the simple logistic counterfactual, which models frequency growth after seeding but does not explicitly model stochastic extinction, importation pressure, or event-level transmission heterogeneity.

The genomic data therefore supports a timing-and-preventability interpretation. Alpha's trajectory in Scotland was set by four interacting factors: (1) social reopening and the autumn B.1.177 wave created a substantial transmission background after phase 3; (2) Alpha's intrinsic transmission advantage allowed growth even under L2; (3) a 5-week cryptic seeding phase in GGC meant detection lagged behind spread; and (4) an early December super-seeding transition appears to have converted local spread into national dissemination before L2 began. Timelier restrictions may have changed the burden and timing of this transition, but the available evidence does not support the stronger claim that they would have changed Alpha's eventual direction toward dominance.

---

## 7. Outputs

### Tables (`part4/tables/`)

| File | Contents |
|---|---|
| `part4_alpha_cluster_chain.csv` | All Alpha clusters ordered by first collection date, with HB/LA/test reason metadata |
| `part4_alpha_chain_overlaps.csv` | Pairwise sequence overlaps between early Alpha clusters (chain evidence) |
| `part4_alpha_phase_demographic_summary.csv` | Early Alpha phase summaries including summed cluster memberships, unique sequences, cluster-weighted profiles, unique-sequence age/sex/SIMD/vaccination profiles, and dates |
| `part4_alpha_clusters_weekly.csv` | Alpha cluster sizes by window |
| `part4_mutation_trajectories.csv` | S:N501Y, S:A222V, S:P681H and other key mutations per window |
| `part4_counterfactual_projections.csv` | Four scenario projections of S:N501Y frequency by date |
| `part4_growth_params.csv` | Fitted logistic growth parameters for F5 and L2 phases |
| `part4_growth_model_sensitivity.csv` | Unweighted, positive-test weighted, and coverage-adjusted S:N501Y GLM sensitivity comparison |
| `part4_scotland_hospital.csv` | Scotland total daily positive cases and hospital occupancy |
| `part4_window_period_map.csv` | Window → policy period mapping |

### Figures (`part4/figures/`)

| File | Contents |
|---|---|
| `fig1_alpha_seeding_chain.png` | Alpha cluster sizes by date with HB colouring; health boards active per week |
| `fig2_n501y_explosive_rise.png` | S:N501Y frequency with logistic fits, Dec 8 annotation, index case marker |
| `fig3_counterfactual_projections.png` | Four scenario frequency traces + hospital occupancy |
| `fig4_growth_rate_comparison.png` | Logit-linear fits for F5/L2; bar chart of growth/decline rates with doubling/halving days |
| `fig5_lineage_displacement.png` | Stacked lineage proportions + hospital occupancy overlay |

---

## 8. Script Architecture

| Script | Function |
|---|---|
| `stage1_data.py` | Loads parquet; builds cluster chain, demographic phase summaries, lineage composition, hospital data, seq ID files |
| `stage2_analysis.py` | Streams Nextclade TSV for mutation trajectories; fits positive-test weighted and sensitivity binomial GLMs for growth rates (with 95% CI); builds counterfactual projections |
| `stage3_figures.py` | Generates all five figures from pre-computed tables; fig2 and fig4 include 95% CI bands/error bars |

Run order: `stage1_data.py → stage2_analysis.py → stage3_figures.py`
