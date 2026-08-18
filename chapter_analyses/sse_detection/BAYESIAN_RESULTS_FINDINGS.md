# Bayesian SSE Detection Results: Findings Reference

This document summarises the consolidated Bayesian model results for the SSE detection analysis. It is intended as a working reference for writing the results section, not as final prose.

## Source Tables

The results below are based on the report-facing summary tables in:

```text
chapter_analyses/sse_detection/results/bayesian_outputs/consolidated_tables/
```

Main summary tables:

- `summary_table_1_diagnostics.csv`: model-level diagnostics.
- `summary_table_2_estimates.csv`: directional estimates, including fixed effects, intercepts, and random intercepts.
- `summary_table_3_random_effect_sds.csv`: directionless random-effect SDs and residual SDs.

Detailed consolidated tables:

- `mixing_logistic_consolidated_results.csv`
- `mixing_linear_consolidated_results.csv`
- `composition_logistic_consolidated_results.csv`
- `composition_linear_consolidated_results.csv`

## Interpretation Rules

Use the expanded models as the main results and the primary models as sensitivity checks. The expanded models include surveillance-context adjusters.

For mixing models, use expanded null-standardised entropy as the main inferential scale. It asks whether cluster composition is more or less mixed than expected under the within-window, cluster-size matched null model. Observed entropy is an absolute heterogeneity companion scale and should be reported as a sensitivity/descriptive result, especially where directions differ from the null-standardised estimates.

For composition models, categorical contrasts are relative to:

- Sex: Male.
- Age group: 25-64.
- SIMD: SIMD Q1.
- Urban/rural: Large Urban.
- Health board: Greater Glasgow and Clyde.

For logistic models, estimates are odds ratios or multiplicative odds. Values above 1 indicate higher candidate odds; values below 1 indicate lower candidate odds.

For linear models, estimates are beta coefficients on the burst-score or burden-score scale (0 to 1). Positive values indicate higher score; negative values indicate lower score.

Direction probabilities are posterior probabilities of the effect direction:

- `P Positive Direction`: `P(beta > 0 | data)` for linear models or `P(OR > 1 | data)` for logistic models.
- `P Negative Direction`: `P(beta < 0 | data)` for linear models or `P(OR < 1 | data)` for logistic models.
- `Direction Probability`: the larger of the two direction probabilities.

Direction probabilities are not p-values. They should be used to support effect-size interpretation, not as the main result.

## Overall Result

The strongest chapter-level finding is that the detection signal is structured by geography, socioeconomic context, demographic composition, surveillance intensity, policy period, and clade. The models suggest that detection is not explained by sampling intensity alone.

The mixing models show that sampled-population heterogeneity is associated with candidate detection, burst score, and burden score. The main mixing interpretation comes from null-standardised entropy, with observed entropy used to describe absolute heterogeneity. The composition models then show which population groups, places, periods, and clades drive those differences.

The most robust findings are:

- Null-standardised mixing entropy shows that several attributes differ from the within-window, cluster-size matched expectation.
- Observed mixing entropy is strongly associated with candidate detection and burst score, but it should be interpreted as absolute heterogeneity rather than the main null-adjusted mixing result.
- Health-board effects are large, especially in the composition models.
- Urban/rural composition is consistently associated with lower candidate odds and lower linear scores outside large urban areas.
- Female samples have lower detection signal than male samples.
- Age 15-24 has higher detection signal than 25-64, while younger children and older age groups generally have lower signal.
- Recent sequencing intensity is positively associated with candidate detection.
- Cumulative incidence is generally negatively associated with detection outcomes.
- Policy period and clade explain substantial residual heterogeneity, especially in the composition models.
- Clade random-effect SDs are generally larger than policy-period SDs, particularly in composition models.

## 1. Model Diagnostics

All models had 0 divergences and 8,000 posterior draws.

The mixing models are diagnostically clean:

- All mixing logistic models have `Diagnostic Status = OK`.
- All mixing linear models have `Diagnostic Status = OK`.
- Mixing models have `Max Rhat = 1.00`.
- Effective sample sizes are generally high.

The composition logistic models are acceptable:

- Primary and expanded composition logistic models have `Diagnostic Status = OK`.
- Both have 0 divergences.
- `Max Rhat` is 1.00 in the primary model and 1.01 in the expanded model.

The composition linear models need a caveat:

- Composition burden-score primary model: `OK`.
- Composition burden-score expanded model: `Warning`.
- Composition burst-score primary and expanded models: `Warning`.
- Warnings are driven by lower ESS and `Max Rhat` up to 1.02, not divergences.

Reporting implication: report composition linear findings, but phrase them more cautiously than mixing and composition logistic findings.

## 2. Health Board And Geographic Structure

### 2.1 Mixing: Health-Board Entropy

The main null-standardised health-board mixing result differs by outcome:

- Candidate odds are lower: OR 0.868, P(lower) = 1.000.
- Burst score is lower: beta -0.030, P(negative) = 1.000.
- Burden score is higher: beta 0.017, P(positive) = 1.000.

Observed health-board mixing is consistently associated with stronger absolute heterogeneity signal.

Expanded observed health-board entropy:

| Outcome | Estimate | 95% HDI | Direction |
| --- | ---: | --- | --- |
| Candidate odds | OR 1.147 | 1.091 to 1.209 | Higher, P = 1.000 |
| Burst score | beta 0.014 | 0.011 to 0.017 | Positive, P = 1.000 |
| Burden score | beta 0.006 | 0.001 to 0.011 | Positive, P = 0.993 |

Interpretation: null-standardised health-board entropy captures whether health-board mixing exceeds the size/window null expectation and should be the main inferential contrast. Observed geographic mixing captures absolute heterogeneity and is associated with higher candidate detection and burst/burden signal, so it should be presented as a companion result rather than collapsed with the null-standardised scale.

### 2.2 Composition: Health Board Contrasts

Health-board composition effects are substantial relative to Greater Glasgow and Clyde.

For candidate odds, the clearest higher signals are:

- Orkney: OR 3.815, 95% HDI 3.004 to 4.953, P(higher) = 1.000.
- Forth Valley: OR 1.129, 95% HDI 1.054 to 1.209, P(higher) = 1.000.
- Western Isles: OR 1.454, 95% HDI 0.992 to 2.075, P(higher) = 0.973.

Candidate odds are lower in:

- Borders: OR 0.833, P(lower) = 0.997.
- Fife: OR 0.917, P(lower) = 0.990.
- Grampian: OR 0.927, P(lower) = 0.989.
- Ayrshire and Arran: OR 0.924, P(lower) = 0.987.

Burst score differs strongly by board:

- Higher burst score: Orkney, Western Isles, Shetland, Lanarkshire, Ayrshire and Arran, and Highland.
- Lower burst score: Grampian, Borders, Fife, Lothian, and Tayside.

Examples:

- Orkney: beta 0.051, P(positive) = 1.000.
- Western Isles: beta 0.032, P(positive) = 0.977.
- Lanarkshire: beta 0.008, P(positive) = 1.000.
- Grampian: beta -0.027, P(negative) = 1.000.
- Borders: beta -0.021, P(negative) = 1.000.

Burden score is lower for most boards relative to Greater Glasgow and Clyde:

- Orkney: beta -0.219, P(negative) = 1.000.
- Grampian: beta -0.055, P(negative) = 1.000.
- Tayside: beta -0.045, P(negative) = 1.000.
- Shetland: beta -0.043, P(negative) = 0.990.
- Dumfries and Galloway: beta -0.036, P(negative) = 1.000.
- Borders: beta -0.033, P(negative) = 1.000.
- Lothian: beta -0.024, P(negative) = 1.000.
- Ayrshire and Arran: beta -0.019, P(negative) = 1.000.

Interpretation: health board is one of the clearest sources of structured variation. Candidate odds, burst score, and burden score do not all move in the same direction, suggesting that different detection summaries capture different aspects of geographic signal.

## 3. Socioeconomic And Urban/Rural Structure

### 3.1 Mixing: SIMD Entropy

The main null-standardised SIMD entropy result is outcome-dependent:

- Candidate odds are slightly lower: OR 0.959, P(lower) = 0.977.
- Burst score is lower: beta -0.0045, P(negative) = 0.999.
- Burden score is higher: beta 0.0096, P(positive) = 1.000.

Observed SIMD entropy is one of the strongest absolute heterogeneity effects.

Expanded observed SIMD entropy:

| Outcome | Estimate | 95% HDI | Direction |
| --- | ---: | --- | --- |
| Candidate odds | OR 1.294 | 1.209 to 1.391 | Higher, P = 1.000 |
| Burst score | beta 0.038 | 0.035 to 0.040 | Positive, P = 1.000 |
| Burden score | beta -0.011 | -0.016 to -0.006 | Negative, P = 1.000 |

This is a strong and coherent observed-scale result for candidate odds and burst score, but burden score moves in the opposite direction on the observed scale. The observed and null-standardised SIMD results should therefore be interpreted as different contrasts: absolute SIMD heterogeneity versus mixing relative to the size/window null expectation.

### 3.2 Mixing: Urban/Rural Entropy

The main null-standardised urban/rural entropy result is directionally similar to the observed scale for all three outcomes:

- Candidate odds: OR 1.069, P(higher) = 0.994.
- Burst score: beta 0.010, P(positive) = 1.000.
- Burden score: beta -0.006, P(negative) = 0.987.

Observed urban/rural entropy is also associated with detection signal.

Expanded observed urban/rural entropy:

| Outcome | Estimate | 95% HDI | Direction |
| --- | ---: | --- | --- |
| Candidate odds | OR 1.083 | 1.030 to 1.139 | Higher, P = 0.999 |
| Burst score | beta 0.009 | 0.007 to 0.012 | Positive, P = 1.000 |
| Burden score | beta -0.007 | -0.012 to -0.002 | Negative, P = 0.997 |

Interpretation: urban/rural mixing is the clearest case where the main null-standardised scale and the observed absolute-heterogeneity scale tell a similar story.

### 3.3 Composition: SIMD Contrasts

SIMD composition effects are weaker than SIMD mixing effects.

For candidate odds relative to SIMD Q1:

- SIMD Q4 is slightly higher: OR 1.047, P(higher) = 0.972.
- SIMD Q5 is slightly lower: OR 0.956, P(lower) = 0.967.
- SIMD Q2 and Q3 are weak or uncertain.

For burst score:

- SIMD Q3: beta 0.003, P(positive) = 0.963.
- SIMD Q4: beta 0.004, P(positive) = 0.985.
- SIMD Q5: beta 0.003, P(positive) = 0.935.

For burden score:

- SIMD Q2: beta -0.0075, P(negative) = 0.998.
- SIMD Q5: beta -0.0085, P(negative) = 0.999.
- SIMD Q3 and Q4 are weak.

Interpretation: the strongest SIMD signal is at the mixing-entropy level rather than as a simple monotonic composition contrast.

### 3.4 Composition: Urban/Rural Contrasts

Urban/rural composition is clearer than SIMD composition.

For candidate odds relative to large urban areas:

- Accessible town: OR 0.883, P(lower) = 1.000.
- Remote rural: OR 0.834, P(lower) = 1.000.
- Remote town: OR 0.795, P(lower) = 1.000.
- Accessible rural: OR 0.944, P(lower) = 0.969.
- Other urban is weak or uncertain.

For burst score, all non-large-urban categories are lower:

- Accessible rural: beta -0.009, P(negative) = 1.000.
- Accessible town: beta -0.011, P(negative) = 1.000.
- Other urban: beta -0.005, P(negative) = 0.998.
- Remote rural: beta -0.026, P(negative) = 1.000.
- Remote town: beta -0.023, P(negative) = 1.000.

For burden score, the strongest reductions are remote rural and remote town:

- Remote rural: beta -0.027, P(negative) = 1.000.
- Remote town: beta -0.026, P(negative) = 1.000.
- Other urban: beta -0.006, P(negative) = 0.993.
- Accessible town is moderately negative: P(negative) = 0.932.
- Accessible rural is weaker: P(negative) = 0.893.

Interpretation: samples from less urban settings generally show lower candidate odds and lower burst/burden scores relative to large urban areas. This is one of the more coherent composition findings.

## 4. Demographic Structure

### 4.1 Mixing: Age Entropy

The main null-standardised age entropy result differs from the observed scale:

- Candidate odds are lower: OR 0.849, P(lower) = 1.000.
- Burst score is lower: beta -0.020, P(negative) = 1.000.
- Burden score is weak or uncertain: beta -0.002, P(negative) = 0.679.

Observed age entropy is strongly associated with absolute heterogeneity in detection outcomes.

Expanded observed age entropy:

| Outcome | Estimate | 95% HDI | Direction |
| --- | ---: | --- | --- |
| Candidate odds | OR 1.106 | 1.052 to 1.162 | Higher, P = 1.000 |
| Burst score | beta 0.022 | 0.020 to 0.025 | Positive, P = 1.000 |
| Burden score | beta -0.013 | -0.018 to -0.008 | Negative, P = 1.000 |

Interpretation: absolute age heterogeneity is higher in candidate and high-burst clusters, but age mixing above the within-window, size-matched null expectation is associated with lower candidate odds and lower burst score.

### 4.2 Mixing: Sex Entropy

The main null-standardised sex entropy result is lower for candidate odds and burst score, with weak burden-score evidence:

- Candidate odds: OR 0.938, P(lower) = 0.996.
- Burst score: beta -0.005, P(negative) = 0.998.
- Burden score is weak or uncertain: beta -0.002, P(negative) = 0.700.

Observed sex entropy has strong linear effects but a weaker candidate-odds effect.

Expanded observed sex entropy:

| Outcome | Estimate | 95% HDI | Direction |
| --- | ---: | --- | --- |
| Candidate odds | OR 1.050 | 0.983 to 1.127 | Higher, P = 0.922 |
| Burst score | beta 0.016 | 0.013 to 0.019 | Positive, P = 1.000 |
| Burden score | beta -0.013 | -0.019 to -0.007 | Negative, P = 1.000 |

Interpretation: sex mixing should be presented carefully because the absolute and null-standardised contrasts point in different directions for candidate odds and burst score.

### 4.3 Composition: Sex

Female samples have consistently lower detection signal than male samples.

Expanded Female vs Male:

| Outcome | Estimate | 95% HDI | Direction |
| --- | ---: | --- | --- |
| Candidate odds | OR 0.956 | 0.929 to 0.983 | Lower, P = 1.000 |
| Burst score | beta -0.0057 | -0.0080 to -0.0034 | Negative, P = 1.000 |
| Burden score | beta -0.0056 | -0.0091 to -0.0022 | Negative, P = 1.000 |

Interpretation: sex composition has a consistent but modest association with all three detection outcomes.

### 4.4 Composition: Age Group

Age-group composition is very coherent.

Relative to age 25-64:

- Age 15-24 has higher detection signal:
  - Candidate odds: OR 1.105, P(higher) = 1.000.
  - Burst score: beta 0.016, P(positive) = 1.000.
  - Burden score: beta 0.006, P(positive) = 0.989.

- Ages 00-04 and 05-14 have lower detection signal:
  - Age 00-04 candidate odds: OR 0.860, P(lower) = 0.997.
  - Age 05-14 candidate odds: OR 0.887, P(lower) = 1.000.
  - Both groups have lower burst and burden scores.

- Age 65-74 is lower:
  - Candidate odds: OR 0.924, P(lower) = 0.983.
  - Burst score: beta -0.012, P(negative) = 1.000.
  - Burden score: beta -0.008, P(negative) = 0.963.

- Age 75+ is lower for linear outcomes, but weaker for candidate odds:
  - Candidate odds: OR 0.949, P(lower) = 0.922.
  - Burst score: beta -0.015, P(negative) = 1.000.
  - Burden score: beta -0.029, P(negative) = 1.000.

Interpretation: age composition is one of the most consistent demographic predictors. The 15-24 group stands out as higher signal, while children and older age groups are generally lower signal.

## 5. Surveillance Context

Surveillance-context effects should be presented after demographic, socioeconomic, and geographic effects. They help clarify whether the main detection signal is partly explained by sequencing intensity, epidemic intensity, policy timing, or lineage structure.

### 5.1 Window Sequencing Proportion

The clearest surveillance-context fixed effect is recent sequencing intensity.

Composition models:

- Candidate odds: OR 1.211, 95% HDI 1.185 to 1.246, P(higher) = 1.000.
- Burst score: beta 0.031, P(positive) = 1.000.
- Burden score: beta 0.038, P(positive) = 1.000.

Mixing models:

- Observed candidate odds: OR 1.156, P(higher) = 0.997.
- Observed burst score: beta 0.006, P(positive) = 0.944.
- Observed burden score is weak: beta 0.003, P(positive) = 0.668.
- Null-standardised burden score is positive: beta 0.015, P(positive) = 0.995.
- Null-standardised burst score is negative: beta -0.018, P(negative) = 1.000.

Interpretation: recent sequencing intensity is positively associated with candidate detection, especially in composition models and observed mixing logistic models. Linear mixing associations vary by outcome and scale.

### 5.2 Cumulative Incidence

Cumulative incidence generally has a negative association with detection outcomes in the composition models.

Composition models:

- Candidate odds: OR 0.843, P(lower) = 1.000.
- Burst score: beta -0.026, P(negative) = 1.000.
- Burden score: beta -0.045, P(negative) = 1.000.

Mixing models:

- Burden score is negative on both scales:
  - Null-standardised: beta -0.076, P(negative) = 1.000.
  - Observed: beta -0.069, P(negative) = 1.000.
- Observed burst score has no clear direction.
- Null-standardised burst score is positive: beta 0.022, P(positive) = 0.998.
- Candidate odds are weak or uncertain in the mixing models.

Interpretation: after accounting for composition and other covariates, cumulative incidence is generally associated with lower candidate odds and lower linear scores. This may indicate that the strongest detection signal is not simply a function of accumulated incidence.

### 5.3 Cumulative Sequencing Proportion

Cumulative sequencing proportion is weaker overall.

Composition models:

- Candidate odds: OR 0.975, P(lower) = 0.986.
- Burden score: beta -0.0056, P(negative) = 1.000.
- Burst score is weak: beta -0.0004, P(negative) = 0.686.

Mixing models:

- Candidate odds are weak or uncertain on both scales.
- Burden score is weak to moderately negative.
- Observed burst score is moderately positive: P(positive) = 0.959.

Interpretation: cumulative sequencing proportion is less central than recent window sequencing proportion.

## 6. Policy Period Effects

Policy-period effects are random intercepts in Table 2. They represent residual period-level deviations after fixed effects and should be interpreted as contextual temporal structure, not as causal policy effects.

### 6.1 Composition Logistic Candidate Outcome

Policy-period effects are strongest in the composition candidate model.

Higher candidate odds:

- T1: multiplicative odds 9.488, P(higher) = 1.000.
- F5: multiplicative odds 4.904, P(higher) = 1.000.
- P3: multiplicative odds 3.127, P(higher) = 0.998.

Lower candidate odds:

- PR: multiplicative odds 0.323, P(lower) = 0.999.
- FE: multiplicative odds 0.383, P(lower) = 0.998.
- L0: multiplicative odds 0.440, P(lower) = 0.992.
- L2: multiplicative odds 0.543, P(lower) = 0.970.
- L21: multiplicative odds 0.543, P(lower) = 0.971.
- SL: multiplicative odds 0.549, P(lower) = 0.968.

### 6.2 Composition Linear Outcomes

Composition burst score:

- Positive period deviations: F5, P3, and T1.
- Negative period deviations: FE, L0, and PR.
- PR is strongly negative: beta -0.246, P(negative) = 1.000.

Composition burden score:

- Positive period deviations: F5, L2, P3, SL, and T1.
- Negative period deviations: FE, L21, and PR.
- PR is strongly negative: beta -0.343, P(negative) = 1.000.
- F5 is strongly positive: beta 0.239, P(positive) = 1.000.

Composition linear period effects should be described with the diagnostic caveat because the relevant linear models include warnings.

### 6.3 Mixing Models

Policy-period effects are weaker in mixing models and mostly appear for burst score.

Observed mixing burst score:

- Positive deviations: F5, PR, and SL.
- Negative deviations: L0, L21, and NN.

There is little strong evidence for policy-period residual structure in mixing candidate odds or burden score.

Interpretation: policy period explains important residual temporal variation in composition models, particularly for candidate odds, but is more secondary in the mixing models.

## 7. Clade Effects

Clade effects are random intercepts in Table 2. They represent lineage-level residual deviations after fixed effects.

### 7.1 Composition Logistic Candidate Outcome

Clade effects are very strong in the composition candidate model.

Higher candidate odds:

- 22D: multiplicative odds 42.521, P(higher) = 1.000.
- 22B: multiplicative odds 12.061, P(higher) = 1.000.
- 21L: multiplicative odds 10.486, P(higher) = 1.000.
- 22A: multiplicative odds 9.679, P(higher) = 1.000.
- 21K: multiplicative odds 6.753, P(higher) = 1.000.
- 21J: multiplicative odds 3.857, P(higher) = 0.995.
- 21I: multiplicative odds 3.743, P(higher) = 0.991.
- 20I: multiplicative odds 2.435, P(higher) = 0.955.

Lower candidate odds:

- 21A: multiplicative odds 0.036, P(lower) = 0.999.
- 20B: multiplicative odds 0.050, P(lower) = 1.000.
- Rec.: multiplicative odds 0.065, P(lower) = 0.992.
- 22C: multiplicative odds 0.088, P(lower) = 0.980.
- 20E: multiplicative odds 0.298, P(lower) = 0.977.

Interpretation: candidate detection is highly lineage-structured after covariate adjustment.

### 7.2 Composition Linear Outcomes

Composition burst score:

- Lower residual burst signal: 19B, 20A, 20B, 20E, 20I, and 21A.
- Higher residual burst signal: 21K, 21L, 22A, 22B, 22D, 23A, and Rec.
- Largest positive burst signal: 23A, beta 0.461, P(positive) = 1.000.
- Strong positive burst signal also appears for 22D, beta 0.360, and 22B, beta 0.239.

Composition burden score:

- Lower residual burden signal: 20A, 20B, 20E, 20I, and 21A.
- Higher residual burden signal: 21L, 22B, 22D, and 22E.
- Largest positive burden signal: 22E, beta 0.710, P(positive) = 1.000.
- Strong negative burden signal: 21A, beta -0.450, P(negative) = 1.000.

Again, composition linear clade effects should be interpreted with the diagnostics caveat.

### 7.3 Mixing Models

Clade effects are less prominent in mixing models.

Mixing logistic candidate models do not show strong clade-specific residual deviations.

Mixing linear burden score shows some clade residual structure:

- Lower burden residuals: 20I and 21J.
- Higher burden residuals: 21L, 22B, and 22D.

Mixing linear burst score has smaller and less coherent clade deviations:

- Observed burst score is positive for 20I.
- Observed burst score is negative for 21K, 21L, and 22C.

Interpretation: clade effects are central for composition candidate detection and composition linear scores, but less central for mixing candidate detection.

## 8. Directionless Variance Components

Random-effect SDs and residual SDs are in Table 3 because they are directionless. They should be interpreted as heterogeneity magnitudes, not as positive or negative effects.

### 8.1 Composition Models

Clade SDs are larger than policy-period SDs in all composition models.

Composition logistic candidate model:

- Clade SD: multiplicative odds SD 9.025, 95% HDI 4.482 to 20.086.
- Policy-period SD: multiplicative odds SD 3.136, 95% HDI 2.096 to 5.474.
- Approximate random-effect variance share: clade 79%, policy period 21%.

Composition linear burst score:

- Clade SD: 0.268, 95% HDI 0.190 to 0.390.
- Policy-period SD: 0.159, 95% HDI 0.100 to 0.260.
- Approximate random-effect variance share: clade 74%, policy period 26%.
- Residual SD: 0.225.

Composition linear burden score:

- Clade SD: 0.371, 95% HDI 0.260 to 0.540.
- Policy-period SD: 0.196, 95% HDI 0.130 to 0.310.
- Approximate random-effect variance share: clade 78%, policy period 22%.
- Residual SD: 0.278.

Interpretation: lineage explains more residual heterogeneity than policy period in the composition models.

### 8.2 Mixing Models

Mixing logistic candidate models also show larger clade SDs than policy-period SDs.

Observed mixing logistic:

- Clade SD: multiplicative odds SD 1.278.
- Policy-period SD: multiplicative odds SD 1.129.
- Approximate random-effect variance share: clade 80%, policy period 20%.

Null-standardised mixing logistic:

- Clade SD: multiplicative odds SD 1.261.
- Policy-period SD: multiplicative odds SD 1.114.
- Approximate random-effect variance share: clade 82%, policy period 18%.

Mixing linear burden score:

- Clade SD dominates policy-period SD.
- Observed scale: clade SD 0.126 vs policy-period SD 0.016.
- Null-standardised scale: clade SD 0.135 vs policy-period SD 0.015.

Mixing linear burst score:

- Policy-period SD is similar to or larger than clade SD.
- Observed scale: policy-period SD 0.053 vs clade SD 0.048.
- Null-standardised scale: policy-period SD 0.031 vs clade SD 0.019.

Interpretation: clade heterogeneity is larger for candidate odds and burden score, while policy-period heterogeneity is relatively more important for burst score in the mixing models.

## 9. Primary vs Expanded Models

The expanded models are the main results because they account for surveillance context. Primary models should be used as sensitivity checks. Within the mixing domain, expanded null-standardised models are the main inferential scale, while expanded observed models are the absolute-heterogeneity companion scale.

Across the main results, the broad findings are stable:

- Mixing entropy effects remain visible after context adjustment, but observed and null-standardised scales answer different questions.
- Demographic composition effects remain coherent.
- Urban/rural composition effects remain coherent.
- Health-board effects remain substantial.
- Clade and policy-period residual structure remains visible after adjustment.

The most important caveat is not instability between primary and expanded models; it is the diagnostics warning in the composition linear models.

## 10. Suggested Results Structure

Recommended ordering:

1. State the modelling framework and outcome definitions.
2. Present model diagnostics using Table 1.
3. Present main fixed effects using Table 2 and forest plots.
4. Start with null-standardised mixing results because they describe broad heterogeneity relative to the size/window null expectation.
5. Then present composition results because they identify which groups drive the signal.
6. Organise substantive effects as:
   - Health board and geography.
   - Socioeconomic and urban/rural structure.
   - Demographic structure.
   - Surveillance context.
7. Present policy-period and clade random intercepts as contextual residual structure.
8. Present Table 3 as variance-component evidence that clade generally contributes more heterogeneity than policy period, especially in composition models.
9. Use observed mixing entropy and primary models as sensitivity/descriptive checks, not as the main inferential mixing result.

## 11. Suggested Figure And Table Strategy

Main text:

- Table 1: diagnostics summary.
- A compact fixed-effect forest plot for mixing entropy effects, showing null-standardised estimates as the main scale and observed estimates as the companion scale.
- Composition forest plots split by:
  - demographic contrasts,
  - socioeconomic/urban-rural contrasts,
  - health-board contrasts,
  - surveillance-context adjusters.
- A short table or panel for clade and policy-period variance components.

Appendix or supplementary material:

- Full Table 2.
- Full Table 3.
- Random-intercept forest plots for policy period and clade.
- Primary-model sensitivity tables.

Direction probabilities should mainly appear in tables and prose. They do not need to be plotted as a primary visual result. The primary plots should show estimates and 95% HDIs.

## 12. Draft Results Summary

The Bayesian SSE detection models indicate that detection signal is structured by both the composition and mixing of sampled populations. The main mixing models use null-standardised entropy, which compares each cluster with a within-window, cluster-size matched expectation. These models show that several attributes have detection associations that differ from their absolute observed-entropy associations, especially age, SIMD, sex, and health board. Observed entropy remains useful as a companion scale: higher observed age, SIMD, urban/rural, and health-board mixing were associated with higher candidate odds and higher burst scores, although burden-score effects sometimes moved in the opposite direction.

Composition models showed coherent demographic and geographic structure. Female samples had consistently lower candidate odds and lower burst and burden scores than male samples. Age 15-24 had higher detection signal than age 25-64, while children and older age groups generally had lower signal. Urban/rural composition was also consistent: samples outside large urban areas, especially remote rural and remote town categories, had lower candidate odds and lower linear scores.

Health-board effects were substantial. Candidate odds were particularly high in Orkney and Forth Valley relative to Greater Glasgow and Clyde, while several boards including Borders, Fife, Grampian, and Ayrshire and Arran had lower candidate odds. Linear outcomes showed additional heterogeneity, with some boards having higher burst scores but lower burden scores. This supports the interpretation that geographic structure is a major component of detection variation.

Surveillance context also mattered. Recent sequencing intensity was positively associated with candidate detection, especially in composition models, while cumulative incidence was generally negatively associated with candidate odds and linear outcomes. Cumulative sequencing proportion was weaker and less central.

Residual temporal and lineage structure remained after fixed-effect adjustment. Policy-period effects were clearest in the composition candidate model, with higher residual candidate odds in F5, P3, and T1, and lower odds in FE, L0, L2, L21, PR, and SL. Clade effects were especially pronounced, with several later clades showing much higher residual candidate odds and some earlier or recombinant clades showing lower odds. Random-effect SDs further support this: clade heterogeneity exceeded policy-period heterogeneity in composition models and in most mixing models.

Overall, the results support a narrative in which SSE detection is shaped jointly by demographic composition, socioeconomic and geographic structure, surveillance intensity, epidemic timing, and viral lineage. The strongest evidence comes from diagnostically clean mixing models and composition logistic models. Composition linear results are useful but should be interpreted with caution because several of those models have diagnostic warnings.
