# Chapter 1 Figure And Table Plan

## Narrative Spine

The manuscript question is:

> Among non-singleton SARS-CoV-2 genomic clusters, are clusters with greater
> excess sociodemographic mixing larger and more geographically dispersed,
> after adjustment for lineage, calendar time, surveillance intensity,
> epidemic intensity, and cluster deprivation?

The figure and table package therefore does four jobs:

1. Define the analysis population (non-singleton clusters, n = 84 067) and the
   observed-minus-expected mixing construct.
2. Present the main pooled association between excess mixing and cluster
   scale.
3. Show whether the association varies by epidemic wave and SIMD domain.
4. Demonstrate that the main pattern survives robustness checks on sliding
   windows, extreme right-tail clusters, predictor-set choice, mixing-null
   definition, and standard-error clustering unit.

## Conventions Inherited From Part 1

These conventions are global to every figure unless an individual spec
overrides them, and are encoded as constants in `make_figures.py`.

- **Plotting helper.** All figures use `utils.style.set_theme`,
  `style.new_figure`, `style.save_figure`, and `style.add_panel_labels`.
  `font_scale=0.85`, `context="paper"`.
- **Save formats.** Every figure is written as PDF, PNG, and TIFF at 600 dpi
  to `chapter1/manuscript/figures/`.
- **Population.** All figures except the descriptive distributions are
  estimated on non-singleton clusters (`cluster_size > 1`); the model files
  carry `component == "positive_count"` only. There is no hurdle component
  for Chapter 1 (singletons are excluded before fitting), so figures and
  captions do not refer to hurdle odds.
- **Wave order.** `["B.1.177", "Alpha", "Delta", "BA.1", "BA.2", "BA.4",
  "BA.5", "BQ.1"]`. XBB and "Other" are not included in chapter 1
  wave-stratified results.
- **SIMD domain order.** `["overall", "income", "employment", "education",
  "health", "access", "crime", "housing"]`, palette `style.SIMD_DOMAIN_PALETTE`.
- **Mixing predictors.** Three primary terms with stable colour mapping:
  - `simd_excess_mixing_z` -> "SIMD excess mixing" -> `palette["crime"]` (yellow)
  - `age_excess_mixing_z` -> "Age excess mixing" -> `palette["housing"]` (brown)
  - `sex_excess_mixing_z` -> "Sex excess mixing" -> `palette["overall"]` (dark)
- **Adjustment terms.** `deprivation_z`, `local_incidence_z`,
  `local_seq_fraction_z`, `window_seq_fraction_z`, `test_positivity_z`,
  `log_cluster_size_z` (size-adjusted spread only).
- **Forest panels.** Log-x axis with `axvline(1.0)` at unity, ticks drawn
  from `[0.5, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0]`, minor ticks
  removed (`NullLocator`/`NullFormatter`), point with CI bar coloured by
  term.
- **Heatmap panels.** `TwoSlopeNorm(vcenter=1.0, vmin, vmax)` for rate
  ratios, `cmap="RdBu_r"`, white gridlines between cells, side colorbar with
  explicit unit label, optional star or hatched outline for cells whose 95%
  CI crosses 1.
- **Panel labels.** `A`, `B`, `C`, ... placed via
  `style.add_panel_labels(axes.ravel(), x=-0.16, y=1.14, size=9)`.
- **SE reporting.** The primary inference uses window-clustered sandwich
  SEs; this is reported in figure captions but not in row labels (every row
  uses the same primary SE except in Figure 4).

## Main Figures

### Figure 1. Non-singleton Cluster Outcomes And Mixing Predictors

**Purpose.** Compact descriptive view of the analysis population: three
cluster-scale outcomes and three primary observed-minus-expected mixing
predictors among non-singleton clusters.

**Layout.** Six-panel 2 x 3 binned-percent figure (`width="double"`,
`height_in=4.9`, `font_scale=0.85`).

- Top row, common grey bars, x-axis binned:
  - Panel A: `cluster_size`, bins `[2, 3, 4-5, 6-10, 11-20, 21-50, >50]`.
  - Panel B: `duration_days`, bins `[0, 1, 2, 3-5, 6-10, 11-15, 15+]`.
  - Panel C: `geographic_spread`, bins `[1, 2, 3-4, 5-9, 10-19, 20-49, >50]`.
- Bottom row, common grey bars, 10 pp histograms, dashed zero line:
  - Panel D: `age_excess_discordance`, x in [-100, 100] pp.
  - Panel E: `sex_excess_discordance`, x in [-100, 100] pp.
  - Panel F: `simd_excess_discordance`, x in [-100, 100] pp.
- Y-axis is "Clusters (%)" on the leftmost column only; share x within each
  row but not across rows.

**Primary inputs.**

- `chapter1/cache/cluster_table.parquet`

**Caption message.** Distributions are restricted to non-singleton clusters
(n = 84 067). Excess mixing is the observed minus lineage-window expected
pair discordance, shown in percentage points before model standardisation.

### Figure 2. Main Pooled Effects Of Excess Mixing On Cluster Scale

**Purpose.** Central result figure. Answers the headline question without
sending the reader to the supplement.

**Layout.** Three-panel forest plot, single row, `width="double"`,
`height_in=3.4`, `font_scale=0.85`. Each panel shares the same y-axis (the
three mixing predictors) and uses a panel-specific log-x range:

- Panel A: cluster size ZTNB count ratio per 1 SD higher predictor
  (`main_effects_results.csv`, `outcome == "cluster_size"`).
- Panel B: geographic spread ZTNB count ratio per 1 SD higher predictor
  (`main_effects_results.csv`, `outcome == "geographic_spread"`).
- Panel C: size-adjusted geographic spread ZTNB count ratio per 1 SD higher
  predictor, adjusted for `log_cluster_size_z`
  (`size_adjusted_spread_results.csv`).

Rows: `simd_excess_mixing_z`, `age_excess_mixing_z`, `sex_excess_mixing_z`
(SIMD first, then age, then sex, matching the head order of part 1).

**Primary inputs.**

- `chapter1/tables/main_effects_results.csv`
- `chapter1/tables/size_adjusted_spread_results.csv`
- `chapter1/tables/main_effects_diagnostics.csv` (for n in caption)
- `chapter1/tables/size_adjusted_spread_diagnostics.csv` (for n in caption)

**Caption message.** Size-adjusted spread is critical because it separates
wider geography from the mechanical tendency for larger clusters to include
more datazones. The adjustment set (deprivation, surveillance intensity,
test positivity, lineage, calendar time) is reported in Supp Table 1; for
clarity Figure 2 plots only the three mixing predictors.

### Figure 3. Epidemic-wave Heterogeneity

**Purpose.** Show whether the pooled mixing-scale association is stable
across variant/wave context.

**Layout.** Two-panel heatmap, `width="double"`, `height_in=4.0`,
`font_scale=0.80`.

- Panel A: cluster size ZTNB count ratio.
- Panel B: geographic spread ZTNB count ratio.
- Rows in each panel: SIMD, age, sex excess mixing (left labels on Panel A
  only).
- Columns: waves in canonical order
  `["B.1.177", "Alpha", "Delta", "BA.1", "BA.2", "BA.4", "BA.5", "BQ.1"]`.
- Cell values: rate ratios per 1 SD higher predictor.
- Colour scale: `TwoSlopeNorm(vcenter=1.0, vmin=min(0.0, observed_min),
  vmax=max(2.0, observed_max))`, `cmap="RdBu_r"`. The colourbar runs along
  the right side and is labelled "ZTNB count ratio per 1 SD higher excess
  mixing".
- Cell annotation: numeric ratio printed in white when |log(ratio)| > 0.3
  and dark grey otherwise. Cells whose 95% CI crosses 1 receive a small
  white open-circle marker so non-significant cells are visually flagged.

**Primary inputs.**

- `chapter1/tables/wave_stratified_results.csv` (primary)
- `chapter1/tables/wave_interaction_results.csv` (companion supplementary
  table, see Supp Table 3)
- `chapter1/tables/wave_stratified_diagnostics.csv` (n by wave)

**Caption message.** The figure should not imply that waves are independent
causal settings. It should show whether the main association concentrates
in waves with deeper sampling and larger reconstructed clusters, and
whether late waves attenuate or destabilise the estimate.

### Figure 4. Robustness Dashboard

**Purpose.** Put the key sensitivity checks in one high-signal figure so
the reader can see at a glance whether the main result is an artefact of
one modelling choice.

**Layout.** Three-column dot-range robustness plot, `width="double"`,
`height_in=6.0`, `font_scale=0.80`.

- Columns:
  - cluster size
  - geographic spread
  - size-adjusted geographic spread
- Rows grouped by sensitivity, each sensitivity contributes three
  point-with-CI rows (SIMD, age, sex), banded by mixing predictor colour:
  - Primary, window-clustered SE (reference row)
  - Health-board clustered SE
  - Non-overlapping windows (`--window-stride 3`)
  - 99% winsorised outcome
  - Top 0.5% outcome excluded
  - Finite-sample standardised mixing
  - Joint-profile adjusted predictor set
- Log-x axis, vertical reference line at 1.0, x-axis tick set from
  `[0.5, 0.7, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0]`, common across the three
  outcomes.
- A thin horizontal divider between sensitivity blocks.

**Primary inputs.**

- `chapter1/tables/main_effects_results.csv`
- `chapter1/tables/size_adjusted_spread_results.csv`
- `chapter1/tables/joint_profile_adjusted_results.csv`
- `chapter1/tables/finite_sample_mixing_sensitivity_results.csv`
- `chapter1/sensitivity/tables_health_board/main_effects_results.csv`
- `chapter1/sensitivity/tables_health_board/size_adjusted_spread_results.csv`
- `chapter1/sensitivity/tables_stride3/main_effects_results.csv`
- `chapter1/sensitivity/tables_stride3/size_adjusted_spread_results.csv`
- `chapter1/sensitivity/tables_winsorise99/main_effects_results.csv`
- `chapter1/sensitivity/tables_winsorise99/size_adjusted_spread_results.csv`
- `chapter1/sensitivity/tables_exclude_tail995/main_effects_results.csv`
- `chapter1/sensitivity/tables_exclude_tail995/size_adjusted_spread_results.csv`

**Caption message.** Health-board clustering is an inference sensitivity,
not a new substantive model. It tests whether uncertainty conclusions are
robust to a coarse spatial clustering unit. The dashboard panels should be
read column by column; the row stack within each panel is purely visual
grouping.

## Main Tables

### Table 1. Analysis Population And Outcome Descriptives

**Purpose.** Give the reader the denominator and the shape of the data
before model estimates appear.

**Recommended columns.**

- sequence rows used
- clusters (all resolutions of the primary cut)
- non-singleton clusters
- windows
- raw Pango lineages
- modelled lineage levels after pooling (with pooling threshold; 30
  non-singleton clusters per lineage)
- for cluster size and distinct datazones: median, IQR, 90th percentile,
  99th percentile, maximum
- fraction of all clusters that are singletons
- fraction of all clusters with one datazone
- fraction of non-singleton clusters with one datazone
- per-wave non-singleton counts (B.1.177 ... BQ.1)

**Primary inputs.**

- `chapter1/tables/dataset_descriptives.csv`
- `chapter1/cache/cluster_table.parquet`

### Table 2. Main Pooled ZTNB Model Estimates

**Purpose.** Numerical companion to Figure 2.

**Recommended columns.**

- outcome
- model (`main` or `main_size_adjusted`)
- predictor
- estimate (log scale)
- robust SE (window-clustered)
- rate ratio
- 95% CI
- p-value
- n observations
- alpha (negative-binomial dispersion)
- convergence status

**Primary inputs.**

- `chapter1/tables/main_effects_results.csv`
- `chapter1/tables/size_adjusted_spread_results.csv`
- `chapter1/tables/main_effects_diagnostics.csv`
- `chapter1/tables/size_adjusted_spread_diagnostics.csv`

**Notes.** Include the adjustment covariates in this table or in a clearly
linked supplementary table. The main printed table can focus on age, sex,
and SIMD mixing if journal space is tight.

## Supplementary Figures

### Supplementary Figure 1. Full Outcome Distributions

**Purpose.** Tail-focused companion to Figure 1 panels A-C with more
granularity than the compact main-text binned figure.

**Panels.** Five-panel layout, `width="double"`, `height_in=5.0`:

- A: all-cluster `cluster_size` histogram on log-y, singleton bar shaded a
  contrasting grey.
- B: non-singleton `cluster_size` histogram, log-y.
- C: all-cluster `geographic_spread` distribution starting at 1, with the
  one-datazone bar shaded.
- D: non-singleton `geographic_spread` distribution starting at 1.
- E: hexbin of `log10(cluster_size)` versus `log10(geographic_spread)` for
  non-singletons, with a 1:1 diagonal for reference.

**Inputs.**

- `chapter1/cache/cluster_table.parquet`

### Supplementary Figure 2. Excess-mixing Predictor Distributions

**Purpose.** Extend Figure 1 panels D-F to the full predictor set used in
the main and joint-profile models, and show the relationship between raw
percentage-point excess mixing and the z-scored model scale.

**Panels.** Six-panel 2 x 3 layout, `width="double"`, `height_in=4.6`:

- A: `age_excess_discordance` (pp).
- B: `sex_excess_discordance` (pp).
- C: `simd_excess_discordance` (pp).
- D: `simd_decile_excess_discordance` (pp).
- E: `demographic_profile_excess_discordance` (pp).
- F: `socio_demographic_profile_excess_discordance` (pp).

Inset on each panel: scatter of raw pp versus the corresponding `_z`
column, with the model-scale standardisation slope drawn as a line.

**Inputs.**

- `chapter1/cache/cluster_table.parquet`
- `chapter1/tables/covariate_scaling.csv`

### Supplementary Figure 3. Excess Mixing Across Cluster-size Bins

**Purpose.** Show whether the three mixing predictors shift monotonically
or non-linearly across the empirical cluster-size distribution.

**Panels.** 2 x 3 grid sharing the cluster-size bins of Figure 1:

- Top row, boxplots of excess mixing in percentage points by cluster-size
  bin, median, IQR, and whiskers, outlying points omitted for readability,
  dashed reference line at zero.
- Bottom row, 100% stacked bars showing the proportion of clusters in each
  bin with negative excess, baseline excess (within +/- 0.5 pp of zero), or
  positive excess.

**Inputs.**

- `chapter1/cache/cluster_table.parquet`

### Supplementary Figure 4. Observed-versus-expected Pair Matrices

**Purpose.** Validate the mixing construction by showing where excess
pairing occurs.

**Panels.** Three heatmaps, `width="double"`, `height_in=3.8`:

- A: SIMD quintile, observed-minus-expected pair probability in pp.
- B: SIMD decile, observed-minus-expected pair probability in pp.
- C: Age band, observed-minus-expected pair probability in pp.

Uses `wave_group == "Overall"` rows of
`observed_expected_mixing_matrices.csv`; categorical order is "1, 2, ..."
for SIMD, and "0-9, 10-19, ..., 75+" for age, sorted with the leading
integer key already used in part 1.

Colour scale: `TwoSlopeNorm(vcenter=0, vmin, vmax)` shared across panels,
`cmap="RdBu_r"`. Colourbar label: "Observed - expected pair probability
(pp)".

**Inputs.**

- `chapter1/tables/observed_expected_mixing_matrices.csv`

### Supplementary Figure 5. Size Adjustment Diagnostics

**Purpose.** Show that the size-adjusted spread inference is not dependent
on the specific size-adjustment functional form.

**Panels.** Two-panel forest plot, `width="onehalf"`, `height_in=3.4`:

- A: pooled main spread vs linear log-size-adjusted spread vs
  spline-log-size-adjusted spread, plotted as three offset rows per
  mixing predictor.
- B: estimated spline contribution of `log_cluster_size` to log-rate, with
  partial residuals (only if the spline term coefficients are recoverable
  from the `size_spline_sensitivity` rows; otherwise this panel is
  dropped).

**Inputs.**

- `chapter1/tables/main_effects_results.csv`
- `chapter1/tables/size_adjusted_spread_results.csv`
- `chapter1/tables/size_spline_sensitivity_results.csv`

### Supplementary Figure 6. SIMD-domain Mixing Predictors

**Purpose.** Show whether the SIMD signal is carried by particular SIMD
domains.

**Layout.** Heatmap with eight domain rows
(`["overall", "income", "employment", "education", "health", "access",
"crime", "housing"]`) and three outcome columns:

- A: cluster size, ZTNB count ratio for the domain-specific
  `_domain_excess_mixing_z` predictor (the row's own "domain quintile").
- B: geographic spread.
- C: size-adjusted geographic spread.

Each domain row is a single estimate per outcome; age and sex mixing are
not repeated here because they are domain-invariant in chapter 1 (see
Figure 2). Colour scale rules match Figure 3.

**Inputs.**

- `chapter1/tables/domain_main_effects_results.csv`
- `chapter1/tables/domain_main_effects_diagnostics.csv`

### Supplementary Figure 7. Profile-predictor Sensitivities

**Purpose.** Keep the high-dimensional profile results available without
making them the main story.

**Panels.** Three-panel forest plot:

- A: demographic profile as single predictor
  (`profile == "demographic_profile"`).
- B: socio-demographic profile as single predictor
  (`profile == "socio_demographic_profile"`).
- C: age/sex/SIMD plus socio-demographic profile in the same model,
  showing all four predictors.

**Inputs.**

- `chapter1/tables/profile_predictor_results.csv`
- `chapter1/tables/joint_profile_adjusted_results.csv`

### Supplementary Figure 8. Model Diagnostics

**Purpose.** Make convergence and overdispersion transparent.

**Panels.** Four-panel grid, `width="double"`, `height_in=5.0`:

- A: alpha by model and outcome (dot-and-bar), main vs sensitivity sources.
- B: log-likelihood by model.
- C: n observations by outcome / model.
- D: tail caps and number excluded for tail-sensitivity runs.

**Inputs.**

- All `*_diagnostics.csv` files in `chapter1/tables/`.
- Diagnostics files under `chapter1/sensitivity/tables_*/`.

## Supplementary Tables

### Supplementary Table 1. Full Main-effects Model Output

**Contents.** All terms from the main pooled models, including adjustment
covariates and lineage/calendar-adjusted estimates.

**Inputs.**

- `chapter1/tables/main_effects_results.csv`
- `chapter1/tables/size_adjusted_spread_results.csv`

### Supplementary Table 2. Wave-stratified Estimates

**Contents.** Wave-specific age, sex, and SIMD mixing effects on cluster
size and geographic spread, ordered by canonical wave order.

**Inputs.**

- `chapter1/tables/wave_stratified_results.csv`
- `chapter1/tables/wave_stratified_diagnostics.csv`

### Supplementary Table 3. Wave-interaction Model

**Contents.** Main and interaction coefficients from the pooled wave
interaction model.

**Inputs.**

- `chapter1/tables/wave_interaction_results.csv`
- `chapter1/tables/wave_interaction_diagnostics.csv`

### Supplementary Table 4. SIMD-decile Sensitivity

**Contents.** Main model rerun with SIMD-decile mixing in place of
SIMD-quintile mixing.

**Inputs.**

- `chapter1/tables/simd_decile_sensitivity_results.csv`
- `chapter1/tables/simd_decile_sensitivity_diagnostics.csv`

### Supplementary Table 5. Finite-sample Mixing Sensitivity

**Contents.** Estimates using finite-sample standardised age, sex, and
SIMD mixing predictors.

**Inputs.**

- `chapter1/tables/finite_sample_mixing_sensitivity_results.csv`
- `chapter1/tables/finite_sample_mixing_sensitivity_diagnostics.csv`

### Supplementary Table 6. Joint-profile Sensitivities

**Contents.** Three rows of models:

- demographic profile as a single predictor
- socio-demographic profile as a single predictor
- age/sex/SIMD plus socio-demographic profile in the same model

**Inputs.**

- `chapter1/tables/profile_predictor_results.csv`
- `chapter1/tables/profile_predictor_diagnostics.csv`
- `chapter1/tables/joint_profile_adjusted_results.csv`
- `chapter1/tables/joint_profile_adjusted_diagnostics.csv`

### Supplementary Table 7. SIMD-domain Estimates

**Contents.** Domain-specific mixing predictor estimates and diagnostics.

**Inputs.**

- `chapter1/tables/domain_main_effects_results.csv`
- `chapter1/tables/domain_main_effects_diagnostics.csv`

### Supplementary Table 8. Robustness And Inference Sensitivities

**Contents.** One harmonised table for the robustness dashboard.

**Recommended sensitivity rows.**

- primary window-clustered SE
- health-board clustered SE
- non-overlapping windows
- 99% winsorised outcome
- top 0.5% outcome excluded
- finite-sample standardised mixing
- joint-profile adjusted predictor set
- null-residual mixing predictor

**Recommended columns.**

- sensitivity label
- outcome
- predictor
- estimate
- rate ratio
- 95% CI
- p-value
- n observations
- alpha
- SE clustering unit
- tail cap or number excluded, where relevant
- notes

**Inputs.**

- primary and size-adjusted results from `chapter1/tables/`
- `chapter1/sensitivity/tables_health_board/`
- `chapter1/sensitivity/tables_stride3/`
- `chapter1/sensitivity/tables_winsorise99/`
- `chapter1/sensitivity/tables_exclude_tail995/`

### Supplementary Table 9. ZTNB Diagnostics Across All Fits

**Contents.** A model-audit table collecting all diagnostics.

**Recommended columns.**

- table source
- model label
- outcome
- component (always `positive_count` for chapter 1)
- n observations
- converged
- log likelihood
- AIC
- alpha
- alpha at or near upper bound, if added
- winsorisation cap or tail-exclusion quantile, where relevant
- warnings or optimizer message, if available

**Inputs.**

- all diagnostics CSVs in `chapter1/tables/`
- all diagnostics CSVs in `chapter1/sensitivity/tables_*`

## Suggested File Naming

Main manuscript figures:

- `fig1_population_measures`
- `fig2_main_pooled_effects`
- `fig3_wave_heterogeneity`
- `fig4_robustness_dashboard`

Supplementary figures:

- `supp_fig1_outcome_distributions`
- `supp_fig2_mixing_distributions`
- `supp_fig3_cluster_size_mixing_boxplots`
- `supp_fig4_observed_expected_matrices`
- `supp_fig5_size_adjustment`
- `supp_fig6_domain_mixing_predictors`
- `supp_fig7_profile_predictors`
- `supp_fig8_model_diagnostics`

Tables:

- `table1_population_descriptives`
- `table2_main_model_estimates`
- `supp_table1_full_main_effects`
- `supp_table2_wave_stratified`
- `supp_table3_wave_interactions`
- `supp_table4_simd_decile`
- `supp_table5_finite_sample_mixing`
- `supp_table6_joint_profile`
- `supp_table7_domain_estimates`
- `supp_table8_robustness`
- `supp_table9_ztnb_diagnostics`

## Figure Priority If Space Is Limited

If the manuscript can only carry three main figures:

1. Keep Figure 1 as the compact six-panel descriptive distribution figure.
2. Keep Figure 2 as the central pooled estimate.
3. Combine wave heterogeneity and robustness into a two-row Figure 3.

If only two main figures are possible:

1. Figure 1: outcome and mixing-predictor distributions.
2. Figure 2: pooled effects plus a small robustness inset.

Wave, domain, profile, and diagnostics material can then move fully into
the supplement.

## Changelog: Improvements From The Part 1 Pass

These changes update the plan to reflect lessons from part 1 and the
narrower chapter 1 scope.

- Added a "Conventions Inherited From Part 1" block so figure scripts and
  readers share the same baseline (wave order, palette, save format, panel
  labels, log-x ticks, share-axis rules, no-hurdle clarification).
- Switched all main-figure language from "hurdle/ZTNB" to "ZTNB only" since
  chapter 1 fits the positive component only on non-singleton clusters.
  Updated component values to the actual `positive_count`.
- Pinned the three mixing predictors and their colour assignments globally
  so Figures 2-4 share visual encoding.
- Specified the wave list (eight waves, no XBB) and stable wave ordering
  for Figure 3 and Supp Table 2 to match chapter 1's
  `wave_stratified_results.csv`.
- Made Figure 3 a heatmap with a non-significance marker, matching the
  successful part 1 fig4 idiom while flagging wide CIs.
- Tightened Figure 4 to three outcome columns and seven sensitivity row
  blocks, all sourced from existing primary plus `tables_*` directories.
- Expanded Figure 1 caption and Supp Fig 1 to be explicit about including
  singletons (Supp Fig 1) versus excluding them (Figure 1).
- Specified TwoSlopeNorm centring and colourbar wording for every heatmap
  figure.
- Replaced ambiguous "inputs" lists with concrete CSV paths and column
  names that match the chapter 1 schema verified against the existing
  tables.
