# Chapter 1 — analysis documentation

## 1. Rationale

Genomic clusters with many sequences or many distinct datazones plausibly
reflect either real transmission expansion (e.g. superspreading) or
sampling-window contraction.  The question this chapter tackles is whether
those larger or more spatially dispersed clusters are also more
*sociodemographically mixed* than we would expect by chance — that is,
whether they are characterised by transmission across sociodemographic
boundaries beyond what marginal composition and cluster size already
imply.

The "beyond what would be expected" part is doing real work and is worth
unpacking before reading the regression tables.

## 2. What "expected" means

Each cluster sits inside an analysis window (three-week sliding window)
and a Pango lineage.  If cluster membership were independent of age,
sex, and area deprivation given the window × lineage stratum, then the
probability that any pair of cases in a cluster differ on a given
variable would equal the corresponding pair-discordance in the full
stratum.  We compute that stratum-level pair discordance, call it
*expected* discordance, and the per-cluster observed minus expected
difference is *excess discordance*.

Positive excess discordance means a cluster is more mixed than the
random-assembly null implies; negative excess discordance means it is
more homogeneous.

This null is conservative: it does *not* hold the cluster's own
marginal composition fixed.  A cluster that happens to draw a peculiar
slice of the stratum can register non-zero excess discordance even if
its members were exchangeable within their own composition.  The
null-regression sensitivity (§5) tightens the null by additionally
adjusting for the cluster's own marginal entropy.

## 3. Outcomes

`cluster_size` and `cluster_n_datazones`, both heavily right-skewed
counts with a structural mass at one.  Both are modelled as
**zero-truncated negative binomial (ZTNB) on the non-singleton sub-
population**.  Cluster size response is the excess count `raw − 1`. This is consistent with the observation that mixing is only defined for non-singleton clusters.

Inference uses analytical first-order scores in the ZTNB log-likelihood,
a numerical Hessian as the sandwich bread, and cluster-robust standard
errors clustered by `window_id`.  The cluster correction follows the
usual `G/(G-1) · (N-1)/(N-K)` finite-sample adjustment.

## 4. Predictors and adjustments

The three excess-mixing predictors are z-scored versions of
observed-minus-expected pair discordance for age band, sex, and SIMD
quintile.  They enter all main models simultaneously so that each
estimate is conditional on the others; this is important because
SIMD-quintile mixing and age-band mixing are weakly correlated.

Adjustment covariates, all z-scored:

- `deprivation_z` (mean cluster SIMD rank, sign flipped so higher means
  more deprived)
- `local_incidence_z` (log1p of mean local cumulative incidence)
- `local_seq_fraction_z` (logit of mean local sequencing fraction)
- `window_seq_fraction_z` (logit of window sequencing proportion)
- `test_positivity_z` (logit of mean local test positivity)

Calendar time is modelled with an 8-df B-spline on `window_idx`.
Lineage is a categorical with rare lineages (<50 clusters) pooled into
"Other rare lineages".  For the wave-interaction model lineage is
replaced by wave dummies because the two are collinear at the broad-
group level.

## 5. Model specifications

### 5.1 Main effects (`fit_main_effects`)

Outcome | Form
---|---
`cluster_size` | `log E[size − 1] = β₀ + β·(excess mixing) + γ·(adjustments) + spline(window_idx) + lineage`
`cluster_n_datazones` | `log datazones = β₀ + β·(excess mixing) + γ·(adjustments) + spline(window_idx) + lineage`

Both equations are fit on non-singleton clusters with ZTNB.

The size-adjusted geographic-spread variant adds `log_cluster_size_z`
to the RHS, which lets us read off whether mixing predicts spread
*beyond what size alone explains*.

### 5.2 Wave interactions (`fit_wave_interactions`)

Same as main effects but for each excess-mixing predictor we add eight
interaction columns (one per wave, with Delta as reference).  Lineage
dummies are dropped and replaced by wave dummies; this is intentional
because the interaction itself absorbs the broad-lineage variation we
care about.

Wave-specific slope estimates are recovered as the main effect plus
the wave interaction, with the corresponding standard error obtained
under the assumption that main and interaction terms are independent
(this is an approximation — the table reports the underlying
coefficients and CIs so the reader can reconstruct the joint
covariance from the diagnostics file if needed).

### 5.3 Sensitivities

**Size spline** (`fit_size_spline_sensitivity`): replaces the linear
`log_cluster_size_z` adjustment with a 4-df B-spline of `log(size)`.
Only the spread ZTNB is fit; the spline coefficients are reported so
the reader can read off the shape of the size-adjustment.

**SIMD decile mixing** (`fit_simd_decile_sensitivity`): refits the
main model with `simd_decile_excess_mixing_z` (observed minus expected
pairwise decile-discordance, z-scored) in place of the
quintile-resolution SIMD predictor.  Age and sex predictors are
unchanged.  The purpose is to check whether the SIMD-mixing →
cluster-scale finding is sensitive to where the quintile cutpoints
fall.  Decile resolution gives finer sensitivity to within-quintile
gradient but exposes the expected-discordance calculation to more
small-cell noise within window × lineage strata, particularly in the
sparser late-pandemic waves; the decile predictor is therefore
reported only as a sensitivity, not as part of the primary model.

**Null-residual mixing** (`build_null_residual_mixing` +
`fit_null_residual_sensitivity`): builds an alternative excess-mixing
predictor as the residual from a per-dimension null regression of
*observed* discordance on `log(size)`, lineage, the calendar spline,
and the cluster's own marginal entropy.  This null additionally
conditions on the cluster's own composition, so it is stricter than
the published observed-minus-expected formulation.  We refit the main
effects with these residual-based predictors in place of the
observed-minus-expected predictors.

**Joint-profile predictors** (`fit_profile_predictor`): refits the main
model using a single joint-profile excess-mixing predictor.  Two
variants are reported side-by-side:

- *Demographic profile* — age × sex.  Speaks to mixing across the
  combined demographic boundary independent of area deprivation.
- *Sociodemographic profile* — age × sex × SIMD quintile.  Speaks to
  mixing across the combined demographic *and* area-deprivation
  boundary.

The two are useful as unidimensional summaries that aggregate boundary
crossing across constituent dimensions.

### 5.4 Stratifications

**SIMD domain** (`fit_domain_main_effects`): refits the main model
seven times, each time replacing overall SIMD mixing with the
quintile-mixing of one SIMD domain (income, employment, education,
health, access, crime, housing).  Age and sex stay in.

**Epidemic wave** (`fit_wave_stratified`): refits the main model
within each wave, using a wave-restricted lineage list.  Waves with
fewer than 50 non-singleton clusters are skipped to avoid unstable
fits, with the omission recorded in the diagnostics file.

## 6. Interpreting the outputs

Coefficients on z-scored predictors are per-SD effects.  The table reports both the raw coefficient and
its exp() rate ratio plus 95% CI; the linear-component column is
populated only by the size-spline rows.

For the wave-interaction model, the rows of interest are the
mixing × wave interactions: a positive interaction coefficient on,
say, `age_excess_mixing_z × wave_BA.2` means that the age-mixing
slope is steeper in BA.2 than in the Delta reference.  The plot
`wave_interaction_slopes.{png,pdf}` reconstructs the implied wave-
specific slope and 95% CI for each predictor.

For the size-adjusted spread model, a positive coefficient on
`age_excess_mixing_z` (etc.) with the size adjustment in place means
the excess-mixing → spread association persists *after* removing the
mechanical correlation that bigger clusters cover more datazones.

For the null-residual sensitivity, agreement with the main-effects
table strengthens the published interpretation; disagreement
suggests the observed-minus-expected formulation is picking up
within-cluster composition effects that the stricter null absorbs.

## 7. Known limitations

1. The expected discordance is computed within window × lineage,
   not within window × lineage × cluster-size.  Very large clusters
   may therefore have systematically lower observed-minus-expected
   discordance simply because their large valid-pair denominator
   stabilises the estimate.  The size-adjusted spread model and the
   null-residual sensitivity are the two principled defences
   against this.
2. The wave-interaction model drops lineage adjustment in favour of
   wave dummies.  Within-wave lineage heterogeneity is therefore
   not adjusted for in that model; if this matters, look at the
   wave-stratified fits (`wave_analysis.py`), which keep lineage
   inside each wave.
3. Mixing is only defined for non-singleton clusters.  Singletons
   are dropped from every fit in this chapter.  Cluster-size
   inference therefore conditions on `size ≥ 2`.
