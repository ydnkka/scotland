# SSE Detection

Code for an alternate-window cluster-transition graph, a demographic-blind superspreading-compatible signal detector, composition summaries, and Bayesian characterisation.

Candidates are graph-derived priorities for review, not confirmed superspreading events or transmission chains.

The genomic-network package owns the underlying window-specific clusters and compatibility networks; this package owns their temporal continuity and downstream characterisation.

## Detection Rationale

The detector operates on window-specific clusters linked across adjacent retained windows when they share sequences. Such an edge records continuity between clustering solutions, not observed transmission.

### Detection versus characterisation

Candidate assignment uses only cluster magnitude, upstream sequence novelty, and downstream sequence burden. Age, sex, deprivation, geography, policy, and mixing entropy are attached later for description and regression. This separation prevents those characteristics from defining the candidate set they are intended to explain.

### Two detection axes

**Local burst** combines context-adjusted cluster size with the fraction of sequences absent from the union of direct-parent memberships. Parentless nodes use size alone because upstream novelty is unobservable, not 100%.

**Onward burden** combines source-size-normalised direct and cumulative reach. Direct burden allocates each successor's new sequences among its parents in proportion to incoming shared-sequence support. Cumulative burden counts unique sequences in all reachable descendants, excluding the source and deduplicating recurrence.

These axes distinguish a large local accumulation from a cluster followed by disproportionate downstream sequence mass. Branching entropy and successor concentration describe the shape of onward spread but do not define candidates.

Burden is applicable only when direct attributable or cumulative burden is positive. Non-propagating nodes remain not applicable rather than receiving a structural-zero rank.

### Calibration and candidates

Within-window percentile composites are calibrated by seeded profile permutation. Profiles are permuted intact to preserve component dependence and missingness. One-sided upper-tail p-values ask whether a score is unusually high.

Discrete scores create exact ties. The operational p-value randomises position within the tied null mass; a conservative version that counts every tie in the upper tail is retained for audit. Candidate status uses the operational p-value:

- high priority at `p <= 0.05` on burst, burden, or both;
- possible review at `p <= 0.10` on either applicable axis;
- size-ineligible below 6 sequences.

Calibration should be examined over all tested nodes and by window, size, component availability, and follow-up. Conservative p-values may pile up near one; false-positive inflation appears near zero.

### Limits

Signals depend on sampling, sequencing coverage, rolling-window construction, EpiLink sparsification, and Leiden resolution. Parentless and sink nodes may reflect observation boundaries. Cumulative burden is right-censored, and missing descendants do not prove containment.

The output is therefore a calibrated prioritisation for external epidemiological review. It cannot verify an SSE without independent evidence.

### Statistical references

- Habiger and Peña (2011), [randomised p-values and nonparametric multiple testing](https://doi.org/10.1080/10485252.2010.482154).
- Hemerik and Goeman (2018), [exact testing with random permutations](https://doi.org/10.1007/s11749-017-0571-1).
- Phipson and Smyth (2010), [non-zero permutation p-values](https://doi.org/10.2202/1544-6115.1585).

## Documentation

- [TECHNICAL.md](TECHNICAL.md): implemented graph, scores, calibration, tiers, and outputs.
- [BAYESIAN_MODELS.md](BAYESIAN_MODELS.md): model frames, formulas, commands, priors, and diagnostics.

## Commands

Run from the repository root:

```bash
# Rebuild detector outputs
python -m analyses.sse_detection.lib.sse.detection

# Build wide cluster-composition tables
python -m analyses.sse_detection.build_composition_tables

# Rebuild detector/Bayesian figures
python -m results.make_figures --domain sse_detection --skip-missing

# Rebuild Bayesian LaTeX table fragments and CSV/parquet companions
python -m results.make_tables --domain sse_detection --skip-missing
```

The detector has no CLI options; its settings are constants in `lib/sse/config.py`. Bayesian fitting commands are listed in [BAYESIAN_MODELS.md](BAYESIAN_MODELS.md).

## Chapter 6 descriptive size analyses

These builders reuse saved detector scores, labels, and entropy values. They do
not fit Bayesian models or recalibrate the detector. They compare cluster and
record weighting, the fixed size bands 6–9, 10–19, 20–49, 50–99, and 100+,
screening routes, saved-p cutoffs 0.01/0.025/0.05/0.10, and restrictions to
primary-labelled clusters of size at least 10 or 20.

```bash
MPLBACKEND=Agg python -m results.make_figures sse_detection:ch6_size_profiles
MPLBACKEND=Agg python -m results.make_tables sse_detection:ch6_size_profiles
python -m unittest analyses.sse_detection.test_characterisation -q
```

`lib/characterisation.py` provides aggregate calculations and reconciliation
checks. `lib/figs/size_profiles.py` writes eight PDF/PNG figures, four LaTeX
tables, and five aggregate CSV/parquet tables named `tab_ch6_*`. These source
tables include category percentages, differences under both weighting schemes,
entropy quartiles, and missing-value denominators. No identifiers are exported.
Burden-only and both-axis comparisons use burden-eligible background clusters.
Intervals in the entropy figures are IQRs, not confidence intervals. These are
unadjusted descriptions; small differences and size overlap do not validate
the screen against confirmed superspreading events.

## Layout

- `lib/sse/`: data loading, transition graph, composition entropy, feature assembly, scoring, and diagnostics.
- `lib/model/`: Bayesian preparation, fitting, and CLI orchestration.
- `lib/figs/`: individual figure builders and publication table builders.
- `results/sse_outputs/`: authoritative detector cluster/edge tables and graph summaries.
- `results/tables/`: composition/source tables and generated CSV/parquet companions.
- `../../results/figures/`: project-level detector and Bayesian figures.
- `../../results/tables/`: project-level publication `.tex` table fragments.
- `results/bayesian_outputs/`: fitted model summaries, diagnostics, metadata, logs, and optional inference data.
