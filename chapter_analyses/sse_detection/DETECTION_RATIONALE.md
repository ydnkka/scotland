# SSE Detection Rationale

The detector operates on window-specific clusters linked across adjacent retained windows when they share sequences. Such an edge records continuity between clustering solutions, not observed transmission.

## Detection versus characterisation

Candidate assignment uses only cluster magnitude, upstream sequence novelty, and downstream sequence burden. Age, sex, deprivation, geography, policy, and mixing entropy are attached later for description and regression. This separation prevents those characteristics from defining the candidate set they are intended to explain.

## Two detection axes

**Local burst** combines context-adjusted cluster size with the fraction of sequences absent from the union of direct-parent memberships. Parentless nodes use size alone because upstream novelty is unobservable, not 100%.

**Onward burden** combines source-size-normalised direct and cumulative reach. Direct burden allocates each successor's new sequences among its parents in proportion to incoming shared-sequence support. Cumulative burden counts unique sequences in all reachable descendants, excluding the source and deduplicating recurrence.

These axes distinguish a large local accumulation from a cluster followed by disproportionate downstream sequence mass. Branching entropy and successor concentration describe the shape of onward spread but do not define candidates.

Burden is applicable only when direct attributable or cumulative burden is positive. Non-propagating nodes remain not applicable rather than receiving a structural-zero rank.

## Calibration and candidates

Within-window percentile composites are calibrated by seeded profile permutation. Profiles are permuted intact to preserve component dependence and missingness. One-sided upper-tail p-values ask whether a score is unusually high.

Discrete scores create exact ties. The operational p-value randomises position within the tied null mass; a conservative version that counts every tie in the upper tail is retained for audit. Candidate status uses the operational p-value:

- high priority at `p <= 0.05` on burst, burden, or both;
- possible review at `p <= 0.10` on either applicable axis;
- size-ineligible below 6 sequences.

Calibration should be examined over all tested nodes and by window, size, component availability, and follow-up. Conservative p-values may pile up near one; false-positive inflation appears near zero.

## Limits

Signals depend on sampling, sequencing coverage, rolling-window construction, EpiLink sparsification, and Leiden resolution. Parentless and sink nodes may reflect observation boundaries. Cumulative burden is right-censored, and missing descendants do not prove containment.

The output is therefore a calibrated prioritisation for external epidemiological review. It cannot verify an SSE without independent evidence.

## Statistical references

- Habiger and Peña (2011), [randomised p-values and nonparametric multiple testing](https://doi.org/10.1080/10485252.2010.482154).
- Hemerik and Goeman (2018), [exact testing with random permutations](https://doi.org/10.1007/s11749-017-0571-1).
- Phipson and Smyth (2010), [non-zero permutation p-values](https://doi.org/10.2202/1544-6115.1585).
