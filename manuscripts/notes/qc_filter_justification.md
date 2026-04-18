# Methods note — justification for the Nextclade QC `good`-only filter

All three manuscripts restrict their analysis dataset to sequences with `nextclade_qc == "good"`. This note documents the empirical basis for that choice so it can be cited in a single line from each paper's Methods section.

## The check

At the primary clustering resolution (Leiden = 0.3) we tabulated how often sequences of each QC status were assigned to a singleton cluster versus a multi-member cluster, and fit a logistic GLM of `is_singleton` on QC status (reference: `good`).

| QC status | Sequence-window rows | P(singleton) | 95% Wilson CI | OR vs. `good` | 95% CI |
|---|---:|---:|---:|---:|---:|
| good | 801,249 | 12.71% | 12.64 – 12.79% | 1.00 (ref.) | — |
| mediocre | 62,894 | 17.68% | 17.38 – 17.98% | 1.47 | 1.44 – 1.51 |
| bad | 60,286 | 23.16% | 22.82 – 23.50% | 2.07 | 2.03 – 2.11 |

A `bad`-QC sequence is roughly twice as likely to be assigned to a singleton cluster as a `good`-QC sequence, and `mediocre` sequences sit at roughly 1.5× the baseline singleton rate. Both log-scale Z-tests give p far below machine precision.

## Why this is the expected direction

The Leiden clustering operates on pairwise single-nucleotide-variant distance. Low-QC genomes carry more ambiguous bases, spurious calls, and missing sites, so they share fewer near-identical variants with anyone else and are systematically pushed to the edges of the sequence-similarity graph. Emerging as a singleton — rather than joining a transmission chain — is therefore the null expectation for a noisy genome, not a biological signal.

## Implication for the analysis

Although the ~13% of sequences classified as `mediocre` or `bad` are a minority, they contribute a disproportionate share (19.8%) of the singleton pool. The inflation is also non-random with respect to the covariates of primary interest: QC is worse earlier in the pandemic when wet-lab protocols were less mature, and QC is plausibly correlated with both sampling intensity (`wn_prop_sequenced`) and regional sequencing capacity, which itself varies by deprivation. Including low-QC sequences therefore risks attributing a QC-compositional artifact to deprivation, demographic, or VOC covariates in the singleton and cluster-size models. Restricting to `good`-QC sequences removes the artifact at a modest cost in sample size (~13% of rows dropped) and keeps the singleton / multi-member split interpretable as a statement about shared genetic signal rather than about data quality.

## Pre-submission sensitivity analysis

Before submission we will refit the `singleton_model` (Paper 2) and Fig 4's quintile ORs (Paper 1) on the unfiltered 924k-row dataset with QC as an additional adjustment covariate (`mediocre` and `bad` dummies vs. `good`). If the deprivation and demographic ORs are stable after QC adjustment, the headline results are not driven by QC composition. If they shift, we report both the QC-filtered and QC-adjusted estimates and discuss.

## Reproducibility — code used to generate the numbers

Save as `manuscripts/notes/qc_filter_check.py` and run from the repository root with the primary analysis parquet in place. All three tables in this note come from a single script; total runtime is a few seconds on a laptop.

```python
"""Reproduce the Nextclade-QC singleton-rate check reported in
manuscripts/notes/qc_filter_justification.md."""
from __future__ import annotations

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm

from manuscripts.common import data

# 1. Load all QC statuses at the primary Leiden resolution (the default
#    loader restricts to qc='good'; pass qc=None to turn that off).
df = data.load_analysis_columns(
    ["sequence_id", "window_id", "cluster_id", "nextclade_qc"],
    resolution=data.PRIMARY_RESOLUTION,
    qc=None,
)

# 2. Label each (window, cluster) as singleton or multi-member, then
#    propagate the flag back to every sequence-window row.
cl = (
    df.groupby(["window_id", "cluster_id"])["sequence_id"]
    .nunique()
    .rename("n_seq")
)
df = df.merge(cl, left_on=["window_id", "cluster_id"], right_index=True)
df["is_singleton"] = (df["n_seq"] == 1).astype(float)


# 3. Descriptive: p(singleton) by QC status with Wilson 95% CIs.
def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float, float]:
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return p, max(0.0, centre - half), min(1.0, centre + half)


print(f"{'QC':<10}{'rows':>10}{'p(singleton)':>16}{'95% Wilson CI':>22}")
for q in ["good", "mediocre", "bad"]:
    sub = df[df["nextclade_qc"] == q]
    k, n = int(sub["is_singleton"].sum()), len(sub)
    p, lo, hi = wilson(k, n)
    print(f"  {q:<8}{n:>10,}{p*100:>14.2f}%{lo*100:>9.2f}%-{hi*100:.2f}%")

# 4. Inferential: logistic GLM of is_singleton ~ QC, 'good' as reference.
d = df.dropna(subset=["nextclade_qc"]).copy()
d["qc_med"] = (d["nextclade_qc"] == "mediocre").astype(float)
d["qc_bad"] = (d["nextclade_qc"] == "bad").astype(float)
X = sm.add_constant(d[["qc_med", "qc_bad"]], has_constant="add")
fit = sm.GLM(d["is_singleton"], X, family=sm.families.Binomial()).fit()

z95 = 1.959964
print(f"\n{'term':<10}{'OR':>8}{'95% CI':>20}{'p':>12}")
for term in ["qc_med", "qc_bad"]:
    est, se = fit.params[term], fit.bse[term]
    OR = float(np.exp(est))
    lo = float(np.exp(est - z95 * se))
    hi = float(np.exp(est + z95 * se))
    p = 2 * (1 - norm.cdf(abs(est / se)))
    print(f"  {term:<8}{OR:>8.3f}  {lo:>7.3f}-{hi:.3f}{p:>12.2e}")
```

Expected output with the dataset used in the three manuscripts:

```
QC              rows    p(singleton)         95% Wilson CI
  good        801,249         12.71%     12.64%-12.79%
  mediocre     62,894         17.68%     17.38%-17.98%
  bad          60,286         23.16%     22.82%-23.50%

term             OR             95% CI           p
  qc_med      1.474   1.443-1.506     0.00e+00
  qc_bad      2.069   2.028-2.111     0.00e+00
```
