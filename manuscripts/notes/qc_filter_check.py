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
