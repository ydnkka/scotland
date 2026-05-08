"""
Part 4 – Stage 2: Mutation trajectories, growth rate estimation,
and counterfactual modelling.

Growth rates are estimated using a binomial GLM (statsmodels) fitted
to sequenced S:N501Y counts, with additional window-level weighting
by confirmed positive-test volume.  Sensitivity analyses compare the
primary positive-test weighted model with an unweighted binomial GLM
and a coverage-adjusted GLM including proportion sequenced.

B.1.177 decline under L2 is kept as an exponential OLS fit
(coverage less problematic in the higher-sequencing L2 period).
"""
import sys
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import linregress
import statsmodels.api as sm
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parent.parent
OUT  = Path(__file__).resolve().parent / "tables"
NEXTCLADE_TSV = REPO / "data/raw/cog_all_scotland_nextclade.tsv"
DATA_PARQUET  = REPO / "data/processed/scotland_clustering_analysis_dataset.parquet"

# ── Load per-window positive-test counts ────────────────────────────────
print("Loading wn_positive_tests …")
win_tests = pq.read_table(str(DATA_PARQUET),
                           columns=["window_id", "wn_positive_tests"]).to_pandas()
win_tests = win_tests.drop_duplicates("window_id")[["window_id", "wn_positive_tests"]]
print(f"  {len(win_tests)} windows with test data")

# ── Load seq→window map ──────────────────────────────────────────────────
print("Loading seq→window map …")
seq_map = pd.read_parquet(OUT / "seq_window_map.parquet")
seq_map["wn_mid_date"] = pd.to_datetime(seq_map["wn_mid_date"])
all_seqids = set(seq_map["sequence_id"])
print(f"  {len(all_seqids):,} sequences, {seq_map['window_id'].nunique()} windows")

wpm = pd.read_csv(OUT / "part4_window_period_map.csv", parse_dates=["wn_mid_date"])
wpm_map = wpm.set_index("window_id")[["window_idx","wn_mid_date","policy_period","policy_intensity"]].to_dict("index")

# ── Key mutations to track ───────────────────────────────────────────────
# Alpha-defining AA mutations + B.1.177 signature for displacement context
TRACK_AA = ["S:N501Y","S:P681H","S:A222V","S:D614G","S:A570D","S:D1118H","N:R203K","N:G204R"]
# Key nuc substitutions for Alpha fingerprint
TRACK_NUC = ["A23063T","C23604A","G28881A","C14408T","C3037T","C22227T","A23403G"]

MIN_WINDOW_SEQS = 5

# ── Stream Nextclade TSV: per-window mutation frequencies ────────────────
print("Streaming Nextclade TSV …")
NC_COLS = ["seqName","substitutions","aaSubstitutions"]

def parse_muts(val):
    if pd.isna(val) or str(val).strip() in ("","nan"):
        return []
    return [m.strip() for m in str(val).split(",") if m.strip()]

wnd_nuc = Counter()   # (window_id, mut) -> n_with_mut
wnd_aa  = Counter()
wnd_tot = Counter()   # window_id -> n_seqs

for chunk in pd.read_csv(str(NEXTCLADE_TSV), sep="\t", usecols=NC_COLS,
                          chunksize=30_000, low_memory=False):
    sub = chunk[chunk["seqName"].isin(all_seqids)].merge(
        seq_map[["sequence_id","window_id"]], left_on="seqName",
        right_on="sequence_id", how="left").dropna(subset=["window_id"])
    for _, row in sub.iterrows():
        wid = row["window_id"]
        wnd_tot[wid] += 1
        for m in parse_muts(row["substitutions"]):
            if m in TRACK_NUC:
                wnd_nuc[(wid, m)] += 1
        for m in parse_muts(row["aaSubstitutions"]):
            if m in TRACK_AA:
                wnd_aa[(wid, m)] += 1

print(f"  Windows with data: {len(wnd_tot)}")

# Build trajectory table
def build_traj(counts, track_set, label):
    rows = []
    for (wid, mut), cnt in counts.items():
        if mut not in track_set:
            continue
        n = wnd_tot.get(wid, 0)
        if n < MIN_WINDOW_SEQS:
            continue
        meta = wpm_map.get(wid, {})
        rows.append({
            "mutation":      mut,
            "window_id":     wid,
            "window_idx":    meta.get("window_idx"),
            "wn_mid_date":   meta.get("wn_mid_date"),
            "policy_period": meta.get("policy_period"),
            "n_seqs":        n,
            "n_with_mut":    cnt,
            "frequency":     cnt / n,
            "mut_type":      label,
        })
    return pd.DataFrame(rows).sort_values(["mutation","window_idx"]) if rows else pd.DataFrame()

traj_nuc = build_traj(wnd_nuc, set(TRACK_NUC), "nucleotide")
traj_aa  = build_traj(wnd_aa,  set(TRACK_AA),  "amino_acid")
traj_all = pd.concat([traj_nuc, traj_aa], ignore_index=True)
traj_all["wn_mid_date"] = pd.to_datetime(traj_all["wn_mid_date"])
# Merge positive-test counts for downstream weighting and context
traj_all = traj_all.merge(win_tests, on="window_id", how="left")
traj_all.to_csv(OUT / "part4_mutation_trajectories.csv", index=False)
print(f"  Trajectories: {traj_all['mutation'].nunique()} mutations × {traj_all['window_id'].nunique()} windows")


# ── Growth rate estimation: binomial GLM ─────────────────────────────────
print("\nFitting growth models (binomial GLM) …")

anchor_date  = pd.Timestamp("2020-11-03 12:00:00")   # first F5 window
anchor_logit = np.log(0.0066 / (1 - 0.0066))          # S:N501Y freq at anchor

# Prepare S:N501Y series with test counts
n501 = traj_aa[traj_aa["mutation"]=="S:N501Y"].sort_values("window_idx").copy()
n501 = n501.merge(win_tests, on="window_id", how="left")
n501["wn_mid_date"] = pd.to_datetime(n501["wn_mid_date"])
n501["n_without"]   = n501["n_seqs"] - n501["n_with_mut"]


def prepare_alpha_glm_data(df_in, anchor):
    d = df_in[
        (df_in["frequency"] > 0.002) & (df_in["frequency"] < 0.998)
    ].dropna(subset=["wn_positive_tests"]).copy()

    d["days"] = (d["wn_mid_date"] - anchor).dt.days.astype(float)
    d["proportion_sequenced"] = d["n_seqs"] / d["wn_positive_tests"]
    return d


def fit_binomial_glm(
    df_in,
    label,
    anchor,
    *,
    model_label="positive_test_weighted",
    use_positive_test_weights=True,
    adjust_for_coverage=False,
    verbose=True,
):
    """Fit binomial GLM for S:N501Y frequency.

    The primary model fits logit(p) ~ days_since_anchor to sequenced
    S:N501Y counts, with per-window positive-test counts used as
    frequency weights.  Sensitivity models omit weights or include
    proportion sequenced as an additional covariate.

    Returns a dictionary with slope, CI, doubling time, pseudo-R2, and
    coverage-coefficient metadata where applicable.
    """
    d = prepare_alpha_glm_data(df_in, anchor)

    covariates = ["days"]
    if adjust_for_coverage:
        covariates.append("proportion_sequenced")
    X = sm.add_constant(d[covariates].values)
    # Binomial response: [successes, failures]
    y = np.column_stack([
        d["n_with_mut"].values.astype(float),
        d["n_without"].values.astype(float),
    ])

    fit_kwargs = {}
    if use_positive_test_weights:
        # Frequency weights: proportional to positive-test count in that window.
        w = d["wn_positive_tests"].values.astype(float)
        fit_kwargs["freq_weights"] = w / w.mean()

    result = sm.GLM(y, X, family=sm.families.Binomial(), **fit_kwargs).fit()

    slope     = result.params[1]
    intercept = result.params[0]
    se        = result.bse[1]
    ci_lo     = slope - 1.96 * se
    ci_hi     = slope + 1.96 * se
    pseudo_r2 = 1 - result.llf / result.llnull

    doubling = np.log(2) / slope if slope > 0 else np.nan
    ci_lo_d  = np.log(2) / ci_hi if ci_hi > 0 else np.nan
    ci_hi_d  = np.log(2) / ci_lo if ci_lo > 0 else np.nan

    coverage_coef = np.nan
    coverage_p = np.nan
    if adjust_for_coverage:
        coverage_coef = result.params[2]
        coverage_p = result.pvalues[2]

    out = {
        "phase": label,
        "model": model_label,
        "formula": "logit(S:N501Y) ~ days + proportion_sequenced"
                   if adjust_for_coverage else "logit(S:N501Y) ~ days",
        "positive_test_weighted": use_positive_test_weights,
        "n_windows": int(len(d)),
        "first_window": d["wn_mid_date"].min(),
        "last_window": d["wn_mid_date"].max(),
        "slope": slope,
        "intercept": intercept,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "doubling_days": doubling,
        "doubling_ci_lo": ci_lo_d,
        "doubling_ci_hi": ci_hi_d,
        "coverage_coef": coverage_coef,
        "coverage_p": coverage_p,
        "positive_tests_min": d["wn_positive_tests"].min(),
        "positive_tests_max": d["wn_positive_tests"].max(),
        "proportion_sequenced_min": d["proportion_sequenced"].min(),
        "proportion_sequenced_max": d["proportion_sequenced"].max(),
        "pseudo_R2": pseudo_r2,
    }
    if verbose:
        print(f"  {label} [{model_label}]: r={slope:.4f}/day  "
              f"(95% CI {ci_lo:.4f}-{ci_hi:.4f}), "
              f"doubling={doubling:.1f}d (CI {ci_lo_d:.1f}-{ci_hi_d:.1f}d), "
              f"pseudo-R2={pseudo_r2:.3f}")
    return out


# Fit F5 phase (3 Nov – 31 Dec, exclude 1 Dec dip anomaly)
f5_fit_data = n501[
    (n501["wn_mid_date"] >= "2020-11-03") &
    (n501["wn_mid_date"] <= "2020-12-31") &
    (n501["wn_mid_date"] != pd.Timestamp("2020-12-01 12:00:00"))
]
f5_primary = fit_binomial_glm(
    f5_fit_data,
    "Alpha under F5",
    anchor_date,
    model_label="positive_test_weighted",
    use_positive_test_weights=True,
    adjust_for_coverage=False,
)
r_f5 = f5_primary["slope"]
i_f5 = f5_primary["intercept"]
t0_f5 = f5_primary["first_window"]
ci_lo_f5 = f5_primary["ci_lo"]
ci_hi_f5 = f5_primary["ci_hi"]
pr2_f5 = f5_primary["pseudo_R2"]

# Fit L2 phase (5 Jan – 16 Mar)
l2_fit_data = n501[
    (n501["wn_mid_date"] >= "2021-01-05") &
    (n501["wn_mid_date"] <= "2021-03-16")
]
l2_primary = fit_binomial_glm(
    l2_fit_data,
    "Alpha under L2",
    anchor_date,
    model_label="positive_test_weighted",
    use_positive_test_weights=True,
    adjust_for_coverage=False,
)
r_l2 = l2_primary["slope"]
i_l2 = l2_primary["intercept"]
t0_l2 = l2_primary["first_window"]
ci_lo_l2 = l2_primary["ci_lo"]
ci_hi_l2 = l2_primary["ci_hi"]
pr2_l2 = l2_primary["pseudo_R2"]

# Sensitivity analysis: compare unweighted, positive-test weighted, and
# coverage-adjusted binomial GLMs.
print("\nGrowth model sensitivity (S:N501Y):")
sensitivity_rows = []
for phase_label, phase_data in [
    ("Alpha under F5", f5_fit_data),
    ("Alpha under L2", l2_fit_data),
]:
    sensitivity_rows.extend([
        fit_binomial_glm(
            phase_data,
            phase_label,
            anchor_date,
            model_label="unweighted_binomial",
            use_positive_test_weights=False,
            adjust_for_coverage=False,
        ),
        fit_binomial_glm(
            phase_data,
            phase_label,
            anchor_date,
            model_label="positive_test_weighted",
            use_positive_test_weights=True,
            adjust_for_coverage=False,
        ),
        fit_binomial_glm(
            phase_data,
            phase_label,
            anchor_date,
            model_label="coverage_adjusted",
            use_positive_test_weights=False,
            adjust_for_coverage=True,
        ),
    ])

sensitivity = pd.DataFrame(sensitivity_rows)
sensitivity["rate_ratio_L2_vs_F5"] = np.nan
sensitivity["pct_slower_L2_vs_F5"] = np.nan
for model_label, grp in sensitivity.groupby("model"):
    slopes = grp.set_index("phase")["slope"]
    if {"Alpha under F5", "Alpha under L2"}.issubset(slopes.index):
        ratio = slopes.loc["Alpha under L2"] / slopes.loc["Alpha under F5"]
        sensitivity.loc[sensitivity["model"] == model_label, "rate_ratio_L2_vs_F5"] = ratio
        sensitivity.loc[sensitivity["model"] == model_label, "pct_slower_L2_vs_F5"] = (1 - ratio) * 100

sensitivity.to_csv(OUT / "part4_growth_model_sensitivity.csv", index=False)

print("\nSensitivity summary (growth rate r/day):")
sens_print = sensitivity.pivot(index="model", columns="phase", values="slope")
sens_print["L2/F5"] = sens_print["Alpha under L2"] / sens_print["Alpha under F5"]
sens_print["L2 slower"] = (1 - sens_print["L2/F5"]) * 100
print(sens_print.to_string(float_format=lambda x: f"{x:.4f}"))

# B.1.177 decline under L2 — kept as exponential OLS (L2 coverage is higher,
# bias less severe; halving is well-determined by 8+ windows)
s222 = traj_aa[traj_aa["mutation"]=="S:A222V"].sort_values("window_idx").copy()
s222["wn_mid_date"] = pd.to_datetime(s222["wn_mid_date"])
l2_b1177 = s222[
    (s222["wn_mid_date"] >= "2021-01-05") &
    (s222["wn_mid_date"] <= "2021-03-09") &
    (s222["frequency"] > 0.005)
]
x_b = (l2_b1177["wn_mid_date"] - l2_b1177["wn_mid_date"].iloc[0]).dt.days.values
y_b = np.log(l2_b1177["frequency"].values)
slope_b, _, r_b, _, _ = linregress(x_b, y_b)
halving = np.log(0.5) / slope_b
print(f"  B.1.177 under L2 (OLS): halving={halving:.1f}d, R²={r_b**2:.3f}")

pct_reduction = (1 - r_l2 / r_f5) * 100
print(f"\n  Rate ratio L2/F5: {r_l2/r_f5:.1%}  →  L2 is {pct_reduction:.0f}% slower than F5")

# Save growth parameters (with 95% CI and pseudo-R²)
growth_params = pd.DataFrame([
    {
        "label":           "Alpha_F5",
        "slope":           r_f5,
        "intercept":       i_f5,
        "t0":              t0_f5,
        "ci_lo":           ci_lo_f5,
        "ci_hi":           ci_hi_f5,
        "doubling_days":   np.log(2) / r_f5,
        "doubling_ci_lo":  np.log(2) / ci_hi_f5,   # inverted
        "doubling_ci_hi":  np.log(2) / ci_lo_f5,
        "model":           "binomial_glm_positive_test_weighted",
        "pseudo_R2":       pr2_f5,
    },
    {
        "label":           "Alpha_L2",
        "slope":           r_l2,
        "intercept":       i_l2,
        "t0":              t0_l2,
        "ci_lo":           ci_lo_l2,
        "ci_hi":           ci_hi_l2,
        "doubling_days":   np.log(2) / r_l2,
        "doubling_ci_lo":  np.log(2) / ci_hi_l2,
        "doubling_ci_hi":  np.log(2) / ci_lo_l2,
        "model":           "binomial_glm_positive_test_weighted",
        "pseudo_R2":       pr2_l2,
    },
    {
        "label":           "B1177_L2_decline",
        "slope":           slope_b,
        "intercept":       None,
        "t0":              l2_b1177["wn_mid_date"].iloc[0],
        "ci_lo":           None,
        "ci_hi":           None,
        "doubling_days":   halving,
        "doubling_ci_lo":  None,
        "doubling_ci_hi":  None,
        "model":           "exponential_decline_ols",
        "pseudo_R2":       r_b**2,
    },
])
growth_params.to_csv(OUT / "part4_growth_params.csv", index=False)


# ── Counterfactual projections ─────────────────────────────────────────────
print("\nBuilding counterfactual projections …")
# anchor_date and anchor_logit already defined above

def logistic_proj(t_days, slope, logit_0):
    return 1 / (1 + np.exp(-(slope * t_days + logit_0)))

dates = pd.date_range("2020-11-03", "2021-03-08", freq="7D")

scenarios = {
    "Actual (F5 → L2 on 5 Jan)": ("2021-01-05", r_f5, r_l2),
    "L2 from 2 Nov (immediate)": ("2020-11-02", r_f5, r_l2),
    "L2 from 2 Dec":             ("2020-12-02", r_f5, r_l2),
    "L2 from 8 Dec (explosion)": ("2020-12-08", r_f5, r_l2),
}

proj_rows = []
for scen_name, (switch_str, r_before, r_after) in scenarios.items():
    switch_date = pd.Timestamp(switch_str)
    freq_at_switch = None
    logit_at_switch = None
    for d in dates:
        if d < switch_date:
            t = (d - anchor_date).days
            f = logistic_proj(t, r_before, anchor_logit)
        else:
            if freq_at_switch is None:
                t_sw = (switch_date - anchor_date).days
                freq_at_switch = logistic_proj(t_sw, r_before, anchor_logit)
                logit_at_switch = np.log(freq_at_switch / (1 - freq_at_switch))
            t2 = (d - switch_date).days
            f = logistic_proj(t2, r_after, logit_at_switch)
        proj_rows.append({
            "scenario":   scen_name,
            "date":       d,
            "frequency":  min(max(f, 0.0), 0.9999),
            "switch_date": switch_str,
        })

proj_df = pd.DataFrame(proj_rows)
proj_df.to_csv(OUT / "part4_counterfactual_projections.csv", index=False)

# Key milestone summary
pivot = proj_df.pivot(index="date", columns="scenario", values="frequency")
key_dates = pd.to_datetime(["2020-12-01","2020-12-08","2020-12-29","2021-01-05","2021-01-19","2021-02-02","2021-02-16"])
print("\nKey dates (Alpha % under each scenario):")
key_date_table = pivot.loc[[d for d in key_dates if d in pivot.index]]
format_pct = lambda x: f"{x:.1%}" if pd.notna(x) else ""
if hasattr(key_date_table, "map"):
    formatted_key_dates = key_date_table.map(format_pct)
else:
    formatted_key_dates = key_date_table.applymap(format_pct)
print(formatted_key_dates.to_string())

# Estimated delay to 50% dominance under each scenario
print("\nEstimated date Alpha reaches 50% under each scenario:")
for col in pivot.columns:
    above50 = pivot[pivot[col] >= 0.50].index
    if len(above50):
        print(f"  {col}: {above50[0].date()}")

print("\nStage 2 complete.")
