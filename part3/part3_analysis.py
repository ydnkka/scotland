"""Part 3: Policy period associations with SARS-CoV-2 genomic cluster structure.

This script characterises how Scottish government COVID-19 policy restriction
periods were associated with genomic cluster size, geographic dispersion, and
demographic mixing.  The analysis is explicitly descriptive; policy periods are
strongly confounded with variant waves and calendar time, so causal inference
is not appropriate.

Analytical components
---------------------
1. Period-level descriptive tables — median cluster size, datazones, mixing
   indices, singleton fraction, and policy intensity for each of the 16 policy
   periods observed in the study data.
2. Weekly aggregate series — ISO-week summaries of cluster outcomes annotated
   with the dominant policy period and its intensity for all downstream figure
   production.
3. Interrupted time-series (ITS) analyses — three pre/post analyses at
   transitions that occur within a relatively stable variant context, where
   the most acute variant-wave confounding is reduced:

     T1-onset  (2020-10-02):  Route-map phase 3 → Pre-tier tightening
                              B.1.177 era, intensity 30 → 55
     L2→SL     (2021-04-02):  Second lockdown → Stay-local Level 3
                              Alpha-dominant period, intensity 95 → 65
     NN-onset  (2021-08-09):  Level 0 → Near-normal (full legal easing)
                              Delta-dominant period, intensity 20 → 10

   For each transition a ±8-week ISO-week window is used. Outcomes are weekly
   medians of log cluster size and log datazones (non-singleton clusters), and
   weekly means of SIMD and age excess-discordance scores (non-singletons with
   valid mixing data). Four OLS segmented-regression (ITS) models are fit per
   transition:

     y_t = β0 + β1·t + β2·D_t + β3·(D_t·t) + ε_t

   where t is the signed week offset from the transition (negative = pre),
   D_t is a 0/1 post-transition indicator, and D_t·t captures slope change.
   Coefficients, 95 % CIs, and p-values are saved for each outcome × transition.

Run from the repository root:

    conda run -n PhD python part3/part3_analysis.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm


# ---------------------------------------------------------------------------
# Bootstrap repo root so utils is importable when run as a script
# ---------------------------------------------------------------------------

def _bootstrap_root() -> Path:
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "config.yaml").exists():
            root = str(cand)
            if root not in sys.path:
                sys.path.insert(0, root)
            return cand
    raise FileNotFoundError("Cannot locate config.yaml.")


ROOT = _bootstrap_root()

from utils.data import load_main_cluster_table          # noqa: E402
from utils.policy import attach_period_pandas, POLICY_PERIODS_PD, PERIOD_ORDER  # noqa: E402


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Study start: P3 (Route-map phase 3) onset — first policy period with data.
STUDY_START = pd.Timestamp("2020-07-10")

# ITS transition definitions: (label, transition_date, pre_period_code, post_period_code)
ITS_TRANSITIONS = [
    (
        "T1_onset",
        pd.Timestamp("2020-10-02"),
        "P3",
        "T1",
        "Route-map phase 3 → Pre-tier tightening",
    ),
    (
        "L2_to_SL",
        pd.Timestamp("2021-04-02"),
        "L2",
        "SL",
        "Second lockdown → Stay-local Level 3",
    ),
    (
        "NN_onset",
        pd.Timestamp("2021-08-09"),
        "L0",
        "NN",
        "Level 0 → Near-normal",
    ),
]

ITS_WINDOW_WEEKS = 8

# Cluster outcomes used in the ITS
ITS_OUTCOMES = {
    "log_cluster_size": {
        "label": "Median log cluster size",
        "agg": "median",
        "filter": "non_singleton",
    },
    "log_datazones": {
        "label": "Median log datazones",
        "agg": "median",
        "filter": "non_singleton",
    },
    "simd_excess_discordance": {
        "label": "Mean SIMD excess discordance",
        "agg": "mean",
        "filter": "non_singleton_mixing",
    },
    "age_excess_discordance": {
        "label": "Mean age excess discordance",
        "agg": "mean",
        "filter": "non_singleton_mixing",
    },
}

OUT_DIR = ROOT / "part3" / "tables"


# ---------------------------------------------------------------------------
# Data loading and preparation
# ---------------------------------------------------------------------------

def load_and_prepare() -> pd.DataFrame:
    """Load the main cluster table and attach policy periods.

    Returns the full cluster-level DataFrame (all clusters, all periods) with
    ``policy_period``, ``policy_period_label``, and ``policy_intensity`` columns
    added, filtered to the study start date.
    """
    df = load_main_cluster_table(root=ROOT)

    # Filter to study window — data begins mid-July 2020 (P3).
    df = df[df["wn_mid_date"] >= STUDY_START].copy()

    # Attach policy periods using window midpoint date.
    df = attach_period_pandas(df, "wn_mid_date")

    # Derived columns used across analyses.
    df["log_datazones"] = np.log(df["cluster_n_datazones"].clip(lower=1))
    # log_cluster_size already present in the main cluster table.

    df["is_non_singleton"] = df["cluster_size"] > 1
    df["week_start"] = df["wn_mid_date"].dt.to_period("W").dt.start_time

    return df


# ---------------------------------------------------------------------------
# Section 1: Period-level descriptive table
# ---------------------------------------------------------------------------

def compute_period_descriptives(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-policy-period summary statistics for all clusters.

    Returns one row per observed policy period, ordered chronologically.
    """
    records = []

    for code in PERIOD_ORDER:
        sub = df[df["policy_period"] == code]
        if sub.empty:
            continue

        ns = sub[sub["is_non_singleton"]]
        n_total = len(sub)
        n_ns = len(ns)
        pct_singleton = 100.0 * (1 - n_ns / n_total) if n_total > 0 else np.nan

        # SIMD/age mixing: non-singleton clusters with valid discordance scores
        ns_simd = ns.dropna(subset=["simd_excess_discordance"])
        ns_age  = ns.dropna(subset=["age_excess_discordance"])

        rec = {
            "period_code":          code,
            "period_label":         sub["policy_period_label"].iloc[0],
            "policy_intensity":     sub["policy_intensity"].iloc[0],
            "n_clusters_total":     n_total,
            "n_clusters_nonsingleton": n_ns,
            "pct_singleton":        round(pct_singleton, 1),
            # cluster size (non-singleton)
            "median_cluster_size":       round(ns["cluster_size"].median(), 1) if n_ns else np.nan,
            "iqr_cluster_size_lo":       round(ns["cluster_size"].quantile(0.25), 1) if n_ns else np.nan,
            "iqr_cluster_size_hi":       round(ns["cluster_size"].quantile(0.75), 1) if n_ns else np.nan,
            # datazones (non-singleton)
            "median_datazones":          round(ns["cluster_n_datazones"].median(), 1) if n_ns else np.nan,
            "iqr_datazones_lo":          round(ns["cluster_n_datazones"].quantile(0.25), 1) if n_ns else np.nan,
            "iqr_datazones_hi":          round(ns["cluster_n_datazones"].quantile(0.75), 1) if n_ns else np.nan,
            # mixing (non-singleton, valid observations)
            "mean_simd_excess_discordance": round(ns_simd["simd_excess_discordance"].mean(), 4) if len(ns_simd) else np.nan,
            "mean_age_excess_discordance":  round(ns_age["age_excess_discordance"].mean(), 4) if len(ns_age) else np.nan,
            "n_simd_valid":          len(ns_simd),
            "n_age_valid":           len(ns_age),
        }
        records.append(rec)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Section 2: Weekly aggregate series
# ---------------------------------------------------------------------------

def compute_weekly_summaries(df: pd.DataFrame) -> pd.DataFrame:
    """Compute ISO-week cluster summaries with policy period annotation.

    Returns one row per ISO week containing the median/mean cluster outcomes
    for that week, together with the dominant policy period and its intensity
    (determined by the modal policy_period value in the week).
    """
    ns = df[df["is_non_singleton"]].copy()

    # Aggregate non-singleton cluster metrics per ISO week.
    agg_ns = (
        ns.groupby("week_start", sort=True)
        .agg(
            n_clusters_nonsingleton=("cluster_id", "count"),
            median_cluster_size=("cluster_size", "median"),
            median_log_cluster_size=("log_cluster_size", "median"),
            median_datazones=("cluster_n_datazones", "median"),
            median_log_datazones=("log_datazones", "median"),
            mean_simd_excess=(
                "simd_excess_discordance",
                lambda x: x.dropna().mean() if x.notna().any() else np.nan,
            ),
            mean_age_excess=(
                "age_excess_discordance",
                lambda x: x.dropna().mean() if x.notna().any() else np.nan,
            ),
        )
        .reset_index()
    )

    # Total cluster count (all clusters, including singletons) per week.
    agg_all = (
        df.groupby("week_start", sort=True)
        .agg(
            n_clusters_total=("cluster_id", "count"),
            pct_singleton=(
                "is_non_singleton",
                lambda x: 100.0 * (1 - x.mean()),
            ),
        )
        .reset_index()
    )

    # Dominant policy period per week (modal value).
    policy_mode = (
        df.groupby("week_start")["policy_period"]
        .agg(lambda s: s.mode().iloc[0] if not s.mode().empty else None)
        .reset_index()
        .rename(columns={"policy_period": "dominant_period_code"})
    )
    policy_mode = policy_mode.merge(
        POLICY_PERIODS_PD[["period_code", "period_label", "intensity"]],
        left_on="dominant_period_code",
        right_on="period_code",
        how="left",
    ).drop(columns=["period_code"]).rename(
        columns={"period_label": "dominant_period_label",
                 "intensity":    "dominant_intensity"}
    )

    result = (
        agg_all
        .merge(agg_ns, on="week_start", how="left")
        .merge(policy_mode, on="week_start", how="left")
        .sort_values("week_start")
        .reset_index(drop=True)
    )
    return result


# ---------------------------------------------------------------------------
# Section 3: Intensity correlation with weekly cluster outcomes
# ---------------------------------------------------------------------------

def compute_intensity_correlations(weekly: pd.DataFrame) -> pd.DataFrame:
    """Compute Spearman correlations between policy intensity and cluster outcomes.

    Correlations are computed (a) pooled across all weeks, and (b) separately
    within each dominant policy period to give a within-period partial picture.
    """
    outcome_cols = [
        "median_log_cluster_size",
        "median_log_datazones",
        "mean_simd_excess",
        "mean_age_excess",
    ]
    records = []
    for col in outcome_cols:
        valid = weekly[["dominant_intensity", col]].dropna()
        if len(valid) < 5:
            continue
        rho = valid["dominant_intensity"].corr(valid[col], method="spearman")
        records.append({
            "outcome": col,
            "spearman_rho_pooled": round(rho, 4),
            "n_weeks": len(valid),
        })
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Section 4: ITS analyses
# ---------------------------------------------------------------------------

def _its_weekly_data(
    df: pd.DataFrame,
    transition_date: pd.Timestamp,
    window_weeks: int,
) -> pd.DataFrame:
    """Extract and aggregate the ITS analysis window for one transition.

    Returns a weekly-level DataFrame with signed week offset (t), post-transition
    indicator (post), and outcome variables.
    """
    half = pd.Timedelta(weeks=window_weeks)
    win_start = transition_date - half
    win_end   = transition_date + half - pd.Timedelta(days=1)

    sub = df[(df["wn_mid_date"] >= win_start) & (df["wn_mid_date"] <= win_end)].copy()
    ns  = sub[sub["is_non_singleton"]].copy()

    agg = (
        ns.groupby("week_start", sort=True)
        .agg(
            n_nonsingleton=("cluster_id", "count"),
            log_cluster_size=("log_cluster_size", "median"),
            log_datazones=("log_datazones", "median"),
            simd_excess_discordance=(
                "simd_excess_discordance",
                lambda x: x.dropna().mean() if x.notna().any() else np.nan,
            ),
            age_excess_discordance=(
                "age_excess_discordance",
                lambda x: x.dropna().mean() if x.notna().any() else np.nan,
            ),
        )
        .reset_index()
    )

    # Signed week offset from transition (0 = first week at or after transition).
    t_dates = agg["week_start"].dt.tz_localize(None)
    tdate = transition_date.tz_localize(None) if transition_date.tzinfo else transition_date

    # Week number relative to transition week.
    agg["t"] = ((t_dates - tdate).dt.days / 7).round().astype(int)
    # Normalise so that pre-transition is negative, first post-week is 0.
    agg["post"] = (agg["t"] >= 0).astype(int)
    # Interaction for slope change after transition.
    agg["t_post"] = agg["t"] * agg["post"]

    return agg.sort_values("t").reset_index(drop=True)


def _fit_its(data: pd.DataFrame, outcome_col: str) -> dict:
    """Fit a segmented OLS ITS model for one outcome.

    Model: y ~ 1 + t + post + t_post

    Returns a dict with coefficient estimates, 95 % CIs, p-values, and R².
    """
    valid = data[["t", "post", "t_post", outcome_col]].dropna()
    if len(valid) < 6:
        return {"n": len(valid), "error": "too few observations"}

    X = sm.add_constant(valid[["t", "post", "t_post"]])
    y = valid[outcome_col]
    model = sm.OLS(y, X).fit()

    result = {"n": len(valid), "r2": round(model.rsquared, 4)}
    for term in ["const", "t", "post", "t_post"]:
        idx = model.params.index.get_loc(term)
        result[f"coef_{term}"]    = round(model.params.iloc[idx], 6)
        result[f"ci_lo_{term}"]   = round(model.conf_int().iloc[idx, 0], 6)
        result[f"ci_hi_{term}"]   = round(model.conf_int().iloc[idx, 1], 6)
        result[f"pval_{term}"]    = round(model.pvalues.iloc[idx], 6)
    return result


def run_its_analyses(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Run ITS analyses for all three transitions and all four outcomes.

    Returns
    -------
    coef_table:
        Long-format DataFrame of ITS coefficients (one row per transition ×
        outcome × parameter).
    weekly_data:
        Dict mapping transition label to the weekly ITS data DataFrame.
    """
    coef_records = []
    weekly_data: dict[str, pd.DataFrame] = {}

    for label, tdate, pre_code, post_code, description in ITS_TRANSITIONS:
        its_df = _its_weekly_data(df, tdate, ITS_WINDOW_WEEKS)
        weekly_data[label] = its_df

        for outcome, meta in ITS_OUTCOMES.items():
            fit = _fit_its(its_df, outcome)
            base = {
                "transition":        label,
                "transition_date":   tdate.date(),
                "pre_period":        pre_code,
                "post_period":       post_code,
                "description":       description,
                "outcome":           outcome,
                "outcome_label":     meta["label"],
            }
            base.update(fit)
            coef_records.append(base)

    coef_table = pd.DataFrame(coef_records)
    return coef_table, weekly_data


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Part 3 analysis — policy period associations")
    print("=" * 55)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Load data ---
    print("Loading main cluster table and attaching policy periods…")
    df = load_and_prepare()
    print(f"  {len(df):,} cluster-window observations (from {STUDY_START.date()})")
    print(f"  {df['is_non_singleton'].sum():,} non-singleton clusters")
    print(f"  Policy periods observed: {sorted(df['policy_period'].dropna().unique())}")

    # --- Section 1: Period descriptives ---
    print("\nComputing period-level descriptive table…")
    period_desc = compute_period_descriptives(df)
    out_path = OUT_DIR / "period_descriptives.csv"
    period_desc.to_csv(out_path, index=False)
    print(f"  Saved → {out_path.relative_to(ROOT)}")
    print(period_desc[["period_code", "policy_intensity", "n_clusters_total",
                         "median_cluster_size", "median_datazones"]].to_string(index=False))

    # --- Section 2: Weekly summaries ---
    print("\nComputing weekly aggregate summaries…")
    weekly = compute_weekly_summaries(df)
    out_path = OUT_DIR / "weekly_summaries.csv"
    weekly.to_csv(out_path, index=False)
    print(f"  {len(weekly)} ISO weeks covered → {out_path.relative_to(ROOT)}")

    # --- Section 3: Intensity correlations ---
    print("\nComputing policy-intensity vs outcome correlations…")
    corr_table = compute_intensity_correlations(weekly)
    out_path = OUT_DIR / "intensity_correlations.csv"
    corr_table.to_csv(out_path, index=False)
    print(corr_table.to_string(index=False))
    print(f"  Saved → {out_path.relative_to(ROOT)}")

    # --- Section 4: ITS analyses ---
    print("\nRunning ITS analyses (3 transitions × 4 outcomes)…")
    coef_table, weekly_data = run_its_analyses(df)

    out_path = OUT_DIR / "its_coefficients.csv"
    coef_table.to_csv(out_path, index=False)
    print(f"  ITS coefficients saved → {out_path.relative_to(ROOT)}")

    # Save per-transition weekly data for figures
    for label, its_df in weekly_data.items():
        out_path = OUT_DIR / f"its_weekly_{label}.csv"
        its_df.to_csv(out_path, index=False)
        print(f"  ITS window data ({label}) → {out_path.relative_to(ROOT)}")

    # Print ITS summary
    print("\nITS level-change estimates (β_post, 95 % CI, p):")
    summary_cols = [
        "transition", "outcome", "coef_post", "ci_lo_post", "ci_hi_post", "pval_post",
    ]
    if all(c in coef_table.columns for c in summary_cols):
        print(coef_table[summary_cols].to_string(index=False))

    print("\nPart 3 analysis complete.")


if __name__ == "__main__":
    main()
