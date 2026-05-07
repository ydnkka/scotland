"""Supplementary Part 3 descriptive question tables.

This script adds secondary policy-context checks that are useful for thesis
interpretation but not part of the primary Part 3 analysis:

1. Do policy-intensity correlations change when intensity is lagged by 1-4
   weeks, including singleton fraction as an outcome?
2. Are ITS level-change estimates sensitive to the width of the pre/post
   window?
3. How close are policy-period starts to lineage dominance/overtake events in
   the surveillance time series?

Run from the repository root:

    conda run -n PhD python part3/supplementary_questions.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _bootstrap_repo_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in [here, *here.parents]:
        if (candidate / "config.yaml").exists():
            root_str = str(candidate)
            if root_str not in sys.path:
                sys.path.insert(0, root_str)
            return candidate
    raise FileNotFoundError("Could not locate config.yaml.")


ROOT = _bootstrap_repo_root()
TABLE_DIR = ROOT / "part3" / "tables"
SURVEILLANCE_TABLE_DIR = ROOT / "surveillance" / "tables"

from part3_analysis import (  # noqa: E402
    ITS_OUTCOMES,
    ITS_TRANSITIONS,
    _fit_its,
    _its_weekly_data,
    compute_weekly_summaries,
    load_and_prepare,
)
from utils.policy import POLICY_PERIODS_PD, PERIOD_ORDER  # noqa: E402


LAGGED_OUTCOMES = {
    "median_log_cluster_size": "Median log cluster size",
    "median_log_datazones": "Median log datazones",
    "pct_singleton": "Singleton fraction",
    "mean_simd_excess": "Mean SIMD excess discordance",
    "mean_age_excess": "Mean age excess discordance",
}


def compute_lagged_intensity_correlations(
    weekly: pd.DataFrame,
    *,
    max_lag_weeks: int = 4,
) -> pd.DataFrame:
    """Spearman correlations between outcomes and current/prior intensity."""
    weekly = weekly.sort_values("week_start").reset_index(drop=True).copy()
    records = []

    for lag in range(max_lag_weeks + 1):
        intensity = weekly["dominant_intensity"].shift(lag)
        for outcome, label in LAGGED_OUTCOMES.items():
            valid = pd.DataFrame(
                {
                    "intensity": intensity,
                    "outcome": weekly[outcome],
                }
            ).dropna()
            if len(valid) < 5:
                continue
            rho = valid["intensity"].corr(valid["outcome"], method="spearman")
            records.append(
                {
                    "outcome": outcome,
                    "outcome_label": label,
                    "intensity_lag_weeks": lag,
                    "spearman_rho": round(float(rho), 4),
                    "n_weeks": int(len(valid)),
                }
            )
    return pd.DataFrame(records)


def compute_its_window_sensitivity(
    df: pd.DataFrame,
    *,
    windows: tuple[int, ...] = (6, 8, 10, 12),
) -> pd.DataFrame:
    """Repeat ITS fits across alternative +/- week windows."""
    records = []

    for window in windows:
        for label, tdate, pre_code, post_code, description in ITS_TRANSITIONS:
            its_df = _its_weekly_data(df, tdate, window)
            for outcome, meta in ITS_OUTCOMES.items():
                fit = _fit_its(its_df, outcome)
                records.append(
                    {
                        "window_weeks_each_side": window,
                        "transition": label,
                        "transition_date": tdate.date(),
                        "pre_period": pre_code,
                        "post_period": post_code,
                        "description": description,
                        "outcome": outcome,
                        "outcome_label": meta["label"],
                        "n": fit.get("n", np.nan),
                        "r2": fit.get("r2", np.nan),
                        "coef_post": fit.get("coef_post", np.nan),
                        "ci_lo_post": fit.get("ci_lo_post", np.nan),
                        "ci_hi_post": fit.get("ci_hi_post", np.nan),
                        "pval_post": fit.get("pval_post", np.nan),
                        "coef_t_post": fit.get("coef_t_post", np.nan),
                        "pval_t_post": fit.get("pval_t_post", np.nan),
                        "error": fit.get("error", ""),
                    }
                )
    return pd.DataFrame(records)


def _nearest_dominance(
    dominance: pd.DataFrame,
    date: pd.Timestamp,
) -> pd.Series | None:
    if dominance.empty:
        return None
    idx = (dominance["time_period"] - date).abs().idxmin()
    return dominance.loc[idx]


def _overtake_context(
    overtakes: pd.DataFrame,
    date: pd.Timestamp,
) -> dict[str, object]:
    before = overtakes[overtakes["time_period"] <= date].tail(1)
    after = overtakes[overtakes["time_period"] >= date].head(1)

    context: dict[str, object] = {
        "previous_overtake_date": pd.NaT,
        "previous_overtake_to": pd.NA,
        "days_since_previous_overtake": np.nan,
        "next_overtake_date": pd.NaT,
        "next_overtake_to": pd.NA,
        "days_until_next_overtake": np.nan,
    }

    if not before.empty:
        row = before.iloc[0]
        context["previous_overtake_date"] = row["time_period"].date()
        context["previous_overtake_to"] = row["dominant_lineage_group"]
        context["days_since_previous_overtake"] = int((date - row["time_period"]).days)

    if not after.empty:
        row = after.iloc[0]
        context["next_overtake_date"] = row["time_period"].date()
        context["next_overtake_to"] = row["dominant_lineage_group"]
        context["days_until_next_overtake"] = int((row["time_period"] - date).days)

    return context


def compute_policy_lineage_context(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate observed policy starts with surveillance lineage context."""
    dominance_path = SURVEILLANCE_TABLE_DIR / "lineage_dominance_by_period.csv"
    overtakes_path = SURVEILLANCE_TABLE_DIR / "lineage_overtake_events.csv"
    if not dominance_path.exists() or not overtakes_path.exists():
        raise FileNotFoundError(
            "Run surveillance/policy_sequences_over_time.py before this script."
        )

    dominance = pd.read_csv(dominance_path, parse_dates=["time_period"])
    overtakes = pd.read_csv(overtakes_path, parse_dates=["time_period"])
    observed_periods = set(df["policy_period"].dropna().unique())
    its_by_post_period = {post: label for label, _, _, post, _ in ITS_TRANSITIONS}

    records = []
    for code in PERIOD_ORDER:
        if code not in observed_periods:
            continue
        row = POLICY_PERIODS_PD[POLICY_PERIODS_PD["period_code"] == code].iloc[0]
        start = pd.Timestamp(row["start_date"])
        nearest = _nearest_dominance(dominance, start)
        base = {
            "period_code": code,
            "period_label": row["period_label"],
            "period_start_date": start.date(),
            "policy_intensity": row["intensity"],
            "is_its_transition": code in its_by_post_period,
            "its_transition": its_by_post_period.get(code, ""),
        }
        if nearest is not None:
            base.update(
                {
                    "nearest_surveillance_week": nearest["time_period"].date(),
                    "dominant_lineage_near_start": nearest["dominant_lineage_group"],
                    "dominant_frequency_near_start": nearest["dominant_frequency"],
                    "days_from_nearest_surveillance_week": int(
                        (start - nearest["time_period"]).days
                    ),
                }
            )
        base.update(_overtake_context(overtakes, start))
        records.append(base)

    return pd.DataFrame(records)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    print("Part 3 supplementary questions")
    print("=" * 34)
    print("Loading Part 3 cluster table...")
    df = load_and_prepare()
    weekly = compute_weekly_summaries(df)

    lagged = compute_lagged_intensity_correlations(weekly)
    lagged_path = TABLE_DIR / "supp_lagged_intensity_correlations.csv"
    lagged.to_csv(lagged_path, index=False)
    print(f"  Saved {lagged_path.relative_to(ROOT)}")

    sensitivity = compute_its_window_sensitivity(df)
    sensitivity_path = TABLE_DIR / "supp_its_window_sensitivity.csv"
    sensitivity.to_csv(sensitivity_path, index=False)
    print(f"  Saved {sensitivity_path.relative_to(ROOT)}")

    lineage_context = compute_policy_lineage_context(df)
    lineage_path = TABLE_DIR / "supp_policy_lineage_context.csv"
    lineage_context.to_csv(lineage_path, index=False)
    print(f"  Saved {lineage_path.relative_to(ROOT)}")

    print("\nQuick checks:")
    singleton_lag0 = lagged[
        (lagged["outcome"] == "pct_singleton")
        & (lagged["intensity_lag_weeks"] == 0)
    ].iloc[0]
    print(
        "  Policy intensity vs singleton fraction, lag 0: "
        f"rho={singleton_lag0['spearman_rho']:.2f} "
        f"(n={int(singleton_lag0['n_weeks'])} weeks)."
    )
    print(
        "  ITS sensitivity table covers "
        f"{sensitivity['window_weeks_each_side'].nunique()} window widths."
    )
    print("Part 3 supplementary tables complete.")


if __name__ == "__main__":
    main()
