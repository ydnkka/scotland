"""Shared policy-index loading and correlation summaries."""

from __future__ import annotations

import numpy as np
import pandas as pd

from utils import load_daily_policy_data


def load_policy_indices(
    start_date: str | pd.Timestamp | None = None,
    end_date: str | pd.Timestamp | None = None,
) -> pd.DataFrame:
    """Load and optionally date-filter Scotland's processed policy table."""
    daily = load_daily_policy_data()
    if start_date is not None:
        daily = daily.loc[daily["date"].ge(pd.Timestamp(start_date))]
    if end_date is not None:
        daily = daily.loc[daily["date"].le(pd.Timestamp(end_date))]
    return daily.sort_values("date", ignore_index=True)


def build_correlation_summary(
    daily: pd.DataFrame,
) -> pd.DataFrame:
    """Summarise complete-day agreement without independence-based p-values."""
    complete = daily.dropna(subset=["stringency_index", "containment_index"]).copy()
    if len(complete) < 2:
        raise ValueError("At least two complete daily index pairs are required.")

    x = complete["stringency_index"].to_numpy(dtype=float)
    y = complete["containment_index"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, deg=1)
    pearson = float(complete["stringency_index"].corr(complete["containment_index"]))
    spearman = float(
        complete["stringency_index"].corr(
            complete["containment_index"],
            method="spearman",
        )
    )
    return pd.DataFrame(
        [
            {
                "start_date": complete["date"].min(),
                "end_date": complete["date"].max(),
                "n_complete_days": len(complete),
                "pearson_r": pearson,
                "spearman_rho": spearman,
                "linear_slope": float(slope),
                "linear_intercept": float(intercept),
                "pearson_r_squared": pearson**2,
            }
        ]
    )
