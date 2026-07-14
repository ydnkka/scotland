from __future__ import annotations

import logging

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #
DIMENSION = ["burst", "burden"]  # expands to <axis>_score, _n, _upper_p

SIZE_COL = "cluster_size"
TIER_COL = "candidate_tier"

HIGH_PRIORITY_TIERS = (
    "high_priority_both_axes",
    "high_priority_burst",
    "high_priority_burden",
)
INELIGIBLE = "size_ineligible"

SPIKE_THRESHOLD = 0.95  # top bin of a 20-bin histogram
EXPECTED_SPIKE = 1 - SPIKE_THRESHOLD
LOGGER = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def axis_cols(axis: str) -> dict[str, str]:
    return {
        "score": f"{axis}_score",
        "n": f"{axis}_score_n",
        "p": f"{axis}_score_upper_p",
    }


def background_frame(df: pd.DataFrame, p_col: str) -> pd.DataFrame:
    """Eligible, non-high-priority nodes with a valid p on the given axis."""
    eligible = df[TIER_COL] != INELIGIBLE
    not_candidate = ~df[TIER_COL].isin(HIGH_PRIORITY_TIERS)
    bg = df.loc[eligible & not_candidate].copy()
    return bg[bg[p_col].notna()]


def spike_fraction(p: pd.Series) -> float:
    if len(p) == 0:
        return np.nan
    return float((p > SPIKE_THRESHOLD).mean())


def summarise_by_component_count(
    bg: pd.DataFrame, *, n_col: str, p_col: str
) -> pd.DataFrame:
    g = bg.groupby(n_col, dropna=False)
    out = g.agg(
        n_nodes=(p_col, "size"),
        median_p=(p_col, "median"),
        mean_p=(p_col, "mean"),
        frac_in_spike=(p_col, spike_fraction),
    ).reset_index()
    out["expected_spike_frac"] = EXPECTED_SPIKE
    out["spike_excess_ratio"] = out["frac_in_spike"] / EXPECTED_SPIKE
    return out.sort_values(n_col)


def summarise_by_size(bg: pd.DataFrame, *, p_col: str) -> pd.DataFrame:
    try:
        bucket = pd.qcut(bg[SIZE_COL], q=6, duplicates="drop")
    except ValueError:
        bucket = pd.cut(bg[SIZE_COL], bins=6)
    g = bg.groupby(bucket, dropna=False, observed=True)
    out = g.agg(
        n_nodes=(p_col, "size"),
        size_min=(SIZE_COL, "min"),
        size_max=(SIZE_COL, "max"),
        median_p=(p_col, "median"),
        frac_in_spike=(p_col, spike_fraction),
    ).reset_index(drop=True)
    out["expected_spike_frac"] = EXPECTED_SPIKE
    out["spike_excess_ratio"] = out["frac_in_spike"] / EXPECTED_SPIKE
    return out


def overall_uniformity(
    bg: pd.DataFrame, *, p_col: str, n_bins: int = 20
) -> pd.DataFrame:
    counts, edges = np.histogram(bg[p_col], bins=n_bins, range=(0, 1))
    expected = len(bg) / n_bins
    return pd.DataFrame(
        {
            "bin_lo": edges[:-1],
            "bin_hi": edges[1:],
            "observed": counts,
            "expected": expected,
            "obs_over_exp": counts / expected,
        }
    )


def report_axis(df: pd.DataFrame, axis: str) -> None:
    cols = axis_cols(axis)
    p_col, n_col = cols["p"], cols["n"]

    bg = background_frame(df, p_col)
    LOGGER.info("\n%s\n# AXIS: %s  (p = %s)\n%s", "#" * 70, axis, p_col, "#" * 70)
    LOGGER.info("Background nodes with valid p: %s", len(bg))
    LOGGER.info(
        "Spike bin = p > %s (expected frac = %.3f)\n",
        SPIKE_THRESHOLD,
        EXPECTED_SPIKE,
    )

    LOGGER.info("=== Overall background p histogram (obs/exp per bin) ===")
    with pd.option_context("display.float_format", lambda v: f"{v:.3f}"):
        LOGGER.info("%s", overall_uniformity(bg, p_col=p_col).to_string(index=False))

    LOGGER.info("\n=== Background p by number of present components ===")
    by_n = summarise_by_component_count(bg, n_col=n_col, p_col=p_col)
    with pd.option_context("display.float_format", lambda v: f"{v:.3f}"):
        LOGGER.info("%s", by_n.to_string(index=False))

    LOGGER.info("\n=== Background p by cluster-size bucket ===")
    with pd.option_context("display.float_format", lambda v: f"{v:.3f}"):
        LOGGER.info("%s", summarise_by_size(bg, p_col=p_col).to_string(index=False))


def report_tier_sizes(df: pd.DataFrame) -> None:
    """Cluster-size distribution within each high-priority tier.

    Checks whether burst/burden candidates pile up at the size floor, where
    low-information calls are most likely.
    """
    LOGGER.info(
        "\n%s\n# Cluster-size distribution within high-priority tiers\n%s",
        "#" * 70,
        "#" * 70,
    )
    sub = df[df[TIER_COL].isin(HIGH_PRIORITY_TIERS)]
    summary = (
        sub.groupby(TIER_COL, observed=True)[SIZE_COL]
        .agg(
            n="size",
            size_min="min",
            size_p25=lambda s: s.quantile(0.25),
            size_median="median",
            size_p75=lambda s: s.quantile(0.75),
            size_max="max",
        )
        .reset_index()
    )
    with pd.option_context("display.float_format", lambda v: f"{v:.1f}"):
        LOGGER.info("%s", summary.to_string(index=False))

    # Fraction of each tier sitting in the bottom size band (floor noise check).
    floor = df.loc[df[TIER_COL] != INELIGIBLE, SIZE_COL].min()
    near_floor = sub[SIZE_COL].le(floor + 2)
    frac = (
        sub.assign(near_floor=near_floor)
        .groupby(TIER_COL, observed=True)["near_floor"]
        .mean()
    )
    LOGGER.info("\nFraction within 2 of the size floor (%.0f):", floor)
    with pd.option_context("display.float_format", lambda v: f"{v:.3f}"):
        LOGGER.info("%s", frac.to_string())


def report_axis_correlation(df: pd.DataFrame) -> None:
    score_cols = [f"{axis}_score" for axis in DIMENSION if f"{axis}_score" in df.columns]
    if len(score_cols) < 2:
        return

    LOGGER.info("\nPairwise score correlations:")
    for i, left in enumerate(score_cols):
        for right in score_cols[i + 1 :]:
            s = df[[left, right]].dropna()
            if len(s) < 2:
                continue
            r = s[left].corr(s[right])
            LOGGER.info("corr(%s, %s) = %.3f  (n=%s)", left, right, r, len(s))


def main(df: pd.DataFrame) -> None:
    for axis in DIMENSION:
        report_axis(df, axis)
    report_tier_sizes(df)
    report_axis_correlation(df)
