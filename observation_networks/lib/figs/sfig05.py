"""Build Chapter 4 Supplementary Figure 5: assortativity confidence intervals."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    ATTRIBUTE_ORDER,
    Paths,
    add_common_args,
    configure_matplotlib,
    paths_from_args,
    read_table,
    save_figure,
    window_idx_from_id,
)


def weighted_mean_ci_from_se(
    values: pd.Series,
    weights: pd.Series,
    standard_errors: pd.Series,
) -> dict[str, float]:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return {
            "weighted_mean": np.nan,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "ci_weight_share": np.nan,
        }

    values = values.loc[mask].astype(float)
    weights = weights.loc[mask].astype(float)
    weighted_mean_value = float(np.average(values, weights=weights))

    se_mask = standard_errors.loc[mask].notna()
    if not se_mask.any():
        return {
            "weighted_mean": weighted_mean_value,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "ci_weight_share": np.nan,
        }

    ci_weights = weights.loc[se_mask]
    ci_standard_errors = standard_errors.loc[mask].loc[se_mask].astype(float)
    normalized = ci_weights / weights.sum()
    combined_se = float(np.sqrt(np.sum((normalized * ci_standard_errors) ** 2)))
    return {
        "weighted_mean": weighted_mean_value,
        "combined_se": combined_se,
        "ci_low": weighted_mean_value - 1.96 * combined_se,
        "ci_high": weighted_mean_value + 1.96 * combined_se,
        "ci_weight_share": float(ci_weights.sum() / weights.sum()),
    }


def compatibility_assortativity_filtered(paths: Paths) -> pd.DataFrame:
    assort = read_table(paths, "compatibility_assortativity")
    assort["window_idx"] = window_idx_from_id(assort["window_id"])
    return assort.loc[
        assort["assortativity"].notna()
        & assort["edge_weight_total"].gt(0)
        & assort["n_categories"].gt(1)
        & assort["n_edge_contributions_used"].ge(20)
    ].copy()


def compatibility_window_assortativity(paths: Paths) -> pd.DataFrame:
    work = compatibility_assortativity_filtered(paths)
    rows = []
    for (window_idx, attribute, label), group in work.groupby(
        ["window_idx", "attribute", "attribute_label"], dropna=False
    ):
        ci = weighted_mean_ci_from_se(
            group["assortativity"],
            group["edge_weight_total"],
            group["assortativity_se"],
        )
        rows.append(
            {
                "window_idx": window_idx,
                "attribute": attribute,
                "attribute_label": label,
                "assortativity": ci["weighted_mean"],
                "assortativity_se": ci["combined_se"],
                "assortativity_ci_low": ci["ci_low"],
                "assortativity_ci_high": ci["ci_high"],
                "ci_weight_share": ci["ci_weight_share"],
                "edge_weight_total": group["edge_weight_total"].sum(),
            }
        )
    return pd.DataFrame(rows)


def compatibility_attribute_summary(paths: Paths) -> pd.DataFrame:
    work = compatibility_assortativity_filtered(paths)
    window_summary = compatibility_window_assortativity(paths)
    rows = []
    for (attribute, label), group in work.groupby(
        ["attribute", "attribute_label"], dropna=False
    ):
        ci = weighted_mean_ci_from_se(
            group["assortativity"],
            group["edge_weight_total"],
            group["assortativity_se"],
        )
        window_group = window_summary.loc[window_summary["attribute"].eq(attribute)]
        rows.append(
            {
                "attribute": attribute,
                "attribute_label": label,
                "weighted_mean": ci["weighted_mean"],
                "combined_se": ci["combined_se"],
                "ci_low": ci["ci_low"],
                "ci_high": ci["ci_high"],
                "window_median": window_group["assortativity"].median(),
                "window_q10": window_group["assortativity"].quantile(0.10),
                "window_q90": window_group["assortativity"].quantile(0.90),
            }
        )
    return pd.DataFrame(rows)


def build(paths: Paths) -> None:
    summary = (
        compatibility_attribute_summary(paths)
        .set_index("attribute_label")
        .reindex(ATTRIBUTE_ORDER)
        .dropna(subset=["weighted_mean"])
        .reset_index()
    )
    y_positions = np.arange(len(summary))
    xerr = [
        summary["weighted_mean"] - summary["ci_low"],
        summary["ci_high"] - summary["weighted_mean"],
    ]

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    ax.errorbar(
        summary["weighted_mean"],
        y_positions,
        xerr=xerr,
        fmt="o",
        ms=4,
        lw=1.1,
        capsize=3,
        color="#1f4e79",
    )
    ax.axvline(0, color="#777777", lw=0.8, ls=":")
    ax.set_yticks(y_positions)
    ax.set_yticklabels(summary["attribute_label"])
    ax.invert_yaxis()
    ax.set_xlabel("Edge-weighted mean compatibility assortativity")
    ax.set_title("Compatibility graph, 95% node-block jackknife CI")
    x_min = min(float(summary["ci_low"].min()), 0.0)
    x_max = float(summary["ci_high"].max())
    padding = max((x_max - x_min) * 0.08, 0.001)
    ax.set_xlim(x_min - padding, x_max + padding)
    save_figure(fig, paths, "sfig05_assortativity_confidence_intervals")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    configure_matplotlib()
    paths = paths_from_args(args)
    build(paths)
    print(f"Wrote sfig05_assortativity_confidence_intervals to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

