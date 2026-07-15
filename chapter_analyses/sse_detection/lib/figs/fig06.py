"""Build Chapter 5 Figure 6: two-axis candidate landscape."""

from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
from .common import (
    Paths,
    add_common_args,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)
from ..sse.io import write_table

FILE_NAME = "ch5_score_landscape"
CANDIDATE_COLORS = {
    "possible_review": "#E6AB02",
    "high_priority_burst": "#D55E00",
    "high_priority_burden": "#0072B2",
    "high_priority_both_axes": "#7B3294",
}

TIER_DISPLAY_NAMES = {
    "possible_review": "Possible Review",
    "high_priority_burst": "Burst",
    "high_priority_burden": "Burden",
    "high_priority_both_axes": "Both Axes",
}


def build_score_landscape(nodes: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "cluster_id",
        "window_id",
        "window_idx",
        "wn_mid_date",
        "policy_period",
        "who_voc",
        "clade",
        "cluster_size",
        "candidate_tier",
        "axes_fired",
        "burden_eligible",
        "burst_score",
        "burden_score",
        "burst_score_null_z",
        "burden_score_null_z",
        "burst_score_upper_p",
        "burden_score_upper_p",
    ]
    return nodes.loc[nodes["sse_tested"].fillna(False), cols].copy()


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_score_landscape(nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")

    x = pd.to_numeric(table["burst_score_null_z"], errors="coerce")
    y = pd.to_numeric(table["burden_score_null_z"], errors="coerce")
    
    # 1. Dynamic offset based on Y-axis data span
    finite_y = y[np.isfinite(y)]
    if len(finite_y) > 0:
        y_span = finite_y.max() - finite_y.min()
        strip_offset = y_span * 0.15  # 15% of the data range
        strip_y = float(finite_y.min() - strip_offset)
        line_y = strip_y + (strip_offset * 0.4)
    else:
        strip_y = -1.0
        line_y = -0.65

    plot_y = y.fillna(strip_y)
    
    fig, ax = styled_new_figure(width="double", height_in=5.2, constrained_layout=True)
    
    # Plot background
    background = table["candidate_tier"].eq("background_or_low_information")
    ax.scatter(
        x[background],
        plot_y[background],
        s=7,
        color="#B8B8B8",
        alpha=0.20,
        linewidth=0,
        rasterized=True,
    )
    
    # Plot priority tiers
    for tier, color in CANDIDATE_COLORS.items():
        tier_data = table[table["candidate_tier"] == tier]
        if not tier_data.empty:
            # 2. Extract specific coordinates
            tx = x[tier_data.index]
            ty = plot_y[tier_data.index]
            sizes = np.clip(np.sqrt(tier_data["cluster_size"]) * 5, 18, 90)
            
            ax.scatter(
                tx,
                ty,
                s=sizes,
                color=color,
                edgecolor="white",
                linewidth=0.6,          # Improved contrast outline
                alpha=0.9,              # Added subtle transparency for overlapping points
                label=TIER_DISPLAY_NAMES.get(tier, tier.title()),
                zorder=3,
            )
            
    # Draw separating line and text
    ax.axhline(line_y, color="#777777", lw=0.7, ls=":")
    ax.text(
        ax.get_xlim()[0], strip_y, "Burden N/A", va="center", ha="left", color="#555555"
    )

    ax.set_xlabel("Local-burst null-standardised score")
    ax.set_ylabel("Onward-burden null-standardised score")
    ax.legend(loc="center left")
    
    outputs = styled_save_figure(fig, paths, f"fig_{FILE_NAME}", tight=False)
    return {"figure": fig, "outputs": outputs}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    paths = paths_from_args(parser.parse_args())
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
