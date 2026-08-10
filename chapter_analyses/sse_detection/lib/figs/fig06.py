"""Build Chapter 5 Figure 6: score landscape and selection funnel."""

from __future__ import annotations

import argparse
from textwrap import fill

import numpy as np
import pandas as pd

from ..sse.config import MIN_CLUSTER_SIZE
from ..sse.io import write_table
from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
    read_table,
    styled_save_figure,
)

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


def build_selection_funnel(nodes: pd.DataFrame) -> pd.DataFrame:
    tier = nodes["candidate_tier"]
    rows = [
        ("All transition-network clusters", len(nodes)),
        (
            f"Size eligible (cluster size >= {MIN_CLUSTER_SIZE})",
            int(nodes["sse_tested"].fillna(False).sum()),
        ),
        (
            "Background or low information",
            int(tier.eq("background_or_low_information").sum()),
        ),
        ("Possible review", int(tier.eq("possible_review").sum())),
        ("High-priority local burst", int(tier.eq("high_priority_burst").sum())),
        ("High-priority onward burden", int(tier.eq("high_priority_burden").sum())),
        ("High priority on both axes", int(tier.eq("high_priority_both_axes").sum())),
    ]
    out = pd.DataFrame(rows, columns=["stage", "n"])
    out["order"] = np.arange(len(out))
    out["pct_all"] = out["n"] / len(nodes)
    return out


def build_score_landscape(nodes: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "cluster_id",
        "window_id",
        "window_idx",
        "wn_mid_date",
        "policy_era",
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


def draw_score_landscape(ax, table: pd.DataFrame) -> None:
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
                linewidth=0.6,  # Improved contrast outline
                alpha=0.9,  # Added subtle transparency for overlapping points
                label=TIER_DISPLAY_NAMES.get(tier, tier.title()),
                zorder=3,
            )

    # Draw separating line and text
    ax.axhline(line_y, color="#777777", lw=0.7, ls=":")
    ax.text(
        ax.get_xlim()[0], strip_y, "Burden N/A", va="center", ha="left", color="#555555"
    )

    ax.set_xlabel("Local burst calibrated score")
    ax.set_ylabel("Onward burden calibrated score")
    ax.legend(loc="center left")


def draw_selection_funnel(ax, table: pd.DataFrame) -> None:
    colors = [
        "#4C78A8",
        "#72A0C1",
        "#B8B8B8",
        "#E6AB02",
        "#D55E00",
        "#0072B2",
        "#7B3294",
    ]
    y = np.arange(len(table))[::-1]
    bars = ax.barh(y, table["n"], color=colors)
    ax.set_xscale("log")
    ax.set_yticks(y, [fill(stage, width=18) for stage in table["stage"]])
    ax.set_xlabel("Clusters (log scale)")
    for bar, value in zip(bars, table["n"]):
        label_x = max(value, 0.8) * 1.12 if value > 0 else 35
        ax.text(
            label_x,
            bar.get_y() + bar.get_height() / 2,
            f"{int(value):,}",
            va="center",
        )


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_score_landscape(nodes)
    funnel = build_selection_funnel(nodes).sort_values("order")
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")
    write_table(funnel, paths.result_table_dir, f"tab_{FILE_NAME}_selection_funnel")

    fig, axes = new_figure(
        ncols=2,
        width="double",
        height_in=4,
        constrained_layout=True,
        gridspec_kw={"width_ratios": [0.4, 0.6]},
    )
    draw_selection_funnel(axes[0], funnel)
    draw_score_landscape(axes[1], table)

    add_panel_labels(axes)
    outputs = styled_save_figure(fig, paths, f"fig_{FILE_NAME}")
    return {"figure": fig, "outputs": outputs}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    paths = paths_from_args(parser.parse_args())
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
