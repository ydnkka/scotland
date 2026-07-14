"""Build the SSE score-distribution and null-calibration figure."""

from __future__ import annotations

import argparse

from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from .common import (
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)


FIGURE_NAME = "fig_ch5_sse_score_null_calibration"

CANDIDATE_COLOR = "#C44E52"
BACKGROUND_COLOR = "#B0B0B0"
DIMENSIONS = ("burst", "burden")
DIMENSION_LABELS = {
    "burst": "Local burst",
    "burden": "Onward burden",
}
DIMENSION_CANDIDATE_TIERS = {
    "burst": ("high_priority_burst", "high_priority_both_axes"),
    "burden": ("high_priority_burden", "high_priority_both_axes"),
}
HIGH_PRIORITY_TIERS = (
    "high_priority_both_axes",
    "high_priority_burst",
    "high_priority_burden",
)
INELIGIBLE_TIER = "size_ineligible"


def _dimension_columns(dimension: str) -> tuple[str, str]:
    return f"{dimension}_score", f"{dimension}_score_upper_p"


def _frames(
    table: pd.DataFrame,
    dimension: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return eligible background and axis-specific candidate nodes."""
    _, p_col = _dimension_columns(dimension)
    eligible = table["candidate_tier"].ne(INELIGIBLE_TIER)
    has_p = table[p_col].notna()
    background = table.loc[
        eligible & has_p & ~table["candidate_tier"].isin(HIGH_PRIORITY_TIERS)
    ]
    candidates = table.loc[
        eligible
        & has_p
        & table["candidate_tier"].isin(DIMENSION_CANDIDATE_TIERS[dimension])
    ]
    return background, candidates


def _legend_handles() -> list[Line2D]:
    return [
        Line2D([], [], color=CANDIDATE_COLOR, lw=5, alpha=0.65, label="Candidate"),
        Line2D([], [], color=BACKGROUND_COLOR, lw=2, label="Background"),
    ]


def build(paths: Paths, *, n_bins: int = 20) -> dict[str, object]:
    """Create a 2x2 figure with null calibration and score distributions."""
    table = read_table(paths, "cluster_table")
    required = {"candidate_tier"}
    for dimension in DIMENSIONS:
        required.update(_dimension_columns(dimension))
    missing = sorted(required.difference(table.columns))
    if missing:
        raise KeyError(f"cluster_table is missing required columns: {missing}")

    fig, axes = styled_new_figure(
        nrows=2,
        ncols=2,
        width="double",
        height_in=6.0,
        sharex=True,
        constrained_layout=True,
    )

    for row, dimension in enumerate(DIMENSIONS):
        score_col, p_col = _dimension_columns(dimension)
        label = DIMENSION_LABELS[dimension]
        background, candidates = _frames(table, dimension)
        p_background = background[p_col].dropna().to_numpy()
        p_candidates = candidates[p_col].dropna().to_numpy()

        calibration_ax = axes[row, 0]
        calibration_ax.hist(
            [p_background, p_candidates],
            bins=n_bins,
            range=(0, 1),
            stacked=True,
            color=[BACKGROUND_COLOR, CANDIDATE_COLOR],
            edgecolor="white",
            linewidth=0.4,
        )
        calibration_ax.axhline(
            len(p_background) / n_bins,
            color="black",
            linestyle="--",
            linewidth=0.8,
            label="Background uniform expectation",
        )
        if row == 0:
            calibration_ax.set_title("Null-model calibration")
        calibration_ax.set_ylabel(f"{label}\nNumber of nodes")
        calibration_ax.tick_params(axis="x", labelbottom=True)
        if row == 0:
            calibration_ax.legend(loc="upper center", frameon=False)

        distribution_ax = axes[row, 1]
        score_arrays = [
            background[score_col].dropna().to_numpy(),
            candidates[score_col].dropna().to_numpy(),
        ]
        finite_scores = np.concatenate(
            [values[np.isfinite(values)] for values in score_arrays if values.size]
        )
        if finite_scores.size:
            score_range = (float(finite_scores.min()), float(finite_scores.max()))
            for values, color, filled in zip(
                score_arrays,
                (BACKGROUND_COLOR, CANDIDATE_COLOR),
                (False, True),
            ):
                if values.size:
                    distribution_ax.hist(
                        values,
                        bins=n_bins,
                        range=score_range,
                        density=True,
                        histtype="stepfilled" if filled else "step",
                        color=color,
                        edgecolor=color,
                        linewidth=1.0,
                        alpha=0.55 if filled else 1.0,
                    )
        if row == 0:
            distribution_ax.set_title("Score distribution")
        distribution_ax.set_ylabel("Density")
        distribution_ax.tick_params(axis="x", labelbottom=True)

    axes[1, 0].set_xlabel(
        "Null-model $p$-value (probability of a score this high)"
    )
    axes[1, 1].set_xlabel("Detection score")

    for ax, label in zip(axes.ravel(), "ABCD"):
        panel_label(ax, label)
    fig.legend(
        handles=_legend_handles(),
        loc="outside lower center",
        ncol=2,
        frameon=False,
    )
    outputs = styled_save_figure(fig, paths, FIGURE_NAME, tight=False)
    return {"figure": fig, "outputs": outputs}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--n-bins", type=int, default=20)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths, n_bins=args.n_bins)
    print(f"Wrote {FIGURE_NAME} to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
