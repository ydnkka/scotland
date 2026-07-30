"""Build the SSE score-distribution and null-calibration figure."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from ..sse.config import DETECTION_RANDOM_SEED
from ..sse.io import write_table
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
SUMMARY_NAME = "tab_ch5_null_calibration_summary"

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
INELIGIBLE_TIER = "size_ineligible"


def _dimension_columns(dimension: str) -> tuple[str, str]:
    return f"{dimension}_score", f"{dimension}_score_upper_p"


def _frames(
    table: pd.DataFrame,
    dimension: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split all nodes tested on an axis by axis-specific significance."""
    _, p_col = _dimension_columns(dimension)
    eligible = table["candidate_tier"].ne(INELIGIBLE_TIER)
    has_p = table[p_col].notna()
    tested = eligible & has_p
    axis_candidate = table["candidate_tier"].isin(
        DIMENSION_CANDIDATE_TIERS[dimension]
    )
    background = table.loc[tested & ~axis_candidate]
    candidates = table.loc[tested & axis_candidate]
    return background, candidates


def _legend_handles() -> list[Line2D]:
    return [
        Line2D([], [], color=CANDIDATE_COLOR, lw=5, alpha=0.65, label="Candidate"),
        Line2D([], [], color=BACKGROUND_COLOR, lw=2, label="Background"),
        Line2D(
            [],
            [],
            color="black",
            lw=0.8,
            linestyle="--",
            label="Uniform expectation",
        ),
    ]


def _uniform_ks_statistic(values: np.ndarray) -> float:
    """Return the one-sample Kolmogorov--Smirnov distance from U(0, 1)."""
    ordered = np.sort(np.asarray(values, dtype=float))
    n = ordered.size
    if n == 0:
        return np.nan
    upper = np.arange(1, n + 1, dtype=float) / n
    lower = np.arange(n, dtype=float) / n
    return float(max(np.max(upper - ordered), np.max(ordered - lower)))


def build_calibration_summary(
    table: pd.DataFrame,
    *,
    n_bins: int = 20,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Summarise conservative and randomized calibration for both axes."""
    rows: list[dict[str, object]] = []
    for axis in DIMENSIONS:
        base_col = f"{axis}_score_upper_p"
        tested = table["candidate_tier"].ne(INELIGIBLE_TIER)
        definitions = {
            "conservative": f"{base_col}_conservative",
            "randomized": f"{base_col}_randomized",
        }
        significant_sets: dict[str, set[object]] = {}
        for definition, p_col in definitions.items():
            values = table.loc[tested & table[p_col].notna(), p_col].to_numpy(float)
            counts = np.histogram(values, bins=n_bins, range=(0, 1))[0]
            expected = values.size / n_bins
            significant = tested & table[p_col].le(alpha)
            significant_sets[definition] = set(table.index[significant])
            rows.append(
                {
                    "axis": axis,
                    "p_value_definition": definition,
                    "n_tested": values.size,
                    "n_significant": int(significant.sum()),
                    "significant_rate": float(significant.sum() / values.size),
                    "n_equal_one": int(np.sum(values == 1)),
                    "ks_distance_uniform": _uniform_ks_statistic(values),
                    "mean_absolute_bin_deviation": float(
                        np.mean(np.abs(counts - expected))
                    ),
                    "max_absolute_bin_deviation": float(
                        np.max(np.abs(counts - expected))
                    ),
                    "expected_per_bin": expected,
                    "n_bins": n_bins,
                    "alpha": alpha,
                    "random_seed": (
                        DETECTION_RANDOM_SEED if definition == "randomized" else np.nan
                    ),
                }
            )

        conservative = significant_sets["conservative"]
        randomized = significant_sets["randomized"]
        union = conservative | randomized
        for row in rows[-2:]:
            row["significant_overlap_n"] = len(conservative & randomized)
            row["randomized_added_n"] = len(randomized - conservative)
            row["randomized_removed_n"] = len(conservative - randomized)
            row["significant_jaccard"] = (
                len(conservative & randomized) / len(union) if union else 1.0
            )
    return pd.DataFrame(rows)


def build(paths: Paths, *, n_bins: int = 20) -> dict[str, object]:
    """Create a 2x2 figure with null calibration and score distributions."""
    table = read_table(paths, "cluster_table")
    required = {"candidate_tier"}
    for dimension in DIMENSIONS:
        required.update(_dimension_columns(dimension))
    missing = sorted(required.difference(table.columns))
    if missing:
        raise KeyError(f"cluster_table is missing required columns: {missing}")

    summary_required = {
        f"{axis}_score_upper_p_{definition}"
        for axis in DIMENSIONS
        for definition in ("conservative", "randomized")
    }
    missing_summary = sorted(summary_required.difference(table.columns))
    if missing_summary:
        raise KeyError(
            f"cluster_table is missing calibration columns: {missing_summary}"
        )
    summary = build_calibration_summary(table, n_bins=n_bins)
    write_table(summary, paths.result_table_dir, SUMMARY_NAME)

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
            (len(p_background) + len(p_candidates)) / n_bins,
            color="black",
            linestyle="--",
            linewidth=0.8,
        )
        if row == 0:
            calibration_ax.set_title("Null-model calibration")
        calibration_ax.set_ylabel(f"{label}\nNumber of nodes")
        calibration_ax.tick_params(axis="x", labelbottom=True)
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
        "Randomized null-model $p$-value\n(probability of a score this high)"
    )
    axes[1, 1].set_xlabel("Detection score")

    for ax, label in zip(axes.ravel(), "ABCD"):
        panel_label(ax, label)
    fig.legend(
        handles=_legend_handles(),
        loc="outside lower center",
        bbox_to_anchor=(0.5, -0.075),
        ncol=3,
        frameon=False,
    )
    outputs = styled_save_figure(fig, paths, FIGURE_NAME, tight=False)
    return {"figure": fig, "outputs": outputs, "summary": summary}


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
