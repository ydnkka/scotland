"""Build descriptive SSE composition differences and score associations."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..model.prep import COMPOSITION_SPECS
from .common import (
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    styled_new_figure,
    styled_save_figure,
)

FIGURE_NAME = "fig_ch5_sse_composition_descriptive"
POINT_COLOR = "#2F6690"
CI_COLOR = "#2F6690"
BOOTSTRAP_SEED = 42


@dataclass(frozen=True)
class Estimate:
    value: float
    low: float
    high: float


def _percentile_interval(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return tuple(np.nanpercentile(values, [2.5, 97.5], axis=0))


def _column_correlations(values: np.ndarray, outcome: np.ndarray) -> np.ndarray:
    """Pearson correlations between every values column and one outcome."""
    centered_values = values - values.mean(axis=0)
    centered_outcome = outcome - outcome.mean()
    numerator = np.sum(centered_values * centered_outcome[:, None], axis=0)
    denominator = np.sqrt(
        np.sum(centered_values**2, axis=0) * np.sum(centered_outcome**2)
    )
    return np.divide(
        numerator,
        denominator,
        out=np.full(values.shape[1], np.nan),
        where=denominator > 0,
    )


def _bootstrap_estimates(
    table: pd.DataFrame,
    level_columns: list[str],
    *,
    n_bootstrap: int,
    random_state: int,
) -> tuple[list[Estimate], list[Estimate], list[Estimate], dict[str, int]]:
    rng = np.random.default_rng(random_state)
    values = table[level_columns].to_numpy(dtype=float)
    candidate = table["sse_status"].eq("candidate").to_numpy()
    background = table["sse_status"].eq("background").to_numpy()
    candidate_values = values[candidate]
    background_values = values[background]

    observed_difference = (
        np.nanmean(candidate_values, axis=0)
        - np.nanmean(background_values, axis=0)
    ) * 100
    difference_draws = np.empty((n_bootstrap, len(level_columns)))
    for draw in range(n_bootstrap):
        sampled_candidate = candidate_values[
            rng.integers(0, len(candidate_values), len(candidate_values))
        ]
        sampled_background = background_values[
            rng.integers(0, len(background_values), len(background_values))
        ]
        difference_draws[draw] = (
            np.nanmean(sampled_candidate, axis=0)
            - np.nanmean(sampled_background, axis=0)
        ) * 100
    difference_low, difference_high = _percentile_interval(difference_draws)

    correlations: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, int]] = {}
    for score_column in ("burst_score", "burden_score"):
        valid = table[score_column].notna().to_numpy() & np.isfinite(values).all(axis=1)
        score_values = values[valid]
        outcome = table.loc[valid, score_column].to_numpy(dtype=float)
        observed = _column_correlations(score_values, outcome)
        draws = np.empty((n_bootstrap, len(level_columns)))
        for draw in range(n_bootstrap):
            sampled = rng.integers(0, len(outcome), len(outcome))
            draws[draw] = _column_correlations(
                score_values[sampled], outcome[sampled]
            )
        low, high = _percentile_interval(draws)
        correlations[score_column] = (observed, low, high, len(outcome))

    def pack(value: np.ndarray, low: np.ndarray, high: np.ndarray) -> list[Estimate]:
        return [Estimate(*map(float, triplet)) for triplet in zip(value, low, high)]

    burst = correlations["burst_score"]
    burden = correlations["burden_score"]
    sample_sizes = {
        "candidate": int(candidate.sum()),
        "background": int(background.sum()),
        "burst": burst[3],
        "burden": burden[3],
    }
    return (
        pack(observed_difference, difference_low, difference_high),
        pack(*burst[:3]),
        pack(*burden[:3]),
        sample_sizes,
    )


def _load_tables(composition_dir: Path) -> list[tuple[str, pd.DataFrame, list[str]]]:
    loaded = []
    metadata_columns = {"cluster_id", "sse_status", "burst_score", "burden_score"}
    for spec in COMPOSITION_SPECS:
        variable = str(spec["column"])
        path = composition_dir / f"cluster_composition_{variable}.parquet"
        if not path.exists():
            raise FileNotFoundError(f"Missing table: {path}")
        table = pd.read_parquet(path)
        missing = metadata_columns.difference(table.columns)
        if missing:
            raise KeyError(f"{path.name} is missing columns: {sorted(missing)}")
        level_columns = [col for col in table.columns if col not in metadata_columns]
        loaded.append((str(spec["label"]), table, level_columns))
    return loaded


def _draw_estimates(ax, estimates: list[Estimate], y: np.ndarray) -> None:
    values = np.asarray([estimate.value for estimate in estimates])
    low = np.asarray([estimate.low for estimate in estimates])
    high = np.asarray([estimate.high for estimate in estimates])
    ax.errorbar(
        values,
        y,
        xerr=np.vstack((values - low, high - values)),
        fmt="o",
        color=POINT_COLOR,
        ecolor=CI_COLOR,
        markersize=3.8,
        elinewidth=0.9,
        capsize=1.8,
        zorder=3,
    )


def build(
    paths: Paths,
    *,
    n_bootstrap: int = 1000,
    random_state: int = BOOTSTRAP_SEED,
) -> dict[str, object]:
    """Create the three-panel descriptive composition figure."""
    if n_bootstrap < 1:
        raise ValueError("n_bootstrap must be at least 1")
    tables = _load_tables(paths.result_table_dir)

    row_labels: list[str] = []
    differences: list[Estimate] = []
    burst_associations: list[Estimate] = []
    burden_associations: list[Estimate] = []
    group_boundaries: list[float | int] = []
    group_ranges: list[tuple[str, int, int]] = []
    sample_sizes: dict[str, int] | None = None

    for table_idx, (variable_label, table, level_columns) in enumerate(tables):
        group_start = len(row_labels)
        diff, burst, burden, sizes = _bootstrap_estimates(
            table,
            level_columns,
            n_bootstrap=n_bootstrap,
            random_state=random_state + table_idx,
        )
        row_labels.extend(str(level) for level in level_columns)
        group_ranges.append((variable_label, group_start, len(row_labels)))
        differences.extend(diff)
        burst_associations.extend(burst)
        burden_associations.extend(burden)
        if table_idx < len(tables) - 1:
            group_boundaries.append(len(row_labels) - 0.5)
        if sample_sizes is None:
            sample_sizes = sizes
        elif sample_sizes != sizes:
            raise ValueError("Composition tables do not have consistent sample sizes")

    assert sample_sizes is not None
    y = np.arange(len(row_labels))[::-1]
    fig, axes = styled_new_figure(
        width="double",
        height_in=max(8.0, 0.25 * len(row_labels) + 1.5),
        nrows=1,
        ncols=3,
        constrained_layout=True,
        sharey=True,
        gridspec_kw={"width_ratios": [1.15, 1, 1]},
    )

    panel_data = (differences, burst_associations, burden_associations)
    titles = (
        "Candidate − background",
        "Burst score",
        "Burden score",
    )
    xlabels = (
        "Mean difference (percentage points)",
        "Correlation ($r$)",
        "Correlation ($r$)",
    )
    subtitles = (
        f"n = {sample_sizes['candidate']:,} vs {sample_sizes['background']:,}",
        f"n = {sample_sizes['burst']:,}",
        f"n = {sample_sizes['burden']:,}",
    )

    for index, (ax, estimates, title, xlabel, subtitle) in enumerate(
        zip(axes, panel_data, titles, xlabels, subtitles)
    ):
        _draw_estimates(ax, estimates, y)
        ax.axvline(0, color="#555555", linestyle="--", linewidth=0.8, zorder=1)
        ax.set_title(f"{title}\n{subtitle}")
        ax.set_xlabel(xlabel)
        ax.set_ylim(-1, len(row_labels))
        ax.grid(axis="x", color="#E5E5E5", linewidth=0.6, zorder=0)
        for boundary in group_boundaries:
            ax.axhline(
                len(row_labels) - 1 - boundary,
                color="#D0D0D0",
                linewidth=0.7,
                zorder=0,
            )
        panel_label(ax, "ABC"[index])

    axes[0].set_yticks(y)
    axes[0].set_yticklabels(row_labels)
    axes[0].tick_params(axis="y", length=0)
    for variable_label, start, stop in group_ranges:
        top = len(row_labels) - 1 - start
        bottom = len(row_labels) - stop
        axes[0].text(
            -0.82,
            (top + bottom) / 2,
            variable_label,
            transform=axes[0].get_yaxis_transform(),
            rotation=90,
            ha="center",
            va="center",
            color="#333333",
            clip_on=False,
        )
    outputs = styled_save_figure(fig, paths, FIGURE_NAME, tight=False)
    return {"figure": fig, "outputs": outputs}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--random-state", type=int, default=BOOTSTRAP_SEED)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(
        paths,
        n_bootstrap=args.n_bootstrap,
        random_state=args.random_state,
    )
    print(f"Wrote {FIGURE_NAME} to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
