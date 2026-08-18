"""Build the parameter sensitivity summary figure."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.ticker import PercentFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
    read_table,
    styled_save_figure,
)

from chapter_analyses.genomic_networks.lib.config import (
    ANALYSIS_RESOLUTION,
    SPARSIFICATION_THRESHOLD,
)

BASELINE_THRESHOLD = SPARSIFICATION_THRESHOLD
FIGURE_NAME = "fig_parameter_sensitivity"
LEIDEN_SUMMARY_TABLE = "leiden_resolution_sensitivity_summary"
SPARSIFICATION_SUMMARY_TABLE = "sparsification_threshold_sensitivity_summary"

LEIDEN_COLOR = "#35618f"
SPARSIFICATION_COLOR = "#b0473c"
REFERENCE_COLOR = "#555555"
BASELINE_EDGE_COLOR = "#222222"

FRAGMENTATION_RATIO_COLS = {
    "median": "median_ratio_clusters_per_1000_sequences_vs_baseline",
    "q25": "q25_ratio_clusters_per_1000_sequences_vs_baseline",
    "q75": "q75_ratio_clusters_per_1000_sequences_vs_baseline",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    return parser.parse_args()


def _line_with_iqr(
    ax: Axes,
    data: pd.DataFrame,
    *,
    x: str,
    median: str,
    q25: str | None = None,
    q75: str | None = None,
    color: str,
    label: str | None = None,
    linestyle: str = "-",
) -> None:
    work = data.sort_values(x)
    x_values = work[x].astype(float).to_numpy()
    median_values = work[median].astype(float).to_numpy()
    ax.plot(
        x_values,
        median_values,
        color=color,
        lw=1.4,
        ls=linestyle,
        label=label,
    )
    if q25 and q75 and q25 in work.columns and q75 in work.columns:
        ax.fill_between(
            x_values,
            work[q25].astype(float).to_numpy(),
            work[q75].astype(float).to_numpy(),
            color=color,
            alpha=0.16,
            lw=0,
        )


def _add_reference_line(ax: Axes, value: float, *, log_axis: bool = False) -> None:
    if log_axis and value <= 0:
        return
    ax.axvline(value, color=REFERENCE_COLOR, lw=0.8, ls=":")


def _format_resolution_axis(ax: Axes) -> None:
    ax.set_xlabel("Leiden resolution")
    ax.set_xlim(0.08, 0.82)


def _baseline_value(
    df: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    baseline: float,
) -> float:
    values = df[x_col].astype(float)
    baseline_rows = df.loc[np.isclose(values, baseline), y_col]
    if baseline_rows.empty:
        raise ValueError(f"Baseline value {baseline:g} not found in {x_col}.")
    return float(baseline_rows.iloc[0])


def _with_fragmentation_ratio(
    leiden: pd.DataFrame,
    *,
    baseline_resolution: float,
) -> pd.DataFrame:
    """Add baseline-relative fragmentation columns if the table predates them."""
    out = leiden.sort_values("resolution").copy()
    if all(col in out.columns for col in FRAGMENTATION_RATIO_COLS.values()):
        return out

    required = {
        "median_clusters_per_1000_sequences",
        "q25_clusters_per_1000_sequences",
        "q75_clusters_per_1000_sequences",
    }
    missing = required - set(out.columns)
    if missing:
        raise KeyError(
            "Leiden sensitivity summary lacks columns needed for fragmentation "
            f"ratios: {sorted(missing)}"
        )

    baseline_median = _baseline_value(
        out,
        x_col="resolution",
        y_col="median_clusters_per_1000_sequences",
        baseline=baseline_resolution,
    )
    for source, target in {
        "median_clusters_per_1000_sequences": FRAGMENTATION_RATIO_COLS["median"],
        "q25_clusters_per_1000_sequences": FRAGMENTATION_RATIO_COLS["q25"],
        "q75_clusters_per_1000_sequences": FRAGMENTATION_RATIO_COLS["q75"],
    }.items():
        out[target] = out[source] / baseline_median
    return out


def _plot_leiden_stability(
    ax: Axes,
    leiden: pd.DataFrame,
    *,
    baseline_resolution: float,
) -> None:
    _line_with_iqr(
        ax,
        leiden,
        x="resolution",
        median="median_ari_vs_baseline",
        q25="q25_ari_vs_baseline",
        q75="q75_ari_vs_baseline",
        color=LEIDEN_COLOR,
    )
    _add_reference_line(ax, baseline_resolution)
    _format_resolution_axis(ax)
    ax.set_title("Partition stability")
    ax.set_ylabel("ARI vs R=0.3")
    ax.set_ylim(-0.02, 1.04)


def _plot_leiden_fragmentation(
    ax: Axes,
    leiden: pd.DataFrame,
    *,
    baseline_resolution: float,
) -> None:
    work = _with_fragmentation_ratio(
        leiden,
        baseline_resolution=baseline_resolution,
    )
    _line_with_iqr(
        ax,
        work,
        x="resolution",
        median=FRAGMENTATION_RATIO_COLS["median"],
        q25=FRAGMENTATION_RATIO_COLS["q25"],
        q75=FRAGMENTATION_RATIO_COLS["q75"],
        color=LEIDEN_COLOR,
    )
    _add_reference_line(ax, baseline_resolution)
    ax.axhline(1.0, color=REFERENCE_COLOR, lw=0.8, ls=":")
    _format_resolution_axis(ax)
    ax.set_title("Cluster fragmentation")
    ax.set_ylabel("Clusters per 1,000 sequences\nrelative to R=0.3")


def _plot_sparsification_tradeoff(
    ax: Axes,
    sparsification: pd.DataFrame,
    *,
    baseline_threshold: float,
) -> None:
    required = {
        "threshold",
        "pooled_retained_edge_fraction",
        "pooled_retained_weight_fraction",
    }
    missing = required - set(sparsification.columns)
    if missing:
        raise KeyError(
            "Sparsification summary lacks columns needed for the trade-off panel: "
            f"{sorted(missing)}"
        )

    work = (
        sparsification.dropna(
            subset=[
                "threshold",
                "pooled_retained_edge_fraction",
                "pooled_retained_weight_fraction",
            ]
        )
        .sort_values("threshold")
        .copy()
    )
    x = work["pooled_retained_edge_fraction"].astype(float)
    y = work["pooled_retained_weight_fraction"].astype(float)
    ax.plot(
        x,
        y,
        color=SPARSIFICATION_COLOR,
        lw=1.4,
        marker="o",
        ms=3.5,
    )

    baseline_rows = work.loc[
        np.isclose(work["threshold"].astype(float), baseline_threshold)
    ]
    if not baseline_rows.empty:
        baseline = baseline_rows.iloc[0]
        baseline_x = float(baseline["pooled_retained_edge_fraction"])
        baseline_y = float(baseline["pooled_retained_weight_fraction"])
        ax.scatter(
            [baseline_x],
            [baseline_y],
            s=36,
            color=SPARSIFICATION_COLOR,
            edgecolor=BASELINE_EDGE_COLOR,
            linewidth=0.8,
            zorder=4,
        )
        ax.annotate(
            f"baseline\n{baseline_threshold:g}",
            xy=(baseline_x, baseline_y),
            xytext=(8, -14),
            textcoords="offset points",
            ha="left",
            va="top",
            arrowprops={
                "arrowstyle": "-",
                "color": REFERENCE_COLOR,
                "lw": 0.7,
            },
        )

    ax.xaxis.set_major_formatter(PercentFormatter(1.0))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_xlabel("Retained pairwise rows")
    ax.set_ylabel("Retained compatibility weight")
    ax.set_title("Sparsification trade-off")
    ax.set_xlim(max(0.0, float(x.min()) - 0.04), min(1.0, float(x.max()) + 0.04))
    ax.set_ylim(max(0.0, float(y.min()) - 0.04), min(1.02, float(y.max()) + 0.02))


def plot_parameter_sensitivity_grid(
    leiden_summary: pd.DataFrame,
    sparsification_summary: pd.DataFrame,
    *,
    paths: Paths,
    baseline_resolution: float = ANALYSIS_RESOLUTION,
    baseline_threshold: float = BASELINE_THRESHOLD,
) -> dict[str, Path]:
    """Plot a compact 1x3 parameter-sensitivity summary."""
    leiden = leiden_summary.sort_values("resolution").copy()
    sparsification = sparsification_summary.sort_values("threshold").copy()

    fig, axes = new_figure(
        width="double",
        height_in=2.75,
        nrows=1,
        ncols=3,
        constrained_layout=True,
    )
    axes = axes.ravel()

    _plot_leiden_stability(
        axes[0],
        leiden,
        baseline_resolution=baseline_resolution,
    )
    _plot_leiden_fragmentation(
        axes[1],
        leiden,
        baseline_resolution=baseline_resolution,
    )
    _plot_sparsification_tradeoff(
        axes[2],
        sparsification,
        baseline_threshold=baseline_threshold,
    )

    for ax in axes:
        ax.tick_params(axis="both", which="major", length=3)

    add_panel_labels(axes)
    return styled_save_figure(fig, paths, FIGURE_NAME)


def build(paths: Paths) -> dict[str, Path]:
    leiden = read_table(paths, LEIDEN_SUMMARY_TABLE)
    sparsification_summary = read_table(paths, SPARSIFICATION_SUMMARY_TABLE)
    return plot_parameter_sensitivity_grid(
        leiden,
        sparsification_summary,
        paths=paths,
    )


def main() -> int:
    args = parse_args()
    paths = paths_from_args(args)
    build(paths)
    print(f"Wrote {FIGURE_NAME} to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
