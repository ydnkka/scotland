"""Build Chapter 4 Supplementary Figure 1: parameter sensitivity summaries."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.ticker import PercentFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (  # noqa: E402
    Paths,
    add_common_args,
    add_panel_labels,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)

from chapter_analyses.genomic_networks.lib.config import (  # noqa: E402
    ANALYSIS_RESOLUTION,
    SPARSIFICATION_THRESHOLD,
)

BASELINE_THRESHOLD = SPARSIFICATION_THRESHOLD
FIGURE_NAME = "fig_ch4_parameter_sensitivity"
LEIDEN_SUMMARY_TABLE = "leiden_resolution_sensitivity_summary"
SPARSIFICATION_SUMMARY_TABLE = "sparsification_threshold_sensitivity_summary"
SPARSIFICATION_DETAIL_TABLE = "sparsification_threshold_sensitivity"

LEIDEN_COLOR = "#35618f"
SPARSIFICATION_COLOR = "#b0473c"
AMI_COLOR = "#7b5ea7"
BACKGROUND_COLOR = "#b8b8b8"
REFERENCE_COLOR = "#555555"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    return parser.parse_args()


def _positive_thresholds(df: pd.DataFrame) -> pd.DataFrame:
    """Drop threshold zero because the figure uses a logarithmic x-axis."""
    return df.loc[df["threshold"].astype(float) > 0].copy()


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


def _format_threshold_axis(ax: Axes) -> None:
    ax.set_xscale("log")
    ax.set_xlabel("Compatibility threshold")


def _add_reference_line(ax: Axes, value: float, *, log_axis: bool = False) -> None:
    if log_axis and value <= 0:
        return
    ax.axvline(value, color=REFERENCE_COLOR, lw=0.8, ls=":")


def _plot_leiden_fragmentation(
    ax: Axes,
    leiden: pd.DataFrame,
    *,
    baseline_resolution: float,
) -> None:
    _line_with_iqr(
        ax,
        leiden,
        x="resolution",
        median="median_clusters_per_1000_sequences",
        q25="q25_clusters_per_1000_sequences",
        q75="q75_clusters_per_1000_sequences",
        color=LEIDEN_COLOR,
    )
    _add_reference_line(ax, baseline_resolution)
    ax.set_title("Leiden sensitivity")
    ax.set_ylabel("Clusters per 1,000 sequences")


def _plot_sparsification_edge_fraction(
    ax: Axes,
    sparsification: pd.DataFrame,
    *,
    baseline_threshold: float,
) -> None:
    work = _positive_thresholds(sparsification)
    y_col = (
        "pooled_retained_weight_fraction"
        if "pooled_retained_weight_fraction" in work.columns
        else "pooled_retained_edge_fraction"
    )
    ylabel = (
        "Retained compatibility weight"
        if y_col == "pooled_retained_weight_fraction"
        else "Retained pairwise rows"
    )
    ax.plot(
        work["threshold"],
        work[y_col],
        color=SPARSIFICATION_COLOR,
        lw=1.4,
    )
    _add_reference_line(ax, baseline_threshold, log_axis=True)
    _format_threshold_axis(ax)
    ax.set_title("Sparsification sensitivity")
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_ylabel(ylabel)


def _plot_leiden_cluster_size(
    ax: Axes,
    leiden: pd.DataFrame,
    *,
    baseline_resolution: float,
) -> None:
    _line_with_iqr(
        ax,
        leiden,
        x="resolution",
        median="median_p90_cluster_size",
        q25="q25_p90_cluster_size",
        q75="q75_p90_cluster_size",
        color=LEIDEN_COLOR,
    )
    _add_reference_line(ax, baseline_resolution)
    ax.set_ylabel("90th percentile cluster size")


def _plot_sparsification_mean_degree(
    ax: Axes,
    sparsification: pd.DataFrame,
    *,
    baseline_threshold: float,
) -> None:
    work = _positive_thresholds(sparsification)
    _line_with_iqr(
        ax,
        work,
        x="threshold",
        median="median_retained_mean_degree",
        q25="q25_retained_mean_degree",
        q75="q75_retained_mean_degree",
        color=SPARSIFICATION_COLOR,
    )
    _add_reference_line(ax, baseline_threshold, log_axis=True)
    _format_threshold_axis(ax)
    ax.set_ylabel("Retained mean degree")


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
        label="ARI",
    )
    if "median_ami_vs_baseline" in leiden.columns:
        _line_with_iqr(
            ax,
            leiden,
            x="resolution",
            median="median_ami_vs_baseline",
            q25="q25_ami_vs_baseline",
            q75="q75_ami_vs_baseline",
            color=AMI_COLOR,
            label="AMI",
            linestyle="--",
        )
    _add_reference_line(ax, baseline_resolution)
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("Agreement with R=0.3")
    ax.set_ylim(-0.02, 1.04)
    ax.legend(loc="lower left")


def _sample_background_groups(
    detail: pd.DataFrame,
    *,
    max_groups: int,
    seed: int,
) -> pd.DataFrame:
    if max_groups <= 0:
        return detail.iloc[0:0].copy()
    stems = pd.Index(detail["pairwise_stem"].dropna().unique())
    if len(stems) > max_groups:
        stems = stems.to_series().sample(max_groups, random_state=seed).to_numpy()
    return detail.loc[detail["pairwise_stem"].isin(stems)].copy()


def _plot_sparsification_group_retention(
    ax: Axes,
    detail: pd.DataFrame,
    *,
    baseline_threshold: float,
    max_background_groups: int,
    background_seed: int,
) -> None:
    work = _positive_thresholds(detail)
    y_col = (
        "retained_weight_fraction"
        if "retained_weight_fraction" in work.columns
        else "retained_edge_fraction"
    )
    ylabel = (
        "Pairwise-group weight retention"
        if y_col == "retained_weight_fraction"
        else "Pairwise-group row retention"
    )
    background = _sample_background_groups(
        work,
        max_groups=max_background_groups,
        seed=background_seed,
    )
    for _, group in background.groupby("pairwise_stem", sort=False):
        group = group.sort_values("threshold")
        ax.plot(
            group["threshold"],
            group[y_col],
            color=BACKGROUND_COLOR,
            lw=0.35,
            alpha=0.12,
            zorder=1,
        )

    summary = (
        work.groupby("threshold", dropna=False)[y_col]
        .agg(
            median="median",
            q25=lambda x: x.quantile(0.25),
            q75=lambda x: x.quantile(0.75),
        )
        .reset_index()
    )
    _line_with_iqr(
        ax,
        summary,
        x="threshold",
        median="median",
        q25="q25",
        q75="q75",
        color=SPARSIFICATION_COLOR,
        label="Median and IQR",
    )
    _add_reference_line(ax, baseline_threshold, log_axis=True)
    _format_threshold_axis(ax)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_ylabel(ylabel)
    ax.legend(loc="lower left")


def plot_parameter_sensitivity_grid(
    leiden_summary: pd.DataFrame,
    sparsification_summary: pd.DataFrame,
    sparsification_detail: pd.DataFrame,
    *,
    paths: Paths,
    baseline_resolution: float = ANALYSIS_RESOLUTION,
    baseline_threshold: float = BASELINE_THRESHOLD,
    max_background_groups: int = 300,
    background_seed: int = 42,
) -> dict[str, Path]:
    """Plot the 3x2 Leiden and sparsification sensitivity grid."""
    leiden = leiden_summary.sort_values("resolution").copy()
    sparsification = sparsification_summary.sort_values("threshold").copy()
    detail = sparsification_detail.sort_values(["threshold", "pairwise_stem"]).copy()

    fig, axes = styled_new_figure(
        width="double",
        height_in=7.0,
        nrows=3,
        ncols=2,
        sharex=False,
        constrained_layout=True,
    )

    _plot_leiden_fragmentation(
        axes[0, 0],
        leiden,
        baseline_resolution=baseline_resolution,
    )
    _plot_sparsification_edge_fraction(
        axes[0, 1],
        sparsification,
        baseline_threshold=baseline_threshold,
    )
    _plot_leiden_cluster_size(
        axes[1, 0],
        leiden,
        baseline_resolution=baseline_resolution,
    )
    _plot_sparsification_mean_degree(
        axes[1, 1],
        sparsification,
        baseline_threshold=baseline_threshold,
    )
    _plot_leiden_stability(
        axes[2, 0],
        leiden,
        baseline_resolution=baseline_resolution,
    )
    _plot_sparsification_group_retention(
        axes[2, 1],
        detail,
        baseline_threshold=baseline_threshold,
        max_background_groups=max_background_groups,
        background_seed=background_seed,
    )

    for ax in axes[:, 0]:
        ax.set_xlim(0.08, 0.82)
    for ax in axes.flat:
        ax.tick_params(axis="both", which="major", length=3)

    add_panel_labels(axes.flat, x=-0.12, y=1.08, size="medium")
    return styled_save_figure(fig, paths, FIGURE_NAME, tight=False)


def build(paths: Paths) -> dict[str, Path]:
    leiden = read_table(paths, LEIDEN_SUMMARY_TABLE)
    sparsification_summary = read_table(paths, SPARSIFICATION_SUMMARY_TABLE)
    sparsification_detail = read_table(paths, SPARSIFICATION_DETAIL_TABLE)
    return plot_parameter_sensitivity_grid(
        leiden,
        sparsification_summary,
        sparsification_detail,
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
