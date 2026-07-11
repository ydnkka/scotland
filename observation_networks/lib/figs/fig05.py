"""Build Chapter 4 Supplementary Figure 6: within-cluster pairwise distances."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

from matplotlib.axes import Axes
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    add_policy_bands,
    date_axis,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)


FIGURE_NAME = "fig_ch4_cluster_pairwise_distances"

MIN_PAIRWISE_ROWS = 10
SNP_COLOR = "#35618f"
TEMPORAL_COLOR = "#b0473c"


def _numeric(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce")


def _weighted_quantile(
    values: pd.Series,
    weights: pd.Series,
    quantile: float,
) -> float:
    work = pd.DataFrame(
        {
            "value": _numeric(values),
            "weight": _numeric(weights),
        }
    ).replace([np.inf, -np.inf], np.nan)
    work = work.dropna()
    work = work.loc[work["weight"].gt(0)]
    if work.empty:
        return np.nan

    work = work.sort_values("value", kind="mergesort")
    values_array = work["value"].to_numpy(dtype=float)
    weights_array = work["weight"].to_numpy(dtype=float)
    total_weight = float(weights_array.sum())
    if total_weight <= 0:
        return np.nan

    positions = (np.cumsum(weights_array) - 0.5 * weights_array) / total_weight
    return float(np.interp(quantile, positions, values_array))


def _window_weighted_quantiles(summary: pd.DataFrame, value_col: str) -> pd.DataFrame:
    rows = []
    for window_idx, group in summary.groupby("window_idx", dropna=False):
        rows.append(
            {
                "window_idx": window_idx,
                "q25": _weighted_quantile(
                    group[value_col], group["n_pairwise_rows"], 0.25
                ),
                "median": _weighted_quantile(
                    group[value_col], group["n_pairwise_rows"], 0.50
                ),
                "q75": _weighted_quantile(
                    group[value_col], group["n_pairwise_rows"], 0.75
                ),
                "n_window_lineages": group["pango_lineage"].nunique(),
                "n_pairwise_rows": _numeric(group["n_pairwise_rows"]).sum(),
            }
        )
    return pd.DataFrame(rows)


def _prepare_inputs(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = read_table(paths, "cluster_pairwise_distance_summary")
    window = read_table(paths, "window_coverage")
    window["wn_mid_date"] = pd.to_datetime(window["wn_mid_date"], errors="coerce")

    summary = summary.loc[
        summary["status"].eq("ok")
        & _numeric(summary["n_pairwise_rows"]).ge(MIN_PAIRWISE_ROWS)
    ].copy()
    for col in (
        "window_idx",
        "n_pairwise_rows",
        "snp_distance_median",
        "temporal_distance_median",
    ):
        summary[col] = _numeric(summary[col])

    summary = summary.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[
            "window_idx",
            "n_pairwise_rows",
            "snp_distance_median",
            "temporal_distance_median",
        ]
    )
    summary = summary.loc[summary["n_pairwise_rows"].gt(0)]
    if summary.empty:
        raise ValueError(
            "No supported pairwise-distance rows remain after filtering "
            f"status == 'ok' and n_pairwise_rows >= {MIN_PAIRWISE_ROWS}."
        )

    return summary, window


def _plot_window_distance(
    ax: Axes,
    window_summary: pd.DataFrame,
    window_coverage: pd.DataFrame,
    *,
    color: str,
    title: str,
    ylabel: str,
    panel: str,
) -> None:
    work = window_summary.merge(
        window_coverage[["window_idx", "wn_mid_date", "policy_period"]],
        on="window_idx",
        how="left",
    ).sort_values("wn_mid_date")
    add_policy_bands(ax, window_coverage)
    x = pd.to_datetime(work["wn_mid_date"], errors="coerce")
    median = _numeric(work["median"]).to_numpy(dtype=float)
    q25 = _numeric(work["q25"]).to_numpy(dtype=float)
    q75 = _numeric(work["q75"]).to_numpy(dtype=float)
    mask = (
        x.notna().to_numpy()
        & np.isfinite(median)
        & np.isfinite(q25)
        & np.isfinite(q75)
    )
    ax.fill_between(
        x.loc[mask],
        q25[mask],
        q75[mask],
        color=color,
        alpha=0.18,
        linewidth=0,
    )
    ax.plot(x.loc[mask], median[mask], color=color, lw=1.35)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    date_axis(ax)
    panel_label(ax, panel)


def _plot_distance_hexbin(ax: Axes, summary: pd.DataFrame, fig) -> None:
    work = summary.loc[
        summary["snp_distance_median"].notna()
        & summary["temporal_distance_median"].notna()
    ]
    hb = ax.hexbin(
        work["snp_distance_median"],
        work["temporal_distance_median"],
        gridsize=36,
        mincnt=1,
        bins="log",
        cmap="viridis",
    )
    fig.colorbar(hb, ax=ax, label="Window-lineage summaries")
    ax.set_title("Within-cluster SNP and temporal distance")
    ax.set_xlabel("Median within-cluster SNP distance")
    ax.set_ylabel("Median days")
    panel_label(ax, "C")


def build(paths: Paths) -> None:
    summary, window = _prepare_inputs(paths)
    snp_window = _window_weighted_quantiles(summary, "snp_distance_median")
    temporal_window = _window_weighted_quantiles(
        summary, "temporal_distance_median"
    )

    fig, axes = styled_new_figure(
        width="double",
        height_in=7.2,
        nrows=3,
        ncols=1,
    )
    _plot_window_distance(
        axes[0],
        snp_window,
        window,
        color=SNP_COLOR,
        title="Within-cluster SNP distance among non-singleton clusters",
        ylabel="Median SNP distance",
        panel="A",
    )
    _plot_window_distance(
        axes[1],
        temporal_window,
        window,
        color=TEMPORAL_COLOR,
        title="Within-cluster temporal distance among non-singleton clusters",
        ylabel="Median days",
        panel="B",
    )
    _plot_distance_hexbin(axes[2], summary, fig)
    fig.subplots_adjust(left=0.11, right=0.96, top=0.95, bottom=0.09, hspace=0.52)
    styled_save_figure(fig, paths, FIGURE_NAME, tight=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    print(f"Wrote {FIGURE_NAME} to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
