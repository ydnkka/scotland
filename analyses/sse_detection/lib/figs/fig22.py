"""Build size cumulative figure."""

from __future__ import annotations

import argparse

import numpy as np
from matplotlib.ticker import FuncFormatter

from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
)
from .size_profiles import COLORS, _finish, load_saved_summaries, prepare

FIGURE_NAME = "fig_size_cumulative"


def build(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    fig, axes = new_figure(width="double", nrows=2, ncols=2, height_in=5.4)
    columns = (
        "cumulative_clusters",
        "cluster_percent",
        "cumulative_records",
        "record_percent",
    )
    labels = (
        "Cumulative clusters",
        "Cumulative clusters (%)",
        "Cumulative sequence records",
        "Cumulative sequence records (%)",
    )
    for ax, col, label in zip(axes.flat, columns, labels):
        for group, block in tables["size_cumulative"].groupby("group"):
            ax.step(
                np.r_[5, block.cluster_size],
                np.r_[0, block[col]],
                where="post",
                color=COLORS[group],
                label=group,
            )
        ax.set_xscale("log")
        ax.set_xticks(
            [6, 10, 20, 50, 100, 500, 1700],
            ["6", "10", "20", "50", "100", "500", "1,700"],
        )
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.18)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
    axes[0, 0].legend()
    fig.supxlabel("Cluster size (number of sequence records)")
    add_panel_labels(axes)
    fig.tight_layout()
    return _finish(fig, paths, FIGURE_NAME)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument(
        "--from-saved-tables",
        action="store_true",
        help="Render saved size aggregates without rebuilding or writing them.",
    )
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths, tables=load_saved_summaries(paths) if args.from_saved_tables else None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
