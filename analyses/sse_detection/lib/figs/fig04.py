"""Build fixed effects mixing forest plot."""

from __future__ import annotations

import argparse

import pandas as pd

from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
    styled_save_figure,
)
from .forest import plot_mixing_forest
from .forest_common import (
    MODEL_OUTCOMES,
    _finish_legend,
    _load_mixing_tables,
)

FIGURE_NAME = "fig04_fixed_effects_mixing"


def build(paths: Paths) -> dict[str, object]:
    """Create the fixed-effect forest plot for mixing models."""
    mixing_logistic, mixing_linear = _load_mixing_tables(paths)

    fig, axes = new_figure(
        nrows=1,
        ncols=3,
        width="double",
        height_in=3.5,
        sharey=True,
    )

    plot_mixing_forest(
        axes[0],
        mixing_logistic,
        plot_scale=["null_standardised", "observed"],
        term_type="mixing_entropy",
    )
    plot_mixing_forest(
        axes[1],
        mixing_linear,
        outcome="burst_score",
        plot_scale=["null_standardised", "observed"],
        term_type="mixing_entropy",
    )
    plot_mixing_forest(
        axes[2],
        mixing_linear,
        outcome="burden_score",
        plot_scale=["null_standardised", "observed"],
        term_type="mixing_entropy",
    )

    for ax, (_, label) in zip(axes, MODEL_OUTCOMES):
        ax.set_title(label)
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)

    add_panel_labels(axes)
    _finish_legend(fig, axes, ncol=4, y=-0.1)
    outputs = styled_save_figure(fig, paths, FIGURE_NAME)
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": pd.concat([mixing_logistic, mixing_linear], ignore_index=True),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
