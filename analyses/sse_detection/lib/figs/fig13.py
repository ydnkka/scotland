"""Build fixed effects composition forest plot."""

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
from .forest import plot_composition_forest
from .forest_common import (
    MODEL_OUTCOMES,
    OUTLIER_LABELS,
    _finish_legend,
    _label_rows,
    _load_composition_tables,
)

FIGURE_NAME = "fig13_fixed_effects_composition"


def build(paths: Paths) -> dict[str, object]:
    """Create the fixed-effect forest plot for composition models."""
    composition_logistic, composition_linear = _load_composition_tables(paths)

    fig, axes = new_figure(
        nrows=3,
        ncols=2,
        width="double",
        height_in=8,
        constrained_layout=True,
    )

    plot_composition_forest(
        axes[0, 0],
        composition_logistic,
        panel=["demographic", "socioeconomic"],
        term_type="categorical_contrast",
        label_col="plot_label",
    )
    plot_composition_forest(
        axes[0, 1],
        composition_logistic,
        panel="geographic",
        term_type="categorical_contrast",
        label_col="plot_label",
        exclude=OUTLIER_LABELS,
    )
    plot_composition_forest(
        axes[1, 0],
        composition_linear,
        outcome="burst_score",
        panel=["demographic", "socioeconomic"],
        term_type="categorical_contrast",
        label_col="plot_label",
    )
    plot_composition_forest(
        axes[1, 1],
        composition_linear,
        outcome="burst_score",
        panel="geographic",
        term_type="categorical_contrast",
        label_col="plot_label",
    )
    plot_composition_forest(
        axes[2, 0],
        composition_linear,
        outcome="burden_score",
        panel=["demographic", "socioeconomic"],
        term_type="categorical_contrast",
        label_col="plot_label",
    )
    plot_composition_forest(
        axes[2, 1],
        composition_linear,
        outcome="burden_score",
        panel="geographic",
        term_type="categorical_contrast",
        label_col="plot_label",
        exclude=OUTLIER_LABELS,
    )

    flat_axes = axes.ravel()
    axes[0, 0].set_title("Sociodemographic")
    axes[0, 1].set_title("Geographic")
    _label_rows(axes[:, 0], tuple(label for _, label in MODEL_OUTCOMES))
    for ax in flat_axes:
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)

    add_panel_labels(flat_axes)
    _finish_legend(fig, flat_axes, ncol=2, y=-0.05)
    outputs = styled_save_figure(fig, paths, FIGURE_NAME)
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": pd.concat(
            [composition_logistic, composition_linear],
            ignore_index=True,
        ),
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
