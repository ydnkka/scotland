"""Build random effects composition forest plot."""

from __future__ import annotations

import argparse

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
    ADJUSTED_TERM_TYPES,
    MODEL_OUTCOMES,
    OUTLIER_LABELS,
    _composition_label_order,
    _finish_legend,
    _load_composition_table,
)

FIGURE_NAME = "fig15_random_effects_composition"


def build(paths: Paths) -> dict[str, object]:
    """Create the forest plot for composition adjusters and random effects."""
    composition_linear = _load_composition_table(paths)
    label_order = _composition_label_order(composition_linear)

    fig, axes = new_figure(
        nrows=1,
        ncols=3,
        width="double",
        height_in=8,
        sharey=True,
    )

    for ax, (outcome, label) in zip(axes, MODEL_OUTCOMES):
        plot_composition_forest(
            ax,
            composition_linear,
            outcome=outcome,
            term_type=ADJUSTED_TERM_TYPES,
            exclude=OUTLIER_LABELS,
            label_order=label_order,
        )
        ax.set_title(label)
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)

    axes[0].set_xscale("log")
    axes[0].set_xlabel("Multiplicative odds (log scale)")
    add_panel_labels(axes)
    _finish_legend(fig, axes, ncol=2, y=0.015)
    outputs = styled_save_figure(fig, paths, FIGURE_NAME)
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": composition_linear.loc[
            composition_linear["plot_term_type"].isin(ADJUSTED_TERM_TYPES)
        ].copy(),
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
