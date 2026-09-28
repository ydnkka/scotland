"""Build random effects mixing forest plot."""

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
from .forest import plot_mixing_forest
from .forest_common import (
    ADJUSTED_TERM_TYPES,
    MIXING_SCALES,
    MODEL_OUTCOMES,
    _finish_legend,
    _load_mixing_table,
    _mixing_label_order,
)

FIGURE_NAME = "fig_random_effects_mixing"


def build(paths: Paths) -> dict[str, object]:
    """Create the forest plot for mixing adjusters and random effects."""
    mixing_linear = _load_mixing_table(paths)
    label_order = _mixing_label_order(mixing_linear)

    fig, axes = new_figure(
        nrows=1,
        ncols=3,
        width="double",
        height_in=8,
        sharey=True,
    )

    for ax, (outcome, label) in zip(axes, MODEL_OUTCOMES):
        plot_mixing_forest(
            ax,
            mixing_linear,
            outcome=outcome,
            plot_scale=MIXING_SCALES,
            term_type=ADJUSTED_TERM_TYPES,
            label_order=label_order,
        )
        ax.set_title(label)
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)

    axes[0].set_xlabel("Multiplicative odds")
    add_panel_labels(axes)
    _finish_legend(fig, axes, ncol=4, y=0.01)
    outputs = styled_save_figure(fig, paths, FIGURE_NAME)
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": mixing_linear.loc[
            mixing_linear["plot_term_type"].isin(ADJUSTED_TERM_TYPES)
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
