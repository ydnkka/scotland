"""Build appendix Bayesian forest plots for random terms and adjusters."""

from __future__ import annotations

import argparse
from collections.abc import Callable

import pandas as pd

from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    paths_from_args,
    styled_new_figure,
    styled_save_figure,
)
from .fig04_main import (
    OUTLIER_LABELS,
    _consolidated_table,
    _finish_legend,
)
from .forest import (
    plot_composition_forest,
    plot_mixing_forest,
)

FIGURE_NAME = {
    "mixing": "fig_ch5_random_effects_mixing",
    "composition": "fig_ch5_random_effects_composition",
}
APP_OUTCOMES = (
    ("burst_score", "Burst score"),
    ("burden_score", "Burden score"),
)
APP_TERM_TYPES = ["continuous_adjuster", "random_intercept"]
MIXING_SCALES = ["null_standardised", "observed"]


def _load_mixing_table(paths: Paths) -> pd.DataFrame:
    return _consolidated_table(paths, "mixing_linear_consolidated_results")


def _load_composition_table(paths: Paths) -> pd.DataFrame:
    return _consolidated_table(paths, "composition_linear_consolidated_results")


def build_mixing(paths: Paths) -> dict[str, object]:
    """Create the appendix forest plot for mixing adjusters and random effects."""
    mixing_linear = _load_mixing_table(paths)

    fig, axes = styled_new_figure(
        nrows=1,
        ncols=2,
        width="double",
        height_in=8,
        sharey=True,
    )
    axes = list(axes)

    for ax, (outcome, label) in zip(axes, APP_OUTCOMES):
        plot_mixing_forest(
            ax,
            mixing_linear,
            outcome=outcome,
            plot_scale=MIXING_SCALES,
            term_type=APP_TERM_TYPES,
        )
        ax.set_title(label)
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)

    add_panel_labels(axes, x=-0.08, y=1.05)
    _finish_legend(fig, axes, ncol=4, y=0.01)
    outputs = styled_save_figure(fig, paths, FIGURE_NAME["mixing"], tight=False)
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": mixing_linear.loc[
            mixing_linear["plot_term_type"].isin(APP_TERM_TYPES)
        ].copy(),
    }


def build_composition(paths: Paths) -> dict[str, object]:
    """Create the appendix forest plot for composition adjusters and random effects."""
    composition_linear = _load_composition_table(paths)

    fig, axes = styled_new_figure(
        nrows=1,
        ncols=2,
        width="double",
        height_in=8,
        sharey=True,
    )
    axes = list(axes)

    for ax, (outcome, label) in zip(axes, APP_OUTCOMES):
        plot_composition_forest(
            ax,
            composition_linear,
            outcome=outcome,
            term_type=APP_TERM_TYPES,
            exclude=OUTLIER_LABELS,
        )
        ax.set_title(label)
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)

    add_panel_labels(axes, x=-0.08, y=1.05)
    _finish_legend(fig, axes, ncol=2, y=0.015)
    outputs = styled_save_figure(fig, paths, FIGURE_NAME["composition"], tight=False)
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": composition_linear.loc[
            composition_linear["plot_term_type"].isin(APP_TERM_TYPES)
        ].copy(),
    }


def build(paths: Paths) -> dict[str, object]:
    """Create both appendix fig04 outputs."""
    return {
        "mixing": build_mixing(paths),
        "composition": build_composition(paths),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument(
        "--domain",
        choices=("all", "mixing", "composition"),
        default="all",
        help="Appendix figure domain to build.",
    )
    args = parser.parse_args()
    paths = paths_from_args(args)
    builders: dict[str, Callable[[Paths], dict[str, object]]] = {
        "mixing": build_mixing,
        "composition": build_composition,
    }
    selected = builders if args.domain == "all" else {args.domain: builders[args.domain]}
    for domain, builder in selected.items():
        builder(paths)
        print(f"Wrote {FIGURE_NAME[domain]} to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
