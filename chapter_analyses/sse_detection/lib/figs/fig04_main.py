"""Build the main-text null-standardised Bayesian mixing forest figure."""

from __future__ import annotations

import argparse

import pandas as pd

from .fig04_app import OUTCOMES, _collect_rows, _sample_size
from .forest import (
    DEFAULT_COLORS,
    DEFAULT_MIXING_FEATURE_ORDER,
    DEFAULT_MODEL_LABELS,
    MIXING_FEATURE_LABELS,
    _add_model_legend,
    _draw_paired_forest_panel,
    _estimate_columns,
    _finish_forest_figure,
    _set_panel_xlim,
    _set_readable_coefficient_ticks,
    _set_readable_or_ticks,
    _y_lookup,
)
from .common import (
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    styled_new_figure,
    styled_save_figure,
)


FIGURE_NAME = "fig_ch5_bayesian_mixing_forest_main"
SCALE = "null"
SCALE_LABEL = "Null-standardised entropy"
SCALE_UNITS = "per 1-SD increase"


def build(paths: Paths) -> dict[str, object]:
    """Create a 3x1 forest plot restricted to null-standardised entropy."""
    collected, _, missing = _collect_rows(paths)
    plot_rows: list[pd.DataFrame] = []

    fig, axes = styled_new_figure(
        width="double",
        height_in=8.0,
        nrows=3,
        ncols=1,
        constrained_layout=True,
        sharey=True,
    )
    y_lookup = _y_lookup(DEFAULT_MIXING_FEATURE_ORDER)
    colors = dict(DEFAULT_COLORS)

    for row_idx, (family, outcome, outcome_label) in enumerate(OUTCOMES):
        ax = axes[row_idx]
        panel = pd.DataFrame(
            collected[(family, outcome)].loc[
            lambda frame: frame["scale"].eq(SCALE)
        ]
        )
        tagged = panel.copy()
        tagged["family"] = family
        tagged["outcome"] = outcome
        plot_rows.append(tagged)

        x_cols = _estimate_columns(family)
        reference = 1.0 if family == "logistic" else 0.0
        _draw_paired_forest_panel(
            ax,
            panel,
            y_lookup=y_lookup,
            row_order=DEFAULT_MIXING_FEATURE_ORDER,
            row_labels=MIXING_FEATURE_LABELS,
            model_labels=DEFAULT_MODEL_LABELS,
            colors=colors,
            x_cols=x_cols,
            reference=reference,
            point_size=4.6,
            interval_lw=1.25,
            dodge=0.12,
        )
        _set_panel_xlim(ax, panel, x_cols=x_cols, reference=reference)
        if family == "logistic":
            ax.set_xscale("log")
            _set_readable_or_ticks(ax)
            measure = "Odds ratio"
        else:
            _set_readable_coefficient_ticks(ax)
            measure = "Coefficient"
        ax.set_xlabel(f"{measure} {SCALE_UNITS}")
        if row_idx == 0:
            ax.set_title(SCALE_LABEL, pad=12)
        n = _sample_size(family, outcome, SCALE, paths.bayesian_result_dir)
        ax.text(
            0.98,
            0.95,
            f"n = {n:,}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            color="#444444",
        )
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)
        ax.tick_params(axis="y", length=0)
        ax.set_ylabel(outcome_label)
        panel_label(ax, chr(ord("A") + row_idx))

    _add_model_legend(fig, axes[0], colors, DEFAULT_MODEL_LABELS, 4.6)
    _finish_forest_figure(fig, axes, title=None)
    outputs = styled_save_figure(fig, paths, FIGURE_NAME, tight=False)
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": pd.concat(plot_rows, ignore_index=True),
        "missing_model_sets": tuple(missing),
    }


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
