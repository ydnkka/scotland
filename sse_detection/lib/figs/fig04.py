"""Build the combined Bayesian mixing-model forest figure."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .forest import (
    DEFAULT_COLORS,
    DEFAULT_MIXING_FEATURE_ORDER,
    DEFAULT_MODEL_LABELS,
    MIXING_FEATURE_LABELS,
    _add_model_legend,
    _collect_mixing_forest_rows,
    _draw_paired_forest_panel,
    _estimate_columns,
    _finish_forest_figure,
    _set_panel_xlim,
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


FIGURE_NAME = "fig_ch5_bayesian_mixing_forest"
OUTCOMES = (
    ("logistic", "candidate", "Candidate status"),
    ("linear", "burst_score", "Burst score"),
    ("linear", "burden_score", "Burden score"),
)
SCALES = (
    ("observed", "Observed entropy", "per 0.1-bit increase"),
    ("null", "Null-standardised entropy", "per 1-SD increase"),
)


def _sample_size(family: str, outcome: str, scale: str, result_dir: Path) -> int:
    path = result_dir / family / "mixing_fit_frame_summary.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing table: {path}")
    summary = pd.read_csv(path)
    rows = summary.loc[
        summary["outcome"].eq(outcome)
        & summary["model_set"].eq(f"{scale}_expanded")
    ]
    if rows.empty:
        raise ValueError(f"No sample-size row for {outcome}/{scale} in {path}")
    return int(rows.iloc[0]["fit_rows"])


def build(
    paths: Paths,
) -> dict[str, object]:
    """Create a 3x2 forest plot for all fitted Bayesian mixing outcomes."""
    collected: dict[tuple[str, str], pd.DataFrame] = {}
    plot_rows: list[pd.DataFrame] = []
    missing: list[str] = []
    for family, outcome, _ in OUTCOMES:
        rows, absent = _collect_mixing_forest_rows(
            paths.bayesian_result_dir / family,
            family=family,
            outcome=outcome,
            feature_order=DEFAULT_MIXING_FEATURE_ORDER,
        )
        if rows.empty:
            raise FileNotFoundError(f"No mixing summaries found for {outcome}")
        collected[(family, outcome)] = rows
        tagged = rows.copy()
        tagged["family"] = family
        tagged["outcome"] = outcome
        plot_rows.append(tagged)
        missing.extend(f"{outcome}/{model_set}" for model_set in absent)

    fig, placeholder_ax = styled_new_figure(
        width="double",
        height_in=8.0,
        constrained_layout=True,
    )
    placeholder_ax.remove()
    axes = fig.subplots(3, 2, sharey=True)
    y_lookup = _y_lookup(DEFAULT_MIXING_FEATURE_ORDER)
    colors = dict(DEFAULT_COLORS)

    for row_idx, (family, outcome, outcome_label) in enumerate(OUTCOMES):
        rows = collected[(family, outcome)]
        for col_idx, (scale, scale_label, scale_units) in enumerate(SCALES):
            ax = axes[row_idx, col_idx]
            panel = rows.loc[rows["scale"].eq(scale)].copy()
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
                measure = "Coefficient"
            ax.set_xlabel(f"{measure} {scale_units}")
            if row_idx == 0:
                ax.set_title(scale_label, pad=12)
            n = _sample_size(family, outcome, scale, paths.bayesian_result_dir)
            ax.text(
                0.98,
                0.95,
                f"n = {n:,}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                # fontsize=8,
                color="#444444",
            )
            ax.grid(axis="x", color="#E6E6E6", lw=0.6)
            ax.tick_params(axis="y", length=0)
            if col_idx == 0:
                ax.set_ylabel(outcome_label, fontweight="bold")
            else:
                ax.tick_params(labelleft=False)
            panel_label(ax, chr(ord("A") + row_idx * 2 + col_idx))

    _add_model_legend(fig, axes[0, 0], colors, DEFAULT_MODEL_LABELS, 4.6)
    _finish_forest_figure(fig, np.asarray(axes), title=None)
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
