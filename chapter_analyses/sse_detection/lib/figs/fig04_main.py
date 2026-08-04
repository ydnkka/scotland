"""Build fixed-effect Bayesian forest plots for SSE detection models."""

from __future__ import annotations

import argparse
from collections.abc import Callable

import pandas as pd

from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
    styled_save_figure,
)
from .forest import (
    plot_composition_forest,
    plot_mixing_forest,
)

FIGURE_NAME = {
    "mixing": "fig_ch5_fixed_effects_mixing",
    "composition": "fig_ch5_fixed_effects_composition",
}
OUTLIER_LABELS = ["Orkney", "Western Isles"]
MIXING_OUTCOMES = (
    ("candidate", "Candidate status"),
    ("burst_score", "Burst score"),
    ("burden_score", "Burden score"),
)
COMPOSITION_OUTCOMES = (
    ("candidate", "Candidate status"),
    ("burst_score", "Burst score"),
    ("burden_score", "Burden score"),
)


def _consolidated_table(paths: Paths, name: str) -> pd.DataFrame:
    table_dir = paths.bayesian_result_dir / "consolidated_tables"
    searched = []
    for suffix, reader in (("parquet", pd.read_parquet), ("csv", pd.read_csv)):
        path = table_dir / f"{name}.{suffix}"
        searched.append(path)
        if path.exists():
            return reader(path)
    locations = ", ".join(str(path) for path in searched)
    raise FileNotFoundError(f"Missing consolidated table {name!r}; searched: {locations}")


def _load_mixing_tables(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame]:
    return (
        _consolidated_table(paths, "mixing_logistic_consolidated_results"),
        _consolidated_table(paths, "mixing_linear_consolidated_results"),
    )


def _load_composition_tables(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame]:
    return (
        _consolidated_table(paths, "composition_logistic_consolidated_results"),
        _consolidated_table(paths, "composition_linear_consolidated_results"),
    )


def _finish_legend(fig, axes, *, ncol: int, y: float) -> None:
    handles, labels = axes[0].get_legend_handles_labels()
    for ax in axes:
        legend = ax.get_legend()
        if legend:
            legend.remove()
    fig.legend(
        handles=handles,
        labels=labels,
        loc="lower center",
        bbox_to_anchor=(0.5, y),
        ncol=ncol,
        frameon=False,
    )


def _label_rows(axes, labels: tuple[str, ...], *, x: float = -0.55) -> None:
    for ax, label in zip(axes, labels):
        ax.annotate(
            label,
            xy=(x, 0.5),
            xycoords="axes fraction",
            ha="center",
            va="center",
            rotation=90,
            fontweight="bold",
        )


def build_mixing(paths: Paths) -> dict[str, object]:
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

    for ax, (_, label) in zip(axes, MIXING_OUTCOMES):
        ax.set_title(label)
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)

    add_panel_labels(axes)
    _finish_legend(fig, axes, ncol=4, y=-0.1)
    outputs = styled_save_figure(fig, paths, FIGURE_NAME["mixing"])
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": pd.concat([mixing_logistic, mixing_linear], ignore_index=True),
    }


def build_composition(paths: Paths) -> dict[str, object]:
    """Create the fixed-effect forest plot for composition models."""
    composition_logistic, composition_linear = _load_composition_tables(paths)

    fig, axes = new_figure(
        nrows=3,
        ncols=2,
        width="double",
        height_in=7,
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
    _label_rows(axes[:, 0], tuple(label for _, label in COMPOSITION_OUTCOMES))
    for ax in flat_axes:
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)

    add_panel_labels(flat_axes)
    _finish_legend(fig, flat_axes, ncol=2, y=-0.05)
    outputs = styled_save_figure(fig, paths, FIGURE_NAME["composition"])
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": pd.concat(
            [composition_logistic, composition_linear],
            ignore_index=True,
        ),
    }


def build(paths: Paths) -> dict[str, object]:
    """Create both fixed-effect fig04 outputs."""
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
        help="Fixed-effect figure domain to build.",
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
