"""Build Chapter 4 Supplementary Figure 5: assortativity variance decomposition."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))

from assortativity_analysis import (
    VARIANCE_REFERENCE_WINDSORISE,
    compatibility_variance_decomposition_long,
    variance_decomposition_summary,
)
from common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
    styled_save_figure,
)

from chapter_analyses.genomic_networks.lib.io import write_table

FIGURE_NAME = "fig_ch4_assortativity_variance_decomposition"


def plot_variance_decomposition(paths: Paths, vd_long: pd.DataFrame) -> None:
    order = (
        vd_long.loc[vd_long["winsorize"].eq(95)]
        .sort_values("additive_model_fraction", ascending=False)["attribute_label"]
        .drop_duplicates()
        .tolist()
    )

    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#1f77b4"])
    palette = {label: cycle[idx % len(cycle)] for idx, label in enumerate(order)}

    fig, axes = new_figure(
        width="double",
        height_in=3.5,
        nrows=1,
        ncols=3,
        constrained_layout=True,
    )
    axes = axes.ravel()

    panels = [
        ("additive_model_fraction", "Additive model", (0, 1)),
        ("window_given_lineage_fraction", "Window | Lineage", (0, 0.5)),
        ("lineage_given_window_fraction", "Lineage | Window", (0, 0.5)),
    ]

    for ax, (col, title, ylim) in zip(axes, panels):
        for label in order:
            group = vd_long.loc[vd_long["attribute_label"].eq(label)].sort_values(
                "winsorize"
            )
            x = group["winsorize"].to_numpy()
            ax.plot(
                x,
                group[col].to_numpy(),
                marker="o",
                color=palette[label],
                label=label,
            )
            if col == "additive_model_fraction":
                ax.plot(
                    x,
                    group["adj_additive_model_fraction"].to_numpy(),
                    color=palette[label],
                    linestyle="--",
                    alpha=0.5,
                    label="_nolegend_",
                )

        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_ylim(*ylim)
        ax.axvline(90, ls="--", color="black", alpha=0.6)
        ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))

    handles, labels = axes[1].get_legend_handles_labels()
    handles.append(Line2D([], [], ls="--", color="black", alpha=0.6))
    labels.append("Reported cap (90)")

    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        bbox_to_anchor=(0.5, -0.15),
    )

    axes[0].text(
        0.03,
        0.05,
        "dashed = adjusted",
        transform=axes[0].transAxes,
        fontsize="small",
        color="grey",
    )

    fig.supxlabel("Inverse-variance weight cap percentile")
    fig.supylabel("Fraction of variance explained")
    add_panel_labels(axes)
    styled_save_figure(fig, paths, FIGURE_NAME)


def build(paths: Paths) -> pd.DataFrame:
    vd_long = compatibility_variance_decomposition_long(paths)
    summary = variance_decomposition_summary(
        vd_long,
        winsorize=VARIANCE_REFERENCE_WINDSORISE,
    )
    write_table(
        vd_long,
        "compatibility_variance_decomposition",
        table_dir=paths.table_dir,
    )
    write_table(
        summary,
        "compatibility_variance_decomposition_summary",
        table_dir=paths.table_dir,
    )
    plot_variance_decomposition(paths, vd_long)
    return summary


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
