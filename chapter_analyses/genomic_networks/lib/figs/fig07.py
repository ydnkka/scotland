"""Build Chapter 4 Figure 3: pooled compatibility assortativity over time."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.dates as mdates
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from assortativity_analysis import (
    compatibility_window_pooled_meta,
    pooled_window_attribute_summary,
)
from common import (
    Paths,
    add_common_args,
    add_panel_labels,
    date_axis,
    new_figure,
    paths_from_args,
    styled_save_figure,
)

from chapter_analyses.genomic_networks.lib.io import write_table

FIGURE_NAME = "fig_ch4_assortativity_pooled_window"


def _date_values(values: pd.Series) -> np.ndarray:
    dates = pd.to_datetime(values, errors="coerce")
    out = np.full(len(dates), np.nan, dtype=float)
    valid = dates.notna().to_numpy()
    if valid.any():
        out[valid] = mdates.date2num(dates.loc[valid].dt.to_pydatetime())
    return out


def plot_pooled_window_meta(
    paths: Paths,
    meta: pd.DataFrame,
    window_lookup: pd.DataFrame,
    *,
    exclude_attrs: list[str] | None = None,
) -> None:
    if exclude_attrs is None:
        exclude_attrs = []

    attributes = meta["attribute_label"].dropna().unique()
    attributes = [attr for attr in attributes if attr not in exclude_attrs]

    n_attrs = len(attributes)
    ncols = 2
    nrows = int(np.ceil(n_attrs / ncols))

    fig, axes = new_figure(
        nrows=nrows,
        ncols=ncols,
        width="double",
        height_in=7,
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )

    axes = axes.ravel()

    for ax, attr in zip(axes, attributes):
        d = meta[meta["attribute_label"] == attr].sort_values("window_idx")
        d = d.merge(
            window_lookup[["window_idx", "wn_mid_date"]],
            on="window_idx",
            how="left",
        )

        x = _date_values(d["wn_mid_date"])

        ax.fill_between(
            x,
            d["pooled_ci_low"].to_numpy(),
            d["pooled_ci_high"].to_numpy(),
            alpha=0.25,
            color="#4C72B0",
            label="95% CI for pooled mean",
        )
        ax.plot(
            x,
            d["pooled_mean"].to_numpy(),
            color="#4C72B0",
            linewidth=2,
            label="Random-effects pooled mean",
        )

        ax.fill_between(
            x,
            d["q25"].to_numpy(),
            d["q75"].to_numpy(),
            alpha=0.20,
            color="#DD8452",
            label="Between-lineage IQR",
        )
        ax.plot(
            x,
            d["median"].to_numpy(),
            color="#DD8452",
            linewidth=1.5,
            linestyle="--",
            label="Median across lineages",
        )

        ax.axhline(0, color="black", linestyle="--", linewidth=1)
        ax.set_title(str(attr))
        ax.set_xlabel("")
        ax.set_ylabel("")
        date_axis(ax)

    for ax in axes[len(attributes) :]:
        ax.remove()

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.075))

    fig.supylabel("Pooled Window Assortativity")
    add_panel_labels(axes)
    styled_save_figure(fig, paths, FIGURE_NAME)


def build(paths: Paths) -> pd.DataFrame:
    window_meta, window_lookup = compatibility_window_pooled_meta(paths)
    summary = pooled_window_attribute_summary(window_meta)
    write_table(
        window_meta,
        "compatibility_window_pooled_meta",
        table_dir=paths.table_dir,
    )
    write_table(
        summary,
        "compatibility_window_pooled_summary",
        table_dir=paths.table_dir,
    )
    plot_pooled_window_meta(
        paths,
        window_meta,
        window_lookup,
        exclude_attrs=["Age band"],
    )
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
