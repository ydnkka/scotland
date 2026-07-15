"""Build Chapter 4 Figure 1: sequence composition by policy period."""

from __future__ import annotations

from pathlib import Path
import argparse
import math
import re
import sys
import textwrap


from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
from matplotlib.axes import Axes
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    POLICY_COLORS,
    POLICY_LABELS,
    POLICY_STRINGENCY_CMAP,
    POLICY_STRINGENCY_NORM,
    Paths,
    add_common_args,
    ordered_policy_values,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)


FIGURE_NAME = "fig_ch4_sequence_composition_by_policy"


def ordered_sequence_categories(
    table: pd.DataFrame,
    attribute: str,
    *,
    sort_by_total: bool = False,
) -> list[str]:
    work = table.loc[table["attribute"].eq(attribute)].copy()
    if work.empty:
        return []
    totals = work.groupby("category", observed=False)["n_sequences"].sum()
    if sort_by_total:
        return totals.sort_values(ascending=False).index.astype(str).tolist()
    if attribute == "sex":
        preferred = ["Female", "Male"]
        return [value for value in preferred if value in totals.index]
    if attribute == "urban_rural":
        preferred = [
            "Large Urban Areas",
            "Other Urban Areas",
            "Accessible Small Towns",
            "Remote Small Towns",
            "Accessible Rural",
            "Remote Rural",
        ]
        return [value for value in preferred if value in totals.index]

    def key(value: object) -> tuple[int, object]:
        text = str(value)
        match = re.search(r"\d+", text)
        if match:
            return (0, int(match.group()))
        return (1, text)

    return sorted(totals.index.astype(str).tolist(), key=key)


def sequence_policy_share_table(
    table: pd.DataFrame,
    attribute: str,
    *,
    sort_by_total: bool = False,
) -> pd.DataFrame:
    work = table.loc[table["attribute"].eq(attribute)].copy()
    policy_columns = ordered_policy_values(work["policy_period"])
    category_order = ordered_sequence_categories(
        work, attribute, sort_by_total=sort_by_total
    )
    pivot = work.pivot_table(
        index="category",
        columns="policy_period",
        values="n_sequences",
        aggfunc="sum",
        fill_value=0,
        observed=False,
    )
    pivot = pivot.reindex(index=category_order, columns=policy_columns, fill_value=0)
    total = pivot.to_numpy().sum()
    if total <= 0:
        return pivot.astype(float)
    return pivot / total


def wrap_labels(labels: pd.Index | list[object], width: int) -> list[str]:
    return [textwrap.fill(str(label), width=width) for label in labels]


def plot_policy_stacked_bars(
    ax: Axes,
    shares: pd.DataFrame,
    *,
    title: str,
    panel: str | None = None,
    wrap_width: int = 18,
    x_max: float | None = None,
) -> None:
    y_positions = np.arange(len(shares.index))
    left = np.zeros(len(shares.index), dtype=float)
    for period in shares.columns:
        values = shares[period].to_numpy(dtype=float)
        ax.barh(
            y_positions,
            values,
            left=left,
            height=0.72,
            color=POLICY_COLORS.get(str(period), "#999999"),
            edgecolor="white",
            linewidth=0.25,
        )
        left += values

    ax.set_title(title)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(wrap_labels(shares.index, wrap_width))
    ax.invert_yaxis()
    ax.xaxis.set_major_formatter(PercentFormatter(1.0))
    ax.grid(axis="x", color="#d9d9d9", lw=0.5, alpha=0.8)
    ax.set_axisbelow(True)
    max_value = float(left.max()) if len(left) else 0.0
    ax.set_xlim(0, x_max if x_max is not None else max_value * 1.12)
    ax.set_xlabel("Share of sequenced records")
    if panel:
        panel_label(ax, panel)


def build(paths: Paths) -> None:
    sequence = read_table(paths, "sequence_composition_by_policy")
    sex = sequence_policy_share_table(sequence, "sex")
    age_group = sequence_policy_share_table(sequence, "age_group")
    simd = sequence_policy_share_table(sequence, "simd_quintile")
    urban_rural = sequence_policy_share_table(sequence, "urban_rural")
    health_board = sequence_policy_share_table(
        sequence, "health_board", sort_by_total=True
    )
    health_split = math.ceil(len(health_board) / 2)
    health_left = health_board.iloc[:health_split]
    health_right = health_board.iloc[health_split:]
    health_x_max = float(health_board.sum(axis=1).max()) * 1.12

    fig, placeholder_ax = styled_new_figure(
        width="double", height_in=8.6, constrained_layout=True
    )
    placeholder_ax.remove()
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[0.85, 1.15, 1.75],
    )
    axes = {
        "sex": fig.add_subplot(grid[0, 0]),
        "age_group": fig.add_subplot(grid[0, 1]),
        "simd": fig.add_subplot(grid[1, 0]),
        "urban_rural": fig.add_subplot(grid[1, 1]),
        "health_left": fig.add_subplot(grid[2, 0]),
        "health_right": fig.add_subplot(grid[2, 1]),
    }

    plot_policy_stacked_bars(
        axes["sex"], sex, title="Sex", panel="A", wrap_width=14
    )
    plot_policy_stacked_bars(
        axes["age_group"],
        age_group,
        title="Age group",
        panel="B",
        wrap_width=14,
    )
    plot_policy_stacked_bars(
        axes["simd"],
        simd,
        title="SIMD quintile",
        panel="C",
        wrap_width=14,
    )
    plot_policy_stacked_bars(
        axes["urban_rural"],
        urban_rural,
        title="Urban/rural class",
        panel="D",
        wrap_width=18,
    )
    plot_policy_stacked_bars(
        axes["health_left"],
        health_left,
        title="Health board",
        panel="E",
        wrap_width=21,
        x_max=health_x_max,
    )
    plot_policy_stacked_bars(
        axes["health_right"],
        health_right,
        title="Health board (continued)",
        wrap_width=21,
        x_max=health_x_max,
    )

    policy_columns = ordered_policy_values(
        sequence.loc[sequence["n_sequences"].fillna(0).gt(0), "policy_period"]
    )
    legend_ncol = 4
    legend_columns = [
        period
        for column in range(legend_ncol)
        for period in policy_columns[column::legend_ncol]
    ]
    handles = [
        Patch(
            facecolor=POLICY_COLORS.get(str(period), "#999999"),
            edgecolor="none",
            label=f"{POLICY_LABELS.get(str(period), str(period))} ({period})",
        )
        for period in legend_columns
    ]
    fig.legend(
        handles=handles,
        loc="outside lower center",
        ncol=legend_ncol,
        title="Policy period",
        columnspacing=1.2,
        handlelength=1.5,
        labelspacing=0.35,
        frameon=False,
    )
    stringency_mappable = ScalarMappable(
        norm=POLICY_STRINGENCY_NORM,
        cmap=POLICY_STRINGENCY_CMAP,
    )
    stringency_colorbar = fig.colorbar(
        stringency_mappable,
        ax=list(axes.values()),
        orientation="horizontal",
        location="bottom",
        ticks=[0, 25, 50, 75, 100],
        shrink=0.42,
        aspect=45,
        pad=0.035,
    )
    stringency_colorbar.set_label("Mean policy stringency index")
    styled_save_figure(fig, paths, FIGURE_NAME, tight=False)


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
