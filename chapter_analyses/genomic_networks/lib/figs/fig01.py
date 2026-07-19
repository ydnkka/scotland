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
    policy_columns = ordered_policy_values(work["policy_era"], column="policy_era")
    category_order = ordered_sequence_categories(
        work, attribute, sort_by_total=sort_by_total
    )
    pivot = work.pivot_table(
        index="category",
        columns="policy_era",
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
    wrap_width: int = 15,
    x_max: float | None = None,
) -> None:
    y_positions = np.arange(len(shares.index))
    left = np.zeros(len(shares.index), dtype=float)
    for era in shares.columns:
        values = shares[era].to_numpy(dtype=float)
        ax.barh(
            y_positions,
            values,
            left=left,
            height=0.72,
            color=POLICY_COLORS.get(str(era), "#999999"),
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
    if panel:
        panel_label(ax, panel)


def build(paths: Paths) -> None:
    sequence = read_table(paths, "sequence_composition_by_policy")
    sex = sequence_policy_share_table(sequence, "sex")
    age_group = sequence_policy_share_table(sequence, "age_group")
    simd = sequence_policy_share_table(sequence, "simd_quintile")
    simd_map1 = {1: "Q1", 2: "Q2", 3: "Q3", 4: "Q4", 5: "Q5"}
    simd_map2 = {"1.0": "Q1", "2.0": "Q2", "3.0": "Q3", "4.0": "Q4", "5.0": "Q5"}
    try:
        simd.index = simd.index.astype(int).map(simd_map1)
    except ValueError:
        simd.index = simd.index.astype(str).map(simd_map2)
    urban_rural = sequence_policy_share_table(sequence, "urban_rural")
    health_board = sequence_policy_share_table(
        sequence, "health_board", sort_by_total=True
    )
    health_split = math.ceil(len(health_board) / 2)
    health_left = health_board.iloc[:health_split]
    health_right = health_board.iloc[health_split:]
    health_x_max = float(health_board.sum(axis=1).max()) * 1.12

    fig, grid = styled_new_figure(
        nrows=3,
        ncols=2,
        width="double", 
        height_in=8, 
        constrained_layout=True,
        gridspec_kw={"height_ratios": [0.2, 0.4, 0.5]},
    )
    axes = {
        "sex": grid[0, 0],
        "simd": grid[0, 1],
        "age_group": grid[1, 0],
        "urban_rural": grid[1, 1],
        "health_left": grid[2, 0],
        "health_right": grid[2, 1],
    }

    plot_policy_stacked_bars(axes["sex"], sex, title="Sex", panel="A", wrap_width=14)
    plot_policy_stacked_bars(
        axes["simd"],
        simd,
        title="SIMD quintile",
        panel="B",
    )
    plot_policy_stacked_bars(
        axes["age_group"],
        age_group,
        title="Age group",
        panel="C",
    )
    plot_policy_stacked_bars(
        axes["urban_rural"],
        urban_rural,
        title="Urban/rural class",
        panel="D",
    )
    plot_policy_stacked_bars(
        axes["health_left"],
        health_left,
        title="Health board",
        panel="E",
        x_max=health_x_max,
    )
    plot_policy_stacked_bars(
        axes["health_right"],
        health_right,
        title="Health board (continued)",
        x_max=health_x_max,
    )

    policy_columns = ordered_policy_values(
        sequence.loc[sequence["n_sequences"].fillna(0).gt(0), "policy_era"],
        column="policy_era",
    )
    legend_ncol = 3
    legend_columns = [
        era
        for column in range(legend_ncol)
        for era in policy_columns[column::legend_ncol]
    ]
    handles = [
        Patch(
            facecolor=POLICY_COLORS.get(str(era), "#999999"),
            edgecolor=POLICY_COLORS.get(str(era), "#999999"),
            label=POLICY_LABELS.get(str(era), str(era).upper().replace("_", " ")),
        )
        for era in legend_columns
    ]
    fig.supxlabel("Proportion of sequences by epidemic era", y=0.075)
    fig.legend(
        handles=handles,
        loc="outside lower center",
        bbox_to_anchor=(0.5, -0.06),
        ncol=legend_ncol,
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
        pad=0.055,
    )
    stringency_colorbar.set_label("Mean restriction stringency")
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
