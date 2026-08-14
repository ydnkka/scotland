"""Build descriptive SSE sequence-composition differences."""

from __future__ import annotations

import argparse
import math
import re
import textwrap
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Patch

from ..model.prep import COMPOSITION_SPECS
from ..sse.detection import load_sequence_data
from ..sse.io import HIGH_PRIORITY_CANDIDATE_TIERS, load_sse_outputs
from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
    styled_save_figure,
)

FIGURE_NAME = "fig_ch5_sse_composition_descriptive"
CANDIDATE_HIGHER_COLOR = "#2F6690"
BACKGROUND_HIGHER_COLOR = "#B75D69"
ZERO_COLOR = "#555555"

CATEGORY_ORDER = {
    "sex": ["Female", "Male"],
    "age_group": ["00-04", "05-14", "15-24", "25-64", "65-74", "75+"],
    "dz_simd_quintile": ["Q1", "Q2", "Q3", "Q4", "Q5"],
    "dz_urban_rural_class": [
        "Large Urban Areas",
        "Other Urban Areas",
        "Accessible Small Towns",
        "Remote Small Towns",
        "Accessible Rural",
        "Remote Rural",
    ],
}


def _composition_column(name: str) -> str:
    matches = [
        str(spec["column"])
        for spec in COMPOSITION_SPECS
        if spec["name"] == name
    ]
    if len(matches) != 1:
        raise KeyError(f"Expected one composition spec named {name!r}")
    return matches[0]


SEX_COLUMN = _composition_column("sex")
AGE_COLUMN = _composition_column("age_group")
SIMD_COLUMN = _composition_column("simd_quintile")
URBAN_RURAL_COLUMN = _composition_column("urban_rural_class")
HEALTH_BOARD_COLUMN = _composition_column("health_board")


def _category_text(value: Any) -> str:
    if pd.isna(value):
        return "Missing"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    text = str(value)
    if re.fullmatch(r"\d+\.0", text):
        return text[:-2]
    return text


def _display_category(column: str, value: Any) -> str:
    text = _category_text(value)
    if column == SIMD_COLUMN and text in {"1", "2", "3", "4", "5"}:
        return f"Q{text}"
    return text


def _wrap_labels(labels: pd.Index | list[object], width: int) -> list[str]:
    return [textwrap.fill(str(label), width=width) for label in labels]


def _ordered_categories(
    data: pd.DataFrame,
    column: str,
    *,
    sort_by_total: bool = False,
) -> list[str]:
    categories = data[column].map(lambda value: _display_category(column, value))
    totals = categories.value_counts(sort=False)
    if totals.empty:
        return []
    if sort_by_total:
        return totals.sort_values(ascending=False).index.astype(str).tolist()
    if column in CATEGORY_ORDER:
        ordered = [value for value in CATEGORY_ORDER[column] if value in totals.index]
        extras = sorted(
            [str(value) for value in totals.index if value not in set(ordered)]
        )
        return [*ordered, *extras]

    def key(value: object) -> tuple[int, object]:
        text = str(value)
        match = re.search(r"\d+", text)
        if match:
            return (0, int(match.group()))
        return (1, text)

    return sorted(totals.index.astype(str).tolist(), key=key)


def _eligible_sequence_data(paths: Paths) -> tuple[pd.DataFrame, dict[str, int]]:
    sse_outputs = load_sse_outputs(paths.table_dir)
    cluster = sse_outputs.cluster_table.copy()
    required_cluster = {"cluster_id", "cluster_size", "candidate_tier"}
    missing_cluster = required_cluster.difference(cluster.columns)
    if missing_cluster:
        raise KeyError(f"cluster_table is missing columns: {sorted(missing_cluster)}")

    cluster["candidate"] = cluster["candidate_tier"].isin(HIGH_PRIORITY_CANDIDATE_TIERS)
    candidate_sizes = cluster.loc[cluster["candidate"], "cluster_size"].dropna()
    if candidate_sizes.empty:
        raise ValueError("No high-priority candidate nodes were found.")
    min_candidate_size = int(candidate_sizes.min())

    sequence = load_sequence_data()
    composition_columns = [str(spec["column"]) for spec in COMPOSITION_SPECS]
    required_sequence = {"cluster_id", *composition_columns}
    missing_sequence = required_sequence.difference(sequence.columns)
    if missing_sequence:
        raise KeyError(f"sequence data is missing columns: {sorted(missing_sequence)}")

    eligible_nodes = cluster.loc[
        cluster["cluster_size"].ge(min_candidate_size),
        ["cluster_id", "candidate"],
    ]
    data = sequence.merge(
        eligible_nodes,
        on="cluster_id",
        how="inner",
        validate="many_to_one",
    )
    if data.empty:
        raise ValueError("No eligible sequence rows found for composition plotting.")
    if data["candidate"].nunique() != 2:
        raise ValueError(
            "Eligible sequence rows must include candidates and background."
        )

    counts = data["candidate"].value_counts().to_dict()
    summary = {
        "candidate": int(counts.get(True, 0)),
        "background": int(counts.get(False, 0)),
        "min_candidate_size": min_candidate_size,
    }
    return data, summary


def composition_difference_table(
    data: pd.DataFrame,
    column: str,
    *,
    sort_by_total: bool = False,
) -> pd.DataFrame:
    """Return candidate-minus-background composition differences.

    Differences are measured in percentage points.
    """
    work = data[[column, "candidate"]].copy()
    work["category"] = work[column].map(lambda value: _display_category(column, value))
    counts = pd.crosstab(work["category"], work["candidate"])
    counts = counts.reindex(columns=[False, True], fill_value=0)
    totals = counts.sum(axis=0).replace(0, np.nan)
    shares = counts.divide(totals, axis=1).fillna(0)

    out = pd.DataFrame(
        {
            "background_n": counts[False].astype(int),
            "candidate_n": counts[True].astype(int),
            "background_share": shares[False],
            "candidate_share": shares[True],
        }
    )
    out["difference_pp"] = (
        out["candidate_share"] - out["background_share"]
    ) * 100
    ordered = _ordered_categories(data, column, sort_by_total=sort_by_total)
    return out.reindex(ordered).dropna(how="all")


def _symmetric_x_limit(tables: list[pd.DataFrame]) -> float:
    absolute_differences = [
        float(table["difference_pp"].abs().max(skipna=True))
        for table in tables
        if not table.empty
    ]
    if not absolute_differences:
        return 1.0
    max_abs = max(absolute_differences)
    if not np.isfinite(max_abs) or max_abs == 0:
        return 1.0
    padded = max_abs * 1.15
    step = 5 if padded > 5 else 1
    return float(math.ceil(padded / step) * step)


def plot_difference_bars(
    ax: Axes,
    differences: pd.DataFrame,
    *,
    title: str,
    panel: str | None = None,
    wrap_width: int = 15,
    x_limit: float | None = None,
) -> None:
    y_positions = np.arange(len(differences.index))
    values = differences["difference_pp"].to_numpy(dtype=float)
    colors = np.where(values >= 0, CANDIDATE_HIGHER_COLOR, BACKGROUND_HIGHER_COLOR)
    ax.barh(
        y_positions,
        values,
        height=0.72,
        color=colors,
        edgecolor="white",
        linewidth=0.25,
    )
    ax.axvline(0, color=ZERO_COLOR, linewidth=0.8, zorder=1)
    ax.set_title(title)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(_wrap_labels(differences.index, wrap_width))
    ax.invert_yaxis()
    ax.grid(axis="x", color="#d9d9d9", lw=0.5, alpha=0.8)
    ax.set_axisbelow(True)
    if x_limit is not None:
        ax.set_xlim(-x_limit, x_limit)
    if panel:
        add_panel_labels(ax, label=panel)


def build(paths: Paths) -> dict[str, object]:
    """Create the six-panel candidate-background sequence composition figure."""
    data, sample_sizes = _eligible_sequence_data(paths)

    sex = composition_difference_table(data, SEX_COLUMN)
    age_group = composition_difference_table(data, AGE_COLUMN)
    simd = composition_difference_table(data, SIMD_COLUMN)
    urban_rural = composition_difference_table(data, URBAN_RURAL_COLUMN)
    health_board = composition_difference_table(
        data,
        HEALTH_BOARD_COLUMN,
        sort_by_total=True,
    )
    health_split = math.ceil(len(health_board) / 2)
    health_left = health_board.iloc[:health_split]
    health_right = health_board.iloc[health_split:]
    x_limit = _symmetric_x_limit(
        [sex, age_group, simd, urban_rural, health_left, health_right]
    )

    fig, grid = new_figure(
        nrows=3,
        ncols=2,
        width="double",
        height_in=8,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [0.2, 0.4, 0.5]},
        sharex=True,
    )
    axes = {
        "sex": grid[0, 0],
        "simd": grid[0, 1],
        "age_group": grid[1, 0],
        "urban_rural": grid[1, 1],
        "health_left": grid[2, 0],
        "health_right": grid[2, 1],
    }

    plot_difference_bars(
        axes["sex"],
        sex,
        title="Sex",
        panel="A",
        wrap_width=14,
        x_limit=x_limit,
    )
    plot_difference_bars(
        axes["simd"],
        simd,
        title="SIMD quintile",
        panel="B",
        x_limit=x_limit,
    )
    plot_difference_bars(
        axes["age_group"],
        age_group,
        title="Age group",
        panel="C",
        x_limit=x_limit,
    )
    plot_difference_bars(
        axes["urban_rural"],
        urban_rural,
        title="Urban/rural class",
        panel="D",
        x_limit=x_limit,
    )
    plot_difference_bars(
        axes["health_left"],
        health_left,
        title="Health board",
        panel="E",
        x_limit=x_limit,
    )
    plot_difference_bars(
        axes["health_right"],
        health_right,
        title="Health board (continued)",
        panel="E",
        x_limit=x_limit,
    )

    xlabel = "Percentage points difference"
    axes["health_left"].set_xlabel(xlabel)
    axes["health_right"].set_xlabel(xlabel)
    fig.legend(
        handles=[
            Patch(
                facecolor=CANDIDATE_HIGHER_COLOR,
                edgecolor=CANDIDATE_HIGHER_COLOR,
                label="Candidate higher",
            ),
            Patch(
                facecolor=BACKGROUND_HIGHER_COLOR,
                edgecolor=BACKGROUND_HIGHER_COLOR,
                label="Background higher",
            ),
        ],
        loc="outside upper center",
        ncol=2,
        columnspacing=1.5,
        handlelength=1.5,
        frameon=False,
    )

    plot_data = pd.concat(
        {
            "sex": sex,
            "simd": simd,
            "age_group": age_group,
            "urban_rural": urban_rural,
            "health_board": health_board,
        },
        names=["panel", "category"],
    ).reset_index()
    outputs = styled_save_figure(fig, paths, FIGURE_NAME)
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": plot_data,
        "sample_sizes": sample_sizes,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--n-bootstrap", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        "--random-state",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    print(f"Wrote {FIGURE_NAME} to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
