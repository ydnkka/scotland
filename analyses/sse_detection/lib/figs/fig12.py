"""Build cluster-level composition profile figures."""

from __future__ import annotations

import argparse
import math
import textwrap
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, PercentFormatter

from ..sse.io import write_table
from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
    read_table,
    styled_save_figure,
)

FILE_NAME = "cluster_composition_profiles"
STATUS_ORDER = ("background", "candidate")
STATUS_LABELS = {
    "background": "Background",
    "candidate": "Candidate",
}
CANDIDATE_HIGHER_COLOR = "#2F6690"
BACKGROUND_HIGHER_COLOR = "#B75D69"
STATUS_COLORS = {
    "background": BACKGROUND_HIGHER_COLOR,
    "candidate": CANDIDATE_HIGHER_COLOR,
}
ZERO_COLOR = "#555555"


@dataclass(frozen=True)
class CompositionFigureSpec:
    """Configuration for one cluster-composition figure."""

    name: str
    variable: str
    label: str
    level_order: tuple[str, ...] | None = None
    level_labels: dict[str, str] = field(default_factory=dict)
    sort_by_total: bool = False

    @property
    def table_name(self) -> str:
        return f"cluster_composition_{self.variable}"

    @property
    def figure_name(self) -> str:
        return f"fig_cluster_composition_{self.name}"

    @property
    def summary_table_name(self) -> str:
        return f"tab_cluster_composition_{self.name}"


COMPOSITION_FIGURES = (
    CompositionFigureSpec(
        name="sex",
        variable="sex",
        label="Sex",
        level_order=("Female", "Male"),
    ),
    CompositionFigureSpec(
        name="age_group",
        variable="age_group",
        label="Age group",
        level_order=("00-04", "05-14", "15-24", "25-64", "65-74", "75+"),
    ),
    CompositionFigureSpec(
        name="simd_quintile",
        variable="dz_simd_quintile",
        label="SIMD quintile",
        level_order=("1", "2", "3", "4", "5"),
        level_labels={str(i): f"Q{i}" for i in range(1, 6)},
    ),
    CompositionFigureSpec(
        name="urban_rural_class",
        variable="dz_urban_rural_class",
        label="Urban/rural class",
        level_order=(
            "Large Urban Areas",
            "Other Urban Areas",
            "Accessible Small Towns",
            "Remote Small Towns",
            "Accessible Rural",
            "Remote Rural",
        ),
    ),
    CompositionFigureSpec(
        name="health_board",
        variable="dz_health_board",
        label="Health board",
        sort_by_total=True,
    ),
)


def _normalise_level(value: object) -> str:
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _display_level(spec: CompositionFigureSpec, value: object) -> str:
    text = _normalise_level(value)
    return spec.level_labels.get(text, text)


def _level_columns(table: pd.DataFrame) -> list[str]:
    metadata = {"cluster_id", "sse_status", "burst_score", "burden_score"}
    return [column for column in table.columns if column not in metadata]


def _ordered_levels(table: pd.DataFrame, spec: CompositionFigureSpec) -> list[str]:
    levels = [_normalise_level(column) for column in _level_columns(table)]
    if not levels:
        raise ValueError(f"{spec.table_name} has no composition-level columns")

    if spec.level_order is not None:
        ordered = [level for level in spec.level_order if level in levels]
        extras = sorted(level for level in levels if level not in set(ordered))
        return [*ordered, *extras]

    if spec.sort_by_total:
        values = table[_level_columns(table)].apply(pd.to_numeric, errors="coerce")
        level_totals = values.mean(axis=0, skipna=True)
        order = level_totals.sort_values(ascending=False).index
        return [_normalise_level(column) for column in order]

    return levels


def build_composition_long(
    table: pd.DataFrame,
    spec: CompositionFigureSpec,
) -> pd.DataFrame:
    """Return one row per cluster, level, and detector status."""
    required = {"cluster_id", "sse_status"}
    missing = sorted(required.difference(table.columns))
    if missing:
        raise KeyError(f"{spec.table_name} is missing columns: {missing}")

    levels = _ordered_levels(table, spec)
    missing_levels = sorted(set(levels).difference(table.columns))
    if missing_levels:
        raise KeyError(f"{spec.table_name} is missing level columns: {missing_levels}")

    work = table[["cluster_id", "sse_status", *levels]].copy()
    work["sse_status"] = work["sse_status"].astype("string").str.lower()
    unexpected_statuses = sorted(set(work["sse_status"].dropna()) - set(STATUS_ORDER))
    if unexpected_statuses:
        raise ValueError(f"Unexpected sse_status values: {unexpected_statuses}")

    for column in levels:
        work[column] = pd.to_numeric(work[column], errors="coerce")

    long = work.melt(
        id_vars=["cluster_id", "sse_status"],
        value_vars=levels,
        var_name="level",
        value_name="proportion",
    )
    long["level"] = pd.Categorical(long["level"], categories=levels, ordered=True)
    long["level_label"] = long["level"].map(lambda value: _display_level(spec, value))
    long["sse_status"] = pd.Categorical(
        long["sse_status"],
        categories=list(STATUS_ORDER),
        ordered=True,
    )
    return long.sort_values(["level", "sse_status", "cluster_id"])


def summarise_composition(long: pd.DataFrame) -> pd.DataFrame:
    """Summarise cluster-level composition by detector status and level."""
    summary = (
        long.groupby(["sse_status", "level", "level_label"], observed=True)[
            "proportion"
        ]
        .agg(
            n_clusters="count",
            mean_proportion="mean",
            median_proportion="median",
            q25=lambda values: values.quantile(0.25),
            q75=lambda values: values.quantile(0.75),
        )
        .reset_index()
    )

    means = summary.pivot(
        index=["level", "level_label"],
        columns="sse_status",
        values="mean_proportion",
    )
    difference = (
        means.get("candidate", pd.Series(dtype=float))
        - means.get("background", pd.Series(dtype=float))
    ) * 100
    difference = difference.rename("mean_difference_pp").reset_index()
    return summary.merge(difference, on=["level", "level_label"], how="left")


def _status_handles() -> list[Patch]:
    return [
        Patch(
            facecolor=STATUS_COLORS[status],
            edgecolor=STATUS_COLORS[status],
            label=STATUS_LABELS[status],
        )
        for status in STATUS_ORDER
    ]


def _wrap_labels(labels: Sequence[object], width: int) -> list[str]:
    return [textwrap.fill(str(label), width=width) for label in labels]


def _draw_distribution_boxplots(
    ax: Axes,
    long: pd.DataFrame,
    spec: CompositionFigureSpec,
) -> None:
    levels = list(long["level"].cat.categories)
    labels = [_display_level(spec, level) for level in levels]
    positions = np.arange(len(levels))
    offsets = {"background": -0.17, "candidate": 0.17}

    for level_index, level in enumerate(levels):
        level_data = long.loc[long["level"].eq(level)]
        for status in STATUS_ORDER:
            values = level_data.loc[
                level_data["sse_status"].eq(status),
                "proportion",
            ].dropna()
            if values.empty:
                continue
            plot = ax.boxplot(
                [values.to_numpy(dtype=float)],
                vert=False,
                positions=[level_index + offsets[status]],
                widths=0.30,
                whis=(5, 95),
                patch_artist=True,
                showfliers=False,
                manage_ticks=False,
            )
            for box in plot["boxes"]:
                box.set(
                    facecolor=STATUS_COLORS[status], edgecolor="#333333", alpha=0.78
                )
            for median in plot["medians"]:
                median.set(color="#222222", linewidth=1.0)
            for item in (*plot["whiskers"], *plot["caps"]):
                item.set(color="#555555", linewidth=0.8)

    ax.set_title("Cluster-level distributions")
    ax.set_xlabel("Within-cluster proportion")
    ax.set_yticks(positions)
    ax.set_yticklabels(_wrap_labels(labels, width=18))
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.02)
    ax.xaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.grid(axis="x", color="#d9d9d9", lw=0.5, alpha=0.8)
    ax.set_axisbelow(True)
    ax.legend(handles=_status_handles(), loc="lower right")


def _symmetric_x_limit(values: pd.Series) -> float:
    max_abs = float(values.abs().max(skipna=True))
    if not np.isfinite(max_abs) or max_abs == 0:
        return 1.0
    padded = max_abs * 1.15
    step = 5 if padded > 5 else 1
    return float(math.ceil(padded / step) * step)


def _draw_mean_difference_bars(
    ax: Axes,
    summary: pd.DataFrame,
    spec: CompositionFigureSpec,
) -> None:
    levels = list(summary["level"].cat.categories)
    level_rows = (
        summary[["level", "level_label", "mean_difference_pp"]]
        .drop_duplicates("level")
        .set_index("level")
        .reindex(levels)
    )
    values = level_rows["mean_difference_pp"].astype(float)
    colors = np.where(values >= 0, CANDIDATE_HIGHER_COLOR, BACKGROUND_HIGHER_COLOR)
    positions = np.arange(len(levels))

    ax.barh(
        positions,
        values.to_numpy(),
        height=0.64,
        color=colors,
        edgecolor="white",
        linewidth=0.25,
    )
    ax.axvline(0, color=ZERO_COLOR, linewidth=0.8, zorder=1)
    ax.set_title("Mean difference")
    ax.set_xlabel("Percentage points\n(candidate minus background)")
    ax.set_yticks(positions)
    ax.tick_params(axis="y", labelleft=False, length=0)
    ax.invert_yaxis()
    ax.set_xlim(-_symmetric_x_limit(values), _symmetric_x_limit(values))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.0f}"))
    ax.grid(axis="x", color="#d9d9d9", lw=0.5, alpha=0.8)
    ax.set_axisbelow(True)


def draw_composition_profile(
    axes: Sequence[Axes],
    long: pd.DataFrame,
    summary: pd.DataFrame,
    spec: CompositionFigureSpec,
) -> None:
    """Draw one composition-profile figure."""
    if len(axes) != 2:
        raise ValueError(f"Expected two axes, got {len(axes)}")
    _draw_distribution_boxplots(axes[0], long, spec)
    _draw_mean_difference_bars(axes[1], summary, spec)


def build_one(paths: Paths, spec: CompositionFigureSpec) -> dict[str, object]:
    """Create one cluster-level composition figure."""
    table = read_table(paths, spec.table_name)
    long = build_composition_long(table, spec)
    summary = summarise_composition(long)
    write_table(summary, paths.result_table_dir, spec.summary_table_name)

    n_levels = int(long["level"].nunique())
    height_in = max(3.2, 1.25 + n_levels * 0.35)
    fig, axes = new_figure(
        nrows=1,
        ncols=2,
        width="double",
        height_in=height_in,
        gridspec_kw={"width_ratios": [1.55, 1.0]},
    )
    draw_composition_profile(axes, long, summary, spec)
    add_panel_labels(axes)
    fig.subplots_adjust(left=0.24, right=0.98, top=0.87, bottom=0.16, wspace=0.18)

    outputs = styled_save_figure(fig, paths, spec.figure_name)
    return {
        "figure": fig,
        "outputs": outputs,
        "plot_data": summary,
    }


def _selected_specs(
    variable_names: Sequence[str] | None,
) -> tuple[CompositionFigureSpec, ...]:
    if not variable_names:
        return COMPOSITION_FIGURES

    by_name: dict[str, CompositionFigureSpec] = {}
    for spec in COMPOSITION_FIGURES:
        by_name[spec.name] = spec
        by_name[spec.variable] = spec

    selected = []
    unknown = []
    for name in variable_names:
        spec = by_name.get(name)
        if spec is None:
            unknown.append(name)
            continue
        if spec not in selected:
            selected.append(spec)
    if unknown:
        available = ", ".join(spec.name for spec in COMPOSITION_FIGURES)
        raise KeyError(
            f"Unknown composition variable(s): {unknown}. Available: {available}"
        )
    return tuple(selected)


def build(
    paths: Paths,
    *,
    variables: Sequence[str] | None = None,
) -> dict[str, dict[str, object]]:
    """Create all configured cluster-level composition figures."""
    return {spec.name: build_one(paths, spec) for spec in _selected_specs(variables)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument(
        "--variables",
        nargs="+",
        help=(
            "Optional subset to build. Use short names such as age_group, "
            "simd_quintile, or source table variables such as dz_simd_quintile."
        ),
    )
    args = parser.parse_args()
    paths = paths_from_args(args)
    results = build(paths, variables=args.variables)
    for name, result in results.items():
        outputs = result["outputs"]
        print(f"Wrote {name}: {outputs['pdf']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
