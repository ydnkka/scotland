"""Build candidate-background entropy-tertile profile differences."""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

from ..sse.io import HIGH_PRIORITY_CANDIDATE_TIERS, write_table
from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
    read_table,
    styled_save_figure,
)

FILE_NAME = "entropy_tertile_profiles"

TERTILE_ORDER = (
    "more_homogeneous",
    "as_expected",
    "more_mixed",
)
TERTILE_LABELS = {
    "more_homogeneous": "Lower tertile",
    "as_expected": "Middle tertile",
    "more_mixed": "Upper tertile",
}
TERTILE_ABBREVIATIONS = {
    "more_homogeneous": "L",
    "as_expected": "M",
    "more_mixed": "U",
}
TERTILE_COLORS = {
    "more_homogeneous": "#1B9E77",
    "as_expected": "#BDBDBD",
    "more_mixed": "#D95F02",
}


@dataclass(frozen=True)
class MixingFeature:
    prefix: str
    label: str


MIXING_FEATURES = (
    MixingFeature("sex", "Sex"),
    MixingFeature("age", "Age group"),
    MixingFeature("simd", "SIMD quintile"),
    MixingFeature("urban_rural", "Urban/rural class"),
    MixingFeature("health_board", "Health board"),
    MixingFeature("local_authority", "Local authority"),
)
ENTROPY_KINDS = (
    ("obs", "Observed"),
    ("z", "Null-standardised"),
)


def _bar_label(feature: MixingFeature, kind_label: str) -> str:
    return f"{feature.label} ({kind_label})"


def _eligible_nodes(nodes: pd.DataFrame) -> pd.DataFrame:
    if "candidate_tier" not in nodes.columns:
        raise KeyError("cluster_table is missing 'candidate_tier'")
    if "cluster_size" not in nodes.columns:
        raise KeyError("cluster_table is missing 'cluster_size'")

    eligible = nodes.copy()
    eligible["candidate"] = eligible["candidate_tier"].isin(
        HIGH_PRIORITY_CANDIDATE_TIERS
    )
    candidate_sizes = eligible.loc[eligible["candidate"], "cluster_size"].dropna()
    if candidate_sizes.empty:
        raise ValueError("No detector candidates found for entropy-tertile profiles.")
    min_candidate_size = int(candidate_sizes.min())
    eligible = eligible.loc[eligible["cluster_size"].ge(min_candidate_size)].copy()
    if eligible["candidate"].nunique() != 2:
        raise ValueError(
            "Eligible cluster rows must include candidates and background."
        )
    eligible["sse_status"] = np.where(
        eligible["candidate"], "candidate", "background"
    )
    eligible["min_candidate_size"] = min_candidate_size

    sort_columns = [
        column
        for column in ("wn_mid_date", "candidate_tier", "cluster_size")
        if column in eligible.columns
    ]
    if sort_columns:
        eligible = eligible.sort_values(sort_columns)
    return eligible


def build_entropy_tertile_profiles(nodes: pd.DataFrame) -> pd.DataFrame:
    """Build candidate-background entropy-tertile profile differences."""
    eligible = _eligible_nodes(nodes)
    rows: list[dict[str, object]] = []
    missing_columns: list[str] = []
    columns_with_missing_values: list[str] = []

    for feature in MIXING_FEATURES:
        for kind, kind_label in ENTROPY_KINDS:
            column = f"{feature.prefix}_entropy_{kind}_tertile"
            if column not in eligible.columns:
                missing_columns.append(column)
                continue
            values = eligible[column].astype("string")
            if values.isna().any():
                columns_with_missing_values.append(column)
                continue

            work = pd.DataFrame(
                {
                    "candidate": eligible["candidate"].to_numpy(),
                    "tertile": values.astype(str).to_numpy(),
                }
            )
            counts = pd.crosstab(work["candidate"], work["tertile"]).reindex(
                index=[False, True],
                columns=TERTILE_ORDER,
                fill_value=0,
            )
            background_total = counts.loc[False].sum()
            candidate_total = counts.loc[True].sum()
            for tertile in TERTILE_ORDER:
                background_n = counts.loc[False, tertile]
                candidate_n = counts.loc[True, tertile]
                background_proportion = (
                    background_n / background_total if background_total else np.nan
                )
                candidate_proportion = (
                    candidate_n / candidate_total if candidate_total else np.nan
                )
                rows.append(
                    {
                        "bar_label": _bar_label(feature, kind_label),
                        "feature": feature.prefix,
                        "feature_label": feature.label,
                        "kind": kind,
                        "kind_label": kind_label,
                        "tertile": tertile,
                        "tertile_label": TERTILE_LABELS[str(tertile)],
                        "background_n": background_n,
                        "candidate_n": candidate_n,
                        "background_total": background_total,
                        "candidate_total": candidate_total,
                        "background_proportion": background_proportion,
                        "candidate_proportion": candidate_proportion,
                        "difference": candidate_proportion - background_proportion,
                        "difference_pp": (
                            candidate_proportion - background_proportion
                        )
                        * 100,
                        "eligible_total": background_total + candidate_total,
                        "min_candidate_size": int(
                            eligible["min_candidate_size"].iloc[0]
                        ),
                    }
                )

    if missing_columns:
        raise KeyError(
            "cluster_table is missing entropy tertile columns: "
            + ", ".join(missing_columns)
        )
    if columns_with_missing_values:
        raise ValueError(
            "Unexpected missing entropy tertile values in: "
            + ", ".join(columns_with_missing_values)
        )
    return pd.DataFrame(rows)


def _feature_order() -> list[str]:
    return [feature.label for feature in MIXING_FEATURES]


def _symmetric_x_limit(table: pd.DataFrame) -> float:
    max_abs = float(table["difference_pp"].abs().max(skipna=True))
    if not np.isfinite(max_abs) or max_abs == 0:
        return 1.0
    padded = max_abs * 1.15
    step = 5 if padded > 5 else 1
    return float(np.ceil(padded / step) * step)


def _legend_handles() -> list[Patch]:
    return [
        Patch(
            facecolor=TERTILE_COLORS[tertile],
            edgecolor=TERTILE_COLORS[tertile],
            label=f"{TERTILE_LABELS[tertile]}",
        )
        for tertile in TERTILE_ORDER
    ]


def _draw_kind_profile(
    ax: Axes,
    table: pd.DataFrame,
    *,
    kind: str,
    kind_label: str,
    x_limit: float,
    show_y_labels: bool,
) -> None:
    y_positions = np.arange(len(MIXING_FEATURES))
    offsets = {
        "more_homogeneous": -0.24,
        "as_expected": 0.0,
        "more_mixed": 0.24,
    }

    for y, feature in zip(y_positions, MIXING_FEATURES):
        data = table.loc[
            table["feature"].eq(feature.prefix) & table["kind"].eq(kind)
        ]
        for tertile in TERTILE_ORDER:
            row = data.loc[data["tertile"].eq(tertile)]
            if row.empty:
                continue
            difference = float(row["difference_pp"].iloc[0])
            ax.barh(
                y + offsets[tertile],
                difference,
                height=0.21,
                color=TERTILE_COLORS[tertile],
                edgecolor="white",
                linewidth=0.45,
            )

    ax.set_title(kind_label)
    ax.set_yticks(y_positions)
    if show_y_labels:
        ax.set_yticklabels(_feature_order())
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.set_xlim(-x_limit, x_limit)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.0f}"))
    ax.axvline(0, color="#555555", lw=0.8)
    ax.grid(axis="x", color="#d9d9d9", lw=0.5, alpha=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)


def draw_entropy_tertile_profiles(axes, table: pd.DataFrame) -> None:
    """Draw candidate-minus-background tertile differences in entropy columns."""
    axes = np.asarray(axes).ravel()
    if len(axes) != len(ENTROPY_KINDS):
        raise ValueError(f"Expected {len(ENTROPY_KINDS)} axes, got {len(axes)}")

    x_limit = _symmetric_x_limit(table)
    for index, (ax, (kind, kind_label)) in enumerate(zip(axes, ENTROPY_KINDS)):
        _draw_kind_profile(
            ax,
            table,
            kind=kind,
            kind_label=kind_label,
            x_limit=x_limit,
            show_y_labels=index == 0,
        )
        ax.set_xlabel("Percentage points\n(candidate minus background share)")

    axes[0].invert_yaxis()

    axes[0].figure.legend(
        handles=_legend_handles(),
        loc="lower center",
        bbox_to_anchor=(0.55, 0.05),
        ncol=len(TERTILE_ORDER),
        columnspacing=1.35,
        handlelength=1.5,
        frameon=False,
    )


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_entropy_tertile_profiles(nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")

    fig, axes = new_figure(
        ncols=2,
        width="double",
        height_in=4.8,
        sharey=True,
    )
    draw_entropy_tertile_profiles(axes, table)
    fig.subplots_adjust(
        left=0.20,
        right=0.98,
        top=0.86,
        bottom=0.24,
        wspace=0.08,
    )
    add_panel_labels(axes)
    outputs = styled_save_figure(fig, paths, f"fig_{FILE_NAME}")
    return {"figure": fig, "outputs": outputs, "plot_data": table}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    paths = paths_from_args(parser.parse_args())
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
