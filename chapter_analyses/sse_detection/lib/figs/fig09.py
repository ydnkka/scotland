"""Build the entropy-tertile profiles among candidates figure."""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter

from ..sse.io import HIGH_PRIORITY_CANDIDATE_TIERS, write_table
from .common import (
    Paths,
    add_common_args,
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
    "more_homogeneous": "More homogeneous",
    "as_expected": "As expected",
    "more_mixed": "More mixed",
}
TERTILE_ABBREVIATIONS = {
    "more_homogeneous": "H",
    "as_expected": "E",
    "more_mixed": "M",
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
    ("z", "Null-adjusted"),
)


def _bar_label(feature: MixingFeature, kind_label: str) -> str:
    return f"{feature.label} ({kind_label})"


def _candidate_nodes(nodes: pd.DataFrame) -> pd.DataFrame:
    if "candidate_tier" not in nodes.columns:
        raise KeyError("cluster_table is missing 'candidate_tier'")
    candidates = nodes.loc[nodes["candidate_tier"].isin(HIGH_PRIORITY_CANDIDATE_TIERS)].copy()
    if candidates.empty:
        raise ValueError("No detector candidates found for entropy-tertile profiles.")
    sort_columns = [
        column
        for column in ("wn_mid_date", "candidate_tier", "cluster_size")
        if column in candidates.columns
    ]
    if sort_columns:
        candidates = candidates.sort_values(sort_columns)
    return candidates


def build_entropy_tertile_profiles(nodes: pd.DataFrame) -> pd.DataFrame:
    """Build profile proportions for observed and null-adjusted entropy tertiles."""
    candidates = _candidate_nodes(nodes)
    rows: list[dict[str, object]] = []
    missing_columns: list[str] = []
    columns_with_missing_values: list[str] = []

    for feature in MIXING_FEATURES:
        for kind, kind_label in ENTROPY_KINDS:
            column = f"{feature.prefix}_entropy_{kind}_tertile"
            if column not in candidates.columns:
                missing_columns.append(column)
                continue
            values = candidates[column].astype("string")
            if values.isna().any():
                columns_with_missing_values.append(column)
                values = values.dropna()
            counts = values.value_counts().reindex(TERTILE_ORDER, fill_value=0)
            total = int(counts.sum())
            for tertile, n in counts.items():
                rows.append(
                    {
                        "bar_label": _bar_label(feature, kind_label),
                        "feature": feature.prefix,
                        "feature_label": feature.label,
                        "kind": kind,
                        "kind_label": kind_label,
                        "tertile": tertile,
                        "tertile_label": TERTILE_LABELS[str(tertile)],
                        "candidate_n": int(n),
                        "candidate_total": total,
                        "proportion": n / total if total else np.nan,
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


def _bar_order() -> list[str]:
    return [
        _bar_label(feature, kind_label)
        for feature in MIXING_FEATURES
        for _, kind_label in ENTROPY_KINDS
    ]


def draw_entropy_tertile_profiles(ax: Axes, table: pd.DataFrame) -> None:
    y_positions = np.arange(len(_bar_order()))
    for y, bar_label in zip(y_positions, _bar_order()):
        data = table.loc[table["bar_label"].eq(bar_label)]
        left = 0.0
        for tertile in TERTILE_ORDER:
            row = data.loc[data["tertile"].eq(tertile)]
            if row.empty:
                continue
            width = float(row["proportion"].iloc[0])
            ax.barh(
                y,
                width,
                left=left,
                height=0.72,
                color=TERTILE_COLORS[tertile],
                edgecolor="white",
                linewidth=0.45,
            )
            if width >= 0.055:
                ax.text(
                    left + width / 2,
                    float(y),
                    f"{TERTILE_ABBREVIATIONS[tertile]}: {width:.0%}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    linespacing=0.85,
                )
            left += width

    ax.set_yticks(y_positions, _bar_order())
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_xlabel("Proportion of detector candidates")
    ax.invert_yaxis()
    ax.grid(axis="x", color="#d9d9d9", lw=0.5, alpha=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.legend(
        handles=[
            Patch(
                facecolor=TERTILE_COLORS[tertile],
                edgecolor=TERTILE_COLORS[tertile],
                label=TERTILE_LABELS[tertile],
            )
            for tertile in TERTILE_ORDER
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        ncol=len(TERTILE_ORDER),
        columnspacing=1.35,
        handlelength=1.5,
        frameon=False,
    )


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_entropy_tertile_profiles(nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")

    fig, ax = new_figure(
        width="double",
        height_in=6.3,
        constrained_layout=True,
    )
    draw_entropy_tertile_profiles(ax, table)
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
