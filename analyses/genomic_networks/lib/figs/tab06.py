"""Build assortativity variance decomposition table."""

from __future__ import annotations

import argparse
from collections.abc import Callable
from typing import Any

import pandas as pd

from utils import write_latex_table

from .assortativity_analysis import (
    compatibility_variance_decomposition_long,
    variance_decomposition_summary,
)
from .common import ATTRIBUTE_ORDER, Paths, add_common_args, paths_from_args, read_table
from .table_common import (
    fmt_float,
    fmt_int,
)

TABLE_NAME = "tab06_assortativity_variance_decomposition"

ATTRIBUTE_COLUMN_LABELS = {
    "SIMD quintile": "SIMD",
    "Urban/rural class": "Urban/rural class",
    "Health board": "HB",
    "Local authority": "LA",
}


def variance_decomposition_summary_table(paths: Paths) -> pd.DataFrame:
    try:
        return read_table(paths, "compatibility_variance_decomposition_summary")
    except FileNotFoundError:
        vd_long = compatibility_variance_decomposition_long(paths)
        return variance_decomposition_summary(vd_long)


def build(paths: Paths) -> None:
    summary = variance_decomposition_summary_table(paths)
    attribute_rows = {}
    for attribute in ATTRIBUTE_ORDER:
        group = summary.loc[summary["attribute_label"].eq(attribute)]
        if not group.empty:
            attribute_rows[attribute] = group.iloc[0]

    metrics: list[tuple[str, str, Callable[[Any], str]]] = [
        ("Rows", "n", fmt_int),
        ("Windows", "n_windows", fmt_int),
        ("Lineages", "n_lineages", fmt_int),
        (
            "Additive model",
            "additive_model_fraction",
            lambda value: fmt_float(value, 3),
        ),
        (
            "Adj. additive",
            "adj_additive_model_fraction",
            lambda value: fmt_float(value, 3),
        ),
        ("Residual", "residual_fraction", lambda value: fmt_float(value, 3)),
        (
            "Window alone",
            "window_alone_fraction",
            lambda value: fmt_float(value, 3),
        ),
        (
            "Lineage alone",
            "lineage_alone_fraction",
            lambda value: fmt_float(value, 3),
        ),
        (
            "Lineage | Window",
            "lineage_given_window_fraction",
            lambda value: fmt_float(value, 3),
        ),
        (
            "Window | Lineage",
            "window_given_lineage_fraction",
            lambda value: fmt_float(value, 3),
        ),
        ("Weight cap", "weight_cap", lambda value: fmt_float(value, 2)),
        ("Capped weights", "n_weights_capped", fmt_int),
        ("Median boot SE", "median_boot_se", lambda value: fmt_float(value, 3)),
        (
            "Median CI width",
            "median_ci_width",
            lambda value: fmt_float(value, 3),
        ),
    ]

    attributes = list(attribute_rows)
    rows = []
    for metric_label, column, formatter in metrics:
        rows.append(
            [
                metric_label,
                *[
                    formatter(attribute_rows[attribute][column])
                    for attribute in attributes
                ],
            ]
        )

    write_latex_table(
        paths.publication_table_dir / f"{TABLE_NAME}.tex",
        caption=(
            "Variance decomposition of weighted compatibility assortativity estimates after "
            "applying a 90th percentile inverse-variance weight cap. "
            "SIMD, population-weighted SIMD quintile; HB, health board; "
            "LA, local authority. Variance components and fractions are unitless proportions; "
            "counts and uncertainty summaries retain their native units."
        ),
        short_caption="Variance decomposition of compatibility assortativity",
        label="tab:assortativity_variance_decomposition",
        columns=[
            "Metric",
            *[
                ATTRIBUTE_COLUMN_LABELS.get(attribute) or attribute
                for attribute in attributes
            ],
        ],
        rows=rows,
        column_spec="l" + "r" * len(attributes),
        addlinespace_after={2, 9},
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
