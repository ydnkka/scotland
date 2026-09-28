"""Build assortativity summary table."""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np
import pandas as pd

from utils import write_latex_table

from .assortativity_analysis import (
    compatibility_window_pooled_meta,
    pooled_window_attribute_summary,
)
from .common import ATTRIBUTE_ORDER, Paths, add_common_args, paths_from_args, read_table
from .table_common import (
    fmt_float,
    fmt_int,
)

TABLE_NAME = "tab_assortativity_summary"


def fmt_ci(low: Any, high: Any, digits: int = 2) -> str:
    if (
        low is None
        or high is None
        or (isinstance(low, float) and np.isnan(low))
        or (isinstance(high, float) and np.isnan(high))
    ):
        return "not estimated"
    return f"{fmt_float(low, digits)}--{fmt_float(high, digits)}"


def pooled_window_attribute_summary_table(paths: Paths) -> pd.DataFrame:
    try:
        summary = read_table(paths, "compatibility_window_pooled_summary")
        if {"gls_mean", "gls_ci_low", "gls_ci_high"}.issubset(summary.columns):
            return summary
    except FileNotFoundError:
        pass
    window_meta, _ = compatibility_window_pooled_meta(paths)
    return pooled_window_attribute_summary(window_meta)


def build(paths: Paths) -> None:
    summary = pooled_window_attribute_summary_table(paths)
    rows = []
    for attribute in ATTRIBUTE_ORDER:
        group = summary.loc[summary["attribute_label"].eq(attribute)]
        if group.empty:
            continue
        row = group.iloc[0]
        mean_col = "gls_mean" if "gls_mean" in row.index else "weighted_mean"
        low_col = "gls_ci_low" if "gls_ci_low" in row.index else "ci_low"
        high_col = "gls_ci_high" if "gls_ci_high" in row.index else "ci_high"
        rows.append(
            [
                attribute,
                fmt_int(row["n_windows"]),
                fmt_int(row["n_estimated_windows"]),
                fmt_float(row["median_n_lineages"], 1),
                (
                    f"{fmt_float(row[mean_col], 4)} "
                    f"({fmt_ci(row[low_col], row[high_col], 4)})"
                ),
                (
                    f"{fmt_float(row['window_median'], 3)} "
                    f"({fmt_float(row['window_q10'], 3)}--"
                    f"{fmt_float(row['window_q90'], 3)})"
                ),
            ]
        )
    write_latex_table(
        paths.publication_table_dir / f"{TABLE_NAME}.tex",
        caption=(
            "Pooled compatibility assortativity estimates calculated across rolling windows. "
            "For each attribute, the table reports the total number of windows, the subset "
            "yielding an estimable coefficient, and the median number of lineages analyzed. "
            "Mean r is the single overall pooled estimate: the GLS mean of the window-level "
            "pooled assortativity estimates, adjusted for covariance induced by overlapping "
            "rolling windows, with its 95% confidence interval. Window r is descriptive: "
            "the median of the same window-level pooled estimates, followed by their "
            "10th--90th percentile range."
        ),
        short_caption="Pooled compatibility assortativity estimates",
        label="tab:assortativity_summary",
        columns=[
            "Attribute",
            "Windows",
            "Est.",
            "Lineages",
            "Mean r (95% CI)",
            "Window r (10th--90th)",
        ],
        rows=rows,
        column_spec="lrrrrll",
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
