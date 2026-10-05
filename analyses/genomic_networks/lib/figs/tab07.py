"""Build simd population weighting table."""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np

from utils import write_latex_table

from .common import Paths, add_common_args, paths_from_args, read_table
from .table_common import (
    fmt_int,
)

TABLE_NAME = "tab07_simd_population_weighting"

SIMD_GROUP_LABELS = {
    5: "quintile",
    10: "decile",
    20: "vigintile",
}


def fmt_percent_points(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{float(value):.1f}%"


def build(paths: Paths) -> None:
    """Write the SIMD validation LaTeX appendix table."""
    group_summary = read_table(paths, "simd_population_weighting_group_summary")
    n_groups = 5
    group_name = SIMD_GROUP_LABELS.get(n_groups, f"{n_groups}-group")
    caption = (
        f"Validation of national SIMD {group_name} boundaries comparing equal-Data-Zone "
        "and population-weighted grouping methodologies. The table details the number of "
        "Data Zones, total population counts, proportional population shares, and the "
        "underlying SIMD rank ranges represented within each constructed group."
    )
    short_caption = f"SIMD {group_name} population-weighting validation"

    methods = ["equal_datazone", "population_weighted"]
    display = group_summary.loc[group_summary["grouping_method"].isin(methods)].copy()
    method_order = {method: idx for idx, method in enumerate(methods)}
    display["_method_sort"] = display["grouping_method"].map(method_order)
    display = display.sort_values(["_method_sort", "simd_group"])
    display["rank_range"] = (
        display["first_simd_rank"].map(fmt_int)
        + "--"
        + display["last_simd_rank"].map(fmt_int)
    )
    rows = []
    addlinespace_after: set[int] = set()
    method_groups = list(display.groupby("_method_sort", sort=False))
    for method_idx, (_, group) in enumerate(method_groups):
        group = group.sort_values("simd_group")
        for row_idx, row in enumerate(group.itertuples(index=False)):
            rows.append(
                [
                    str(row.grouping_method_label) if row_idx == 0 else "",
                    row.simd_group,
                    fmt_int(row.n_datazones),
                    fmt_int(row.total_population),
                    fmt_percent_points(row.pct_population),
                    str(row.rank_range),
                ]
            )
        if method_idx < len(method_groups) - 1 and rows:
            addlinespace_after.add(len(rows) - 1)

    write_latex_table(
        paths.publication_table_dir / f"{TABLE_NAME}.tex",
        caption=caption,
        short_caption=short_caption,
        label="tab:simd_population_weighting",
        columns=[
            "Method",
            "Group",
            "Zones",
            "Pop.",
            "Share",
            "Ranks",
        ],
        rows=rows,
        column_spec="lrrrrl",
        addlinespace_after=addlinespace_after,
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
