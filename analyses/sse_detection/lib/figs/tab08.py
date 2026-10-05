"""Build sensitivity sizes table."""

from __future__ import annotations

import argparse

from .common import Paths, add_common_args, paths_from_args
from .size_profiles import (
    SCENARIO_LABELS,
    SCENARIO_ORDER,
    load_saved_summaries,
    prepare,
)
from .table_common import write_size_table

TABLE_NAME = "tab08_sensitivity_sizes"


def build(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    size = tables["size_summary"]
    output_path = paths.publication_table_dir / f"{TABLE_NAME}.tex"
    rows = []
    for scenario, label in zip(SCENARIO_ORDER, SCENARIO_LABELS):
        for group in ("CSC", "Background"):
            r = size.loc[
                size.scenario.eq(scenario)
                & size.size_band.eq("All")
                & size.group.eq(group)
            ].iloc[0]
            rows.append(
                [
                    label.replace("≤", " <= ").replace("≥", " >= "),
                    group,
                    f"{r.n_clusters:,}",
                    f"{r.n_records:,}",
                    f"{r.size_median:g}",
                ]
            )
    write_size_table(
        output_path,
        "Descriptive sensitivity comparisons using saved scores. Primary uses $p\\leq0.05$ and size $\\geq6$. Alternative p cutoffs retain the same testing floor; size restrictions retain primary labels without recalibration.",
        "tab:sensitivity_sizes",
        ["Comparison", "Group", "Clusters", "Records", "Median size"],
        rows,
    )
    return {"tex": output_path}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument(
        "--from-saved-tables",
        action="store_true",
        help="Render saved size aggregates without rebuilding or writing them.",
    )
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths, tables=load_saved_summaries(paths) if args.from_saved_tables else None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
