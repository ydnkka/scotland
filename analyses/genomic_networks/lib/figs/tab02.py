"""Build policy denominators table."""

from __future__ import annotations

import argparse

import pandas as pd

from utils import addlinespace_after_group_changes, write_latex_table

from .common import Paths, add_common_args, paths_from_args, read_table, sort_by_policy
from .table_common import (
    fmt_int,
    fmt_percent,
)

TABLE_NAME = "tab_policy_denominators"


def build(paths: Paths) -> None:
    denominators = read_table(paths, "window_denominator_contrasts")
    has_sequences = pd.to_numeric(
        denominators["median_window_sequences"], errors="coerce"
    ).gt(0)
    denominators = denominators.loc[
        has_sequences | denominators["policy_period"].astype(str).eq("P2")
    ]
    denominators = sort_by_policy(denominators, column="policy_period")
    rows = []
    for row in denominators.itertuples(index=False):
        rows.append(
            [
                str(row.policy_era).capitalize().replace("_", " "),
                row.policy_period,
                fmt_int(row.n_windows),
                fmt_int(row.median_window_positive_tests),
                fmt_int(row.median_window_sequences),
                f"{fmt_percent(row.median_window_prop_sequenced)} ({fmt_percent(row.min_window_prop_sequenced)}--{fmt_percent(row.max_window_prop_sequenced)})",
            ]
        )
    write_latex_table(
        paths.publication_table_dir / f"{TABLE_NAME}.tex",
        caption=(
            "Rolling-window observation denominators stratified by epidemic era and "
            "policy period. For each period, the table outlines the total number of "
            "analytical windows alongside the median window-level counts of confirmed positive "
            "tests and sequenced genomes. Sequencing coverage (the percentage of "
            "confirmed positive tests successfully sequenced) is reported as a median, with "
            "absolute minimum and maximum extremes presented in parentheses."
        ),
        short_caption="Rolling-window observation denominators by policy period",
        label="tab:policy_denominators",
        columns=[
            "Epidemic era",
            "Period",
            "Windows",
            "Tests",
            "Genomes",
            "Coverage",
        ],
        rows=rows,
        column_spec="llrrrl",
        addlinespace_after=addlinespace_after_group_changes(denominators["policy_era"]),
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
