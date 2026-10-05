"""Build bayesian fixed effects full table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from utils import write_latex_longtable

from ..sse.config import BAYESIAN_OUTPUT_DIR
from .common import Paths, add_common_args, latex_table_path, paths_from_args
from .table_common import (
    _format_direction,
    _format_effect_interval,
    _format_probability,
    _read_summary_table,
    _sort_for_report,
    _write_data_table,
)

TABLE_NAME = "tab05_bayesian_fixed_effects_full"


def build_fixed_effects_full_table(
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
) -> pd.DataFrame:
    """Return the complete estimate summary table."""
    return _sort_for_report(
        _read_summary_table(result_dir, "estimates"),
        extra_columns=["Parameter"],
    )


def build(paths: Paths) -> dict[str, Path]:
    table = build_fixed_effects_full_table(paths.bayesian_result_dir)
    name = TABLE_NAME
    _write_data_table(table, paths.result_table_dir, name)
    rows = []
    for _, row in table.iterrows():
        domain_family = f"{row['Domain']}; {row['Family']}"
        rows.append(
            [
                domain_family,
                row["Outcome"],
                row["Scale"],
                row["Model"],
                row["Parameter"],
                row["Term Type"],
                _format_effect_interval(row),
                row["Effect Scale"],
                _format_probability(row["P Positive Direction"]),
                _format_probability(row["P Negative Direction"]),
                _format_direction(row),
                row["Direction Band"],
                row["Diagnostic Status"],
            ]
        )
    tex_path = latex_table_path(paths, name)
    write_latex_longtable(
        tex_path,
        caption=(
            "Complete Bayesian characterisation estimate summary. This table "
            "contains all rows from summary table 2, including intercepts, "
            "focal fixed effects, continuous adjusters, and random-intercept "
            "level summaries."
        ),
        short_caption="Complete Bayesian characterisation estimate summary.",
        label="tab:bayesian_fixed_effects_full",
        columns=[
            "Domain; family",
            "Outcome",
            "Scale",
            "Model",
            "Parameter",
            "Term type",
            "Effect (95% HDI)",
            "Effect scale",
            "P(+)",
            "P(-)",
            "Direction",
            "Band",
            "Status",
        ],
        rows=rows,
        column_spec=(
            r"P{0.070\linewidth}P{0.060\linewidth}P{0.065\linewidth}"
            r"P{0.050\linewidth}P{0.185\linewidth}P{0.070\linewidth}"
            r"P{0.095\linewidth}P{0.065\linewidth}P{0.035\linewidth}"
            r"P{0.035\linewidth}P{0.075\linewidth}P{0.050\linewidth}"
            r"P{0.040\linewidth}"
        ),
        tiny=True,
    )
    return {
        "csv": paths.result_table_dir / f"{name}.csv",
        "parquet": paths.result_table_dir / f"{name}.parquet",
        "tex": tex_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
