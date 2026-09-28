"""Build bayesian random effect sds table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from utils import write_latex_longtable

from ..sse.config import BAYESIAN_OUTPUT_DIR
from .common import Paths, add_common_args, latex_table_path, paths_from_args
from .table_common import (
    _format_effect_interval,
    _format_int,
    _format_percent,
    _read_summary_table,
    _sort_for_report,
    _write_data_table,
)

TABLE_NAME = "tab_bayesian_random_effect_sds"


def build_random_effects_table(result_dir: Path = BAYESIAN_OUTPUT_DIR) -> pd.DataFrame:
    """Return random-effect and residual-SD summaries."""
    return _sort_for_report(
        _read_summary_table(result_dir, "random_effects"),
        extra_columns=["Grouping Factor", "Parameter"],
    )


def build(paths: Paths) -> dict[str, Path]:
    table = build_random_effects_table(paths.bayesian_result_dir)
    name = TABLE_NAME
    _write_data_table(table, paths.result_table_dir, name)
    rows = [
        [
            row["Domain"],
            row["Family"],
            row["Outcome"],
            row["Scale"],
            row["Model"],
            row["Component Type"],
            row["Grouping Factor"],
            _format_effect_interval(row),
            row["Effect Scale"],
            _format_percent(row["Random Effect Variance Share"]),
            _format_int(row["Random Effect SD Rank"]),
            row["Diagnostic Status"],
        ]
        for _, row in table.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    write_latex_longtable(
        tex_path,
        caption=(
            "Random-effect and residual standard-deviation summaries from "
            "Bayesian characterisation models. Variance share and rank are shown "
            "for random-effect standard deviations only."
        ),
        short_caption="Bayesian random-effect and residual standard deviations.",
        label="tab:bayesian_random_effect_sds",
        columns=[
            "Domain",
            "Family",
            "Outcome",
            "Scale",
            "Model",
            "Component",
            "Group",
            "SD (95% HDI)",
            "Effect scale",
            "Variance share",
            "Rank",
            "Diagnostic",
        ],
        rows=rows,
        column_spec=(
            r"P{0.065\linewidth}P{0.055\linewidth}P{0.075\linewidth}"
            r"P{0.090\linewidth}P{0.060\linewidth}P{0.080\linewidth}"
            r"P{0.075\linewidth}P{0.105\linewidth}P{0.115\linewidth}"
            r"P{0.065\linewidth}P{0.045\linewidth}P{0.060\linewidth}"
        ),
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
