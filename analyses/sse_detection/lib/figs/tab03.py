"""Build bayesian fixed effects focal table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from utils import table_caption, write_latex_longtable

from ..sse.config import BAYESIAN_OUTPUT_DIR
from .common import Paths, add_common_args, latex_table_path, paths_from_args
from .table_common import (
    _format_direction,
    _format_effect_interval,
    _read_summary_table,
    _sort_for_report,
    _write_data_table,
)

TABLE_NAME = "tab03_bayesian_fixed_effects_focal"


def build_fixed_effects_focal_table(
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
) -> pd.DataFrame:
    """Return focal fixed effects from expanded interpretation models."""
    table = _read_summary_table(result_dir, "estimates")
    scale = table["Scale"].fillna("")
    focal = table["Term Type"].isin(["Categorical Contrast", "Mixing Entropy"])
    interpretation_model = table["Model"].eq("Expanded") & (
        table["Domain"].eq("Composition")
        | (table["Domain"].eq("Mixing") & scale.eq("Null Standardised"))
    )
    out = table.loc[focal & interpretation_model].copy()
    return _sort_for_report(out, extra_columns=["Parameter"])


def build(paths: Paths) -> dict[str, Path]:
    table = build_fixed_effects_focal_table(paths.bayesian_result_dir)
    name = TABLE_NAME
    _write_data_table(table, paths.result_table_dir, name)
    rows = [
        [
            row["Domain"],
            row["Family"],
            row["Outcome"],
            row["Scale"],
            row["Model"],
            row["Parameter"],
            _format_effect_interval(row),
            _format_direction(row),
            row["Diagnostic Status"],
        ]
        for _, row in table.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    write_latex_longtable(
        tex_path,
        **table_caption("tab:bayesian_fixed_effects_focal"),
        columns=[
            "Domain",
            "Family",
            "Outcome",
            "Scale",
            "Model",
            "Term",
            "Effect (95% HDI)",
            "Direction",
            "Diagnostic",
        ],
        rows=rows,
        column_spec=(
            r"P{0.065\linewidth}P{0.055\linewidth}P{0.075\linewidth}"
            r"P{0.090\linewidth}P{0.060\linewidth}P{0.255\linewidth}"
            r"P{0.120\linewidth}P{0.120\linewidth}P{0.070\linewidth}"
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
