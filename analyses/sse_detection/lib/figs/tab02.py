"""Build bayesian model diagnostics table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from utils import write_latex_longtable

from ..sse.config import BAYESIAN_OUTPUT_DIR
from .common import Paths, add_common_args, latex_table_path, paths_from_args
from .table_common import (
    _format_float,
    _format_int,
    _read_summary_table,
    _sort_for_report,
    _write_data_table,
)

TABLE_NAME = "tab02_bayesian_model_diagnostics"


def build_model_diagnostics_table(
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
) -> pd.DataFrame:
    """Return the all-domain diagnostic summary table."""
    return _sort_for_report(_read_summary_table(result_dir, "diagnostics"))


def build(paths: Paths) -> dict[str, Path]:
    table = build_model_diagnostics_table(paths.bayesian_result_dir)
    name = TABLE_NAME
    _write_data_table(table, paths.result_table_dir, name)
    rows = [
        [
            row["Domain"],
            row["Family"],
            row["Outcome"],
            row["Scale"],
            row["Model"],
            row["Diagnostic Status"],
            _format_int(row["Divergences"]),
            _format_float(row["Min BFMI"], 3),
            _format_float(row["Max Rhat"], 3),
            _format_int(row["Min Bulk ESS"]),
            _format_int(row["Min Tail ESS"]),
            _format_int(row["Max Tree Depth"]),
        ]
        for _, row in table.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    write_latex_longtable(
        tex_path,
        caption=(
            "Sampling diagnostics for all 18 Bayesian characterisation models. "
            "Candidate denotes CSC status. OK indicates that all specified sampling "
            "criteria were met; Warning indicates at least one failed criterion. "
            "Div. counts post-tuning divergent transitions; BFMI is Bayesian "
            "fraction of missing information, Rhat is $\\widehat R$, and ESS is "
            "effective sample size. Definitions and thresholds are given in "
            "Table~\\ref{tab:app_bayes_diagnostic_guide}; model specifications are "
            "defined in Table~\\ref{tab:ch6_model_terminology}."
        ),
        caption_is_latex=True,
        short_caption=(
            "Bayesian characterisation sampling diagnostics across composition and "
            "mixing models."
        ),
        label="tab:bayesian_model_diagnostics",
        columns=[
            "Domain",
            "Family",
            "Outcome",
            "Scale",
            "Model",
            "Status",
            "Div.",
            "Min BFMI",
            "Max Rhat",
            "Min bulk ESS",
            "Min tail ESS",
            "Max tree depth",
        ],
        rows=rows,
        column_spec=(
            r"P{0.070\linewidth}P{0.055\linewidth}P{0.075\linewidth}"
            r"P{0.090\linewidth}P{0.060\linewidth}P{0.060\linewidth}"
            r"P{0.045\linewidth}P{0.055\linewidth}P{0.055\linewidth}"
            r"P{0.065\linewidth}P{0.065\linewidth}P{0.060\linewidth}"
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
