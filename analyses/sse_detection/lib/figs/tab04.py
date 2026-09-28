"""Build bayesian fixed effects intercepts table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from utils import write_latex_grouped_column_table

from ..sse.config import BAYESIAN_OUTPUT_DIR
from .common import Paths, add_common_args, latex_table_path, paths_from_args
from .table_common import (
    _display_text,
    _format_float,
    _read_summary_table,
    _write_data_table,
)

TABLE_NAME = "tab_bayesian_fixed_effects_intercepts"


def _intercept_digits(row: pd.Series) -> int:
    values = [
        float(value)
        for value in (
            row.get("Estimate"),
            row.get("HDI 95 Low"),
            row.get("HDI 95 High"),
        )
        if not pd.isna(value)
    ]
    if (
        "odds" in _display_text(row.get("Effect Scale")).lower()
        and max(
            (abs(value) for value in values),
            default=0.0,
        )
        < 0.01
    ):
        return 4
    return 3


def _format_intercept_interval(row: pd.Series) -> str:
    digits = _intercept_digits(row)
    estimate = _format_float(row.get("Estimate"), digits)
    low = _format_float(row.get("HDI 95 Low"), digits)
    high = _format_float(row.get("HDI 95 High"), digits)
    if not estimate or not low or not high:
        return estimate
    return f"{estimate} [{low}, {high}]"


def _intercept_lookup_key(row: pd.Series) -> tuple[str, str, str, str]:
    return (
        _display_text(row.get("Domain")),
        _display_text(row.get("Outcome")),
        _display_text(row.get("Scale")),
        _display_text(row.get("Model")),
    )


def _intercept_lookup(
    table: pd.DataFrame,
) -> dict[tuple[str, str, str, str], pd.Series]:
    lookup = {}
    for _, row in table.iterrows():
        key = _intercept_lookup_key(row)
        if key in lookup:
            raise ValueError(f"Duplicate intercept rows for {key}.")
        lookup[key] = row
    return lookup


def _intercept_cell(
    lookup: dict[tuple[str, str, str, str], pd.Series],
    *,
    domain: str,
    outcome: str,
    scale: str,
    model: str,
) -> str:
    key = (domain, outcome, scale, model)
    row = lookup.get(key)
    if row is None:
        raise ValueError(f"Missing intercept row for {key}.")
    return _format_intercept_interval(row)


def build_fixed_effects_intercepts_table(
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
) -> pd.DataFrame:
    """Return compact intercept estimates for primary and expanded models."""
    table = _read_summary_table(result_dir, "estimates")
    intercepts = table.loc[table["Term Type"].eq("Intercept")].copy()
    lookup = _intercept_lookup(intercepts)
    rows = []
    for outcome in ("Candidate", "Burst Score", "Burden Score"):
        for scale_idx, scale in enumerate(("Observed", "Null Standardised")):
            show_composition = scale_idx == 0
            rows.append(
                {
                    "Outcome": outcome if scale_idx == 0 else "",
                    "Mixing Scale": scale,
                    "Mixing Primary": _intercept_cell(
                        lookup,
                        domain="Mixing",
                        outcome=outcome,
                        scale=scale,
                        model="Primary",
                    ),
                    "Mixing Expanded": _intercept_cell(
                        lookup,
                        domain="Mixing",
                        outcome=outcome,
                        scale=scale,
                        model="Expanded",
                    ),
                    "Composition Primary": (
                        _intercept_cell(
                            lookup,
                            domain="Composition",
                            outcome=outcome,
                            scale="",
                            model="Primary",
                        )
                        if show_composition
                        else ""
                    ),
                    "Composition Expanded": (
                        _intercept_cell(
                            lookup,
                            domain="Composition",
                            outcome=outcome,
                            scale="",
                            model="Expanded",
                        )
                        if show_composition
                        else ""
                    ),
                }
            )
    return pd.DataFrame(rows)


def build(paths: Paths) -> dict[str, Path]:
    table = build_fixed_effects_intercepts_table(paths.bayesian_result_dir)
    name = TABLE_NAME
    _write_data_table(table, paths.result_table_dir, name)
    value_columns = [
        "Mixing Primary",
        "Mixing Expanded",
        "Composition Primary",
        "Composition Expanded",
    ]
    rows = [
        [
            row["Outcome"],
            row["Mixing Scale"],
            *[row[column] for column in value_columns],
        ]
        for _, row in table.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    write_latex_grouped_column_table(
        tex_path,
        caption=(
            "Intercept estimates from primary and expanded Bayesian "
            "characterisation models. Values are posterior means with 95% HDIs. CSC "
            "intercepts are baseline odds; intercepts for the local burst and onward "
            "burden models are expected scores. These baselines correspond to zero "
            "continuous predictors and group deviations, with categorical predictors at their "
            "reference levels. Mixing intercepts are shown separately for observed "
            "and null-standardised entropy specifications; composition intercepts "
            "are shown once per outcome."
        ),
        short_caption=(
            "Intercept estimates from primary and expanded Bayesian models."
        ),
        label="tab:bayesian_fixed_effects_intercepts",
        row_columns=["Outcome", "Mixing scale"],
        column_groups=[
            ("Mixing Model Intercept", ["Primary", "Expanded"]),
            ("Composition Model Intercept", ["Primary", "Expanded"]),
        ],
        rows=rows,
        column_spec=(
            r"P{0.105\linewidth}P{0.120\linewidth}"
            r"P{0.135\linewidth}P{0.135\linewidth}"
            r"P{0.135\linewidth}P{0.135\linewidth}"
        ),
        addlinespace_after={idx for idx in range(1, len(rows) - 1, 2)},
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
