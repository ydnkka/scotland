"""Build sequence composition by policy table."""

from __future__ import annotations

import argparse
from typing import Any

import pandas as pd

from utils import write_latex_grouped_column_table

from .common import (
    POLICY_LABELS,
    Paths,
    add_common_args,
    paths_from_args,
    read_table,
    sort_by_policy,
)
from .table_common import (
    fmt_int,
    fmt_percent,
)

TABLE_NAME = "tab01_sequence_composition_by_policy"

SEQUENCE_COMPOSITION_ATTRIBUTES = (
    "sex",
    "simd_quintile",
    "age_group",
    "urban_rural",
    "health_board",
)

SEQUENCE_COMPOSITION_CATEGORY_LABELS = {
    "simd_quintile": {
        "1.0": "Q1",
        "2.0": "Q2",
        "3.0": "Q3",
        "4.0": "Q4",
        "5.0": "Q5",
        "1": "Q1",
        "2": "Q2",
        "3": "Q3",
        "4": "Q4",
        "5": "Q5",
    },
}


def sequence_composition_category_order(
    table: pd.DataFrame,
    attribute: str,
) -> list[str]:
    work = table.loc[table["attribute"].eq(attribute)]
    totals = work.groupby("category", observed=False)["n_sequences"].sum()
    if totals.empty:
        return []

    if attribute == "sex":
        preferred = ["Female", "Male"]
        return [value for value in preferred if value in totals.index.astype(str)]

    if attribute == "simd_quintile":
        return sorted(totals.index.astype(str), key=lambda value: float(value))

    if attribute == "age_group":
        preferred = ["00-04", "05-14", "15-24", "25-64", "65-74", "75+"]
        return [value for value in preferred if value in totals.index.astype(str)]

    if attribute == "urban_rural":
        preferred = [
            "Large Urban Areas",
            "Other Urban Areas",
            "Accessible Small Towns",
            "Remote Small Towns",
            "Accessible Rural",
            "Remote Rural",
        ]
        return [value for value in preferred if value in totals.index.astype(str)]

    if attribute == "health_board":
        return totals.sort_values(ascending=False).index.astype(str).tolist()

    return totals.sort_values(ascending=False).index.astype(str).tolist()


def display_sequence_category(attribute: str, category: object) -> str:
    text = str(category)
    return SEQUENCE_COMPOSITION_CATEGORY_LABELS.get(attribute, {}).get(text, text)


def format_sequence_composition_cell(row: pd.Series | None) -> str:
    if row is None:
        return "0 (0.0%)"
    n_sequences = pd.to_numeric(pd.Series([row["n_sequences"]]), errors="coerce").iloc[
        0
    ]
    if pd.notna(n_sequences) and float(n_sequences) == 0:
        return "0 (0.0%)"
    if bool(row.get("small_cell", False)):
        return "<5"
    return f"{fmt_int(n_sequences)} ({fmt_percent(row['proportion'])})"


def sequence_composition_policy_groups(
    composition: pd.DataFrame,
) -> tuple[pd.DataFrame, list[tuple[str, list[str]]]]:
    periods = composition[["policy_era", "policy_period"]].drop_duplicates()
    periods = sort_by_policy(periods, column="policy_period")

    column_groups: list[tuple[str, list[str]]] = []
    for era, group in periods.groupby("policy_era", sort=False, dropna=False):
        label = POLICY_LABELS.get(str(era), str(era).capitalize().replace("_", " "))
        column_groups.append((label, group["policy_period"].astype(str).tolist()))

    return periods, column_groups


def build(paths: Paths) -> None:
    composition = read_table(paths, "sequence_composition_by_policy")
    composition = composition.loc[
        composition["attribute"].isin(SEQUENCE_COMPOSITION_ATTRIBUTES)
    ].copy()
    if composition.empty:
        raise FileNotFoundError(
            "sequence_composition_by_policy contains no figure attributes"
        )

    periods, column_groups = sequence_composition_policy_groups(composition)
    period_codes = periods["policy_period"].astype(str).tolist()

    lookup = {
        (str(row["attribute"]), str(row["category"]), str(row["policy_period"])): row
        for row in composition.to_dict(orient="records")
    }

    rows: list[list[Any]] = []
    addlinespace_after: set[int] = set()
    for attribute_idx, attribute in enumerate(SEQUENCE_COMPOSITION_ATTRIBUTES):
        attr_rows = composition.loc[composition["attribute"].eq(attribute)]
        if attr_rows.empty:
            continue
        attribute_label = str(attr_rows["attribute_label"].dropna().iloc[0])
        categories = sequence_composition_category_order(composition, attribute)
        for category_idx, category in enumerate(categories):
            rows.append(
                [
                    attribute_label if category_idx == 0 else "",
                    display_sequence_category(attribute, category),
                    *[
                        format_sequence_composition_cell(
                            pd.Series(lookup.get((attribute, str(category), period)))
                            if lookup.get((attribute, str(category), period))
                            is not None
                            else None
                        )
                        for period in period_codes
                    ],
                ]
            )
        if attribute_idx < len(SEQUENCE_COMPOSITION_ATTRIBUTES) - 1 and rows:
            addlinespace_after.add(len(rows) - 1)

    write_latex_grouped_column_table(
        paths.publication_table_dir / f"{TABLE_NAME}.tex",
        caption=(
            "Sequence composition by epidemic era and policy period for demographic, socioeconomic, or geographic variables. "
            "Columns are grouped by epidemic era, with policy-period codes shown as subcolumns. "
            "Cells report the number and percentage of unique sequences in that policy "
            "period with the corresponding attribute category. Counts from one to four are "
            "suppressed and shown as <5."
        ),
        short_caption="Sequence composition by epidemic era and policy period",
        label="tab:sequence_composition_by_policy",
        row_columns=["Attribute", "Category"],
        column_groups=column_groups,
        rows=rows,
        column_spec=f"ll*{{{len(period_codes)}}}{{r}}",
        addlinespace_after=addlinespace_after,
        landscape=True,
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
