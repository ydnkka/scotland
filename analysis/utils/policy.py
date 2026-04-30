"""Scotland COVID-19 policy period helpers."""

from __future__ import annotations

import polars as pl

_PERIOD_CODES = [
    "E0", "L1", "P1", "P2", "P3", "T1", "F5", "L2",
    "SL", "L3", "L21", "L0", "NN", "OM", "FE", "PR",
]

POLICY_PERIODS = pl.DataFrame({
    "period_code": [
        "E0", "L1", "P1", "P2", "P3", "T1", "F5", "L2",
        "SL", "L3", "L21", "L0", "NN", "OM", "FE", "PR",
    ],
    "period_label": [
        "Emergence", "First lockdown", "Route map phase 1", "Route map phase 2",
        "Route map phase 3", "Pre-tier tightening", "Five-tier framework",
        "Second lockdown", "Stay local — Level 3", "Level 3", "Level 2 / Level 1",
        "Level 0", "Near-normal", "Omicron wave", "Final easing", "Post-restriction",
    ],
    "start_date": [
        "2020-03-01", "2020-03-24", "2020-05-29", "2020-06-19", "2020-07-10",
        "2020-10-02", "2020-11-02", "2021-01-05", "2021-04-02", "2021-04-26",
        "2021-05-17", "2021-07-19", "2021-08-09", "2021-11-29", "2022-01-24", "2022-04-18",
    ],
    "end_date": [
        "2020-03-23", "2020-05-28", "2020-06-18", "2020-07-09", "2020-10-01",
        "2020-11-01", "2021-01-04", "2021-04-01", "2021-04-25", "2021-05-16",
        "2021-07-18", "2021-08-08", "2021-11-28", "2022-01-23", "2022-04-17", "2023-05-05",
    ],
    "intensity": [15, 100, 72, 52, 30, 55, 65, 95, 65, 55, 38, 20, 10, 42, 15, 3],
}).with_columns([
    pl.col("start_date").str.to_date(),
    pl.col("end_date").str.to_date(),
    pl.col("period_code").cast(pl.Enum(_PERIOD_CODES)),
])


def assign_period(dates: pl.Series) -> pl.Series:
    """Assign a policy period code to each date in the series."""
    date_rows = pl.DataFrame({
        "_row_id": pl.int_range(0, len(dates), eager=True),
        "_date": dates.cast(pl.Date),
    })

    matches = (
        date_rows
        .join(POLICY_PERIODS.select(["period_code", "start_date", "end_date"]), how="cross")
        .filter(pl.col("_date").is_between(pl.col("start_date"), pl.col("end_date")))
        .select(["_row_id", "period_code"])
    )

    return (
        date_rows
        .select("_row_id")
        .join(matches, on="_row_id", how="left")
        .sort("_row_id")
        .get_column("period_code")
    )


def attach_period(df: pl.DataFrame, date_col: str) -> pl.DataFrame:
    """Attach policy period code, label, and intensity using the given date column."""
    period_lookup = POLICY_PERIODS.select([
        pl.col("period_code").alias("policy_period"),
        pl.col("period_label").alias("policy_period_label"),
        pl.col("intensity").alias("policy_intensity"),
    ])

    return (
        df
        .with_columns(assign_period(df[date_col]).alias("policy_period"))
        .join(period_lookup, on="policy_period", how="left")
    )
