"""Policy-period definitions and helpers for Scottish COVID-19 restrictions.

The exported period table assigns ordered restriction phases to calendar dates
and exposes compact labels/intensity scores for plotting and modelling.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


__all__ = [
    "PERIOD_ORDER",
    "POLICY_PERIODS",
    "PERIOD_LABELS",
    "PERIOD_INTENSITY",
    "assign_period",
    "attach_period",
]


_PERIOD_CODES = [
    "E0",
    "L1",
    "P1",
    "P2",
    "P3",
    "T1",
    "F5",
    "L2",
    "SL",
    "L3",
    "L21",
    "L0",
    "NN",
    "OM",
    "FE",
    "PR",
]

PERIOD_ORDER: list[str] = _PERIOD_CODES

POLICY_PERIODS: pd.DataFrame = pd.DataFrame(
    {
        "period_code": [
            "E0",
            "L1",
            "P1",
            "P2",
            "P3",
            "T1",
            "F5",
            "L2",
            "SL",
            "L3",
            "L21",
            "L0",
            "NN",
            "OM",
            "FE",
            "PR",
        ],
        "period_label": [
            "Emergence",
            "First lockdown",
            "Route map phase 1",
            "Route map phase 2",
            "Route map phase 3",
            "Pre-tier tightening",
            "Five-tier framework",
            "Second lockdown",
            "Stay local — Level 3",
            "Level 3",
            "Level 2 / Level 1",
            "Level 0",
            "Near-normal",
            "Omicron wave",
            "Final easing",
            "Post-restriction",
        ],
        "start_date": pd.to_datetime(
            [
                "2020-03-01",
                "2020-03-24",
                "2020-05-29",
                "2020-06-19",
                "2020-07-10",
                "2020-10-02",
                "2020-11-02",
                "2021-01-05",
                "2021-04-02",
                "2021-04-26",
                "2021-05-17",
                "2021-07-19",
                "2021-08-09",
                "2021-11-29",
                "2022-01-24",
                "2022-04-18",
            ]
        ).normalize(),
        "end_date": pd.to_datetime(
            [
                "2020-03-23",
                "2020-05-28",
                "2020-06-18",
                "2020-07-09",
                "2020-10-01",
                "2020-11-01",
                "2021-01-04",
                "2021-04-01",
                "2021-04-25",
                "2021-05-16",
                "2021-07-18",
                "2021-08-08",
                "2021-11-28",
                "2022-01-23",
                "2022-04-17",
                "2023-05-05",
            ]
        ).normalize(),
        "intensity": [15, 100, 72, 52, 30, 55, 65, 95, 65, 55, 38, 20, 10, 42, 15, 3],
    }
)

POLICY_PERIODS["period_code"] = pd.Categorical(
    POLICY_PERIODS["period_code"],
    categories=_PERIOD_CODES,
    ordered=True,
)

PERIOD_LABELS: dict[str, str] = dict(
    zip(
        POLICY_PERIODS["period_code"].astype(str),
        POLICY_PERIODS["period_label"],
    )
)

PERIOD_INTENSITY: dict[str, int] = dict(
    zip(
        POLICY_PERIODS["period_code"].astype(str),
        POLICY_PERIODS["intensity"],
    )
)


def assign_period(dates: pd.Series) -> pd.Series:
    """Assign a policy period code to each date in the series.

    Parameters
    ----------
    dates:
        Series containing date-like values.

    Returns
    -------
    pandas.Series
        Series of policy period codes. Rows whose date falls outside all
        defined periods receive ``None``.
    """
    dates = pd.to_datetime(dates).dt.normalize()

    codes = np.full(len(dates), None, dtype=object)

    for _, row in POLICY_PERIODS.iterrows():
        mask = dates.between(row["start_date"], row["end_date"], inclusive="both")
        codes[mask.to_numpy()] = str(row["period_code"])

    return pd.Series(
        pd.Categorical(codes, categories=_PERIOD_CODES, ordered=True),
        index=dates.index,
        name="policy_period",
    )


def attach_period(df: pd.DataFrame, date_col: str) -> pd.DataFrame:
    """Attach policy period code, label, and intensity using the given date column.

    Parameters
    ----------
    df:
        Input DataFrame containing a date or datetime column.
    date_col:
        Name of the column to use for period assignment. Only the date
        component is used.

    Returns
    -------
    pandas.DataFrame
        A copy of ``df`` with three new columns:
        ``policy_period``, ``policy_period_label``, and ``policy_intensity``.
    """
    result = df.copy()
    result["policy_period"] = assign_period(result[date_col])

    period_lookup = POLICY_PERIODS[
        [
            "period_code",
            "period_label",
            "intensity",
        ]
    ].rename(
        columns={
            "period_code": "policy_period",
            "period_label": "policy_period_label",
            "intensity": "policy_intensity",
        }
    )

    result = result.merge(
        period_lookup,
        on="policy_period",
        how="left",
    )

    return result
