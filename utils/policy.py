"""Policy-period definitions and helpers for Scottish COVID-19 restrictions.

The exported period table assigns ordered restriction phases to calendar dates
and exposes compact labels/stringency scores for plotting and modelling.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


__all__ = [
    "PERIOD_ORDER",
    "POLICY_PERIODS",
    "PERIOD_LABELS",
    "PERIOD_STRINGENCY",
    "assign_period",
    "attach_period",
]

#  Mean stringency index values for each policy period, derived from the
#  OxCGRT Stringency Index (<https://github.com/OxCGRT/covid-policy-dataset/>),
#  File: `timeseries_indices/OxCGRT_timeseries_StringencyIndex_v1.csv`

stringency_index = {
    "E0": np.float64(23.911739130434785),
    "L1": np.float64(79.63000000000008),
    "P1": np.float64(75.13238095238096),
    "P2": np.float64(76.45619047619047),
    "P3": np.float64(67.64571428571442),
    "T1": np.float64(64.80999999999996),
    "F5": np.float64(70.88968750000001),
    "L2": np.float64(85.53586206896566),
    "SL": np.float64(69.90375000000003),
    "L3": np.float64(58.330000000000005),
    "L21": np.float64(56.13063492063501),
    "L0": np.float64(52.77999999999999),
    "NN": np.float64(31.637857142857182),
    "OM": np.float64(34.83857142857142),
    "FE": np.float64(19.786666666666626),
    "PR": np.float64(8.184418604651126),
}


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
    }
)

PERIOD_ORDER: list[str] = POLICY_PERIODS["period_code"].tolist()
POLICY_PERIODS["policy_stringency"] = POLICY_PERIODS["period_code"].map(stringency_index)

POLICY_PERIODS["period_code"] = pd.Categorical(
    POLICY_PERIODS["period_code"],
    categories=PERIOD_ORDER,
    ordered=True,
)

PERIOD_LABELS: dict[str, str] = dict(
    zip(
        POLICY_PERIODS["period_code"].astype(str),
        POLICY_PERIODS["period_label"],
    )
)

PERIOD_STRINGENCY: dict[str, int] = dict(
    zip(
        POLICY_PERIODS["period_code"].astype(str),
        POLICY_PERIODS["policy_stringency"],
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
        pd.Categorical(codes, categories=PERIOD_ORDER, ordered=True),
        index=dates.index,
        name="policy_period",
    )


def attach_period(df: pd.DataFrame, date_col: str) -> pd.DataFrame:
    """Attach policy period code, label, and stringency using the given date column.

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
        ``policy_period``, ``policy_period_label``, and ``policy_stringency``.
    """
    result = df.copy()
    result["policy_period"] = assign_period(result[date_col])

    period_lookup = POLICY_PERIODS[
        [
            "period_code",
            "period_label",
            "policy_stringency",
        ]
    ].rename(
        columns={
            "period_code": "policy_period",
            "period_label": "policy_period_label",
            "policy_stringency": "policy_stringency",
        }
    )

    result = result.merge(
        period_lookup,
        on="policy_period",
        how="left",
    )

    return result
