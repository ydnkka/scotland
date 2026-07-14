"""Policy-period definitions and helpers for Scottish COVID-19 restrictions.

The exported period table assigns ordered restriction phases to calendar dates
and exposes compact labels plus OxCGRT policy indices for plotting and modelling.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


__all__ = [
    "PERIOD_ORDER",
    "POLICY_PERIODS",
    "PERIOD_LABELS",
    "PERIOD_STRINGENCY",
    "PERIOD_CONTAINMENT",
    "OXCGRT_STRINGENCY_PATH",
    "OXCGRT_CONTAINMENT_PATH",
    "load_oxcgrt_stringency",
    "load_oxcgrt_containment",
    "derive_period_stringency",
    "derive_period_containment",
    "assign_period",
    "attach_period",
]

OXCGRT_STRINGENCY_PATH = (
    Path(__file__).resolve().parents[1]
    / "data/raw/oxcgrt/OxCGRT_timeseries_StringencyIndex_v1.csv"
)
OXCGRT_CONTAINMENT_PATH = (
    Path(__file__).resolve().parents[1]
    / "data/raw/oxcgrt/OxCGRT_timeseries_ContainmentHealthIndex_v1.csv"
)


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


def _load_oxcgrt_index(
    path: Path | str,
    *,
    value_col: str,
    index_label: str,
    region_name: str = "Scotland",
) -> pd.DataFrame:
    """Load one daily OxCGRT index for a region from a wide timeseries table."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"OxCGRT {index_label} table not found: {path}")

    table = pd.read_csv(path)
    if "RegionName" not in table.columns:
        raise KeyError(f"OxCGRT {index_label} table needs 'RegionName'.")
    selected = table.loc[table["RegionName"].eq(region_name)]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one OxCGRT row for {region_name!r}; found {len(selected)}."
        )

    parsed_dates = pd.to_datetime(table.columns, format="%d%b%Y", errors="coerce")
    date_mask = parsed_dates.notna()
    date_columns = table.columns[date_mask]
    values = pd.to_numeric(selected.iloc[0][date_columns], errors="coerce")
    return pd.DataFrame(
        {
            "date": parsed_dates[date_mask],
            value_col: values.to_numpy(dtype=float),
        }
    ).sort_values("date", ignore_index=True)


def load_oxcgrt_stringency(
    path: Path | str = OXCGRT_STRINGENCY_PATH,
    *,
    region_name: str = "Scotland",
) -> pd.DataFrame:
    """Load the daily OxCGRT Stringency Index for one region."""
    return _load_oxcgrt_index(
        path,
        value_col="stringency_index",
        index_label="Stringency Index",
        region_name=region_name,
    )


def load_oxcgrt_containment(
    path: Path | str = OXCGRT_CONTAINMENT_PATH,
    *,
    region_name: str = "Scotland",
) -> pd.DataFrame:
    """Load the daily OxCGRT Containment and Health Index for one region."""
    return _load_oxcgrt_index(
        path,
        value_col="containment_index",
        index_label="Containment and Health Index",
        region_name=region_name,
    )


def _derive_period_index(
    policy_periods: pd.DataFrame,
    daily_index: pd.DataFrame,
    *,
    value_col: str,
) -> pd.Series:
    """Return an inclusive daily-index mean for each policy-period interval."""
    required_period = {"period_code", "start_date", "end_date"}
    missing_period = sorted(required_period - set(policy_periods.columns))
    if missing_period:
        raise KeyError(f"Missing policy-period columns: {missing_period}")
    required_daily = {"date", value_col}
    missing_daily = sorted(required_daily - set(daily_index.columns))
    if missing_daily:
        raise KeyError(f"Missing daily index columns: {missing_daily}")

    dates = pd.to_datetime(daily_index["date"], errors="coerce").dt.normalize()
    values = pd.to_numeric(daily_index[value_col], errors="coerce")
    means = {
        str(row.period_code): float(
            values.loc[dates.between(row.start_date, row.end_date, inclusive="both")].mean()
        )
        for row in policy_periods.itertuples(index=False)
    }
    return policy_periods["period_code"].astype(str).map(means).astype(float)


def derive_period_stringency(
    policy_periods: pd.DataFrame,
    daily_stringency: pd.DataFrame,
) -> pd.Series:
    """Return mean daily stringency for each inclusive policy-period interval."""
    return _derive_period_index(
        policy_periods,
        daily_stringency,
        value_col="stringency_index",
    )


def derive_period_containment(
    policy_periods: pd.DataFrame,
    daily_containment: pd.DataFrame,
) -> pd.Series:
    """Return mean daily containment for each inclusive policy-period interval."""
    return _derive_period_index(
        policy_periods,
        daily_containment,
        value_col="containment_index",
    )


POLICY_PERIODS["policy_stringency"] = derive_period_stringency(
    POLICY_PERIODS,
    load_oxcgrt_stringency(),
)
POLICY_PERIODS["policy_containment"] = derive_period_containment(
    POLICY_PERIODS,
    load_oxcgrt_containment(),
)

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

PERIOD_STRINGENCY: dict[str, float] = dict(
    zip(
        POLICY_PERIODS["period_code"].astype(str),
        POLICY_PERIODS["policy_stringency"],
    )
)
PERIOD_CONTAINMENT: dict[str, float] = dict(
    zip(
        POLICY_PERIODS["period_code"].astype(str),
        POLICY_PERIODS["policy_containment"],
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
    """Attach policy period code, label, and OxCGRT indices using a date column.

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
        A copy of ``df`` with ``policy_period``, ``policy_period_label``,
        ``policy_stringency``, and ``policy_containment``.
    """
    result = df.copy()
    result["policy_period"] = assign_period(result[date_col])

    period_lookup = POLICY_PERIODS[
        [
            "period_code",
            "period_label",
            "policy_stringency",
            "policy_containment",
        ]
    ].rename(
        columns={
            "period_code": "policy_period",
            "period_label": "policy_period_label",
            "policy_stringency": "policy_stringency",
            "policy_containment": "policy_containment",
        }
    )

    result = result.merge(
        period_lookup,
        on="policy_period",
        how="left",
    )

    return result
