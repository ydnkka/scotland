"""Shared policy and lineage lookup helpers.

The processed daily policy calendar is the source of truth for period and era
labels. This module normalises that calendar, provides window-to-policy mappers
based on the window midpoint date, and exposes lineage/clade lookup helpers for
the analysis packages.
"""

from __future__ import annotations

from collections.abc import Iterable as IterableABC
import re
from typing import Iterable

import numpy as np
import pandas as pd

from .data import CLADES, Paths, pango_lineages_for_clades

__all__ = [
    "load_policy_calendar",
    "policy_descriptors",
    "policy_order",
    "policy_period_labels",
    "policy_era_labels",
    "attach_policy_calendar",
    "window_idx_from_id",
    "window_id_from_idx",
    "window_policy_lookup",
    "window_id_to_policy_period_map",
    "window_id_to_policy_era_map",
    "window_id_to_policy_period",
    "window_id_to_policy_era",
    "lineage_clade_lookup",
    "lineage_to_clade_map",
    "clades_for_pango_lineages",
    "pango_lineages_for_clades",
]


_POLICY_COLUMN_ALIASES = {
    "policy_period": "period_code",
    "policy_period_label": "period_label",
    "policy_period_start_date": "period_start_date",
    "policy_period_end_date": "period_end_date",
    "policy_period_order": "period_order",
}
_POLICY_SOURCE_TO_NORMALISED = {
    source: normalized for normalized, source in _POLICY_COLUMN_ALIASES.items()
}


def _as_list(value: object | Iterable[object] | None) -> list[object] | None:
    """Return scalar or iterable input as a list, treating strings as scalars."""
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    if isinstance(value, IterableABC):
        return list(value)
    return [value]


def _normalise_str_values(
    values: object | Iterable[object] | None,
    *,
    name: str,
) -> list[str] | None:
    """Return string values with duplicates removed in input order."""
    raw_values = _as_list(values)
    if raw_values is None:
        return None

    out: list[str] = []
    for value in raw_values:
        if value is None:
            raise ValueError(f"{name} cannot contain None")
        out.append(str(value))

    return list(dict.fromkeys(out))


def _requested_policy_columns(
    columns: str | Iterable[str] | None,
) -> list[str] | None:
    """Translate normalised policy column names to the stored parquet schema."""
    raw_columns = _as_list(columns)
    if raw_columns is None:
        return None

    selected = [
        _POLICY_COLUMN_ALIASES.get(str(column), str(column)) for column in raw_columns
    ]
    if "date" not in selected:
        selected.insert(0, "date")
    return list(dict.fromkeys(selected))


def _normalise_window_ids(
    windows: object | Iterable[object] | None,
) -> list[str] | None:
    """Normalise window identifiers to the processed ``W###`` format."""
    raw_windows = _as_list(windows)
    if raw_windows is None:
        return None

    out: list[str] = []
    for window in raw_windows:
        if window is None:
            raise ValueError("windows cannot contain None")

        if isinstance(window, (int, np.integer)):
            idx = int(window)
            if idx < 1:
                raise ValueError("window indices must be positive")
            out.append(f"W{idx:03d}")
            continue

        if isinstance(window, (float, np.floating)) and float(window).is_integer():
            idx = int(window)
            if idx < 1:
                raise ValueError("window indices must be positive")
            out.append(f"W{idx:03d}")
            continue

        value = str(window).strip()
        if not value:
            raise ValueError("windows cannot contain empty strings")

        upper_value = value.upper()
        if upper_value.startswith("W") and upper_value[1:].isdigit():
            idx = int(upper_value[1:])
            if idx < 1:
                raise ValueError("window indices must be positive")
            out.append(f"W{idx:03d}")
        elif value.isdigit():
            idx = int(value)
            if idx < 1:
                raise ValueError("window indices must be positive")
            out.append(f"W{idx:03d}")
        else:
            out.append(upper_value)

    return list(dict.fromkeys(out))


def load_policy_calendar(
    columns: str | Iterable[str] | None = None,
    *,
    start_date: str | pd.Timestamp | None = None,
    end_date: str | pd.Timestamp | None = None,
    paths: Paths | None = None,
) -> pd.DataFrame:
    """Load the processed Scotland policy calendar with normalised column names."""
    paths = paths or Paths.from_config()
    selected = _requested_policy_columns(columns)

    policy = pd.read_parquet(paths.policy, columns=selected)
    policy["date"] = pd.to_datetime(policy["date"], errors="raise").dt.normalize()
    if policy["date"].duplicated().any():
        raise ValueError("Processed policy table contains duplicate dates.")

    if "period_start_date" in policy.columns:
        policy["period_start_date"] = pd.to_datetime(
            policy["period_start_date"], errors="raise"
        ).dt.normalize()
    if "period_end_date" in policy.columns:
        policy["period_end_date"] = pd.to_datetime(
            policy["period_end_date"], errors="raise"
        ).dt.normalize()
    if "period_order" in policy.columns:
        policy["period_order"] = pd.to_numeric(
            policy["period_order"], errors="raise"
        ).astype(int)

    if start_date is not None:
        policy = policy.loc[policy["date"].ge(pd.Timestamp(start_date).normalize())]
    if end_date is not None:
        policy = policy.loc[policy["date"].le(pd.Timestamp(end_date).normalize())]

    policy = policy.rename(columns=_POLICY_SOURCE_TO_NORMALISED)
    return policy.sort_values("date", ignore_index=True)


def policy_descriptors(*, paths: Paths | None = None) -> pd.DataFrame:
    """Return the unique policy periods in chronological order."""
    policy = load_policy_calendar(
        [
            "policy_period",
            "policy_period_label",
            "policy_period_start_date",
            "policy_period_end_date",
            "policy_period_order",
            "policy_era",
        ],
        paths=paths,
    )
    return (
        policy[
            [
                "policy_period",
                "policy_period_label",
                "policy_period_start_date",
                "policy_period_end_date",
                "policy_period_order",
                "policy_era",
            ]
        ]
        .drop_duplicates(subset=["policy_period"], keep="first")
        .sort_values("policy_period_order", kind="mergesort")
        .reset_index(drop=True)
    )


def policy_order(column: str, *, paths: Paths | None = None) -> list[str]:
    """Return policy categories in their canonical chronological order."""
    descriptors = policy_descriptors(paths=paths)

    if column in {"policy_period", "period_code"}:
        return descriptors["policy_period"].astype(str).tolist()
    if column in {"policy_period_label", "period_label"}:
        return descriptors["policy_period_label"].astype(str).tolist()
    if column == "policy_era":
        return list(dict.fromkeys(descriptors["policy_era"].astype(str)))

    raise KeyError(
        f"Unsupported policy column {column!r}; expected policy_period, "
        "policy_period_label, or policy_era."
    )


def policy_period_labels(*, paths: Paths | None = None) -> dict[str, str]:
    """Return a mapping from policy period code to human-readable label."""
    descriptors = policy_descriptors(paths=paths)
    return dict(
        zip(
            descriptors["policy_period"].astype(str),
            descriptors["policy_period_label"].astype(str),
        )
    )


def policy_era_labels(*, paths: Paths | None = None) -> dict[str, str]:
    """Return a display-label mapping for policy eras."""
    return {
        era: era.upper().replace("_", " ")
        for era in policy_order("policy_era", paths=paths)
    }


def attach_policy_calendar(
    df: pd.DataFrame,
    date_col: str,
    *,
    paths: Paths | None = None,
    drop_existing: bool = True,
) -> pd.DataFrame:
    """Attach the processed policy calendar to a frame using a date column."""
    if date_col not in df.columns:
        raise KeyError(f"{date_col!r} is required for the policy-calendar join")

    calendar = load_policy_calendar(paths=paths).rename(columns={"date": "_policy_date"})
    policy_cols = [col for col in calendar.columns if col != "_policy_date"]

    out = df.copy()
    if drop_existing:
        out = out.drop(columns=policy_cols, errors="ignore")

    out["_policy_date"] = pd.to_datetime(out[date_col], errors="coerce").dt.normalize()
    merged = out.merge(
        calendar,
        on="_policy_date",
        how="left",
        validate="many_to_one",
    )
    return merged.drop(columns="_policy_date")


def window_idx_from_id(value: object | pd.Series) -> object:
    """Extract a numeric window index from a window ID or series of IDs."""
    if isinstance(value, pd.Series):
        extracted = value.astype("string").str.extract(r"(\d+)")[0]
        return pd.to_numeric(extracted, errors="coerce")

    if value is None or pd.isna(value):
        return None

    match = re.search(r"(\d+)", str(value))
    if not match:
        return None
    return int(match.group(1))


def window_id_from_idx(value: object | pd.Series) -> object:
    """Convert a numeric window index or series of indices to ``W###`` IDs."""
    if isinstance(value, pd.Series):
        return value.map(window_id_from_idx)

    if value is None or pd.isna(value):
        return None

    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        match = re.search(r"(\d+)", stripped)
        if not match:
            raise ValueError(f"Cannot convert {value!r} to a window ID")
        idx = int(match.group(1))
    elif isinstance(value, (int, np.integer)):
        idx = int(value)
    elif isinstance(value, (float, np.floating)):
        if not float(value).is_integer():
            raise ValueError(f"Window indices must be whole numbers, got {value!r}")
        idx = int(value)
    else:
        match = re.search(r"(\d+)", str(value))
        if not match:
            raise ValueError(f"Cannot convert {value!r} to a window ID")
        idx = int(match.group(1))

    if idx < 1:
        raise ValueError("window indices must be positive")
    return f"W{idx:03d}"


def window_policy_lookup(
    windows: str | int | Iterable[str | int] | None = None,
    *,
    paths: Paths | None = None,
) -> pd.DataFrame:
    """Return one row per window with policy-period and policy-era annotations."""
    paths = paths or Paths.from_config()
    window_ids = _normalise_window_ids(windows)

    window = pd.read_parquet(
        paths.analysis_dataset,
        columns=["window_id", "window_idx", "wn_mid_date"],
    )
    window["window_id"] = window["window_id"].astype(str)
    window["window_idx"] = pd.to_numeric(window["window_idx"], errors="coerce")
    window["wn_mid_date"] = pd.to_datetime(window["wn_mid_date"], errors="coerce")
    window = window.drop_duplicates(["window_id", "window_idx"]).sort_values(
        "window_idx"
    )
    if window_ids is not None:
        window = window.loc[window["window_id"].isin(window_ids)]
    if window.empty:
        detail = f"windows={window_ids}" if window_ids is not None else "the dataset"
        raise ValueError(f"No window-policy rows found for {detail}.")

    return attach_policy_calendar(window, "wn_mid_date", paths=paths)


def window_id_to_policy_period_map(
    windows: str | int | Iterable[str | int] | None = None,
    *,
    paths: Paths | None = None,
) -> dict[str, str]:
    """Return a mapping from window ID to policy-period code."""
    lookup = window_policy_lookup(windows=windows, paths=paths)
    return dict(
        zip(
            lookup["window_id"].astype(str),
            lookup["policy_period"].astype(str),
        )
    )


def window_id_to_policy_era_map(
    windows: str | int | Iterable[str | int] | None = None,
    *,
    paths: Paths | None = None,
) -> dict[str, str]:
    """Return a mapping from window ID to policy-era label."""
    lookup = window_policy_lookup(windows=windows, paths=paths)
    return dict(
        zip(
            lookup["window_id"].astype(str),
            lookup["policy_era"].astype(str),
        )
    )


def window_id_to_policy_period(
    window_id: str | int,
    *,
    paths: Paths | None = None,
) -> str | None:
    """Return the policy period for a single window ID or window index."""
    normalised = window_id_from_idx(window_id)
    if normalised is None:
        return None
    return window_id_to_policy_period_map(paths=paths).get(str(normalised))


def window_id_to_policy_era(
    window_id: str | int,
    *,
    paths: Paths | None = None,
) -> str | None:
    """Return the policy era for a single window ID or window index."""
    normalised = window_id_from_idx(window_id)
    if normalised is None:
        return None
    return window_id_to_policy_era_map(paths=paths).get(str(normalised))


def lineage_clade_lookup(
    windows: str | int | Iterable[str | int] | None = None,
    *,
    paths: Paths | None = None,
    display_labels: bool = False,
) -> pd.DataFrame:
    """Return unique Pango-lineage to Nextclade-clade pairs.

    When *windows* is provided, the lookup is restricted to those rolling
    windows before duplicates are removed. If *display_labels* is True, an
    additional ``clade_label`` column is included with the display labels from
    ``utils.data.CLADES``.
    """
    paths = paths or Paths.from_config()
    window_ids = _normalise_window_ids(windows)

    columns = ["pango_lineage", "clade"]
    if window_ids is not None:
        columns = ["window_id", *columns]

    df = pd.read_parquet(paths.analysis_dataset, columns=columns)
    df = df.dropna(subset=["pango_lineage", "clade"]).copy()
    df["pango_lineage"] = df["pango_lineage"].astype(str)
    df["clade"] = df["clade"].astype(str)

    if window_ids is not None:
        df["window_id"] = df["window_id"].astype(str)
        df = df.loc[df["window_id"].isin(window_ids)]

    df = (
        df.drop_duplicates(subset=["pango_lineage", "clade"])
        .sort_values(["pango_lineage", "clade"], kind="mergesort")
        .reset_index(drop=True)
    )

    if df.empty:
        detail = f"windows={window_ids}" if window_ids is not None else "the dataset"
        raise ValueError(f"No lineage/clade pairs found for {detail}.")

    if display_labels:
        df["clade_label"] = df["clade"].map(CLADES).fillna(df["clade"])

    return df.drop(columns=["window_id"], errors="ignore")


def lineage_to_clade_map(
    windows: str | int | Iterable[str | int] | None = None,
    *,
    paths: Paths | None = None,
    display_labels: bool = False,
) -> dict[str, tuple[str, ...]]:
    """Return a lineage-to-clade mapping with deterministic tuple values."""
    lookup = lineage_clade_lookup(
        windows=windows,
        paths=paths,
        display_labels=display_labels,
    )
    clade_col = "clade_label" if display_labels else "clade"
    grouped = (
        lookup.groupby("pango_lineage", sort=False)[clade_col]
        .agg(lambda values: tuple(dict.fromkeys(values.astype(str))))
        .to_dict()
    )
    return grouped


def clades_for_pango_lineages(
    lineages: str | Iterable[str],
    *,
    windows: str | int | Iterable[str | int] | None = None,
    paths: Paths | None = None,
    display_labels: bool = False,
) -> list[str]:
    """Resolve one or more Pango lineages to their associated clade labels."""
    lineage_values = _normalise_str_values(lineages, name="lineages")
    if lineage_values is None:
        raise ValueError("lineages must contain at least one value")

    lineage_map = lineage_to_clade_map(
        windows=windows,
        paths=paths,
        display_labels=display_labels,
    )

    clades: list[str] = []
    for lineage in lineage_values:
        clades.extend(lineage_map.get(lineage, ()))

    clades = list(dict.fromkeys(clades))
    if not clades:
        detail = f"lineages={lineage_values}"
        if windows is not None:
            detail += f", windows={_normalise_window_ids(windows)}"
        raise ValueError(f"No clades found for {detail}.")
    return clades
