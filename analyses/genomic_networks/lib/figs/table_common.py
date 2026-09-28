"""Shared data loading and formatting for publication tables."""

from __future__ import annotations

from typing import Any

import numpy as np


def fmt_int(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    try:
        value_float = float(value)
    except (TypeError, ValueError):
        return str(value)
    if value_float.is_integer():
        return f"{int(value_float):,}"
    return f"{value_float:,.0f}"


def fmt_percent(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{100 * float(value):.1f}%"


def fmt_float(value: Any, digits: int = 2) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{float(value):.{digits}f}"


def fmt_iqr(
    median: Any,
    q25: Any,
    q75: Any,
    digits: int = 0,
) -> str:
    if median is None or (isinstance(median, float) and np.isnan(median)):
        return ""
    return (
        f"{fmt_float(median, digits)} "
        f"({fmt_float(q25, digits)}--{fmt_float(q75, digits)})"
    )
