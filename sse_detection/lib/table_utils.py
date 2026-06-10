"""Shared table and label helpers for SSE report figures."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np
import pandas as pd


def read_table(
    table: pd.DataFrame | str | Path | Any,
    *,
    skipinitialspace: bool = True,
) -> pd.DataFrame:
    """Return a dataframe from an in-memory table or a CSV/parquet path."""
    if isinstance(table, pd.DataFrame):
        return table.copy()
    path = Path(table)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, skipinitialspace=skipinitialspace)


def clean_strings(df: pd.DataFrame) -> pd.DataFrame:
    """Strip whitespace from column names and string-like values."""
    out = df.copy()
    out.columns = [str(col).strip() for col in out.columns]

    for col in out.select_dtypes(include=["object", "string"]).columns:
        out[col] = out[col].str.strip()

    return out


def pretty_text(value: Any, label_map: Mapping[str, str] | None = None) -> str:
    """Return a display label for a model term or category value."""
    if pd.isna(value):
        return ""
    text = str(value)
    if label_map and text in label_map:
        return label_map[text]
    return text.replace("_", " ").strip().capitalize()


def term_level(term: Any) -> str:
    """Extract the level from a patsy categorical term string."""
    if pd.isna(term):
        return ""
    match = re.search(r"\[T\.(.*)\]$", str(term))
    return match.group(1) if match else str(term)


def forest_xlim(panel: pd.DataFrame) -> tuple[float, float]:
    """Return a padded positive x-axis range for odds-ratio forest panels."""
    values = pd.to_numeric(
        pd.concat([panel["or_low"], panel["or_high"], pd.Series([1.0])]),
        errors="coerce",
    )
    values = values[np.isfinite(values) & values.gt(0)]
    if values.empty:
        return (0.75, 1.35)
    lo = float(values.min())
    hi = float(values.max())
    log_lo = np.log(lo)
    log_hi = np.log(hi)
    pad = max((log_hi - log_lo) * 0.12, 0.08)
    return float(np.exp(log_lo - pad)), float(np.exp(log_hi + pad))
