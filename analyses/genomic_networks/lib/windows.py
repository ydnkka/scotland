"""Rolling-window identifier helpers for genomic-network analyses."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any


def normalise_window(value: Any) -> str:
    """Return a canonical ``W000``-style window identifier when possible."""
    text = str(value).strip()
    upper = text.upper()
    if upper.startswith("W") and upper[1:].isdigit():
        return f"W{int(upper[1:]):03d}"
    if upper.isdigit():
        return f"W{int(upper):03d}"
    return text


def normalise_windows(values: Iterable[Any] | None) -> list[str] | None:
    """Return canonical window identifiers, preserving ``None`` as no filter."""
    if values is None:
        return None
    return [normalise_window(value) for value in values]
