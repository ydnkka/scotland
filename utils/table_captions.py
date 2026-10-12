"""Load publication captions and column descriptions by LaTeX label."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml

TABLE_CAPTIONS_PATH = (
    Path(__file__).resolve().parents[1] / "results" / "table_captions.yaml"
)
_PLACEHOLDER = re.compile(r"\{\{\s*([A-Za-z_]\w*)\s*\}\}")


def load_table_metadata(label: str) -> dict[str, Any]:
    """Read one entry, including editorial column descriptions and notes.

    Resolve the registry relative to this module, independent of the working
    directory. Read it afresh so edits also take effect in interactive sessions.
    """
    with TABLE_CAPTIONS_PATH.open(encoding="utf-8") as handle:
        registry = yaml.safe_load(handle)
    if not isinstance(registry, dict):
        raise TypeError(f"Expected a label-keyed mapping in {TABLE_CAPTIONS_PATH}")
    if label not in registry:
        raise KeyError(f"No caption for {label!r} in {TABLE_CAPTIONS_PATH}")
    entry = registry[label]
    if not isinstance(entry, dict):
        raise TypeError(f"Expected a metadata mapping for {label!r}")
    for field in ("caption", "short_caption"):
        if not isinstance(entry.get(field), str) or not entry[field].strip():
            raise ValueError(f"{label!r}: {field} must be a non-empty string")
    if not isinstance(entry.get("caption_is_latex", False), bool):
        raise TypeError(f"{label!r}: caption_is_latex must be a boolean")
    columns = entry.get("columns", {})
    if not isinstance(columns, dict) or any(
        not isinstance(key, str) or not isinstance(value, str)
        for key, value in columns.items()
    ):
        raise ValueError(f"{label!r}: columns must map headings to descriptions")
    return entry


def table_caption(label: str, **values: Any) -> dict[str, Any]:
    """Return caption keyword arguments for the shared LaTeX table writers.

    Only ``{{name}}`` placeholders are substituted, leaving LaTeX braces and
    math untouched. Supply preformatted values, e.g. ``n_pairs="1,234"``.
    Column descriptions and notes are metadata, not rendered caption text.
    """
    entry = load_table_metadata(label)

    def substitute(match: re.Match[str]) -> str:
        name = match.group(1)
        if name not in values:
            raise ValueError(f"{label!r}: missing caption value {name!r}")
        return str(values[name])

    return {
        "label": label,
        "caption": _PLACEHOLDER.sub(substitute, entry["caption"].strip()),
        "short_caption": _PLACEHOLDER.sub(substitute, entry["short_caption"].strip()),
        "caption_is_latex": entry.get("caption_is_latex", False),
    }
