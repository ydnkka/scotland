"""Consistent attribute names for publication figures and tables."""

import re


URBAN_RURAL_CLASS_LABEL = "Urban/rural class"


def normalise_attribute_label(label: str) -> str:
    """Standardise display text, including labels in older saved summaries."""
    return re.sub(
        r"\b(?:urban[/–-]+rural(?:\s+class)?|settlement type)\b",
        URBAN_RURAL_CLASS_LABEL,
        label,
        flags=re.IGNORECASE,
    )
