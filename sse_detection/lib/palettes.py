"""Categorical colour maps for the SSE-detection label space.

Colour decisions
----------------
The label space has two principal axes: ``sse_role`` (origin / continuation
character of the node) and ``sse_onward_dynamic`` (what happens downstream).
There are six roles and eight onward-dynamics, so the joint ``sse_category``
space has up to 48 cells — too many for a flat categorical palette.

The strategy used here is:

* ``ROLE_PALETTE`` assigns a colourblind-friendly base hue per role.
* ``DYNAMIC_PALETTE`` assigns stable colours for plots that encode only
onward dynamics.
* ``sse_category`` colours are derived from the role's hue, with lightness
modulated by the onward-dynamic so categories that share a role read as a
family. ``not_sse_like`` is always neutral grey.

Anything that needs the role colour directly (e.g. the role x onward-dynamic
heatmap, the metric-space scatter) reads from ``ROLE_PALETTE``. Anything that
needs the full category granularity (e.g. the meta-cluster subgraph) reads
from ``SSE_CATEGORY_PALETTE`` or calls
:func:`sse_category_palette_from` to build a palette for the subset present.
"""

from __future__ import annotations

from typing import Iterable
from typing import Callable

from matplotlib.colors import to_rgb

__all__ = [
    "sse_category_palette_from",
    "SSE_CATEGORY_PALETTE",
    "ROLE_ORDER",
    "ROLE_PALETTE",
    "DYNAMIC_ORDER",
    "DYNAMIC_PALETTE",
    "LIFECYCLE_PALETTE",
    "WAVE_GROUPS",
    "WAVE_GROUP_PALETTE",
]

NOT_SSE_COLOR = "#D0D3D8"

# ---------------------------------------------------------------------------
# Role / onward-dynamic ordering and base hues
# ---------------------------------------------------------------------------

ROLE_ORDER: list[str] = [
    "putative_birth",
    "relay_amplifier",
    "merged_relay",
    "terminal_sink",
    "isolated_burst",
    "unclear_origin",
]

ROLE_PALETTE: dict[str, str] = {
    "putative_birth": "#0072B2",
    "relay_amplifier": "#009E73",
    "merged_relay": "#D55E00",
    "terminal_sink": "#CC79A7",
    "isolated_burst": "#E69F00",
    "unclear_origin": "#6C6F73",
}

DYNAMIC_ORDER: list[str] = [
    "no_observed_onward_spread",
    "contained_burst",
    "single_successor_chain",
    "dominant_branch",
    "high_volume_onward_spread",
    "multi_branch_seeder",
    "multi_branch_expander",
    "diverse_population_broadcaster",
    "weak_or_ambiguous_onward_spread",
]

DYNAMIC_PALETTE: dict[str, str] = {
    "no_observed_onward_spread": "#6C6F73",
    "contained_burst": "#8A7B3F",
    "single_successor_chain": "#0072B2",
    "dominant_branch": "#332288",
    "high_volume_onward_spread": "#D55E00",
    "multi_branch_seeder": "#E69F00",
    "multi_branch_expander": "#CC79A7",
    "diverse_population_broadcaster": "#009E73",
    "weak_or_ambiguous_onward_spread": "#A6ACB3",
}

# Lightness multipliers per onward-dynamic. 0 = base hue, 1 = white.
# Order roughly: most "muted" outcomes are darker, most "active" outcomes
# are kept close to the base hue.
_DYNAMIC_LIGHTNESS: dict[str, float] = {
    "no_observed_onward_spread": 0.50,
    "contained_burst": 0.40,
    "single_successor_chain": 0.30,
    "dominant_branch": 0.22,
    "high_volume_onward_spread": 0.12,
    "multi_branch_seeder": 0.18,
    "multi_branch_expander": 0.05,
    "diverse_population_broadcaster": 0.00,
    "weak_or_ambiguous_onward_spread": 0.58,
}

LIFECYCLE_PALETTE: dict[str, str] = {
    "birth": "#0072B2",
    "continuation": "#009E73",
    "death": "#CC79A7",
    "unknown": "#6C6F73",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _lighten(color: str, amount: float) -> tuple[float, float, float]:
    """Blend ``color`` toward white by ``amount`` in [0, 1]."""
    amount = max(0.0, min(1.0, float(amount)))
    r, g, b = to_rgb(color)
    return r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount


def _split_category(category: str) -> tuple[str, str] | None:
    if not isinstance(category, str) or "__" not in category:
        return None
    role, dynamic = category.split("__", 1)
    return role, dynamic


def sse_category_palette_from(
    categories: Iterable[str],
    *,
    role_palette: dict[str, str] | None = None,
    dynamic_lightness: dict[str, float] | None = None,
) -> dict[str, tuple[float, float, float] | str]:
    """Build an ``sse_category`` colour map for the categories provided.

    Categories follow the ``role__dynamic`` schema. The role contributes the
    hue and the dynamic contributes the lightness modulation, so two
    categories that share a role read as a colour family.

    Parameters
    ----------
    categories
        Iterable of category labels (including any of ``"not_sse_like"``,
        ``None``, or unknown strings, which all map to the neutral grey).
    role_palette
        Optional override of the role hue map.
    dynamic_lightness
        Optional override of the per-dynamic lightness multipliers.
    """
    roles = role_palette or ROLE_PALETTE
    lightness = dynamic_lightness or _DYNAMIC_LIGHTNESS

    palette: dict[str, tuple[float, float, float] | str] = {}
    for cat in categories:
        if cat is None or cat == "not_sse_like" or cat != cat:  # NaN check
            palette[cat] = NOT_SSE_COLOR
            continue
        parts = _split_category(cat)
        if parts is None:
            palette[cat] = NOT_SSE_COLOR
            continue
        role, dynamic = parts
        base = roles.get(role, "#8C8C8C")
        palette[cat] = _lighten(base, lightness.get(dynamic, 0.35))
    return palette


def _build_full_category_palette() -> dict[str, tuple[float, float, float] | str]:
    cats = {"not_sse_like"}
    for role in ROLE_ORDER:
        for dyn in DYNAMIC_ORDER:
            cats.add(f"{role}__{dyn}")
    return sse_category_palette_from(sorted(cats))


SSE_CATEGORY_PALETTE: dict[str, tuple[float, float, float] | str] = (
    _build_full_category_palette()
)


WAVE_GROUPS: dict[str, Callable[[str], bool]] = {
    "B.1.177": lambda lineage: lineage.startswith("B.1.177"),
    "Alpha": lambda lineage: lineage == "B.1.1.7" or lineage.startswith("B.1.1.7."),
    "Delta": lambda lineage: lineage.startswith("AY.") or lineage == "B.1.617.2",
    "BA.1": lambda lineage: lineage.startswith("BA.1"),
    "BA.2": lambda lineage: lineage.startswith("BA.2"),
    "BA.4": lambda lineage: lineage.startswith("BA.4"),
    "BA.5": lambda lineage: lineage.startswith("BA.5") or lineage.startswith("BE."),
    "BQ.1": lambda lineage: lineage.startswith("BQ."),
    "XBB": lambda lineage: lineage.startswith("XBB"),
}

WAVE_GROUP_PALETTE: dict[str, str] = {
    "B.1.177": "#4e79a7",
    "Alpha": "#f28e2b",
    "Delta": "#e15759",
    "BA.1": "#76b7b2",
    "BA.2": "#59a14f",
    "BA.4": "#edc948",
    "BA.5": "#b07aa1",
    "BQ.1": "#ff9da7",
    "XBB": "#9c755f",
}
