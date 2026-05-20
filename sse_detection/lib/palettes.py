"""Categorical colour maps for the SSE-detection label space.

Colour decisions
----------------
The label space has two principal axes: ``sse_role`` (origin / continuation
character of the node) and ``sse_onward_dynamic`` (what happens downstream).
There are six roles and eight onward-dynamics, so the joint ``sse_category``
space has up to 48 cells — too many for a flat categorical palette.

The strategy used here is:

* ``ROLE_PALETTE`` assigns a base hue per role.
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

from matplotlib.colors import to_rgb

NOT_SSE_COLOR = "#bfbfbf"

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
    "putative_birth": "#3A6EA5",
    "relay_amplifier": "#1F7A75",
    "merged_relay": "#C75C2C",
    "terminal_sink": "#7C8A43",
    "isolated_burst": "#6B5F7A",
    "unclear_origin": "#8C8C8C",
}

DYNAMIC_ORDER: list[str] = [
    "no_observed_onward_spread",
    "contained_burst",
    "single_dominant_chain",
    "high_volume_onward_spread",
    "multi_branch_seeder",
    "multi_branch_expander",
    "diffuse_spatial_broadcaster",
    "weak_or_ambiguous_onward_spread",
]

# Lightness multipliers per onward-dynamic. 0 = base hue, 1 = white.
# Order roughly: most "muted" outcomes are darker, most "active" outcomes
# are kept close to the base hue.
_DYNAMIC_LIGHTNESS: dict[str, float] = {
    "no_observed_onward_spread": 0.55,
    "contained_burst": 0.45,
    "single_dominant_chain": 0.30,
    "high_volume_onward_spread": 0.15,
    "multi_branch_seeder": 0.20,
    "multi_branch_expander": 0.05,
    "diffuse_spatial_broadcaster": 0.00,
    "weak_or_ambiguous_onward_spread": 0.65,
}

LIFECYCLE_PALETTE: dict[str, str] = {
    "birth": "#3A6EA5",
    "continuation": "#1F7A75",
    "death": "#7C8A43",
    "unknown": "#8C8C8C",
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
