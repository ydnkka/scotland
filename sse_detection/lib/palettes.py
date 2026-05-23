"""Categorical colour maps for the SSE-detection label space.

``sse_category`` is the compact epidemiological transmission-phenotype label.
The graph-diagnostic ``sse_graph_category`` preserves the more granular
``sse_role`` x ``sse_onward_dynamic`` composite.
"""

from __future__ import annotations

from typing import Iterable, Any

from matplotlib.colors import to_rgb

__all__ = [
    "sse_category_palette_from",
    "sse_graph_category_palette_from",
    "SSE_CATEGORY_ORDER",
    "SSE_CATEGORY_PALETTE",
    "SSE_GRAPH_CATEGORY_PALETTE",
    "ROLE_ORDER",
    "ROLE_PALETTE",
    "DYNAMIC_ORDER",
    "DYNAMIC_PALETTE",
    "LIFECYCLE_PALETTE",
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

SSE_CATEGORY_ORDER: list[str] = [
    "not_sse_like",
    "mixed_population_dissemination",
    "putative_introduction_burst",
    "secondary_relay_amplification",
    "diffuse_branching_transmission",
    "focused_branching_transmission",
    "sustained_single_chain",
    "contained_local_burst",
    "high_volume_onward_transmission",
    "ambiguous_amplification_signal",
]

SSE_CATEGORY_PALETTE: dict[str, str] = {
    "not_sse_like": NOT_SSE_COLOR,
    "mixed_population_dissemination": "#009E73",
    "putative_introduction_burst": "#56B4E9",
    "secondary_relay_amplification": "#D55E00",
    "diffuse_branching_transmission": "#E69F00",
    "focused_branching_transmission": "#332288",
    "sustained_single_chain": "#0072B2",
    "contained_local_burst": "#8A7B3F",
    "high_volume_onward_transmission": "#CC79A7",
    "ambiguous_amplification_signal": "#A6ACB3",
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
    category_palette: dict[str, str] | None = None,
) -> Any:
    """Build an ``sse_category`` colour map for the categories provided.

    Unknown labels map to neutral grey. Legacy ``role__dynamic`` labels are
    accepted as a compatibility fallback and coloured as graph categories.

    Parameters
    ----------
    categories
        Iterable of category labels (including any of ``"not_sse_like"``,
        ``None``, or unknown strings, which all map to the neutral grey).
    category_palette
        Optional override of the compact category colour map.
    """
    category_colours = category_palette or SSE_CATEGORY_PALETTE

    palette: dict[str, tuple[float, float, float] | str] = {}
    graph_palette = sse_graph_category_palette_from(categories)
    for cat in categories:
        if cat is None or cat == "not_sse_like" or cat != cat:  # NaN check
            palette[cat] = NOT_SSE_COLOR
        elif cat in category_colours:
            palette[cat] = category_colours[cat]
        elif isinstance(cat, str) and "__" in cat:
            palette[cat] = graph_palette.get(cat, NOT_SSE_COLOR)
        else:
            palette[cat] = NOT_SSE_COLOR
    return palette


def sse_graph_category_palette_from(
    categories: Iterable[str],
    *,
    role_palette: dict[str, str] | None = None,
    dynamic_lightness: dict[str, float] | None = None,
) -> Any:
    """Build a palette for ``sse_graph_category`` composite labels."""
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


def _build_full_graph_category_palette() -> dict[str, tuple[float, float, float] | str]:
    cats = {"not_sse_like"}
    for role in ROLE_ORDER:
        for dyn in DYNAMIC_ORDER:
            cats.add(f"{role}__{dyn}")
    return sse_graph_category_palette_from(sorted(cats))


SSE_GRAPH_CATEGORY_PALETTE: dict[str, tuple[float, float, float] | str] = (
    _build_full_graph_category_palette()
)
