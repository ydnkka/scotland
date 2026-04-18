"""Publication-grade plotting defaults for the Scotland clustering papers.

The style is deliberately conservative: a sans-serif body font at ~8-9 pt,
thin axes, no top/right spines, and a single-column width tuned to *Lancet*
and *Virus Evolution* templates. Palettes are colour-blind friendly and
use perceptually uniform ramps where a quantitative variable is mapped.

Usage
-----
>>> from manuscripts.common import style
>>> style.set_theme()
>>> fig, ax = style.new_figure("single")
"""

from __future__ import annotations

from typing import Literal, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb

# ---------------------------------------------------------------------------
# Figure size presets (inches). Aligned with common journal guidance.
# ---------------------------------------------------------------------------

FIG_WIDTHS_IN = {
    "single": 3.5,       # single column
    "onehalf": 5.2,      # 1.5-column
    "double": 7.2,       # full-page double column
    "slide": 10.0,       # for talks
}

DEFAULT_HEIGHT_IN = 2.6

# ---------------------------------------------------------------------------
# Categorical palettes
# ---------------------------------------------------------------------------

# SIMD quintile: 1 = most deprived (deep red) -> 5 = least (desaturated blue).
# Hand-picked from ColorBrewer RdBu to remain legible in greyscale.
SIMD_QUINTILE_PALETTE: dict[int, str] = {
    1: "#b2182b",
    2: "#ef8a62",
    3: "#cccccc",
    4: "#67a9cf",
    5: "#2166ac",
}

# WHO variants of concern. Order follows emergence in Scotland.
WHO_VOC_PALETTE: dict[str, str] = {
    "Alpha":        "#4e79a7",
    "Beta":         "#f28e2b",
    "Gamma":        "#e15759",
    "Delta":        "#59a14f",
    "Eta":          "#edc948",
    "Kappa":        "#b07aa1",
    "Mu":           "#ff9da7",
    "Theta":        "#9c755f",
    "Omicron":      "#af2d2d",
    "recombinant":  "#bab0ac",
    "None":         "#cccccc",
}

# Leiden resolution: monotonic palette from low -> high.
RESOLUTION_CMAP = LinearSegmentedColormap.from_list(
    "resolution", ["#4b0082", "#9e2f8f", "#d94b76", "#ee7b4c", "#f0c300"]
)

# SIMD-domain palette (seven domains plus overall). Qualitative but stable.
SIMD_DOMAIN_PALETTE: dict[str, str] = {
    "overall":    "#2b2b2b",
    "income":     "#4e79a7",
    "employment": "#f28e2b",
    "education":  "#59a14f",
    "health":     "#e15759",
    "access":     "#b07aa1",
    "crime":      "#edc948",
    "housing":    "#9c755f",
}

# ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------


def set_theme() -> None:
    """Apply the paper-wide matplotlib rcParams."""
    mpl.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 400,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "pdf.fonttype": 42,           # embed TrueType — editable in Illustrator
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.titlesize": 9.0,
            "axes.titleweight": "bold",
            "axes.labelsize": 8.5,
            "axes.labelweight": "regular",
            "axes.linewidth": 0.6,
            "axes.edgecolor": "#333333",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "legend.fontsize": 7.5,
            "legend.frameon": False,
            "lines.linewidth": 1.1,
            "lines.markersize": 3.0,
            "image.cmap": "viridis",
        }
    )


def new_figure(
    width: Literal["single", "onehalf", "double", "slide"] = "single",
    height_in: float | None = None,
    nrows: int = 1,
    ncols: int = 1,
    **subplots_kwargs,
) -> tuple[plt.Figure, plt.Axes]:
    """Create a figure with one of the paper-width presets."""
    set_theme()
    w = FIG_WIDTHS_IN[width]
    h = height_in if height_in is not None else DEFAULT_HEIGHT_IN * nrows
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(w, h), **subplots_kwargs)
    return fig, ax


def save_figure(fig: plt.Figure, out_path, *, formats: Sequence[str] = ("pdf", "png")) -> list:
    """Save a figure in every requested format. Returns list of paths written."""
    from pathlib import Path

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in formats:
        p = out_path.with_suffix(f".{fmt}")
        fig.savefig(p, format=fmt)
        written.append(p)
    return written


def lighten(color: str, amount: float = 0.4) -> tuple[float, float, float]:
    """Blend `color` toward white by `amount` in [0, 1]."""
    r, g, b = to_rgb(color)
    return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)
