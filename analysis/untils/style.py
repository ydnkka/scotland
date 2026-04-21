"""Publication-grade plotting defaults for the Scotland clustering papers.

The style is deliberately conservative: a sans-serif body font at ~8-9 pt,
thin axes, no top/right spines, and a single-column width tuned to *Lancet*
and *Virus Evolution* templates. Palettes are colour-blind friendly and
use perceptually uniform ramps where a quantitative variable is mapped.

Usage
-----
>>> from analysis.untils import style
>>> style.set_theme()
>>> fig, ax = style.new_figure("single")
"""

from __future__ import annotations

from typing import Literal, Sequence
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
import seaborn as sns

# ---------------------------------------------------------------------------
# Figure size presets (inches). Aligned with untils journal guidance.
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

# SIMD quintile: 1 = most deprived (red) -> 5 = least (blue).
SIMD_QUINTILE_PALETTE: dict[int, str] = {
    1: "#d7191c",
    2: "#fdae61",
    3: "#ffffbf",
    4: "#abdda4",
    5: "#2b83ba",
}

# WHO variants of concern. Order follows emergence in Scotland.
WHO_VOC_PALETTE: dict[str, str] = {
    "None":        "#cccccc",
    "Alpha":       "#4e79a7",
    "Beta":        "#f28e2b",
    "Eta":         "#e15759",
    "Gamma":       "#59a14f",
    "Omicron":     "#edc948",
    "Kappa":       "#b07aa1",
    "Delta":       "#ff9da7",
    "Theta":       "#9c755f",
    "Mu":          "#af2d2d",
    "recombinant": "#bab0ac",
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


def set_theme(
    context: Literal["paper", "talk", "poster"] = "paper",
    font: str = "Arial",
    font_scale: float = 1.0,
) -> None:
    sns.set_theme(
        style="white",
        context=context,
        font_scale=font_scale,
        rc={
            "font.family": "sans-serif",
            "font.sans-serif": [font, "Arial", "Liberation Sans", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.03,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.fontsize": 8,
            "legend.title_fontsize": 8,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.5,
            "lines.markersize": 5,
            "patch.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "axes.grid": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "image.cmap": "viridis",
        },
    )


def save_figure(
    fig: plt.Figure,
    out_path: Path,
    width: Literal["single", "onehalf", "double", "slide"] = "single",
    *,
    height_in: float | None = None,
    dpi: int = 600,
    save_pdf: bool = True,
    save_png: bool = False,
    save_tiff: bool = False,
    save_eps: bool = False,
) -> dict[str, Path]:
    if not (300 <= dpi <= 600):
        raise ValueError("dpi should usually be between 300 and 600 for PLOS scripts.")
    if height_in is not None and height_in <= 0:
        raise ValueError("height_in must be positive when provided.")

    current_w, current_h = fig.get_size_inches()
    w = FIG_WIDTHS_IN[width]
    if height_in is None:
        aspect = current_h / current_w
        h = w * aspect
    else:
        h = height_in

    fig.set_size_inches(w, h, forward=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    saved_paths: dict[str, Path] = {}

    if save_pdf:
        pdf_path = out_path.with_suffix(".pdf")
        fig.savefig(pdf_path, dpi=dpi, transparent=False)
        saved_paths["pdf"] = pdf_path

    if save_png:
        png_path =out_path.with_suffix(".png")
        fig.savefig(png_path, dpi=dpi, transparent=False)
        saved_paths["png"] = png_path

    if save_eps:
        eps_path = out_path.with_suffix(".eps")
        fig.savefig(eps_path, format="eps", dpi=dpi, transparent=False)
        saved_paths["eps"] = eps_path

    if save_tiff:
        tif_path = out_path.with_suffix(".tif")
        fig.savefig(
            tif_path,
            format="tiff",
            dpi=dpi,
            transparent=False,
            pil_kwargs={"compression": "tiff_lzw"},
        )
        saved_paths["tiff"] = tif_path

    plt.close(fig)

    return saved_paths


def add_panel_labels(
    axes: Sequence[plt.Axes],
    *,
    x: float = 0,
    y: float = 1.1,
    size: float | str = "medium",
) -> None:
    import string
    labels = list(string.ascii_uppercase)
    for ax, label in zip(axes, labels):
        ax.text(
            x,
            y,
            label,
            transform=ax.transAxes,
            fontsize=size,
            fontweight="bold",
            va="top",
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


def lighten(color: str, amount: float = 0.4) -> tuple[float, float, float]:
    """Blend `color` toward white by `amount` in [0, 1]."""
    r, g, b = to_rgb(color)
    return r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount
