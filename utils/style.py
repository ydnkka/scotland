"""Publication-grade plotting defaults.

The style is deliberately conservative: a sans-serif body font at ~8-9 pt,
thin axes, no top/right spines, and a single-column width tuned to *Lancet*
and *Virus Evolution* templates. Palettes are colour-blind friendly and
use perceptually uniform ramps where a quantitative variable is mapped.

Usage
-----
>>> from analysis.utils import style
>>> style.set_theme()
>>> fig, ax = style.new_figure("single")
"""

from __future__ import annotations

from typing import Literal, Sequence, Any
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import to_rgb
import seaborn as sns

__all__ = [
    "set_theme",
    "save_figure",
    "add_panel_labels",
    "new_figure",
    "lighten",
    "FIG_WIDTHS_IN",
    "WIDTHS",
    "CONTEXTS",
    "DEFAULT_HEIGHT_IN",
]

# ---------------------------------------------------------------------------
# Figure size presets (inches). Aligned with utils journal guidance.
# ---------------------------------------------------------------------------

FIG_WIDTHS_IN = {
    "single": 3.5,  # single column
    "onehalf": 5.2,  # 1.5-column
    "double": 7.2,  # full-page double column
    "slide": 10.0,  # for talks
}

WIDTHS = Literal["single", "onehalf", "double", "slide"]
CONTEXTS = Literal["paper", "talk", "poster", "notebook"]

DEFAULT_HEIGHT_IN = 2.6

# ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------


def set_theme(
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> None:
    """Apply the project matplotlib/seaborn theme for figures."""
    # Seaborn's context scaling factors relative to "paper"
    _context_scales = {"paper": 1.0, "notebook": 1.0, "talk": 1.3, "poster": 1.6}
    scale = _context_scales.get(context, 1.0) * font_scale

    sns.set_theme(
        style="white",
        context=context,  # still sets line/marker base scaling via seaborn internals
        font_scale=font_scale,
        rc={
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Arial", "Liberation Sans", "DejaVu Sans"],
            # --- everything below now scales ---
            "font.size": 9 * scale,
            "axes.titlesize": 10 * scale,
            "axes.labelsize": 9 * scale,
            "xtick.labelsize": 8 * scale,
            "ytick.labelsize": 8 * scale,
            "legend.fontsize": 8 * scale,
            "legend.title_fontsize": 8 * scale,
            "axes.linewidth": 0.8 * scale,
            "lines.linewidth": 1.5 * scale,
            "lines.markersize": 5 * scale,
            "patch.linewidth": 0.8 * scale,
            "xtick.major.width": 0.8 * scale,
            "ytick.major.width": 0.8 * scale,
            "xtick.major.size": 3.5 * scale,
            "ytick.major.size": 3.5 * scale,
            # --- these are booleans/strings, don't scale ---
            "mathtext.fontset": "dejavusans",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.03,
            "axes.spines.top": False,
            "axes.spines.right": False,
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
    fig: Figure,
    out_path: Path,
    width: WIDTHS = "single",
    *,
    width_in: float | None = None,
    height_in: float | None = None,
    dpi: int = 600,
    save_pdf: bool = True,
    save_png: bool = False,
    save_tiff: bool = False,
    save_eps: bool = False,
    save_svg: bool = False,
) -> dict[str, Path]:
    """Save a figure to the requested raster/vector formats."""
    if not (300 <= dpi <= 600):
        raise ValueError("dpi should usually be between 300 and 600 for PLOS scripts.")
    if height_in is not None and height_in <= 0:
        raise ValueError("height_in must be positive when provided.")

    current_w, current_h = fig.get_size_inches()

    w = width_in if width_in is not None else FIG_WIDTHS_IN[width]
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
        png_path = out_path.with_suffix(".png")
        fig.savefig(png_path, dpi=dpi, transparent=False)
        saved_paths["png"] = png_path

    if save_eps:
        eps_path = out_path.with_suffix(".eps")
        fig.savefig(eps_path, format="eps", dpi=dpi, transparent=False)
        saved_paths["eps"] = eps_path

    if save_svg:
        svg_path = out_path.with_suffix(".svg")
        fig.savefig(svg_path, format="eps", dpi=dpi, transparent=False)
        saved_paths["svg"] = svg_path

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
    axes: Sequence[Axes],
    *,
    x: float = 0,
    y: float = 1.1,
    size: float | str = "medium",
) -> None:
    """Place sequential panel labels on a list of axes."""
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
    width: WIDTHS = "single",
    height_in: float | None = None,
    width_in: float | None = None,
    nrows: int = 1,
    ncols: int = 1,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    **subplots_kwargs,
) -> tuple[Figure, Any]:
    """Create a figure with one of the paper-width presets.

    Parameters
    ----------
    width:
        Named width preset (ignored when *width_in* is given).
    height_in:
        Explicit figure height in inches.  Defaults to
        ``DEFAULT_HEIGHT_IN * nrows``.
    width_in:
        Explicit figure width in inches.  Overrides the *width* preset,
        useful for large multi-panel figures that need more room than the
        standard ``"double"`` (7.2 in) column allows.
    nrows, ncols:
        Passed to ``plt.subplots``.  The default is a single panel, but this
        is often overridden for multi-panel figures.
    context:
        Seaborn context to set for the figure.
        This controls the base scaling of fonts and lines.
        The default is ``"paper"``, but for talks or posters, you may want
        to use ``"talk"`` or ``"poster"`` for larger text.
    font_scale:
        Additional scaling factor for fonts and lines.  This is applied on top
        of the scaling from the *context* parameter, allowing for fine-tuning of
        text size without changing the overall context.
        The default is 1.0 (no additional scaling), but you can increase this for
        larger figures or decrease it for smaller ones.
    **subplots_kwargs:
        Additional keyword arguments passed to ``plt.subplots``.  This allows you to customize the subplots further, for example by setting ``sharex=True`` or ``gridspec_kw={"width            _ratios": [1, 2]}`` for uneven panels.  Note that the figure size is determined by the *width* and *height_in* parameters, so you should not set the *figsize* argument in *subplots_kwargs* when using this function.
    """
    set_theme(context=context, font_scale=font_scale)
    w = width_in if width_in is not None else FIG_WIDTHS_IN[width]
    h = height_in if height_in is not None else DEFAULT_HEIGHT_IN * nrows
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(w, h), **subplots_kwargs)
    return fig, ax


def lighten(color: str, amount: float = 0.4) -> tuple[float, float, float]:
    """Blend `color` toward white by `amount` in [0, 1]."""
    r, g, b = to_rgb(color)
    return r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount
