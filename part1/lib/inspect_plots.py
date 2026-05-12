"""Quick-look summary plots used during model checking.

These figures are intended for inspecting results just after a fitting run —
they are *not* the publication-ready manuscript figures (those live in
``manuscript/make_figures.py``).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import (
    MIXING_PREDICTOR_TERMS,
    MIXING_VARIABLES,
    PRIMARY_TERMS,
    TERM_LABELS,
)
from .data_prep import repo_root


# ---------------------------------------------------------------------------
# Environment + style helpers
# ---------------------------------------------------------------------------


def setup_matplotlib_cache() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)


def load_plot_style():
    setup_matplotlib_cache()
    import matplotlib

    matplotlib.use("Agg")
    root = repo_root()
    root_str = str(root)
    if root_str not in sys.path:
        sys.path.insert(0, root_str)
    from utils import style  # type: ignore[import]

    return style


def term_colours(style) -> dict[str, str]:
    palette = style.SIMD_DOMAIN_PALETTE
    return {
        "deprivation_z": palette["overall"],
        "index_deprivation_z": palette["overall"],
        "local_incidence_z": palette["income"],
        "local_seq_fraction_z": palette["employment"],
        "window_seq_fraction_z": palette["education"],
        "test_positivity_z": palette["health"],
        "log_cluster_size_z": palette["access"],
        "simd_excess_mixing_z": palette["crime"],
        "age_excess_mixing_z": palette["housing"],
        "sex_excess_mixing_z": palette["overall"],
        "profile_excess_mixing_z": palette["income"],
    }


# ---------------------------------------------------------------------------
# Quick-look figures
# ---------------------------------------------------------------------------


def plot_count_effects(
    results: pd.DataFrame,
    out_base: Path,
    primary_terms: list[str] | None = None,
) -> None:
    """Forest plot of count-model coefficients for the two main outcomes."""
    style = load_plot_style()
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator

    primary_terms = list(primary_terms or PRIMARY_TERMS)
    outcomes = ["cluster_size", "geographic_dispersion"]
    primary = results[
        results["outcome"].isin(outcomes) & results["term"].isin(primary_terms)
    ].copy()
    if primary.empty:
        return
    components = ["hurdle_binary", "positive_zero_truncated_count"]
    colours = term_colours(style)
    ratio_ticks = [0.8, 1.0, 1.5, 2.0, 3.0, 4.0]

    fig, axes = style.new_figure(
        width="double",
        height_in=4.4,
        nrows=len(outcomes),
        ncols=2,
        sharex=True,
        font_scale=0.85,
        layout="constrained",
    )
    term_offsets = np.linspace(-0.3, 0.3, len(primary_terms))
    term_positions = dict(zip(primary_terms, term_offsets))

    for i, outcome in enumerate(outcomes):
        for j, component in enumerate(components):
            ax = axes[i, j]
            sub = primary[(primary["outcome"] == outcome) & (primary["component"] == component)]
            for _, row in sub.iterrows():
                y = term_positions[row["term"]]
                ax.plot(
                    [row["ratio_ci_low"], row["ratio_ci_high"]],
                    [y, y],
                    color=colours[row["term"]],
                    linewidth=1.1,
                    solid_capstyle="round",
                )
                ax.scatter(
                    row["ratio"], y,
                    color=colours[row["term"]],
                    s=18, zorder=3,
                    label=TERM_LABELS[row["term"]],
                )
            ax.axvline(1, color="#666666", linewidth=0.8, linestyle="--")
            ax.set_xscale("log")
            ax.set_xlim(0.75, 4.5)
            ax.xaxis.set_major_locator(FixedLocator(ratio_ticks))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
            ax.xaxis.set_minor_locator(NullLocator())
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_yticks([])
            ax.grid(axis="x", color="#dddddd", linewidth=0.6)
            if i == 0:
                ax.set_title(
                    "Hurdle: any excess" if component == "hurdle_binary" else "Positive: ZTNB",
                )
            if j == 0:
                ax.set_ylabel(
                    {"cluster_size": "Size", "geographic_dispersion": "Datazones"}[outcome],
                    rotation=0, ha="right", va="center", labelpad=20,
                )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    fig.legend(
        unique.values(), unique.keys(),
        loc="lower center",
        bbox_to_anchor=(0.6, -0.065),
        ncol=3, frameon=False,
        columnspacing=1.4, handlelength=1.2,
    )
    fig.supxlabel("Adjusted ratio per 1 SD higher cluster-level covariate", x=0.6)
    style.save_figure(
        fig, out_base, width="double", height_in=4.8,
        dpi=600, save_pdf=True, save_png=True,
    )


def plot_mixing_predictor_count_effects(
    results: pd.DataFrame,
    out_base: Path,
) -> None:
    """Forest plot of Line-2 mixing-predictor count-model coefficients."""
    style = load_plot_style()
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator

    outcomes = ["cluster_size", "geographic_dispersion"]
    data = results[
        results["outcome"].isin(outcomes)
        & results["component"].isin(["hurdle_binary", "positive_zero_truncated_count"])
        & results["term"].isin(MIXING_PREDICTOR_TERMS)
    ].copy()
    if data.empty:
        return

    components = ["hurdle_binary", "positive_zero_truncated_count"]
    colours = term_colours(style)
    ci_min = float(data["ratio_ci_low"].min())
    ci_max = float(data["ratio_ci_high"].max())
    x_min = max(0.25, np.floor(ci_min * 10.0) / 10.0)
    x_max = min(5.0, np.ceil(ci_max * 10.0) / 10.0)
    x_min = min(x_min, 0.8)
    x_max = max(x_max, 1.2)
    tick_candidates = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
    ratio_ticks = [tick for tick in tick_candidates if x_min <= tick <= x_max]

    fig, axes = style.new_figure(
        width="double",
        height_in=4.4,
        nrows=len(outcomes),
        ncols=2,
        sharex=True,
        font_scale=0.85,
        layout="constrained",
    )
    term_offsets = np.linspace(-0.27, 0.27, len(MIXING_PREDICTOR_TERMS))
    term_positions = dict(zip(MIXING_PREDICTOR_TERMS, term_offsets))

    for i, outcome in enumerate(outcomes):
        for j, component in enumerate(components):
            ax = axes[i, j]
            sub = data[(data["outcome"] == outcome) & (data["component"] == component)]
            if not sub.empty:
                for _, row in sub.iterrows():
                    y = term_positions[row["term"]]
                    ax.plot(
                        [row["ratio_ci_low"], row["ratio_ci_high"]],
                        [y, y],
                        color=colours[row["term"]],
                        linewidth=1.1,
                        solid_capstyle="round",
                    )
                    ax.scatter(
                        row["ratio"], y,
                        color=colours[row["term"]],
                        s=18, zorder=3,
                        label=TERM_LABELS[row["term"]],
                    )
                ax.axvline(1, color="#666666", linewidth=0.8, linestyle="--")
                ax.set_xscale("log")
                ax.set_xlim(x_min, x_max)
                ax.xaxis.set_major_locator(FixedLocator(ratio_ticks))
                ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
                ax.xaxis.set_minor_locator(NullLocator())
                ax.xaxis.set_minor_formatter(NullFormatter())
                ax.set_yticks([])
                ax.grid(axis="x", color="#dddddd", linewidth=0.6)
            else:
                ax.text(
                    2.0, 0.5, "Not estimable",
                    ha="center", va="center",
                    fontsize=8, color="#666666",
                )
                ax.set_yticks([])
            if i == 0:
                ax.set_title(
                    "Hurdle: any excess" if component == "hurdle_binary" else "Positive: ZTNB",
                )
            if j == 0:
                ax.set_ylabel(
                    {"cluster_size": "Size", "geographic_dispersion": "Datazones"}[outcome],
                    rotation=0, ha="right", va="center", labelpad=20,
                )

    handles, labels = axes[1, 0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    fig.legend(
        unique.values(), unique.keys(),
        loc="lower center",
        bbox_to_anchor=(0.6, -0.05),
        ncol=4, frameon=False,
        columnspacing=1.1, handlelength=1.2,
    )
    fig.supxlabel("Adjusted ratio per 1 SD higher excess-mixing predictor", x=0.6)
    style.save_figure(
        fig, out_base, width="double", height_in=4.8,
        dpi=600, save_pdf=True, save_png=True,
    )


def plot_mixing_effects(
    results: pd.DataFrame,
    out_base: Path,
    primary_terms: list[str] | None = None,
) -> None:
    """Forest plot of Line-1 mixing-outcome coefficients (pp scale)."""
    import math as _math

    style = load_plot_style()
    from matplotlib.ticker import FuncFormatter, MultipleLocator

    primary_terms = list(primary_terms or PRIMARY_TERMS)
    models = ["simd", "age", "sex", "profile"]
    model_positions = {model: i for i, model in enumerate(models)}
    terms = primary_terms + ["log_cluster_size_z"]
    term_offsets = np.linspace(-0.32, 0.32, len(terms))
    term_positions = dict(zip(terms, term_offsets))
    colours = term_colours(style)

    max_abs = float(
        np.nanmax(
            np.abs(
                results[
                    [
                        "ci_low_percentage_points",
                        "coefficient_percentage_points",
                        "ci_high_percentage_points",
                    ]
                ].to_numpy()
            )
        )
    )
    x_limit = max(8, int(_math.ceil(max_abs / 2.0) * 2))

    fig, ax = style.new_figure(width="double", height_in=5, font_scale=0.85)
    for _, row in results.iterrows():
        y = model_positions[row["outcome"]] + term_positions[row["term"]]
        ax.plot(
            [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
            [y, y],
            color=colours[row["term"]],
            linewidth=1.2,
            solid_capstyle="round",
        )
        ax.scatter(
            row["coefficient_percentage_points"], y,
            color=colours[row["term"]],
            s=18, zorder=3,
            label=TERM_LABELS[row["term"]],
        )

    ax.axvline(0, color="#666666", linewidth=0.8, linestyle="--")
    ax.set_xlim(-x_limit, x_limit)
    ax.xaxis.set_major_locator(MultipleLocator(2 if x_limit <= 10 else 5))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.set_xlabel(
        "Change in excess pairwise discordance (pp per 1 SD higher covariate)",
        labelpad=5,
    )
    ax.set_yticks(list(model_positions.values()))
    ax.set_yticklabels([MIXING_VARIABLES[m]["short_label"] for m in models])
    ax.set_ylim(-0.6, len(models) - 0.4)
    ax.invert_yaxis()
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(
        unique.values(), unique.keys(),
        loc="upper center",
        bbox_to_anchor=(0.5, -0.17),
        ncol=3, frameon=False,
        columnspacing=1.4, handlelength=1.2,
    )
    fig.subplots_adjust(bottom=0.38, left=0.2, right=0.98)
    style.save_figure(
        fig, out_base, width="double", height_in=4.6,
        dpi=600, save_pdf=True, save_png=True,
    )
