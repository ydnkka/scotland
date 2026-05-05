"""Create publication-ready figures for the Part 1 main analysis.

Outputs are written to ``part1/main/manuscript/figures`` as PDF, PNG, and TIFF.
The script uses the shared project plotting module at ``utils/style.py``.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


OUTCOME_LABELS = {
    "cluster_size": "Cluster size",
    "duration": "Duration",
    "geographic_dispersion": "Geographic spread",
    "duration_size_adjusted": "Duration, size-adjusted",
    "geographic_dispersion_size_adjusted": "Geographic spread, size-adjusted",
}

TERM_LABELS = {
    "deprivation_z": "SIMD deprivation",
    "local_incidence_z": "Local incidence",
    "local_seq_fraction_z": "Local sequencing",
    "window_seq_fraction_z": "Window sequencing",
    "test_positivity_z": "Test positivity",
    "log_cluster_size_z": "Cluster size",
}

PRIMARY_TERMS = [
    "deprivation_z",
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
]

MIXING_TERMS = PRIMARY_TERMS + ["log_cluster_size_z"]

MIXING_LABELS = {
    "simd": "SIMD",
    "age": "Age",
    "sex": "Sex",
    "profile": "Joint profile",
    "age_sex": "Joint age-sex",
}

COMPONENT_LABELS = {
    "hurdle_binary": "Hurdle odds",
    "positive_zero_truncated_count": "ZTNB count ratio",
    "log_linear": "Log-linear",
}

DOMAIN_ORDER = [
    "overall",
    "income",
    "employment",
    "education",
    "health",
    "access",
    "crime",
    "housing",
]

DOMAIN_LABELS = {
    "overall": "Overall",
    "income": "Income",
    "employment": "Employment",
    "education": "Education",
    "health": "Health",
    "access": "Access",
    "crime": "Crime",
    "housing": "Housing",
}

WAVE_ORDER = ["B.1.177", "Alpha", "Delta", "BA.1", "BA.2", "BA.4", "BA.5", "BQ.1"]


def repo_root(start: Path | None = None) -> Path:
    p = (start or Path(__file__)).resolve()
    for candidate in [p, *p.parents]:
        if (candidate / "config.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate config.yaml.")


def load_style(root: Path):
    sys.path.insert(0, str(root))
    from utils import style

    return style


def setup_environment() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)


def term_colours(style) -> dict[str, str]:
    palette = style.SIMD_DOMAIN_PALETTE
    return {
        "deprivation_z": palette["overall"],
        "local_incidence_z": palette["income"],
        "local_seq_fraction_z": palette["employment"],
        "window_seq_fraction_z": palette["education"],
        "test_positivity_z": palette["health"],
        "log_cluster_size_z": palette["access"],
    }


def save_all(style, fig, out_base: Path, width: str, height_in: float) -> None:
    style.save_figure(
        fig,
        out_base,
        width=width,
        height_in=height_in,
        dpi=600,
        save_pdf=True,
        save_png=True,
        save_tiff=True,
    )


def remove_stale_main_duplicates(out_dir: Path) -> None:
    """Drop old main-figure exports now represented only as supplements."""
    for stem in ["fig4_wave_specific_domain_demographic_mixing"]:
        for suffix in [".pdf", ".png", ".tif"]:
            path = out_dir / f"{stem}{suffix}"
            if path.exists():
                path.unlink()


def draw_ratio_panel(
    ax,
    df: pd.DataFrame,
    terms: list[str],
    colours: dict[str, str],
    *,
    title: str,
    show_ylabels: bool,
    xlim: tuple[float, float],
    xlabel: str | None = None,
) -> None:
    from matplotlib.ticker import NullFormatter, NullLocator

    y_positions = np.arange(len(terms))[::-1]
    position = dict(zip(terms, y_positions))
    for term in terms:
        row = df[df["term"] == term]
        if row.empty:
            continue
        row = row.iloc[0]
        y = position[term]
        ax.plot(
            [row["ratio_ci_low"], row["ratio_ci_high"]],
            [y, y],
            color=colours[term],
            linewidth=1.1,
            solid_capstyle="round",
        )
        ax.scatter(
            row["ratio"],
            y,
            color=colours[term],
            edgecolor="white",
            linewidth=0.3,
            s=18,
            zorder=3,
        )

    ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    if xlim[1] >= 3.9:
        ticks = [0.8, 1.0, 1.5, 2.0, 3.0, 4.0]
    else:
        ticks = [0.9, 1.0, 1.5, 2.0, 3.0]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{tick:.1f}" for tick in ticks])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_title(title, pad=4)
    ax.set_ylim(-0.7, len(terms) - 0.3)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([TERM_LABELS[t] for t in terms] if show_ylabels else [])
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)


def draw_difference_panel(
    ax,
    df: pd.DataFrame,
    terms: list[str],
    colours: dict[str, str],
    *,
    title: str,
    show_ylabels: bool,
    xlim: tuple[float, float],
    xlabel: str | None = None,
) -> None:
    y_positions = np.arange(len(terms))[::-1]
    position = dict(zip(terms, y_positions))
    for term in terms:
        row = df[df["term"] == term]
        if row.empty:
            continue
        row = row.iloc[0]
        y = position[term]
        ax.plot(
            [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
            [y, y],
            color=colours[term],
            linewidth=1.1,
            solid_capstyle="round",
        )
        ax.scatter(
            row["coefficient_percentage_points"],
            y,
            color=colours[term],
            edgecolor="white",
            linewidth=0.3,
            s=18,
            zorder=3,
        )

    ax.axvline(0.0, color="#666666", linestyle="--", linewidth=0.7)
    ax.set_xlim(*xlim)
    ax.set_title(title, pad=4)
    ax.set_ylim(-0.7, len(terms) - 0.3)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([TERM_LABELS[t] for t in terms] if show_ylabels else [])
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)


def plot_main_count_results(style, count_results: pd.DataFrame, out_dir: Path) -> None:
    colours = term_colours(style)
    outcomes = ["cluster_size", "duration", "geographic_dispersion"]
    components = ["hurdle_binary", "positive_zero_truncated_count"]

    fig, axes = style.new_figure(
        width="double",
        height_in=6.0,
        nrows=3,
        ncols=2,
        sharex=True,
        font_scale=0.85,
    )

    for i, outcome in enumerate(outcomes):
        for j, component in enumerate(components):
            ax = axes[i, j]
            sub = count_results[
                (count_results["outcome"] == outcome)
                & (count_results["component"] == component)
                & (count_results["term"].isin(PRIMARY_TERMS))
            ]
            title = f"{OUTCOME_LABELS[outcome]}: {COMPONENT_LABELS[component]}"
            draw_ratio_panel(
                ax,
                sub,
                PRIMARY_TERMS,
                colours,
                title=title,
                show_ylabels=j == 0,
                xlim=(0.75, 4.0),
                xlabel=(
                    "Odds ratio per 1 SD higher covariate"
                    if component == "hurdle_binary" and i == 2
                    else "ZTNB count ratio per 1 SD higher covariate"
                    if component == "positive_zero_truncated_count" and i == 2
                    else None
                ),
            )

    style.add_panel_labels(axes.ravel(), x=-0.18, y=1.14, size=9)
    fig.subplots_adjust(left=0.21, right=0.98, top=0.94, bottom=0.10, hspace=0.42, wspace=0.14)
    save_all(style, fig, out_dir / "fig1_main_cluster_outcomes", "double", 6.0)


def plot_main_mixing_results(style, mixing_results: pd.DataFrame, out_dir: Path) -> None:
    colours = term_colours(style)
    outcomes = ["simd", "age", "sex", "profile"]
    fig, axes = style.new_figure(
        width="double",
        height_in=5.0,
        nrows=2,
        ncols=2,
        sharex=True,
        font_scale=0.85,
    )

    for idx, outcome in enumerate(outcomes):
        ax = axes.ravel()[idx]
        sub = mixing_results[mixing_results["outcome"] == outcome]
        draw_difference_panel(
            ax,
            sub,
            MIXING_TERMS,
            colours,
            title=f"{MIXING_LABELS[outcome]} excess mixing",
            show_ylabels=idx % 2 == 0,
            xlim=(-8.5, 8.5),
        )

    style.add_panel_labels(axes.ravel(), x=-0.18, y=1.13, size=9)
    fig.supxlabel("Change in excess mixing (pp per 1 SD higher covariate)", y=0.04, fontsize=8)
    fig.subplots_adjust(left=0.20, right=0.98, top=0.93, bottom=0.14, hspace=0.36, wspace=0.14)
    save_all(style, fig, out_dir / "fig2_main_cluster_mixing", "double", 5.0)


def binned_percent(values: pd.Series, bins: list[float], labels: list[str]) -> pd.DataFrame:
    cats = pd.cut(values, bins=bins, labels=labels, include_lowest=True, right=True)
    pct = cats.value_counts(sort=False, normalize=True).mul(100)
    out = pct.rename("percent").reset_index()
    return out.rename(columns={out.columns[0]: "bin"})


def plot_outcome_distributions(style, cluster_table: pd.DataFrame, out_dir: Path) -> None:
    grey = "#6f6f6f"
    specs = [
        (
            "cluster_size",
            "Cluster size",
            [-np.inf, 1.5, 2.5, 3.5, 5.5, 10.5, 20.5, 50.5, np.inf],
            ["1", "2", "3", "4-5", "6-10", "11-20", "21-50", ">50"],
        ),
        (
            "duration_days",
            "Duration (days)",
            [-np.inf, 0.5, 1.5, 2.5, 5.5, 10.5, 15.5, np.inf],
            ["0", "1", "2", "3-5", "6-10", "11-15", ">15"],
        ),
        (
            "cluster_n_datazones",
            "Distinct datazones",
            [-np.inf, 1.5, 2.5, 3.5, 5.5, 10.5, 20.5, 50.5, np.inf],
            ["1", "2", "3", "4-5", "6-10", "11-20", "21-50", ">50"],
        ),
    ]

    fig, axes = style.new_figure(
        width="double",
        height_in=2.45,
        nrows=1,
        ncols=3,
        font_scale=0.85,
    )

    for ax, (col, title, bins, labels) in zip(axes, specs):
        data = binned_percent(cluster_table[col], bins, labels)
        ax.bar(data["bin"].astype(str), data["percent"], color=grey, width=0.78)
        ax.set_title(title, pad=4)
        ax.set_ylabel("Clusters (%)" if ax is axes[0] else "")
        ax.set_ylim(0, max(65, data["percent"].max() * 1.15))
        ax.tick_params(axis="x", rotation=45)
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes, x=-0.18, y=1.14, size=9)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.86, bottom=0.28, wspace=0.26)
    save_all(style, fig, out_dir / "supp_fig1_outcome_distributions", "double", 2.45)


def plot_size_adjusted_sensitivity(style, count_results: pd.DataFrame, out_dir: Path) -> None:
    colours = term_colours(style)
    outcomes = ["duration_size_adjusted", "geographic_dispersion_size_adjusted"]
    fig, axes = style.new_figure(
        width="double",
        height_in=3.0,
        nrows=1,
        ncols=2,
        sharex=True,
        font_scale=0.85,
    )
    for idx, outcome in enumerate(outcomes):
        sub = count_results[
            (count_results["outcome"] == outcome)
            & (count_results["component"] == "positive_zero_truncated_count")
        ]
        draw_ratio_panel(
            axes[idx],
            sub,
            MIXING_TERMS,
            colours,
            title=OUTCOME_LABELS[outcome],
            show_ylabels=idx == 0,
            xlim=(0.85, 3.2),
            xlabel="ZTNB count ratio per 1 SD higher covariate",
        )
    style.add_panel_labels(axes, x=-0.18, y=1.14, size=9)
    fig.subplots_adjust(left=0.20, right=0.98, top=0.86, bottom=0.17, wspace=0.12)
    save_all(style, fig, out_dir / "supp_fig2_size_adjusted_positive_counts", "double", 3.0)


def plot_loglinear_comparison(
    style,
    count_results: pd.DataFrame,
    loglinear_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    from matplotlib.ticker import NullFormatter, NullLocator

    outcomes = ["cluster_size", "duration", "geographic_dispersion"]
    model_colours = {
        "Log-linear": "#666666",
        "Hurdle": "#4e79a7",
        "ZTNB positive": "#f28e2b",
    }
    markers = {"Log-linear": "o", "Hurdle": "s", "ZTNB positive": "^"}

    log_map = {
        "cluster_size": "cluster_size",
        "duration": "duration",
        "geographic_dispersion": "geographic_dispersion",
    }
    pieces = []
    log = loglinear_results[
        loglinear_results["model"].isin(log_map)
        & loglinear_results["term"].isin(PRIMARY_TERMS)
    ].copy()
    log["outcome"] = log["model"].map(log_map)
    log["model_type"] = "Log-linear"
    log = log.rename(
        columns={
            "geometric_mean_ratio": "ratio",
            "ci_low": "ratio_ci_low",
            "ci_high": "ratio_ci_high",
        }
    )
    pieces.append(log[["outcome", "term", "model_type", "ratio", "ratio_ci_low", "ratio_ci_high"]])

    hurdle = count_results[
        (count_results["outcome"].isin(outcomes))
        & (count_results["component"] == "hurdle_binary")
        & (count_results["term"].isin(PRIMARY_TERMS))
    ].copy()
    hurdle["model_type"] = "Hurdle"
    pieces.append(hurdle[["outcome", "term", "model_type", "ratio", "ratio_ci_low", "ratio_ci_high"]])

    ztnb = count_results[
        (count_results["outcome"].isin(outcomes))
        & (count_results["component"] == "positive_zero_truncated_count")
        & (count_results["term"].isin(PRIMARY_TERMS))
    ].copy()
    ztnb["model_type"] = "ZTNB positive"
    pieces.append(ztnb[["outcome", "term", "model_type", "ratio", "ratio_ci_low", "ratio_ci_high"]])
    comp = pd.concat(pieces, ignore_index=True)

    fig, axes = style.new_figure(
        width="double",
        height_in=4.2,
        nrows=1,
        ncols=3,
        sharex=True,
        font_scale=0.85,
    )
    offsets = {"Log-linear": -0.18, "Hurdle": 0.0, "ZTNB positive": 0.18}
    y_positions = np.arange(len(PRIMARY_TERMS))[::-1]
    pos = dict(zip(PRIMARY_TERMS, y_positions))
    for idx, outcome in enumerate(outcomes):
        ax = axes[idx]
        sub = comp[comp["outcome"] == outcome]
        for model_type in ["Log-linear", "Hurdle", "ZTNB positive"]:
            model_sub = sub[sub["model_type"] == model_type]
            for _, row in model_sub.iterrows():
                y = pos[row["term"]] + offsets[model_type]
                ax.plot(
                    [row["ratio_ci_low"], row["ratio_ci_high"]],
                    [y, y],
                    color=model_colours[model_type],
                    linewidth=0.9,
                    solid_capstyle="round",
                )
                ax.scatter(
                    row["ratio"],
                    y,
                    color=model_colours[model_type],
                    marker=markers[model_type],
                    s=17,
                    zorder=3,
                    label=model_type,
                )
        ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_xscale("log")
        # Extend the right xlim slightly so the "4.0" tick label doesn't clip
        ax.set_xlim(0.75, 4.5)
        ticks = [0.8, 1.0, 1.5, 2.0, 3.0, 4.0]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{tick:.1f}" for tick in ticks])
        # Rotate tick labels to prevent crowding on narrow panels
        ax.tick_params(axis="x", rotation=35)
        for label in ax.get_xticklabels():
            label.set_ha("right")
            label.set_rotation_mode("anchor")
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_title(OUTCOME_LABELS[outcome], pad=4)
        ax.set_ylim(-0.8, len(PRIMARY_TERMS) - 0.2)
        ax.set_yticks(y_positions)
        ax.set_yticklabels([TERM_LABELS[t] for t in PRIMARY_TERMS] if idx == 0 else [])
        # No per-panel xlabel — use shared supxlabel below
        ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    # Single shared x-axis label (all panels share the same x axis)
    fig.supxlabel(
        "Model-specific ratio per 1 SD higher covariate",
        y=0.01,
        fontsize=8,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    fig.legend(
        unique.values(),
        unique.keys(),
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.54, 0.07),
    )
    style.add_panel_labels(axes, x=-0.18, y=1.12, size=9)
    # Extra bottom space for rotated tick labels + shared supxlabel + legend
    fig.subplots_adjust(left=0.18, right=0.99, top=0.87, bottom=0.28, wspace=0.18)
    save_all(style, fig, out_dir / "supp_fig3_loglinear_vs_hurdle_ztnb", "double", 4.2)


def domain_effect_rows(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["term"].astype(str).eq(df["domain"].astype(str) + "_deprivation_z")].copy()


def plot_simd_domain_outcomes(
    style,
    domain_outcomes: pd.DataFrame,
    out_dir: Path,
) -> None:
    from matplotlib.ticker import NullFormatter, NullLocator

    colours = style.SIMD_DOMAIN_PALETTE
    outcomes = ["cluster_size", "duration", "geographic_dispersion"]
    components = ["hurdle_binary", "positive_zero_truncated_count"]
    data = domain_effect_rows(domain_outcomes)
    ci_min = float(data["ratio_ci_low"].min())
    ci_max = float(data["ratio_ci_high"].max())
    xlim = (max(0.5, ci_min * 0.95), min(1.5, ci_max * 1.05))
    if xlim[0] > 0.9:
        xlim = (0.9, xlim[1])
    if xlim[1] < 1.1:
        xlim = (xlim[0], 1.1)

    fig, axes = style.new_figure(
        width="double",
        height_in=6.0,
        nrows=3,
        ncols=2,
        sharex=True,
        font_scale=0.85,
    )
    y_positions = np.arange(len(DOMAIN_ORDER))[::-1]
    pos = dict(zip(DOMAIN_ORDER, y_positions))
    for idx, outcome in enumerate(outcomes):
        for jdx, component in enumerate(components):
            ax = axes[idx, jdx]
            sub = data[(data["outcome"] == outcome) & (data["component"] == component)]
            for domain in DOMAIN_ORDER:
                row = sub[sub["domain"] == domain]
                if row.empty:
                    continue
                row = row.iloc[0]
                y = pos[domain]
                ax.plot(
                    [row["ratio_ci_low"], row["ratio_ci_high"]],
                    [y, y],
                    color=colours[domain],
                    linewidth=1.0,
                    solid_capstyle="round",
                )
                ax.scatter(
                    row["ratio"],
                    y,
                    color=colours[domain],
                    edgecolor="white",
                    linewidth=0.3,
                    s=18,
                    zorder=3,
                )
            ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
            ax.set_xscale("log")
            ax.set_xlim(*xlim)
            ticks = [tick for tick in [0.8, 0.9, 1.0, 1.1, 1.2] if xlim[0] <= tick <= xlim[1]]
            ax.set_xticks(ticks)
            ax.set_xticklabels([f"{tick:.1f}" for tick in ticks])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_title(f"{OUTCOME_LABELS[outcome]}: {COMPONENT_LABELS[component]}", pad=4)
            ax.set_yticks(y_positions)
            ax.set_yticklabels([DOMAIN_LABELS[d] for d in DOMAIN_ORDER] if jdx == 0 else [])
            ax.set_xlabel(
                (
                    "Odds ratio per 1 SD higher domain deprivation"
                    if component == "hurdle_binary"
                    else "ZTNB count ratio per 1 SD higher domain deprivation"
                )
                if idx == len(outcomes) - 1
                else ""
            )
            ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes.ravel(), x=-0.18, y=1.12, size=9)
    fig.subplots_adjust(left=0.18, right=0.99, top=0.93, bottom=0.09, hspace=0.42, wspace=0.12)
    save_all(style, fig, out_dir / "supp_fig4_simd_domain_cluster_outcomes", "double", 6.0)


def plot_main_domain_mixing_results(
    style,
    domain_mixing: pd.DataFrame,
    domain_demo: pd.DataFrame,
    out_dir: Path,
) -> None:
    colours = style.SIMD_DOMAIN_PALETTE
    domain_quintile = domain_effect_rows(domain_mixing).copy()
    domain_quintile["panel"] = "domain_quintile"
    demo = domain_effect_rows(domain_demo).copy()
    demo["panel"] = demo["mixing"]
    data = pd.concat([domain_quintile, demo], ignore_index=True)
    panels = [
        ("domain_quintile", "Domain-quintile mixing"),
        ("age", "Age mixing"),
        ("sex", "Sex mixing"),
        ("age_sex", "Joint age-sex mixing"),
    ]
    ci_min = float(data["ci_low_percentage_points"].min())
    ci_max = float(data["ci_high_percentage_points"].max())
    limit = max(2.5, np.ceil(max(abs(ci_min), abs(ci_max)) * 2) / 2)
    xlim = (-limit, limit)

    fig, axes = style.new_figure(
        width="double",
        height_in=5.0,
        nrows=2,
        ncols=2,
        sharex=True,
        font_scale=0.85,
    )
    y_positions = np.arange(len(DOMAIN_ORDER))[::-1]
    pos = dict(zip(DOMAIN_ORDER, y_positions))

    for idx, (panel, title) in enumerate(panels):
        ax = axes.ravel()[idx]
        sub = data[data["panel"] == panel]
        for domain in DOMAIN_ORDER:
            row = sub[sub["domain"] == domain]
            if row.empty:
                continue
            row = row.iloc[0]
            y = pos[domain]
            ax.plot(
                [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
                [y, y],
                color=colours[domain],
                linewidth=1.1,
                solid_capstyle="round",
            )
            ax.scatter(
                row["coefficient_percentage_points"],
                y,
                color=colours[domain],
                edgecolor="white",
                linewidth=0.3,
                s=19,
                zorder=3,
            )
        ax.axvline(0.0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_xlim(*xlim)
        ax.set_title(title, pad=4)
        ax.set_yticks(y_positions)
        ax.set_yticklabels([DOMAIN_LABELS[d] for d in DOMAIN_ORDER] if idx % 2 == 0 else [])
        ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes.ravel(), x=-0.16, y=1.12, size=9)
    fig.supxlabel("Change in excess mixing (pp per 1 SD higher domain deprivation)", y=0.04, fontsize=8)
    fig.subplots_adjust(left=0.18, right=0.99, top=0.90, bottom=0.14, hspace=0.32, wspace=0.12)
    save_all(style, fig, out_dir / "fig3_simd_domain_mixing", "double", 5.0)


def plot_wave_cluster_outcomes(
    style,
    wave_count_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    from matplotlib.ticker import NullFormatter, NullLocator

    colour = style.SIMD_DOMAIN_PALETTE["overall"]
    outcomes = ["cluster_size", "duration", "geographic_dispersion"]
    components = ["hurdle_binary", "positive_zero_truncated_count"]
    data = wave_count_results[
        wave_count_results["term"].eq("deprivation_z")
        & wave_count_results["outcome"].isin(outcomes)
        & wave_count_results["component"].isin(components)
    ].copy()
    waves = [wave for wave in WAVE_ORDER if wave in set(data["wave_group"])]
    if data.empty or not waves:
        return

    xlims: dict[str, tuple[float, float]] = {}
    ticks: dict[str, list[float]] = {}
    for component in components:
        sub = data[data["component"] == component]
        ci_min = float(sub["ratio_ci_low"].min())
        ci_max = float(sub["ratio_ci_high"].max())
        lower = max(0.5, np.floor(ci_min * 10) / 10)
        upper = min(3.0, np.ceil(ci_max * 10) / 10)
        if component == "hurdle_binary":
            lower = min(lower, 0.8)
            upper = max(upper, 1.2)
            tick_candidates = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
        else:
            lower = min(lower, 0.6)
            upper = max(upper, 2.7)
            tick_candidates = [0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
        xlims[component] = (lower, upper)
        ticks[component] = [tick for tick in tick_candidates if lower <= tick <= upper]

    fig, axes = style.new_figure(
        width="double",
        height_in=6.0,
        nrows=3,
        ncols=2,
        sharex=False,
        font_scale=0.85,
    )
    y_positions = np.arange(len(waves))[::-1]
    pos = dict(zip(waves, y_positions))

    for idx, outcome in enumerate(outcomes):
        for jdx, component in enumerate(components):
            ax = axes[idx, jdx]
            sub = data[(data["outcome"] == outcome) & (data["component"] == component)]
            for wave in waves:
                row = sub[sub["wave_group"] == wave]
                if row.empty:
                    continue
                row = row.iloc[0]
                y = pos[wave]
                ax.plot(
                    [row["ratio_ci_low"], row["ratio_ci_high"]],
                    [y, y],
                    color=colour,
                    linewidth=1.1,
                    solid_capstyle="round",
                )
                ax.scatter(
                    row["ratio"],
                    y,
                    color=colour,
                    edgecolor="white",
                    linewidth=0.3,
                    s=20,
                    zorder=3,
                )
            ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
            ax.set_xscale("log")
            ax.set_xlim(*xlims[component])
            ax.set_xticks(ticks[component])
            ax.set_xticklabels([f"{tick:g}" for tick in ticks[component]])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_title(f"{OUTCOME_LABELS[outcome]}: {COMPONENT_LABELS[component]}", pad=4)
            ax.set_yticks(y_positions)
            ax.set_yticklabels(waves if jdx == 0 else [])
            if idx == len(outcomes) - 1:
                ax.set_xlabel(
                    "Odds ratio per 1 SD higher SIMD deprivation"
                    if component == "hurdle_binary"
                    else "ZTNB count ratio per 1 SD higher SIMD deprivation"
                )
            ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes.ravel(), x=-0.16, y=1.12, size=9)
    fig.subplots_adjust(left=0.13, right=0.99, top=0.93, bottom=0.10, hspace=0.42, wspace=0.13)
    save_all(style, fig, out_dir / "fig4_wave_specific_cluster_outcomes", "double", 6.0)


def plot_simd_domain_quintile_mixing(
    style,
    domain_mixing: pd.DataFrame,
    out_dir: Path,
) -> None:
    colours = style.SIMD_DOMAIN_PALETTE
    data = domain_effect_rows(domain_mixing)
    fig, ax = style.new_figure(width="onehalf", height_in=3.4, font_scale=0.9)
    y_positions = np.arange(len(DOMAIN_ORDER))[::-1]
    pos = dict(zip(DOMAIN_ORDER, y_positions))
    for domain in DOMAIN_ORDER:
        row = data[data["domain"] == domain]
        if row.empty:
            continue
        row = row.iloc[0]
        y = pos[domain]
        ax.plot(
            [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
            [y, y],
            color=colours[domain],
            linewidth=1.2,
            solid_capstyle="round",
        )
        ax.scatter(
            row["coefficient_percentage_points"],
            y,
            color=colours[domain],
            edgecolor="white",
            linewidth=0.3,
            s=24,
            zorder=3,
        )
    ax.axvline(0.0, color="#666666", linestyle="--", linewidth=0.7)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([DOMAIN_LABELS[d] for d in DOMAIN_ORDER])
    # Wrap the long xlabel so it fits within the 1.5-column (5.2 inch) figure
    ax.set_xlabel(
        "Change in domain-quintile excess mixing\n(pp per 1 SD higher domain deprivation)",
        labelpad=6,
    )
    ax.set_xlim(-3.2, 0.8)
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)
    # Increase bottom margin to accommodate the two-line xlabel
    fig.subplots_adjust(left=0.26, right=0.98, top=0.96, bottom=0.22)
    save_all(style, fig, out_dir / "supp_fig5_simd_domain_quintile_mixing", "onehalf", 3.4)


def plot_domain_demographic_mixing(
    style,
    domain_demo: pd.DataFrame,
    out_dir: Path,
) -> None:
    colours = style.SIMD_DOMAIN_PALETTE
    data = domain_effect_rows(domain_demo)
    mixings = ["age", "sex", "age_sex"]
    fig, axes = style.new_figure(
        width="double",
        height_in=3.8,
        nrows=1,
        ncols=3,
        sharex=True,
        font_scale=0.85,
    )
    y_positions = np.arange(len(DOMAIN_ORDER))[::-1]
    pos = dict(zip(DOMAIN_ORDER, y_positions))
    for idx, mixing in enumerate(mixings):
        ax = axes[idx]
        sub = data[data["mixing"] == mixing]
        for domain in DOMAIN_ORDER:
            row = sub[sub["domain"] == domain]
            if row.empty:
                continue
            row = row.iloc[0]
            y = pos[domain]
            ax.plot(
                [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
                [y, y],
                color=colours[domain],
                linewidth=1.1,
                solid_capstyle="round",
            )
            ax.scatter(
                row["coefficient_percentage_points"],
                y,
                color=colours[domain],
                edgecolor="white",
                linewidth=0.3,
                s=19,
                zorder=3,
            )
        ax.axvline(0.0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_xlim(-2.0, 2.4)
        ax.set_title(f"{MIXING_LABELS[mixing]} mixing", pad=4)
        ax.set_yticks(y_positions)
        ax.set_yticklabels([DOMAIN_LABELS[d] for d in DOMAIN_ORDER] if idx == 0 else [])
        ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes, x=-0.18, y=1.12, size=9)
    fig.supxlabel("Change in excess mixing (pp per 1 SD higher domain deprivation)", y=0.04, fontsize=8)
    fig.subplots_adjust(left=0.18, right=0.99, top=0.86, bottom=0.18, wspace=0.13)
    save_all(style, fig, out_dir / "supp_fig6_simd_domain_demographic_mixing", "double", 3.8)


def plot_wave_domain_demographic_mixing(
    style,
    wave_domain_demo: pd.DataFrame,
    out_dir: Path,
    out_name: str = "supp_fig7_wave_specific_domain_demographic_mixing",
) -> None:
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.pyplot as plt

    data = domain_effect_rows(wave_domain_demo)
    mixings = ["age", "sex", "age_sex"]
    waves = [wave for wave in WAVE_ORDER if wave in set(data["wave_group"])]
    fig, axes = style.new_figure(
        width="double",
        height_in=6.2,
        nrows=3,
        ncols=1,
        font_scale=0.82,
    )
    vmax = max(5.5, np.nanmax(np.abs(data["coefficient_percentage_points"])))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)

    image = None
    for idx, mixing in enumerate(mixings):
        ax = axes[idx]
        sub = data[data["mixing"] == mixing]
        matrix = (
            sub.pivot_table(
                index="domain",
                columns="wave_group",
                values="coefficient_percentage_points",
                aggfunc="first",
            )
            .reindex(index=DOMAIN_ORDER, columns=waves)
        )
        image = ax.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", norm=norm, aspect="auto")
        ax.set_title(f"{MIXING_LABELS[mixing]} mixing", pad=4)
        ax.set_yticks(np.arange(len(DOMAIN_ORDER)))
        ax.set_yticklabels([DOMAIN_LABELS[d] for d in DOMAIN_ORDER])
        ax.set_xticks(np.arange(len(waves)))
        ax.set_xticklabels(waves, rotation=35, ha="right")
        ax.tick_params(length=0)
        for y in np.arange(len(DOMAIN_ORDER) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.6)
        for x in np.arange(len(waves) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.6)

    assert image is not None
    fig.subplots_adjust(left=0.17, right=0.84, top=0.93, bottom=0.12, hspace=0.46)
    cbar_ax = fig.add_axes([0.875, 0.20, 0.022, 0.62])
    cbar = fig.colorbar(image, cax=cbar_ax)
    cbar.set_label("pp per 1 SD higher domain deprivation")
    style.add_panel_labels(axes, x=-0.14, y=1.10, size=9)
    save_all(style, fig, out_dir / out_name, "double", 6.2)
    plt.close("all")


def plot_observed_expected_matrices(
    style,
    matrices: pd.DataFrame,
    out_dir: Path,
) -> None:
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.pyplot as plt

    def category_key(value: object) -> tuple[int, str]:
        text = str(value)
        if text == "75+":
            return (75, text)
        digits = "".join(ch if ch.isdigit() else " " for ch in text).split()
        return (int(digits[0]) if digits else 999, text)

    specs = [
        ("simd", "SIMD quintile"),
        ("age", "Age band"),
    ]
    fig, axes = style.new_figure(
        width="double",
        height_in=3.7,
        nrows=1,
        ncols=2,
        font_scale=0.8,
    )
    overall = matrices[matrices["wave_group"] == "Overall"].copy()
    vmax = max(2.5, np.nanmax(np.abs(overall["excess_percentage_points"])))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    image = None
    for idx, (variable, title) in enumerate(specs):
        ax = axes[idx]
        sub = overall[overall["variable"] == variable]
        row_order = sorted(sub["category_i"].astype(str).unique(), key=category_key)
        col_order = sorted(sub["category_j"].astype(str).unique(), key=category_key)
        matrix = (
            sub.assign(category_i=sub["category_i"].astype(str), category_j=sub["category_j"].astype(str))
            .pivot_table(
                index="category_i",
                columns="category_j",
                values="excess_percentage_points",
                aggfunc="first",
            )
            .reindex(index=row_order, columns=col_order)
        )
        image = ax.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", norm=norm, aspect="equal")
        ax.set_title(f"{title}: observed - expected", pad=4)
        ax.set_xticks(np.arange(len(col_order)))
        ax.set_yticks(np.arange(len(row_order)))
        ax.set_xticklabels(col_order, rotation=45, ha="right")
        ax.set_yticklabels(row_order)
        ax.set_xlabel(title)
        ax.set_ylabel(title if idx == 0 else "")
        ax.tick_params(length=0)
        for y in np.arange(len(row_order) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.5)
        for x in np.arange(len(col_order) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.5)

    assert image is not None
    fig.subplots_adjust(left=0.08, right=0.84, top=0.87, bottom=0.21, wspace=0.18)
    cbar_ax = fig.add_axes([0.875, 0.28, 0.024, 0.50])
    cbar = fig.colorbar(image, cax=cbar_ax)
    cbar.set_label("Observed - expected pair probability (pp)")
    style.add_panel_labels(axes, x=-0.14, y=1.12, size=9)
    save_all(style, fig, out_dir / "supp_fig8_observed_expected_mixing_matrices", "double", 3.7)
    plt.close("all")


def run(
    root: Path,
    tables_dir: Path | None = None,
    out_dir: Path | None = None,
    cache_dir: Path | None = None,
) -> None:
    setup_environment()
    style = load_style(root)

    main_dir = root / "part1" / "main"
    if tables_dir is None:
        tables_dir = main_dir / "tables"
    if out_dir is None:
        out_dir = main_dir / "manuscript" / "figures"
    if cache_dir is None:
        cache_dir = main_dir / "cache"

    out_dir.mkdir(parents=True, exist_ok=True)
    # Only remove stale duplicates when writing to the canonical figures directory
    if out_dir == main_dir / "manuscript" / "figures":
        remove_stale_main_duplicates(out_dir)

    count_results = pd.read_csv(tables_dir / "main_hurdle_count_model_results.csv")
    mixing_results = pd.read_csv(tables_dir / "main_mixing_model_results.csv")
    cluster_table = pd.read_parquet(cache_dir / "main_cluster_table.parquet")

    plot_main_count_results(style, count_results, out_dir)
    plot_main_mixing_results(style, mixing_results, out_dir)
    plot_outcome_distributions(style, cluster_table, out_dir)
    plot_size_adjusted_sensitivity(style, count_results, out_dir)

    domain_mixing = None
    domain_demo = None

    domain_outcome_path = tables_dir / "main_simd_domain_hurdle_count_model_results.csv"
    if domain_outcome_path.exists():
        plot_simd_domain_outcomes(style, pd.read_csv(domain_outcome_path), out_dir)

    domain_mixing_path = tables_dir / "main_simd_domain_quintile_mixing_model_results.csv"
    if domain_mixing_path.exists():
        domain_mixing = pd.read_csv(domain_mixing_path)
        plot_simd_domain_quintile_mixing(style, domain_mixing, out_dir)

    domain_demo_path = tables_dir / "main_simd_domain_demographic_mixing_model_results.csv"
    if domain_demo_path.exists():
        domain_demo = pd.read_csv(domain_demo_path)
        plot_domain_demographic_mixing(style, domain_demo, out_dir)

    if domain_mixing is not None and domain_demo is not None:
        plot_main_domain_mixing_results(style, domain_mixing, domain_demo, out_dir)

    wave_count_path = tables_dir / "main_wave_specific_hurdle_count_model_results.csv"
    if wave_count_path.exists():
        plot_wave_cluster_outcomes(style, pd.read_csv(wave_count_path), out_dir)

    wave_domain_demo_path = tables_dir / "main_wave_specific_domain_demographic_mixing_model_results.csv"
    if wave_domain_demo_path.exists():
        wave_domain_demo = pd.read_csv(wave_domain_demo_path)
        plot_wave_domain_demographic_mixing(style, wave_domain_demo, out_dir)

    observed_expected_path = tables_dir / "main_observed_expected_mixing_matrices.csv"
    if observed_expected_path.exists():
        plot_observed_expected_matrices(style, pd.read_csv(observed_expected_path), out_dir)

    loglinear_path = tables_dir / "main_loglinear_count_model_results.csv"
    if loglinear_path.exists():
        loglinear_results = pd.read_csv(loglinear_path)
        plot_loglinear_comparison(style, count_results, loglinear_results, out_dir)

    print(f"Wrote manuscript figures to {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing the model result CSV tables to plot. "
            "Defaults to part1/main/tables (primary results). "
            "Pass the --tables-dir used for a sensitivity run to plot those results, "
            "e.g. --tables-dir part1/main/tables_health_board."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help=(
            "Directory to write figures into. "
            "Defaults to part1/main/manuscript/figures. "
            "Set a different path for sensitivity figures to avoid overwriting "
            "primary figures, e.g. --out-dir part1/main/manuscript/figures_health_board."
        ),
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing main_cluster_table.parquet. "
            "Defaults to part1/main/cache. Must match the --cache-dir used by "
            "main_analysis.py for the same sensitivity run."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    run(
        root,
        tables_dir=args.tables_dir.resolve() if args.tables_dir else None,
        out_dir=args.out_dir.resolve() if args.out_dir else None,
        cache_dir=args.cache_dir.resolve() if args.cache_dir else None,
    )


if __name__ == "__main__":
    main()
