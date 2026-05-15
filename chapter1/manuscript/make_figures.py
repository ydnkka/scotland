"""Create manuscript figures for Chapter 1.

Main figures
------------
fig1_population_measures
    Binned distributions of non-singleton cluster-scale outcomes and
    observed-minus-expected excess-mixing predictors.
fig2_main_pooled_effects
    Forest plot of the three primary mixing predictors on cluster size,
    geographic spread, and size-adjusted geographic spread.
fig3_wave_heterogeneity
    Heatmap of wave-stratified mixing effects on cluster size and
    geographic spread.
fig4_robustness_dashboard
    Dot-range robustness dashboard across primary plus six sensitivities.

Supplementary figures
---------------------
supp_fig1_outcome_distributions
    Cluster-size and datazone distributions, all vs non-singleton.
supp_fig2_mixing_distributions
    Distributions of all six excess-mixing predictors.
supp_fig3_cluster_size_mixing_boxplots
    Boxplots and sign-category proportions for excess mixing across
    cluster-size bins.
supp_fig4_observed_expected_matrices
    Observed-minus-expected pair matrices for SIMD quintile, SIMD decile,
    and age band.
supp_fig5_size_adjustment
    Forest comparing main, linear-size-adjusted, and spline-size-adjusted
    spread models.
supp_fig6_domain_mixing_predictors
    Heatmap of SIMD-domain-specific quintile mixing on the three outcomes.
supp_fig7_profile_predictors
    Forest plot of demographic and socio-demographic profile sensitivities.
supp_fig8_model_diagnostics
    Alpha, log-likelihood, n, and tail-cap diagnostics across model fits.

Outputs are written to ``chapter1/manuscript/figures`` as PDF, PNG, and TIFF.
The script uses the shared project plotting module at ``utils/style.py``.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import tempfile

_TMP_ROOT = (
    "/private/tmp" if Path("/private/tmp").exists() else tempfile.gettempdir()
)
os.environ.setdefault("MPLCONFIGDIR", f"{_TMP_ROOT}/scotland-mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", f"{_TMP_ROOT}/scotland-xdg-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Constants (shared with the figure_table_plan.md "Conventions" block)
# ---------------------------------------------------------------------------

SIZE_BIN_EDGES = [-np.inf, 2.5, 3.5, 5.5, 10.5, 20.5, 50.5, np.inf]
SIZE_BIN_LABELS = ["2", "3", "4-5", "6-10", "11-20", "21-50", ">50"]
SIZE_BIN_TICK_LABELS = ["2", "3", "4-5", "6-10", "11-\n20", "21-\n50", ">50"]

MIXING_SPECS = [
    ("age_excess_discordance", "Age"),
    ("sex_excess_discordance", "Sex"),
    ("simd_excess_discordance", "SIMD"),
]

BASELINE_TOLERANCE_PP = 0.5
BASELINE_TOLERANCE = BASELINE_TOLERANCE_PP / 100

MIXING_CATEGORY = ["Negative", "Baseline", "Positive"]
MIXING_CATEGORY_LABELS = {
    "Negative": "Negative excess",
    "Baseline": "Baseline (+/-0.5 pp)",
    "Positive": "Positive excess",
}
MIXING_CATEGORY_COLORS = {
    "Negative": "#4e79a7",
    "Baseline": "#bdbdbd",
    "Positive": "#e15759",
}

# Three primary mixing predictors in the order shown across all model
# figures. SIMD first because it is the headline finding.
PRIMARY_MIXING_TERMS = [
    "simd_excess_mixing_z",
    "age_excess_mixing_z",
    "sex_excess_mixing_z",
]
PRIMARY_MIXING_LABELS = {
    "simd_excess_mixing_z": "SIMD excess mixing",
    "age_excess_mixing_z": "Age excess mixing",
    "sex_excess_mixing_z": "Sex excess mixing",
}

# Finite-sample standardised equivalents used by the finite-sample
# sensitivity model.
FINITE_SAMPLE_TERMS = {
    "simd_excess_mixing_z": "simd_finite_sample_mixing_z",
    "age_excess_mixing_z": "age_finite_sample_mixing_z",
    "sex_excess_mixing_z": "sex_finite_sample_mixing_z",
}

# Wave order verified against chapter1/tables/wave_stratified_results.csv
# (no XBB in chapter 1).
WAVE_ORDER = ["B.1.177", "Alpha", "Delta", "BA.1", "BA.2", "BA.4", "BA.5", "BQ.1"]

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

OUTCOME_LABELS = {
    "cluster_size": "Cluster size",
    "geographic_spread": "Geographic spread",
    "geographic_spread_size_adjusted": "Geographic spread,\nsize-adjusted",
}


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


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


def primary_mixing_palette(style) -> dict[str, str]:
    palette = style.SIMD_DOMAIN_PALETTE
    return {
        "simd_excess_mixing_z": palette["crime"],
        "age_excess_mixing_z": palette["housing"],
        "sex_excess_mixing_z": palette["overall"],
    }


def save_all(style, fig, out_base: Path, width: str, height_in: float) -> dict[str, Path]:
    return style.save_figure(
        fig,
        out_base,
        width=width,
        height_in=height_in,
        dpi=600,
        save_pdf=True,
        save_png=True,
        save_tiff=True,
    )


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------


def load_non_singleton_clusters(root: Path) -> pd.DataFrame:
    cluster_path = root / "chapter1" / "cache" / "cluster_table.parquet"
    clusters = pd.read_parquet(cluster_path)
    non_singletons = clusters.loc[clusters["cluster_size"] > 1].copy()
    if non_singletons.empty:
        raise ValueError("No non-singleton clusters available for manuscript figures.")
    return non_singletons


def load_all_clusters(root: Path) -> pd.DataFrame:
    return pd.read_parquet(root / "chapter1" / "cache" / "cluster_table.parquet")


def read_csv_if_exists(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
# Drawing primitives
# ---------------------------------------------------------------------------


def binned_percent(values: pd.Series, bins: list[float], labels: list[str]) -> pd.DataFrame:
    cats = pd.cut(values, bins=bins, labels=labels, include_lowest=True, right=True)
    pct = cats.value_counts(sort=False, normalize=True).mul(100)
    out = pct.rename("percent").reset_index()
    return out.rename(columns={out.columns[0]: "bin"})


def histogram_percent(
    values: pd.Series,
    *,
    bins: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    clean = values.dropna().to_numpy(dtype=float)
    counts, edges = np.histogram(clean, bins=bins)
    total = counts.sum()
    percent = counts / total * 100 if total else counts.astype(float)
    centres = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)
    return centres, widths, percent


def configure_log_axis(ax, xlim: tuple[float, float], ticks: list[float]) -> None:
    from matplotlib.ticker import NullFormatter, NullLocator

    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{tick:g}" for tick in ticks])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())


def add_panel_labels(style, axes) -> None:
    style.add_panel_labels(np.asarray(axes).ravel(), x=-0.16, y=1.14, size=9)


def pick_log_ticks(
    xlim: tuple[float, float],
    candidates: list[float] | None = None,
) -> list[float]:
    """Pick log-axis tick values, with denser candidates for narrow ranges."""
    if candidates is None:
        span = xlim[1] / max(xlim[0], 1e-6)
        if span <= 1.6:
            candidates = [0.9, 1.0, 1.1, 1.2, 1.3, 1.5]
        elif span <= 3.0:
            candidates = [0.7, 0.85, 1.0, 1.25, 1.5, 2.0]
        else:
            candidates = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
    return [c for c in candidates if xlim[0] <= c <= xlim[1]] or [1.0]


# ---------------------------------------------------------------------------
# FIGURE 1 — descriptive distributions
# ---------------------------------------------------------------------------


def make_figure_1(root: Path, out_dir: Path) -> dict[str, Path]:
    """Create a six-panel Figure 1 using the Part 1 supplementary style."""
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.85)

    non_singletons = load_non_singleton_clusters(root)

    grey = "#6f6f6f"
    outcome_specs = [
        (
            "cluster_size",
            "Cluster size",
            "Sequences per cluster",
            SIZE_BIN_EDGES,
            SIZE_BIN_LABELS,
        ),
        (
            "duration_days",
            "Duration",
            "Days",
            [-np.inf, 0.5, 1.5, 2.5, 5.5, 10.5, 15.5, np.inf],
            ["0", "1", "2", "3-5", "6-10", "11-15", "15+"],
        ),
        (
            "geographic_spread",
            "Geographic spread",
            "Unique datazones present",
            [-np.inf, 1.5, 2.5, 4.5, 9.5, 19.5, 49.5, np.inf],
            ["1", "2", "3-4", "5-9", "10-19", "20-49", ">50"],
        ),
    ]
    mixing_bins = np.arange(-100, 101, 10)

    fig, axes = style.new_figure(
        width="double",
        height_in=4.9,
        nrows=2,
        ncols=3,
        font_scale=0.85,
    )

    for ax, (col, title, xlabel, bins, labels) in zip(axes[0], outcome_specs):
        data = binned_percent(non_singletons[col], bins, labels)
        ax.bar(data["bin"].astype(str), data["percent"], color=grey, width=0.78)
        ax.set_title(title, pad=4)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Clusters (%)" if ax is axes[0, 0] else "")
        ax.set_ylim(0, max(30, data["percent"].max() * 1.15))
        ax.tick_params(axis="x", rotation=45, length=0)
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    for ax, (col, title) in zip(axes[1], MIXING_SPECS):
        centres, widths, percent = histogram_percent(
            non_singletons[col] * 100,
            bins=mixing_bins,
        )
        ax.bar(centres, percent, width=widths * 0.92, color=grey, align="center")
        ax.axvline(0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_title(title, pad=4)
        ax.set_xlabel("Excess mixing (pp)")
        ax.set_ylabel("Clusters (%)" if ax is axes[1, 0] else "")
        ax.set_xlim(mixing_bins[0], mixing_bins[-1])
        ax.set_xticks([-100, -50, 0, 50, 100])
        ax.set_ylim(0, max(30, percent.max() * 1.15))
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    add_panel_labels(style, axes)
    fig.subplots_adjust(left=0.085, right=0.99, top=0.91, bottom=0.13, wspace=0.32, hspace=0.52)
    return save_all(style, fig, out_dir / "fig1_population_measures", "double", 4.9)


# ---------------------------------------------------------------------------
# FIGURE 2 — main pooled forest plot
# ---------------------------------------------------------------------------


def draw_forest_panel(
    ax,
    df: pd.DataFrame,
    terms: list[str],
    labels: dict[str, str],
    colours: dict[str, str],
    *,
    title: str,
    show_ylabels: bool,
    xlim: tuple[float, float],
    xlabel: str,
) -> None:
    y_positions = np.arange(len(terms))[::-1]
    pos = dict(zip(terms, y_positions))
    for term in terms:
        row = df[df["term"] == term]
        if row.empty:
            continue
        row = row.iloc[0]
        y = pos[term]
        ax.plot(
            [row["ratio_lower"], row["ratio_upper"]],
            [y, y],
            color=colours[term],
            linewidth=1.3,
            solid_capstyle="round",
        )
        ax.scatter(
            row["ratio"],
            y,
            color=colours[term],
            edgecolor="white",
            linewidth=0.3,
            s=22,
            zorder=3,
        )

    ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
    configure_log_axis(ax, xlim, pick_log_ticks(xlim))
    ax.set_title(title, pad=4)
    ax.set_ylim(-0.7, len(terms) - 0.3)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([labels[t] for t in terms] if show_ylabels else [])
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)


def make_figure_2(root: Path, out_dir: Path) -> dict[str, Path]:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.85)

    tables = root / "chapter1" / "tables"
    main = pd.read_csv(tables / "main_effects_results.csv")
    sas = pd.read_csv(tables / "size_adjusted_spread_results.csv")

    colours = primary_mixing_palette(style)

    panels = [
        (
            "Cluster size\nZTNB count ratio",
            main[(main["outcome"] == "cluster_size") & (main["model"] == "main")],
        ),
        (
            "Geographic spread\nZTNB count ratio",
            main[(main["outcome"] == "geographic_spread") & (main["model"] == "main")],
        ),
        (
            OUTCOME_LABELS["geographic_spread_size_adjusted"] + "\nZTNB count ratio",
            sas[sas["model"] == "main_size_adjusted"],
        ),
    ]

    # Compute a panel-specific xlim using just the three primary mixing terms.
    xlims = []
    for _, sub in panels:
        sub_mix = sub[sub["term"].isin(PRIMARY_MIXING_TERMS)]
        if sub_mix.empty:
            xlims.append((0.8, 2.0))
            continue
        lo = max(0.3, float(sub_mix["ratio_lower"].min()) * 0.92)
        hi = min(8.0, float(sub_mix["ratio_upper"].max()) * 1.08)
        lo = min(lo, 0.9)
        hi = max(hi, 1.2)
        xlims.append((lo, hi))

    fig, axes = style.new_figure(
        width="double",
        height_in=3.4,
        nrows=1,
        ncols=3,
        font_scale=0.85,
    )

    for idx, ((title, sub), xlim) in enumerate(zip(panels, xlims)):
        draw_forest_panel(
            axes[idx],
            sub,
            PRIMARY_MIXING_TERMS,
            PRIMARY_MIXING_LABELS,
            colours,
            title=title,
            show_ylabels=(idx == 0),
            xlim=xlim,
            xlabel="Count ratio per 1 SD\nhigher excess mixing",
        )

    style.add_panel_labels(axes, x=-0.10, y=1.18, size=9)
    fig.subplots_adjust(left=0.205, right=0.985, top=0.85, bottom=0.22, wspace=0.20)
    return save_all(style, fig, out_dir / "fig2_main_pooled_effects", "double", 3.4)


# ---------------------------------------------------------------------------
# FIGURE 3 — wave heatmap
# ---------------------------------------------------------------------------


def make_figure_3(root: Path, out_dir: Path) -> dict[str, Path]:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.80)

    tables = root / "chapter1" / "tables"
    wave = pd.read_csv(tables / "wave_stratified_results.csv")
    wave = wave[wave["term"].isin(PRIMARY_MIXING_TERMS)].copy()
    waves = [w for w in WAVE_ORDER if w in set(wave["wave"].astype(str))]
    if wave.empty or not waves:
        return {}

    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.pyplot as plt

    panels = [("cluster_size", "Cluster size"), ("geographic_spread", "Geographic spread")]
    values = wave[wave["outcome"].isin([p[0] for p in panels])]
    vmax = max(2.0, float(np.nanmax(values["ratio"])))
    vmin = min(0.5, float(np.nanmin(values["ratio"])))
    if vmax <= 1.0:
        vmax = 1.5
    if vmin >= 1.0:
        vmin = 0.5
    norm = TwoSlopeNorm(vcenter=1.0, vmin=vmin, vmax=vmax)

    fig, axes = style.new_figure(
        width="double",
        height_in=3.6,
        nrows=1,
        ncols=2,
        font_scale=0.80,
    )

    image = None
    for ax_idx, (outcome, title) in enumerate(panels):
        ax = axes[ax_idx]
        sub = wave[wave["outcome"] == outcome]
        matrix = (
            sub.pivot_table(
                index="term",
                columns="wave",
                values="ratio",
                aggfunc="first",
            )
            .reindex(index=PRIMARY_MIXING_TERMS, columns=waves)
        )
        ci_lower = (
            sub.pivot_table(
                index="term",
                columns="wave",
                values="ratio_lower",
                aggfunc="first",
            )
            .reindex(index=PRIMARY_MIXING_TERMS, columns=waves)
        )
        ci_upper = (
            sub.pivot_table(
                index="term",
                columns="wave",
                values="ratio_upper",
                aggfunc="first",
            )
            .reindex(index=PRIMARY_MIXING_TERMS, columns=waves)
        )

        arr = matrix.to_numpy(dtype=float)
        image = ax.imshow(arr, cmap="RdBu_r", norm=norm, aspect="auto")
        ax.set_title(title, pad=6)
        ax.set_yticks(np.arange(len(PRIMARY_MIXING_TERMS)))
        ax.set_yticklabels(
            [PRIMARY_MIXING_LABELS[t] for t in PRIMARY_MIXING_TERMS]
            if ax_idx == 0
            else []
        )
        ax.set_xticks(np.arange(len(waves)))
        ax.set_xticklabels(waves, rotation=35, ha="right")
        ax.tick_params(length=0)
        for y in np.arange(len(PRIMARY_MIXING_TERMS) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.6)
        for x in np.arange(len(waves) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.6)

        # Annotate each cell with the ratio; mark non-significant cells.
        for i, term in enumerate(PRIMARY_MIXING_TERMS):
            for j, w in enumerate(waves):
                val = arr[i, j]
                if np.isnan(val):
                    continue
                lower = ci_lower.loc[term, w]
                upper = ci_upper.loc[term, w]
                text_colour = (
                    "white" if abs(np.log(val)) > 0.30 else "#222222"
                )
                ax.text(
                    j,
                    i,
                    f"{val:.2f}",
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color=text_colour,
                )
                if not (np.isnan(lower) or np.isnan(upper)) and lower <= 1.0 <= upper:
                    # Open white circle marker for non-significant cells.
                    ax.scatter(
                        j,
                        i + 0.32,
                        s=20,
                        facecolors="none",
                        edgecolors="white",
                        linewidths=0.9,
                    )

    assert image is not None
    fig.subplots_adjust(left=0.14, right=0.86, top=0.88, bottom=0.22, wspace=0.06)
    cbar_ax = fig.add_axes([0.885, 0.24, 0.020, 0.58])
    cbar = fig.colorbar(image, cax=cbar_ax, extend="both")
    cbar.set_label("ZTNB count ratio per 1 SD higher excess mixing", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)
    style.add_panel_labels(axes, x=-0.08, y=1.10, size=9)
    plt.close("all") if False else None
    return save_all(style, fig, out_dir / "fig3_wave_heterogeneity", "double", 3.6)


# ---------------------------------------------------------------------------
# FIGURE 4 — robustness dashboard
# ---------------------------------------------------------------------------


def collect_robustness_rows(root: Path) -> pd.DataFrame:
    tables = root / "chapter1" / "tables"
    sens = root / "chapter1" / "sensitivity"

    def standardise(df: pd.DataFrame, sensitivity: str, outcome_override: str | None = None) -> pd.DataFrame:
        keep = df[["outcome", "term", "ratio", "ratio_lower", "ratio_upper", "p_value"]].copy()
        keep["sensitivity"] = sensitivity
        if outcome_override is not None:
            keep["outcome"] = outcome_override
        return keep

    rows: list[pd.DataFrame] = []

    main = pd.read_csv(tables / "main_effects_results.csv")
    main = main[main["model"] == "main"]
    rows.append(standardise(main, "Primary (window SE)"))

    sas = pd.read_csv(tables / "size_adjusted_spread_results.csv")
    sas = sas[sas["model"] == "main_size_adjusted"]
    rows.append(standardise(sas, "Primary (window SE)", outcome_override="geographic_spread_size_adjusted"))

    for label, subdir in [
        ("Health-board SE", "tables_health_board"),
        ("Non-overlapping windows (stride 3)", "tables_stride3"),
        ("99% winsorised", "tables_winsorise99"),
        ("Top 0.5% excluded", "tables_exclude_tail995"),
    ]:
        m_path = sens / subdir / "main_effects_results.csv"
        s_path = sens / subdir / "size_adjusted_spread_results.csv"
        if m_path.exists():
            df = pd.read_csv(m_path)
            df = df[df["model"] == "main"]
            rows.append(standardise(df, label))
        if s_path.exists():
            df = pd.read_csv(s_path)
            df = df[df["model"] == "main_size_adjusted"]
            rows.append(standardise(df, label, outcome_override="geographic_spread_size_adjusted"))

    fs_path = tables / "finite_sample_mixing_sensitivity_results.csv"
    if fs_path.exists():
        fs = pd.read_csv(fs_path)
        # Map finite-sample terms back to their primary equivalents.
        fs_map = {v: k for k, v in FINITE_SAMPLE_TERMS.items()}
        fs = fs[fs["term"].isin(fs_map)].copy()
        fs["term"] = fs["term"].map(fs_map)
        rows.append(standardise(fs, "Finite-sample mixing"))

    jp_path = tables / "joint_profile_adjusted_results.csv"
    if jp_path.exists():
        jp = pd.read_csv(jp_path)
        rows.append(standardise(jp, "Joint-profile adjusted"))

    out = pd.concat(rows, ignore_index=True)
    out = out[out["term"].isin(PRIMARY_MIXING_TERMS)].copy()
    return out


def make_figure_4(root: Path, out_dir: Path) -> dict[str, Path]:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.80)

    rows = collect_robustness_rows(root)
    if rows.empty:
        return {}

    sens_order = [
        "Primary (window SE)",
        "Health-board SE",
        "Non-overlapping windows (stride 3)",
        "99% winsorised",
        "Top 0.5% excluded",
        "Finite-sample mixing",
        "Joint-profile adjusted",
    ]
    sens_order = [s for s in sens_order if s in set(rows["sensitivity"])]

    outcomes = [
        ("cluster_size", "Cluster size"),
        ("geographic_spread", "Geographic spread"),
        ("geographic_spread_size_adjusted", "Spread (size-adj)"),
    ]

    colours = primary_mixing_palette(style)
    offsets = {
        "simd_excess_mixing_z": +0.24,
        "age_excess_mixing_z": 0.0,
        "sex_excess_mixing_z": -0.24,
    }

    n_rows_per_sens = 1.0
    y_centres = {sens: idx * n_rows_per_sens for idx, sens in enumerate(sens_order[::-1])}

    # Calculate a shared xlim across all outcome columns separately to keep
    # comparisons fair within each outcome.
    xlim_by_outcome: dict[str, tuple[float, float]] = {}
    for outcome, _ in outcomes:
        sub = rows[rows["outcome"] == outcome]
        if sub.empty:
            xlim_by_outcome[outcome] = (0.8, 1.2)
            continue
        lo = max(0.3, float(sub["ratio_lower"].min()) * 0.92)
        hi = min(6.0, float(sub["ratio_upper"].max()) * 1.08)
        lo = min(lo, 0.9)
        hi = max(hi, 1.2)
        xlim_by_outcome[outcome] = (lo, hi)

    import matplotlib.pyplot as plt
    fig, axes = style.new_figure(
        width="double",
        height_in=5.6,
        nrows=1,
        ncols=3,
        font_scale=0.80,
    )

    for col_idx, (outcome, title) in enumerate(outcomes):
        ax = axes[col_idx]
        sub = rows[rows["outcome"] == outcome]
        for sens in sens_order:
            y0 = y_centres[sens]
            sens_sub = sub[sub["sensitivity"] == sens]
            for term in PRIMARY_MIXING_TERMS:
                row = sens_sub[sens_sub["term"] == term]
                if row.empty:
                    continue
                row = row.iloc[0]
                y = y0 + offsets[term]
                ax.plot(
                    [row["ratio_lower"], row["ratio_upper"]],
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
        # Sensitivity dividers between blocks.
        for idx in range(len(sens_order) - 1):
            y_mid = (idx * n_rows_per_sens + (idx + 1) * n_rows_per_sens) / 2.0
            ax.axhline(y_mid, color="#dddddd", linewidth=0.5)
        ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
        xlim = xlim_by_outcome[outcome]
        configure_log_axis(ax, xlim, pick_log_ticks(xlim))
        ax.set_title(title, pad=4)
        y_ticks = [y_centres[s] for s in sens_order]
        ax.set_yticks(y_ticks)
        if col_idx == 0:
            ax.set_yticklabels(sens_order)
        else:
            ax.set_yticklabels([])
        ax.set_ylim(min(y_ticks) - 0.6, max(y_ticks) + 0.6)
        ax.set_xlabel("ZTNB count ratio per 1 SD")
        ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    # Legend across the three mixing predictors.
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="-",
            color=colours[t],
            markeredgecolor="white",
            markeredgewidth=0.3,
            label=PRIMARY_MIXING_LABELS[t],
        )
        for t in PRIMARY_MIXING_TERMS
    ]
    fig.legend(
        handles=handles,
        labels=[PRIMARY_MIXING_LABELS[t] for t in PRIMARY_MIXING_TERMS],
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.55, -0.005),
    )

    style.add_panel_labels(axes, x=-0.08, y=1.08, size=9)
    fig.subplots_adjust(left=0.27, right=0.99, top=0.92, bottom=0.13, wspace=0.10)
    return save_all(style, fig, out_dir / "fig4_robustness_dashboard", "double", 5.6)


# ---------------------------------------------------------------------------
# SUPP FIG 1 — Full outcome distributions
# ---------------------------------------------------------------------------


def make_supp_fig_1(root: Path, out_dir: Path) -> dict[str, Path]:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.80)

    clusters = load_all_clusters(root)
    non_singletons = clusters.loc[clusters["cluster_size"] > 1]

    grey = "#6f6f6f"
    accent = "#4e79a7"

    import matplotlib.pyplot as plt

    fig, axes = style.new_figure(
        width="double",
        height_in=5.0,
        nrows=2,
        ncols=3,
        font_scale=0.80,
    )

    # Pre-compute log-spaced bins for cluster sizes and datazones.
    size_bins = np.logspace(0, np.log10(max(clusters["cluster_size"].max(), 2.0)) + 0.05, 30)
    dz_bins = np.logspace(0, np.log10(max(clusters["cluster_n_datazones"].max(), 2.0)) + 0.05, 30)

    def log_hist(ax, values, bins, *, mark_first: bool, xlabel: str, title: str):
        counts, edges = np.histogram(values.dropna(), bins=bins)
        centres = np.sqrt(edges[:-1] * edges[1:])
        widths = np.diff(edges)
        colours = [accent if mark_first and i == 0 else grey for i in range(len(centres))]
        ax.bar(centres, counts, width=widths * 0.95, color=colours, align="center")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_title(title, pad=4)
        ax.grid(axis="y", color="#dddddd", linewidth=0.5, which="both")

    # Panel A — all clusters cluster_size.
    log_hist(
        axes[0, 0],
        clusters["cluster_size"],
        size_bins,
        mark_first=True,
        xlabel="Cluster size",
        title="All clusters",
    )
    axes[0, 0].set_ylabel("Clusters (count, log)")

    # Panel B — non-singleton cluster_size.
    log_hist(
        axes[0, 1],
        non_singletons["cluster_size"],
        size_bins,
        mark_first=False,
        xlabel="Cluster size",
        title="Non-singletons",
    )

    # Panel C — all clusters datazones.
    log_hist(
        axes[0, 2],
        clusters["cluster_n_datazones"],
        dz_bins,
        mark_first=True,
        xlabel="Distinct datazones",
        title="All clusters",
    )

    # Panel D — non-singleton datazones.
    log_hist(
        axes[1, 0],
        non_singletons["cluster_n_datazones"],
        dz_bins,
        mark_first=False,
        xlabel="Distinct datazones",
        title="Non-singletons",
    )
    axes[1, 0].set_ylabel("Clusters (count, log)")

    # Panel E — hexbin of log-size vs log-datazones.
    ax = axes[1, 1]
    log_size = np.log10(non_singletons["cluster_size"].astype(float))
    log_dz = np.log10(non_singletons["cluster_n_datazones"].clip(lower=1).astype(float))
    hb = ax.hexbin(
        log_size,
        log_dz,
        gridsize=24,
        bins="log",
        cmap="Greys",
        mincnt=1,
    )
    diag_max = float(max(log_size.max(), log_dz.max()))
    ax.plot([0, diag_max], [0, diag_max], color="#e15759", linewidth=0.8, linestyle="--")
    ax.set_xlabel("log10 cluster size")
    ax.set_ylabel("log10 datazones")
    ax.set_title("Size vs datazones (non-singletons)", pad=4)
    cb = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("log10 count", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)

    # Panel F — annotation, used to display key descriptive stats.
    ax = axes[1, 2]
    ax.axis("off")
    desc = clusters["cluster_size"].describe(percentiles=[0.5, 0.9, 0.99]).to_dict()
    dz_desc = clusters["cluster_n_datazones"].describe(percentiles=[0.5, 0.9, 0.99]).to_dict()
    lines = [
        f"All clusters: n = {len(clusters):,}",
        f"Singletons: {(clusters['cluster_size'] == 1).sum():,}",
        f"Non-singletons: {(clusters['cluster_size'] > 1).sum():,}",
        "",
        f"Cluster size median (IQR): {desc['50%']:.0f} "
        f"(25 to 75 pct: {clusters['cluster_size'].quantile(0.25):.0f}-"
        f"{clusters['cluster_size'].quantile(0.75):.0f})",
        f"Cluster size 90th pct: {desc['90%']:.0f}",
        f"Cluster size 99th pct: {desc['99%']:.0f}",
        f"Cluster size max: {desc['max']:.0f}",
        "",
        f"Datazones median: {dz_desc['50%']:.0f}",
        f"Datazones 90th pct: {dz_desc['90%']:.0f}",
        f"Datazones 99th pct: {dz_desc['99%']:.0f}",
    ]
    ax.text(0, 0.95, "\n".join(lines), va="top", ha="left", fontsize=7, family="monospace")

    style.add_panel_labels(axes.ravel(), x=-0.14, y=1.12, size=9)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.10, wspace=0.40, hspace=0.45)
    return save_all(style, fig, out_dir / "supp_fig1_outcome_distributions", "double", 5.0)


# ---------------------------------------------------------------------------
# SUPP FIG 2 — Mixing predictor distributions
# ---------------------------------------------------------------------------


def make_supp_fig_2(root: Path, out_dir: Path) -> dict[str, Path]:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.80)

    non_singletons = load_non_singleton_clusters(root)
    grey = "#6f6f6f"

    panels = [
        ("age_excess_discordance", "age_excess_mixing_z", "Age"),
        ("sex_excess_discordance", "sex_excess_mixing_z", "Sex"),
        ("simd_excess_discordance", "simd_excess_mixing_z", "SIMD quintile"),
        ("simd_decile_excess_discordance", "simd_decile_excess_mixing_z", "SIMD decile"),
        ("demographic_profile_excess_discordance", "demographic_profile_excess_mixing_z", "Demographic"),
        (
            "socio_demographic_profile_excess_discordance",
            "socio_demographic_profile_excess_mixing_z",
            "Socio-demographic",
        ),
    ]

    bins = np.arange(-100, 101, 10)

    fig, axes = style.new_figure(
        width="double",
        height_in=4.6,
        nrows=2,
        ncols=3,
        font_scale=0.80,
    )

    for ax, (raw_col, z_col, title) in zip(axes.ravel(), panels):
        if raw_col not in non_singletons.columns:
            ax.axis("off")
            continue
        centres, widths, percent = histogram_percent(non_singletons[raw_col] * 100, bins=bins)
        ax.bar(centres, percent, width=widths * 0.92, color=grey, align="center")
        ax.axvline(0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_title(title, pad=4)
        ax.set_xlabel("Excess mixing (pp)")
        ax.set_xlim(-100, 100)
        ax.set_xticks([-100, -50, 0, 50, 100])
        ax.set_ylim(0, max(30, percent.max() * 1.15))
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    axes[0, 0].set_ylabel("Clusters (%)")
    axes[1, 0].set_ylabel("Clusters (%)")

    style.add_panel_labels(axes.ravel(), x=-0.18, y=1.18, size=9)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.92, bottom=0.10, wspace=0.30, hspace=0.55)
    return save_all(style, fig, out_dir / "supp_fig2_mixing_distributions", "double", 4.6)


# ---------------------------------------------------------------------------
# SUPP FIG 3 — Boxplots and stacked bars across cluster-size bins (existing)
# ---------------------------------------------------------------------------


def make_cluster_size_mixing_boxplots(root: Path, out_dir: Path) -> dict[str, Path]:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.85)

    non_singletons = load_non_singleton_clusters(root)
    non_singletons["cluster_size_bin"] = pd.cut(
        non_singletons["cluster_size"],
        bins=SIZE_BIN_EDGES,
        labels=SIZE_BIN_LABELS,
        include_lowest=True,
        right=True,
    )

    fig, axes = style.new_figure(
        width="double",
        height_in=4.8,
        nrows=2,
        ncols=3,
        font_scale=0.85,
        gridspec_kw={"height_ratios": [2.0, 1.1]},
    )
    top_axes = axes[0]
    bottom_axes = axes[1]

    boxprops = dict(facecolor="#d9d9d9", edgecolor="#4f4f4f", linewidth=0.8)
    medianprops = dict(color="#222222", linewidth=1.1)
    whiskerprops = dict(color="#4f4f4f", linewidth=0.8)
    capprops = dict(color="#4f4f4f", linewidth=0.8)

    for ax, (col, title) in zip(top_axes, MIXING_SPECS):
        grouped = [
            non_singletons.loc[non_singletons["cluster_size_bin"] == label, col]
            .dropna()
            .mul(100)
            .to_numpy()
            for label in SIZE_BIN_LABELS
        ]
        ax.boxplot(
            grouped,
            patch_artist=True,
            widths=0.62,
            showfliers=False,
            boxprops=boxprops,
            medianprops=medianprops,
            whiskerprops=whiskerprops,
            capprops=capprops,
        )
        ax.axhline(0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_title(title, pad=4)
        ax.set_xlabel("")
        ax.set_ylabel("Excess mixing (pp)" if ax is top_axes[0] else "")
        ax.set_xticks(range(1, len(SIZE_BIN_LABELS) + 1))
        ax.set_xticklabels([])
        ax.tick_params(axis="x", length=0)
        ax.set_ylim(-100, 100)
        ax.set_yticks([-100, -50, 0, 50, 100])
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    x = np.arange(len(SIZE_BIN_LABELS))
    for ax, (col, _title) in zip(bottom_axes, MIXING_SPECS):
        category = pd.Series(
            np.select(
                [
                    non_singletons[col] < -BASELINE_TOLERANCE,
                    non_singletons[col].abs() <= BASELINE_TOLERANCE,
                    non_singletons[col] > BASELINE_TOLERANCE,
                ],
                MIXING_CATEGORY,
                default="Missing",
            ),
            index=non_singletons.index,
        )
        props = (
            pd.crosstab(non_singletons["cluster_size_bin"], category, normalize="index")
            .reindex(index=SIZE_BIN_LABELS, columns=MIXING_CATEGORY, fill_value=0)
            .mul(100)
        )

        bottom = np.zeros(len(SIZE_BIN_LABELS))
        for label in MIXING_CATEGORY:
            heights = props[label].to_numpy(dtype=float)
            ax.bar(
                x,
                heights,
                bottom=bottom,
                color=MIXING_CATEGORY_COLORS[label],
                edgecolor="white",
                linewidth=0.4,
                width=0.78,
                label=MIXING_CATEGORY_LABELS[label],
            )
            bottom += heights

        ax.set_xlabel("Cluster size bin")
        ax.set_ylabel("Clusters (%)" if ax is bottom_axes[0] else "")
        ax.set_xticks(x)
        ax.set_xticklabels(SIZE_BIN_TICK_LABELS)
        ax.tick_params(axis="x", length=0)
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 50, 100])
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    for ax, label in zip(top_axes, list("ABC")):
        ax.text(
            -0.16,
            1.15,
            label,
            transform=ax.transAxes,
            fontsize=9,
            fontweight="bold",
            va="top",
        )
    handles, labels = bottom_axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.54, 0.005),
        handlelength=1.4,
        columnspacing=1.4,
    )
    fig.subplots_adjust(
        left=0.105,
        right=0.99,
        top=0.9,
        bottom=0.19,
        wspace=0.32,
        hspace=0.18,
    )

    return save_all(style, fig, out_dir / "supp_fig3_cluster_size_mixing_boxplots", "double", 4.8)


# ---------------------------------------------------------------------------
# SUPP FIG 4 — Observed-expected pair matrices
# ---------------------------------------------------------------------------


def make_supp_fig_4(root: Path, out_dir: Path) -> dict[str, Path] | None:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.80)

    path = root / "chapter1" / "tables" / "observed_expected_mixing_matrices.csv"
    if not path.exists():
        return None
    matrices = pd.read_csv(path)
    overall = matrices[matrices["wave_group"] == "Overall"].copy()
    if overall.empty:
        return None

    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.pyplot as plt

    def category_key(value: object) -> tuple[int, str]:
        text = str(value)
        if text == "75+":
            return (75, text)
        digits = "".join(ch if ch.isdigit() else " " for ch in text).split()
        return (int(digits[0]) if digits else 999, text)

    specs = [
        ("simd_quintile", "SIMD quintile"),
        ("simd_decile", "SIMD decile"),
        ("age", "Age band"),
    ]
    specs = [s for s in specs if s[0] in set(overall["variable"])]

    fig, axes = style.new_figure(
        width="double",
        height_in=3.8,
        nrows=1,
        ncols=len(specs),
        font_scale=0.80,
    )
    if len(specs) == 1:
        axes = np.array([axes])

    vmax = max(0.5, float(np.nanmax(np.abs(overall["excess_percentage_points"]))))
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)

    image = None
    for idx, (variable, title) in enumerate(specs):
        ax = axes[idx]
        sub = overall[overall["variable"] == variable]
        row_order = sorted(sub["category_i"].astype(str).unique(), key=category_key)
        col_order = sorted(sub["category_j"].astype(str).unique(), key=category_key)
        matrix = (
            sub.assign(
                category_i=sub["category_i"].astype(str),
                category_j=sub["category_j"].astype(str),
            )
            .pivot_table(
                index="category_i",
                columns="category_j",
                values="excess_percentage_points",
                aggfunc="first",
            )
            .reindex(index=row_order, columns=col_order)
        )
        image = ax.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", norm=norm, aspect="equal")
        ax.set_title(title, pad=4)
        ax.set_xticks(np.arange(len(col_order)))
        ax.set_yticks(np.arange(len(row_order)))
        ax.set_xticklabels(col_order, rotation=45 if "age" in variable or "decile" in variable else 0, ha="right")
        ax.set_yticklabels(row_order)
        ax.tick_params(length=0)
        for y in np.arange(len(row_order) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.4)
        for x in np.arange(len(col_order) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.4)
        ax.set_xlabel(title)
        if idx == 0:
            ax.set_ylabel(title)

    assert image is not None
    fig.subplots_adjust(left=0.08, right=0.86, top=0.86, bottom=0.20, wspace=0.30)
    cbar_ax = fig.add_axes([0.885, 0.24, 0.020, 0.58])
    cbar = fig.colorbar(image, cax=cbar_ax, extend="both")
    cbar.set_label("Observed - expected pair probability (pp)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)
    style.add_panel_labels(axes, x=-0.10, y=1.10, size=9)
    return save_all(style, fig, out_dir / "supp_fig4_observed_expected_matrices", "double", 3.8)


# ---------------------------------------------------------------------------
# SUPP FIG 5 — Size adjustment forest
# ---------------------------------------------------------------------------


def make_supp_fig_5(root: Path, out_dir: Path) -> dict[str, Path] | None:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.85)

    tables = root / "chapter1" / "tables"
    main = pd.read_csv(tables / "main_effects_results.csv")
    sas = pd.read_csv(tables / "size_adjusted_spread_results.csv")
    spline_path = tables / "size_spline_sensitivity_results.csv"
    spline = pd.read_csv(spline_path) if spline_path.exists() else None

    def select(df: pd.DataFrame, model: str | None) -> pd.DataFrame:
        sub = df[df["term"].isin(PRIMARY_MIXING_TERMS)]
        if model is not None:
            sub = sub[sub["model"] == model]
        return sub

    pieces: list[tuple[str, pd.DataFrame]] = [
        ("Main spread", select(main[main["outcome"] == "geographic_spread"], "main")),
        ("Linear log-size adjusted", select(sas, "main_size_adjusted")),
    ]
    if spline is not None and not spline.empty:
        pieces.append(("Spline log-size adjusted", select(spline, "size_spline")))

    colours = primary_mixing_palette(style)
    model_offsets = {label: (i - (len(pieces) - 1) / 2) * 0.22 for i, (label, _) in enumerate(pieces)}
    markers = {0: "o", 1: "s", 2: "^"}

    import matplotlib.pyplot as plt
    fig, ax = style.new_figure(width="onehalf", height_in=3.4, font_scale=0.85)

    y_positions = np.arange(len(PRIMARY_MIXING_TERMS))[::-1]
    pos = dict(zip(PRIMARY_MIXING_TERMS, y_positions))
    all_lower = []
    all_upper = []
    for idx, (label, sub) in enumerate(pieces):
        for term in PRIMARY_MIXING_TERMS:
            row = sub[sub["term"] == term]
            if row.empty:
                continue
            row = row.iloc[0]
            y = pos[term] + model_offsets[label]
            all_lower.append(float(row["ratio_lower"]))
            all_upper.append(float(row["ratio_upper"]))
            ax.plot(
                [row["ratio_lower"], row["ratio_upper"]],
                [y, y],
                color=colours[term],
                linewidth=1.0,
                solid_capstyle="round",
            )
            ax.scatter(
                row["ratio"],
                y,
                color=colours[term],
                edgecolor="white",
                linewidth=0.3,
                marker=markers[idx],
                s=22,
                zorder=3,
                label=label if term == PRIMARY_MIXING_TERMS[0] else None,
            )

    ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
    lo = max(0.3, min(all_lower) * 0.92) if all_lower else 0.5
    hi = min(8.0, max(all_upper) * 1.08) if all_upper else 3.0
    lo = min(lo, 0.8)
    hi = max(hi, 1.5)
    configure_log_axis(ax, (lo, hi), pick_log_ticks((lo, hi)))
    ax.set_yticks(y_positions)
    ax.set_yticklabels([PRIMARY_MIXING_LABELS[t] for t in PRIMARY_MIXING_TERMS])
    ax.set_ylim(-0.7, len(PRIMARY_MIXING_TERMS) - 0.3)
    ax.set_xlabel("ZTNB count ratio per 1 SD higher excess mixing")
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        seen = {}
        for h, l in zip(handles, labels):
            seen.setdefault(l, h)
        fig.legend(
            seen.values(),
            seen.keys(),
            loc="lower center",
            ncol=len(pieces),
            frameon=False,
            bbox_to_anchor=(0.55, -0.02),
        )

    fig.subplots_adjust(left=0.30, right=0.98, top=0.94, bottom=0.20)
    return save_all(style, fig, out_dir / "supp_fig5_size_adjustment", "onehalf", 3.4)


# ---------------------------------------------------------------------------
# SUPP FIG 6 — Domain mixing predictors heatmap
# ---------------------------------------------------------------------------


def make_supp_fig_6(root: Path, out_dir: Path) -> dict[str, Path] | None:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.80)

    path = root / "chapter1" / "tables" / "domain_main_effects_results.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)

    # Identify each domain's own domain-quintile predictor row.
    df["is_domain_term"] = df.apply(
        lambda r: str(r["term"]) == f"{r['domain']}_domain_excess_mixing_z", axis=1
    )
    domain_quintile = df[df["is_domain_term"]].copy()
    if domain_quintile.empty:
        return None

    # Chapter 1 domain results currently include only cluster_size and
    # geographic_spread; size-adjusted spread is supplied via the
    # size_adjusted_spread table if domain-stratified runs exist.
    outcomes = [
        ("cluster_size", "Cluster size"),
        ("geographic_spread", "Geographic spread"),
    ]

    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.pyplot as plt

    fig, axes = style.new_figure(
        width="double",
        height_in=3.8,
        nrows=1,
        ncols=2,
        font_scale=0.80,
    )

    values = domain_quintile[domain_quintile["outcome"].isin([o for o, _ in outcomes])]
    vmax = max(2.0, float(np.nanmax(values["ratio"])))
    vmin = min(0.5, float(np.nanmin(values["ratio"])))
    if vmax <= 1.0:
        vmax = 1.5
    if vmin >= 1.0:
        vmin = 0.5
    norm = TwoSlopeNorm(vcenter=1.0, vmin=vmin, vmax=vmax)

    image = None
    for idx, (outcome, title) in enumerate(outcomes):
        ax = axes[idx]
        sub = domain_quintile[domain_quintile["outcome"] == outcome]
        matrix = (
            sub.pivot_table(index="domain", columns="outcome", values="ratio", aggfunc="first")
            .reindex(DOMAIN_ORDER)
        )
        arr = matrix.to_numpy(dtype=float).reshape(-1, 1)
        image = ax.imshow(arr, cmap="RdBu_r", norm=norm, aspect="auto")
        ax.set_title(title, pad=10)
        ax.set_yticks(np.arange(len(DOMAIN_ORDER)))
        ax.set_yticklabels([DOMAIN_LABELS[d] for d in DOMAIN_ORDER] if idx == 0 else [])
        ax.set_xticks([0])
        ax.set_xticklabels(["Domain quintile"], rotation=0, ha="center")
        ax.tick_params(length=0)
        for y in np.arange(len(DOMAIN_ORDER) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.6)
        ax.axvline(0.5, color="white", linewidth=0.6)
        ax.axvline(-0.5, color="white", linewidth=0.6)
        for i in range(len(DOMAIN_ORDER)):
            v = arr[i, 0]
            if np.isnan(v):
                continue
            colour = "white" if abs(np.log(v)) > 0.3 else "#222222"
            ax.text(0, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color=colour)

    assert image is not None
    fig.subplots_adjust(left=0.20, right=0.84, top=0.82, bottom=0.15, wspace=0.30)
    cbar_ax = fig.add_axes([0.865, 0.20, 0.025, 0.62])
    cbar = fig.colorbar(image, cax=cbar_ax, extend="both")
    cbar.set_label("ZTNB count ratio per 1 SD higher domain excess mixing", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)
    style.add_panel_labels(axes, x=-0.20, y=1.16, size=9)
    return save_all(style, fig, out_dir / "supp_fig6_domain_mixing_predictors", "double", 3.8)


# ---------------------------------------------------------------------------
# SUPP FIG 7 — Profile predictor sensitivities
# ---------------------------------------------------------------------------


def make_supp_fig_7(root: Path, out_dir: Path) -> dict[str, Path] | None:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.80)

    tables = root / "chapter1" / "tables"
    profile_path = tables / "profile_predictor_results.csv"
    joint_path = tables / "joint_profile_adjusted_results.csv"
    if not (profile_path.exists() and joint_path.exists()):
        return None
    profile = pd.read_csv(profile_path)
    joint = pd.read_csv(joint_path)

    demo_terms = ["demographic_profile_excess_mixing_z"]
    socio_terms = ["socio_demographic_profile_excess_mixing_z"]
    joint_terms = PRIMARY_MIXING_TERMS + ["socio_demographic_profile_excess_mixing_z"]
    label_map = {
        "demographic_profile_excess_mixing_z": "Demographic profile",
        "socio_demographic_profile_excess_mixing_z": "Socio-demographic profile",
        **PRIMARY_MIXING_LABELS,
    }
    colours = primary_mixing_palette(style)
    palette = style.SIMD_DOMAIN_PALETTE
    colours["demographic_profile_excess_mixing_z"] = palette["income"]
    colours["socio_demographic_profile_excess_mixing_z"] = palette["employment"]

    panels = [
        ("Demographic profile alone", profile[profile["profile"] == "demographic_profile"], demo_terms),
        ("Socio-demographic profile alone", profile[profile["profile"] == "socio_demographic_profile"], socio_terms),
        ("Age/Sex/SIMD + Socio profile", joint, joint_terms),
    ]

    import matplotlib.pyplot as plt
    fig, axes = style.new_figure(
        width="double",
        height_in=3.8,
        nrows=1,
        ncols=3,
        font_scale=0.80,
    )

    # Two outcomes only, plotted as offsets within the same panel row.
    outcome_offsets = {"cluster_size": +0.15, "geographic_spread": -0.15}
    outcome_markers = {"cluster_size": "o", "geographic_spread": "^"}

    for idx, (title, sub, terms) in enumerate(panels):
        ax = axes[idx]
        y_positions = np.arange(len(terms))[::-1]
        pos = dict(zip(terms, y_positions))
        all_lower, all_upper = [], []
        for outcome, marker in outcome_markers.items():
            outcome_sub = sub[sub["outcome"] == outcome]
            for term in terms:
                row = outcome_sub[outcome_sub["term"] == term]
                if row.empty:
                    continue
                row = row.iloc[0]
                y = pos[term] + outcome_offsets[outcome]
                all_lower.append(float(row["ratio_lower"]))
                all_upper.append(float(row["ratio_upper"]))
                ax.plot(
                    [row["ratio_lower"], row["ratio_upper"]],
                    [y, y],
                    color=colours.get(term, "#444444"),
                    linewidth=1.0,
                    solid_capstyle="round",
                )
                ax.scatter(
                    row["ratio"],
                    y,
                    color=colours.get(term, "#444444"),
                    edgecolor="white",
                    linewidth=0.3,
                    marker=marker,
                    s=22,
                    zorder=3,
                    label=outcome if (term == terms[0] and idx == 0) else None,
                )
        ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
        lo = max(0.3, min(all_lower) * 0.92) if all_lower else 0.5
        hi = min(6.0, max(all_upper) * 1.08) if all_upper else 2.0
        lo = min(lo, 0.8)
        hi = max(hi, 1.5)
        configure_log_axis(ax, (lo, hi), pick_log_ticks((lo, hi)))
        ax.set_yticks(y_positions)
        ax.set_yticklabels([label_map[t] for t in terms])
        ax.set_ylim(-0.7, len(terms) - 0.3)
        ax.set_xlabel("ZTNB count ratio per 1 SD")
        ax.set_title(title, pad=4)
        ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        seen = {}
        for h, l in zip(handles, labels):
            seen.setdefault(l, h)
        fig.legend(
            seen.values(),
            [{"cluster_size": "Cluster size", "geographic_spread": "Geographic spread"}.get(l, l) for l in seen],
            loc="lower center",
            ncol=2,
            frameon=False,
            bbox_to_anchor=(0.55, -0.02),
        )

    style.add_panel_labels(axes, x=-0.18, y=1.10, size=9)
    fig.subplots_adjust(left=0.25, right=0.99, top=0.90, bottom=0.20, wspace=0.65)
    return save_all(style, fig, out_dir / "supp_fig7_profile_predictors", "double", 3.8)


# ---------------------------------------------------------------------------
# SUPP FIG 8 — Model diagnostics
# ---------------------------------------------------------------------------


def gather_diagnostics(root: Path) -> pd.DataFrame:
    tables = root / "chapter1" / "tables"
    sens = root / "chapter1" / "sensitivity"
    pieces: list[pd.DataFrame] = []
    for path in tables.glob("*_diagnostics.csv"):
        df = pd.read_csv(path)
        df["source"] = path.stem
        df["bucket"] = "primary"
        pieces.append(df)
    for d in sens.iterdir() if sens.exists() else []:
        if not d.is_dir() or not d.name.startswith("tables_"):
            continue
        for path in d.glob("*_diagnostics.csv"):
            df = pd.read_csv(path)
            df["source"] = f"{d.name}/{path.stem}"
            df["bucket"] = d.name.replace("tables_", "")
            pieces.append(df)
    if not pieces:
        return pd.DataFrame()
    return pd.concat(pieces, ignore_index=True)


SOURCE_ABBREV = {
    "main_effects_diagnostics": "main",
    "size_adjusted_spread_diagnostics": "size-adjusted",
    "wave_stratified_diagnostics": "wave-strat",
    "wave_interaction_diagnostics": "wave-int",
    "domain_main_effects_diagnostics": "domain",
    "simd_decile_sensitivity_diagnostics": "simd-decile",
    "finite_sample_mixing_sensitivity_diagnostics": "finite-sample",
    "joint_profile_adjusted_diagnostics": "joint-profile",
    "null_residual_sensitivity_diagnostics": "null-residual",
    "profile_predictor_diagnostics": "profile",
    "size_spline_sensitivity_diagnostics": "size-spline",
}


def _abbreviate_source(name: str) -> str:
    bucket = ""
    if "/" in name:
        bucket_part, stem = name.split("/", 1)
        bucket = bucket_part.replace("tables_", "") + ":"
    else:
        stem = name
    return bucket + SOURCE_ABBREV.get(stem, stem.replace("_diagnostics", ""))


def make_supp_fig_8(root: Path, out_dir: Path) -> dict[str, Path] | None:
    style = load_style(root)
    style.set_theme(context="paper", font_scale=0.78)

    diag = gather_diagnostics(root)
    if diag.empty:
        return None
    diag = diag.copy()
    diag["short_source"] = diag["source"].apply(_abbreviate_source)

    import matplotlib.pyplot as plt
    fig, axes = style.new_figure(
        width="double",
        height_in=6.4,
        nrows=2,
        ncols=2,
        font_scale=0.78,
    )

    # Panel A — alpha by outcome, jittered points coloured by bucket.
    ax = axes[0, 0]
    sub = diag.dropna(subset=["alpha"]).copy()
    if not sub.empty:
        outcomes = sorted(sub["outcome"].dropna().unique())
        outcome_pos = {o: i for i, o in enumerate(outcomes)}
        buckets = sorted(sub["bucket"].dropna().unique())
        bucket_colours = dict(zip(buckets, plt.get_cmap("tab10").colors[: len(buckets)]))
        for _, row in sub.iterrows():
            x = outcome_pos.get(row["outcome"], -1)
            if x < 0:
                continue
            jitter = (hash(row["source"]) % 21 - 10) / 60.0
            ax.scatter(
                x + jitter,
                row["alpha"],
                s=18,
                color=bucket_colours.get(row["bucket"], "#777777"),
                edgecolor="white",
                linewidth=0.3,
                alpha=0.85,
                label=row["bucket"],
            )
        ax.set_xticks(list(outcome_pos.values()))
        ax.set_xticklabels(
            [OUTCOME_LABELS.get(o, o).replace("\n", " ") for o in outcomes],
            rotation=10,
        )
        ax.set_yscale("log")
        ax.set_ylabel("Dispersion alpha (log)")
        ax.set_title("Negative-binomial dispersion", pad=4)
        ax.grid(axis="y", color="#dddddd", linewidth=0.5, which="both")
        handles, labels = ax.get_legend_handles_labels()
        seen: dict[str, object] = {}
        for h, l in zip(handles, labels):
            seen.setdefault(l, h)
        ax.legend(seen.values(), seen.keys(), fontsize=6, frameon=False, loc="upper right")

    # Panel B — log-likelihood by bucket (median across model fits within bucket).
    ax = axes[0, 1]
    if "llf" in diag.columns:
        sub = diag.dropna(subset=["llf"]).copy()
        sub["label"] = sub["bucket"] + " / " + sub["outcome"].astype(str)
        agg = sub.groupby("label", as_index=False)["llf"].median().sort_values("llf")
        if not agg.empty:
            ax.barh(range(len(agg)), agg["llf"], color="#888888")
            ax.set_yticks(range(len(agg)))
            ax.set_yticklabels(agg["label"], fontsize=6)
            ax.set_xlabel("Median log likelihood")
            ax.set_title("Log-likelihood (median by bucket/outcome)", pad=4)
            ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    # Panel C — median n observations by bucket/outcome (one bar per bucket).
    ax = axes[1, 0]
    if "n_obs" in diag.columns:
        sub = diag.dropna(subset=["n_obs"]).copy()
        sub["label"] = sub["bucket"] + " / " + sub["outcome"].astype(str)
        agg = sub.groupby("label", as_index=False)["n_obs"].median().sort_values("n_obs")
        if not agg.empty:
            ax.barh(range(len(agg)), agg["n_obs"], color="#444444")
            ax.set_yticks(range(len(agg)))
            ax.set_yticklabels(agg["label"], fontsize=6)
            ax.set_xlabel("Median n observations")
            ax.set_title("Observations by bucket/outcome", pad=4)
            ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    # Panel D — tail caps for the canonical main/size-adjusted fits in the
    # tail-sensitivity buckets only.
    ax = axes[1, 1]
    tail_cols = [
        c
        for c in ["n_tail_excluded", "winsorise_cap", "tail_exclude_cap"]
        if c in diag.columns
    ]
    if tail_cols:
        canonical_stems = {
            "main_effects_diagnostics",
            "size_adjusted_spread_diagnostics",
        }

        def _stem(src: str) -> str:
            return src.split("/", 1)[1] if "/" in src else src

        sub = diag[
            diag["bucket"].isin(["winsorise99", "exclude_tail995"])
            & diag["source"].apply(lambda s: _stem(s) in canonical_stems)
        ][["bucket", "source", "outcome", *tail_cols]].copy()
        rows = []
        for _, row in sub.iterrows():
            for col in tail_cols:
                val = row.get(col)
                if pd.isna(val):
                    continue
                try:
                    val_f = float(val)
                except (TypeError, ValueError):
                    continue
                if val_f == 0:
                    continue
                short_outcome = {"cluster_size": "size", "geographic_spread": "spread"}.get(
                    str(row["outcome"]), str(row["outcome"])
                )
                lbl = f"{row['bucket']} ({short_outcome})"
                rows.append({"label": lbl, "metric": col, "value": val_f})
        if rows:
            tail_df = pd.DataFrame(rows)
            metrics = sorted(tail_df["metric"].unique())
            metric_colours = dict(zip(metrics, ["#4e79a7", "#f28e2b", "#e15759", "#59a14f"]))
            ordered_labels = sorted(tail_df["label"].unique())
            for i, lbl in enumerate(ordered_labels):
                for m in metrics:
                    val = tail_df[(tail_df["label"] == lbl) & (tail_df["metric"] == m)]["value"]
                    if val.empty:
                        continue
                    ax.scatter(
                        val.iloc[0],
                        i,
                        s=26,
                        color=metric_colours[m],
                        edgecolor="white",
                        linewidth=0.3,
                        label=m,
                    )
            ax.set_yticks(range(len(ordered_labels)))
            ax.set_yticklabels(ordered_labels, fontsize=6.5)
            ax.set_xscale("symlog")
            ax.set_xlabel("Cap or n excluded (symlog)")
            ax.set_title("Tail-rule metadata", pad=4)
            ax.grid(axis="x", color="#dddddd", linewidth=0.5, which="both")
            handles, labels = ax.get_legend_handles_labels()
            seen = {}
            for h, l in zip(handles, labels):
                seen.setdefault(l, h)
            ax.legend(seen.values(), seen.keys(), fontsize=6, frameon=False, loc="lower right")
        else:
            ax.axis("off")
    else:
        ax.axis("off")

    style.add_panel_labels(axes.ravel(), x=-0.20, y=1.10, size=9)
    fig.subplots_adjust(left=0.22, right=0.99, top=0.93, bottom=0.10, wspace=0.55, hspace=0.40)
    return save_all(style, fig, out_dir / "supp_fig8_model_diagnostics", "double", 6.4)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def run(root: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, dict[str, Path] | None] = {}
    builders = [
        ("fig1_population_measures", make_figure_1),
        ("fig2_main_pooled_effects", make_figure_2),
        ("fig3_wave_heterogeneity", make_figure_3),
        ("fig4_robustness_dashboard", make_figure_4),
        ("supp_fig1_outcome_distributions", make_supp_fig_1),
        ("supp_fig2_mixing_distributions", make_supp_fig_2),
        ("supp_fig3_cluster_size_mixing_boxplots", make_cluster_size_mixing_boxplots),
        ("supp_fig4_observed_expected_matrices", make_supp_fig_4),
        ("supp_fig5_size_adjustment", make_supp_fig_5),
        ("supp_fig6_domain_mixing_predictors", make_supp_fig_6),
        ("supp_fig7_profile_predictors", make_supp_fig_7),
        ("supp_fig8_model_diagnostics", make_supp_fig_8),
    ]
    for name, builder in builders:
        try:
            outputs[name] = builder(root, out_dir)
        except Exception as exc:  # noqa: BLE001 - log and continue
            print(f"[chapter1.manuscript] FAILED {name}: {exc}")
            outputs[name] = None

    for name, saved in outputs.items():
        if not saved:
            print(f"[chapter1.manuscript] skipped {name}")
            continue
        for fmt, path in saved.items():
            print(f"[chapter1.manuscript] wrote {name} {fmt}: {path}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to chapter1/manuscript/figures.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    root = args.root.resolve()
    out_dir = (
        args.out_dir.resolve()
        if args.out_dir is not None
        else root / "chapter1" / "manuscript" / "figures"
    )
    run(root=root, out_dir=out_dir)


if __name__ == "__main__":
    main()
