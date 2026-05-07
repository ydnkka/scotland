"""Create publication-ready figures for the Part 2 vaccination characterisation analysis.

Outputs are written to ``part2/manuscript/figures`` as PDF, PNG, and TIFF.
The script uses the shared project plotting module at ``utils/style.py``.

Main figures
------------
fig1  Vaccinated-case proportion over time, by JCVI rollout age group and SIMD quintile
fig2  Cluster vaccination profile and mean proportion vaccinated by wave and cluster size
fig3  Vaccination-status mixing category proportions by wave
fig4  Demographic mixing categories (SIMD, age, sex, joint-profile) by wave
fig5  Geographic dispersion category by wave
fig6  Booster coverage and days since last dose by wave and SIMD deprivation quintile

Supplementary figures
---------------------
supp_fig1  Weekly evolution of vaccination-status mixing category fractions
supp_fig2  Deprivation gradient in booster coverage and dose recency by SIMD domain and wave
supp_fig3  Cross-category heatmaps: fraction "more mix" by cluster size and SIMD quintile

Figures 4 and supp_fig3 require pyarrow (available in the PhD conda environment) to read
the cluster-level parquet caches.  They are skipped gracefully if pyarrow is absent.

Run from the repository root:

    conda run -n PhD python part2/manuscript/make_figures.py
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Ordering / label / colour constants
# ---------------------------------------------------------------------------

WAVE_ORDER = [
    "B.1.177", "Alpha", "Delta",
    "BA.1", "BA.2", "BA.4", "BA.5", "BQ.1", "XBB", "Other",
]
MAIN_WAVES = [w for w in WAVE_ORDER if w != "Other"]

JCVI_AGE_ORDER = [
    "00-14", "15-19", "20-29", "30-39", "40-49",
    "50-54", "55-59", "60-64", "65-69", "70-74", "75+",
]

SIMD_QUINTILE_INT = [1, 2, 3, 4, 5]
SIMD_QUINTILE_STR_MAP = {
    "1 most deprived": 1,
    "2": 2,
    "3": 3,
    "4": 4,
    "5 least deprived": 5,
}

# Vaccination-status mixing categories (non-singleton clusters only)
MIXING_CAT_ORDER = ["homogeneous", "baseline", "mixed"]
MIXING_CAT_LABELS = {
    "homogeneous": "Homogeneous",
    "baseline":    "Baseline",
    "mixed":       "Mixed",
}
MIXING_CAT_COLOURS = {
    "homogeneous": "#2b83ba",
    "baseline":    "#b0b0b0",
    "mixed":       "#d7191c",
}

# Cluster vaccination profile
PROFILE_ORDER = ["none vaccinated", "mixed vaccination", "all vaccinated"]
PROFILE_LABELS = {
    "none vaccinated": "None vaccinated",
    "mixed vaccination": "Mixed",
    "all vaccinated":  "All vaccinated",
}
PROFILE_COLOURS = {
    "none vaccinated": "#cccccc",
    "mixed vaccination": "#f28e2b",
    "all vaccinated":  "#4e79a7",
}

# Cluster size categories
SIZE_ORDER  = ["small/moderate", "large", "very large"]
SIZE_LABELS = {
    "small/moderate": "Small/moderate",
    "large":          "Large",
    "very large":     "Very large",
}
SIZE_COLOURS = {
    "small/moderate": "#abdda4",
    "large":          "#fdae61",
    "very large":     "#d7191c",
}

# SIMD domain order (matches Part 1)
DOMAIN_ORDER = [
    "overall", "income", "employment", "education",
    "health", "access", "crime", "housing",
]
DOMAIN_LABELS = {
    "overall":    "Overall",
    "income":     "Income",
    "employment": "Employment",
    "education":  "Education",
    "health":     "Health",
    "access":     "Access",
    "crime":      "Crime",
    "housing":    "Housing",
}

# Minimum cases per weekly stratum before masking proportion as NaN
MIN_N_WEEKLY = 20
SMOOTH_WINDOW = 4   # weeks, centred rolling mean


# ---------------------------------------------------------------------------
# Shared helpers  (mirror Part 1 make_figures.py patterns)
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


def setup_environment() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)


def save_all(style, fig, out_base: Path, width: str, height_in: float) -> None:
    style.save_figure(
        fig, out_base,
        width=width, height_in=height_in, dpi=600,
        save_pdf=True, save_png=True, save_tiff=True,
    )


def waves_present(df: pd.DataFrame, col: str = "wave_group") -> list[str]:
    """Return WAVE_ORDER entries that exist in *df[col]*."""
    present = set(df[col].dropna().unique())
    return [w for w in WAVE_ORDER if w in present]


# ---------------------------------------------------------------------------
# Figure 1 — Vaccinated-case proportion over time
# ---------------------------------------------------------------------------

def plot_vaccinated_cases_over_time(
    style, case_weekly: pd.DataFrame, out_dir: Path
) -> None:
    """Double-column, 2-panel: (A) JCVI age group, (B) SIMD quintile."""

    case_weekly = case_weekly.copy()
    case_weekly["case_week"] = pd.to_datetime(case_weekly["case_week"])

    # Build age-group palette: plasma from youngest (purple) → oldest (yellow)
    n_age = len(JCVI_AGE_ORDER)
    age_palette = {
        g: plt.cm.plasma(i / (n_age - 1))
        for i, g in enumerate(JCVI_AGE_ORDER)
    }

    fig, axes = style.new_figure(
        width="double", height_in=3.6, nrows=1, ncols=2,
        font_scale=0.85,
    )

    # --- Panel A: JCVI rollout age group ---
    ax = axes[0]
    age_df = (
        case_weekly[case_weekly["stratum_type"] == "vaccination_age_group"]
        .sort_values(["stratum", "case_week"])
        .copy()
    )
    age_df.loc[age_df["n_vaccination_known"] < MIN_N_WEEKLY, "prop_vaccinated"] = np.nan

    for age_group in JCVI_AGE_ORDER:
        sub = age_df[age_df["stratum"] == age_group]
        if sub.empty:
            continue
        smooth = (
            sub["prop_vaccinated"]
            .rolling(SMOOTH_WINDOW, center=True, min_periods=2)
            .mean()
        )
        ax.plot(
            sub["case_week"], smooth * 100,
            color=age_palette[age_group],
            linewidth=1.1,
            label=age_group,
        )

    ax.set_xlim(age_df["case_week"].min(), age_df["case_week"].max())
    ax.set_ylim(0, 105)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.set_ylabel("Vaccinated sequenced cases (%)")
    ax.set_title("By JCVI rollout age group", pad=4)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    # Compact 2-column legend inside panel
    ax.legend(
        title="Age group",
        fontsize=6,
        title_fontsize=6.5,
        loc="upper left",
        frameon=False,
        handlelength=1.0,
        handletextpad=0.4,
        columnspacing=0.6,
        borderpad=0.2,
    )

    # --- Panel B: SIMD quintile ---
    ax = axes[1]
    simd_df = (
        case_weekly[case_weekly["stratum_type"] == "simd_quintile"]
        .sort_values(["stratum", "case_week"])
        .copy()
    )
    simd_df.loc[simd_df["n_vaccination_known"] < MIN_N_WEEKLY, "prop_vaccinated"] = np.nan
    simd_df["quintile_int"] = simd_df["stratum"].map(SIMD_QUINTILE_STR_MAP)

    for q in SIMD_QUINTILE_INT:
        label_str = {v: k for k, v in SIMD_QUINTILE_STR_MAP.items()}[q]
        sub = simd_df[simd_df["quintile_int"] == q]
        if sub.empty:
            continue
        smooth = (
            sub["prop_vaccinated"]
            .rolling(SMOOTH_WINDOW, center=True, min_periods=2)
            .mean()
        )
        # leg_label = (
        #     f"Q{q} (most deprived)" if q == 1
        #     else f"Q{q} (least deprived)" if q == 5
        #     else f"Q{q}"
        # )
        ax.plot(
            sub["case_week"], smooth * 100,
            color=style.SIMD_QUINTILE_PALETTE[q],
            linewidth=1.3,
            label=label_str,
        )

    ax.set_xlim(simd_df["case_week"].min(), simd_df["case_week"].max())
    ax.set_ylim(0, 105)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.set_ylabel("")
    ax.set_title("By overall SIMD deprivation quintile", pad=4)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.legend(
        title="SIMD quintile",
        fontsize=6.5,
        title_fontsize=7,
        loc="upper left",
        frameon=False,
        handlelength=1.2,
        borderpad=0.2,
    )

    style.add_panel_labels(axes, x=-0.07, y=1.12, size=9)
    fig.text(0.5, 0.01, "Calendar week", ha="center", fontsize=8)
    fig.subplots_adjust(
        left=0.09, right=0.99, top=0.90, bottom=0.12, wspace=0.22,
    )
    save_all(style, fig, out_dir / "fig1_vaccinated_cases_over_time", "double", 3.6)


# ---------------------------------------------------------------------------
# Figure 2 — Vaccination-status mixing categories by wave
# ---------------------------------------------------------------------------

def plot_vaccination_mixing_by_wave(
    style, cluster_wave_cat: pd.DataFrame, out_dir: Path
) -> None:
    """Onehalf-column stacked proportional bar chart, non-singleton clusters only."""

    data = cluster_wave_cat[
        cluster_wave_cat["category_variable"] == "vaccination_mixing_category"
    ].copy()
    # Drop singletons (not available)
    data = data[data["category"].isin(MIXING_CAT_ORDER)]

    waves = waves_present(data)

    # Compute proportions per wave
    totals = data.groupby("wave_group")["n_clusters"].sum().rename("total")
    data = data.join(totals, on="wave_group")
    data["fraction"] = data["n_clusters"] / data["total"]

    pivot = (
        data.pivot_table(
            index="wave_group", columns="category",
            values="fraction", aggfunc="first",
        )
        .reindex(index=waves, columns=MIXING_CAT_ORDER)
        .fillna(0.0)
    )
    counts = (
        data.groupby("wave_group")["n_clusters"].sum()
        .reindex(waves)
    )

    fig, ax = style.new_figure(
        width="onehalf", height_in=3.8, font_scale=0.85,
    )

    x = np.arange(len(waves))
    bottoms = np.zeros(len(waves))

    for cat in MIXING_CAT_ORDER:
        vals = pivot[cat].to_numpy()
        ax.bar(
            x, vals * 100,
            bottom=bottoms * 100,
            color=MIXING_CAT_COLOURS[cat],
            label=MIXING_CAT_LABELS[cat],
            width=0.72,
            edgecolor="white",
            linewidth=0.4,
        )
        bottoms += vals

    # Annotate total cluster count above each bar
    for i, wave in enumerate(waves):
        n = counts[wave]
        ax.text(
            i, 101.5,
            f"n={n:,}" if n < 1000 else f"n={n/1000:.1f}k",
            ha="center", va="bottom", fontsize=5.5, color="#444444",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(waves, rotation=40, ha="right", fontsize=7)
    ax.set_ylim(0, 110)
    ax.set_ylabel("Non-singleton clusters (%)")
    ax.set_title("Vaccination-status mixing by wave", pad=4)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.legend(
        title="Mixing category",
        fontsize=7,
        title_fontsize=7.5,
        loc="lower right",
        frameon=False,
        handlelength=1.0,
    )

    fig.subplots_adjust(left=0.14, right=0.97, top=0.92, bottom=0.20)
    save_all(style, fig, out_dir / "fig3_vaccination_mixing_by_wave", "onehalf", 3.8)


# ---------------------------------------------------------------------------
# Figure 3 — Cluster vaccination profile and mean proportion vaccinated
#            by wave and cluster size category
# ---------------------------------------------------------------------------

def plot_cluster_vaccination_by_wave_and_category(
    style, cluster_wave_cat: pd.DataFrame, out_dir: Path
) -> None:
    """Double-column, 2-panel figure.

    Panel A — stacked proportional bar: vaccination profile by wave.
    Panel B — mean cluster proportion vaccinated by wave, coloured by size category.
    """

    # --- Panel A data: vaccination profile by wave ---
    prof_data = cluster_wave_cat[
        cluster_wave_cat["category_variable"] == "cluster_vaccination_profile"
    ].copy()
    prof_data = prof_data[prof_data["category"].isin(PROFILE_ORDER)]
    waves_a = waves_present(prof_data)

    totals_a = prof_data.groupby("wave_group")["n_clusters"].sum().rename("total")
    prof_data = prof_data.join(totals_a, on="wave_group")
    prof_data["fraction"] = prof_data["n_clusters"] / prof_data["total"]

    pivot_a = (
        prof_data.pivot_table(
            index="wave_group", columns="category",
            values="fraction", aggfunc="first",
        )
        .reindex(index=waves_a, columns=PROFILE_ORDER)
        .fillna(0.0)
    )

    # --- Panel B data: mean cluster proportion vaccinated by size × wave ---
    size_data = cluster_wave_cat[
        cluster_wave_cat["category_variable"] == "cluster_size_category"
    ].copy()
    size_data = size_data[size_data["category"].isin(SIZE_ORDER)]
    waves_b = waves_present(size_data)

    fig, axes = style.new_figure(
        width="double", height_in=3.8, nrows=1, ncols=2,
        font_scale=0.85,
    )

    # --- Panel A: stacked bar ---
    ax = axes[0]
    x = np.arange(len(waves_a))
    bottoms = np.zeros(len(waves_a))

    for prof in PROFILE_ORDER:
        vals = pivot_a[prof].to_numpy()
        ax.bar(
            x, vals * 100,
            bottom=bottoms * 100,
            color=PROFILE_COLOURS[prof],
            label=PROFILE_LABELS[prof],
            width=0.72,
            edgecolor="white",
            linewidth=0.4,
        )
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels(waves_a, rotation=40, ha="right", fontsize=7)
    ax.set_ylim(0, 103)
    ax.set_ylabel("Clusters (%)")
    ax.set_title("Cluster vaccination profile by wave", pad=4)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.legend(
        title="Profile",
        fontsize=7,
        title_fontsize=7.5,
        loc="center left",
        frameon=False,
        handlelength=1.0,
    )

    # --- Panel B: dot plot mean proportion vaccinated by wave × size ---
    ax = axes[1]
    x = np.arange(len(waves_b))
    offsets = {"small/moderate": -0.22, "large": 0.0, "very large": 0.22}

    for size_cat in SIZE_ORDER:
        sub = size_data[size_data["category"] == size_cat].copy()
        sub = sub.set_index("wave_group").reindex(waves_b)
        y = sub["mean_prop_vaccinated"].to_numpy() * 100
        ax.scatter(
            x + offsets[size_cat], y,
            color=SIZE_COLOURS[size_cat],
            s=22,
            zorder=3,
            label=SIZE_LABELS[size_cat],
            edgecolor="white",
            linewidth=0.3,
        )
        ax.plot(
            x + offsets[size_cat], y,
            color=SIZE_COLOURS[size_cat],
            linewidth=0.7,
            linestyle="--",
            alpha=0.6,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(waves_b, rotation=40, ha="right", fontsize=7)
    ax.set_ylim(0, 105)
    ax.set_ylabel("Mean cluster proportion vaccinated (%)")
    ax.set_title("Mean vaccination by wave and cluster size", pad=4)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.legend(
        title="Cluster size",
        fontsize=7,
        title_fontsize=7.5,
        loc="upper left",
        frameon=False,
        handlelength=1.2,
    )

    style.add_panel_labels(axes, x=-0.07, y=1.12, size=9)
    fig.subplots_adjust(
        left=0.09, right=0.99, top=0.90, bottom=0.22, wspace=0.28,
    )
    save_all(
        style, fig,
        out_dir / "fig2_cluster_vaccination_by_wave_and_category",
        "double", 3.8,
    )


# ---------------------------------------------------------------------------
# Figure 4 — Booster coverage and dose recency by wave × SIMD quintile
# ---------------------------------------------------------------------------

def plot_dose_recency_by_simd(
    style, cluster_wave_simd: pd.DataFrame, out_dir: Path
) -> None:
    """Double-column 2-panel heatmap.

    Panel A — mean booster coverage among vaccinated cluster members (%).
    Panel B — mean days since last prior vaccination.

    Cells where mean_prop_vaccinated < 0.05 are masked (insufficient vaccination).
    """

    overall = cluster_wave_simd[
        cluster_wave_simd["simd_domain"] == "overall"
    ].copy()

    waves = [w for w in WAVE_ORDER if w in set(overall["wave_group"])]
    quintiles = [1, 2, 3, 4, 5]

    def _pivot(metric: str) -> np.ndarray:
        piv = (
            overall.pivot_table(
                index="simd_quintile", columns="wave_group",
                values=metric, aggfunc="first",
            )
            .reindex(index=quintiles, columns=waves)
            .to_numpy(dtype=float)
            .copy()
        )
        # Mask cells where vaccination was negligible
        mask_piv = (
            overall.pivot_table(
                index="simd_quintile", columns="wave_group",
                values="mean_prop_vaccinated", aggfunc="first",
            )
            .reindex(index=quintiles, columns=waves)
            .to_numpy(dtype=float)
        )
        piv[mask_piv < 0.05] = np.nan
        return piv

    booster_mat  = _pivot("mean_prop_boosted_vaccinated_members") * 100
    days_mat     = _pivot("mean_days_since_vaccination")

    fig, axes = style.new_figure(
        width="double", height_in=3.2, nrows=1, ncols=2,
        font_scale=0.85,
    )

    ytick_labels = [
        "Q1 (most\ndeprived)", "Q2", "Q3", "Q4", "Q5 (least\ndeprived)",
    ]

    def _draw_heatmap(ax, matrix, cmap, vmin, vmax, title, cbarlabel, fmt=".0f"):
        # Grey background for masked cells
        ax.set_facecolor("#e8e8e8")
        img = ax.imshow(
            matrix, cmap=cmap, aspect="auto",
            vmin=vmin, vmax=vmax,
            interpolation="nearest",
        )
        ax.set_xticks(np.arange(len(waves)))
        ax.set_xticklabels(waves, rotation=45, ha="right", fontsize=6.5)
        ax.set_yticks(np.arange(len(quintiles)))
        ax.set_yticklabels(ytick_labels, fontsize=7)
        ax.tick_params(length=0)

        # Cell gridlines
        for y in np.arange(len(quintiles) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.5)
        for x in np.arange(len(waves) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.5)

        # Annotate cell values
        for r in range(len(quintiles)):
            for c in range(len(waves)):
                val = matrix[r, c]
                if not np.isfinite(val):
                    continue
                text_col = "white" if (val - vmin) / (vmax - vmin) > 0.6 else "#333333"
                ax.text(
                    c, r, f"{val:{fmt}}",
                    ha="center", va="center",
                    fontsize=5, color=text_col,
                )

        ax.set_title(title, pad=4)
        return img

    img_a = _draw_heatmap(
        axes[0], booster_mat,
        cmap="YlGnBu",
        vmin=0, vmax=100,
        title="Booster coverage among\nvaccinated cluster members (%)",
        cbarlabel="Booster coverage (%)",
    )
    img_b = _draw_heatmap(
        axes[1], days_mat,
        cmap="YlOrBr",
        vmin=0, vmax=300,
        title="Mean days since last\nprior dose (vaccinated members)",
        cbarlabel="Days since last dose",
        fmt=".0f",
    )

    # Colorbars
    fig.subplots_adjust(
        left=0.13, right=0.87, top=0.85, bottom=0.22, wspace=0.38,
    )
    for img, ax in [(img_a, axes[0]), (img_b, axes[1])]:
        pos = ax.get_position()
        cbar_ax = fig.add_axes([pos.x1 + 0.01, pos.y0, 0.012, pos.height])
        cb = fig.colorbar(img, cax=cbar_ax)
        cb.ax.tick_params(labelsize=6)

    style.add_panel_labels(axes, x=-0.14, y=1.15, size=9)
    save_all(
        style, fig,
        out_dir / "fig6_dose_recency_by_simd",
        "double", 3.2,
    )


# ---------------------------------------------------------------------------
# Supplementary Figure 1 — Weekly vaccination-status mixing fractions
# ---------------------------------------------------------------------------

def plot_weekly_mixing_evolution(
    style, cluster_weekly_cat: pd.DataFrame, out_dir: Path
) -> None:
    """Double-column stacked-area chart showing mixing category fractions over time."""

    cluster_weekly_cat = cluster_weekly_cat.copy()
    cluster_weekly_cat["cluster_week"] = pd.to_datetime(cluster_weekly_cat["cluster_week"])

    data = cluster_weekly_cat[
        cluster_weekly_cat["category_variable"] == "vaccination_mixing_category"
    ].copy()
    data = data[data["category"].isin(MIXING_CAT_ORDER)]

    # Compute total non-singleton clusters per week
    totals = (
        data.groupby("cluster_week")["n_clusters"].sum()
        .rename("total")
    )
    data = data.join(totals, on="cluster_week")
    data["fraction"] = data["n_clusters"] / data["total"]

    pivot = (
        data.pivot_table(
            index="cluster_week", columns="category",
            values="fraction", aggfunc="first",
        )
        .reindex(columns=MIXING_CAT_ORDER)
        .fillna(0.0)
        .sort_index()
    )

    # 4-week rolling smooth
    smooth = pivot.rolling(SMOOTH_WINDOW, center=True, min_periods=2).mean()

    fig, ax = style.new_figure(
        width="double", height_in=3.0, font_scale=0.85,
    )

    bottoms = np.zeros(len(smooth))
    for cat in MIXING_CAT_ORDER:
        vals = smooth[cat].to_numpy()
        # Replace NaN with 0 for stacking (early/late edge artefacts)
        vals = np.where(np.isfinite(vals), vals, 0.0)
        ax.fill_between(
            smooth.index, bottoms * 100, (bottoms + vals) * 100,
            color=MIXING_CAT_COLOURS[cat],
            label=MIXING_CAT_LABELS[cat],
            alpha=0.90,
            linewidth=0,
        )
        bottoms += vals

    ax.set_xlim(smooth.index.min(), smooth.index.max())
    ax.set_ylim(0, 100)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.set_xlabel("Calendar week")
    ax.set_ylabel("Non-singleton clusters (%)")
    ax.set_title(
        "Weekly vaccination-status mixing category fractions\n"
        "(4-week rolling mean, non-singleton clusters)",
        pad=4,
    )
    ax.grid(axis="y", color="white", linewidth=0.7, alpha=0.5)
    ax.legend(
        title="Mixing category",
        fontsize=7.5,
        title_fontsize=8,
        loc="upper left",
        frameon=False,
        handlelength=1.2,
    )

    fig.subplots_adjust(left=0.09, right=0.99, top=0.86, bottom=0.14)
    save_all(
        style, fig,
        out_dir / "supp_fig1_weekly_mixing_evolution",
        "double", 3.0,
    )


# ---------------------------------------------------------------------------
# Supplementary Figure 2 — Deprivation gradient in dose metrics by domain
# ---------------------------------------------------------------------------

def plot_domain_dose_gradient(
    style, cluster_wave_simd: pd.DataFrame, out_dir: Path
) -> None:
    """Double-column 2-panel heatmap.

    Rows: SIMD domains.  Columns: waves.
    Colour: Q5 minus Q1 gradient in (A) booster coverage (pp) and
    (B) mean days since last dose.
    Positive booster gradient = less-deprived clusters have higher booster coverage.
    Positive days gradient = less-deprived clusters have more days since last dose.
    """

    waves = [w for w in WAVE_ORDER if w in set(cluster_wave_simd["wave_group"])]
    domains = [d for d in DOMAIN_ORDER if d in set(cluster_wave_simd["simd_domain"])]

    def _gradient_matrix(metric: str) -> np.ndarray:
        q1 = cluster_wave_simd[cluster_wave_simd["simd_quintile"] == 1].copy()
        q5 = cluster_wave_simd[cluster_wave_simd["simd_quintile"] == 5].copy()
        q1 = q1.pivot_table(
            index="simd_domain", columns="wave_group",
            values=metric, aggfunc="first",
        ).reindex(index=domains, columns=waves)
        q5 = q5.pivot_table(
            index="simd_domain", columns="wave_group",
            values=metric, aggfunc="first",
        ).reindex(index=domains, columns=waves)
        # Mask where vaccination was negligible in Q1
        q1_vacc = cluster_wave_simd[cluster_wave_simd["simd_quintile"] == 1].pivot_table(
            index="simd_domain", columns="wave_group",
            values="mean_prop_vaccinated", aggfunc="first",
        ).reindex(index=domains, columns=waves)
        grad = (q5 - q1).to_numpy(dtype=float).copy()
        mask = q1_vacc.to_numpy(dtype=float) < 0.05
        grad[mask] = np.nan
        return grad

    booster_grad = _gradient_matrix("mean_prop_boosted_vaccinated_members") * 100
    days_grad    = _gradient_matrix("mean_days_since_vaccination")

    fig, axes = style.new_figure(
        width="double", height_in=4.2, nrows=1, ncols=2,
        font_scale=0.82,
    )

    domain_tick_labels = [DOMAIN_LABELS[d] for d in domains]

    def _draw(ax, matrix, cmap, vlim, title, fmt="+.0f"):
        ax.set_facecolor("#e8e8e8")
        img = ax.imshow(
            matrix, cmap=cmap, aspect="auto",
            vmin=-vlim, vmax=vlim,
            interpolation="nearest",
        )
        ax.set_xticks(np.arange(len(waves)))
        ax.set_xticklabels(waves, rotation=45, ha="right", fontsize=6.5)
        ax.set_yticks(np.arange(len(domains)))
        ax.set_yticklabels(domain_tick_labels, fontsize=7)
        ax.tick_params(length=0)
        for y in np.arange(len(domains) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.5)
        for x in np.arange(len(waves) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.5)
        for r in range(len(domains)):
            for c in range(len(waves)):
                val = matrix[r, c]
                if not np.isfinite(val):
                    continue
                text_col = "white" if abs(val) / vlim > 0.6 else "#333333"
                ax.text(
                    c, r, f"{val:{fmt}}",
                    ha="center", va="center",
                    fontsize=5.5, color=text_col,
                )
        ax.set_title(title, pad=4)
        return img

    img_a = _draw(
        axes[0], booster_grad,
        cmap="RdBu", vlim=20,
        title="Booster coverage gradient\n(Q5 − Q1, pp)",
    )
    img_b = _draw(
        axes[1], days_grad,
        cmap="RdBu_r", vlim=40,
        title="Days since last dose gradient\n(Q5 − Q1, days)",
    )

    fig.subplots_adjust(
        left=0.13, right=0.87, top=0.87, bottom=0.24, wspace=0.52,
    )
    for img, ax in [(img_a, axes[0]), (img_b, axes[1])]:
        pos = ax.get_position()
        cbar_ax = fig.add_axes([pos.x1 + 0.01, pos.y0, 0.012, pos.height])
        cb = fig.colorbar(img, cax=cbar_ax)
        cb.ax.tick_params(labelsize=6)

    style.add_panel_labels(axes, x=-0.14, y=1.12, size=9)
    save_all(
        style, fig,
        out_dir / "supp_fig2_domain_dose_gradient",
        "double", 4.2,
    )


# ---------------------------------------------------------------------------
# Shared helper — load and merge parquet caches for cluster-level analyses
# ---------------------------------------------------------------------------

DEMO_MIX_CATS = [
    "simd_mixing_category",
    "age_mixing_category",
    "sex_mixing_category",
    "profile_mixing_category",
]
DEMO_MIX_LABELS = {
    "simd_mixing_category":    "SIMD deprivation mixing",
    "age_mixing_category":     "Age mixing",
    "sex_mixing_category":     "Sex mixing",
    "profile_mixing_category": "Joint profile mixing",
}
MIX_CAT_POOLED_ORDER = ["less mix", "baseline", "more mix"]
MIX_CAT_POOLED_COLOURS = {
    "less mix": "#2b83ba",
    "baseline": "#b0b0b0",
    "more mix": "#d7191c",
}
MIX_CAT_POOLED_LABELS = {
    "less mix": "Less mix",
    "baseline": "Baseline",
    "more mix": "More mix",
}


def load_cluster_demo_mix(cache_dir: Path) -> pd.DataFrame | None:
    """Return non-singleton clusters with wave + all demographic mixing categories.

    Merges ``vaccination_cluster_table.parquet`` (wave info) with
    ``cluster_categories.parquet`` (demographic mixing categories).
    Returns None if pyarrow is unavailable.
    """
    try:
        vct = pd.read_parquet(
            cache_dir / "vaccination_cluster_table.parquet",
            columns=["cluster_id", "wave_group", "cluster_size_category",
                     "geographic_dispersion_category"],
        )
        cats = pd.read_parquet(
            cache_dir / "cluster_categories.parquet",
            columns=["cluster_id", "simd_quintile", "simd_quintile_label",
                     "cluster_size_category", "geographic_dispersion_category",
                     *DEMO_MIX_CATS],
        )
    except Exception as exc:
        print(f"  [skip] Could not load parquet caches ({exc}). "
              "Install pyarrow in the active environment.")
        return None
    return vct.merge(cats, on="cluster_id", how="inner", suffixes=("", "_cat"))


# ---------------------------------------------------------------------------
# Figure 4 — Demographic mixing categories by wave (4 panels)
# ---------------------------------------------------------------------------

def plot_demographic_mixing_by_wave(
    style, cluster_demo: pd.DataFrame, out_dir: Path
) -> None:
    """Double-column 2×2 figure: one stacked bar panel per mixing type (wave on x-axis)."""

    waves = waves_present(cluster_demo, col="wave_group")

    fig, axes = style.new_figure(
        width="double", height_in=5.6, nrows=2, ncols=2,
        font_scale=0.85, sharey=True,
    )

    for idx, (mix_col, mix_label) in enumerate(DEMO_MIX_LABELS.items()):
        ax = axes.ravel()[idx]

        # Count per wave × category, restrict to valid mixing categories
        grp = (
            cluster_demo[cluster_demo[mix_col].isin(MIX_CAT_POOLED_ORDER)]
            .groupby(["wave_group", mix_col], observed=True)["cluster_id"]
            .count()
            .reset_index()
            .rename(columns={"cluster_id": "n_clusters", mix_col: "category"})
        )
        totals = grp.groupby("wave_group")["n_clusters"].transform("sum")
        grp["fraction"] = grp["n_clusters"] / totals

        pivot = (
            grp.pivot_table(
                index="wave_group", columns="category",
                values="fraction", aggfunc="first",
            )
            .reindex(index=waves, columns=MIX_CAT_POOLED_ORDER)
            .fillna(0.0)
        )
        counts = grp.groupby("wave_group")["n_clusters"].sum().reindex(waves)

        x = np.arange(len(waves))
        bottoms = np.zeros(len(waves))
        for cat in MIX_CAT_POOLED_ORDER:
            vals = pivot[cat].to_numpy()
            ax.bar(
                x, vals * 100,
                bottom=bottoms * 100,
                color=MIX_CAT_POOLED_COLOURS[cat],
                label=MIX_CAT_POOLED_LABELS[cat],
                width=0.72,
                edgecolor="white",
                linewidth=0.4,
            )
            bottoms += vals

        # Annotate total N
        for i, wave in enumerate(waves):
            n = counts.get(wave, 0)
            ax.text(
                i, 101.5,
                f"n={int(n/1000):.0f}k" if n >= 1000 else f"n={int(n)}",
                ha="center", va="bottom", fontsize=5.5, color="#444444",
            )

        ax.set_xticks(x)
        ax.set_xticklabels(waves, rotation=40, ha="right", fontsize=7)
        ax.set_ylim(0, 110)
        ax.set_title(mix_label, pad=4)
        ax.set_ylabel("Non-singleton clusters (%)" if idx % 2 == 0 else "")
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

        if idx == 1:
            ax.legend(
                title="Mixing category",
                fontsize=7, title_fontsize=7.5,
                loc="lower right", frameon=False, handlelength=1.0,
            )

    style.add_panel_labels(axes.ravel(), x=-0.06, y=1.12, size=9)
    fig.subplots_adjust(
        left=0.09, right=0.99, top=0.93, bottom=0.16, hspace=0.52, wspace=0.12,
    )
    save_all(style, fig, out_dir / "fig4_demographic_mixing_by_wave", "double", 5.6)


# ---------------------------------------------------------------------------
# Figure 5 — Geographic dispersion category by wave
# ---------------------------------------------------------------------------

GEO_ORDER = ["low/moderate dispersion", "large dispersion", "very large dispersion"]
GEO_LABELS = {
    "low/moderate dispersion": "Low/moderate",
    "large dispersion":        "Large",
    "very large dispersion":   "Very large",
}
GEO_COLOURS = {
    "low/moderate dispersion": "#abdda4",
    "large dispersion":        "#fdae61",
    "very large dispersion":   "#d7191c",
}


def plot_geographic_dispersion_by_wave(
    style, cluster_wave_cat: pd.DataFrame, out_dir: Path
) -> None:
    """Double-column 2-panel figure.

    Panel A — Stacked proportional bar: geographic dispersion category by wave.
    Panel B — Mean proportion vaccinated by wave × dispersion category (dot plot).
    """

    geo_data = cluster_wave_cat[
        cluster_wave_cat["category_variable"] == "geographic_dispersion_category"
    ].copy()
    geo_data = geo_data[geo_data["category"].isin(GEO_ORDER)]
    waves = waves_present(geo_data)

    # Panel A proportions
    totals_a = geo_data.groupby("wave_group")["n_clusters"].sum().rename("total")
    geo_data = geo_data.join(totals_a, on="wave_group")
    geo_data["fraction"] = geo_data["n_clusters"] / geo_data["total"]
    pivot_a = (
        geo_data.pivot_table(
            index="wave_group", columns="category",
            values="fraction", aggfunc="first",
        )
        .reindex(index=waves, columns=GEO_ORDER)
        .fillna(0.0)
    )

    fig, axes = style.new_figure(
        width="double", height_in=3.8, nrows=1, ncols=2,
        font_scale=0.85,
    )

    # --- Panel A: stacked bar ---
    ax = axes[0]
    x = np.arange(len(waves))
    bottoms = np.zeros(len(waves))
    for geo in GEO_ORDER:
        vals = pivot_a[geo].to_numpy()
        ax.bar(
            x, vals * 100,
            bottom=bottoms * 100,
            color=GEO_COLOURS[geo],
            label=GEO_LABELS[geo],
            width=0.72,
            edgecolor="white",
            linewidth=0.4,
        )
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels(waves, rotation=40, ha="right", fontsize=7)
    ax.set_ylim(0, 103)
    ax.set_ylabel("Clusters (%)")
    ax.set_title("Geographic dispersion category by wave", pad=4)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.legend(
        title="Dispersion",
        fontsize=7, title_fontsize=7.5,
        loc="lower right", frameon=False, handlelength=1.0,
    )

    # --- Panel B: mean prop vaccinated by wave × dispersion ---
    ax = axes[1]
    offsets = {
        "low/moderate dispersion": -0.22,
        "large dispersion":         0.0,
        "very large dispersion":    0.22,
    }
    x = np.arange(len(waves))
    for geo in GEO_ORDER:
        sub = (
            geo_data[geo_data["category"] == geo]
            .set_index("wave_group")
            .reindex(waves)
        )
        y = sub["mean_prop_vaccinated"].to_numpy() * 100
        ax.scatter(
            x + offsets[geo], y,
            color=GEO_COLOURS[geo],
            s=22, zorder=3,
            label=GEO_LABELS[geo],
            edgecolor="white", linewidth=0.3,
        )
        ax.plot(
            x + offsets[geo], y,
            color=GEO_COLOURS[geo],
            linewidth=0.7, linestyle="--", alpha=0.6,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(waves, rotation=40, ha="right", fontsize=7)
    ax.set_ylim(0, 105)
    ax.set_ylabel("Mean cluster proportion vaccinated (%)")
    ax.set_title("Mean vaccination by wave and dispersion", pad=4)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.legend(
        title="Dispersion",
        fontsize=7, title_fontsize=7.5,
        loc="upper left", frameon=False, handlelength=1.2,
    )

    style.add_panel_labels(axes, x=-0.07, y=1.12, size=9)
    fig.subplots_adjust(
        left=0.09, right=0.99, top=0.90, bottom=0.22, wspace=0.30,
    )
    save_all(style, fig, out_dir / "fig5_geographic_dispersion_by_wave", "double", 3.8)


# ---------------------------------------------------------------------------
# Supplementary Figures 3 & 4 — Cross-category heatmaps
# ---------------------------------------------------------------------------

def _cross_category_heatmap(
    style,
    cluster_demo: pd.DataFrame,
    out_dir: Path,
    cat_value: str,
    cmap: str,
    cbar_label: str,
    out_stem: str,
) -> None:
    """Shared 2×2 heatmap builder for cross-category mixing fraction figures.

    Rows = SIMD quintile (1–5), columns = cluster size category.
    Shared x/y tick labels (outer panels only) and a single shared colourbar.
    """
    size_order = SIZE_ORDER
    quintiles  = [1, 2, 3, 4, 5]
    qlabels    = ["Q1(most\ndeprived)", "Q2", "Q3", "Q4", "Q5(least\ndeprived)"]
    slabels    = [SIZE_LABELS[s] for s in size_order]

    fig, axes = style.new_figure(
        width="double", height_in=5.0, nrows=2, ncols=2,
        font_scale=0.85,
    )

    last_img = None

    for idx, (mix_col, mix_label) in enumerate(DEMO_MIX_LABELS.items()):
        ax = axes.ravel()[idx]
        row, col = divmod(idx, 2)

        valid = cluster_demo[cluster_demo[mix_col].isin(MIX_CAT_POOLED_ORDER)].copy()
        grp = (
            valid.groupby(["simd_quintile", "cluster_size_category", mix_col],
                          observed=True)["cluster_id"]
            .count()
            .reset_index()
            .rename(columns={"cluster_id": "n", mix_col: "cat"})
        )
        totals = grp.groupby(["simd_quintile", "cluster_size_category"],
                             observed=True)["n"].transform("sum")
        grp["frac"] = grp["n"] / totals

        mat = (
            grp[grp["cat"] == cat_value]
            .pivot_table(
                index="simd_quintile", columns="cluster_size_category",
                values="frac", aggfunc="first",
            )
            .reindex(index=quintiles, columns=size_order)
            .to_numpy(dtype=float)
            .copy()
        ) * 100   # 0–100 scale for both imshow and annotation

        ax.set_facecolor("#e8e8e8")
        img = ax.imshow(
            mat, cmap=cmap, aspect="auto",
            vmin=0, vmax=100, interpolation="nearest",
        )
        last_img = img

        # Tick positions on every panel; labels only on outer edges
        ax.set_xticks(np.arange(len(size_order)))
        ax.set_yticks(np.arange(len(quintiles)))
        if row == 1:
            ax.set_xticklabels(slabels, fontsize=7, rotation=20, ha="right")
        else:
            ax.set_xticklabels([])
        if col == 0:
            ax.set_yticklabels(qlabels, fontsize=6.5)
        else:
            ax.set_yticklabels([])
        ax.tick_params(length=0)

        for y in np.arange(len(quintiles) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.6)
        for x in np.arange(len(size_order) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.6)

        for r in range(len(quintiles)):
            for c in range(len(size_order)):
                val = mat[r, c]
                if not np.isfinite(val):
                    ax.text(c, r, "—", ha="center", va="center",
                            fontsize=7, color="#888888")
                    continue
                text_col = "white" if val > 60 else "#333333"
                ax.text(c, r, f"{val:.0f}%",
                        ha="center", va="center", fontsize=6.5, color=text_col)

        ax.set_title(f'{mix_label}\n({cbar_label})', pad=4)

    style.add_panel_labels(axes.ravel(), x=-0.12, y=1.10, size=9)
    fig.subplots_adjust(
        left=0.13, right=0.85, top=0.92, bottom=0.12, hspace=0.44, wspace=0.08,
    )

    # Single shared colourbar spanning the full subplot-area height
    top_pos = axes[0, 1].get_position()
    bot_pos = axes[1, 1].get_position()
    cbar_ax = fig.add_axes([
        top_pos.x1 + 0.025,
        bot_pos.y0,
        0.022,
        top_pos.y1 - bot_pos.y0,
    ])
    cb = fig.colorbar(last_img, cax=cbar_ax)
    cb.set_label(cbar_label, fontsize=7)
    cb.ax.tick_params(labelsize=6.5)

    save_all(style, fig, out_dir / out_stem, "double", 5.0)


def plot_cross_category_heatmap(
    style, cluster_demo: pd.DataFrame, out_dir: Path
) -> None:
    """Supp fig 3 — fraction 'more mix' by SIMD quintile × cluster size."""
    _cross_category_heatmap(
        style, cluster_demo, out_dir,
        cat_value="more mix",
        cmap="Reds",
        cbar_label='% "more mix"',
        out_stem="supp_fig3_cross_category_heatmap",
    )


def plot_cross_category_less_mix_heatmap(
    style, cluster_demo: pd.DataFrame, out_dir: Path
) -> None:
    """Supp fig 4 — fraction 'less mix' by SIMD quintile × cluster size."""
    _cross_category_heatmap(
        style, cluster_demo, out_dir,
        cat_value="less mix",
        cmap="Blues",
        cbar_label='% "less mix"',
        out_stem="supp_fig4_cross_category_less_mix_heatmap",
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def remove_stale_figures(out_dir: Path) -> None:
    """Delete output files whose stems are no longer produced by this script."""
    current_stems = {
        "fig1_vaccinated_cases_over_time",
        "fig2_cluster_vaccination_by_wave_and_category",
        "fig3_vaccination_mixing_by_wave",
        "fig4_demographic_mixing_by_wave",
        "fig5_geographic_dispersion_by_wave",
        "fig6_dose_recency_by_simd",
        "supp_fig1_weekly_mixing_evolution",
        "supp_fig2_domain_dose_gradient",
        "supp_fig3_cross_category_heatmap",
        "supp_fig4_cross_category_less_mix_heatmap",
    }
    for f in out_dir.iterdir():
        if f.suffix in {".pdf", ".png", ".tif"} and f.stem not in current_stems:
            try:
                f.unlink()
            except OSError:
                pass  # skip if filesystem does not allow deletion


def run(
    root: Path,
    tables_dir: Path | None = None,
    cache_dir: Path | None = None,
    out_dir: Path | None = None,
) -> None:
    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)

    setup_environment()
    style = load_style(root)

    tables_dir = tables_dir or root / "part2" / "tables"
    cache_dir  = cache_dir  or root / "part2" / "cache"
    out_dir    = out_dir    or root / "part2" / "manuscript" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    remove_stale_figures(out_dir)

    # --- CSV tables (always available) ---
    case_weekly = pd.read_csv(
        tables_dir / "vaccination_case_weekly_summary.csv",
        parse_dates=["case_week"],
    )
    cluster_wave_cat = pd.read_csv(
        tables_dir / "vaccination_cluster_wave_category_summary.csv",
    )
    cluster_wave_simd = pd.read_csv(
        tables_dir / "vaccination_cluster_wave_simd_domain_summary.csv",
    )
    cluster_weekly_cat = pd.read_csv(
        tables_dir / "vaccination_cluster_weekly_category_summary.csv",
        parse_dates=["cluster_week"],
    )

    # --- Parquet cache (requires pyarrow; gracefully skipped if absent) ---
    cluster_demo = load_cluster_demo_mix(cache_dir)

    print("Generating fig1 — vaccinated cases over time …")
    plot_vaccinated_cases_over_time(style, case_weekly, out_dir)

    print("Generating fig2 — cluster vaccination by wave and category …")
    plot_cluster_vaccination_by_wave_and_category(style, cluster_wave_cat, out_dir)

    print("Generating fig3 — vaccination mixing by wave …")
    plot_vaccination_mixing_by_wave(style, cluster_wave_cat, out_dir)

    if cluster_demo is not None:
        print("Generating fig4 — demographic mixing by wave …")
        plot_demographic_mixing_by_wave(style, cluster_demo, out_dir)
    else:
        print("Skipping fig4 — demographic mixing (parquets unavailable).")

    print("Generating fig5 — geographic dispersion by wave …")
    plot_geographic_dispersion_by_wave(style, cluster_wave_cat, out_dir)

    print("Generating fig6 — dose recency by SIMD …")
    plot_dose_recency_by_simd(style, cluster_wave_simd, out_dir)

    print("Generating supp_fig1 — weekly mixing evolution …")
    plot_weekly_mixing_evolution(style, cluster_weekly_cat, out_dir)

    print("Generating supp_fig2 — domain dose gradient …")
    plot_domain_dose_gradient(style, cluster_wave_simd, out_dir)

    if cluster_demo is not None:
        print("Generating supp_fig3 — cross-category heatmap (more mix) …")
        plot_cross_category_heatmap(style, cluster_demo, out_dir)
        print("Generating supp_fig4 — cross-category heatmap (less mix) …")
        plot_cross_category_less_mix_heatmap(style, cluster_demo, out_dir)
    else:
        print("Skipping supp_fig3/4 — cross-category heatmaps (parquets unavailable).")

    print(f"\nDone. Figures written to:\n  {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument(
        "--tables-dir", type=Path, default=None,
        help="Part 2 tables directory. Defaults to part2/tables.",
    )
    parser.add_argument(
        "--cache-dir", type=Path, default=None,
        help="Part 2 cache directory (parquets). Defaults to part2/cache.",
    )
    parser.add_argument(
        "--out-dir", type=Path, default=None,
        help="Figure output directory. Defaults to part2/manuscript/figures.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = (args.root or repo_root()).resolve()
    run(
        root,
        tables_dir=args.tables_dir.resolve() if args.tables_dir else None,
        cache_dir=args.cache_dir.resolve() if args.cache_dir else None,
        out_dir=args.out_dir.resolve() if args.out_dir else None,
    )


if __name__ == "__main__":
    main()
