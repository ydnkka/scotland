"""Age-conditional vaccination-mixing figures for SSE candidates."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any

from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, NullFormatter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils import WIDTHS, CONTEXTS, new_figure, add_panel_labels

from .association_pipeline import VACCINATION_MIXING_TERTILE_ORDER
from .palettes import (
    BORDER,
    CANDIDATE_COLOR,
    CANDIDATE_DARK,
    GRAY,
    GRAY_LIGHT,
    GRID,
    INK,
    ORANGE_DARK,
    SSE_CATEGORY_ORDER,
    SSE_CATEGORY_PALETTE,
    TEAL_DARK,
)

__all__ = [
    "make_vaccination_figures",
    "plot_vaccination_report",
]


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = (
    PROJECT_ROOT / "sse_detection" / "results" / "vaccination_outputs"
)
REFERENCE_COLOR = "#FF0000"

TERTILE_LABELS = {
    "more_homogeneous": "More\nhomogeneous",
    "as_expected": "As\nexpected",
    "more_mixed": "More\nmixed",
}

SSE_CATEGORY_LABELS = {
    "mixed_population_dissemination": "Mixed-population",
    "putative_introduction_burst": "Introduction burst",
    "secondary_relay_amplification": "Secondary relay",
    "diffuse_branching_transmission": "Diffuse branching",
    "focused_branching_transmission": "Focused branching",
    "sustained_single_chain": "Single chain",
    "contained_local_burst": "Contained burst",
    "high_volume_onward_transmission": "High-volume onward",
    "ambiguous_amplification_signal": "Ambiguous",
}


def _read_table(table: pd.DataFrame | str | Path | Any) -> pd.DataFrame:
    if isinstance(table, pd.DataFrame):
        return table.copy()
    path = Path(table)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, skipinitialspace=True)


def _clean_strings(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(col).strip() for col in out.columns]
    for col in out.select_dtypes(include=["object", "string", "category"]).columns:
        present = out[col].notna()
        out.loc[present, col] = out.loc[present, col].astype(str).str.strip()
    return out


def _pretty_text(value: Any, label_map: dict[str, str] | None = None) -> str:
    if pd.isna(value):
        return ""
    text = str(value)
    if label_map and text in label_map:
        return label_map[text]
    return text.replace("_", " ").strip().capitalize()


def _term_level(term: Any) -> str:
    if pd.isna(term):
        return ""
    match = re.search(r"\[T\.(.*)\]$", str(term))
    return match.group(1) if match else str(term)


def _prepare_node_features(node_features: pd.DataFrame) -> pd.DataFrame:
    df = _clean_strings(node_features)
    df["candidate"] = pd.to_numeric(df["candidate"], errors="coerce").fillna(0).astype(int)
    df["vaccination_mix_entropy_z"] = pd.to_numeric(
        df["vaccination_mix_entropy_z"],
        errors="coerce",
    )
    df["candidate_label"] = np.where(df["candidate"].eq(1), "Candidate", "Background")
    return df.dropna(subset=["vaccination_mix_entropy_z"])


def _prepare_summary(summary: pd.DataFrame) -> pd.DataFrame:
    df = _clean_strings(summary)
    df["vaccination_mix_tertile"] = pd.Categorical(
        df["vaccination_mix_tertile"],
        categories=VACCINATION_MIXING_TERTILE_ORDER,
        ordered=True,
    )
    if "candidate_rate" not in df.columns:
        df["candidate_rate"] = df["n_candidates"] / df["n_nodes"]
    df["tertile_label"] = df["vaccination_mix_tertile"].map(TERTILE_LABELS).astype(str)
    return df.sort_values("vaccination_mix_tertile")


def _prepare_category_mix(category_summary: pd.DataFrame) -> pd.DataFrame:
    df = _clean_strings(category_summary)
    if "candidate_category_share" not in df.columns:
        totals = df.groupby("vaccination_mix_tertile", dropna=False)[
            "n_candidates"
        ].transform("sum")
        df["candidate_category_share"] = df["n_candidates"] / totals
    categories = [
        category
        for category in SSE_CATEGORY_ORDER
        if category != "not_sse_like" and category in set(df["sse_category"])
    ]
    extras = sorted(set(df["sse_category"]) - set(categories))
    categories.extend(extras)
    pivot = df.pivot_table(
        index="vaccination_mix_tertile",
        columns="sse_category",
        values="candidate_category_share",
        aggfunc="sum",
        fill_value=0,
    )
    pivot = pivot.reindex(
        index=VACCINATION_MIXING_TERTILE_ORDER,
        columns=categories,
        fill_value=0,
    )
    return pivot.loc[pivot.sum(axis=1).gt(0)]


def _prepare_or_forest(odds_ratios: pd.DataFrame) -> pd.DataFrame:
    df = _clean_strings(odds_ratios)
    if "model_set" in df.columns:
        df = df.loc[df["model_set"].astype(str).eq("primary")].copy()
    if "predictor_set" in df.columns:
        df = df.loc[df["predictor_set"].astype(str).eq("single")].copy()

    rows: list[dict[str, Any]] = []
    continuous = df.loc[df["predictor"].astype(str).eq("vaccination_mix_entropy_z")]
    if not continuous.empty:
        row = continuous.iloc[0]
        rows.append(
            {
                "group": "Continuous",
                "label": "Mixing z-score",
                "odds_ratio": row.odds_ratio,
                "or_low": row.or_low,
                "or_high": row.or_high,
                "is_reference": False,
                "order": 0,
            }
        )

    tertile = df.loc[df["predictor"].astype(str).eq("vaccination_mix_tertile")]
    tertile_lookup = {
        _term_level(row.term): row
        for row in tertile.itertuples(index=False)
    }
    for order, level in enumerate(VACCINATION_MIXING_TERTILE_ORDER, start=1):
        is_reference = level == "as_expected"
        if is_reference:
            rows.append(
                {
                    "group": "Tertile",
                    "label": f"{_pretty_text(level, TERTILE_LABELS)} (ref)",
                    "odds_ratio": 1.0,
                    "or_low": 1.0,
                    "or_high": 1.0,
                    "is_reference": True,
                    "order": order,
                }
            )
            continue
        if level not in tertile_lookup:
            continue
        row = tertile_lookup[level]
        rows.append(
            {
                "group": "Tertile",
                "label": _pretty_text(level, TERTILE_LABELS),
                "odds_ratio": row.odds_ratio,
                "or_low": row.or_low,
                "or_high": row.or_high,
                "is_reference": False,
                "order": order,
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No vaccination-mixing odds-ratio rows were available.")
    for col in ["odds_ratio", "or_low", "or_high"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return (
        out.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["odds_ratio"])
        .sort_values("order")
    )


def _forest_xlim(panel: pd.DataFrame) -> tuple[float, float]:
    values = pd.to_numeric(
        pd.concat([panel["or_low"], panel["or_high"], pd.Series([1.0])]),
        errors="coerce",
    )
    values = values[np.isfinite(values) & values.gt(0)]
    if values.empty:
        return (0.75, 1.35)
    lo = float(values.min())
    hi = float(values.max())
    pad = max((np.log(hi) - np.log(lo)) * 0.12, 0.08)
    return float(np.exp(np.log(lo) - pad)), float(np.exp(np.log(hi) + pad))


def _plot_z_distribution_panel(ax: Any, node_features: pd.DataFrame) -> None:
    labels = ["Background", "Candidate"]
    data = [
        node_features.loc[
            node_features["candidate_label"].eq(label),
            "vaccination_mix_entropy_z",
        ].dropna()
        for label in labels
    ]
    box = ax.boxplot(
        data,
        vert=False,
        labels=labels,
        patch_artist=True,
        widths=0.55,
        showfliers=False,
    )
    for patch, color in zip(box["boxes"], [GRAY_LIGHT, CANDIDATE_COLOR]):
        patch.set_facecolor(color)
        patch.set_edgecolor(INK)
        patch.set_alpha(0.85)
    for key in ["whiskers", "caps", "medians"]:
        for artist in box[key]:
            artist.set_color(INK)
    ax.axvline(0, color=GRAY, lw=1.0, ls="--", zorder=0)
    ax.set_xlabel("Vaccination-mixing z-score")
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def _plot_candidate_rate_panel(ax: Any, summary: pd.DataFrame) -> None:
    y = np.arange(len(summary))
    ax.barh(
        y,
        summary["candidate_rate"],
        color=CANDIDATE_COLOR,
        edgecolor=CANDIDATE_DARK,
        linewidth=0.8,
        alpha=0.9,
    )
    for ypos, row in zip(y, summary.itertuples(index=False)):
        ax.annotate(
            f"{int(row.n_candidates):,}/{int(row.n_nodes):,}",
            (row.candidate_rate, ypos),
            xytext=(5, 0),
            textcoords="offset points",
            va="center",
            ha="left",
            fontsize="small",
            color=INK,
        )
    xmax = max(0.05, float(summary["candidate_rate"].max())) * 1.23
    ax.set_xlim(0, min(xmax, 1.0))
    ax.set_yticks(y)
    ax.set_yticklabels(summary["tertile_label"])
    ax.invert_yaxis()
    ax.set_xlabel("Candidate-node rate")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{100 * value:.0f}%"))
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def _plot_category_mix_panel(ax: Any, category_mix: pd.DataFrame) -> None:
    y = np.arange(len(category_mix))
    left = np.zeros(len(category_mix), dtype=float)
    handles = []
    labels = []
    for category in category_mix.columns:
        values = category_mix[category].to_numpy(dtype=float)
        if not np.any(values > 0):
            continue
        color = SSE_CATEGORY_PALETTE.get(category, GRAY)
        handle = ax.barh(
            y,
            values,
            left=left,
            color=color,
            edgecolor="white",
            linewidth=0.5,
            alpha=0.95,
        )
        handles.append(handle[0])
        labels.append(_pretty_text(category, SSE_CATEGORY_LABELS))
        left += values
    ax.set_xlim(0, 1)
    ax.set_yticks(y)
    ax.set_yticklabels([TERTILE_LABELS.get(str(idx), str(idx)) for idx in category_mix.index])
    ax.invert_yaxis()
    ax.set_xlabel("Candidate category share")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{100 * value:.0f}%"))
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    if handles:
        ax.legend(
            handles,
            labels,
            loc="lower left",
            bbox_to_anchor=(0, -0.60),
            ncol=2,
            frameon=False,
            fontsize="small",
            columnspacing=1.0,
            handlelength=1.1,
        )


def _plot_or_forest_panel(ax: Any, panel: pd.DataFrame) -> None:
    panel = panel.copy().reset_index(drop=True)
    y = np.arange(len(panel))[::-1]
    group_boundaries = []
    previous_group = None
    for idx, row in enumerate(panel.itertuples(index=False)):
        ypos = y[idx]
        low = row.or_low if pd.notna(row.or_low) else row.odds_ratio
        high = row.or_high if pd.notna(row.or_high) else row.odds_ratio
        overlaps_one = min(low, high) <= 1 <= max(low, high)
        if row.is_reference:
            color = REFERENCE_COLOR
        elif overlaps_one:
            color = GRAY
        else:
            color = ORANGE_DARK if row.odds_ratio > 1 else TEAL_DARK
        if row.is_reference:
            ax.plot(1.0, ypos, marker="o", ms=6.5, color=color, zorder=3)
        else:
            ax.plot([low, high], [ypos, ypos], color=color, lw=2.1, alpha=0.45)
            ax.plot(row.odds_ratio, ypos, marker="o", ms=6.5, color=color, zorder=3)
        ax.annotate(
            f"{row.odds_ratio:.2f}",
            (row.odds_ratio, ypos),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize="small",
            color=color,
            fontweight="bold" if not row.is_reference and not overlaps_one else "normal",
        )
        if previous_group is not None and row.group != previous_group:
            group_boundaries.append(ypos + 0.5)
        previous_group = row.group
    for boundary in group_boundaries:
        ax.axhline(boundary, color=BORDER, lw=0.8)
    ax.axvline(1.0, color=GRAY_LIGHT, lw=1.0, ls="--", zorder=0)
    ax.set_xscale("log")
    xlim = _forest_xlim(panel)
    ax.set_xlim(*xlim)
    ticks = [
        tick
        for tick in [0.5, 0.6, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
        if xlim[0] <= tick <= xlim[1]
    ]
    if ticks:
        ax.set_xticks(ticks)
    ax.set_yticks(y)
    ax.set_yticklabels(panel["label"].tolist())
    ax.set_ylim(-0.6, len(panel) - 0.05)
    ax.set_xlabel("Adjusted candidate odds ratio")
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())


def plot_vaccination_report(
    node_features: pd.DataFrame | str | Path | Any,
    mixing_summary: pd.DataFrame | str | Path | Any,
    category_summary: pd.DataFrame | str | Path | Any,
    odds_ratios: pd.DataFrame | str | Path | Any,
    *,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 6.3,
    width_ratios: tuple[float, float] | None = None,
    height_ratios: tuple[float, float] | None = None,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Plot age-conditional vaccination-mixing context for SSE candidates."""
    node = _prepare_node_features(_read_table(node_features))
    summary = _prepare_summary(_read_table(mixing_summary))
    category_mix = _prepare_category_mix(_read_table(category_summary))
    forest = _prepare_or_forest(_read_table(odds_ratios))

    fig, axes = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        nrows=2,
        ncols=2,
        context=context,
        font_scale=font_scale,
        gridspec_kw={
            "width_ratios": width_ratios,
            "height_ratios": height_ratios,
            "wspace": 0.10,
        },
        layout="constrained",
    )
    ax_dist = axes[0, 0]
    ax_rate = axes[0, 1]
    ax_mix = axes[1, 0]
    ax_forest = axes[1, 1]

    _plot_z_distribution_panel(ax_dist, node)
    _plot_candidate_rate_panel(ax_rate, summary)
    _plot_category_mix_panel(ax_mix, category_mix)
    _plot_or_forest_panel(ax_forest, forest)

    ax_dist.set_title("Mixing against age-window null")
    ax_rate.set_title("Candidate frequency")
    ax_mix.set_title("Candidate phenotype mix")
    ax_forest.set_title("Adjusted ORs")
    add_panel_labels(axes.ravel())
    plt.close(fig)
    return fig


def make_vaccination_figures(
    *,
    results_dir: Path | str = DEFAULT_RESULTS_DIR,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 6.3,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> dict[str, Figure]:
    """Build vaccination-mixing figures from the exported analysis CSVs."""
    results_path = Path(results_dir)
    fig = plot_vaccination_report(
        results_path / "vaccination_mixing_age_conditional_node_features.csv",
        results_path / "vaccination_mixing_age_conditional_summary.csv",
        results_path / "vaccination_mixing_age_conditional_category_summary.csv",
        results_path / "vaccination_mixing_age_conditional_odds_ratios.csv",
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
    )
    return {"vaccination_report": fig}
