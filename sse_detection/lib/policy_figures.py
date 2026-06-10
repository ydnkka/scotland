"""Compact policy-era context figures for SSE candidates."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping

from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, NullFormatter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils import WIDTHS, CONTEXTS, new_figure, add_panel_labels

from .association_pipeline import POLICY_ERA_ORDER
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
from .table_utils import (
    clean_strings as _clean_strings,
    forest_xlim as _forest_xlim,
    pretty_text as _pretty_text,
    read_table as _read_table,
    term_level as _term_level,
)

__all__ = [
    "make_policy_figures",
    "plot_policy_report",
]


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = (
    PROJECT_ROOT / "sse_detection" / "results" / "policy_outputs"
)

POLICY_ERA_LABELS = {
    "early_restriction_easing": "Early easing",
    "autumn_winter_restrictions": "Autumn/winter\nrestrictions",
    "spring_summer_2021_easing": "Spring/summer\n2021 easing",
    "near_normal_delta": "Near-normal\nDelta",
    "omicron_response": "Omicron\nresponse",
    "post_restriction": "Post-\nrestriction",
}

SSE_CATEGORY_LABELS = {
    "burst": "Burst",
    "burden": "Burden",
    "burst+burden": "Burst+burden",
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

REFERENCE_COLOR = "#FF0000"


def _filter_model_single(
    df: pd.DataFrame,
    *,
    model_set: str = "primary",
    predictor_set: str = "single",
) -> pd.DataFrame:
    out = df.copy()
    if "model_set" in out.columns:
        out = out.loc[out["model_set"].astype(str).eq(model_set)].copy()
    if "predictor_set" in out.columns:
        out = out.loc[out["predictor_set"].astype(str).eq(predictor_set)].copy()
    return out


def _prepare_policy_summary(policy_summary: pd.DataFrame) -> pd.DataFrame:
    df = _clean_strings(policy_summary)
    if "candidate_rate" not in df.columns:
        df["candidate_rate"] = df["n_candidates"] / df["n_nodes"]
    df["policy_era"] = pd.Categorical(
        df["policy_era"],
        categories=POLICY_ERA_ORDER,
        ordered=True,
    )
    df = df.sort_values("policy_era")
    df["policy_label"] = df["policy_era"].map(POLICY_ERA_LABELS).astype(str)
    return df


def _prepare_category_mix(category_summary: pd.DataFrame) -> pd.DataFrame:
    df = _clean_strings(category_summary)
    category_col = "sse_signature" if "sse_signature" in df.columns else "sse_category"
    share_col = (
        "candidate_signature_share"
        if "candidate_signature_share" in df.columns
        else "candidate_category_share"
    )
    if category_col not in df.columns:
        raise ValueError(
            "Policy category/signature summary needs `sse_signature` or "
            "`sse_category`."
        )
    if share_col not in df.columns:
        totals = df.groupby("policy_era", dropna=False)["n_candidates"].transform("sum")
        df[share_col] = df["n_candidates"] / totals
    categories = [
        category
        for category in SSE_CATEGORY_ORDER
        if category not in {"none", "not_sse_like"} and category in set(df[category_col])
    ]
    extras = sorted(set(df[category_col]) - set(categories) - {"none", "not_sse_like"})
    categories.extend(extras)

    pivot = df.pivot_table(
        index="policy_era",
        columns=category_col,
        values=share_col,
        aggfunc="sum",
        fill_value=0,
    )
    pivot = pivot.reindex(index=POLICY_ERA_ORDER, columns=categories, fill_value=0)
    pivot = pivot.loc[pivot.sum(axis=1).gt(0)]
    return pivot


def _policy_era_or_rows(
    table: pd.DataFrame,
    *,
    model_set: str = "primary",
    label_map: Mapping[str, str] | None = None,
    include_reference: bool = False,
) -> list[dict[str, Any]]:
    df = _clean_strings(table)
    if "predictor" not in df.columns:
        return []
    df = df.loc[df["predictor"].astype(str).eq("policy_era")].copy()
    df = _filter_model_single(df, model_set=model_set)
    if df.empty:
        return []
    reference = "post_restriction"
    if "reference" in df.columns:
        observed_references = [
            str(value).strip()
            for value in df["reference"]
            if pd.notna(value) and str(value).strip()
        ]
        if observed_references:
            reference = observed_references[0]

    rows = []
    if include_reference:
        reference_label = _pretty_text(reference, label_map or POLICY_ERA_LABELS)
        rows.append(
            {
                "group": "Policy era",
                "label": f"{reference_label} (ref)",
                "policy_era": reference,
                "odds_ratio": 1.0,
                "or_low": 1.0,
                "or_high": 1.0,
                "is_reference": True,
            }
        )
    for row in df.itertuples(index=False):
        level = _term_level(getattr(row, "term"))
        rows.append(
            {
                "group": "Policy era",
                "label": _pretty_text(level, label_map or POLICY_ERA_LABELS),
                "policy_era": level,
                "odds_ratio": getattr(row, "odds_ratio"),
                "or_low": getattr(row, "or_low"),
                "or_high": getattr(row, "or_high"),
                "is_reference": False,
            }
        )
    order = {era: idx for idx, era in enumerate(POLICY_ERA_ORDER)}
    return sorted(rows, key=lambda row: order.get(str(row["policy_era"]), 999))


def _prepare_policy_era_forest(
    policy_odds_ratios: pd.DataFrame,
    *,
    model_set: str = "primary",
) -> pd.DataFrame:
    rows = _policy_era_or_rows(
        policy_odds_ratios,
        model_set=model_set,
        include_reference=True,
    )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError(
            f"No {model_set} single policy-era odds-ratio rows were available."
        )
    for col in ["odds_ratio", "or_low", "or_high"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out.replace([np.inf, -np.inf], np.nan).dropna(subset=["odds_ratio"])


def _plot_candidate_rate_panel(
    ax: Any,
    summary: pd.DataFrame,
    *,
    show_y_labels: bool = True,
    invert_y_axis: bool = True,
) -> None:
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
            f"{int(row.n_candidates):,}/{int(row.n_nodes):,}", # type: ignore
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
    if show_y_labels:
        ax.set_yticklabels(summary["policy_label"])
    if invert_y_axis:
        ax.invert_yaxis()
    ax.set_xlabel("Candidate-node rate")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{100 * value:.0f}%"))
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelleft=show_y_labels)


def _plot_category_mix_panel(
    ax: Any,
    category_mix: pd.DataFrame,
    policy_labels: Iterable[str],
    *,
    show_y_labels: bool = True,
    invert_y_axis: bool = True,
) -> None:
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
    if show_y_labels:
        ax.set_yticklabels(list(policy_labels))
    if invert_y_axis:
        ax.invert_yaxis()
    ax.set_xlabel("Candidate category share")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{100 * value:.0f}%"))
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelleft=show_y_labels)
    if handles:
        ax.legend(
            handles,
            labels,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            ncol=1,
            frameon=False,
            fontsize="small",
            handlelength=1.1,
            borderaxespad=0.0,
        )


def _or_tick_label(value: float) -> str:
    if value < 1:
        return f"{value:.2g}"
    return f"{value:g}"


def _plot_or_forest_panel(
    ax: Any,
    panel: pd.DataFrame,
    *,
    show_y_labels: bool = True,
) -> None:
    panel = panel.copy().reset_index(drop=True)
    if "is_reference" not in panel.columns:
        panel["is_reference"] = False
    else:
        panel["is_reference"] = panel["is_reference"].where(
            panel["is_reference"].notna(),
            False,
        ).astype(bool)
    y = np.arange(len(panel))[::-1]
    group_boundaries = []
    previous_group = None
    for idx, row in enumerate(panel.itertuples(index=False)):
        ypos = y[idx]
        low = row.or_low if pd.notna(row.or_low) else row.odds_ratio
        high = row.or_high if pd.notna(row.or_high) else row.odds_ratio
        overlaps_one = min(low, high) <= 1 <= max(low, high) # type: ignore
        is_reference = bool(row.is_reference)
        if is_reference:
            color = REFERENCE_COLOR
        elif overlaps_one:
            color = GRAY
        else:
            color = ORANGE_DARK if row.odds_ratio > 1 else TEAL_DARK # type: ignore
        if is_reference:
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
            fontweight="bold" if not is_reference and not overlaps_one else "normal",
        )
        if previous_group is not None and row.group != previous_group:
            group_boundaries.append(ypos + 0.5)
        previous_group = row.group

    for boundary in group_boundaries:
        ax.axhline(boundary, color=BORDER, lw=0.8)
    ax.axvline(1, color=GRAY_LIGHT, lw=1.0, ls="--", zorder=0)
    ax.set_xscale("log")
    xlim = _forest_xlim(panel)
    ax.set_xlim(*xlim)
    # ticks = [
    #     tick
    #     for tick in [-10, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
    #     if xlim[0] <= tick <= xlim[1]
    # ]

    ticks = [0.15, 0.5, 1.0, 2.0, 4.0, 8.0]
    
    if ticks:
        ax.set_xticks(ticks)
    ax.set_yticks(y)
    if show_y_labels:
        ax.set_yticklabels(panel["label"].tolist())
    ax.set_ylim(-0.6, len(panel) - 0.05)
    ax.set_xlabel("Adjusted candidate odds ratio")
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelleft=show_y_labels)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: _or_tick_label(value)))
    ax.xaxis.set_minor_formatter(NullFormatter())


def plot_policy_report(
    policy_summary: pd.DataFrame | str | Path | Any,
    policy_signature_summary: pd.DataFrame | str | Path | Any,
    policy_odds_ratios: pd.DataFrame | str | Path | Any,
    *,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 6.3,
    width_ratios: tuple[float, float] | None = None,
    height_ratios: tuple[float, float] | None = None,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Plot a compact policy-era report figure.

    The figure contains policy-era candidate rates, candidate-category mix by
    policy era, and primary plus expanded policy-era odds-ratio forests.
    """
    policy = _prepare_policy_summary(_read_table(policy_summary))
    category_mix = _prepare_category_mix(_read_table(policy_signature_summary))
    category_mix = category_mix.reindex(policy["policy_era"].astype(str), fill_value=0)
    policy_odds = _read_table(policy_odds_ratios)
    primary_policy_forest = _prepare_policy_era_forest(
        policy_odds,
        model_set="primary",
    )
    expanded_policy_forest = _prepare_policy_era_forest(
        policy_odds,
        model_set="expanded",
    )

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
    axes = axes
    ax_rate = axes[0, 0]
    ax_mix = axes[0, 1]
    ax_primary_policy_forest = axes[1, 0]
    ax_expanded_policy_forest = axes[1, 1]

    ax_mix.sharey(ax_rate)
    ax_expanded_policy_forest.sharey(ax_primary_policy_forest)
    ax_expanded_policy_forest.sharex(ax_primary_policy_forest)

    _plot_candidate_rate_panel(ax_rate, policy)
    _plot_category_mix_panel(
        ax_mix,
        category_mix,
        policy["policy_label"],
        show_y_labels=False,
        invert_y_axis=False,
    )
    _plot_or_forest_panel(ax_primary_policy_forest, primary_policy_forest)
    _plot_or_forest_panel(
        ax_expanded_policy_forest,
        expanded_policy_forest,
        show_y_labels=False,
    )

    ax_primary_policy_forest.set_xlabel("Primary OR vs post-restriction")
    ax_expanded_policy_forest.set_xlabel("Expanded OR vs post-restriction")
    add_panel_labels(axes.ravel())
    plt.close(fig)
    return fig


def make_policy_figures(
    *,
    results_dir: Path | str = DEFAULT_RESULTS_DIR,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 6.3,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> dict[str, Figure]:
    """Build policy-era figures from the exported analysis CSVs."""
    results_path = Path(results_dir)
    fig = plot_policy_report(
        results_path / "policy_era_candidate_summary.csv",
        results_path / "policy_era_signature_summary.csv",
        results_path / "policy_odds_ratios.csv",
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
    )
    return {"policy_report": fig}
