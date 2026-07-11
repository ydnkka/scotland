"""Shared compatibility-assortativity time-series panels."""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Literal

from matplotlib.axes import Axes
import matplotlib.dates as mdates
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (  # noqa: E402
    Paths,
    add_policy_bands,
    date_axis,
    panel_label,
    read_table,
    styled_new_figure,
    styled_save_figure,
    window_idx_from_id,
)


ASSORTATIVITY_ATTRIBUTES: tuple[tuple[str, str], ...] = (
    ("age_band", "Age band"),
    ("age_group", "Age group"),
    ("simd_quintile", "SIMD quintile"),
    ("health_board", "Health board"),
    ("urban_rural", "Urban/rural class"),
    ("local_authority", "Local authority"),
)
ASSORTATIVITY_ATTRIBUTE_ORDER = [attribute for attribute, _ in ASSORTATIVITY_ATTRIBUTES]
ASSORTATIVITY_ATTRIBUTE_LABELS = dict(ASSORTATIVITY_ATTRIBUTES)

MIN_EDGE_CONTRIBUTIONS = 20
LINE_COLOR = "#1f4e79"
RIBBON_COLOR = "#4e79a7"


def _weighted_mean_ci_from_se(
    values: pd.Series,
    weights: pd.Series,
    standard_errors: pd.Series,
) -> dict[str, float]:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return {
            "weighted_mean": np.nan,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "ci_weight_share": np.nan,
        }

    values = values.loc[mask].astype(float)
    weights = weights.loc[mask].astype(float)
    weighted_mean = float(np.average(values, weights=weights))

    se_mask = standard_errors.loc[mask].notna()
    if not se_mask.any():
        return {
            "weighted_mean": weighted_mean,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "ci_weight_share": np.nan,
        }

    ci_weights = weights.loc[se_mask]
    ci_standard_errors = standard_errors.loc[mask].loc[se_mask].astype(float)
    normalized = ci_weights / weights.sum()
    combined_se = float(np.sqrt(np.sum((normalized * ci_standard_errors) ** 2)))
    return {
        "weighted_mean": weighted_mean,
        "combined_se": combined_se,
        "ci_low": weighted_mean - 1.96 * combined_se,
        "ci_high": weighted_mean + 1.96 * combined_se,
        "ci_weight_share": float(ci_weights.sum() / weights.sum()),
    }


def _weighted_quantile(
    values: pd.Series,
    weights: pd.Series,
    quantile: float,
) -> float:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return np.nan

    value_array = values.loc[mask].astype(float).to_numpy()
    weight_array = weights.loc[mask].astype(float).to_numpy()
    order = np.argsort(value_array)
    value_array = value_array[order]
    weight_array = weight_array[order]
    midpoint = np.cumsum(weight_array) - 0.5 * weight_array
    midpoint = midpoint / weight_array.sum()
    return float(
        np.interp(
            quantile,
            midpoint,
            value_array,
            left=value_array[0],
            right=value_array[-1],
        )
    )


def compatibility_assortativity_filtered(paths: Paths) -> pd.DataFrame:
    assort = read_table(paths, "compatibility_assortativity")
    assort["window_idx"] = window_idx_from_id(assort["window_id"])
    return assort.loc[
        assort["attribute"].isin(ASSORTATIVITY_ATTRIBUTE_ORDER)
        & assort["assortativity"].notna()
        & assort["edge_weight_total"].gt(0)
        & assort["n_categories"].gt(1)
        & assort["n_edge_contributions_used"].ge(MIN_EDGE_CONTRIBUTIONS)
    ].copy()


def compatibility_window_assortativity(paths: Paths) -> pd.DataFrame:
    work = compatibility_assortativity_filtered(paths)
    rows = []
    for (window_idx, attribute, label), group in work.groupby(
        ["window_idx", "attribute", "attribute_label"], dropna=False
    ):
        ci = _weighted_mean_ci_from_se(
            group["assortativity"],
            group["edge_weight_total"],
            group["assortativity_se"],
        )
        rows.append(
            {
                "window_idx": window_idx,
                "attribute": attribute,
                "attribute_label": label,
                "assortativity": ci["weighted_mean"],
                "assortativity_se": ci["combined_se"],
                "assortativity_ci_low": ci["ci_low"],
                "assortativity_ci_high": ci["ci_high"],
                "ci_weight_share": ci["ci_weight_share"],
                "assortativity_q25": _weighted_quantile(
                    group["assortativity"], group["edge_weight_total"], 0.25
                ),
                "assortativity_q75": _weighted_quantile(
                    group["assortativity"], group["edge_weight_total"], 0.75
                ),
                "edge_weight_total": group["edge_weight_total"].sum(),
                "eligible_networks": group["pairwise_stem"].nunique()
                if "pairwise_stem" in group.columns
                else len(group),
            }
        )
    return pd.DataFrame(rows)


def window_coverage_for_assortativity(paths: Paths) -> pd.DataFrame:
    window = read_table(paths, "window_coverage").copy()
    window["wn_mid_date"] = pd.to_datetime(window["wn_mid_date"], errors="coerce")
    return window.sort_values("window_idx")


def compatibility_assortativity_timeseries(
    paths: Paths,
    *,
    uncertainty: Literal["iqr", "ci"],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = compatibility_window_assortativity(paths)
    if summary.empty:
        return summary, window_coverage_for_assortativity(paths)

    if uncertainty == "iqr":
        summary["uncertainty_low"] = summary["assortativity_q25"]
        summary["uncertainty_high"] = summary["assortativity_q75"]
    elif uncertainty == "ci":
        summary["uncertainty_low"] = summary["assortativity_ci_low"]
        summary["uncertainty_high"] = summary["assortativity_ci_high"]
    else:
        raise ValueError("uncertainty must be 'iqr' or 'ci'")

    window = window_coverage_for_assortativity(paths)
    summary = summary.merge(
        window[["window_idx", "wn_mid_date", "policy_period"]],
        on="window_idx",
        how="left",
    )
    summary["attribute"] = pd.Categorical(
        summary["attribute"],
        categories=ASSORTATIVITY_ATTRIBUTE_ORDER,
        ordered=True,
    )
    return summary.sort_values(["attribute", "window_idx"]), window


def _date_values(values: pd.Series) -> np.ndarray:
    dates = pd.to_datetime(values, errors="coerce")
    out = np.full(len(dates), np.nan, dtype=float)
    valid = dates.notna().to_numpy()
    if valid.any():
        out[valid] = mdates.date2num(dates.loc[valid].dt.to_pydatetime())
    return out


def _numeric_values(values: pd.Series) -> np.ndarray:
    return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)


def _panel_table(
    summary: pd.DataFrame, window: pd.DataFrame, attribute: str
) -> pd.DataFrame:
    base = window[["window_idx", "wn_mid_date"]].copy()
    panel = summary.loc[
        summary["attribute"].eq(attribute),
        [
            "window_idx",
            "assortativity",
            "uncertainty_low",
            "uncertainty_high",
        ],
    ]
    return base.merge(panel, on="window_idx", how="left").sort_values("window_idx")


def _plot_attribute_panel(
    ax: Axes,
    panel: pd.DataFrame,
    window: pd.DataFrame,
    *,
    label: str,
) -> None:
    add_policy_bands(ax, window)
    x = _date_values(panel["wn_mid_date"])
    y = _numeric_values(panel["assortativity"])
    lower = _numeric_values(panel["uncertainty_low"])
    upper = _numeric_values(panel["uncertainty_high"])
    ribbon_mask = np.isfinite(x) & np.isfinite(lower) & np.isfinite(upper)

    ax.fill_between(
        x,
        lower,
        upper,
        where=ribbon_mask,
        color=RIBBON_COLOR,
        alpha=0.22,
        linewidth=0,
        interpolate=False,
        zorder=1,
    )
    ax.plot(x, y, color=LINE_COLOR, lw=1.35, zorder=2)
    ax.axhline(0, color="#777777", lw=0.65, ls=":", zorder=1.5)
    ax.xaxis_date()
    ax.set_title(label)
    date_axis(ax)


def _shared_y_limits(summary: pd.DataFrame) -> tuple[float, float]:
    values = pd.concat(
        [
            summary["assortativity"],
            summary["uncertainty_low"],
            summary["uncertainty_high"],
        ],
        ignore_index=True,
    )
    values = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan)
    values = values.dropna()
    if values.empty:
        return -0.1, 1.0
    y_min = min(float(values.min()), 0.0)
    y_max = max(float(values.max()), 0.0)
    padding = max((y_max - y_min) * 0.08, 0.05)
    return y_min - padding, y_max + padding


def plot_compatibility_assortativity_grid(
    paths: Paths,
    *,
    figure_name: str,
    uncertainty: Literal["iqr", "ci"],
) -> None:
    summary, window = compatibility_assortativity_timeseries(
        paths, uncertainty=uncertainty
    )

    fig, axes = styled_new_figure(
        width="double",
        height_in=7.4,
        nrows=3,
        ncols=2,
        sharex=True,
        sharey=True,
    )
    axes_array = np.asarray(axes).reshape(3, 2)
    y_limits = _shared_y_limits(summary)

    for idx, (attribute, label) in enumerate(ASSORTATIVITY_ATTRIBUTES):
        ax = axes_array.flat[idx]
        panel = _panel_table(summary, window, attribute)
        _plot_attribute_panel(
            ax,
            panel,
            window,
            label=label,
        )
        ax.set_ylim(*y_limits)
        panel_label(ax, chr(ord("A") + idx))
        if idx % 2 == 0:
            ax.set_ylabel("Compatibility assortativity")
        if idx < 4:
            ax.tick_params(labelbottom=False)

    fig.subplots_adjust(left=0.08, right=0.985, top=0.94, bottom=0.09, hspace=0.38)
    styled_save_figure(fig, paths, figure_name, tight=False)
