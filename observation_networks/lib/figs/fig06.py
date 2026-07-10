"""Build Chapter 4 Figure 6: vaccination context in the observed cohort."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

import matplotlib.dates as mdates
from matplotlib.axes import Axes
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    Paths,
    add_common_args,
    add_policy_bands,
    date_axis,
    new_blank_figure,
    panel_label,
    paths_from_args,
    read_table,
    save_figure,
    sort_by_policy,
)


FIGURE_NAME = "fig06_vaccination_context"
MIN_DAYS_SERIES_COUNT = 20

DOSE_GROUPS = (
    ("Unvaccinated", "n_unvaccinated", "prop_unvaccinated"),
    ("One dose", "n_one_dose", "prop_one_dose"),
    ("Two doses", "n_two_doses", "prop_two_doses"),
    (
        "Booster/3+ doses",
        "n_booster_or_three_plus",
        "prop_booster_or_three_plus",
    ),
    (
        "Dose unknown",
        "n_vaccinated_dose_unknown",
        "prop_vaccinated_dose_unknown",
    ),
)

DOSE_COLORS = {
    "Unvaccinated": "#747474",
    "One dose": "#59a14f",
    "Two doses": "#4e79a7",
    "Booster/3+ doses": "#b07aa1",
    "Dose unknown": "#bab0ac",
}


def _date_values(values: pd.Series) -> np.ndarray:
    dates = pd.to_datetime(values, errors="coerce")
    return mdates.date2num(dates.dt.to_pydatetime())


def _numeric(values: pd.Series) -> np.ndarray:
    return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)


def _active_dose_groups(table: pd.DataFrame, value_kind: str) -> list[tuple[str, str]]:
    if value_kind not in {"count", "proportion"}:
        raise ValueError("value_kind must be 'count' or 'proportion'")

    groups = []
    for label, count_col, prop_col in DOSE_GROUPS:
        col = count_col if value_kind == "count" else prop_col
        total = pd.to_numeric(table[col], errors="coerce").sum() if col in table else 0
        if total > 0:
            groups.append((label, col))
    return groups


def _plot_window_counts(ax: Axes, window: pd.DataFrame) -> None:
    add_policy_bands(ax, window)
    x = pd.to_datetime(window["wn_mid_date"], errors="coerce")
    groups = _active_dose_groups(window, "count")
    values = [_numeric(window[col]) for _, col in groups]
    colors = [DOSE_COLORS[label] for label, _ in groups]
    ax.stackplot(x, values, colors=colors, linewidth=0, alpha=0.94)
    ax.set_title("Rolling-window vaccination context")
    ax.set_ylabel("Sequences")
    date_axis(ax)
    panel_label(ax, "A")


def _plot_window_proportions(ax: Axes, window: pd.DataFrame) -> None:
    add_policy_bands(ax, window)
    x = pd.to_datetime(window["wn_mid_date"], errors="coerce")
    vaccinated = ax.plot(
        x,
        window["prop_vaccinated"],
        color="#1f4e79",
        lw=1.35,
        label="Vaccinated",
    )[0]
    booster = ax.plot(
        x,
        window["prop_booster"],
        color="#8e3b8a",
        lw=1.35,
        label="Booster recorded",
    )[0]
    ax.set_title("Vaccinated and booster share")
    ax.set_ylabel("Share of sequences")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.legend(
        [vaccinated, booster],
        ["Vaccinated", "Booster recorded"],
        loc="upper left",
    )
    date_axis(ax)
    panel_label(ax, "B")


def _plot_days_since_vaccination(ax: Axes, window: pd.DataFrame) -> None:
    add_policy_bands(ax, window)
    x = _date_values(window["wn_mid_date"])
    median = _numeric(window["median_days_since_vaccination"])
    lower = _numeric(window["q25_days_since_vaccination"])
    upper = _numeric(window["q75_days_since_vaccination"])
    enough_vaccinated = (
        pd.to_numeric(window["n_days_since_vaccination"], errors="coerce")
        .ge(MIN_DAYS_SERIES_COUNT)
        .to_numpy(dtype=bool)
    )
    ribbon_mask = (
        enough_vaccinated & np.isfinite(x) & np.isfinite(lower) & np.isfinite(upper)
    )
    line_mask = enough_vaccinated & np.isfinite(x) & np.isfinite(median)
    ax.fill_between(
        x[ribbon_mask],
        lower[ribbon_mask],
        upper[ribbon_mask],
        color="#f28e2b",
        alpha=0.24,
        linewidth=0,
    )
    ax.plot(x[line_mask], median[line_mask], color="#b85c00", lw=1.35)
    ax.xaxis_date()
    ax.set_title("Time since vaccination among vaccinated sequences")
    ax.set_ylabel("Days")
    date_axis(ax)
    panel_label(ax, "C")


def _plot_policy_composition(ax: Axes, policy: pd.DataFrame) -> None:
    policy = sort_by_policy(policy)
    labels = policy["policy_period"].astype(str).tolist()
    x = np.arange(len(policy))
    left = np.zeros(len(policy), dtype=float)
    groups = _active_dose_groups(policy, "proportion")
    for label, col in groups:
        values = _numeric(policy[col])
        ax.bar(
            x,
            values,
            bottom=left,
            width=0.78,
            color=DOSE_COLORS[label],
            edgecolor="white",
            linewidth=0.35,
        )
        left += np.nan_to_num(values, nan=0.0)
    ax.set_title("Dose composition by policy period")
    ax.set_ylabel("Share of sequences")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=0)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.grid(axis="y", color="#d9d9d9", lw=0.5, alpha=0.8)
    ax.set_axisbelow(True)
    panel_label(ax, "D")


def build(paths: Paths) -> None:
    window = read_table(paths, "vaccination_window_context")
    policy = read_table(paths, "vaccination_context_by_policy")
    for col in ("wn_start_date", "wn_mid_date", "wn_end_date"):
        if col in window.columns:
            window[col] = pd.to_datetime(window[col], errors="coerce")
    window = window.sort_values("window_idx")

    fig = new_blank_figure(
        "double",
        width_in=8.8,
        height_in=7.3,
    )
    grid = fig.add_gridspec(2, 2, hspace=0.54, wspace=0.30)
    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
    ]

    _plot_window_counts(axes[0], window)
    _plot_window_proportions(axes[1], window)
    _plot_days_since_vaccination(axes[2], window)
    _plot_policy_composition(axes[3], policy)

    dose_handles = [
        Patch(facecolor=DOSE_COLORS[label], edgecolor="none", label=label)
        for label, _ in _active_dose_groups(policy, "proportion")
    ]
    fig.legend(
        handles=dose_handles,
        loc="lower center",
        ncol=len(dose_handles),
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=1.1,
        handlelength=1.5,
        frameon=False,
    )
    fig.subplots_adjust(left=0.08, right=0.985, top=0.94, bottom=0.15)
    save_figure(fig, paths, FIGURE_NAME, tight=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    print(f"Wrote {FIGURE_NAME} to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
