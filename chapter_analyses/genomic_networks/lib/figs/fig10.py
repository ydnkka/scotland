"""Build Chapter 4 Figure 10: test-reason composition by epidemic era."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    POLICY_LABELS,
    Paths,
    add_common_args,
    ordered_policy_values,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)

FIGURE_NAME = "fig_ch4_test_reason_by_policy_era"


def _format_count(value: float) -> str:
    if not np.isfinite(value):
        return ""
    return f"{round(value):,}"


def _pretty_reason(value: object) -> str:
    return str(value).replace("_", " ").title()


def _prepare_heatmap(summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if summary.empty:
        raise ValueError("test_reason_by_policy_era is empty")

    work = summary.copy()
    work["test_reason"] = work["test_reason"].astype("string").fillna("missing")
    work["policy_era"] = work["policy_era"].astype("string").fillna("missing")

    counts = work.pivot_table(
        index="test_reason",
        columns="policy_era",
        values="n_sequences",
        aggfunc="sum",
        fill_value=0,
    )

    policy_order = ordered_policy_values(work["policy_era"], column="policy_era")
    extra_policy_eras = [era for era in counts.columns if era not in policy_order]
    extra_policy_eras = sorted(
        extra_policy_eras,
        key=lambda value: (str(value) == "missing", str(value)),
    )
    counts = counts.reindex(columns=[*policy_order, *extra_policy_eras], fill_value=0)

    row_totals = counts.sum(axis=1).sort_values(ascending=False, kind="mergesort")
    counts = counts.loc[row_totals.index]
    column_totals = counts.sum(axis=0)
    shares = counts.div(column_totals.replace(0, np.nan), axis=1).fillna(0.0)

    return shares, counts


def build(paths: Paths) -> None:
    summary = read_table(paths, "test_reason_by_policy_era")
    shares, counts = _prepare_heatmap(summary)

    x_labels = [
        fill(
            POLICY_LABELS.get(str(era), str(era).replace("_", " ").upper()),
            width=14,
        )
        for era in shares.columns
    ]
    y_labels = [fill(_pretty_reason(reason), width=22) for reason in shares.index]
    values = shares.to_numpy(dtype=float)
    count_values = counts.to_numpy(dtype=float)
    vmax = float(np.nanmax(values))
    if not np.isfinite(vmax) or vmax <= 0:
        raise ValueError("No positive counts were available to plot")

    fig, ax = styled_new_figure(
        width="double",
        height_in=5.8,
        constrained_layout=True,
    )
    image = ax.imshow(
        values,
        aspect="auto",
        cmap="Blues",
        vmin=0.0,
        vmax=vmax,
    )

    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_xticklabels(x_labels)
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("Epidemic era")
    ax.set_ylabel("Test reason")
    plt.setp(
        ax.get_xticklabels(),
        rotation=30,
        ha="right",
        rotation_mode="anchor",
    )
    ax.set_xticks(np.arange(-0.5, len(x_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(y_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=8)

    text_threshold = vmax * 0.55
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            ax.text(
                j,
                i,
                _format_count(float(count_values[i, j])),
                ha="center",
                va="center",
                fontsize=7,
                color="white" if value >= text_threshold else "#1f1f1f",
            )

    cbar = fig.colorbar(image, ax=ax, shrink=0.92, pad=0.02)
    cbar.set_label("Share of sequences within epidemic era")
    cbar.ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))

    styled_save_figure(fig, paths, FIGURE_NAME, tight=False)


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
