"""Build Chapter 5 Figure 8: detector eligibility and selection flow."""

from __future__ import annotations
import argparse
import numpy as np
import pandas as pd


from .common import (
    Paths,
    add_common_args,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)
from ..sse.config import MIN_CLUSTER_SIZE
from ..sse.io import write_table

FILE_NAME = "ch5_selection_funnel"


def build_selection_funnel(nodes: pd.DataFrame) -> pd.DataFrame:
    tier = nodes["candidate_tier"]
    rows = [
        ("All transition-network clusters", len(nodes)),
        (
            f"Size eligible (cluster size ≥ {MIN_CLUSTER_SIZE})",
            int(nodes["sse_tested"].fillna(False).sum()),
        ),
        (
            "Background or low information",
            int(tier.eq("background_or_low_information").sum()),
        ),
        ("Possible review", int(tier.eq("possible_review").sum())),
        ("High-priority local burst", int(tier.eq("high_priority_burst").sum())),
        ("High-priority onward burden", int(tier.eq("high_priority_burden").sum())),
        ("High priority on both axes", int(tier.eq("high_priority_both_axes").sum())),
    ]
    out = pd.DataFrame(rows, columns=["stage", "n"])
    out["order"] = np.arange(len(out))
    out["pct_all"] = out["n"] / len(nodes)
    return out


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_selection_funnel(nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")
    table = table.sort_values("order")
    fig, ax = styled_new_figure(width="double", height_in=4.8, constrained_layout=True)
    colors = [
        "#4C78A8",
        "#72A0C1",
        "#B8B8B8",
        "#E6AB02",
        "#D55E00",
        "#0072B2",
        "#7B3294",
    ]
    y = np.arange(len(table))[::-1]
    bars = ax.barh(y, table["n"], color=colors)
    ax.set_xscale("log")
    ax.set_yticks(y, table["stage"])
    ax.set_xlabel("Clusters (log scale)")
    # ax.set_title("Detector eligibility and candidate-selection flow")
    for bar, value in zip(bars, table["n"]):
        if value > 0:
            ax.text(
                max(value, 0.8) * 1.12,
                bar.get_y() + bar.get_height() / 2,
                f"{int(value):,}",
                va="center",
            )
        else:
            ax.text(
                35,
                bar.get_y() + bar.get_height() / 2,
                f"{int(value):,}",
                va="center",
            )
    outputs = styled_save_figure(fig, paths, f"fig_{FILE_NAME}", tight=False)
    return {"figure": fig, "outputs": outputs}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    paths = paths_from_args(parser.parse_args())
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
