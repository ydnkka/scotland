"""Build Chapter 5 Figure 11: stratified null-model calibration."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from ..sse.io import write_table
from .common import (
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)

FILE_NAME = "ch5_stratified_calibration"


def build_stratified_calibration(nodes: pd.DataFrame) -> pd.DataFrame:
    tested = nodes.loc[nodes["sse_tested"].fillna(False)].copy()
    tested["size_stratum"] = pd.cut(
        tested["cluster_size"],
        [5, 7, 11, 24, np.inf],
        labels=["6-7", "8-11", "12-24", "25+"],
        include_lowest=True,
    )
    tested["coverage_stratum"] = pd.qcut(
        tested["wn_prop_sequenced"],
        4,
        labels=["Q1", "Q2", "Q3", "Q4"],
        duplicates="drop",
    )
    rows = []
    for axis in ("burst", "burden"):
        axis_frame = tested.loc[tested[f"{axis}_score_upper_p"].notna()]
        if axis == "burden":
            axis_frame = axis_frame.loc[axis_frame["burden_eligible"].fillna(False)]
        for kind, column in (
            ("Cluster size", "size_stratum"),
            ("Sequencing coverage", "coverage_stratum"),
        ):
            for stratum, group in axis_frame.groupby(column, observed=True):
                values = group[f"{axis}_score_upper_p"].to_numpy(float)
                for threshold in np.linspace(0.05, 1.0, 20):
                    rows.append(
                        {
                            "axis": axis,
                            "stratifier": kind,
                            "stratum": str(stratum),
                            "threshold": threshold,
                            "empirical_cdf": float(np.mean(values <= threshold)),
                            "n": len(values),
                        }
                    )
    return pd.DataFrame(rows)


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_stratified_calibration(nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")
    fig, axes = styled_new_figure(
        width="double",
        height_in=6.2,
        nrows=2,
        ncols=2,
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    for row, axis_name in enumerate(("burst", "burden")):
        for col, stratifier in enumerate(("Cluster size", "Sequencing coverage")):
            ax = axes[row, col]
            data = table.loc[
                table["axis"].eq(axis_name) & table["stratifier"].eq(stratifier)
            ]
            ax.plot(
                [0, 1],
                [0, 1],
                color="black",
                ls="--",
                lw=0.8,
                label="Uniform expectation",
            )
            for stratum, group in data.groupby("stratum", sort=False):
                group = group.sort_values("threshold")
                ax.plot(
                    group["threshold"],
                    group["empirical_cdf"],
                    label=f"{stratum} (n={int(group['n'].iloc[0]):,})",
                )
            ax.set_title(f"{axis_name.title()} | {stratifier.lower()}")
            ax.legend(
                loc="upper left", frameon=True, facecolor="#ffffff7b", edgecolor="#ffffff7b"
            )
            panel_label(ax, chr(ord("A") + row * 2 + col))

    fig.supxlabel("Randomized null-model p-value")
    fig.supylabel("Empirical cumulative proportion")
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
