"""Build Chapter 5 Figure 5: high-priority candidates over time."""

from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
from .common import (
    Paths,
    add_common_args,
    add_policy_bands,
    date_axis,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
    wilson,
)
from ..sse.io import write_table


FILE_NAME = "ch5_candidate_timeline"


def build_candidate_timeline(nodes: pd.DataFrame) -> pd.DataFrame:
    tested = nodes.loc[nodes["sse_tested"].fillna(False)].copy()
    tested["burst_candidate"] = tested["candidate_tier"].isin(
        {"high_priority_burst", "high_priority_both_axes"}
    )
    tested["burden_candidate"] = tested["candidate_tier"].isin(
        {"high_priority_burden", "high_priority_both_axes"}
    )
    rows = []
    keys = ["window_id", "window_idx", "wn_mid_date", "policy_period", "who_voc"]
    for key, group in tested.groupby(keys, dropna=False, observed=True):
        for axis in ("burst", "burden"):
            eligible = (
                group
                if axis == "burst"
                else group.loc[group["burden_eligible"].fillna(False)]
            )
            n = len(eligible)
            k = int(eligible[f"{axis}_candidate"].sum()) if n else 0
            rows.append(
                dict(zip(keys, key)) | {"axis": axis, "eligible_n": n, "candidate_n": k}
            )
    out = pd.DataFrame(rows)
    out["candidate_rate"] = out["candidate_n"].div(out["eligible_n"].replace(0, np.nan))
    low, high = wilson(out["candidate_n"], out["eligible_n"].replace(0, np.nan))
    out["ci95_low"], out["ci95_high"] = low, high
    return out


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_candidate_timeline(nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")
    table["wn_mid_date"] = pd.to_datetime(table["wn_mid_date"])
    summary = (
        table.groupby(["wn_mid_date", "policy_period", "axis"], observed=True)[
            ["eligible_n", "candidate_n"]
        ]
        .sum()
        .reset_index()
    )
    summary["candidate_rate"] = summary["candidate_n"] / summary["eligible_n"].replace(
        0, np.nan
    )
    fig, axes = styled_new_figure(
        width="double",
        height_in=5.6,
        nrows=2,
        ncols=1,
        sharex=True,
        constrained_layout=True,
    )
    colors = {"burst": "#D55E00", "burden": "#0072B2"}
    labels = {"burst": "Local burst", "burden": "Onward burden"}
    context = summary.drop_duplicates(["wn_mid_date", "policy_period"])
    for ax in axes:
        add_policy_bands(ax, context)
    for axis, group in summary.groupby("axis"):
        group = group.sort_values("wn_mid_date")
        axes[0].plot(
            group["wn_mid_date"],
            group["candidate_n"],
            color=colors[str(axis)],
            label=labels[str(axis)],
        )
        axes[1].plot(
            group["wn_mid_date"],
            100 * group["candidate_rate"],
            color=colors[str(axis)],
            label=labels[str(axis)],
        )
    axes[0].set_title("High-priority candidates over epidemic time")
    axes[0].set_ylabel("Candidates")
    axes[0].legend(loc="upper left")
    axes[1].set_ylabel("Candidates per 100\neligible clusters")
    axes[1].set_xlabel("Window midpoint")
    date_axis(axes[1])
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
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
