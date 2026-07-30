"""Build Chapter 5 Figure 9: candidate rates by epidemic context."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from ..sse.io import write_table
from .common import (
    CLADES,
    HIGH_PRIORITY,
    POLICY_LABELS,
    POLICY_ORDER,
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    read_table,
    sort_by_policy,
    styled_new_figure,
    styled_save_figure,
    wilson,
)

FILE_NAME = "ch5_candidate_context_rates"

EPOCHS = [
    ("Epidemic era", "policy_era"),
    ("Policy period", "policy_period"),
    ("Clade", "clade"),
]

FACTOR = 100


def build_group_rates(nodes: pd.DataFrame) -> pd.DataFrame:
    tested = nodes.loc[nodes["sse_tested"].fillna(False)].copy()
    tested["clade"] = tested["clade"].map(CLADES).fillna("Other")
    tested["candidate"] = tested["candidate_tier"].isin(HIGH_PRIORITY)
    rows = []
    for dimension, column in EPOCHS:
        for value, group in tested.groupby(column, dropna=False, observed=True):
            n, k = len(group), int(group["candidate"].sum())
            rows.append(
                {
                    "dimension": dimension,
                    column: str(value),
                    "eligible_n": n,
                    "candidate_n": k,
                }
            )
    out = pd.DataFrame(rows)
    out["candidate_rate"] = out["candidate_n"] / out["eligible_n"]
    out["ci95_low"], out["ci95_high"] = wilson(out["candidate_n"], out["eligible_n"])
    return out


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_group_rates(nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")
    fig, axes = styled_new_figure(
        width="double",
        height_in=6.2,
        nrows=1,
        ncols=3,
        constrained_layout=True,
    )
    for ax, (dimension, column), label in zip(axes, EPOCHS, "ABC"):
        data = table.loc[table["dimension"].eq(dimension)].copy()
        if column in POLICY_ORDER:
            data = sort_by_policy(data, column)
        else:
            data = data.sort_values(column, ignore_index=True)
        if column == "policy_era":
            data[column] = data[column].map(POLICY_LABELS)
        rate = (FACTOR * data["candidate_rate"]).to_numpy(dtype=float)
        y = np.arange(len(data))
        # Clip tiny negative values introduced by floating-point roundoff.
        xerr = np.vstack(
            (
                np.maximum(rate - FACTOR * data["ci95_low"].to_numpy(dtype=float), 0.0),
                np.maximum(
                    FACTOR * data["ci95_high"].to_numpy(dtype=float) - rate, 0.0
                ),
            )
        )
        ax.errorbar(
            rate,
            y,
            xerr=xerr,
            fmt="o",
            color="#2F6690",
            capsize=2,
        )

        ax.set_yticks(y, data[column], rotation=0)
        ax.set_title(f"{dimension}")
        panel_label(ax, label)
    fig.supxlabel(f"Candidates per {FACTOR} eligible clusters")
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
