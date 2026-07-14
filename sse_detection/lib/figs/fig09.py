"""Build Chapter 5 Figure 9: candidate rates by epidemic context."""

from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
from .common import (
    HIGH_PRIORITY,
    POLICY_ORDER,
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
    wilson,
)
from ..sse.io import write_table

FILE_NAME = "ch5_candidate_context_rates"


def build_group_rates(nodes: pd.DataFrame) -> pd.DataFrame:
    tested = nodes.loc[nodes["sse_tested"].fillna(False)].copy()
    tested["candidate"] = tested["candidate_tier"].isin(HIGH_PRIORITY)
    rows = []
    for dimension, column in (
        ("Policy period", "policy_period"),
        ("Variant", "who_voc"),
    ):
        for value, group in tested.groupby(column, dropna=False, observed=True):
            n, k = len(group), int(group["candidate"].sum())
            rows.append(
                {
                    "dimension": dimension,
                    "group": str(value),
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
        width="double", height_in=6.2, nrows=2, ncols=1, constrained_layout=True
    )
    for ax, dimension, label in zip(axes, ("Policy period", "Variant"), "AB"):
        data = table.loc[table["dimension"].eq(dimension)].copy()
        if dimension == "Policy period":
            rank = {value: idx for idx, value in enumerate(POLICY_ORDER)}
            data = data.sort_values("group", key=lambda s: s.map(rank).fillna(999))
        else:
            data = data.sort_values("eligible_n", ascending=False)
        x = np.arange(len(data))
        rate = 100 * data["candidate_rate"]
        ax.errorbar(
            x,
            rate,
            yerr=[rate - 100 * data["ci95_low"], 100 * data["ci95_high"] - rate],
            fmt="o",
            color="#2F6690",
            capsize=2,
        )
        rotation = 35 if dimension == "Variant" else 0
        ax.set_xticks(
            x, data["group"], rotation=rotation, ha="right" if rotation else "center"
        )
        ax.set_ylabel("Candidates per 100 eligible clusters")
        ax.set_title(f"Candidate rate by {dimension.lower()}")
        panel_label(ax, label)
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
