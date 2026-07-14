"""Build Chapter 5 Figure 10: detector-threshold robustness."""

from __future__ import annotations
import argparse
import pandas as pd
from .common import (
    HIGH_PRIORITY,
    Paths,
    add_common_args,
    panel_label,
    paths_from_args,
    read_table,
    styled_new_figure,
    styled_save_figure,
)
from ..sse.io import write_table

FILE_NAME = "ch5_threshold_robustness"


def build_threshold_robustness(nodes: pd.DataFrame) -> pd.DataFrame:
    tested = nodes.loc[nodes["sse_tested"].fillna(False)].copy()
    baseline = set(
        tested.loc[tested["candidate_tier"].isin(HIGH_PRIORITY), "cluster_id"]
    )
    rows = []
    for min_size in (6, 8, 10, 15, 20):
        size_ok = tested["cluster_size"].ge(min_size)
        for alpha in (0.01, 0.025, 0.05, 0.10):
            selected = tested.loc[
                size_ok
                & (
                    tested["burst_score_upper_p"].le(alpha)
                    | (
                        tested["burden_eligible"].fillna(False)
                        & tested["burden_score_upper_p"].le(alpha)
                    )
                ),
                "cluster_id",
            ]
            chosen = set(selected)
            union = chosen | baseline
            rows.append(
                {
                    "min_cluster_size": min_size,
                    "alpha": alpha,
                    "eligible_n": int(size_ok.sum()),
                    "candidate_n": len(chosen),
                    "baseline_overlap_n": len(chosen & baseline),
                    "baseline_jaccard": len(chosen & baseline) / len(union)
                    if union
                    else 1.0,
                }
            )
    return pd.DataFrame(rows)


def build(paths: Paths) -> dict[str, object]:
    nodes = read_table(paths, "cluster_table")
    table = build_threshold_robustness(nodes)
    write_table(table, paths.result_table_dir, f"tab_{FILE_NAME}")
    sizes = sorted(table["min_cluster_size"].unique())
    alphas = sorted(table["alpha"].unique())
    fig, axes = styled_new_figure(
        width="double", height_in=4.2, nrows=1, ncols=2, constrained_layout=True
    )
    specs = (
        (axes[0], "candidate_n", "Candidates retained", "Blues", "A"),
        (axes[1], "baseline_jaccard", "Agreement with primary set", "Greens", "B"),
    )
    for ax, value, title, cmap, label in specs:
        matrix = table.pivot(
            index="min_cluster_size", columns="alpha", values=value
        ).reindex(index=sizes, columns=alphas)
        image = ax.imshow(
            matrix,
            aspect="auto",
            cmap=cmap,
            vmin=0,
            vmax=1 if value == "baseline_jaccard" else None,
        )
        ax.set_xticks(range(len(alphas)), [f"{a:.3g}" for a in alphas])
        ax.set_yticks(range(len(sizes)), sizes)
        ax.set_xlabel("Upper-tail significance threshold")
        ax.set_ylabel("Minimum cluster size")
        ax.set_title(title)
        for i in range(len(sizes)):
            for j in range(len(alphas)):
                cell = matrix.iloc[i, j]
                ax.text(
                    j,
                    i,
                    f"{cell:.2f}" if value == "baseline_jaccard" else f"{cell:.0f}",
                    ha="center",
                    va="center",
                )
        fig.colorbar(
            image,
            ax=ax,
            label="Jaccard similarity" if value == "baseline_jaccard" else "Candidates",
        )
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
