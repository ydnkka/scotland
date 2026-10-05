"""Build cluster pairwise distance summary table."""

from __future__ import annotations

import argparse

from utils import write_latex_table

from .common import Paths, add_common_args, paths_from_args, read_table
from .table_common import (
    fmt_int,
    fmt_iqr,
)

TABLE_NAME = "tab04_cluster_pairwise_distance_summary"


def build(paths: Paths) -> None:
    summary = read_table(paths, "cluster_pairwise_distance_summary")
    metric_labels = {
        "snp_distance": "Genetic (SNPs)",
        "temporal_distance": "Temporal (days)",
    }
    metric_order = {name: idx for idx, name in enumerate(metric_labels)}
    weighting_labels = {
        "unweighted": "Unweighted",
        "pair_count_weighted": "Pair-weighted",
    }
    weighting_order = {name: idx for idx, name in enumerate(weighting_labels)}
    summary = summary.copy()
    summary["_metric_sort"] = summary["distance_metric"].map(metric_order).fillna(999)
    summary["_weighting_sort"] = summary["weighting"].map(weighting_order).fillna(999)
    summary = summary.sort_values(
        ["_metric_sort", "_weighting_sort", "distance_metric", "weighting"],
        kind="mergesort",
    )

    rows = []
    for row in summary.itertuples(index=False):
        rows.append(
            [
                metric_labels.get(str(row.distance_metric), str(row.distance_metric)),
                weighting_labels.get(str(row.weighting), str(row.weighting)),
                fmt_int(row.n_window_lineage_summaries),
                fmt_int(row.n_windows),
                fmt_int(row.n_pairwise_rows),
                fmt_iqr(row.median, row.q25, row.q75, digits=1),
            ]
        )

    write_latex_table(
        paths.publication_table_dir / f"{TABLE_NAME}.tex",
        caption=(
            "Distributions of within-cluster genetic (SNPs) and temporal (days) "
            "distances. For each metric, the table contrasts unweighted and "
            "pair-weighted quartiles of window-lineage median pairwise distances. "
            "Analyses are restricted to summaries containing a minimum of ten observed "
            "within-cluster pairs. Median pairwise distances are reported alongside "
            "corresponding interquartile ranges enclosed in parentheses."
        ),
        short_caption="Within-cluster pairwise distance summaries",
        label="tab:cluster_pairwise_distance_summary",
        columns=[
            "Metric",
            "Weighting",
            "Summaries",
            "Windows",
            "Pairs",
            "Median",
        ],
        rows=rows,
        column_spec="llrrrl",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
