"""Build cluster period summary table."""

from __future__ import annotations

import argparse

from utils import addlinespace_after_group_changes, write_latex_table

from .common import Paths, add_common_args, paths_from_args, read_table, sort_by_policy
from .table_common import (
    fmt_float,
    fmt_int,
    fmt_iqr,
)

TABLE_NAME = "tab03_cluster_period_summary"


def build(paths: Paths) -> None:
    clusters = read_table(paths, "cluster_period_summary")
    clusters = sort_by_policy(clusters, column="policy_period")
    rows = []
    for row in clusters.itertuples(index=False):
        cells = [
            str(row.policy_era).capitalize().replace("_", " "),
            row.policy_period,
            fmt_int(row.n_sequence_memberships),
            f"{fmt_int(row.n_clusters)} ({fmt_int(row.n_non_singleton_clusters)})",
        ]
        if fmt_int(row.n_non_singleton_clusters) != "0":
            cells.append(
                f"{fmt_float(row.median_non_singleton_cluster_size, 0)} "
                f"({fmt_float(row.p90_non_singleton_cluster_size, 0)}; "
                f"{fmt_int(row.max_non_singleton_cluster_size)})"
            )
            cells.append(
                f"{fmt_float(row.median_non_singleton_datazones, 0)} "
                f"({fmt_float(row.p90_non_singleton_datazones, 0)}; "
                f"{fmt_int(row.max_non_singleton_datazones)})"
            )
            cells.append(
                fmt_iqr(
                    row.median_non_singleton_spatial_distance_km,
                    row.q25_non_singleton_spatial_distance_km,
                    row.q75_non_singleton_spatial_distance_km,
                    digits=0,
                )
            )
            cells.append(
                fmt_iqr(
                    row.median_non_singleton_duration_days,
                    row.q25_non_singleton_duration_days,
                    row.q75_non_singleton_duration_days,
                    digits=0,
                )
            )
        else:
            cells.extend(["NA"] * 4)
        rows.append(cells)
    write_latex_table(
        paths.publication_table_dir / f"{TABLE_NAME}.tex",
        caption=(
            "Summary of Scottish EpiLink clusters stratified by epidemic era and policy "
            "period. The table reports sequence-window memberships and total distinct clusters, "
            "with non-singleton cluster counts provided in parentheses. Further metrics "
            "for non-singleton clusters summarise distributions [median (90th percentile; "
            "maximum)] for both cluster size and affected Data Zones. Residential spatial distance "
            "(in kilometers) and temporal span (in days) are reported as medians alongside "
            "their corresponding interquartile ranges."
        ),
        short_caption="Scottish EpiLink cluster summary by policy period",
        label="tab:cluster_period_summary",
        columns=[
            "Epidemic era",
            "Period",
            "Members",
            "Clusters",
            "Size",
            "Zones",
            "Spatial (km)",
            "Span (days)",
        ],
        rows=rows,
        addlinespace_after=addlinespace_after_group_changes(clusters["policy_era"]),
        column_spec="llrrrrrr",
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
