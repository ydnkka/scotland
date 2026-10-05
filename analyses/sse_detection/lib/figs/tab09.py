"""Build route sizes table."""

from __future__ import annotations

import argparse

from .common import Paths, add_common_args, paths_from_args
from .size_profiles import load_saved_summaries, prepare
from .table_common import write_size_table

TABLE_NAME = "tab09_route_sizes"


def build(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    size = tables["size_summary"]
    output_path = paths.publication_table_dir / f"{TABLE_NAME}.tex"
    rows = []
    for r in size.loc[size.scenario.eq("route")].itertuples():
        rows.append(
            [
                r.group.replace("Burden-eligible background", "Burden background"),
                f"{r.n_clusters:,}",
                f"{r.n_records:,}",
                f"{r.size_median:g} ({r.size_q25:g}--{r.size_q75:g})",
            ]
        )
    write_size_table(
        output_path,
        "Cluster sizes by screening route. Burden background is the subset of background clusters with positive observed downstream burden; it is the comparator for burden-only and both-axis candidates. Rows therefore do not form a disjoint partition.",
        "tab:route_sizes",
        ["Group", "Clusters", "Records", "Median (IQR)"],
        rows,
    )
    return {"tex": output_path}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument(
        "--from-saved-tables",
        action="store_true",
        help="Render saved size aggregates without rebuilding or writing them.",
    )
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths, tables=load_saved_summaries(paths) if args.from_saved_tables else None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
