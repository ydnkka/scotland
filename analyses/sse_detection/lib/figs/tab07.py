"""Build size overview table."""

from __future__ import annotations

import argparse

from ..characterisation import SIZE_BANDS
from .common import Paths, add_common_args, paths_from_args
from .size_profiles import load_saved_summaries, prepare
from .table_common import write_size_table

TABLE_NAME = "tab_size_overview"


def build(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    size = tables["size_summary"]
    output_path = paths.publication_table_dir / f"{TABLE_NAME}.tex"
    rows = []
    for band in ("All", *SIZE_BANDS):
        for group in ("CSC", "Background"):
            r = size.loc[
                size.scenario.eq("primary")
                & size.size_band.eq(band)
                & size.group.eq(group)
            ].iloc[0]
            rows.append(
                [
                    band,
                    group,
                    f"{r.n_clusters:,}",
                    f"{r.n_records:,}",
                    f"{r.size_median:g} ({r.size_q25:g}--{r.size_q75:g})",
                ]
            )
    write_size_table(
        output_path,
        "Cluster sizes and sequence records in the primary comparison. Records count sequence appearances in retained windows, not unique people. The interquartile range (IQR) describes the middle half of cluster sizes.",
        "tab:size_overview",
        ["Size band", "Group", "Clusters", "Records", "Median (IQR)"],
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
