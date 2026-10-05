"""Build composition by size cluster weighted figure."""

from __future__ import annotations

import argparse

from ..characterisation import SIZE_BANDS
from .common import Paths, add_common_args, paths_from_args
from .size_profiles import _finish, draw_composition, load_saved_summaries, prepare

FIGURE_NAME = "fig25_composition_by_size_cluster_weighted"


def build(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    data = tables["composition_contrasts"]
    data = data.loc[
        data.weighting.eq("clusters")
        & (data.scenario.eq("primary") & data.size_band.ne("All"))
    ].copy()
    fig = draw_composition(
        data,
        xcol="size_band",
        order=SIZE_BANDS,
        labels=SIZE_BANDS,
    )
    return _finish(fig, paths, FIGURE_NAME)


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
