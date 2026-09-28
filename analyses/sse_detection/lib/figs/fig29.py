"""Build sensitivity entropy figure."""

from __future__ import annotations

import argparse

from .common import Paths, add_common_args, paths_from_args
from .size_profiles import (
    SCENARIO_LABELS,
    SCENARIO_ORDER,
    _finish,
    draw_entropy,
    load_saved_summaries,
    prepare,
)

FIGURE_NAME = "fig_sensitivity_entropy"


def build(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    data = tables["entropy"]
    data = data.loc[data.scenario.isin(SCENARIO_ORDER) & data.size_band.eq("All")]
    fig = draw_entropy(
        data,
        xcol="scenario",
        order=SCENARIO_ORDER,
        labels=SCENARIO_LABELS,
        sensitivity=True,
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
