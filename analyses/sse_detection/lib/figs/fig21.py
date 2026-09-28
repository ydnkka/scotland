"""Build the health board composition overview."""

from __future__ import annotations

import argparse

from .common import Paths, add_common_args, paths_from_args, styled_save_figure
from .composition_profiles import COMPOSITION_FIGURES, draw_overview

FIGURE_NAME = "fig_cluster_composition_overview_health_board"


def build(paths: Paths) -> dict[str, object]:
    fig, summaries = draw_overview(
        paths,
        COMPOSITION_FIGURES[4:],
        height=6.8,
        ratios=[1.0],
        panel_labels=("I", "J"),
    )
    outputs = styled_save_figure(fig, paths, FIGURE_NAME)
    return {"figure": fig, "outputs": outputs, "plot_data": summaries}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
