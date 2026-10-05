"""Build route entropy figure."""

from __future__ import annotations

import argparse

from ..characterisation import ENTROPY_ATTRIBUTES
from .common import (
    Paths,
    add_common_args,
    add_panel_labels,
    new_figure,
    paths_from_args,
)
from .size_profiles import _finish, load_saved_summaries, prepare

FIGURE_NAME = "fig27_route_entropy"


def build(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    data = tables["entropy"].loc[lambda t: t.scenario.eq("route")]
    groups = [
        "Background",
        "Burden-eligible background",
        "Burst only",
        "Burden only",
        "Both axes",
    ]
    counts = (
        tables["size_summary"]
        .loc[lambda t: t.scenario.eq("route")]
        .set_index("group")
        .n_clusters
    )
    short = [
        f"{g.replace('Burden-eligible background', 'Burden background')} ({counts[g]:,})"
        for g in groups
    ]
    colors = ["#888888", "#444444", "#D55E00", "#0072B2", "#009E73"]
    fig, axes = new_figure(width="double", nrows=6, ncols=2, height_in=9.4)
    for row, (key, _, label) in enumerate(ENTROPY_ATTRIBUTES):
        for col, scale in enumerate(("obs", "z")):
            ax = axes[row, col]
            block = (
                data.loc[data.attribute.eq(key) & data.scale.eq(scale)]
                .set_index("group")
                .reindex(groups)
            )
            for y, (group, color) in enumerate(zip(groups, colors)):
                b = block.loc[group]
                ax.errorbar(
                    b["median"],
                    y,
                    xerr=[[b["median"] - b.q25], [b.q75 - b["median"]]],
                    fmt="o",
                    markersize=3,
                    capsize=2,
                    color=color,
                )
            ax.set_yticks(range(5), short if col == 0 else [""] * 5, fontsize=8)
            ax.invert_yaxis()
            ax.set_title(
                label + (" · observed" if scale == "obs" else " · null-standardised"),
                fontsize=8,
                loc="left",
            )
            ax.tick_params(axis="x", labelsize=9)
            if scale == "obs":
                ax.set_xlim(-0.03, 1.03)
            else:
                ax.axvline(0, color="#777777", ls=":", lw=0.7)
    add_panel_labels(axes)
    fig.tight_layout(h_pad=0.8)
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
