"""Build Chapter 4 Supplementary Figure 5: assortativity confidence intervals."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from assortativity_panels import plot_compatibility_assortativity_grid  # noqa: E402
from common import (  # noqa: E402
    Paths,
    add_common_args,
    paths_from_args,
)


FIGURE_NAME = "fig_ch4_assortativity_confidence_intervals"


def build(paths: Paths) -> None:
    plot_compatibility_assortativity_grid(
        paths,
        figure_name=FIGURE_NAME,
        uncertainty="ci",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    print(f"Wrote {FIGURE_NAME} to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
