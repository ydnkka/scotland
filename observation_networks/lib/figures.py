"""Figure orchestration for Chapter 4 thesis artifacts.

The plotting implementations live in :mod:`observation_networks.lib.figs`.
This module keeps a stable programmatic entry point for rebuilding the full
figure set from saved tables.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Callable, Iterable, Any

from .config import FIGURES_DIR, TABLES_DIR
from .figs import (
    fig01,
    fig02,
    fig03,
    fig04,
    fig05,
    fig06,
    sfig01,
    sfig02,
    sfig03,
    sfig05,
    tables,
)
from .figs.common import Paths


LOGGER = logging.getLogger(__name__)

BuildResult = Any
BuildFunction = Callable[[Paths], BuildResult]


@dataclass(frozen=True)
class ArtifactBuilder:
    name: str
    build: BuildFunction


FIGURE_BUILDERS: tuple[ArtifactBuilder, ...] = (
    ArtifactBuilder("fig01_sequence_composition_by_policy", fig01.build),
    ArtifactBuilder("fig02_cluster_landscape", fig02.build),
    ArtifactBuilder("fig03_assortativity_baseline", fig03.build),
    ArtifactBuilder("fig04_mixing_matrices", fig04.build),
    ArtifactBuilder("fig05_transition_graph_baseline", fig05.build),
    ArtifactBuilder(fig06.FIGURE_NAME, fig06.build),
    ArtifactBuilder(sfig01.FIGURE_NAME, sfig01.build),
    ArtifactBuilder("sfig02_compatibility_topology", sfig02.build),
    ArtifactBuilder("sfig03_simd_population_weighting", sfig03.build),
    ArtifactBuilder("sfig05_assortativity_confidence_intervals", sfig05.build),
)

TABLE_BUILDERS: tuple[ArtifactBuilder, ...] = tuple(
    ArtifactBuilder(name, build) for name, build in tables.TABLE_WRITERS
)


def make_paths(
    *,
    table_dir: Path = TABLES_DIR,
    figure_dir: Path = FIGURES_DIR,
) -> Paths:
    return Paths(table_dir=table_dir, figure_dir=figure_dir)


def _iter_builders(
    builders: tuple[ArtifactBuilder, ...],
    names: Iterable[str] | None = None,
    *,
    kind: str,
) -> tuple[ArtifactBuilder, ...]:
    if names is None:
        return builders

    requested = (names,) if isinstance(names, str) else tuple(names)
    by_name = {builder.name: builder for builder in builders}
    missing = [name for name in requested if name not in by_name]
    if missing:
        available = ", ".join(by_name)
        raise KeyError(
            f"Unknown Chapter 4 {kind} builder(s): {', '.join(missing)}. "
            f"Available builders: {available}"
        )
    return tuple(by_name[name] for name in requested)


def iter_figure_builders(names: Iterable[str] | None = None) -> tuple[ArtifactBuilder, ...]:
    return _iter_builders(FIGURE_BUILDERS, names, kind="figure")


def iter_table_builders(names: Iterable[str] | None = None) -> tuple[ArtifactBuilder, ...]:
    return _iter_builders(TABLE_BUILDERS, names, kind="table")


def build_figure(
    name: str,
    *,
    table_dir: Path = TABLES_DIR,
    figure_dir: Path = FIGURES_DIR,
) -> BuildResult:
    paths = make_paths(table_dir=table_dir, figure_dir=figure_dir)
    builder = iter_figure_builders([name])[0]
    return builder.build(paths)


def build_table(
    name: str,
    *,
    table_dir: Path = TABLES_DIR,
    figure_dir: Path = FIGURES_DIR,
) -> BuildResult:
    paths = make_paths(table_dir=table_dir, figure_dir=figure_dir)
    builder = iter_table_builders([name])[0]
    return builder.build(paths)


def _build_all(
    builders: tuple[ArtifactBuilder, ...],
    *,
    table_dir: Path,
    figure_dir: Path,
    names: Iterable[str] | None,
    skip_missing: bool,
    logger: logging.Logger | None,
    kind: str,
) -> dict[str, BuildResult]:
    log = logger or LOGGER
    paths = make_paths(table_dir=table_dir, figure_dir=figure_dir)
    results: dict[str, BuildResult] = {}

    for builder in _iter_builders(builders, names, kind=kind):
        try:
            log.info("Writing %s %s", kind, builder.name)
            results[builder.name] = builder.build(paths)
        except FileNotFoundError as exc:
            if not skip_missing:
                raise
            log.warning("Skipping %s %s: %s", kind, builder.name, exc)

    return results


def build_all_figures(
    *,
    table_dir: Path = TABLES_DIR,
    figure_dir: Path = FIGURES_DIR,
    names: Iterable[str] | None = None,
    skip_missing: bool = False,
    logger: logging.Logger | None = None,
) -> dict[str, BuildResult]:
    return _build_all(
        FIGURE_BUILDERS,
        table_dir=table_dir,
        figure_dir=figure_dir,
        names=names,
        skip_missing=skip_missing,
        logger=logger,
        kind="figure",
    )


def build_all_tables(
    *,
    table_dir: Path = TABLES_DIR,
    figure_dir: Path = FIGURES_DIR,
    names: Iterable[str] | None = None,
    skip_missing: bool = False,
    logger: logging.Logger | None = None,
) -> dict[str, BuildResult]:
    return _build_all(
        TABLE_BUILDERS,
        table_dir=table_dir,
        figure_dir=figure_dir,
        names=names,
        skip_missing=skip_missing,
        logger=logger,
        kind="table",
    )
