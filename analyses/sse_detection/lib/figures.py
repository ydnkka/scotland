"""Figure orchestration for SSE detection artifacts."""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .figs import (
    fig01,
    fig02,
    fig03,
    fig04_app,
    fig04_main,
    fig05,
    fig06,
    fig07,
    fig08,
    fig09,
    fig10,
    fig11,
    tables,
)
from .figs.common import (
    DEFAULT_RESULT_TABLE_DIR,
    DEFAULT_TABLE_DIR,
    FIGURE_DIR,
    Paths,
)
from .sse.config import BAYESIAN_OUTPUT_DIR

LOGGER = logging.getLogger(__name__)

BuildResult = Any
BuildFunction = Callable[[Paths], BuildResult]


@dataclass(frozen=True)
class ArtifactBuilder:
    name: str
    build: BuildFunction


FIGURE_BUILDERS: tuple[ArtifactBuilder, ...] = (
    ArtifactBuilder(fig01.FIGURE_NAME, fig01.build),
    ArtifactBuilder(fig02.FIGURE_NAME, fig02.build),
    ArtifactBuilder(fig03.FIGURE_NAME, fig03.build),
    ArtifactBuilder(fig04_main.FIGURE_NAME["mixing"], fig04_main.build_mixing),
    ArtifactBuilder(
        fig04_main.FIGURE_NAME["composition"],
        fig04_main.build_composition,
    ),
    ArtifactBuilder(fig04_app.FIGURE_NAME["mixing"], fig04_app.build_mixing),
    ArtifactBuilder(
        fig04_app.FIGURE_NAME["composition"],
        fig04_app.build_composition,
    ),
    ArtifactBuilder(fig05.FILE_NAME, fig05.build),
    ArtifactBuilder(fig06.FILE_NAME, fig06.build),
    ArtifactBuilder(fig07.FILE_NAME, fig07.build),
    ArtifactBuilder(fig08.FILE_NAME, fig08.build),
    ArtifactBuilder(fig09.FILE_NAME, fig09.build),
    ArtifactBuilder(fig10.FILE_NAME, fig10.build),
    ArtifactBuilder(fig11.FILE_NAME, fig11.build),
)


TABLE_BUILDERS: tuple[ArtifactBuilder, ...] = tuple(
    ArtifactBuilder(name, writer) for name, writer in tables.TABLE_WRITERS
)


def make_paths(
    *,
    table_dir: Path = DEFAULT_TABLE_DIR,
    figure_dir: Path = FIGURE_DIR,
    bayesian_result_dir: Path = BAYESIAN_OUTPUT_DIR,
    result_table_dir: Path = DEFAULT_RESULT_TABLE_DIR,
) -> Paths:
    return Paths(
        table_dir=table_dir,
        figure_dir=figure_dir,
        bayesian_result_dir=bayesian_result_dir,
        result_table_dir=result_table_dir,
    )


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
            f"Unknown SSE {kind} builder(s): {', '.join(missing)}. "
            f"Available builders: {available}"
        )
    return tuple(by_name[name] for name in requested)


def iter_figure_builders(
    names: Iterable[str] | None = None,
) -> tuple[ArtifactBuilder, ...]:
    return _iter_builders(FIGURE_BUILDERS, names, kind="figure")


def iter_table_builders(
    names: Iterable[str] | None = None,
) -> tuple[ArtifactBuilder, ...]:
    return _iter_builders(TABLE_BUILDERS, names, kind="table")


def build_figure(
    name: str,
    *,
    table_dir: Path = DEFAULT_TABLE_DIR,
    figure_dir: Path = FIGURE_DIR,
    bayesian_result_dir: Path = BAYESIAN_OUTPUT_DIR,
    result_table_dir: Path = DEFAULT_RESULT_TABLE_DIR,
) -> BuildResult:
    paths = make_paths(
        table_dir=table_dir,
        figure_dir=figure_dir,
        bayesian_result_dir=bayesian_result_dir,
        result_table_dir=result_table_dir,
    )
    builder = iter_figure_builders([name])[0]
    return builder.build(paths)


def build_table(
    name: str,
    *,
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
    output_dir: Path = DEFAULT_RESULT_TABLE_DIR,
    figure_dir: Path = FIGURE_DIR,
) -> BuildResult:
    paths = make_paths(
        figure_dir=figure_dir,
        bayesian_result_dir=result_dir,
        result_table_dir=output_dir,
    )
    builder = iter_table_builders([name])[0]
    return builder.build(paths)


def _build_all(
    builders: tuple[ArtifactBuilder, ...],
    *,
    paths: Paths,
    names: Iterable[str] | None,
    skip_missing: bool,
    logger: logging.Logger | None,
    kind: str,
) -> dict[str, BuildResult]:
    log = logger or LOGGER
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
    table_dir: Path = DEFAULT_TABLE_DIR,
    figure_dir: Path = FIGURE_DIR,
    bayesian_result_dir: Path = BAYESIAN_OUTPUT_DIR,
    result_table_dir: Path = DEFAULT_RESULT_TABLE_DIR,
    names: Iterable[str] | None = None,
    skip_missing: bool = False,
    logger: logging.Logger | None = None,
) -> dict[str, BuildResult]:
    paths = make_paths(
        table_dir=table_dir,
        figure_dir=figure_dir,
        bayesian_result_dir=bayesian_result_dir,
        result_table_dir=result_table_dir,
    )
    return _build_all(
        FIGURE_BUILDERS,
        paths=paths,
        names=names,
        skip_missing=skip_missing,
        logger=logger,
        kind="figure",
    )


def build_all_tables(
    *,
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
    output_dir: Path = DEFAULT_RESULT_TABLE_DIR,
    figure_dir: Path = FIGURE_DIR,
    names: Iterable[str] | None = None,
    skip_missing: bool = False,
    logger: logging.Logger | None = None,
) -> dict[str, BuildResult]:
    paths = make_paths(
        figure_dir=figure_dir,
        bayesian_result_dir=result_dir,
        result_table_dir=output_dir,
    )
    return _build_all(
        TABLE_BUILDERS,
        paths=paths,
        names=names,
        skip_missing=skip_missing,
        logger=logger,
        kind="table",
    )
