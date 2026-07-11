"""Figure orchestration for SSE detection artifacts."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Any, Callable, Iterable

from .figs import fig05
from .figs.common import DEFAULT_TABLE_DIR, FIGURE_DIR, Paths


LOGGER = logging.getLogger(__name__)

BuildResult = Any
BuildFunction = Callable[[Paths], BuildResult]


@dataclass(frozen=True)
class ArtifactBuilder:
    name: str
    build: BuildFunction


FIGURE_BUILDERS: tuple[ArtifactBuilder, ...] = (
    ArtifactBuilder(fig05.FIGURE_NAME, fig05.build),
)


def make_paths(
    *,
    table_dir: Path = DEFAULT_TABLE_DIR,
    figure_dir: Path = FIGURE_DIR,
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
            f"Unknown SSE {kind} builder(s): {', '.join(missing)}. "
            f"Available builders: {available}"
        )
    return tuple(by_name[name] for name in requested)


def iter_figure_builders(
    names: Iterable[str] | None = None,
) -> tuple[ArtifactBuilder, ...]:
    return _iter_builders(FIGURE_BUILDERS, names, kind="figure")


def build_figure(
    name: str,
    *,
    table_dir: Path = DEFAULT_TABLE_DIR,
    figure_dir: Path = FIGURE_DIR,
) -> BuildResult:
    paths = make_paths(table_dir=table_dir, figure_dir=figure_dir)
    builder = iter_figure_builders([name])[0]
    return builder.build(paths)


def build_all_figures(
    *,
    table_dir: Path = DEFAULT_TABLE_DIR,
    figure_dir: Path = FIGURE_DIR,
    names: Iterable[str] | None = None,
    skip_missing: bool = False,
    logger: logging.Logger | None = None,
) -> dict[str, BuildResult]:
    log = logger or LOGGER
    paths = make_paths(table_dir=table_dir, figure_dir=figure_dir)
    results: dict[str, BuildResult] = {}

    for builder in iter_figure_builders(names):
        try:
            log.info("Writing figure %s", builder.name)
            results[builder.name] = builder.build(paths)
        except FileNotFoundError as exc:
            if not skip_missing:
                raise
            log.warning("Skipping figure %s: %s", builder.name, exc)

    return results
