"""Central registry for publication figure and table builders."""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from analyses.genomic_networks.lib.config import TABLES_DIR as GENOMIC_TABLES_DIR
from analyses.genomic_networks.lib.figs import (
    fig01 as genomic_fig01,
    fig02 as genomic_fig02,
    fig03 as genomic_fig03,
    fig04 as genomic_fig04,
    fig05 as genomic_fig05,
    fig06 as genomic_fig06,
    fig07 as genomic_fig07,
    fig08 as genomic_fig08,
    fig09 as genomic_fig09,
    fig10 as genomic_fig10,
    tables as genomic_tables,
)
from analyses.genomic_networks.lib.figs.common import Paths as GenomicPaths
from analyses.sse_detection.lib.figs import (
    fig01 as sse_fig01,
    fig02 as sse_fig02,
    fig03 as sse_fig03,
    fig04_app as sse_fig04_app,
    fig04_main as sse_fig04_main,
    fig05 as sse_fig05,
    fig06 as sse_fig06,
    fig07 as sse_fig07,
    fig08 as sse_fig08,
    fig09 as sse_fig09,
    fig10 as sse_fig10,
    fig11 as sse_fig11,
    tables as sse_tables,
)
from analyses.sse_detection.lib.figs.common import (
    DEFAULT_RESULT_TABLE_DIR as SSE_RESULT_TABLE_DIR,
    DEFAULT_TABLE_DIR as SSE_TABLE_DIR,
    Paths as SSEPaths,
)
from analyses.sse_detection.lib.sse.config import BAYESIAN_OUTPUT_DIR
from analyses.surveillance.lib.config import TABLES_DIR as SURVEILLANCE_TABLES_DIR
from analyses.surveillance.lib.figs import (
    fig01 as surveillance_fig01,
    fig02 as surveillance_fig02,
)

from .config import FIGURES_DIR, TABLES_DIR

LOGGER = logging.getLogger(__name__)
BuilderKind = Literal["figure", "table"]
BuildFunction = Callable[["BuildContext"], Any]


@dataclass(frozen=True)
class BuildContext:
    """Output directories for top-level figures and LaTeX tables."""

    figure_dir: Path = FIGURES_DIR
    table_dir: Path = TABLES_DIR


@dataclass(frozen=True)
class ArtifactBuilder:
    """A named figure/table builder in one analysis domain."""

    domain: str
    name: str
    kind: BuilderKind
    build: BuildFunction

    @property
    def key(self) -> str:
        return f"{self.domain}:{self.name}"


DOMAINS: tuple[str, ...] = (
    "surveillance",
    "genomic_networks",
    "sse_detection",
)


def _surveillance_figure_builder(
    name: str,
    build_func: Callable[..., Any],
) -> ArtifactBuilder:
    def build(context: BuildContext) -> Any:
        return build_func(
            figure_dir=context.figure_dir,
            table_dir=SURVEILLANCE_TABLES_DIR,
            write_figure=True,
            write_tables=False,
        )

    return ArtifactBuilder("surveillance", name, "figure", build)


def _surveillance_table_builder(
    name: str,
    build_func: Callable[..., Any],
) -> ArtifactBuilder:
    def build(context: BuildContext) -> Any:
        return build_func(
            figure_dir=context.figure_dir,
            table_dir=SURVEILLANCE_TABLES_DIR,
            write_figure=False,
            write_tables=True,
        )

    return ArtifactBuilder("surveillance", name, "table", build)


def _genomic_figure_builder(
    name: str,
    build_func: Callable[[GenomicPaths], Any],
) -> ArtifactBuilder:
    def build(context: BuildContext) -> Any:
        paths = GenomicPaths(
            table_dir=GENOMIC_TABLES_DIR,
            figure_dir=context.figure_dir,
        )
        return build_func(paths)

    return ArtifactBuilder("genomic_networks", name, "figure", build)


def _genomic_table_builder(
    name: str,
    build_func: Callable[[GenomicPaths], Any],
) -> ArtifactBuilder:
    def build(context: BuildContext) -> Any:
        paths = GenomicPaths(
            table_dir=GENOMIC_TABLES_DIR,
            figure_dir=context.table_dir,
        )
        return build_func(paths)

    return ArtifactBuilder("genomic_networks", name, "table", build)


def _sse_figure_builder(
    name: str,
    build_func: Callable[[SSEPaths], Any],
) -> ArtifactBuilder:
    def build(context: BuildContext) -> Any:
        paths = SSEPaths(
            table_dir=SSE_TABLE_DIR,
            figure_dir=context.figure_dir,
            bayesian_result_dir=BAYESIAN_OUTPUT_DIR,
            result_table_dir=SSE_RESULT_TABLE_DIR,
        )
        return build_func(paths)

    return ArtifactBuilder("sse_detection", name, "figure", build)


def _sse_table_builder(
    name: str,
    build_func: Callable[[SSEPaths], Any],
) -> ArtifactBuilder:
    def build(context: BuildContext) -> Any:
        paths = SSEPaths(
            table_dir=SSE_TABLE_DIR,
            figure_dir=context.table_dir,
            bayesian_result_dir=BAYESIAN_OUTPUT_DIR,
            result_table_dir=SSE_RESULT_TABLE_DIR,
        )
        return build_func(paths)

    return ArtifactBuilder("sse_detection", name, "table", build)


def figure_builders() -> tuple[ArtifactBuilder, ...]:
    """Return all top-level figure builders in build order."""
    return (
        _surveillance_figure_builder(
            surveillance_fig01.FIGURE_NAME,
            surveillance_fig01.build,
        ),
        _surveillance_figure_builder(
            surveillance_fig02.FIGURE_NAME,
            surveillance_fig02.build,
        ),
        _genomic_figure_builder(genomic_fig02.FIGURE_NAME, genomic_fig02.build),
        _genomic_figure_builder(genomic_fig05.FIGURE_NAME, genomic_fig05.build),
        _genomic_figure_builder(genomic_fig07.FIGURE_NAME, genomic_fig07.build),
        _genomic_figure_builder(genomic_fig08.FIGURE_NAME, genomic_fig08.build),
        _genomic_figure_builder(genomic_fig03.FIGURE_NAME, genomic_fig03.build),
        _genomic_figure_builder(genomic_fig06.FIGURE_NAME, genomic_fig06.build),
        _genomic_figure_builder(genomic_fig09.FIGURE_NAME, genomic_fig09.build),
        _genomic_figure_builder(genomic_fig10.FIGURE_NAME, genomic_fig10.build),
        _genomic_figure_builder(genomic_fig01.FIGURE_NAME, genomic_fig01.build),
        _genomic_figure_builder(genomic_fig04.FIGURE_NAME, genomic_fig04.build),
        _sse_figure_builder(sse_fig01.FIGURE_NAME, sse_fig01.build),
        _sse_figure_builder(sse_fig02.FIGURE_NAME, sse_fig02.build),
        _sse_figure_builder(sse_fig03.FIGURE_NAME, sse_fig03.build),
        _sse_figure_builder(
            sse_fig04_main.FIGURE_NAME["mixing"],
            sse_fig04_main.build_mixing,
        ),
        _sse_figure_builder(
            sse_fig04_main.FIGURE_NAME["composition"],
            sse_fig04_main.build_composition,
        ),
        _sse_figure_builder(
            sse_fig04_app.FIGURE_NAME["mixing"],
            sse_fig04_app.build_mixing,
        ),
        _sse_figure_builder(
            sse_fig04_app.FIGURE_NAME["composition"],
            sse_fig04_app.build_composition,
        ),
        _sse_figure_builder(sse_fig05.FILE_NAME, sse_fig05.build),
        _sse_figure_builder(sse_fig06.FILE_NAME, sse_fig06.build),
        _sse_figure_builder(sse_fig07.FILE_NAME, sse_fig07.build),
        _sse_figure_builder(sse_fig08.FILE_NAME, sse_fig08.build),
        _sse_figure_builder(sse_fig09.FILE_NAME, sse_fig09.build),
        _sse_figure_builder(sse_fig10.FILE_NAME, sse_fig10.build),
        _sse_figure_builder(sse_fig11.FILE_NAME, sse_fig11.build),
    )


def table_builders() -> tuple[ArtifactBuilder, ...]:
    """Return all top-level table builders in build order."""
    return (
        _surveillance_table_builder(
            surveillance_fig01.FIGURE_NAME,
            surveillance_fig01.build,
        ),
        _surveillance_table_builder(
            surveillance_fig02.FIGURE_NAME,
            surveillance_fig02.build,
        ),
        *(
            _genomic_table_builder(name, build_func)
            for name, build_func in genomic_tables.TABLE_WRITERS
        ),
        *(
            _sse_table_builder(name, build_func)
            for name, build_func in sse_tables.TABLE_WRITERS
        ),
    )


def _normalise_requested(values: Iterable[str] | None) -> tuple[str, ...] | None:
    if values is None:
        return None
    requested = tuple(value for value in values if value)
    return requested or None


def _filter_by_domain(
    builders: Sequence[ArtifactBuilder],
    domains: Iterable[str] | None,
) -> tuple[ArtifactBuilder, ...]:
    selected_domains = _normalise_requested(domains)
    if selected_domains is None:
        return tuple(builders)

    unknown = sorted(set(selected_domains) - set(DOMAINS))
    if unknown:
        raise KeyError(
            "Unknown domain(s): "
            + ", ".join(unknown)
            + ". Available domains: "
            + ", ".join(DOMAINS)
        )
    return tuple(builder for builder in builders if builder.domain in selected_domains)


def select_builders(
    builders: Sequence[ArtifactBuilder],
    names: Iterable[str] | None = None,
    *,
    domains: Iterable[str] | None = None,
) -> tuple[ArtifactBuilder, ...]:
    """Select builders by optional domain and builder names.

    Names may be fully qualified (``domain:name``) or unqualified when unique
    in the filtered builder set.
    """
    filtered = _filter_by_domain(builders, domains)
    requested = _normalise_requested(names)
    if requested is None:
        return filtered

    by_key = {builder.key: builder for builder in filtered}
    selected = []
    for name in requested:
        if name in by_key:
            selected.append(by_key[name])
            continue

        matches = [builder for builder in filtered if builder.name == name]
        if not matches:
            available = ", ".join(builder.key for builder in filtered)
            raise KeyError(
                f"Unknown builder {name!r}. Available builders: {available}"
            )
        if len(matches) > 1:
            options = ", ".join(builder.key for builder in matches)
            raise KeyError(
                f"Ambiguous builder {name!r}; qualify it as one of: {options}"
            )
        selected.append(matches[0])
    return tuple(selected)


def list_builders(
    builders: Sequence[ArtifactBuilder],
    *,
    domains: Iterable[str] | None = None,
) -> str:
    """Return a newline-delimited builder listing."""
    filtered = select_builders(builders, domains=domains)
    return "\n".join(builder.key for builder in filtered)


def build_selected(
    builders: Sequence[ArtifactBuilder],
    *,
    names: Iterable[str] | None = None,
    domains: Iterable[str] | None = None,
    figure_dir: Path = FIGURES_DIR,
    table_dir: Path = TABLES_DIR,
    skip_missing: bool = False,
    logger: logging.Logger | None = None,
) -> dict[str, Any]:
    """Build selected artifacts and return results keyed by ``domain:name``."""
    log = logger or LOGGER
    context = BuildContext(figure_dir=figure_dir, table_dir=table_dir)
    context.figure_dir.mkdir(parents=True, exist_ok=True)
    context.table_dir.mkdir(parents=True, exist_ok=True)

    outputs: dict[str, Any] = {}
    for builder in select_builders(builders, names, domains=domains):
        try:
            log.info("Building %s %s", builder.kind, builder.key)
            outputs[builder.key] = builder.build(context)
        except FileNotFoundError as exc:
            if not skip_missing:
                raise
            log.warning("Skipping %s %s: %s", builder.kind, builder.key, exc)
    return outputs


def build_figures(
    *,
    names: Iterable[str] | None = None,
    domains: Iterable[str] | None = None,
    figure_dir: Path = FIGURES_DIR,
    table_dir: Path = TABLES_DIR,
    skip_missing: bool = False,
    logger: logging.Logger | None = None,
) -> dict[str, Any]:
    return build_selected(
        figure_builders(),
        names=names,
        domains=domains,
        figure_dir=figure_dir,
        table_dir=table_dir,
        skip_missing=skip_missing,
        logger=logger,
    )


def build_tables(
    *,
    names: Iterable[str] | None = None,
    domains: Iterable[str] | None = None,
    figure_dir: Path = FIGURES_DIR,
    table_dir: Path = TABLES_DIR,
    skip_missing: bool = False,
    logger: logging.Logger | None = None,
) -> dict[str, Any]:
    return build_selected(
        table_builders(),
        names=names,
        domains=domains,
        figure_dir=figure_dir,
        table_dir=table_dir,
        skip_missing=skip_missing,
        logger=logger,
    )
