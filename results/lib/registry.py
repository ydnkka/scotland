"""Central registry for publication figure and table builders."""

from __future__ import annotations

import json
import logging
import re
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from analyses.genomic_networks.lib.config import TABLES_DIR as GENOMIC_TABLES_DIR
from analyses.genomic_networks.lib.figs import fig01 as genomic_fig01
from analyses.genomic_networks.lib.figs import fig02 as genomic_fig02
from analyses.genomic_networks.lib.figs import fig03 as genomic_fig03
from analyses.genomic_networks.lib.figs import fig04 as genomic_fig04
from analyses.genomic_networks.lib.figs import fig05 as genomic_fig05
from analyses.genomic_networks.lib.figs import fig06 as genomic_fig06
from analyses.genomic_networks.lib.figs import fig07 as genomic_fig07
from analyses.genomic_networks.lib.figs import fig08 as genomic_fig08
from analyses.genomic_networks.lib.figs import fig09 as genomic_fig09
from analyses.genomic_networks.lib.figs import fig10 as genomic_fig10
from analyses.genomic_networks.lib.figs import tab01 as genomic_tab01
from analyses.genomic_networks.lib.figs import tab02 as genomic_tab02
from analyses.genomic_networks.lib.figs import tab03 as genomic_tab03
from analyses.genomic_networks.lib.figs import tab04 as genomic_tab04
from analyses.genomic_networks.lib.figs import tab05 as genomic_tab05
from analyses.genomic_networks.lib.figs import tab06 as genomic_tab06
from analyses.genomic_networks.lib.figs import tab07 as genomic_tab07
from analyses.genomic_networks.lib.figs.common import Paths as GenomicPaths
from analyses.sse_detection.lib.figs import fig01 as sse_fig01
from analyses.sse_detection.lib.figs import fig02 as sse_fig02
from analyses.sse_detection.lib.figs import fig03 as sse_fig03
from analyses.sse_detection.lib.figs import fig04 as sse_fig04
from analyses.sse_detection.lib.figs import fig05 as sse_fig05
from analyses.sse_detection.lib.figs import fig06 as sse_fig06
from analyses.sse_detection.lib.figs import fig07 as sse_fig07
from analyses.sse_detection.lib.figs import fig08 as sse_fig08
from analyses.sse_detection.lib.figs import fig09 as sse_fig09
from analyses.sse_detection.lib.figs import fig10 as sse_fig10
from analyses.sse_detection.lib.figs import fig11 as sse_fig11
from analyses.sse_detection.lib.figs import fig13 as sse_fig13
from analyses.sse_detection.lib.figs import fig14 as sse_fig14
from analyses.sse_detection.lib.figs import fig15 as sse_fig15
from analyses.sse_detection.lib.figs import fig20 as sse_fig20
from analyses.sse_detection.lib.figs import fig21 as sse_fig21
from analyses.sse_detection.lib.figs import fig22 as sse_fig22
from analyses.sse_detection.lib.figs import fig23 as sse_fig23
from analyses.sse_detection.lib.figs import fig24 as sse_fig24
from analyses.sse_detection.lib.figs import fig25 as sse_fig25
from analyses.sse_detection.lib.figs import fig26 as sse_fig26
from analyses.sse_detection.lib.figs import fig27 as sse_fig27
from analyses.sse_detection.lib.figs import fig28 as sse_fig28
from analyses.sse_detection.lib.figs import fig29 as sse_fig29
from analyses.sse_detection.lib.figs import tab01 as sse_tab01
from analyses.sse_detection.lib.figs import tab02 as sse_tab02
from analyses.sse_detection.lib.figs import tab03 as sse_tab03
from analyses.sse_detection.lib.figs import tab04 as sse_tab04
from analyses.sse_detection.lib.figs import tab05 as sse_tab05
from analyses.sse_detection.lib.figs import tab06 as sse_tab06
from analyses.sse_detection.lib.figs import tab07 as sse_tab07
from analyses.sse_detection.lib.figs import tab08 as sse_tab08
from analyses.sse_detection.lib.figs import tab09 as sse_tab09
from analyses.sse_detection.lib.figs import tab10 as sse_tab10
from analyses.sse_detection.lib.figs.common import (
    DEFAULT_RESULT_TABLE_DIR as SSE_RESULT_TABLE_DIR,
)
from analyses.sse_detection.lib.figs.common import DEFAULT_TABLE_DIR as SSE_TABLE_DIR
from analyses.sse_detection.lib.figs.common import Paths as SSEPaths
from analyses.sse_detection.lib.sse.config import BAYESIAN_OUTPUT_DIR
from analyses.surveillance.lib.config import TABLES_DIR as SURVEILLANCE_TABLES_DIR
from analyses.surveillance.lib.figs import fig01 as surveillance_fig01
from analyses.surveillance.lib.figs import fig02 as surveillance_fig02
from analyses.surveillance.lib.figs import tab01 as surveillance_tab01
from analyses.surveillance.lib.figs import tab02 as surveillance_tab02
from analyses.surveillance.lib.figs import tab03 as surveillance_tab03
from analyses.surveillance.lib.figs import tab04 as surveillance_tab04
from analyses.surveillance.lib.figs import tab05 as surveillance_tab05
from analyses.surveillance.lib.figs import tab06 as surveillance_tab06
from analyses.surveillance.lib.figs import tab07 as surveillance_tab07
from analyses.surveillance.lib.figs import tab08 as surveillance_tab08

from .config import FIGURES_DIR, RESULTS_DIR, TABLES_DIR

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
    module: str = ""

    @property
    def key(self) -> str:
        return f"{self.domain}:{self.name}"

    @property
    def script(self) -> str:
        return self.module.rsplit(".", 1)[-1]

    @property
    def legacy_name(self) -> str:
        prefix = "fig" if self.kind == "figure" else "tab"
        return f"{prefix}_{self.name.removeprefix(self.script + '_')}"


DOMAINS: tuple[str, ...] = (
    "surveillance",
    "genomic_networks",
    "sse_detection",
)


def _validate_numbered_name(name: str, build_func: Callable[..., Any]) -> str:
    module = build_func.__module__
    script = module.rsplit(".", 1)[-1]
    if not re.fullmatch(r"(?:fig|tab)\d{2}", script) or not name.startswith(script + "_"):
        raise ValueError(f"Export stem {name!r} must begin with its producer script {script!r}")
    return module


def _surveillance_figure_builder(
    name: str,
    build_func: Callable[..., Any],
) -> ArtifactBuilder:
    module = _validate_numbered_name(name, build_func)
    def build(context: BuildContext) -> Any:
        return build_func(
            figure_dir=context.figure_dir,
        )

    return ArtifactBuilder("surveillance", name, "figure", build, module)


def _surveillance_table_builder(
    name: str,
    build_func: Callable[..., Any],
) -> ArtifactBuilder:
    module = _validate_numbered_name(name, build_func)
    def build(context: BuildContext) -> Any:
        return build_func(
            table_dir=SURVEILLANCE_TABLES_DIR,
        )

    return ArtifactBuilder("surveillance", name, "table", build, module)


def _genomic_figure_builder(
    name: str,
    build_func: Callable[[GenomicPaths], Any],
) -> ArtifactBuilder:
    module = _validate_numbered_name(name, build_func)
    def build(context: BuildContext) -> Any:
        paths = GenomicPaths(
            table_dir=GENOMIC_TABLES_DIR,
            figure_dir=context.figure_dir,
        )
        return build_func(paths)

    return ArtifactBuilder("genomic_networks", name, "figure", build, module)


def _genomic_table_builder(
    name: str,
    build_func: Callable[[GenomicPaths], Any],
) -> ArtifactBuilder:
    module = _validate_numbered_name(name, build_func)
    def build(context: BuildContext) -> Any:
        paths = GenomicPaths(
            table_dir=GENOMIC_TABLES_DIR,
            figure_dir=context.figure_dir,
            publication_table_dir=context.table_dir,
        )
        return build_func(paths)

    return ArtifactBuilder("genomic_networks", name, "table", build, module)


def _sse_figure_builder(
    name: str,
    build_func: Callable[[SSEPaths], Any],
) -> ArtifactBuilder:
    module = _validate_numbered_name(name, build_func)
    def build(context: BuildContext) -> Any:
        paths = SSEPaths(
            table_dir=SSE_TABLE_DIR,
            figure_dir=context.figure_dir,
            bayesian_result_dir=BAYESIAN_OUTPUT_DIR,
            result_table_dir=SSE_RESULT_TABLE_DIR,
        )
        return build_func(paths)

    return ArtifactBuilder("sse_detection", name, "figure", build, module)


def _sse_table_builder(
    name: str,
    build_func: Callable[[SSEPaths], Any],
) -> ArtifactBuilder:
    module = _validate_numbered_name(name, build_func)
    def build(context: BuildContext) -> Any:
        paths = SSEPaths(
            table_dir=SSE_TABLE_DIR,
            figure_dir=context.figure_dir,
            publication_table_dir=context.table_dir,
            bayesian_result_dir=BAYESIAN_OUTPUT_DIR,
            result_table_dir=SSE_RESULT_TABLE_DIR,
        )
        return build_func(paths)

    return ArtifactBuilder("sse_detection", name, "table", build, module)


def figure_builders() -> tuple[ArtifactBuilder, ...]:
    """Return each individual figure builder in build order."""
    return (
        _surveillance_figure_builder(
            surveillance_fig01.FIGURE_NAME, surveillance_fig01.build
        ),
        _surveillance_figure_builder(
            surveillance_fig02.FIGURE_NAME, surveillance_fig02.build
        ),
        _genomic_figure_builder(genomic_fig01.FIGURE_NAME, genomic_fig01.build),
        _genomic_figure_builder(genomic_fig02.FIGURE_NAME, genomic_fig02.build),
        _genomic_figure_builder(genomic_fig03.FIGURE_NAME, genomic_fig03.build),
        _genomic_figure_builder(genomic_fig04.FIGURE_NAME, genomic_fig04.build),
        _genomic_figure_builder(genomic_fig05.FIGURE_NAME, genomic_fig05.build),
        _genomic_figure_builder(genomic_fig06.FIGURE_NAME, genomic_fig06.build),
        _genomic_figure_builder(genomic_fig07.FIGURE_NAME, genomic_fig07.build),
        _genomic_figure_builder(genomic_fig08.FIGURE_NAME, genomic_fig08.build),
        _genomic_figure_builder(genomic_fig09.FIGURE_NAME, genomic_fig09.build),
        _genomic_figure_builder(genomic_fig10.FIGURE_NAME, genomic_fig10.build),
        _sse_figure_builder(sse_fig01.FIGURE_NAME, sse_fig01.build),
        _sse_figure_builder(sse_fig02.FIGURE_NAME, sse_fig02.build),
        _sse_figure_builder(sse_fig03.FIGURE_NAME, sse_fig03.build),
        _sse_figure_builder(sse_fig04.FIGURE_NAME, sse_fig04.build),
        _sse_figure_builder(sse_fig05.FIGURE_NAME, sse_fig05.build),
        _sse_figure_builder(sse_fig06.FIGURE_NAME, sse_fig06.build),
        _sse_figure_builder(sse_fig07.FIGURE_NAME, sse_fig07.build),
        _sse_figure_builder(sse_fig08.FIGURE_NAME, sse_fig08.build),
        _sse_figure_builder(sse_fig09.FIGURE_NAME, sse_fig09.build),
        _sse_figure_builder(sse_fig10.FIGURE_NAME, sse_fig10.build),
        _sse_figure_builder(sse_fig11.FIGURE_NAME, sse_fig11.build),
        _sse_figure_builder(sse_fig13.FIGURE_NAME, sse_fig13.build),
        _sse_figure_builder(sse_fig14.FIGURE_NAME, sse_fig14.build),
        _sse_figure_builder(sse_fig15.FIGURE_NAME, sse_fig15.build),
        _sse_figure_builder(sse_fig20.FIGURE_NAME, sse_fig20.build),
        _sse_figure_builder(sse_fig21.FIGURE_NAME, sse_fig21.build),
        _sse_figure_builder(sse_fig22.FIGURE_NAME, sse_fig22.build),
        _sse_figure_builder(sse_fig23.FIGURE_NAME, sse_fig23.build),
        _sse_figure_builder(sse_fig24.FIGURE_NAME, sse_fig24.build),
        _sse_figure_builder(sse_fig25.FIGURE_NAME, sse_fig25.build),
        _sse_figure_builder(sse_fig26.FIGURE_NAME, sse_fig26.build),
        _sse_figure_builder(sse_fig27.FIGURE_NAME, sse_fig27.build),
        _sse_figure_builder(sse_fig28.FIGURE_NAME, sse_fig28.build),
        _sse_figure_builder(sse_fig29.FIGURE_NAME, sse_fig29.build),
    )


def table_builders() -> tuple[ArtifactBuilder, ...]:
    """Return each individual table builder in build order."""
    return (
        _surveillance_table_builder(
            surveillance_tab01.TABLE_NAME, surveillance_tab01.build
        ),
        _surveillance_table_builder(
            surveillance_tab02.TABLE_NAME, surveillance_tab02.build
        ),
        _surveillance_table_builder(
            surveillance_tab03.TABLE_NAME, surveillance_tab03.build
        ),
        _surveillance_table_builder(
            surveillance_tab04.TABLE_NAME, surveillance_tab04.build
        ),
        _surveillance_table_builder(
            surveillance_tab05.TABLE_NAME, surveillance_tab05.build
        ),
        _surveillance_table_builder(
            surveillance_tab06.TABLE_NAME, surveillance_tab06.build
        ),
        _surveillance_table_builder(
            surveillance_tab07.TABLE_NAME, surveillance_tab07.build
        ),
        _surveillance_table_builder(
            surveillance_tab08.TABLE_NAME, surveillance_tab08.build
        ),
        _genomic_table_builder(genomic_tab01.TABLE_NAME, genomic_tab01.build),
        _genomic_table_builder(genomic_tab02.TABLE_NAME, genomic_tab02.build),
        _genomic_table_builder(genomic_tab03.TABLE_NAME, genomic_tab03.build),
        _genomic_table_builder(genomic_tab04.TABLE_NAME, genomic_tab04.build),
        _genomic_table_builder(genomic_tab05.TABLE_NAME, genomic_tab05.build),
        _genomic_table_builder(genomic_tab06.TABLE_NAME, genomic_tab06.build),
        _genomic_table_builder(genomic_tab07.TABLE_NAME, genomic_tab07.build),
        _sse_table_builder(sse_tab01.TABLE_NAME, sse_tab01.build),
        _sse_table_builder(sse_tab02.TABLE_NAME, sse_tab02.build),
        _sse_table_builder(sse_tab03.TABLE_NAME, sse_tab03.build),
        _sse_table_builder(sse_tab04.TABLE_NAME, sse_tab04.build),
        _sse_table_builder(sse_tab05.TABLE_NAME, sse_tab05.build),
        _sse_table_builder(sse_tab06.TABLE_NAME, sse_tab06.build),
        _sse_table_builder(sse_tab07.TABLE_NAME, sse_tab07.build),
        _sse_table_builder(sse_tab08.TABLE_NAME, sse_tab08.build),
        _sse_table_builder(sse_tab09.TABLE_NAME, sse_tab09.build),
        _sse_table_builder(sse_tab10.TABLE_NAME, sse_tab10.build),
    )


def publication_assets() -> list[dict[str, Any]]:
    """Describe the numbered PDF/PNG and LaTeX exports for migration and thesis sync."""
    assets = []
    for builder in (*figure_builders(), *table_builders()):
        # Surveillance table scripts produce analysis-local CSV/parquet data.
        if builder.kind == "table" and builder.domain == "surveillance":
            continue
        assets.append({
            "domain": builder.domain,
            "kind": builder.kind,
            "module": builder.module,
            "script": builder.script,
            "stem": builder.name,
            "legacy_stem": builder.legacy_name,
            "directory": "figures" if builder.kind == "figure" else "tables",
            "formats": ["pdf", "png"] if builder.kind == "figure" else ["tex"],
        })
    if len({row["stem"] for row in assets}) != len(assets):
        raise ValueError("Publication export stems must be unique across analysis domains")
    return assets


def write_asset_manifest(path: Path = RESULTS_DIR / "asset_manifest.json") -> Path:
    """Record producer identities independently of manuscript figure/table numbering."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"schema_version": 1, "assets": publication_assets()}, indent=2) + "\n")
    return path


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
            raise KeyError(f"Unknown builder {name!r}. Available builders: {available}")
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
