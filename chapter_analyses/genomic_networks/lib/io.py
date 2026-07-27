"""Input/output helpers for Chapter 4 analysis products."""

from __future__ import annotations

import sys
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

from .config import (
    ANALYSIS_COLUMNS,
    ANALYSIS_RESOLUTION,
    FIGURES_DIR,
    INTERMEDIATE_DIR,
    PROJECT_ROOT,
    RESULTS_DIR,
    SPARSIFICATION_THRESHOLD,
    TABLES_DIR,
)

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import load_analysis_columns, load_pairwise_edges


def ensure_results_dirs() -> None:
    """Create the standard Chapter 4 output directories."""
    for path in (RESULTS_DIR, TABLES_DIR, FIGURES_DIR, INTERMEDIATE_DIR):
        path.mkdir(parents=True, exist_ok=True)


def load_sequence_data(
    columns: Sequence[str] | None = None,
    *,
    resolution: float | None = ANALYSIS_RESOLUTION,
    window_stride: int | None = None,
    window_offset: int = 0,
    renumber_windows: bool = True,
    weighted_simd: bool = True,
) -> pd.DataFrame:
    """Load sequence-window rows for Chapter 4.

    The defaults match the intended main analysis: good-QC genomes at Leiden
    resolution 0.3 with population-weighted SIMD groups. Policy fields are
    persisted in the analysis dataset. Passing ``window_stride=2`` gives the
    alternate-window graph input used by the Chapter 5 detector.
    """
    return load_analysis_columns(
        columns or ANALYSIS_COLUMNS,
        resolution=resolution,
        window_stride=window_stride,
        window_offset=window_offset,
        renumber_windows=renumber_windows,
        weighted_simd=weighted_simd,
    )


def load_pairwise_compatibility_edges(
    *,
    windows: str | int | Iterable[str | int] | None = None,
    clades: str | Iterable[str] | None = None,
    pango_lineages: str | Iterable[str] | None = None,
    sequence_ids: str | Iterable[str] | None = None,
    compatibility_threshold: float | None = SPARSIFICATION_THRESHOLD,
    columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Load sparse pairwise compatibility edges for network mixing summaries."""
    selected = list(
        columns
        or (
            "window_id",
            "pango_lineage",
            "id1",
            "id2",
            "epilink_compatibility",
        )
    )
    return load_pairwise_edges(
        selected,
        windows=windows,
        clades=clades,
        pango_lineages=pango_lineages,
        sequence_ids=sequence_ids,
        compatibility_threshold=compatibility_threshold,
    )


def table_path(
    name: str, *, suffix: str = ".parquet", table_dir: Path = TABLES_DIR
) -> Path:
    """Return a standard table path under ``results/tables``."""
    clean = name.removesuffix(".parquet").removesuffix(".csv")
    return table_dir / f"{clean}{suffix}"


def write_table(
    df: pd.DataFrame,
    name: str,
    *,
    table_dir: Path = TABLES_DIR,
    formats: Sequence[str] = ("parquet", "csv"),
    index: bool = False,
) -> dict[str, Path]:
    """Write a table in one or more standard formats."""
    table_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    for fmt in formats:
        fmt = fmt.lower().lstrip(".")
        path = table_path(name, suffix=f".{fmt}", table_dir=table_dir)
        if fmt == "parquet":
            df.to_parquet(path, index=index)
        elif fmt == "csv":
            df.to_csv(path, index=index)
        else:
            raise ValueError(f"Unsupported table format: {fmt}")
        written[fmt] = path

    return written


def read_table(name: str, *, table_dir: Path = TABLES_DIR) -> pd.DataFrame:
    """Read a standard Chapter 4 output table, preferring parquet."""
    parquet_path = table_path(name, suffix=".parquet", table_dir=table_dir)
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)

    csv_path = table_path(name, suffix=".csv", table_dir=table_dir)
    if csv_path.exists():
        return pd.read_csv(csv_path)

    raise FileNotFoundError(f"No table named {name!r} found in {table_dir}")
