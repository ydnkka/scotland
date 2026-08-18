"""Standard table I/O for surveillance outputs."""

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from .config import FIGURES_DIR, TABLES_DIR


def ensure_results_dirs() -> None:
    """Create the standard surveillance output directories."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)


def write_table(
    table: pd.DataFrame,
    name: str,
    *,
    table_dir: Path = TABLES_DIR,
    formats: Sequence[str] = ("csv", "parquet"),
) -> dict[str, Path]:
    """Write an index-free table in the requested standard formats."""
    output = table.copy()
    if output.index.name is not None:
        output = output.reset_index()
    elif not isinstance(output.index, pd.RangeIndex):
        output = output.reset_index(drop=True)

    table_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for fmt in formats:
        clean = fmt.lower().lstrip(".")
        path = table_dir / f"{name}.{clean}"
        if clean == "csv":
            output.to_csv(path, index=False)
        elif clean == "parquet":
            output.to_parquet(path, index=False)
        else:
            raise ValueError(f"Unsupported table format: {fmt}")
        paths[clean] = path
    return paths
