"""Loaders for the parquet outputs produced by ``sse_detection.ipynb``."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


HIGH_PRIORITY_CANDIDATE_TIERS = frozenset(
    {
        "high_priority_both_axes",
        "high_priority_burst",
        "high_priority_burden",
    }
)

_DEFAULT_TABLE_NAMES = (
    "cluster_table",
    "edge_table",
)


@dataclass
class SseOutputs:
    """Container for the static parquet outputs of the SSE pipeline."""

    cluster_table: pd.DataFrame
    edge_table: pd.DataFrame
    output_dir: Path

    @property
    def node_stats(self) -> pd.DataFrame:
        """Compatibility alias for association code written against node stats."""
        return self.cluster_table

    @property
    def candidates(self) -> pd.DataFrame:
        """Subset of high-priority burst/burden candidate nodes."""
        if "candidate_tier" not in self.cluster_table.columns:
            raise KeyError("cluster_table needs 'candidate_tier'")
        return self.cluster_table.loc[
            self.cluster_table["candidate_tier"].isin(HIGH_PRIORITY_CANDIDATE_TIERS)
        ].copy()


def _coerce_datetime(df: pd.DataFrame, cols: Iterable[str]) -> None:
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_datetime(df[col], errors="coerce")


def load_sse_outputs(
    output_dir: Path | str,
    *,
    tables: Iterable[str] = _DEFAULT_TABLE_NAMES,
) -> SseOutputs:
    """Load the SSE-pipeline parquet outputs from ``output_dir``.

    Parameters
    ----------
    output_dir
        Directory containing the parquet files. Typically
        ``<project_root>/sse_detection/results/sse_outputs``.
    tables
        Override which tables to load (mostly useful for testing).
    """
    output_dir = Path(output_dir)
    if not output_dir.is_dir():
        raise FileNotFoundError(f"SSE output directory not found: {output_dir}")

    found: dict[str, pd.DataFrame] = {}
    missing: list[str] = []

    for name in tables:
        path = output_dir / f"{name}.parquet"
        if not path.exists():
            missing.append(name)
            continue
        found[name] = pd.read_parquet(path)

    if missing:
        raise FileNotFoundError(
            "Missing expected SSE output parquet files: " + ", ".join(missing)
        )

    cluster_table = found["cluster_table"]
    edge_table = found["edge_table"]

    _coerce_datetime(
        cluster_table,
        [
            "wn_start_date",
            "wn_mid_date",
            "wn_end_date",
            "first_collection_date",
            "last_collection_date",
        ],
    )
    if (
        "who_voc" in cluster_table.columns
        and "who_voc_plot" not in cluster_table.columns
    ):
        cluster_table["who_voc_plot"] = cluster_table["who_voc"].fillna("None")

    return SseOutputs(
        cluster_table=cluster_table,
        edge_table=edge_table,
        output_dir=output_dir,
    )
