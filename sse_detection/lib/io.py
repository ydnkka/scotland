"""Loaders for the parquet outputs produced by ``sse_detection.ipynb``."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


_DEFAULT_TABLE_NAMES = (
    "node_stats",
    "edge_table",
    "candidate_review",
    "meta_summary",
)

_WEEKLY_COLUMNS = [
    "week",
    "new_sequences",
    "meta_cluster_id",
    "clade",
    "cc_size",
    "cc_size_prev",
    "norm_change",
    "is_sse",
]


@dataclass
class SseOutputs:
    """Container for the static parquet outputs of the SSE pipeline.

    ``weekly_growth`` is deliberately omitted from the bulk container because
    it has ~7M rows; load it separately via :func:`load_weekly_growth`.
    """

    node_stats: pd.DataFrame
    edge_table: pd.DataFrame
    candidate_review: pd.DataFrame
    meta_summary: pd.DataFrame
    output_dir: Path

    @property
    def candidates(self) -> pd.DataFrame:
        """Subset of ``node_stats`` flagged as ``sse_candidate``."""
        return self.node_stats.loc[self.node_stats["sse_candidate"]].copy()


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

    node_stats = found["node_stats"]
    edge_table = found["edge_table"]
    candidate_review = found["candidate_review"]
    meta_summary = found["meta_summary"]

    _coerce_datetime(
        node_stats,
        ["wn_mid_date", "first_collection_date", "last_collection_date"],
    )
    _coerce_datetime(
        candidate_review,
        ["wn_mid_date", "first_collection_date", "last_collection_date"],
    )
    _coerce_datetime(
        meta_summary,
        [
            "first_window_mid_date",
            "last_window_mid_date",
            "first_collection_date",
            "last_collection_date",
        ],
    )

    if "who_voc" in node_stats.columns:
        node_stats["who_voc_plot"] = node_stats["who_voc"].fillna("None")
    if "who_voc" in candidate_review.columns:
        candidate_review["who_voc_plot"] = candidate_review["who_voc"].fillna("None")
    if "who_voc" in meta_summary.columns:
        meta_summary["who_voc_plot"] = meta_summary["who_voc"].fillna("None")

    return SseOutputs(
        node_stats=node_stats,
        edge_table=edge_table,
        candidate_review=candidate_review,
        meta_summary=meta_summary,
        output_dir=output_dir,
    )


def load_weekly_growth(
    output_dir: Path | str,
    *,
    only_sse: bool = True,
    columns: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Load the weekly meta-cluster growth table.

    By default returns only the SSE-flagged weeks, which is what the
    plotting figures actually need. Pass ``only_sse=False`` to load the
    full ~7M-row table (use with care).
    """
    output_dir = Path(output_dir)
    path = output_dir / "weekly_growth.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Weekly growth parquet not found: {path}")

    cols = list(columns) if columns else _WEEKLY_COLUMNS

    if only_sse:
        try:
            weekly = pd.read_parquet(
                path, columns=cols, filters=[("is_sse", "=", True)]
            )
        except Exception:
            weekly = pd.read_parquet(path, columns=cols).query("is_sse")
    else:
        weekly = pd.read_parquet(path, columns=cols)

    weekly["week"] = pd.to_datetime(weekly["week"])
    return weekly
