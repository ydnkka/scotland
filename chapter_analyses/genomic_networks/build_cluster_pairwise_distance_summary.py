"""Summarise pairwise distances among non-singleton cluster sequences.

For each selected rolling-window/Pango-lineage group, this script finds
non-singleton clusters at the Chapter 4 default Leiden resolution, combines the
selected cluster members, loads the corresponding physical pairwise parquet
file, and summarises SNP and temporal pairwise distances.

Example development run:

    python -m chapter_analyses.genomic_networks.build_cluster_pairwise_distance_summary \
        --windows W080 W081 --max-clusters-per-window-lineage 25
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import load_analysis_columns, load_pairwise_edges

from .lib.config import (
    ANALYSIS_RESOLUTION,
    PROJECT_ROOT,
    TABLES_DIR,
)
from .lib.io import ensure_results_dirs, write_table

LOGGER = logging.getLogger(__name__)

PAIRWISE_DATASET_DIR = PROJECT_ROOT / "data/processed/pairwise_distances_dataset"
DEFAULT_TABLE_NAME = "cluster_pairwise_distance_summary"
SEQUENCE_COLUMNS = [
    "window_id",
    "window_idx",
    "sequence_id",
    "cluster_id",
    "cluster_size",
    "pango_lineage",
]
PAIRWISE_COLUMNS = [
    "id1",
    "id2",
    "snp_distance",
    "temporal_distance",
]
OUTPUT_COLUMNS = [
    "window_idx",
    "window_id",
    "pango_lineage",
    "pairwise_stem",
    "status",
    "n_selected_clusters",
    "n_selected_sequences",
    "n_selected_possible_pairs",
    "n_pairwise_rows",
    "snp_distance_median",
    "snp_distance_q25",
    "snp_distance_q75",
    "snp_distance_iqr",
    "temporal_distance_median",
    "temporal_distance_q25",
    "temporal_distance_q75",
    "temporal_distance_iqr",
]
GROUP_COLUMNS = ["window_id", "window_idx", "pango_lineage"]
CLUSTER_COLUMNS = [*GROUP_COLUMNS, "cluster_id"]


@dataclass(frozen=True)
class PairwiseMetadata:
    window_id: str
    window_idx: int | float
    pango_lineage: str
    pairwise_stem: str
    nunique_sequences: int


def _normalise_window(value: Any) -> str:
    text = str(value).strip()
    upper = text.upper()
    if upper.startswith("W") and upper[1:].isdigit():
        return f"W{int(upper[1:]):03d}"
    if upper.isdigit():
        return f"W{int(upper):03d}"
    return text


def _normalise_windows(values: Iterable[Any] | None) -> list[str] | None:
    if values is None:
        return None
    return [_normalise_window(value) for value in values]


def _pairwise_metadata_from_path(path: Path) -> PairwiseMetadata:
    parts = path.stem.split("_")
    if len(parts) < 3 or not parts[-1].isdigit():
        raise ValueError(
            "Pairwise filename must look like {window}_{lineage}_{count}: "
            f"{path.name}"
        )
    window_id = _normalise_window(parts[0])
    window_idx = int(window_id[1:]) if window_id[1:].isdigit() else np.nan
    return PairwiseMetadata(
        window_id=window_id,
        window_idx=window_idx,
        pango_lineage="_".join(parts[1:-1]),
        pairwise_stem=path.stem,
        nunique_sequences=int(parts[-1]),
    )


def _pairwise_file_lookup(pairwise_dir: Path) -> dict[tuple[str, str], Path]:
    files = sorted(pairwise_dir.glob("*.parquet"), key=lambda path: path.stem)
    if not files:
        raise FileNotFoundError(f"No pairwise parquet files found in {pairwise_dir}")

    lookup: dict[tuple[str, str], Path] = {}
    for path in files:
        metadata = _pairwise_metadata_from_path(path)
        lookup[(metadata.window_id, metadata.pango_lineage)] = path
    return lookup


def _load_sequence_rows(
    *,
    resolution: float,
    windows: Sequence[str] | None,
    pango_lineages: Sequence[str] | None,
    max_windows: int | None,
) -> pd.DataFrame:
    df = load_analysis_columns(
        SEQUENCE_COLUMNS,
        resolution=resolution,
    )
    df = df.dropna(subset=["window_id", "sequence_id", "cluster_id", "pango_lineage"])
    df["window_id"] = df["window_id"].map(_normalise_window)

    if windows is not None:
        df = df.loc[df["window_id"].isin(set(windows))]
    if pango_lineages:
        df = df.loc[df["pango_lineage"].isin(set(pango_lineages))]
    if max_windows is not None:
        retained_windows = sorted(df["window_id"].dropna().unique())[:max_windows]
        df = df.loc[df["window_id"].isin(retained_windows)]

    if df.empty:
        raise ValueError("No sequence rows remain after applying filters.")
    return df.reset_index(drop=True)


def _cluster_counts(sequence_rows: pd.DataFrame) -> pd.DataFrame:
    out = (
        sequence_rows.groupby(CLUSTER_COLUMNS, dropna=False)
        .agg(
            n_cluster_sequences=("sequence_id", "nunique"),
            cluster_size=("cluster_size", "first"),
        )
        .reset_index()
    )
    out["cluster_size"] = pd.to_numeric(out["cluster_size"], errors="coerce")
    return out


def _select_non_singleton_clusters(
    cluster_counts: pd.DataFrame,
    *,
    min_cluster_size: int,
    max_clusters_per_window_lineage: int | None,
    selection: str,
) -> pd.DataFrame:
    eligible = cluster_counts.loc[
        cluster_counts["n_cluster_sequences"].ge(min_cluster_size)
    ].copy()
    if eligible.empty:
        raise ValueError("No non-singleton clusters remain after applying filters.")

    ascending_size = selection == "smallest"
    if selection in {"largest", "smallest"}:
        eligible = eligible.sort_values(
            [
                "window_idx",
                "window_id",
                "pango_lineage",
                "n_cluster_sequences",
                "cluster_id",
            ],
            ascending=[True, True, True, ascending_size, True],
            kind="mergesort",
        )
    elif selection == "cluster_id":
        eligible = eligible.sort_values(
            ["window_idx", "window_id", "pango_lineage", "cluster_id"],
            kind="mergesort",
        )
    else:
        raise ValueError(f"Unsupported cluster selection method: {selection}")

    if max_clusters_per_window_lineage is None:
        return eligible.reset_index(drop=True)

    return (
        eligible.groupby(GROUP_COLUMNS, dropna=False, group_keys=False)
        .head(max_clusters_per_window_lineage)
        .reset_index(drop=True)
    )


def _group_selection_table(
    sequence_rows: pd.DataFrame,
    selected_clusters: pd.DataFrame,
    *,
    min_cluster_size: int,
    max_groups: int | None,
) -> tuple[pd.DataFrame, dict[tuple[str, str], dict[str, str]]]:
    eligible_counts = (
        _cluster_counts(sequence_rows)
        .loc[lambda x: x["n_cluster_sequences"].ge(min_cluster_size)]
        .groupby(GROUP_COLUMNS, dropna=False)
        .agg(
            n_eligible_clusters=("cluster_id", "nunique"),
            n_eligible_sequences=("n_cluster_sequences", "sum"),
        )
        .reset_index()
    )

    selected_sequences = sequence_rows.merge(
        selected_clusters[CLUSTER_COLUMNS],
        on=CLUSTER_COLUMNS,
        how="inner",
    )
    selected_counts = (
        selected_sequences.groupby(GROUP_COLUMNS, dropna=False)
        .agg(
            n_selected_clusters=("cluster_id", "nunique"),
            n_selected_sequences=("sequence_id", "nunique"),
        )
        .reset_index()
    )
    selected_possible_pairs = (
        selected_sequences.groupby(CLUSTER_COLUMNS, dropna=False)["sequence_id"]
        .nunique()
        .rename("n_cluster_sequences")
        .reset_index()
    )
    selected_possible_pairs["n_cluster_possible_pairs"] = (
        selected_possible_pairs["n_cluster_sequences"]
        * (selected_possible_pairs["n_cluster_sequences"] - 1)
        // 2
    )
    selected_possible_pairs = (
        selected_possible_pairs.groupby(GROUP_COLUMNS, dropna=False)[
            "n_cluster_possible_pairs"
        ]
        .sum()
        .rename("n_selected_possible_pairs")
        .reset_index()
    )
    selected_counts = selected_counts.merge(
        selected_possible_pairs,
        on=GROUP_COLUMNS,
        how="left",
    )
    groups = eligible_counts.merge(
        selected_counts,
        on=GROUP_COLUMNS,
        how="inner",
    ).sort_values(["window_idx", "window_id", "pango_lineage"], kind="mergesort")

    if max_groups is not None:
        groups = groups.head(max_groups).copy()
        keep = groups[["window_id", "pango_lineage"]].drop_duplicates()
        selected_sequences = selected_sequences.merge(
            keep,
            on=["window_id", "pango_lineage"],
            how="inner",
        )

    cluster_by_sequence = {
        (str(window_id), str(pango_lineage)): dict(
            zip(
                group["sequence_id"].astype(str),
                group["cluster_id"].astype(str),
                strict=True,
            )
        )
        for (window_id, pango_lineage), group in selected_sequences.groupby(
            ["window_id", "pango_lineage"],
            sort=False,
        )
    }
    return groups.reset_index(drop=True), cluster_by_sequence


def _distance_summary(values: pd.Series, prefix: str) -> dict[str, float | int]:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        return {
            f"{prefix}_distance_median": np.nan,
            f"{prefix}_distance_q25": np.nan,
            f"{prefix}_distance_q75": np.nan,
            f"{prefix}_distance_iqr": np.nan,
        }

    q25 = float(numeric.quantile(0.25))
    q75 = float(numeric.quantile(0.75))
    return {
        f"{prefix}_distance_median": float(numeric.median()),
        f"{prefix}_distance_q25": q25,
        f"{prefix}_distance_q75": q75,
        f"{prefix}_distance_iqr": q75 - q25,
    }


def _empty_distance_summary(prefix: str) -> dict[str, float | int]:
    return _distance_summary(pd.Series(dtype=float), prefix)


def _to_int(value: Any) -> int:
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return int(value)
    return int(str(value))


def _summarise_one_group(
    row: dict[Any, Any],
    *,
    selected_cluster_by_sequence: dict[str, str],
    pairwise_path: Path | None,
) -> dict[str, Any]:
    n_selected_clusters = _to_int(row["n_selected_clusters"])
    n_selected_sequences = _to_int(row["n_selected_sequences"])
    n_selected_possible_pairs = _to_int(row["n_selected_possible_pairs"])
    base = {
        "window_idx": row["window_idx"],
        "window_id": row["window_id"],
        "pango_lineage": row["pango_lineage"],
        "n_selected_clusters": n_selected_clusters,
        "n_selected_sequences": n_selected_sequences,
        "n_selected_possible_pairs": n_selected_possible_pairs,
        "pairwise_stem": np.nan,
        "n_pairwise_rows": 0,
    }

    if pairwise_path is None:
        return {
            **base,
            **_empty_distance_summary("snp"),
            **_empty_distance_summary("temporal"),
            "status": "missing_pairwise_file",
        }

    metadata = _pairwise_metadata_from_path(pairwise_path)
    edges = load_pairwise_edges(
        PAIRWISE_COLUMNS,
        compatibility_threshold=None,
        pairwise_dataset=pairwise_path,
    )
    base["pairwise_stem"] = metadata.pairwise_stem

    if edges.empty:
        return {
            **base,
            **_empty_distance_summary("snp"),
            **_empty_distance_summary("temporal"),
            "status": "empty_pairwise_file",
        }

    id1 = edges["id1"].astype(str)
    id2 = edges["id2"].astype(str)
    id1_cluster = id1.map(selected_cluster_by_sequence)
    id2_cluster = id2.map(selected_cluster_by_sequence)
    same_cluster = (
        id1_cluster.notna() & id2_cluster.notna() & id1_cluster.eq(id2_cluster)
    )
    selected_edges = edges.loc[same_cluster]
    base["n_pairwise_rows"] = int(len(selected_edges))

    if selected_edges.empty:
        return {
            **base,
            **_empty_distance_summary("snp"),
            **_empty_distance_summary("temporal"),
            "status": "no_selected_pairwise_rows",
        }

    return {
        **base,
        **_distance_summary(selected_edges["snp_distance"], "snp"),
        **_distance_summary(selected_edges["temporal_distance"], "temporal"),
        "status": "ok",
    }


def build_cluster_pairwise_distance_summary(
    *,
    windows: Sequence[str] | None,
    pango_lineages: Sequence[str] | None = None,
    max_windows: int | None = None,
    max_groups: int | None = None,
    min_cluster_size: int = 2,
    max_clusters_per_window_lineage: int | None = None,
    cluster_selection: str = "largest",
    resolution: float = ANALYSIS_RESOLUTION,
    pairwise_dir: Path = PAIRWISE_DATASET_DIR,
) -> pd.DataFrame:
    """Build the window-lineage pairwise distance summary table."""
    if min_cluster_size < 2:
        raise ValueError("min_cluster_size must be at least 2.")
    if (
        max_clusters_per_window_lineage is not None
        and max_clusters_per_window_lineage < 1
    ):
        raise ValueError("max_clusters_per_window_lineage must be at least 1.")
    if max_windows is not None and max_windows < 1:
        raise ValueError("max_windows must be at least 1.")
    if max_groups is not None and max_groups < 1:
        raise ValueError("max_groups must be at least 1.")

    sequence_rows = _load_sequence_rows(
        resolution=resolution,
        windows=windows,
        pango_lineages=pango_lineages,
        max_windows=max_windows,
    )
    cluster_counts = _cluster_counts(sequence_rows)
    selected_clusters = _select_non_singleton_clusters(
        cluster_counts,
        min_cluster_size=min_cluster_size,
        max_clusters_per_window_lineage=max_clusters_per_window_lineage,
        selection=cluster_selection,
    )
    groups, cluster_by_sequence = _group_selection_table(
        sequence_rows,
        selected_clusters,
        min_cluster_size=min_cluster_size,
        max_groups=max_groups,
    )
    pairwise_lookup = _pairwise_file_lookup(pairwise_dir)

    LOGGER.info(
        "Summarising pairwise distances for %s window-lineage groups",
        f"{len(groups):,}",
    )
    rows: list[dict[str, Any]] = []
    for idx, group_row in enumerate(groups.to_dict("records"), start=1):
        key = (str(group_row["window_id"]), str(group_row["pango_lineage"]))
        if idx == 1 or idx % 100 == 0 or idx == len(groups):
            LOGGER.info(
                "Processing %s/%s: %s %s",
                f"{idx:,}",
                f"{len(groups):,}",
                key[0],
                key[1],
            )
        rows.append(
            _summarise_one_group(
                group_row,
                selected_cluster_by_sequence=cluster_by_sequence[key],
                pairwise_path=pairwise_lookup.get(key),
            )
        )

    out = pd.DataFrame(rows)
    return out[OUTPUT_COLUMNS].sort_values(
        ["window_idx", "window_id", "pango_lineage"],
        kind="mergesort",
    ).reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--windows",
        nargs="*",
        help="Specific windows to process, e.g. W080 W081 or 80 81.",
    )
    parser.add_argument(
        "--all-windows",
        action="store_true",
        help="Process all retained analysis windows.",
    )
    parser.add_argument(
        "--pango-lineages",
        nargs="*",
        default=None,
        help="Optional Pango lineage filter.",
    )
    parser.add_argument(
        "--max-windows",
        type=int,
        default=None,
        help="Development cap on selected windows after filtering.",
    )
    parser.add_argument(
        "--max-groups",
        type=int,
        default=None,
        help="Development cap on selected window-lineage groups.",
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=2,
        help="Minimum within-window cluster sequence count. Default: 2.",
    )
    parser.add_argument(
        "--max-clusters-per-window-lineage",
        type=int,
        default=None,
        help=(
            "Optional cap on selected non-singleton clusters per window-lineage. "
            "Default keeps all eligible clusters."
        ),
    )
    parser.add_argument(
        "--cluster-selection",
        choices=("largest", "smallest", "cluster_id"),
        default="largest",
        help=(
            "How to choose clusters when --max-clusters-per-window-lineage is set. "
            "Default: largest."
        ),
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=ANALYSIS_RESOLUTION,
        help=f"Leiden resolution filter. Default: {ANALYSIS_RESOLUTION}.",
    )
    parser.add_argument(
        "--pairwise-dir",
        type=Path,
        default=PAIRWISE_DATASET_DIR,
        help=(
            "Directory of per-window-lineage pairwise parquet files. "
            f"Default: {PAIRWISE_DATASET_DIR}"
        ),
    )
    parser.add_argument(
        "--output-name",
        default=DEFAULT_TABLE_NAME,
        help=f"Output table stem under --table-dir. Default: {DEFAULT_TABLE_NAME}.",
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=TABLES_DIR,
        help=f"Output table directory. Default: {TABLES_DIR}.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=("parquet", "csv"),
        help="Output formats passed to write_table. Default: parquet csv.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print selected group counts and exit without loading pairwise files.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=args.log_level, format="%(levelname)s: %(message)s")

    if not args.all_windows and not args.windows:
        raise SystemExit("Specify --windows or --all-windows.")
    if args.all_windows and args.windows:
        raise SystemExit("Use either --windows or --all-windows, not both.")

    windows = None if args.all_windows else _normalise_windows(args.windows)

    if args.dry_run:
        sequence_rows = _load_sequence_rows(
            resolution=args.resolution,
            windows=windows,
            pango_lineages=args.pango_lineages,
            max_windows=args.max_windows,
        )
        selected_clusters = _select_non_singleton_clusters(
            _cluster_counts(sequence_rows),
            min_cluster_size=args.min_cluster_size,
            max_clusters_per_window_lineage=args.max_clusters_per_window_lineage,
            selection=args.cluster_selection,
        )
        groups, _ = _group_selection_table(
            sequence_rows,
            selected_clusters,
            min_cluster_size=args.min_cluster_size,
            max_groups=args.max_groups,
        )
        print(f"Selected window-lineage groups: {len(groups):,}")
        print(f"Selected clusters: {groups['n_selected_clusters'].sum():,}")
        print(f"Selected sequences: {groups['n_selected_sequences'].sum():,}")
        preview_cols = [
            "window_id",
            "pango_lineage",
            "n_eligible_clusters",
            "n_selected_clusters",
            "n_selected_sequences",
        ]
        print(groups[preview_cols].head(20).to_string(index=False))
        return 0

    ensure_results_dirs()
    summary = build_cluster_pairwise_distance_summary(
        windows=windows,
        pango_lineages=args.pango_lineages,
        max_windows=args.max_windows,
        max_groups=args.max_groups,
        min_cluster_size=args.min_cluster_size,
        max_clusters_per_window_lineage=args.max_clusters_per_window_lineage,
        cluster_selection=args.cluster_selection,
        resolution=args.resolution,
        pairwise_dir=args.pairwise_dir,
    )
    written = write_table(
        summary,
        args.output_name,
        table_dir=args.table_dir,
        formats=args.formats,
    )
    for fmt, path in written.items():
        LOGGER.info("Wrote %s output to %s", fmt, path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
