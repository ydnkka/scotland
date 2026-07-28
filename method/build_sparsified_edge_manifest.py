"""Build the sparse pairwise edge-count manifest.

Run from the Scotland repository root:

    python3 method/build_sparsified_edge_manifest.py

The output is used by ``chapter_analyses.genomic_networks.build_mixing`` to schedule large
pairwise scans by the number of edges retained after sparsification.
"""

from __future__ import annotations

import argparse
import logging
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq
import yaml

LOGGER = logging.getLogger(__name__)


def load_config(path: Path) -> dict:
    """Load the YAML pipeline configuration file."""
    with open(path) as f:
        return yaml.safe_load(f)


def resolve_path(root: Path, override: Path | None, configured: str | Path) -> Path:
    """Resolve an optional CLI override or config path against the project root."""
    path = override if override is not None else Path(configured)
    return path if path.is_absolute() else root / path


def _normalise_window(value: str) -> str:
    value = str(value).strip()
    upper = value.upper()
    if upper.startswith("W") and upper[1:].isdigit():
        return f"W{int(upper[1:]):03d}"
    if upper.isdigit():
        return f"W{int(upper):03d}"
    return value


def _metadata_from_pairwise_path(path: Path) -> dict[str, object]:
    parts = path.stem.split("_")
    if len(parts) < 3 or not parts[-1].isdigit():
        raise ValueError(
            "Pairwise filename must look like {window}_{lineage}_{count}: "
            f"{path.name}"
        )

    window_id = _normalise_window(parts[0])
    lineage = "_".join(parts[1:-1])
    window_match = re.search(r"(\d+)$", window_id)
    window_idx = int(window_match.group(1)) if window_match else math.nan

    return {
        "window_idx": window_idx,
        "window_id": window_id,
        "pango_lineage": lineage,
        "pairwise_stem": path.stem,
        "nunique_sequences": int(parts[-1]),
    }


def _count_sparse_edges(
    path: Path,
    *,
    threshold: float,
    score_column: str,
    batch_size: int,
) -> tuple[int, int]:
    parquet_file = pq.ParquetFile(path)
    total_rows = int(parquet_file.metadata.num_rows)
    sparse_edges = 0

    for batch in parquet_file.iter_batches(
        columns=[score_column],
        batch_size=batch_size,
    ):
        scores = batch.column(0)
        retained = pc.greater(scores, threshold)  # type: ignore
        retained_count = pc.sum(  # type: ignore
            pc.fill_null(retained, False),
        ).as_py()
        sparse_edges += int(retained_count or 0)

    return total_rows, sparse_edges


def build_manifest(
    *,
    pairwise_dir: Path,
    threshold: float,
    score_column: str,
    batch_size: int,
    max_files: int | None = None,
) -> pd.DataFrame:
    """Return one sparse edge-count row per pairwise parquet file."""
    files = sorted(pairwise_dir.glob("*.parquet"), key=lambda path: path.stem)
    if max_files is not None:
        files = files[:max_files]
    if not files:
        raise FileNotFoundError(f"No pairwise parquet files found in {pairwise_dir}")

    rows: list[dict[str, object]] = []
    for idx, path in enumerate(files, start=1):
        if idx == 1 or idx % 100 == 0 or idx == len(files):
            LOGGER.info("Scanning pairwise file %s/%s: %s", idx, len(files), path.name)

        row = _metadata_from_pairwise_path(path)
        total_edges, sparse_edges = _count_sparse_edges(
            path,
            threshold=threshold,
            score_column=score_column,
            batch_size=batch_size,
        )
        row.update(
            {
                "total_pairwise_edges": total_edges,
                "sparse_edges": sparse_edges,
                "sparse_edge_fraction": (
                    sparse_edges / total_edges if total_edges else np.nan
                ),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        ["window_idx", "pango_lineage", "pairwise_stem"],
        kind="mergesort",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--pairwise-dir",
        type=Path,
        default=None,
        help="Override pairwise parquet directory from config.",
    )
    parser.add_argument(
        "--out-path",
        type=Path,
        default=None,
        help="Override sparse edge manifest parquet path from config.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Override pipeline.sparsification from config.",
    )
    parser.add_argument(
        "--score-column",
        default="epilink_compatibility",
        help="Pairwise score column used for sparsification.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1_000_000,
        help="Rows per parquet scan batch. Default: 1,000,000.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Development cap on the number of pairwise files to scan.",
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

    cfg = load_config(args.root / args.config)
    pipe = cfg["pipeline"]
    proc = cfg["data"]["processed"]

    pairwise_dir = resolve_path(
        args.root,
        args.pairwise_dir,
        proc.get(
            "pairwise_distances_dataset",
            "data/processed/pairwise_distances_dataset",
        ),
    )
    out_path = resolve_path(
        args.root,
        args.out_path,
        proc.get(
            "sparsified_edge_manifest",
            "data/processed/sparsified_edge_counts_by_window_lineage.parquet",
        ),
    )
    threshold = args.threshold if args.threshold is not None else pipe["sparsification"]

    if threshold < 0:
        raise SystemExit("--threshold must be non-negative.")
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be positive.")

    LOGGER.info("Pairwise input dir: %s", pairwise_dir)
    LOGGER.info("Manifest output path: %s", out_path)
    LOGGER.info("Sparsification threshold: %g", threshold)

    manifest = build_manifest(
        pairwise_dir=pairwise_dir,
        threshold=threshold,
        score_column=args.score_column,
        batch_size=args.batch_size,
        max_files=args.max_files,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_parquet(out_path, index=False)
    
    LOGGER.info(
        "Wrote %s rows at threshold %g to %s",
        f"{len(manifest):,}",
        threshold,
        str(out_path)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
