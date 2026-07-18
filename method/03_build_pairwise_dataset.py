#!/usr/bin/env python3
"""
Process pairwise genetic distances into an extended, distributed Parquet dataset.

Transforms raw TN93 CSVs into a Zstd-compressed Parquet dataset, enriched with 
calculated SNP distances, temporal distances (days), and EpiLink compatibility weights 
to facilitate large-scale epidemiological network analysis.

Input:
    - CSVs: `<WindowID>_<Lineage>_<N_Unique>.csv` (columns: ID1, ID2, Distance).
    - Metadata: Parquet file with `sequence_id` and `collection_date`.

Output Schema (Parquet Dataset):
    window_id, pango_lineage, nunique_sequences, id1, id2, tn93_distance, 
    snp_distance, temporal_distance, epilink_compatibility

Usage:
    $ python3 method/03_build_pairwise_dataset.py
    $ python3 method/03_build_pairwise_dataset.py --config config.yaml --root /path/to/repo

Example dataset usage with Pandas and igraph for clustering:
    Pandas + igraph (Memory-efficient graph construction & clustering)
    >>> import pandas as pd
    >>> import igraph as ig
    >>> edges = pd.read_parquet(
    ...     "data/dataset_dir/",
    ...     columns=["id1", "id2", "epilink_compatibility"],      # Projection pushdown
    ...     filters=[("epilink_compatibility", ">", 0.001)]       # Predicate pushdown
    ... )
    >>> g = ig.Graph.TupleList(edges.itertuples(index=False), directed=False, edge_attrs=['weight'])
    >>> communities = g.community_leiden(weights="weight", resolution=0.2, n_iterations=10)
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from epilink_wrapper import estimate_epilink_compatibility


def load_config(path: Path) -> dict:
    """Load the YAML pipeline configuration file."""
    with open(path) as f:
        return yaml.safe_load(f)


def load_metadata_dates(path: Path) -> pd.Series:
    """Load sequence collection dates indexed by sequence ID."""
    df = pd.read_parquet(path, columns=["sequence_id", "collection_date"])
    df = df.drop_duplicates("sequence_id").set_index("sequence_id")
    ser = pd.to_datetime(df["collection_date"], errors="coerce").dropna()
    return ser


def parse_group_label(stem: str) -> tuple:
    """Extract Surveillance/Window ID, Lineage, and N-unique from filename."""
    # Matches format: "W001_BA.2_124"
    m = re.match(r"^([^_]+)_([^_]+)_(\d+)$", stem)
    if m:
        return m.group(1), m.group(2), int(m.group(3))
    return None, None, None


def process_pairwise_dataset(
    input_dir: Path,
    metadata_path: Path,
    output_dataset_dir: Path,
    alignment_length: int,
    force: bool = False,
):
    """Converts a directory of TN93 CSVs into a Parquet dataset of pairwise edges."""
    if not input_dir.exists():
        raise FileNotFoundError(f"TN93 results directory not found: {input_dir}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata parquet not found: {metadata_path}")

    output_dataset_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Loading metadata dates...")
    date_map = load_metadata_dates(metadata_path)
    csv_files = sorted(input_dir.glob("*.csv"))

    if not csv_files:
        raise ValueError(f"No TN93 CSV files found in {input_dir}")

    logging.info(f"Found {len(csv_files)} CSV files to process.")

    for tn93_csv in csv_files:
        stem = tn93_csv.stem
        window_id, lineage, n_unique = parse_group_label(stem)

        out_path = output_dataset_dir / f"{stem}.parquet"

        # Skip if already processed (useful if the script crashes and you need to resume)
        if out_path.exists() and not force:
            continue

        # 1. Load TN93 distances
        df = pd.read_csv(tn93_csv)
        df = df.dropna(subset=["ID1", "ID2", "Distance"]).copy()
        df["Distance"] = pd.to_numeric(df["Distance"], errors="coerce")
        df = df[df["Distance"].between(0, 1, inclusive="both")]
        df = df[df["ID1"] != df["ID2"]]

        # Deduplicate undirected pairs
        a, b = df["ID1"].astype(str), df["ID2"].astype(str)
        key = pd.DataFrame({"a": np.minimum(a, b), "b": np.maximum(a, b)})
        df = df.loc[~key.duplicated(keep="first")].copy()

        # 2. Calculate SNP Distance
        df["snp_distance"] = np.rint(df["Distance"] * alignment_length).astype(float)

        # 3. Calculate Temporal Distance
        df["date1"] = pd.to_datetime(df["ID1"].map(date_map), errors="coerce")
        df["date2"] = pd.to_datetime(df["ID2"].map(date_map), errors="coerce")

        # Drop rows missing dates so we can compute temporal distance
        df = df.dropna(subset=["date1", "date2"])
        df["temporal_distance"] = (df["date1"] - df["date2"]).abs().dt.days.astype(float)

        if df.empty:
            logging.info(f"{stem}: No valid pairs left after date filtering. Skipping.")
            continue

        # 4. Calculate EpiLink Compatibility
        try:
            df["epilink_compatibility"] = estimate_epilink_compatibility(
                genetic_distance=df["snp_distance"].to_numpy(dtype=float),
                temporal_distance=df["temporal_distance"].to_numpy(dtype=float),
            )
        except Exception as exc:
            logging.warning(f"{stem}: epilink failed ({exc}). Filling with NaNs.")
            df["epilink_compatibility"] = np.nan

        # 5. Format to final schema
        df = df.rename(columns={
            "ID1": "id1",
            "ID2": "id2",
            "Distance": "tn93_distance"
        })

        df["window_id"] = window_id
        df["pango_lineage"] = lineage
        df["nunique_sequences"] = n_unique

        # Select and reorder columns
        final_cols = [
            "window_id", "pango_lineage", "nunique_sequences",
            "id1", "id2", "tn93_distance", "snp_distance",
            "temporal_distance", "epilink_compatibility"
        ]
        df = df[final_cols]

        # 6. Save directly to the Parquet Dataset directory
        df.to_parquet(
            out_path,
            engine="pyarrow",
            compression="zstd",
            index=False
        )
        logging.info(f"Processed {stem}: Wrote {len(df):,} edges.")


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Override TN93 results directory from config.",
    )
    ap.add_argument(
        "--metadata",
        type=Path,
        default=None,
        help="Override metadata parquet from config.",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Override pairwise parquet output directory from config.",
    )
    ap.add_argument(
        "--alignment-length",
        type=int,
        default=None,
        help="Override pipeline.alignment_length from config.",
    )
    ap.add_argument("--force", action="store_true", help="Rebuild existing output parquets.")
    ap.add_argument("--log-level", default="INFO")
    return ap.parse_args()


def resolve_path(root: Path, override: Path | None, configured: str | Path) -> Path:
    """Resolve an optional CLI override or config path against the project root."""
    path = override if override is not None else Path(configured)
    return path if path.is_absolute() else root / path


def main() -> int:
    """Resolve config paths and build the pairwise parquet dataset."""
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    cfg = load_config(args.root / args.config)
    pipe = cfg["pipeline"]
    proc = cfg["data"]["processed"]

    pairwise_default = proc.get(
        "pairwise_distances_dataset",
        "data/processed/pairwise_distances_dataset",
    )
    input_dir = resolve_path(args.root, args.input_dir, proc["tn93_results_dir"])
    metadata_path = resolve_path(args.root, args.metadata, proc["metadata"])
    output_dataset_dir = resolve_path(args.root, args.output_dir, pairwise_default)
    alignment_length = (
        args.alignment_length
        if args.alignment_length is not None
        else pipe["alignment_length"]
    )

    logging.info("TN93 input dir: %s", input_dir)
    logging.info("Metadata path: %s", metadata_path)
    logging.info("Pairwise output dir: %s", output_dataset_dir)
    logging.info("Alignment length: %s", alignment_length)

    try:
        process_pairwise_dataset(
            input_dir,
            metadata_path,
            output_dataset_dir,
            alignment_length=alignment_length,
            force=args.force,
        )
    except Exception as exc:
        logging.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
