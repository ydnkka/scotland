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
    $ python3 method/create_pairwise_dataset.py

Example dataset usage with Pandas and igraph for clustering:
    Pandas + igraph (Memory-efficient graph construction & clustering)
    >>> import pandas as pd
    >>> import igraph as ig
    >>> edges = pd.read_parquet(
    ...     "data/dataset_dir/",
    ...     columns=["id1", "id2", "epilink_compatibility"],      # Projection pushdown
    ...     filters=[("epilink_compatibility", ">", 0.0001)]      # Predicate pushdown
    ... )
    >>> g = ig.Graph.TupleList(edges.itertuples(index=False), directed=False, edge_attrs=['weight'])
    >>> communities = g.community_leiden(weights="weight", resolution=0.2, n_iterations=10)
"""

from __future__ import annotations


import re
import logging
from pathlib import Path
import numpy as np
import pandas as pd


from epilink_wrapper import estimate_epilink_compatibility

ALIGNMENT_LENGTH_DEFAULT = 29903

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
    alignment_length: int = ALIGNMENT_LENGTH_DEFAULT
):
    """Converts a directory of TN93 CSVs into a Parquet dataset of pairwise edges."""
    
    # Create the dataset directory
    output_dataset_dir.mkdir(parents=True, exist_ok=True)
    
    logging.info("Loading metadata dates...")
    date_map = load_metadata_dates(metadata_path)
    
    csv_files = list(input_dir.glob("*.csv"))
    logging.info(f"Found {len(csv_files)} CSV files to process.")

    for tn93_csv in csv_files:
        stem = tn93_csv.stem
        window_id, lineage, n_unique = parse_group_label(stem)
        
        if not window_id:
            logging.warning(f"Skipping {tn93_csv.name}: Filename format mismatch.")
            continue
            
        out_path = output_dataset_dir / f"{stem}.parquet"
        
        # Skip if already processed (useful if the script crashes and you need to resume)
        if out_path.exists():
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

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    INPUT_DIR = Path("data/processed/group_fastas/tn93_results")
    METADATA_PATH = Path("data/processed/scotland_sequence_metadata.parquet")
    OUTPUT_DATASET_DIR = Path("data/processed/pairwise_distances_dataset")
    
    process_pairwise_dataset(INPUT_DIR, METADATA_PATH, OUTPUT_DATASET_DIR)
