"""
wave_dates.py
-------------
Identifies epidemic wave start, peak, and end dates from a SARS-CoV-2
sequencing dataset by grouping pango_lineages into named variant waves
and computing weekly sequence counts.

A wave is considered:
  - Started  : first week where the variant exceeds THRESHOLD of all sequences
  - Peak     : week with the highest absolute sequence count
  - Ended    : last week where the variant exceeds THRESHOLD of all sequences

Usage:
    python wave_dates.py --input <path/to/dataset.parquet> [--threshold 0.05]
"""

import argparse
import pandas as pd


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

WAVE_GROUPS = {
    "B.1.177":  lambda l: l.startswith("B.1.177"),
    "Alpha":    lambda l: l == "B.1.1.7" or l.startswith("B.1.1.7."),
    "Delta":    lambda l: l.startswith("AY.") or l == "B.1.617.2",
    "BA.1":     lambda l: l.startswith("BA.1"),
    "BA.2":     lambda l: l.startswith("BA.2"),
    "BA.4":     lambda l: l.startswith("BA.4"),
    "BA.5":     lambda l: l.startswith("BA.5") or l.startswith("BE."),
    "BQ.1":     lambda l: l.startswith("BQ."),
    "XBB":      lambda l: l.startswith("XBB"),
}

WAVE_LABELS = {
    "B.1.177": "B.1.177 (pre-Alpha)",
    "Alpha":   "Alpha (B.1.1.7)",
    "Delta":   "Delta (AY.*/B.1.617.2)",
    "BA.1":    "Omicron BA.1",
    "BA.2":    "Omicron BA.2",
    "BA.4":    "Omicron BA.4",
    "BA.5":    "Omicron BA.5 / BE.*",
    "BQ.1":    "Omicron BQ.1",
    "XBB":     "Omicron XBB",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def assign_wave(lineage: str) -> str:
    if not isinstance(lineage, str):
        return "Other"
    for name, predicate in WAVE_GROUPS.items():
        if predicate(lineage):
            return name
    return "Other"


def compute_wave_table(parquet_path: str, threshold: float = 0.05) -> pd.DataFrame:
    """
    Read the parquet file, compute weekly counts per wave group,
    and return a summary DataFrame with start / peak / end dates.
    """
    # Read only the columns we need
    df = pd.read_parquet(parquet_path, columns=["collection_date", "pango_lineage"])
    df.dropna(subset=["collection_date"], inplace=True)

    df["collection_date"] = pd.to_datetime(df["collection_date"])
    # Week anchor = Monday of the ISO week
    df["week"] = df["collection_date"].dt.to_period("W").apply(lambda p: p.start_time)
    df["wave_group"] = df["pango_lineage"].apply(assign_wave)

    # Weekly counts per wave group
    weekly = (
        df.groupby(["week", "wave_group"])
        .size()
        .unstack(fill_value=0)
    )
    weekly["total"] = weekly.sum(axis=1)

    rows = []
    for wave_key in WAVE_GROUPS:
        if wave_key not in weekly.columns:
            continue

        counts = weekly[wave_key]
        pct    = counts / weekly["total"]

        dominant = pct[pct > threshold]
        if dominant.empty:
            continue

        start    = dominant.index[0]
        end      = dominant.index[-1]
        peak_idx = counts.idxmax()
        peak_val = int(counts.max())

        rows.append({
            "Wave":          WAVE_LABELS[wave_key],
            "Start (>5%)":   start.date(),
            "Peak week":     peak_idx.date(),
            "Peak seqs":     peak_val,
            "End (<5%)":     end.date(),
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute SARS-CoV-2 wave dates.")
    parser.add_argument(
        "--input",
        default="../data/processed/scotland_clustering_analysis_dataset.parquet",
        help="Path to the input parquet file.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="Fraction of weekly sequences to define wave start/end (default: 0.05).",
    )
    args = parser.parse_args()

    print(f"Reading: {args.input}")
    print(f"Dominance threshold: {args.threshold:.0%}\n")

    table = compute_wave_table(args.input, threshold=args.threshold)

    # Pretty-print
    col_widths = {col: max(len(col), table[col].astype(str).str.len().max())
                  for col in table.columns}
    header = "  ".join(col.ljust(col_widths[col]) for col in table.columns)
    sep    = "  ".join("-" * col_widths[col] for col in table.columns)
    print(header)
    print(sep)
    for _, row in table.iterrows():
        print("  ".join(str(row[col]).ljust(col_widths[col]) for col in table.columns))

    # Also save to CSV alongside the script
    table.to_parquet("tables/wave_dates.parquet", index=False)
    print("\nSaved to: tables/wave_dates.parquet")


if __name__ == "__main__":
    main()