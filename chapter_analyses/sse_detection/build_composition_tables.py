"""Build cluster-level socio-demographic composition tables for SSE analysis.

Each output contains one row per cluster meeting the detector's minimum size.
Level columns contain within-cluster proportions among non-missing observations,
followed by the detector status and burst/burden scores.

Run from the repository root::

    python -m chapter_analyses.sse_detection.build_composition_tables
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from .lib.model.prep import COMPOSITION_SPECS
from .lib.sse.config import MIN_CLUSTER_SIZE, RESULTS_DIR, SSE_OUTPUT_DIR
from .lib.sse.detection import load_sequence_data
from .lib.sse.io import HIGH_PRIORITY_CANDIDATE_TIERS


LOGGER = logging.getLogger(__name__)
DEFAULT_OUTPUT_DIR = RESULTS_DIR / "tables"
DEFAULT_VARIABLES = tuple(str(spec["column"]) for spec in COMPOSITION_SPECS)
DETECTION_COLUMNS = ("candidate_tier", "burst_score", "burden_score")


def _level_name(value: object) -> str:
    """Return a stable, parquet-safe column name for a variable level."""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _validate_inputs(
    sequence_data: pd.DataFrame,
    cluster_table: pd.DataFrame,
    variables: Sequence[str],
) -> None:
    sequence_required = {"cluster_id", *variables}
    cluster_required = {"cluster_id", "cluster_size", *DETECTION_COLUMNS}
    missing_sequence = sorted(sequence_required.difference(sequence_data.columns))
    missing_cluster = sorted(cluster_required.difference(cluster_table.columns))
    if missing_sequence:
        raise KeyError(f"Sequence data is missing columns: {missing_sequence}")
    if missing_cluster:
        raise KeyError(f"Cluster table is missing columns: {missing_cluster}")
    if cluster_table["cluster_id"].duplicated().any():
        raise ValueError("cluster_table must contain one row per cluster_id")


def build_variable_table(
    sequence_data: pd.DataFrame,
    cluster_table: pd.DataFrame,
    variable: str,
    *,
    min_cluster_size: int = MIN_CLUSTER_SIZE,
) -> pd.DataFrame:
    """Build one wide cluster-composition table for ``variable``.

    Missing values do not contribute to the denominator. Clusters with no
    observed value for the variable are retained, with all level proportions
    set to missing.
    """
    eligible = cluster_table.loc[
        cluster_table["cluster_size"].ge(min_cluster_size),
        ["cluster_id", *DETECTION_COLUMNS],
    ].copy()
    eligible["sse_status"] = eligible["candidate_tier"].isin(
        HIGH_PRIORITY_CANDIDATE_TIERS
    ).map({True: "candidate", False: "background"})

    observed = sequence_data.loc[
        sequence_data["cluster_id"].isin(eligible["cluster_id"]),
        ["cluster_id", variable],
    ].dropna(subset=[variable])

    if observed.empty:
        proportions = pd.DataFrame(index=pd.Index([], name="cluster_id"))
    else:
        counts = pd.crosstab(observed["cluster_id"], observed[variable])
        proportions = counts.div(counts.sum(axis=1), axis=0)
        proportions = proportions.reindex(
            sorted(proportions.columns, key=lambda value: str(value)), axis=1
        )
        level_names = [_level_name(value) for value in proportions.columns]
        if len(level_names) != len(set(level_names)):
            raise ValueError(
                f"String conversion produces duplicate levels for {variable}: "
                f"{level_names}"
            )
        proportions.columns = level_names

    out = eligible[["cluster_id"]].merge(
        proportions.reset_index(), on="cluster_id", how="left"
    )
    metadata = eligible[
        ["cluster_id", "sse_status", "burst_score", "burden_score"]
    ]
    out = out.merge(metadata, on="cluster_id", how="left", validate="one_to_one")
    return out


def build_composition_tables(
    sequence_data: pd.DataFrame,
    cluster_table: pd.DataFrame,
    *,
    variables: Iterable[str] = DEFAULT_VARIABLES,
    min_cluster_size: int = MIN_CLUSTER_SIZE,
) -> dict[str, pd.DataFrame]:
    """Build composition tables keyed by source variable name."""
    variables = tuple(dict.fromkeys(variables))
    if not variables:
        raise ValueError("At least one composition variable is required")
    _validate_inputs(sequence_data, cluster_table, variables)
    return {
        variable: build_variable_table(
            sequence_data,
            cluster_table,
            variable,
            min_cluster_size=min_cluster_size,
        )
        for variable in variables
    }


def write_composition_tables(
    tables: dict[str, pd.DataFrame], output_dir: Path
) -> dict[str, tuple[Path, Path]]:
    """Write each table as parquet and CSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, tuple[Path, Path]] = {}
    for variable, table in tables.items():
        stem = f"cluster_composition_{variable}"
        parquet_path = output_dir / f"{stem}.parquet"
        csv_path = output_dir / f"{stem}.csv"
        table.to_parquet(parquet_path, index=False)
        table.to_csv(csv_path, index=False)
        outputs[variable] = (parquet_path, csv_path)
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cluster-table",
        type=Path,
        default=SSE_OUTPUT_DIR / "cluster_table.parquet",
        help="Detector cluster table.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for the composition tables.",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_VARIABLES),
        help="Sequence-level categorical variables to tabulate.",
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=MIN_CLUSTER_SIZE,
        help="Minimum cluster size retained in every output.",
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
    if args.min_cluster_size < 1:
        raise ValueError("--min-cluster-size must be at least 1")
    if not args.cluster_table.exists():
        raise FileNotFoundError(f"Missing cluster table: {args.cluster_table}")

    LOGGER.info("Loading detector cluster table")
    cluster_table = pd.read_parquet(args.cluster_table)
    LOGGER.info("Loading sequence-level analysis data")
    sequence_data = load_sequence_data()
    tables = build_composition_tables(
        sequence_data,
        cluster_table,
        variables=args.variables,
        min_cluster_size=args.min_cluster_size,
    )
    outputs = write_composition_tables(tables, args.output_dir)
    for variable, paths in outputs.items():
        LOGGER.info(
            "Wrote %s (%d clusters, %d levels): %s and %s",
            variable,
            len(tables[variable]),
            len(tables[variable].columns) - 4,
            paths[0],
            paths[1],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
