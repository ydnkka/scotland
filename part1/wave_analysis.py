"""Wave-stratified formulation cluster outcome models.

Fits the hurdle/ZTNB count models separately within epidemic wave
groups for cluster size and geographic spread.  The figure generated from
these tables focuses on the wave-specific SIMD-deprivation coefficient; the
tables retain all covariates. A companion sensitivity adds the
excess-mixing metrics as outcome predictors (Line 2 supplement).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

# Allow ``python part1/wave_analysis.py`` to import lib/.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.constants import MIXING_PREDICTOR_TERMS  # noqa: E402
from lib.data_prep import (  # noqa: E402
    assign_wave,
    ensure_mixing_predictor_columns,
    load_cluster_table,
    repo_root,
)
from lib.fit_models import (  # noqa: E402
    fit_wave_outcome_models,
    summarise_wave_outcomes,
)


# ---------------------------------------------------------------------------
# Top-level run
# ---------------------------------------------------------------------------


def run(
    root: Path,
    *,
    maxiter: int,
    min_clusters: int,
    min_windows: int,
    min_positive: int,
    min_events: int,
    tables_dir: Path | None = None,
    cache_dir: Path | None = None,
) -> None:
    main_dir = root / "part1"
    tables_dir = tables_dir or main_dir / "tables"
    cache_dir = cache_dir or main_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)

    clusters = ensure_mixing_predictor_columns(
        load_cluster_table(root=root, cache_dir=cache_dir)
    )
    clusters["wave_group"] = clusters["pango_lineage"].astype(str).map(assign_wave)

    print("Fitting wave-specific cluster outcome models", flush=True)
    results, diagnostics = fit_wave_outcome_models(
        clusters,
        maxiter=maxiter,
        min_clusters=min_clusters,
        min_windows=min_windows,
        min_positive=min_positive,
        min_events=min_events,
    )
    results.to_csv(tables_dir / "wave_specific_hurdle_count_model_results.csv", index=False)
    diagnostics.to_csv(
        tables_dir / "wave_specific_hurdle_count_model_diagnostics.csv",
        index=False,
    )

    print("Fitting wave-specific cluster outcome models with mixing predictors", flush=True)
    mixing_results, mixing_diagnostics = fit_wave_outcome_models(
        clusters,
        maxiter=maxiter,
        min_clusters=min_clusters,
        min_windows=min_windows,
        min_positive=min_positive,
        min_events=min_events,
        extra_terms=MIXING_PREDICTOR_TERMS,
        predictor_set="primary_plus_mixing",
        skip_cluster_size_binary=True,
    )
    mixing_results.to_csv(
        tables_dir / "wave_specific_mixing_predictor_hurdle_count_model_results.csv",
        index=False,
    )
    mixing_diagnostics.to_csv(
        tables_dir / "wave_specific_mixing_predictor_hurdle_count_model_diagnostics.csv",
        index=False,
    )
    summarise_wave_outcomes(clusters).to_csv(
        tables_dir / "wave_cluster_outcome_descriptives.csv",
        index=False,
    )
    print(
        f"Wrote {tables_dir / 'wave_specific_hurdle_count_model_results.csv'}",
        flush=True,
    )
    print(
        f"Wrote {tables_dir / 'wave_specific_mixing_predictor_hurdle_count_model_results.csv'}",
        flush=True,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument("--maxiter", type=int, default=1000)
    parser.add_argument("--min-clusters", type=int, default=1000)
    parser.add_argument("--min-windows", type=int, default=4)
    parser.add_argument("--min-positive", type=int, default=500)
    parser.add_argument("--min-events", type=int, default=50)
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=None,
        help=(
            "Directory for output CSV tables. Defaults to part1/tables. "
            "Set to match the --tables-dir used by overall_analysis.py for the "
            "same sensitivity run."
        ),
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing the cluster_table.parquet cache. "
            "Defaults to part1/cache. Must match the --cache-dir used "
            "by overall_analysis.py for the same sensitivity run."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    run(
        root=args.root.resolve(),
        maxiter=args.maxiter,
        min_clusters=args.min_clusters,
        min_windows=args.min_windows,
        min_positive=args.min_positive,
        min_events=args.min_events,
        tables_dir=args.tables_dir.resolve() if args.tables_dir else None,
        cache_dir=args.cache_dir.resolve() if args.cache_dir else None,
    )


if __name__ == "__main__":
    main()
