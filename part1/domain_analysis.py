"""Part 1 supplements — SIMD subdomain and wave-stratified mixing analyses.

This script extends ``overall_analysis.py``. It uses the same primary Leiden
resolution, QC filter, lineage pooling, calendar spline, and window-clustered
standard errors as the Part 1 analysis.

Outputs include:

* SIMD-domain hurdle/ZTNB count models (Line 1 supplement).
* SIMD-domain hurdle/ZTNB count models with mixing predictors
  (Line 2 supplement).
* SIMD-domain quintile mixing models.
* SIMD-domain demographic (age, sex, age × sex) mixing models.
* Wave-specific SIMD-domain demographic mixing models.
* Primary-resolution observed-vs-expected pair-probability matrices.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

# Allow ``python part1/domain_analysis.py`` to import lib/.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.constants import (  # noqa: E402
    CALENDAR_SPLINE_DF,
    DEMOGRAPHIC_MIXING,
    DEMOGRAPHIC_MIXING_PREDICTOR_TERMS,
    DOMAINS,
    LINEAGE_MIN_CLUSTERS,
    MATRIX_VARIABLES,
    PRIMARY_RESOLUTION,
    QC_DEFAULT,
    SHARED_COUNT_TERMS,
    SHARED_MIXING_TERMS,
    WAVE_LABELS,
    WAVE_ORDER,
)
from lib.data_prep import (  # noqa: E402
    assign_wave,
    build_domain_wave_cluster_table,
    build_matrix_for_variable,
    expected_stratum_discordance,
    observed_cluster_discordance,
    read_domain_wave_sequence_rows,
    repo_root,
    summarise_domain_wave_dataset,
)
from lib.estimators import lineage_levels  # noqa: E402
from lib.fit_models import (  # noqa: E402
    fit_domain_count_models,
    fit_domain_demographic_mixing_models,
    fit_domain_quintile_mixing_models,
    fit_wave_domain_demographic_mixing_models,
)


__all__ = [
    # constants surfaced for downstream scripts
    "DEMOGRAPHIC_MIXING", "DEMOGRAPHIC_MIXING_PREDICTOR_TERMS", "DOMAINS",
    "MATRIX_VARIABLES", "SHARED_COUNT_TERMS", "SHARED_MIXING_TERMS",
    "WAVE_LABELS", "WAVE_ORDER",
    # data prep helpers
    "assign_wave", "build_domain_wave_cluster_table", "build_matrix_for_variable",
    "expected_stratum_discordance", "observed_cluster_discordance",
    "read_domain_wave_sequence_rows",
    # fit pipelines
    "fit_domain_count_models", "fit_domain_demographic_mixing_models",
    "fit_domain_quintile_mixing_models", "fit_wave_domain_demographic_mixing_models",
    # entry points
    "main", "parse_args", "run",
]


# ---------------------------------------------------------------------------
# Top-level run
# ---------------------------------------------------------------------------


def run(
    root: Path,
    qc: str | None,
    primary_resolution: float,
    lineage_min_clusters: int,
    calendar_spline_df: int,
    maxiter: int,
    min_wave_clusters: int,
    min_wave_windows: int,
) -> None:
    out_dir = root / "part1"
    tables_dir = out_dir / "tables"
    cache_dir = out_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    print("Reading main-formulation domain/wave sequence rows", flush=True)
    seq = read_domain_wave_sequence_rows(qc=qc, primary_resolution=primary_resolution)
    print(f"Building domain/wave cluster table from {len(seq):,} sequence rows", flush=True)
    clusters, scaling, dropped = build_domain_wave_cluster_table(
        seq,
        lineage_min_clusters=lineage_min_clusters,
        calendar_spline_df=calendar_spline_df,
    )
    calendar_cols = [col for col in clusters.columns if col.startswith("calendar_spline_")]
    lineage_levels_all = lineage_levels(clusters)

    clusters.to_parquet(cache_dir / "domain_wave_cluster_table.parquet", index=False)
    scaling.to_csv(tables_dir / "domain_wave_covariate_scaling.csv", index=False)
    summarise_domain_wave_dataset(seq, clusters, qc, primary_resolution, dropped).to_csv(
        tables_dir / "domain_wave_dataset_descriptives.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # Line 1 supplement — domain deprivation effects on count outcomes
    # ------------------------------------------------------------------
    print("Fitting SIMD-domain hurdle/ZTNB count models", flush=True)
    count_results, count_diagnostics = fit_domain_count_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
        maxiter=maxiter,
    )
    count_results.to_csv(
        tables_dir / "simd_domain_hurdle_count_model_results.csv",
        index=False,
    )
    count_diagnostics.to_csv(
        tables_dir / "simd_domain_hurdle_count_model_diagnostics.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # Line 2 supplement — domain deprivation + mixing predictors
    # ------------------------------------------------------------------
    print("Fitting SIMD-domain hurdle/ZTNB count models with mixing predictors", flush=True)
    mixing_predictor_count_results, mixing_predictor_count_diagnostics = (
        fit_domain_count_models(
            clusters,
            lineage_levels_all=lineage_levels_all,
            calendar_cols=calendar_cols,
            maxiter=maxiter,
            include_mixing_predictors=True,
        )
    )
    mixing_predictor_count_results.to_csv(
        tables_dir / "simd_domain_mixing_predictor_hurdle_count_model_results.csv",
        index=False,
    )
    mixing_predictor_count_diagnostics.to_csv(
        tables_dir / "simd_domain_mixing_predictor_hurdle_count_model_diagnostics.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # Line 1 supplement — domain deprivation effects on mixing outcomes
    # ------------------------------------------------------------------
    print("Fitting SIMD-domain quintile mixing models", flush=True)
    domain_mixing, domain_mixing_diag = fit_domain_quintile_mixing_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
    )
    domain_mixing.to_csv(
        tables_dir / "simd_domain_quintile_mixing_model_results.csv",
        index=False,
    )
    domain_mixing_diag.to_csv(
        tables_dir / "simd_domain_quintile_mixing_model_diagnostics.csv",
        index=False,
    )

    print("Fitting SIMD-domain demographic mixing models", flush=True)
    domain_demo, domain_demo_diag = fit_domain_demographic_mixing_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
    )
    domain_demo.to_csv(
        tables_dir / "simd_domain_demographic_mixing_model_results.csv",
        index=False,
    )
    domain_demo_diag.to_csv(
        tables_dir / "simd_domain_demographic_mixing_model_diagnostics.csv",
        index=False,
    )

    print("Fitting wave-specific SIMD-domain demographic mixing models", flush=True)
    wave_demo, wave_demo_diag = fit_wave_domain_demographic_mixing_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
        min_clusters=min_wave_clusters,
        min_windows=min_wave_windows,
    )
    wave_demo.to_csv(
        tables_dir / "wave_specific_domain_demographic_mixing_model_results.csv",
        index=False,
    )
    wave_demo_diag.to_csv(
        tables_dir / "wave_specific_domain_demographic_mixing_model_diagnostics.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # Observed-vs-expected pair-probability matrices
    # ------------------------------------------------------------------
    print("Building primary-resolution observed-vs-expected mixing matrices", flush=True)
    matrices = pd.concat(
        [build_matrix_for_variable(seq, variable) for variable in MATRIX_VARIABLES],
        ignore_index=True,
    )
    matrices.to_csv(tables_dir / "observed_expected_mixing_matrices.csv", index=False)

    print(f"Wrote main domain/wave tables to {tables_dir}", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument("--qc", default=QC_DEFAULT)
    parser.add_argument("--primary-resolution", type=float, default=PRIMARY_RESOLUTION)
    parser.add_argument("--lineage-min-clusters", type=int, default=LINEAGE_MIN_CLUSTERS)
    parser.add_argument("--calendar-spline-df", type=int, default=CALENDAR_SPLINE_DF)
    parser.add_argument("--maxiter", type=int, default=1000)
    parser.add_argument("--min-wave-clusters", type=int, default=1000)
    parser.add_argument("--min-wave-windows", type=int, default=4)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    qc = None if str(args.qc).lower() == "none" else str(args.qc)
    run(
        root=args.root.resolve(),
        qc=qc,
        primary_resolution=args.primary_resolution,
        lineage_min_clusters=args.lineage_min_clusters,
        calendar_spline_df=args.calendar_spline_df,
        maxiter=args.maxiter,
        min_wave_clusters=args.min_wave_clusters,
        min_wave_windows=args.min_wave_windows,
    )


if __name__ == "__main__":
    main()
