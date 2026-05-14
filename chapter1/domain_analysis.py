"""Chapter 1 — SIMD-domain stratification.

Refits the main-effects model with each SIMD domain's quintile mixing in
place of overall SIMD mixing.  Age and sex excess mixing remain in every
model so the SIMD-channel estimate is directly comparable across domains.

Output tables and figures end up in ``chapter1/tables/`` and
``chapter1/figures/`` alongside the overall-analysis outputs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.constants import (  # noqa: E402
    CALENDAR_SPLINE_DF,
    LINEAGE_MIN_CLUSTERS,
    PRIMARY_RESOLUTION,
    MATRIX_VARIABLES,
    QC_DEFAULT,
)
from lib.data_prep import (  # noqa: E402
    build_domain_cluster_table,
    read_sequence_rows,
    repo_root,
    build_matrix_for_variable
)
from lib.fit_models import fit_domain_main_effects  # noqa: E402
from lib.plots import plot_stratified_forest  # noqa: E402


def run(
    root: Path,
    *,
    qc: str | None = QC_DEFAULT,
    primary_resolution: float = PRIMARY_RESOLUTION,
    lineage_min_clusters: int = LINEAGE_MIN_CLUSTERS,
    calendar_spline_df: int = CALENDAR_SPLINE_DF,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    tables_dir: Path | None = None,
    figures_dir: Path | None = None,
    cache_dir: Path | None = None,
) -> None:
    out_dir = root / "chapter1"
    tables_dir = tables_dir or out_dir / "tables"
    figures_dir = figures_dir or out_dir / "figures"
    cache_dir = cache_dir or out_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    print("[chapter1.domain] loading sequence rows (with domain ranks)", flush=True)
    seq = read_sequence_rows(
        qc=qc, primary_resolution=primary_resolution, include_domains=True,
    )
    print(
        f"[chapter1.domain] building cluster table from {len(seq):,} rows",
        flush=True,
    )
    clusters, scaling, _ = build_domain_cluster_table(
        seq,
        lineage_min_clusters=lineage_min_clusters,
        calendar_spline_df=calendar_spline_df,
    )
    clusters.to_parquet(cache_dir / "domain_cluster_table.parquet", index=False)
    scaling.to_csv(tables_dir / "domain_covariate_scaling.csv", index=False)

    print("[chapter1.domain] fitting domain-stratified main effects", flush=True)
    domain_results, domain_diag = fit_domain_main_effects(
        clusters, cluster_by=cluster_by, maxiter=maxiter,
    )
    domain_results.to_csv(
        tables_dir / "domain_main_effects_results.csv", index=False,
    )
    domain_diag.to_csv(
        tables_dir / "domain_main_effects_diagnostics.csv", index=False,
    )

    print("[chapter1.domain] Building primary-resolution observed-vs-expected mixing matrices", flush=True)
    matrices = pd.concat(
        [build_matrix_for_variable(seq, variable) for variable in MATRIX_VARIABLES],
        ignore_index=True,
    )
    matrices.to_csv(tables_dir / "observed_expected_mixing_matrices.csv", index=False)

    plot_stratified_forest(
        domain_results,
        figures_dir / "domain_forest",
        stratum_col="domain",
        terms=[
            "age_excess_mixing_z",
            "sex_excess_mixing_z",
            "overall_domain_excess_mixing_z",
            "income_domain_excess_mixing_z",
            "employment_domain_excess_mixing_z",
            "education_domain_excess_mixing_z",
            "health_domain_excess_mixing_z",
            "access_domain_excess_mixing_z",
            "crime_domain_excess_mixing_z",
            "housing_domain_excess_mixing_z",
        ],
        title="Excess mixing → cluster scale, by SIMD domain",
    )
    print(f"[chapter1.domain] tables written under {tables_dir}", flush=True)
    print(f"[chapter1.domain] figures written under {figures_dir}", flush=True)


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument("--qc", default=QC_DEFAULT)
    parser.add_argument("--primary-resolution", type=float, default=PRIMARY_RESOLUTION)
    parser.add_argument("--lineage-min-clusters", type=int, default=LINEAGE_MIN_CLUSTERS)
    parser.add_argument("--calendar-spline-df", type=int, default=CALENDAR_SPLINE_DF)
    parser.add_argument("--cluster-by", default="window_id",
                        choices=["window_id", "health_board"])
    parser.add_argument("--maxiter", type=int, default=1000)
    parser.add_argument("--tables-dir", type=Path, default=None)
    parser.add_argument("--figures-dir", type=Path, default=None)
    parser.add_argument("--cache-dir", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    qc = None if str(args.qc).lower() == "none" else str(args.qc)
    run(
        root=args.root.resolve(),
        qc=qc,
        primary_resolution=args.primary_resolution,
        lineage_min_clusters=args.lineage_min_clusters,
        calendar_spline_df=args.calendar_spline_df,
        cluster_by=args.cluster_by,
        maxiter=args.maxiter,
        tables_dir=args.tables_dir.resolve() if args.tables_dir else None,
        figures_dir=args.figures_dir.resolve() if args.figures_dir else None,
        cache_dir=args.cache_dir.resolve() if args.cache_dir else None,
    )


if __name__ == "__main__":
    main()
