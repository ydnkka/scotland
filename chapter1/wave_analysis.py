"""Chapter 1 — wave stratification.

Refits the main-effects model separately within each epidemic wave.
Waves with fewer than ``--min-clusters-per-wave`` non-singleton clusters
are skipped to avoid unstable fits, with the omission noted
in the diagnostics table.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.constants import (  # noqa: E402
    CALENDAR_SPLINE_DF,
    LINEAGE_MIN_CLUSTERS,
    PRIMARY_RESOLUTION,
    QC_DEFAULT,
)
from lib.data_prep import (  # noqa: E402
    build_cluster_table,
    read_sequence_rows,
    repo_root,
)
from lib.fit_models import fit_wave_stratified  # noqa: E402
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
    min_clusters_per_wave: int = 50,
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

    print("[chapter1.wave] loading sequence rows", flush=True)
    seq = read_sequence_rows(qc=qc, primary_resolution=primary_resolution)
    print(
        f"[chapter1.wave] building cluster table from {len(seq):,} rows",
        flush=True,
    )
    clusters, scaling, _ = build_cluster_table(
        seq,
        lineage_min_clusters=lineage_min_clusters,
        calendar_spline_df=calendar_spline_df,
    )
    clusters.to_parquet(cache_dir / "wave_cluster_table.parquet", index=False)

    print("[chapter1.wave] fitting wave-stratified main effects", flush=True)
    wave_results, wave_diag = fit_wave_stratified(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        min_clusters_per_wave=min_clusters_per_wave,
    )
    wave_results.to_csv(
        tables_dir / "wave_stratified_results.csv", index=False,
    )
    wave_diag.to_csv(
        tables_dir / "wave_stratified_diagnostics.csv", index=False,
    )

    plot_stratified_forest(
        wave_results,
        figures_dir / "wave_stratified_forest",
        stratum_col="wave",
        title="Excess mixing → cluster scale, by epidemic wave",
    )
    print(f"[chapter1.wave] tables written under {tables_dir}", flush=True)
    print(f"[chapter1.wave] figures written under {figures_dir}", flush=True)


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
    parser.add_argument("--min-clusters-per-wave", type=int, default=200)
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
        min_clusters_per_wave=args.min_clusters_per_wave,
        tables_dir=args.tables_dir.resolve() if args.tables_dir else None,
        figures_dir=args.figures_dir.resolve() if args.figures_dir else None,
        cache_dir=args.cache_dir.resolve() if args.cache_dir else None,
    )


if __name__ == "__main__":
    main()
