"""Chapter 1 — overall analysis.

Runs the core fits on the full primary-resolution cluster table:

1. **Main effects** — cluster size (ZTNB) and geographic spread (ZTNB)
   regressed on excess age, sex, and SIMD mixing plus adjustments, with
   lineage and an 8-df calendar B-spline as nuisance.  Both outcomes are
   ZTNB-only on the non-singleton sub-population.
2. **Wave interactions** — same outcomes; mixing × wave interactions added,
   lineage adjustment replaced by wave dummies.
3. **Size-adjusted spread (linear log size)** — adds log(cluster_size) as
   covariate to the spread ZTNB.
4. **Sensitivity: spline on log(cluster_size)** in the spread model.
5. **Sensitivity: SIMD-decile mixing** — swap the SIMD-quintile excess
   mixing predictor for an SIMD-decile version; age and sex remain.
6. **Sensitivity: finite-sample mixing** — replace raw observed-minus-
   expected excess mixing with finite-sample standardised excess mixing.
7. **Sensitivity: joint-profile adjusted predictor set** — add the joint
   age × sex × SIMD profile term to the three main mixing predictors.
8. **Sensitivity: null-residual mixing** — refit main effects using
   residuals from a per-dimension null regression as the predictor.
9. **Supplementary: joint profile predictors** — one fit per joint
   profile: demographic (age × sex) and sociodemographic (age × sex ×
   SIMD).

Run from the repo root::

    conda run -n PhD python chapter1/overall_analysis.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Make ``import lib.…`` work when invoked as a script.
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
    summarise_dataset,
)
from lib.constants import PROFILE_PREDICTORS  # noqa: E402
from lib.fit_models import (  # noqa: E402
    build_null_residual_mixing,
    fit_finite_sample_mixing_sensitivity,
    fit_joint_profile_adjusted_sensitivity,
    fit_main_effects,
    fit_null_residual_sensitivity,
    fit_profile_predictor,
    fit_simd_decile_sensitivity,
    fit_size_spline_sensitivity,
    fit_wave_interactions,
)
from lib.plots import (  # noqa: E402
    plot_main_effects_forest,
    plot_wave_interaction_slopes,
)


def run(
    root: Path,
    *,
    qc: str | None = QC_DEFAULT,
    primary_resolution: float = PRIMARY_RESOLUTION,
    lineage_min_clusters: int = LINEAGE_MIN_CLUSTERS,
    calendar_spline_df: int = CALENDAR_SPLINE_DF,
    cluster_by: str = "window_id",
    maxiter: int = 1000,
    winsorise_quantile: float = 0.0,
    exclude_tail_quantile: float = 0.0,
    window_stride: int = 1,
    tables_dir: Path | None = None,
    figures_dir: Path | None = None,
    cache_dir: Path | None = None,
    sample_clusters: int | None = None,
) -> None:
    if winsorise_quantile > 0.0 and exclude_tail_quantile > 0.0:
        raise ValueError(
            "Use either --winsorise-quantile or --exclude-tail-quantile, not both."
        )
    for name, value in {
        "winsorise_quantile": winsorise_quantile,
        "exclude_tail_quantile": exclude_tail_quantile,
    }.items():
        if value < 0.0 or value >= 1.0:
            raise ValueError(f"{name} must be in [0, 1); got {value}.")
    if window_stride < 1:
        raise ValueError(f"window_stride must be >= 1; got {window_stride}.")

    out_dir = root / "chapter1"
    tables_dir = tables_dir or out_dir / "tables"
    figures_dir = figures_dir or out_dir / "figures"
    cache_dir = cache_dir or out_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"[chapter1] loading sequence rows (qc={qc}, "
        f"resolution={primary_resolution})", flush=True,
    )
    seq = read_sequence_rows(qc=qc, primary_resolution=primary_resolution)
    print(f"[chapter1] building cluster table from {len(seq):,} rows", flush=True)
    clusters, scaling, dropped = build_cluster_table(
        seq,
        lineage_min_clusters=lineage_min_clusters,
        calendar_spline_df=calendar_spline_df,
    )
    if window_stride > 1:
        keep_idx = clusters["window_idx"] % window_stride == 0
        n_before = len(clusters)
        clusters = clusters.loc[keep_idx].copy()
        print(
            f"[chapter1] non-overlapping window filter "
            f"(stride={window_stride}): {n_before:,} -> {len(clusters):,} "
            "clusters",
            flush=True,
        )
    if sample_clusters is not None and sample_clusters < len(clusters):
        clusters = clusters.sample(sample_clusters, random_state=42).reset_index(drop=True)
        print(
            f"[chapter1] sampled cluster table to {len(clusters):,} clusters "
            "(quick-test mode)", flush=True,
        )

    clusters.to_parquet(cache_dir / "cluster_table.parquet", index=False)
    scaling.to_csv(tables_dir / "covariate_scaling.csv", index=False)
    summarise_dataset(
        seq, clusters, qc, primary_resolution, dropped,
    ).to_csv(tables_dir / "dataset_descriptives.csv", index=False)

    # ------------------------------------------------------------------
    # 1. Main effects
    # ------------------------------------------------------------------
    print("[chapter1] fitting main-effects model", flush=True)
    main_results, main_diag = fit_main_effects(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )
    main_results.to_csv(tables_dir / "main_effects_results.csv", index=False)
    main_diag.to_csv(tables_dir / "main_effects_diagnostics.csv", index=False)

    # 1b. Size-adjusted geographic spread
    print("[chapter1] fitting size-adjusted geographic-spread model", flush=True)
    size_adj_results, size_adj_diag = fit_main_effects(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        include_log_size=True,
        model_label="main",
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )
    size_adj_results.to_csv(
        tables_dir / "size_adjusted_spread_results.csv", index=False,
    )
    size_adj_diag.to_csv(
        tables_dir / "size_adjusted_spread_diagnostics.csv", index=False,
    )

    # ------------------------------------------------------------------
    # 2. Wave interactions
    # ------------------------------------------------------------------
    print("[chapter1] fitting wave-interaction model", flush=True)
    wave_results, wave_diag = fit_wave_interactions(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )
    wave_results.to_csv(
        tables_dir / "wave_interaction_results.csv", index=False,
    )
    wave_diag.to_csv(
        tables_dir / "wave_interaction_diagnostics.csv", index=False,
    )

    # ------------------------------------------------------------------
    # 3. Size-spline sensitivity
    # ------------------------------------------------------------------
    print("[chapter1] fitting size-spline sensitivity", flush=True)
    spline_results, spline_diag = fit_size_spline_sensitivity(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )
    spline_results.to_csv(
        tables_dir / "size_spline_sensitivity_results.csv", index=False,
    )
    spline_diag.to_csv(
        tables_dir / "size_spline_sensitivity_diagnostics.csv", index=False,
    )

    # ------------------------------------------------------------------
    # 3b. SIMD-decile sensitivity (swap quintile mixing for decile)
    # ------------------------------------------------------------------
    print("[chapter1] fitting SIMD-decile sensitivity", flush=True)
    decile_results, decile_diag = fit_simd_decile_sensitivity(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )
    decile_results.to_csv(
        tables_dir / "simd_decile_sensitivity_results.csv", index=False,
    )
    decile_diag.to_csv(
        tables_dir / "simd_decile_sensitivity_diagnostics.csv", index=False,
    )

    # ------------------------------------------------------------------
    # 3c. Finite-sample standardised mixing sensitivity
    # ------------------------------------------------------------------
    print("[chapter1] fitting finite-sample mixing sensitivity", flush=True)
    finite_results, finite_diag = fit_finite_sample_mixing_sensitivity(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )
    finite_results.to_csv(
        tables_dir / "finite_sample_mixing_sensitivity_results.csv",
        index=False,
    )
    finite_diag.to_csv(
        tables_dir / "finite_sample_mixing_sensitivity_diagnostics.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # 3d. Predictor-set sensitivity (main predictors + joint profile)
    # ------------------------------------------------------------------
    print("[chapter1] fitting joint-profile adjusted sensitivity", flush=True)
    joint_results, joint_diag = fit_joint_profile_adjusted_sensitivity(
        clusters,
        cluster_by=cluster_by,
        maxiter=maxiter,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )
    joint_results.to_csv(
        tables_dir / "joint_profile_adjusted_results.csv", index=False,
    )
    joint_diag.to_csv(
        tables_dir / "joint_profile_adjusted_diagnostics.csv", index=False,
    )

    # ------------------------------------------------------------------
    # 4. Null-residual mixing sensitivity
    # ------------------------------------------------------------------
    print("[chapter1] building null-residual mixing", flush=True)
    clusters_with_null = build_null_residual_mixing(clusters)
    print("[chapter1] fitting null-residual sensitivity", flush=True)
    null_results, null_diag = fit_null_residual_sensitivity(
        clusters_with_null,
        cluster_by=cluster_by,
        maxiter=maxiter,
        winsorise_quantile=winsorise_quantile,
        exclude_tail_quantile=exclude_tail_quantile,
    )
    null_results.to_csv(
        tables_dir / "null_residual_sensitivity_results.csv", index=False,
    )
    null_diag.to_csv(
        tables_dir / "null_residual_sensitivity_diagnostics.csv", index=False,
    )

    # ------------------------------------------------------------------
    # 5. Joint-profile predictors (demographic + sociodemographic)
    # ------------------------------------------------------------------
    profile_results_frames = []
    profile_diag_frames = []
    for profile_name in PROFILE_PREDICTORS:
        print(
            f"[chapter1] fitting joint-profile predictor: {profile_name}",
            flush=True,
        )
        profile_results, profile_diag = fit_profile_predictor(
            clusters,
            profile_name=profile_name,
            cluster_by=cluster_by,
            maxiter=maxiter,
            winsorise_quantile=winsorise_quantile,
            exclude_tail_quantile=exclude_tail_quantile,
        )
        if not profile_results.empty:
            profile_results["profile"] = profile_name
        profile_diag["profile"] = profile_name
        profile_results_frames.append(profile_results)
        profile_diag_frames.append(profile_diag)

    profile_results_all = (
        pd.concat(profile_results_frames, ignore_index=True)
        if profile_results_frames else pd.DataFrame()
    )
    profile_diag_all = (
        pd.concat(profile_diag_frames, ignore_index=True)
        if profile_diag_frames else pd.DataFrame()
    )
    profile_results_all.to_csv(
        tables_dir / "profile_predictor_results.csv", index=False,
    )
    profile_diag_all.to_csv(
        tables_dir / "profile_predictor_diagnostics.csv", index=False,
    )

    # ------------------------------------------------------------------
    # Diagnostic plots
    # ------------------------------------------------------------------
    print("[chapter1] writing diagnostic plots", flush=True)
    plot_main_effects_forest(
        main_results,
        figures_dir / "main_effects_forest",
        title="Excess sociodemographic mixing → cluster scale (main effects)",
    )
    plot_wave_interaction_slopes(
        wave_results,
        figures_dir / "wave_interaction_slopes",
        title="Excess-mixing slopes on cluster scale, by epidemic wave",
    )

    print(f"[chapter1] tables written under {tables_dir}", flush=True)
    print(f"[chapter1] figures written under {figures_dir}", flush=True)


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
    parser.add_argument(
        "--winsorise-quantile",
        type=float,
        default=0.0,
        metavar="Q",
        help=(
            "Winsorise each ZTNB outcome at the Q-th quantile before fitting. "
            "Use 0 (default) to disable."
        ),
    )
    parser.add_argument(
        "--exclude-tail-quantile",
        type=float,
        default=0.0,
        metavar="Q",
        help=(
            "Exclude rows above the Q-th outcome quantile before each ZTNB fit. "
            "For example, 0.995 excludes the top 0.5%% tail. Use 0 to disable."
        ),
    )
    parser.add_argument(
        "--window-stride",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Retain only windows where window_idx %% N == 0. "
            "Default 1 keeps all sliding windows."
        ),
    )
    parser.add_argument("--tables-dir", type=Path, default=None)
    parser.add_argument("--figures-dir", type=Path, default=None)
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument(
        "--sample-clusters", type=int, default=None,
        help="Subsample the cluster table to this many rows (quick-test).",
    )
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
        winsorise_quantile=args.winsorise_quantile,
        exclude_tail_quantile=args.exclude_tail_quantile,
        window_stride=args.window_stride,
        tables_dir=args.tables_dir.resolve() if args.tables_dir else None,
        figures_dir=args.figures_dir.resolve() if args.figures_dir else None,
        cache_dir=args.cache_dir.resolve() if args.cache_dir else None,
        sample_clusters=args.sample_clusters,
    )


if __name__ == "__main__":
    main()
