"""Primary Part 1 analysis — Line 1 (deprivation) + Line 2 (mixing) overall.

This script orchestrates the main two-part Part 1 analysis.  It uses one
Leiden resolution to avoid treating repeated cluster resolutions as
independent, then fits:

* **Line 1** — area-level SIMD deprivation (plus surveillance covariates) as
  predictors of cluster outcomes:

  * Hurdle (binomial) and zero-truncated negative-binomial (positive count)
    components of cluster size and geographic dispersion
    (:func:`lib.fit_models.fit_count_models`).
  * Linear models for observed-minus-expected within-cluster excess mixing
    (:func:`lib.fit_models.fit_mixing_models`).

* **Line 2** — cluster-level excess-mixing scores as predictors of cluster
  outcomes, controlling for the Line-1 covariates
  (:func:`lib.fit_models.fit_mixing_predictor_count_models`).

The actual fitting routines and data-prep helpers live in the ``lib``
subpackage; this script is a thin orchestrator that mirrors the structure
of ``manuscript/make_figures.py``.

Run from the repository root with::

    conda run -n PhD python part1/overall_analysis.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

# Allow ``python part1/overall_analysis.py`` to import the local lib package
# even when this script is run from elsewhere.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.constants import (  # noqa: E402
    CALENDAR_SPLINE_DF,
    COUNT_MODEL_SPECS,
    CountModelSpec,
    LINEAGE_MIN_CLUSTERS,
    MIXING_PREDICTOR_TERMS,
    MIXING_VARIABLES,
    PRIMARY_RESOLUTION,
    PRIMARY_TERMS,
    QC_DEFAULT,
    SEQUENCE_COLUMNS,
    TERM_LABELS,
)
from lib.data_prep import (  # noqa: E402
    analysis_dataset_path,
    build_cluster_table,
    ensure_mixing_predictor_columns,
    expected_stratum_discordance,
    load_analysis_columns_pandas,
    load_cluster_table,
    load_simd_columns_pandas,
    observed_cluster_discordance,
    read_sequence_rows,
    repo_root,
    summarise_dataset,
)
from lib.estimators import (  # noqa: E402
    build_exog,
    fit_ztnb,
    lineage_levels,
    logit_clipped,
    zscore,
)
from lib.fit_models import (  # noqa: E402
    fit_binary_component,
    fit_count_models,
    fit_mixing_models,
    fit_mixing_predictor_count_models,
    fit_positive_component,
)
from lib.inspect_plots import (  # noqa: E402
    load_plot_style,
    plot_count_effects,
    plot_mixing_effects,
    plot_mixing_predictor_count_effects,
    setup_matplotlib_cache,
    term_colours,
)


# ---------------------------------------------------------------------------
# Backward-compatible re-exports
# ---------------------------------------------------------------------------
#
# Other scripts (``loglinear_sensitivity.py``,
# ``wave_analysis.py``, ``domain_analysis.py``) and
# the manuscript figure script import names from this module.  Keep those
# names available here so callers do not need to update their imports.

__all__ = [
    # constants
    "CALENDAR_SPLINE_DF", "COUNT_MODEL_SPECS", "CountModelSpec",
    "LINEAGE_MIN_CLUSTERS", "MIXING_PREDICTOR_TERMS", "MIXING_VARIABLES",
    "PRIMARY_RESOLUTION", "PRIMARY_TERMS", "QC_DEFAULT", "SEQUENCE_COLUMNS",
    "TERM_LABELS",
    # data prep
    "analysis_dataset_path", "build_cluster_table",
    "ensure_mixing_predictor_columns", "expected_stratum_discordance",
    "load_analysis_columns_pandas", "load_cluster_table",
    "load_simd_columns_pandas", "observed_cluster_discordance",
    "read_sequence_rows", "repo_root", "summarise_dataset",
    # estimators
    "build_exog", "fit_ztnb", "lineage_levels", "logit_clipped", "zscore",
    # fit pipelines
    "fit_binary_component", "fit_count_models", "fit_mixing_models",
    "fit_mixing_predictor_count_models", "fit_positive_component",
    # quick-look plots
    "load_plot_style", "plot_count_effects", "plot_mixing_effects",
    "plot_mixing_predictor_count_effects", "setup_matplotlib_cache",
    "term_colours",
    # entry point
    "main", "parse_args", "run",
]


# ---------------------------------------------------------------------------
# Top-level run
# ---------------------------------------------------------------------------


def _resolve_primary_terms(use_index_simd: bool) -> list[str]:
    """Return the active primary covariate set without mutating any global."""
    if not use_index_simd:
        return list(PRIMARY_TERMS)
    return [
        "index_deprivation_z" if t == "deprivation_z" else t
        for t in PRIMARY_TERMS
    ]


def run(
    root: Path,
    qc: str | None,
    primary_resolution: float,
    lineage_min_clusters: int,
    calendar_spline_df: int,
    maxiter: int,
    cluster_by: str = "window_id",
    use_size_offset: bool = False,
    winsorise_quantile: float = 0.0,
    use_index_simd: bool = False,
    window_stride: int = 1,
    tables_dir: Path | None = None,
    figures_dir: Path | None = None,
    cache_dir: Path | None = None,
) -> None:
    out_dir = root / "part1"
    tables_dir = tables_dir or out_dir / "tables"
    figures_dir = figures_dir or out_dir / "figures"
    cache_dir = cache_dir or out_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    path = analysis_dataset_path(root)
    print(f"Reading primary-resolution sequence rows from {path}", flush=True)
    seq = read_sequence_rows(qc=qc, primary_resolution=primary_resolution)
    print(f"Building cluster table from {len(seq):,} sequence rows", flush=True)
    clusters, scaling, dropped = build_cluster_table(
        seq,
        lineage_min_clusters=lineage_min_clusters,
        calendar_spline_df=calendar_spline_df,
    )

    # Non-overlapping window sensitivity: retain only every Nth window so
    # successive analysis windows do not share sequences.  ``window_stride=1``
    # (the default) keeps all windows and is the primary analysis.
    if window_stride > 1:
        keep_idx = clusters["window_idx"] % window_stride == 0
        n_before = len(clusters)
        clusters = clusters.loc[keep_idx].copy()
        print(
            f"Non-overlapping window filter (stride={window_stride}): "
            f"{n_before:,} → {len(clusters):,} clusters",
            flush=True,
        )

    # Index-case SIMD sensitivity: swap deprivation_z for index_deprivation_z
    # in the primary covariate list — threaded through the fit functions so
    # the module-level PRIMARY_TERMS is never mutated.
    primary_terms = _resolve_primary_terms(use_index_simd)
    if use_index_simd:
        print(
            "Using index-case SIMD (index_deprivation_z) instead of mean cluster SIMD.",
            flush=True,
        )

    calendar_cols = [col for col in clusters.columns if col.startswith("calendar_spline_")]
    lineage_levels_all = lineage_levels(clusters)

    clusters.to_parquet(cache_dir / "cluster_table.parquet", index=False)
    scaling.to_csv(tables_dir / "covariate_scaling.csv", index=False)
    descriptives = summarise_dataset(seq, clusters, qc, primary_resolution, dropped, scaling)
    descriptives.to_csv(tables_dir / "dataset_descriptives.csv", index=False)
    print(
        f"Fitting count models for {len(clusters):,} clusters, "
        f"{len(lineage_levels_all)} lineage model levels, "
        f"{len(calendar_cols)} calendar spline terms",
        flush=True,
    )

    # ------------------------------------------------------------------
    # Line 1 — deprivation as exposure
    # ------------------------------------------------------------------
    count_results, count_diagnostics = fit_count_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
        maxiter=maxiter,
        cluster_by=cluster_by,
        use_size_offset=use_size_offset,
        winsorise_quantile=winsorise_quantile,
        primary_terms=primary_terms,
    )
    count_results.to_csv(tables_dir / "hurdle_count_model_results.csv", index=False)
    count_diagnostics.to_csv(tables_dir / "hurdle_count_model_diagnostics.csv", index=False)

    print("Fitting mixing models", flush=True)
    mixing_results, mixing_diagnostics = fit_mixing_models(
        clusters,
        lineage_levels_all=lineage_levels_all,
        calendar_cols=calendar_cols,
        cluster_by=cluster_by,
        primary_terms=primary_terms,
    )
    mixing_results.to_csv(tables_dir / "mixing_model_results.csv", index=False)
    mixing_diagnostics.to_csv(tables_dir / "mixing_model_diagnostics.csv", index=False)

    # ------------------------------------------------------------------
    # Line 2 — excess mixing as predictor
    # ------------------------------------------------------------------
    print("Fitting count models with mixing predictors", flush=True)
    mixing_predictor_count_results, mixing_predictor_count_diagnostics = (
        fit_mixing_predictor_count_models(
            clusters,
            lineage_levels_all=lineage_levels_all,
            calendar_cols=calendar_cols,
            maxiter=maxiter,
            cluster_by=cluster_by,
            use_size_offset=use_size_offset,
            winsorise_quantile=winsorise_quantile,
            primary_terms=primary_terms,
        )
    )
    mixing_predictor_count_results.to_csv(
        tables_dir / "mixing_predictor_hurdle_count_model_results.csv",
        index=False,
    )
    mixing_predictor_count_diagnostics.to_csv(
        tables_dir / "mixing_predictor_hurdle_count_model_diagnostics.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # Quick-look summary plots (inspection, not publication)
    # ------------------------------------------------------------------
    plot_count_effects(
        count_results, figures_dir / "hurdle_count_effects",
        primary_terms=primary_terms,
    )
    plot_mixing_predictor_count_effects(
        mixing_predictor_count_results,
        figures_dir / "mixing_predictor_hurdle_count_effects",
    )
    plot_mixing_effects(
        mixing_results, figures_dir / "mixing_effects",
        primary_terms=primary_terms,
    )

    print(f"Wrote {tables_dir / 'hurdle_count_model_results.csv'}", flush=True)
    print(
        f"Wrote {tables_dir / 'mixing_predictor_hurdle_count_model_results.csv'}",
        flush=True,
    )
    print(f"Wrote {tables_dir / 'mixing_model_results.csv'}", flush=True)
    print(f"Wrote {figures_dir / 'hurdle_count_effects.png'}", flush=True)
    print(
        f"Wrote {figures_dir / 'mixing_predictor_hurdle_count_effects.png'}",
        flush=True,
    )
    print(f"Wrote {figures_dir / 'mixing_effects.png'}", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument(
        "--qc",
        default=QC_DEFAULT,
        help="Nextclade QC status to retain. Use 'none' to disable QC filtering.",
    )
    parser.add_argument("--primary-resolution", type=float, default=PRIMARY_RESOLUTION)
    parser.add_argument("--lineage-min-clusters", type=int, default=LINEAGE_MIN_CLUSTERS)
    parser.add_argument("--calendar-spline-df", type=int, default=CALENDAR_SPLINE_DF)
    parser.add_argument("--maxiter", type=int, default=1000)
    parser.add_argument(
        "--cluster-by",
        default="window_id",
        choices=["window_id", "health_board"],
        help=(
            "Column used to define clusters for the sandwich standard-error estimator. "
            "'window_id' (default, primary analysis) clusters by sliding analysis window "
            "to account for temporal dependency. 'health_board' clusters by NHS Health "
            "Board (14 groups) to account for spatial dependency, as pre-specified in the SAP."
        ),
    )
    parser.add_argument(
        "--use-size-offset",
        action="store_true",
        default=False,
        help=(
            "Include log(wn_no_sequences) as an offset in the cluster-size hurdle "
            "and positive count model."
        ),
    )
    parser.add_argument(
        "--winsorise-quantile",
        type=float,
        default=0.0,
        metavar="Q",
        help=(
            "Winsorise the positive count outcome at the Q-th quantile before fitting "
            "the ZTNB. Set to 0 (default) to disable."
        ),
    )
    parser.add_argument(
        "--use-index-simd",
        action="store_true",
        default=False,
        help=(
            "Sensitivity S2: replace mean cluster SIMD rank with the SIMD rank of the "
            "index case (earliest collection date) as the deprivation exposure."
        ),
    )
    parser.add_argument(
        "--window-stride",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Retain only windows where window_idx %% N == 0. Default 1 keeps all windows."
        ),
    )
    parser.add_argument("--tables-dir", type=Path, default=None)
    parser.add_argument("--figures-dir", type=Path, default=None)
    parser.add_argument("--cache-dir", type=Path, default=None)
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
        cluster_by=args.cluster_by,
        use_size_offset=args.use_size_offset,
        winsorise_quantile=args.winsorise_quantile,
        use_index_simd=args.use_index_simd,
        window_stride=args.window_stride,
        tables_dir=args.tables_dir.resolve() if args.tables_dir else None,
        figures_dir=args.figures_dir.resolve() if args.figures_dir else None,
        cache_dir=args.cache_dir.resolve() if args.cache_dir else None,
    )


if __name__ == "__main__":
    main()
