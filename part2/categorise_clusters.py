"""Shared Part 2 cluster categorisation helpers and cached-table CLI.

The script reads the cached main cluster table and adds:

* cluster size categories: small/moderate, large, very large
* geographic dispersion categories using distinct-datazone counts
* mean-cluster SIMD quintile categories
* SIMD, age, sex, and joint-profile excess-mixing categories

Run from the repository root with:

    conda run -n PhD python part2/categorise_clusters.py
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def _bootstrap_repo_root_for_utils() -> Path:
    """Ensure the repository root is importable for direct script execution."""
    here = Path(__file__).resolve()
    for candidate in [here, *here.parents]:
        if (candidate / "config.yaml").exists():
            root_str = str(candidate)
            if root_str not in sys.path:
                sys.path.insert(0, root_str)
            return candidate
    raise FileNotFoundError("Could not locate config.yaml.")


_bootstrap_repo_root_for_utils()

from utils.data import (  # noqa: E402
    load_main_cluster_table,
    main_cluster_table_path,
    repo_root as data_repo_root,
)


DEFAULT_SIMD_RANK_MAX = 6976.0

SIZE_LABELS = ("small/moderate", "large", "very large")
GEOGRAPHY_LABELS = (
    "low/moderate dispersion",
    "large dispersion",
    "very large dispersion",
)
MIXING_LABELS = ("less mix", "baseline", "more mix", "not available")

MIXING_SPECS = {
    "simd": {
        "source": "simd_excess_discordance",
        "category": "simd_mixing_category",
        "label": "SIMD quintile mixing",
    },
    "age": {
        "source": "age_excess_discordance",
        "category": "age_mixing_category",
        "label": "Age-band mixing",
    },
    "sex": {
        "source": "sex_excess_discordance",
        "category": "sex_mixing_category",
        "label": "Sex mixing",
    },
    "profile": {
        "source": "profile_excess_discordance",
        "category": "profile_mixing_category",
        "label": "Joint SIMD-age-sex profile mixing",
    },
}

REQUIRED_COLUMNS = {
    "cluster_id",
    "cluster_size",
    "cluster_n_datazones",
    "mean_simd_rank",
    *(spec["source"] for spec in MIXING_SPECS.values()),
}

PREFERRED_COMPACT_COLUMNS = [
    "cluster_id",
    "window_id",
    "window_idx",
    "wn_mid_date",
    "pango_lineage",
    "lineage_model",
    "cluster_size",
    "cluster_size_category",
    "cluster_n_datazones",
    "geographic_dispersion_category",
    "mean_simd_rank",
    "simd_quintile",
    "simd_quintile_label",
    "simd_excess_discordance",
    "simd_mixing_category",
    "age_excess_discordance",
    "age_mixing_category",
    "sex_excess_discordance",
    "sex_mixing_category",
    "profile_excess_discordance",
    "profile_mixing_category",
]


def resolve_repo_path(root: Path, path: Path | None, default: Path) -> Path:
    """Return an absolute path, resolving relative paths under ``root``."""
    target = path if path is not None else default
    return target if target.is_absolute() else root / target


def validate_probability(value: float, name: str) -> float:
    if not 0.0 < value < 1.0:
        raise ValueError(f"{name} must be between 0 and 1, got {value!r}.")
    return value


def validate_columns(clusters: pd.DataFrame) -> None:
    missing = sorted(REQUIRED_COLUMNS - set(clusters.columns))
    if missing:
        joined = ", ".join(missing)
        raise KeyError(f"Cluster table is missing required columns: {joined}")


def threshold_from_quantile(values: pd.Series, quantile: float) -> int:
    finite = pd.to_numeric(values, errors="coerce").dropna()
    if finite.empty:
        raise ValueError(f"Cannot compute quantile threshold for empty {values.name!r}.")
    return max(1, int(math.ceil(float(finite.quantile(quantile)))))


def resolve_thresholds(
    values: pd.Series,
    *,
    large_min: int | None,
    very_large_min: int | None,
    large_quantile: float,
    very_large_quantile: float,
) -> tuple[int, int]:
    """Resolve integer thresholds for large and very-large categories."""
    large_quantile = validate_probability(large_quantile, "large_quantile")
    very_large_quantile = validate_probability(very_large_quantile, "very_large_quantile")
    if large_quantile >= very_large_quantile:
        raise ValueError("large_quantile must be smaller than very_large_quantile.")

    resolved_large = (
        int(large_min)
        if large_min is not None
        else threshold_from_quantile(values, large_quantile)
    )
    resolved_very_large = (
        int(very_large_min)
        if very_large_min is not None
        else threshold_from_quantile(values, very_large_quantile)
    )

    if resolved_large < 1 or resolved_very_large < 1:
        raise ValueError("Category thresholds must be positive integers.")
    if resolved_large >= resolved_very_large:
        raise ValueError(
            "Large-category threshold must be smaller than the very-large threshold "
            f"(got {resolved_large} and {resolved_very_large})."
        )
    return resolved_large, resolved_very_large


def categorise_count(
    values: pd.Series,
    *,
    large_min: int,
    very_large_min: int,
    labels: Iterable[str],
) -> pd.Categorical:
    ordered_labels = list(labels)
    numeric = pd.to_numeric(values, errors="coerce")
    out = pd.Series(pd.NA, index=values.index, dtype="object")
    out[numeric < large_min] = ordered_labels[0]
    out[(numeric >= large_min) & (numeric < very_large_min)] = ordered_labels[1]
    out[numeric >= very_large_min] = ordered_labels[2]
    return pd.Categorical(out, categories=ordered_labels, ordered=True)


def simd_quintile_from_rank(rank: pd.Series, rank_max: float) -> pd.Series:
    if rank_max <= 0:
        raise ValueError("simd_rank_max must be positive.")
    numeric = pd.to_numeric(rank, errors="coerce")
    quintile = np.ceil(numeric / (rank_max / 5.0))
    quintile = pd.Series(quintile, index=rank.index).clip(lower=1, upper=5)
    return quintile.astype("Int64")


def categorise_mixing(values: pd.Series, baseline_band: float) -> pd.Categorical:
    """Classify excess discordance relative to zero expected excess mixing."""
    if baseline_band < 0:
        raise ValueError("mixing_baseline_band must be non-negative.")

    numeric = pd.to_numeric(values, errors="coerce")
    out = pd.Series("not available", index=values.index, dtype="object")
    out[numeric < -baseline_band] = "less mix"
    out[numeric.abs() <= baseline_band] = "baseline"
    out[numeric > baseline_band] = "more mix"
    return pd.Categorical(out, categories=list(MIXING_LABELS), ordered=True)


def add_categories(
    clusters: pd.DataFrame,
    *,
    large_size_min: int,
    very_large_size_min: int,
    large_geography_min: int,
    very_large_geography_min: int,
    simd_rank_max: float,
    mixing_baseline_band: float,
) -> pd.DataFrame:
    categorised = clusters.copy()
    categorised["cluster_size_category"] = categorise_count(
        categorised["cluster_size"],
        large_min=large_size_min,
        very_large_min=very_large_size_min,
        labels=SIZE_LABELS,
    )
    categorised["geographic_dispersion_category"] = categorise_count(
        categorised["cluster_n_datazones"],
        large_min=large_geography_min,
        very_large_min=very_large_geography_min,
        labels=GEOGRAPHY_LABELS,
    )

    categorised["simd_quintile"] = simd_quintile_from_rank(
        categorised["mean_simd_rank"],
        rank_max=simd_rank_max,
    )
    simd_labels = {
        1: "1 most deprived",
        2: "2",
        3: "3",
        4: "4",
        5: "5 least deprived",
    }
    categorised["simd_quintile_label"] = pd.Categorical(
        categorised["simd_quintile"].map(simd_labels),
        categories=[simd_labels[i] for i in range(1, 6)],
        ordered=True,
    )

    for spec in MIXING_SPECS.values():
        categorised[spec["category"]] = categorise_mixing(
            categorised[spec["source"]],
            baseline_band=mixing_baseline_band,
        )
    return categorised


def summarise_categories(clusters: pd.DataFrame, category_columns: list[str]) -> pd.DataFrame:
    rows = []
    total = len(clusters)
    for col in category_columns:
        counts = clusters[col].value_counts(dropna=False, sort=False)
        for category, count in counts.items():
            rows.append(
                {
                    "category_variable": col,
                    "category": category,
                    "n_clusters": int(count),
                    "fraction_clusters": float(count / total) if total else np.nan,
                }
            )
    return pd.DataFrame(rows)


def build_threshold_table(
    *,
    large_size_min: int,
    very_large_size_min: int,
    large_size_quantile: float,
    very_large_size_quantile: float,
    large_geography_min: int,
    very_large_geography_min: int,
    large_geography_quantile: float,
    very_large_geography_quantile: float,
    simd_rank_max: float,
    mixing_baseline_band: float,
) -> pd.DataFrame:
    rows = [
        {
            "category_variable": "cluster_size_category",
            "source_column": "cluster_size",
            "category": "small/moderate",
            "rule": f"cluster_size < {large_size_min}",
            "default_basis": (
                f"non-singleton {large_size_quantile:.2f} and "
                f"{very_large_size_quantile:.2f} quantiles"
            ),
        },
        {
            "category_variable": "cluster_size_category",
            "source_column": "cluster_size",
            "category": "large",
            "rule": f"{large_size_min} <= cluster_size < {very_large_size_min}",
            "default_basis": (
                f"non-singleton {large_size_quantile:.2f} and "
                f"{very_large_size_quantile:.2f} quantiles"
            ),
        },
        {
            "category_variable": "cluster_size_category",
            "source_column": "cluster_size",
            "category": "very large",
            "rule": f"cluster_size >= {very_large_size_min}",
            "default_basis": (
                f"non-singleton {large_size_quantile:.2f} and "
                f"{very_large_size_quantile:.2f} quantiles"
            ),
        },
        {
            "category_variable": "geographic_dispersion_category",
            "source_column": "cluster_n_datazones",
            "category": "low/moderate dispersion",
            "rule": f"cluster_n_datazones < {large_geography_min}",
            "default_basis": (
                f"non-singleton {large_geography_quantile:.2f} and "
                f"{very_large_geography_quantile:.2f} quantiles"
            ),
        },
        {
            "category_variable": "geographic_dispersion_category",
            "source_column": "cluster_n_datazones",
            "category": "large dispersion",
            "rule": (
                f"{large_geography_min} <= cluster_n_datazones < "
                f"{very_large_geography_min}"
            ),
            "default_basis": (
                f"non-singleton {large_geography_quantile:.2f} and "
                f"{very_large_geography_quantile:.2f} quantiles"
            ),
        },
        {
            "category_variable": "geographic_dispersion_category",
            "source_column": "cluster_n_datazones",
            "category": "very large dispersion",
            "rule": f"cluster_n_datazones >= {very_large_geography_min}",
            "default_basis": (
                f"non-singleton {large_geography_quantile:.2f} and "
                f"{very_large_geography_quantile:.2f} quantiles"
            ),
        },
    ]

    for quintile in range(1, 6):
        lower = (quintile - 1) * simd_rank_max / 5.0
        upper = quintile * simd_rank_max / 5.0
        rule = (
            f"1 <= mean_simd_rank <= {upper:.1f}"
            if quintile == 1
            else f"{lower:.1f} < mean_simd_rank <= {upper:.1f}"
        )
        rows.append(
            {
                "category_variable": "simd_quintile",
                "source_column": "mean_simd_rank",
                "category": quintile,
                "rule": rule,
                "default_basis": f"rank max {simd_rank_max:.0f}; quintile 1 is most deprived",
            }
        )

    for name, spec in MIXING_SPECS.items():
        rows.extend(
            [
                {
                    "category_variable": spec["category"],
                    "source_column": spec["source"],
                    "category": "less mix",
                    "rule": f"{spec['source']} < -{mixing_baseline_band}",
                    "default_basis": (
                        "observed-minus-expected excess discordance; "
                        "negative means less mixed than expected"
                    ),
                },
                {
                    "category_variable": spec["category"],
                    "source_column": spec["source"],
                    "category": "baseline",
                    "rule": (
                        f"abs({spec['source']}) <= {mixing_baseline_band}"
                    ),
                    "default_basis": (
                        "observed-minus-expected excess discordance; "
                        f"baseline band is +/-{mixing_baseline_band}"
                    ),
                },
                {
                    "category_variable": spec["category"],
                    "source_column": spec["source"],
                    "category": "more mix",
                    "rule": f"{spec['source']} > {mixing_baseline_band}",
                    "default_basis": (
                        "observed-minus-expected excess discordance; "
                        "positive means more mixed than expected"
                    ),
                },
                {
                    "category_variable": spec["category"],
                    "source_column": spec["source"],
                    "category": "not available",
                    "rule": f"{spec['source']} is missing",
                    "default_basis": f"{spec['label']} is undefined for singleton/invalid clusters",
                },
            ]
        )
        rows[-1]["mixing_dimension"] = name
        rows[-2]["mixing_dimension"] = name
        rows[-3]["mixing_dimension"] = name
        rows[-4]["mixing_dimension"] = name

    return pd.DataFrame(rows)


def build_combination_summary(clusters: pd.DataFrame) -> pd.DataFrame:
    group_cols = [
        "cluster_size_category",
        "geographic_dispersion_category",
        "simd_quintile",
        "simd_quintile_label",
        "simd_mixing_category",
        "age_mixing_category",
        "sex_mixing_category",
        "profile_mixing_category",
    ]
    summary = (
        clusters.groupby(group_cols, observed=True, dropna=False)
        .agg(
            n_clusters=("cluster_id", "size"),
            mean_cluster_size=("cluster_size", "mean"),
            median_cluster_size=("cluster_size", "median"),
            mean_datazones=("cluster_n_datazones", "mean"),
            median_datazones=("cluster_n_datazones", "median"),
        )
        .reset_index()
        .sort_values("n_clusters", ascending=False)
    )
    summary["fraction_clusters"] = summary["n_clusters"] / len(clusters)
    return summary


def compact_category_table(clusters: pd.DataFrame) -> pd.DataFrame:
    columns = [col for col in PREFERRED_COMPACT_COLUMNS if col in clusters.columns]
    return clusters.loc[:, columns].copy()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Categorise cached Part 1 clusters for Part 2.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing main_cluster_table.parquet. Relative paths are "
            "resolved from the repository root. Default: part1/main/cache."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Default: part2.",
    )
    parser.add_argument(
        "--large-size-min",
        type=int,
        default=None,
        help="Minimum cluster_size for the large category. Default: ceil(size 90th percentile).",
    )
    parser.add_argument(
        "--very-large-size-min",
        type=int,
        default=None,
        help=(
            "Minimum cluster_size for the very large category. "
            "Default: ceil(size 99th percentile)."
        ),
    )
    parser.add_argument(
        "--large-size-quantile",
        type=float,
        default=0.90,
        help="Quantile used for default large-size threshold.",
    )
    parser.add_argument(
        "--very-large-size-quantile",
        type=float,
        default=0.99,
        help="Quantile used for default very-large-size threshold.",
    )
    parser.add_argument(
        "--large-geography-min",
        type=int,
        default=None,
        help=(
            "Minimum cluster_n_datazones for the large geographic-dispersion "
            "category. Default: ceil(datazone-count 90th percentile)."
        ),
    )
    parser.add_argument(
        "--very-large-geography-min",
        type=int,
        default=None,
        help=(
            "Minimum cluster_n_datazones for the very-large geographic-dispersion "
            "category. Default: ceil(datazone-count 99th percentile)."
        ),
    )
    parser.add_argument(
        "--large-geography-quantile",
        type=float,
        default=0.90,
        help="Quantile used for default large geographic-dispersion threshold.",
    )
    parser.add_argument(
        "--very-large-geography-quantile",
        type=float,
        default=0.99,
        help="Quantile used for default very-large geographic-dispersion threshold.",
    )
    parser.add_argument(
        "--simd-rank-max",
        type=float,
        default=DEFAULT_SIMD_RANK_MAX,
        help=(
            "Maximum SIMD rank used to convert mean cluster rank to quintiles. "
            "Default: 6976."
        ),
    )
    parser.add_argument(
        "--mixing-baseline-band",
        type=float,
        default=0.01,
        help=(
            "Absolute excess-discordance band classified as baseline. "
            "Default: 0.01, i.e. +/-1 percentage point."
        ),
    )
    parser.add_argument(
        "--no-combinations",
        action="store_true",
        help="Skip the full cross-category combination summary.",
    )
    return parser.parse_args()


def run(args: argparse.Namespace) -> None:
    root = data_repo_root()
    input_cache_dir = (
        resolve_repo_path(root, args.cache_dir, Path("part1/main/cache"))
        if args.cache_dir
        else None
    )
    cluster_cache = main_cluster_table_path(root=root, cache_dir=input_cache_dir)
    out_dir = resolve_repo_path(root, args.out_dir, Path("part2"))
    cache_dir = out_dir / "cache"
    tables_dir = out_dir / "tables"
    cache_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading cached cluster table: {cluster_cache}", flush=True)
    clusters = load_main_cluster_table(root=root, cache_dir=input_cache_dir)
    validate_columns(clusters)
    n_before = len(clusters)
    clusters = clusters.loc[clusters["cluster_size"] > 1].copy()
    if clusters.empty:
        raise ValueError("No non-singleton clusters available for categorisation.")
    print(
        f"Retained {len(clusters):,}/{n_before:,} non-singleton clusters before categorising.",
        flush=True,
    )

    large_size_min, very_large_size_min = resolve_thresholds(
        clusters["cluster_size"],
        large_min=args.large_size_min,
        very_large_min=args.very_large_size_min,
        large_quantile=args.large_size_quantile,
        very_large_quantile=args.very_large_size_quantile,
    )
    large_geography_min, very_large_geography_min = resolve_thresholds(
        clusters["cluster_n_datazones"],
        large_min=args.large_geography_min,
        very_large_min=args.very_large_geography_min,
        large_quantile=args.large_geography_quantile,
        very_large_quantile=args.very_large_geography_quantile,
    )

    print(
        "Using thresholds: "
        f"cluster_size >= {large_size_min} large, >= {very_large_size_min} very large; "
        f"cluster_n_datazones >= {large_geography_min} large dispersion, "
        f">= {very_large_geography_min} very large dispersion; "
        f"mixing baseline +/-{args.mixing_baseline_band}.",
        flush=True,
    )

    categorised = add_categories(
        clusters,
        large_size_min=large_size_min,
        very_large_size_min=very_large_size_min,
        large_geography_min=large_geography_min,
        very_large_geography_min=very_large_geography_min,
        simd_rank_max=args.simd_rank_max,
        mixing_baseline_band=args.mixing_baseline_band,
    )

    category_columns = [
        "cluster_size_category",
        "geographic_dispersion_category",
        "simd_quintile_label",
        *(spec["category"] for spec in MIXING_SPECS.values()),
    ]

    categorised.to_parquet(cache_dir / "cluster_categories.parquet", index=False)
    compact_category_table(categorised).to_csv(
        tables_dir / "cluster_categories.csv",
        index=False,
    )
    summarise_categories(categorised, category_columns).to_csv(
        tables_dir / "cluster_category_summary.csv",
        index=False,
    )
    build_threshold_table(
        large_size_min=large_size_min,
        very_large_size_min=very_large_size_min,
        large_size_quantile=args.large_size_quantile,
        very_large_size_quantile=args.very_large_size_quantile,
        large_geography_min=large_geography_min,
        very_large_geography_min=very_large_geography_min,
        large_geography_quantile=args.large_geography_quantile,
        very_large_geography_quantile=args.very_large_geography_quantile,
        simd_rank_max=args.simd_rank_max,
        mixing_baseline_band=args.mixing_baseline_band,
    ).to_csv(tables_dir / "cluster_category_thresholds.csv", index=False)

    if not args.no_combinations:
        build_combination_summary(categorised).to_csv(
            tables_dir / "cluster_category_combinations.csv",
            index=False,
        )

    print(
        f"Wrote categorised clusters and summaries under {out_dir} "
        f"({len(categorised):,} clusters).",
        flush=True,
    )


if __name__ == "__main__":
    run(parse_args())
