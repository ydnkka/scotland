"""Build Leiden-resolution and sparsification sensitivity tables.

Run from the Scotland repository root:

    python -m chapter_analyses.genomic_networks.build_sensitivity_tables

For a quick development check:

    python -m chapter_analyses.genomic_networks.build_sensitivity_tables --max-windows 3 \
        --max-files 10
"""

from __future__ import annotations

import argparse
import logging
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any, TypedDict, Union, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow.parquet as pq
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

from .lib.config import PROJECT_ROOT, QC_FILTER, SPARSIFICATION_THRESHOLD, TABLES_DIR
from .lib.io import ensure_results_dirs, write_table


LOGGER = logging.getLogger(__name__)

AggregationFunc = Union[str, Callable[[pd.Series], float]]
NamedAggregation = tuple[str, AggregationFunc]


class PairwiseScan(TypedDict):
    total_rows: int
    scan_rows: int
    scan_weight_sum: float
    retained_scan_edges: npt.NDArray[np.int64]
    retained_scan_weight_sums: npt.NDArray[np.float64]


ANALYSIS_DATASET_PATH = (
    PROJECT_ROOT / "data/processed/scotland_clustering_analysis_dataset.parquet"
)
PAIRWISE_DATASET_DIR = PROJECT_ROOT / "data/processed/pairwise_distances_dataset"
EDGE_MANIFEST_PATH = (
    PROJECT_ROOT / "data/processed/sparsified_edge_counts_by_window_lineage.parquet"
)

DEFAULT_SPARSIFICATION_THRESHOLDS = (
    0.0,
    1e-6,
    1e-5,
    1e-4,
    1e-3,
    1e-2,
    5e-2,
    1e-1,
)

LEIDEN_COLUMNS = [
    "resolution",
    "nextclade_qc",
    "window_id",
    "window_idx",
    "sequence_id",
    "cluster_id",
    "cluster_size",
    "cluster_n_datazones",
    "cluster_duration_days",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only",
        choices=("all", "leiden", "sparsification"),
        default="all",
        help="Build only one sensitivity family. Default: all.",
    )
    parser.add_argument(
        "--analysis-dataset",
        type=Path,
        default=ANALYSIS_DATASET_PATH,
        help="Multi-resolution clustering analysis parquet.",
    )
    parser.add_argument(
        "--pairwise-dir",
        type=Path,
        default=PAIRWISE_DATASET_DIR,
        help="Directory of window-lineage pairwise parquet files.",
    )
    parser.add_argument(
        "--edge-manifest",
        type=Path,
        default=EDGE_MANIFEST_PATH,
        help="Saved pairwise edge-count manifest used for metadata.",
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=TABLES_DIR,
        help="Output table directory. Default: chapter_analyses/genomic_networks/results/tables.",
    )
    parser.add_argument(
        "--baseline-resolution",
        type=float,
        default=0.3,
        help="Reference Leiden resolution. Default: EpiLink manuscript main analysis resolution.",
    )
    parser.add_argument(
        "--include-ami",
        action="store_true",
        help=(
            "Add exact adjusted mutual information columns for Leiden sensitivity. "
            "This is substantially slower than ARI on large high-resolution windows."
        ),
    )
    parser.add_argument(
        "--baseline-threshold",
        type=float,
        default=SPARSIFICATION_THRESHOLD,
        help="Reference EpiLink compatibility threshold. Default: Chapter 4 baseline.",
    )
    parser.add_argument(
        "--thresholds",
        nargs="+",
        type=float,
        default=list(DEFAULT_SPARSIFICATION_THRESHOLDS),
        help="Sparsification thresholds to evaluate. Baseline is added if absent.",
    )
    parser.add_argument(
        "--qc",
        default=QC_FILTER,
        help="Nextclade QC value(s), comma-separated. Use 'all' to skip QC filtering.",
    )
    parser.add_argument(
        "--windows",
        nargs="*",
        default=None,
        help="Optional window IDs or indices to retain, e.g. W080 81.",
    )
    parser.add_argument(
        "--pango-lineages",
        nargs="*",
        default=None,
        help="Optional Pango lineage filters for sparsification scans.",
    )
    parser.add_argument(
        "--max-windows",
        type=int,
        default=None,
        help="Development cap on the first N selected windows.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Development cap on the first N selected pairwise files.",
    )
    parser.add_argument(
        "--max-row-groups-per-file",
        type=int,
        default=None,
        help="Approximate sparsification scan by reading only first N row groups.",
    )
    parser.add_argument(
        "--max-rows-per-file",
        type=int,
        default=None,
        help="Approximate sparsification scan by reading only first N rows per file.",
    )
    parser.add_argument(
        "--score-column",
        default="epilink_compatibility",
        help="Pairwise score column used for sparsification. Default: epilink_compatibility.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    )
    return parser.parse_args()


def _normalise_window(value: str | int) -> str:
    text = str(value).strip()
    upper = text.upper()
    if upper.startswith("W") and upper[1:].isdigit():
        return f"W{int(upper[1:]):03d}"
    if upper.isdigit():
        return f"W{int(upper):03d}"
    return text


def _normalise_windows(values: Sequence[str] | None) -> list[str] | None:
    if values is None:
        return None
    return [_normalise_window(value) for value in values]


def _normalise_qc(qc: str | None) -> set[str] | None:
    if qc is None:
        return None
    text = str(qc).strip()
    if not text or text.lower() in {"all", "none"}:
        return None
    return {part.strip() for part in text.split(",") if part.strip()}


def _float_key(value: float) -> float:
    return float(f"{float(value):.12g}")


def _sorted_unique_float(values: Iterable[float]) -> list[float]:
    out = sorted({_float_key(value) for value in values})
    if any(value < 0 for value in out):
        raise ValueError("Sparsification thresholds must be non-negative.")
    return out


def _threshold_label(value: float) -> str:
    return f"{value:g}"


def _as_float(value: object) -> float:
    return float(cast(Any, value))


def _quantile(probability: float) -> Callable[[pd.Series], float]:
    def inner(values: pd.Series) -> float:
        return float(values.quantile(probability))

    inner.__name__ = f"q{int(probability * 100):02d}"
    return inner


def _apply_window_filters(
    df: pd.DataFrame,
    *,
    windows: Sequence[str] | None,
    max_windows: int | None,
) -> pd.DataFrame:
    out = df
    if windows is not None:
        wanted = set(windows)
        out = out.loc[out["window_id"].isin(wanted)]
    if max_windows is not None:
        keep = sorted(out["window_id"].dropna().unique())[:max_windows]
        out = out.loc[out["window_id"].isin(keep)]
    return out.copy()


def build_leiden_sensitivity_tables(
    *,
    analysis_dataset: Path,
    baseline_resolution: float,
    include_ami: bool,
    qc: set[str] | None,
    windows: Sequence[str] | None,
    max_windows: int | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return window-level and resolution-level Leiden sensitivity tables."""
    LOGGER.info("Loading multi-resolution clustering rows from %s", analysis_dataset)
    df = pd.read_parquet(analysis_dataset, columns=LEIDEN_COLUMNS)
    if qc is not None:
        df = df.loc[df["nextclade_qc"].isin(qc)]
    df = _apply_window_filters(df, windows=windows, max_windows=max_windows)
    if df.empty:
        raise ValueError("No clustering rows remain after sensitivity filters.")

    LOGGER.info(
        "Summarising %s rows across %s resolutions and %s windows",
        f"{len(df):,}",
        df["resolution"].nunique(),
        df["window_id"].nunique(),
    )
    cluster_rows = df.drop_duplicates(
        ["resolution", "window_id", "window_idx", "cluster_id"]
    )
    window = (
        cluster_rows.groupby(["resolution", "window_id", "window_idx"], dropna=False)
        .agg(
            n_clusters=("cluster_id", "nunique"),
            n_sequence_memberships=("cluster_size", "sum"),
            median_cluster_size=("cluster_size", "median"),
            p90_cluster_size=("cluster_size", _quantile(0.90)),
            max_cluster_size=("cluster_size", "max"),
            singleton_clusters=("cluster_size", lambda x: int((x == 1).sum())),
            median_duration_days=("cluster_duration_days", "median"),
            median_datazones=("cluster_n_datazones", "median"),
            max_datazones=("cluster_n_datazones", "max"),
        )
        .reset_index()
    )
    sequence_counts = (
        df.groupby(["resolution", "window_id", "window_idx"], dropna=False)[
            "sequence_id"
        ]
        .nunique()
        .rename("n_sequences")
        .reset_index()
    )
    window = window.merge(
        sequence_counts,
        on=["resolution", "window_id", "window_idx"],
        how="left",
    )
    window["clusters_per_1000_sequences"] = (
        1000 * window["n_clusters"] / window["n_sequences"].replace(0, np.nan)
    )
    window["singleton_cluster_share"] = window["singleton_clusters"] / window[
        "n_clusters"
    ].replace(0, np.nan)
    window["singleton_sequence_share"] = window["singleton_clusters"] / window[
        "n_sequences"
    ].replace(0, np.nan)

    ari = _build_leiden_ari_table(
        df,
        baseline_resolution=baseline_resolution,
        include_ami=include_ami,
    )
    window = window.merge(
        ari,
        on=["resolution", "window_id"],
        how="left",
    )
    window = _add_leiden_baseline_contrasts(
        window,
        baseline_resolution=baseline_resolution,
    )
    window = window.sort_values(["resolution", "window_idx"]).reset_index(drop=True)
    summary = _build_leiden_summary(window, include_ami=include_ami)
    return window, summary


def _build_leiden_ari_table(
    df: pd.DataFrame,
    *,
    baseline_resolution: float,
    include_ami: bool,
) -> pd.DataFrame:
    labels = df[
        ["resolution", "window_id", "sequence_id", "cluster_id"]
    ].drop_duplicates()

    baseline_mask = np.isclose(labels["resolution"], baseline_resolution)
    baseline = labels.loc[baseline_mask, ["window_id", "sequence_id", "cluster_id"]]
    baseline = baseline.rename(columns={"cluster_id": "baseline_cluster_id"})  # type: ignore

    rows: list[dict[str, object]] = []
    baseline_by_window: dict[str, pd.DataFrame] = {
        str(window_id): group.drop(columns="window_id")
        for window_id, group in baseline.groupby("window_id", sort=False)
    }

    for (resolution, window_id), group in labels.groupby(
        ["resolution", "window_id"],
        sort=False,
    ):
        window_id_text = str(window_id)
        resolution_value = _as_float(resolution)
        base = baseline_by_window.get(window_id_text)
        if base is None or base.empty:
            ari = np.nan
            ami = np.nan if include_ami else None
            n_compared = 0
        else:
            merged = base.merge(
                group[["sequence_id", "cluster_id"]],
                on="sequence_id",
                how="inner",
            )
            n_compared = int(len(merged))
            if n_compared < 2:
                ari = np.nan
                ami = np.nan if include_ami else None
            elif np.isclose(resolution_value, baseline_resolution):
                ari = 1.0
                ami = 1.0 if include_ami else None
            else:
                baseline_codes = pd.factorize(
                    merged["baseline_cluster_id"],
                    sort=False,
                )[0]
                cluster_codes = pd.factorize(merged["cluster_id"], sort=False)[0]
                ari = float(
                    adjusted_rand_score(
                        baseline_codes,
                        cluster_codes,
                    )
                )
                if include_ami:
                    ami = float(
                        adjusted_mutual_info_score(
                            baseline_codes,
                            cluster_codes,
                        )
                    )
                else:
                    ami = None
        row: dict[str, object] = {
            "resolution": resolution_value,
            "window_id": window_id_text,
            "ari_vs_baseline_resolution": ari,
            "n_sequences_compared_to_baseline": n_compared,
        }
        if include_ami:
            row["ami_vs_baseline_resolution"] = ami
        rows.append(row)

    return pd.DataFrame(rows)


def _add_leiden_baseline_contrasts(
    window: pd.DataFrame,
    *,
    baseline_resolution: float,
) -> pd.DataFrame:
    metric_cols = [
        "clusters_per_1000_sequences",
        "p90_cluster_size",
        "max_cluster_size",
        "singleton_cluster_share",
        "singleton_sequence_share",
    ]
    baseline_cols = ["window_id", *metric_cols]

    baseline = window.loc[
        np.isclose(window["resolution"], baseline_resolution),
        baseline_cols,
    ].set_axis(
        ["window_id", *(f"baseline_{c}" for c in metric_cols)],
        axis="columns",
    )

    out = window.merge(baseline, on="window_id", how="left", validate="m:1")

    for col in metric_cols:
        base_col = f"baseline_{col}"
        out[f"delta_{col}_vs_baseline"] = out[col] - out[base_col]
        out[f"ratio_{col}_vs_baseline"] = out[col] / out[base_col].replace(0, np.nan)

    return out


def _build_leiden_summary(
    window: pd.DataFrame,
    *,
    include_ami: bool,
) -> pd.DataFrame:
    aggregations: dict[str, NamedAggregation] = {
        "n_windows": ("window_id", "nunique"),
        "median_ari_vs_baseline": ("ari_vs_baseline_resolution", "median"),
        "q25_ari_vs_baseline": ("ari_vs_baseline_resolution", _quantile(0.25)),
        "q75_ari_vs_baseline": ("ari_vs_baseline_resolution", _quantile(0.75)),
        "min_ari_vs_baseline": ("ari_vs_baseline_resolution", "min"),
    }
    if include_ami:
        aggregations.update(
            {
                "median_ami_vs_baseline": ("ami_vs_baseline_resolution", "median"),
                "q25_ami_vs_baseline": (
                    "ami_vs_baseline_resolution",
                    _quantile(0.25),
                ),
                "q75_ami_vs_baseline": (
                    "ami_vs_baseline_resolution",
                    _quantile(0.75),
                ),
                "min_ami_vs_baseline": ("ami_vs_baseline_resolution", "min"),
            }
        )
    aggregations.update(
        {
            "median_clusters_per_1000_sequences": (
                "clusters_per_1000_sequences",
                "median",
            ),
            "q25_clusters_per_1000_sequences": (
                "clusters_per_1000_sequences",
                _quantile(0.25),
            ),
            "q75_clusters_per_1000_sequences": (
                "clusters_per_1000_sequences",
                _quantile(0.75),
            ),
            "median_p90_cluster_size": ("p90_cluster_size", "median"),
            "q25_p90_cluster_size": ("p90_cluster_size", _quantile(0.25)),
            "q75_p90_cluster_size": ("p90_cluster_size", _quantile(0.75)),
            "median_max_cluster_size": ("max_cluster_size", "median"),
            "median_singleton_cluster_share": ("singleton_cluster_share", "median"),
            "q25_singleton_cluster_share": (
                "singleton_cluster_share",
                _quantile(0.25),
            ),
            "q75_singleton_cluster_share": (
                "singleton_cluster_share",
                _quantile(0.75),
            ),
            "median_singleton_sequence_share": (
                "singleton_sequence_share",
                "median",
            ),
            "median_duration_days": ("median_duration_days", "median"),
            "median_datazones": ("median_datazones", "median"),
        }
    )
    summary = (
        window.groupby("resolution", dropna=False)
        .agg(**aggregations)
        .reset_index()
        .sort_values("resolution")
    )
    return summary


def _metadata_from_pairwise_path(path: Path) -> dict[str, object]:
    parts = path.stem.split("_")
    if len(parts) < 3 or not parts[-1].isdigit():
        raise ValueError(
            "Pairwise filename must look like {window}_{lineage}_{count}: "
            f"{path.name}"
        )
    window_id = _normalise_window(parts[0])
    return {
        "window_id": window_id,
        "window_idx": int(window_id[1:]) if window_id[1:].isdigit() else np.nan,
        "pango_lineage": "_".join(parts[1:-1]),
        "pairwise_stem": path.stem,
        "nunique_sequences": int(parts[-1]),
    }


def _load_edge_manifest(edge_manifest: Path) -> pd.DataFrame:
    if not edge_manifest.exists():
        LOGGER.warning("Edge manifest not found: %s", edge_manifest)
        return pd.DataFrame()
    manifest = pd.read_parquet(edge_manifest)
    if "pairwise_stem" not in manifest.columns:
        raise ValueError(f"Edge manifest lacks pairwise_stem: {edge_manifest}")
    return manifest


def _select_pairwise_files(
    pairwise_dir: Path,
    *,
    windows: Sequence[str] | None,
    pango_lineages: Sequence[str] | None,
    max_windows: int | None,
    max_files: int | None,
) -> list[Path]:
    files = sorted(pairwise_dir.glob("*.parquet"), key=lambda path: path.stem)
    if not files:
        raise FileNotFoundError(f"No pairwise parquet files found in {pairwise_dir}")

    meta = pd.DataFrame(_metadata_from_pairwise_path(path) for path in files)

    if windows is not None:
        meta = meta.loc[meta["window_id"].isin(set(windows))]
    if pango_lineages:
        meta = meta.loc[meta["pango_lineage"].isin(set(pango_lineages))]
    if max_windows is not None:
        keep = sorted(meta["window_id"].dropna().unique())[:max_windows]
        meta = meta.loc[meta["window_id"].isin(keep)]

    selected_stems = set(meta["pairwise_stem"])
    selected = [path for path in files if path.stem in selected_stems]
    if max_files is not None:
        selected = selected[:max_files]
    if not selected:
        raise ValueError("No pairwise files remain after sensitivity filters.")
    return selected


def build_sparsification_sensitivity_tables(
    *,
    pairwise_dir: Path,
    edge_manifest: Path,
    thresholds: Sequence[float],
    baseline_threshold: float,
    windows: Sequence[str] | None,
    pango_lineages: Sequence[str] | None,
    max_windows: int | None,
    max_files: int | None,
    max_row_groups_per_file: int | None,
    max_rows_per_file: int | None,
    score_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return pairwise-file and threshold-level sparsification sensitivity tables."""
    thresholds = _sorted_unique_float([*thresholds, baseline_threshold])
    manifest = _load_edge_manifest(edge_manifest)
    files = _select_pairwise_files(
        pairwise_dir,
        windows=windows,
        pango_lineages=pango_lineages,
        max_windows=max_windows,
        max_files=max_files,
    )
    manifest_by_stem = (
        manifest.set_index("pairwise_stem").to_dict("index")
        if not manifest.empty
        else {}
    )

    LOGGER.info(
        "Scanning %s pairwise files across %s thresholds",
        f"{len(files):,}",
        len(thresholds),
    )
    rows: list[dict[str, object]] = []
    for idx, path in enumerate(files, start=1):
        if idx == 1 or idx % 100 == 0 or idx == len(files):
            LOGGER.info("Scanning pairwise file %s/%s: %s", idx, len(files), path.name)
        metadata = _metadata_from_pairwise_path(path)
        metadata.update(manifest_by_stem.get(path.stem, {}))  # type: ignore
        scanned = _scan_pairwise_file(
            path,
            thresholds=thresholds,
            score_column=score_column,
            max_row_groups=max_row_groups_per_file,
            max_rows=max_rows_per_file,
        )
        total_pairwise_edges = int(
            metadata.get("total_pairwise_edges", scanned["total_rows"])  # type: ignore
        )
        n_sequences = int(metadata["nunique_sequences"])  # type: ignore
        scan_rows = int(scanned["scan_rows"])  # type: ignore
        scan_fraction = (
            scan_rows / scanned["total_rows"] if scanned["total_rows"] else np.nan  # type: ignore
        )
        estimated = scan_rows != scanned["total_rows"]
        scan_weight_sum = float(scanned["scan_weight_sum"])
        total_weight_sum = (
            scan_weight_sum
            if not estimated
            else scan_weight_sum / scan_fraction
            if scan_fraction and not pd.isna(scan_fraction)
            else np.nan
        )

        for threshold, retained_scan_edges, retained_scan_weight_sum in zip(
            thresholds,
            scanned["retained_scan_edges"],  # type: ignore
            scanned["retained_scan_weight_sums"],  # type: ignore
            strict=True,
        ):
            retained_scan_weight_sum = float(retained_scan_weight_sum)
            retained_fraction = retained_scan_edges / scan_rows if scan_rows else np.nan
            retained_weight_fraction = (
                retained_scan_weight_sum / scan_weight_sum
                if scan_weight_sum
                else np.nan
            )
            retained_edges = (
                int(retained_scan_edges)
                if not estimated
                else int(round(retained_fraction * total_pairwise_edges))
            )
            retained_weight_sum = (
                retained_scan_weight_sum
                if not estimated
                else retained_weight_fraction * total_weight_sum
            )
            rows.append(
                {
                    "threshold": threshold,
                    "threshold_label": _threshold_label(threshold),
                    "window_idx": metadata.get("window_idx"),
                    "window_id": metadata.get("window_id"),
                    "pango_lineage": metadata.get("pango_lineage"),
                    "pairwise_stem": path.stem,
                    "nunique_sequences": n_sequences,
                    "total_pairwise_edges": total_pairwise_edges,
                    "total_weight_sum": total_weight_sum,
                    "retained_edges": retained_edges,
                    "retained_scan_edges": int(retained_scan_edges),
                    "retained_weight_sum": retained_weight_sum,
                    "retained_scan_weight_sum": retained_scan_weight_sum,
                    "retained_edge_fraction": retained_fraction,
                    "retained_weight_fraction": retained_weight_fraction,
                    "retained_mean_degree": (
                        2 * retained_edges / n_sequences if n_sequences else np.nan
                    ),
                    "scan_rows": scan_rows,
                    "scan_weight_sum": scan_weight_sum,
                    "scan_fraction": scan_fraction,
                    "estimated_from_partial_scan": estimated,
                }
            )

    detail = pd.DataFrame(rows)
    detail = _add_sparsification_baseline_contrasts(
        detail,
        baseline_threshold=baseline_threshold,
    )
    detail = detail.sort_values(
        ["threshold", "window_idx", "pango_lineage", "pairwise_stem"],
        kind="mergesort",
    ).reset_index(drop=True)
    summary = _build_sparsification_summary(detail)
    return detail, summary


def _scan_pairwise_file(
    path: Path,
    *,
    thresholds: Sequence[float],
    score_column: str,
    max_row_groups: int | None,
    max_rows: int | None,
) -> PairwiseScan:
    parquet_file = pq.ParquetFile(path)
    total_rows = int(parquet_file.metadata.num_rows)
    row_groups = range(parquet_file.num_row_groups)
    if max_row_groups is not None:
        row_groups = range(min(parquet_file.num_row_groups, max_row_groups))

    counts = np.zeros(len(thresholds), dtype=np.int64)
    weight_sums = np.zeros(len(thresholds), dtype=np.float64)
    scan_rows = 0
    scan_weight_sum = 0.0
    for row_group in row_groups:
        if max_rows is not None and scan_rows >= max_rows:
            break
        table = parquet_file.read_row_group(row_group, columns=[score_column])
        values = table.column(score_column).to_numpy(zero_copy_only=False).astype(float)
        if max_rows is not None and scan_rows + len(values) > max_rows:
            values = values[: max_rows - scan_rows]
        scan_rows += len(values)
        scan_weight_sum += float(np.nansum(values))
        for idx, threshold in enumerate(thresholds):
            retained = values > threshold
            counts[idx] += int(np.count_nonzero(retained))
            weight_sums[idx] += float(np.nansum(values[retained]))

    return {
        "total_rows": total_rows,
        "scan_rows": scan_rows,
        "scan_weight_sum": scan_weight_sum,
        "retained_scan_edges": counts,
        "retained_scan_weight_sums": weight_sums,
    }


def _add_sparsification_baseline_contrasts(
    detail: pd.DataFrame,
    *,
    baseline_threshold: float,
) -> pd.DataFrame:
    baseline = detail.loc[
        np.isclose(detail["threshold"].astype(float), baseline_threshold),
        [
            "pairwise_stem",
            "retained_edges",
            "retained_edge_fraction",
            "retained_weight_sum",
            "retained_weight_fraction",
            "retained_mean_degree",
        ],
    ].rename(
        columns={
            "retained_edges": "baseline_retained_edges",
            "retained_edge_fraction": "baseline_retained_edge_fraction",
            "retained_weight_sum": "baseline_retained_weight_sum",
            "retained_weight_fraction": "baseline_retained_weight_fraction",
            "retained_mean_degree": "baseline_retained_mean_degree",
        }
    )  # type: ignore
    out = detail.merge(baseline, on="pairwise_stem", how="left")
    out["retained_edge_ratio_vs_baseline"] = out["retained_edges"] / out[
        "baseline_retained_edges"
    ].replace(0, np.nan)
    out["delta_retained_edge_fraction_vs_baseline"] = (
        out["retained_edge_fraction"] - out["baseline_retained_edge_fraction"]
    )
    out["retained_weight_ratio_vs_baseline"] = out["retained_weight_sum"] / out[
        "baseline_retained_weight_sum"
    ].replace(0, np.nan)
    out["delta_retained_weight_fraction_vs_baseline"] = (
        out["retained_weight_fraction"] - out["baseline_retained_weight_fraction"]
    )
    out["retained_mean_degree_ratio_vs_baseline"] = out["retained_mean_degree"] / out[
        "baseline_retained_mean_degree"
    ].replace(0, np.nan)
    return out


def _build_sparsification_summary(detail: pd.DataFrame) -> pd.DataFrame:
    summary = (
        detail.groupby(["threshold", "threshold_label"], dropna=False)
        .agg(
            n_pairwise_groups=("pairwise_stem", "nunique"),
            n_windows=("window_id", "nunique"),
            n_lineages=("pango_lineage", "nunique"),
            total_pairwise_edges=("total_pairwise_edges", "sum"),
            total_retained_edges=("retained_edges", "sum"),
            total_weight_sum=("total_weight_sum", "sum"),
            total_retained_weight_sum=("retained_weight_sum", "sum"),
            median_retained_edge_fraction=("retained_edge_fraction", "median"),
            q25_retained_edge_fraction=("retained_edge_fraction", _quantile(0.25)),
            q75_retained_edge_fraction=("retained_edge_fraction", _quantile(0.75)),
            q10_retained_edge_fraction=("retained_edge_fraction", _quantile(0.10)),
            q90_retained_edge_fraction=("retained_edge_fraction", _quantile(0.90)),
            median_retained_weight_fraction=("retained_weight_fraction", "median"),
            q25_retained_weight_fraction=("retained_weight_fraction", _quantile(0.25)),
            q75_retained_weight_fraction=("retained_weight_fraction", _quantile(0.75)),
            q10_retained_weight_fraction=("retained_weight_fraction", _quantile(0.10)),
            q90_retained_weight_fraction=("retained_weight_fraction", _quantile(0.90)),
            median_retained_mean_degree=("retained_mean_degree", "median"),
            q25_retained_mean_degree=("retained_mean_degree", _quantile(0.25)),
            q75_retained_mean_degree=("retained_mean_degree", _quantile(0.75)),
            median_retained_edge_ratio_vs_baseline=(
                "retained_edge_ratio_vs_baseline",
                "median",
            ),
            q25_retained_edge_ratio_vs_baseline=(
                "retained_edge_ratio_vs_baseline",
                _quantile(0.25),
            ),
            q75_retained_edge_ratio_vs_baseline=(
                "retained_edge_ratio_vs_baseline",
                _quantile(0.75),
            ),
            median_retained_weight_ratio_vs_baseline=(
                "retained_weight_ratio_vs_baseline",
                "median",
            ),
            q25_retained_weight_ratio_vs_baseline=(
                "retained_weight_ratio_vs_baseline",
                _quantile(0.25),
            ),
            q75_retained_weight_ratio_vs_baseline=(
                "retained_weight_ratio_vs_baseline",
                _quantile(0.75),
            ),
            min_scan_fraction=("scan_fraction", "min"),
            estimated_from_partial_scan=("estimated_from_partial_scan", "any"),
        )
        .reset_index()
        .sort_values("threshold")
    )
    summary["pooled_retained_edge_fraction"] = summary[
        "total_retained_edges"
    ] / summary["total_pairwise_edges"].replace(0, np.nan)
    summary["pooled_retained_weight_fraction"] = summary[
        "total_retained_weight_sum"
    ] / summary["total_weight_sum"].replace(0, np.nan)
    return summary


def _write_table(df: pd.DataFrame, name: str, *, table_dir: Path) -> None:
    LOGGER.info("Writing %s (%s rows)", name, f"{len(df):,}")
    write_table(df, name, table_dir=table_dir, formats=("parquet", "csv"))


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=args.log_level, format="%(levelname)s: %(message)s")
    ensure_results_dirs()
    args.table_dir.mkdir(parents=True, exist_ok=True)

    windows = _normalise_windows(args.windows)
    qc = _normalise_qc(args.qc)
    thresholds = _sorted_unique_float(args.thresholds)

    if args.only in {"all", "leiden"}:
        detail, summary = build_leiden_sensitivity_tables(
            analysis_dataset=args.analysis_dataset,
            baseline_resolution=args.baseline_resolution,
            include_ami=args.include_ami,
            qc=qc,
            windows=windows,
            max_windows=args.max_windows,
        )
        _write_table(
            detail, "leiden_resolution_window_sensitivity", table_dir=args.table_dir
        )
        _write_table(
            summary, "leiden_resolution_sensitivity_summary", table_dir=args.table_dir
        )

    if args.only in {"all", "sparsification"}:
        detail, summary = build_sparsification_sensitivity_tables(
            pairwise_dir=args.pairwise_dir,
            edge_manifest=args.edge_manifest,
            thresholds=thresholds,
            baseline_threshold=args.baseline_threshold,
            windows=windows,
            pango_lineages=args.pango_lineages,
            max_windows=args.max_windows,
            max_files=args.max_files,
            max_row_groups_per_file=args.max_row_groups_per_file,
            max_rows_per_file=args.max_rows_per_file,
            score_column=args.score_column,
        )
        _write_table(
            detail, "sparsification_threshold_sensitivity", table_dir=args.table_dir
        )
        _write_table(
            summary,
            "sparsification_threshold_sensitivity_summary",
            table_dir=args.table_dir,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
