"""Build compatibility-network mixing matrices and assortativity summaries.

Pairwise compatibility edge scans can be large. For development, pass a small
window set, for example:

    python -m chapter_analyses.genomic_networks.build_mixing --windows W080 W081

For the full Chapter 4 run:

    python -m chapter_analyses.genomic_networks.build_mixing --all-windows --workers 4 \
        --include-giants --giant-workers 1

Use --dry-run first to inspect the small-file and giant-file schedule.
Giant files are skipped by default.
"""

from __future__ import annotations

import argparse
import logging
import math
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, cast

import numpy as np
import pandas as pd

from .lib.config import (
    DEFAULT_MIXING_ATTRIBUTES,
    INTERMEDIATE_DIR,
    PROJECT_ROOT,
    SPARSIFICATION_THRESHOLD,
)
from .lib.io import (
    ensure_results_dirs,
    load_sequence_data,
    write_table,
)
from .lib.mixing import (
    build_degree_assortativity_for_edge_table,
    build_mixing_for_edge_table,
    specs_by_name,
)
from utils import load_pairwise_edges

LOGGER = logging.getLogger(__name__)

PAIRWISE_DATASET_DIR = PROJECT_ROOT / "data/processed/pairwise_distances_dataset"
EDGE_MANIFEST_PATH = (
    PROJECT_ROOT / "data/processed/sparsified_edge_counts_by_window_lineage.parquet"
)
INTERMEDIATE_TABLE_DIRS = {
    "matrix": INTERMEDIATE_DIR / "mixing_matrix",
    "summary": INTERMEDIATE_DIR / "comp_assortativity",
    "topology": INTERMEDIATE_DIR / "deg_assortativity",
}
FINAL_TABLES = {
    "matrix": "compatibility_mixing_matrix",
    "summary": "compatibility_assortativity",
    "topology": "compatibility_degree_assortativity",
}
PAIRWISE_COLUMNS = [
    "id1",
    "id2",
    "epilink_compatibility",
]
JACKKNIFE_LEAVE_ONE_NODE_LIMIT = 1_000
JACKKNIFE_MIN_BLOCKS = 50
JACKKNIFE_TARGET_BLOCK_SIZE = 1_000
JACKKNIFE_MIN_FINITE_REPLICATES = 5


@dataclass(frozen=True)
class PairwiseMetadata:
    window_id: str
    pango_lineage: str
    sequence_count: int


@dataclass(frozen=True)
class ScheduledTask:
    args: tuple
    edge_cost: int | None
    edge_cost_source: str

    @property
    def pairwise_path(self) -> Path:
        return self.args[0]

    @property
    def stem(self) -> str:
        return self.pairwise_path.stem


def _format_task_progress(processed: int, total: int) -> str:
    if total <= 0:
        return "0/0 tasks (0.0%)"
    return f"{processed:,}/{total:,} tasks ({processed / total:.1%})"


def _format_status_counts(statuses: dict[str, int]) -> str:
    return (
        f"{statuses.get('ok', 0):,} ok, "
        f"{statuses.get('skipped', 0):,} skipped, "
        f"{statuses.get('no_edges', 0):,} no_edges"
    )


def _normalise_window(value: str) -> str:
    value = str(value).strip()
    upper = value.upper()
    if upper.startswith("W") and upper[1:].isdigit():
        return f"W{int(upper[1:]):03d}"
    if upper.isdigit():
        return f"W{int(upper):03d}"
    return value


def _pairwise_metadata_from_path(path: Path) -> PairwiseMetadata:
    # Stems are {window}_{lineage}_{count}; Pango lineages use dots, not
    # underscores, but join the middle fields so the parser is robust.
    parts = path.stem.split("_")
    if len(parts) < 3 or not parts[-1].isdigit():
        raise ValueError(
            f"Pairwise filename must look like {{window}}_{{lineage}}_{{count}}: "
            f"{path.name}"
        )
    lineage = "_".join(parts[1:-1])
    if not lineage:
        raise ValueError(
            f"Could not parse Pango lineage from pairwise filename: {path.name}"
        )
    return PairwiseMetadata(
        window_id=_normalise_window(parts[0]),
        pango_lineage=lineage,
        sequence_count=int(parts[-1]),
    )


def _window_from_pairwise_path(path: Path) -> str:
    return _pairwise_metadata_from_path(path).window_id


def _select_pairwise_files(
    pairwise_dir: Path,
    *,
    windows: Sequence[str],
    max_windows: int | None,
) -> list[Path]:
    files = sorted(pairwise_dir.glob("*.parquet"), key=lambda path: path.stem)
    if windows:
        wanted = set(windows)
        files = [path for path in files if _window_from_pairwise_path(path) in wanted]

    if max_windows is not None:
        retained_windows = sorted({_window_from_pairwise_path(path) for path in files})
        retained_windows = retained_windows[:max_windows]
        wanted = set(retained_windows)
        files = [path for path in files if _window_from_pairwise_path(path) in wanted]

    return files


def _output_paths_for_stem(stem: str) -> dict[str, Path]:
    return {
        kind: table_dir / f"{stem}.parquet"
        for kind, table_dir in INTERMEDIATE_TABLE_DIRS.items()
    }


def _all_outputs_exist(paths: dict[str, Path]) -> bool:
    return all(path.exists() for path in paths.values())


def _ensure_intermediate_dirs() -> None:
    for table_dir in INTERMEDIATE_TABLE_DIRS.values():
        table_dir.mkdir(parents=True, exist_ok=True)


def _read_pairwise_edges(
    pairwise_path: Path,
    *,
    compatibility_threshold: float | None,
) -> pd.DataFrame:
    return load_pairwise_edges(
        PAIRWISE_COLUMNS,
        compatibility_threshold=compatibility_threshold,
        pairwise_dataset=pairwise_path,
    )


def _with_pairwise_metadata(
    df: pd.DataFrame,
    metadata: PairwiseMetadata,
    stem: str,
    *,
    kind: str,
) -> pd.DataFrame:
    out = df.copy()
    out["window_id"] = metadata.window_id
    out["pango_lineage"] = metadata.pango_lineage
    out["pairwise_stem"] = stem

    metadata_cols = ["window_id", "pango_lineage", "pairwise_stem"]
    non_metadata = [col for col in out.columns if col not in metadata_cols]
    if kind == "topology":
        ordered = ["window_id", "pango_lineage", *non_metadata, "pairwise_stem"]
    elif kind == "summary" and "attribute_label" in non_metadata:
        insert_at = non_metadata.index("attribute_label") + 1
        ordered = [
            *non_metadata[:insert_at],
            "window_id",
            "pango_lineage",
            *non_metadata[insert_at:],
            "pairwise_stem",
        ]
    else:
        ordered = [*non_metadata, "window_id", "pango_lineage", "pairwise_stem"]
    return out[ordered]


def _as_finite_weight_array(
    edges: pd.DataFrame,
    *,
    weight_col: str,
) -> np.ndarray:
    return (
        pd.to_numeric(edges[weight_col], errors="coerce")
        .fillna(0.0)
        .to_numpy(dtype=float)
    )


def _edge_arrays_for_jackknife(
    edges: pd.DataFrame,
    *,
    source_col: str,
    target_col: str,
    weight_col: str,
) -> dict[str, np.ndarray | pd.Index]:
    work = edges[[source_col, target_col]].copy()
    work["_edge_weight"] = _as_finite_weight_array(edges, weight_col=weight_col)
    work = work.dropna(subset=[source_col, target_col])
    work = work.loc[work["_edge_weight"].gt(0)]

    if work.empty:
        return {
            "vertex_names": pd.Index([], dtype="object"),
            "source_vertices": np.array([], dtype=int),
            "target_vertices": np.array([], dtype=int),
            "edge_weights": np.array([], dtype=float),
        }

    endpoints = work[[source_col, target_col]].to_numpy().ravel()
    vertex_codes, vertex_names = pd.factorize(endpoints, sort=False)
    vertex_names = pd.Index(vertex_names)
    return {
        "vertex_names": vertex_names,
        "source_vertices": vertex_codes[0::2].astype(int),
        "target_vertices": vertex_codes[1::2].astype(int),
        "edge_weights": work["_edge_weight"].to_numpy(dtype=float),
    }


def _jackknife_block_count(n_vertices: int, requested_blocks: int) -> int:
    if requested_blocks <= 0 or n_vertices < 2:
        return 0
    if n_vertices <= JACKKNIFE_LEAVE_ONE_NODE_LIMIT:
        return n_vertices

    adaptive_blocks = max(
        JACKKNIFE_MIN_BLOCKS,
        math.ceil(math.sqrt(n_vertices)),
        math.ceil(n_vertices / JACKKNIFE_TARGET_BLOCK_SIZE),
    )
    return min(requested_blocks, adaptive_blocks, n_vertices)


def _assign_jackknife_blocks(
    vertex_names: pd.Index,
    *,
    requested_blocks: int,
    seed: int,
) -> np.ndarray:
    n_vertices = len(vertex_names)
    n_blocks = _jackknife_block_count(n_vertices, requested_blocks)
    if n_blocks == 0:
        return np.array([], dtype=int)

    keys = pd.Series(vertex_names, dtype="object").astype("string") + f"|{seed}"
    hashes = pd.util.hash_pandas_object(keys, index=False).to_numpy(dtype=np.uint64)
    order = np.argsort(hashes, kind="mergesort")
    blocks = np.empty(n_vertices, dtype=int)
    blocks[order] = np.arange(n_vertices, dtype=int) % n_blocks
    return blocks


def _assortativity_from_dense(mixing_matrix: np.ndarray) -> float:
    total = float(mixing_matrix.sum())
    if total <= 0:
        return np.nan

    e = mixing_matrix / total
    observed = float(np.trace(e))
    expected = float(np.dot(e.sum(axis=1), e.sum(axis=0)))
    denominator = 1.0 - expected
    r = np.nan if np.isclose(denominator, 0.0) else (observed - expected) / denominator
    return float(r) if not pd.isna(r) else np.nan


def _jackknife_estimates_for_categories(
    source_categories: np.ndarray,
    target_categories: np.ndarray,
    source_blocks: np.ndarray,
    target_blocks: np.ndarray,
    edge_weights: np.ndarray,
    *,
    n_categories: int,
    n_blocks: int,
) -> np.ndarray:
    n_cells = n_categories * n_categories
    flat = source_categories * n_categories + target_categories
    reverse_flat = target_categories * n_categories + source_categories

    full = np.bincount(flat, weights=edge_weights, minlength=n_cells)
    full += np.bincount(reverse_flat, weights=edge_weights, minlength=n_cells)
    full = full.reshape(n_categories, n_categories)

    removal = np.bincount(
        source_blocks * n_cells + flat,
        weights=edge_weights,
        minlength=n_blocks * n_cells,
    )
    removal += np.bincount(
        source_blocks * n_cells + reverse_flat,
        weights=edge_weights,
        minlength=n_blocks * n_cells,
    )

    different_blocks = source_blocks != target_blocks
    if different_blocks.any():
        removal += np.bincount(
            target_blocks[different_blocks] * n_cells + flat[different_blocks],
            weights=edge_weights[different_blocks],
            minlength=n_blocks * n_cells,
        )
        removal += np.bincount(
            target_blocks[different_blocks] * n_cells + reverse_flat[different_blocks],
            weights=edge_weights[different_blocks],
            minlength=n_blocks * n_cells,
        )

    removal = removal.reshape(n_blocks, n_categories, n_categories)
    return np.array(
        [
            _assortativity_from_dense(full - removal[block_idx])
            for block_idx in range(n_blocks)
        ],
        dtype=float,
    )


def _jackknife_uncertainty_for_attribute(
    labels: pd.Series,
    source_vertices: np.ndarray,
    target_vertices: np.ndarray,
    edge_weights: np.ndarray,
    vertex_names: pd.Index,
    *,
    missing_label: str | None,
    requested_blocks: int,
    seed: int,
) -> dict[str, float | int | str]:
    if requested_blocks <= 0:
        return _empty_jackknife_uncertainty(0)

    if missing_label is not None:
        labels = labels.fillna(missing_label)
    valid_vertices = ~labels.isna()
    labels = labels.astype("string")

    edge_mask = (
        valid_vertices.to_numpy()[source_vertices]
        & valid_vertices.to_numpy()[target_vertices]
        & (edge_weights > 0)
    )
    if not edge_mask.any():
        return _empty_jackknife_uncertainty(0)

    source_vertices = source_vertices[edge_mask]
    target_vertices = target_vertices[edge_mask]
    edge_weights = edge_weights[edge_mask]

    used_vertices = np.unique(np.concatenate([source_vertices, target_vertices]))
    remap = np.full(len(labels), -1, dtype=int)
    remap[used_vertices] = np.arange(used_vertices.size)
    source_vertices = remap[source_vertices]
    target_vertices = remap[target_vertices]

    used_vertex_names = pd.Index(vertex_names.take(used_vertices))
    vertex_blocks = _assign_jackknife_blocks(
        used_vertex_names,
        requested_blocks=requested_blocks,
        seed=seed,
    )
    if vertex_blocks.size == 0:
        return _empty_jackknife_uncertainty(0, n_vertices=used_vertices.size)

    used_labels = labels.iloc[used_vertices].astype(str).to_numpy()
    unique_labels = np.array(sorted(pd.unique(used_labels).tolist()), dtype=object)
    label_to_index = {label: idx for idx, label in enumerate(unique_labels)}
    category_indices = np.fromiter(
        (label_to_index[label] for label in used_labels),
        dtype=int,
        count=len(used_labels),
    )

    source_categories = category_indices[source_vertices]
    target_categories = category_indices[target_vertices]
    source_blocks = vertex_blocks[source_vertices]
    target_blocks = vertex_blocks[target_vertices]
    n_blocks = int(vertex_blocks.max()) + 1

    estimates = _jackknife_estimates_for_categories(
        source_categories,
        target_categories,
        source_blocks,
        target_blocks,
        edge_weights,
        n_categories=len(unique_labels),
        n_blocks=n_blocks,
    )
    finite = estimates[np.isfinite(estimates)]
    if finite.size < JACKKNIFE_MIN_FINITE_REPLICATES:
        return _empty_jackknife_uncertainty(
            n_blocks,
            n_vertices=used_vertices.size,
            replicates=finite.size,
        )

    estimate_mean = float(finite.mean())
    se = float(
        np.sqrt(
            (finite.size - 1) / finite.size * np.square(finite - estimate_mean).sum()
        )
    )
    return {
        "uncertainty_method": _jackknife_uncertainty_method(
            n_blocks,
            n_vertices=used_vertices.size,
        ),
        "jackknife_vertices_used": int(used_vertices.size),
        "jackknife_blocks_used": int(n_blocks),
        "jackknife_replicates": int(finite.size),
        "jackknife_assortativity_mean": estimate_mean,
        "assortativity_se": se,
        "assortativity_ci_low": np.nan,
        "assortativity_ci_high": np.nan,
    }


def _jackknife_uncertainty_method(n_blocks: int, *, n_vertices: int) -> str:
    if n_vertices >= 2 and n_blocks == n_vertices:
        return "leave_one_node_jackknife"
    return "node_block_jackknife"


def _empty_jackknife_uncertainty(
    n_blocks: int,
    *,
    n_vertices: int = 0,
    replicates: int = 0,
) -> dict[str, float | int | str]:
    return {
        "uncertainty_method": _jackknife_uncertainty_method(
            n_blocks,
            n_vertices=n_vertices,
        ),
        "jackknife_vertices_used": int(n_vertices),
        "jackknife_blocks_used": int(n_blocks),
        "jackknife_replicates": int(replicates),
        "jackknife_assortativity_mean": np.nan,
        "assortativity_se": np.nan,
        "assortativity_ci_low": np.nan,
        "assortativity_ci_high": np.nan,
    }


def _add_jackknife_uncertainty(
    summary: pd.DataFrame,
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    attributes,
    *,
    node_id_col: str,
    source_col: str,
    target_col: str,
    weight_col: str,
    missing_label: str | None,
    requested_blocks: int,
    seed: int,
) -> pd.DataFrame:
    if summary.empty or requested_blocks <= 0:
        return summary

    edge_arrays = _edge_arrays_for_jackknife(
        edges,
        source_col=source_col,
        target_col=target_col,
        weight_col=weight_col,
    )
    vertex_names = cast(pd.Index, edge_arrays["vertex_names"])
    if vertex_names.empty:
        return summary

    source_vertices = cast(np.ndarray, edge_arrays["source_vertices"])
    target_vertices = cast(np.ndarray, edge_arrays["target_vertices"])
    edge_weights = cast(np.ndarray, edge_arrays["edge_weights"])

    lookup = node_attribute_lookup(
        nodes,
        node_id_col=node_id_col,
        attributes=attributes,
    )
    attr_cols = [spec.column for spec in attributes if spec.column in lookup.columns]
    attr_lookup = (
        lookup.set_index(node_id_col)[attr_cols] if attr_cols else pd.DataFrame()
    )

    uncertainty_rows = []
    for spec in attributes:
        if spec.column not in attr_lookup.columns:
            continue
        labels = pd.Series(vertex_names, dtype="object").map(attr_lookup[spec.column])
        uncertainty = _jackknife_uncertainty_for_attribute(
            labels,
            source_vertices,
            target_vertices,
            edge_weights,
            vertex_names,
            missing_label=missing_label,
            requested_blocks=requested_blocks,
            seed=seed,
        )
        uncertainty_rows.append({"attribute": spec.name, **uncertainty})

    if not uncertainty_rows:
        return summary

    out = summary.merge(pd.DataFrame(uncertainty_rows), on="attribute", how="left")
    finite_se = out["assortativity_se"].notna()
    out.loc[finite_se, "assortativity_ci_low"] = (
        out.loc[finite_se, "assortativity"]
        - out.loc[finite_se, "assortativity_se"] * 1.96
    )
    out.loc[finite_se, "assortativity_ci_high"] = (
        out.loc[finite_se, "assortativity"]
        + out.loc[finite_se, "assortativity_se"] * 1.96
    )
    return out


def _write_parquet_atomic(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    tmp_file = tempfile.NamedTemporaryFile(
        dir=path.parent,
        prefix=f".{path.stem}.",
        suffix=f".tmp{path.suffix}",
        delete=False,
    )
    try:
        tmp_path = Path(tmp_file.name)
        tmp_file.close()
        df.to_parquet(tmp_path, index=False)
        os.replace(tmp_path, path)
    finally:
        tmp_file.close()
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink()


def _write_missing_outputs(
    outputs: dict[str, pd.DataFrame],
    paths: dict[str, Path],
    *,
    force: bool,
) -> None:
    for kind, df in outputs.items():
        path = paths[kind]
        if path.exists() and not force:
            continue
        _write_parquet_atomic(df, path)


def _process_pairwise_file(
    pairwise_path: Path,
    nodes: pd.DataFrame,
    attributes,
    compatibility_threshold: float,
    missing_label: str | None,
    jackknife_blocks: int,
    jackknife_seed: int,
    force: bool,
) -> tuple[str, str]:
    """Build and write the mixing outputs for a single pairwise parquet file.

    Runs inside a worker process, so it must only use picklable arguments and
    module-level imports. Returns the pairwise file stem and a status string so
    the parent process can do all the logging with the configured handlers.
    """
    stem = pairwise_path.stem
    metadata = _pairwise_metadata_from_path(pairwise_path)
    output_paths = _output_paths_for_stem(stem)
    if not force and _all_outputs_exist(output_paths):
        return stem, "skipped"

    edges = _read_pairwise_edges(
        pairwise_path,
        compatibility_threshold=compatibility_threshold,
    )
    if edges.empty:
        _write_missing_outputs(
            {
                "matrix": _with_pairwise_metadata(
                    pd.DataFrame(), metadata, stem, kind="matrix"
                ),
                "summary": _with_pairwise_metadata(
                    pd.DataFrame(), metadata, stem, kind="summary"
                ),
                "topology": _with_pairwise_metadata(
                    pd.DataFrame(), metadata, stem, kind="topology"
                ),
            },
            output_paths,
            force=force,
        )
        return stem, "no_edges"

    topology = build_degree_assortativity_for_edge_table(
        edges,
        source_col="id1",
        target_col="id2",
        weight_col="epilink_compatibility",
        group_cols=None,
    )

    matrix, summary = build_mixing_for_edge_table(
        edges,
        nodes,
        attributes=attributes,
        node_id_col="sequence_id",
        source_col="id1",
        target_col="id2",
        weight_col="epilink_compatibility",
        group_cols=None,
        symmetric=True,
        missing_label=missing_label,
    )
    summary = _add_jackknife_uncertainty(
        summary,
        edges,
        nodes,
        attributes,
        node_id_col="sequence_id",
        source_col="id1",
        target_col="id2",
        weight_col="epilink_compatibility",
        missing_label=missing_label,
        requested_blocks=jackknife_blocks,
        seed=jackknife_seed,
    )
    _write_missing_outputs(
        {
            "matrix": _with_pairwise_metadata(matrix, metadata, stem, kind="matrix"),
            "summary": _with_pairwise_metadata(summary, metadata, stem, kind="summary"),
            "topology": _with_pairwise_metadata(
                topology, metadata, stem, kind="topology"
            ),
        },
        output_paths,
        force=force,
    )
    return stem, "ok"


def _sort_output_table(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    preferred = [
        "window_id",
        "pango_lineage",
        "pairwise_stem",
        "attribute",
        "source_category",
        "target_category",
    ]
    sort_cols = [col for col in preferred if col in df.columns]
    if not sort_cols:
        return df.reset_index(drop=True)
    return df.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)


def _concat_intermediate_table(kind: str, stems: Sequence[str]) -> pd.DataFrame:
    table_dir = INTERMEDIATE_TABLE_DIRS[kind]
    frames: list[pd.DataFrame] = []
    for stem in stems:
        path = table_dir / f"{stem}.parquet"
        if not path.exists():
            continue
        df = pd.read_parquet(path)
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return _sort_output_table(pd.concat(frames, ignore_index=True, sort=False))


def _load_edge_costs(edge_manifest: Path) -> dict[str, int]:
    if not edge_manifest.exists():
        LOGGER.warning(
            "Edge manifest not found at %s; falling back to parquet file sizes",
            edge_manifest,
        )
        return {}
    try:
        manifest = pd.read_parquet(
            edge_manifest,
            columns=["pairwise_stem", "sparse_edges"],
        )
    except ValueError as exc:
        raise SystemExit(
            "Edge manifest must contain pairwise_stem and sparse_edges: "
            f"{edge_manifest}"
        ) from exc

    costs: dict[str, int] = {}
    for row in manifest.itertuples(index=False):
        sparse_edges = row.sparse_edges
        if pd.isna(sparse_edges):
            continue
        costs[str(row.pairwise_stem)] = int(cast(int, sparse_edges))
    return costs


def _edge_cost_for_pairwise_path(
    pairwise_path: Path,
    manifest_costs: dict[str, int],
) -> tuple[int | None, str]:
    manifest_cost = manifest_costs.get(pairwise_path.stem)
    if manifest_cost is not None:
        return manifest_cost, "manifest"

    # Missing manifest rows still need a deterministic cost; file size is a
    # conservative local proxy, and stat failures are treated as giant later.
    try:
        return pairwise_path.stat().st_size, "file_size"
    except OSError:
        return None, "unknown"


def _is_giant_task(task: ScheduledTask, threshold: int) -> bool:
    return task.edge_cost is None or task.edge_cost >= threshold


def _descending_edge_cost(task: ScheduledTask) -> float:
    return float("inf") if task.edge_cost is None else float(task.edge_cost)


def _format_edge_cost(edge_cost: int | None) -> str:
    return "unknown" if edge_cost is None else f"{edge_cost:,}"


def _phase_edge_summary(tasks: Sequence[ScheduledTask]) -> str:
    known_total = sum(task.edge_cost for task in tasks if task.edge_cost is not None)
    unknown = sum(1 for task in tasks if task.edge_cost is None)
    suffix = f", {unknown:,} unknown" if unknown else ""
    return f"{known_total:,} summed edges/cost{suffix}"


def _print_dry_run_schedule(
    small_tasks: Sequence[ScheduledTask],
    giant_tasks: Sequence[ScheduledTask],
    *,
    workers: int,
    giant_workers: int,
    giant_threshold: int,
    include_giants: bool,
) -> None:
    total = len(small_tasks) + len(giant_tasks)
    processing_total = len(small_tasks) + (len(giant_tasks) if include_giants else 0)
    print("Dry run: no pairwise files will be processed.")
    print(f"Total selected files: {total:,}")
    print(f"Files that would be processed: {processing_total:,}")
    print(
        "Small phase: "
        f"{len(small_tasks):,} files, {_phase_edge_summary(small_tasks)}, "
        f"workers={workers:,}"
    )
    if include_giants:
        print(
            "Giant phase: "
            f"{len(giant_tasks):,} files, {_phase_edge_summary(giant_tasks)}, "
            f"workers={giant_workers:,}, threshold={giant_threshold:,}"
        )
    else:
        print(
            "Giant phase: "
            f"{len(giant_tasks):,} files skipped by default, "
            f"{_phase_edge_summary(giant_tasks)}, threshold={giant_threshold:,}; "
            "pass --include-giants to process them"
        )
    print("Giants in processing order:")
    if not giant_tasks:
        print("  (none)")
        return
    for task in giant_tasks:
        print(
            "  "
            f"{task.stem}\t{_format_edge_cost(task.edge_cost)}\t"
            f"{task.edge_cost_source}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--windows",
        nargs="*",
        help="Specific windows to process, e.g. W080 W081 or 80 81.",
    )
    parser.add_argument(
        "--all-windows",
        action="store_true",
        help="Process every available retained analysis window.",
    )
    parser.add_argument(
        "--max-windows",
        type=int,
        default=None,
        help="Optional cap for development runs after window selection.",
    )
    parser.add_argument(
        "--attributes",
        nargs="*",
        default=None,
        help=(
            "Attribute names to process. Defaults to all: "
            + ", ".join(spec.name for spec in DEFAULT_MIXING_ATTRIBUTES)
        ),
    )
    parser.add_argument(
        "--compatibility-threshold",
        type=float,
        default=SPARSIFICATION_THRESHOLD,
        help="Retain pairwise edges with epilink_compatibility above this threshold.",
    )
    parser.add_argument(
        "--pairwise-dir",
        type=Path,
        default=PAIRWISE_DATASET_DIR,
        help=(
            "Directory of per-window-lineage pairwise parquet files. "
            f"Default: {PAIRWISE_DATASET_DIR}"
        ),
    )
    parser.add_argument(
        "--missing-label",
        default=None,
        help="Optional label for missing node attributes. Default drops missing pairs.",
    )
    parser.add_argument(
        "--jackknife-blocks",
        type=int,
        default=1_000,
        help=(
            "Maximum deterministic node blocks for assortativity uncertainty. "
            "Attributes with up to 1,000 contributing vertices use standard "
            "leave-one-node jackknife; larger attributes use adaptive balanced "
            "blocks capped by this value. Use 0 to skip uncertainty columns. "
            "Default: 1,000."
        ),
    )
    parser.add_argument(
        "--jackknife-seed",
        type=int,
        default=42,
        help="Seed mixed into deterministic jackknife assignment. Default: 42.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help=(
            "Worker processes for the small-file phase. Default use 1. "
            "Use 1 to run serially (easier debugging/profiling)."
        ),
    )
    parser.add_argument(
        "--giant-threshold",
        type=int,
        default=50_000_000,
        help=(
            "Files with cost >= this many sparse edges are processed in the "
            "giant phase. Manifest misses fall back to file size; unknown costs "
            "are treated as giant. Default: 50,000,000."
        ),
    )
    parser.add_argument(
        "--giant-workers",
        type=int,
        default=1,
        help=(
            "Worker processes for the giant-file phase when --include-giants "
            "is passed. Default: 1."
        ),
    )
    parser.add_argument(
        "--include-giants",
        action="store_true",
        help=(
            "Process files at or above --giant-threshold. By default these "
            "memory-heavy files are skipped."
        ),
    )
    parser.add_argument(
        "--edge-manifest",
        type=Path,
        default=EDGE_MANIFEST_PATH,
        help=(
            "Parquet with pairwise_stem and sparse_edges columns used for scheduling. "
            f"Default: {EDGE_MANIFEST_PATH}"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the edge-count-aware schedule and exit without processing files.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute per-file intermediate outputs even when they already exist.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help=(
            "Log aggregate INFO progress every N completed pairwise files. "
            "Per-file messages are DEBUG. Default: 100."
        ),
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

    if not args.all_windows and not args.windows:
        raise SystemExit("Specify --windows or --all-windows.")
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1.")
    if args.giant_workers < 1:
        raise SystemExit("--giant-workers must be at least 1.")
    if args.giant_threshold < 0:
        raise SystemExit("--giant-threshold must be non-negative.")
    if args.progress_every < 1:
        raise SystemExit("--progress-every must be at least 1.")
    if args.jackknife_blocks < 0:
        raise SystemExit("--jackknife-blocks must be non-negative.")
    if 0 < args.jackknife_blocks < JACKKNIFE_MIN_FINITE_REPLICATES:
        raise SystemExit(
            "--jackknife-blocks must be 0 or at least "
            f"{JACKKNIFE_MIN_FINITE_REPLICATES}."
        )
    pairwise_dir = args.pairwise_dir
    if not pairwise_dir.exists():
        raise SystemExit(f"Pairwise dataset directory not found: {pairwise_dir}")
    if args.jackknife_blocks > 0:
        LOGGER.info(
            "Computing jackknife uncertainty with leave-one-node up to %s "
            "vertices and adaptive blocks capped at %s",
            f"{JACKKNIFE_LEAVE_ONE_NODE_LIMIT:,}",
            f"{args.jackknife_blocks:,}",
        )

    attributes = specs_by_name(args.attributes)
    attr_cols = [spec.column for spec in attributes]
    columns = ["window_id", "window_idx", "sequence_id", *attr_cols]
    sequence_df = load_sequence_data(columns=columns)

    if args.all_windows:
        windows = sorted(sequence_df["window_id"].dropna().unique())
    else:
        windows = [_normalise_window(window) for window in args.windows]

    pairwise_files = _select_pairwise_files(
        pairwise_dir,
        windows=windows,
        max_windows=args.max_windows,
    )
    if not pairwise_files:
        raise SystemExit("No pairwise parquet files matched the requested windows.")

    node_by_window = {
        window_id: group.copy()
        for window_id, group in sequence_df.groupby("window_id", sort=False)
    }
    edge_costs = _load_edge_costs(args.edge_manifest)

    # Build the per-file tasks up front. Each task corresponds to one
    # window-lineage pairwise parquet file.
    tasks: list[ScheduledTask] = []
    for pairwise_path in pairwise_files:
        window_id = _window_from_pairwise_path(pairwise_path)
        nodes = node_by_window.get(window_id)
        if nodes is None or nodes.empty:
            LOGGER.warning("Skipping %s: no node rows found", pairwise_path.stem)
            continue
        edge_cost, edge_cost_source = _edge_cost_for_pairwise_path(
            pairwise_path,
            edge_costs,
        )
        tasks.append(
            ScheduledTask(
                args=(
                    pairwise_path,
                    nodes,
                    attributes,
                    args.compatibility_threshold,
                    args.missing_label,
                    args.jackknife_blocks,
                    args.jackknife_seed,
                    args.force,
                ),
                edge_cost=edge_cost,
                edge_cost_source=edge_cost_source,
            )
        )

    total_tasks = len(tasks)
    if total_tasks == 0:
        raise SystemExit("No pairwise tasks could be built from the requested inputs.")

    small_tasks = [
        task for task in tasks if not _is_giant_task(task, args.giant_threshold)
    ]
    giant_tasks = sorted(
        (task for task in tasks if _is_giant_task(task, args.giant_threshold)),
        key=_descending_edge_cost,
        reverse=True,
    )

    if args.dry_run:
        _print_dry_run_schedule(
            small_tasks,
            giant_tasks,
            workers=args.workers,
            giant_workers=args.giant_workers,
            giant_threshold=args.giant_threshold,
            include_giants=args.include_giants,
        )
        return 0

    if giant_tasks and not args.include_giants:
        LOGGER.info(
            "Skipping %s giant pairwise files at threshold %s; pass "
            "--include-giants to process them",
            f"{len(giant_tasks):,}",
            f"{args.giant_threshold:,}",
        )
    process_tasks = [*small_tasks, *giant_tasks] if args.include_giants else small_tasks
    total_process_tasks = len(process_tasks)
    if total_process_tasks == 0:
        raise SystemExit(
            "No pairwise tasks left after skipping giant files. "
            "Pass --include-giants to process the selected giant files."
        )

    ensure_results_dirs()
    _ensure_intermediate_dirs()

    LOGGER.info(
        "Processing compatibility mixing for %s pairwise files (%s)",
        f"{total_process_tasks:,}",
        _format_task_progress(0, total_process_tasks),
    )

    statuses: dict[str, int] = {"ok": 0, "skipped": 0, "no_edges": 0}

    def _record(stem: str, status: str) -> None:
        statuses[status] = statuses.get(status, 0) + 1
        if status == "no_edges":
            LOGGER.warning("Skipping %s: no compatibility edges found", stem)

    def _log_progress(processed: int) -> None:
        is_progress_interval = processed % args.progress_every == 0
        if not is_progress_interval and processed != total_process_tasks:
            return
        LOGGER.info(
            "Processed %s pairwise files (%s)",
            _format_task_progress(processed, total_process_tasks),
            _format_status_counts(statuses),
        )

    processed_total = 0

    def _run_phase(
        phase_name: str,
        phase_tasks: Sequence[ScheduledTask],
        *,
        workers: int,
    ) -> None:
        nonlocal processed_total
        if not phase_tasks:
            LOGGER.info("Skipping %s phase: no pairwise files", phase_name)
            return

        LOGGER.info(
            "Starting %s phase: %s pairwise files, %s",
            phase_name,
            f"{len(phase_tasks):,}",
            f"workers={workers:,}",
        )

        if workers == 1:
            for task in phase_tasks:
                pairwise_path = task.pairwise_path
                LOGGER.debug(
                    "Processing compatibility mixing for %s (%s)",
                    pairwise_path.stem,
                    _format_task_progress(processed_total, total_process_tasks),
                )
                _record(*_process_pairwise_file(*task.args))
                processed_total += 1
                LOGGER.debug(
                    "Finished compatibility mixing for %s (%s)",
                    pairwise_path.stem,
                    _format_task_progress(processed_total, total_process_tasks),
                )
                _log_progress(processed_total)
            return

        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_stem = {
                executor.submit(_process_pairwise_file, *task.args): task.stem
                for task in phase_tasks
            }
            for future in as_completed(future_to_stem):
                stem = future_to_stem[future]
                _record(*future.result())
                processed_total += 1
                LOGGER.debug(
                    "Finished compatibility mixing for %s (%s)",
                    stem,
                    _format_task_progress(processed_total, total_process_tasks),
                )
                _log_progress(processed_total)

    # Run small files first, then the memory-heavy files. This avoids clumping
    # the largest BA.2-era parquets on 16 GB machines using macOS spawn.
    _run_phase("small", small_tasks, workers=args.workers)
    if args.include_giants:
        _run_phase("giant", giant_tasks, workers=args.giant_workers)

    LOGGER.info(
        "Pairwise mixing chunks complete: %s",
        _format_status_counts(statuses),
    )

    LOGGER.info("Writing compatibility mixing outputs")
    stems = [task.stem for task in process_tasks]
    matrix_table = _concat_intermediate_table("matrix", stems)
    summary_table = _concat_intermediate_table("summary", stems)
    topology_table = _concat_intermediate_table("topology", stems)

    write_table(matrix_table, FINAL_TABLES["matrix"], formats=("parquet",))
    write_table(
        summary_table,
        FINAL_TABLES["summary"],
        formats=("csv", "parquet"),
    )
    write_table(
        topology_table,
        FINAL_TABLES["topology"],
        formats=("csv", "parquet"),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
