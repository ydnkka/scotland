"""Weighted categorical mixing matrices and assortativity summaries."""

from __future__ import annotations

from typing import Iterable, Sequence, Any

import numpy as np
import pandas as pd

from .config import DEFAULT_MIXING_ATTRIBUTES, AttributeSpec


def specs_by_name(
    names: Sequence[str] | None,
    *,
    specs: Iterable[AttributeSpec] = DEFAULT_MIXING_ATTRIBUTES,
) -> tuple[AttributeSpec, ...]:
    """Return attribute specs matching ``names`` while preserving spec order."""
    available = tuple(specs)
    if names is None:
        return available
    wanted = set(names)
    unknown = wanted - {spec.name for spec in available}
    if unknown:
        raise ValueError(f"Unknown attribute spec(s): {sorted(unknown)}")
    return tuple(spec for spec in available if spec.name in wanted)


def _first_non_missing(values: pd.Series) -> Any:
    values = values.dropna()
    if values.empty:
        return np.nan
    return values.iloc[0]


def node_attribute_lookup(
    nodes: pd.DataFrame,
    *,
    node_id_col: str,
    attributes: Iterable[AttributeSpec],
) -> pd.DataFrame:
    """Return one attribute row per node."""
    attr_cols = [spec.column for spec in attributes if spec.column in nodes.columns]
    if node_id_col not in nodes.columns:
        raise KeyError(f"Missing node id column: {node_id_col}")
    if not attr_cols:
        return nodes[[node_id_col]].drop_duplicates().copy()

    work = nodes[[node_id_col, *attr_cols]].copy()
    if work[node_id_col].duplicated().any():
        work = (
            work.groupby(node_id_col, dropna=False)[attr_cols]
            .agg(_first_non_missing)
            .reset_index()
        )
    else:
        work = work.drop_duplicates(node_id_col)
    return work


def _as_finite_weight_array(
    edges: pd.DataFrame,
    *,
    weight_col: str | None,
) -> np.ndarray:
    if weight_col is None or weight_col not in edges.columns:
        return np.ones(len(edges), dtype=float)
    return (
        pd.to_numeric(edges[weight_col], errors="coerce")
        .fillna(0.0)
        .to_numpy(dtype=float)
    )


def _assortativity_components_from_dense(
    mixing_matrix: np.ndarray,
) -> dict[str, float]:
    total = float(mixing_matrix.sum())
    if total <= 0:
        return {
            "assortativity": np.nan,
            "observed_same_category_weight": np.nan,
            "expected_same_category_weight": np.nan,
        }

    e = mixing_matrix / total
    observed = float(np.trace(e))
    expected = float(np.dot(e.sum(axis=1), e.sum(axis=0)))
    denominator = 1.0 - expected
    r = np.nan if np.isclose(denominator, 0.0) else (observed - expected) / denominator
    return {
        "assortativity": float(r) if not pd.isna(r) else np.nan,
        "observed_same_category_weight": observed,
        "expected_same_category_weight": expected,
    }


def _dense_mixing_for_categories(
    category_indices: np.ndarray,
    source_vertices: np.ndarray,
    target_vertices: np.ndarray,
    edge_weights: np.ndarray,
    *,
    n_categories: int,
    symmetric: bool,
) -> np.ndarray:
    matrix = np.zeros((n_categories, n_categories), dtype=float)
    if source_vertices.size == 0:
        return matrix

    source_groups = category_indices[source_vertices]
    target_groups = category_indices[target_vertices]
    flat = source_groups * n_categories + target_groups
    matrix = np.bincount(
        flat, weights=edge_weights, minlength=n_categories * n_categories
    ).reshape(n_categories, n_categories)
    if symmetric:
        matrix = matrix + matrix.T
    return matrix


def _edge_arrays_from_edge_table(
    edges: pd.DataFrame,
    *,
    source_col: str,
    target_col: str,
    weight_col: str | None,
) -> dict[str, Any]:
    required = {source_col, target_col}
    missing = required - set(edges.columns)
    if missing:
        raise KeyError(f"Missing edge columns: {sorted(missing)}")

    work = edges[[source_col, target_col]].copy()
    work["_edge_weight"] = _as_finite_weight_array(edges, weight_col=weight_col)
    work = work.dropna(subset=[source_col, target_col])
    work = work.loc[work["_edge_weight"].gt(0)]

    if work.empty:
        return {
            "n_edges_observed": len(edges),
            "vertex_names": pd.Index([], dtype="object"),
            "source_vertices": np.array([], dtype=int),
            "target_vertices": np.array([], dtype=int),
            "edge_weights": np.array([], dtype=float),
        }

    endpoints = work[[source_col, target_col]].to_numpy().ravel()
    vertex_codes, vertex_names = pd.factorize(endpoints, sort=False)
    return {
        "n_edges_observed": len(edges),
        "vertex_names": vertex_names,
        "source_vertices": vertex_codes[0::2].astype(int),
        "target_vertices": vertex_codes[1::2].astype(int),
        "edge_weights": work["_edge_weight"].to_numpy(dtype=float),
    }


def _labels_for_vertices(
    vertex_names: pd.Index,
    attr_lookup: pd.DataFrame,
    *,
    attribute: AttributeSpec,
    missing_label: str | None,
) -> pd.Series:
    labels = pd.Series(vertex_names, dtype="object").map(attr_lookup[attribute.column])
    if missing_label is not None:
        labels = labels.fillna(missing_label)
    return labels


def _weighted_pearson_correlation(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    weights: Sequence[float] | np.ndarray | None = None,
) -> float:
    """Return a weighted Pearson correlation, or NaN if undefined."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if weights is None:
        weights = np.ones_like(x, dtype=float)
    else:
        weights = np.asarray(weights, dtype=float)

    if not (x.shape == y.shape == weights.shape):
        raise ValueError("x, y, and weights must align")

    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0)
    if not mask.any():
        return np.nan

    x = x[mask]
    y = y[mask]
    weights = weights[mask]
    weight_total = weights.sum()
    if weight_total <= 0:
        return np.nan

    x_mean = np.average(x, weights=weights)
    y_mean = np.average(y, weights=weights)
    x_centered = x - x_mean
    y_centered = y - y_mean
    covariance = np.average(x_centered * y_centered, weights=weights)
    x_variance = np.average(x_centered**2, weights=weights)
    y_variance = np.average(y_centered**2, weights=weights)
    denominator = np.sqrt(x_variance * y_variance)
    if np.isclose(denominator, 0.0):
        return np.nan
    return float(covariance / denominator)


def weighted_nominal_assortativity(
    labels: pd.Series[Any] | Sequence[Any],
    source_vertices: Sequence[int] | np.ndarray,
    target_vertices: Sequence[int] | np.ndarray,
    edge_weights: Sequence[float] | np.ndarray,
    *,
    symmetric: bool = True,
    missing_label: str | None = None,
) -> dict[str, Any]:
    """Compute weighted nominal assortativity from compact edge arrays.

    The inputs are compact NumPy-style edge arrays: each edge references source
    and target vertex indices, and ``labels`` gives one categorical label per
    vertex. This avoids repeated pandas joins/groupbys during mixing summaries.
    """
    labels = pd.Series(labels, dtype="object")
    source_vertices = np.asarray(source_vertices, dtype=int)
    target_vertices = np.asarray(target_vertices, dtype=int)
    edge_weights = np.asarray(edge_weights, dtype=float)

    if not (
        source_vertices.shape[0] == target_vertices.shape[0] == edge_weights.shape[0]
    ):
        raise ValueError(
            "source_vertices, target_vertices, and edge_weights must align"
        )
    if not np.all(np.isfinite(edge_weights)):
        raise ValueError("edge weights must be finite numbers")
    if source_vertices.size:
        max_vertex = int(max(source_vertices.max(), target_vertices.max()))
        if max_vertex >= len(labels):
            raise ValueError("edge vertex index exceeds labels length")

    if missing_label is not None:
        labels = labels.fillna(missing_label)
    valid_vertices = ~labels.isna()
    labels = labels.astype("string")

    if source_vertices.size:
        edge_mask = (
            valid_vertices.to_numpy()[source_vertices]
            & valid_vertices.to_numpy()[target_vertices]
            & (edge_weights > 0)
        )
        source_vertices = source_vertices[edge_mask]
        target_vertices = target_vertices[edge_mask]
        edge_weights = edge_weights[edge_mask]

    if source_vertices.size:
        used_vertices = np.unique(np.concatenate([source_vertices, target_vertices]))
    else:
        used_vertices = np.array([], dtype=int)

    if used_vertices.size:
        remap = np.full(len(labels), -1, dtype=int)
        remap[used_vertices] = np.arange(used_vertices.size)
        source_vertices = remap[source_vertices]
        target_vertices = remap[target_vertices]
        used_labels = labels.iloc[used_vertices].astype(str).to_numpy()
        unique_labels = np.array(
            sorted(pd.unique(used_labels).tolist()),
            dtype=object,
        )
        label_to_index = {label: idx for idx, label in enumerate(unique_labels)}
        category_indices = np.fromiter(
            (label_to_index[label] for label in used_labels),
            dtype=int,
            count=len(used_labels),
        )
    else:
        unique_labels = np.array([], dtype=object)
        category_indices = np.array([], dtype=int)

    n_categories = len(unique_labels)
    mixing_matrix = _dense_mixing_for_categories(
        category_indices,
        source_vertices,
        target_vertices,
        edge_weights,
        n_categories=n_categories,
        symmetric=symmetric,
    )
    contribution_matrix = _dense_mixing_for_categories(
        category_indices,
        source_vertices,
        target_vertices,
        np.ones_like(edge_weights, dtype=float),
        n_categories=n_categories,
        symmetric=symmetric,
    )

    components = _assortativity_components_from_dense(mixing_matrix)
    r_obs = components["assortativity"]

    total_weight = mixing_matrix.sum()
    normalized_mixing = (
        mixing_matrix / total_weight if total_weight > 0 else mixing_matrix
    )
    return {
        "observed_r": r_obs,
        "assortativity": r_obs,
        "mixing_matrix": normalized_mixing,
        "mixing_matrix_raw_weights": mixing_matrix,
        "mixing_matrix_contributions": contribution_matrix.astype(int),
        "groups": unique_labels.tolist(),
        **components,
    }


def _dense_result_to_long_matrix(
    result: dict[str, Any],
    *,
    attribute: AttributeSpec,
) -> pd.DataFrame:
    groups = result["groups"]
    weights = np.asarray(result["mixing_matrix_raw_weights"], dtype=float)
    contributions = np.asarray(result["mixing_matrix_contributions"], dtype=int)
    rows: list[dict[str, Any]] = []
    total = float(weights.sum())
    columns = [
        "attribute",
        "attribute_label",
        "source_category",
        "target_category",
        "edge_weight",
        "edge_contributions",
        "edge_weight_proportion",
    ]
    for source_idx, source_category in enumerate(groups):
        for target_idx, target_category in enumerate(groups):
            edge_weight = float(weights[source_idx, target_idx])
            edge_contributions = int(contributions[source_idx, target_idx])
            if edge_weight <= 0 and edge_contributions <= 0:
                continue
            rows.append(
                {
                    "attribute": attribute.name,
                    "attribute_label": attribute.label,
                    "source_category": str(source_category),
                    "target_category": str(target_category),
                    "edge_weight": edge_weight,
                    "edge_contributions": edge_contributions,
                    "edge_weight_proportion": (
                        edge_weight / total if total > 0 else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows, columns=columns)


def degree_strength_assortativity_from_edge_arrays(
    source_vertices: Sequence[int] | np.ndarray,
    target_vertices: Sequence[int] | np.ndarray,
    edge_weights: Sequence[float] | np.ndarray,
    *,
    n_vertices: int,
) -> dict[str, float]:
    """Compute degree/strength assortativity diagnostics for an edge array.

    ``degree_assortativity`` is the ordinary edge-level Pearson correlation of
    endpoint degrees. ``weighted_degree_assortativity`` uses edge weights in
    that correlation. ``strength_assortativity`` correlates endpoint strengths
    using edge weights.
    """
    source_vertices = np.asarray(source_vertices, dtype=int)
    target_vertices = np.asarray(target_vertices, dtype=int)
    edge_weights = np.asarray(edge_weights, dtype=float)
    if not (
        source_vertices.shape[0] == target_vertices.shape[0] == edge_weights.shape[0]
    ):
        raise ValueError(
            "source_vertices, target_vertices, and edge_weights must align"
        )
    if n_vertices < 0:
        raise ValueError("n_vertices must be non-negative")
    if n_vertices == 0 or source_vertices.size == 0:
        return {
            "n_nodes": int(n_vertices),
            "n_edges_used": int(source_vertices.size),
            "edge_weight_total": float(edge_weights.sum()),
            "mean_degree": np.nan,
            "max_degree": np.nan,
            "mean_strength": np.nan,
            "max_strength": np.nan,
            "degree_assortativity": np.nan,
            "weighted_degree_assortativity": np.nan,
            "strength_assortativity": np.nan,
        }

    degree = (
        np.bincount(source_vertices, minlength=n_vertices)
        + np.bincount(target_vertices, minlength=n_vertices)
    ).astype(float)

    strength = np.bincount(
        source_vertices, weights=edge_weights, minlength=n_vertices
    ) + np.bincount(target_vertices, weights=edge_weights, minlength=n_vertices)

    source_degree = degree[source_vertices]
    target_degree = degree[target_vertices]
    source_strength = strength[source_vertices]
    target_strength = strength[target_vertices]
    return {
        "n_nodes": int(n_vertices),
        "n_edges_used": int(source_vertices.size),
        "edge_weight_total": float(edge_weights.sum()),
        "mean_degree": float(degree.mean()),
        "max_degree": float(degree.max()),
        "mean_strength": float(strength.mean()),
        "max_strength": float(strength.max()),
        "degree_assortativity": _weighted_pearson_correlation(
            source_degree,
            target_degree,
        ),
        "weighted_degree_assortativity": _weighted_pearson_correlation(
            source_degree,
            target_degree,
            edge_weights,
        ),
        "strength_assortativity": _weighted_pearson_correlation(
            source_strength,
            target_strength,
            edge_weights,
        ),
    }


def _iter_edge_groups(
    edges: pd.DataFrame,
    group_cols: Sequence[str] | None,
) -> Iterable[tuple[dict[str, Any], pd.DataFrame]]:
    group_cols = [col for col in (group_cols or []) if col]
    if not group_cols:
        yield {}, edges
        return

    missing = [col for col in group_cols if col not in edges.columns]
    if missing:
        raise KeyError(f"Missing edge grouping columns: {missing}")

    grouped = edges.groupby(group_cols, dropna=False, sort=True)
    for key, group in grouped:
        if not isinstance(key, tuple):
            key = (key,)
        yield dict(zip(group_cols, key)), group


def build_degree_assortativity_for_edge_table(
    edges: pd.DataFrame,
    *,
    source_col: str = "id1",
    target_col: str = "id2",
    weight_col: str | None = "epilink_compatibility",
    group_cols: Sequence[str] | None = ("window_id",),
) -> pd.DataFrame:
    """Build degree/strength assortativity diagnostics for each edge group."""
    rows: list[dict[str, Any]] = []
    for group_values, edge_group in _iter_edge_groups(edges, group_cols):
        edge_arrays = _edge_arrays_from_edge_table(
            edge_group,
            source_col=source_col,
            target_col=target_col,
            weight_col=weight_col,
        )
        diagnostics = degree_strength_assortativity_from_edge_arrays(
            edge_arrays["source_vertices"],
            edge_arrays["target_vertices"],
            edge_arrays["edge_weights"],
            n_vertices=len(edge_arrays["vertex_names"]),
        )
        rows.append(
            {
                **group_values,
                "n_edges_observed": edge_arrays["n_edges_observed"],
                **diagnostics,
            }
        )
    return pd.DataFrame(rows)


def _build_mixing_for_edge_table_numpy(
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    *,
    attributes: Iterable[AttributeSpec],
    node_id_col: str,
    source_col: str,
    target_col: str,
    weight_col: str | None,
    group_cols: Sequence[str] | None,
    symmetric: bool,
    missing_label: str | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build mixing outputs using pandas preparation and NumPy arrays."""
    matrix_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    attributes = tuple(attributes)

    lookup = node_attribute_lookup(
        nodes,
        node_id_col=node_id_col,
        attributes=attributes,
    )
    attr_cols = [spec.column for spec in attributes if spec.column in lookup.columns]
    attr_lookup = (
        lookup.set_index(node_id_col)[attr_cols] if attr_cols else pd.DataFrame()
    )

    for group_values, edge_group in _iter_edge_groups(edges, group_cols):
        edge_arrays = _edge_arrays_from_edge_table(
            edge_group,
            source_col=source_col,
            target_col=target_col,
            weight_col=weight_col,
        )

        for spec in attributes:
            if spec.column not in attr_lookup.columns:
                continue

            labels = _labels_for_vertices(
                edge_arrays["vertex_names"],
                attr_lookup,
                attribute=spec,
                missing_label=missing_label,
            )
            result = weighted_nominal_assortativity(
                labels,
                edge_arrays["source_vertices"],
                edge_arrays["target_vertices"],
                edge_arrays["edge_weights"],
                symmetric=symmetric,
                missing_label=missing_label,
            )

            matrix = _dense_result_to_long_matrix(result, attribute=spec)
            for key, value in group_values.items():
                matrix[key] = value
            matrix_parts.append(matrix)

            summary = {
                "attribute": spec.name,
                "attribute_label": spec.label,
                **group_values,
                "n_edges_observed": edge_arrays["n_edges_observed"],
                "n_edge_contributions_used": int(
                    np.asarray(result["mixing_matrix_contributions"]).sum()
                ),
                "edge_weight_total": float(
                    np.asarray(result["mixing_matrix_raw_weights"]).sum()
                ),
                "n_categories": len(result["groups"]),
                "assortativity": result["assortativity"],
                "observed_same_category_weight": result[
                    "observed_same_category_weight"
                ],
                "expected_same_category_weight": result[
                    "expected_same_category_weight"
                ],
            }
            summary_rows.append(summary)

    matrix_table = (
        pd.concat(matrix_parts, ignore_index=True, sort=False)
        if matrix_parts
        else pd.DataFrame()
    )
    summary_table = pd.DataFrame(summary_rows)
    return matrix_table, summary_table


def build_mixing_for_edge_table(
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    *,
    attributes: Iterable[AttributeSpec] = DEFAULT_MIXING_ATTRIBUTES,
    node_id_col: str = "sequence_id",
    source_col: str = "id1",
    target_col: str = "id2",
    weight_col: str | None = "epilink_compatibility",
    group_cols: Sequence[str] | None = ("window_id",),
    symmetric: bool = True,
    missing_label: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build mixing matrices and assortativity summaries for an edge table."""
    return _build_mixing_for_edge_table_numpy(
        edges,
        nodes,
        attributes=attributes,
        node_id_col=node_id_col,
        source_col=source_col,
        target_col=target_col,
        weight_col=weight_col,
        group_cols=group_cols,
        symmetric=symmetric,
        missing_label=missing_label,
    )

