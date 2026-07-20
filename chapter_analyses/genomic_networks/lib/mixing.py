"""Weighted mixing matrices and assortativity summaries.

This module implements weighted categorical assortativity, numeric
assortativity for degree/strength diagnostics, and multiplier-bootstrap
uncertainty for edge-weighted statistics.
"""

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


def prepare_graph_arrays(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    node_id_col: str = "sequence_id",
    source_col: str = "id1",
    target_col: str = "id2",
    weight_col: str | None = "epilink_compatibility",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, pd.Series]:
    """
    Convert node and edge dataframes into efficient NumPy arrays.

    Returns
    -------
    node_ids
        Node IDs in internal ordering.
    i, j
        Internal integer node indices for edge endpoints.
    w
        Edge weights.
    id_to_idx
        Mapping from original node ID to internal integer index.
    """
    if node_id_col not in nodes.columns:
        raise KeyError(f"Missing node id column: {node_id_col}")

    required_edge_cols = {source_col, target_col}
    if weight_col is not None:
        required_edge_cols.add(weight_col)

    missing = required_edge_cols - set(edges.columns)
    if missing:
        raise KeyError(f"Missing edge columns: {sorted(missing)}")

    if nodes[node_id_col].duplicated().any():
        raise ValueError("Duplicate node IDs found in nodes dataframe.")

    node_ids = nodes[node_id_col].to_numpy()
    id_to_idx = pd.Series(np.arange(len(nodes), dtype=np.int32), index=node_ids)

    edge_work = edges[[source_col, target_col]].copy()

    if weight_col is None:
        edge_work["_edge_weight"] = 1.0
    else:
        edge_work["_edge_weight"] = (
            pd.to_numeric(edges[weight_col], errors="coerce")
            .fillna(0.0)
            .to_numpy(dtype=np.float64)
        )

    edge_work = edge_work.dropna(subset=[source_col, target_col])
    edge_work = edge_work.loc[edge_work["_edge_weight"].gt(0)]

    i = edge_work[source_col].map(id_to_idx).to_numpy()
    j = edge_work[target_col].map(id_to_idx).to_numpy()

    if pd.isna(i).any() or pd.isna(j).any():
        raise ValueError("Some edges reference node IDs not present in nodes.")

    i = i.astype(np.int32)
    j = j.astype(np.int32)
    w = edge_work["_edge_weight"].to_numpy(dtype=np.float64)

    mask = i != j
    i = i[mask]
    j = j[mask]
    w = w[mask]

    if np.any(w < 0):
        raise ValueError(
            "Negative edge weights found. Assortativity usually assumes non-negative weights."
        )

    return node_ids, i, j, w, id_to_idx


def complete_mixing_matrix(
    matrix: np.ndarray,
    group_labels: pd.Index,
    all_categories: pd.Index,
    fill_value: float = 0.0,
    return_dataframe: bool = True,
) -> tuple[pd.DataFrame | np.ndarray, pd.Index]:
    """Expand a mixing matrix so that it includes all requested categories."""
    matrix = np.asarray(matrix)
    group_labels = pd.Index(group_labels)
    all_categories = pd.Index(all_categories)

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("`matrix` must be a square 2D array.")

    if matrix.shape[0] != len(group_labels):
        raise ValueError("`matrix` shape does not match length of `group_labels`.")

    if len(pd.Index(group_labels)) != len(pd.Index(group_labels).unique()):
        raise ValueError("`group_labels` contains duplicates.")

    if len(pd.Index(all_categories)) != len(pd.Index(all_categories).unique()):
        raise ValueError("`all_categories` contains duplicates.")

    extra_labels = set(group_labels) - set(all_categories)
    if extra_labels:
        raise ValueError(
            "`group_labels` contains labels not present in `all_categories`: "
            f"{sorted(extra_labels)}"
        )

    completed = np.full(
        (len(all_categories), len(all_categories)),
        fill_value,
        dtype=matrix.dtype,
    )

    label_to_pos = {label: pos for pos, label in enumerate(all_categories)}
    old_positions = np.array(
        [label_to_pos[label] for label in group_labels],
        dtype=np.int64,
    )

    completed[np.ix_(old_positions, old_positions)] = matrix

    if return_dataframe:
        return pd.DataFrame(
            completed, index=all_categories, columns=all_categories
        ), all_categories

    return completed, all_categories


def weighted_numeric_assortativity(
    i: np.ndarray,
    j: np.ndarray,
    w: np.ndarray,
    x: np.ndarray,
) -> float:
    """Weighted assortativity for a scalar node attribute x."""
    x = np.asarray(x, dtype=np.float64)

    if len(w) == 0 or np.sum(w) <= 0:
        return np.nan

    total_weight_directed = 2.0 * np.sum(w)

    sx = np.sum(w * (x[i] + x[j]))
    sxx = np.sum(w * (x[i] ** 2 + x[j] ** 2))
    sxy = 2.0 * np.sum(w * x[i] * x[j])

    mean_x = sx / total_weight_directed
    var_x = sxx / total_weight_directed - mean_x**2
    cov_xy = sxy / total_weight_directed - mean_x**2

    if var_x <= 0 or not np.isfinite(var_x):
        return np.nan

    return float(cov_xy / var_x)


def weighted_categorical_assortativity(
    i: np.ndarray,
    j: np.ndarray,
    w: np.ndarray,
    categories: pd.Index,
    levels: pd.Index | None = None,
) -> tuple[float, pd.DataFrame | np.ndarray, pd.Index]:
    """
    Weighted assortativity for categorical node attributes.

    Returns
    -------
    r
        Weighted categorical assortativity.
    mixing_matrix
        Normalized weighted mixing matrix as a DataFrame.
    labels
        Category labels corresponding to rows/columns.
    """
    codes, labels = pd.factorize(categories, sort=True)

    if np.any(codes < 0):
        raise ValueError(
            "Missing categorical values found. Please impute or filter first."
        )

    if len(w) == 0 or np.sum(w) <= 0:
        all_levels = levels if levels is not None else labels
        empty = np.zeros((len(all_levels), len(all_levels)), dtype=float)
        matrix, labels = complete_mixing_matrix(
            empty,
            group_labels=labels,
            all_categories=all_levels,
            fill_value=0.0,
            return_dataframe=True,
        )
        return np.nan, matrix, labels

    k = len(labels)
    flat_index = codes[i] * k + codes[j]

    e = np.bincount(
        flat_index,
        weights=w,
        minlength=k * k,
    ).reshape(k, k)

    # Undirected graph: count both directions.
    e = e + e.T

    total = e.sum()
    if total <= 0:
        r = np.nan
        e = e.astype(float)
    else:
        e = e / total
        a = e.sum(axis=1)
        expected_same = np.sum(a**2)
        observed_same = np.trace(e)
        denominator = 1.0 - expected_same
        r = (
            np.nan
            if denominator <= 0
            else float((observed_same - expected_same) / denominator)
        )

    e, labels = complete_mixing_matrix(
        e,
        group_labels=labels,
        all_categories=levels if levels is not None else labels,
        fill_value=0.0,
        return_dataframe=True,
    )

    return r, e, labels


def raw_categorical_mixing_matrix(
    i: np.ndarray,
    j: np.ndarray,
    w: np.ndarray,
    categories: pd.Index,
    levels: pd.Index | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return normalized weighted matrix and contribution-count matrix."""

    codes, labels = pd.factorize(categories, sort=True)

    if np.any(codes < 0):
        raise ValueError("Missing categorical values found.")

    k = len(labels)

    if len(w) == 0:
        weighted = np.zeros((k, k), dtype=float)
        counts = np.zeros((k, k), dtype=int)
    else:
        flat_index = codes[i] * k + codes[j]

        weighted = np.bincount(
            flat_index,
            weights=w,
            minlength=k * k,
        ).reshape(k, k)
        weighted = weighted + weighted.T

        counts = np.bincount(
            flat_index,
            weights=np.ones_like(w, dtype=float),
            minlength=k * k,
        ).reshape(k, k)
        counts = counts + counts.T
        counts = counts.astype(int)

    weighted_completed, _ = complete_mixing_matrix(
        weighted,
        group_labels=labels,
        all_categories=levels if levels is not None else labels,
        fill_value=0.0,
        return_dataframe=True,
    )
    counts_completed, _ = complete_mixing_matrix(
        counts,
        group_labels=labels,
        all_categories=levels if levels is not None else labels,
        fill_value=0,
        return_dataframe=True,
    )

    weighted_completed = pd.DataFrame(weighted_completed)
    counts_completed = pd.DataFrame(counts_completed)

    total = weighted_completed.to_numpy(float).sum()
    normalized = (
        weighted_completed / total if total > 0 else weighted_completed.astype(float)
    )

    return normalized, counts_completed


def multiplier_bootstrap(
    stat_fn,
    w: np.ndarray,
    B: int = 500,
    alpha: float = 0.05,
    seed: int = 123,
    chunk_size=None,
) -> tuple[float, tuple[float, float], float, np.ndarray]:
    """Multiplier/Bayesian bootstrap for an edge-weighted statistic."""
    rng = np.random.default_rng(seed)

    point = stat_fn(w)
    boot = np.empty(B, dtype=np.float64)

    for b in range(B):
        g = rng.exponential(scale=1.0, size=len(w))
        boot[b] = stat_fn(w * g)

    finite = boot[np.isfinite(boot)]
    if finite.size == 0:
        return point, (np.nan, np.nan), np.nan, boot

    lo, hi = np.quantile(finite, [alpha / 2, 1 - alpha / 2])
    se = finite.std(ddof=1) if finite.size > 1 else np.nan

    return point, (float(lo), float(hi)), float(se), boot


def node_strengths(
    i: np.ndarray,
    j: np.ndarray,
    w: np.ndarray,
    n_nodes: int,
) -> np.ndarray:
    """Weighted degree / node strength for an undirected graph."""
    strength = np.zeros(n_nodes, dtype=np.float64)
    np.add.at(strength, i, w)
    np.add.at(strength, j, w)
    return strength


def node_degrees(
    i: np.ndarray,
    j: np.ndarray,
    n_nodes: int,
) -> np.ndarray:
    """Unweighted degree for an undirected graph."""
    degree = np.zeros(n_nodes, dtype=np.float64)
    np.add.at(degree, i, 1.0)
    np.add.at(degree, j, 1.0)
    return degree


def _iter_edge_groups(
    edges: pd.DataFrame,
    group_cols: Sequence[str] | None,
):
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


def _dense_result_to_long_matrix(
    *,
    matrix: pd.DataFrame,
    contributions: pd.DataFrame,
    attribute: AttributeSpec,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    columns = [
        "attribute",
        "attribute_label",
        "source_category",
        "target_category",
        "edge_weight",
        "edge_contributions",
        "edge_weight_proportion",
    ]

    for source_category in matrix.index:
        for target_category in matrix.columns:
            edge_weight_proportion = matrix.loc[source_category, target_category]
            edge_contributions = contributions.loc[source_category, target_category]

            if edge_weight_proportion <= 0 and edge_contributions <= 0:
                continue

            rows.append(
                {
                    "attribute": attribute.name,
                    "attribute_label": attribute.label,
                    "source_category": str(source_category),
                    "target_category": str(target_category),
                    # Since matrix is normalized, this column is now a proportion.
                    # Kept for backward compatibility with the old schema.
                    "edge_weight": edge_weight_proportion,
                    "edge_contributions": edge_contributions,
                    "edge_weight_proportion": edge_weight_proportion,
                }
            )

    return pd.DataFrame(rows, columns=columns)


def _attribute_levels(
    nodes: pd.DataFrame,
    column: str,
    *,
    missing_label: str | None,
) -> list[Any]:
    values = nodes[column]
    if missing_label is not None:
        values = values.fillna(missing_label)
    else:
        values = values.dropna()

    return sorted(pd.unique(values.astype(str)).tolist())


def _labels_for_nodes(
    nodes: pd.DataFrame,
    column: str,
    *,
    missing_label: str | None,
) -> pd.Series:
    labels = nodes[column].copy()
    if missing_label is not None:
        labels = labels.fillna(missing_label)
    return labels.astype("string")


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
    bootstrap_replicates: int = 500,
    bootstrap_alpha: float = 0.05,
    bootstrap_seed: int = 123,
    min_edges: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build mixing matrices and assortativity summaries for an edge table."""
    if not symmetric:
        raise NotImplementedError(
            "The new implementation currently assumes undirected graphs."
        )

    matrix_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []

    attributes = tuple(attributes)

    for group_idx, (group_values, edge_group) in enumerate(
        _iter_edge_groups(edges, group_cols)
    ):
        n_edges_observed = len(edge_group)

        node_ids, i, j, w, _ = prepare_graph_arrays(
            nodes,
            edge_group,
            node_id_col=node_id_col,
            source_col=source_col,
            target_col=target_col,
            weight_col=weight_col,
        )

        n_edges_used = int(len(w))
        edge_weight_total = float(w.sum())

        for attr_idx, spec in enumerate(attributes):
            if spec.column not in nodes.columns:
                continue

            labels = _labels_for_nodes(
                nodes,
                spec.column,
                missing_label=missing_label,
            )

            valid_nodes = ~labels.isna()
            if len(w) > 0:
                edge_mask = valid_nodes.to_numpy()[i] & valid_nodes.to_numpy()[j]
            else:
                edge_mask = np.array([], dtype=bool)

            i_attr = i[edge_mask]
            j_attr = j[edge_mask]
            w_attr = w[edge_mask]

            if len(w_attr) < min_edges:
                summary_rows.append(
                    {
                        "attribute": spec.name,
                        "attribute_label": spec.label,
                        **group_values,
                        "n_edges_observed": n_edges_observed,
                        "n_edges_used": int(len(w_attr)),
                        "edge_weight_total": float(w_attr.sum()),
                        "n_categories": np.nan,
                        "assortativity": np.nan,
                        "observed_same_category_weight": np.nan,
                        "expected_same_category_weight": np.nan,
                        "uncertainty_method": "multiplier_bootstrap",
                        "bootstrap_replicates": int(bootstrap_replicates),
                        "bootstrap_finite_replicates": 0,
                        "assortativity_se": np.nan,
                        "assortativity_ci_low": np.nan,
                        "assortativity_ci_high": np.nan,
                        "skipped_reason": f"fewer_than_min_edges:{min_edges}",
                    }
                )
                continue

            attr_labels = labels.to_numpy()
            levels = _attribute_levels(
                nodes,
                spec.column,
                missing_label=missing_label,
            )

            r, normalized_matrix, completed_labels = weighted_categorical_assortativity(
                i_attr,
                j_attr,
                w_attr,
                attr_labels,
                levels=levels,
            )

            normalized_raw, contribution_matrix = raw_categorical_mixing_matrix(
                i_attr,
                j_attr,
                w_attr,
                attr_labels,
                levels=levels,
            )

            # Bootstrap uncertainty.
            if bootstrap_replicates > 0 and len(w_attr) > 0:
                seed = bootstrap_seed + group_idx * 10_000 + attr_idx

                def stat_fn(w_boot):
                    return weighted_categorical_assortativity(
                        i_attr,
                        j_attr,
                        w_boot,
                        attr_labels,
                        levels=levels,
                    )[0]

                point, ci, se, boot = multiplier_bootstrap(
                    stat_fn,
                    w_attr,
                    B=bootstrap_replicates,
                    alpha=bootstrap_alpha,
                    seed=seed,
                )
                finite_boot = int(np.isfinite(boot).sum())
            else:
                point = r
                ci = (np.nan, np.nan)
                se = np.nan
                finite_boot = 0

            matrix = _dense_result_to_long_matrix(
                matrix=normalized_raw,
                contributions=contribution_matrix,
                attribute=spec,
            )
            for key, value in group_values.items():
                matrix[key] = value
            matrix_parts.append(matrix)

            e_np = normalized_matrix.to_numpy(dtype=float)
            a = e_np.sum(axis=1)
            observed_same = float(np.trace(e_np)) if e_np.sum() > 0 else np.nan
            expected_same = float(np.sum(a**2)) if e_np.sum() > 0 else np.nan

            summary_rows.append(
                {
                    "attribute": spec.name,
                    "attribute_label": spec.label,
                    **group_values,
                    "n_edges_observed": n_edges_observed,
                    "n_edges_used": int(len(w_attr)),
                    "edge_weight_total": float(w_attr.sum()),
                    "n_categories": int(len(completed_labels)),
                    "assortativity": float(point) if np.isfinite(point) else np.nan,
                    "observed_same_category_weight": observed_same,
                    "expected_same_category_weight": expected_same,
                    "uncertainty_method": "multiplier_bootstrap",
                    "bootstrap_replicates": int(bootstrap_replicates),
                    "bootstrap_finite_replicates": finite_boot,
                    "assortativity_se": se,
                    "assortativity_ci_low": ci[0],
                    "assortativity_ci_high": ci[1],
                    "skipped_reason": None,
                }
            )

    matrix_table = (
        pd.concat(matrix_parts, ignore_index=True, sort=False)
        if matrix_parts
        else pd.DataFrame()
    )
    summary_table = pd.DataFrame(summary_rows)

    return matrix_table, summary_table


def build_degree_assortativity_for_edge_table(
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    *,
    node_id_col: str = "sequence_id",
    source_col: str = "id1",
    target_col: str = "id2",
    weight_col: str | None = "epilink_compatibility",
    group_cols: Sequence[str] | None = ("window_id",),
    bootstrap_replicates: int = 500,
    bootstrap_alpha: float = 0.05,
    bootstrap_seed: int = 123,
    min_edges: int = 0,
) -> pd.DataFrame:
    """Build degree/strength assortativity diagnostics for each edge group."""
    rows: list[dict[str, Any]] = []

    for group_idx, (group_values, edge_group) in enumerate(
        _iter_edge_groups(edges, group_cols)
    ):
        n_edges_observed = len(edge_group)

        node_ids, i, j, w, _ = prepare_graph_arrays(
            nodes,
            edge_group,
            node_id_col=node_id_col,
            source_col=source_col,
            target_col=target_col,
            weight_col=weight_col,
        )

        n_nodes = len(node_ids)
        n_edges_used = len(w)

        if n_nodes == 0 or n_edges_used == 0 or n_edges_used < min_edges:
            rows.append(
                {
                    **group_values,
                    "n_nodes": int(n_nodes),
                    "n_edges_observed": int(n_edges_observed),
                    "n_edges_used": int(n_edges_used),
                    "edge_weight_total": float(w.sum()),
                    "mean_degree": np.nan,
                    "max_degree": np.nan,
                    "mean_strength": np.nan,
                    "max_strength": np.nan,
                    "degree_assortativity": np.nan,
                    "weighted_degree_assortativity": np.nan,
                    "strength_assortativity": np.nan,
                    "strength_assortativity_se": np.nan,
                    "strength_assortativity_ci_low": np.nan,
                    "strength_assortativity_ci_high": np.nan,
                    "uncertainty_method": "multiplier_bootstrap",
                    "bootstrap_replicates": int(bootstrap_replicates),
                    "bootstrap_finite_replicates": 0,
                    "skipped_reason": (
                        f"fewer_than_min_edges:{min_edges}"
                        if n_edges_used < min_edges
                        else "no_edges"
                    ),
                }
            )
            continue

        degree = node_degrees(i, j, n_nodes)
        strength = node_strengths(i, j, w, n_nodes)

        degree_r = weighted_numeric_assortativity(i, j, np.ones_like(w), degree)
        weighted_degree_r = weighted_numeric_assortativity(i, j, w, degree)
        strength_r = weighted_numeric_assortativity(i, j, w, strength)

        if bootstrap_replicates > 0:
            seed = bootstrap_seed + group_idx * 10_000

            def stat_fn(w_boot):
                strength_boot = node_strengths(i, j, w_boot, n_nodes)
                return weighted_numeric_assortativity(i, j, w_boot, strength_boot)

            point, ci, se, boot = multiplier_bootstrap(
                stat_fn,
                w,
                B=bootstrap_replicates,
                alpha=bootstrap_alpha,
                seed=seed,
            )
            finite_boot = int(np.isfinite(boot).sum())
        else:
            point = strength_r
            ci = (np.nan, np.nan)
            se = np.nan
            finite_boot = 0

        rows.append(
            {
                **group_values,
                "n_nodes": int(n_nodes),
                "n_edges_observed": int(n_edges_observed),
                "n_edges_used": int(n_edges_used),
                "edge_weight_total": float(w.sum()),
                "mean_degree": float(degree.mean()),
                "max_degree": float(degree.max()),
                "mean_strength": float(strength.mean()),
                "max_strength": float(strength.max()),
                "degree_assortativity": degree_r,
                "weighted_degree_assortativity": weighted_degree_r,
                "strength_assortativity": float(point)
                if np.isfinite(point)
                else np.nan,
                "strength_assortativity_se": se,
                "strength_assortativity_ci_low": ci[0],
                "strength_assortativity_ci_high": ci[1],
                "uncertainty_method": "multiplier_bootstrap",
                "bootstrap_replicates": int(bootstrap_replicates),
                "bootstrap_finite_replicates": finite_boot,
                "skipped_reason": None,
            }
        )

    return pd.DataFrame(rows)
