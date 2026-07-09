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


def categorical_mixing_matrix(
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    *,
    attribute: AttributeSpec,
    node_id_col: str,
    source_col: str,
    target_col: str,
    weight_col: str | None = None,
    symmetric: bool = True,
    missing_label: str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build a long weighted categorical mixing matrix for one attribute.

    When ``symmetric=True``, every edge contributes in both directions. This is
    the usual representation for an undirected compatibility network and makes
    row and column marginals comparable. Directed transition graphs can pass
    ``symmetric=False`` to retain source-to-target orientation.
    """
    required = {source_col, target_col}
    missing = required - set(edges.columns)
    if missing:
        raise KeyError(f"Missing edge columns: {sorted(missing)}")
    if attribute.column not in nodes.columns:
        raise KeyError(f"Missing node attribute column: {attribute.column}")

    lookup = node_attribute_lookup(
        nodes,
        node_id_col=node_id_col,
        attributes=(attribute,),
    )
    source_lookup = lookup.rename(
        columns={
            node_id_col: source_col,
            attribute.column: "source_category",
        }
    )
    target_lookup = lookup.rename(
        columns={
            node_id_col: target_col,
            attribute.column: "target_category",
        }
    )

    edge_cols = [source_col, target_col]
    if weight_col is not None and weight_col in edges.columns:
        edge_cols.append(weight_col)
    work = edges[edge_cols].copy()
    if weight_col is None or weight_col not in work.columns:
        work["_edge_weight"] = 1.0
        weight_col = "_edge_weight"

    work = work.merge(source_lookup, on=source_col, how="left").merge(
        target_lookup,
        on=target_col,
        how="left",
    )
    n_edges_observed = len(work)

    if missing_label is None:
        work = work.dropna(subset=["source_category", "target_category"])
    else:
        work["source_category"] = work["source_category"].fillna(missing_label)
        work["target_category"] = work["target_category"].fillna(missing_label)

    work["source_category"] = work["source_category"].astype("string")
    work["target_category"] = work["target_category"].astype("string")
    work["_weight"] = pd.to_numeric(work[weight_col], errors="coerce").fillna(0.0)
    work = work.loc[work["_weight"].gt(0)].copy()

    if symmetric and not work.empty:
        rev = work.rename(
            columns={
                "source_category": "target_category",
                "target_category": "source_category",
            }
        )
        work = pd.concat([work, rev], ignore_index=True, sort=False)

    matrix = (
        work.groupby(["source_category", "target_category"], dropna=False)
        .agg(
            edge_weight=("_weight", "sum"),
            edge_contributions=("_weight", "size"),
        )
        .reset_index()
    )
    total = matrix["edge_weight"].sum()
    matrix["edge_weight_proportion"] = (
        matrix["edge_weight"] / total if total > 0 else np.nan
    )

    diagnostics = {
        "n_edges_observed": n_edges_observed,
        "n_edge_contributions_used": int(matrix["edge_contributions"].sum()),
        "edge_weight_total": float(total),
        "n_categories": int(
            pd.concat([matrix["source_category"], matrix["target_category"]])
            .dropna()
            .nunique()
        )
        if not matrix.empty
        else 0,
    }
    return matrix, diagnostics


def assortativity_from_matrix(matrix: pd.DataFrame) -> dict[str, float]:
    """Compute categorical assortativity from a long mixing matrix."""
    if matrix.empty or matrix["edge_weight"].sum() <= 0:
        return {
            "assortativity": np.nan,
            "observed_same_category_weight": np.nan,
            "expected_same_category_weight": np.nan,
        }

    categories = sorted(
        set(matrix["source_category"].dropna().astype(str))
        | set(matrix["target_category"].dropna().astype(str))
    )
    wide = (
        matrix.pivot_table(
            index="source_category",
            columns="target_category",
            values="edge_weight",
            aggfunc="sum",
            fill_value=0.0,
        )
        .reindex(index=categories, columns=categories, fill_value=0.0)
        .astype(float)
    )
    e = wide / wide.to_numpy().sum()
    observed = float(np.trace(e.to_numpy()))
    expected = float((e.sum(axis=1).to_numpy() * e.sum(axis=0).to_numpy()).sum())
    denom = 1.0 - expected
    r = np.nan if denom <= 0 else (observed - expected) / denom
    return {
        "assortativity": float(r) if not pd.isna(r) else np.nan,
        "observed_same_category_weight": observed,
        "expected_same_category_weight": expected,
    }


def _iter_edge_groups(
    edges: pd.DataFrame,
    group_cols: Sequence[str] | None,
)-> Iterable[tuple[dict[str, Any], pd.DataFrame]]:
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
    matrix_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    attributes = tuple(attributes)

    for group_values, edge_group in _iter_edge_groups(edges, group_cols):
        for spec in attributes:
            if spec.column not in nodes.columns:
                continue
            matrix, diagnostics = categorical_mixing_matrix(
                edge_group,
                nodes,
                attribute=spec,
                node_id_col=node_id_col,
                source_col=source_col,
                target_col=target_col,
                weight_col=weight_col,
                symmetric=symmetric,
                missing_label=missing_label,
            )
            for key, value in group_values.items():
                matrix[key] = value
            matrix.insert(0, "attribute", spec.name)
            matrix.insert(1, "attribute_label", spec.label)
            matrix_parts.append(matrix)

            summary = {
                "attribute": spec.name,
                "attribute_label": spec.label,
                **group_values,
                **diagnostics,
                **assortativity_from_matrix(matrix),
            }
            summary_rows.append(summary)

    matrix_table = (
        pd.concat(matrix_parts, ignore_index=True, sort=False)
        if matrix_parts
        else pd.DataFrame()
    )
    summary_table = pd.DataFrame(summary_rows)
    return matrix_table, summary_table


def transition_attribute_specs(
    cluster_table: pd.DataFrame,
    *,
    attributes: Iterable[AttributeSpec] = DEFAULT_MIXING_ATTRIBUTES,
) -> tuple[AttributeSpec, ...]:
    """Return specs pointing at modal cluster-attribute columns."""
    out: list[AttributeSpec] = []
    for spec in attributes:
        modal_col = f"modal_{spec.name}"
        if modal_col in cluster_table.columns:
            out.append(
                AttributeSpec(
                    name=spec.name,
                    column=modal_col,
                    label=spec.label,
                    ordered=spec.ordered,
                )
            )
    return tuple(out)


def build_transition_mixing(
    edge_table: pd.DataFrame,
    cluster_table: pd.DataFrame,
    *,
    attributes: Iterable[AttributeSpec] = DEFAULT_MIXING_ATTRIBUTES,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build directed cluster-level mixing for the temporal transition graph."""
    specs = transition_attribute_specs(cluster_table, attributes=attributes)
    return build_mixing_for_edge_table(
        edge_table,
        cluster_table,
        attributes=specs,
        node_id_col="cluster_id",
        source_col="source",
        target_col="target",
        weight_col="n_shared_sequences",
        group_cols=("source_window_id",),
        symmetric=False,
        missing_label=None,
    )

