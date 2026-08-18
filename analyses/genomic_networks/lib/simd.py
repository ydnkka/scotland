"""SIMD population-weighting validation tables."""

from __future__ import annotations

import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .config import PROJECT_ROOT

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import load_datazone_info

SOURCE_GROUP_COLUMNS = {
    5: "dz_simd_quintile",
    10: "dz_simd_decile",
    20: "dz_simd_vigintile",
}

@dataclass(frozen=True)
class SimdValidationTables:
    """Container for generated SIMD validation outputs."""

    datazone_assignments: pd.DataFrame
    group_summary: pd.DataFrame
    movement_table: pd.DataFrame
    change_summary: pd.DataFrame
    diagnostics: pd.DataFrame


def load_simd_datazone_lookup(n_groups: int = 5) -> pd.DataFrame:
    """Load Data Zone SIMD rank, population, and source grouping columns."""
    if n_groups not in SOURCE_GROUP_COLUMNS:
        raise ValueError(f"Unsupported SIMD group count: {n_groups}")

    source_col = SOURCE_GROUP_COLUMNS[n_groups]
    lookup = load_datazone_info(
        [
            "dz_population",
            "dz_simd_rank",
            source_col,
        ],
        weighted_simd=False,
    )
    if "datazone" not in lookup.columns:
        lookup = lookup.reset_index()

    out = (
        lookup.drop(columns="geometry", errors="ignore")
        .rename(
            columns={
                "dz_population": "population",
                "dz_simd_rank": "simd_rank",
                source_col: "source_group",
            }
        )
        .dropna(subset=["datazone", "population", "simd_rank"])
        .copy()
    )
    out["population"] = pd.to_numeric(out["population"], errors="raise").astype(int)
    out["simd_rank"] = pd.to_numeric(out["simd_rank"], errors="raise").astype(int)
    out["source_group"] = pd.to_numeric(out["source_group"], errors="raise").astype(int)
    return out.sort_values("simd_rank").reset_index(drop=True)


def assign_equal_datazone_groups(
    df: pd.DataFrame,
    *,
    n_groups: int,
    rank_col: str = "simd_rank",
    group_col: str = "equal_datazone_group",
) -> pd.DataFrame:
    """Assign groups by equal numbers of Data Zones ordered by SIMD rank."""
    if n_groups < 2:
        raise ValueError("n_groups must be at least 2")
    if df[rank_col].isna().any():
        raise ValueError(f"Missing values found in {rank_col}")

    out = df.sort_values(rank_col, ascending=True).copy()
    positions = np.arange(1, len(out) + 1, dtype=float)
    out[group_col] = np.ceil(positions * n_groups / len(out)).astype(int)
    out[group_col] = out[group_col].clip(1, n_groups)
    return out.sort_index()


def assign_population_weighted_groups(
    df: pd.DataFrame,
    *,
    n_groups: int,
    rank_col: str = "simd_rank",
    population_col: str = "population",
    group_col: str = "population_weighted_group",
) -> pd.DataFrame:
    """Assign groups by cumulative population ordered by SIMD rank."""
    if n_groups < 2:
        raise ValueError("n_groups must be at least 2")
    if df[rank_col].isna().any():
        raise ValueError(f"Missing values found in {rank_col}")
    if df[population_col].isna().any():
        raise ValueError(f"Missing values found in {population_col}")
    if (df[population_col] < 0).any():
        raise ValueError(f"Negative values found in {population_col}")

    out = df.sort_values(rank_col, ascending=True).copy()
    total_population = out[population_col].sum()
    if total_population <= 0:
        raise ValueError("Total population must be greater than zero")

    out["cumulative_population"] = out[population_col].cumsum()
    out["cumulative_population_share"] = (
        out["cumulative_population"] / total_population
    )
    out[group_col] = np.ceil(out["cumulative_population_share"] * n_groups).astype(int)
    out[group_col] = out[group_col].clip(1, n_groups)
    return out.sort_index()


def build_datazone_assignments(n_groups: int = 5) -> pd.DataFrame:
    """Return one row per Data Zone with source/equal/population-weighted groups."""
    out = load_simd_datazone_lookup(n_groups=n_groups)
    out = assign_equal_datazone_groups(out, n_groups=n_groups)
    out = assign_population_weighted_groups(out, n_groups=n_groups)
    out["source_minus_population_weighted"] = (
        out["source_group"] - out["population_weighted_group"]
    )
    out["equal_minus_population_weighted"] = (
        out["equal_datazone_group"] - out["population_weighted_group"]
    )
    return out.sort_values("simd_rank").reset_index(drop=True)


def _summarise_one_grouping(
    assignments: pd.DataFrame,
    *,
    method: str,
    method_label: str,
    group_col: str,
    n_groups: int,
) -> pd.DataFrame:
    total_population = assignments["population"].sum()
    total_datazones = len(assignments)
    summary = (
        assignments.groupby(group_col, dropna=False)
        .agg(
            n_datazones=("datazone", "nunique"),
            total_population=("population", "sum"),
            first_simd_rank=("simd_rank", "min"),
            last_simd_rank=("simd_rank", "max"),
        )
        .rename_axis("simd_group")
        .reset_index()
    )
    summary["pct_datazones"] = 100 * summary["n_datazones"] / total_datazones
    summary["pct_population"] = 100 * summary["total_population"] / total_population
    summary.insert(0, "grouping_method", method)
    summary.insert(1, "grouping_method_label", method_label)
    summary["n_groups"] = n_groups
    return summary


def build_group_summary(
    assignments: pd.DataFrame,
    *,
    n_groups: int = 5,
) -> pd.DataFrame:
    """Summarise population and rank ranges for source/equal/weighted groups."""
    parts = [
        _summarise_one_grouping(
            assignments,
            method="source",
            method_label="Stored SIMD grouping",
            group_col="source_group",
            n_groups=n_groups,
        ),
        _summarise_one_grouping(
            assignments,
            method="equal_datazone",
            method_label="Equal Data Zone grouping",
            group_col="equal_datazone_group",
            n_groups=n_groups,
        ),
        _summarise_one_grouping(
            assignments,
            method="population_weighted",
            method_label="Population-weighted grouping",
            group_col="population_weighted_group",
            n_groups=n_groups,
        ),
    ]
    return pd.concat(parts, ignore_index=True)


def build_movement_table(assignments: pd.DataFrame) -> pd.DataFrame:
    """Cross-tabulate source/equal-Data-Zone groups against weighted groups."""
    rows: list[pd.DataFrame] = []
    comparisons = {
        "source": "source_group",
        "equal_datazone": "equal_datazone_group",
    }
    for method, source_col in comparisons.items():
        counts = (
            assignments.groupby([source_col, "population_weighted_group"])
            .agg(
                n_datazones=("datazone", "nunique"),
                total_population=("population", "sum"),
            )
            .reset_index()
            .rename(columns={source_col: "comparison_group"})
        )
        counts.insert(0, "comparison_method", method)
        counts["pct_datazones"] = 100 * counts["n_datazones"] / len(assignments)
        counts["pct_population"] = (
            100 * counts["total_population"] / assignments["population"].sum()
        )
        rows.append(counts)
    return pd.concat(rows, ignore_index=True)


def build_change_summary(assignments: pd.DataFrame) -> pd.DataFrame:
    """Summarise how many Data Zones move after population weighting."""
    rows: list[pd.DataFrame] = []
    change_cols = {
        "source": "source_minus_population_weighted",
        "equal_datazone": "equal_minus_population_weighted",
    }
    for method, change_col in change_cols.items():
        summary = (
            assignments.groupby(change_col)
            .agg(
                n_datazones=("datazone", "nunique"),
                total_population=("population", "sum"),
            )
            .rename_axis("group_difference")
            .reset_index()
        )
        summary.insert(0, "comparison_method", method)
        summary["pct_datazones"] = 100 * summary["n_datazones"] / len(assignments)
        summary["pct_population"] = (
            100 * summary["total_population"] / assignments["population"].sum()
        )
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def build_diagnostics(assignments: pd.DataFrame, group_summary: pd.DataFrame) -> pd.DataFrame:
    """Return reproducibility diagnostics for the SIMD grouping implementation."""
    source_equal_match = bool(
        assignments["source_group"].equals(assignments["equal_datazone_group"])
    )

    pop_summary = group_summary.loc[
        group_summary["grouping_method"].eq("population_weighted")
    ].copy()
    equal_summary = group_summary.loc[
        group_summary["grouping_method"].eq("equal_datazone")
    ].copy()

    raw_rows = [
        ("n_datazones", len(assignments)),
        ("total_population", int(assignments["population"].sum())),
        ("source_matches_equal_datazone", source_equal_match),
        (
            "max_abs_population_share_deviation_weighted_pct_points",
            float((pop_summary["pct_population"] - 100 / len(pop_summary)).abs().max()),
        ),
        (
            "max_abs_population_share_deviation_equal_datazone_pct_points",
            float(
                (equal_summary["pct_population"] - 100 / len(equal_summary)).abs().max()
            ),
        ),
        (
            "n_datazones_moved_source_vs_weighted",
            int(
                assignments["source_group"]
                .ne(assignments["population_weighted_group"])
                .sum()
            ),
        ),
    ]
    rows = []
    for metric, value in raw_rows:
        numeric_value = np.nan
        if isinstance(value, (int, float, np.integer, np.floating, bool, np.bool_)):
            numeric_value = float(value)
        rows.append(
            {
                "metric": metric,
                "value": str(value),
                "numeric_value": numeric_value,
            }
        )
    return pd.DataFrame(rows)


def build_simd_validation_tables(n_groups: int = 5) -> SimdValidationTables:
    """Build all SIMD population-weighting validation outputs."""
    assignments = build_datazone_assignments(n_groups=n_groups)
    group_summary = build_group_summary(assignments, n_groups=n_groups)
    return SimdValidationTables(
        datazone_assignments=assignments,
        group_summary=group_summary,
        movement_table=build_movement_table(assignments),
        change_summary=build_change_summary(assignments),
        diagnostics=build_diagnostics(assignments, group_summary),
    )
