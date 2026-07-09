"""Observed cohort, denominator, and coverage summaries for Chapter 4."""

from __future__ import annotations

from typing import Iterable, Sequence, Any
import sys

import numpy as np
import pandas as pd

from .config import (
    DEFAULT_MIXING_ATTRIBUTES,
    DISCLOSURE_MIN_CELL,
    PROJECT_ROOT,
    AttributeSpec,
)


if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import attach_period  # noqa: E402


def sequence_level_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Return one row per sequence from sequence-window input."""
    if "sequence_id" not in df.columns:
        raise KeyError("'sequence_id' is required")

    sort_cols = [col for col in ("collection_date", "window_idx") if col in df.columns]
    work = df.sort_values(sort_cols) if sort_cols else df
    return work.drop_duplicates("sequence_id").reset_index(drop=True)


def _metric_row(metric: str, value: Any) -> dict[str, Any]:
    numeric_value = np.nan
    if isinstance(value, (int, float, np.integer, np.floating)) and not pd.isna(value):
        numeric_value = float(value)
    return {
        "metric": metric,
        "value": "" if pd.isna(value) else str(value),
        "numeric_value": numeric_value,
    }


def build_cohort_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Build high-level counts for the observed sequenced record."""
    seq = sequence_level_frame(df)

    rows = [
        _metric_row("sequence_window_rows", len(df)),
        _metric_row("unique_sequences", seq["sequence_id"].nunique()),
    ]

    if "patient_id" in seq.columns:
        rows.append(_metric_row("unique_patients", seq["patient_id"].nunique()))
    if "collection_date" in seq.columns:
        dates = pd.to_datetime(seq["collection_date"], errors="coerce")
        rows.extend(
            [
                _metric_row("first_collection_date", dates.min().date()),
                _metric_row("last_collection_date", dates.max().date()),
            ]
        )
    if "window_id" in df.columns:
        rows.append(_metric_row("windows", df["window_id"].nunique()))
    if "cluster_id" in df.columns:
        rows.append(_metric_row("clusters", df["cluster_id"].nunique()))
    if "clade" in seq.columns:
        rows.append(_metric_row("clades", seq["clade"].nunique()))
    if "pango_lineage" in seq.columns:
        rows.append(_metric_row("pango_lineages", seq["pango_lineage"].nunique()))

    return pd.DataFrame(rows)


def build_window_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """Return one row per rolling window with sequence and denominator counts."""
    required = {
        "window_id",
        "window_idx",
        "wn_start_date",
        "wn_mid_date",
        "wn_end_date",
        "wn_no_sequences",
        "wn_positive_tests",
        "wn_prop_sequenced",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing window coverage columns: {sorted(missing)}")

    out = (
        df[list(required)]
        .drop_duplicates(["window_id", "window_idx"])
        .sort_values("window_idx")
        .reset_index(drop=True)
    )
    for col in ("wn_start_date", "wn_mid_date", "wn_end_date"):
        out[col] = pd.to_datetime(out[col], errors="coerce")

    out = attach_period(out, "wn_mid_date")
    out["sequences_per_positive_test"] = out["wn_no_sequences"] / out[
        "wn_positive_tests"
    ].replace(0, np.nan)
    return out


def build_clade_window_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Count observed sequences by window and clade."""
    required = {"window_id", "window_idx", "sequence_id", "clade"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing clade-window columns: {sorted(missing)}")

    work = df[["window_id", "window_idx", "sequence_id", "clade"]].drop_duplicates()
    counts = (
        work.groupby(["window_id", "window_idx", "clade"], dropna=False)["sequence_id"]
        .nunique()
        .rename("n_sequences")
        .reset_index()
    )
    totals = (
        counts.groupby(["window_id", "window_idx"], dropna=False)["n_sequences"]
        .sum()
        .rename("window_sequences")
        .reset_index()
    )
    out = counts.merge(totals, on=["window_id", "window_idx"], how="left")
    out["proportion"] = out["n_sequences"] / out["window_sequences"].replace(0, np.nan)
    return out.sort_values(["window_idx", "clade"]).reset_index(drop=True)


def _normalise_group_cols(group_cols: Sequence[str] | None) -> list[str]:
    return [col for col in (group_cols or []) if col]


def build_sequence_composition(
    df: pd.DataFrame,
    *,
    attributes: Iterable[AttributeSpec] = DEFAULT_MIXING_ATTRIBUTES,
    group_cols: Sequence[str] | None = ("policy_period",),
    min_cell: int = DISCLOSURE_MIN_CELL,
) -> pd.DataFrame:
    """Build long sequence-composition tables for selected categorical variables."""
    seq = sequence_level_frame(df)
    group_cols = _normalise_group_cols(group_cols)
    missing_groups = [col for col in group_cols if col not in seq.columns]
    if missing_groups:
        raise KeyError(f"Missing grouping columns: {missing_groups}")

    rows: list[pd.DataFrame] = []
    for spec in attributes:
        if spec.column not in seq.columns:
            continue

        cols = [*group_cols, spec.column, "sequence_id"]
        work = seq[cols].copy()
        work[spec.column] = work[spec.column].astype("string").fillna("Missing")

        if group_cols:
            counts = (
                work.groupby([*group_cols, spec.column], dropna=False)["sequence_id"]
                .nunique()
                .rename("n_sequences")
                .reset_index()
            )
            totals = (
                counts.groupby(group_cols, dropna=False)["n_sequences"]
                .sum()
                .rename("group_sequences")
                .reset_index()
            )
            counts = counts.merge(totals, on=group_cols, how="left")
        else:
            counts = (
                work.groupby(spec.column, dropna=False)["sequence_id"]
                .nunique()
                .rename("n_sequences")
                .reset_index()
            )
            counts["group_sequences"] = counts["n_sequences"].sum()

        counts = counts.rename(columns={spec.column: "category"})
        counts.insert(0, "attribute", spec.name)
        counts.insert(1, "attribute_label", spec.label)
        counts["proportion"] = counts["n_sequences"] / counts[
            "group_sequences"
        ].replace(0, np.nan)
        counts["small_cell"] = counts["n_sequences"].between(1, min_cell - 1)
        rows.append(counts)

    if not rows:
        return pd.DataFrame(
            columns=[
                "attribute",
                "attribute_label",
                "category",
                *group_cols,
                "n_sequences",
                "group_sequences",
                "proportion",
                "small_cell",
            ]
        )

    return pd.concat(rows, ignore_index=True, sort=False)


def build_denominator_contrasts(window_coverage: pd.DataFrame) -> pd.DataFrame:
    """Summarise rolling-window sequence denominators by policy period."""
    required = {
        "window_id",
        "policy_period",
        "wn_no_sequences",
        "wn_positive_tests",
    }
    missing = required - set(window_coverage.columns)
    if missing:
        raise KeyError(f"Missing denominator columns: {sorted(missing)}")

    out = (
        window_coverage.groupby("policy_period", dropna=False)
        .agg(
            n_windows=("window_id", "nunique"),
            median_window_sequences=("wn_no_sequences", "median"),
            median_window_positive_tests=("wn_positive_tests", "median"),
            median_window_prop_sequenced=("wn_prop_sequenced", "median"),
            min_window_prop_sequenced=("wn_prop_sequenced", "min"),
            max_window_prop_sequenced=("wn_prop_sequenced", "max"),
        )
        .reset_index()
    )
    return out
