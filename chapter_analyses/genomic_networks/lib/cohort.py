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

from utils import load_daily_policy_data, policy_order  # noqa: E402


VACCINATION_DOSE_GROUPS = (
    "Unvaccinated",
    "One dose",
    "Two doses",
    "Booster/3+ doses",
    "Vaccinated dose unknown",
)


def _attach_policy_calendar(df: pd.DataFrame, date_col: str) -> pd.DataFrame:
    """Attach policy fields by directly joining the processed daily calendar."""
    if date_col not in df.columns:
        raise KeyError(f"{date_col!r} is required for the policy-calendar join")

    lookup = load_daily_policy_data(
        ["period_code", "period_label", "policy_era", "period_order"]
    ).rename(
        columns={
            "date": "_policy_date",
            "period_code": "policy_period",
            "period_label": "policy_period_label",
        }
    )
    policy_cols = [col for col in lookup.columns if col != "_policy_date"]
    out = df.drop(columns=policy_cols, errors="ignore").copy()
    out["_policy_date"] = pd.to_datetime(out[date_col], errors="coerce").dt.normalize()
    return out.merge(
        lookup,
        on="_policy_date",
        how="left",
        validate="many_to_one",
    ).drop(columns="_policy_date")


def _policy_period_order() -> list[str]:
    """Return policy codes in the ordering stored in the policy calendar."""
    policy = load_daily_policy_data(["period_code", "period_order", "policy_era"])
    return (
        policy.drop_duplicates()
        .sort_values("period_order")["period_code"]
        .astype(str)
        .tolist()
    )


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
        cluster_columns = ["cluster_id"]
        if "cluster_size" in df.columns:
            cluster_columns.append("cluster_size")
        clusters = df[cluster_columns].drop_duplicates("cluster_id")
        rows.append(_metric_row("clusters", clusters["cluster_id"].nunique()))
        if "cluster_size" in clusters.columns:
            cluster_size = pd.to_numeric(clusters["cluster_size"], errors="coerce")
            rows.extend(
                [
                    _metric_row("singleton_clusters", cluster_size.eq(1).sum()),
                    _metric_row(
                        "non_singleton_clusters",
                        cluster_size.gt(1).sum(),
                    ),
                ]
            )
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

    out = _attach_policy_calendar(out, "wn_mid_date")
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
    group_cols: Sequence[str] | None = ("policy_period", "policy_era"),
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


def build_test_reason_by_policy_era(df: pd.DataFrame) -> pd.DataFrame:
    """Count unique sequences by test reason and epidemic era."""
    seq = sequence_level_frame(df)
    required = {"sequence_id", "test_reason", "policy_era"}
    missing = required - set(seq.columns)
    if missing:
        raise KeyError(f"Missing test-reason summary columns: {sorted(missing)}")

    work = seq[["sequence_id", "test_reason", "policy_era"]].copy()
    work["test_reason"] = work["test_reason"].astype("string").fillna("missing")
    work["policy_era"] = work["policy_era"].astype("string").fillna("missing")

    out = (
        work.groupby(["test_reason", "policy_era"], dropna=False)["sequence_id"]
        .nunique()
        .rename("n_sequences")
        .reset_index()
    )

    era_order = {era: idx for idx, era in enumerate(policy_order("policy_era"))}
    out["_policy_sort"] = out["policy_era"].map(era_order).fillna(999)
    return (
        out.sort_values(["_policy_sort", "test_reason"], kind="mergesort")
        .drop(columns="_policy_sort")
        .reset_index(drop=True)
    )


def _vaccination_context_frame(df: pd.DataFrame) -> pd.DataFrame:
    required = {
        "sequence_id",
        "is_vaccinated",
        "vacc_dose_number",
        "vacc_booster",
        "days_since_vaccination",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing vaccination context columns: {sorted(missing)}")

    out = df.copy()
    out["_vaccinated"] = pd.to_numeric(out["is_vaccinated"], errors="coerce").eq(1)
    out["_booster"] = pd.to_numeric(out["vacc_booster"], errors="coerce").eq(1)
    out["_dose_number"] = pd.to_numeric(out["vacc_dose_number"], errors="coerce")
    out["_days_since_vaccination"] = pd.to_numeric(
        out["days_since_vaccination"], errors="coerce"
    )

    out["_dose_group"] = np.select(
        [
            ~out["_vaccinated"],
            out["_vaccinated"] & (out["_booster"] | out["_dose_number"].ge(3)),
            out["_vaccinated"] & out["_dose_number"].eq(2),
            out["_vaccinated"] & out["_dose_number"].eq(1),
        ],
        [
            "Unvaccinated",
            "Booster/3+ doses",
            "Two doses",
            "One dose",
        ],
        default="Vaccinated dose unknown",
    )
    return out


def _vaccination_summary_rows(
    df: pd.DataFrame,
    group_cols: Sequence[str],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for keys, group in df.groupby(list(group_cols), dropna=False, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        n_sequences = int(group["sequence_id"].nunique())
        vaccinated = group["_vaccinated"]
        booster = group["_booster"]
        dose = group["_dose_number"]
        days = group["_days_since_vaccination"]
        valid_days = group.loc[vaccinated & days.ge(0), "_days_since_vaccination"]
        valid_dose = group.loc[vaccinated & dose.gt(0), "_dose_number"]
        dose_counts = group["_dose_group"].value_counts()

        row: dict[str, Any] = dict(zip(group_cols, keys, strict=True))
        row.update(
            {
                "n_sequences": n_sequences,
                "n_vaccinated": int(vaccinated.sum()),
                "n_unvaccinated": int(dose_counts.get("Unvaccinated", 0)),
                "n_one_dose": int(dose_counts.get("One dose", 0)),
                "n_two_doses": int(dose_counts.get("Two doses", 0)),
                "n_booster_or_three_plus": int(dose_counts.get("Booster/3+ doses", 0)),
                "n_vaccinated_dose_unknown": int(
                    dose_counts.get("Vaccinated dose unknown", 0)
                ),
                "n_booster": int(booster.sum()),
                "n_days_since_vaccination": int(valid_days.notna().sum()),
                "median_dose_number_vaccinated": valid_dose.median(),
                "median_days_since_vaccination": valid_days.median(),
                "q25_days_since_vaccination": valid_days.quantile(0.25),
                "q75_days_since_vaccination": valid_days.quantile(0.75),
            }
        )
        rows.append(row)

    out = pd.DataFrame(rows)
    for count_col in (
        "n_vaccinated",
        "n_unvaccinated",
        "n_one_dose",
        "n_two_doses",
        "n_booster_or_three_plus",
        "n_vaccinated_dose_unknown",
        "n_booster",
    ):
        prop_col = count_col.removeprefix("n_")
        out[f"prop_{prop_col}"] = out[count_col] / out["n_sequences"].replace(0, np.nan)

    out["prop_booster_among_vaccinated"] = out["n_booster"] / out[
        "n_vaccinated"
    ].replace(0, np.nan)
    return out


def build_vaccination_context_by_policy(df: pd.DataFrame) -> pd.DataFrame:
    """Summarise sequence-level vaccination context by policy period."""
    seq = sequence_level_frame(df)
    if "policy_period" not in seq.columns or "policy_era" not in seq.columns:
        raise KeyError(
            "'policy_period' and 'policy_era' are required in the persisted sequence metadata"
        )

    work = _vaccination_context_frame(seq)
    out = _vaccination_summary_rows(work, ("policy_period", "policy_era"))
    policy_order = {period: idx for idx, period in enumerate(_policy_period_order())}
    out["_policy_sort"] = out["policy_period"].astype(str).map(policy_order).fillna(999)
    return out.sort_values(["_policy_sort", "policy_period"]).drop(
        columns="_policy_sort"
    )


def build_vaccination_window_context(df: pd.DataFrame) -> pd.DataFrame:
    """Summarise rolling-window vaccination context with one row per window."""
    required = {
        "window_id",
        "window_idx",
        "wn_start_date",
        "wn_mid_date",
        "wn_end_date",
        "sequence_id",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing vaccination window columns: {sorted(missing)}")

    work = df.drop_duplicates(["window_id", "sequence_id"]).copy()
    for col in ("wn_start_date", "wn_mid_date", "wn_end_date"):
        work[col] = pd.to_datetime(work[col], errors="coerce")
    work = _attach_policy_calendar(work, "wn_mid_date")

    work = _vaccination_context_frame(work)
    group_cols = (
        "window_id",
        "window_idx",
        "wn_start_date",
        "wn_mid_date",
        "wn_end_date",
        "policy_period",
        "policy_era",
    )
    out = _vaccination_summary_rows(work, group_cols)
    return out.sort_values("window_idx").reset_index(drop=True)


def build_denominator_contrasts(window_coverage: pd.DataFrame) -> pd.DataFrame:
    """Summarise rolling-window sequence denominators by policy period."""
    required = {
        "window_id",
        "policy_period",
        "policy_era",
        "wn_no_sequences",
        "wn_positive_tests",
    }
    missing = required - set(window_coverage.columns)
    if missing:
        raise KeyError(f"Missing denominator columns: {sorted(missing)}")

    out = (
        window_coverage.groupby(["policy_period", "policy_era"], dropna=False)
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
