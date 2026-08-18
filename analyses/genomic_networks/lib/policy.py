"""Policy-period labelling and ordering helpers."""

from __future__ import annotations

import sys

import pandas as pd

from .config import PROJECT_ROOT

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import policy_descriptors


def canonicalise_policy_fields(frame: pd.DataFrame) -> pd.DataFrame:
    """Make policy labels and eras agree with the selected policy period."""
    if "policy_period" not in frame.columns:
        return frame

    descriptors = policy_descriptors().drop_duplicates("policy_period")
    lookup = descriptors.set_index(descriptors["policy_period"].astype(str))
    period = frame["policy_period"].astype(str)
    out = frame.copy()
    for source_col, target_col in (
        ("policy_era", "policy_era"),
        ("policy_period_label", "policy_period_label"),
    ):
        if source_col not in lookup.columns:
            continue
        mapped = period.map(lookup[source_col])
        if target_col in out.columns:
            out[target_col] = mapped.where(mapped.notna(), out[target_col])
        else:
            out[target_col] = mapped
    return out


def sort_by_policy_period(frame: pd.DataFrame) -> pd.DataFrame:
    """Sort rows by the canonical policy-period order when available."""
    if "policy_period" not in frame.columns:
        return frame

    descriptors = policy_descriptors().drop_duplicates("policy_period")
    order = dict(
        zip(
            descriptors["policy_period"].astype(str),
            descriptors["policy_period_order"],
            strict=True,
        )
    )
    out = frame.copy()
    out["_policy_sort"] = out["policy_period"].astype(str).map(order).fillna(999)
    sort_cols = ["_policy_sort", "policy_period"]
    if "policy_era" in out.columns:
        sort_cols.append("policy_era")
    return (
        out.sort_values(sort_cols, kind="mergesort")
        .drop(columns="_policy_sort")
        .reset_index(drop=True)
    )
