"""Build Chapter 4 LaTeX table fragments."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from assortativity_analysis import (
    compatibility_variance_decomposition_long,
    compatibility_window_pooled_meta,
    pooled_window_attribute_summary,
    variance_decomposition_summary,
)
from common import (
    ATTRIBUTE_ORDER,
    POLICY_LABELS,
    Paths,
    add_common_args,
    paths_from_args,
    read_table,
    sort_by_policy,
    window_idx_from_id,
)

from utils import (
    addlinespace_after_group_changes,
    write_latex_grouped_column_table,
    write_latex_table,
)

SIMD_GROUP_LABELS = {
    5: "quintile",
    10: "decile",
    20: "vigintile",
}
TABLE_NAMES = {
    "sequence_composition_by_policy": "tab_ch4_sequence_composition_by_policy",
    "policy_denominators": "tab_ch4_policy_denominators",
    "cluster_period_summary": "tab_ch4_cluster_period_summary",
    "cluster_pairwise_distance_summary": "tab_ch4_cluster_pairwise_distance_summary",
    "assortativity_summary": "tab_ch4_assortativity_summary",
    "variance_decomposition_summary": "tab_ch4_assortativity_variance_decomposition",
    "simd_population_weighting": "tab_ch4_simd_population_weighting",
}
ATTRIBUTE_COLUMN_LABELS = {
    "SIMD quintile": "SIMD",
    "Urban/rural class": "Urban/rural",
    "Health board": "HB",
    "Local authority": "LA",
}
SEQUENCE_COMPOSITION_ATTRIBUTES = (
    "sex",
    "simd_quintile",
    "age_group",
    "urban_rural",
    "health_board",
)
SEQUENCE_COMPOSITION_CATEGORY_LABELS = {
    "simd_quintile": {
        "1.0": "Q1",
        "2.0": "Q2",
        "3.0": "Q3",
        "4.0": "Q4",
        "5.0": "Q5",
        "1": "Q1",
        "2": "Q2",
        "3": "Q3",
        "4": "Q4",
        "5": "Q5",
    },
}


def fmt_int(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    try:
        value_float = float(value)
    except (TypeError, ValueError):
        return str(value)
    if value_float.is_integer():
        return f"{int(value_float):,}"
    return f"{value_float:,.1f}"


def fmt_percent(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{100 * float(value):.1f}%"


def fmt_percent_points(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{float(value):.1f}%"


def fmt_float(value: Any, digits: int = 2) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{float(value):.{digits}f}"


def fmt_ci(low: Any, high: Any, digits: int = 2) -> str:
    if (
        low is None
        or high is None
        or (isinstance(low, float) and np.isnan(low))
        or (isinstance(high, float) and np.isnan(high))
    ):
        return "not estimated"
    return f"{fmt_float(low, digits)}--{fmt_float(high, digits)}"


def fmt_iqr(
    median: Any,
    q25: Any,
    q75: Any,
    digits: int = 0,
) -> str:
    if median is None or (isinstance(median, float) and np.isnan(median)):
        return ""
    return (
        f"{fmt_float(median, digits)} "
        f"({fmt_float(q25, digits)}--{fmt_float(q75, digits)})"
    )


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return np.nan
    return float(np.average(values.loc[mask], weights=weights.loc[mask]))


def weighted_mean_ci_from_se(
    values: pd.Series,
    weights: pd.Series,
    standard_errors: pd.Series,
) -> dict[str, float]:
    mask = values.notna() & weights.notna() & weights.gt(0)
    if not mask.any():
        return {
            "weighted_mean": np.nan,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "ci_weight_share": np.nan,
        }

    values = values.loc[mask].astype(float)
    weights = weights.loc[mask].astype(float)
    weighted_mean_value = float(np.average(values, weights=weights))

    se_mask = standard_errors.loc[mask].notna()
    if not se_mask.any():
        return {
            "weighted_mean": weighted_mean_value,
            "combined_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "ci_weight_share": np.nan,
        }

    ci_weights = weights.loc[se_mask]
    ci_standard_errors = standard_errors.loc[mask].loc[se_mask].astype(float)
    normalized = ci_weights / weights.sum()
    combined_se = float(np.sqrt(np.sum((normalized * ci_standard_errors) ** 2)))
    return {
        "weighted_mean": weighted_mean_value,
        "combined_se": combined_se,
        "ci_low": weighted_mean_value - 1.96 * combined_se,
        "ci_high": weighted_mean_value + 1.96 * combined_se,
        "ci_weight_share": float(ci_weights.sum() / weights.sum()),
    }


def compatibility_assortativity_filtered(paths: Paths) -> pd.DataFrame:
    assort = read_table(paths, "compatibility_assortativity")
    assort["window_idx"] = window_idx_from_id(assort["window_id"])
    return assort.loc[
        assort["assortativity"].notna()
        & assort["edge_weight_total"].gt(0)
        & assort["n_categories"].gt(1)
        & assort["n_edge_contributions_used"].ge(20)
    ].copy()


def compatibility_window_assortativity(paths: Paths) -> pd.DataFrame:
    work = compatibility_assortativity_filtered(paths)
    rows = []
    for (window_idx, attribute, label), group in work.groupby(
        ["window_idx", "attribute", "attribute_label"], dropna=False
    ):
        ci = weighted_mean_ci_from_se(
            group["assortativity"],
            group["edge_weight_total"],
            group["assortativity_se"],
        )
        rows.append(
            {
                "window_idx": window_idx,
                "attribute": attribute,
                "attribute_label": label,
                "assortativity": ci["weighted_mean"],
                "assortativity_se": ci["combined_se"],
                "assortativity_ci_low": ci["ci_low"],
                "assortativity_ci_high": ci["ci_high"],
                "ci_weight_share": ci["ci_weight_share"],
                "edge_weight_total": group["edge_weight_total"].sum(),
                "eligible_networks": group["pairwise_stem"].nunique()
                if "pairwise_stem" in group.columns
                else len(group),
            }
        )
    return pd.DataFrame(rows)


def compatibility_attribute_summary(paths: Paths) -> pd.DataFrame:
    work = compatibility_assortativity_filtered(paths)
    window_summary = compatibility_window_assortativity(paths)
    rows = []
    for (attribute, label), group in work.groupby(
        ["attribute", "attribute_label"], dropna=False
    ):
        ci = weighted_mean_ci_from_se(
            group["assortativity"],
            group["edge_weight_total"],
            group["assortativity_se"],
        )
        window_group = window_summary.loc[
            window_summary["attribute"].eq(str(attribute))
        ]
        rows.append(
            {
                "attribute": attribute,
                "attribute_label": label,
                "n_windows": window_group["window_idx"].nunique(),
                "n_networks": group["pairwise_stem"].nunique()
                if "pairwise_stem" in group.columns
                else len(group),
                "weighted_mean": ci["weighted_mean"],
                "combined_se": ci["combined_se"],
                "ci_low": ci["ci_low"],
                "ci_high": ci["ci_high"],
                "ci_weight_share": ci["ci_weight_share"],
                "window_median": window_group["assortativity"].median(),
                "window_q10": window_group["assortativity"].quantile(0.10),
                "window_q90": window_group["assortativity"].quantile(0.90),
            }
        )
    return pd.DataFrame(rows)


def pooled_window_attribute_summary_table(paths: Paths) -> pd.DataFrame:
    try:
        summary = read_table(paths, "compatibility_window_pooled_summary")
        if {"gls_mean", "gls_ci_low", "gls_ci_high"}.issubset(summary.columns):
            return summary
    except FileNotFoundError:
        pass
    window_meta, _ = compatibility_window_pooled_meta(paths)
    return pooled_window_attribute_summary(window_meta)


def variance_decomposition_summary_table(paths: Paths) -> pd.DataFrame:
    try:
        return read_table(paths, "compatibility_variance_decomposition_summary")
    except FileNotFoundError:
        vd_long = compatibility_variance_decomposition_long(paths)
        return variance_decomposition_summary(vd_long)


def sequence_composition_category_order(
    table: pd.DataFrame,
    attribute: str,
) -> list[str]:
    work = table.loc[table["attribute"].eq(attribute)]
    totals = work.groupby("category", observed=False)["n_sequences"].sum()
    if totals.empty:
        return []

    if attribute == "sex":
        preferred = ["Female", "Male"]
        return [value for value in preferred if value in totals.index.astype(str)]

    if attribute == "simd_quintile":
        return sorted(totals.index.astype(str), key=lambda value: float(value))

    if attribute == "age_group":
        preferred = ["00-04", "05-14", "15-24", "25-64", "65-74", "75+"]
        return [value for value in preferred if value in totals.index.astype(str)]

    if attribute == "urban_rural":
        preferred = [
            "Large Urban Areas",
            "Other Urban Areas",
            "Accessible Small Towns",
            "Remote Small Towns",
            "Accessible Rural",
            "Remote Rural",
        ]
        return [value for value in preferred if value in totals.index.astype(str)]

    if attribute == "health_board":
        return totals.sort_values(ascending=False).index.astype(str).tolist()

    return totals.sort_values(ascending=False).index.astype(str).tolist()


def display_sequence_category(attribute: str, category: object) -> str:
    text = str(category)
    return SEQUENCE_COMPOSITION_CATEGORY_LABELS.get(attribute, {}).get(text, text)


def format_sequence_composition_cell(row: pd.Series | None) -> str:
    if row is None:
        return "0 (0.0%)"
    n_sequences = pd.to_numeric(pd.Series([row["n_sequences"]]), errors="coerce").iloc[
        0
    ]
    if pd.notna(n_sequences) and float(n_sequences) == 0:
        return "0 (0.0%)"
    if bool(row.get("small_cell", False)):
        return "<5"
    return f"{fmt_int(n_sequences)} ({fmt_percent(row['proportion'])})"


def sequence_composition_policy_groups(
    composition: pd.DataFrame,
) -> tuple[pd.DataFrame, list[tuple[str, list[str]]]]:
    periods = composition[["policy_era", "policy_period"]].drop_duplicates()
    periods = sort_by_policy(periods, column="policy_period")

    column_groups: list[tuple[str, list[str]]] = []
    for era, group in periods.groupby("policy_era", sort=False, dropna=False):
        label = POLICY_LABELS.get(str(era), str(era).capitalize().replace("_", " "))
        column_groups.append((label, group["policy_period"].astype(str).tolist()))

    return periods, column_groups


def write_sequence_composition_by_policy_table(paths: Paths) -> None:
    composition = read_table(paths, "sequence_composition_by_policy")
    composition = composition.loc[
        composition["attribute"].isin(SEQUENCE_COMPOSITION_ATTRIBUTES)
    ].copy()
    if composition.empty:
        raise FileNotFoundError(
            "sequence_composition_by_policy contains no figure attributes"
        )

    periods, column_groups = sequence_composition_policy_groups(composition)
    period_codes = periods["policy_period"].astype(str).tolist()

    lookup = {
        (str(row["attribute"]), str(row["category"]), str(row["policy_period"])): row
        for row in composition.to_dict(orient="records")
    }

    rows: list[list[Any]] = []
    addlinespace_after: set[int] = set()
    for attribute_idx, attribute in enumerate(SEQUENCE_COMPOSITION_ATTRIBUTES):
        attr_rows = composition.loc[composition["attribute"].eq(attribute)]
        if attr_rows.empty:
            continue
        attribute_label = str(attr_rows["attribute_label"].dropna().iloc[0])
        categories = sequence_composition_category_order(composition, attribute)
        for category_idx, category in enumerate(categories):
            rows.append(
                [
                    attribute_label if category_idx == 0 else "",
                    display_sequence_category(attribute, category),
                    *[
                        format_sequence_composition_cell(
                            pd.Series(lookup.get((attribute, str(category), period)))
                            if lookup.get((attribute, str(category), period))
                            is not None
                            else None
                        )
                        for period in period_codes
                    ],
                ]
            )
        if attribute_idx < len(SEQUENCE_COMPOSITION_ATTRIBUTES) - 1 and rows:
            addlinespace_after.add(len(rows) - 1)

    write_latex_grouped_column_table(
        paths,
        TABLE_NAMES["sequence_composition_by_policy"],
        caption=(
            "Sequence composition by epidemic era and policy period for demographic, socioeconomic, or geographic variables. "
            "Columns are grouped by epidemic era, with policy-period codes shown as subcolumns. "
            "Cells report the number and percentage of unique sequences in that policy "
            "period with the corresponding attribute category. Counts from one to four are "
            "suppressed and shown as <5."
        ),
        short_caption="Sequence composition by epidemic era and policy period",
        label=TABLE_NAMES["sequence_composition_by_policy"],
        row_columns=["Attribute", "Category"],
        column_groups=column_groups,
        rows=rows,
        column_spec=f"ll*{{{len(period_codes)}}}{{r}}",
        addlinespace_after=addlinespace_after,
        landscape=True,
    )


def write_policy_denominator_table(paths: Paths) -> None:
    denominators = read_table(paths, "window_denominator_contrasts")
    has_sequences = pd.to_numeric(
        denominators["median_window_sequences"], errors="coerce"
    ).gt(0)
    denominators = denominators.loc[
        has_sequences | denominators["policy_period"].astype(str).eq("P2")
    ]
    denominators = sort_by_policy(denominators, column="policy_period")
    rows = []
    for row in denominators.itertuples(index=False):
        rows.append(
            [
                str(row.policy_era).capitalize().replace("_", " "),
                row.policy_period,
                fmt_int(row.n_windows),
                fmt_int(row.median_window_sequences),
                fmt_int(row.median_window_positive_tests),
                f"{fmt_percent(row.median_window_prop_sequenced)} ({fmt_percent(row.min_window_prop_sequenced)}--{fmt_percent(row.max_window_prop_sequenced)})",
            ]
        )
    write_latex_table(
        paths,
        TABLE_NAMES["policy_denominators"],
        caption=(
            "Rolling-window observation denominators stratified by epidemic era and "
            "policy period. For each period, the table outlines the total number of "
            "analytical windows alongside the median window-level counts of sequenced "
            "genomes and positive PCR tests. Genomic sampling coverage (the percentage of "
            "positive tests successfully sequenced) is reported as a median, with "
            "absolute minimum and maximum extremes presented in parentheses."
        ),
        short_caption="Rolling-window observation denominators by policy period",
        label=TABLE_NAMES["policy_denominators"],
        columns=[
            "Epidemic era",
            "Period",
            "Windows",
            "Genomes",
            "Tests",
            "Coverage",
        ],
        rows=rows,
        column_spec="llrrrl",
        addlinespace_after=addlinespace_after_group_changes(denominators["policy_era"]),
    )


def write_cluster_period_table(paths: Paths) -> None:
    clusters = read_table(paths, "cluster_period_summary")
    clusters = sort_by_policy(clusters, column="policy_period")
    rows = []
    for row in clusters.itertuples(index=False):
        cells = [
            str(row.policy_era).capitalize().replace("_", " "),
            row.policy_period,
            fmt_int(row.n_sequence_memberships),
            f"{fmt_int(row.n_clusters)} ({fmt_int(row.n_non_singleton_clusters)})",
        ]
        if fmt_int(row.n_non_singleton_clusters) != "0":
            cells.append(
                f"{fmt_float(row.median_non_singleton_cluster_size, 0)}; "
                f"{fmt_float(row.p90_non_singleton_cluster_size, 0)}; "
                f"{fmt_int(row.max_non_singleton_cluster_size)}"
            )
            cells.append(
                f"{fmt_float(row.median_non_singleton_datazones, 0)}; "
                f"{fmt_float(row.p90_non_singleton_datazones, 0)}; "
                f"{fmt_int(row.max_non_singleton_datazones)}"
            )
            cells.append(
                fmt_iqr(
                    row.median_non_singleton_spatial_distance_km,
                    row.q25_non_singleton_spatial_distance_km,
                    row.q75_non_singleton_spatial_distance_km,
                    digits=1,
                )
            )
            cells.append(
                fmt_iqr(
                    row.median_non_singleton_duration_days,
                    row.q25_non_singleton_duration_days,
                    row.q75_non_singleton_duration_days,
                    digits=1,
                )
            )
        else:
            cells.extend(["NA"] * 4)
        rows.append(cells)
    write_latex_table(
        paths,
        TABLE_NAMES["cluster_period_summary"],
        caption=(
            "Summary of Scottish EpiLink clusters stratified by epidemic era and policy "
            "period. The table reports sequence-window memberships and total distinct clusters, "
            "with non-singleton cluster counts provided in parentheses. Further metrics "
            "for non-singleton clusters summarise distributions (median; 90th percentile; "
            "and maximum) for both cluster size and affected Data Zones. Residential reach "
            "(in kilometers) and temporal span (in days) are reported as medians alongside "
            "their corresponding interquartile ranges."
        ),
        short_caption="Scottish EpiLink cluster summary by policy period",
        label=TABLE_NAMES["cluster_period_summary"],
        columns=[
            "Epidemic era",
            "Period",
            "Members",
            "Clusters",
            "Size",
            "Zones",
            "Reach (km)",
            "Span (d)",
        ],
        rows=rows,
        addlinespace_after=addlinespace_after_group_changes(clusters["policy_era"]),
        column_spec="llrrrrrr",
    )


def write_cluster_pairwise_distance_summary_table(paths: Paths) -> None:
    summary = read_table(paths, "cluster_pairwise_distance_summary")
    metric_labels = {
        "snp_distance": "Genetic (SNPs)",
        "temporal_distance": "Temporal (days)",
    }
    metric_order = {name: idx for idx, name in enumerate(metric_labels)}
    weighting_labels = {
        "unweighted": "Unweighted",
        "pair_count_weighted": "Pair-weighted",
    }
    weighting_order = {name: idx for idx, name in enumerate(weighting_labels)}
    summary = summary.copy()
    summary["_metric_sort"] = summary["distance_metric"].map(metric_order).fillna(999)
    summary["_weighting_sort"] = summary["weighting"].map(weighting_order).fillna(999)
    summary = summary.sort_values(
        ["_metric_sort", "_weighting_sort", "distance_metric", "weighting"],
        kind="mergesort",
    )

    rows = []
    for row in summary.itertuples(index=False):
        rows.append(
            [
                metric_labels.get(str(row.distance_metric), str(row.distance_metric)),
                weighting_labels.get(str(row.weighting), str(row.weighting)),
                fmt_int(row.n_window_lineage_summaries),
                fmt_int(row.n_windows),
                fmt_int(row.n_pairwise_rows),
                fmt_iqr(row.median, row.q25, row.q75, digits=1),
            ]
        )

    write_latex_table(
        paths,
        TABLE_NAMES["cluster_pairwise_distance_summary"],
        caption=(
            "Distributions of within-cluster genetic (SNPs) and temporal (days) "
            "distances. For each metric, the table contrasts unweighted and "
            "pair-weighted quartiles of window-lineage median pairwise distances. "
            "Analyses are restricted to summaries containing a minimum of ten observed "
            "within-cluster pairs. Median pairwise distances are reported alongside "
            "corresponding interquartile ranges enclosed in parentheses."
        ),
        short_caption="Within-cluster pairwise distance summaries",
        label=TABLE_NAMES["cluster_pairwise_distance_summary"],
        columns=[
            "Metric",
            "Weighting",
            "Summaries",
            "Windows",
            "Pairs",
            "Median",
        ],
        rows=rows,
        column_spec="llrrrl",
    )


def write_assortativity_summary_table(paths: Paths) -> None:
    summary = pooled_window_attribute_summary_table(paths)
    rows = []
    for attribute in ATTRIBUTE_ORDER:
        group = summary.loc[summary["attribute_label"].eq(attribute)]
        if group.empty:
            continue
        row = group.iloc[0]
        mean_col = "gls_mean" if "gls_mean" in row.index else "weighted_mean"
        low_col = "gls_ci_low" if "gls_ci_low" in row.index else "ci_low"
        high_col = "gls_ci_high" if "gls_ci_high" in row.index else "ci_high"
        rows.append(
            [
                attribute,
                fmt_int(row["n_windows"]),
                fmt_int(row["n_estimated_windows"]),
                fmt_float(row["median_n_lineages"], 1),
                (
                    f"{fmt_float(row[mean_col], 4)} "
                    f"({fmt_ci(row[low_col], row[high_col], 4)})"
                ),
                (
                    f"{fmt_float(row['window_median'], 3)} "
                    f"({fmt_float(row['window_q10'], 3)}--"
                    f"{fmt_float(row['window_q90'], 3)})"
                ),
            ]
        )
    write_latex_table(
        paths,
        TABLE_NAMES["assortativity_summary"],
        caption=(
            "Pooled compatibility assortativity estimates calculated across rolling windows. "
            "For each attribute, the table reports the total number of windows, the subset "
            "yielding an estimable coefficient, and the median number of lineages analyzed. "
            "Mean r is the single overall pooled estimate: the GLS mean of the window-level "
            "pooled assortativity estimates, adjusted for covariance induced by overlapping "
            "rolling windows, with its 95% confidence interval. Window r is descriptive: "
            "the median of the same window-level pooled estimates, followed by their "
            "10th--90th percentile range."
        ),
        short_caption="Pooled compatibility assortativity estimates",
        label=TABLE_NAMES["assortativity_summary"],
        columns=[
            "Attribute",
            "Windows",
            "Est.",
            "Lineages",
            "Mean r (95% CI)",
            "Window r (10th--90th)",
        ],
        rows=rows,
        column_spec="lrrrrll",
    )


def write_variance_decomposition_summary_table(paths: Paths) -> None:
    summary = variance_decomposition_summary_table(paths)
    attribute_rows = {}
    for attribute in ATTRIBUTE_ORDER:
        group = summary.loc[summary["attribute_label"].eq(attribute)]
        if not group.empty:
            attribute_rows[attribute] = group.iloc[0]

    metrics: list[tuple[str, str, Callable[[Any], str]]] = [
        ("Rows", "n", fmt_int),
        ("Windows", "n_windows", fmt_int),
        ("Lineages", "n_lineages", fmt_int),
        (
            "Additive model",
            "additive_model_fraction",
            lambda value: fmt_float(value, 3),
        ),
        (
            "Adj. additive",
            "adj_additive_model_fraction",
            lambda value: fmt_float(value, 3),
        ),
        ("Residual", "residual_fraction", lambda value: fmt_float(value, 3)),
        (
            "Window alone",
            "window_alone_fraction",
            lambda value: fmt_float(value, 3),
        ),
        (
            "Lineage alone",
            "lineage_alone_fraction",
            lambda value: fmt_float(value, 3),
        ),
        (
            "Lineage | Window",
            "lineage_given_window_fraction",
            lambda value: fmt_float(value, 3),
        ),
        (
            "Window | Lineage",
            "window_given_lineage_fraction",
            lambda value: fmt_float(value, 3),
        ),
        ("Weight cap", "weight_cap", lambda value: fmt_float(value, 2)),
        ("Capped weights", "n_weights_capped", fmt_int),
        ("Median boot SE", "median_boot_se", lambda value: fmt_float(value, 3)),
        (
            "Median CI width",
            "median_ci_width",
            lambda value: fmt_float(value, 3),
        ),
    ]

    attributes = list(attribute_rows)
    rows = []
    for metric_label, column, formatter in metrics:
        rows.append(
            [
                metric_label,
                *[
                    formatter(attribute_rows[attribute][column])
                    for attribute in attributes
                ],
            ]
        )

    write_latex_table(
        paths,
        TABLE_NAMES["variance_decomposition_summary"],
        caption=(
            "Variance decomposition of weighted compatibility assortativity estimates "
            "following the application of a 90th percentile inverse-variance weight cap. "
            "Target attributes are evaluated in columns using abbreviated labels (SIMD "
            "denotes SIMD quintile; HB denotes Health board; LA denotes Local authority). "
            "Variance components and fractional metrics are presented as unitless "
            "proportions, whereas underlying count and measurement uncertainty summaries "
            "are reported in their native units."
        ),
        short_caption="Variance decomposition of compatibility assortativity",
        label=TABLE_NAMES["variance_decomposition_summary"],
        columns=[
            "Metric",
            *[
                ATTRIBUTE_COLUMN_LABELS.get(attribute) or attribute
                for attribute in attributes
            ],
        ],
        rows=rows,
        column_spec="l" + "r" * len(attributes),
        addlinespace_after={2, 9},
    )


def write_simd_population_weighting_table(paths: Paths) -> None:
    """Write the SIMD validation LaTeX appendix table."""
    group_summary = read_table(paths, "simd_population_weighting_group_summary")
    n_groups = 5
    group_name = SIMD_GROUP_LABELS.get(n_groups, f"{n_groups}-group")
    caption = (
        f"Validation of national SIMD {group_name} boundaries comparing equal-Data-Zone "
        "and population-weighted grouping methodologies. The table details the number of "
        "Data Zones, total population counts, proportional population shares, and the "
        "underlying SIMD rank ranges represented within each constructed group."
    )
    short_caption = f"SIMD {group_name} population-weighting validation"

    methods = ["equal_datazone", "population_weighted"]
    display = group_summary.loc[group_summary["grouping_method"].isin(methods)].copy()
    method_order = {method: idx for idx, method in enumerate(methods)}
    display["_method_sort"] = display["grouping_method"].map(method_order)
    display = display.sort_values(["_method_sort", "simd_group"])
    display["rank_range"] = (
        display["first_simd_rank"].map(fmt_int)
        + "--"
        + display["last_simd_rank"].map(fmt_int)
    )
    rows = []
    addlinespace_after: set[int] = set()
    method_groups = list(display.groupby("_method_sort", sort=False))
    for method_idx, (_, group) in enumerate(method_groups):
        group = group.sort_values("simd_group")
        for row_idx, row in enumerate(group.itertuples(index=False)):
            rows.append(
                [
                    str(row.grouping_method_label) if row_idx == 0 else "",
                    row.simd_group,
                    fmt_int(row.n_datazones),
                    fmt_int(row.total_population),
                    fmt_percent_points(row.pct_population),
                    str(row.rank_range),
                ]
            )
        if method_idx < len(method_groups) - 1 and rows:
            addlinespace_after.add(len(rows) - 1)

    write_latex_table(
        paths,
        TABLE_NAMES["simd_population_weighting"],
        caption=caption,
        short_caption=short_caption,
        label=TABLE_NAMES["simd_population_weighting"],
        columns=[
            "Method",
            "Group",
            "Zones",
            "Pop.",
            "Share",
            "Ranks",
        ],
        rows=rows,
        column_spec="lrrrrl",
        addlinespace_after=addlinespace_after,
    )


TABLE_WRITERS: tuple[tuple[str, Callable[[Paths], object]], ...] = (
    (
        TABLE_NAMES["sequence_composition_by_policy"],
        write_sequence_composition_by_policy_table,
    ),
    (TABLE_NAMES["policy_denominators"], write_policy_denominator_table),
    (TABLE_NAMES["cluster_period_summary"], write_cluster_period_table),
    (
        TABLE_NAMES["cluster_pairwise_distance_summary"],
        write_cluster_pairwise_distance_summary_table,
    ),
    (TABLE_NAMES["assortativity_summary"], write_assortativity_summary_table),
    (
        TABLE_NAMES["variance_decomposition_summary"],
        write_variance_decomposition_summary_table,
    ),
    (TABLE_NAMES["simd_population_weighting"], write_simd_population_weighting_table),
)


def write_tables(paths: Paths) -> None:
    paths.figure_dir.mkdir(parents=True, exist_ok=True)
    skipped: list[str] = []
    for _, write_table in TABLE_WRITERS:
        try:
            write_table(paths)
        except FileNotFoundError as err:
            skipped.append(str(err))
            print(f"Skipping table writer: {err}")
    if skipped:
        print(
            f"Skipped {len(skipped)} table writer(s) because required source tables were missing."
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    write_tables(paths)
    print(f"Wrote Chapter 4 LaTeX tables to {paths.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
