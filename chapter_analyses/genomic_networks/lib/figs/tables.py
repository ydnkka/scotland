"""Build Chapter 4 observation-network LaTeX table fragments."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys
from typing import Callable, Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

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
from chapter_analyses.genomic_networks.lib.config import TABLES_DIR  # noqa: E402
from assortativity_analysis import (  # noqa: E402
    compatibility_variance_decomposition_long,
    compatibility_window_pooled_meta,
    pooled_window_attribute_summary,
    variance_decomposition_summary,
)


SIMD_GROUP_LABELS = {
    5: "quintile",
    10: "decile",
    20: "vigintile",
}
TABLE_NAMES = {
    "cohort_objects": "tab_ch4_cohort_objects",
    "policy_denominators": "tab_ch4_policy_denominators",
    "vaccination_context": "tab_ch4_vaccination_context",
    "test_reason_by_policy_era": "tab_ch4_test_reason_by_policy_era",
    "cluster_period_summary": "tab_ch4_cluster_period_summary",
    "assortativity_summary": "tab_ch4_assortativity_summary",
    "variance_decomposition_summary": "tab_ch4_assortativity_variance_decomposition",
    "simd_population_weighting": "tab_ch4_simd_population_weighting",
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


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def latex_column_spec(column_spec: str | None, n_columns: int) -> str:
    spec = (column_spec or ("l" * n_columns)).strip()
    if spec.startswith("@{"):
        return spec
    return f"@{{}}{spec}@{{}}"


def render_latex_table(
    *,
    caption: str,
    label: str,
    columns: list[str],
    rows: list[list[object]],
    column_spec: str | None = None,
    small: bool = True,
    addlinespace_after: set[int] | None = None,
    tabcolsep: str = "4pt",
    arraystretch: str = "1.12",
) -> str:
    column_spec = latex_column_spec(column_spec, len(columns))
    addlinespace_after = addlinespace_after or set()

    header = " & ".join(f"\\textbf{{{latex_escape(col)}}}" for col in columns)
    body_lines = []
    for row_idx, row in enumerate(rows):
        if len(row) != len(columns):
            raise ValueError(
                f"Table row {row_idx} has {len(row)} cells; expected {len(columns)}."
            )
        body_lines.append(
            "    " + " & ".join(latex_escape(cell) for cell in row) + r" \\"
        )
        if row_idx in addlinespace_after and row_idx < len(rows) - 1:
            body_lines.append(r"    \addlinespace[0.35em]")

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        f"\\caption[{latex_escape(caption)}]{{{latex_escape(caption)}}}\\label{{{label}}}",
        r"\begingroup",
    ]
    if small:
        lines.append(r"\small")
    lines.extend(
        [
            f"\\setlength{{\\tabcolsep}}{{{tabcolsep}}}",
            f"\\renewcommand{{\\arraystretch}}{{{arraystretch}}}",
            r"\begin{adjustbox}{max width=\textwidth,center}",
            f"\\begin{{tabular}}{{{column_spec}}}",
            r"\toprule",
            f"{header} \\\\",
            r"\midrule",
            *body_lines,
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{adjustbox}",
            r"\endgroup",
            r"\end{table}",
        ]
    )
    return "\n".join(lines)


def write_latex_table(
    paths: Paths,
    name: str,
    *,
    caption: str,
    label: str,
    columns: list[str],
    rows: list[list[Any]],
    column_spec: str | None = None,
    small: bool = True,
    addlinespace_after: set[int] | None = None,
) -> None:
    paths.figure_dir.mkdir(parents=True, exist_ok=True)
    content = render_latex_table(
        caption=caption,
        label=label,
        columns=columns,
        rows=rows,
        column_spec=column_spec,
        small=small,
        addlinespace_after=addlinespace_after,
    )
    (paths.figure_dir / f"{name}.tex").write_text(content + "\n")


def simd_appendix_table_tex(
    group_summary: pd.DataFrame,
    *,
    n_groups: int = 5,
    caption: str | None = None,
    label: str = TABLE_NAMES["simd_population_weighting"],
) -> str:
    """Render a compact LaTeX table for the SIMD validation appendix."""
    group_name = SIMD_GROUP_LABELS.get(n_groups, f"{n_groups}-group")
    caption = caption or (
        f"Validation of national population-weighted SIMD {group_name} groupings."
    )

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

    return render_latex_table(
        caption=caption,
        label=label,
        columns=[
            "Grouping",
            "SIMD group",
            "Data zones",
            "Population",
            "Population share",
            "SIMD rank range",
        ],
        rows=rows,
        column_spec="lrrrrl",
        addlinespace_after=addlinespace_after,
    )


def write_simd_appendix_table(
    group_summary: pd.DataFrame,
    *,
    n_groups: int = 5,
    path: Path | None = None,
) -> Path:
    """Write the SIMD validation LaTeX appendix table and return its path."""
    path = path or TABLES_DIR / f"{TABLE_NAMES['simd_population_weighting']}.tex"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(simd_appendix_table_tex(group_summary, n_groups=n_groups))
    return path


def write_simd_population_weighting_table(paths: Paths) -> Path:
    group_summary = read_table(paths, "simd_population_weighting_group_summary")
    return write_simd_appendix_table(
        group_summary,
        path=paths.figure_dir / f"{TABLE_NAMES['simd_population_weighting']}.tex",
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


def write_cohort_table(paths: Paths) -> None:
    cohort = read_table(paths, "cohort_summary")
    cohort_map = cohort.set_index("metric")["value"].to_dict()
    rows = [
        ["Unique sequences", fmt_int(cohort_map.get("unique_sequences"))],
        ["Unique patients", fmt_int(cohort_map.get("unique_patients"))],
        [
            "Collection dates",
            f"{cohort_map.get('first_collection_date')} to "
            f"{cohort_map.get('last_collection_date')}",
        ],
        ["Rolling windows", fmt_int(cohort_map.get("windows"))],
        ["Window-level clusters (all)", fmt_int(cohort_map.get("clusters"))],
        [
            "Window-level singleton clusters",
            fmt_int(cohort_map.get("singleton_clusters")),
        ],
        [
            "Window-level non-singleton clusters",
            fmt_int(cohort_map.get("non_singleton_clusters")),
        ],
        ["Nextclade clades", fmt_int(cohort_map.get("clades"))],
        ["Pango lineages", fmt_int(cohort_map.get("pango_lineages"))],
    ]
    write_latex_table(
        paths,
        TABLE_NAMES["cohort_objects"],
        caption="Chapter 4 analysis cohort and window-cluster object counts",
        label=TABLE_NAMES["cohort_objects"],
        columns=["Quantity", "Value"],
        rows=rows,
        column_spec="lr",
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
                row.policy_period,
                str(row.policy_era).upper().replace("_", " "),
                fmt_int(row.n_windows),
                fmt_int(row.median_window_sequences),
                fmt_int(row.median_window_positive_tests),
                fmt_percent(row.median_window_prop_sequenced),
                f"{fmt_percent(row.min_window_prop_sequenced)}--"
                f"{fmt_percent(row.max_window_prop_sequenced)}",
            ]
        )
    write_latex_table(
        paths,
        TABLE_NAMES["policy_denominators"],
        caption="Rolling-window observation denominators by policy period",
        label=TABLE_NAMES["policy_denominators"],
        columns=[
            "Policy period",
            "Epidemic era",
            "Windows",
            "Median sequences",
            "Median positives",
            "Median sequenced",
            "Coverage range",
        ],
        rows=rows,
        column_spec="lrrrrr",
    )


def write_vaccination_context_table(paths: Paths) -> None:
    vaccination = read_table(paths, "vaccination_context_by_policy")
    vaccination = sort_by_policy(vaccination, column="policy_period")
    rows = []
    for row in vaccination.itertuples(index=False):
        rows.append(
            [
                row.policy_period,
                str(row.policy_era).upper().replace("_", " "),
                fmt_int(row.n_sequences),
                fmt_percent(row.prop_vaccinated),
                fmt_percent(row.prop_booster),
                fmt_float(row.median_dose_number_vaccinated, 1),
                fmt_iqr(
                    row.median_days_since_vaccination,
                    row.q25_days_since_vaccination,
                    row.q75_days_since_vaccination,
                    0,
                ),
            ]
        )
    write_latex_table(
        paths,
        TABLE_NAMES["vaccination_context"],
        caption="Vaccination context of sequenced records by policy period",
        label=TABLE_NAMES["vaccination_context"],
        columns=[
            "Policy period",
            "Epidemic era",
            "Sequences",
            "Vaccinated",
            "Booster recorded",
            "Median dose",
            "Days since vaccination",
        ],
        rows=rows,
        column_spec="lrrrrl",
    )


def write_test_reason_by_policy_era_table(paths: Paths) -> None:
    summary = read_table(paths, "test_reason_by_policy_era").copy()
    if summary.empty:
        pivot = pd.DataFrame(columns=["test_reason"])
    else:
        summary["test_reason"] = summary["test_reason"].astype("string").fillna("missing")
        summary["policy_era"] = summary["policy_era"].astype("string").fillna("missing")
        pivot = summary.pivot_table(
            index="test_reason",
            columns="policy_era",
            values="n_sequences",
            aggfunc="sum",
            fill_value=0,
        )
        era_order = list(POLICY_LABELS)
        ordered_eras = list(era_order)
        extra_eras = [era for era in pivot.columns if era not in era_order]
        extra_eras = sorted(extra_eras, key=lambda value: (str(value) == "missing", str(value)))
        pivot = pivot.reindex(columns=[*ordered_eras, *extra_eras], fill_value=0)
        pivot = pivot.reset_index()
        value_columns = [*ordered_eras, *extra_eras]
        pivot["__total__"] = pivot[value_columns].sum(axis=1)
        pivot = pivot.sort_values(
            ["__total__", "test_reason"],
            ascending=[False, True],
            kind="mergesort",
        ).drop(columns="__total__")

    value_columns = [column for column in pivot.columns if column != "test_reason"]
    rows = []
    for row in pivot.itertuples(index=False):
        rows.append(
            [
                str(row.test_reason).replace("_", " ").title(),
                *[fmt_int(getattr(row, column)) for column in value_columns],
            ]
        )

    columns = [
        "Test reason",
        *[
            POLICY_LABELS.get(column, str(column).replace("_", " ").title())
            for column in value_columns
        ],
    ]
    write_latex_table(
        paths,
        TABLE_NAMES["test_reason_by_policy_era"],
        caption="Unique-sequence test-reason counts by epidemic era",
        label=TABLE_NAMES["test_reason_by_policy_era"],
        columns=columns,
        rows=rows,
        column_spec="l" + "r" * len(value_columns),
    )


def write_cluster_period_table(paths: Paths) -> None:
    clusters = read_table(paths, "cluster_period_summary")
    clusters = sort_by_policy(clusters, column="policy_period")
    rows = []
    for row in clusters.itertuples(index=False):
        rows.append(
            [
                row.policy_period,
                str(row.policy_era).upper().replace("_", " "),
                fmt_int(row.n_clusters),
                fmt_int(row.n_singleton_clusters),
                fmt_int(row.n_non_singleton_clusters),
                fmt_float(row.median_non_singleton_cluster_size, 1),
                fmt_float(row.p90_non_singleton_cluster_size, 1),
                fmt_int(row.max_non_singleton_cluster_size),
                fmt_float(row.median_non_singleton_duration_days, 1),
                fmt_float(row.median_non_singleton_datazones, 1),
            ]
        )
    write_latex_table(
        paths,
        TABLE_NAMES["cluster_period_summary"],
        caption=(
            "EpiLink cluster counts by policy period; size, duration, and "
            "Data Zone summaries are restricted to non-singleton clusters"
        ),
        label=TABLE_NAMES["cluster_period_summary"],
        columns=[
            "Policy period",
            "Epidemic era",
            "All clusters",
            "Singletons",
            "Non-singletons",
            "Median size",
            "P90 size",
            "Max size",
            "Median duration (days)",
            "Median Data Zones",
        ],
        rows=rows,
        column_spec="lrrrrrrrr",
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
                fmt_float(row[mean_col], 4),
                fmt_ci(row[low_col], row[high_col], 4),
                f"{fmt_float(row['window_median'], 3)} "
                f"({fmt_float(row['window_q10'], 3)}--"
                f"{fmt_float(row['window_q90'], 3)})",
            ]
        )
    write_latex_table(
        paths,
        TABLE_NAMES["assortativity_summary"],
        caption=(
            "Pooled assortativity summaries with overlap-based correlated GLS "
            "confidence intervals"
        ),
        label=TABLE_NAMES["assortativity_summary"],
        columns=[
            "Attribute",
            "Windows",
            "Estimated windows",
            "Median lineages",
            "GLS mean r",
            "95% CI",
            "Window median (10--90%)",
        ],
        rows=rows,
        column_spec="lrrrrll",
    )


def write_variance_decomposition_summary_table(paths: Paths) -> None:
    summary = variance_decomposition_summary_table(paths)
    rows = []
    for attribute in ATTRIBUTE_ORDER:
        group = summary.loc[summary["attribute_label"].eq(attribute)]
        if group.empty:
            continue
        row = group.iloc[0]
        rows.append(
            [
                attribute,
                fmt_int(row["n"]),
                fmt_int(row["n_windows"]),
                fmt_int(row["n_lineages"]),
                fmt_float(row["additive_model_fraction"], 3),
                fmt_float(row["adj_additive_model_fraction"], 3),
                fmt_float(row["residual_fraction"], 3),
                fmt_float(row["window_alone_fraction"], 3),
                fmt_float(row["lineage_alone_fraction"], 3),
                fmt_float(row["lineage_given_window_fraction"], 3),
                fmt_float(row["window_given_lineage_fraction"], 3),
                fmt_float(row["weight_cap"], 2),
                fmt_int(row["n_weights_capped"]),
                fmt_float(row["median_boot_se"], 3),
                fmt_float(row["median_ci_width"], 3),
            ]
        )

    write_latex_table(
        paths,
        TABLE_NAMES["variance_decomposition_summary"],
        caption=(
            "Variance decomposition of weighted assortativity at the 90th "
            "percentile inverse-variance weight cap"
        ),
        label=TABLE_NAMES["variance_decomposition_summary"],
        columns=[
            "Attribute",
            "Rows",
            "Windows",
            "Lineages",
            "Additive model",
            "Adj. additive",
            "Residual",
            "Window alone",
            "Lineage alone",
            "Lineage | Window",
            "Window | Lineage",
            "Weight cap",
            "Capped weights",
            "Median boot SE",
            "Median CI width",
        ],
        rows=rows,
        column_spec="l" + "r" * 14,
    )


TABLE_WRITERS: tuple[tuple[str, Callable[[Paths], object]], ...] = (
    (TABLE_NAMES["cohort_objects"], write_cohort_table),
    (TABLE_NAMES["policy_denominators"], write_policy_denominator_table),
    (TABLE_NAMES["vaccination_context"], write_vaccination_context_table),
    (TABLE_NAMES["test_reason_by_policy_era"], write_test_reason_by_policy_era_table),
    (TABLE_NAMES["cluster_period_summary"], write_cluster_period_table),
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
        print(f"Skipped {len(skipped)} table writer(s) because required source tables were missing.")


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
