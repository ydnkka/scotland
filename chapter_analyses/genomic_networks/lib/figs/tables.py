"""Build Chapter 4 observation-network LaTeX table fragments."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys
from typing import Callable

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (
    ATTRIBUTE_ORDER,
    Paths,
    add_common_args,
    paths_from_args,
    read_table,
    sort_by_policy,
    window_idx_from_id,
)
from chapter_analyses.genomic_networks.lib.config import TABLES_DIR  # noqa: E402


SIMD_GROUP_LABELS = {
    5: "quintile",
    10: "decile",
    20: "vigintile",
}
TABLE_NAMES = {
    "cohort_objects": "tab_ch4_cohort_objects",
    "policy_denominators": "tab_ch4_policy_denominators",
    "vaccination_context": "tab_ch4_vaccination_context",
    "cluster_period_summary": "tab_ch4_cluster_period_summary",
    "assortativity_summary": "tab_ch4_assortativity_summary",
    "simd_population_weighting": "tab_ch4_simd_population_weighting",
}
TABLE_LABELS = {
    "cohort_objects": "tab:ch4_cohort_objects",
    "policy_denominators": "tab:ch4_policy_denominators",
    "vaccination_context": "tab:ch4_vaccination_context",
    "cluster_period_summary": "tab:ch4_cluster_period_summary",
    "assortativity_summary": "tab:ch4_assortativity_summary",
    "simd_population_weighting": "tab:ch4_simd_population_weighting",
}


def fmt_int(value: float | int | str | None) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    try:
        value_float = float(value)
    except (TypeError, ValueError):
        return str(value)
    if value_float.is_integer():
        return f"{int(value_float):,}"
    return f"{value_float:,.1f}"


def fmt_percent(value: float | int | str | None) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{100 * float(value):.1f}%"


def fmt_percent_points(value: float | int | str | None) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{float(value):.1f}%"


def fmt_float(value: float | int | str | None, digits: int = 2) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{float(value):.{digits}f}"


def fmt_ci(low: float | int | None, high: float | int | None, digits: int = 2) -> str:
    if (
        low is None
        or high is None
        or (isinstance(low, float) and np.isnan(low))
        or (isinstance(high, float) and np.isnan(high))
    ):
        return "not estimated"
    return f"{fmt_float(low, digits)}--{fmt_float(high, digits)}"


def fmt_iqr(
    median: float | int | None,
    q25: float | int | None,
    q75: float | int | None,
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
    rows: list[list[object]],
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
    label: str = TABLE_LABELS["simd_population_weighting"],
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
                    str(int(row.simd_group)),
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
        window_group = window_summary.loc[window_summary["attribute"].eq(attribute)]
        rows.append(
            {
                "graph": "Compatibility",
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
        ["Window-level clusters", fmt_int(cohort_map.get("clusters"))],
        ["Nextclade clades", fmt_int(cohort_map.get("clades"))],
        ["Pango lineages", fmt_int(cohort_map.get("pango_lineages"))],
    ]
    write_latex_table(
        paths,
        TABLE_NAMES["cohort_objects"],
        caption="Chapter 4 analysis cohort and window-cluster object counts",
        label=TABLE_LABELS["cohort_objects"],
        columns=["Quantity", "Value"],
        rows=rows,
        column_spec="lr",
    )


def write_policy_denominator_table(paths: Paths) -> None:
    denominators = read_table(paths, "window_denominator_contrasts")
    denominators = sort_by_policy(denominators)
    rows = []
    for row in denominators.itertuples(index=False):
        rows.append(
            [
                row.policy_period,
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
        label=TABLE_LABELS["policy_denominators"],
        columns=[
            "Period",
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
    vaccination = sort_by_policy(vaccination)
    rows = []
    for row in vaccination.itertuples(index=False):
        rows.append(
            [
                row.policy_period,
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
        label=TABLE_LABELS["vaccination_context"],
        columns=[
            "Period",
            "Sequences",
            "Vaccinated",
            "Booster recorded",
            "Median dose",
            "Days since vaccination",
        ],
        rows=rows,
        column_spec="lrrrrl",
    )


def write_cluster_period_table(paths: Paths) -> None:
    clusters = read_table(paths, "cluster_period_summary")
    clusters = sort_by_policy(clusters)
    rows = []
    for row in clusters.itertuples(index=False):
        rows.append(
            [
                row.policy_period,
                fmt_int(row.n_clusters),
                fmt_float(row.median_cluster_size, 1),
                fmt_float(row.p90_cluster_size, 1),
                fmt_int(row.max_cluster_size),
                fmt_float(row.median_duration_days, 1),
                fmt_float(row.median_datazones, 1),
            ]
        )
    write_latex_table(
        paths,
        TABLE_NAMES["cluster_period_summary"],
        caption="Window-level EpiLink cluster summaries by policy period",
        label=TABLE_LABELS["cluster_period_summary"],
        columns=[
            "Period",
            "Clusters",
            "Median size",
            "P90 size",
            "Max size",
            "Median duration",
            "Median Data Zones",
        ],
        rows=rows,
        column_spec="lrrrrrr",
    )


def write_assortativity_summary_table(paths: Paths) -> None:
    summary = compatibility_attribute_summary(paths)
    rows = []
    addlinespace_after: set[int] = set()
    graph_groups = list(summary.groupby("graph", sort=False))
    for graph_idx, (graph_label, df) in enumerate(graph_groups):
        graph_row_start = len(rows)
        for attribute in ATTRIBUTE_ORDER:
            group = df.loc[df["attribute_label"].eq(attribute)]
            if group.empty:
                continue
            row = group.iloc[0]
            rows.append(
                [
                    graph_label if len(rows) == graph_row_start else "",
                    attribute,
                    fmt_int(row["n_windows"]),
                    fmt_int(row["n_networks"])
                    if not pd.isna(row["n_networks"])
                    else "",
                    fmt_float(row["weighted_mean"], 4),
                    fmt_ci(row["ci_low"], row["ci_high"], 4),
                    f"{fmt_float(row['window_median'], 3)} "
                    f"({fmt_float(row['window_q10'], 3)}--"
                    f"{fmt_float(row['window_q90'], 3)})",
                ]
            )
        if graph_idx < len(graph_groups) - 1 and len(rows) > graph_row_start:
            addlinespace_after.add(len(rows) - 1)
    write_latex_table(
        paths,
        TABLE_NAMES["assortativity_summary"],
        caption="Weighted assortativity summaries with compatibility confidence intervals",
        label=TABLE_LABELS["assortativity_summary"],
        columns=[
            "Graph",
            "Attribute",
            "Windows",
            "Networks",
            "Weighted mean r",
            "95% CI",
            "Window median (10--90%)",
        ],
        rows=rows,
        column_spec="llrrrll",
        addlinespace_after=addlinespace_after,
    )


TABLE_WRITERS: tuple[tuple[str, Callable[[Paths], object]], ...] = (
    (TABLE_NAMES["cohort_objects"], write_cohort_table),
    (TABLE_NAMES["policy_denominators"], write_policy_denominator_table),
    (TABLE_NAMES["vaccination_context"], write_vaccination_context_table),
    (TABLE_NAMES["cluster_period_summary"], write_cluster_period_table),
    (TABLE_NAMES["assortativity_summary"], write_assortativity_summary_table),
    (TABLE_NAMES["simd_population_weighting"], write_simd_population_weighting_table),
)


def write_tables(paths: Paths) -> None:
    paths.figure_dir.mkdir(parents=True, exist_ok=True)
    for _, write_table in TABLE_WRITERS:
        write_table(paths)


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
