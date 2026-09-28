"""Shared data loading and formatting for publication tables."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from utils.display_labels import normalise_attribute_label
from utils.latex_tables import latex_escape

SUMMARY_TABLE_STEMS = {
    "diagnostics": "summary_table_1_diagnostics",
    "estimates": "summary_table_2_estimates",
    "random_effects": "summary_table_3_random_effect_sds",
}

DOMAIN_ORDER = {"Mixing": 0, "Composition": 1}

FAMILY_ORDER = {"Logistic": 0, "Linear": 1}

OUTCOME_ORDER = {
    "Candidate": 0,
    "Candidate Status": 0,
    "Burst Score": 1,
    "Burden Score": 2,
}

SCALE_ORDER = {"": -1, "Observed": 0, "Null Standardised": 1}

MODEL_ORDER = {"Primary": 0, "Expanded": 1}

TERM_TYPE_ORDER = {
    "Mixing Entropy": 0,
    "Categorical Contrast": 0,
    "Continuous Adjuster": 1,
    "Intercept": 2,
    "Random Intercept": 3,
}


def _write_data_table(table: pd.DataFrame, output_dir: Path, name: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / f"{name}.csv", index=False)
    table.to_parquet(output_dir / f"{name}.parquet", index=False)


def _consolidated_table_dir(result_dir: Path) -> Path:
    return result_dir / "consolidated_tables"


def _read_summary_table(result_dir: Path, table: str) -> pd.DataFrame:
    path = _consolidated_table_dir(result_dir) / f"{SUMMARY_TABLE_STEMS[table]}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing table: {path}")
    frame = pd.read_csv(path)
    for column in ("Parameter", "Grouping Factor"):
        if column in frame:
            frame[column] = frame[column].map(
                lambda value: (
                    normalise_attribute_label(value)
                    if isinstance(value, str)
                    else value
                )
            )
    return frame


def _display_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value)


def _format_int(value: Any) -> str:
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return f"{int(float(value)):,}"


def _format_float(value: Any, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return f"{float(value):.{digits}f}"


def _format_probability(value: Any) -> str:
    return _format_float(value, 3)


def _format_percent(value: Any, digits: int = 1) -> str:
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return f"{100 * float(value):.{digits}f}%"


def _estimate_digits(effect_scale: Any) -> int:
    text = _display_text(effect_scale).lower()
    return 2 if "odds" in text else 3


def _format_effect_interval(row: pd.Series) -> str:
    digits = _estimate_digits(row.get("Effect Scale", ""))
    estimate = _format_float(row.get("Estimate"), digits)
    low = _format_float(row.get("HDI 95 Low"), digits)
    high = _format_float(row.get("HDI 95 High"), digits)
    if not estimate or not low or not high:
        return estimate
    return f"{estimate} [{low}, {high}]"


def _format_direction(row: pd.Series) -> str:
    direction = _display_text(row.get("Favoured Direction"))
    probability = _format_probability(row.get("Direction Probability"))
    if not direction or not probability:
        return ""
    return f"P({direction}) = {probability}"


def _sort_for_report(
    table: pd.DataFrame,
    *,
    extra_columns: list[str] | None = None,
) -> pd.DataFrame:
    extra_columns = extra_columns or []
    out = table.copy()
    domain = out.get("Domain", pd.Series("", index=out.index)).fillna("")
    out["_domain_order"] = domain.map(DOMAIN_ORDER).fillna(99)
    family = out.get("Family", pd.Series("", index=out.index)).fillna("")
    out["_family_order"] = family.map(FAMILY_ORDER).fillna(99)
    outcome = out.get("Outcome", pd.Series("", index=out.index)).fillna("")
    out["_outcome_order"] = outcome.map(OUTCOME_ORDER).fillna(99)
    scale = out.get("Scale", pd.Series("", index=out.index)).fillna("")
    out["_scale_order"] = scale.map(SCALE_ORDER).fillna(99)
    model = out.get("Model", pd.Series("", index=out.index)).fillna("")
    out["_model_order"] = model.map(MODEL_ORDER).fillna(99)
    if "Term Type" in out:
        out["_term_type_order"] = out["Term Type"].map(TERM_TYPE_ORDER).fillna(99)
    sort_columns = [
        "_domain_order",
        "_family_order",
        "_outcome_order",
        "_scale_order",
        "_model_order",
    ]
    if "_term_type_order" in out:
        sort_columns.append("_term_type_order")
    sort_columns.extend(extra_columns)
    return (
        out.sort_values(sort_columns, kind="stable")
        .drop(columns=[column for column in out if str(column).startswith("_")])
        .reset_index(drop=True)
    )


def write_size_table(
    path: Path,
    caption: str,
    label: str,
    headers: list[str],
    rows: list[list],
    *,
    long=False,
):
    path.parent.mkdir(parents=True, exist_ok=True)
    spec = (
        "@{}"
        + "l" * min(2, len(headers))
        + "r" * (len(headers) - min(2, len(headers)))
        + "@{}"
    )
    header = " & ".join(latex_escape(x) for x in headers) + r" \\"
    body = [" & ".join(latex_escape(x) for x in row) + r" \\" for row in rows]
    if long:
        lines = [
            r"\begingroup\small",
            r"\setlength{\tabcolsep}{4pt}",
            r"\begin{longtable}{" + spec + "}",
            r"\caption{" + caption + r"}\label{" + label + r"}\\",
            r"\toprule",
            header,
            r"\midrule\endfirsthead",
            r"\multicolumn{"
            + str(len(headers))
            + r"}{l}{\tablename\ \thetable\ continued}\\",
            r"\toprule",
            header,
            r"\midrule\endhead",
            r"\bottomrule\endfoot",
            *body,
            r"\end{longtable}\endgroup",
        ]
    else:
        lines = [
            r"\begin{table}[htbp]\centering",
            r"\caption{" + caption + r"}\label{" + label + "}",
            r"\begin{thesistablebody}{" + spec + "}",
            r"\toprule",
            header,
            r"\midrule",
            *body,
            r"\bottomrule\end{thesistablebody}",
            r"\end{table}",
        ]
    path.write_text("\n".join(lines) + "\n")
