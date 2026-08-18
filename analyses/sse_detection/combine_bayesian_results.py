"""Consolidate saved Bayesian SSE regression model summaries.

The model-fitting commands write one directory per model specification. This
script joins each saved ``summary.csv`` with its model grid row, fit-frame
summary, metadata, and headline diagnostics. It writes four consolidated tables:

- mixing logistic
- mixing linear
- composition logistic
- composition linear

It also writes three report-facing summary tables for diagnostics, directional
estimates, and directionless variance components.

Run from the repository root::

    python -m analyses.sse_detection.combine_bayesian_results
"""

from __future__ import annotations

import re
from collections.abc import Hashable, Sequence
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from .lib.concurrent_io import (
    atomic_write_csv,
    atomic_write_parquet,
)
from .lib.sse.config import BAYESIAN_OUTPUT_DIR

Family = Literal["logistic", "linear"]
Domain = Literal["mixing", "composition"]

FAMILIES: tuple[Family, ...] = ("logistic", "linear")
DOMAINS: tuple[Domain, ...] = ("mixing", "composition")

MIXING_FEATURE_LABELS: dict[str, str] = {
    "sex_entropy": "Sex",
    "age_entropy": "Age",
    "simd_entropy": "SIMD",
    "urban_rural_entropy": "Urban/rural",
    "health_board_entropy": "Health board",
}
DEFAULT_MIXING_FEATURE_ORDER: tuple[str, ...] = (
    "sex_entropy",
    "age_entropy",
    "simd_entropy",
    "urban_rural_entropy",
    "health_board_entropy",
)
COMPOSITION_VARIABLE_LABELS: dict[str, str] = {
    "sex": "Sex",
    "age_group": "Age group",
    "dz_simd_quintile": "SIMD",
    "urban_rural_class": "Urban/rural",
    "dz_urban_rural_class": "Urban/rural",
    "health_board": "Health board",
    "dz_health_board": "Health board",
}
COMPOSITION_REFERENCE_LABELS: dict[str, str] = {
    "sex": "Male",
    "age_group": "25-64",
    "dz_simd_quintile": "Q1",
    "urban_rural_class": "Large Urban Areas",
    "dz_urban_rural_class": "Large Urban Areas",
    "health_board": "Glasgow and Clyde",
    "dz_health_board": "Glasgow and Clyde",
}
DEFAULT_COMPOSITION_VARIABLE_ORDER: tuple[str, ...] = (
    "sex",
    "age_group",
    "dz_simd_quintile",
    "dz_urban_rural_class",
    "dz_health_board",
)
COMPOSITION_PANEL_GROUPS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    (
        "demographic",
        "Demographic",
        ("sex", "age_group"),
    ),
    (
        "socioeconomic",
        "Sociodemographic",
        (
            "dz_simd_quintile",
            "urban_rural_class",
            "dz_urban_rural_class",
        ),
    ),
    (
        "geographic",
        "Geographic",
        (
            "health_board",
            "dz_health_board",
        ),
    ),
)

MODEL_PREFIX_COLUMNS = [
    "family",
    "domain",
    "outcome",
    "model_set",
    "predictor",
    "parameter",
    "parameter_type",
]
MIXING_PLOT_COLUMNS = [
    "plot_model_kind",
    "plot_scale",
    "plot_term_type",
    "plot_variable",
    "plot_label",
    "plot_panel",
    "plot_variable_order",
    "plot_estimate",
    "plot_hdi95_low",
    "plot_hdi95_high",
    "plot_reference_value",
]
COMPOSITION_PLOT_COLUMNS = [
    "plot_model_kind",
    "plot_term_type",
    "plot_variable",
    "plot_level",
    "plot_reference",
    "plot_level_label",
    "plot_reference_label",
    "plot_contrast_id",
    "plot_label",
    "plot_contrast_label",
    "plot_panel",
    "plot_variable_order",
    "plot_level_order",
    "plot_estimate",
    "plot_hdi95_low",
    "plot_hdi95_high",
    "plot_reference_value",
]
MODEL_SUFFIX_COMMON_COLUMNS = [
    "formula",
    "model_dir",
    "summary_path",
    "diagnostics_path",
    "metadata_path",
    "n_rows",
    "use_sample",
    "full_rows",
    "fit_rows",
    "fit_fraction",
]
LOGISTIC_SUFFIX_COLUMNS = [
    "candidate_rate",
    "full_candidates",
    "fit_candidates",
]
LINEAR_SUFFIX_COLUMNS = [
    "outcome_mean",
    "outcome_sd",
    "full_outcome_mean",
    "fit_outcome_mean",
    "full_outcome_sd",
    "fit_outcome_sd",
]
FIT_SUMMARY_COMMON_COLUMNS = {
    "use_sample",
    "full_rows",
    "fit_rows",
    "fit_fraction",
}
LOGISTIC_FIT_SUMMARY_COLUMNS = {
    "full_candidates",
    "fit_candidates",
}
LINEAR_FIT_SUMMARY_COLUMNS = {
    "full_outcome_mean",
    "fit_outcome_mean",
    "full_outcome_sd",
    "fit_outcome_sd",
}
DIAGNOSTIC_COLUMNS = [
    "diagnostic_status",
    "divergences",
    "draws",
    "divergence_percent",
    "min_bfmi",
    "max_rhat",
    "min_bulk_ess",
    "min_tail_ess",
    "max_tree_depth",
]
POSTERIOR_COLUMNS = [
    "mean",
    "sd",
    "hdi95_lb",
    "hdi95_ub",
    "ess_bulk",
    "ess_tail",
    "r_hat",
    "mcse_mean",
    "mcse_sd",
]
LOGISTIC_COLUMNS = [
    "OR_mean",
    "OR_sd",
    "OR_hdi95_lb",
    "OR_hdi95_ub",
    "OR_ess_bulk",
    "OR_ess_tail",
    "OR_r_hat",
    "OR_mcse_mean",
    "OR_mcse_sd",
    "P(OR > 1 | data)",
    "P(OR < 1 | data)",
]
LINEAR_COLUMNS = ["P(beta > 0 | data)", "P(beta < 0 | data)"]
PLOT_CONTEXT_LABELS: dict[str, str] = {
    "wn_prop_sequenced_z": "Window seq. proportion",
    "dz_cum_incidence_per_capita_z": "Cumulative incidence",
    "dz_cum_prop_sequenced_z": "Cumulative seq. proportion",
}
COMPOSITION_PANEL_LOOKUP = {
    variable: panel_id
    for panel_id, _, variables in COMPOSITION_PANEL_GROUPS
    for variable in variables
}
COMPOSITION_VARIABLE_ORDER = {
    variable: index for index, variable in enumerate(DEFAULT_COMPOSITION_VARIABLE_ORDER)
}
MIXING_FEATURE_ORDER = {
    variable: index for index, variable in enumerate(DEFAULT_MIXING_FEATURE_ORDER)
}

_DIVERGENCE_RE = re.compile(r"(\d+)\s*/\s*(\d+).*?\(([0-9.]+)%\)")
_BFMI_RE = re.compile(r"min=([0-9.]+)")
_CATEGORICAL_RE = re.compile(r"^C\((?P<inside>.+?)\)\[(?P<level>.+)\]$")
_BARE_CATEGORICAL_RE = re.compile(
    r"^(?P<variable>[A-Za-z_][A-Za-z0-9_]*)\[(?P<level>.+)\]$"
)
_TREATMENT_VARIABLE_RE = re.compile(
    r"^(?P<variable>[^,]+),\s*Treatment\(reference=(?P<reference>.+)\)$"
)
_RANDOM_INTERCEPT_RE = re.compile(r"^1\|(?P<group>[^\[]+)\[(?P<level>.+)\]$")

GROUP_LABELS: dict[str, str] = {
    "policy_period": "Policy period",
    "clade": "Clade",
}
GROUP_ORDER: dict[str, int] = {
    "policy_period": 100,
    "clade": 101,
}


def build_consolidated_tables(
    result_dir: Path | str = BAYESIAN_OUTPUT_DIR,
    *,
    families: Sequence[Family] = FAMILIES,
    domains: Sequence[Domain] = DOMAINS,
    strict: bool = False,
) -> tuple[dict[tuple[Domain, Family], pd.DataFrame], pd.DataFrame]:
    """Build requested family/domain consolidated tables.

    Missing or still-running models are skipped by default and reported in the
    returned missing-output table. Set ``strict=True`` to raise instead.
    """
    result_path = Path(result_dir)
    tables: dict[tuple[Domain, Family], pd.DataFrame] = {}
    missing_rows: list[dict[Hashable, Any]] = []
    for domain in domains:
        for family in families:
            table, missing = build_consolidated_table(
                result_path,
                family=family,
                domain=domain,
                strict=strict,
            )
            tables[(domain, family)] = table
            missing_rows.extend(missing)
    return tables, pd.DataFrame(missing_rows)


def build_consolidated_table(
    result_dir: Path | str,
    *,
    family: Family,
    domain: Domain,
    strict: bool = False,
) -> tuple[pd.DataFrame, list[dict[Hashable, Any]]]:
    """Return one consolidated table for a family/domain pair."""
    result_path = Path(result_dir)
    grid = _load_model_grid(result_path, family, domain)
    fit_summary = _load_fit_summary(result_path, family, domain)
    table_parts: list[pd.DataFrame] = []
    missing_rows: list[dict[Hashable, Any]] = []

    for row in grid.to_dict("records"):
        model_dir = _local_model_dir(result_path, row, family=family, domain=domain)
        summary_path = model_dir / "summary.csv"
        diagnostics_path = model_dir / "diagnostics.csv"
        metadata_path = model_dir / "metadata.csv"
        missing_standard = [
            name
            for name, path in (
                ("summary", summary_path),
                ("diagnostics", diagnostics_path),
                ("metadata", metadata_path),
            )
            if not path.exists()
        ]
        if missing_standard:
            missing_rows.append(
                {
                    "family": family,
                    "domain": domain,
                    "outcome": row.get("outcome"),
                    "model_set": row.get("model_set"),
                    "model_dir": str(model_dir),
                    "missing_outputs": ";".join(missing_standard),
                }
            )
            if strict:
                raise FileNotFoundError(
                    "Missing Bayesian model output(s) for "
                    f"{family}/{domain}/{row.get('outcome')}/{row.get('model_set')}: "
                    + ", ".join(missing_standard)
                )
            if "summary" in missing_standard:
                continue

        summary = pd.read_csv(summary_path)
        if "parameter" not in summary.columns:
            summary.insert(0, "parameter", summary.index.astype(str))

        context = _model_context(
            result_path,
            row,
            family=family,
            domain=domain,
            model_dir=model_dir,
            metadata_path=metadata_path,
            diagnostics_path=diagnostics_path,
            fit_summary=fit_summary,
        )
        summary = summary.copy()
        summary.insert(
            1,
            "parameter_type",
            summary["parameter"].map(_parameter_type),
        )
        combined = pd.concat(
            [_repeat_context(context, len(summary)), summary],
            axis=1,
        )
        table_parts.append(_add_plot_columns(combined, family=family, domain=domain))

    if not table_parts:
        return _empty_table(family, domain), missing_rows

    table = pd.concat(table_parts, ignore_index=True)
    return _ordered_table(table, family, domain), missing_rows


def write_consolidated_tables(
    tables: dict[tuple[Domain, Family], pd.DataFrame],
    output_dir: Path | str,
    *,
    write_parquet: bool = True,
) -> dict[tuple[Domain, Family], dict[str, Path]]:
    """Write consolidated tables to CSV and optionally parquet."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    outputs: dict[tuple[Domain, Family], dict[str, Path]] = {}
    for (domain, family), table in tables.items():
        stem = consolidated_table_stem(domain, family)
        csv_path = output_path / f"{stem}.csv"
        atomic_write_csv(table, csv_path, index=False)
        paths = {"csv": csv_path}
        if write_parquet:
            parquet_path = output_path / f"{stem}.parquet"
            atomic_write_parquet(table, parquet_path, index=False)
            paths["parquet"] = parquet_path
        outputs[(domain, family)] = paths
    return outputs


REPORT_TABLE_STEMS: dict[str, str] = {
    "diagnostics": "summary_table_1_diagnostics",
    "estimates": "summary_table_2_estimates",
    "random_effects": "summary_table_3_random_effect_sds",
}

REPORT_COLUMN_LABELS: dict[str, str] = {
    "domain": "Domain",
    "family": "Family",
    "outcome": "Outcome",
    "scale": "Scale",
    "model": "Model",
    "parameter": "Parameter",
    "term_type": "Term Type",
    "component_type": "Component Type",
    "grouping_factor": "Grouping Factor",
    "estimate": "Estimate",
    "hdi95_low": "HDI 95 Low",
    "hdi95_high": "HDI 95 High",
    "effect_scale": "Effect Scale",
    "p_positive_direction": "P Positive Direction",
    "p_negative_direction": "P Negative Direction",
    "favoured_direction": "Favoured Direction",
    "direction_probability": "Direction Probability",
    "direction_band": "Direction Band",
    "random_effect_variance_share": "Random Effect Variance Share",
    "random_effect_sd_rank": "Random Effect SD Rank",
    "diagnostic_status": "Diagnostic Status",
    "divergences": "Divergences",
    "min_bfmi": "Min BFMI",
    "max_rhat": "Max Rhat",
    "min_bulk_ess": "Min Bulk ESS",
    "min_tail_ess": "Min Tail ESS",
    "max_tree_depth": "Max Tree Depth",
}


def build_report_summary_tables(
    tables: dict[tuple[Domain, Family], pd.DataFrame],
) -> dict[str, pd.DataFrame]:
    """Build report-facing summary tables from the consolidated outputs."""
    return {
        "diagnostics": _build_report_diagnostics_summary(tables),
        "estimates": _build_report_estimate_summary(tables),
        "random_effects": _build_report_random_effects_summary(tables),
    }


def write_report_summary_tables(
    tables: dict[tuple[Domain, Family], pd.DataFrame],
    output_dir: Path | str,
    *,
    write_parquet: bool = True,
) -> dict[str, dict[str, Path]]:
    """Write the three report-facing summary tables."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    summaries = build_report_summary_tables(tables)
    outputs: dict[str, dict[str, Path]] = {}
    for name, summary in summaries.items():
        if summary.empty:
            continue

        stem = REPORT_TABLE_STEMS[name]
        csv_path = output_path / f"{stem}.csv"
        atomic_write_csv(summary, csv_path, index=False)
        paths = {"csv": csv_path}
        if write_parquet:
            parquet_path = output_path / f"{stem}.parquet"
            atomic_write_parquet(summary, parquet_path, index=False)
            paths["parquet"] = parquet_path
        outputs[name] = paths
    return outputs


def _build_report_diagnostics_summary(
    tables: dict[tuple[Domain, Family], pd.DataFrame],
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    metric_columns = [
        "diagnostic_status",
        "divergences",
        "min_bfmi",
        "max_rhat",
        "min_bulk_ess",
        "min_tail_ess",
        "max_tree_depth",
    ]
    for (domain, family), table in tables.items():
        if table.empty:
            continue

        available_metrics = [column for column in metric_columns if column in table]
        model_columns = [
            column
            for column in ["outcome", "model_set", "plot_model_kind", "plot_scale"]
            if column in table
        ]
        sub = table.loc[:, model_columns + available_metrics].drop_duplicates()
        sub.insert(0, "family", family)
        sub.insert(0, "domain", domain)
        sub["scale"] = sub.apply(_report_scale_label, axis=1)
        sub["model"] = sub.apply(_report_model_label, axis=1)
        parts.append(sub)

    if not parts:
        return pd.DataFrame()

    summary = pd.concat(parts, ignore_index=True)
    summary["domain"] = summary["domain"].map(_pretty_value)
    summary["family"] = summary["family"].map(_pretty_value)
    summary["outcome"] = summary["outcome"].map(_pretty_value)
    summary["diagnostic_status"] = summary["diagnostic_status"].map(
        _diagnostic_status_label
    )
    columns = [
        "domain",
        "family",
        "outcome",
        "scale",
        "model",
        *metric_columns,
    ]
    columns = [column for column in columns if column in summary]
    return _finalise_report_table(summary.loc[:, columns], columns).sort_values(
        ["Domain", "Family", "Outcome", "Scale", "Model"],
        ignore_index=True,
    )


def _build_report_estimate_summary(
    tables: dict[tuple[Domain, Family], pd.DataFrame],
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    direction_types = {"fixed_effect", "intercept", "random_intercept"}
    for (domain, family), table in tables.items():
        if table.empty or "parameter_type" not in table:
            continue

        sub = table.loc[table["parameter_type"].isin(direction_types)].copy()
        if sub.empty:
            continue

        out = _base_report_rows(sub, domain=domain, family=family)
        out["parameter"] = sub.apply(_report_parameter_label, axis=1)
        out["term_type"] = sub["plot_term_type"].map(_pretty_value)
        out["effect_scale"] = sub.apply(
            lambda row: _effect_scale_label(family, row["parameter_type"]), # type: ignore  # noqa: B023
            axis=1,
        )
        estimates = _estimate_columns(sub, family=family)
        out["estimate"] = estimates["estimate"]
        out["hdi95_low"] = estimates["hdi95_low"]
        out["hdi95_high"] = estimates["hdi95_high"]
        direction = _direction_columns(sub, family=family)
        out = pd.concat([out, direction], axis=1)
        out["diagnostic_status"] = sub["diagnostic_status"].map(
            _diagnostic_status_label
        )
        parts.append(out)

    if not parts:
        return pd.DataFrame()

    columns = [
        "domain",
        "family",
        "outcome",
        "scale",
        "model",
        "parameter",
        "term_type",
        "estimate",
        "hdi95_low",
        "hdi95_high",
        "effect_scale",
        "p_positive_direction",
        "p_negative_direction",
        "favoured_direction",
        "direction_probability",
        "direction_band",
        "diagnostic_status",
    ]
    summary = pd.concat(parts, ignore_index=True)
    return _finalise_report_table(summary.loc[:, columns], columns).sort_values(
        ["Domain", "Family", "Outcome", "Scale", "Model", "Term Type", "Parameter"],
        ignore_index=True,
    )


def _build_report_random_effects_summary(
    tables: dict[tuple[Domain, Family], pd.DataFrame],
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    component_types = {"random_effect_sd", "residual_sd"}
    for (domain, family), table in tables.items():
        if table.empty or "parameter_type" not in table:
            continue

        sub = table.loc[table["parameter_type"].isin(component_types)].copy()
        if sub.empty:
            continue

        out = _base_report_rows(sub, domain=domain, family=family)
        out["component_type"] = sub["parameter_type"].map(_pretty_value)
        out["grouping_factor"] = sub.apply(_grouping_factor_label, axis=1)
        out["parameter"] = sub.apply(_report_parameter_label, axis=1)
        out["effect_scale"] = sub.apply(
            lambda row: _effect_scale_label(family, row["parameter_type"]), # type: ignore  # noqa: B023
            axis=1,
        )
        estimates = _estimate_columns(sub, family=family)
        out["estimate"] = estimates["estimate"]
        out["hdi95_low"] = estimates["hdi95_low"]
        out["hdi95_high"] = estimates["hdi95_high"]
        out["random_effect_variance_share"] = np.nan
        out["random_effect_sd_rank"] = np.nan
        out["diagnostic_status"] = sub["diagnostic_status"].map(
            _diagnostic_status_label
        )

        sd_mask = sub["parameter_type"].eq("random_effect_sd")
        if sd_mask.any():
            model_key = [
                "domain",
                "family",
                "outcome",
                "scale",
                "model",
            ]
            sd_rows = out.loc[sd_mask, model_key].copy()
            sd_rows["_sd_variance"] = (
                pd.to_numeric(sub.loc[sd_mask, "mean"], errors="coerce") ** 2
            )
            variance_total = sd_rows.groupby(model_key, dropna=False)[
                "_sd_variance"
            ].transform("sum")
            out.loc[sd_mask, "random_effect_variance_share"] = (
                sd_rows["_sd_variance"] / variance_total
            ).to_numpy()
            out.loc[sd_mask, "random_effect_sd_rank"] = (
                sd_rows.groupby(model_key, dropna=False)["_sd_variance"]
                .rank(ascending=False, method="dense")
                .astype(float)
                .to_numpy()
            )
        parts.append(out)

    if not parts:
        return pd.DataFrame()

    columns = [
        "domain",
        "family",
        "outcome",
        "scale",
        "model",
        "component_type",
        "grouping_factor",
        "parameter",
        "estimate",
        "hdi95_low",
        "hdi95_high",
        "effect_scale",
        "random_effect_variance_share",
        "random_effect_sd_rank",
        "diagnostic_status",
    ]
    summary = pd.concat(parts, ignore_index=True)
    return _finalise_report_table(summary.loc[:, columns], columns).sort_values(
        [
            "Domain",
            "Family",
            "Outcome",
            "Scale",
            "Model",
            "Component Type",
            "Grouping Factor",
            "Parameter",
        ],
        ignore_index=True,
    )


def _base_report_rows(
    frame: pd.DataFrame, *, domain: Domain, family: Family
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "domain": _pretty_value(domain),
            "family": _pretty_value(family),
            "outcome": frame["outcome"].map(_pretty_value),
            "scale": frame.apply(_report_scale_label, axis=1),
            "model": frame.apply(_report_model_label, axis=1),
        },
        index=frame.index,
    )


def _estimate_columns(frame: pd.DataFrame, *, family: Family) -> pd.DataFrame:
    mean = pd.to_numeric(frame["mean"], errors="coerce")
    low = pd.to_numeric(frame["hdi95_lb"], errors="coerce")
    high = pd.to_numeric(frame["hdi95_ub"], errors="coerce")
    if family == "logistic":
        mean = np.exp(mean)
        low = np.exp(low)
        high = np.exp(high)
    return pd.DataFrame(
        {"estimate": mean, "hdi95_low": low, "hdi95_high": high},
        index=frame.index,
    )


def _direction_columns(frame: pd.DataFrame, *, family: Family) -> pd.DataFrame:
    if family == "logistic":
        p_positive = pd.to_numeric(frame["P(OR > 1 | data)"], errors="coerce")
        p_negative = pd.to_numeric(frame["P(OR < 1 | data)"], errors="coerce")
        positive_label = "Higher Odds"
        negative_label = "Lower Odds"
    else:
        p_positive = pd.to_numeric(frame["P(beta > 0 | data)"], errors="coerce")
        p_negative = pd.to_numeric(frame["P(beta < 0 | data)"], errors="coerce")
        positive_label = "Positive"
        negative_label = "Negative"

    direction_probability = pd.concat([p_positive, p_negative], axis=1).max(axis=1)
    favoured_direction = pd.Series(
        np.where(p_positive >= p_negative, positive_label, negative_label),
        index=frame.index,
    )
    return pd.DataFrame(
        {
            "p_positive_direction": p_positive,
            "p_negative_direction": p_negative,
            "favoured_direction": favoured_direction,
            "direction_probability": direction_probability,
            "direction_band": direction_probability.map(_direction_band_label),
        },
        index=frame.index,
    )


def _direction_band_label(probability: float) -> str:
    if pd.isna(probability):
        return ""
    if probability >= 0.975:
        return "Very Strong"
    if probability >= 0.95:
        return "Strong"
    if probability >= 0.90:
        return "Moderate"
    return "Weak Or Uncertain"


def _effect_scale_label(family: Family, parameter_type: str) -> str:
    if family == "logistic":
        if parameter_type == "intercept":
            return "Odds"
        if parameter_type == "random_effect_sd":
            return "Multiplicative Odds SD"
        if parameter_type == "random_intercept":
            return "Multiplicative Odds"
        return "Odds Ratio"
    if parameter_type in {"random_effect_sd", "residual_sd"}:
        return "Outcome-Scale SD"
    return "Beta"


def _grouping_factor_label(row: pd.Series) -> str:
    parameter_type = str(row.get("parameter_type", ""))
    if parameter_type == "intercept":
        return ""
    if parameter_type == "residual_sd":
        return "Residual"

    variable = row.get("plot_variable", "")
    if pd.isna(variable):
        return ""

    variable = str(variable)
    label = (
        MIXING_FEATURE_LABELS.get(variable)
        or COMPOSITION_VARIABLE_LABELS.get(variable)
        or {
            "policy_period": "Policy Period",
            "clade": "Clade",
        }.get(variable)
    )
    return _pretty_value(label or variable)


def _report_parameter_label(row: pd.Series) -> str:
    if str(row.get("plot_term_type", "")) == "categorical_contrast":
        contrast_label = row.get("plot_contrast_label")
        if _has_value(contrast_label):
            return _pretty_parameter_label(contrast_label)

    plot_label = row.get("plot_label")
    if _has_value(plot_label):
        return _pretty_parameter_label(plot_label)

    parameter = row.get("parameter")
    return "" if not _has_value(parameter) else _pretty_parameter_label(parameter)


def _pretty_parameter_label(value: object) -> str:
    text = "" if not _has_value(value) else str(value).strip()
    replacements = {
        "Age group:": "Age Group:",
        "Health board:": "Health Board:",
        "Policy period:": "Policy Period:",
        "Policy period SD": "Policy Period SD",
        "Urban/rural:": "Urban/Rural:",
        "Urban/rural": "Urban/Rural",
        "Accessible rural": "Accessible Rural",
        "Accessible town": "Accessible Town",
        "Large urban": "Large Urban",
        "Other urban": "Other Urban",
        "Remote rural": "Remote Rural",
        "Remote town": "Remote Town",
        "Window seq. proportion": "Window Sequencing Proportion",
        "Cumulative seq. proportion": "Cumulative Sequencing Proportion",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def _report_model_label(row: pd.Series) -> str:
    model_kind = row.get("plot_model_kind")
    if _has_value(model_kind):
        return _pretty_value(model_kind)

    model_set = str(row.get("model_set", ""))
    if "expanded" in model_set:
        return "Expanded"
    if "primary" in model_set:
        return "Primary"
    return _pretty_value(model_set)


def _report_scale_label(row: pd.Series) -> str:
    scale = row.get("plot_scale", "")
    return "" if not _has_value(scale) else _pretty_value(scale)


def _diagnostic_status_label(value: object) -> str:
    if not _has_value(value):
        return ""
    text = str(value).strip()
    if text.upper() == "OK":
        return "OK"
    return _pretty_value(text)


def _pretty_value(value: object) -> str:
    if not _has_value(value):
        return ""

    text = str(value).strip()
    special = {
        "null_standardised": "Null Standardised",
        "observed": "Observed",
        "ok": "OK",
        "warning": "Warning",
    }
    if text.lower() in special:
        return special[text.lower()]

    text = re.sub(r"[_-]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip().title()
    replacements = {
        " Simd": " SIMD",
        "Simd ": "SIMD ",
        "Simd": "SIMD",
        " Sd": " SD",
        " Or ": " OR ",
        "Hdi": "HDI",
        "Bfmi": "BFMI",
        " Ess": " ESS",
        "Rhat": "Rhat",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def _finalise_report_table(table: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    out = table.copy()
    out = out.rename(
        columns={column: REPORT_COLUMN_LABELS[column] for column in columns}
    )
    for column in out.columns:
        if pd.api.types.is_object_dtype(out[column]) or pd.api.types.is_string_dtype(
            out[column]
        ):
            out[column] = out[column].fillna("")
    return out


def _has_value(value: Any) -> bool:
    if value is None:
        return False
    try:
        if pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    return str(value).strip() != ""


def consolidated_table_stem(domain: Domain, family: Family) -> str:
    """Return the output table stem for one family/domain pair."""
    return f"{domain}_{family}_consolidated_results"


def _load_model_grid(result_dir: Path, family: Family, domain: Domain) -> pd.DataFrame:
    family_dir = result_dir / family
    path = family_dir / f"{domain}_model_grid.csv"
    if path.exists():
        return _filter_grid(pd.read_csv(path), family=family, domain=domain)

    combined_path = family_dir / "model_grid.csv"
    if combined_path.exists():
        return _filter_grid(pd.read_csv(combined_path), family=family, domain=domain)
    return _empty_grid(family, domain)


def _filter_grid(
    frame: pd.DataFrame, *, family: Family, domain: Domain
) -> pd.DataFrame:
    out = frame.copy()
    if "family" in out:
        out = out.loc[out["family"].astype(str).eq(family)]
    if "domain" in out:
        out = out.loc[out["domain"].astype(str).eq(domain)]
    return out.reset_index(drop=True)


def _empty_grid(family: Family, domain: Domain) -> pd.DataFrame:
    return pd.DataFrame(
        columns=["family", "domain", "outcome", "model_set", "predictor", "formula"],
    ).assign(family=family, domain=domain)


def _load_fit_summary(
    result_dir: Path,
    family: Family,
    domain: Domain,
) -> pd.DataFrame:
    family_dir = result_dir / family
    path = family_dir / f"{domain}_fit_frame_summary.csv"
    if path.exists():
        return _filter_fit_summary(pd.read_csv(path), family=family, domain=domain)

    combined_path = family_dir / "fit_frame_summary.csv"
    if combined_path.exists():
        return _filter_fit_summary(
            pd.read_csv(combined_path),
            family=family,
            domain=domain,
        )
    return pd.DataFrame()


def _filter_fit_summary(
    frame: pd.DataFrame,
    *,
    family: Family,
    domain: Domain,
) -> pd.DataFrame:
    out = frame.copy()
    if "family" in out:
        out = out.loc[out["family"].astype(str).eq(family)]
    if "domain" in out:
        out = out.loc[out["domain"].astype(str).eq(domain)]
    return out.reset_index(drop=True)


def _local_model_dir(
    result_dir: Path,
    row: dict[Hashable, Any],
    *,
    family: Family,
    domain: Domain,
) -> Path:
    outcome = str(row.get("outcome", "candidate"))
    model_set = str(row.get("model_set", ""))
    if family == "linear":
        canonical = result_dir / family / domain / outcome / model_set
    else:
        canonical = result_dir / family / domain / model_set
    if canonical.exists():
        return canonical

    raw_path = row.get("model_dir")
    if raw_path is not None and pd.notna(raw_path):
        candidate = Path(str(raw_path))
        if candidate.exists():
            return candidate
    return canonical


def _model_context(
    result_dir: Path,
    row: dict[Hashable, Any],
    *,
    family: Family,
    domain: Domain,
    model_dir: Path,
    metadata_path: Path,
    diagnostics_path: Path,
    fit_summary: pd.DataFrame,
) -> dict[Hashable, Any]:
    metadata = _read_first_row(metadata_path)
    sample = _matching_fit_summary_row(fit_summary, row)
    diagnostics = _diagnostic_context(diagnostics_path)
    formula = metadata.get("formula") or row.get("formula")
    context: dict[Hashable, Any] = {
        "family": family,
        "domain": domain,
        "outcome": row.get("outcome") or metadata.get("outcome"),
        "model_set": row.get("model_set") or metadata.get("model_set"),
        "predictor": row.get("predictor"),
        "formula": formula,
        "model_dir": _display_path(model_dir, result_dir),
        "summary_path": _display_path(model_dir / "summary.csv", result_dir),
        "diagnostics_path": _display_path(diagnostics_path, result_dir),
        "metadata_path": _display_path(metadata_path, result_dir),
        "n_rows": metadata.get("n_rows"),
        "use_sample": metadata.get("use_sample", sample.get("use_sample")),
    }
    if family == "logistic":
        context["candidate_rate"] = metadata.get("candidate_rate")
    else:
        context["outcome_mean"] = metadata.get("outcome_mean")
        context["outcome_sd"] = metadata.get("outcome_sd")
    context.update(sample)
    context.update(diagnostics)
    return context


def _read_first_row(path: Path) -> dict[Hashable, Any]:
    if not path.exists():
        return {}
    frame = pd.read_csv(path)
    if frame.empty:
        return {}
    return frame.iloc[0].to_dict()


def _matching_fit_summary_row(
    fit_summary: pd.DataFrame,
    grid_row: dict[Hashable, Any],
) -> dict[Hashable, Any]:
    if fit_summary.empty:
        return {}
    mask = pd.Series(True, index=fit_summary.index)
    for column in ("family", "domain", "outcome", "model_set"):
        if column in fit_summary.columns and column in grid_row:
            mask &= fit_summary[column].astype(str).eq(str(grid_row[column]))
    matches = fit_summary.loc[mask]
    if matches.empty:
        return {}
    row = matches.iloc[0].to_dict()
    family = str(grid_row.get("family", ""))
    keep = set(FIT_SUMMARY_COMMON_COLUMNS)
    if family == "logistic":
        keep.update(LOGISTIC_FIT_SUMMARY_COLUMNS)
    elif family == "linear":
        keep.update(LINEAR_FIT_SUMMARY_COLUMNS)
    return {key: value for key, value in row.items() if key in keep}


def _add_plot_columns(
    table: pd.DataFrame, *, family: Family, domain: Domain
) -> pd.DataFrame:
    """Add plotting-ready labels and effect columns."""
    out = table.copy()
    out["plot_model_kind"] = out["model_set"].map(_plot_model_kind)
    if domain == "mixing":
        out["plot_scale"] = out["model_set"].map(_plot_scale)
    estimate_col, lower_col, upper_col = _plot_effect_columns(family)
    out["plot_estimate"] = out.get(estimate_col, pd.NA)
    out["plot_hdi95_low"] = out.get(lower_col, pd.NA)
    out["plot_hdi95_high"] = out.get(upper_col, pd.NA)
    if family == "logistic":
        _set_multiplicative_odds_plot_values(out)
    out["plot_reference_value"] = 1.0 if family == "logistic" else 0.0

    for column in _plot_metadata_columns(domain):
        out[column] = pd.NA

    columns = list(_plot_metadata_columns(domain))
    out.loc[:, columns] = out.apply(
        lambda row: _plot_row(row, domain=domain),
        axis=1,
        result_type="expand",
    )
    return out


def _composition_plot_columns() -> tuple[str, ...]:
    return (
        "plot_term_type",
        "plot_variable",
        "plot_level",
        "plot_reference",
        "plot_level_label",
        "plot_reference_label",
        "plot_contrast_id",
        "plot_label",
        "plot_contrast_label",
        "plot_panel",
        "plot_variable_order",
        "plot_level_order",
    )


def _mixing_plot_columns() -> tuple[str, ...]:
    return (
        "plot_term_type",
        "plot_variable",
        "plot_label",
        "plot_panel",
        "plot_variable_order",
    )


def _plot_metadata_columns(domain: Domain) -> tuple[str, ...]:
    return (
        _composition_plot_columns()
        if domain == "composition"
        else _mixing_plot_columns()
    )


def _plot_row(row: pd.Series, *, domain: Domain) -> pd.Series:
    if row["parameter_type"] == "fixed_effect":
        return (
            _composition_plot_row(row)
            if domain == "composition"
            else _mixing_plot_row(row)
        )
    return _generic_plot_row(row, domain=domain)


def _composition_plot_row(row: pd.Series) -> pd.Series:
    parsed = parse_categorical_parameter(str(row["parameter"]))
    variable = parsed.get("variable")
    if variable:
        level = parsed.get("level")
        reference = parsed.get("reference") or COMPOSITION_REFERENCE_LABELS.get(
            variable,
            "reference",
        )
        variable_label = COMPOSITION_VARIABLE_LABELS.get(variable, variable)
        level_label = _composition_display_label(variable, level)
        reference_label = _composition_display_label(variable, reference)
        contrast_id = f"{variable}={level}"
        return pd.Series(
            {
                "plot_term_type": "categorical_contrast",
                "plot_variable": variable,
                "plot_level": level,
                "plot_reference": reference,
                "plot_level_label": level_label,
                "plot_reference_label": reference_label,
                "plot_contrast_id": contrast_id,
                "plot_label": level_label,
                "plot_contrast_label": (
                    f"{variable_label}: {level_label} vs {reference_label}"
                ),
                "plot_panel": COMPOSITION_PANEL_LOOKUP.get(variable, "other"),
                "plot_variable_order": COMPOSITION_VARIABLE_ORDER.get(
                    variable,
                    len(COMPOSITION_VARIABLE_ORDER),
                ),
                "plot_level_order": _plot_level_order(level),
            }
        )

    variable = str(row["parameter"])
    variable_label = PLOT_CONTEXT_LABELS.get(variable, variable)
    return pd.Series(
        {
            "plot_term_type": "continuous_adjuster",
            "plot_variable": variable,
            "plot_level": pd.NA,
            "plot_reference": pd.NA,
            "plot_level_label": pd.NA,
            "plot_reference_label": pd.NA,
            "plot_contrast_id": variable,
            "plot_label": variable_label,
            "plot_contrast_label": pd.NA,
            "plot_panel": "context" if variable in PLOT_CONTEXT_LABELS else "other",
            "plot_variable_order": (
                len(COMPOSITION_VARIABLE_ORDER)
                + list(PLOT_CONTEXT_LABELS).index(variable)
                if variable in PLOT_CONTEXT_LABELS
                else len(COMPOSITION_VARIABLE_ORDER) + len(PLOT_CONTEXT_LABELS)
            ),
            "plot_level_order": "0",
        }
    )


def _mixing_plot_row(row: pd.Series) -> pd.Series:
    parameter = str(row["parameter"])
    variable = (
        parameter
        if parameter in PLOT_CONTEXT_LABELS
        else _mixing_feature_stem(parameter)
    )
    variable_label = MIXING_FEATURE_LABELS.get(
        variable,
        PLOT_CONTEXT_LABELS.get(variable, variable),
    )
    if variable in MIXING_FEATURE_ORDER:
        panel = (
            str(row.get("plot_scale")) if pd.notna(row.get("plot_scale")) else "mixing"
        )
        variable_order = MIXING_FEATURE_ORDER[variable]
    else:
        panel = "context" if variable in PLOT_CONTEXT_LABELS else "other"
        variable_order = (
            len(MIXING_FEATURE_ORDER) + list(PLOT_CONTEXT_LABELS).index(variable)
            if variable in PLOT_CONTEXT_LABELS
            else len(MIXING_FEATURE_ORDER) + len(PLOT_CONTEXT_LABELS)
        )
    return pd.Series(
        {
            "plot_term_type": (
                "mixing_entropy"
                if variable in MIXING_FEATURE_ORDER
                else "continuous_adjuster"
            ),
            "plot_variable": variable,
            "plot_label": variable_label,
            "plot_panel": panel,
            "plot_variable_order": variable_order,
        }
    )


def _generic_plot_row(row: pd.Series, *, domain: Domain) -> pd.Series:
    parameter = str(row["parameter"])
    parameter_type = str(row["parameter_type"])
    group, level = _random_effect_group_level(parameter, parameter_type)

    if parameter_type == "intercept":
        variable = "intercept"
        label = "Intercept"
        panel = "intercept"
        variable_order = -1
    elif parameter_type == "residual_sd":
        variable = "residual_sd"
        label = "Residual SD"
        panel = "scale"
        variable_order = 200
    elif parameter_type == "random_effect_sd":
        variable = group or parameter.removeprefix("1|").removesuffix("_sigma")
        label = f"{_group_label(variable)} SD"
        panel = "random_effect_sd"
        variable_order = GROUP_ORDER.get(variable, 199)
    elif parameter_type == "random_intercept":
        variable = group or "random_effect"
        level_label = _random_effect_level_label(variable, level)
        label = f"{_group_label(variable)}: {level_label}" if level_label else parameter
        panel = "random_intercept"
        variable_order = GROUP_ORDER.get(variable, 150)
    else:
        variable = parameter
        label = parameter
        panel = parameter_type
        variable_order = 999

    values: dict[str, Any] = {
        "plot_term_type": parameter_type,
        "plot_variable": variable,
        "plot_label": label,
        "plot_panel": panel,
        "plot_variable_order": variable_order,
    }
    if domain == "composition":
        values.update(
            {
                "plot_level": level if level else pd.NA,
                "plot_reference": pd.NA,
                "plot_level_label": (
                    _random_effect_level_label(variable, level) if level else pd.NA
                ),
                "plot_reference_label": pd.NA,
                "plot_contrast_id": (
                    f"{variable}={level}" if variable and level else variable
                ),
                "plot_contrast_label": pd.NA,
                "plot_level_order": _plot_level_order(level) if level else pd.NA,
            }
        )
    return pd.Series(values)


def _plot_effect_columns(family: Family) -> tuple[str, str, str]:
    if family == "logistic":
        return "OR_mean", "OR_hdi95_lb", "OR_hdi95_ub"
    return "mean", "hdi95_lb", "hdi95_ub"


def _set_multiplicative_odds_plot_values(out: pd.DataFrame) -> None:
    """Use exponentiated log-scale summaries for logistic random-effect plots."""
    mask = out["parameter_type"].isin(
        {"intercept", "random_intercept", "random_effect_sd"}
    )
    if not mask.any():
        return
    for source, target in (
        ("mean", "plot_estimate"),
        ("hdi95_lb", "plot_hdi95_low"),
        ("hdi95_ub", "plot_hdi95_high"),
    ):
        out.loc[mask, target] = np.exp(
            pd.to_numeric(out.loc[mask, source], errors="coerce")
        )


def _plot_model_kind(model_set: Any) -> Any:
    text = str(model_set)
    if text in {"primary", "expanded"}:
        return text
    if text.endswith("_primary"):
        return "primary"
    if text.endswith("_expanded"):
        return "expanded"
    return pd.NA


def _plot_scale(model_set: Any) -> Any:
    text = str(model_set)
    if text.startswith("observed_"):
        return "observed"
    if text.startswith("null_"):
        return "null_standardised"
    return pd.NA


def _random_effect_group_level(
    parameter: str,
    parameter_type: str,
) -> tuple[str | None, str | None]:
    if parameter_type == "random_intercept":
        match = _RANDOM_INTERCEPT_RE.match(parameter)
        if match:
            return match.group("group"), match.group("level")
    if parameter_type == "random_effect_sd" and parameter.startswith("1|"):
        group = parameter.removeprefix("1|").removesuffix("_sigma")
        return group, None
    return None, None


def _group_label(group: str | None) -> str:
    if not group:
        return "Group"
    return GROUP_LABELS.get(group, group.replace("_", " ").title())


def _random_effect_level_label(group: str | None, level: str | None) -> str | None:
    if level is None:
        return None
    if group == "clade" and level == "recombinant":
        return "Rec."
    return level


def _mixing_feature_stem(parameter: str) -> str:
    for suffix in ("_obs_x10", "_z"):
        if str(parameter).endswith(suffix):
            return str(parameter)[: -len(suffix)]
    return str(parameter)


def _plot_level_order(value: Any) -> str:
    group, key = _level_sort_key(str(value))
    if isinstance(key, float):
        return f"{group}:{key:012.6f}"
    return f"{group}:{key}"


def parse_categorical_parameter(parameter: str) -> dict[str, str | None]:
    """Parse Bambi/formulae categorical parameter names."""
    parameter = str(parameter)
    match = _CATEGORICAL_RE.match(parameter)
    if match:
        inside = match.group("inside")
        level = _clean_level(match.group("level"))
        treatment = _TREATMENT_VARIABLE_RE.match(inside)
        if treatment:
            variable = treatment.group("variable").strip()
            reference = _clean_reference(treatment.group("reference"))
        else:
            variable = inside.strip()
            reference = None
        return {"variable": variable, "level": level, "reference": reference}

    match = _BARE_CATEGORICAL_RE.match(parameter)
    if match:
        return {
            "variable": match.group("variable").strip(),
            "level": _clean_level(match.group("level")),
            "reference": None,
        }

    return {"variable": None, "level": None, "reference": None}


def _clean_reference(value: str) -> str:
    value = str(value).strip()
    if (value.startswith("'") and value.endswith("'")) or (
        value.startswith('"') and value.endswith('"')
    ):
        value = value[1:-1]
    return _clean_level(value)


def _clean_level(value: str) -> str:
    value = str(value).strip()
    value = value.removeprefix("T.")
    if (value.startswith("'") and value.endswith("'")) or (
        value.startswith('"') and value.endswith('"')
    ):
        value = value[1:-1]
    return value


def _composition_display_label(variable: str, value: Any) -> str:
    text = str(value)
    if variable == "dz_simd_quintile":
        return f"SIMD Q{int(eval(text))}"
    if variable in {"urban_rural_class", "dz_urban_rural_class"}:
        return _shorten_urban_rural_label(text)
    if variable in {"health_board", "dz_health_board"}:
        return _shorten_health_board_label(text)
    return text


def _shorten_urban_rural_label(value: str) -> str:
    replacements = {
        "Large Urban Areas": "Large urban",
        "Other Urban Areas": "Other urban",
        "Accessible Small Towns": "Accessible town",
        "Remote Small Towns": "Remote town",
        "Very Remote Small Towns": "Very remote town",
        "Accessible Rural": "Accessible rural",
        "Remote Rural": "Remote rural",
        "Very Remote Rural": "Very remote rural",
    }
    return replacements.get(value, value.replace(" Areas", "").replace(" Small", ""))


def _shorten_health_board_label(value: str) -> str:
    text = value.removeprefix("NHS ").removeprefix("Nhs ")
    return text.replace(" & ", " and ")


def _level_sort_key(value: str) -> tuple[int, float | str]:
    text = str(value)
    try:
        return (0, float(text))
    except ValueError:
        return (1, text)


def _diagnostic_context(path: Path) -> dict[str, Any]:
    out: dict[str, Any] = {column: pd.NA for column in DIAGNOSTIC_COLUMNS}
    if not path.exists():
        return out

    diagnostics = pd.read_csv(path)
    if diagnostics.empty:
        return out
    statuses = set(diagnostics.get("Status", pd.Series(dtype=str)).astype(str))
    out["diagnostic_status"] = "WARNING" if "WARNING" in statuses else "OK"
    by_name = diagnostics.set_index("Diagnostic")["Value"].astype(str).to_dict()

    divergence = _DIVERGENCE_RE.search(by_name.get("Divergences", ""))
    if divergence:
        out["divergences"] = int(divergence.group(1))
        out["draws"] = int(divergence.group(2))
        out["divergence_percent"] = float(divergence.group(3))

    bfmi = _BFMI_RE.search(by_name.get("BFMI", ""))
    if bfmi:
        out["min_bfmi"] = float(bfmi.group(1))
    for diagnostic, column in (
        ("Max R-hat", "max_rhat"),
        ("Min bulk ESS", "min_bulk_ess"),
        ("Min tail ESS", "min_tail_ess"),
        ("Max tree depth", "max_tree_depth"),
    ):
        if diagnostic in by_name:
            out[column] = _float_or_na(by_name[diagnostic])
    return out


def _float_or_na(value: Any) -> Any:
    try:
        return float(value)
    except (TypeError, ValueError):
        return pd.NA


def _repeat_context(context: dict[Hashable, Any], rows: int) -> pd.DataFrame:
    return pd.DataFrame([context] * rows)


def _parameter_type(parameter: Any) -> str:
    name = str(parameter)
    if name == "Intercept":
        return "intercept"
    if name == "sigma":
        return "residual_sd"
    if name.startswith("1|") and name.endswith("_sigma"):
        return "random_effect_sd"
    if name.startswith("1|"):
        return "random_intercept"
    return "fixed_effect"


def _display_path(path: Path, result_dir: Path) -> str:
    try:
        return str(path.relative_to(result_dir))
    except ValueError:
        return str(path)


def _empty_table(family: Family, domain: Domain) -> pd.DataFrame:
    return pd.DataFrame(columns=_ordered_columns(pd.Index([]), family, domain))


def _ordered_table(table: pd.DataFrame, family: Family, domain: Domain) -> pd.DataFrame:
    columns = _ordered_columns(table.columns, family, domain)
    for column in columns:
        if column not in table.columns:
            table[column] = pd.NA
    return table.loc[:, columns]


def _ordered_columns(
    existing: pd.Index[Any],
    family: Family,
    domain: Domain,
) -> pd.Index[Any]:
    family_columns = LOGISTIC_COLUMNS if family == "logistic" else LINEAR_COLUMNS
    suffix_columns = [
        *MODEL_SUFFIX_COMMON_COLUMNS,
        *(LOGISTIC_SUFFIX_COLUMNS if family == "logistic" else LINEAR_SUFFIX_COLUMNS),
    ]
    plot_columns = (
        COMPOSITION_PLOT_COLUMNS if domain == "composition" else MIXING_PLOT_COLUMNS
    )
    preferred = [
        *MODEL_PREFIX_COLUMNS,
        *plot_columns,
        *suffix_columns,
        *DIAGNOSTIC_COLUMNS,
        *POSTERIOR_COLUMNS,
        *family_columns,
    ]
    seen = set()
    ordered = []
    for column in preferred:
        if column not in seen:
            ordered.append(column)
            seen.add(column)
    for column in existing:
        if column not in seen:
            ordered.append(column)
            seen.add(column)
    return pd.Index(ordered)


def main() -> int:
    """Build all consolidated and report-facing Bayesian result tables."""
    result_dir = BAYESIAN_OUTPUT_DIR.resolve()
    if not result_dir.exists():
        raise FileNotFoundError(f"Missing Bayesian output directory: {result_dir}")

    output_dir = result_dir / "consolidated_tables"
    tables, missing = build_consolidated_tables(
        result_dir,
        families=FAMILIES,
        domains=DOMAINS,
        strict=False,
    )
    outputs = write_consolidated_tables(
        tables,
        output_dir,
        write_parquet=True,
    )
    summary_outputs = write_report_summary_tables(
        tables,
        output_dir,
        write_parquet=True,
    )
    for (domain, family), paths in outputs.items():
        table = tables[(domain, family)]
        print(
            f"Wrote {domain}/{family}: {len(table):,} rows, "
            f"{len(table.columns):,} columns -> {paths['csv']}"
        )
    for name, paths in summary_outputs.items():
        print(f"Wrote report {name} summary -> {paths['csv']}")
    if not missing.empty:
        print(
            f"Skipped {len(missing):,} model(s) with incomplete outputs. "
            "Set strict=True in build_consolidated_tables(...) to fail on these."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
