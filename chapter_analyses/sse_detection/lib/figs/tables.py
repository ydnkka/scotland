"""Build publication and supplementary tables for Bayesian characterisation models."""

from __future__ import annotations

import argparse
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from utils import write_latex_longtable, write_latex_table

from ..sse.config import BAYESIAN_OUTPUT_DIR
from .common import Paths, add_common_args, latex_table_path, paths_from_args
from .forest import (
    DEFAULT_MIXING_FEATURE_ORDER,
    MIXING_FEATURE_LABELS,
    _collect_mixing_forest_rows,
)

OUTCOMES = (
    ("logistic", "candidate", "Candidate status"),
    ("linear", "burst_score", "Burst score"),
    ("linear", "burden_score", "Burden score"),
)
TABLE_NAMES = {
    "model_specifications": "tab_ch5_bayesian_model_specifications",
    "model_diagnostics": "tab_ch5_bayesian_model_diagnostics",
    "fixed_effects_main": "tab_ch5_bayesian_fixed_effects_main",
    "fixed_effects_full": "tab_app_ch5_bayesian_fixed_effects_full",
    "random_effects": "tab_app_ch5_bayesian_random_effect_sds",
}
SUMMARY_TABLE_STEMS = {
    "diagnostics": "summary_table_1_diagnostics",
    "estimates": "summary_table_2_estimates",
    "random_effects": "summary_table_3_random_effect_sds",
}
CONSOLIDATED_RESULT_STEMS = (
    "mixing_logistic_consolidated_results",
    "mixing_linear_consolidated_results",
    "composition_logistic_consolidated_results",
    "composition_linear_consolidated_results",
)
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
    return pd.read_csv(path)


def _read_consolidated_results(result_dir: Path) -> pd.DataFrame:
    frames = []
    base = _consolidated_table_dir(result_dir)
    for stem in CONSOLIDATED_RESULT_STEMS:
        path = base / f"{stem}.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing table: {path}")
        frames.append(pd.read_csv(path))
    return pd.concat(frames, ignore_index=True, sort=False)


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


def _sort_for_thesis(
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


def _sample_summaries(result_dir: Path) -> pd.DataFrame:
    frames = []
    for family in ("logistic", "linear"):
        path = result_dir / family / "mixing_fit_frame_summary.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing table: {path}")
        frame = pd.read_csv(path)
        frame["family"] = family
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def build_full_results_table(result_dir: Path = BAYESIAN_OUTPUT_DIR) -> pd.DataFrame:
    """Return all focal primary/expanded mixing estimates on both scales."""
    samples = _sample_summaries(result_dir)
    records: list[dict[str, Any]] = []
    for family, outcome, outcome_label in OUTCOMES:
        rows, missing = _collect_mixing_forest_rows(
            result_dir / family,
            family=family,
            outcome=outcome,
            feature_order=DEFAULT_MIXING_FEATURE_ORDER,
        )
        if missing:
            raise FileNotFoundError(
                f"Missing {outcome} model summaries: {', '.join(missing)}"
            )
        for row in rows.itertuples(index=False):
            if family == "logistic":
                estimate = cast(float, row.OR_mean)
                low = cast(float, row.OR_hdi95_lb)
                high = cast(float, row.OR_hdi95_ub)
                summary = pd.read_csv(str(row.summary_path)).set_index("parameter")
                source = summary.loc[row.parameter]
                p_positive = cast(float, source["P(OR > 1 | data)"])
                p_negative = cast(float, source["P(OR < 1 | data)"])
                effect_type = "odds ratio"
                reference = 1.0
                positive_label, negative_label = "OR > 1", "OR < 1"
            else:
                estimate = cast(float, row.mean)
                low = cast(float, row.hdi95_lb)
                high = cast(float, row.hdi95_ub)
                summary = pd.read_csv(str(row.summary_path)).set_index("parameter")
                source = summary.loc[row.parameter]
                p_positive = cast(float, source["P(beta > 0 | data)"])
                p_negative = cast(float, source["P(beta < 0 | data)"])
                effect_type = "coefficient"
                reference = 0.0
                positive_label, negative_label = "beta > 0", "beta < 0"
            if cast(float, estimate) >= reference:
                direction = positive_label
                direction_probability = p_positive
            else:
                direction = negative_label
                direction_probability = p_negative
            sample = samples.loc[
                samples["family"].eq(family)
                & samples["outcome"].eq(outcome)
                & samples["model_set"].eq(str(row.model_set))
            ]
            if sample.empty:
                raise ValueError(
                    f"Missing sample summary for {outcome}/{row.model_set}"
                )
            records.append(
                {
                    "outcome": outcome,
                    "outcome_label": outcome_label,
                    "mixing_variable": str(row.feature),
                    "mixing_variable_label": MIXING_FEATURE_LABELS[str(row.feature)],
                    "scale": str(row.scale),
                    "scaling": (
                        "per 0.1-bit increase"
                        if str(row.scale) == "observed"
                        else "per 1-SD increase"
                    ),
                    "model": str(row.model_kind),
                    "effect_type": effect_type,
                    "estimate": estimate,
                    "hdi95_low": low,
                    "hdi95_high": high,
                    "direction": direction,
                    "direction_probability": direction_probability,
                    "n": int(sample.iloc[0]["fit_rows"]),
                }
            )
    out = pd.DataFrame(records)
    outcome_order = {value[1]: idx for idx, value in enumerate(OUTCOMES)}
    feature_order = {
        value: idx for idx, value in enumerate(DEFAULT_MIXING_FEATURE_ORDER)
    }
    scale_order = {"observed": 0, "null": 1}
    model_order = {"primary": 0, "expanded": 1}
    return (
        out.assign(
            _outcome=out["outcome"].map(outcome_order),
            _feature=out["mixing_variable"].map(feature_order),
            _scale=out["scale"].map(scale_order),
            _model=out["model"].map(model_order),
        )
        .sort_values(["_outcome", "_scale", "_feature", "_model"])
        .drop(columns=["_outcome", "_feature", "_scale", "_model"])
        .reset_index(drop=True)
    )


def build_main_results_table(result_dir: Path = BAYESIAN_OUTPUT_DIR) -> pd.DataFrame:
    """Return expanded null-standardised estimates used for main interpretation."""
    full = build_full_results_table(result_dir)
    return full.loc[
        full["scale"].eq("null") & full["model"].eq("expanded")
    ].reset_index(drop=True)


def build_sample_table(result_dir: Path = BAYESIAN_OUTPUT_DIR) -> pd.DataFrame:
    samples = _sample_summaries(result_dir)
    eligible_n = int(
        samples.loc[samples["model_set"].str.startswith("observed"), "full_rows"].max()
    )
    rows = []
    for family, outcome, outcome_label in OUTCOMES:
        for scale in ("observed", "null"):
            source = samples.loc[
                samples["family"].eq(family)
                & samples["outcome"].eq(outcome)
                & samples["model_set"].eq(f"{scale}_expanded")
            ]
            if source.empty:
                raise ValueError(f"Missing sample summary for {outcome}/{scale}")
            row = source.iloc[0]
            outcome_available = int(
                samples.loc[samples["outcome"].eq(outcome), "full_rows"].max()
            )
            rows.append(
                {
                    "outcome": outcome_label,
                    "scale": scale,
                    "eligible_clusters": eligible_n,
                    "outcome_available": outcome_available,
                    "outcome_missing": eligible_n - outcome_available,
                    "model_n": int(row["fit_rows"]),
                    "mixing_complete_case_exclusions": outcome_available
                    - int(row["fit_rows"]),
                    "candidate_n": (
                        int(row["fit_candidates"])
                        if "fit_candidates" in row.index
                        and pd.notna(row["fit_candidates"])
                        else np.nan
                    ),
                    "candidate_rate": (
                        float(row["fit_candidate_rate"])
                        if "fit_candidate_rate" in row.index
                        and pd.notna(row["fit_candidate_rate"])
                        else np.nan
                    ),
                    "outcome_mean": (
                        float(row["fit_outcome_mean"])
                        if "fit_outcome_mean" in row.index
                        and pd.notna(row["fit_outcome_mean"])
                        else np.nan
                    ),
                    "outcome_sd": (
                        float(row["fit_outcome_sd"])
                        if "fit_outcome_sd" in row.index
                        and pd.notna(row["fit_outcome_sd"])
                        else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def _diagnostic_value(table: pd.DataFrame, diagnostic: str) -> str:
    rows = table.loc[table["Diagnostic"].eq(diagnostic), "Value"]
    return "" if rows.empty else str(rows.iloc[0])


def build_diagnostics_table(result_dir: Path = BAYESIAN_OUTPUT_DIR) -> pd.DataFrame:
    records = []
    for family, outcome, outcome_label in OUTCOMES:
        base = result_dir / family / "mixing"
        if family == "linear":
            base = base / outcome
        for model_set in (
            "observed_primary",
            "observed_expanded",
            "null_primary",
            "null_expanded",
        ):
            path = base / model_set / "diagnostics.csv"
            if not path.exists():
                raise FileNotFoundError(f"Missing table: {path}")
            table = pd.read_csv(path)
            divergence_text = _diagnostic_value(table, "Divergences")
            bfmi_text = _diagnostic_value(table, "BFMI")
            divergence_match = re.search(
                r"(\d+)\s*/\s*(\d+).*\(([0-9.]+)%\)", divergence_text
            )
            bfmi_match = re.search(r"min=([0-9.]+)", bfmi_text)
            statuses = set(table["Status"].astype(str))
            overall = "WARNING" if "WARNING" in statuses else "OK"
            records.append(
                {
                    "outcome": outcome_label,
                    "scale": model_set.split("_")[0],
                    "model": model_set.split("_")[1],
                    "divergences": int(divergence_match.group(1))
                    if divergence_match
                    else np.nan,
                    "draws": int(divergence_match.group(2))
                    if divergence_match
                    else np.nan,
                    "divergence_percent": float(divergence_match.group(3))
                    if divergence_match
                    else np.nan,
                    "min_bfmi": float(bfmi_match.group(1)) if bfmi_match else np.nan,
                    "max_rhat": float(_diagnostic_value(table, "Max R-hat")),
                    "min_bulk_ess": float(_diagnostic_value(table, "Min bulk ESS")),
                    "min_tail_ess": float(_diagnostic_value(table, "Min tail ESS")),
                    "max_tree_depth": float(_diagnostic_value(table, "Max tree depth")),
                    "status": overall,
                }
            )
    return pd.DataFrame(records)


def build_specification_table(result_dir: Path = BAYESIAN_OUTPUT_DIR) -> pd.DataFrame:
    frames = []
    for family in ("logistic", "linear"):
        path = result_dir / family / "mixing_model_grid.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing table: {path}")
        frame = pd.read_csv(path)
        frame["scale"] = frame["model_set"].str.split("_").str[0]
        frame["model"] = frame["model_set"].str.split("_").str[1]
        frame["predictor_scaling"] = np.where(
            frame["scale"].eq("observed"),
            "0.1-bit entropy increase",
            "1-SD null-standardised entropy increase",
        )
        frame["context_adjusters"] = np.where(
            frame["model"].eq("expanded"),
            "window sequencing proportion; cumulative incidence; cumulative sequencing proportion",
            "none",
        )
        frame["varying_intercepts"] = "policy period; clade"
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)[
        [
            "family",
            "outcome",
            "scale",
            "model",
            "predictor_scaling",
            "context_adjusters",
            "varying_intercepts",
            "formula",
        ]
    ]


def _effect_text(row: pd.Series) -> str:
    digits = 3 if row["effect_type"] == "coefficient" else 2
    return (
        f"{row['estimate']:.{digits}f} "
        f"[{row['hdi95_low']:.{digits}f}, {row['hdi95_high']:.{digits}f}]"
    )


def _model_display_fields(table: pd.DataFrame) -> pd.DataFrame:
    work = table.copy()
    work["Domain"] = work["domain"].str.replace("_", " ").str.title()
    work["Family"] = work["family"].str.replace("_", " ").str.title()
    work["Outcome"] = work["outcome"].str.replace("_", " ").str.title()
    work["Outcome"] = work["Outcome"].replace({"Candidate": "Candidate"})
    work["Scale"] = np.where(
        work["domain"].eq("mixing"),
        np.where(
            work["model_set"].astype(str).str.startswith("observed"),
            "Observed",
            "Null Standardised",
        ),
        "",
    )
    work["Model"] = np.where(
        work["model_set"].astype(str).str.contains("expanded"),
        "Expanded",
        "Primary",
    )
    return work


def _predictor_scale_label(row: pd.Series) -> str:
    domain = _display_text(row.get("domain"))
    scale = _display_text(row.get("Scale"))
    if domain == "composition":
        return "Categorical sequence composition"
    if scale == "Observed":
        return "Observed entropy per 0.1"
    if scale == "Null Standardised":
        return "Null-standardised entropy per 1 SD"
    return ""


def _primary_predictor_label(row: pd.Series) -> str:
    if _display_text(row.get("domain")) == "composition":
        return (
            "Sex; age group; SIMD quintile; urban/rural class; health board "
            "categorical contrasts"
        )
    return "Sex, age-group, SIMD, urban/rural, and health-board entropy"


def _additional_fixed_effects_label(row: pd.Series) -> str:
    if _display_text(row.get("Model")) != "Expanded":
        return "None"
    return (
        "Window sequencing proportion; cumulative incidence; cumulative "
        "sequencing proportion"
    )


def _fitted_outcome_summary(row: pd.Series) -> str:
    family = _display_text(row.get("family"))
    if family == "logistic":
        candidate_count = row.get("fit_candidates", row.get("full_candidates", np.nan))
        rate = row.get("candidate_rate", row.get("fit_candidate_rate", np.nan))
        count = _format_int(candidate_count)
        percentage = _format_percent(rate)
        if count and percentage:
            return f"{count} ({percentage})"
        return count

    mean = row.get("fit_outcome_mean", row.get("outcome_mean", np.nan))
    sd = row.get("fit_outcome_sd", row.get("outcome_sd", np.nan))
    mean_text = _format_float(mean, 3)
    sd_text = _format_float(sd, 3)
    if mean_text and sd_text:
        return f"{mean_text} ({sd_text})"
    return mean_text


def build_model_sample_specification_table(
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
) -> pd.DataFrame:
    """Return fitted sample accounting and compact model specification rows."""
    consolidated = _read_consolidated_results(result_dir)
    context_columns = [
        "domain",
        "family",
        "outcome",
        "model_set",
        "predictor",
        "formula",
        "full_rows",
        "fit_rows",
        "fit_fraction",
        "candidate_rate",
        "fit_candidate_rate",
        "full_candidates",
        "fit_candidates",
        "fit_outcome_mean",
        "fit_outcome_sd",
        "outcome_mean",
        "outcome_sd",
    ]
    available_columns = [
        column for column in context_columns if column in consolidated.columns
    ]
    models = (
        consolidated.loc[:, available_columns].drop_duplicates().reset_index(drop=True)
    )
    models = _model_display_fields(models)
    out = pd.DataFrame(
        {
            "Domain": models["Domain"],
            "Family": models["Family"],
            "Outcome": models["Outcome"],
            "Scale": models["Scale"],
            "Model": models["Model"],
            "Predictor Scale": models.apply(_predictor_scale_label, axis=1),
            "Available Rows": models["full_rows"],
            "Fit Rows": models["fit_rows"],
            "Fit Fraction": models["fit_fraction"],
            "Fitted Outcome": models.apply(_fitted_outcome_summary, axis=1),
            "Primary Predictors": models.apply(_primary_predictor_label, axis=1),
            "Additional Fixed Effects": models.apply(
                _additional_fixed_effects_label, axis=1
            ),
            "Varying Intercepts": "Policy period; clade",
            "Formula": models["formula"],
        }
    )
    return _sort_for_thesis(out)


def build_model_diagnostics_table(
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
) -> pd.DataFrame:
    """Return the all-domain diagnostic summary table."""
    return _sort_for_thesis(_read_summary_table(result_dir, "diagnostics"))


def build_fixed_effects_main_table(
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
) -> pd.DataFrame:
    """Return focal fixed effects from expanded interpretation models."""
    table = _read_summary_table(result_dir, "estimates")
    scale = table["Scale"].fillna("")
    focal = table["Term Type"].isin(["Categorical Contrast", "Mixing Entropy"])
    interpretation_model = table["Model"].eq("Expanded") & (
        table["Domain"].eq("Composition")
        | (table["Domain"].eq("Mixing") & scale.eq("Null Standardised"))
    )
    out = table.loc[focal & interpretation_model].copy()
    return _sort_for_thesis(out, extra_columns=["Parameter"])


def build_fixed_effects_full_table(
    result_dir: Path = BAYESIAN_OUTPUT_DIR,
) -> pd.DataFrame:
    """Return the complete thesis estimate summary table."""
    return _sort_for_thesis(
        _read_summary_table(result_dir, "estimates"),
        extra_columns=["Parameter"],
    )


def build_random_effects_table(result_dir: Path = BAYESIAN_OUTPUT_DIR) -> pd.DataFrame:
    """Return random-effect and residual-SD summaries."""
    return _sort_for_thesis(
        _read_summary_table(result_dir, "random_effects"),
        extra_columns=["Grouping Factor", "Parameter"],
    )


def write_model_sample_specification_table(paths: Paths) -> dict[str, Path]:
    table = build_model_sample_specification_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["model_specifications"]
    _write_data_table(table, paths.result_table_dir, name)
    rows = [
        [
            row["Domain"],
            row["Family"],
            row["Outcome"],
            row["Predictor Scale"],
            row["Model"],
            _format_int(row["Fit Rows"]),
            row["Fitted Outcome"],
        ]
        for _, row in table.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    write_latex_table(
        tex_path,
        caption=(
            "Bayesian characterisation fitted samples and model specifications. "
            "Fitted outcome gives candidate count and rate for logistic models, "
            "and posterior-outcome mean and standard deviation for Gaussian models. "
            "Primary predictors vary by domain: composition models use sex, age "
            "group, SIMD quintile, urban/rural class, and health-board categorical "
            "contrasts; mixing models use sex, age-group, SIMD, urban/rural, and "
            "health-board entropy predictors. Additional fixed effects vary by "
            "model: primary models include none, while expanded models include "
            "window sequencing proportion, cumulative incidence, and cumulative "
            "sequencing proportion. "
            "All specifications include varying intercepts for policy period and clade."
        ),
        short_caption="Bayesian characterisation fitted samples and model specifications.",
        label="tab:ch5_bayesian_model_specifications",
        columns=[
            "Domain",
            "Family",
            "Outcome",
            "Predictor scale",
            "Model",
            "Fitted rows",
            "Fitted outcome",
        ],
        rows=rows,
        column_spec=(
            r"P{0.105\linewidth}P{0.080\linewidth}P{0.110\linewidth}"
            r"P{0.230\linewidth}P{0.095\linewidth}P{0.105\linewidth}"
            r"P{0.150\linewidth}"
        ),
    )
    return {
        "csv": paths.result_table_dir / f"{name}.csv",
        "parquet": paths.result_table_dir / f"{name}.parquet",
        "tex": tex_path,
    }


def write_model_diagnostics_table(paths: Paths) -> dict[str, Path]:
    table = build_model_diagnostics_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["model_diagnostics"]
    _write_data_table(table, paths.result_table_dir, name)
    rows = [
        [
            row["Domain"],
            row["Family"],
            row["Outcome"],
            row["Scale"],
            row["Model"],
            row["Diagnostic Status"],
            _format_int(row["Divergences"]),
            _format_float(row["Min BFMI"], 3),
            _format_float(row["Max Rhat"], 3),
            _format_int(row["Min Bulk ESS"]),
            _format_int(row["Min Tail ESS"]),
            _format_int(row["Max Tree Depth"]),
        ]
        for _, row in table.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    write_latex_longtable(
        tex_path,
        caption=(
            "Bayesian characterisation sampling diagnostics across composition "
            "and mixing models."
        ),
        label="tab:ch5_bayesian_model_diagnostics",
        columns=[
            "Domain",
            "Family",
            "Outcome",
            "Scale",
            "Model",
            "Status",
            "Div.",
            "Min BFMI",
            "Max Rhat",
            "Min bulk ESS",
            "Min tail ESS",
            "Max tree depth",
        ],
        rows=rows,
        column_spec=(
            r"P{0.070\linewidth}P{0.055\linewidth}P{0.075\linewidth}"
            r"P{0.090\linewidth}P{0.060\linewidth}P{0.060\linewidth}"
            r"P{0.045\linewidth}P{0.055\linewidth}P{0.055\linewidth}"
            r"P{0.065\linewidth}P{0.065\linewidth}P{0.060\linewidth}"
        ),
    )
    return {
        "csv": paths.result_table_dir / f"{name}.csv",
        "parquet": paths.result_table_dir / f"{name}.parquet",
        "tex": tex_path,
    }


def write_fixed_effects_main_table(paths: Paths) -> dict[str, Path]:
    table = build_fixed_effects_main_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["fixed_effects_main"]
    _write_data_table(table, paths.result_table_dir, name)
    rows = [
        [
            row["Domain"],
            row["Family"],
            row["Outcome"],
            row["Scale"],
            row["Model"],
            row["Parameter"],
            _format_effect_interval(row),
            _format_direction(row),
            row["Diagnostic Status"],
        ]
        for _, row in table.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    write_latex_longtable(
        tex_path,
        caption=(
            "Focal fixed-effect estimates from expanded Bayesian "
            "characterisation models. Values are posterior means with 95 percent "
            "highest-density intervals."
        ),
        short_caption="Focal fixed-effect estimates from expanded Bayesian models.",
        label="tab:ch5_bayesian_fixed_effects_main",
        columns=[
            "Domain",
            "Family",
            "Outcome",
            "Scale",
            "Model",
            "Term",
            "Effect (95% HDI)",
            "Direction",
            "Diagnostic",
        ],
        rows=rows,
        column_spec=(
            r"P{0.065\linewidth}P{0.055\linewidth}P{0.075\linewidth}"
            r"P{0.090\linewidth}P{0.060\linewidth}P{0.255\linewidth}"
            r"P{0.120\linewidth}P{0.120\linewidth}P{0.070\linewidth}"
        ),
    )
    return {
        "csv": paths.result_table_dir / f"{name}.csv",
        "parquet": paths.result_table_dir / f"{name}.parquet",
        "tex": tex_path,
    }


def write_fixed_effects_full_table(paths: Paths) -> dict[str, Path]:
    table = build_fixed_effects_full_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["fixed_effects_full"]
    _write_data_table(table, paths.result_table_dir, name)
    rows = []
    for _, row in table.iterrows():
        domain_family = f"{row['Domain']}; {row['Family']}"
        rows.append(
            [
                domain_family,
                row["Outcome"],
                row["Scale"],
                row["Model"],
                row["Parameter"],
                row["Term Type"],
                _format_effect_interval(row),
                row["Effect Scale"],
                _format_probability(row["P Positive Direction"]),
                _format_probability(row["P Negative Direction"]),
                _format_direction(row),
                row["Direction Band"],
                row["Diagnostic Status"],
            ]
        )
    tex_path = latex_table_path(paths, name)
    write_latex_longtable(
        tex_path,
        caption=(
            "Complete Bayesian characterisation estimate summary. This table "
            "contains all rows from thesis summary table 2, including intercepts, "
            "focal fixed effects, continuous adjusters, and random-intercept "
            "level summaries."
        ),
        short_caption="Complete Bayesian characterisation estimate summary.",
        label="tab:app_ch5_bayesian_fixed_effects_full",
        columns=[
            "Domain; family",
            "Outcome",
            "Scale",
            "Model",
            "Parameter",
            "Term type",
            "Effect (95% HDI)",
            "Effect scale",
            "P(+)",
            "P(-)",
            "Direction",
            "Band",
            "Diagnostic",
        ],
        rows=rows,
        column_spec=(
            r"P{0.070\linewidth}P{0.060\linewidth}P{0.065\linewidth}"
            r"P{0.050\linewidth}P{0.185\linewidth}P{0.070\linewidth}"
            r"P{0.095\linewidth}P{0.065\linewidth}P{0.035\linewidth}"
            r"P{0.035\linewidth}P{0.075\linewidth}P{0.050\linewidth}"
            r"P{0.040\linewidth}"
        ),
        tiny=True,
    )
    return {
        "csv": paths.result_table_dir / f"{name}.csv",
        "parquet": paths.result_table_dir / f"{name}.parquet",
        "tex": tex_path,
    }


def write_random_effects_table(paths: Paths) -> dict[str, Path]:
    table = build_random_effects_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["random_effects"]
    _write_data_table(table, paths.result_table_dir, name)
    rows = [
        [
            row["Domain"],
            row["Family"],
            row["Outcome"],
            row["Scale"],
            row["Model"],
            row["Component Type"],
            row["Grouping Factor"],
            _format_effect_interval(row),
            row["Effect Scale"],
            _format_percent(row["Random Effect Variance Share"]),
            _format_int(row["Random Effect SD Rank"]),
            row["Diagnostic Status"],
        ]
        for _, row in table.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    write_latex_longtable(
        tex_path,
        caption=(
            "Random-effect and residual standard-deviation summaries from "
            "Bayesian characterisation models. Variance share and rank are shown "
            "for random-effect standard deviations only."
        ),
        short_caption="Bayesian random-effect and residual standard deviations.",
        label="tab:app_ch5_bayesian_random_effect_sds",
        columns=[
            "Domain",
            "Family",
            "Outcome",
            "Scale",
            "Model",
            "Component",
            "Group",
            "SD (95% HDI)",
            "Effect scale",
            "Variance share",
            "Rank",
            "Diagnostic",
        ],
        rows=rows,
        column_spec=(
            r"P{0.065\linewidth}P{0.055\linewidth}P{0.075\linewidth}"
            r"P{0.090\linewidth}P{0.060\linewidth}P{0.080\linewidth}"
            r"P{0.075\linewidth}P{0.105\linewidth}P{0.115\linewidth}"
            r"P{0.065\linewidth}P{0.045\linewidth}P{0.060\linewidth}"
        ),
    )
    return {
        "csv": paths.result_table_dir / f"{name}.csv",
        "parquet": paths.result_table_dir / f"{name}.parquet",
        "tex": tex_path,
    }


TABLE_WRITERS: tuple[tuple[str, Callable[[Paths], dict[str, Path]]], ...] = (
    (TABLE_NAMES["model_specifications"], write_model_sample_specification_table),
    (TABLE_NAMES["model_diagnostics"], write_model_diagnostics_table),
    (TABLE_NAMES["fixed_effects_main"], write_fixed_effects_main_table),
    (TABLE_NAMES["fixed_effects_full"], write_fixed_effects_full_table),
    (TABLE_NAMES["random_effects"], write_random_effects_table),
)


def write_tables(
    paths: Paths,
) -> dict[str, dict[str, Path]]:
    paths.result_table_dir.mkdir(parents=True, exist_ok=True)
    return {name: writer(paths) for name, writer in TABLE_WRITERS}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    outputs = write_tables(paths)
    for name in outputs:
        print(f"Wrote {name} to {paths.result_table_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
