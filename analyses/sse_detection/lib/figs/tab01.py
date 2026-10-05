"""Build bayesian model specifications table."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from utils import write_latex_table

from ..sse.config import BAYESIAN_OUTPUT_DIR
from .common import Paths, add_common_args, latex_table_path, paths_from_args
from .table_common import (
    _consolidated_table_dir,
    _display_text,
    _format_float,
    _format_int,
    _format_percent,
    _sort_for_report,
    _write_data_table,
)

TABLE_NAME = "tab01_bayesian_model_specifications"

CONSOLIDATED_RESULT_STEMS = (
    "mixing_logistic_consolidated_results",
    "mixing_linear_consolidated_results",
    "composition_logistic_consolidated_results",
    "composition_linear_consolidated_results",
)


def _read_consolidated_results(result_dir: Path) -> pd.DataFrame:
    frames = []
    base = _consolidated_table_dir(result_dir)
    for stem in CONSOLIDATED_RESULT_STEMS:
        path = base / f"{stem}.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing table: {path}")
        frames.append(pd.read_csv(path))
    return pd.concat(frames, ignore_index=True, sort=False)


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
    return _sort_for_report(out)


def build(paths: Paths) -> dict[str, Path]:
    table = build_model_sample_specification_table(paths.bayesian_result_dir)
    name = TABLE_NAME
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
            "Fitted samples for the association models. Candidate outcomes give "
            "CSC counts and percentages; score outcomes give observed means and "
            "standard deviations. Expanded models add the three standardised "
            "surveillance measures to the primary predictors. All models include "
            "policy-period and clade varying intercepts."
        ),
        short_caption="Bayesian characterisation fitted samples and model specifications.",
        label="tab:bayesian_model_specifications",
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
            r"P{0.127\linewidth}P{0.080\linewidth}P{0.103\linewidth}"
            r"P{0.195\linewidth}P{0.110\linewidth}P{0.105\linewidth}"
            r"P{0.155\linewidth}"
        ),
    )
    return {
        "csv": paths.result_table_dir / f"{name}.csv",
        "parquet": paths.result_table_dir / f"{name}.parquet",
        "tex": tex_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
