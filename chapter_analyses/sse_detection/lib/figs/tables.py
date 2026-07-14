"""Build publication and supplementary tables for Bayesian mixing models."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
from typing import Callable

import numpy as np
import pandas as pd

from .forest import (
    DEFAULT_MIXING_FEATURE_ORDER,
    MIXING_FEATURE_LABELS,
    _collect_mixing_forest_rows,
)
from ..sse.config import BAYESIAN_OUTPUT_DIR
from .common import Paths, add_common_args, latex_table_path, paths_from_args


OUTCOMES = (
    ("logistic", "candidate", "Candidate status"),
    ("linear", "burst_score", "Burst score"),
    ("linear", "burden_score", "Burden score"),
)
TABLE_NAMES = {
    "main": "tab_ch5_bayesian_mixing_main",
    "full": "tab_ch5_bayesian_mixing_full",
    "samples": "tab_ch5_bayesian_mixing_samples",
    "diagnostics": "tab_ch5_bayesian_mixing_diagnostics",
    "specifications": "tab_ch5_bayesian_mixing_specifications",
}


def _latex_escape(value: object) -> str:
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


def _render_latex_table(
    *,
    caption: str,
    label: str,
    columns: list[str],
    rows: list[list[object]],
    column_spec: str,
    addlinespace_after: set[int] | None = None,
    landscape: bool = False,
) -> str:
    addlinespace_after = addlinespace_after or set()
    body = []
    for row_idx, row in enumerate(rows):
        if len(row) != len(columns):
            raise ValueError(f"Row {row_idx} has the wrong number of cells")
        body.append("    " + " & ".join(_latex_escape(cell) for cell in row) + r" \\")
        if row_idx in addlinespace_after and row_idx < len(rows) - 1:
            body.append(r"    \addlinespace[0.35em]")
    header = " & ".join(f"\\textbf{{{_latex_escape(col)}}}" for col in columns)
    lines = []
    if landscape:
        lines.append(r"\begin{landscape}")
    lines.extend(
        [
            r"\begin{table}[htbp]",
            r"\centering",
            f"\\caption[{_latex_escape(caption)}]{{{_latex_escape(caption)}}}"
            f"\\label{{{label}}}",
            r"\begingroup",
            r"\small",
            r"\setlength{\tabcolsep}{4pt}",
            r"\renewcommand{\arraystretch}{1.12}",
            r"\begin{adjustbox}{max width=\textwidth,center}",
            f"\\begin{{tabular}}{{@{{}}{column_spec}@{{}}}}",
            r"\toprule",
            f"{header} " + r"\\",
            r"\midrule",
            *body,
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{adjustbox}",
            r"\endgroup",
            r"\end{table}",
        ]
    )
    if landscape:
        lines.append(r"\end{landscape}")
    return "\n".join(lines) + "\n"


def _write_data_table(table: pd.DataFrame, output_dir: Path, name: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / f"{name}.csv", index=False)
    table.to_parquet(output_dir / f"{name}.parquet", index=False)


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
    records: list[dict[str, object]] = []
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
                estimate = float(row.OR_mean)
                low = float(row.OR_hdi95_lb)
                high = float(row.OR_hdi95_ub)
                summary = pd.read_csv(row.summary_path).set_index("parameter")
                source = summary.loc[row.parameter]
                p_positive = float(source["P(OR > 1 | data)"])
                p_negative = float(source["P(OR < 1 | data)"])
                effect_type = "odds ratio"
                reference = 1.0
                positive_label, negative_label = "OR > 1", "OR < 1"
            else:
                estimate = float(row.mean)
                low = float(row.hdi95_lb)
                high = float(row.hdi95_ub)
                summary = pd.read_csv(row.summary_path).set_index("parameter")
                source = summary.loc[row.parameter]
                p_positive = float(source["P(beta > 0 | data)"])
                p_negative = float(source["P(beta < 0 | data)"])
                effect_type = "coefficient"
                reference = 0.0
                positive_label, negative_label = "beta > 0", "beta < 0"
            if estimate >= reference:
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
                raise ValueError(f"Missing sample summary for {outcome}/{row.model_set}")
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
    return full.loc[full["scale"].eq("null") & full["model"].eq("expanded")].reset_index(
        drop=True
    )


def build_sample_table(result_dir: Path = BAYESIAN_OUTPUT_DIR) -> pd.DataFrame:
    samples = _sample_summaries(result_dir)
    eligible_n = int(samples.loc[samples["model_set"].str.startswith("observed"), "full_rows"].max())
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
            divergence_match = re.search(r"(\d+)\s*/\s*(\d+).*\(([0-9.]+)%\)", divergence_text)
            bfmi_match = re.search(r"min=([0-9.]+)", bfmi_text)
            statuses = set(table["Status"].astype(str))
            overall = "WARNING" if "WARNING" in statuses else "OK"
            records.append(
                {
                    "outcome": outcome_label,
                    "scale": model_set.split("_")[0],
                    "model": model_set.split("_")[1],
                    "divergences": int(divergence_match.group(1)) if divergence_match else np.nan,
                    "draws": int(divergence_match.group(2)) if divergence_match else np.nan,
                    "divergence_percent": float(divergence_match.group(3)) if divergence_match else np.nan,
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


def write_main_table(paths: Paths) -> dict[str, Path]:
    table = build_main_results_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["main"]
    _write_data_table(table, paths.result_table_dir, name)
    rows = []
    boundaries = set()
    for outcome_idx, (_, group) in enumerate(table.groupby("outcome", sort=False)):
        start = len(rows)
        for row_idx, (_, row) in enumerate(group.iterrows()):
            rows.append(
                [
                    row["outcome_label"] if row_idx == 0 else "",
                    row["mixing_variable_label"],
                    _effect_text(row),
                    f"P({row['direction']}) = {row['direction_probability']:.3f}",
                    f"{int(row['n']):,}",
                ]
            )
        if outcome_idx < table["outcome"].nunique() - 1:
            boundaries.add(len(rows) - 1)
    tex_path = latex_table_path(paths, name)
    tex_path.write_text(
        _render_latex_table(
            caption=(
                "Expanded null-standardised Bayesian mixing-model estimates. "
                "Values are posterior means with 95 percent highest-density intervals."
            ),
            label="tab:ch5_bayesian_mixing_main",
            columns=["Outcome", "Mixing measure", "Effect (95% HDI)", "Direction", "n"],
            rows=rows,
            column_spec="lllrr",
            addlinespace_after=boundaries,
        )
    )
    return {"csv": paths.result_table_dir / f"{name}.csv", "parquet": paths.result_table_dir / f"{name}.parquet", "tex": tex_path}


def write_full_table(paths: Paths) -> dict[str, Path]:
    table = build_full_results_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["full"]
    _write_data_table(table, paths.result_table_dir, name)
    # Separate LaTeX fragments keep each supplementary table page-sized.
    outputs: dict[str, Path] = {
        "csv": paths.result_table_dir / f"{name}.csv",
        "parquet": paths.result_table_dir / f"{name}.parquet",
    }
    for outcome, group in table.groupby("outcome", sort=False):
        rows = [
            [
                row["scale"].title(),
                row["model"].title(),
                row["mixing_variable_label"],
                _effect_text(row),
                f"P({row['direction']}) = {row['direction_probability']:.3f}",
                f"{int(row['n']):,}",
            ]
            for _, row in group.iterrows()
        ]
        tex_path = latex_table_path(paths, f"{name}_{outcome}")
        tex_path.write_text(
            _render_latex_table(
                caption=f"Full primary and expanded Bayesian mixing estimates for {group.iloc[0]['outcome_label'].lower()}.",
                label=f"tab:ch5_bayesian_mixing_full_{outcome}",
                columns=["Scale", "Model", "Mixing measure", "Effect (95% HDI)", "Direction", "n"],
                rows=rows,
                column_spec="llllrr",
            )
        )
        outputs[f"tex_{outcome}"] = tex_path
    return outputs


def write_sample_table(paths: Paths) -> dict[str, Path]:
    table = build_sample_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["samples"]
    _write_data_table(table, paths.result_table_dir, name)
    rows = []
    for _, row in table.iterrows():
        candidate = "" if pd.isna(row["candidate_n"]) else f"{int(row['candidate_n']):,} ({100 * row['candidate_rate']:.1f}%)"
        outcome_summary = "" if pd.isna(row["outcome_mean"]) else f"{row['outcome_mean']:.3f} ({row['outcome_sd']:.3f})"
        rows.append(
            [
                row["outcome"],
                row["scale"].title(),
                f"{int(row['model_n']):,}",
                f"{int(row['outcome_missing']):,}",
                f"{int(row['mixing_complete_case_exclusions']):,}",
                candidate,
                outcome_summary,
            ]
        )
    tex_path = latex_table_path(paths, name)
    tex_path.write_text(
        _render_latex_table(
            caption="Bayesian mixing-model sample sizes and outcome availability.",
            label="tab:ch5_bayesian_mixing_samples",
            columns=["Outcome", "Scale", "Model n", "Outcome missing", "Mixing excluded", "Candidates", "Mean (SD)"],
            rows=rows,
            column_spec="llrrrrr",
        )
    )
    return {"csv": paths.result_table_dir / f"{name}.csv", "parquet": paths.result_table_dir / f"{name}.parquet", "tex": tex_path}


def write_diagnostics_table(paths: Paths) -> dict[str, Path]:
    table = build_diagnostics_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["diagnostics"]
    _write_data_table(table, paths.result_table_dir, name)
    rows = [
        [
            row["outcome"], row["scale"].title(), row["model"].title(),
            f"{int(row['divergences'])}", f"{row['min_bfmi']:.3f}",
            f"{row['max_rhat']:.3f}", f"{int(row['min_bulk_ess'])}",
            f"{int(row['min_tail_ess'])}", row["status"],
        ]
        for _, row in table.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    tex_path.write_text(
        _render_latex_table(
            caption="Bayesian mixing-model sampling diagnostics.",
            label="tab:ch5_bayesian_mixing_diagnostics",
            columns=["Outcome", "Scale", "Model", "Divergences", "Min BFMI", "Max R-hat", "Min bulk ESS", "Min tail ESS", "Status"],
            rows=rows,
            column_spec="lllrrrrrl",
        )
    )
    return {"csv": paths.result_table_dir / f"{name}.csv", "parquet": paths.result_table_dir / f"{name}.parquet", "tex": tex_path}


def write_specification_table(paths: Paths) -> dict[str, Path]:
    table = build_specification_table(paths.bayesian_result_dir)
    name = TABLE_NAMES["specifications"]
    _write_data_table(table, paths.result_table_dir, name)
    compact = table.drop_duplicates(["scale", "model"])
    rows = [
        [
            row["scale"].title(), row["model"].title(), row["predictor_scaling"],
            (
                "None"
                if row["context_adjusters"] == "none"
                else "Sequencing coverage, incidence, cumulative coverage"
            ),
        ]
        for _, row in compact.iterrows()
    ]
    tex_path = latex_table_path(paths, name)
    tex_path.write_text(
        _render_latex_table(
            caption=(
                "Bayesian mixing-model specifications and predictor scaling. "
                "All models include varying intercepts for policy period and clade; "
                "complete formulas are provided in the accompanying data table."
            ),
            label="tab:ch5_bayesian_mixing_specifications",
            columns=["Scale", "Model", "Predictor scaling", "Additional fixed effects"],
            rows=rows,
            column_spec="llll",
        )
    )
    return {"csv": paths.result_table_dir / f"{name}.csv", "parquet": paths.result_table_dir / f"{name}.parquet", "tex": tex_path}


TABLE_WRITERS: tuple[tuple[str, Callable[[Paths], dict[str, Path]]], ...] = (
    (TABLE_NAMES["main"], write_main_table),
    (TABLE_NAMES["full"], write_full_table),
    (TABLE_NAMES["samples"], write_sample_table),
    (TABLE_NAMES["diagnostics"], write_diagnostics_table),
    (TABLE_NAMES["specifications"], write_specification_table),
)


def write_tables(
    paths: Paths,
) -> dict[str, dict[str, Path]]:
    paths.result_table_dir.mkdir(parents=True, exist_ok=True)
    return {
        name: writer(paths) for name, writer in TABLE_WRITERS
    }


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
