"""Prepare SSE regression data, formulas, and output paths.

It prepares:

- node-level mixing models, where the row is a cluster/node
- sequence-level composition models, where the row is a sequence-window record
- logistic candidate outcomes and linear score outcomes
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Sequence, Any

import numpy as np
import pandas as pd

from ..sse.config import PROJECT_ROOT, RESULTS_DIR
from ..concurrent_io import atomic_write_csv, atomic_write_parquet, exclusive_file_lock
from ..sse.entropy import (
    DEFAULT_MIXING_FEATURES,
    OBSERVED_MIXING_FEATURES_X10,
    add_observed_mixing_entropy_scales,
)
from ..sse.io import HIGH_PRIORITY_CANDIDATE_TIERS, SseOutputs


RANDOM_SEED = 123
CLUSTER_ID_COL = "cluster_id"
GROUP_VARS = ("policy_period", "clade")
SCORE_OUTCOMES = ("burst_score", "burden_score")

EXCLUDE_NULL_FEATURES = frozenset(
    {
        "datazone_entropy_z",
        "local_authority_entropy_z",
    }
)
EXCLUDE_OBSERVED_FEATURES = frozenset(
    {
        "datazone_entropy_obs_x10",
        "local_authority_entropy_obs_x10",
    }
)

NULL_MIXING_FEATURES = [
    feature
    for feature in DEFAULT_MIXING_FEATURES
    if feature not in EXCLUDE_NULL_FEATURES
]
OBSERVED_MIXING_FEATURES = [
    feature
    for feature in OBSERVED_MIXING_FEATURES_X10
    if feature not in EXCLUDE_OBSERVED_FEATURES
]

COMPOSITION_SPECS = (
    {
        "name": "sex",
        "column": "sex",
        "reference": "Male",
        "label": "Sex",
    },
    {
        "name": "age_band",
        "column": "age_band",
        "reference": "20-24",
        "label": "Age band",
    },
    {
        "name": "simd_quintile",
        "column": "dz_simd_quintile",
        "reference": 1,
        "label": "SIMD quintile",
    },
    {
        "name": "urban_rural_class",
        "column": "dz_urban_rural_class",
        "reference": "Large Urban Areas",
        "label": "Urban/rural class",
    },
    {
        "name": "health_board",
        "column": "dz_health_board",
        "reference": "Greater Glasgow and Clyde",
        "label": "Health board",
    },
)
COMPOSITION_FEATURES = {
    str(spec["column"]): spec["reference"] for spec in COMPOSITION_SPECS
}

STANDARDISE_SPECS = {
    "wn_prop_sequenced_z": "wn_prop_sequenced",
    "dz_cum_incidence_per_capita_z": "dz_cum_incidence_per_capita",
    "dz_cum_prop_sequenced_z": "dz_cum_prop_sequenced",
}
EPIDEMIC_CONTEXT_ADJUSTERS = [
    "wn_prop_sequenced_z",
    "dz_cum_incidence_per_capita_z",
    "dz_cum_prop_sequenced_z",
]

Domain = Literal["mixing", "composition"]
RegressionFamily = Literal["logistic", "linear"]


def required_regression_cluster_columns(
    score_outcomes: Sequence[str] = SCORE_OUTCOMES,
    cluster_id_col: str = CLUSTER_ID_COL,
) -> tuple[str, ...]:
    """Return cluster-table columns required by saved regression specifications."""
    observed_entropy_sources = [
        feature.removesuffix("_x10") for feature in OBSERVED_MIXING_FEATURES
    ]
    return tuple(
        _unique_preserve_order(
            [
                cluster_id_col,
                "candidate_tier",
                "cluster_size",
                *score_outcomes,
                *GROUP_VARS,
                *NULL_MIXING_FEATURES,
                *observed_entropy_sources,
                *STANDARDISE_SPECS.values(),
            ]
        )
    )


def validate_regression_cluster_columns(
    cluster_table: pd.DataFrame,
    *,
    score_outcomes: Sequence[str] = SCORE_OUTCOMES,
    cluster_id_col: str = CLUSTER_ID_COL,
) -> None:
    """Raise if an SSE cluster table cannot support the regression workflow."""
    required = required_regression_cluster_columns(
        score_outcomes=score_outcomes,
        cluster_id_col=cluster_id_col,
    )
    missing = [col for col in required if col not in cluster_table.columns]
    if missing:
        raise ValueError(
            "SSE cluster table is missing regression-required columns: "
            + ", ".join(missing)
        )


@dataclass(frozen=True)
class SampleSpec:
    """Sampling settings for a model frame.

    Use ``rows`` for a fixed maximum row count or ``fraction`` for a fraction of
    the complete-case frame.  Logistic frames may also request a target
    ``positive_fraction`` for case/control development samples.
    """

    rows: int | None = None
    fraction: float | None = None
    positive_fraction: float | None = None
    random_state: int = RANDOM_SEED

    def target_rows(self, n_rows: int) -> int | None:
        if self.rows is not None and self.fraction is not None:
            raise ValueError("SampleSpec accepts either rows or fraction, not both.")
        if self.rows is not None:
            if self.rows < 1:
                raise ValueError("SampleSpec.rows must be at least 1.")
            return min(int(self.rows), n_rows)
        if self.fraction is not None:
            if not 0 < self.fraction <= 1:
                raise ValueError("SampleSpec.fraction must be in (0, 1].")
            return max(1, min(n_rows, int(round(n_rows * self.fraction))))
        return None


@dataclass(frozen=True)
class ModelSpec:
    """One candidate model definition before fitting."""

    domain: Domain
    model_set: str
    predictor: str
    predictors: tuple[str, ...]
    terms: tuple[str, ...]


@dataclass
class RegressionDataBundle:
    """Aligned eligible data for node-level and sequence-level regressions."""

    eligible_nodes: pd.DataFrame
    eligible_sequence_data: pd.DataFrame
    eligibility_summary: pd.DataFrame
    min_candidate_size: int


@dataclass
class PreparedModelFrame:
    """Complete-case and sampled frame for one model configuration."""

    family: RegressionFamily
    domain: Domain
    outcome: str
    model_set: str
    predictor: str
    formula: str
    full_df: pd.DataFrame
    fit_df: pd.DataFrame
    output_dir: Path

    @property
    def result_key(self) -> str:
        if self.family == "logistic" and self.outcome == "candidate":
            return self.model_set
        return f"{self.outcome}_{self.model_set}"

    def grid_row(self) -> dict[str, object]:
        return {
            "family": self.family,
            "domain": self.domain,
            "outcome": self.outcome,
            "model_set": self.model_set,
            "predictor": self.predictor,
            "formula": self.formula,
            "model_dir": str(self.output_dir),
        }

    def summary_row(self) -> dict[str, object]:
        return model_frame_summary_row(
            self.result_key,
            self.full_df,
            self.fit_df,
            outcome=self.outcome,
            family=self.family,
            domain=self.domain,
            model_set=self.model_set,
            model_dir=self.output_dir,
        )


@dataclass
class PreparedRegressionRun:
    """Prepared model frames and run-level tables."""

    family: RegressionFamily
    result_dir: Path
    frames: dict[tuple[Domain, str, str], PreparedModelFrame]
    model_grid: pd.DataFrame
    fit_frame_summary: pd.DataFrame
    run_config: pd.DataFrame

    def select(
        self,
        *,
        domain: Domain,
        model_set: str,
        outcome: str | None = None,
    ) -> PreparedModelFrame:
        """Return a prepared model frame for a domain/model/outcome."""
        if outcome is None:
            outcome = "candidate" if self.family == "logistic" else None
        if outcome is None:
            matches = [
                frame
                for (frame_domain, _, frame_model_set), frame in self.frames.items()
                if frame_domain == domain and frame_model_set == model_set
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Please specify outcome; selection matched {len(matches)} frames."
                )
            return matches[0]
        return self.frames[(domain, outcome, model_set)]


def add_standardised_adjusters(data: pd.DataFrame) -> pd.DataFrame:
    """Add observed entropy scales and standardised surveillance adjusters."""
    out = add_observed_mixing_entropy_scales(data)
    for target, source in STANDARDISE_SPECS.items():
        if source not in out.columns:
            continue
        values = out[source].astype(float)
        sd = values.std(skipna=True)
        out[target] = (
            np.nan
            if pd.isna(sd) or sd == 0
            else (values - values.mean(skipna=True)) / sd
        )
    return out


def prepare_regression_data(
    sse_outputs: SseOutputs,
    *,
    sequence_data: pd.DataFrame | None = None,
    score_outcomes: Sequence[str] = SCORE_OUTCOMES,
    cluster_id_col: str = CLUSTER_ID_COL,
) -> RegressionDataBundle:
    """Load and align eligible node and sequence-level regression rows."""
    validate_regression_cluster_columns(
        sse_outputs.cluster_table,
        score_outcomes=score_outcomes,
        cluster_id_col=cluster_id_col,
    )
    cluster_data = add_standardised_adjusters(sse_outputs.cluster_table.copy())
    if sequence_data is None:
        from ..sse.detection import load_sequence_data as _load_sequence_data

        sequence_data = _load_sequence_data()
    sequence_data = add_standardised_adjusters(sequence_data)

    cluster_data["candidate"] = cluster_data["candidate_tier"].isin(
        HIGH_PRIORITY_CANDIDATE_TIERS
    )
    candidate_sizes = cluster_data.loc[
        cluster_data["candidate"], "cluster_size"
    ].dropna()
    if candidate_sizes.empty:
        raise ValueError("No high-priority candidate nodes were found.")

    missing_scores = [col for col in score_outcomes if col not in cluster_data.columns]
    if missing_scores:
        raise ValueError(f"Missing score columns in cluster table: {missing_scores}")

    min_candidate_size = int(candidate_sizes.min())
    eligible_nodes = cluster_data.loc[
        cluster_data["cluster_size"].ge(min_candidate_size)
    ].copy()

    merge_cols = _unique_preserve_order(
        [
            cluster_id_col,
            "candidate",
            *score_outcomes,
            *[
                col
                for col in ("burden_eligible", "candidate_tier")
                if col in eligible_nodes.columns
            ],
        ]
    )
    eligible_sequence_data = sequence_data.merge(
        eligible_nodes.loc[:, merge_cols],
        on=cluster_id_col,
        how="inner",
    )

    summary_rows = []
    for label, df in (
        ("eligible_nodes", eligible_nodes),
        ("eligible_sequence_data", eligible_sequence_data),
    ):
        row: dict[str, object] = {
            "dataset": label,
            "rows": len(df),
            "candidate_rate": float(df["candidate"].mean()),
            "candidates": int(df["candidate"].sum()),
        }
        for outcome in score_outcomes:
            row[f"{outcome}_nonmissing"] = int(df[outcome].notna().sum())
        summary_rows.append(row)

    return RegressionDataBundle(
        eligible_nodes=eligible_nodes,
        eligible_sequence_data=eligible_sequence_data,
        eligibility_summary=pd.DataFrame(summary_rows),
        min_candidate_size=min_candidate_size,
    )


def prepare_regression_run(
    data: RegressionDataBundle,
    *,
    family: RegressionFamily,
    result_dir: Path | str,
    outcomes: Sequence[str] | None = None,
    domains: Sequence[Domain] = ("mixing", "composition"),
    group_vars: Sequence[str] = GROUP_VARS,
    mixing_sample: SampleSpec | None = None,
    composition_sample: SampleSpec | None = None,
    write_tables: bool = False,
) -> PreparedRegressionRun:
    """Build all requested complete-case frames, sampled frames, and formulas."""
    result_dir = Path(result_dir)
    if outcomes is None:
        outcomes = ("candidate",) if family == "logistic" else SCORE_OUTCOMES

    frames: dict[tuple[Domain, str, str], PreparedModelFrame] = {}
    for domain in domains:
        sample = mixing_sample if domain == "mixing" else composition_sample
        for outcome in outcomes:
            for spec in default_model_specs(domain):
                base_df = (
                    data.eligible_nodes
                    if domain == "mixing"
                    else data.eligible_sequence_data
                )
                categorical_vars = (
                    (*group_vars, *COMPOSITION_FEATURES)
                    if domain == "composition"
                    else tuple(group_vars)
                )
                full_df = get_complete_case_data(
                    base_df,
                    outcome=outcome,
                    predictors=spec.predictors,
                    group_vars=group_vars,
                    categorical_vars=categorical_vars,
                    outcome_kind="binary" if family == "logistic" else "continuous",
                )
                fit_df = make_model_fit_data(
                    full_df,
                    family=family,
                    outcome=outcome,
                    sample=sample,
                    categorical_vars=tuple(COMPOSITION_FEATURES)
                    if domain == "composition"
                    else (),
                    compact_category_vars=group_vars,
                )
                formula = formula_with_varying_intercepts(
                    outcome=outcome,
                    terms=spec.terms,
                    group_vars=group_vars,
                )
                model_dir = model_output_dir(
                    result_dir,
                    domain=domain,
                    outcome=outcome,
                    model_set=spec.model_set,
                    include_outcome=family == "linear",
                )
                frame = PreparedModelFrame(
                    family=family,
                    domain=domain,
                    outcome=outcome,
                    model_set=spec.model_set,
                    predictor=spec.predictor,
                    formula=formula,
                    full_df=full_df,
                    fit_df=fit_df,
                    output_dir=model_dir,
                )
                frames[(domain, outcome, spec.model_set)] = frame

    model_grid = pd.DataFrame([frame.grid_row() for frame in frames.values()])
    fit_summary = pd.DataFrame([frame.summary_row() for frame in frames.values()])
    run_config = run_config_table(
        family=family,
        result_dir=result_dir,
        outcomes=outcomes,
        domains=domains,
        group_vars=group_vars,
        mixing_sample=mixing_sample,
        composition_sample=composition_sample,
    )

    prepared = PreparedRegressionRun(
        family=family,
        result_dir=result_dir,
        frames=frames,
        model_grid=model_grid,
        fit_frame_summary=fit_summary,
        run_config=run_config,
    )
    if write_tables:
        write_prepared_run_tables(prepared)
    return prepared


def default_model_specs(domain: Domain) -> tuple[ModelSpec, ...]:
    """Return primary and expanded model specs for a regression domain."""
    if domain == "mixing":
        return (
            ModelSpec(
                domain="mixing",
                model_set="null_primary",
                predictor="null_predictors",
                predictors=tuple(NULL_MIXING_FEATURES),
                terms=tuple(NULL_MIXING_FEATURES),
            ),
            ModelSpec(
                domain="mixing",
                model_set="null_expanded",
                predictor="null_predictors_plus_context",
                predictors=tuple([*NULL_MIXING_FEATURES, *EPIDEMIC_CONTEXT_ADJUSTERS]),
                terms=tuple([*NULL_MIXING_FEATURES, *EPIDEMIC_CONTEXT_ADJUSTERS]),
            ),
            ModelSpec(
                domain="mixing",
                model_set="observed_primary",
                predictor="observed_predictors",
                predictors=tuple(OBSERVED_MIXING_FEATURES),
                terms=tuple(OBSERVED_MIXING_FEATURES),
            ),
            ModelSpec(
                domain="mixing",
                model_set="observed_expanded",
                predictor="observed_predictors_plus_context",
                predictors=tuple(
                    [*OBSERVED_MIXING_FEATURES, *EPIDEMIC_CONTEXT_ADJUSTERS]
                ),
                terms=tuple([*OBSERVED_MIXING_FEATURES, *EPIDEMIC_CONTEXT_ADJUSTERS]),
            ),
        )

    composition_terms = tuple(
        treatment_term(column, reference)
        for column, reference in COMPOSITION_FEATURES.items()
    )
    return (
        ModelSpec(
            domain="composition",
            model_set="primary",
            predictor="composition_predictors",
            predictors=tuple(COMPOSITION_FEATURES),
            terms=composition_terms,
        ),
        ModelSpec(
            domain="composition",
            model_set="expanded",
            predictor="composition_predictors_plus_context",
            predictors=tuple([*COMPOSITION_FEATURES, *EPIDEMIC_CONTEXT_ADJUSTERS]),
            terms=tuple([*composition_terms, *EPIDEMIC_CONTEXT_ADJUSTERS]),
        ),
    )


def get_complete_case_data(
    df: pd.DataFrame,
    *,
    outcome: str,
    predictors: Iterable[str],
    group_vars: Sequence[str] = GROUP_VARS,
    categorical_vars: Sequence[str] = (),
    id_cols: Sequence[str] = (CLUSTER_ID_COL,),
    outcome_kind: Literal["binary", "continuous"] = "continuous",
) -> pd.DataFrame:
    """Create a complete-case model frame with stable categorical dtypes."""
    required_cols = _unique_preserve_order(
        [outcome, *id_cols, *predictors, *group_vars]
    )
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in dataframe: {missing_cols}")

    out = df.loc[:, required_cols].dropna().copy()
    if outcome_kind == "binary":
        out[outcome] = out[outcome].astype(int)
    else:
        out[outcome] = out[outcome].astype(float)

    categorical_cols = _unique_preserve_order(
        [*categorical_vars, *group_vars, *id_cols]
    )
    for col in categorical_cols:
        if col in out.columns:
            out[col] = pd.Categorical(out[col], categories=_categories_from(df[col]))
    return out


def make_model_fit_data(
    data: pd.DataFrame,
    *,
    family: RegressionFamily,
    outcome: str,
    sample: SampleSpec | None = None,
    categorical_vars: Sequence[str] = (),
    compact_category_vars: Sequence[str] = (),
) -> pd.DataFrame:
    """Return either a full complete-case frame or a sampled fit frame."""
    if sample is None:
        return data.copy()

    max_rows = sample.target_rows(len(data))
    if max_rows is None or max_rows >= len(data):
        return data.copy().reset_index(drop=True)

    if family == "logistic":
        sampled = sample_binary_outcome_data(
            data,
            outcome=outcome,
            max_rows=max_rows,
            positive_fraction=sample.positive_fraction,
            random_state=sample.random_state,
            categorical_vars=categorical_vars,
        )
    else:
        sampled = sample_rows_preserving_categories(
            data,
            max_rows=max_rows,
            random_state=sample.random_state,
            categorical_vars=categorical_vars,
        )
    return compact_categorical_levels(sampled, compact_category_vars)


def sample_rows_preserving_categories(
    data: pd.DataFrame,
    *,
    max_rows: int,
    random_state: int = RANDOM_SEED,
    categorical_vars: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Sample rows while seeding at least one row per observed category level."""
    if len(data) <= max_rows:
        return data.copy().reset_index(drop=True)

    category_cols = [
        col for col in _unique_preserve_order(categorical_vars or ()) if col in data
    ]
    selected_indices: set[object] = set()
    seed_indices = _category_seed_indices(
        data, category_cols, random_state=random_state
    )
    if len(seed_indices) > max_rows:
        raise ValueError(
            "max_rows is too small to include all observed categorical levels "
            f"({len(seed_indices)} seed rows needed; max_rows={max_rows})."
        )
    selected_indices.update(seed_indices)

    remaining = data.drop(index=list(selected_indices), errors="ignore")
    sampled_parts = [data.loc[seed_indices]] if seed_indices else []
    sampled_parts.append(
        remaining.sample(
            n=max_rows - len(seed_indices),
            random_state=random_state,
            replace=False,
        )
    )
    out = pd.concat(sampled_parts, axis=0).sample(frac=1, random_state=random_state)
    return _restore_categories(out, data, category_cols)


def sample_binary_outcome_data(
    data: pd.DataFrame,
    *,
    outcome: str,
    max_rows: int,
    positive_fraction: float | None = None,
    random_state: int = RANDOM_SEED,
    categorical_vars: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Sample binary data with optional target positive fraction."""
    positives = data.loc[data[outcome] == 1]
    negatives = data.loc[data[outcome] == 0]
    if positives.empty or negatives.empty:
        raise ValueError(f"'{outcome}' must contain at least one 0 and one 1.")

    if positive_fraction is None:
        positive_fraction = float(data[outcome].mean())
    positive_fraction = float(np.clip(positive_fraction, 0.01, 0.99))

    category_cols = [
        col for col in _unique_preserve_order(categorical_vars or ()) if col in data
    ]
    seed_indices = _category_seed_indices(
        data, category_cols, random_state=random_state
    )
    if len(seed_indices) > max_rows:
        raise ValueError(
            "max_rows is too small to include all observed categorical levels "
            f"({len(seed_indices)} seed rows needed; max_rows={max_rows})."
        )

    selected_indices: set[object] = set(seed_indices)
    seed_outcomes = (
        data.loc[seed_indices, outcome] if seed_indices else pd.Series(dtype=int)
    )
    seed_pos = int((seed_outcomes == 1).sum())
    seed_neg = int((seed_outcomes == 0).sum())

    n_pos = min(len(positives), max(1, int(round(max_rows * positive_fraction))))
    n_neg = min(len(negatives), max_rows - n_pos)
    if n_neg <= 0:
        raise ValueError("Sampling settings left no room for negative controls.")

    remaining_capacity = max_rows - len(seed_indices)
    n_pos_fill = max(0, n_pos - seed_pos)
    n_neg_fill = max(0, n_neg - seed_neg)
    if n_pos_fill + n_neg_fill > remaining_capacity:
        scale = remaining_capacity / (n_pos_fill + n_neg_fill)
        n_pos_fill = int(round(n_pos_fill * scale))
        n_neg_fill = remaining_capacity - n_pos_fill

    def sample_available(frame: pd.DataFrame, n: int, seed_offset: int) -> pd.DataFrame:
        available = frame.drop(index=list(selected_indices), errors="ignore")
        n = min(n, len(available))
        if n <= 0:
            return available.iloc[0:0]
        out = available.sample(
            n=n,
            random_state=random_state + seed_offset,
            replace=False,
        )
        selected_indices.update(out.index.tolist())
        return out

    sampled_parts = [data.loc[seed_indices]] if seed_indices else []
    sampled_parts.append(sample_available(positives, n_pos_fill, 10_000))
    sampled_parts.append(sample_available(negatives, n_neg_fill, 20_000))
    sampled = pd.concat(sampled_parts, axis=0)

    if len(sampled) < max_rows:
        sampled = pd.concat(
            [
                sampled,
                sample_available(data, max_rows - len(sampled), 30_000),
            ],
            axis=0,
        )

    sampled = sampled.sample(frac=1, random_state=random_state)
    return _restore_categories(sampled, data, category_cols)


def treatment_term(variable: str, reference: str | int | float | None) -> str:
    """Return a formulae categorical term with an optional reference level."""
    if reference is None:
        return f"C({variable})"
    return f"C({variable}, Treatment(reference={reference!r}))"


def formula_with_varying_intercepts(
    outcome: str,
    terms: str | Sequence[str],
    group_vars: Sequence[str],
) -> str:
    """Build a fixed-effect formula plus group-level intercepts."""
    if isinstance(terms, str):
        terms = [terms]
    rhs_terms = [*terms, *(f"(1|{group})" for group in group_vars)]
    return f"{outcome} ~ " + " + ".join(rhs_terms)


def compact_categorical_levels(
    data: pd.DataFrame,
    columns: Sequence[str],
) -> pd.DataFrame:
    """Drop unobserved categories from categorical columns."""
    out = data.copy()
    for col in columns:
        if col not in out:
            continue
        if isinstance(out[col].dtype, pd.CategoricalDtype):
            out[col] = out[col].cat.remove_unused_categories()
        else:
            out[col] = pd.Categorical(out[col])
    return out


def model_frame_summary_row(
    label: str,
    full_df: pd.DataFrame,
    fit_df: pd.DataFrame,
    *,
    outcome: str,
    family: RegressionFamily,
    domain: Domain,
    model_set: str,
    model_dir: Path,
) -> dict[str, object]:
    """Summarise the full and fitting frames side by side."""
    row: dict[str, object] = {
        "frame": label,
        "family": family,
        "domain": domain,
        "outcome": outcome,
        "model_set": model_set,
        "model_dir": str(model_dir),
        "full_rows": len(full_df),
        "fit_rows": len(fit_df),
        "fit_fraction": len(fit_df) / len(full_df) if len(full_df) else np.nan,
        "use_sample": len(fit_df) != len(full_df),
    }
    if family == "logistic":
        row.update(
            {
                "full_candidates": int(full_df[outcome].sum()),
                "fit_candidates": int(fit_df[outcome].sum()),
                "full_candidate_rate": float(full_df[outcome].mean()),
                "fit_candidate_rate": float(fit_df[outcome].mean()),
            }
        )
    else:
        row.update(
            {
                "full_outcome_mean": float(full_df[outcome].mean()),
                "fit_outcome_mean": float(fit_df[outcome].mean()),
                "full_outcome_sd": float(full_df[outcome].std()),
                "fit_outcome_sd": float(fit_df[outcome].std()),
            }
        )
    return row


def default_result_dir(project_root: Path | str, *, family: RegressionFamily) -> Path:
    """Return the conventional result directory for a regression family."""
    name = (
        "bayesian_socio_geo_demo_logistic_regression"
        if family == "logistic"
        else "bayesian_socio_geo_demo_linear_regression"
    )
    relative_results = RESULTS_DIR.relative_to(PROJECT_ROOT)
    return Path(project_root) / relative_results / name


def model_output_dir(
    result_dir: Path | str,
    *,
    domain: Domain,
    outcome: str,
    model_set: str,
    include_outcome: bool,
) -> Path:
    """Build a deterministic output directory from model configuration."""
    parts = [Path(result_dir), Path(domain)]
    if include_outcome:
        parts.append(Path(_slug(outcome)))
    parts.append(Path(_slug(model_set)))
    return Path(*parts)


def model_output_files(model_dir: Path | str) -> dict[str, Path]:
    """Return conventional output file names for one fitted model."""
    model_dir = Path(model_dir)
    return {
        "summary": model_dir / "summary.csv",
        "diagnostics": model_dir / "diagnostics.csv",
        "metadata": model_dir / "metadata.csv",
        "idata": model_dir / "idata.nc",
    }


def write_prepared_run_tables(prepared: PreparedRegressionRun) -> None:
    """Write run config, model grid, and fit-frame summary tables."""
    prepared.result_dir.mkdir(parents=True, exist_ok=True)
    with exclusive_file_lock(prepared.result_dir / ".prepared_run_tables.lock"):
        atomic_write_csv(
            prepared.run_config,
            prepared.result_dir / "run_config.csv",
            index=False,
        )
        atomic_write_csv(
            prepared.model_grid,
            prepared.result_dir / "model_grid.csv",
            index=False,
        )
        atomic_write_csv(
            prepared.fit_frame_summary,
            prepared.result_dir / "fit_frame_summary.csv",
            index=False,
        )
        for domain, table in prepared.model_grid.groupby("domain", sort=False):
            atomic_write_csv(
                table,
                prepared.result_dir / f"{domain}_model_grid.csv",
                index=False,
            )
        for domain, table in prepared.fit_frame_summary.groupby("domain", sort=False):
            atomic_write_csv(
                table,
                prepared.result_dir / f"{domain}_fit_frame_summary.csv",
                index=False,
            )


def write_fit_frames(
    prepared: PreparedRegressionRun,
    *,
    file_name: str = "fit_frame.parquet",
) -> pd.DataFrame:
    """Write each prepared fitting frame below its model output directory."""
    manifest_rows = []
    with exclusive_file_lock(prepared.result_dir / ".fit_frames.lock"):
        for frame in prepared.frames.values():
            frame.output_dir.mkdir(parents=True, exist_ok=True)
            path = frame.output_dir / file_name
            atomic_write_parquet(frame.fit_df, path, index=False)
            manifest_rows.append(
                {
                    "family": frame.family,
                    "domain": frame.domain,
                    "outcome": frame.outcome,
                    "model_set": frame.model_set,
                    "fit_frame_path": str(path),
                    "fit_rows": len(frame.fit_df),
                }
            )
        manifest = pd.DataFrame(manifest_rows)
        atomic_write_csv(
            manifest,
            prepared.result_dir / "fit_frame_manifest.csv",
            index=False,
        )
    return manifest


def run_config_table(
    *,
    family: RegressionFamily,
    result_dir: Path,
    outcomes: Sequence[str],
    domains: Sequence[Domain],
    group_vars: Sequence[str],
    mixing_sample: SampleSpec | None,
    composition_sample: SampleSpec | None,
) -> pd.DataFrame:
    """Represent run-level settings as a one-row table."""
    return pd.DataFrame(
        [
            {
                "family": family,
                "result_dir": str(result_dir),
                "outcomes": ",".join(outcomes),
                "domains": ",".join(domains),
                "group_vars": ",".join(group_vars),
                "composition_predictors": ",".join(COMPOSITION_FEATURES),
                "null_mixing_predictors": ",".join(NULL_MIXING_FEATURES),
                "observed_mixing_predictors": ",".join(OBSERVED_MIXING_FEATURES),
                **_sample_config("mixing", mixing_sample),
                **_sample_config("composition", composition_sample),
            }
        ]
    )


def _sample_config(prefix: str, sample: SampleSpec | None) -> dict[str, object]:
    if sample is None:
        return {
            f"use_{prefix}_sample": False,
            f"{prefix}_sample_rows": np.nan,
            f"{prefix}_sample_fraction": np.nan,
            f"{prefix}_positive_fraction": np.nan,
        }
    return {
        f"use_{prefix}_sample": True,
        f"{prefix}_sample_rows": sample.rows,
        f"{prefix}_sample_fraction": sample.fraction,
        f"{prefix}_positive_fraction": sample.positive_fraction,
    }


def _unique_preserve_order(items: Iterable[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for item in items:
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


def _categories_from(series: pd.Series) -> list:
    if isinstance(series.dtype, pd.CategoricalDtype):
        return list(series.cat.categories)
    return series.dropna().drop_duplicates().tolist()


def _category_seed_indices(
    data: pd.DataFrame,
    category_cols: Sequence[str],
    *,
    random_state: int,
) -> list[Any]:
    missing_by_col = {
        col: {
            level for level in _categories_from(data[col]) if data[col].eq(level).any()
        }
        for col in category_cols
    }
    seed_indices: list[Any] = []
    selected_indices: set[Any] = set()

    for col, levels in missing_by_col.items():
        for level in list(levels):
            if level not in missing_by_col[col]:
                continue
            candidates = data.loc[data[col].eq(level)]
            coverage = pd.Series(0, index=candidates.index, dtype=int)
            for cover_col, missing_levels in missing_by_col.items():
                if missing_levels:
                    coverage += candidates[cover_col].isin(missing_levels).astype(int)
            best_indices = coverage.loc[coverage.eq(coverage.max())].index
            idx = (
                data.loc[best_indices]
                .sample(
                    n=1,
                    random_state=random_state + len(seed_indices),
                    replace=False,
                )
                .index[0]
            )
            if idx not in selected_indices:
                seed_indices.append(idx)
                selected_indices.add(idx)
            row = data.loc[idx]
            for cover_col in category_cols:
                missing_by_col[cover_col].discard(row[cover_col])
    return seed_indices


def _restore_categories(
    sampled: pd.DataFrame,
    original: pd.DataFrame,
    category_cols: Sequence[str],
) -> pd.DataFrame:
    out = sampled.reset_index(drop=True)
    for col in category_cols:
        out[col] = pd.Categorical(out[col], categories=_categories_from(original[col]))
    return out


def _slug(value: object) -> str:
    text = str(value).strip().lower()
    chars = [char if char.isalnum() else "_" for char in text]
    slug = "_".join("".join(chars).split("_"))
    return slug or "model"
