"""Parameterized association-regression runner for SSE sensitivity notebooks."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from utils.data import CLADES, load_analysis_columns

from .io import load_sse_outputs
from .regression import (
    AssociationModel,
    bh_adjust,
    categorical_term,
    cluster_se_diagnostics,
    fit_binomial_glm,
    fit_conditional_logit,
    fit_firth_logit,
    make_formula,
    model_fit_stats,
    model_variables_from_terms,
    parameter_names_for_term,
    robust_wald_for_prefix,
    tidy_odds_ratios,
    tidy_single_parameter_wald,
)


__all__ = [
    "COMPOSITION_SPECS",
    "DEFAULT_MIXING_FEATURES",
    "OBSERVED_MIXING_FEATURES",
    "OBSERVED_MIXING_FEATURES_X10",
    "OBSERVED_MIXING_REFERENCE_X10",
    "WINDOW_SURVEILLANCE_ADJUSTERS",
    "EXPANDED_CONTEXT_ADJUSTERS",
    "TABLE_DISPLAY_COLUMNS",
    "AssociationFrames",
    "add_clade_group",
    "add_standardised_adjusters",
    "default_model_sets",
    "load_association_frames",
    "run_association_pipeline",
    "run_main_association_analysis",
    "select_table_columns",
]


WINDOW_SURVEILLANCE_ADJUSTERS = [
    "z_wn_prop_sequenced",
    "z_log1p_wn_positive_tests",
]

EXPANDED_CONTEXT_ADJUSTERS = [
    "z_dz_cum_prop_sequenced",
    "z_dz_cum_incidence_per_capita",
    "z_dz_7d_test_positivity",
    "z_log1p_dz_cum_positive_tests",
]

COMPOSITION_SPECS = [
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
        "fallback_references": ["25-29", "30-34"],
        "label": "Age band",
    },
    {
        "name": "simd_quintile",
        "column": "dz_simd_quintile",
        "reference": "1",
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
]

DEFAULT_MIXING_FEATURES = [
    "sex_entropy_z",
    "age_entropy_z",
    "simd_entropy_z",
    "urban_rural_entropy_z",
    "health_board_entropy_z",
]

OBSERVED_MIXING_FEATURES = [
    "sex_entropy_obs",
    "age_entropy_obs",
    "simd_entropy_obs",
    "urban_rural_entropy_obs",
    "health_board_entropy_obs",
]
OBSERVED_MIXING_FEATURES_X10 = [
    f"{feature}_x10" for feature in OBSERVED_MIXING_FEATURES
]
OBSERVED_MIXING_REFERENCE_X10 = "per 0.1 increase in observed normalised entropy"

STANDARDISE_SPECS = {
    "z_wn_prop_sequenced": "wn_prop_sequenced",
    "z_log1p_wn_positive_tests": "log1p_wn_positive_tests",
    "z_dz_cum_prop_sequenced": "dz_cum_prop_sequenced",
    "z_dz_cum_incidence_per_capita": "dz_cum_incidence_per_capita",
    "z_dz_7d_test_positivity": "dz_7d_test_positivity",
    "z_log1p_dz_cum_positive_tests": "log1p_dz_cum_positive_tests",
}

TABLE_DISPLAY_COLUMNS = {
    "wald": [
        "domain",
        "model_set",
        "predictor_set",
        "predictor",
        "label",
        "reference",
        "term",
        "chi2",
        "df",
        "P>chi2",
        "p_adj_bh",
        "n_model_rows",
        "n_sequences",
        "n_nodes",
        "dropped_nonvarying_rows",
        "dropped_nonvarying_strata",
        "dropped_nonvarying_detail",
    ],
    "odds_ratios": [
        "domain",
        "model_set",
        "predictor_set",
        "predictor",
        "label",
        "reference",
        "term",
        "estimate",
        "std_error",
        "p_value",
        "odds_ratio",
        "or_low",
        "or_high",
    ],
    "fit_stats": [
        "domain",
        "model_set",
        "predictor_set",
        "predictor",
        "r2_mcfadden",
        "converged",
        "aic",
        "bic_llf",
        "log_likelihood",
        "ll_null",
        "n_model_rows",
        "n_sequences",
        "n_nodes",
    ],
}


@dataclass
class AssociationFrames:
    """Prepared node and sequence-level model frames."""

    node_stats: pd.DataFrame
    eligible_nodes: pd.DataFrame
    node_model_base: pd.DataFrame
    composition_base: pd.DataFrame
    cluster_diagnostics: pd.DataFrame
    min_candidate_size: int


def _project_root(project_root: Path | str | None = None) -> Path:
    root = Path(project_root or Path.cwd()).resolve()
    candidates = [root, *root.parents]
    for candidate in candidates:
        if (candidate / "config.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate config.yaml from project_root.")


def _dedupe(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))


def _slug(value: object) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_").lower()
    return slug or "all"


def _concat_or_empty(tables: list[pd.DataFrame]) -> pd.DataFrame:
    tables = [table for table in tables if table is not None and not table.empty]
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def _has_categorical_adjuster(adjusters: Sequence[str], variable: str) -> bool:
    pattern = re.compile(rf"^C\(\s*{re.escape(variable)}\s*(?:,|\))")
    return any(pattern.search(term) for term in adjusters)


def _model_name(base: str, group_label: object | None) -> str:
    if group_label is None:
        return base
    return f"{base}__{_slug(group_label)}"


def _group_metadata(
    group_col: str | None, group_label: object | None
) -> dict[str, object]:
    if group_col is None:
        return {}
    return {"analysis_group_col": group_col, group_col: group_label}


def _iter_groups(
    data: pd.DataFrame,
    group_col: str | None,
    group_values: Sequence[object] | None,
    *,
    outcome: str = "candidate",
) -> Iterable[tuple[object | None, pd.DataFrame]]:
    if group_col is None:
        if data[outcome].nunique(dropna=True) >= 2:
            yield None, data
        return

    observed = set(data[group_col].dropna().unique())
    ordered = [value for value in (group_values or []) if value in observed]
    extras = sorted(observed - set(ordered), key=str)

    for value in [*ordered, *extras]:
        d = data.loc[data[group_col].eq(value)].copy()
        if d.empty:
            continue
        if d[outcome].nunique(dropna=True) < 2:
            print(f"Skipping {group_col}={value!r}: outcome does not vary", flush=True)
            continue
        yield value, d


def add_clade_group(
    data: pd.DataFrame,
    *,
    source_col: str = "clade",
    target_col: str = "clade_group",
) -> pd.DataFrame:
    """Map Nextclade clades onto the curated top-clade labels plus Other."""
    if source_col not in data.columns:
        raise KeyError(f"{source_col!r} is required to build clade groups.")
    out = data.copy()
    out[target_col] = out[source_col].map(CLADES).fillna("Other")
    return out


def add_standardised_adjusters(data: pd.DataFrame) -> pd.DataFrame:
    """Add standardised surveillance and context adjusters used by notebooks."""
    out = data.copy()
    for feature, scaled_feature in zip(
        OBSERVED_MIXING_FEATURES,
        OBSERVED_MIXING_FEATURES_X10,
    ):
        if feature in out.columns:
            out[scaled_feature] = out[feature].astype(float) * 10
    if "wn_positive_tests" in out.columns:
        out["log1p_wn_positive_tests"] = np.log1p(out["wn_positive_tests"])
    if "dz_cum_positive_tests" in out.columns:
        out["log1p_dz_cum_positive_tests"] = np.log1p(out["dz_cum_positive_tests"])
    for target, source in STANDARDISE_SPECS.items():
        if source not in out.columns:
            continue
        values = out[source].astype(float)
        sd = values.std(skipna=True)
        if pd.isna(sd) or sd == 0:
            out[target] = np.nan
        else:
            out[target] = (values - values.mean(skipna=True)) / sd
    return out


def default_model_sets(
    *,
    variant_adjuster: str | None = "clade",
    window_adjustment: str = "fixed_effects",
) -> dict[str, list[str]]:
    """Build primary and expanded adjustment sets for association models."""
    if window_adjustment == "fixed_effects":
        window_terms = ["C(window_idx)"]
    elif window_adjustment == "surveillance":
        window_terms = list(WINDOW_SURVEILLANCE_ADJUSTERS)
    else:
        raise ValueError("window_adjustment must be 'fixed_effects' or 'surveillance'.")

    variant_terms = [f"C({variant_adjuster})"] if variant_adjuster else []
    primary = [*window_terms, *variant_terms]
    return {
        "primary": primary,
        "expanded": [*primary, *EXPANDED_CONTEXT_ADJUSTERS],
    }


def load_association_frames(
    *,
    project_root: Path | str | None = None,
    output_dir: Path | str | None = None,
    cluster_se: str = "cluster_id",
    variant_adjuster: str | None = "clade",
    group_by_clade: bool = False,
    clade_group_col: str = "clade_group",
    window_stride: int = 2,
    run_composition: bool = True,
) -> AssociationFrames:
    """Load SSE outputs and construct complete model base frames."""
    root = _project_root(project_root)
    output_path = (
        Path(output_dir)
        if output_dir
        else root / "sse_detection" / "results" / "sse_outputs"
    )
    outs = load_sse_outputs(output_path)
    node_stats = outs.node_stats.copy()

    if group_by_clade:
        node_stats = add_clade_group(node_stats, target_col=clade_group_col)
        if variant_adjuster == "clade":
            variant_adjuster = None

    min_candidate_size = int(
        node_stats.loc[node_stats["sse_candidate"], "cluster_size"].min()
    )
    eligible_nodes = node_stats.loc[
        node_stats["cluster_size"].ge(min_candidate_size)
    ].copy()

    if cluster_se not in eligible_nodes.columns:
        raise KeyError(f"{cluster_se!r} is not present in node_stats.")
    if variant_adjuster and variant_adjuster not in eligible_nodes.columns:
        raise KeyError(f"{variant_adjuster!r} is not present in node_stats.")

    node_model_base = eligible_nodes.copy()
    node_model_base["candidate"] = node_model_base["sse_candidate"].astype(int)
    if variant_adjuster:
        node_model_base[variant_adjuster] = (
            node_model_base[variant_adjuster].fillna("Missing").astype(str)
        )
    node_model_base[cluster_se] = (
        node_model_base[cluster_se].fillna(node_model_base["cluster_id"]).astype(str)
    )
    node_model_base = add_standardised_adjusters(node_model_base)

    composition_base = pd.DataFrame()
    diagnostics = [
        cluster_se_diagnostics(
            node_model_base,
            cluster_se,
            outcome="candidate",
        ).assign(analysis_frame="node_mixing")
    ]

    if run_composition:
        node_key_cols = [
            "cluster_id",
            "meta_cluster_id",
            "sse_candidate",
            cluster_se,
            "cluster_size",
        ]
        if group_by_clade:
            node_key_cols.append(clade_group_col)
        node_key = eligible_nodes[_dedupe(node_key_cols)].drop_duplicates("cluster_id")

        sequence_columns: list[str] = list(
            {
                "window_id",
                "window_idx",
                "cluster_id",
                "sequence_id",
                "wn_prop_sequenced",
                "wn_positive_tests",
                "dz_cum_prop_sequenced",
                "dz_cum_incidence_per_capita",
                "dz_7d_test_positivity",
                "dz_cum_positive_tests",
                *(spec["column"] for spec in COMPOSITION_SPECS),
                *([variant_adjuster] if variant_adjuster else []),
            }
        )

        sequence_raw = load_analysis_columns(
            sequence_columns,
            add_policy=False,
            window_stride=window_stride,
        )
        composition_base = sequence_raw.merge(node_key, on="cluster_id", how="inner")
        composition_base["candidate"] = composition_base["sse_candidate"].astype(int)
        if variant_adjuster:
            composition_base[variant_adjuster] = (
                composition_base[variant_adjuster].fillna("Missing").astype(str)
            )
        composition_base[cluster_se] = (
            composition_base[cluster_se]
            .fillna(composition_base["cluster_id"])
            .astype(str)
        )
        composition_base = add_standardised_adjusters(composition_base)
        diagnostics.insert(
            0,
            cluster_se_diagnostics(
                composition_base,
                cluster_se,
                outcome="candidate",
            ).assign(analysis_frame="composition"),
        )

    return AssociationFrames(
        node_stats=node_stats,
        eligible_nodes=eligible_nodes,
        node_model_base=node_model_base,
        composition_base=composition_base,
        cluster_diagnostics=pd.concat(diagnostics, ignore_index=True),
        min_candidate_size=min_candidate_size,
    )


def fit_association_result(
    data: pd.DataFrame,
    formula: str,
    *,
    model_method: str,
    cluster_se: str,
    window_strata: str,
):
    if model_method == "conditional_logit_by_window":
        return fit_conditional_logit(data, formula, strata_col=window_strata)
    if model_method == "firth_glm":
        return fit_firth_logit(data, formula)
    if model_method == "glm_clustered":
        return fit_binomial_glm(data, formula, cluster_col=cluster_se)
    raise ValueError(f"Unknown model_method={model_method!r}")


def fit_exposure_association(
    data: pd.DataFrame,
    *,
    outcome: str,
    exposure: str,
    adjusters: list[str],
    model_name: str,
    model_method: str,
    cluster_se: str,
    window_strata: str,
    reference=None,
    categorical: bool = True,
) -> AssociationModel:
    exposure_term = (
        categorical_term(exposure, reference=reference) if categorical else exposure
    )
    formula = make_formula(outcome, exposure_term, adjusters)
    result = fit_association_result(
        data,
        formula,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
    )
    odds = tidy_odds_ratios(
        result,
        model_name=model_name,
        term_filter=exposure_term,
    )
    if categorical:
        wald = robust_wald_for_prefix(
            result,
            exposure_term,
            model_name=model_name,
            term=exposure,
        )
    else:
        wald = tidy_single_parameter_wald(
            result,
            [exposure],
            model_name=model_name,
        )
    wald["formula"] = formula
    return AssociationModel(result=result, odds_ratios=odds, wald=wald, formula=formula)


def resolve_reference(data: pd.DataFrame, column: str, preferred, fallbacks=None):
    levels = set(data[column].dropna().astype(str))
    candidates = [preferred, *(fallbacks or [])]
    for ref in candidates:
        if ref is None:
            continue
        ref_str = str(ref)
        if ref_str in levels:
            return ref_str
        for level in levels:
            try:
                if float(level) == float(ref_str):
                    return level
            except ValueError:
                pass
    counts = data[column].dropna().astype(str).value_counts()
    if counts.empty:
        raise ValueError(f"No observed levels for {column!r}.")
    return str(counts.index[0])


def complete_case(data: pd.DataFrame, required: list[str]) -> pd.DataFrame:
    missing = [col for col in required if col not in data.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")
    return data.dropna(subset=required).copy()


def drop_nonvarying_levels(
    data: pd.DataFrame,
    columns: list[str],
    *,
    outcome: str = "candidate",
) -> tuple[pd.DataFrame, dict[str, int], dict[str, int]]:
    d = data.copy()
    dropped_rows: dict[str, int] = {}
    dropped_strata: dict[str, int] = {}
    changed = True
    while changed:
        changed = False
        for col in columns:
            if col not in d.columns:
                continue
            strata_nunique = d.groupby(col, dropna=False)[outcome].transform("nunique")
            varies = strata_nunique.gt(1)
            n_drop = int((~varies).sum())
            if n_drop:
                dropped_rows[col] = dropped_rows.get(col, 0) + n_drop
                dropped_strata[col] = dropped_strata.get(col, 0) + int(
                    d.loc[~varies, col].nunique(dropna=False)
                )
                d = d.loc[varies].copy()
                changed = True
    return d, dropped_rows, dropped_strata


def dropped_metadata(
    dropped_rows: dict[str, int], dropped_strata: dict[str, int]
) -> dict[str, object]:
    return {
        "dropped_nonvarying_rows": sum(dropped_rows.values()),
        "dropped_nonvarying_strata": sum(dropped_strata.values()),
        "dropped_nonvarying_detail": repr(
            {
                "rows": dropped_rows,
                "strata": dropped_strata,
            }
        ),
    }


def add_model_metadata(table: pd.DataFrame, **metadata) -> pd.DataFrame:
    out = table.copy()
    for key, value in reversed(list(metadata.items())):
        out.insert(0, key, value)
    return out


def add_fit_metadata(fit_stats: pd.DataFrame, **metadata) -> pd.DataFrame:
    out = fit_stats.copy()
    for key, value in reversed(list(metadata.items())):
        out.insert(0, key, value)
    return out


def bh_adjust_by(
    table: pd.DataFrame,
    group_cols: list[str],
    p_col: str = "P>chi2",
) -> pd.DataFrame:
    if table.empty:
        return table.copy()
    out = table.copy()
    out["p_adj_bh"] = np.nan
    for _, idx in out.groupby(group_cols, dropna=False).groups.items():
        adjusted = bh_adjust(out.loc[idx], p_col=p_col)
        out.loc[idx, "p_adj_bh"] = adjusted["p_adj_bh"].to_numpy()
    return out


def _prepare_model_frame(
    source: pd.DataFrame,
    *,
    predictors: list[str],
    adjusters: list[str],
    required_base: list[str],
    window_strata: str,
    drop_window_nonvarying: bool,
    categorical_predictors: bool = False,
) -> tuple[pd.DataFrame, dict[str, int], dict[str, int]]:
    required = required_base + predictors + model_variables_from_terms(adjusters)
    required = _dedupe(required)
    d = complete_case(source, required)
    if categorical_predictors:
        for col in predictors:
            d[col] = d[col].astype(str)
    strata = [window_strata] if drop_window_nonvarying else []
    d, dropped_rows, dropped_strata = drop_nonvarying_levels(d, strata)
    if d.empty:
        raise ValueError("No complete-case rows remain after filtering.")
    if d["candidate"].nunique(dropna=True) < 2:
        raise ValueError("Outcome does not vary after filtering.")
    return d, dropped_rows, dropped_strata


def _record_failure(
    failures: list[dict[str, object]],
    *,
    domain: str,
    model_set: str,
    predictor_set: str,
    predictor: str,
    group_col: str | None,
    group_label: object | None,
    error: Exception,
) -> None:
    row = {
        "domain": domain,
        "model_set": model_set,
        "predictor_set": predictor_set,
        "predictor": predictor,
        "error_type": type(error).__name__,
        "error": str(error),
    }
    row |= _group_metadata(group_col, group_label)
    failures.append(row)
    where = f", {group_col}={group_label!r}" if group_col else ""
    print(
        f"Failed {domain} {model_set} {predictor_set} {predictor}{where}: {error}",
        flush=True,
    )


def fit_single_composition_models(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    wald_tables = []
    or_tables = []
    fit_tables = []
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )

    for spec in COMPOSITION_SPECS:
        predictor = spec["column"]
        try:
            d, dropped_rows, dropped_strata = _prepare_model_frame(
                source,
                predictors=[predictor],
                adjusters=adjusters,
                required_base=[
                    "candidate",
                    "cluster_id",
                    "sequence_id",
                    cluster_se,
                    window_strata,
                ],
                window_strata=window_strata,
                drop_window_nonvarying=drop_window,
                categorical_predictors=True,
            )
            reference = resolve_reference(
                d,
                predictor,
                spec.get("reference"),
                spec.get("fallback_references"),
            )
            base_name = f"composition__{model_set}__single__{spec['name']}"
            model_name = _model_name(base_name, group_label)
            fit = fit_exposure_association(
                d,
                outcome="candidate",
                exposure=predictor,
                adjusters=adjusters,
                model_name=model_name,
                model_method=model_method,
                cluster_se=cluster_se,
                window_strata=window_strata,
                reference=reference,
                categorical=True,
            )
            meta = {
                "domain": "composition",
                "model_set": model_set,
                "predictor_set": "single",
                "predictor": spec["name"],
                "label": spec["label"],
                "reference": reference,
                "n_model_rows": len(d),
                "n_sequences": d["sequence_id"].nunique(),
                "n_nodes": d["cluster_id"].nunique(),
                **_group_metadata(group_col, group_label),
                **dropped_metadata(dropped_rows, dropped_strata),
            }
            wald_tables.append(add_model_metadata(fit.wald, **meta))
            or_tables.append(add_model_metadata(fit.odds_ratios, **meta))
            fit_tables.append(
                add_fit_metadata(
                    model_fit_stats(
                        fit.result, model_name=model_name, formula=fit.formula
                    ),
                    **meta,
                )
            )
            print(f"Fitted {model_name}: {len(d):,} rows", flush=True)
        except Exception as exc:  # keep other clade/predictor fits moving
            _record_failure(
                failures,
                domain="composition",
                model_set=model_set,
                predictor_set="single",
                predictor=spec["name"],
                group_col=group_col,
                group_label=group_label,
                error=exc,
            )

    return {
        "wald": _concat_or_empty(wald_tables),
        "odds": _concat_or_empty(or_tables),
        "fit_stats": _concat_or_empty(fit_tables),
    }


def fit_joint_composition_model(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )
    predictors = [spec["column"] for spec in COMPOSITION_SPECS]
    try:
        d, dropped_rows, dropped_strata = _prepare_model_frame(
            source,
            predictors=predictors,
            adjusters=adjusters,
            required_base=[
                "candidate",
                "cluster_id",
                "sequence_id",
                cluster_se,
                window_strata,
            ],
            window_strata=window_strata,
            drop_window_nonvarying=drop_window,
            categorical_predictors=True,
        )

        terms = []
        references = {}
        for spec in COMPOSITION_SPECS:
            reference = resolve_reference(
                d,
                spec["column"],
                spec.get("reference"),
                spec.get("fallback_references"),
            )
            references[spec["name"]] = reference
            terms.append(categorical_term(spec["column"], reference))

        formula = "candidate ~ " + " + ".join(terms + adjusters)
        model_name = _model_name(f"composition__{model_set}__joint", group_label)
        result = fit_association_result(
            d,
            formula,
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
        )

        wald_tables = []
        for spec, term in zip(COMPOSITION_SPECS, terms):
            wald = robust_wald_for_prefix(
                result,
                term,
                model_name=model_name,
                term=spec["name"],
            )
            meta = {
                "domain": "composition",
                "model_set": model_set,
                "predictor_set": "joint",
                "predictor": spec["name"],
                "label": spec["label"],
                "reference": references[spec["name"]],
                "n_model_rows": len(d),
                "n_sequences": d["sequence_id"].nunique(),
                "n_nodes": d["cluster_id"].nunique(),
                **_group_metadata(group_col, group_label),
                **dropped_metadata(dropped_rows, dropped_strata),
            }
            wald_tables.append(add_model_metadata(wald, **meta))

        term_names = [
            name for term in terms for name in parameter_names_for_term(result, term)
        ]
        odds = tidy_odds_ratios(result, model_name=model_name)
        odds = odds.loc[odds["term"].isin(term_names)].copy()
        common_meta = {
            "domain": "composition",
            "model_set": model_set,
            "predictor_set": "joint",
            "predictor": "all_composition",
            "label": "All composition predictors",
            "reference": repr(references),
            "n_model_rows": len(d),
            "n_sequences": d["sequence_id"].nunique(),
            "n_nodes": d["cluster_id"].nunique(),
            **_group_metadata(group_col, group_label),
            **dropped_metadata(dropped_rows, dropped_strata),
        }
        fit_stats = model_fit_stats(result, model_name=model_name, formula=formula)
        print(f"Fitted {model_name}: {len(d):,} rows", flush=True)
        return {
            "wald": _concat_or_empty(wald_tables),
            "odds": add_model_metadata(odds, **common_meta),
            "fit_stats": add_fit_metadata(fit_stats, **common_meta),
        }
    except Exception as exc:
        _record_failure(
            failures,
            domain="composition",
            model_set=model_set,
            predictor_set="joint",
            predictor="all_composition",
            group_col=group_col,
            group_label=group_label,
            error=exc,
        )
        return {
            "wald": pd.DataFrame(),
            "odds": pd.DataFrame(),
            "fit_stats": pd.DataFrame(),
        }


def run_composition_model_set(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    single = fit_single_composition_models(
        source=source,
        model_set=model_set,
        adjusters=adjusters,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
        group_col=group_col,
        group_label=group_label,
        failures=failures,
    )
    joint = fit_joint_composition_model(
        source=source,
        model_set=model_set,
        adjusters=adjusters,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
        group_col=group_col,
        group_label=group_label,
        failures=failures,
    )
    return {
        "wald": _concat_or_empty([single["wald"], joint["wald"]]),
        "odds": _concat_or_empty([single["odds"], joint["odds"]]),
        "fit_stats": _concat_or_empty([single["fit_stats"], joint["fit_stats"]]),
    }


def fit_single_mixing_models(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    mixing_features: Sequence[str],
    mixing_reference: str,
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    wald_tables = []
    or_tables = []
    fit_tables = []
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )

    for feature in mixing_features:
        if feature not in source.columns:
            print(f"Skipping {feature}: not found", flush=True)
            continue
        try:
            d, dropped_rows, dropped_strata = _prepare_model_frame(
                source,
                predictors=[feature],
                adjusters=adjusters,
                required_base=["candidate", "cluster_id", cluster_se, window_strata],
                window_strata=window_strata,
                drop_window_nonvarying=drop_window,
            )
            base_name = f"mixing__{model_set}__single__{feature}"
            model_name = _model_name(base_name, group_label)
            fit = fit_exposure_association(
                d,
                outcome="candidate",
                exposure=feature,
                adjusters=adjusters,
                model_name=model_name,
                model_method=model_method,
                cluster_se=cluster_se,
                window_strata=window_strata,
                categorical=False,
            )
            meta = {
                "domain": "node_mixing",
                "model_set": model_set,
                "predictor_set": "single",
                "predictor": feature,
                "label": feature.replace("_", " "),
                "reference": mixing_reference,
                "n_model_rows": len(d),
                "n_nodes": d["cluster_id"].nunique(),
                **_group_metadata(group_col, group_label),
                **dropped_metadata(dropped_rows, dropped_strata),
            }
            wald_tables.append(add_model_metadata(fit.wald, **meta))
            or_tables.append(
                add_model_metadata(
                    fit.odds_ratios.loc[fit.odds_ratios["term"].eq(feature)].copy(),
                    **meta,
                )
            )
            fit_tables.append(
                add_fit_metadata(
                    model_fit_stats(
                        fit.result, model_name=model_name, formula=fit.formula
                    ),
                    **meta,
                )
            )
            print(f"Fitted {model_name}: {len(d):,} nodes", flush=True)
        except Exception as exc:
            _record_failure(
                failures,
                domain="node_mixing",
                model_set=model_set,
                predictor_set="single",
                predictor=feature,
                group_col=group_col,
                group_label=group_label,
                error=exc,
            )

    return {
        "wald": _concat_or_empty(wald_tables),
        "odds": _concat_or_empty(or_tables),
        "fit_stats": _concat_or_empty(fit_tables),
    }


def fit_joint_mixing_model(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    mixing_features: Sequence[str],
    mixing_reference: str,
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    features = [feature for feature in mixing_features if feature in source.columns]
    drop_window = (
        _has_categorical_adjuster(adjusters, window_strata)
        or model_method == "conditional_logit_by_window"
    )
    try:
        d, dropped_rows, dropped_strata = _prepare_model_frame(
            source,
            predictors=features,
            adjusters=adjusters,
            required_base=["candidate", "cluster_id", cluster_se, window_strata],
            window_strata=window_strata,
            drop_window_nonvarying=drop_window,
        )
        formula = "candidate ~ " + " + ".join(features + adjusters)
        model_name = _model_name(f"mixing__{model_set}__joint", group_label)
        result = fit_association_result(
            d,
            formula,
            model_method=model_method,
            cluster_se=cluster_se,
            window_strata=window_strata,
        )

        common_meta = {
            "domain": "node_mixing",
            "model_set": model_set,
            "predictor_set": "joint",
            "predictor": "all_mixing",
            "label": "All mixing predictors",
            "reference": mixing_reference,
            "n_model_rows": len(d),
            "n_nodes": d["cluster_id"].nunique(),
            **_group_metadata(group_col, group_label),
            **dropped_metadata(dropped_rows, dropped_strata),
        }
        wald = tidy_single_parameter_wald(result, features, model_name=model_name)
        odds = tidy_odds_ratios(result, model_name=model_name)
        odds = odds.loc[odds["term"].isin(features)].copy()
        fit_stats = model_fit_stats(result, model_name=model_name, formula=formula)
        print(f"Fitted {model_name}: {len(d):,} nodes", flush=True)
        return {
            "wald": add_model_metadata(wald, **common_meta),
            "odds": add_model_metadata(odds, **common_meta),
            "fit_stats": add_fit_metadata(fit_stats, **common_meta),
        }
    except Exception as exc:
        _record_failure(
            failures,
            domain="node_mixing",
            model_set=model_set,
            predictor_set="joint",
            predictor="all_mixing",
            group_col=group_col,
            group_label=group_label,
            error=exc,
        )
        return {
            "wald": pd.DataFrame(),
            "odds": pd.DataFrame(),
            "fit_stats": pd.DataFrame(),
        }


def run_mixing_model_set(
    *,
    source: pd.DataFrame,
    model_set: str,
    adjusters: list[str],
    mixing_features: Sequence[str],
    mixing_reference: str,
    model_method: str,
    cluster_se: str,
    window_strata: str,
    group_col: str | None,
    group_label: object | None,
    failures: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    single = fit_single_mixing_models(
        source=source,
        model_set=model_set,
        adjusters=adjusters,
        mixing_features=mixing_features,
        mixing_reference=mixing_reference,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
        group_col=group_col,
        group_label=group_label,
        failures=failures,
    )
    joint = fit_joint_mixing_model(
        source=source,
        model_set=model_set,
        adjusters=adjusters,
        mixing_features=mixing_features,
        mixing_reference=mixing_reference,
        model_method=model_method,
        cluster_se=cluster_se,
        window_strata=window_strata,
        group_col=group_col,
        group_label=group_label,
        failures=failures,
    )
    return {
        "wald": _concat_or_empty([single["wald"], joint["wald"]]),
        "odds": _concat_or_empty([single["odds"], joint["odds"]]),
        "fit_stats": _concat_or_empty([single["fit_stats"], joint["fit_stats"]]),
    }


def clean_export_table(table: pd.DataFrame) -> pd.DataFrame:
    out = table.copy()
    out.columns = [str(col).strip() for col in out.columns]
    for col in out.select_dtypes(include=["object", "string"]).columns:
        present = out[col].notna()
        out.loc[present, col] = out.loc[present, col].astype(str).str.strip()
    return out


def select_table_columns(
    table: pd.DataFrame,
    kind: str,
    *,
    sort_by: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Return a notebook-display view with available standard columns only."""
    if kind not in TABLE_DISPLAY_COLUMNS:
        raise KeyError(
            f"Unknown table kind {kind!r}. "
            f"Expected one of: {sorted(TABLE_DISPLAY_COLUMNS)}."
        )
    columns = [col for col in TABLE_DISPLAY_COLUMNS[kind] if col in table.columns]
    out = table.loc[:, columns].copy()
    if sort_by:
        present_sort = [col for col in sort_by if col in out.columns]
        if present_sort:
            out = out.sort_values(present_sort)
    return out


def _summarise_tables(
    *,
    composition_results: list[dict[str, pd.DataFrame]],
    mixing_results: list[dict[str, pd.DataFrame]],
) -> dict[str, pd.DataFrame]:
    summary: dict[str, pd.DataFrame] = {}

    if composition_results:
        composition_wald = _concat_or_empty(
            [result["wald"] for result in composition_results]
        )
        composition_or = _concat_or_empty(
            [result["odds"] for result in composition_results]
        )
        composition_fit_stats = _concat_or_empty(
            [result["fit_stats"] for result in composition_results]
        )
        if not composition_wald.empty:
            group_cols = ["domain", "model_set", "predictor_set"]
            if "analysis_group_col" in composition_wald.columns:
                group_cols.extend(
                    col
                    for col in composition_wald["analysis_group_col"]
                    .dropna()
                    .astype(str)
                    .unique()
                    if col in composition_wald.columns
                )
            composition_wald = bh_adjust_by(composition_wald, group_cols)
        summary.update(
            {
                "composition_wald.csv": composition_wald,
                "composition_odds_ratios.csv": composition_or,
                "composition_fit_stats.csv": composition_fit_stats,
            }
        )

    if mixing_results:
        mixing_wald = _concat_or_empty([result["wald"] for result in mixing_results])
        mixing_or = _concat_or_empty([result["odds"] for result in mixing_results])
        mixing_fit_stats = _concat_or_empty(
            [result["fit_stats"] for result in mixing_results]
        )
        if not mixing_wald.empty:
            group_cols = ["domain", "model_set", "predictor_set"]
            if "analysis_group_col" in mixing_wald.columns:
                group_cols.extend(
                    col
                    for col in mixing_wald["analysis_group_col"]
                    .dropna()
                    .astype(str)
                    .unique()
                    if col in mixing_wald.columns
                )
            mixing_wald = bh_adjust_by(mixing_wald, group_cols)
        summary.update(
            {
                "mixing_wald.csv": mixing_wald,
                "mixing_odds_ratios.csv": mixing_or,
                "mixing_fit_stats.csv": mixing_fit_stats,
            }
        )

    return summary


def run_association_pipeline(
    *,
    project_root: Path | str | None = None,
    output_dir: Path | str | None = None,
    result_dir: Path | str | None = None,
    result_subdir: str = "association_outputs",
    cluster_se: str = "cluster_id",
    window_strata: str = "window_idx",
    model_method: str = "firth_glm",
    variant_adjuster: str | None = "clade",
    window_adjustment: str = "fixed_effects",
    composition_model_sets: Mapping[str, Sequence[str]] | None = None,
    mixing_model_sets: Mapping[str, Sequence[str]] | None = None,
    mixing_features: Sequence[str] = DEFAULT_MIXING_FEATURES,
    mixing_reference: str = "per 1 null-model SD increase in entropy",
    group_by_clade: bool = False,
    clade_group_col: str = "clade_group",
    clade_group_values: Sequence[object] | None = None,
    window_stride: int = 2,
    run_composition: bool = True,
    run_mixing: bool = True,
) -> dict[str, Any]:
    """Run composition and/or mixing association models and export CSV tables."""
    root = _project_root(project_root)
    output_path = (
        Path(output_dir)
        if output_dir
        else root / "sse_detection" / "results" / "sse_outputs"
    )
    result_path = (
        Path(result_dir)
        if result_dir
        else root / "sse_detection" / "results" / result_subdir
    )
    result_path.mkdir(parents=True, exist_ok=True)

    if group_by_clade and variant_adjuster == "clade":
        variant_adjuster = None

    if composition_model_sets is None:
        composition_model_sets = default_model_sets(
            variant_adjuster=variant_adjuster,
            window_adjustment=window_adjustment,
        )
    if mixing_model_sets is None:
        mixing_model_sets = default_model_sets(
            variant_adjuster=variant_adjuster,
            window_adjustment=window_adjustment,
        )

    frames = load_association_frames(
        project_root=root,
        output_dir=output_path,
        cluster_se=cluster_se,
        variant_adjuster=variant_adjuster,
        group_by_clade=group_by_clade,
        clade_group_col=clade_group_col,
        window_stride=window_stride,
        run_composition=run_composition,
    )

    group_col = clade_group_col if group_by_clade else None
    if clade_group_values is None:
        clade_group_values = [*CLADES.values(), "Other"]

    failures: list[dict[str, object]] = []
    composition_results: list[dict[str, pd.DataFrame]] = []
    mixing_results: list[dict[str, pd.DataFrame]] = []

    if run_composition:
        for model_set, adjusters in composition_model_sets.items():
            for group_label, source in _iter_groups(
                frames.composition_base,
                group_col,
                clade_group_values,
            ):
                composition_results.append(
                    run_composition_model_set(
                        source=source,
                        model_set=model_set,
                        adjusters=list(adjusters),
                        model_method=model_method,
                        cluster_se=cluster_se,
                        window_strata=window_strata,
                        group_col=group_col,
                        group_label=group_label,
                        failures=failures,
                    )
                )

    if run_mixing:
        for model_set, adjusters in mixing_model_sets.items():
            for group_label, source in _iter_groups(
                frames.node_model_base,
                group_col,
                clade_group_values,
            ):
                mixing_results.append(
                    run_mixing_model_set(
                        source=source,
                        model_set=model_set,
                        adjusters=list(adjusters),
                        mixing_features=mixing_features,
                        mixing_reference=mixing_reference,
                        model_method=model_method,
                        cluster_se=cluster_se,
                        window_strata=window_strata,
                        group_col=group_col,
                        group_label=group_label,
                        failures=failures,
                    )
                )

    summary_tables = _summarise_tables(
        composition_results=composition_results,
        mixing_results=mixing_results,
    )

    for filename, table in summary_tables.items():
        path = result_path / filename
        clean_export_table(table).to_csv(path, index=False)
        print(f"saved {filename}: {len(table):,} rows", flush=True)

    failures_df = pd.DataFrame(failures)
    if not failures_df.empty:
        clean_export_table(failures_df).to_csv(
            result_path / "model_failures.csv", index=False
        )
        print(f"saved model_failures.csv: {len(failures_df):,} rows", flush=True)

    frames.cluster_diagnostics.to_csv(
        result_path / "cluster_diagnostics.csv", index=False
    )

    return {
        "result_dir": result_path,
        "frames": frames,
        "summary_tables": summary_tables,
        "failures": failures_df,
        "cluster_diagnostics": frames.cluster_diagnostics,
    }


def run_main_association_analysis(
    *,
    project_root: Path | str | None = None,
    result_dir: Path | str | None = None,
    model_method: str = "firth_glm",
    variant_adjuster: str | None = "clade",
    window_adjustment: str = "fixed_effects",
    run_composition: bool = True,
    run_mixing: bool = True,
    **kwargs,
) -> dict[str, Any]:
    """Run the primary overall socio-geodemographic association analysis.

    This is a thin preset around :func:`run_association_pipeline` for the main
    notebook: no clade stratification, primary and expanded model sets,
    sequence-level composition models, and node-level entropy z-score mixing
    models saved under ``sse_detection/results/association_outputs`` by default.
    """
    return run_association_pipeline(
        project_root=project_root,
        result_dir=result_dir,
        result_subdir="association_outputs",
        model_method=model_method,
        variant_adjuster=variant_adjuster,
        window_adjustment=window_adjustment,
        group_by_clade=False,
        run_composition=run_composition,
        run_mixing=run_mixing,
        **kwargs,
    )
