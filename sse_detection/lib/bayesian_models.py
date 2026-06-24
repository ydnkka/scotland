"""Bambi fitting, posterior summaries, diagnostics, and result writing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import xarray as xr

from .concurrent_io import (
    atomic_write_csv,
    atomic_write_netcdf,
    exclusive_file_lock,
)
from .regression_prep import GROUP_VARS, PreparedModelFrame, PreparedRegressionRun


RANDOM_SEED = 123
_RANDOM_EFFECT_RE = re.compile(r"\([^|()]+\|([^()]+)\)")


@dataclass(frozen=True)
class BayesianFitConfig:
    """Sampling and prior settings shared by Bambi regression models."""

    draws: int = 2_000
    tune: int = 2_000
    chains: int = 4
    cores: int = 4
    target_accept: float = 0.99
    inference_method: str = "pymc"
    random_seed: int = RANDOM_SEED
    fixed_prior_sigma: float | None = None
    intercept_prior_sigma: float | None = None
    random_effect_sigma: float | None = None
    residual_sigma: float | None = None
    log_likelihood: bool = True
    noncentered: bool = True
    progressbar: bool = False
    quiet: bool = True


@dataclass
class BayesianModelResult:
    """Container for one fitted Bambi model and derived summaries."""

    family: str
    formula: str
    model: Any
    idata: Any
    summary: pd.DataFrame
    diagnostics: pd.DataFrame
    var_names: list[str]
    n_rows: int
    outcome: str
    outcome_mean: float
    outcome_sd: float | None = None

    @property
    def candidate_rate(self) -> float | None:
        return self.outcome_mean if self.family == "logistic" else None

    def as_dict(self) -> dict[str, object]:
        """Return a notebook-compatible dictionary."""
        out: dict[str, object] = {
            "family": self.family,
            "formula": self.formula,
            "model": self.model,
            "idata": self.idata,
            "summary": self.summary,
            "diagnostics": self.diagnostics,
            "var_names": self.var_names,
            "n_rows": self.n_rows,
            "outcome": self.outcome,
            "outcome_mean": self.outcome_mean,
        }
        if self.family == "logistic":
            out["candidate_rate"] = self.outcome_mean
        if self.outcome_sd is not None:
            out["outcome_sd"] = self.outcome_sd
        return out


def fit_prepared_model(
    frame: PreparedModelFrame,
    *,
    config: BayesianFitConfig | None = None,
    var_names: Sequence[str] | None = None,
    categorical: Sequence[str] = GROUP_VARS,
    display_tables: bool = False,
    print_diagnostics: bool = True,
    include_auxiliary: bool = True,
) -> BayesianModelResult:
    """Fit and summarise one ``PreparedModelFrame``."""
    return fit_and_summarise_model(
        data=frame.fit_df,
        formula=frame.formula,
        family=frame.family,
        config=config,
        var_names=var_names,
        categorical=categorical,
        display_tables=display_tables,
        print_diagnostics=print_diagnostics,
        include_auxiliary=include_auxiliary,
    )


def fit_prepared_run(
    prepared: PreparedRegressionRun,
    *,
    domains: Sequence[str] | None = None,
    outcomes: Sequence[str] | None = None,
    model_sets: Sequence[str] | None = None,
    config: BayesianFitConfig | None = None,
    var_names: Sequence[str] | None = None,
    categorical: Sequence[str] = GROUP_VARS,
    display_tables: bool = False,
    print_diagnostics: bool = True,
    include_auxiliary: bool = True,
) -> dict[str, BayesianModelResult]:
    """Fit selected frames from a prepared regression run."""
    domain_set = set(domains) if domains is not None else None
    outcome_set = set(outcomes) if outcomes is not None else None
    model_set_filter = set(model_sets) if model_sets is not None else None

    results: dict[str, BayesianModelResult] = {}
    for frame in prepared.frames.values():
        if domain_set is not None and frame.domain not in domain_set:
            continue
        if outcome_set is not None and frame.outcome not in outcome_set:
            continue
        if model_set_filter is not None and frame.model_set not in model_set_filter:
            continue
        results[frame.result_key] = fit_prepared_model(
            frame,
            config=config,
            var_names=var_names,
            categorical=categorical,
            display_tables=display_tables,
            print_diagnostics=print_diagnostics,
            include_auxiliary=include_auxiliary,
        )
    return results


def fit_and_summarise_model(
    data: pd.DataFrame,
    formula: str,
    *,
    family: str,
    config: BayesianFitConfig | None = None,
    var_names: Sequence[str] | None = None,
    categorical: Sequence[str] = GROUP_VARS,
    display_tables: bool = False,
    print_diagnostics: bool = True,
    include_auxiliary: bool = True,
) -> BayesianModelResult:
    """Fit a Bambi model and return posterior summaries plus diagnostics."""
    config = BayesianFitConfig() if config is None else config
    model, idata = fit_bayesian_model(
        data=data,
        formula=formula,
        family=family,
        categorical=categorical,
        config=config,
    )
    summary, diagnostics = summarise_bambi_idata(
        idata,
        family=family,
        var_names=var_names,
        formula=formula,
        print_diagnostics=print_diagnostics,
        display_tables=display_tables,
        include_auxiliary=include_auxiliary,
    )
    response = response_from_formula(formula)
    selected_vars = available_posterior_vars(
        idata,
        var_names,
        formula=formula,
        include_auxiliary=include_auxiliary,
    )
    return BayesianModelResult(
        family=family,
        formula=formula,
        model=model,
        idata=idata,
        summary=summary,
        diagnostics=diagnostics,
        var_names=selected_vars,
        n_rows=len(data),
        outcome=response,
        outcome_mean=float(data[response].mean()),
        outcome_sd=None if family == "logistic" else float(data[response].std()),
    )


def fit_bayesian_model(
    data: pd.DataFrame,
    formula: str,
    *,
    family: str,
    categorical: Sequence[str] | None = None,
    config: BayesianFitConfig | None = None,
) -> tuple[Any, Any]:
    """Fit a hierarchical Bambi model for logistic or linear regression."""
    bmb, _ = _bayesian_deps()
    config = BayesianFitConfig() if config is None else config
    response = response_from_formula(formula)
    if response not in data.columns:
        raise ValueError(f"Response variable '{response}' not found in data.")

    family = _normalise_family(family)
    if family == "logistic":
        bambi_family = "bernoulli"
        outcome_mean = float(np.clip(data[response].mean(), 1e-6, 1 - 1e-6))
        priors = {
            "Intercept": bmb.Prior(
                "Normal",
                mu=_logit(outcome_mean),
                sigma=_coalesce(config.intercept_prior_sigma, 1.5),
            ),
            "common": bmb.Prior(
                "Normal",
                mu=0,
                sigma=_coalesce(config.fixed_prior_sigma, 1.0),
            ),
            "group_specific": bmb.Prior(
                "Normal",
                mu=0,
                sigma=bmb.Prior(
                    "HalfNormal",
                    sigma=_coalesce(config.random_effect_sigma, 1.0),
                ),
            ),
        }
    else:
        bambi_family = "gaussian"
        priors = {
            "Intercept": bmb.Prior(
                "Normal",
                mu=float(data[response].mean()),
                sigma=_coalesce(config.intercept_prior_sigma, 1.0),
            ),
            "common": bmb.Prior(
                "Normal",
                mu=0,
                sigma=_coalesce(config.fixed_prior_sigma, 0.5),
            ),
            "group_specific": bmb.Prior(
                "Normal",
                mu=0,
                sigma=bmb.Prior(
                    "HalfNormal",
                    sigma=_coalesce(config.random_effect_sigma, 0.5),
                ),
            ),
            "sigma": bmb.Prior(
                "HalfNormal",
                sigma=_coalesce(config.residual_sigma, 0.5),
            ),
        }

    group_cols = random_effect_groups(formula)
    categorical_cols = _unique_preserve_order([*(categorical or ()), *group_cols])
    inference_method = config.inference_method
    if inference_method == "nutpie" and has_categorical_formula_terms(formula):
        print(
            "Using PyMC sampler for categorical formula terms; "
            "nutpie can mis-convert categorical coefficient dimensions."
        )
        inference_method = "pymc"

    model = bmb.Model(
        formula=formula,
        data=data,
        family=bambi_family,
        priors=priors,
        categorical=categorical_cols or None,
        noncentered=config.noncentered,
    )
    fit_kwargs = {
        "draws": config.draws,
        "tune": config.tune,
        "chains": config.chains,
        "cores": config.cores,
        "inference_method": inference_method,
        "target_accept": config.target_accept,
        "random_seed": config.random_seed,
        "progressbar": config.progressbar,
        "quiet": config.quiet,
    }
    try:
        idata = model.fit(**fit_kwargs)
    except (EOFError, OSError) as exc:
        if config.cores <= 1 or (
            isinstance(exc, OSError) and "end of file" not in str(exc).lower()
        ):
            raise
        print(
            f"PyMC multiprocessing failed with {type(exc).__name__}; "
            "retrying with cores=1."
        )
        fit_kwargs["cores"] = 1
        idata = model.fit(**fit_kwargs)

    if config.log_likelihood:
        try:
            model.compute_log_likelihood(idata)
        except Exception as exc:
            print(f"Error computing log likelihood: {exc}")
    return model, idata


def summarise_bambi_idata(
    idata: Any,
    *,
    family: str,
    var_names: Sequence[str] | None = None,
    formula: str | None = None,
    ci_prob: float = 0.95,
    print_diagnostics: bool = True,
    rhat_threshold: float = 1.01,
    ess_threshold: int = 400,
    display_tables: bool = False,
    float_digits: int = 4,
    include_auxiliary: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarise posterior variables and print useful sampler diagnostics."""
    _, az = _bayesian_deps()
    family = _normalise_family(family)
    selected_vars = available_posterior_vars(
        idata,
        var_names,
        formula=formula,
        include_auxiliary=include_auxiliary,
    )
    summary = az.summary(
        idata,
        var_names=selected_vars,
        ci_prob=ci_prob,
        ci_kind="hdi",
    )
    posterior_ds = posterior_dataset(idata)
    beta_ds = posterior_ds[selected_vars]

    if family == "logistic":
        scale_summary = pd.DataFrame(
            az.summary(np.exp(beta_ds), ci_prob=ci_prob, ci_kind="hdi")
        ).add_prefix("OR_")
        probability_summary = posterior_probability_summary(
            beta_ds,
            positive_label="P(OR > 1 | data)",
            negative_label="P(OR < 1 | data)",
        )
        combined_summary = pd.concat(
            [pd.DataFrame(summary), scale_summary, probability_summary],
            axis=1,
        )
    else:
        direction_summary = posterior_probability_summary(
            beta_ds,
            positive_label="P(beta > 0 | data)",
            negative_label="P(beta < 0 | data)",
        )
        combined_summary = pd.concat(
            [pd.DataFrame(summary), direction_summary],
            axis=1,
        )

    diagnostics = model_diagnostics(
        idata,
        ci_prob=ci_prob,
        rhat_threshold=rhat_threshold,
        ess_threshold=ess_threshold,
    )
    if print_diagnostics:
        print_diagnostic_report(
            diagnostics,
            combined_summary,
            display_tables=display_tables,
            float_digits=float_digits,
        )
    return combined_summary, diagnostics


def model_diagnostics(
    idata: Any,
    *,
    ci_prob: float = 0.95,
    rhat_threshold: float = 1.01,
    ess_threshold: int = 400,
) -> pd.DataFrame:
    """Return divergences, BFMI, R-hat, ESS, and tree-depth checks."""
    _, az = _bayesian_deps()
    diagnostic_rows = []

    if hasattr(idata, "sample_stats") and "diverging" in idata.sample_stats:
        n_div = float(idata.sample_stats["diverging"].sum().item())
        n_total = float(idata.sample_stats["diverging"].size)
        div_rate = n_div / n_total
        diagnostic_rows.append(
            {
                "Diagnostic": "Divergences",
                "Value": f"{n_div:.0f} / {n_total:.0f} ({div_rate:.2%})",
                "Status": "OK" if n_div == 0 else "WARNING",
                "Interpretation": (
                    "No divergent transitions."
                    if n_div == 0
                    else "Investigate divergent transitions."
                ),
            }
        )

    try:
        energy = idata.sample_stats["energy"]
        bfmi_result = az.bfmi(energy)
        bfmi = (
            np.asarray(
                bfmi_result.values if hasattr(bfmi_result, "values") else bfmi_result
            )
            .astype(float)
            .ravel()
        )
        bfmi = bfmi[np.isfinite(bfmi)]
        if bfmi.size == 0:
            raise ValueError("BFMI returned no finite values.")
        min_bfmi = float(np.nanmin(bfmi))
        bfmi_by_chain = ", ".join(f"{value:.3f}" for value in bfmi)
        diagnostic_rows.append(
            {
                "Diagnostic": "BFMI",
                "Value": f"min={min_bfmi:.3f}; chains=[{bfmi_by_chain}]",
                "Status": "OK" if min_bfmi >= 0.3 else "WARNING",
                "Interpretation": (
                    "Energy exploration looks acceptable."
                    if min_bfmi >= 0.3
                    else "One or more chains have BFMI < 0.3."
                ),
            }
        )
    except Exception as exc:
        diagnostic_rows.append(
            {
                "Diagnostic": "BFMI",
                "Value": "Could not compute",
                "Status": "NA",
                "Interpretation": str(exc),
            }
        )

    diagnostic_summary = pd.DataFrame(az.summary(idata, ci_prob=ci_prob, ci_kind="hdi"))
    if "r_hat" in diagnostic_summary.columns:
        max_rhat = float(diagnostic_summary["r_hat"].max(skipna=True))
        diagnostic_rows.append(
            {
                "Diagnostic": "Max R-hat",
                "Value": f"{max_rhat:.4f}",
                "Status": "OK" if max_rhat <= rhat_threshold else "WARNING",
                "Interpretation": (
                    f"All posterior variables are at or below {rhat_threshold}."
                    if max_rhat <= rhat_threshold
                    else f"Some posterior variables exceed {rhat_threshold}."
                ),
            }
        )

    if "ess_bulk" in diagnostic_summary.columns:
        min_bulk_ess = float(diagnostic_summary["ess_bulk"].min(skipna=True))
        diagnostic_rows.append(
            {
                "Diagnostic": "Min bulk ESS",
                "Value": f"{min_bulk_ess:.1f}",
                "Status": "OK" if min_bulk_ess >= ess_threshold else "WARNING",
                "Interpretation": (
                    f"All bulk ESS values are at least {ess_threshold}."
                    if min_bulk_ess >= ess_threshold
                    else f"Some bulk ESS values are below {ess_threshold}."
                ),
            }
        )

    if "ess_tail" in diagnostic_summary.columns:
        min_tail_ess = float(diagnostic_summary["ess_tail"].min(skipna=True))
        diagnostic_rows.append(
            {
                "Diagnostic": "Min tail ESS",
                "Value": f"{min_tail_ess:.1f}",
                "Status": "OK" if min_tail_ess >= ess_threshold else "WARNING",
                "Interpretation": (
                    f"All tail ESS values are at least {ess_threshold}."
                    if min_tail_ess >= ess_threshold
                    else f"Some tail ESS values are below {ess_threshold}."
                ),
            }
        )

    if hasattr(idata, "sample_stats") and "tree_depth" in idata.sample_stats:
        max_tree_depth = float(idata.sample_stats["tree_depth"].max().item())
        diagnostic_rows.append(
            {
                "Diagnostic": "Max tree depth",
                "Value": str(max_tree_depth),
                "Status": "INFO",
                "Interpretation": "Maximum observed tree depth.",
            }
        )

    return pd.DataFrame(diagnostic_rows)


def posterior_probability_summary(
    beta_ds: xr.Dataset,
    *,
    positive_label: str,
    negative_label: str,
) -> pd.DataFrame:
    """Return posterior sign probabilities for each scalar or indexed variable."""
    rows = []
    for var_name, beta in beta_ds.data_vars.items():
        prob_positive = (beta > 0).mean(dim=("chain", "draw"))
        prob_negative = (beta < 0).mean(dim=("chain", "draw"))

        if prob_positive.ndim == 0:
            rows.append(
                {
                    "Variable": var_name,
                    positive_label: float(prob_positive.item()),
                    negative_label: float(prob_negative.item()),
                }
            )
            continue

        prob_positive_df = prob_positive.to_dataframe(name=positive_label)
        prob_negative_df = prob_negative.to_dataframe(name=negative_label)
        prob_df = pd.concat([prob_positive_df, prob_negative_df], axis=1).reset_index()
        coord_cols = [
            col
            for col in prob_df.columns
            if col not in {positive_label, negative_label}
        ]
        for _, row in prob_df.iterrows():
            coords = ", ".join(str(row[col]) for col in coord_cols)
            rows.append(
                {
                    "Variable": f"{var_name}[{coords}]",
                    positive_label: float(row[positive_label]),
                    negative_label: float(row[negative_label]),
                }
            )

    return pd.DataFrame(rows).set_index("Variable")


def save_model_result(
    result: BayesianModelResult | MappingResult,
    *,
    model_dir: Path | str,
    domain: str | None = None,
    model_set: str | None = None,
    use_sample: bool | None = None,
    save_idata: bool = False,
) -> dict[str, object]:
    """Write one model result and return a manifest row."""
    result_dict = _result_as_dict(result)
    model_dir = Path(model_dir)
    model_dir.mkdir(parents=True, exist_ok=True)

    for key in ("summary", "diagnostics"):
        value = result_dict.get(key)
        if isinstance(value, pd.DataFrame):
            write_summary_table(value, model_dir / f"{key}.csv")
        else:
            print(f"Warning: '{key}' not found or not a DataFrame.")

    if save_idata and result_dict.get("idata") is not None:
        atomic_write_netcdf(result_dict["idata"], model_dir / "idata.nc")

    metadata = pd.DataFrame(
        [
            {
                "family": result_dict.get("family"),
                "domain": domain,
                "outcome": result_dict.get("outcome"),
                "model_set": model_set,
                "formula": result_dict.get("formula"),
                "n_rows": result_dict.get("n_rows"),
                "candidate_rate": result_dict.get("candidate_rate"),
                "outcome_mean": result_dict.get("outcome_mean"),
                "outcome_sd": result_dict.get("outcome_sd"),
                "use_sample": use_sample,
            }
        ]
    )
    atomic_write_csv(metadata, model_dir / "metadata.csv", index=False)

    return {
        "family": result_dict.get("family"),
        "domain": domain,
        "outcome": result_dict.get("outcome"),
        "model_set": model_set,
        "model_dir": str(model_dir),
        "n_rows": result_dict.get("n_rows"),
        "candidate_rate": result_dict.get("candidate_rate"),
        "outcome_mean": result_dict.get("outcome_mean"),
        "outcome_sd": result_dict.get("outcome_sd"),
        "use_sample": use_sample,
    }


def save_prepared_model_result(
    result: BayesianModelResult | MappingResult,
    frame: PreparedModelFrame,
    *,
    save_idata: bool = False,
) -> dict[str, object]:
    """Write a result using the output directory from a prepared frame."""
    return save_model_result(
        result,
        model_dir=frame.output_dir,
        domain=frame.domain,
        model_set=frame.model_set,
        use_sample=len(frame.fit_df) != len(frame.full_df),
        save_idata=save_idata,
    )


def save_prepared_run_results(
    results: dict[str, BayesianModelResult | MappingResult],
    prepared: PreparedRegressionRun,
    *,
    save_idata: bool = False,
) -> pd.DataFrame:
    """Save all fitted results from a prepared run and return a manifest."""
    manifest_rows = []
    for frame in prepared.frames.values():
        result = results.get(frame.result_key)
        if result is None:
            continue
        manifest_rows.append(
            save_prepared_model_result(result, frame, save_idata=save_idata)
        )
    manifest = pd.DataFrame(manifest_rows)
    prepared.result_dir.mkdir(parents=True, exist_ok=True)
    with exclusive_file_lock(prepared.result_dir / ".saved_model_manifest.lock"):
        atomic_write_csv(
            manifest,
            prepared.result_dir / "saved_model_manifest.csv",
            index=False,
        )
        for domain, table in manifest.groupby("domain", sort=False):
            atomic_write_csv(
                table,
                prepared.result_dir / f"{domain}_saved_model_manifest.csv",
                index=False,
            )
    return manifest


def write_summary_table(summary: pd.DataFrame, path: Path | str) -> None:
    """Write a summary table while preserving parameter names from the index."""
    out = summary.copy()
    if "parameter" not in out.columns:
        out.insert(0, "parameter", out.index.astype(str))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    atomic_write_csv(out, path, index=False)


def print_diagnostic_report(
    diagnostics: pd.DataFrame,
    summary: pd.DataFrame | None = None,
    *,
    display_tables: bool = False,
    float_digits: int = 4,
) -> None:
    """Print diagnostics and, optionally, the posterior summary."""
    _print_section("Bayesian model diagnostics")
    _show_table(diagnostics, display_tables=display_tables, float_digits=float_digits)
    if summary is not None:
        _print_section("Posterior summary")
        _show_table(summary, display_tables=display_tables, float_digits=float_digits)


def posterior_dataset(idata: Any) -> xr.Dataset:
    """Return posterior samples as an xarray Dataset across ArviZ versions."""
    posterior = idata.posterior
    if isinstance(posterior, xr.Dataset):
        return posterior
    if hasattr(posterior, "to_dataset"):
        return posterior.to_dataset()
    raise TypeError("idata.posterior cannot be converted to an xarray Dataset.")


def available_posterior_vars(
    idata: Any,
    var_names: Sequence[str] | None,
    *,
    formula: str | None = None,
    include_auxiliary: bool = True,
) -> list[str]:
    """Resolve requested posterior variables against variables in ``idata``."""
    available = [str(var) for var in posterior_dataset(idata).data_vars]
    formula_vars = (
        posterior_vars_from_formula(formula) if formula is not None else available
    )
    requested = list(formula_vars if var_names is None else var_names)
    if include_auxiliary and formula is not None:
        requested = _unique_preserve_order(
            [*requested, *random_effect_sigma_vars_from_formula(formula), "sigma"]
        )

    available_set = set(available)
    selected = [var for var in requested if var in available_set]
    missing = [var for var in requested if var not in available_set]
    if missing:
        print("Skipping unavailable posterior variables:")
        for var in missing:
            print(f"  - {var}")
    if not selected:
        raise KeyError("None of the requested posterior variables were found.")
    return selected


def posterior_vars_from_formula(formula: str) -> list[str]:
    """Return Bambi posterior variable names implied by formula terms."""
    if "~" not in formula:
        raise ValueError("Formula must contain '~'.")
    rhs_terms = split_formula_terms(formula.split("~", 1)[1])
    include_intercept = True
    posterior_vars: list[str] = []
    for term in rhs_terms:
        if term in {"0", "-1"}:
            include_intercept = False
            continue
        if term == "1":
            include_intercept = True
            continue
        random_effect = re.fullmatch(r"\(([^|()]+)\|([^()]+)\)", term)
        if random_effect:
            posterior_vars.append(
                f"{random_effect.group(1).strip()}|{random_effect.group(2).strip()}"
            )
        else:
            posterior_vars.append(term)
    if include_intercept:
        posterior_vars.insert(0, "Intercept")
    return _unique_preserve_order(posterior_vars)


def random_effect_sigma_vars_from_formula(formula: str) -> list[str]:
    """Return Bambi random-effect SD variable names implied by formula terms."""
    if "~" not in formula:
        raise ValueError("Formula must contain '~'.")
    sigma_vars = []
    for term in split_formula_terms(formula.split("~", 1)[1]):
        random_effect = re.fullmatch(r"\(([^|()]+)\|([^()]+)\)", term)
        if random_effect:
            effect = random_effect.group(1).strip()
            group = random_effect.group(2).strip()
            sigma_vars.append(f"{effect}|{group}_sigma")
    return _unique_preserve_order(sigma_vars)


def split_formula_terms(rhs: str) -> list[str]:
    """Split a formula RHS on top-level plus signs."""
    terms: list[str] = []
    start = 0
    depth = 0
    quote: str | None = None
    for idx, char in enumerate(rhs):
        if quote is not None:
            if char == quote:
                quote = None
            continue
        if char in {"'", '"'}:
            quote = char
        elif char == "(":
            depth += 1
        elif char == ")":
            depth = max(0, depth - 1)
        elif char == "+" and depth == 0:
            term = rhs[start:idx].strip()
            if term:
                terms.append(term)
            start = idx + 1
    term = rhs[start:].strip()
    if term:
        terms.append(term)
    return terms


def response_from_formula(formula: str) -> str:
    """Return the response variable from a model formula."""
    if "~" not in formula:
        raise ValueError("Formula must contain '~'.")
    return formula.split("~", 1)[0].strip()


def random_effect_groups(formula: str) -> list[str]:
    """Extract grouping variables from terms such as ``(1|clade)``."""
    return [match.strip() for match in _RANDOM_EFFECT_RE.findall(formula)]


def has_categorical_formula_terms(formula: str) -> bool:
    """Return True when the formula uses explicit formulae ``C(...)`` terms."""
    return "C(" in formula


def _bayesian_deps() -> tuple[Any, Any]:
    try:
        import arviz as az
        import bambi as bmb
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Bayesian model fitting requires optional dependencies 'bambi' and "
            "'arviz'. Install them in the notebook environment before fitting."
        ) from exc
    return bmb, az


def _normalise_family(family: str) -> str:
    if family in {"logistic", "bernoulli", "binary"}:
        return "logistic"
    if family in {"linear", "gaussian", "normal"}:
        return "linear"
    raise ValueError("family must be 'logistic'/'bernoulli' or 'linear'/'gaussian'.")


def _logit(p: float) -> float:
    return float(np.log(p / (1 - p)))


def _coalesce(value: float | None, default: float) -> float:
    return default if value is None else value


def _unique_preserve_order(items: Iterable[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for item in items:
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


def _show_table(
    df: pd.DataFrame,
    *,
    display_tables: bool,
    float_digits: int,
) -> None:
    if display_tables:
        try:
            from IPython.display import display

            display(df.style.format(precision=float_digits))
            return
        except Exception:
            pass
    print(_format_df_for_print(df, float_digits=float_digits))


def _format_df_for_print(
    df: pd.DataFrame,
    *,
    float_digits: int = 4,
    width: int = 160,
    max_colwidth: int = 100,
) -> str:
    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.width",
        width,
        "display.max_colwidth",
        max_colwidth,
        "display.float_format",
        lambda x: f"{x:,.{float_digits}f}",
    ):
        return df.to_string()


def _print_section(title: str, char: str = "=") -> None:
    print(f"\n{title}")
    print(char * len(title))


MappingResult = dict[str, object]


def _result_as_dict(result: BayesianModelResult | MappingResult) -> dict[str, object]:
    if isinstance(result, BayesianModelResult):
        return result.as_dict()
    return dict(result)


__all__ = [
    "BayesianFitConfig",
    "BayesianModelResult",
    "available_posterior_vars",
    "fit_and_summarise_model",
    "fit_bayesian_model",
    "fit_prepared_model",
    "fit_prepared_run",
    "has_categorical_formula_terms",
    "model_diagnostics",
    "posterior_dataset",
    "posterior_probability_summary",
    "posterior_vars_from_formula",
    "print_diagnostic_report",
    "random_effect_groups",
    "random_effect_sigma_vars_from_formula",
    "response_from_formula",
    "save_model_result",
    "save_prepared_model_result",
    "save_prepared_run_results",
    "split_formula_terms",
    "summarise_bambi_idata",
    "write_summary_table",
]
