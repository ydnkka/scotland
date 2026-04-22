import json
import os
from typing import Literal

import pandas as pd
import bambi as bmb
import arviz as az


def fit_bambi_model(
    data: pd.DataFrame,
     run_id: int | str,
    dependent: str,
    *,
    family: Literal["negativebinomial", "binomial"] = "negativebinomial",
    fixed_effects: list[str] = None,
    interaction_effects: list[str] = None,
    random_effects: list[str] = None,
    offset: str = "log_seq_prop",
    formula: str = None,
    save_dir: str = 'bambi_outputs',
    draws: int = 1000,
    tune: int = 1000,
    target_accept: float = 0.95,
    chains: int = 4,
    random_seed: int = 42,
    force_refit: bool = False,
) -> tuple[bmb.Model, az.InferenceData]:
    """
    Fit a Bambi hierarchical model and cache results to disk.

    This function automates the construction of a Bambi model formula, fits
    the model using NUTS (via numpyro if specified), and saves the resulting
    trace and summary statistics for reproducibility.

    Args:
        data: The pandas DataFrame containing the variables.
        run_id: A unique identifier for the model run (used for file naming).
        dependent: The target variable (y-axis).
        family: The likelihood family. Supported: 'negativebinomial' or 'binomial'.
        fixed_effects: Main effect terms. Defaults to SIMD quintile and Epoch.
        interaction_effects: Terms specifying relationships between variables.
            Use ':' for a specific interaction (e.g., 'var1:var2').
            Use '*' for main effects plus interaction (e.g., 'var1*var2').
        random_effects: Hierarchical/Grouping terms. Default is (1 | window_id).
        offset: The name of the column to be used as an offset (log-scale).
        formula: Optional full formula string. If provided, overrides effect lists.
        save_dir: Directory where outputs (.nc and .csv) will be stored.
        draws: Number of posterior samples to draw.
        tune: Number of burn-in/tuning iterations.
        target_accept: Target acceptance rate for the NUTS sampler.
        chains: Number of independent MCMC chains.
        random_seed: Seed for reproducibility.
        force_refit: If True, ignores cached files and re-runs the model.

    Returns:
        A tuple containing the initialized (bmb.Model) and the (az.InferenceData) trace.

    Example:
        # To specify an interaction between two categorical variables:
        fit_bambi_model(
            data=df,
            run_id="interaction_test",
            interaction_effects=["C(epoch):C(simd_quintile_mode)"]
        )
    """

    family = family.lower()
    supported_families = {"negativebinomial", "binomial"}
    if family not in supported_families:
        raise ValueError(
            f"Unsupported family '{family}'. Use one of: {sorted(supported_families)}"
        )

    # --- Logic for Formula Construction ---
    if formula is None:
        if fixed_effects is None:
            fixed_effects = ["C(simd_quintile_mode, Treatment(3))"]

        if random_effects is None:
            random_effects = ["(1 | window_id)"]

        terms = list(fixed_effects)

        if interaction_effects:
            terms.extend(interaction_effects)

        terms.extend(random_effects)

        if offset:
            terms.append(f"offset({offset})")

        formula = f"{dependent} ~ {' + '.join(terms)}"

    # --- Caching and Execution Logic ---
    os.makedirs(save_dir, exist_ok=True)
    trace_path = os.path.join(save_dir, f"{run_id}_{family}_bambi_trace.nc")
    formula_path = os.path.join(save_dir, f"{run_id}_{family}_bambi_formula.json")
    summary_path = os.path.join(save_dir, f"{run_id}_{family}_bambi_summary.csv")

    if not force_refit and os.path.exists(trace_path):
        print(f"Loading cached results for run_id: {run_id}")
        model = bmb.Model(formula, data, family=family)
        trace = az.from_netcdf(trace_path)
        return model, trace

    print(f"Fitting model: {run_id}\nFormula: {formula}")
    model = bmb.Model(formula, data, family=family)

    trace = model.fit(
        draws=draws, tune=tune, target_accept=target_accept,
        chains=chains, random_seed=random_seed,
        inference_method="nuts_numpyro"
    )

    summary = az.summary(trace)
    summary.insert(0, 'terms', summary.index)
    summary.insert(1, 'run_id', run_id)
    summary.insert(2, 'family', family)
    summary.reset_index(drop=True, inplace=True)

    az.to_netcdf(trace, trace_path)
    summary.to_csv(summary_path, index=False)
    with open(formula_path, 'w') as f:
        json.dump({'formula': formula}, f)

    return model, trace