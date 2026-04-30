import json
from pathlib import Path
from typing import Literal

import arviz as az
import bambi as bmb
import pandas as pd
import polars as pl


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def fit_bambi_model(
    data: pd.DataFrame | pl.DataFrame,
    run_id: int | str,
    dependent: str,
    fixed_effects: list[str],
    family: Literal["negativebinomial", "binomial"],
    *,
    trials: str = None,
    interaction_effects: list[str] = None,
    random_effects: list[str] = None,
    offset: str  = None,
    formula: str = None,
    save_dir: str | Path = "bambi_outputs",
    draws: int = 1000,
    tune: int = 1000,
    target_accept: float = 0.95,
    chains: int = 4,
    random_seed: int = 42,
    force_refit: bool = False,
) -> tuple[bmb.Model, az.InferenceData]:
    """Fit a Bambi hierarchical model and cache results to disk.

    Parameters
    ----------
    data:
        Source data. Polars DataFrames are converted to pandas automatically.
    run_id:
        Unique identifier for the model run, used as a file-name prefix.
    dependent:
        Response variable name.
    family:
        Likelihood family — ``"negativebinomial"`` or ``"binomial"``.
    fixed_effects:
        Main-effect terms.
    interaction_effects:
        Interaction terms. Use ``":"`` for a pure interaction
        (e.g. ``"var1:var2"``), ``"*"`` to include main effects too.
    random_effects:
        Grouping / hierarchical terms.
    trials:
        Column of trial counts for ``binomial`` proportion models. When
        supplied the LHS becomes ``proportion(dependent, trials)``.
    offset:
        Column (on the log scale) to include as a model offset. Pass
        ``None`` to omit entirely.
    formula:
        Full formula string. When provided it overrides all effect lists,
        *trials*, and *offset* arguments.
    save_dir:
        Directory for cached ``.nc`` trace, ``.csv`` summary, and
        ``.json`` formula files.
    draws, tune, target_accept, chains, random_seed:
        PyMC / numpyro sampler settings.
    force_refit:
        When ``True``, ignore any cached trace and re-fit from scratch.

    Returns
    -------
    (model, trace)
        The initialised ``bmb.Model`` and the ``az.InferenceData`` trace.
    """
    if isinstance(data, pl.DataFrame):
        data = data.to_pandas()

    family = family.lower()
    if family not in {"negativebinomial", "binomial"}:
        raise ValueError(
            f"Unsupported family '{family}'. Choose 'negativebinomial' or 'binomial'."
        )
    if dependent not in data.columns:
        raise ValueError(f"Dependent column '{dependent}' not found in data.")
    if trials is not None and trials not in data.columns:
        raise ValueError(f"Trials column '{trials}' not found in data.")
    if offset is not None and offset not in data.columns:
        raise ValueError(f"Offset column '{offset}' not found in data.")

    # --- Formula construction ---
    if formula is None:
        effects = fixed_effects or []

        if interaction_effects:
            effects.extend(interaction_effects)

        if random_effects is not None:
            effects.extend(random_effects)

        if offset is not None:
            effects.append(f"offset({offset})")

        lhs = (
            f"proportion({dependent}, {trials})"
            if family == "binomial" and trials is not None
            else dependent
        )
        formula = f"{lhs} ~ {' + '.join(effects)}"

    # --- Caching ---
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    stem        = f"{run_id}_{family}"
    trace_path   = save_dir / f"{stem}_trace.nc"
    summary_path = save_dir / f"{stem}_summary.csv"
    formula_path = save_dir / f"{stem}_formula.json"

    if not force_refit and trace_path.exists():
        print(f"[{run_id}] Loading cached trace from {trace_path}")
        model = bmb.Model(formula, data, family=family)
        trace = az.from_netcdf(trace_path)
        return model, trace

    print(f"[{run_id}] Fitting  : {formula}")
    model = bmb.Model(formula, data, family=family, dropna=True)
    trace = model.fit(
        draws=draws,
        tune=tune,
        target_accept=target_accept,
        chains=chains,
        random_seed=random_seed,
        inference_method="numpyro",  # comment out for Ubuntu
    )

    summary = az.summary(trace)
    summary.insert(0, "terms",  summary.index)
    summary.insert(1, "run_id", run_id)
    summary.insert(2, "family", family)
    summary.reset_index(drop=True, inplace=True)

    az.to_netcdf(trace, str(trace_path))
    summary.to_csv(summary_path, index=False)
    formula_path.write_text(json.dumps({"formula": formula}, indent=2))

    return model, trace


# ---------------------------------------------------------------------------
# Domain wrappers
# ---------------------------------------------------------------------------

def fit_cluster_model(
    data: pd.DataFrame | pl.DataFrame,
    run_id: int | str,
    dependent: str,
    *,
    fixed_effects: list[str] = None,
    interaction_effects: list[str] = None,
    random_effects: list[str] = None,
    offset: str = "log_seq_prop",
    **kwargs,
) -> tuple[bmb.Model, az.InferenceData]:
    """Fit a negative-binomial model on cluster-level data.
    """
    fixed_effects = fixed_effects or []
    random_effects = random_effects or []
    return fit_bambi_model(
        data,
        run_id,
        dependent,
        family="negativebinomial",
        fixed_effects=fixed_effects,
        interaction_effects=interaction_effects,
        random_effects=random_effects,
        offset=offset,
        **kwargs,
    )


def fit_individual_model(
    data: pd.DataFrame | pl.DataFrame,
    run_id: int | str,
    success: str,
    trials: str,
    *,
    fixed_effects: list[str] | None = None,
    interaction_effects: list[str] | None = None,
    random_effects: list[str] | None = None,
    offset: str = None,
    **kwargs,
) -> tuple[bmb.Model, az.InferenceData]:
    """Fit a binomial proportion model on individual / patient-level data.
    """
    fixed_effects = fixed_effects or []
    random_effects = random_effects or []
    return fit_bambi_model(
        data,
        run_id,
        success,
        family="binomial",
        trials=trials,
        fixed_effects=fixed_effects,
        interaction_effects=interaction_effects,
        random_effects=random_effects,
        offset=offset,
        **kwargs,
    )