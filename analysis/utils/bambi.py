import os

import arviz as az
import json
import bambi as bmb


def fit_bambi_model(
    cluster_data,
    run_id: int | str,
    *,
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
    Fit a Bambi negative-binomial hierarchical model and cache results to disk.

    Fits the model using NUTS via PyMC. On subsequent calls with the same
    ``run_id``, the saved model and trace are loaded from ``save_dir`` instead
    of refitting, unless ``force_refit=True``.

    Parameters
    ----------
    cluster_data : pd.DataFrame
        Input data. Must contain all variables referenced in ``formula``,
        including ``epoch``, ``simd_quintile_mode``, ``window_id``,
        ``n_sequences_shifted``, and ``log_seq_prop`` (the offset).
    run_id : int or str
        Unique identifier for this run, used to name the cached output files.
        Use distinct values to cache results for different formulae or data
        subsets without overwriting each other.
    formula : str, optional
        R-style Bambi formula. Categorical reference levels must be specified
        using ``Treatment()``, e.g. ``C(epoch, Treatment('Delta'))``.
        Defaults to::

            n_sequences_shifted
                ~ C(simd_quintile_mode, Treatment(1))
                + C(epoch, Treatment('Delta'))
                + (1 | window_id)
                + offset(log_seq_prop)

    save_dir : str, default 'bambi_outputs'
        Directory in which to save (or look for) the fitted model and trace.
        Created automatically if it does not exist.
    draws : int, default 1000
        Number of posterior samples drawn per chain after tuning.
    tune : int, default 1000
        Number of tuning (warm-up) steps per chain, discarded after sampling.
    target_accept : float, default 0.9
        Target acceptance rate for the NUTS sampler. Increase towards 1.0
        if divergences are observed; values above 0.95 may slow sampling.
    chains : int, default 4
        Number of independent Markov chains. Must be ≥ 2 for reliable
        R-hat convergence diagnostics.
    random_seed : int, default 42
        Seed passed to PyMC for reproducibility.
    force_refit : bool, default False
        If ``True``, refit and overwrite any existing cached results for
        this ``run_id``. If ``False``, return cached results when available.

    Returns
    -------
    model : bmb.Model
        The compiled Bambi model object.
    trace : az.InferenceData
        ArviZ ``InferenceData`` object containing the posterior samples,
        suitable for use with ``az.summary()``, ``az.plot_trace()``, etc.

    Raises
    ------
    FileNotFoundError
        If ``force_refit=False`` and expected cache files are missing or
        incomplete.
    ValueError
        If ``formula`` references columns not present in ``cluster_data``.
    """
    if formula is None:
        formula = (
            "n_sequences_shifted ~ C(simd_quintile_mode, Treatment(3))"
            " + C(epoch, Treatment('Delta')) + (1|window_id) + offset(log_seq_prop)"
        )

    os.makedirs(save_dir, exist_ok=True)
    trace_path = os.path.join(save_dir, f'bambi_trace_{run_id}.nc')
    formula_path = os.path.join(save_dir, f'bambi_formula_{run_id}.json')
    summary_path = os.path.join(save_dir, f'bambi_summary_{run_id}.csv')

    # --- Load from disk if available ---
    if not force_refit and os.path.exists(trace_path) and os.path.exists(formula_path):
        print("Found saved trace — loading from disk.")
        trace = az.from_netcdf(trace_path)
        with open(formula_path) as f:
            saved = json.load(f)
        # Rebuild model (fast — no sampling)
        cluster_data = cluster_data.copy()
        cluster_data['n_sequences_shifted'] = cluster_data['n_sequences'] - 1
        model = bmb.Model(saved['formula'], data=cluster_data, family='negativebinomial')
        model.build()
        print(f"  Trace : {trace_path}")
        print(f"  Formula: {saved['formula']}")
        return model, trace

    # --- Prepare data ---
    cluster_data = cluster_data.copy()
    cluster_data['n_sequences_shifted'] = cluster_data['n_sequences'] - 1

    print(f"Formula  : {formula}")
    print(f"Obs      : {len(cluster_data):,}")
    print(f"Windows  : {cluster_data['window_id'].nunique()}")
    print(f"Epochs   : {cluster_data['epoch'].nunique()}")

    # --- Build and sample ---
    model = bmb.Model(formula, data=cluster_data, family='negativebinomial')
    model.build()

    trace = model.fit(
        draws=draws,
        tune=tune,
        target_accept=target_accept,
        inference_method='numpyro',
        chains=chains,
        random_seed=random_seed,
    )
    summary = az.summary(trace)
    summary['parameter'] = summary.index
    summary['run_id'] = run_id
    summary.reset_index(drop=True, inplace=True)


    # --- Save trace + summary + formula ---
    print(f"\nSaving trace   → {trace_path}")
    az.to_netcdf(trace, trace_path)

    print(f"Saving summary → {summary_path}")
    summary.to_csv(summary_path, index=False)

    print(f"Saving formula → {formula_path}")
    with open(formula_path, 'w') as f:
        json.dump({'formula': formula}, f)

    print("Done.")
    return model, trace
