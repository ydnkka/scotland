"""Thin wrapper around the EpiLink scorer used by the Scotland clustering pipeline.

Exposes :func:`estimate_epilink_compatibility`, which scores pairwise
(genetic distance, temporal distance) observations using a module-level,
parameterised :class:`epilink.EpiLink` instance so the scorer is built once
per process and reused across every call.

The default natural-history / inference parameters match the baseline
configuration used in the EpiLink manuscript evaluation
(see ``projects/epilink/src/evaluation/models.py``). Call
:func:`configure_epilink` once at program start to override defaults
(e.g. mutation process, MC samples, seed) before the first scoring call.
"""

from __future__ import annotations

from collections.abc import Mapping
from threading import Lock
from typing import Any

import numpy as np

from epilink import EpiLink, InfectiousnessToTransmission, NaturalHistoryParameters


# ---------------------------------------------------------------------------
# Defaults (match projects/epilink/src/evaluation baseline)
# ---------------------------------------------------------------------------

DEFAULT_NATURAL_HISTORY: dict[str, Any] = {
    "incubation_shape": 5.807,
    "incubation_scale": 0.948,
    "latent_shape": 3.38,
    "symptomatic_rate": 0.37,
    "symptomatic_shape": 1.0,
    "transmission_rate_ratio": 2.29,
    "testing_delay_shape": 1.0,
    "testing_delay_scale": 1.0,
    "substitution_rate": 1e-3,
    "relaxation": 0.33,
    "genome_length": 29903,
}

DEFAULT_MUTATION_PROCESS = "stochastic"
DEFAULT_MAXIMUM_DEPTH = 0
DEFAULT_MC_SAMPLES = 10_000
DEFAULT_TARGET: tuple[str, ...] = ("ad(0)", "ca(0,0)")
DEFAULT_RNG_SEED = 42


# ---------------------------------------------------------------------------
# Shared scorer (built lazily, rebuilt on configure_epilink)
# ---------------------------------------------------------------------------

_EPILINK_INSTANCE: EpiLink | None = None
_BUILD_LOCK = Lock()


def _build_natural_history_parameters(
    parameters: Mapping[str, Any] | None,
) -> NaturalHistoryParameters:
    merged = {**DEFAULT_NATURAL_HISTORY, **dict(parameters or {})}
    return NaturalHistoryParameters(
        incubation_shape=float(merged["incubation_shape"]),
        incubation_scale=float(merged["incubation_scale"]),
        latent_shape=float(merged["latent_shape"]),
        symptomatic_rate=float(merged["symptomatic_rate"]),
        symptomatic_shape=float(merged["symptomatic_shape"]),
        transmission_rate_ratio=float(merged["transmission_rate_ratio"]),
        testing_delay_shape=float(merged["testing_delay_shape"]),
        testing_delay_scale=float(merged["testing_delay_scale"]),
        substitution_rate=float(merged["substitution_rate"]),
        relaxation=float(merged["relaxation"]),
        genome_length=int(merged["genome_length"]),
    )


def configure_epilink(
    natural_history: Mapping[str, Any] | None = None,
    *,
    mutation_process: str = DEFAULT_MUTATION_PROCESS,
    maximum_depth: int = DEFAULT_MAXIMUM_DEPTH,
    mc_samples: int = DEFAULT_MC_SAMPLES,
    target: tuple[str, ...] = DEFAULT_TARGET,
    rng_seed: int = DEFAULT_RNG_SEED,
) -> EpiLink:
    """Build (or rebuild) the shared EpiLink instance and return it.

    Parameters
    ----------
    natural_history : Mapping, optional
        Overrides applied on top of :data:`DEFAULT_NATURAL_HISTORY`.
    mutation_process : str
        ``"deterministic"`` or ``"stochastic"``.
    maximum_depth, mc_samples : int
        Forwarded to :class:`epilink.EpiLink`.
    target : tuple[str, ...]
        EpiLink target labels (defaults to ``("ad(0)", "ca(0,0)")``).
    rng_seed : int
        Seed for the underlying ``InfectiousnessToTransmission`` sampler.
    """
    global _EPILINK_INSTANCE
    nhp = _build_natural_history_parameters(natural_history)
    transmission_profile = InfectiousnessToTransmission(parameters=nhp, rng_seed=rng_seed)
    with _BUILD_LOCK:
        _EPILINK_INSTANCE = EpiLink(
            mutation_process=mutation_process,
            transmission_profile=transmission_profile,
            maximum_depth=int(maximum_depth),
            mc_samples=int(mc_samples),
            target=tuple(target),
        )
    return _EPILINK_INSTANCE


def get_epilink() -> EpiLink:
    """Return the shared EpiLink instance, building it with defaults on first use."""
    if _EPILINK_INSTANCE is None:
        with _BUILD_LOCK:
            if _EPILINK_INSTANCE is None:
                configure_epilink()
    return _EPILINK_INSTANCE  # type: ignore[return-value]


def estimate_epilink_compatibility(
    genetic_distance,
    temporal_distance,
) -> np.ndarray:
    """Score pair compatibility with the shared EpiLink instance.

    Parameters
    ----------
    genetic_distance : array-like of float
        Pairwise SNP distances (integer counts cast to float, not proportions).
    temporal_distance : array-like of float
        Absolute sampling-date differences in days.

    Returns
    -------
    numpy.ndarray
        EpiLink compatibility scores in ``[0, 1]``, one per input pair.
    """
    scorer = get_epilink()
    genetic = np.asarray(genetic_distance, dtype=float)
    temporal = np.asarray(temporal_distance, dtype=float)
    scores = scorer.score_target(
        sample_time_difference=temporal,
        genetic_distance=genetic,
    )
    return np.asarray(scores, dtype=float)
