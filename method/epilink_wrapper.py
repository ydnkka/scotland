"""Thin wrapper around the EpiLink scorer used by the Scotland clustering pipeline."""

from __future__ import annotations

from typing import Any

import numpy as np

from epilink import EpiLink, InfectiousnessToTransmission, NaturalHistoryParameters


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

NATURAL_HISTORY: dict[str, Any] = {
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
MUTATION_PROCESS = "stochastic"
MAXIMUM_DEPTH = 0
MC_SAMPLES = 10_000
TARGET: tuple[str, ...] = ("ad(0)", "ca(0,0)")
RNG_SEED = 42


def _get_epilink() -> EpiLink:
    """Return EpiLink instance."""
    profile = InfectiousnessToTransmission(
        parameters=NaturalHistoryParameters(**NATURAL_HISTORY),
        rng_seed=RNG_SEED
    )
    epilink = EpiLink(
        mutation_process=MUTATION_PROCESS,
        transmission_profile=profile,
        maximum_depth=MAXIMUM_DEPTH,
        mc_samples=MC_SAMPLES,
        target=TARGET,
    )
    return epilink


def estimate_epilink_compatibility(
    genetic_distance,
    temporal_distance,
) -> np.ndarray:
    """Score pair compatibility with the with EpiLink.

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
    epilink = _get_epilink()
    genetic = np.asarray(genetic_distance, dtype=float)
    temporal = np.asarray(temporal_distance, dtype=float)
    scores = epilink.score_target(
        sample_time_difference=temporal,
        genetic_distance=genetic,
    )
    return np.asarray(scores, dtype=float)

