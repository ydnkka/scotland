"""Thin wrapper around the EpiLink scorer used by the Scotland clustering pipeline.

Design note
-----------
``estimate_epilink_compatibility`` is called once per (window, lineage) group from a
separate child process launched by GNU parallel.  Because each invocation is its own
Python process the module-level singleton ``_EPILINK_INSTANCE`` is initialised exactly
once per process — avoiding the cost of rebuilding the Monte-Carlo lookup tables
(MC_SAMPLES draws) on every call while keeping the parallel-process isolation intact.
"""

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
    # Symptomatic rate of 37% reflects population-level ascertainment in the
    # Scottish surveillance context; adjust if using a different ascertainment model.
    "symptomatic_rate": 0.37,
    "symptomatic_shape": 1.0,
    "transmission_rate_ratio": 2.29,
    "testing_delay_shape": 1.0,
    "testing_delay_scale": 1.0,
    # Standard SARS-CoV-2 clock rate (~1e-3 substitutions/site/year).
    # Must match the alignment_length used to convert TN93 proportional distances
    # to SNP counts in 03_process_group.py (default 29903 nt).
    "substitution_rate": 1e-3,
    "relaxation": 0.33,
    "genome_length": 29903,  # keep in sync with pipeline.alignment_length in config.yaml
}
# MAXIMUM_DEPTH = 0 : only direct or co-primary transmission links are scored.
# Increase to consider transmission chains of greater depth at the cost of
# wider compatibility intervals and slower scoring.
MUTATION_PROCESS = "stochastic"
MAXIMUM_DEPTH = 0
MC_SAMPLES = 10_000
TARGET: tuple[str, ...] = ("ad(0)", "ca(0,0)")
RNG_SEED = 42

# ---------------------------------------------------------------------------
# Module-level singleton — built once per process, reused across calls.
# ---------------------------------------------------------------------------
_EPILINK_INSTANCE: EpiLink | None = None


def _get_epilink() -> EpiLink:
    """Return the module-level EpiLink instance, constructing it on first call."""
    global _EPILINK_INSTANCE
    if _EPILINK_INSTANCE is None:
        profile = InfectiousnessToTransmission(
            parameters=NaturalHistoryParameters(**NATURAL_HISTORY),
            rng_seed=RNG_SEED,
        )
        _EPILINK_INSTANCE = EpiLink(
            mutation_process=MUTATION_PROCESS,
            transmission_profile=profile,
            maximum_depth=MAXIMUM_DEPTH,
            mc_samples=MC_SAMPLES,
            target=TARGET,
        )
    return _EPILINK_INSTANCE


def estimate_epilink_compatibility(
    genetic_distance,
    temporal_distance,
) -> np.ndarray:
    """Score pair compatibility with EpiLink.

    Parameters
    ----------
    genetic_distance : array-like of float
        Pairwise SNP distances (integer counts cast to float, not proportions).
        Must be computed with the same genome length as ``NATURAL_HISTORY["genome_length"]``.
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
