"""Thin wrapper around the EpiLink scorer used by the Scotland clustering pipeline."""

from __future__ import annotations

import numpy as np
from epilink import EpiLink, InfectiousnessToTransmission, NaturalHistoryParameters

RNG_SEED = 42


def estimate_epilink_compatibility(
    genetic_distance,
    temporal_distance,
) -> np.ndarray:
    profile = InfectiousnessToTransmission(
        parameters=NaturalHistoryParameters(),
        rng_seed=RNG_SEED,
    )
    epilink = EpiLink(
        transmission_profile=profile,
    )
    genetic = np.asarray(genetic_distance, dtype=float)
    temporal = np.asarray(temporal_distance, dtype=float)
    scores = epilink.score_target(
        sample_time_difference=temporal,
        genetic_distance=genetic,
    )
    return np.asarray(scores, dtype=float)
