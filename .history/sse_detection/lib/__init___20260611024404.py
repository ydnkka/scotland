"""Public helpers for the SSE detection workflow."""

from .entropy import (
    MIXING_TERTILE_FEATURES,
    MIXING_TERTILE_ORDER,
    DEFAULT_MIXING_FEATURES,
    add_observed_mixing_entropy_scales,
    add_mixing_tertiles,
    cluster_age_conditional_binary_entropy,
    cluster_socio_demo_entropy,
    observed_mixing_entropy_scales,
    onward_edge_entropy,
    vaccination_mixing_features,
)
from .io import HIGH_PRIORITY_CANDIDATE_TIERS, SseOutputs, load_sse_outputs

from .sse_detection import load_sequence_data

from .palettes import (
    BACKGROUND_COLOR,
    BACKGROUND_DARK,
    BORDER,
    CANDIDATE_COLOR,
    CANDIDATE_DARK,
    GRAY,
    GRAY_LIGHT,
    GRID,
    INK,
    INK_SOFT,
    ORANGE_DARK,
    SSE_SIGNATURE_ORDER,
    SSE_SIGNATURE_PALETTE,
    TEAL_DARK,
)
    __all__ = [
    "BACKGROUND_COLOR",
    "BACKGROUND_DARK",
    "BORDER",
    "CANDIDATE_COLOR",
    "CANDIDATE_DARK",
    "GRAY",
    "GRAY_LIGHT",
    "GRID",
    "INK",
    "INK_SOFT",
    "ORANGE_DARK",
    "SSE_SIGNATURE_ORDER",
    "SSE_SIGNATURE_PALETTE",
    "TEAL_DARK",
]

__all__ = [
    "HIGH_PRIORITY_CANDIDATE_TIERS",
    "MIXING_TERTILE_FEATURES",
    "DEFAULT_MIXING_FEATURES",
    "MIXING_TERTILE_ORDER",
    "add_observed_mixing_entropy_scales",
    "SseOutputs",
    "add_mixing_tertiles",
    "cluster_age_conditional_binary_entropy",
    "cluster_socio_demo_entropy",
    "load_sse_outputs",
    "observed_mixing_entropy_scales",
    "onward_edge_entropy",
    "vaccination_mixing_features",
    "load_sequence_data",
]
