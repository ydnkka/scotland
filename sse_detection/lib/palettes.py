#!/usr/bin/env python3
from __future__ import annotations

__all__ = [
    "BACKGROUND_COLOR",
    "BACKGROUND_DARK",
    "CANDIDATE_COLOR",
    "CANDIDATE_DARK",
    "ROLE_ORDER",
    "ROLE_PALETTE",
]


CANDIDATE_COLOR = "#C44E52"
CANDIDATE_DARK = "#9C0207"
BACKGROUND_COLOR = "#B0B0B0"
BACKGROUND_DARK = "#444343"

ROLE_ORDER: list[str] = [
    "putative_birth",
    "relay_amplifier",
    "merged_relay",
    "terminal_sink",
    "isolated_burst",
    "unclear_origin",
]

ROLE_PALETTE: dict[str, str] = {
    "putative_birth": "#0072B2",
    "relay_amplifier": "#009E73",
    "merged_relay": "#D55E00",
    "terminal_sink": "#CC79A7",
    "isolated_burst": "#E69F00",
    "unclear_origin": "#6C6F73",
}
