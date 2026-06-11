#!/usr/bin/env python3
from __future__ import annotations

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
    "ROLE_ORDER",
    "ROLE_PALETTE",
    "SSE_CATEGORY_ORDER",
    "SSE_CATEGORY_PALETTE",
    "SSE_SIGNATURE_ORDER",
    "SSE_SIGNATURE_PALETTE",
    "TEAL_DARK",
]


CANDIDATE_COLOR = "#C44E52"
CANDIDATE_DARK = "#9C0207"
BACKGROUND_COLOR = "#B0B0B0"
BACKGROUND_DARK = "#444343"
BORDER = "#D0D7DE"
GRAY = "#808080"
GRAY_LIGHT = "#6C757D"
GRID = "#E9ECEF"
INK = "#212529"
INK_SOFT = "#495057"
ORANGE_DARK = "#D95F02"
TEAL_DARK = "#007C89"

SSE_SIGNATURE_ORDER: list[str] = ["burst", "burden", "burst+burden"]
SSE_SIGNATURE_PALETTE: dict[str, str] = {
    "burst": CANDIDATE_COLOR,
    "burden": TEAL_DARK,
    "burst+burden": ORANGE_DARK,
}

