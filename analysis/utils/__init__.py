"""Shared utilities for Scotland clustering chapters.

Submodules
----------
data
    Data loading and preprocessing helpers that use `config.yaml` paths and
    read the processed parquet file in a memory-conscious way.
style
    Publication-grade matplotlib theme, figure size presets, and consistent
    categorical colour palettes for SIMD quintiles, WHO variants, Leiden
    resolutions, and Nextclade clades.
stats
    Light statistical helpers: negative-binomial cluster-size regression,
    bootstrapped confidence intervals, singleton odds ratios, and tidy
    coefficient tables for forest plots.
"""

from . import data, stats, style  # noqa: F401
