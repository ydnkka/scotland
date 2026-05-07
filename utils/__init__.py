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
policy
    Ordered policy-period helpers for Scotland COVID-19 restriction phases.
"""

from importlib import import_module

__all__ = ["data", "policy", "style"]


def __getattr__(name: str):
    if name in __all__:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
