"""Shared building blocks for the Part 1 analysis.

This package separates the Part 1 codebase along the same lines as the
manuscript figure script ``manuscript/make_figures.py``:

* :mod:`constants`     — labels, term sets, model specifications.
* :mod:`data_prep`     — sequence loading and cluster-table construction.
* :mod:`estimators`    — low-level fitting primitives (ZTNB, sandwich SEs).
* :mod:`fit_models`    — high-level model fits orchestrated for each
                         outcome / line of inquiry.
* :mod:`inspect_plots` — quick-look plots used during model checking.

The top-level ``*.py`` scripts are thin orchestrators that compose
these building blocks for the two lines of inquiry in Part 1:

* **Line 1** — area-level SIMD deprivation as the exposure for cluster
  outcomes (cluster size, geographic spread) and for within-cluster
  excess mixing.
* **Line 2** — within-cluster excess mixing as a predictor of cluster
  scale (size and geographic spread).
"""
