"""Run Bayesian SSE sequence-composition regression models.

Examples
--------
Run all composition models:
    python -m chapter_analyses.sse_detection.lib.model.composition

Run one model version:
    python -m chapter_analyses.sse_detection.lib.model.composition --family linear --outcome burden_score --model-set expanded
"""

from __future__ import annotations

import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
    from chapter_analyses.sse_detection.lib.model.runner import main_for_domain
else:
    from .runner import main_for_domain


def main() -> int:
    """Run the composition-model CLI."""
    return main_for_domain("composition")


if __name__ == "__main__":
    raise SystemExit(main())
