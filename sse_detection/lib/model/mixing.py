"""Run Bayesian SSE mixing regression models.

Examples
--------
Run all mixing models:
    python -m sse_detection.lib.model.mixing

Run one model version:
    python -m sse_detection.lib.model.mixing --family logistic --model-set null_primary
"""

from __future__ import annotations

from pathlib import Path
import sys


if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from sse_detection.lib.model.runner import main_for_domain
else:
    from .runner import main_for_domain


def main() -> int:
    """Run the mixing-model CLI."""
    return main_for_domain("mixing")


if __name__ == "__main__":
    raise SystemExit(main())
