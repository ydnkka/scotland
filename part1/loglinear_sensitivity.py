"""Log-linear sensitivity for Part 1 count outcomes.

Single-component log-linear (Poisson-style) OLS fits on log-cluster size
and log-distinct-datazones, contrasted in the manuscript supplement with
the corresponding hurdle and ZTNB positive-count components from the
two-part main model.  Two runs are produced: the primary covariate set
(Line 1) and the primary set augmented with the excess-mixing predictors
(Line 2).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

# Allow ``python part1/loglinear_sensitivity.py`` to import lib/.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.constants import MIXING_PREDICTOR_TERMS  # noqa: E402
from lib.data_prep import (  # noqa: E402
    ensure_mixing_predictor_columns,
    load_cluster_table,
    repo_root,
)
from lib.fit_models import fit_loglinear_models  # noqa: E402


def run(root: Path) -> None:
    main_dir = root / "part1"
    tables_dir = main_dir / "tables"
    clusters = ensure_mixing_predictor_columns(
        load_cluster_table(root=root, cache_dir=main_dir / "cache")
    )

    # Line 1 — primary covariates only.
    results = fit_loglinear_models(clusters)
    out = tables_dir / "loglinear_count_model_results.csv"
    results.to_csv(out, index=False)
    print(f"Wrote {out}", flush=True)

    # Line 2 — primary covariates plus excess-mixing predictors.
    mixing_predictor_results = fit_loglinear_models(
        clusters,
        extra_terms=MIXING_PREDICTOR_TERMS,
        predictor_set="primary_plus_mixing",
    )
    mixing_out = tables_dir / "mixing_predictor_loglinear_count_model_results.csv"
    mixing_predictor_results.to_csv(mixing_out, index=False)
    print(f"Wrote {mixing_out}", flush=True)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    run(args.root.resolve())


if __name__ == "__main__":
    main()
