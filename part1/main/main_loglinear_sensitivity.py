"""Main-formulation log-linear sensitivity models for Part 1 count outcomes."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm

try:
    from .main_analysis import PRIMARY_TERMS, build_exog, lineage_levels, repo_root
except ImportError:
    from main_analysis import PRIMARY_TERMS, build_exog, lineage_levels, repo_root


OUTCOMES = {
    "cluster_size": {
        "label": "Cluster size",
        "source": "cluster_size",
        "log_plus": 0,
    },
    "duration": {
        "label": "Duration",
        "source": "duration_days",
        "log_plus": 1,
    },
    "geographic_dispersion": {
        "label": "Geographic dispersion",
        "source": "cluster_n_datazones",
        "log_plus": 0,
    },
}

TERM_LABELS = {
    "deprivation_z": "Mean SIMD deprivation",
    "local_incidence_z": "Local cumulative incidence",
    "local_seq_fraction_z": "Local sequencing fraction",
    "window_seq_fraction_z": "Window sequencing proportion",
    "test_positivity_z": "Local test positivity",
}


def fit_loglinear_models(clusters: pd.DataFrame) -> pd.DataFrame:
    calendar_cols = [col for col in clusters.columns if col.startswith("calendar_spline_")]
    lineage_levels_all = lineage_levels(clusters)
    frames = []
    for outcome, spec in OUTCOMES.items():
        terms = PRIMARY_TERMS.copy()
        use = clusters.dropna(subset=[spec["source"], *terms, *calendar_cols, "lineage_model"]).copy()
        y_raw = use[spec["source"]].astype(float) + float(spec["log_plus"])
        y = np.log(y_raw)
        x = build_exog(use, terms, calendar_cols, lineage_levels_all)
        groups = use["window_id"].astype(str).to_numpy()
        result = sm.OLS(y, x).fit(cov_type="cluster", cov_kwds={"groups": groups})

        names = list(result.model.exog_names)
        params = np.asarray(result.params, dtype=float)
        bse = np.asarray(result.bse, dtype=float)
        pvalues = np.asarray(result.pvalues, dtype=float)
        idx = {name: i for i, name in enumerate(names)}

        rows = []
        for term in terms:
            i = idx[term]
            coef = float(params[i])
            stderr = float(bse[i])
            rows.append(
                {
                    "model": outcome,
                    "model_label": spec["label"],
                    "outcome": spec["source"],
                    "term": term,
                    "term_label": TERM_LABELS[term],
                    "coefficient_log_ratio": coef,
                    "std_error_clustered_by_window": stderr,
                    "z": coef / stderr if stderr > 0 else np.nan,
                    "p_value": float(pvalues[i]),
                    "geometric_mean_ratio": float(np.exp(coef)),
                    "ci_low": float(np.exp(coef - 1.96 * stderr)),
                    "ci_high": float(np.exp(coef + 1.96 * stderr)),
                    "n_observations": int(len(use)),
                    "r2": float(result.rsquared),
                }
            )
        frames.append(pd.DataFrame(rows))
    return pd.concat(frames, ignore_index=True)


def run(root: Path) -> None:
    main_dir = root / "part1" / "main"
    tables_dir = main_dir / "tables"
    clusters = pd.read_parquet(main_dir / "cache" / "main_cluster_table.parquet")
    results = fit_loglinear_models(clusters)
    out = tables_dir / "main_loglinear_count_model_results.csv"
    results.to_csv(out, index=False)
    print(f"Wrote {out}", flush=True)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    run(args.root.resolve())


if __name__ == "__main__":
    main()
