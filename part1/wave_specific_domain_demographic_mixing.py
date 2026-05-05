"""Wave-specific SIMD domain effects on demographic cluster mixing."""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from .simd_domain_analysis import DOMAINS
    from .simd_domain_demographic_mixing import (
        MIXING_OUTCOMES,
        DomainMixingSpec,
        build_cluster_table,
        fit_model,
        read_sequence_rows,
    )
    from .cluster_outcome_models import categorical_levels, repo_root
except ImportError:
    from simd_domain_analysis import DOMAINS
    from simd_domain_demographic_mixing import (
        MIXING_OUTCOMES,
        DomainMixingSpec,
        build_cluster_table,
        fit_model,
        read_sequence_rows,
    )
    from cluster_outcome_models import categorical_levels, repo_root


QC_DEFAULT = "good"

WAVE_ORDER = [
    "B.1.177",
    "Alpha",
    "Delta",
    "BA.1",
    "BA.2",
    "BA.4",
    "BA.5",
    "BQ.1",
    "XBB",
]

WAVE_LABELS = {
    "B.1.177": "B.1.177",
    "Alpha": "Alpha",
    "Delta": "Delta",
    "BA.1": "BA.1",
    "BA.2": "BA.2",
    "BA.4": "BA.4",
    "BA.5": "BA.5",
    "BQ.1": "BQ.1",
    "XBB": "XBB",
    "Other": "Other",
}


def assign_wave(lineage: str) -> str:
    if not isinstance(lineage, str):
        return "Other"
    if lineage.startswith("B.1.177"):
        return "B.1.177"
    if lineage == "B.1.1.7" or lineage.startswith("B.1.1.7."):
        return "Alpha"
    if lineage.startswith("AY.") or lineage == "B.1.617.2":
        return "Delta"
    if lineage.startswith("BA.1"):
        return "BA.1"
    if lineage.startswith("BA.2"):
        return "BA.2"
    if lineage.startswith("BA.4"):
        return "BA.4"
    if lineage.startswith("BA.5") or lineage.startswith("BE."):
        return "BA.5"
    if lineage.startswith("BQ."):
        return "BQ.1"
    if lineage.startswith("XBB"):
        return "XBB"
    return "Other"


def add_wave(seq: pd.DataFrame) -> pd.DataFrame:
    seq = seq.copy()
    seq["wave_group"] = seq["pango_lineage"].astype(str).map(assign_wave).astype("category")
    return seq


def fit_wave_specific_models(
    clusters: pd.DataFrame,
    *,
    min_clusters: int,
    min_windows: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    multi = clusters[
        (clusters["cluster_size"] >= 2)
        & (clusters["wave_group"].isin(WAVE_ORDER))
    ].copy()

    results = []
    diagnostics = []
    for wave in WAVE_ORDER:
        wave_df = multi[multi["wave_group"] == wave].copy()
        n_windows = wave_df["window_id"].nunique()
        if len(wave_df) < min_clusters or n_windows < min_windows:
            diagnostics.append(
                {
                    "wave_group": wave,
                    "wave_label": WAVE_LABELS[wave],
                    "skipped": True,
                    "reason": "below minimum clusters/windows",
                    "n_clusters": int(len(wave_df)),
                    "n_windows": int(n_windows),
                }
            )
            continue

        levels = categorical_levels(wave_df)
        for domain in DOMAINS:
            for mixing in MIXING_OUTCOMES:
                spec = DomainMixingSpec(
                    domain=domain,
                    mixing=mixing,
                    outcome=f"{mixing}_excess_discordance",
                )
                result, diag = fit_model(wave_df, spec, levels)
                result.insert(0, "wave_group", wave)
                result.insert(1, "wave_label", WAVE_LABELS[wave])
                diag.update(
                    {
                        "wave_group": wave,
                        "wave_label": WAVE_LABELS[wave],
                        "skipped": False,
                        "reason": "",
                    }
                )
                diagnostics.append(diag)
                results.append(result)

    return pd.concat(results, ignore_index=True), pd.DataFrame(diagnostics)


def plot_wave_heatmaps(results: pd.DataFrame, out_base: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    domain_terms = {
        domain: f"{domain}_deprivation_z"
        for domain in DOMAINS
    }
    plot_df = pd.concat(
        [
            results[(results["domain"] == domain) & (results["term"] == term)]
            for domain, term in domain_terms.items()
        ],
        ignore_index=True,
    )

    domains = list(DOMAINS)
    waves = [wave for wave in WAVE_ORDER if wave in set(plot_df["wave_group"])]
    vmax = max(0.5, float(np.nanmax(np.abs(plot_df["coefficient_percentage_points"]))))

    fig, axes = plt.subplots(1, len(MIXING_OUTCOMES), figsize=(10.5, 4.2), sharey=True)
    if len(MIXING_OUTCOMES) == 1:
        axes = [axes]

    for ax, mixing in zip(axes, MIXING_OUTCOMES):
        matrix = np.full((len(domains), len(waves)), np.nan)
        sub = plot_df[plot_df["mixing"] == mixing]
        for i, domain in enumerate(domains):
            for j, wave in enumerate(waves):
                row = sub[(sub["domain"] == domain) & (sub["wave_group"] == wave)]
                if not row.empty:
                    matrix[i, j] = row.iloc[0]["coefficient_percentage_points"]

        im = ax.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_title(MIXING_OUTCOMES[mixing]["short_label"])
        ax.set_xticks(range(len(waves)))
        ax.set_xticklabels([WAVE_LABELS[w] for w in waves], rotation=45, ha="right")
        ax.set_yticks(range(len(domains)))
        ax.set_yticklabels([DOMAINS[d]["label"] for d in domains])

    cbar = fig.colorbar(im, ax=axes, shrink=0.72, pad=0.02)
    cbar.set_label("pp per 1 SD", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.subplots_adjust(left=0.12, right=0.9, bottom=0.24, wspace=0.08)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def run(root: Path, qc: str | None, min_clusters: int, min_windows: int) -> None:
    tables_dir = root / "part1" / "tables"
    figures_dir = root / "part1" / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("Reading sequence rows", flush=True)
    seq = add_wave(read_sequence_rows(root, qc))
    print(f"Building wave-specific cluster table from {len(seq):,} rows", flush=True)
    clusters, scaling = build_cluster_table(seq)

    print("Fitting wave-specific domain-demographic mixing models", flush=True)
    results, diagnostics = fit_wave_specific_models(
        clusters,
        min_clusters=min_clusters,
        min_windows=min_windows,
    )

    results.to_csv(
        tables_dir / "wave_specific_domain_demographic_mixing_model_results.csv",
        index=False,
    )
    diagnostics.to_csv(
        tables_dir / "wave_specific_domain_demographic_mixing_model_diagnostics.csv",
        index=False,
    )
    scaling.to_csv(
        tables_dir / "wave_specific_domain_demographic_mixing_covariate_scaling.csv",
        index=False,
    )

    del seq, clusters
    gc.collect()

    plot_wave_heatmaps(
        results,
        figures_dir / "wave_specific_domain_demographic_mixing_effects",
    )

    print(
        f"Wrote {tables_dir / 'wave_specific_domain_demographic_mixing_model_results.csv'}",
        flush=True,
    )
    print(
        f"Wrote {figures_dir / 'wave_specific_domain_demographic_mixing_effects.png'}",
        flush=True,
    )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument("--qc", default=QC_DEFAULT)
    parser.add_argument("--min-clusters", type=int, default=5000)
    parser.add_argument("--min-windows", type=int, default=3)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    qc = None if str(args.qc).lower() == "none" else str(args.qc)
    run(args.root.resolve(), qc, args.min_clusters, args.min_windows)


if __name__ == "__main__":
    main()
