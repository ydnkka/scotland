"""SIMD domain deprivation as predictors of demographic cluster mixing.

The outcomes are excess pairwise discordance for age band, sex, and joint
age-sex profile. Each SIMD domain is modelled one at a time as the exposure,
with lineage, window, Leiden resolution, local incidence, sequencing intensity,
test positivity, and cluster size adjusted.
"""

from __future__ import annotations

import argparse
import gc
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.linalg import pinvh
from scipy.stats import norm

try:
    from .cluster_mixing_analysis import (
        expected_stratum_discordance,
        observed_cluster_discordance,
    )
    from .cluster_outcome_models import (
        analysis_dataset_path,
        build_design_matrix,
        categorical_levels,
        clustered_ols_covariance,
        logit_clipped,
        repo_root,
        zscore,
    )
    from .simd_domain_analysis import DOMAINS
except ImportError:
    from cluster_mixing_analysis import (
        expected_stratum_discordance,
        observed_cluster_discordance,
    )
    from cluster_outcome_models import (
        analysis_dataset_path,
        build_design_matrix,
        categorical_levels,
        clustered_ols_covariance,
        logit_clipped,
        repo_root,
        zscore,
    )
    from simd_domain_analysis import DOMAINS


QC_DEFAULT = "good"

SEQUENCE_COLUMNS = [
    "cluster_id",
    "sequence_id",
    "resolution",
    "window_id",
    "window_idx",
    "wn_mid_date",
    "wn_prop_sequenced",
    "collection_date",
    "datazone",
    "pango_lineage",
    "nextclade_qc",
    "age_band",
    "sex",
    "dz_cum_incidence_per_capita",
    "dz_cum_prop_sequenced",
    "dz_7d_test_positivity",
]

MIXING_OUTCOMES = {
    "age": {
        "column": "age_band",
        "label": "Age-band mixing",
        "short_label": "Age",
    },
    "sex": {
        "column": "sex",
        "label": "Sex mixing",
        "short_label": "Sex",
    },
    "age_sex": {
        "column": "age_sex_profile",
        "label": "Joint age-sex profile mixing",
        "short_label": "Age-sex",
    },
}

SHARED_TERMS = [
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
    "log_cluster_size_z",
]


@dataclass(frozen=True)
class DomainMixingSpec:
    domain: str
    mixing: str
    outcome: str


def read_sequence_rows(root: Path, qc: str | None) -> pd.DataFrame:
    rank_cols = [spec["rank_col"] for spec in DOMAINS.values()]
    columns = list(dict.fromkeys([*SEQUENCE_COLUMNS, *rank_cols]))
    filters = None if qc is None else [("nextclade_qc", "==", qc)]
    df = pd.read_parquet(
        analysis_dataset_path(root),
        columns=columns,
        filters=filters,
        engine="pyarrow",
    )

    for col in [
        "cluster_id",
        "sequence_id",
        "window_id",
        "pango_lineage",
        "nextclade_qc",
        "datazone",
        "age_band",
        "sex",
    ]:
        df[col] = df[col].astype("category")

    df["collection_date"] = pd.to_datetime(df["collection_date"])
    df["wn_mid_date"] = pd.to_datetime(df["wn_mid_date"])
    df["resolution_label"] = df["resolution"].map(lambda x: f"{x:.1f}")
    df["age_sex_profile"] = (
        df["age_band"].astype(str) + "|" + df["sex"].astype(str)
    ).astype("category")
    return df


def build_cluster_table(seq: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rank_cols = [spec["rank_col"] for spec in DOMAINS.values()]
    required = [
        "cluster_id",
        "sequence_id",
        "resolution",
        "window_id",
        "window_idx",
        "pango_lineage",
        "dz_cum_incidence_per_capita",
        "dz_cum_prop_sequenced",
        "wn_prop_sequenced",
        "dz_7d_test_positivity",
        *rank_cols,
    ]
    before = len(seq)
    seq = seq.dropna(subset=required).copy()
    dropped = before - len(seq)

    agg = {
        "cluster_size": ("sequence_id", "nunique"),
        "resolution": ("resolution", "first"),
        "resolution_label": ("resolution_label", "first"),
        "window_id": ("window_id", "first"),
        "window_idx": ("window_idx", "first"),
        "wn_mid_date": ("wn_mid_date", "first"),
        "pango_lineage": ("pango_lineage", "first"),
        "mean_local_incidence_per_capita": ("dz_cum_incidence_per_capita", "mean"),
        "mean_local_seq_fraction": ("dz_cum_prop_sequenced", "mean"),
        "mean_window_seq_fraction": ("wn_prop_sequenced", "mean"),
        "mean_test_positivity": ("dz_7d_test_positivity", "mean"),
    }
    if "wave_group" in seq.columns:
        agg["wave_group"] = ("wave_group", "first")
    for domain, spec in DOMAINS.items():
        agg[f"{domain}_mean_rank"] = (spec["rank_col"], "mean")

    clusters = seq.groupby("cluster_id", observed=True, sort=False).agg(**agg).reset_index()

    for prefix, spec in MIXING_OUTCOMES.items():
        obs = observed_cluster_discordance(seq, spec["column"], prefix)
        exp = expected_stratum_discordance(seq, spec["column"], prefix)
        clusters = clusters.merge(obs, on="cluster_id", how="left")
        clusters = clusters.merge(
            exp,
            on=["window_id", "pango_lineage", "resolution_label"],
            how="left",
        )
        clusters[f"{prefix}_excess_discordance"] = (
            clusters[f"{prefix}_discordance"] - clusters[f"{prefix}_expected_discordance"]
        )

    clusters["local_incidence_log"] = np.log1p(
        clusters["mean_local_incidence_per_capita"].clip(lower=0) * 1000
    )
    clusters["local_seq_fraction_logit"] = logit_clipped(clusters["mean_local_seq_fraction"])
    clusters["window_seq_fraction_logit"] = logit_clipped(clusters["mean_window_seq_fraction"])
    clusters["test_positivity_logit"] = logit_clipped(clusters["mean_test_positivity"].fillna(0))
    clusters["log_cluster_size"] = np.log(clusters["cluster_size"])

    scaling_rows = []
    transforms = {
        "local_incidence_z": "local_incidence_log",
        "local_seq_fraction_z": "local_seq_fraction_logit",
        "window_seq_fraction_z": "window_seq_fraction_logit",
        "test_positivity_z": "test_positivity_logit",
        "log_cluster_size_z": "log_cluster_size",
    }
    for z_col, raw_col in transforms.items():
        clusters[z_col], mean, sd = zscore(clusters[raw_col])
        scaling_rows.append(
            {
                "standardised_column": z_col,
                "source_column": raw_col,
                "source_mean": mean,
                "source_sd": sd,
            }
        )

    for domain in DOMAINS:
        raw_col = f"{domain}_deprivation_raw"
        z_col = f"{domain}_deprivation_z"
        clusters[raw_col] = -clusters[f"{domain}_mean_rank"]
        clusters[z_col], mean, sd = zscore(clusters[raw_col])
        scaling_rows.append(
            {
                "standardised_column": z_col,
                "source_column": raw_col,
                "source_mean": mean,
                "source_sd": sd,
            }
        )

    scaling = pd.DataFrame(scaling_rows)
    scaling.attrs["dropped_sequence_rows_missing_model_fields"] = dropped
    return clusters, scaling


def fit_model(
    df: pd.DataFrame,
    spec: DomainMixingSpec,
    levels: dict[str, list],
) -> tuple[pd.DataFrame, dict]:
    domain_term = f"{spec.domain}_deprivation_z"
    terms = [domain_term, *SHARED_TERMS]
    use = df.dropna(subset=[spec.outcome, *terms]).copy()
    X, feature_names = build_design_matrix(use, terms, levels)
    y = use[spec.outcome].to_numpy(dtype=np.float64)

    xtx = (X.T @ X).toarray()
    xty = X.T @ y
    beta = pinvh(xtx, rtol=1e-10) @ xty
    fitted = X @ beta
    residuals = y - fitted
    cov, residual_sd = clustered_ols_covariance(
        X,
        residuals,
        groups=use["window_id"].astype(str).to_numpy(),
    )
    se = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
    idx = {name: i for i, name in enumerate(feature_names)}

    rows = []
    for term in terms:
        i = idx[term]
        coef = float(beta[i])
        stderr = float(se[i])
        z = coef / stderr if stderr > 0 else np.nan
        rows.append(
            {
                "domain": spec.domain,
                "domain_label": DOMAINS[spec.domain]["label"],
                "mixing": spec.mixing,
                "mixing_label": MIXING_OUTCOMES[spec.mixing]["label"],
                "term": term,
                "term_label": (
                    f"{DOMAINS[spec.domain]['label']} deprivation"
                    if term == domain_term
                    else term
                ),
                "coefficient_excess_discordance": coef,
                "coefficient_percentage_points": coef * 100,
                "std_error_clustered_by_window": stderr,
                "z": z,
                "p_value": float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan,
                "ci_low": coef - 1.96 * stderr,
                "ci_high": coef + 1.96 * stderr,
                "ci_low_percentage_points": (coef - 1.96 * stderr) * 100,
                "ci_high_percentage_points": (coef + 1.96 * stderr) * 100,
            }
        )

    ssr = float(np.sum(residuals**2))
    sst = float(np.sum((y - y.mean()) ** 2))
    diagnostics = {
        "domain": spec.domain,
        "domain_label": DOMAINS[spec.domain]["label"],
        "mixing": spec.mixing,
        "mixing_label": MIXING_OUTCOMES[spec.mixing]["label"],
        "n_clusters": int(len(use)),
        "n_features": int(X.shape[1]),
        "residual_sd": residual_sd,
        "r2": 1 - ssr / sst,
    }
    return pd.DataFrame(rows), diagnostics


def fit_domain_demographic_models(clusters: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    multi = clusters[clusters["cluster_size"] >= 2].copy()
    levels = categorical_levels(multi)
    results = []
    diagnostics = []
    for domain in DOMAINS:
        for mixing in MIXING_OUTCOMES:
            spec = DomainMixingSpec(
                domain=domain,
                mixing=mixing,
                outcome=f"{mixing}_excess_discordance",
            )
            result, diag = fit_model(multi, spec, levels)
            results.append(result)
            diagnostics.append(diag)
    return pd.concat(results, ignore_index=True), pd.DataFrame(diagnostics)


def plot_effects(results: pd.DataFrame, out_base: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    domains = list(DOMAINS)
    outcomes = list(MIXING_OUTCOMES)
    colours = {
        "age": "#4e79a7",
        "sex": "#f28e2b",
        "age_sex": "#59a14f",
    }
    offsets = {"age": -0.22, "sex": 0.0, "age_sex": 0.22}
    y_pos = {domain: i for i, domain in enumerate(domains)}

    plot_rows = []
    for domain in domains:
        term = f"{domain}_deprivation_z"
        plot_rows.append(results[(results["domain"] == domain) & (results["term"] == term)])
    plot_df = pd.concat(plot_rows, ignore_index=True)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for _, row in plot_df.iterrows():
        y = y_pos[row["domain"]] + offsets[row["mixing"]]
        ax.plot(
            [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
            [y, y],
            color=colours[row["mixing"]],
            linewidth=1.2,
        )
        ax.scatter(
            row["coefficient_percentage_points"],
            y,
            color=colours[row["mixing"]],
            s=18,
            label=MIXING_OUTCOMES[row["mixing"]]["short_label"],
        )

    ax.axvline(0, color="#666666", linewidth=0.8, linestyle="--")
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels([DOMAINS[d]["label"] for d in domains])
    ax.set_xlabel(
        "Change in demographic excess mixing, percentage points\n"
        "per 1 SD higher domain deprivation"
    )
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(
        unique.values(),
        unique.keys(),
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=3,
        frameon=False,
    )

    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(bottom=0.3, left=0.2, right=0.98)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def summarise_demographic_mixing(clusters: pd.DataFrame, qc: str | None, dropped: int) -> pd.DataFrame:
    rows = [
        {"measure": "clusters_total", "value": len(clusters)},
        {"measure": "clusters_size_ge_2", "value": int((clusters["cluster_size"] >= 2).sum())},
        {"measure": "sequence_rows_dropped_missing_model_fields", "value": dropped},
        {"measure": "qc_filter", "value": qc or "none"},
    ]
    for prefix, spec in MIXING_OUTCOMES.items():
        values = clusters.loc[
            clusters["cluster_size"] >= 2,
            f"{prefix}_excess_discordance",
        ].dropna()
        desc = values.describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9])
        rows.extend(
            {
                "measure": f"{spec['short_label']} excess_discordance",
                "statistic": key,
                "value": value,
            }
            for key, value in desc.items()
        )
    return pd.DataFrame(rows)


def run(root: Path, qc: str | None) -> None:
    tables_dir = root / "part1" / "tables"
    figures_dir = root / "part1" / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("Reading sequence rows", flush=True)
    seq = read_sequence_rows(root, qc)
    print(f"Building cluster demographic-mixing table from {len(seq):,} rows", flush=True)
    clusters, scaling = build_cluster_table(seq)
    dropped = int(scaling.attrs.get("dropped_sequence_rows_missing_model_fields", 0))

    print(f"Fitting domain-demographic mixing models for {len(clusters):,} clusters", flush=True)
    results, diagnostics = fit_domain_demographic_models(clusters)
    descriptives = summarise_demographic_mixing(clusters, qc, dropped)

    results.to_csv(tables_dir / "simd_domain_demographic_mixing_model_results.csv", index=False)
    diagnostics.to_csv(tables_dir / "simd_domain_demographic_mixing_model_diagnostics.csv", index=False)
    descriptives.to_csv(tables_dir / "simd_domain_demographic_mixing_descriptives.csv", index=False)
    scaling.to_csv(tables_dir / "simd_domain_demographic_mixing_covariate_scaling.csv", index=False)

    del seq, clusters
    gc.collect()

    plot_effects(results, figures_dir / "simd_domain_demographic_mixing_effects")

    print(f"Wrote {tables_dir / 'simd_domain_demographic_mixing_model_results.csv'}", flush=True)
    print(f"Wrote {figures_dir / 'simd_domain_demographic_mixing_effects.png'}", flush=True)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument(
        "--qc",
        default=QC_DEFAULT,
        help="Nextclade QC status to retain. Use 'none' to disable QC filtering.",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    qc = None if str(args.qc).lower() == "none" else str(args.qc)
    run(args.root.resolve(), qc)


if __name__ == "__main__":
    main()
