"""Compare SIMD domains in cluster outcome and mixing models."""

from __future__ import annotations

import argparse
import gc
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml
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


QC_DEFAULT = "good"

DOMAINS = {
    "overall": {
        "label": "Overall",
        "rank_col": "dz_simd_rank",
        "quintile_col": "dz_simd_quintile",
    },
    "income": {
        "label": "Income",
        "rank_col": "dz_simd_income_rank",
    },
    "employment": {
        "label": "Employment",
        "rank_col": "dz_simd_employment_rank",
    },
    "education": {
        "label": "Education",
        "rank_col": "dz_simd_education_rank",
    },
    "health": {
        "label": "Health",
        "rank_col": "dz_simd_health_rank",
    },
    "access": {
        "label": "Access",
        "rank_col": "dz_simd_access_rank",
    },
    "crime": {
        "label": "Crime",
        "rank_col": "dz_simd_crime_rank",
    },
    "housing": {
        "label": "Housing",
        "rank_col": "dz_simd_housing_rank",
    },
}

BASE_SEQUENCE_COLUMNS = [
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
    "dz_cum_incidence_per_capita",
    "dz_cum_prop_sequenced",
    "dz_7d_test_positivity",
]

SHARED_TERMS = [
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
]

OUTCOME_SPECS = {
    "cluster_size": {
        "label": "Cluster size",
        "outcome": "cluster_size",
        "scale": "log",
    },
    "duration": {
        "label": "Duration",
        "outcome": "duration_days_plus1",
        "scale": "log",
    },
    "geographic_dispersion": {
        "label": "Geographic dispersion",
        "outcome": "cluster_n_datazones",
        "scale": "log",
    },
}


@dataclass(frozen=True)
class DomainModelSpec:
    domain: str
    outcome_name: str
    outcome: str
    label: str
    scale: str
    include_size: bool = False


def processed_simd_path(root: Path) -> Path:
    with open(root / "config.yaml", "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return root / cfg["data"]["processed"]["simd"]


def domain_rank_maxima(root: Path) -> dict[str, float]:
    cols = [spec["rank_col"] for spec in DOMAINS.values()]
    simd = pd.read_parquet(processed_simd_path(root), columns=cols)
    return {
        domain: float(simd[spec["rank_col"]].max())
        for domain, spec in DOMAINS.items()
    }


def rank_to_quintile(rank: pd.Series, max_rank: float) -> pd.Series:
    quintile = np.ceil(rank.astype(float) / (max_rank / 5.0))
    return quintile.clip(1, 5).astype("Int64").astype("category")


def read_sequence_rows(root: Path, qc: str | None) -> pd.DataFrame:
    rank_cols = [spec["rank_col"] for spec in DOMAINS.values()]
    columns = BASE_SEQUENCE_COLUMNS + ["dz_simd_quintile", *rank_cols]
    columns = list(dict.fromkeys(columns))
    filters = None if qc is None else [("nextclade_qc", "==", qc)]
    df = pd.read_parquet(
        analysis_dataset_path(root),
        columns=columns,
        filters=filters,
        engine="pyarrow",
    )

    maxima = domain_rank_maxima(root)
    for domain, spec in DOMAINS.items():
        q_col = f"{domain}_domain_quintile"
        if domain == "overall":
            df[q_col] = df["dz_simd_quintile"].astype("category")
        else:
            df[q_col] = rank_to_quintile(df[spec["rank_col"]], maxima[domain])

    for col in ["cluster_id", "sequence_id", "window_id", "pango_lineage", "nextclade_qc", "datazone"]:
        df[col] = df[col].astype("category")

    df["collection_date"] = pd.to_datetime(df["collection_date"])
    df["wn_mid_date"] = pd.to_datetime(df["wn_mid_date"])
    df["resolution_label"] = df["resolution"].map(lambda x: f"{x:.1f}")
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
        "cluster_n_datazones": ("datazone", "nunique"),
        "cluster_start_date": ("collection_date", "min"),
        "cluster_end_date": ("collection_date", "max"),
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
    for domain, spec in DOMAINS.items():
        agg[f"{domain}_mean_rank"] = (spec["rank_col"], "mean")

    clusters = seq.groupby("cluster_id", observed=True, sort=False).agg(**agg).reset_index()
    clusters["duration_days"] = (
        clusters["cluster_end_date"] - clusters["cluster_start_date"]
    ).dt.days.astype(int)
    clusters["duration_days_plus1"] = clusters["duration_days"] + 1

    for domain in DOMAINS:
        q_col = f"{domain}_domain_quintile"
        obs = observed_cluster_discordance(seq, q_col, domain)
        exp = expected_stratum_discordance(seq, q_col, domain)
        clusters = clusters.merge(obs, on="cluster_id", how="left")
        clusters = clusters.merge(
            exp,
            on=["window_id", "pango_lineage", "resolution_label"],
            how="left",
        )
        clusters[f"{domain}_excess_discordance"] = (
            clusters[f"{domain}_discordance"] - clusters[f"{domain}_expected_discordance"]
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


def fit_log_linear(
    df: pd.DataFrame,
    spec: DomainModelSpec,
    levels: dict[str, list],
) -> tuple[pd.DataFrame, dict]:
    domain_term = f"{spec.domain}_deprivation_z"
    numeric_terms = [domain_term, *SHARED_TERMS]
    if spec.include_size:
        numeric_terms.append("log_cluster_size_z")

    use = df.dropna(subset=[spec.outcome, *numeric_terms]).copy()
    X, feature_names = build_design_matrix(use, numeric_terms, levels)
    y_raw = use[spec.outcome].to_numpy(dtype=np.float64)
    y = np.log(y_raw) if spec.scale == "log" else y_raw

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
    for term in numeric_terms:
        i = idx[term]
        coef = float(beta[i])
        stderr = float(se[i])
        z = coef / stderr if stderr > 0 else np.nan
        row = {
            "domain": spec.domain,
            "domain_label": DOMAINS[spec.domain]["label"],
            "model": spec.outcome_name,
            "model_label": spec.label,
            "term": term,
            "term_label": (
                f"{DOMAINS[spec.domain]['label']} deprivation"
                if term == domain_term
                else term
            ),
            "coefficient": coef,
            "std_error_clustered_by_window": stderr,
            "z": z,
            "p_value": float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan,
        }
        if spec.scale == "log":
            row.update(
                {
                    "geometric_mean_ratio": float(np.exp(coef)),
                    "ci_low": float(np.exp(coef - 1.96 * stderr)),
                    "ci_high": float(np.exp(coef + 1.96 * stderr)),
                }
            )
        else:
            row.update(
                {
                    "coefficient_percentage_points": coef * 100,
                    "ci_low": coef - 1.96 * stderr,
                    "ci_high": coef + 1.96 * stderr,
                    "ci_low_percentage_points": (coef - 1.96 * stderr) * 100,
                    "ci_high_percentage_points": (coef + 1.96 * stderr) * 100,
                }
            )
        rows.append(row)

    ssr = float(np.sum(residuals**2))
    sst = float(np.sum((y - y.mean()) ** 2))
    diag = {
        "domain": spec.domain,
        "domain_label": DOMAINS[spec.domain]["label"],
        "model": spec.outcome_name,
        "model_label": spec.label,
        "n_clusters": int(len(use)),
        "n_features": int(X.shape[1]),
        "residual_sd": residual_sd,
        "r2": 1 - ssr / sst,
    }
    return pd.DataFrame(rows), diag


def fit_domain_analyses(clusters: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    levels_all = categorical_levels(clusters)
    multi = clusters[clusters["cluster_size"] >= 2].copy()
    levels_multi = categorical_levels(multi)

    outcome_results = []
    outcome_diag = []
    for domain in DOMAINS:
        for outcome_name, outcome in OUTCOME_SPECS.items():
            spec = DomainModelSpec(
                domain=domain,
                outcome_name=outcome_name,
                outcome=outcome["outcome"],
                label=outcome["label"],
                scale=outcome["scale"],
            )
            result, diag = fit_log_linear(clusters, spec, levels_all)
            outcome_results.append(result)
            outcome_diag.append(diag)

    mixing_results = []
    mixing_diag = []
    for domain in DOMAINS:
        spec = DomainModelSpec(
            domain=domain,
            outcome_name=f"{domain}_mixing",
            outcome=f"{domain}_excess_discordance",
            label=f"{DOMAINS[domain]['label']} quintile mixing",
            scale="identity",
            include_size=True,
        )
        result, diag = fit_log_linear(multi, spec, levels_multi)
        mixing_results.append(result)
        mixing_diag.append(diag)

    return (
        pd.concat(outcome_results, ignore_index=True),
        pd.DataFrame(outcome_diag),
        pd.concat(mixing_results, ignore_index=True),
        pd.DataFrame(mixing_diag),
    )


def plot_domain_outcomes(results: pd.DataFrame, out_base: Path) -> None:
    domain_rows = []
    for domain in DOMAINS:
        term = f"{domain}_deprivation_z"
        domain_rows.append(results[(results["domain"] == domain) & (results["term"] == term)])
    plot_df = pd.concat(domain_rows, ignore_index=True)

    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    models = list(OUTCOME_SPECS)
    colours = {
        "cluster_size": "#4e79a7",
        "duration": "#f28e2b",
        "geographic_dispersion": "#59a14f",
    }
    domains = list(DOMAINS)
    y_pos = {domain: i for i, domain in enumerate(domains)}
    offsets = {"cluster_size": -0.22, "duration": 0.0, "geographic_dispersion": 0.22}

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for _, row in plot_df.iterrows():
        y = y_pos[row["domain"]] + offsets[row["model"]]
        ax.plot([row["ci_low"], row["ci_high"]], [y, y], color=colours[row["model"]], linewidth=1.2)
        ax.scatter(row["geometric_mean_ratio"], y, color=colours[row["model"]], s=18, label=row["model_label"])

    ax.axvline(1, color="#666666", linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels([DOMAINS[d]["label"] for d in domains])
    ax.set_xlabel("Adjusted geometric mean ratio per 1 SD higher domain deprivation")
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=False)
    fig.subplots_adjust(bottom=0.25, left=0.2, right=0.98)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_domain_mixing(results: pd.DataFrame, out_base: Path) -> None:
    plot_df = []
    for domain in DOMAINS:
        term = f"{domain}_deprivation_z"
        plot_df.append(results[(results["domain"] == domain) & (results["term"] == term)])
    plot_df = pd.concat(plot_df, ignore_index=True)

    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    domains = list(DOMAINS)
    y_pos = {domain: i for i, domain in enumerate(domains)}
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for _, row in plot_df.iterrows():
        y = y_pos[row["domain"]]
        ax.plot(
            [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
            [y, y],
            color="#2b2b2b",
            linewidth=1.2,
        )
        ax.scatter(row["coefficient_percentage_points"], y, color="#2b2b2b", s=18)

    ax.axvline(0, color="#666666", linewidth=0.8, linestyle="--")
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels([DOMAINS[d]["label"] for d in domains])
    ax.set_xlabel("Change in excess domain-quintile mixing, percentage points per 1 SD higher domain deprivation")
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    fig.subplots_adjust(left=0.2, right=0.98)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def run(root: Path, qc: str | None) -> None:
    tables_dir = root / "part1" / "tables"
    figures_dir = root / "part1" / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("Reading sequence rows", flush=True)
    seq = read_sequence_rows(root, qc)
    print(f"Building domain cluster table from {len(seq):,} rows", flush=True)
    clusters, scaling = build_cluster_table(seq)

    print(f"Fitting domain models for {len(clusters):,} clusters", flush=True)
    outcome_results, outcome_diag, mixing_results, mixing_diag = fit_domain_analyses(clusters)

    outcome_results.to_csv(tables_dir / "simd_domain_outcome_model_results.csv", index=False)
    outcome_diag.to_csv(tables_dir / "simd_domain_outcome_model_diagnostics.csv", index=False)
    mixing_results.to_csv(tables_dir / "simd_domain_mixing_model_results.csv", index=False)
    mixing_diag.to_csv(tables_dir / "simd_domain_mixing_model_diagnostics.csv", index=False)
    scaling.to_csv(tables_dir / "simd_domain_covariate_scaling.csv", index=False)

    del seq, clusters
    gc.collect()

    plot_domain_outcomes(outcome_results, figures_dir / "simd_domain_outcome_effects")
    plot_domain_mixing(mixing_results, figures_dir / "simd_domain_mixing_effects")

    print(f"Wrote {tables_dir / 'simd_domain_outcome_model_results.csv'}", flush=True)
    print(f"Wrote {tables_dir / 'simd_domain_mixing_model_results.csv'}", flush=True)
    print(f"Wrote {figures_dir / 'simd_domain_outcome_effects.png'}", flush=True)
    print(f"Wrote {figures_dir / 'simd_domain_mixing_effects.png'}", flush=True)


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
