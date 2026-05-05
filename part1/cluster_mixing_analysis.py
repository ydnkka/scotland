"""Socioeconomic and demographic mixing within genomic clusters.

This analysis treats mixing as within-cluster pairwise discordance. For example,
SIMD mixing is the probability that two non-identical cases drawn from the same
cluster have different SIMD quintiles. Observed discordance is compared with the
expected discordance among all sampled cases from the same lineage, window, and
Leiden resolution; the model outcome is observed minus expected discordance.

Positive excess discordance means clusters are more mixed than expected given
the social/demographic composition of the sampled stratum. Negative values mean
clusters are more homogeneous, or assortative, than expected.
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

SEQUENCE_COLUMNS = [
    "cluster_id",
    "sequence_id",
    "resolution",
    "window_id",
    "window_idx",
    "wn_mid_date",
    "wn_prop_sequenced",
    "collection_date",
    "pango_lineage",
    "nextclade_qc",
    "age_band",
    "sex",
    "dz_simd_rank",
    "dz_simd_quintile",
    "dz_cum_incidence_per_capita",
    "dz_cum_prop_sequenced",
    "dz_7d_test_positivity",
]

MIXING_VARIABLES = {
    "simd": {
        "column": "dz_simd_quintile",
        "label": "SIMD quintile mixing",
        "short_label": "SIMD",
    },
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
    "profile": {
        "column": "socio_demographic_profile",
        "label": "Joint SIMD-age-sex profile mixing",
        "short_label": "Joint profile",
    },
}

PRIMARY_TERMS = [
    "deprivation_z",
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
    "log_cluster_size_z",
]

TERM_LABELS = {
    "deprivation_z": "Mean SIMD deprivation",
    "local_incidence_z": "Local cumulative incidence",
    "local_seq_fraction_z": "Local sequencing fraction",
    "window_seq_fraction_z": "Window sequencing proportion",
    "test_positivity_z": "Local test positivity",
    "log_cluster_size_z": "Cluster size",
}


@dataclass(frozen=True)
class MixingModel:
    name: str
    outcome: str
    label: str


def read_sequence_rows(path: Path, qc: str | None) -> pd.DataFrame:
    filters = None if qc is None else [("nextclade_qc", "==", qc)]
    df = pd.read_parquet(
        path,
        columns=SEQUENCE_COLUMNS,
        filters=filters,
        engine="pyarrow",
    )

    categorical = [
        "cluster_id",
        "sequence_id",
        "window_id",
        "pango_lineage",
        "nextclade_qc",
        "age_band",
        "sex",
        "dz_simd_quintile",
    ]
    for col in categorical:
        df[col] = df[col].astype("category")

    df["wn_mid_date"] = pd.to_datetime(df["wn_mid_date"])
    df["resolution_label"] = df["resolution"].map(lambda x: f"{x:.1f}")
    df["socio_demographic_profile"] = (
        df["dz_simd_quintile"].astype(str)
        + "|"
        + df["age_band"].astype(str)
        + "|"
        + df["sex"].astype(str)
    ).astype("category")
    return df


def pairwise_discordance_from_counts(counts: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """Return pairwise discordance from long counts by group and category."""
    totals = counts.groupby(group_cols, observed=True)["n"].sum().rename("n_valid")
    same_pairs = (
        counts.assign(same_pairs=counts["n"] * (counts["n"] - 1))
        .groupby(group_cols, observed=True)["same_pairs"]
        .sum()
    )
    out = pd.concat([totals, same_pairs], axis=1).reset_index()
    denom = out["n_valid"] * (out["n_valid"] - 1)
    out["discordance"] = np.where(
        denom > 0,
        1 - out["same_pairs"] / denom,
        np.nan,
    )
    return out.drop(columns=["same_pairs"])


def observed_cluster_discordance(
    seq: pd.DataFrame,
    variable: str,
    prefix: str,
) -> pd.DataFrame:
    counts = (
        seq.dropna(subset=[variable])
        .groupby(["cluster_id", variable], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = pairwise_discordance_from_counts(counts, ["cluster_id"])
    return out.rename(
        columns={
            "n_valid": f"{prefix}_n_valid",
            "discordance": f"{prefix}_discordance",
        }
    )


def expected_stratum_discordance(
    seq: pd.DataFrame,
    variable: str,
    prefix: str,
) -> pd.DataFrame:
    strata = ["window_id", "pango_lineage", "resolution_label"]
    counts = (
        seq.dropna(subset=[variable])
        .groupby(strata + [variable], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = pairwise_discordance_from_counts(counts, strata)
    return out.rename(
        columns={
            "n_valid": f"{prefix}_stratum_n_valid",
            "discordance": f"{prefix}_expected_discordance",
        }
    )


def build_cluster_table(seq: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = [
        "cluster_id",
        "sequence_id",
        "resolution",
        "window_id",
        "window_idx",
        "pango_lineage",
        "dz_simd_rank",
        "dz_cum_incidence_per_capita",
        "dz_cum_prop_sequenced",
        "wn_prop_sequenced",
        "dz_7d_test_positivity",
    ]
    before = len(seq)
    seq = seq.dropna(subset=required).copy()
    dropped = before - len(seq)

    clusters = (
        seq.groupby("cluster_id", observed=True, sort=False)
        .agg(
            cluster_size=("sequence_id", "nunique"),
            resolution=("resolution", "first"),
            resolution_label=("resolution_label", "first"),
            window_id=("window_id", "first"),
            window_idx=("window_idx", "first"),
            wn_mid_date=("wn_mid_date", "first"),
            pango_lineage=("pango_lineage", "first"),
            mean_simd_rank=("dz_simd_rank", "mean"),
            mean_local_incidence_per_capita=("dz_cum_incidence_per_capita", "mean"),
            mean_local_seq_fraction=("dz_cum_prop_sequenced", "mean"),
            mean_window_seq_fraction=("wn_prop_sequenced", "mean"),
            mean_test_positivity=("dz_7d_test_positivity", "mean"),
        )
        .reset_index()
    )

    for prefix, spec in MIXING_VARIABLES.items():
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

    clusters["deprivation_raw"] = -clusters["mean_simd_rank"]
    clusters["local_incidence_log"] = np.log1p(
        clusters["mean_local_incidence_per_capita"].clip(lower=0) * 1000
    )
    clusters["local_seq_fraction_logit"] = logit_clipped(clusters["mean_local_seq_fraction"])
    clusters["window_seq_fraction_logit"] = logit_clipped(clusters["mean_window_seq_fraction"])
    clusters["test_positivity_logit"] = logit_clipped(clusters["mean_test_positivity"].fillna(0))
    clusters["log_cluster_size"] = np.log(clusters["cluster_size"])

    scaling_rows = []
    transforms = {
        "deprivation_z": "deprivation_raw",
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

    scaling = pd.DataFrame(scaling_rows)
    scaling.attrs["dropped_sequence_rows_missing_model_fields"] = dropped
    return clusters, scaling


def fit_ols(
    df: pd.DataFrame,
    model: MixingModel,
    levels: dict[str, list],
) -> tuple[pd.DataFrame, dict]:
    use = df.dropna(subset=[model.outcome, *PRIMARY_TERMS]).copy()
    X, feature_names = build_design_matrix(use, PRIMARY_TERMS, levels)
    y = use[model.outcome].to_numpy(dtype=np.float64)

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
    feature_to_idx = {name: i for i, name in enumerate(feature_names)}

    rows = []
    for term in PRIMARY_TERMS:
        idx = feature_to_idx[term]
        coef = float(beta[idx])
        stderr = float(se[idx])
        z = coef / stderr if stderr > 0 else np.nan
        rows.append(
            {
                "model": model.name,
                "model_label": model.label,
                "outcome": model.outcome,
                "term": term,
                "term_label": TERM_LABELS[term],
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
        "model": model.name,
        "model_label": model.label,
        "outcome": model.outcome,
        "n_clusters": int(len(use)),
        "n_features": int(X.shape[1]),
        "n_lineages": int(len(levels["pango_lineage"])),
        "n_windows": int(len(levels["window_id"])),
        "n_resolutions": int(len(levels["resolution_label"])),
        "mean_outcome": float(np.mean(y)),
        "sd_outcome": float(np.std(y, ddof=0)),
        "residual_sd": residual_sd,
        "r2": 1 - ssr / sst,
    }
    return pd.DataFrame(rows), diagnostics


def summarise_mixing(clusters: pd.DataFrame, qc: str | None, dropped: int) -> pd.DataFrame:
    rows = [
        {"measure": "clusters_total", "value": len(clusters)},
        {"measure": "clusters_size_ge_2", "value": int((clusters["cluster_size"] >= 2).sum())},
        {"measure": "sequence_rows_dropped_missing_model_fields", "value": dropped},
        {"measure": "windows", "value": clusters["window_id"].nunique()},
        {"measure": "pango_lineages", "value": clusters["pango_lineage"].nunique()},
        {"measure": "leiden_resolutions", "value": clusters["resolution"].nunique()},
        {"measure": "qc_filter", "value": qc or "none"},
    ]

    for prefix, spec in MIXING_VARIABLES.items():
        for col in [
            f"{prefix}_discordance",
            f"{prefix}_expected_discordance",
            f"{prefix}_excess_discordance",
        ]:
            values = clusters.loc[clusters["cluster_size"] >= 2, col].dropna()
            desc = values.describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9])
            rows.extend(
                {
                    "measure": f"{spec['short_label']} {col}",
                    "statistic": key,
                    "value": value,
                }
                for key, value in desc.items()
            )
    return pd.DataFrame(rows)


def plot_mixing_effects(results: pd.DataFrame, out_base: Path) -> None:
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
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 8,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    models = ["simd", "age", "sex", "profile"]
    model_positions = {model: i for i, model in enumerate(models)}
    term_offsets = np.linspace(-0.32, 0.32, len(PRIMARY_TERMS))
    term_positions = dict(zip(PRIMARY_TERMS, term_offsets))
    colours = {
        "deprivation_z": "#2b2b2b",
        "local_incidence_z": "#4e79a7",
        "local_seq_fraction_z": "#59a14f",
        "window_seq_fraction_z": "#f28e2b",
        "test_positivity_z": "#b07aa1",
        "log_cluster_size_z": "#7f7f7f",
    }

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for _, row in results.iterrows():
        y = model_positions[row["model"]] + term_positions[row["term"]]
        ax.plot(
            [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
            [y, y],
            color=colours[row["term"]],
            linewidth=1.2,
            solid_capstyle="round",
        )
        ax.scatter(
            row["coefficient_percentage_points"],
            y,
            color=colours[row["term"]],
            s=18,
            zorder=3,
            label=TERM_LABELS[row["term"]],
        )

    ax.axvline(0, color="#666666", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Change in excess pairwise discordance, percentage points per 1 SD higher covariate")
    ax.set_yticks(list(model_positions.values()))
    ax.set_yticklabels([MIXING_VARIABLES[m]["short_label"] for m in models])
    ax.set_ylim(-0.6, len(models) - 0.4)
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
        bbox_to_anchor=(0.5, -0.22),
        ncol=2,
        frameon=False,
        columnspacing=1.4,
        handlelength=1.2,
    )

    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(bottom=0.34, left=0.2, right=0.98)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def run(root: Path, qc: str | None) -> None:
    tables_dir = root / "part1" / "tables"
    figures_dir = root / "part1" / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    path = analysis_dataset_path(root)
    print(f"Reading sequence rows from {path}", flush=True)
    seq = read_sequence_rows(path, qc)
    print(f"Building mixing table from {len(seq):,} sequence-window-resolution rows", flush=True)
    clusters, scaling = build_cluster_table(seq)
    dropped = int(scaling.attrs.get("dropped_sequence_rows_missing_model_fields", 0))
    multi = clusters[clusters["cluster_size"] >= 2].copy()
    levels = categorical_levels(multi)

    models = [
        MixingModel(prefix, f"{prefix}_excess_discordance", spec["label"])
        for prefix, spec in MIXING_VARIABLES.items()
    ]

    print(
        f"Fitting mixing models for {len(multi):,} non-singleton clusters, "
        f"{len(levels['pango_lineage'])} lineages, {len(levels['window_id'])} windows, "
        f"{len(levels['resolution_label'])} resolutions",
        flush=True,
    )
    result_frames = []
    diagnostics = []
    for model in models:
        print(f"  - {model.name}", flush=True)
        result, diag = fit_ols(multi, model, levels)
        result_frames.append(result)
        diagnostics.append(diag)

    results = pd.concat(result_frames, ignore_index=True)
    diagnostics_df = pd.DataFrame(diagnostics)
    descriptives = summarise_mixing(clusters, qc, dropped)

    results.to_csv(tables_dir / "cluster_mixing_model_results.csv", index=False)
    diagnostics_df.to_csv(tables_dir / "cluster_mixing_model_diagnostics.csv", index=False)
    descriptives.to_csv(tables_dir / "cluster_mixing_descriptives.csv", index=False)
    scaling.to_csv(tables_dir / "cluster_mixing_covariate_scaling.csv", index=False)

    del seq, clusters, multi, levels, result_frames
    gc.collect()

    plot_mixing_effects(results, figures_dir / "cluster_mixing_model_effects")

    print(f"Wrote {tables_dir / 'cluster_mixing_model_results.csv'}", flush=True)
    print(f"Wrote {tables_dir / 'cluster_mixing_model_diagnostics.csv'}", flush=True)
    print(f"Wrote {figures_dir / 'cluster_mixing_model_effects.png'}", flush=True)


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
    run(args.root.resolve(), qc=qc)


if __name__ == "__main__":
    main()
