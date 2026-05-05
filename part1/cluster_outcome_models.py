"""Adjusted models for SARS-CoV-2 genomic cluster outcomes in Scotland.

This script answers the Part 1 question:

    After accounting for lineage, calendar time, Leiden resolution, local
    incidence, and sequencing intensity, are socioeconomic deprivation and local
    surveillance conditions associated with larger, longer-lasting, or more
    geographically dispersed genomic clusters?

The unit of analysis is an inferred cluster within a sliding window and Leiden
resolution. Member-level SIMD and surveillance variables are averaged within
each cluster before modelling. Outcomes are modelled on the log scale, so
exponentiated coefficients are adjusted geometric mean ratios.
"""

from __future__ import annotations

import argparse
import gc
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.linalg import pinvh
from scipy.stats import norm


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
    "dz_simd_rank",
    "dz_cum_incidence_per_capita",
    "dz_cum_prop_sequenced",
    "dz_7d_test_positivity",
]

PRIMARY_TERMS = [
    "deprivation_z",
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
]

TERM_LABELS = {
    "deprivation_z": "Mean SIMD deprivation",
    "local_incidence_z": "Local cumulative incidence",
    "local_seq_fraction_z": "Local sequencing fraction",
    "window_seq_fraction_z": "Window sequencing proportion",
    "test_positivity_z": "Local test positivity",
    "log_cluster_size_z": "Cluster size",
}

MODEL_LABELS = {
    "cluster_size": "Cluster size",
    "duration": "Duration",
    "geographic_dispersion": "Geographic dispersion",
    "duration_size_adjusted": "Duration, size-adjusted",
    "geographic_dispersion_size_adjusted": "Geographic dispersion, size-adjusted",
}


@dataclass(frozen=True)
class ModelSpec:
    name: str
    outcome: str
    include_size: bool = False


MODEL_SPECS = [
    ModelSpec("cluster_size", "cluster_size"),
    ModelSpec("duration", "duration_days_plus1"),
    ModelSpec("geographic_dispersion", "cluster_n_datazones"),
    ModelSpec("duration_size_adjusted", "duration_days_plus1", include_size=True),
    ModelSpec("geographic_dispersion_size_adjusted", "cluster_n_datazones", include_size=True),
]


def repo_root(start: Path | None = None) -> Path:
    """Walk upward until config.yaml is found."""
    p = (start or Path(__file__)).resolve()
    for candidate in [p, *p.parents]:
        if (candidate / "config.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate config.yaml.")


def analysis_dataset_path(root: Path) -> Path:
    with open(root / "config.yaml", "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return root / cfg["data"]["processed"]["analysis_dataset"]


def zscore(values: pd.Series) -> tuple[pd.Series, float, float]:
    mean = float(values.mean())
    sd = float(values.std(ddof=0))
    if not math.isfinite(sd) or sd == 0:
        raise ValueError(f"Cannot standardise {values.name!r}: zero or invalid SD.")
    return (values - mean) / sd, mean, sd


def logit_clipped(values: pd.Series, eps: float = 1e-5) -> pd.Series:
    clipped = values.clip(lower=eps, upper=1 - eps)
    return np.log(clipped / (1 - clipped))


def read_sequence_rows(path: Path, qc: str | None) -> pd.DataFrame:
    """Read only columns needed for the cluster-level analysis."""
    filters = None if qc is None else [("nextclade_qc", "==", qc)]
    df = pd.read_parquet(
        path,
        columns=SEQUENCE_COLUMNS,
        filters=filters,
        engine="pyarrow",
    )

    for col in ["cluster_id", "sequence_id", "window_id", "datazone", "pango_lineage", "nextclade_qc"]:
        df[col] = df[col].astype("category")

    df["collection_date"] = pd.to_datetime(df["collection_date"])
    df["wn_mid_date"] = pd.to_datetime(df["wn_mid_date"])
    return df


def build_cluster_dataset(seq: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse sequence rows to one row per inferred cluster."""
    required = [
        "cluster_id",
        "sequence_id",
        "resolution",
        "window_id",
        "window_idx",
        "collection_date",
        "datazone",
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

    grouped = seq.groupby("cluster_id", observed=True, sort=False)
    clusters = grouped.agg(
        cluster_size=("sequence_id", "nunique"),
        cluster_n_datazones=("datazone", "nunique"),
        cluster_start_date=("collection_date", "min"),
        cluster_end_date=("collection_date", "max"),
        resolution=("resolution", "first"),
        window_id=("window_id", "first"),
        window_idx=("window_idx", "first"),
        wn_mid_date=("wn_mid_date", "first"),
        pango_lineage=("pango_lineage", "first"),
        mean_simd_rank=("dz_simd_rank", "mean"),
        mean_local_incidence_per_capita=("dz_cum_incidence_per_capita", "mean"),
        mean_local_seq_fraction=("dz_cum_prop_sequenced", "mean"),
        mean_window_seq_fraction=("wn_prop_sequenced", "mean"),
        mean_test_positivity=("dz_7d_test_positivity", "mean"),
    ).reset_index()

    clusters["duration_days"] = (
        clusters["cluster_end_date"] - clusters["cluster_start_date"]
    ).dt.days.astype(int)
    clusters["duration_days_plus1"] = clusters["duration_days"] + 1
    clusters["resolution_label"] = clusters["resolution"].map(lambda x: f"{x:.1f}")

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


def categorical_levels(df: pd.DataFrame) -> dict[str, list]:
    lineage_levels = (
        df["pango_lineage"]
        .astype(str)
        .value_counts()
        .sort_values(ascending=False)
        .index.tolist()
    )
    window_levels = (
        df.sort_values(["window_idx", "window_id"])["window_id"]
        .astype(str)
        .drop_duplicates()
        .tolist()
    )
    resolution_levels = (
        df.sort_values("resolution")["resolution_label"]
        .astype(str)
        .drop_duplicates()
        .tolist()
    )
    return {
        "pango_lineage": lineage_levels,
        "window_id": window_levels,
        "resolution_label": resolution_levels,
    }


def build_design_matrix(
    df: pd.DataFrame,
    numeric_terms: list[str],
    levels: dict[str, list],
) -> tuple[sparse.csr_matrix, list[str]]:
    """Build a sparse design matrix with dropped-baseline categorical effects."""
    n = len(df)
    blocks: list[sparse.csr_matrix] = []
    names: list[str] = []

    intercept = sparse.csr_matrix(np.ones((n, 1), dtype=np.float64))
    blocks.append(intercept)
    names.append("Intercept")

    numeric = sparse.csr_matrix(df[numeric_terms].to_numpy(dtype=np.float64))
    blocks.append(numeric)
    names.extend(numeric_terms)

    for col, cats in levels.items():
        cat = pd.Categorical(df[col].astype(str), categories=cats)
        codes = cat.codes
        if np.any(codes < 0):
            missing = sorted(set(df.loc[codes < 0, col].astype(str)))
            raise ValueError(f"{col} has values outside configured levels: {missing[:5]}")

        keep = codes > 0
        row = np.flatnonzero(keep)
        col_idx = codes[keep] - 1
        mat = sparse.csr_matrix(
            (np.ones(len(row), dtype=np.float64), (row, col_idx)),
            shape=(n, len(cats) - 1),
        )
        blocks.append(mat)
        names.extend([f"{col}={level}" for level in cats[1:]])

    return sparse.hstack(blocks, format="csr"), names


def clustered_ols_covariance(
    X: sparse.csr_matrix,
    residuals: np.ndarray,
    groups: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Window-clustered sandwich covariance for a log-linear fixed-effect model."""
    n, p = X.shape
    xtx = (X.T @ X).toarray()
    xtx_inv = pinvh(xtx, rtol=1e-10)

    group_codes, inverse = np.unique(groups, return_inverse=True)
    scores = np.zeros((len(group_codes), p), dtype=np.float64)
    for group_idx in range(len(group_codes)):
        idx = np.flatnonzero(inverse == group_idx)
        scores[group_idx, :] = X[idx, :].T @ residuals[idx]

    meat = scores.T @ scores
    correction = (len(group_codes) / (len(group_codes) - 1)) * ((n - 1) / (n - p))
    cov = correction * (xtx_inv @ meat @ xtx_inv)
    residual_sd = float(np.sqrt(np.sum(residuals**2) / (n - p)))
    return cov, residual_sd


def fit_one_model(
    df: pd.DataFrame,
    spec: ModelSpec,
    levels: dict[str, list],
) -> tuple[pd.DataFrame, dict]:
    numeric_terms = PRIMARY_TERMS.copy()
    if spec.include_size:
        numeric_terms.append("log_cluster_size_z")

    X, feature_names = build_design_matrix(df, numeric_terms, levels)
    y_raw = df[spec.outcome].to_numpy(dtype=np.float64)
    y = np.log(y_raw)

    xtx = (X.T @ X).toarray()
    xty = X.T @ y
    beta = pinvh(xtx, rtol=1e-10) @ xty
    fitted = X @ beta
    residuals = y - fitted

    cov, residual_sd = clustered_ols_covariance(
        X,
        residuals,
        groups=df["window_id"].astype(str).to_numpy(),
    )
    se = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
    feature_to_idx = {name: i for i, name in enumerate(feature_names)}

    rows = []
    for term in numeric_terms:
        idx = feature_to_idx[term]
        coef = float(beta[idx])
        stderr = float(se[idx])
        z = coef / stderr if stderr > 0 else np.nan
        rows.append(
            {
                "model": spec.name,
                "model_label": MODEL_LABELS[spec.name],
                "outcome": spec.outcome,
                "term": term,
                "term_label": TERM_LABELS[term],
                "coefficient_log_ratio": coef,
                "std_error_clustered_by_window": stderr,
                "z": z,
                "p_value": float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan,
                "geometric_mean_ratio": float(np.exp(coef)),
                "ci_low": float(np.exp(coef - 1.96 * stderr)),
                "ci_high": float(np.exp(coef + 1.96 * stderr)),
            }
        )

    diagnostics = {
        "model": spec.name,
        "model_label": MODEL_LABELS[spec.name],
        "outcome": spec.outcome,
        "n_clusters": int(len(df)),
        "n_features": int(X.shape[1]),
        "n_lineages": int(len(levels["pango_lineage"])),
        "n_windows": int(len(levels["window_id"])),
        "n_resolutions": int(len(levels["resolution_label"])),
        "mean_outcome": float(np.mean(y_raw)),
        "max_outcome": float(np.max(y_raw)),
        "residual_sd_log_scale": residual_sd,
    }
    ssr = float(np.sum(residuals**2))
    sst = float(np.sum((y - y.mean()) ** 2))
    diagnostics["r2_log_scale"] = 1 - ssr / sst
    diagnostics["adjusted_r2_log_scale"] = 1 - (ssr / (len(y) - X.shape[1])) / (
        sst / (len(y) - 1)
    )
    return pd.DataFrame(rows), diagnostics


def describe_clusters(
    seq: pd.DataFrame,
    clusters: pd.DataFrame,
    qc: str | None,
    dropped_missing: int,
) -> pd.DataFrame:
    outcome_summary = clusters[
        ["cluster_size", "duration_days", "cluster_n_datazones"]
    ].describe(percentiles=[0.5, 0.9, 0.95, 0.99]).T.reset_index(names="measure")

    overview = pd.DataFrame(
        [
            {"measure": "sequence_rows_used", "count": len(seq)},
            {"measure": "sequence_rows_dropped_missing_model_fields", "count": dropped_missing},
            {"measure": "clusters", "count": len(clusters)},
            {"measure": "windows", "count": clusters["window_id"].nunique()},
            {"measure": "pango_lineages", "count": clusters["pango_lineage"].nunique()},
            {"measure": "leiden_resolutions", "count": clusters["resolution"].nunique()},
            {
                "measure": "singleton_cluster_fraction",
                "count": float((clusters["cluster_size"] == 1).mean()),
            },
            {"measure": "qc_filter", "count": qc or "none"},
        ]
    )
    return pd.concat([overview, outcome_summary], ignore_index=True, sort=False)


def plot_effects(results: pd.DataFrame, out_base: Path) -> None:
    """Plot primary terms for the three unadjusted-for-size outcome models."""
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

    primary_models = ["cluster_size", "duration", "geographic_dispersion"]
    terms = PRIMARY_TERMS
    plot_df = results[
        results["model"].isin(primary_models) & results["term"].isin(terms)
    ].copy()

    model_positions = {model: i for i, model in enumerate(primary_models)}
    term_offsets = np.linspace(-0.28, 0.28, len(terms))
    term_positions = dict(zip(terms, term_offsets))
    colours = {
        "deprivation_z": "#2b2b2b",
        "local_incidence_z": "#4e79a7",
        "local_seq_fraction_z": "#59a14f",
        "window_seq_fraction_z": "#f28e2b",
        "test_positivity_z": "#b07aa1",
    }

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for _, row in plot_df.iterrows():
        y = model_positions[row["model"]] + term_positions[row["term"]]
        ax.plot(
            [row["ci_low"], row["ci_high"]],
            [y, y],
            color=colours[row["term"]],
            linewidth=1.4,
            solid_capstyle="round",
        )
        ax.scatter(
            row["geometric_mean_ratio"],
            y,
            color=colours[row["term"]],
            s=22,
            zorder=3,
            label=TERM_LABELS[row["term"]],
        )

    ax.axvline(1, color="#666666", linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlabel("Adjusted geometric mean ratio per 1 SD higher cluster-level covariate")
    ax.set_yticks(list(model_positions.values()))
    ax.set_yticklabels([MODEL_LABELS[m] for m in primary_models])
    ax.set_ylim(-0.55, len(primary_models) - 0.45)
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
        bbox_to_anchor=(0.5, -0.24),
        ncol=2,
        frameon=False,
        columnspacing=1.4,
        handlelength=1.2,
    )

    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(bottom=0.34, left=0.23, right=0.98)
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
    print(
        f"Building cluster-level dataset from {len(seq):,} sequence-window-resolution rows",
        flush=True,
    )
    clusters, scaling = build_cluster_dataset(seq)
    dropped_missing = int(scaling.attrs.get("dropped_sequence_rows_missing_model_fields", 0))

    levels = categorical_levels(clusters)
    print(
        "Fitting models for "
        f"{len(clusters):,} clusters, {len(levels['pango_lineage'])} lineages, "
        f"{len(levels['window_id'])} windows, {len(levels['resolution_label'])} resolutions",
        flush=True,
    )

    all_results = []
    diagnostics = []
    for spec in MODEL_SPECS:
        print(f"  - {spec.name}", flush=True)
        result, diag = fit_one_model(clusters, spec, levels)
        all_results.append(result)
        diagnostics.append(diag)

    results = pd.concat(all_results, ignore_index=True)
    diagnostics_df = pd.DataFrame(diagnostics)
    descriptives = describe_clusters(seq, clusters, qc, dropped_missing)

    results.to_csv(tables_dir / "cluster_outcome_model_results.csv", index=False)
    diagnostics_df.to_csv(tables_dir / "cluster_outcome_model_diagnostics.csv", index=False)
    descriptives.to_csv(tables_dir / "cluster_outcome_descriptives.csv", index=False)
    scaling.to_csv(tables_dir / "cluster_outcome_covariate_scaling.csv", index=False)

    del seq, clusters, levels, all_results
    gc.collect()

    plot_effects(results, figures_dir / "cluster_outcome_model_effects")

    print(f"Wrote {tables_dir / 'cluster_outcome_model_results.csv'}", flush=True)
    print(f"Wrote {tables_dir / 'cluster_outcome_model_diagnostics.csv'}", flush=True)
    print(f"Wrote {figures_dir / 'cluster_outcome_model_effects.png'}", flush=True)


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
