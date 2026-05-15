from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .data import FeatureSpec


ENTROPY_FEATURES = [
    "age_entropy_norm",
    "sex_entropy_norm",
    "simd_entropy_norm",
    "demographic_entropy_norm",
    "sociodemographic_entropy_norm",
]

MIXING_ADJUSTMENT_NUMERIC_FEATURES = [
    "cluster_mid_days_since_epoch",
    "log_wn_no_sequences",
    "wn_prop_sequenced",
    "dz_cum_prop_sequenced_mean",
    "log_dz_cum_incidence_per_capita_mean",
    "dz_7d_test_positivity_mean",
    "log_dz_population_density_mean",
]

MIXING_ADJUSTMENT_CATEGORICAL_FEATURES = [
    "who_voc_model",
]


@dataclass(slots=True)
class MixingFeatureSet:
    name: str
    numeric_features: list[str]
    categorical_features: list[str]

    @property
    def all_features(self) -> list[str]:
        return self.numeric_features + self.categorical_features


def add_cluster_uid(df: pd.DataFrame) -> pd.DataFrame:
    """Add a stable cluster key for the resolution/window/cluster combination."""
    out = df.copy()
    key_cols = [c for c in ["resolution", "window_id", "cluster_id"] if c in out.columns]
    if not key_cols:
        raise ValueError("At least one cluster identifier column is required.")
    key = _clean_category(out[key_cols[0]]).fillna("Missing")
    for col in key_cols[1:]:
        key = key.str.cat(_clean_category(out[col]).fillna("Missing"), sep=" | ")
    out["cluster_uid"] = key
    return out


def _clean_category(series: pd.Series) -> pd.Series:
    clean = series.astype("string").str.strip()
    return clean.mask(clean.isin(["", "nan", "NaN", "<NA>", "None"]))


def _make_joint_category(df: pd.DataFrame, columns: list[str], name: str) -> pd.Series:
    clean = pd.DataFrame({col: _clean_category(df[col]) for col in columns}, index=df.index)
    complete = clean.notna().all(axis=1)
    out = pd.Series(pd.NA, index=df.index, dtype="string", name=name)
    if complete.any():
        joined = clean.loc[complete, columns[0]]
        for col in columns[1:]:
            joined = joined.str.cat(clean.loc[complete, col], sep=" x ")
        out.loc[complete] = joined
    return out


def _category_count(series: pd.Series) -> int:
    return int(_clean_category(series).dropna().nunique())


def entropy_by_cluster(
    df: pd.DataFrame,
    value_col: str,
    prefix: str,
    cluster_col: str = "cluster_uid",
    n_possible_categories: int | None = None,
    finite_sample_normalisation: bool = True,
) -> pd.DataFrame:
    """Compute raw and normalised Shannon entropy for one categorical variable.

    Normalisation uses log(min(K, n_valid)) by default, where K is the number of
    possible categories in the analysis set and n_valid is the cluster's valid
    count for that variable. This keeps the maximum at 1 even for small clusters
    that cannot possibly contain every global category.
    """
    all_clusters = df[[cluster_col]].drop_duplicates().copy()
    valid = df[[cluster_col, value_col]].copy()
    valid[value_col] = _clean_category(valid[value_col])
    valid = valid.dropna(subset=[value_col])

    if valid.empty:
        return all_clusters.assign(
            **{
                f"{prefix}_entropy": np.nan,
                f"{prefix}_entropy_norm": np.nan,
                f"{prefix}_n_valid": 0,
                f"{prefix}_n_categories": 0,
            }
        )

    if n_possible_categories is None:
        n_possible_categories = int(valid[value_col].nunique())

    counts = valid.groupby([cluster_col, value_col], observed=True).size().rename("n")
    totals = counts.groupby(level=0).sum().rename(f"{prefix}_n_valid")
    proportions = counts / counts.groupby(level=0).transform("sum")
    entropy = (-(proportions * np.log(proportions))).groupby(level=0).sum()
    entropy = entropy.rename(f"{prefix}_entropy")
    n_categories = counts.groupby(level=0).size().rename(f"{prefix}_n_categories")

    out = pd.concat([entropy, totals, n_categories], axis=1).reset_index()
    denominator_categories = np.full(len(out), n_possible_categories, dtype=float)
    if finite_sample_normalisation:
        denominator_categories = np.minimum(denominator_categories, out[f"{prefix}_n_valid"])
    denominator = np.log(denominator_categories)

    out[f"{prefix}_entropy_norm"] = np.nan
    positive_denom = denominator > 0
    out.loc[positive_denom, f"{prefix}_entropy_norm"] = (
        out.loc[positive_denom, f"{prefix}_entropy"] / denominator[positive_denom]
    )
    out.loc[
        out[f"{prefix}_n_valid"].gt(0) & ~positive_denom,
        f"{prefix}_entropy_norm",
    ] = 0.0

    out = all_clusters.merge(out, on=cluster_col, how="left")
    out[f"{prefix}_n_valid"] = out[f"{prefix}_n_valid"].fillna(0).astype(int)
    out[f"{prefix}_n_categories"] = out[f"{prefix}_n_categories"].fillna(0).astype(int)
    return out


def _mode_or_missing(series: pd.Series) -> str:
    mode = _clean_category(series).dropna().mode()
    if mode.empty:
        return "Missing"
    return str(mode.iloc[0])


def make_cluster_mixing_table(
    model_df: pd.DataFrame,
    feature_spec: FeatureSpec,
    finite_sample_normalisation: bool = True,
) -> pd.DataFrame:
    """Build one row per non-singleton cluster with composition entropy features."""
    seq = add_cluster_uid(model_df)
    simd_feature = feature_spec.simd_overall_feature

    seq["age_entropy_category"] = _clean_category(seq["age_band"])
    seq["sex_entropy_category"] = _clean_category(seq["sex"])
    seq["simd_entropy_category"] = _clean_category(seq[simd_feature])
    seq["demographic_entropy_category"] = _make_joint_category(
        seq, ["age_entropy_category", "sex_entropy_category"], "demographic_entropy_category"
    )
    seq["sociodemographic_entropy_category"] = _make_joint_category(
        seq,
        ["age_entropy_category", "sex_entropy_category", "simd_entropy_category"],
        "sociodemographic_entropy_category",
    )

    grouped = seq.groupby("cluster_uid", observed=True)
    cluster = grouped.agg(
        cluster_id=("cluster_id", "first"),
        window_id=("window_id", "first"),
        window_idx=("window_idx", "first"),
        resolution=("resolution", "first"),
        pango_lineage=("pango_lineage", _mode_or_missing),
        who_voc_model=("who_voc_model", _mode_or_missing),
        cluster_type=("cluster_type", "first"),
        is_large_cluster=("is_large_cluster", "max"),
        reported_cluster_size=("cluster_size", "max"),
        observed_cluster_size=("sequence_id", "nunique"),
        n_patients=("patient_id", "nunique"),
        cluster_start_date=("collection_date", "min"),
        cluster_end_date=("collection_date", "max"),
        cluster_mid_days_since_epoch=("days_since_epoch", "median"),
        log_wn_no_sequences=("log_wn_no_sequences", "median"),
        wn_prop_sequenced=("wn_prop_sequenced", "median"),
        dz_cum_prop_sequenced_mean=("dz_cum_prop_sequenced", "mean"),
        log_dz_cum_incidence_per_capita_mean=("log_dz_cum_incidence_per_capita", "mean"),
        dz_7d_test_positivity_mean=("dz_7d_test_positivity", "mean"),
        log_dz_population_density_mean=("log_dz_population_density", "mean"),
    ).reset_index()

    entropy_specs = [
        ("age_entropy_category", "age"),
        ("sex_entropy_category", "sex"),
        ("simd_entropy_category", "simd"),
        ("demographic_entropy_category", "demographic"),
        ("sociodemographic_entropy_category", "sociodemographic"),
    ]
    for value_col, prefix in entropy_specs:
        cluster = cluster.merge(
            entropy_by_cluster(
                seq,
                value_col=value_col,
                prefix=prefix,
                n_possible_categories=_category_count(seq[value_col]),
                finite_sample_normalisation=finite_sample_normalisation,
            ),
            on="cluster_uid",
            how="left",
        )

    cluster["window_lineage_group"] = (
        cluster["window_id"].astype("string") + " | " + cluster["pango_lineage"].astype("string")
    )
    cluster["observed_fraction_reported"] = (
        cluster["observed_cluster_size"] / cluster["reported_cluster_size"].replace(0, np.nan)
    )
    cluster = cluster.loc[cluster["cluster_type"].isin(["small_cluster", "large_cluster"])].copy()
    cluster["is_large_cluster"] = cluster["is_large_cluster"].astype("int8")
    return cluster.reset_index(drop=True)


def make_mixing_feature_sets(adjusted: bool = True) -> dict[str, MixingFeatureSet]:
    """Return entropy-only and context-adjusted feature sets."""
    feature_sets = {
        "entropy_only": MixingFeatureSet(
            name="entropy_only",
            numeric_features=list(ENTROPY_FEATURES),
            categorical_features=[],
        )
    }
    if adjusted:
        feature_sets["context_adjusted"] = MixingFeatureSet(
            name="context_adjusted",
            numeric_features=list(ENTROPY_FEATURES) + list(MIXING_ADJUSTMENT_NUMERIC_FEATURES),
            categorical_features=list(MIXING_ADJUSTMENT_CATEGORICAL_FEATURES),
        )
    return feature_sets


def entropy_summary_by_cluster_type(cluster_table: pd.DataFrame) -> pd.DataFrame:
    """Summarise entropy distributions by small versus large clusters."""
    summary = (
        cluster_table.groupby("cluster_type", observed=True)[ENTROPY_FEATURES]
        .agg(["count", "mean", "std", "median", "min", "max"])
        .reset_index()
    )
    summary.columns = [
        "_".join(str(part) for part in col if part).strip("_")
        if isinstance(col, tuple)
        else col
        for col in summary.columns
    ]
    return summary


def plot_entropy_by_cluster_type(cluster_table: pd.DataFrame):
    """Boxplots of normalised entropy features by cluster type."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    long = cluster_table.melt(
        id_vars=["cluster_type"],
        value_vars=ENTROPY_FEATURES,
        var_name="entropy_feature",
        value_name="normalised_entropy",
    )
    fig, ax = plt.subplots(figsize=(11, 5))
    sns.boxplot(
        data=long,
        x="entropy_feature",
        y="normalised_entropy",
        hue="cluster_type",
        showfliers=False,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Normalised entropy")
    ax.set_title("Cluster composition mixing by cluster type")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(title="")
    fig.tight_layout()
    return fig


def plot_entropy_correlation(cluster_table: pd.DataFrame):
    """Correlation heatmap for the entropy predictors."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(7, 6))
    corr = cluster_table[ENTROPY_FEATURES].corr(method="spearman")
    sns.heatmap(corr, vmin=-1, vmax=1, cmap="vlag", annot=True, fmt=".2f", square=True, ax=ax)
    ax.set_title("Spearman correlation between entropy features")
    fig.tight_layout()
    return fig


def run_mixing_entropy_workflow(
    label: str,
    model_df: pd.DataFrame,
    feature_spec: FeatureSpec,
    config,
    out_dir: str | Path,
    membership_definition: str,
    save_models: bool = True,
) -> dict:
    """Run and save the cluster mixing entropy workflow for one membership definition."""
    import matplotlib.pyplot as plt

    from .models import (
        fit_binary_model,
        make_split_data,
        permutation_importance_table,
        sample_for_importance,
    )
    from .outputs import save_binary_result, save_run_metadata
    from .plots import (
        accumulated_local_effects_curves,
        partial_dependence_curves,
        plot_ale_curves,
        plot_binary_diagnostics,
        plot_importance,
        plot_partial_dependence_curves,
    )

    out = Path(out_dir)
    fig_dir = out / "figures"
    data_dir = out / "data"
    fig_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    cluster_table = make_cluster_mixing_table(
        model_df,
        feature_spec,
        finite_sample_normalisation=True,
    )
    cluster_table.to_parquet(data_dir / f"{label}_cluster_mixing_entropy_table.parquet", index=False)
    cluster_table.to_csv(data_dir / f"{label}_cluster_mixing_entropy_table.csv", index=False)

    cluster_counts = (
        cluster_table["cluster_type"]
        .value_counts(sort=False)
        .rename_axis("cluster_type")
        .reset_index(name="n_clusters")
    )
    cluster_counts = cluster_counts.loc[cluster_counts["n_clusters"].gt(0)].copy()
    cluster_counts["prop_clusters"] = cluster_counts["n_clusters"] / cluster_counts["n_clusters"].sum()
    cluster_counts.to_csv(out / f"{label}_outcome_counts.csv", index=False)

    entropy_summary = entropy_summary_by_cluster_type(cluster_table)
    entropy_summary.to_csv(out / f"{label}_entropy_summary_by_type.csv", index=False)

    fig = plot_entropy_by_cluster_type(cluster_table)
    fig.savefig(fig_dir / f"{label}_entropy_by_cluster_type.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / f"{label}_entropy_by_cluster_type.pdf", bbox_inches="tight")
    plt.close(fig)

    fig = plot_entropy_correlation(cluster_table)
    fig.savefig(fig_dir / f"{label}_entropy_correlation.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / f"{label}_entropy_correlation.pdf", bbox_inches="tight")
    plt.close(fig)

    feature_sets = make_mixing_feature_sets(adjusted=True)
    results = {}
    splits = {}
    metrics_rows = []

    for model_name, feature_set in feature_sets.items():
        X = cluster_table[feature_set.all_features].copy()
        y = cluster_table["is_large_cluster"].copy()
        groups = cluster_table["window_lineage_group"].copy()

        split = make_split_data(X, y, groups, config)
        result = fit_binary_model(
            name=f"{label}_{model_name}",
            X_train=split.X_train,
            y_train=split.y_train,
            X_test=split.X_test,
            y_test=split.y_test,
            config=config,
            numeric_features=feature_set.numeric_features,
            categorical_features=feature_set.categorical_features,
            negative_name="small_cluster",
            positive_name="large_cluster",
        )
        result.importance = permutation_importance_table(
            result.model,
            split.X_test,
            split.y_test,
            feature_set.all_features,
            config,
            scoring="average_precision",
        )
        save_binary_result(result, out, save_model=save_models)

        prediction_table = cluster_table.loc[
            split.X_test.index,
            [
                "cluster_uid",
                "cluster_type",
                "reported_cluster_size",
                "observed_cluster_size",
                "window_id",
                "pango_lineage",
                "who_voc_model",
            ],
        ].copy()
        prediction_table["true_is_large_cluster"] = split.y_test.to_numpy()
        prediction_table["pred_is_large_cluster"] = result.predictions
        prediction_table["pred_large_cluster_probability"] = result.probabilities
        prediction_table.to_csv(out / f"{label}_{model_name}_holdout_predictions.csv", index=False)

        fig = plot_binary_diagnostics(result, split.y_test, title_prefix=f"{label}: {model_name}")
        fig.savefig(fig_dir / f"{label}_{model_name}_diagnostics.png", dpi=300, bbox_inches="tight")
        fig.savefig(fig_dir / f"{label}_{model_name}_diagnostics.pdf", bbox_inches="tight")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(result.importance))))
        plot_importance(result.importance, f"{label}: {model_name} permutation importance", ax=ax)
        fig.savefig(
            fig_dir / f"{label}_{model_name}_permutation_importance.png",
            dpi=300,
            bbox_inches="tight",
        )
        fig.savefig(fig_dir / f"{label}_{model_name}_permutation_importance.pdf", bbox_inches="tight")
        plt.close(fig)

        results[model_name] = result
        splits[model_name] = split
        metrics_rows.append(
            {
                "analysis": label,
                "membership_definition": membership_definition,
                "model": model_name,
                "splitter": split.splitter_name,
                **result.metrics,
            }
        )

    metrics = pd.DataFrame(metrics_rows)
    metrics.to_csv(out / f"{label}_model_metrics.csv", index=False)

    effect_name = "context_adjusted"
    effect_result = results[effect_name]
    effect_split = splits[effect_name]
    X_effect, _ = sample_for_importance(effect_split.X_test, effect_split.y_test, config)

    pdp_curves = partial_dependence_curves(
        effect_result.model,
        X_effect,
        ENTROPY_FEATURES,
        grid_resolution=25,
        positive_class_index=1,
    )
    ale_curves = accumulated_local_effects_curves(
        effect_result.model,
        X_effect,
        ENTROPY_FEATURES,
        bins=20,
        positive_class_index=1,
    )
    pdp_curves.to_csv(out / f"{label}_context_adjusted_partial_dependence_curves.csv", index=False)
    ale_curves.to_csv(out / f"{label}_context_adjusted_ale_curves.csv", index=False)

    fig = plot_partial_dependence_curves(pdp_curves)
    fig.savefig(
        fig_dir / f"{label}_context_adjusted_partial_dependence_curves.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        fig_dir / f"{label}_context_adjusted_partial_dependence_curves.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    fig = plot_ale_curves(ale_curves)
    fig.savefig(fig_dir / f"{label}_context_adjusted_ale_curves.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / f"{label}_context_adjusted_ale_curves.pdf", bbox_inches="tight")
    plt.close(fig)

    save_run_metadata(
        config,
        feature_spec,
        out,
        extra={
            "analysis": "cluster_composition_mixing_entropy",
            "analysis_label": label,
            "membership_definition": membership_definition,
            "outcome": "large_cluster_vs_small_cluster_non_singletons",
            "entropy_features": ENTROPY_FEATURES,
            "entropy_normalisation": "log(min(K_global, n_valid_in_cluster))",
            "split_group": "window_id x pango_lineage",
            "n_non_singleton_clusters": int(len(cluster_table)),
            "output_dir": str(out),
        },
    )

    return {
        "cluster_table": cluster_table,
        "cluster_counts": cluster_counts,
        "entropy_summary": entropy_summary,
        "metrics": metrics,
        "results": results,
        "splits": splits,
        "pdp_curves": pdp_curves,
        "ale_curves": ale_curves,
    }
