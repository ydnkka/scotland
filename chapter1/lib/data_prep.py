"""Data preparation for Chapter 1.

Builds the cluster-level analysis table used by all Chapter 1 fits.

The two structurally important quantities are:

* per-cluster *observed* pairwise discordance, computed from the sequence
  rows belonging to that cluster.
* per-(window, lineage) *expected* pairwise discordance, computed from all
  sequence rows in the same analysis window and Pango lineage.

Excess mixing is observed minus expected; positive values mean clusters
are more sociodemographically mixed than would be expected if cases were
drawn at random from the same lineage and analysis window.

Two builders are exposed:

* :func:`build_cluster_table` — primary table used by the overall and wave
  analyses.  Discordance is computed for age, sex, SIMD quintile, and the
  joint age × sex profile and SIMD × age × sex profile.
* :func:`build_domain_cluster_table` — adds per-SIMD-domain quintile mixing
  for the domain analysis.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from patsy import dmatrix

from .constants import (
    CALENDAR_SPLINE_DF,
    DOMAINS,
    LINEAGE_MIN_CLUSTERS,
    MIXING_VARIABLES,
    MATRIX_VARIABLES,
    WAVE_LABELS,
    PRIMARY_RESOLUTION,
    QC_DEFAULT,
    SEQUENCE_COLUMNS,
    WAVE_ORDER,
)
from .estimators import logit_clipped, zscore


# ---------------------------------------------------------------------------
# Repo bootstrap so utils.data can be imported when this script is run direct
# ---------------------------------------------------------------------------


def _bootstrap_repo_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in [here, *here.parents]:
        if (candidate / "config.yaml").exists():
            root_str = str(candidate)
            if root_str not in sys.path:
                sys.path.insert(0, root_str)
            return candidate
    raise FileNotFoundError("Could not locate config.yaml.")


REPO_ROOT = _bootstrap_repo_root()

from utils.data import (  # noqa: E402
    Paths as DataPaths,
    load_analysis_columns_pandas,
    load_simd_columns_pandas,
    repo_root as _repo_root,
)


def repo_root(start: Path | None = None) -> Path:
    return _repo_root(start)


def analysis_dataset_path(root: Path) -> Path:
    return DataPaths.from_config(root).analysis_dataset


# ---------------------------------------------------------------------------
# Wave assignment
# ---------------------------------------------------------------------------


def assign_wave(lineage: str) -> str:
    """Bucket a Pango lineage into the published wave grouping."""
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


# ---------------------------------------------------------------------------
# Sequence loading
# ---------------------------------------------------------------------------


def read_sequence_rows(
    qc: str | None = QC_DEFAULT,
    primary_resolution: float = PRIMARY_RESOLUTION,
    include_domains: bool = False,
) -> pd.DataFrame:
    """Load sequence rows for the Chapter 1 cluster table."""
    columns = list(SEQUENCE_COLUMNS)
    if include_domains:
        domain_ranks = [spec["rank_col"] for spec in DOMAINS.values()]
        columns = list(dict.fromkeys([*columns, *domain_ranks]))

    seq = load_analysis_columns_pandas(
        columns=columns,
        resolution=primary_resolution,
        qc=qc,
    )

    categorical = [
        "cluster_id", "sequence_id", "window_id", "datazone", "pango_lineage",
        "nextclade_qc", "age_band", "sex", "dz_simd_quintile", "dz_simd_decile",
    ]
    for col in categorical:
        if col in seq.columns:
            seq[col] = seq[col].astype("category")

    seq["collection_date"] = pd.to_datetime(seq["collection_date"])
    seq["wn_mid_date"] = pd.to_datetime(seq["wn_mid_date"])
    seq["wave_group"] = (
        seq["pango_lineage"].astype(str).map(assign_wave).astype("category")
    )

    # Demographic profile (age × sex)
    complete_demo = seq[["age_band", "sex"]].notna().all(axis=1)
    seq["demographic_profile"] = pd.NA
    seq.loc[complete_demo, "demographic_profile"] = (
        seq.loc[complete_demo, "age_band"].astype(str)
        + "|"
        + seq.loc[complete_demo, "sex"].astype(str)
    )
    seq["demographic_profile"] = seq["demographic_profile"].astype("category")

    # Sociodemographic profile (SIMD quintile × age × sex)
    complete_sd = seq[["dz_simd_quintile", "age_band", "sex"]].notna().all(axis=1)
    seq["socio_demographic_profile"] = pd.NA
    seq.loc[complete_sd, "socio_demographic_profile"] = (
        seq.loc[complete_sd, "dz_simd_quintile"].astype(str)
        + "|"
        + seq.loc[complete_sd, "age_band"].astype(str)
        + "|"
        + seq.loc[complete_sd, "sex"].astype(str)
    )
    seq["socio_demographic_profile"] = seq["socio_demographic_profile"].astype(
        "category"
    )

    if include_domains:
        maxima = _domain_rank_maxima()
        for domain, spec in DOMAINS.items():
            q_col = f"{domain}_domain_quintile"
            if domain == "overall":
                seq[q_col] = seq["dz_simd_quintile"].astype("category")
            else:
                seq[q_col] = _rank_to_quintile(seq[spec["rank_col"]], maxima[domain])

    return seq


def _domain_rank_maxima() -> dict[str, float]:
    cols = [spec["rank_col"] for spec in DOMAINS.values()]
    simd = load_simd_columns_pandas(columns=cols)
    return {
        domain: float(simd[spec["rank_col"]].max())
        for domain, spec in DOMAINS.items()
    }


def _rank_to_quintile(rank: pd.Series, max_rank: float) -> pd.Series:
    quintile = np.ceil(rank.astype(float) / (max_rank / 5.0))
    return quintile.clip(1, 5).astype("Int64").astype("category")


# ---------------------------------------------------------------------------
# Pairwise discordance
# ---------------------------------------------------------------------------


def _pairwise_discordance_from_counts(
    counts: pd.DataFrame,
    group_cols: list[str],
) -> pd.DataFrame:
    totals = (
        counts.groupby(group_cols, observed=True)["n"]
        .sum()
        .rename("n_valid")
    )
    same_pairs = (
        counts.assign(same_pairs=counts["n"] * (counts["n"] - 1))
        .groupby(group_cols, observed=True)["same_pairs"]
        .sum()
    )
    out = pd.concat([totals, same_pairs], axis=1).reset_index()
    denom = out["n_valid"] * (out["n_valid"] - 1)
    out["discordance"] = np.nan
    mask = denom > 0
    out.loc[mask, "discordance"] = (
        1 - out.loc[mask, "same_pairs"] / denom.loc[mask]
    )
    return out.drop(columns=["same_pairs"])


def observed_cluster_discordance(
    seq: pd.DataFrame, variable: str, prefix: str,
) -> pd.DataFrame:
    """Per-cluster pairwise discordance for ``variable``."""
    counts = (
        seq.dropna(subset=[variable])
        .groupby(["cluster_id", variable], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = _pairwise_discordance_from_counts(counts, ["cluster_id"])
    return out.rename(columns={
        "n_valid": f"{prefix}_n_valid",
        "discordance": f"{prefix}_discordance",
    })


def expected_stratum_discordance(
    seq: pd.DataFrame, variable: str, prefix: str,
) -> pd.DataFrame:
    """Per-(window, lineage) discordance — the random-assembly expectation."""
    strata = ["window_id", "pango_lineage"]
    counts = (
        seq.dropna(subset=[variable])
        .groupby(strata + [variable], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = _pairwise_discordance_from_counts(counts, strata)
    return out.rename(columns={
        "n_valid": f"{prefix}_stratum_n_valid",
        "discordance": f"{prefix}_expected_discordance",
    })


# ---------------------------------------------------------------------------
# Marginal cluster entropy (for the null-residual sensitivity)
# ---------------------------------------------------------------------------


def _cluster_entropy(
    seq: pd.DataFrame, variable: str, prefix: str,
) -> pd.DataFrame:
    """Shannon entropy of cluster composition for ``variable``."""
    counts = (
        seq.dropna(subset=[variable])
        .groupby(["cluster_id", variable], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    totals = counts.groupby("cluster_id", observed=True)["n"].transform("sum")
    p = counts["n"] / totals.replace(0, np.nan)
    counts["plogp"] = np.where(p > 0, -p * np.log(p), 0.0)
    entropy = (
        counts.groupby("cluster_id", observed=True)["plogp"]
        .sum()
        .rename(f"{prefix}_entropy")
        .reset_index()
    )
    return entropy


# ---------------------------------------------------------------------------
# Calendar spline + lineage pooling
# ---------------------------------------------------------------------------


def _attach_calendar_spline(
    clusters: pd.DataFrame, calendar_spline_df: int,
) -> pd.DataFrame:
    calendar = dmatrix(
        f"bs(window_idx, df={calendar_spline_df}, degree=3, "
        "include_intercept=False) - 1",
        clusters,
        return_type="dataframe",
    )
    calendar.columns = [
        f"calendar_spline_{i + 1}" for i in range(calendar.shape[1])
    ]
    return pd.concat(
        [clusters.reset_index(drop=True), calendar.reset_index(drop=True)],
        axis=1,
    )


def _pool_lineages(
    clusters: pd.DataFrame, lineage_min_clusters: int,
) -> tuple[pd.DataFrame, int, int]:
    counts = clusters["pango_lineage"].astype(str).value_counts()
    common = set(counts[counts >= lineage_min_clusters].index)
    clusters = clusters.copy()
    clusters["lineage_model"] = np.where(
        clusters["pango_lineage"].astype(str).isin(common),
        clusters["pango_lineage"].astype(str),
        "Other rare lineages",
    )
    return clusters, len(counts), len(common)


# ---------------------------------------------------------------------------
# Cluster table — primary
# ---------------------------------------------------------------------------


def build_cluster_table(
    seq: pd.DataFrame,
    lineage_min_clusters: int = LINEAGE_MIN_CLUSTERS,
    calendar_spline_df: int = CALENDAR_SPLINE_DF,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Build the primary Chapter 1 cluster table."""
    required = [
        "cluster_id", "sequence_id", "window_id", "window_idx",
        "collection_date", "datazone", "pango_lineage",
        "dz_simd_rank", "dz_cum_incidence_per_capita",
        "dz_cum_prop_sequenced", "wn_prop_sequenced", "dz_7d_test_positivity",
    ]
    before = len(seq)
    seq = seq.dropna(subset=required).copy()
    dropped = before - len(seq)

    clusters = (
        seq.groupby("cluster_id", observed=True, sort=False)
        .agg(
            cluster_size=("sequence_id", "nunique"),
            cluster_n_datazones=("datazone", "nunique"),
            cluster_start_date=("collection_date", "min"),
            cluster_end_date=("collection_date", "max"),
            resolution=("resolution", "first"),
            window_id=("window_id", "first"),
            window_idx=("window_idx", "first"),
            wn_mid_date=("wn_mid_date", "first"),
            pango_lineage=("pango_lineage", "first"),
            wave_group=("wave_group", "first"),
            mean_simd_rank=("dz_simd_rank", "mean"),
            mean_local_incidence_per_capita=("dz_cum_incidence_per_capita", "mean"),
            mean_local_seq_fraction=("dz_cum_prop_sequenced", "mean"),
            mean_window_seq_fraction=("wn_prop_sequenced", "mean"),
            mean_test_positivity=("dz_7d_test_positivity", "mean"),
            wn_no_sequences=("wn_no_sequences", "first"),
            health_board=(
                "dz_health_board_code",
                lambda x: x.mode().iloc[0] if not x.mode().empty else pd.NA,
            ),
        )
        .reset_index()
    )

    clusters["cluster_size_gt1"] = (clusters["cluster_size"] > 1).astype(int)
    clusters["datazones_gt1"] = (clusters["cluster_n_datazones"] > 1).astype(int)
    clusters["cluster_size_excess"] = clusters["cluster_size"] - 1
    clusters["geographic_spread"] = clusters["cluster_n_datazones"]
    clusters["duration_days"] = (
        clusters["cluster_end_date"] - clusters["cluster_start_date"]
    ).dt.days.astype(int)

    # Mixing: observed + expected discordance for each variable + cluster entropy
    for prefix, spec in MIXING_VARIABLES.items():
        obs = observed_cluster_discordance(seq, spec["column"], prefix)
        exp = expected_stratum_discordance(seq, spec["column"], prefix)
        ent = _cluster_entropy(seq, spec["column"], prefix)
        clusters = clusters.merge(obs, on="cluster_id", how="left")
        clusters = clusters.merge(
            exp, on=["window_id", "pango_lineage"], how="left",
        )
        clusters = clusters.merge(ent, on="cluster_id", how="left")
        clusters[f"{prefix}_excess_discordance"] = (
            clusters[f"{prefix}_discordance"]
            - clusters[f"{prefix}_expected_discordance"]
        )

    # Transforms
    clusters["deprivation_raw"] = -clusters["mean_simd_rank"]
    clusters["local_incidence_log"] = np.log1p(
        clusters["mean_local_incidence_per_capita"].clip(lower=0) * 1000
    )
    clusters["local_seq_fraction_logit"] = logit_clipped(
        clusters["mean_local_seq_fraction"]
    )
    clusters["window_seq_fraction_logit"] = logit_clipped(
        clusters["mean_window_seq_fraction"]
    )
    clusters["test_positivity_logit"] = logit_clipped(
        clusters["mean_test_positivity"].fillna(0)
    )
    clusters["log_cluster_size"] = np.log(clusters["cluster_size"])

    scaling_rows: list[dict[str, object]] = []
    transforms = {
        "deprivation_z":          "deprivation_raw",
        "local_incidence_z":      "local_incidence_log",
        "local_seq_fraction_z":   "local_seq_fraction_logit",
        "window_seq_fraction_z":  "window_seq_fraction_logit",
        "test_positivity_z":      "test_positivity_logit",
        "log_cluster_size_z":     "log_cluster_size",
    }
    transforms.update({
        f"{prefix}_excess_mixing_z": f"{prefix}_excess_discordance"
        for prefix in MIXING_VARIABLES
    })
    for z_col, raw_col in transforms.items():
        clusters[z_col], mean, sd = zscore(clusters[raw_col])
        scaling_rows.append({
            "standardised_column": z_col,
            "source_column": raw_col,
            "source_mean": mean,
            "source_sd": sd,
        })

    clusters, n_lin_total, n_lin_common = _pool_lineages(
        clusters, lineage_min_clusters,
    )
    clusters = _attach_calendar_spline(clusters, calendar_spline_df)

    # Reference-coded wave_group with the chosen reference dropped via patsy
    clusters["wave_group"] = pd.Categorical(
        clusters["wave_group"].astype(str),
        categories=list(WAVE_ORDER) + ["Other"],
    )

    scaling = pd.DataFrame(scaling_rows)
    scaling.attrs["dropped_sequence_rows_missing_model_fields"] = dropped
    scaling.attrs["lineages_total"] = n_lin_total
    scaling.attrs["lineages_modelled"] = n_lin_common + int(
        n_lin_common < n_lin_total
    )
    scaling.attrs["lineage_min_clusters"] = lineage_min_clusters
    return clusters, scaling, dropped


# ---------------------------------------------------------------------------
# Cluster table — SIMD-domain extension
# ---------------------------------------------------------------------------


def build_domain_cluster_table(
    seq: pd.DataFrame,
    lineage_min_clusters: int = LINEAGE_MIN_CLUSTERS,
    calendar_spline_df: int = CALENDAR_SPLINE_DF,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Extend the cluster table with per-domain SIMD-quintile mixing."""
    clusters, scaling, dropped = build_cluster_table(
        seq,
        lineage_min_clusters=lineage_min_clusters,
        calendar_spline_df=calendar_spline_df,
    )

    scaling_rows: list[dict[str, object]] = []
    for domain in DOMAINS:
        variable = f"{domain}_domain_quintile"
        if variable not in seq.columns:
            continue
        prefix = f"{domain}_domain"
        obs = observed_cluster_discordance(seq, variable, prefix)
        exp = expected_stratum_discordance(seq, variable, prefix)
        clusters = clusters.merge(obs, on="cluster_id", how="left")
        clusters = clusters.merge(
            exp, on=["window_id", "pango_lineage"], how="left",
        )
        clusters[f"{prefix}_excess_discordance"] = (
            clusters[f"{prefix}_discordance"]
            - clusters[f"{prefix}_expected_discordance"]
        )
        z_col = f"{prefix}_excess_mixing_z"
        clusters[z_col], mean, sd = zscore(
            clusters[f"{prefix}_excess_discordance"]
        )
        scaling_rows.append({
            "standardised_column": z_col,
            "source_column": f"{prefix}_excess_discordance",
            "source_mean": mean,
            "source_sd": sd,
        })

    domain_scaling = pd.DataFrame(scaling_rows)
    combined_scaling = pd.concat(
        [scaling, domain_scaling], ignore_index=True,
    )
    for key, value in scaling.attrs.items():
        combined_scaling.attrs[key] = value
    return clusters, combined_scaling, dropped


# ---------------------------------------------------------------------------
# Descriptives
# ---------------------------------------------------------------------------


def summarise_dataset(
    seq: pd.DataFrame,
    clusters: pd.DataFrame,
    qc: str | None,
    primary_resolution: float,
    dropped: int,
) -> pd.DataFrame:
    """Compact descriptive summary table written to ``tables/``."""
    non_singleton = int((clusters["cluster_size"] >= 2).sum())
    rows: list[dict[str, object]] = [
        {"measure": "sequence_rows_used", "statistic": "count", "value": len(seq)},
        {"measure": "sequence_rows_dropped_missing_model_fields", "statistic": "count", "value": dropped},
        {"measure": "clusters", "statistic": "count", "value": len(clusters)},
        {"measure": "non_singleton_clusters", "statistic": "count", "value": non_singleton},
        {"measure": "primary_leiden_resolution", "statistic": "value", "value": primary_resolution},
        {"measure": "windows", "statistic": "count", "value": clusters["window_id"].nunique()},
        {"measure": "pango_lineages_raw", "statistic": "count", "value": clusters["pango_lineage"].nunique()},
        {"measure": "pango_lineage_model_levels", "statistic": "count", "value": clusters["lineage_model"].nunique()},
        {"measure": "qc_filter", "statistic": "value", "value": qc or "none"},
    ]

    outcomes = [
        ("cluster_size", "Cluster size"),
        ("cluster_n_datazones", "Distinct datazones"),
        ("duration_days", "Duration days"),
    ]
    percentiles = [0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
    for col, label in outcomes:
        desc = clusters[col].describe(percentiles=percentiles)
        rows.extend(
            {"measure": label, "statistic": str(stat), "value": float(value)}
            for stat, value in desc.items()
        )

    for prefix, spec in MIXING_VARIABLES.items():
        col = f"{prefix}_excess_discordance"
        values = clusters.loc[clusters["cluster_size"] >= 2, col].dropna()
        desc = values.describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9])
        rows.extend(
            {
                "measure": f"{spec['short_label']} excess mixing",
                "statistic": str(stat),
                "value": float(value),
            }
            for stat, value in desc.items()
        )

    wave_counts = (
        clusters.loc[clusters["cluster_size"] >= 2, "wave_group"]
        .astype(str)
        .value_counts()
        .sort_index()
    )
    rows.extend(
        {
            "measure": f"non_singleton_clusters_wave_{wave}",
            "statistic": "count",
            "value": int(value),
        }
        for wave, value in wave_counts.items()
    )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Observed-vs-expected pair-probability matrices
# ---------------------------------------------------------------------------


def _observed_ordered_pairs(
    cluster_counts: pd.DataFrame,
    levels: list,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide = (
        cluster_counts.pivot_table(
            index=["cluster_id", "wave_group", "window_id", "pango_lineage"],
            columns="category",
            values="n",
            aggfunc="sum",
            fill_value=0,
            observed=True,
        )
        .reindex(columns=levels, fill_value=0)
        .astype(np.int32)
        .reset_index()
    )
    wide["n_valid"] = wide[levels].sum(axis=1)
    wide = wide[wide["n_valid"] >= 2].copy()

    rows: list[dict[str, object]] = []
    for left in levels:
        for right in levels:
            values = wide[left].astype(np.int64) * wide[right].astype(np.int64)
            if left == right:
                values = wide[left].astype(np.int64) * (wide[left].astype(np.int64) - 1)
            by_wave = values.groupby(wide["wave_group"], observed=True).sum()
            for wave, n_pairs in by_wave.items():
                rows.append(
                    {
                        "wave_group": wave,
                        "category_i": left,
                        "category_j": right,
                        "observed_pairs": float(n_pairs),
                    }
                )
    return pd.DataFrame(rows), wide


def _expected_ordered_pairs(
    cluster_wide: pd.DataFrame,
    stratum_counts: pd.DataFrame,
    levels: list,
) -> pd.DataFrame:
    cluster_wide = cluster_wide.copy()
    cluster_wide["ordered_pairs"] = cluster_wide["n_valid"] * (cluster_wide["n_valid"] - 1)
    stratum_cols = ["wave_group", "window_id", "pango_lineage"]
    stratum_pair_totals = (
        cluster_wide.groupby(stratum_cols, observed=True)["ordered_pairs"]
        .sum()
        .rename("cluster_ordered_pairs")
        .reset_index()
    )

    stratum_wide = (
        stratum_counts.pivot_table(
            index=stratum_cols,
            columns="category",
            values="n",
            aggfunc="sum",
            fill_value=0,
            observed=True,
        )
        .reindex(columns=levels, fill_value=0)
        .astype(np.int64)
        .reset_index()
    )
    stratum_wide["stratum_n"] = stratum_wide[levels].sum(axis=1)
    stratum_wide = stratum_wide.merge(stratum_pair_totals, on=stratum_cols, how="inner")
    denom = stratum_wide["stratum_n"] * (stratum_wide["stratum_n"] - 1)

    rows: list[dict[str, object]] = []
    for left in levels:
        for right in levels:
            numerator = stratum_wide[left].astype(np.float64) * stratum_wide[right].astype(np.float64)
            if left == right:
                numerator = stratum_wide[left].astype(np.float64) * (
                    stratum_wide[left].astype(np.float64) - 1
                )
            expected = stratum_wide["cluster_ordered_pairs"] * numerator / denom
            expected = expected.replace([np.inf, -np.inf], np.nan).fillna(0)
            by_wave = expected.groupby(stratum_wide["wave_group"], observed=True).sum()
            for wave, n_pairs in by_wave.items():
                rows.append(
                    {
                        "wave_group": wave,
                        "category_i": left,
                        "category_j": right,
                        "expected_pairs": float(n_pairs),
                    }
                )
    return pd.DataFrame(rows)


def build_matrix_for_variable(seq: pd.DataFrame, variable: str) -> pd.DataFrame:
    """Build the observed-vs-expected pair-probability matrix for one variable."""
    spec = MATRIX_VARIABLES[variable]
    levels = list(spec["levels"])
    work = seq.dropna(subset=[spec["column"]]).copy()
    work = work[work["wave_group"].isin(WAVE_ORDER)].copy()
    work["category"] = work[spec["column"]]

    cluster_counts = (
        work.groupby(
            ["cluster_id", "wave_group", "window_id", "pango_lineage", "category"],
            observed=True,
        )
        .size()
        .rename("n")
        .reset_index()
    )
    observed, cluster_wide = _observed_ordered_pairs(cluster_counts, levels)

    stratum_counts = (
        work.groupby(["wave_group", "window_id", "pango_lineage", "category"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    expected = _expected_ordered_pairs(cluster_wide, stratum_counts, levels)

    matrix = observed.merge(
        expected,
        on=["wave_group", "category_i", "category_j"],
        how="outer",
    ).fillna({"observed_pairs": 0, "expected_pairs": 0})
    matrix["variable"] = variable
    matrix["variable_label"] = spec["label"]

    overall = (
        matrix.groupby(["variable", "variable_label", "category_i", "category_j"], observed=True)[
            ["observed_pairs", "expected_pairs"]
        ]
        .sum()
        .reset_index()
    )
    overall["wave_group"] = "Overall"
    matrix = pd.concat([matrix, overall], ignore_index=True)

    totals = (
        matrix.groupby(["variable", "wave_group"], observed=True)[
            ["observed_pairs", "expected_pairs"]
        ]
        .sum()
        .rename(
            columns={
                "observed_pairs": "total_observed_pairs",
                "expected_pairs": "total_expected_pairs",
            }
        )
        .reset_index()
    )
    matrix = matrix.merge(totals, on=["variable", "wave_group"], how="left")
    matrix["observed_probability"] = matrix["observed_pairs"] / matrix["total_observed_pairs"]
    matrix["expected_probability"] = matrix["expected_pairs"] / matrix["total_expected_pairs"]
    matrix["excess_probability"] = matrix["observed_probability"] - matrix["expected_probability"]
    matrix["excess_percentage_points"] = matrix["excess_probability"] * 100
    matrix["observed_expected_ratio"] = matrix["observed_probability"] / matrix["expected_probability"]
    matrix["wave_label"] = matrix["wave_group"].map(
        lambda w: "Overall" if w == "Overall" else WAVE_LABELS.get(w, w)
    )
    return matrix[
        [
            "variable", "variable_label", "wave_group", "wave_label",
            "category_i", "category_j",
            "observed_pairs", "expected_pairs",
            "observed_probability", "expected_probability",
            "excess_probability", "excess_percentage_points",
            "observed_expected_ratio",
        ]
    ]


__all__ = [
    "REPO_ROOT",
    "analysis_dataset_path",
    "assign_wave",
    "build_cluster_table",
    "build_domain_cluster_table",
    "expected_stratum_discordance",
    "observed_cluster_discordance",
    "read_sequence_rows",
    "repo_root",
    "summarise_dataset",
    "build_matrix_for_variable"
]
