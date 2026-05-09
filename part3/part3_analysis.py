"""Build Part 3 policy-period analyses for the Scotland clustering project.

The analysis is descriptive by design.  Policy periods are used as epidemic
context and as anchors for selected interrupted time-series summaries, not as
causal interventions.  The Alpha case study is rebuilt directly from the
sequence-level processed dataset and the raw Nextclade mutation table.

Run from the repository root:

    conda run -n PhD python part3/part3_analysis.py
"""

from __future__ import annotations

import os
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils import data as data_utils
from utils import policy


TABLE_DIR = ROOT / "part3" / "tables"
CACHE_DIR = ROOT / "part3" / "cache"

PRIMARY_RESOLUTION = data_utils.PRIMARY_RESOLUTION
PRIMARY_QC = "good"

SELECTED_PERIODS = ["P3", "T1", "F5", "L2", "SL", "L0", "NN"]
CONTEXT_PERIODS = ["OM", "FE", "PR"]

TRANSITIONS = {
    "t1_onset": {
        "label": "T1 onset",
        "from_to": "P3 -> T1",
        "date": pd.Timestamp("2020-10-02"),
    },
    "l2_to_sl": {
        "label": "L2 to SL",
        "from_to": "L2 -> SL",
        "date": pd.Timestamp("2021-04-02"),
    },
    "nn_onset": {
        "label": "NN onset",
        "from_to": "L0 -> NN",
        "date": pd.Timestamp("2021-08-09"),
    },
}

ITS_OUTCOMES = {
    "median_log_cluster_size": "Median log cluster size",
    "median_log_datazones": "Median log datazones",
    "mean_simd_excess_discordance": "Mean SIMD excess discordance",
    "mean_age_excess_discordance": "Mean age excess discordance",
}

MUTATION_MARKERS = {
    "S:N501Y": "s_n501y",
    "S:A222V": "s_a222v",
    "S:P681H": "s_p681h",
    "S:A570D": "s_a570d",
    "S:D1118H": "s_d1118h",
    "N:R203K": "n_r203k",
    "N:G204R": "n_g204r",
}


@dataclass(frozen=True)
class GrowthFit:
    analysis: str
    marker: str
    marker_slug: str
    period_code: str
    period_label: str
    weight_scheme: str
    n_windows: int
    n_success: int
    n_total: int
    start_date: pd.Timestamp
    end_date: pd.Timestamp
    origin_date: pd.Timestamp
    intercept: float
    slope_per_day: float
    intercept_se: float
    slope_se_per_day: float
    intercept_pvalue: float
    slope_pvalue: float
    slope_ci_low_per_day: float
    slope_ci_high_per_day: float
    aic: float
    pseudo_r2: float


def setup_environment() -> None:
    """Use writable cache paths for matplotlib and friends."""
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)


def read_config() -> dict:
    with open(ROOT / "config.yaml") as f:
        return yaml.safe_load(f)


def nextclade_tsv_path() -> Path:
    """Return the Nextclade TSV path from the config."""
    cfg = read_config()
    return ROOT / cfg["data"]["raw"]["nextclade_tsv"]


def safe_log(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    return np.log(values.where(values > 0))


def attach_policy(df: pd.DataFrame, date_col: str) -> pd.DataFrame:
    out = df.copy()
    out[date_col] = pd.to_datetime(out[date_col])
    return policy.attach_period_pandas(out, date_col)


def ordered_policy_table() -> pd.DataFrame:
    periods = policy.POLICY_PERIODS_PD.copy()
    periods["period_order"] = np.arange(len(periods))
    return periods


def load_policy_cluster_table() -> pd.DataFrame:
    """Load the Part 1 primary cluster cache and add Part 3 policy fields."""
    cluster = data_utils.load_main_cluster_table(root=ROOT).copy()
    for col in ["wn_mid_date", "cluster_start_date", "cluster_end_date"]:
        if col in cluster:
            cluster[col] = pd.to_datetime(cluster[col])

    if "log_cluster_size" not in cluster:
        cluster["log_cluster_size"] = safe_log(cluster["cluster_size"])
    cluster["log_datazones"] = safe_log(cluster["cluster_n_datazones"])
    cluster["is_non_singleton"] = cluster["cluster_size"] > 1

    cluster = attach_policy(cluster, "wn_mid_date")
    cluster["is_selected_policy_phase"] = cluster["policy_period"].isin(SELECTED_PERIODS)
    cluster["is_context_policy_phase"] = cluster["policy_period"].isin(CONTEXT_PERIODS)

    out_path = CACHE_DIR / "policy_cluster_table.parquet"
    cluster.to_parquet(out_path, index=False)
    return cluster


def summarise_periods(cluster: pd.DataFrame) -> pd.DataFrame:
    """Create a full policy-period descriptive table."""

    def q25(x: pd.Series) -> float:
        return x.quantile(0.25)

    def q75(x: pd.Series) -> float:
        return x.quantile(0.75)

    non_singleton = cluster[cluster["is_non_singleton"]].copy()
    grouped_all = (
        cluster.groupby("policy_period", dropna=False)
        .agg(
            observed_start=("wn_mid_date", "min"),
            observed_end=("wn_mid_date", "max"),
            total_cluster_rows=("cluster_id", "size"),
            total_sequences_represented=("cluster_size", "sum"),
            singleton_clusters=("is_non_singleton", lambda s: int((~s).sum())),
            non_singleton_clusters=("is_non_singleton", "sum"),
            mean_cluster_size=("cluster_size", "mean"),
            median_cluster_size=("cluster_size", "median"),
            q25_cluster_size=("cluster_size", q25),
            q75_cluster_size=("cluster_size", q75),
            median_log_cluster_size=("log_cluster_size", "median"),
            mean_cluster_n_datazones=("cluster_n_datazones", "mean"),
            median_cluster_n_datazones=("cluster_n_datazones", "median"),
            median_log_datazones=("log_datazones", "median"),
            n_windows=("window_id", "nunique"),
        )
        .reset_index()
    )

    grouped_ns = (
        non_singleton.groupby("policy_period", dropna=False)
        .agg(
            non_singleton_median_cluster_size=("cluster_size", "median"),
            non_singleton_median_log_cluster_size=("log_cluster_size", "median"),
            non_singleton_median_datazones=("cluster_n_datazones", "median"),
            non_singleton_median_log_datazones=("log_datazones", "median"),
            mean_simd_excess_discordance=("simd_excess_discordance", "mean"),
            mean_age_excess_discordance=("age_excess_discordance", "mean"),
            mean_sex_excess_discordance=("sex_excess_discordance", "mean"),
            mean_profile_excess_discordance=("profile_excess_discordance", "mean"),
        )
        .reset_index()
    )

    periods = ordered_policy_table().rename(
        columns={
            "period_code": "policy_period",
            "period_label": "policy_period_label",
            "intensity": "policy_intensity",
        }
    )

    out = (
        periods.merge(grouped_all, on="policy_period", how="left")
        .merge(grouped_ns, on="policy_period", how="left")
        .sort_values("period_order")
    )

    count_cols = [
        "total_cluster_rows",
        "total_sequences_represented",
        "singleton_clusters",
        "non_singleton_clusters",
        "n_windows",
    ]
    for col in count_cols:
        out[col] = out[col].fillna(0).astype(int)
    out["singleton_fraction"] = np.where(
        out["total_cluster_rows"] > 0,
        out["singleton_clusters"] / out["total_cluster_rows"],
        np.nan,
    )
    out["chapter_role"] = np.select(
        [
            out["policy_period"].isin(SELECTED_PERIODS),
            out["policy_period"].isin(CONTEXT_PERIODS),
        ],
        ["selected phase", "context/supplement"],
        default="supplement",
    )

    out.to_csv(TABLE_DIR / "period_descriptives.csv", index=False)
    return out


def summarise_weekly(cluster: pd.DataFrame) -> pd.DataFrame:
    non_singleton = cluster[cluster["is_non_singleton"]].copy()

    all_counts = (
        cluster.groupby(["window_id", "window_idx", "wn_mid_date"], as_index=False)
        .agg(
            total_clusters=("cluster_id", "size"),
            singleton_clusters=("is_non_singleton", lambda s: int((~s).sum())),
            non_singleton_clusters=("is_non_singleton", "sum"),
            total_sequences_represented=("cluster_size", "sum"),
            wn_no_sequences=("wn_no_sequences", "max"),
            mean_window_seq_fraction=("mean_window_seq_fraction", "mean"),
        )
    )

    ns_summary = (
        non_singleton.groupby(["window_id", "window_idx", "wn_mid_date"], as_index=False)
        .agg(
            median_cluster_size=("cluster_size", "median"),
            median_log_cluster_size=("log_cluster_size", "median"),
            median_datazones=("cluster_n_datazones", "median"),
            median_log_datazones=("log_datazones", "median"),
            mean_simd_excess_discordance=("simd_excess_discordance", "mean"),
            mean_age_excess_discordance=("age_excess_discordance", "mean"),
            mean_sex_excess_discordance=("sex_excess_discordance", "mean"),
            mean_profile_excess_discordance=("profile_excess_discordance", "mean"),
        )
    )

    weekly = all_counts.merge(
        ns_summary, on=["window_id", "window_idx", "wn_mid_date"], how="left"
    )
    weekly["singleton_fraction"] = weekly["singleton_clusters"] / weekly["total_clusters"]
    weekly = attach_policy(weekly, "wn_mid_date")
    weekly = weekly.sort_values("window_idx")
    weekly.to_csv(TABLE_DIR / "weekly_summaries.csv", index=False)
    return weekly


def intensity_correlations(weekly: pd.DataFrame) -> pd.DataFrame:
    outcomes = [
        "total_clusters",
        "non_singleton_clusters",
        "singleton_fraction",
        "median_log_cluster_size",
        "median_log_datazones",
        "mean_simd_excess_discordance",
        "mean_age_excess_discordance",
        "mean_sex_excess_discordance",
        "mean_profile_excess_discordance",
    ]
    rows: list[dict] = []
    for outcome in outcomes:
        dat = weekly[["policy_intensity", outcome]].dropna()
        if len(dat) < 3 or dat["policy_intensity"].nunique() < 2:
            rho = np.nan
            pvalue = np.nan
        else:
            rho, pvalue = stats.spearmanr(dat["policy_intensity"], dat[outcome])
        rows.append(
            {
                "outcome": outcome,
                "n_weeks": len(dat),
                "spearman_rho": rho,
                "pvalue": pvalue,
                "interpretation": "descriptive/confounded by variant phase, surveillance, immunity, and calendar time",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "intensity_correlations.csv", index=False)
    return out


def fit_segmented_ols(
    weekly: pd.DataFrame,
    transition_slug: str,
    transition_date: pd.Timestamp,
    outcome: str,
    window_weeks: int,
) -> tuple[list[dict], pd.DataFrame]:
    lo = transition_date - pd.Timedelta(weeks=window_weeks)
    hi = transition_date + pd.Timedelta(weeks=window_weeks)
    dat = weekly.loc[
        weekly["wn_mid_date"].between(lo, hi),
        [
            "window_id",
            "window_idx",
            "wn_mid_date",
            "policy_period",
            "policy_intensity",
            outcome,
        ],
    ].copy()
    dat = dat.dropna(subset=[outcome])
    dat["t_weeks"] = (dat["wn_mid_date"] - transition_date).dt.days / 7.0
    dat["post"] = (dat["wn_mid_date"] >= transition_date).astype(int)
    dat["post_t_weeks"] = dat["post"] * dat["t_weeks"]

    rows: list[dict] = []
    dat[f"fitted_{outcome}"] = np.nan
    if len(dat) < 7 or dat["post"].nunique() < 2:
        for term in ["const", "t_weeks", "post", "post_t_weeks"]:
            rows.append(
                {
                    "transition": transition_slug,
                    "transition_label": TRANSITIONS[transition_slug]["label"],
                    "from_to": TRANSITIONS[transition_slug]["from_to"],
                    "transition_date": transition_date.date().isoformat(),
                    "window_weeks": window_weeks,
                    "outcome": outcome,
                    "outcome_label": ITS_OUTCOMES[outcome],
                    "term": term,
                    "estimate": np.nan,
                    "std_error": np.nan,
                    "pvalue": np.nan,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "n_weeks": len(dat),
                    "adj_r2": np.nan,
                }
            )
        return rows, dat

    x = sm.add_constant(dat[["t_weeks", "post", "post_t_weeks"]], has_constant="add")
    try:
        model = sm.OLS(dat[outcome].astype(float), x.astype(float)).fit(cov_type="HC1")
        dat[f"fitted_{outcome}"] = model.predict(x)
        conf = model.conf_int()
        for term in ["const", "t_weeks", "post", "post_t_weeks"]:
            rows.append(
                {
                    "transition": transition_slug,
                    "transition_label": TRANSITIONS[transition_slug]["label"],
                    "from_to": TRANSITIONS[transition_slug]["from_to"],
                    "transition_date": transition_date.date().isoformat(),
                    "window_weeks": window_weeks,
                    "outcome": outcome,
                    "outcome_label": ITS_OUTCOMES[outcome],
                    "term": term,
                    "estimate": model.params.get(term, np.nan),
                    "std_error": model.bse.get(term, np.nan),
                    "pvalue": model.pvalues.get(term, np.nan),
                    "ci_low": conf.loc[term, 0] if term in conf.index else np.nan,
                    "ci_high": conf.loc[term, 1] if term in conf.index else np.nan,
                    "n_weeks": int(model.nobs),
                    "adj_r2": model.rsquared_adj,
                }
            )
    except Exception as exc:  # pragma: no cover - retained for robust reruns
        warnings.warn(f"ITS model failed for {transition_slug} {outcome}: {exc}")
        for term in ["const", "t_weeks", "post", "post_t_weeks"]:
            rows.append(
                {
                    "transition": transition_slug,
                    "transition_label": TRANSITIONS[transition_slug]["label"],
                    "from_to": TRANSITIONS[transition_slug]["from_to"],
                    "transition_date": transition_date.date().isoformat(),
                    "window_weeks": window_weeks,
                    "outcome": outcome,
                    "outcome_label": ITS_OUTCOMES[outcome],
                    "term": term,
                    "estimate": np.nan,
                    "std_error": np.nan,
                    "pvalue": np.nan,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "n_weeks": len(dat),
                    "adj_r2": np.nan,
                }
            )
    return rows, dat


def run_its(weekly: pd.DataFrame) -> pd.DataFrame:
    coeff_rows: list[dict] = []
    for transition_slug, meta in TRANSITIONS.items():
        primary_weekly: pd.DataFrame | None = None
        for window_weeks in [8, 6, 10, 12]:
            for outcome in ITS_OUTCOMES:
                rows, fitted = fit_segmented_ols(
                    weekly=weekly,
                    transition_slug=transition_slug,
                    transition_date=meta["date"],
                    outcome=outcome,
                    window_weeks=window_weeks,
                )
                coeff_rows.extend(rows)
                if window_weeks == 8:
                    cols = [
                        "window_id",
                        "window_idx",
                        "wn_mid_date",
                        "policy_period",
                        "policy_intensity",
                        "t_weeks",
                        "post",
                        outcome,
                        f"fitted_{outcome}",
                    ]
                    fitted = fitted[cols].copy()
                    if primary_weekly is None:
                        primary_weekly = fitted
                    else:
                        primary_weekly = primary_weekly.merge(
                            fitted[
                                [
                                    "window_id",
                                    outcome,
                                    f"fitted_{outcome}",
                                ]
                            ],
                            on="window_id",
                            how="left",
                        )

        if primary_weekly is not None:
            primary_weekly = primary_weekly.sort_values("window_idx")
            primary_weekly.to_csv(
                TABLE_DIR / f"its_weekly_{transition_slug}.csv", index=False
            )

    coeffs = pd.DataFrame(coeff_rows)
    coeffs.to_csv(TABLE_DIR / "its_coefficients.csv", index=False)
    return coeffs


def load_sequence_table() -> pd.DataFrame:
    """Read the primary sequence-level table needed for the Alpha case study."""
    paths = data_utils.Paths.from_config(ROOT)
    columns = [
        "sequence_id",
        "window_id",
        "window_idx",
        "wn_start_date",
        "wn_mid_date",
        "wn_end_date",
        "wn_no_sequences",
        "wn_positive_tests",
        "wn_prop_sequenced",
        "resolution",
        "nextclade_qc",
        "cluster_id",
        "cluster_size",
        "cluster_n_datazones",
        "collection_date",
        "pango_lineage",
        "dz_health_board",
        "dz_local_authority",
        "hb_hospital_occupancy",
    ]
    try:
        seq = pd.read_parquet(
            paths.analysis_dataset,
            columns=columns,
            filters=[
                ("resolution", "==", PRIMARY_RESOLUTION),
                ("nextclade_qc", "==", PRIMARY_QC),
            ],
        )
    except Exception:
        seq = pd.read_parquet(paths.analysis_dataset, columns=columns)
        seq = seq[
            (seq["resolution"] == PRIMARY_RESOLUTION)
            & (seq["nextclade_qc"] == PRIMARY_QC)
        ].copy()

    for col in ["wn_start_date", "wn_mid_date", "wn_end_date", "collection_date"]:
        seq[col] = pd.to_datetime(seq[col])

    seq["pango_lineage"] = seq["pango_lineage"].fillna("unknown").astype(str)
    seq["is_alpha_pango"] = seq["pango_lineage"].str.startswith("B.1.1.7")
    seq["is_b1177_pango"] = seq["pango_lineage"].str.startswith("B.1.177")
    seq = attach_policy(seq, "wn_mid_date")
    seq.to_parquet(CACHE_DIR / "alpha_sequence_table.parquet", index=False)
    return seq


def alpha_phase_definitions(seq: pd.DataFrame) -> list[dict]:
    """Define phase windows by date and record the observed window labels."""
    window_table = (
        seq[["window_id", "window_idx", "wn_mid_date"]]
        .drop_duplicates()
        .sort_values("window_idx")
    )
    phase_specs = [
        (
            "cryptic_early",
            "Cryptic/early Alpha phase",
            pd.Timestamp("2020-10-27"),
            pd.Timestamp("2020-12-02"),
        ),
        (
            "multi_region_expansion",
            "Multi-region Alpha expansion",
            pd.Timestamp("2020-12-08"),
            pd.Timestamp("2020-12-23"),
        ),
        (
            "f5_l2_bridge",
            "F5/L2 bridge",
            pd.Timestamp("2020-12-23"),
            pd.Timestamp("2020-12-30"),
        ),
    ]

    phases: list[dict] = []
    for slug, label, start, end in phase_specs:
        wins = window_table[window_table["wn_mid_date"].between(start, end)].copy()
        if wins.empty:
            start_idx = (window_table["wn_mid_date"] - start).abs().idxmin()
            end_idx = (window_table["wn_mid_date"] - end).abs().idxmin()
            lo = min(window_table.loc[start_idx, "window_idx"], window_table.loc[end_idx, "window_idx"])
            hi = max(window_table.loc[start_idx, "window_idx"], window_table.loc[end_idx, "window_idx"])
            wins = window_table[window_table["window_idx"].between(lo, hi)].copy()
        phases.append(
            {
                "phase": slug,
                "phase_label": label,
                "target_start_date": start,
                "target_end_date": end,
                "start_window_idx": int(wins["window_idx"].min()),
                "end_window_idx": int(wins["window_idx"].max()),
                "start_window_id": wins.sort_values("window_idx")["window_id"].iloc[0],
                "end_window_id": wins.sort_values("window_idx")["window_id"].iloc[-1],
                "observed_start_mid_date": wins["wn_mid_date"].min(),
                "observed_end_mid_date": wins["wn_mid_date"].max(),
            }
        )
    return phases


def summarise_alpha_phases(seq: pd.DataFrame) -> pd.DataFrame:
    alpha = seq[seq["is_alpha_pango"]].copy()
    phases = alpha_phase_definitions(seq)
    rows = []
    for phase in phases:
        mask_phase = seq["window_idx"].between(
            phase["start_window_idx"], phase["end_window_idx"]
        )
        phase_all = seq[mask_phase].copy()
        phase_alpha = alpha[
            alpha["window_idx"].between(
                phase["start_window_idx"], phase["end_window_idx"]
            )
        ].copy()
        alpha_clusters = (
            phase_alpha.groupby(["window_id", "cluster_id"], as_index=False)
            .agg(
                n_alpha_sequences=("sequence_id", "nunique"),
                cluster_size=("cluster_size", "max"),
                cluster_n_datazones=("cluster_n_datazones", "max"),
                n_health_boards=("dz_health_board", "nunique"),
                n_local_authorities=("dz_local_authority", "nunique"),
            )
        )
        hb_counts = (
            phase_alpha["dz_health_board"].dropna().value_counts().head(5)
        )
        la_counts = (
            phase_alpha["dz_local_authority"].dropna().value_counts().head(5)
        )
        rows.append(
            {
                **phase,
                "n_sequences_all": int(phase_all["sequence_id"].nunique()),
                "n_alpha_sequences": int(phase_alpha["sequence_id"].nunique()),
                "alpha_sequence_fraction": (
                    phase_alpha["sequence_id"].nunique()
                    / phase_all["sequence_id"].nunique()
                    if phase_all["sequence_id"].nunique()
                    else np.nan
                ),
                "n_alpha_clusters": int(alpha_clusters["cluster_id"].nunique()),
                "median_alpha_cluster_size": alpha_clusters["cluster_size"].median(),
                "max_alpha_cluster_size": alpha_clusters["cluster_size"].max(),
                "median_alpha_cluster_datazones": alpha_clusters[
                    "cluster_n_datazones"
                ].median(),
                "max_alpha_cluster_datazones": alpha_clusters[
                    "cluster_n_datazones"
                ].max(),
                "n_health_boards": int(phase_alpha["dz_health_board"].nunique()),
                "n_local_authorities": int(phase_alpha["dz_local_authority"].nunique()),
                "top_health_boards": "; ".join(
                    f"{idx} ({val})" for idx, val in hb_counts.items()
                ),
                "top_local_authorities": "; ".join(
                    f"{idx} ({val})" for idx, val in la_counts.items()
                ),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "alpha_phase_summary.csv", index=False)
    return out


def summarise_alpha_emergence(seq: pd.DataFrame) -> None:
    alpha = seq[seq["is_alpha_pango"]].copy()
    cluster_emergence = (
        alpha.groupby(["window_id", "window_idx", "wn_mid_date", "cluster_id"], as_index=False)
        .agg(
            n_alpha_sequences=("sequence_id", "nunique"),
            cluster_size=("cluster_size", "max"),
            cluster_n_datazones=("cluster_n_datazones", "max"),
            n_health_boards=("dz_health_board", "nunique"),
            n_local_authorities=("dz_local_authority", "nunique"),
            health_boards=("dz_health_board", lambda s: "; ".join(sorted(set(s.dropna()))[:8])),
            local_authorities=("dz_local_authority", lambda s: "; ".join(sorted(set(s.dropna()))[:8])),
        )
        .sort_values(["window_idx", "cluster_size", "n_alpha_sequences"], ascending=[True, False, False])
    )
    cluster_emergence.to_csv(TABLE_DIR / "alpha_cluster_emergence.csv", index=False)

    hb_weekly = (
        alpha.groupby(["window_id", "window_idx", "wn_mid_date", "dz_health_board"], as_index=False)
        .agg(
            n_alpha_sequences=("sequence_id", "nunique"),
            n_alpha_clusters=("cluster_id", "nunique"),
            median_alpha_cluster_size=("cluster_size", "median"),
            max_alpha_cluster_size=("cluster_size", "max"),
        )
        .rename(columns={"dz_health_board": "health_board"})
        .sort_values(["window_idx", "health_board"])
    )
    hb_weekly.to_csv(TABLE_DIR / "alpha_health_board_weekly.csv", index=False)

    la_weekly = (
        alpha.groupby(["window_id", "window_idx", "wn_mid_date", "dz_local_authority"], as_index=False)
        .agg(
            n_alpha_sequences=("sequence_id", "nunique"),
            n_alpha_clusters=("cluster_id", "nunique"),
            median_alpha_cluster_size=("cluster_size", "median"),
            max_alpha_cluster_size=("cluster_size", "max"),
        )
        .rename(columns={"dz_local_authority": "local_authority"})
        .sort_values(["window_idx", "local_authority"])
    )
    la_weekly.to_csv(TABLE_DIR / "alpha_local_authority_weekly.csv", index=False)


def read_nextclade_mutation_flags(sequence_ids: set[str]) -> pd.DataFrame:
    """Read mutation flags from Nextclade in chunks and keep project sequences."""
    path = nextclade_tsv_path()
    chunks: list[pd.DataFrame] = []
    usecols = ["seqName", "aaSubstitutions"]
    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype=str,
        chunksize=100_000,
        low_memory=False,
    ):
        chunk = chunk[chunk["seqName"].isin(sequence_ids)].copy()
        if chunk.empty:
            continue
        chunk = chunk.reset_index(drop=True)
        aa = chunk["aaSubstitutions"].fillna("")
        out = pd.DataFrame({"sequence_id": chunk["seqName"].to_numpy()})
        for marker, slug in MUTATION_MARKERS.items():
            out[slug] = aa.str.contains(marker, regex=False).astype("int8").to_numpy()
        chunks.append(out)

    if not chunks:
        raise RuntimeError(
            f"No sequence IDs from the processed data matched {path}. "
            "Check sequence_id/seqName naming."
        )

    flags = pd.concat(chunks, ignore_index=True)
    agg = {slug: "max" for slug in MUTATION_MARKERS.values()}
    flags = flags.groupby("sequence_id", as_index=False).agg(agg)
    flags.to_parquet(CACHE_DIR / "nextclade_mutation_flags.parquet", index=False)
    return flags


def build_mutation_trajectories(seq: pd.DataFrame) -> pd.DataFrame:
    mapping_cols = [
        "sequence_id",
        "window_id",
        "window_idx",
        "wn_mid_date",
        "wn_positive_tests",
        "wn_prop_sequenced",
        "wn_no_sequences",
        "is_alpha_pango",
        "is_b1177_pango",
    ]
    mapping = seq[mapping_cols].drop_duplicates("sequence_id").copy()
    flags = read_nextclade_mutation_flags(set(mapping["sequence_id"]))
    joined = mapping.merge(flags, on="sequence_id", how="inner")

    grouped = (
        joined.groupby(["window_id", "window_idx", "wn_mid_date"], as_index=False)
        .agg(
            mutation_records=("sequence_id", "nunique"),
            wn_positive_tests=("wn_positive_tests", "max"),
            wn_prop_sequenced=("wn_prop_sequenced", "max"),
            wn_no_sequences=("wn_no_sequences", "max"),
            alpha_pango_sequences=("is_alpha_pango", "sum"),
            b1177_pango_sequences=("is_b1177_pango", "sum"),
            **{f"n_{slug}": (slug, "sum") for slug in MUTATION_MARKERS.values()},
        )
        .sort_values("window_idx")
    )

    for slug in MUTATION_MARKERS.values():
        grouped[f"freq_{slug}"] = grouped[f"n_{slug}"] / grouped["mutation_records"]

    alpha_cluster_counts = (
        seq[seq["is_alpha_pango"]]
        .groupby(["window_id"], as_index=False)
        .agg(alpha_pango_clusters=("cluster_id", "nunique"))
    )
    grouped = grouped.merge(alpha_cluster_counts, on="window_id", how="left")
    grouped["alpha_pango_clusters"] = grouped["alpha_pango_clusters"].fillna(0).astype(int)

    hb_occupancy = (
        seq.dropna(subset=["dz_health_board"])
        .groupby(["window_id", "dz_health_board"], as_index=False)
        .agg(hb_hospital_occupancy=("hb_hospital_occupancy", "max"))
        .groupby("window_id", as_index=False)
        .agg(
            hb_hospital_occupancy_total=("hb_hospital_occupancy", "sum"),
            hb_hospital_boards_with_data=("hb_hospital_occupancy", "count"),
        )
    )
    grouped = grouped.merge(hb_occupancy, on="window_id", how="left")
    grouped = attach_policy(grouped, "wn_mid_date")
    grouped.to_csv(TABLE_DIR / "alpha_mutation_trajectories.csv", index=False)
    return grouped


def glm_weight_vector(dat: pd.DataFrame, scheme: str) -> pd.Series:
    if scheme == "positive_test_weighted":
        weights = dat["wn_positive_tests"].replace(0, np.nan)
    elif scheme == "sequence_count_weighted":
        weights = dat["mutation_records"].replace(0, np.nan)
    elif scheme == "coverage_adjusted":
        coverage = dat["wn_prop_sequenced"].replace(0, np.nan)
        weights = dat["mutation_records"] / coverage
    elif scheme == "unweighted_weekly":
        weights = pd.Series(1.0, index=dat.index)
    else:
        raise ValueError(f"Unknown weight scheme: {scheme}")
    return weights.fillna(dat["mutation_records"]).astype(float).clip(lower=1.0)


def fit_growth_model(
    traj: pd.DataFrame,
    *,
    analysis: str,
    marker: str,
    marker_slug: str,
    period_code: str,
    weight_scheme: str,
) -> GrowthFit | None:
    periods = ordered_policy_table().set_index("period_code")
    period_row = periods.loc[period_code]
    start = pd.Timestamp(period_row["start_date"])
    end = pd.Timestamp(period_row["end_date"])
    dat = traj[
        traj["wn_mid_date"].between(start, end)
        & traj["mutation_records"].gt(0)
    ].copy()
    if dat.empty:
        return None

    success_col = f"n_{marker_slug}"
    dat["success"] = dat[success_col].astype(float)
    dat["total"] = dat["mutation_records"].astype(float)
    dat["prop"] = dat["success"] / dat["total"]
    dat["days"] = (dat["wn_mid_date"] - start).dt.days.astype(float)
    dat = dat.dropna(subset=["prop", "days"])
    if len(dat) < 4 or dat["prop"].nunique() < 2:
        return None

    x = sm.add_constant(dat[["days"]], has_constant="add")
    weights = glm_weight_vector(dat, weight_scheme)
    try:
        model = sm.GLM(
            dat["prop"].astype(float),
            x.astype(float),
            family=sm.families.Binomial(),
            freq_weights=weights,
        ).fit()
    except Exception as exc:
        warnings.warn(
            f"Growth model failed for {analysis}, {period_code}, {weight_scheme}: {exc}"
        )
        return None

    null = getattr(model, "null_deviance", np.nan)
    dev = getattr(model, "deviance", np.nan)
    pseudo_r2 = 1 - dev / null if null and np.isfinite(null) and null > 0 else np.nan
    conf = model.conf_int()
    return GrowthFit(
        analysis=analysis,
        marker=marker,
        marker_slug=marker_slug,
        period_code=period_code,
        period_label=period_row["period_label"],
        weight_scheme=weight_scheme,
        n_windows=int(len(dat)),
        n_success=int(dat["success"].sum()),
        n_total=int(dat["total"].sum()),
        start_date=dat["wn_mid_date"].min(),
        end_date=dat["wn_mid_date"].max(),
        origin_date=start,
        intercept=float(model.params["const"]),
        slope_per_day=float(model.params["days"]),
        intercept_se=float(model.bse["const"]),
        slope_se_per_day=float(model.bse["days"]),
        intercept_pvalue=float(model.pvalues["const"]),
        slope_pvalue=float(model.pvalues["days"]),
        slope_ci_low_per_day=float(conf.loc["days", 0]),
        slope_ci_high_per_day=float(conf.loc["days", 1]),
        aic=float(model.aic),
        pseudo_r2=float(pseudo_r2) if np.isfinite(pseudo_r2) else np.nan,
    )


def growth_fit_to_row(fit: GrowthFit) -> dict:
    slope_week = fit.slope_per_day * 7.0
    low_week = fit.slope_ci_low_per_day * 7.0
    high_week = fit.slope_ci_high_per_day * 7.0
    doubling_days = np.log(2) / fit.slope_per_day if fit.slope_per_day > 0 else np.nan
    halving_days = np.log(2) / abs(fit.slope_per_day) if fit.slope_per_day < 0 else np.nan
    return {
        **fit.__dict__,
        "start_date": fit.start_date.date().isoformat(),
        "end_date": fit.end_date.date().isoformat(),
        "origin_date": fit.origin_date.date().isoformat(),
        "slope_per_week": slope_week,
        "slope_ci_low_per_week": low_week,
        "slope_ci_high_per_week": high_week,
        "odds_ratio_per_week": np.exp(slope_week),
        "odds_ratio_ci_low_per_week": np.exp(low_week),
        "odds_ratio_ci_high_per_week": np.exp(high_week),
        "doubling_time_days_if_growth": doubling_days,
        "halving_time_days_if_decline": halving_days,
    }


def run_growth_models(traj: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    specs = [
        ("alpha_f5_n501y", "S:N501Y", "s_n501y", "F5"),
        ("alpha_l2_n501y", "S:N501Y", "s_n501y", "L2"),
        ("b1177_l2_a222v", "S:A222V", "s_a222v", "L2"),
    ]
    weight_schemes = [
        "positive_test_weighted",
        "sequence_count_weighted",
        "unweighted_weekly",
        "coverage_adjusted",
    ]
    fits: list[GrowthFit] = []
    for analysis, marker, slug, period_code in specs:
        for scheme in weight_schemes:
            fit = fit_growth_model(
                traj,
                analysis=analysis,
                marker=marker,
                marker_slug=slug,
                period_code=period_code,
                weight_scheme=scheme,
            )
            if fit is not None:
                fits.append(fit)

    rows = [growth_fit_to_row(fit) for fit in fits]
    all_fits = pd.DataFrame(rows)
    sensitivity = all_fits.copy()
    sensitivity.to_csv(TABLE_DIR / "alpha_growth_model_sensitivity.csv", index=False)

    primary = all_fits[all_fits["weight_scheme"] == "positive_test_weighted"].copy()
    primary.to_csv(TABLE_DIR / "alpha_growth_params.csv", index=False)
    return primary, sensitivity


def inv_logit(x: np.ndarray | float) -> np.ndarray | float:
    return 1 / (1 + np.exp(-np.asarray(x)))


def nearest_window(traj: pd.DataFrame, date: pd.Timestamp) -> pd.Series:
    idx = (traj["wn_mid_date"] - date).abs().idxmin()
    return traj.loc[idx]


def build_counterfactuals(
    traj: pd.DataFrame,
    growth_params: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = growth_params.set_index("analysis")
    if "alpha_f5_n501y" not in primary.index or "alpha_l2_n501y" not in primary.index:
        raise RuntimeError("Counterfactuals require primary F5 and L2 N501Y fits.")

    f5 = primary.loc["alpha_f5_n501y"]
    l2 = primary.loc["alpha_l2_n501y"]
    f5_origin = pd.Timestamp(f5["origin_date"])
    l2_slope = float(l2["slope_per_day"])

    observed_crossing = traj.loc[traj["freq_s_n501y"] >= 0.5, "wn_mid_date"].min()
    if pd.isna(observed_crossing):
        observed_crossing = pd.NaT

    scenarios = [
        ("actual_l2_start", "Actual L2 start", pd.Timestamp("2021-01-05")),
        ("expansion_date_2020_12_08", "Earlier L2 from 2020-12-08", pd.Timestamp("2020-12-08")),
        ("nearest_w021_2020_12_02", "Earlier L2 from 2020-12-02", pd.Timestamp("2020-12-02")),
        ("f5_start_2020_11_02", "L2 from F5 start", pd.Timestamp("2020-11-02")),
    ]

    projection_rows: list[dict] = []
    trajectory_rows: list[dict] = []
    grid = pd.date_range("2020-11-02", "2021-04-15", freq="D")

    actual_reach_date: pd.Timestamp | None = None
    for slug, label, switch_date in scenarios:
        near = nearest_window(traj, switch_date)
        switch_logit = float(f5["intercept"]) + float(f5["slope_per_day"]) * (
            switch_date - f5_origin
        ).days
        switch_freq = float(inv_logit(switch_logit))
        if l2_slope <= 0:
            reach_date = pd.NaT
        else:
            days_after_switch = (0 - switch_logit) / l2_slope
            reach_date = switch_date + pd.Timedelta(days=float(days_after_switch))
        if slug == "actual_l2_start":
            actual_reach_date = reach_date
        projection_rows.append(
            {
                "scenario": slug,
                "scenario_label": label,
                "requested_switch_date": switch_date.date().isoformat(),
                "nearest_observed_window_id": near["window_id"],
                "nearest_observed_window_mid_date": near["wn_mid_date"].date().isoformat(),
                "nearest_observed_n501y_frequency": near["freq_s_n501y"],
                "projected_switch_frequency": switch_freq,
                "projected_50pct_date": (
                    reach_date.date().isoformat() if pd.notna(reach_date) else ""
                ),
                "observed_50pct_window_mid_date": (
                    observed_crossing.date().isoformat()
                    if pd.notna(observed_crossing)
                    else ""
                ),
                "days_vs_actual_l2_projection": np.nan,
            }
        )

        for date in grid:
            if date <= switch_date:
                logit = float(f5["intercept"]) + float(f5["slope_per_day"]) * (
                    date - f5_origin
                ).days
            else:
                logit = switch_logit + l2_slope * (date - switch_date).days
            trajectory_rows.append(
                {
                    "date": date,
                    "scenario": slug,
                    "scenario_label": label,
                    "requested_switch_date": switch_date,
                    "projected_n501y_frequency": float(inv_logit(logit)),
                }
            )

    projections = pd.DataFrame(projection_rows)
    if actual_reach_date is not None and pd.notna(actual_reach_date):
        actual = pd.Timestamp(actual_reach_date).normalize()
        projections["days_vs_actual_l2_projection"] = projections["projected_50pct_date"].apply(
            lambda x: (pd.Timestamp(x).normalize() - actual).days if x else np.nan
        )

    trajectories = pd.DataFrame(trajectory_rows)
    projections.to_csv(TABLE_DIR / "alpha_counterfactual_projections.csv", index=False)
    trajectories.to_csv(TABLE_DIR / "alpha_counterfactual_trajectories.csv", index=False)
    return projections, trajectories


def main() -> None:
    setup_environment()

    print("Loading and summarising policy cluster table...")
    cluster = load_policy_cluster_table()
    summarise_periods(cluster)
    weekly = summarise_weekly(cluster)
    intensity_correlations(weekly)
    run_its(weekly)

    print("Loading sequence-level data for Alpha case study...")
    seq = load_sequence_table()
    summarise_alpha_phases(seq)
    summarise_alpha_emergence(seq)

    print("Reading Nextclade mutation flags and fitting Alpha/B.1.177 growth models...")
    traj = build_mutation_trajectories(seq)
    growth_params, _ = run_growth_models(traj)
    build_counterfactuals(traj, growth_params)

    print("Part 3 analysis complete.")
    print(f"Wrote tables to {TABLE_DIR}")
    print(f"Wrote caches to {CACHE_DIR}")


if __name__ == "__main__":
    main()
