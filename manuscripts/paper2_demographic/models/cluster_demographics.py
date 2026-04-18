"""Cluster-level regression frames for the demographic paper.

Joins the pre-computed `cluster_demographic_features.parquet` (from
`analysis/demographic.py`) to the window-level `wn_prop_sequenced` and
`who_voc` that are needed as adjustment variables but do not live in the
cluster-level parquet.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from manuscripts.common import data, stats


@lru_cache(maxsize=4)
def build_cluster_regression_frame(resolution: float = data.PRIMARY_RESOLUTION) -> pd.DataFrame:
    feat = data.load_cluster_demographic_features()
    feat = feat[feat["resolution"] == resolution].copy()

    # Window-level attributes (prop sequenced, dominant VOC per window) from the master parquet.
    wn_cols = ["window_id", "wn_prop_sequenced", "who_voc", "sequence_id"]
    wn = data.load_analysis_columns(wn_cols, resolution=resolution)

    def _mode(s: pd.Series):
        m = s.mode()
        return m.iloc[0] if len(m) > 0 else np.nan

    window_attrs = (
        wn.groupby("window_id", observed=True)
        .agg(
            wn_prop_sequenced=("wn_prop_sequenced", "first"),
            who_voc=("who_voc", _mode),
            wn_n_sequences=("sequence_id", "nunique"),
        )
        .reset_index()
    )
    out = feat.merge(window_attrs, on="window_id", how="left")
    return out


def cluster_size_model(frame: pd.DataFrame):
    """NB GLM: n_sequences ~ demographics + VOC (offset=prop_seq)."""
    df = frame.dropna(
        subset=["n_sequences", "median_age", "frac_female", "frac_vaccinated",
                "mean_vacc_dose", "wn_prop_sequenced"]
    ).copy()
    voc = stats.one_hot(df["who_voc"].fillna("None"), reference="Omicron", prefix="voc")
    preds_cont = df[["median_age", "frac_female", "frac_vaccinated", "mean_vacc_dose"]]
    # Standardise continuous predictors for comparability.
    preds_cont_z = (preds_cont - preds_cont.mean()) / preds_cont.std()
    X = pd.concat([preds_cont_z, voc], axis=1)
    df_model = pd.concat([df[["n_sequences", "wn_prop_sequenced"]], X], axis=1)
    return stats.negbin_cluster_size(
        df_model, outcome="n_sequences",
        predictors=X.columns.tolist(), offset="wn_prop_sequenced",
    )


def singleton_model(frame: pd.DataFrame, *, min_prop_seq: float = 1e-3):
    """Logistic GLM for `is_singleton`.

    Notes
    -----
    Windows with essentially no sequencing (`wn_prop_sequenced < min_prop_seq`)
    are dropped before fitting; they are structurally uninformative and produce
    log-offset values so extreme that the IRLS inner sigmoid overflows. The
    surviving log-prop-sequenced column is standardised to put every predictor
    on the same ~O(1) scale, which is what keeps the `exp()` step in the
    logistic link numerically safe.
    """
    df = frame.dropna(
        subset=["is_singleton", "median_age", "frac_female", "frac_vaccinated",
                "mean_vacc_dose", "wn_prop_sequenced"]
    ).copy()
    df = df[df["wn_prop_sequenced"] >= min_prop_seq]
    voc = stats.one_hot(df["who_voc"].fillna("None"), reference="Omicron", prefix="voc")
    preds_cont = df[["median_age", "frac_female", "frac_vaccinated", "mean_vacc_dose"]]
    preds_cont_z = (preds_cont - preds_cont.mean()) / preds_cont.std()
    log_prop = np.log(df["wn_prop_sequenced"].astype(float))
    log_prop_z = (log_prop - log_prop.mean()) / log_prop.std()
    X = pd.concat([preds_cont_z, voc], axis=1)
    X["log_prop_seq_z"] = log_prop_z.values
    df_model = pd.concat([df[["is_singleton"]].reset_index(drop=True),
                          X.reset_index(drop=True)], axis=1)
    return stats.logit_singleton(df_model, predictors=X.columns.tolist())
