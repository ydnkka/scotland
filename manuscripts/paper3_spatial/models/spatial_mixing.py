"""Per-cluster and pair-level spatial feature builders.

Distances use British National Grid Easting/Northing (metres). We divide by
1000 to report kilometres throughout.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Iterable

import numpy as np
import pandas as pd

from manuscripts.common import data


@lru_cache(maxsize=4)
def build_cluster_spatial_features(
    resolution: float = data.PRIMARY_RESOLUTION,
) -> pd.DataFrame:
    """One row per cluster with spatial summary stats.

    Columns
    -------
    n_sequences, n_datazones, centroid_x, centroid_y,
    mean_pairwise_km, max_pairwise_km, bbox_diag_km,
    simd_quintile_mode, who_voc, wn_mid_date
    """
    cols = [
        "window_id", "cluster_id", "resolution", "sequence_id",
        "dz_xcoord", "dz_ycoord", "datazone",
        "wn_mid_date", "who_voc", "dz_simd_quintile",
    ]
    df = data.load_analysis_columns(cols, resolution=resolution)
    df = df.dropna(subset=["dz_xcoord", "dz_ycoord"])

    def _mode(s: pd.Series):
        m = s.mode()
        return m.iloc[0] if len(m) > 0 else np.nan

    def _pair_summary(g: pd.DataFrame) -> pd.Series:
        xy = g[["dz_xcoord", "dz_ycoord"]].drop_duplicates().to_numpy()
        if xy.shape[0] <= 1:
            return pd.Series({
                "mean_pairwise_km": 0.0, "max_pairwise_km": 0.0, "bbox_diag_km": 0.0,
            })
        # Compute pairwise distances via a scipy helper if available; fall back to manual.
        try:
            from scipy.spatial.distance import pdist
            d = pdist(xy)
        except ImportError:  # pragma: no cover
            diff = xy[:, None, :] - xy[None, :, :]
            d = np.sqrt((diff ** 2).sum(-1))
            d = d[np.triu_indices_from(d, k=1)]
        bbox = np.sqrt(
            (xy[:, 0].max() - xy[:, 0].min()) ** 2
            + (xy[:, 1].max() - xy[:, 1].min()) ** 2
        )
        return pd.Series({
            "mean_pairwise_km": float(d.mean()) / 1000.0,
            "max_pairwise_km": float(d.max()) / 1000.0,
            "bbox_diag_km": float(bbox) / 1000.0,
        })

    grp = df.groupby(["window_id", "cluster_id"], observed=True)
    base = grp.agg(
        n_sequences=("sequence_id", "nunique"),
        n_datazones=("datazone", "nunique"),
        centroid_x=("dz_xcoord", "mean"),
        centroid_y=("dz_ycoord", "mean"),
        wn_mid_date=("wn_mid_date", "first"),
        who_voc=("who_voc", _mode),
        simd_quintile_mode=("dz_simd_quintile", _mode),
    ).reset_index()
    pair = grp.apply(_pair_summary, include_groups=False).reset_index()
    out = base.merge(pair, on=["window_id", "cluster_id"], how="left")
    return out


def sample_pairs_within_windows(
    n_per_window: int = 2000,
    resolution: float = data.PRIMARY_RESOLUTION,
    rng_seed: int = 42,
) -> pd.DataFrame:
    """Sample random pairs of sequences *within the same window* and record
    whether they fall in the same cluster, plus their BNG distance.

    Used for distance-decay figure (Fig. 3). We sample to keep pair counts
    manageable (full pairwise would be ~10^11 pairs for Omicron windows).
    """
    rng = np.random.default_rng(rng_seed)
    cols = [
        "window_id", "resolution", "cluster_id", "sequence_id",
        "dz_xcoord", "dz_ycoord", "who_voc", "wn_mid_date",
    ]
    df = data.load_analysis_columns(cols, resolution=resolution)
    df = df.dropna(subset=["dz_xcoord", "dz_ycoord"]).drop_duplicates(
        subset=["window_id", "sequence_id"]
    )
    rows = []
    for wid, g in df.groupby("window_id", observed=True):
        if len(g) < 2:
            continue
        idx = rng.integers(0, len(g), size=(n_per_window, 2))
        idx = idx[idx[:, 0] != idx[:, 1]]
        pair_a = g.iloc[idx[:, 0]].reset_index(drop=True)
        pair_b = g.iloc[idx[:, 1]].reset_index(drop=True)
        d_m = np.sqrt(
            (pair_a["dz_xcoord"] - pair_b["dz_xcoord"]) ** 2
            + (pair_a["dz_ycoord"] - pair_b["dz_ycoord"]) ** 2
        )
        rows.append(pd.DataFrame({
            "window_id": wid,
            "distance_km": d_m.to_numpy() / 1000.0,
            "co_cluster": (pair_a["cluster_id"].to_numpy() == pair_b["cluster_id"].to_numpy()),
            "who_voc": pair_a["who_voc"].to_numpy(),
            "wn_mid_date": pair_a["wn_mid_date"].to_numpy(),
        }))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
