"""Figure 4 — Urban vs. rural cluster footprint.

Rural / urban is proxied from SIMD access-domain rank: clusters whose
modal access rank falls in the bottom 20% (highest access-deprivation =
least accessible = rural) are classed as *rural*; top 20% urban. Remaining
60% are labelled *mixed*.

Three-panel row:
    (A) Cluster size distribution (log1p)
    (B) Mean pairwise distance distribution (km, log1p)
    (C) Bounding-box diagonal (km, log1p)

Each panel contains three violins (rural/mixed/urban) with embedded medians.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, style
from manuscripts.paper3_spatial.models.spatial_mixing import build_cluster_spatial_features


def _attach_access_rank() -> pd.Series:
    """Compute modal access-domain rank per cluster, returned as a Series
    indexed by (window_id, cluster_id)."""
    df = data.load_analysis_columns(
        ["window_id", "cluster_id", "dz_simd_access_rank"],
        resolution=data.PRIMARY_RESOLUTION,
    )
    return (
        df.groupby(["window_id", "cluster_id"], observed=True)["dz_simd_access_rank"]
        .mean()
    )


def build_data() -> pd.DataFrame:
    feats = build_cluster_spatial_features()
    feats = feats[feats["n_sequences"] >= 2].copy()
    access = _attach_access_rank().rename("access_rank_mean")
    feats = feats.set_index(["window_id", "cluster_id"]).join(access).reset_index()
    feats = feats.dropna(subset=["access_rank_mean"])
    # Low access rank => most access-deprived (rural).
    q20, q80 = feats["access_rank_mean"].quantile([0.2, 0.8])
    labels = np.where(
        feats["access_rank_mean"] <= q20, "Rural",
        np.where(feats["access_rank_mean"] >= q80, "Urban", "Mixed"),
    )
    feats["locale"] = pd.Categorical(labels, categories=["Urban", "Mixed", "Rural"], ordered=True)
    return feats


def _violin(ax, values: list[np.ndarray], labels: list[str], colors: list[str]) -> None:
    parts = ax.violinplot(values, positions=range(len(values)), widths=0.75, showextrema=False)
    for body, c in zip(parts["bodies"], colors):
        body.set_facecolor(c); body.set_alpha(0.55); body.set_edgecolor("none")
    for i, v in enumerate(values):
        ax.scatter([i], [np.nanmedian(v)], zorder=3, color="black", s=8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)


def make_figure(df: pd.DataFrame) -> plt.Figure:
    fig, axes = style.new_figure(
        width="double", height_in=3.2, nrows=1, ncols=3,
        gridspec_kw={"wspace": 0.28},
    )
    colors = {"Urban": "#4e79a7", "Mixed": "#999999", "Rural": "#59a14f"}
    locales = ["Urban", "Mixed", "Rural"]

    metrics = [
        ("n_sequences",      "Cluster size  (log1p)"),
        ("mean_pairwise_km", "Mean pairwise km  (log1p)"),
        ("bbox_diag_km",     "Bounding-box diagonal km  (log1p)"),
    ]
    for ax, (col, label) in zip(axes, metrics):
        values = [np.log1p(df.loc[df["locale"] == loc, col].to_numpy()) for loc in locales]
        _violin(ax, values, locales, [colors[l] for l in locales])
        ax.set_ylabel(label)

    fig.suptitle("Urban vs. rural cluster footprint (rural/urban proxied from SIMD access rank)",
                 x=0.02, ha="left", y=1.02, fontsize=9.5, fontweight="bold")
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper3_spatial/output"
    df = build_data()
    fig = make_figure(df)
    paths_out = style.save_figure(fig, out_dir / "fig4_urban_rural_footprint")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
