"""Figure 2 — Within-cluster geographic spread (mean pairwise km), by SIMD quintile × epoch.

Twin panels:
    (A) Ridgeline of log(mean_pairwise_km + 1) per SIMD quintile.
    (B) Small boxplot matrix: rows = VOC epochs, columns = SIMD quintiles.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, style
from manuscripts.paper3_spatial.models.spatial_mixing import build_cluster_spatial_features


def build_data() -> pd.DataFrame:
    df = build_cluster_spatial_features()
    df = df[df["n_sequences"] >= 2]  # only multi-member clusters
    df["epoch"] = data.assign_epoch(df["wn_mid_date"])
    df = df.dropna(subset=["simd_quintile_mode", "epoch", "mean_pairwise_km"])
    df["simd_quintile_mode"] = df["simd_quintile_mode"].astype(int)
    return df


def _kde(x: np.ndarray, grid: np.ndarray) -> np.ndarray:
    from scipy.stats import gaussian_kde
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return np.zeros_like(grid)
    try:
        return gaussian_kde(x, bw_method="silverman")(grid)
    except Exception:
        return np.zeros_like(grid)


def make_figure(df: pd.DataFrame) -> plt.Figure:
    fig, (ax_a, ax_b) = style.new_figure(
        width="double", height_in=4.0, nrows=1, ncols=2,
        gridspec_kw={"width_ratios": [3, 4], "wspace": 0.22},
    )

    # Panel A: ridgeline of log(mean_pairwise_km + 1) per quintile
    grid = np.linspace(0, np.log1p(df["mean_pairwise_km"].quantile(0.995)), 400)
    for i, q in enumerate([1, 2, 3, 4, 5]):
        vals = np.log1p(df.loc[df["simd_quintile_mode"] == q, "mean_pairwise_km"].to_numpy())
        dens = _kde(vals, grid)
        dens /= (dens.max() + 1e-9)
        y = i * 1.0
        c = style.SIMD_QUINTILE_PALETTE[q]
        ax_a.fill_between(grid, y, y + dens * 0.9, color=c, alpha=0.45, lw=0)
        ax_a.plot(grid, y + dens * 0.9, color=c, lw=1.0)
        ax_a.text(-0.1, y + 0.15, f"Q{q}", va="center", ha="right", fontsize=8)
    ax_a.set_yticks([])
    ax_a.set_xlabel("log(1 + mean pairwise distance, km)")
    ax_a.set_title("Spread distribution by SIMD quintile", fontsize=9.0)

    # Panel B: mini boxplot matrix
    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS if lbl in set(df["epoch"].astype(str).unique())]
    quintiles = [1, 2, 3, 4, 5]
    width = 0.8 / len(quintiles)
    for i, epoch in enumerate(epochs):
        for j, q in enumerate(quintiles):
            sub = df[(df["epoch"] == epoch) & (df["simd_quintile_mode"] == q)]
            v = np.log1p(sub["mean_pairwise_km"].to_numpy())
            if len(v) < 3:
                continue
            pos = i + (j - 2) * width
            parts = ax_b.boxplot(
                v, positions=[pos], widths=width * 0.85, patch_artist=True,
                showfliers=False, showcaps=False,
                boxprops=dict(facecolor=style.SIMD_QUINTILE_PALETTE[q],
                              edgecolor="black", lw=0.3),
                medianprops=dict(color="black", lw=0.8),
                whiskerprops=dict(color="black", lw=0.5),
            )
    ax_b.set_xticks(range(len(epochs)))
    ax_b.set_xticklabels(epochs, rotation=20, ha="right", fontsize=7.5)
    ax_b.set_ylabel("log(1 + mean pairwise km)")
    ax_b.set_title("By VOC epoch × SIMD quintile", fontsize=9.0)

    # Custom legend for quintile colours
    from matplotlib.patches import Patch
    handles = [Patch(color=style.SIMD_QUINTILE_PALETTE[q], label=f"Q{q}") for q in quintiles]
    ax_b.legend(handles=handles, title="SIMD quintile", loc="upper right",
                ncol=5, fontsize=6.5, handlelength=0.9, frameon=False)

    fig.suptitle("Within-cluster geographic spread across deprivation and VOC",
                 x=0.02, ha="left", y=1.01, fontsize=9.5, fontweight="bold")
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper3_spatial/output"
    df = build_data()
    fig = make_figure(df)
    paths_out = style.save_figure(fig, out_dir / "fig2_within_cluster_spread")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
