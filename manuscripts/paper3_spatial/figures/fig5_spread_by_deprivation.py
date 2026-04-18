"""Figure 5 — Hexbin of (cluster size, footprint) with deprivation colour overlay.

We plot log(1 + mean_pairwise_km) against log(1 + n_sequences). The hexbin
density is in greyscale; on top, per-hex mean modal SIMD quintile is shown
as a scatter colour overlay only for hexes with ≥ 10 clusters. This keeps
the density context while layering the deprivation signal.
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
    df = df[df["n_sequences"] >= 2]
    df = df.dropna(subset=["mean_pairwise_km", "simd_quintile_mode"])
    df["log_size"] = np.log1p(df["n_sequences"])
    df["log_km"] = np.log1p(df["mean_pairwise_km"])
    df["simd_quintile_mode"] = df["simd_quintile_mode"].astype(int)
    return df


def make_figure(df: pd.DataFrame) -> plt.Figure:
    fig, (ax_hex, ax_quint) = style.new_figure(
        width="double", height_in=3.6, nrows=1, ncols=2,
        gridspec_kw={"wspace": 0.18},
    )

    # Density hexbin (greyscale)
    hb = ax_hex.hexbin(
        df["log_size"], df["log_km"], gridsize=40,
        cmap="Greys", mincnt=1, linewidths=0.1,
    )
    cb = fig.colorbar(hb, ax=ax_hex, pad=0.02, shrink=0.8)
    cb.set_label("Clusters per hex", fontsize=7)
    cb.ax.tick_params(labelsize=7)
    ax_hex.set_xlabel("log(1 + cluster size)")
    ax_hex.set_ylabel("log(1 + mean pairwise km)")
    ax_hex.set_title("Size × footprint density (all clusters)", fontsize=9.0)

    # Mean deprivation per hex, overlaid on same axes
    gridsize = 40
    xedges = np.linspace(df["log_size"].min(), df["log_size"].max() + 1e-6, gridsize + 1)
    yedges = np.linspace(df["log_km"].min(), df["log_km"].max() + 1e-6, gridsize + 1)
    ix = np.clip(np.searchsorted(xedges, df["log_size"]) - 1, 0, gridsize - 1)
    iy = np.clip(np.searchsorted(yedges, df["log_km"]) - 1, 0, gridsize - 1)
    df = df.assign(ix=ix, iy=iy)
    hex_stats = (
        df.groupby(["ix", "iy"]).agg(
            mean_quint=("simd_quintile_mode", "mean"),
            n=("simd_quintile_mode", "size"),
        ).reset_index()
    )
    hex_stats = hex_stats[hex_stats["n"] >= 10]

    cx = 0.5 * (xedges[hex_stats["ix"]] + xedges[hex_stats["ix"] + 1])
    cy = 0.5 * (yedges[hex_stats["iy"]] + yedges[hex_stats["iy"] + 1])

    sc = ax_quint.scatter(
        cx, cy, c=hex_stats["mean_quint"], s=10 + hex_stats["n"] ** 0.5,
        cmap="RdBu", vmin=1, vmax=5, edgecolor="black", linewidth=0.2, alpha=0.9,
    )
    cb2 = fig.colorbar(sc, ax=ax_quint, pad=0.02, shrink=0.8)
    cb2.set_label("Mean modal SIMD quintile", fontsize=7)
    cb2.ax.tick_params(labelsize=7)
    ax_quint.set_xlabel("log(1 + cluster size)")
    ax_quint.set_ylabel("log(1 + mean pairwise km)")
    ax_quint.set_title("Deprivation gradient (≥10 clusters per bin)", fontsize=9.0)

    fig.suptitle("Where do deprived clusters sit in the size–footprint plane?",
                 x=0.02, ha="left", y=1.02, fontsize=9.5, fontweight="bold")
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper3_spatial/output"
    df = build_data()
    fig = make_figure(df)
    paths_out = style.save_figure(fig, out_dir / "fig5_spread_by_deprivation")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
