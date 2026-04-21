"""Figure 1 — Map of Scotland showing cluster centroids by VOC.

One point per cluster, placed at the mean BNG easting/northing of its
members' data zones. Point size is proportional to cluster size
(log scaled); colour encodes the dominant WHO VOC. Small multiples by
VOC make each layer legible.

If the Scottish Government DZ boundary shapefile is present at
`data/raw/datazone/sg_datazone_bdry_2011.shp`, we overlay a simplified
dissolved coastline; otherwise we draw axes only with BNG ticks.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, style
from manuscripts.paper3_spatial.models.spatial_mixing import build_cluster_spatial_features

KEEP_VOCS = ("Alpha", "Delta", "Omicron")


def _scotland_outline(ax, shp_path: Path | None) -> None:
    """Draw a simplified Scotland outline if shapefile is available."""
    if shp_path is None or not shp_path.exists():
        # Rough bounding box of Scotland in BNG (metres).
        ax.set_xlim(0, 5e5)
        ax.set_ylim(5.3e5, 1.25e6)
        ax.set_aspect("equal")
        return
    try:
        import geopandas as gpd
        gdf = gpd.read_file(shp_path)
        gdf = gdf.to_crs(epsg=27700)
        outline = gdf.dissolve().boundary
        outline.plot(ax=ax, color="#555555", linewidth=0.4)
    except Exception as e:  # pragma: no cover
        print(f"  (skipping coastline overlay: {e})")
    ax.set_aspect("equal")


def make_figure(cluster_spatial: pd.DataFrame, shp_path: Path | None) -> plt.Figure:
    df = cluster_spatial.copy()
    df = df[df["who_voc"].isin(KEEP_VOCS)]
    fig, axes = style.new_figure(
        width="double", height_in=4.8, nrows=1, ncols=len(KEEP_VOCS),
        sharex=True, sharey=True, gridspec_kw={"wspace": 0.04},
    )
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    for ax, voc in zip(axes, KEEP_VOCS):
        sub = df[df["who_voc"] == voc]
        _scotland_outline(ax, shp_path)
        sizes = 6 + np.log1p(sub["n_sequences"]) * 6
        ax.scatter(
            sub["centroid_x"], sub["centroid_y"],
            s=sizes,
            c=style.WHO_VOC_PALETTE.get(voc, "#333333"),
            alpha=0.40, edgecolor="none", linewidth=0,
        )
        ax.set_title(f"{voc}\n{len(sub):,} clusters", fontsize=8.8)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle("Geographic distribution of SARS-CoV-2 transmission cluster centroids by VOC",
                 x=0.02, ha="left", y=1.01, fontsize=9.5, fontweight="bold")
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper3_spatial/figures"
    shp = paths.root / "data/raw/datazone/sg_datazone_bdry_2011.shp"
    cluster_spatial = build_cluster_spatial_features()
    fig = make_figure(cluster_spatial, shp if shp.exists() else None)
    paths_out = style.save_figure(fig, out_dir / "fig1_cluster_centroid_map")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
