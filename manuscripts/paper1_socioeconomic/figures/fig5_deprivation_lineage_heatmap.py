"""Figure 5 — Median SIMD rank of clusters over time × top lineages.

A heat-map whose rows are the K most abundant PANGO lineages across the
study and columns are sliding analysis windows (`wn_mid_date`). Cell colour
encodes the median SIMD rank across sequences contributing to that
(lineage, window) cell, re-expressed as a rank percentile so values are
comparable over time despite small-N cells.

A paired bar along the right margin shows total sequence count per lineage
on a log axis.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from manuscripts.common import data, style

TOP_K_LINEAGES = 20


def build_matrix() -> tuple[pd.DataFrame, pd.Series]:
    df = data.load_analysis_columns(
        ["sequence_id", "pango_lineage", "wn_mid_date", "dz_simd_rank", "collection_date"],
        resolution=data.PRIMARY_RESOLUTION,
        qc=None,
    ).drop_duplicates("sequence_id")

    top = df["pango_lineage"].value_counts().head(TOP_K_LINEAGES).index.tolist()
    df = df[df["pango_lineage"].isin(top)]

    # Monthly bucketing is more readable than weekly windows.
    df["month"] = df["collection_date"].dt.to_period("M").dt.start_time
    piv = df.groupby(["pango_lineage", "month"])["dz_simd_rank"].median().reset_index()
    mat = piv.pivot(index="pango_lineage", columns="month", values="dz_simd_rank")
    # Re-order by cumulative abundance, most common on top. NB:
    # ``value_counts()`` is already descending; reindexing to ``mat.index``
    # here would silently re-sort it to pandas' default alphabetical order
    # (the bug we had before) — so take ``value_counts()`` directly and
    # push that order onto ``mat``.
    counts = df["pango_lineage"].value_counts()
    mat = mat.loc[counts.index]
    # Normalize to rank-percentile across the whole dataset (1 = least deprived).
    max_rank = df["dz_simd_rank"].max()
    mat_pct = mat / max_rank
    return mat_pct, counts


def make_figure(mat_pct: pd.DataFrame, counts: pd.Series) -> plt.Figure:
    # NB: we do NOT share the y-axis here. Sharing made ``ax_bar.set_yticks([])``
    # wipe ticks from the heatmap too, and ``ax_bar.invert_yaxis()`` flipped
    # the shared axis so the most-abundant lineage rendered at the bottom.
    # We instead align the two panels by manually syncing ``ylim``.
    fig, (ax_heat, ax_bar) = style.new_figure(
        width="double", height_in=5.0, nrows=1, ncols=2,
        gridspec_kw={"width_ratios": [5, 1], "wspace": 0.02},
    )

    # Heatmap: diverging palette centred on 0.5 (median SIMD rank across
    # Scotland). The classic RdBu midpoint (``#f7f7f7``) is near-white and
    # disappears into the figure background; swap it for a visible light
    # grey so "middle-deprivation" cells remain distinguishable from empty
    # (NaN) cells on a white canvas.
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "simd_div", ["#b2182b", "#ef8a62", "#bdbdbd", "#67a9cf", "#2166ac"]
    )
    im = ax_heat.imshow(
        mat_pct.values, aspect="auto", cmap=cmap,
        vmin=0.1, vmax=0.9, interpolation="nearest",
    )
    ax_heat.set_yticks(range(len(mat_pct.index)))
    ax_heat.set_yticklabels(mat_pct.index, fontsize=7)
    xticks_every = max(1, len(mat_pct.columns) // 10)
    ax_heat.set_xticks(range(0, len(mat_pct.columns), xticks_every))
    ax_heat.set_xticklabels(
        [d.strftime("%Y-%m") for d in mat_pct.columns[::xticks_every]],
        rotation=45, ha="right", fontsize=7,
    )
    ax_heat.set_xlabel("Collection month")

    # Colorbar sits above the heatmap; set the figure suptitle above the
    # colorbar so nothing overlaps the heatmap's title row.
    cbar = fig.colorbar(im, ax=ax_heat, pad=0.02, shrink=0.55, location="top")
    cbar.set_label("SIMD rank percentile  (0 = most deprived, 1 = least)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)
    fig.suptitle(
        "Median SIMD rank (percentile) of clusters by lineage over time",
        x=0.02, ha="left", y=1.02, fontsize=9.5, fontweight="bold",
    )

    # Right-margin: lineage abundance. ``counts`` is in the same row order as
    # ``mat_pct`` (most-abundant first), so plotting at y = 0..N-1 and copying
    # the heatmap's ylim puts each bar next to its matching row.
    ax_bar.barh(
        np.arange(len(counts)), counts.values,
        color="#555555", edgecolor="none", log=True,
    )
    ax_bar.set_xlabel("Sequences (log)")
    ax_bar.set_ylim(ax_heat.get_ylim())
    ax_bar.tick_params(axis="y", left=False, labelleft=False)

    # The ``location="top"`` colorbar steals vertical space from ``ax_heat``
    # but NOT ``ax_bar``, leaving ``ax_bar`` physically taller than the
    # heatmap. Matching ylim is not enough — equal data ranges on
    # different-height axes place row *i* and bar *i* at different figure-y
    # pixels, which is the visible misalignment. Force a draw so colorbar
    # layout has run, then re-anchor ``ax_bar`` to the heatmap's final bbox.
    fig.canvas.draw()
    pos_heat = ax_heat.get_position()
    pos_bar = ax_bar.get_position()
    ax_bar.set_position(
        [pos_bar.x0, pos_heat.y0, pos_bar.width, pos_heat.height]
    )

    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper1_socioeconomic/output"
    mat, counts = build_matrix()
    fig = make_figure(mat, counts)
    paths_out = style.save_figure(fig, out_dir / "fig5_deprivation_lineage_heatmap")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
