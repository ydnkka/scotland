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
from functools import lru_cache

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from manuscripts.common import data, style

TOP_K_LINEAGES = 20

@lru_cache(maxsize=1)
def _load_top_lineage_frame() -> tuple[pd.DataFrame, pd.Series, float]:
    """Load the sequence-level slice restricted to the top-K PANGO lineages.

    Returns the frame (with a ``month`` column added), the lineage-order
    ``counts`` series (descending), and the dataset-wide maximum SIMD rank
    used for percentile normalisation.
    """
    df = data.load_analysis_columns(
        ["sequence_id", "pango_lineage", "wn_mid_date", "dz_simd_rank", "collection_date"],
        resolution=data.PRIMARY_RESOLUTION,
        qc=None,
    ).drop_duplicates("sequence_id")

    top = df["pango_lineage"].value_counts().head(TOP_K_LINEAGES).index.tolist()
    df = df[df["pango_lineage"].isin(top)].copy()
    df["month"] = df["collection_date"].dt.to_period("M").dt.start_time
    counts = df["pango_lineage"].value_counts()  # descending
    max_rank = float(df["dz_simd_rank"].max())
    return df, counts, max_rank


def build_matrix() -> tuple[pd.DataFrame, pd.Series]:
    df, counts, max_rank = _load_top_lineage_frame()
    piv = df.groupby(["pango_lineage", "month"])["dz_simd_rank"].median().reset_index()
    mat = piv.pivot(index="pango_lineage", columns="month", values="dz_simd_rank")
    # Re-order by cumulative abundance, most untils on top. NB:
    # ``value_counts()`` is already descending; reindexing to ``mat.index``
    # here would silently re-sort it to pandas' default alphabetical order
    # (the bug we had before) — so take ``value_counts()`` directly and
    # push that order onto ``mat``.
    mat = mat.loc[counts.index]
    # Normalize to rank-percentile across the whole dataset (1 = least deprived).
    mat_pct = mat / max_rank
    return mat_pct, counts


def build_cell_long_table() -> pd.DataFrame:
    """Long-format table underlying the heatmap.

    One row per non-empty ``(lineage, month)`` cell. Columns:
    ``lineage``, ``month``, ``n_sequences``, ``median_simd_rank``,
    ``median_simd_rank_percentile``. Percentile is normalised by the
    dataset-wide maximum SIMD rank (same convention the heatmap uses).
    """
    df, counts, max_rank = _load_top_lineage_frame()
    agg = (
        df.groupby(["pango_lineage", "month"])
        .agg(
            n_sequences=("sequence_id", "nunique"),
            median_simd_rank=("dz_simd_rank", "median"),
        )
        .reset_index()
        .rename(columns={"pango_lineage": "lineage"})
    )
    agg["median_simd_rank_percentile"] = agg["median_simd_rank"] / max_rank
    # Keep lineage order consistent with the heatmap (most abundant first).
    order = {lin: i for i, lin in enumerate(counts.index)}
    agg["_ord"] = agg["lineage"].map(order)
    agg = agg.sort_values(["_ord", "month"]).drop(columns="_ord").reset_index(drop=True)
    return agg


def build_lineage_summary() -> pd.DataFrame:
    """Per-lineage summary across all sequences of the top-K lineages.

    One row per lineage. Columns: ``lineage``, ``n_sequences``,
    ``median_simd_rank_percentile``, ``q1_percentile``, ``q3_percentile``,
    ``first_month``, ``last_month``, ``n_months_present``. Rows are sorted
    by ``n_sequences`` descending to match heatmap row order.
    """
    df, counts, max_rank = _load_top_lineage_frame()
    df = df.dropna(subset=["dz_simd_rank"]).copy()
    df["pct"] = df["dz_simd_rank"] / max_rank

    rows = []
    for lineage in counts.index:
        sub = df[df["pango_lineage"] == lineage]
        if sub.empty:
            rows.append({
                "lineage": lineage, "n_sequences": 0,
                "median_simd_rank_percentile": np.nan,
                "q1_percentile": np.nan, "q3_percentile": np.nan,
                "first_month": pd.NaT, "last_month": pd.NaT,
                "n_months_present": 0,
            })
            continue
        months_present = sub["month"].dropna().unique()
        rows.append({
            "lineage": lineage,
            "n_sequences": int(sub["sequence_id"].nunique()),
            "median_simd_rank_percentile": float(sub["pct"].median()),
            "q1_percentile": float(sub["pct"].quantile(0.25)),
            "q3_percentile": float(sub["pct"].quantile(0.75)),
            "first_month": pd.Timestamp(months_present.min()) if len(months_present) else pd.NaT,
            "last_month": pd.Timestamp(months_present.max()) if len(months_present) else pd.NaT,
            "n_months_present": int(len(months_present)),
        })
    return pd.DataFrame(rows)


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
        "simd_div", list(style.SIMD_QUINTILE_PALETTE.values())
    )
    im = ax_heat.imshow(
        mat_pct.values, aspect="auto", cmap=cmap,
        vmin=0.1, vmax=0.9, interpolation="nearest",
    )
    ax_heat.set_yticks(range(len(mat_pct.index)))
    ax_heat.set_yticklabels(mat_pct.index)
    xticks_every = max(1, len(mat_pct.columns) // 10)
    ax_heat.set_xticks(range(0, len(mat_pct.columns), xticks_every))
    ax_heat.set_xticklabels(
        [d.strftime("%Y-%m") for d in mat_pct.columns[::xticks_every]],
        rotation=45, ha="right", fontsize=7,
    )
    ax_heat.set_xlabel("Collection month")

    # Colorbar sits above the heatmap; set the figure suptitle above the
    # colorbar so nothing overlaps the heatmap's title row.
    cbar = fig.colorbar(im, ax=ax_heat, pad=0.05, shrink=0.75, location="top")
    cbar.set_label("SIMD rank percentile  (0 = most deprived, 1 = least)")
    cbar.ax.tick_params(labelsize=7)
    # fig.suptitle(
    #     "Median SIMD rank (percentile) of clusters by lineage over time",
    #     x=0.02, ha="left", y=1.02, fontsize=9.5, fontweight="bold",
    # )

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


def main(out_dir: Path | str = None) -> dict[str, Path]:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir: Path = Path(out_dir) if out_dir else paths.root / "analysis/ch1_socioeconomic/figures"

    tables_dir = out_dir.parent / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    build_lineage_summary().to_csv(
        tables_dir / "fig5_lineage_summary.csv", index=False
    )
    build_cell_long_table().to_csv(
        tables_dir / "fig5_cells_long.csv", index=False
    )

    mat, counts = build_matrix()
    fig = make_figure(mat, counts)
    paths_out = style.save_figure(
        fig, out_dir / "fig5_deprivation_lineage_heatmap",
        width="double", save_png=True, save_pdf=True,
    )
    plt.close(fig)
    return paths_out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.figures)
    print(f"Wrote:\n   " + "\n   ".join(f"{k}: {v}" for k, v in p.items()))
