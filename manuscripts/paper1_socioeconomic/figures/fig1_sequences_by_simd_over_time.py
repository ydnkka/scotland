"""Figure 1 — Sequences per week by SIMD quintile, Scotland 2020-2023.

Two stacked panels share an x-axis:
    (A) Weekly sequenced case counts, stratified by SIMD quintile of the
        patient's data zone (1 = most deprived, 5 = least).
    (B) Proportion of positive tests sequenced in that week
        (`wn_prop_sequenced`), which is the critical surveillance-intensity
        covariate used in every downstream model.

VOC epochs are shaded across both panels so readers can orient immediately.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from manuscripts.common import data, style


def _weekly_counts_by_simd(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate sequence counts per ISO week × SIMD quintile."""
    df = df.copy()
    df["week"] = df["collection_date"].dt.to_period("W-SUN").dt.start_time
    counts = (
        df.dropna(subset=["dz_simd_quintile"])
        .groupby(["week", "dz_simd_quintile"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    return counts


def _epoch_spans():
    from manuscripts.common.data import VOC_EPOCHS

    return [(lbl, pd.Timestamp(s), pd.Timestamp(e)) for lbl, s, e in VOC_EPOCHS]


def _shade_epochs(ax):
    colors = ["#ffffff", "#f5e6e6", "#e6efdd", "#e6e6f5", "#f5efe0"]
    for (lbl, s, e), c in zip(_epoch_spans(), colors):
        ax.axvspan(s, e, color=c, alpha=0.55, zorder=0)


def _epoch_labels(ax):
    y = ax.get_ylim()[1]
    for lbl, s, e in _epoch_spans():
        ax.text(
            s + (e - s) / 2, y * 0.985, lbl,
            ha="center", va="top", fontsize=6.5, color="#555555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.65, pad=0.8),
        )


def make_figure(df_seq: pd.DataFrame, df_prop: pd.DataFrame) -> plt.Figure:
    fig, (ax_a, ax_b) = style.new_figure(
        width="double", height_in=4.4, nrows=2, ncols=1, sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
    )

    # Panel A: stacked lines per SIMD quintile
    pivot = (
        df_seq.pivot_table(
            index="week", columns="dz_simd_quintile", values="n", aggfunc="sum"
        ).fillna(0).sort_index()
    )
    pivot = pivot.reindex(columns=[1, 2, 3, 4, 5])

    _shade_epochs(ax_a)
    for q in [1, 2, 3, 4, 5]:
        if q not in pivot.columns:
            continue
        ax_a.plot(
            pivot.index, pivot[q].values,
            color=style.SIMD_QUINTILE_PALETTE[q],
            lw=1.2, label=f"Q{q}" + (" (most deprived)" if q == 1 else (" (least)" if q == 5 else "")),
        )
    ax_a.set_ylabel("Sequenced cases per week")
    ax_a.legend(
        title="SIMD quintile", loc="upper center", ncol=5, columnspacing=0.8,
        handlelength=1.2, bbox_to_anchor=(0.5, 1.2), frameon=False,
    )
    _epoch_labels(ax_a)

    # Panel B: wn_prop_sequenced over time
    _shade_epochs(ax_b)
    ax_b.plot(
        df_prop["wn_mid_date"], df_prop["wn_prop_sequenced"],
        color="#333333", lw=0.9,
    )
    ax_b.fill_between(
        df_prop["wn_mid_date"], 0, df_prop["wn_prop_sequenced"],
        color="#999999", alpha=0.25,
    )
    ax_b.set_ylim(0, max(0.3, df_prop["wn_prop_sequenced"].max() * 1.05))
    ax_b.set_ylabel("Prop. sequenced")
    ax_b.set_xlabel("")

    for ax in (ax_a, ax_b):
        ax.margins(x=0.005)

    # fig.suptitle(
    #     "Weekly sequenced cases by SIMD quintile and surveillance intensity",
    #     x=0.08, ha="left", y=0.995, fontsize=9.5, fontweight="bold",
    # )
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper1_socioeconomic/output"

    # Sequence-level slice (deduplicated on sequence_id to avoid window-overlap double count).
    seq = data.load_analysis_columns(
        ["sequence_id", "collection_date", "dz_simd_quintile"],
        resolution=data.PRIMARY_RESOLUTION,
    ).drop_duplicates("sequence_id")

    # Window-level prop sequenced (one value per window).
    wn = data.load_analysis_columns(
        ["window_id", "wn_mid_date", "wn_prop_sequenced"],
        resolution=data.PRIMARY_RESOLUTION,
    ).drop_duplicates("window_id").sort_values("wn_mid_date")

    weekly = _weekly_counts_by_simd(seq)
    fig = make_figure(weekly, wn)
    paths_out = style.save_figure(fig, out_dir / "fig1_sequences_by_simd_over_time")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
