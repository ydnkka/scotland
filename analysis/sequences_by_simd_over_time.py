
"""Sequences per week by SIMD quintile, Scotland 2020-2023.

Two stacked panels share an x-axis:
    (A) Weekly sequenced case counts, stratified by SIMD quintile of the
        patient's data zone (1 = most deprived, 5 = least).
    (B) Proportion of positive tests sequenced in that week
        (`wn_prop_sequenced`), which is the critical surveillance-intensity
        covariate used in every downstream model.

VOC epochs are shaded across both panels so readers can orient immediately.
"""

import pandas as pd

import matplotlib.pyplot as plt
from utils import data, style

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
    return [
        (lbl, pd.Timestamp(s), pd.Timestamp(e))
        for lbl, s, e in data.VOC_EPOCHS
    ]


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


def main():
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = paths.root / "analysis/figures"

    # Sequence-level slice (deduplicated on sequence_id to avoid window-overlap double count).
    seq = data.load_analysis_columns(
        ["sequence_id", "collection_date", "dz_simd_quintile"],
        resolution=data.PRIMARY_RESOLUTION, qc=None
    ).drop_duplicates("sequence_id")

    seq = _weekly_counts_by_simd(seq)

    # Window-level prop sequenced (one value per window).
    wn = data.load_analysis_columns(
        ["window_id", "wn_mid_date", "wn_prop_sequenced"],
        resolution=data.PRIMARY_RESOLUTION, qc=None
    ).drop_duplicates("window_id").sort_values("wn_mid_date")

    fig, (ax_a, ax_b) = style.new_figure(
        width="double", height_in=4.4, nrows=2, ncols=1, sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    # Panel A: stacked lines per SIMD quintile
    pivot = (
        seq.pivot_table(
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
    ax_a.set_ylim(0, max(0.3, seq["n"].max() * 1.1))
    ax_a.set_ylabel("Sequenced cases per week")
    ax_a.legend(
        title="Neighbourhood overall deprivation quintile", loc="upper center", ncol=5, columnspacing=0.8,
        handlelength=1.2, bbox_to_anchor=(0.5, 1.2), frameon=False,
    )
    _epoch_labels(ax_a)

    # Panel B: wn_prop_sequenced over time
    _shade_epochs(ax_b)
    ax_b.plot(
        wn["wn_mid_date"], wn["wn_prop_sequenced"],
        color="#333333", lw=0.9,
    )
    ax_b.fill_between(
        wn["wn_mid_date"], 0, wn["wn_prop_sequenced"],
        color="#999999", alpha=0.25,
    )
    ax_b.set_ylim(0, max(0.3, wn["wn_prop_sequenced"].max() * 1.05))
    ax_b.set_ylabel("Prop. sequenced")
    ax_b.set_xlabel("")

    for ax in (ax_a, ax_b):
        ax.margins(x=0.005)

    style.add_panel_labels([ax_a, ax_b], y=1.15)
    _ = style.save_figure(
        fig, out_dir / "sequences_by_simd_over_time",
        width="double", save_png=True, save_pdf=True
    )
    plt.close(fig)


if __name__ == "__main__":
    main()