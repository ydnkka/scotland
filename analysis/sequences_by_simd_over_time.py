"""Sequences per week by SIMD quintile, Scotland 2020-2023.

Two stacked panels share an x-axis:
    (A) Weekly sequenced case counts, stratified by SIMD quintile of the
        patient's data zone (1 = most deprived, 5 = least).
    (B) Proportion of positive tests sequenced in that week
        (`wn_prop_sequenced`), which is the critical surveillance-intensity
        covariate used in every downstream model.

Waves are shaded across both panels so readers can orient immediately.
"""

import polars as pl
import matplotlib.pyplot as plt
from utils import data, style


def _weekly_counts_by_simd(df: pl.DataFrame) -> pl.DataFrame:
    """Aggregate sequence counts per ISO week × SIMD quintile."""
    return (
        df.filter(pl.col("dz_simd_quintile").is_not_null())
        .with_columns(
            pl.col("collection_date")
            # Match the previous pandas W-SUN start_time labelling: Monday-start weeks.
            .dt.truncate("1w")
            .alias("week")
        )
        .group_by(["week", "dz_simd_quintile"])
        .agg(pl.len().alias("n"))
        .sort("week")
    )


def _wave_spans():
    return [(lbl, s, e) for lbl, s, e in data.WAVES]


def _shade_waves(ax):
    colors = ["#ffffff", "#f5e6e6", "#e6efdd", "#e6e6f5", "#f5efe0", "#e6f5f5"]
    for (lbl, s, e), c in zip(_wave_spans(), colors):
        ax.axvspan(s, e, color=c, alpha=0.55, zorder=0)
        ax.axvline(s, color="black", lw=0.5, ls="--", zorder=1)


def _wave_labels(ax):
    y = ax.get_ylim()[1]
    for lbl, s, e in _wave_spans():
        lbl = lbl.split('_')
        mid = s + (e - s) / 2
        ax.text(
            mid, y * 0.985,
            f"Wave {lbl[0][2]}\n{lbl[1]}\nn={int(lbl[2][1:]):,}",
            ha="center", va="top", fontsize=6.5, color="#555555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.65, pad=0.8),
        )


def main():
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = paths.root / "analysis/figures"

    # Sequence-level slice (deduplicated on sequence_id to avoid window-overlap double count).
    seq = (
        data.load_analysis_columns(
            ["sequence_id", "collection_date", "dz_simd_quintile"],
            resolution=data.PRIMARY_RESOLUTION, qc=None
        )
        .unique("sequence_id")
    )

    seq = _weekly_counts_by_simd(seq)

    # Window-level prop sequenced (one value per window).
    wn = (
        data.load_analysis_columns(
            ["window_id", "wn_mid_date", "wn_prop_sequenced"],
            resolution=data.PRIMARY_RESOLUTION, qc=None
        )
        .unique("window_id")
        .sort("wn_mid_date")
    )

    fig, (ax_a, ax_b) = style.new_figure(
        width="double", height_in=4.4, nrows=2, ncols=1, sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    # Panel A: stacked lines per SIMD quintile
    pivot = (
        seq.pivot(index="week", on="dz_simd_quintile", values="n", aggregate_function="sum")
        .sort("week")
    )
    # Ensure all quintile columns exist
    for q in [1, 2, 3, 4, 5]:
        col = str(q)
        if col not in pivot.columns:
            pivot = pivot.with_columns(pl.lit(0.0).alias(col))

    weeks = pivot["week"].to_list()

    _shade_waves(ax_a)
    for q in [1, 2, 3, 4, 5]:
        col = str(q)
        if col not in pivot.columns:
            continue
        ax_a.plot(
            weeks,
            pivot[col].fill_null(0).to_list(),
            color=style.SIMD_QUINTILE_PALETTE[q],
            lw=1.2,
            label=f"Q{q}" + (" (most deprived)" if q == 1 else (" (least)" if q == 5 else "")),
        )

    ax_a.set_ylim(0, max(0.3, seq["n"].max() * 1.2))
    ax_a.set_ylabel("Sequenced cases per week")
    ax_a.legend(
        title="Neighbourhood overall deprivation quintile", loc="upper center", ncol=5,
        columnspacing=0.8, handlelength=1.2, bbox_to_anchor=(0.5, 1.2), frameon=False,
    )
    _wave_labels(ax_a)

    # Panel B: wn_prop_sequenced over time
    mid_dates = wn["wn_mid_date"].to_list()
    prop_seq  = wn["wn_prop_sequenced"].to_list()

    _shade_waves(ax_b)
    ax_b.plot(mid_dates, prop_seq, color="#333333", lw=0.9)
    ax_b.fill_between(mid_dates, 0, prop_seq, color="#999999", alpha=0.25)
    ax_b.set_ylim(0, max(0.3, wn["wn_prop_sequenced"].max() * 1.05))
    ax_b.set_ylabel("Prop. sequenced")
    ax_b.set_xlabel("")

    for ax in (ax_a, ax_b):
        ax.margins(x=0.005)

    style.add_panel_labels([ax_a, ax_b], y=1.15)
    _ = style.save_figure(
        fig, out_dir / "sequences_by_simd_over_time",
        width="double", save_png=True, save_pdf=True,
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
