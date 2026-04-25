"""5 × 1 surveillance overview figure.

    (A) Weekly sequenced case counts (scaled to per-100 max)
    (B) Proportion of positive tests sequenced (surveillance intensity)
    (C) Proportion in non-singleton clusters – by age group
    (D) Proportion in non-singleton clusters – by sex
    (E) Proportion in non-singleton clusters – by SIMD deprivation quintile
"""

import polars as pl
import matplotlib.pyplot as plt
from utils import data, style


# ── helpers ──────────────────────────────────────────────────────────────────

def _cluster_trend_pl(df: pl.DataFrame, group_col: str) -> pl.DataFrame:
    return (
        df
        .group_by(["wn_mid_date", group_col])
        .agg(pl.col("in_non_singleton").mean())
        .sort("wn_mid_date")
        .pivot(on=group_col, index="wn_mid_date", values="in_non_singleton")
    )


def _weekly_counts(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.filter(pl.col("dz_simd_quintile").is_not_null())
        .with_columns(pl.col("collection_date").dt.truncate("1w").alias("week"))
        .group_by("week")
        .agg(pl.len().alias("n"))
        .sort("week")
    )


def _shade_waves(ax):
    colors = ["#ffffff", "#f5e6e6", "#e6efdd", "#e6e6f5", "#f5efe0", "#e6f5f5"]
    for (lbl, s, e), c in zip(data.WAVES, colors):
        ax.axvspan(s, e, color=c, alpha=0.55, zorder=0)
        ax.axvline(s, color="black", lw=0.5, ls="--", zorder=1)


def _wave_labels(ax):
    y = ax.get_ylim()[1]
    for lbl, s, e in data.WAVES:
        parts = lbl.split("_")
        mid = s + (e - s) / 2
        ax.text(
            mid, y * 0.985,
            f"{parts[1]}\nn={int(parts[2][1:]):,}",
            ha="center", va="top", color="#555555", fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.65, pad=0.8),
        )


def _plot_stratified(ax, trend_df: pl.DataFrame, palette: dict,
                     label_map: dict | None = None):
    dates = trend_df["wn_mid_date"].to_list()
    cols  = [c for c in trend_df.columns if c != "wn_mid_date"]
    _shade_waves(ax)
    for col in cols:
        label = label_map.get(col, col) if label_map else col
        ax.plot(
            dates,
            trend_df[col].to_list(),
            color=palette[col],
            lw=1.2, marker="o", markersize=2.5, alpha=0.9,
            label=label,
        )
    ax.set_ylim(0.2, 1.0)
    ax.set_ylabel("Prop. in cluster (n>1)")


def _outside_legend(ax, title):
    ax.legend(
        title=title, loc="center left",
        bbox_to_anchor=(1.01, 0.5), borderaxespad=0,
        frameon=True,
    )


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    style.set_theme(context="talk")
    paths   = data.Paths.from_config()
    out_dir = paths.root / "analysis/figures"

    # ── load data ────────────────────────────────────────────────────────────
    seq_raw = (
        data.load_analysis_columns(
            ["sequence_id", "collection_date", "dz_simd_quintile"],
            resolution=data.PRIMARY_RESOLUTION, qc=None,
        )
        .unique("sequence_id")
    )
    weekly_n = _weekly_counts(seq_raw)

    wn = (
        data.load_analysis_columns(
            ["window_id", "wn_mid_date", "wn_prop_sequenced"],
            resolution=data.PRIMARY_RESOLUTION, qc=None,
        )
        .unique("window_id")
        .sort("wn_mid_date")
    )

    cluster_df = (
        data.load_individual_features(format="long")
        .with_columns(
            pl.when(pl.col("is_female") == 1).then(pl.lit("Female"))
              .when(pl.col("is_female") == 0).then(pl.lit("Male"))
              .otherwise(pl.lit("Unknown"))
              .alias("sex")
        )
    )

    age_trend  = _cluster_trend_pl(cluster_df, "age_group")
    sex_trend  = _cluster_trend_pl(cluster_df, "sex")
    simd_trend = _cluster_trend_pl(cluster_df, "dz_simd_quintile")

    # ── figure: 5 rows × 1 col ───────────────────────────────────────────────
    fig, (ax_A, ax_B, ax_C, ax_D, ax_E) = style.new_figure(
        width="slide", height_in=10,
        nrows=5, ncols=1,
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2, 1.1, 1.1, 1.1, 1.1], "hspace": 0.08},
    )

    # ── A: weekly sequences (scaled to per-100 max) ───────────────────────────
    weeks  = weekly_n["week"].to_list()
    counts = (weekly_n["n"] / weekly_n["n"].max() * 100)

    _shade_waves(ax_A)
    ax_A.plot(weeks, counts, color="#2c7bb6", lw=1.3)
    ax_A.fill_between(weeks, 0, counts, color="#2c7bb6", alpha=0.18)
    ax_A.set_ylim(0, 125)
    ax_A.set_ylabel("Sequences ($\\times$100)")
    _wave_labels(ax_A)

    # ── B: proportion sequenced ───────────────────────────────────────────────
    mid_dates = wn["wn_mid_date"].to_list()
    prop_seq  = wn["wn_prop_sequenced"].to_list()

    _shade_waves(ax_B)
    ax_B.plot(mid_dates, prop_seq, color="#333333", lw=0.9)
    ax_B.fill_between(mid_dates, 0, prop_seq, color="#999999", alpha=0.25)
    ax_B.set_ylim(0, max(0.3, wn["wn_prop_sequenced"].max() * 1.05))
    ax_B.set_ylabel("Prop. sequenced")

    # ── C: cluster proportion by age group ───────────────────────────────────
    AGE_PALETTE = {
        "00-09": "#1f77b4", "10-19": "#ff7f0e", "20-39": "#2ca02c",
        "40–59": "#d62728", "60-74": "#9467bd", "elderly": "#8c564b",
    }
    _plot_stratified(ax_C, age_trend, AGE_PALETTE)
    _outside_legend(ax_C, "Age group")

    # ── D: cluster proportion by sex ─────────────────────────────────────────
    SEX_PALETTE = {"Male": "#1565C0", "Female": "#E91E63"}
    _plot_stratified(ax_D, sex_trend, SEX_PALETTE)
    _outside_legend(ax_D, "Sex")

    # ── E: cluster proportion by SIMD quintile ───────────────────────────────
    SIMD_PALETTE = {str(q): style.SIMD_QUINTILE_PALETTE[q] for q in range(1, 6)}
    SIMD_LABELMAP = {
        "1": "Q1 (most deprived)", "2": "Q2", "3": "Q3", "4": "Q4",
        "5": "Q5 (least deprived)",
    }
    _plot_stratified(ax_E, simd_trend, SIMD_PALETTE, label_map=SIMD_LABELMAP)
    _outside_legend(ax_E, "SIMD quintile")
    fig.supxlabel("Date (from 2020-2023)")

    # ── shared formatting ─────────────────────────────────────────────────────
    for ax in (ax_A, ax_B, ax_C, ax_D, ax_E):
        ax.margins(x=0.005)
        ax.set_facecolor("#fafafa")

    # style.add_panel_labels([ax_A, ax_B, ax_C, ax_D, ax_E], y=1.02)

    _ = style.save_figure(
        fig, out_dir / "surveillance_overview_5x1",
        width="slide", save_png=True, save_pdf=True,
    )
    plt.close(fig)


if __name__ == "__main__":
    main()