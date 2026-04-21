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
import numpy as np
import pandas as pd

from analysis.untils import data, style


def build_panel_a_summary(df_seq: pd.DataFrame) -> pd.DataFrame:
    """Sequences per (epoch, SIMD quintile) with within-epoch / within-quintile shares.

    Uses the same sequence-level slice that feeds Panel A. Quintiles are labelled
    1 (most deprived) to 5 (least). Epoch is assigned from ``collection_date``.
    """
    df = df_seq.dropna(subset=["dz_simd_quintile"]).copy()
    df["dz_simd_quintile"] = df["dz_simd_quintile"].astype(int)
    df["epoch"] = data.assign_epoch(df["collection_date"])
    df = df.dropna(subset=["epoch"])

    tab = (
        df.groupby(["epoch", "dz_simd_quintile"], observed=True)
        .size()
        .rename("n_sequences")
        .reset_index()
        .rename(columns={"dz_simd_quintile": "simd_quintile"})
    )

    # Fill missing (epoch, quintile) combos with 0 so readers see gaps explicitly.
    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS]
    full = pd.MultiIndex.from_product(
        [epochs, [1, 2, 3, 4, 5]], names=["epoch", "simd_quintile"]
    )
    tab = (
        tab.set_index(["epoch", "simd_quintile"])
        .reindex(full, fill_value=0)
        .reset_index()
    )

    epoch_totals = tab.groupby("epoch")["n_sequences"].transform("sum")
    quintile_totals = tab.groupby("simd_quintile")["n_sequences"].transform("sum")
    tab["pct_within_epoch"] = np.where(
        epoch_totals > 0, 100.0 * tab["n_sequences"] / epoch_totals, np.nan
    )
    tab["pct_within_quintile"] = np.where(
        quintile_totals > 0, 100.0 * tab["n_sequences"] / quintile_totals, np.nan
    )
    return tab


def build_panel_b_summary(df_prop: pd.DataFrame) -> pd.DataFrame:
    """Per-epoch summary of surveillance intensity (``wn_prop_sequenced``).

    One row per epoch. ``n_weeks`` is the number of windows contributing to the
    summary (one value per window).
    """
    df = df_prop.dropna(subset=["wn_prop_sequenced"]).copy()
    df["epoch"] = data.assign_epoch(df["wn_mid_date"])
    df = df.dropna(subset=["epoch"])

    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS]
    rows = []
    for epoch in epochs:
        v = df.loc[df["epoch"] == epoch, "wn_prop_sequenced"].to_numpy()
        if v.size == 0:
            rows.append({
                "epoch": epoch, "n_weeks": 0,
                "median": np.nan, "q1": np.nan, "q3": np.nan,
                "min": np.nan, "max": np.nan, "mean": np.nan,
            })
            continue
        rows.append({
            "epoch": epoch,
            "n_weeks": int(v.size),
            "median": float(np.median(v)),
            "q1": float(np.quantile(v, 0.25)),
            "q3": float(np.quantile(v, 0.75)),
            "min": float(v.min()),
            "max": float(v.max()),
            "mean": float(v.mean()),
        })
    return pd.DataFrame(rows)


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


def make_figure(df_seq: pd.DataFrame, df_prop: pd.DataFrame) -> plt.Figure:
    fig, (ax_a, ax_b) = style.new_figure(
        width="double", height_in=4.4, nrows=2, ncols=1, sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1]},
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
    ax_a.set_ylim(0, max(0.3, df_seq["n"].max() * 1.1))
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

    style.add_panel_labels([ax_a, ax_b], y=1.15)
    return fig


def main(out_dir: Path | str = None) -> dict[str, Path]:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir: Path = Path(out_dir) if out_dir else paths.root / "analysis/ch1_socioeconomic/figures"

    # Sequence-level slice (deduplicated on sequence_id to avoid window-overlap double count).
    seq = data.load_analysis_columns(
        ["sequence_id", "collection_date", "dz_simd_quintile"],
        resolution=data.PRIMARY_RESOLUTION, qc=None
    ).drop_duplicates("sequence_id")

    # Window-level prop sequenced (one value per window).
    wn = data.load_analysis_columns(
        ["window_id", "wn_mid_date", "wn_prop_sequenced"],
        resolution=data.PRIMARY_RESOLUTION, qc=None
    ).drop_duplicates("window_id").sort_values("wn_mid_date")

    tables_dir = out_dir.parent / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    build_panel_a_summary(seq).to_csv(
        tables_dir / "fig1_sequences_by_simd.csv", index=False
    )
    build_panel_b_summary(wn).to_csv(
        tables_dir / "fig1_prop_sequenced_by_epoch.csv", index=False
    )

    weekly = _weekly_counts_by_simd(seq)
    fig = make_figure(weekly, wn)
    paths_out = style.save_figure(
        fig, out_dir / "fig1_sequences_by_simd_over_time",
        width="double", save_png=True, save_pdf=True
    )
    plt.close(fig)
    return paths_out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.figures)
    print(f"Wrote:\n   " + "\n   ".join(f"{k}: {v}" for k, v in p.items()))
