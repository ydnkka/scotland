"""Figure 2 — Cluster size distribution by SIMD quintile, faceted by VOC epoch.

Cluster-level dataset: one row per `(window_id, cluster_id)`. Cluster SIMD
quintile is defined as the modal quintile among its members
(`simd_quintile_mode`). Cluster size is plotted on a log1p axis because the
distribution is heavy-tailed (modal size = 1; maxima reach thousands).

Each panel shows the five SIMD quintiles side-by-side; a notched box plot is
overlaid on a violin to convey both the density and the median/IQR. A
Kruskal-Wallis p-value annotates each panel.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis.utils import data, style


def _panel_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["epoch"] = data.assign_epoch(df["wn_mid_date"])
    df = df.dropna(subset=["simd_quintile_mode", "epoch"])
    df["simd_quintile_mode"] = df["simd_quintile_mode"].astype(int)
    return df


def _kw_pvalue(group_values: list[np.ndarray]) -> float:
    from scipy.stats import kruskal

    try:
        return float(kruskal(*group_values).pvalue)
    except ValueError:
        return float("nan")


def build_summary_table(cluster_simd: pd.DataFrame) -> pd.DataFrame:
    """Per-(epoch, SIMD quintile) summary statistics for cluster size.

    One row per (epoch, quintile). The Kruskal-Wallis p-value is computed once
    per epoch across the five quintiles and repeated across that epoch's rows
    for convenience.
    """
    df = _panel_data(cluster_simd)
    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS]

    rows = []
    for epoch in epochs:
        sub = df[df["epoch"] == epoch]
        values_by_q = {
            q: sub.loc[sub["simd_quintile_mode"] == q, "n_sequences"].to_numpy()
            for q in range(1, 6)
        }
        p = _kw_pvalue([v for v in values_by_q.values() if len(v) > 1])
        n_epoch = sum(len(v) for v in values_by_q.values())
        for q, v in values_by_q.items():
            if len(v) == 0:
                rows.append({
                    "epoch": epoch, "simd_quintile": q,
                    "n_clusters": 0, "n_sequences": 0,
                    "median": np.nan, "q1": np.nan, "q3": np.nan,
                    "mean": np.nan, "min": np.nan, "max": np.nan,
                    "kw_pvalue": p, "n_epoch": n_epoch,
                })
                continue
            rows.append({
                "epoch": epoch,
                "simd_quintile": q,
                "n_clusters": int(v.size),
                "n_sequences": int(v.sum()),
                "median": float(np.median(v)),
                "q1": float(np.quantile(v, 0.25)),
                "q3": float(np.quantile(v, 0.75)),
                "mean": float(v.mean()),
                "variances": float(v.var()),
                "stddev": float(v.std()),
                "min": int(v.min()),
                "max": int(v.max()),
                "kw_pvalue": p,
                "n_epoch": n_epoch,
            })
    return pd.DataFrame(rows)


def _violin_box(ax, values_by_q: dict[int, np.ndarray]) -> None:
    positions = list(range(1, 6))
    data_list = [np.log1p(values_by_q.get(q, np.array([], dtype=float))) for q in positions]
    if all(len(d) == 0 for d in data_list):
        ax.set_axis_off()
        return

    parts = ax.violinplot(
        data_list, positions=positions, widths=0.8, showextrema=False, showmedians=False,
    )
    for body, q in zip(parts["bodies"], positions):
        body.set_facecolor(style.SIMD_QUINTILE_PALETTE[q])
        body.set_edgecolor("none")
        body.set_alpha(0.55)

    ax.boxplot(
        data_list, positions=positions, widths=0.18, notch=True, showcaps=False,
        patch_artist=True, medianprops=dict(color="black", lw=0.9),
        boxprops=dict(facecolor="white", edgecolor="black", lw=0.5),
        whiskerprops=dict(color="black", lw=0.5),
        flierprops=dict(marker="", linestyle=""),
    )
    ax.set_xticks(positions)
    ax.set_xticklabels([f"Q{q}" for q in positions])


def make_figure(cluster_simd: pd.DataFrame) -> plt.Figure:
    df = _panel_data(cluster_simd)
    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS]

    fig, axes = style.new_figure(
        width="double", height_in=3.6, nrows=1, ncols=len(epochs),
        sharey=True, gridspec_kw={"wspace": 0.1},
    )
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])

    for ax, epoch in zip(axes, epochs):
        sub = df[df["epoch"] == epoch]
        values_by_q = {
            q: sub.loc[sub["simd_quintile_mode"] == q, "n_sequences"].to_numpy()
            for q in range(1, 6)
        }
        _violin_box(ax, values_by_q)
        n = sum(len(v) for v in values_by_q.values())
        p = _kw_pvalue([v for v in values_by_q.values() if len(v) > 1])
        ptxt = "p < 0.001" if p < 1e-3 else (f"p = {p:.3f}" if np.isfinite(p) else "p = NA")
        ax.set_title(f"{epoch}\nn = {n:,}\n{ptxt}")

    axes[0].set_ylabel("Cluster size  (log1p scale)")
    for ax in axes:
        ax.set_ylim(bottom=0)

    fig.supxlabel("SIMD quintile (cluster mode)", y=-0.03)
    return fig


def main(out_dir: Path | str = None) -> dict[str, Path]:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir: Path = Path(out_dir) if out_dir else paths.root / "analysis/ch1_socioeconomic/figures"

    cs = data.load_cluster_simd_features()

    fig = make_figure(cs)
    paths_out = style.save_figure(
        fig, out_dir / "fig2_cluster_size_by_simd",
        width="double",
        save_png=True, save_pdf=True,
    )
    plt.close(fig)

    summary = build_summary_table(cs)
    tables_dir = out_dir.parent / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(tables_dir / "fig2_cluster_size_by_simd.csv", index=False)
    paths_out["csv"] = tables_dir / "fig2_cluster_size_by_simd.csv"

    return paths_out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.figures)
    print(f"Wrote:\n   " + "\n   ".join(f"{k}: {v}" for k, v in p.items()))
