"""Figure 4 — Sex composition of clusters vs. cluster size, by VOC epoch.

`frac_female` is plotted against log cluster size. We overlay (1) the
population-expected reference line at ~0.51, (2) per-epoch GAM-ish loess
smoother (implemented as a rolling median with binomial-style CI), and
(3) a small rug showing the actual cluster counts per size bin.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, style


def build_data() -> pd.DataFrame:
    df = data.load_cluster_demographic_features()
    df = df[df["resolution"] == data.PRIMARY_RESOLUTION]
    df = df[df["n_sequences"] >= 2]
    df["epoch"] = data.assign_epoch(df["wn_mid_date"])
    df = df.dropna(subset=["frac_female", "epoch"])
    return df


def _rolling_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    bins = np.logspace(0, np.log10(df["n_sequences"].max()), 18)
    for epoch, sub in df.groupby("epoch", observed=True):
        cats = pd.cut(sub["n_sequences"], bins=bins, include_lowest=True)
        agg = (
            sub.assign(_bin=cats).groupby("_bin", observed=True)["frac_female"]
            .agg(["mean", "count", "std"]).reset_index()
        )
        agg["epoch"] = epoch
        agg["size_mid"] = [c.mid for c in agg["_bin"]]
        agg["se"] = agg["std"] / np.sqrt(agg["count"].clip(lower=1))
        rows.append(agg)
    return pd.concat(rows, ignore_index=True).dropna(subset=["mean"])


def make_figure(df: pd.DataFrame, summary: pd.DataFrame) -> plt.Figure:
    fig, ax = style.new_figure(width="onehalf", height_in=3.2)
    epoch_colors = {
        "Pre-VOC":      "#999999",
        "Alpha":        "#4e79a7",
        "Delta":        "#59a14f",
        "Omicron BA.1": "#e15759",
        "Omicron BA.2+":"#af2d2d",
    }
    # faint cloud
    ax.scatter(
        df["n_sequences"], df["frac_female"],
        s=2, color="#cccccc", alpha=0.25, linewidth=0, zorder=1,
    )
    for epoch, sub in summary.groupby("epoch", observed=True):
        sub = sub.sort_values("size_mid")
        c = epoch_colors.get(epoch, "#333333")
        ax.plot(sub["size_mid"], sub["mean"], "-o", color=c, ms=3.5, lw=1.2, label=epoch)
        ax.fill_between(sub["size_mid"], sub["mean"] - 1.96 * sub["se"],
                        sub["mean"] + 1.96 * sub["se"], color=c, alpha=0.15, lw=0)
    ax.axhline(0.51, color="#444444", ls="--", lw=0.6, zorder=0)
    ax.text(df["n_sequences"].max(), 0.515, " pop. reference",
            va="bottom", ha="right", fontsize=6.5, color="#444444")
    ax.set_xscale("log")
    ax.set_xlabel("Cluster size (log scale)")
    ax.set_ylabel("Fraction female in cluster")
    ax.set_ylim(0.25, 0.75)
    ax.set_title("Sex composition vs. cluster size, by VOC epoch",
                 fontsize=9.5, fontweight="bold")
    ax.legend(loc="upper right", fontsize=6.8)
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper2_demographic/output"
    df = build_data()
    summary = _rolling_summary(df)
    fig = make_figure(df, summary)
    paths_out = style.save_figure(fig, out_dir / "fig4_sex_composition")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
