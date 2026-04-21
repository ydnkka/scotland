"""Figure 2 — Vaccination prevalence vs. cluster size, by VOC epoch.

Bin clusters by size quintile within each epoch. Plot mean `frac_vaccinated`
and bootstrap 95% CI per bin, with separate line per epoch.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, stats, style


def _binned_summary(df: pd.DataFrame, n_bins: int = 6) -> pd.DataFrame:
    """Within each epoch, split clusters into log-spaced size bins and
    compute mean + bootstrap CI of `frac_vaccinated`."""
    rows = []
    for epoch, sub in df.groupby("epoch", observed=True):
        sub = sub[sub["n_sequences"] > 0]
        if sub.empty:
            continue
        log_size = np.log1p(sub["n_sequences"])
        bins = np.linspace(log_size.min(), log_size.max(), n_bins + 1)
        cats = pd.cut(log_size, bins=bins, include_lowest=True)
        for cat, g in sub.groupby(cats, observed=True):
            v = g["frac_vaccinated"].dropna().to_numpy()
            if len(v) < 5:
                continue
            point, lo, hi = stats.bootstrap_ci(v, np.mean, n_boot=800)
            rows.append({
                "epoch": epoch,
                "size_mid": np.expm1((cat.left + cat.right) / 2),
                "mean": point, "lo": lo, "hi": hi, "n": len(v),
            })
    return pd.DataFrame(rows)


def build_data() -> pd.DataFrame:
    df = data.load_cluster_demographic_features()
    df = df[df["resolution"] == data.PRIMARY_RESOLUTION]
    df["epoch"] = data.assign_epoch(df["wn_mid_date"])
    df = df.dropna(subset=["frac_vaccinated", "epoch"])
    return df


def make_figure(summary: pd.DataFrame) -> plt.Figure:
    fig, ax = style.new_figure(width="onehalf", height_in=3.2)
    epoch_colors = {
        "Pre-VOC":      "#999999",
        "Alpha":        "#4e79a7",
        "Delta":        "#59a14f",
        "Omicron BA.1": "#e15759",
        "Omicron BA.2+":"#af2d2d",
    }
    for epoch, sub in summary.groupby("epoch", observed=True):
        sub = sub.sort_values("size_mid")
        c = epoch_colors.get(epoch, "#666666")
        ax.plot(sub["size_mid"], sub["mean"], "-o", color=c, label=epoch,
                ms=3.5, lw=1.2)
        ax.fill_between(sub["size_mid"], sub["lo"], sub["hi"], color=c, alpha=0.2, lw=0)
    ax.set_xscale("log")
    ax.set_xlabel("Cluster size (log scale)")
    ax.set_ylabel("Fraction vaccinated at sampling")
    ax.set_title("Vaccination prevalence vs. cluster size, by VOC epoch",
                 fontsize=9.5, fontweight="bold")
    ax.legend(loc="lower right", fontsize=6.8)
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper2_demographic/figures"
    df = build_data()
    summary = _binned_summary(df)
    summary.to_csv(out_dir.parent / "tables" / "fig2_vacc_vs_size.csv", index=False)
    fig = make_figure(summary)
    paths_out = style.save_figure(fig, out_dir / "fig2_vaccination_vs_cluster_size")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
