"""Figure 3 — Distance-decay of cluster co-membership, per VOC epoch.

For pairs of sequences sampled within the same window, plot the probability
they belong to the same cluster as a function of pairwise DZ-centroid
distance. Pairs are binned in log-spaced distance bins; Wilson 95% CIs
are shown. One line per VOC epoch.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, style
from manuscripts.paper3_spatial.models.spatial_mixing import sample_pairs_within_windows


def _wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def build_data(n_per_window: int = 1500) -> pd.DataFrame:
    return sample_pairs_within_windows(n_per_window=n_per_window)


def _summarise(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["epoch"] = data.assign_epoch(pd.to_datetime(df["wn_mid_date"]))
    df = df.dropna(subset=["epoch"])
    bins = np.logspace(-1, 3, 20)  # 0.1 km to 1000 km
    rows = []
    for epoch, sub in df.groupby("epoch", observed=True):
        sub = sub[sub["distance_km"] > 0]
        cats = pd.cut(sub["distance_km"], bins=bins, include_lowest=True)
        for cat, g in sub.groupby(cats, observed=True):
            n = len(g)
            k = int(g["co_cluster"].sum())
            lo, hi = _wilson_ci(k, n)
            rows.append({
                "epoch": epoch,
                "distance_km_mid": (cat.left + cat.right) / 2,
                "p_co_cluster": k / n if n > 0 else np.nan,
                "ci_low": lo, "ci_high": hi, "n_pairs": n,
            })
    return pd.DataFrame(rows).dropna(subset=["p_co_cluster"])


def make_figure(summary: pd.DataFrame) -> plt.Figure:
    fig, ax = style.new_figure(width="onehalf", height_in=3.3)
    epoch_colors = {
        "Pre-VOC":      "#999999",
        "Alpha":        "#4e79a7",
        "Delta":        "#59a14f",
        "Omicron BA.1": "#e15759",
        "Omicron BA.2+":"#af2d2d",
    }
    for epoch, sub in summary.groupby("epoch", observed=True):
        sub = sub.sort_values("distance_km_mid")
        c = epoch_colors.get(epoch, "#333333")
        ax.plot(sub["distance_km_mid"], sub["p_co_cluster"], "-o",
                color=c, ms=3.5, lw=1.1, label=epoch)
        ax.fill_between(sub["distance_km_mid"], sub["ci_low"], sub["ci_high"],
                        color=c, alpha=0.15, lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Pairwise DZ-centroid distance  (km, log)")
    ax.set_ylabel("Pr(same cluster)  (log)")
    ax.set_title("Distance-decay of cluster co-membership by VOC epoch",
                 fontsize=9.5, fontweight="bold")
    ax.legend(loc="lower left", fontsize=7)
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper3_spatial/output"
    pair_df = build_data()
    summary = _summarise(pair_df)
    summary.to_csv(out_dir.parent / "tables" / "fig3_distance_decay.csv", index=False)
    fig = make_figure(summary)
    paths_out = style.save_figure(fig, out_dir / "fig3_distance_decay")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
