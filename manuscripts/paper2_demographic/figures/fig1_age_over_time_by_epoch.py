"""Figure 1 — Median cluster age distribution by VOC epoch.

A ridgeline (stacked KDE) showing how the median age of sequenced cluster
members shifted across the Alpha, Delta, and Omicron sub-epochs. One ridge
per epoch, shared x-axis (years of age).

Excludes singletons (n_sequences == 1) to focus on clusters where a median
over multiple members is meaningful.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, style


def _kde(x: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Gaussian KDE with Silverman's rule, NaN-safe."""
    from scipy.stats import gaussian_kde

    x = x[np.isfinite(x)]
    if len(x) < 5:
        return np.zeros_like(grid)
    try:
        kde = gaussian_kde(x, bw_method="silverman")
        return kde(grid)
    except Exception:
        return np.zeros_like(grid)


def build_panel_data() -> pd.DataFrame:
    df = data.load_cluster_demographic_features()
    df = df[df["resolution"] == data.PRIMARY_RESOLUTION]
    df = df[df["is_singleton"] == 0]
    df = df.dropna(subset=["median_age", "wn_mid_date"])
    df["epoch"] = data.assign_epoch(df["wn_mid_date"])
    df = df.dropna(subset=["epoch"])
    return df


def make_figure(df: pd.DataFrame) -> plt.Figure:
    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS]
    fig, ax = style.new_figure(width="onehalf", height_in=3.8)

    grid = np.linspace(0, 100, 400)
    y_step = 1.0
    for i, epoch in enumerate(epochs):
        sub = df.loc[df["epoch"] == epoch, "median_age"].to_numpy()
        dens = _kde(sub, grid)
        dens = dens / (dens.max() + 1e-9) * 0.9  # normalise each ridge
        y = i * y_step
        ax.fill_between(grid, y, y + dens, color="#4e79a7", alpha=0.35, lw=0)
        ax.plot(grid, y + dens, color="#2b4f72", lw=1.0)
        median = float(np.nanmedian(sub))
        ax.plot([median, median], [y, y + 0.2], color="black", lw=0.8)
        ax.text(
            101, y + 0.05,
            f"  {epoch}\n  n = {len(sub):,}  median={median:.0f}y",
            va="bottom", ha="left", fontsize=7,
        )
    ax.set_xlim(0, 100)
    ax.set_xlabel("Median age of cluster members (years)")
    ax.set_ylim(-0.1, len(epochs) * y_step)
    ax.set_yticks([])
    ax.set_title("Median cluster age shifts across VOC epochs",
                 fontsize=9.5, fontweight="bold")
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper2_demographic/figures"
    df = build_panel_data()
    fig = make_figure(df)
    paths_out = style.save_figure(fig, out_dir / "fig1_age_over_time_by_epoch")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
