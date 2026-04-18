"""Figure 3 — Within-cluster age diversity vs. cluster size.

Hexbin density of clusters on (log cluster size, age std) with a running
median (and 10th / 90th percentile band) overlaid. A null-expectation curve
is added: std of age midpoints under random sampling without replacement
from the age distribution of the relevant VOC epoch, simulated once.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, style

N_NULL_DRAWS_PER_SIZE = 500


def build_data() -> pd.DataFrame:
    df = data.load_cluster_demographic_features()
    df = df[df["resolution"] == data.PRIMARY_RESOLUTION]
    df = df[df["n_sequences"] >= 2]
    df = df.dropna(subset=["age_diversity", "n_sequences"])
    return df


def _null_curve(age_pool: np.ndarray, sizes: np.ndarray, rng) -> np.ndarray:
    """Expected std of age midpoints under random draws without replacement."""
    out = []
    age_pool = age_pool[np.isfinite(age_pool)]
    for s in sizes:
        if len(age_pool) < s:
            out.append(np.nan); continue
        draws = np.stack(
            [rng.choice(age_pool, size=int(s), replace=False) for _ in range(N_NULL_DRAWS_PER_SIZE)]
        )
        out.append(np.mean(draws.std(axis=1)))
    return np.asarray(out)


def make_figure(df: pd.DataFrame) -> plt.Figure:
    fig, ax = style.new_figure(width="onehalf", height_in=3.2)
    x = np.log1p(df["n_sequences"].to_numpy())
    y = df["age_diversity"].to_numpy()

    hb = ax.hexbin(x, y, gridsize=35, cmap="viridis", mincnt=1, linewidths=0.1)
    cb = fig.colorbar(hb, ax=ax, pad=0.02, shrink=0.8)
    cb.set_label("Clusters per hex", fontsize=7)
    cb.ax.tick_params(labelsize=7)

    # Rolling median and band
    df_sorted = df.sort_values("n_sequences")
    bins = np.logspace(0, np.log10(df["n_sequences"].max()), 25)
    cats = pd.cut(df_sorted["n_sequences"], bins=bins, include_lowest=True)
    rolled = (
        df_sorted.groupby(cats, observed=True)["age_diversity"]
        .agg(["median", lambda v: v.quantile(0.1), lambda v: v.quantile(0.9), "count"])
    )
    rolled.columns = ["median", "p10", "p90", "count"]
    rolled = rolled[rolled["count"] >= 10]
    centres = np.array([c.mid for c in rolled.index])
    ax.plot(np.log1p(centres), rolled["median"], color="#ffffff", lw=2.2, zorder=3)
    ax.plot(np.log1p(centres), rolled["median"], color="#d62728", lw=1.4, zorder=4,
            label="Observed median")
    ax.fill_between(np.log1p(centres), rolled["p10"], rolled["p90"],
                    color="#d62728", alpha=0.15, zorder=2)

    # Null expectation (single sample of the age pool over all clusters)
    age_pool = data.load_analysis_columns(
        ["age_midpoint", "sequence_id"], resolution=data.PRIMARY_RESOLUTION
    ).drop_duplicates("sequence_id")["age_midpoint"].to_numpy()
    rng = np.random.default_rng(42)
    null = _null_curve(age_pool, centres.astype(int), rng)
    ax.plot(np.log1p(centres), null, color="#333333", lw=1.0, ls="--",
            label="Random-draw null")

    ax.set_xlabel("log(1 + cluster size)")
    ax.set_ylabel("Within-cluster age std (years)")
    ax.set_title("Age homogeneity shrinks below random expectation in larger clusters",
                 fontsize=9.5, fontweight="bold", loc="left")
    ax.legend(loc="lower right", fontsize=7)
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper2_demographic/output"
    df = build_data()
    fig = make_figure(df)
    paths_out = style.save_figure(fig, out_dir / "fig3_age_homogeneity")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
