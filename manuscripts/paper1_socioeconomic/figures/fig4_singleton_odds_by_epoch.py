"""Figure 4 — Odds of singleton vs. non-singleton cluster by SIMD quintile, by VOC epoch.

One panel per VOC epoch. Each panel shows the adjusted OR (logistic GLM)
for cluster singleton-ness given SIMD quintile, using Q5 (least deprived)
as the reference. Errors bars are 95% Wald CIs.

Interpretation: OR < 1 in Q1 means members of the most deprived quintile
are *less* likely to be in a singleton cluster - i.e. *more* likely to be
part of a genetically linked onward chain.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, style
from manuscripts.paper1_socioeconomic.models import simd_models


def make_figure(tab: pd.DataFrame) -> plt.Figure:
    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS if lbl in set(tab["epoch"])]
    fig, axes = style.new_figure(
        width="double", height_in=3.4, nrows=1, ncols=len(epochs),
        sharey=True, gridspec_kw={"wspace": 0.1},
    )
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])

    quintiles_plotted = [1, 2, 3, 4]  # 5 is reference
    for ax, epoch in zip(axes, epochs):
        sub = tab[tab["epoch"] == epoch].sort_values("quintile")
        if sub.empty:
            ax.set_axis_off()
            continue
        est = sub.set_index("quintile").reindex(quintiles_plotted)
        for q, row in est.iterrows():
            values = row[["estimate", "conf_low", "conf_high"]].to_numpy(dtype=float)
            if not np.isfinite(values).all():
                continue
            estimate, conf_low, conf_high = values
            ax.errorbar(
                estimate,
                float(q),
                xerr=np.array([[estimate - conf_low], [conf_high - estimate]]),
                fmt="o",
                markersize=5,
                color="#333333",
                ecolor="#888888",
                elinewidth=0.8,
                capsize=2,
                markerfacecolor=style.SIMD_QUINTILE_PALETTE[int(q)],
                markeredgecolor="black",
                markeredgewidth=0.4,
            )
        ax.axvline(1.0, color="#bbbbbb", lw=0.7, ls="--", zorder=0)
        ax.set_yticks(quintiles_plotted)
        ax.set_yticklabels([f"Q{q}" for q in quintiles_plotted])
        ax.set_title(epoch, fontsize=8.8)

    axes[0].set_ylabel("SIMD quintile")
    fig.supxlabel("OR vs. Q5 (least deprived)")
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper1_socioeconomic/output"
    cluster_df = simd_models.build_cluster_regression_frame()
    tab = simd_models.build_singleton_epoch_table(cluster_df)
    table_dir = out_dir.parent / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    tab.to_csv(table_dir / "fig4_singleton_ors.csv", index=False)
    fig = make_figure(tab)
    paths_out = style.save_figure(fig, out_dir / "fig4_singleton_odds_by_epoch")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
