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

from analysis.utils import data, style
from analysis.ch1_socioeconomic.models import simd_models


def _load_cached_table(table: Path) -> pd.DataFrame | None:
    if not table.exists():
        return None
    tab = pd.read_csv(table)
    if not simd_models.is_current_model_output(tab):
        return None
    return tab.drop(columns=["model_version"])


def make_figure(tab: pd.DataFrame) -> plt.Figure:
    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS if lbl in set(tab["epoch"])]
    n_epochs = len(epochs)
    quintiles_plotted = [1, 2, 3, 4]   # 5 is reference

    fig, axes = style.new_figure(
        width="double",
        height_in=3.5,
        nrows=1,
        ncols=n_epochs,
        sharey=True,
        gridspec_kw={"wspace": 0.12},
    )
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])

    # y-offsets so primary and QC-adjusted markers don't overlap
    Y_OFFSET = {"primary": -0.15, "qc_adjusted": +0.15}
    MARKER = {"primary": "o",   "qc_adjusted": "D"}
    MFILL = {"primary": True,  "qc_adjusted": False}
    LABEL = {"primary": "Good-QC only (primary)",
             "qc_adjusted": "All QC + covariates (sensitivity)"}

    legend_handles: dict[str, object] = {}

    for ax, epoch in zip(axes, epochs):
        sub = tab[tab["epoch"] == epoch]
        if sub.empty:
            ax.set_axis_off()
            continue

        for model in ("primary", "qc_adjusted"):
            msub = sub[sub["model"] == model].set_index("quintile").reindex(quintiles_plotted)
            for q, row in msub.iterrows():
                vals = row[["estimate", "conf_low", "conf_high"]].to_numpy(dtype=float)
                if not np.isfinite(vals).all():
                    continue
                est, clo, chi = vals
                color = style.SIMD_QUINTILE_PALETTE[int(q)]
                mfc = color if MFILL[model] else "white"
                y_pos = float(q) + Y_OFFSET[model]
                h = ax.errorbar(
                    est, y_pos,
                    xerr=np.array([[est - clo], [chi - est]]),
                    fmt=MARKER[model],
                    markersize=4.5,
                    color=color,
                    ecolor=color,
                    elinewidth=0.7,
                    capsize=2.0,
                    markerfacecolor=mfc,
                    markeredgecolor=color,
                    markeredgewidth=0.8,
                    label=LABEL[model] if q == 1 else "_nolegend_",
                )
                if model not in legend_handles:
                    legend_handles[model] = h

        ax.axvline(1.0, color="#bbbbbb", lw=0.7, ls="--", zorder=0)
        ax.set_yticks(quintiles_plotted)
        ax.set_yticklabels([f"Q{q}" for q in quintiles_plotted])
        ax.set_title(epoch)

    axes[0].set_ylabel("SIMD quintile (vs Q5)")
    fig.supxlabel("OR vs. Q5 (least deprived)", y=-0.02)

    # Legend below the figure using the first two handles
    handles = [legend_handles.get("primary"), legend_handles.get("qc_adjusted")]
    labels = [LABEL["primary"], LABEL["qc_adjusted"]]
    handles_clean = [(h, l) for h, l in zip(handles, labels) if h is not None]
    if handles_clean:
        hh, ll = zip(*handles_clean)
        fig.legend(
            hh, ll,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.15),
            ncol=2,
            frameon=False,
        )
    return fig


def main(out_dir: Path | str = None) -> dict[str, Path]:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir: Path = Path(out_dir) if out_dir else paths.root / "analysis/ch1_socioeconomic/figures"
    table_dir = out_dir.parent / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    table = table_dir / "fig4_singleton_ors.csv"
    tab = _load_cached_table(table)
    if tab is None:
        cluster_df = simd_models.build_cluster_regression_frame()
        tab = simd_models.build_singleton_epoch_table(cluster_df)
        simd_models.tag_model_output(tab).to_csv(table, index=False)
    fig = make_figure(tab)
    paths_out = style.save_figure(
        fig, out_dir / "fig4_singleton_odds_by_epoch",
        width="double", save_png=True, save_pdf=True
    )
    plt.close(fig)
    paths_out["csv"] = table
    return paths_out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.figures)
    print(f"Wrote:\n   " + "\n   ".join(f"{k}: {v}" for k, v in p.items()))
