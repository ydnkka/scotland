"""Figure 6 — VOC-stratified shifts in the demographic-size association.

Refit the cluster-size NB GLM **within each VOC epoch** and plot, for each
standardised predictor, its IRR and 95% CI across epochs. This visualises
how the strength of each demographic effect changes as dominant lineage and
population immunity evolve.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, stats, style
from manuscripts.paper2_demographic.models import cluster_demographics

PREDICTORS = ["median_age", "frac_female", "frac_vaccinated", "mean_vacc_dose"]
LABELS = {
    "median_age":      "Median age",
    "frac_female":     "Frac. female",
    "frac_vaccinated": "Frac. vaccinated",
    "mean_vacc_dose":  "Mean vacc. dose",
}


def _fit_epoch(epoch_df: pd.DataFrame) -> pd.DataFrame:
    if len(epoch_df) < 100:
        return pd.DataFrame()
    fit = cluster_demographics.cluster_size_model(epoch_df)
    tidy = stats.tidy_glm(fit)
    tidy = tidy[tidy["term"].isin(PREDICTORS)].copy()
    return tidy


def build_table() -> pd.DataFrame:
    frame = cluster_demographics.build_cluster_regression_frame()
    frame["epoch"] = data.assign_epoch(frame["wn_mid_date"])
    rows = []
    for epoch in [lbl for lbl, *_ in data.VOC_EPOCHS]:
        sub = frame[frame["epoch"] == epoch]
        tab = _fit_epoch(sub)
        if tab.empty:
            continue
        tab["epoch"] = epoch
        tab["n_obs"] = int(len(sub.dropna(subset=["median_age", "frac_female",
                                                  "frac_vaccinated", "mean_vacc_dose",
                                                  "wn_prop_sequenced", "n_sequences"])))
        rows.append(tab)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def make_figure(tab: pd.DataFrame) -> plt.Figure:
    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS if lbl in set(tab["epoch"])]
    fig, axes = style.new_figure(
        width="double", height_in=3.2, nrows=1, ncols=len(PREDICTORS), sharey=True,
        gridspec_kw={"wspace": 0.08},
    )
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])

    y = np.arange(len(epochs))
    for ax, pred in zip(axes, PREDICTORS):
        sub = tab[tab["term"] == pred].set_index("epoch").reindex(epochs)
        ax.errorbar(
            sub["estimate"], y,
            xerr=[sub["estimate"] - sub["conf_low"], sub["conf_high"] - sub["estimate"]],
            fmt="o", ms=5, color="#333333", ecolor="#888", elinewidth=0.8, capsize=2,
            markerfacecolor="#59a14f", markeredgecolor="black", markeredgewidth=0.4,
        )
        ax.axvline(1.0, color="#bbbbbb", lw=0.7, ls="--", zorder=0)
        ax.set_xscale("log")
        ax.set_title(LABELS[pred], fontsize=8.8)
        ax.set_xlabel("IRR (log)")

    axes[0].set_yticks(y)
    axes[0].set_yticklabels(epochs)

    fig.suptitle("Within-epoch IRRs for demographic predictors of cluster size",
                 x=0.02, ha="left", y=1.02, fontsize=9.5, fontweight="bold")
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper2_demographic/output"
    tab = build_table()
    tab.to_csv(out_dir.parent / "tables" / "fig6_voc_stratified_irrs.csv", index=False)
    fig = make_figure(tab)
    paths_out = style.save_figure(fig, out_dir / "fig6_voc_stratified_shifts")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
