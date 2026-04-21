"""Figure 5 — Forest plot of demographic IRRs (cluster size) and ORs (singleton).

Two side-by-side forest plots from the mutually-adjusted models in
`models.cluster_demographics`. Predictors are standardised so IRRs/ORs are
per 1 SD of the predictor and directly comparable.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, stats, style
from manuscripts.paper2_demographic.models import cluster_demographics


DEMOG_TERMS = ["median_age", "frac_female", "frac_vaccinated", "mean_vacc_dose"]
LABELS = {
    "median_age":      "Median age  (+1 SD)",
    "frac_female":     "Frac. female  (+1 SD)",
    "frac_vaccinated": "Frac. vaccinated  (+1 SD)",
    "mean_vacc_dose":  "Mean vacc. dose  (+1 SD)",
}


def _prepare_tables(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    fit_size = cluster_demographics.cluster_size_model(frame)
    tab_size = stats.tidy_glm(fit_size).query("term in @DEMOG_TERMS").copy()

    fit_single = cluster_demographics.singleton_model(frame)
    tab_single = stats.tidy_glm(fit_single).query("term in @DEMOG_TERMS").copy()

    return tab_size, tab_single


def _forest(ax, tab: pd.DataFrame, title: str, x_label: str) -> None:
    tab = tab.set_index("term").loc[DEMOG_TERMS]
    y = np.arange(len(tab))
    ax.errorbar(
        tab["estimate"], y,
        xerr=[tab["estimate"] - tab["conf_low"], tab["conf_high"] - tab["estimate"]],
        fmt="o", ms=5, color="#333333", ecolor="#888", elinewidth=0.8, capsize=2,
        markerfacecolor="#4e79a7", markeredgecolor="black", markeredgewidth=0.4,
    )
    ax.axvline(1.0, color="#bbbbbb", lw=0.7, ls="--", zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS[t] for t in tab.index])
    ax.set_xscale("log")
    ax.set_xlabel(x_label)
    ax.set_title(title, fontsize=9.0)


def make_figure(tab_size: pd.DataFrame, tab_single: pd.DataFrame) -> plt.Figure:
    fig, (ax1, ax2) = style.new_figure(
        width="double", height_in=2.8, nrows=1, ncols=2, sharey=True,
        gridspec_kw={"wspace": 0.08},
    )
    _forest(ax1, tab_size,   "Cluster size (IRR per 1 SD)", "IRR  (log scale)")
    _forest(ax2, tab_single, "Singleton status (OR per 1 SD)", "OR  (log scale)")
    fig.suptitle("Demographic predictors of cluster size and singleton status",
                 x=0.02, ha="left", y=1.03, fontsize=9.5, fontweight="bold")
    return fig


def main(out_dir: Path | None = None) -> Path:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir = Path(out_dir) if out_dir else paths.root / "manuscripts/paper2_demographic/figures"
    frame = cluster_demographics.build_cluster_regression_frame()
    tab_size, tab_single = _prepare_tables(frame)
    tab_size.to_csv(out_dir.parent / "tables" / "fig5_irr_size.csv", index=False)
    tab_single.to_csv(out_dir.parent / "tables" / "fig5_or_singleton.csv", index=False)
    fig = make_figure(tab_size, tab_single)
    paths_out = style.save_figure(fig, out_dir / "fig5_demographic_forest")
    plt.close(fig)
    return paths_out[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.output)
    print(f"Wrote {p}")
