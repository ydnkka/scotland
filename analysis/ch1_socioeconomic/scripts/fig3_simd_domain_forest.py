"""Figure 3 — Forest plot of IRRs for cluster size by SIMD domain.

For each SIMD domain (overall, income, employment, education, health, access,
crime, housing), fit a separate negative-binomial GLM predicting cluster size
from the standardised domain rank, adjusting for VOC, a cr() time spline on
``wn_mid_date``, and ``log(wn_prop_sequenced)`` offset. Plot the IRR for a
1-SD decrease in rank (i.e. *more* deprivation) with its 95% CI.

Interpretation: IRR > 1 means more deprivation is associated with larger clusters.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis.untils import data, style
from analysis.ch1_socioeconomic.models import simd_models


def make_figure(tab: pd.DataFrame) -> plt.Figure:
    tab = tab.sort_values("estimate")
    fig, ax = style.new_figure(width="onehalf", height_in=3.2)

    y = np.arange(len(tab))
    ax.errorbar(
        tab["estimate"], y,
        xerr=[tab["estimate"] - tab["conf_low"], tab["conf_high"] - tab["estimate"]],
        fmt="none", ecolor="#777777", elinewidth=0.8, capsize=2, zorder=1,
    )
    ax.scatter(
        tab["estimate"], y,
        c=[style.SIMD_DOMAIN_PALETTE[d] for d in tab.index],
        edgecolors="black", linewidths=0.4, s=25, zorder=2,
    )
    ax.axvline(1.0, color="#bbbbbb", lw=0.7, ls="--", zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([d.capitalize() for d in tab.index])
    ax.set_xlabel("IRR for a 1-SD *increase* in deprivation  (95% CI)")
    # ax.set_title("Cluster-size effect of SIMD domain, adjusted for VOC and surveillance")
    # Annotate n and p on the right margin.
    for yi, (name, row) in zip(y, tab.iterrows()):
        ax.text(
            ax.get_xlim()[1], yi,
            f"p={row['p_value']:.2g}",
            va="center", fontsize=6.8, color="#555555",
        )
    return fig


def main(out_dir: Path | str = None) -> dict[str, Path]:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir: Path = Path(out_dir) if out_dir else paths.root / "analysis/ch1_socioeconomic/figures"
    table = out_dir.parent / "tables" / "fig3_domain_irrs.csv"

    if table.exists():
        tab = pd.read_csv(table, index_col=0)
        fig = make_figure(tab)
        paths_out = style.save_figure(
            fig, out_dir / "fig3_simd_domain_forest",
            width="onehalf", save_png=True, save_pdf=True,
        )
    else:
        cluster_df = simd_models.build_cluster_regression_frame()
        tab = simd_models.build_domain_forest_table(cluster_df)
        tab.to_csv(out_dir.parent / "tables" / "fig3_domain_irrs.csv")
        fig = make_figure(tab)
        paths_out = style.save_figure(
            fig, out_dir / "fig3_simd_domain_forest",
            width="onehalf", save_png=True, save_pdf=True,
        )
    plt.close(fig)
    return paths_out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", type=Path, default=None)
    args = ap.parse_args()
    p = main(args.figures)
    print(f"Wrote:\n   " + "\n   ".join(f"{k}: {v}" for k, v in p.items()))
