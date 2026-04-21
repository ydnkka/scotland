"""Figure 6 — SIMD-domain decomposition of the cluster-size effect.

Fit a mutually-adjusted negative-binomial GLM that includes *all seven*
SIMD domain ranks simultaneously plus VOC dummies and a cr() time
spline on ``wn_mid_date``. Each domain's contribution is plotted as its
Shapley share of the total domain-attributable log-likelihood gain over
the covariate-only model.

Answers: when every SIMD domain is allowed to compete, which domains
carry the signal?
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
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
    fig, ax = style.new_figure(width="onehalf", height_in=3.0)

    colors = [style.SIMD_DOMAIN_PALETTE[d] for d in tab["domain"]]
    ax.barh(
        tab["domain"], tab["share"] * 100,
        color=colors, edgecolor="black", lw=0.4,
    )
    for y, (_, row) in enumerate(tab.iterrows()):
        sign = "+" if row["estimate"] > 0 else "−"
        ax.text(
            row["share"] * 100 + 0.6, y,
            f"{sign}{abs(row['estimate']):.2f}  (p={row['p_value']:.2g})",
            va="center", fontsize=6.8, color="#333333",
        )
    # ax.set_xlabel("Share of total domain-attributable log-likelihood gain  (%)")
    # ax.set_title(
    #     "SIMD-domain decomposition in a mutually-adjusted model",
    #     fontsize=9.0, fontweight="bold",
    # )
    ax.set_xlim(0, max(35, tab["share"].max() * 100 * 1.3))
    return fig


def main(out_dir: Path | None = None) -> dict[str, Path]:
    style.set_theme()
    paths = data.Paths.from_config()
    out_dir: Path = Path(out_dir) if out_dir else paths.root / "analysis/ch1_socioeconomic/figures"

    table_dir = out_dir.parent / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    table = table_dir / "fig6_domain_decomposition.csv"
    tab = _load_cached_table(table)
    if tab is None:
        cluster_df = simd_models.build_cluster_regression_frame()
        tab = simd_models.build_domain_decomposition_table_cv_shapley(cluster_df)
        simd_models.tag_model_output(tab).to_csv(table, index=False)
    fig = make_figure(tab)
    paths_out = style.save_figure(
        fig, out_dir / "fig6_domain_decomposition",
        width="onehalf", save_png=True, save_pdf=True,
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
