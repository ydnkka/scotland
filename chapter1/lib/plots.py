"""Diagnostic plots for Chapter 1.

Three plot families, written to PNG + PDF:

* :func:`plot_main_effects_forest` — forest plot of coefficients (rate ratios
  for count components; log-odds for the spread hurdle), grouped by outcome
  × component.
* :func:`plot_wave_interaction_slopes` — for each excess-mixing predictor,
  the implied slope (main effect + interaction) for each wave, with 95% CIs.
* :func:`plot_domain_or_wave_grid` — a generic grid of forest panels keyed
  by ``domain`` or ``wave`` for the stratified analyses.

All plots are diagnostic, not publication-grade: simple matplotlib with the
project's plot style if available.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .constants import EXCESS_MIXING_TERMS, TERM_LABELS


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------


def _set_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.grid.axis": "x",
        "grid.alpha": 0.25,
        "legend.frameon": False,
    })


def _save(fig: plt.Figure, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def _component_x(component: str, results: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series, float]:
    """Pick the right effect-scale series for plotting."""
    if component == "positive_count":
        # rate ratio on the log mean
        return (
            results["ratio"],
            results["ratio_lower"],
            results["ratio_upper"],
            1.0,
        )
    if component == "hurdle":
        return (
            results["ratio"],
            results["ratio_lower"],
            results["ratio_upper"],
            1.0,
        )
    return (
        results["estimate"],
        results["estimate_lower"],
        results["estimate_upper"],
        0.0,
    )


# ---------------------------------------------------------------------------
# Main-effects forest plot
# ---------------------------------------------------------------------------


def plot_main_effects_forest(
    results: pd.DataFrame,
    out_path: Path,
    *,
    terms: Iterable[str] | None = None,
    title: str = "Excess mixing → cluster scale",
) -> None:
    _set_style()
    terms = list(terms) if terms is not None else list(EXCESS_MIXING_TERMS)
    sub = results[results["term"].isin(terms)].copy()
    if sub.empty:
        return

    groups = list(sub.groupby(["outcome", "component"]).groups.keys())
    n_panels = len(groups)
    n_cols = min(2, n_panels)
    n_rows = math.ceil(n_panels / n_cols)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(5.4 * n_cols, max(2.6, 1.2 * len(terms)) * n_rows),
        squeeze=False,
    )

    for ax, (outcome, component) in zip(axes.flat, groups):
        panel = (
            sub[(sub["outcome"] == outcome) & (sub["component"] == component)]
            .set_index("term")
            .reindex(terms)
            .reset_index()
        )
        est, lo, hi, ref = _component_x(component, panel)
        y_pos = np.arange(len(panel))[::-1]
        labels = [TERM_LABELS.get(t, t) for t in panel["term"]]
        ax.errorbar(
            est, y_pos,
            xerr=[est - lo, hi - est],
            fmt="o", capsize=3, color="#1f3a5f",
        )
        ax.axvline(ref, color="#888", linestyle=":", linewidth=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_xlabel(
            "Rate ratio" if component != "linear" else "Estimate"
        )
        ax.set_title(f"{outcome} — {component}", fontsize=10)

    for ax in axes.flat[n_panels:]:
        ax.set_visible(False)

    fig.suptitle(title, fontsize=12, y=1.02)
    fig.tight_layout()
    _save(fig, out_path)


# ---------------------------------------------------------------------------
# Wave-interaction slope plot
# ---------------------------------------------------------------------------


def plot_wave_interaction_slopes(
    results: pd.DataFrame,
    out_path: Path,
    *,
    title: str = "Excess-mixing slopes by wave",
) -> None:
    """For each excess-mixing predictor and each (outcome, component),
    plot the wave-specific slope as (main effect + interaction)."""
    _set_style()
    main_terms = list(EXCESS_MIXING_TERMS)

    sub = results.copy()
    # split rows into main effect (term == predictor) and interaction
    # rows (term starts with predictor + "__x__wave_").
    panels = list(sub.groupby(["outcome", "component"]).groups.keys())
    n_panels = len(panels)
    if n_panels == 0:
        return
    n_cols = min(2, n_panels)
    n_rows = math.ceil(n_panels / n_cols)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(6.4 * n_cols, max(3.0, 1.4 * len(main_terms)) * n_rows),
        squeeze=False,
    )

    colors = {
        "age_excess_mixing_z":  "#3a5fbf",
        "sex_excess_mixing_z":  "#bf553a",
        "simd_excess_mixing_z": "#3aa05f",
    }

    for ax, (outcome, component) in zip(axes.flat, panels):
        panel = sub[(sub["outcome"] == outcome) & (sub["component"] == component)]
        # Collect waves that appear in interaction columns
        waves = sorted({
            t.split("__x__wave_")[1]
            for t in panel["term"]
            if "__x__wave_" in t
        })
        if not waves:
            ax.set_visible(False)
            continue
        x_idx = np.arange(len(waves))
        for j, predictor in enumerate(main_terms):
            main = panel[panel["term"] == predictor]
            if main.empty:
                continue
            main_est = float(main["estimate"].iloc[0])
            main_se = float(main["std_error"].iloc[0])
            ests, los, his = [], [], []
            for wave in waves:
                inter = panel[panel["term"] == f"{predictor}__x__wave_{wave}"]
                if inter.empty:
                    ests.append(main_est)
                    los.append(main_est - 1.96 * main_se)
                    his.append(main_est + 1.96 * main_se)
                else:
                    inter_est = float(inter["estimate"].iloc[0])
                    inter_se = float(inter["std_error"].iloc[0])
                    est = main_est + inter_est
                    se = math.sqrt(main_se ** 2 + inter_se ** 2)
                    ests.append(est)
                    los.append(est - 1.96 * se)
                    his.append(est + 1.96 * se)
            offset = (j - 1) * 0.2
            ax.errorbar(
                x_idx + offset, ests,
                yerr=[np.array(ests) - np.array(los),
                      np.array(his) - np.array(ests)],
                fmt="o", capsize=2.5,
                color=colors.get(predictor, "#444"),
                label=TERM_LABELS.get(predictor, predictor),
            )
        ax.axhline(0, color="#888", linestyle=":", linewidth=0.8)
        ax.set_xticks(x_idx)
        ax.set_xticklabels(waves, rotation=45, ha="right")
        ax.set_ylabel("Slope (log-scale)")
        ax.set_title(f"{outcome} — {component}", fontsize=10)
        ax.legend(loc="best", fontsize=8)

    for ax in axes.flat[n_panels:]:
        ax.set_visible(False)

    fig.suptitle(title, fontsize=12, y=1.02)
    fig.tight_layout()
    _save(fig, out_path)


# ---------------------------------------------------------------------------
# Domain / wave-stratified grid
# ---------------------------------------------------------------------------


def plot_stratified_forest(
    results: pd.DataFrame,
    out_path: Path,
    *,
    stratum_col: str,
    terms: Iterable[str] | None = None,
    title: str = "",
) -> None:
    _set_style()
    terms = list(terms) if terms is not None else list(EXCESS_MIXING_TERMS)
    sub = results[results["term"].isin(terms)].copy()
    if sub.empty or stratum_col not in sub.columns:
        return

    panels = list(sub.groupby(["outcome", "component"]).groups.keys())
    n_panels = len(panels)
    n_cols = min(2, n_panels)
    n_rows = math.ceil(n_panels / n_cols)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(6.4 * n_cols, 0.45 * sub[stratum_col].nunique() * len(terms) * n_rows),
        squeeze=False,
    )

    strata = list(sub[stratum_col].dropna().unique())
    colors = {
        "age_excess_mixing_z":  "#3a5fbf",
        "sex_excess_mixing_z":  "#bf553a",
        "simd_excess_mixing_z": "#3aa05f",
        "income_domain_excess_mixing_z":     "#3aa05f",
        "employment_domain_excess_mixing_z": "#3aa05f",
        "education_domain_excess_mixing_z":  "#3aa05f",
        "health_domain_excess_mixing_z":     "#3aa05f",
        "access_domain_excess_mixing_z":     "#3aa05f",
        "crime_domain_excess_mixing_z":      "#3aa05f",
        "housing_domain_excess_mixing_z":    "#3aa05f",
        "overall_domain_excess_mixing_z":    "#3aa05f",
    }

    for ax, (outcome, component) in zip(axes.flat, panels):
        panel = sub[
            (sub["outcome"] == outcome) & (sub["component"] == component)
        ].copy()
        labels = []
        positions = []
        pos = 0
        for stratum in strata:
            block = panel[panel[stratum_col] == stratum]
            for term in terms:
                row = block[block["term"] == term]
                if row.empty:
                    continue
                est, lo, hi, ref = _component_x(component, row)
                ax.errorbar(
                    float(est.iloc[0]), pos,
                    xerr=[[float(est.iloc[0]) - float(lo.iloc[0])],
                          [float(hi.iloc[0]) - float(est.iloc[0])]],
                    fmt="o", capsize=2.5,
                    color=colors.get(term, "#444"),
                )
                labels.append(f"{stratum} · {TERM_LABELS.get(term, term)}")
                positions.append(pos)
                pos += 1
            pos += 0.5
        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.axvline(
            1.0 if component != "linear" else 0.0,
            color="#888", linestyle=":", linewidth=0.8,
        )
        ax.set_xlabel("Rate ratio" if component != "linear" else "Estimate")
        ax.set_title(f"{outcome} — {component}", fontsize=10)

    for ax in axes.flat[n_panels:]:
        ax.set_visible(False)

    if title:
        fig.suptitle(title, fontsize=12, y=1.02)
    fig.tight_layout()
    _save(fig, out_path)
