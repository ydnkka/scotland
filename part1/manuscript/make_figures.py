"""Create publication-ready figures for the Part 1 main analysis.

The manuscript is organised around two complementary lines of inquiry:

* Line 1 — deprivation as exposure for cluster outcomes (size, spread, mixing).
* Line 2 — within-cluster excess mixing as predictor of cluster scale.

For each line there are two main figures (overall and wave-specific). Domain
extensions, outcome-distribution descriptives, mixing distributions, and other
sensitivities live in the supplement.

Main figures:
    fig1_deprivation_overall         — Line 1 overall (counts + mixing)
    fig2_deprivation_wave_specific   — Line 1 by epidemic wave
    fig3_mixing_overall              — Line 2 overall (mixing predictors)
    fig4_mixing_wave_specific        — Line 2 by epidemic wave

Supplementary figures:
    supp_fig1_outcome_distributions
    supp_fig2_mixing_distributions
    supp_fig3_observed_expected_matrices
    supp_fig4_deprivation_domain_outcomes
    supp_fig5_deprivation_domain_mixing
    supp_fig6_deprivation_domain_wave_mixing
    supp_fig7_mixing_domain_outcomes
    supp_fig8_deprivation_size_adjusted
    supp_fig9_deprivation_loglinear
    supp_fig10_mixing_loglinear

Outputs are written to ``part1/manuscript/figures`` as PDF, PNG, and TIFF.
The script uses the shared project plotting module at ``utils/style.py``.

Note on mixing-predictor models: excess mixing is undefined for singletons, so
the cluster-size hurdle component is not estimable (its comparison group would
be the excluded singletons). The cluster-size figure for Line 2 therefore shows
only the positive ZTNB component, while geographic spread shows both the
hurdle (multi- vs single-datazone among non-singletons) and the positive ZTNB
component, as in the underlying results files.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Constants and labels
# ---------------------------------------------------------------------------

OUTCOME_LABELS = {
    "cluster_size": "Cluster size",
    "geographic_dispersion": "Geographic spread",
    "geographic_dispersion_size_adjusted": "Geographic spread, size-adjusted",
}

TERM_LABELS = {
    "deprivation_z": "SIMD deprivation",
    "index_deprivation_z": "Index-case SIMD deprivation",
    "local_incidence_z": "Local incidence",
    "local_seq_fraction_z": "Local sequencing",
    "window_seq_fraction_z": "Window sequencing",
    "test_positivity_z": "Test positivity",
    "log_cluster_size_z": "Cluster size",
    "simd_excess_mixing_z": "SIMD excess mixing",
    "age_excess_mixing_z": "Age excess mixing",
    "sex_excess_mixing_z": "Sex excess mixing",
    "profile_excess_mixing_z": "Joint-profile excess mixing",
    "age_sex_excess_mixing_z": "Age-sex excess mixing",
}

SURVEILLANCE_TERMS = [
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
]

PRIMARY_TERMS = ["deprivation_z", *SURVEILLANCE_TERMS]
MIXING_TERMS = PRIMARY_TERMS + ["log_cluster_size_z"]
MIXING_PREDICTOR_TERMS = [
    "simd_excess_mixing_z",
    "age_excess_mixing_z",
    "sex_excess_mixing_z",
    "profile_excess_mixing_z",
]

DOMAIN_MIXING_PREDICTOR_ORDER = ["domain_quintile", "age", "sex", "age_sex"]
DOMAIN_MIXING_PREDICTOR_LABELS = {
    "domain_quintile": "Domain quintile",
    "age": "Age",
    "sex": "Sex",
    "age_sex": "Age-sex",
}

MIXING_LABELS = {
    "simd": "SIMD",
    "age": "Age",
    "sex": "Sex",
    "profile": "Joint profile",
    "age_sex": "Joint age-sex",
}

COMPONENT_LABELS = {
    "hurdle_binary": "Hurdle odds",
    "positive_zero_truncated_count": "ZTNB count ratio",
    "log_linear": "Log-linear",
}

DOMAIN_ORDER = [
    "overall",
    "income",
    "employment",
    "education",
    "health",
    "access",
    "crime",
    "housing",
]

DOMAIN_LABELS = {
    "overall": "Overall",
    "income": "Income",
    "employment": "Employment",
    "education": "Education",
    "health": "Health",
    "access": "Access",
    "crime": "Crime",
    "housing": "Housing",
}

WAVE_ORDER = ["B.1.177", "Alpha", "Delta", "BA.1", "BA.2", "BA.4", "BA.5", "BQ.1"]
COUNT_OUTCOMES = ["cluster_size", "geographic_dispersion"]
COUNT_COMPONENTS = ["hurdle_binary", "positive_zero_truncated_count"]
SIZE_ADJUSTED_OUTCOMES = ["geographic_dispersion_size_adjusted"]

# Mixing-predictor components that are estimable. Cluster-size hurdle is
# excluded because its comparison group (singletons) has undefined mixing
# scores; see manuscript Methods and the diagnostics CSV.
MIXING_PREDICTOR_COMPONENTS_BY_OUTCOME = {
    "cluster_size": ["positive_zero_truncated_count"],
    "geographic_dispersion": ["hurdle_binary", "positive_zero_truncated_count"],
}


# ---------------------------------------------------------------------------
# Path / environment helpers
# ---------------------------------------------------------------------------


def repo_root(start: Path | None = None) -> Path:
    p = (start or Path(__file__)).resolve()
    for candidate in [p, *p.parents]:
        if (candidate / "config.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate config.yaml.")


def load_style(root: Path):
    sys.path.insert(0, str(root))
    from utils import style

    return style


def setup_environment() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)


def term_colours(style) -> dict[str, str]:
    palette = style.SIMD_DOMAIN_PALETTE
    return {
        "deprivation_z": palette["overall"],
        "index_deprivation_z": palette["overall"],
        "local_incidence_z": palette["income"],
        "local_seq_fraction_z": palette["employment"],
        "window_seq_fraction_z": palette["education"],
        "test_positivity_z": palette["health"],
        "log_cluster_size_z": palette["access"],
        "simd_excess_mixing_z": palette["crime"],
        "age_excess_mixing_z": palette["housing"],
        "sex_excess_mixing_z": palette["overall"],
        "profile_excess_mixing_z": palette["income"],
        "age_sex_excess_mixing_z": palette["income"],
    }


def primary_terms_for_results(*frames: pd.DataFrame) -> list[str]:
    observed_terms: set[str] = set()
    for frame in frames:
        if "term" in frame:
            observed_terms.update(frame["term"].dropna().astype(str))
    exposure = (
        "index_deprivation_z"
        if "index_deprivation_z" in observed_terms and "deprivation_z" not in observed_terms
        else "deprivation_z"
    )
    return [exposure, *SURVEILLANCE_TERMS]


def mixing_terms_for_results(*frames: pd.DataFrame) -> list[str]:
    return [*primary_terms_for_results(*frames), "log_cluster_size_z"]


def save_all(style, fig, out_base: Path, width: str, height_in: float) -> None:
    style.save_figure(
        fig,
        out_base,
        width=width,
        height_in=height_in,
        dpi=600,
        save_pdf=True,
        save_png=True,
        save_tiff=True,
    )


def write_hurdle_geographic_spread_table(
    results: pd.DataFrame,
    out_path: Path,
    group_cols: list[str],
) -> None:
    """Write a tidy CSV of the dropped hurdle / geographic-spread results.

    The hurdle component of geographic spread is no longer plotted in
    Fig 4 / Supp Fig 7 (SIMD coefficients blow up on the odds scale,
    making the heatmap uninformative). It is reported instead as a
    supplementary table with this function.

    A ``notes`` column flags rows where the window-clustered sandwich
    standard error failed numerically (Hessian effectively singular,
    yielding NaN SE/CI/p-value). The maximum-likelihood point estimates
    (coefficient, ratio) on those rows are still valid.

    The function handles both data schemas: wave-specific results, where the
    relevant mixing predictors are identified by `term`, and SIMD-domain
    results, where each row's "domain-quintile" predictor is named
    `<domain>_domain_excess_mixing_z` and must be matched via the
    `domain_mixing_predictor_key` helper.
    """
    sub = results[
        (results["outcome"] == "geographic_dispersion")
        & (results["component"] == "hurdle_binary")
    ].copy()
    if sub.empty:
        return

    if "domain" in sub.columns:
        sub["predictor"] = sub.apply(domain_mixing_predictor_key, axis=1)
        sub = sub[sub["predictor"].isin(DOMAIN_MIXING_PREDICTOR_ORDER)].copy()
        ordering_col = "predictor"
    else:
        sub = sub[sub["term"].isin(MIXING_PREDICTOR_TERMS)].copy()
        ordering_col = "term"
    if sub.empty:
        return

    # Flag rows where cluster-robust inference failed.
    se_col = "std_error_clustered_by_window"
    se_missing = sub[se_col].isna() if se_col in sub.columns else pd.Series(False, index=sub.index)
    sub["notes"] = np.where(
        se_missing,
        "cluster-robust SE unavailable (Hessian singular); point estimate only",
        "",
    )

    keep_cols = [
        *group_cols,
        "predictor",
        "term",
        "coefficient",
        "std_error_clustered_by_window",
        "ratio",
        "ratio_ci_low",
        "ratio_ci_high",
        "p_value",
        "n_observations",
        "n_events",
        "notes",
    ]
    keep_cols = [c for c in keep_cols if c in sub.columns]
    sort_keys = group_cols + [ordering_col]
    sub = sub[keep_cols].sort_values(sort_keys).reset_index(drop=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(out_path, index=False)


# ---------------------------------------------------------------------------
# SUPPLEMENTARY FIGURE CAPTIONS + TABLES DOCUMENT
# ---------------------------------------------------------------------------


SUPPLEMENTARY_FIGURE_CAPTIONS: list[tuple[str, str, str]] = [
    (
        "supp_fig1_outcome_distributions",
        "Cluster-outcome and excess-mixing distributions among non-singleton clusters",
        "Two-row distributional summary for the 84,067 non-singleton clusters. "
        "Top row: histograms of cluster size, cluster duration, and number of "
        "distinct datazones. Bottom row: histograms of observed-minus-expected "
        "excess mixing for age, sex, and SIMD-deprivation quintile composition. "
        "Modal peaks sit at the structural minima (size 2, duration 0 days, "
        "datazones 1) and medians at 3, 4 days, and 3 respectively. Age and sex "
        "excess mixing centre slightly above zero; SIMD-quintile excess mixing "
        "centres slightly below zero, motivating the regression mixing analyses.",
    ),
    (
        "supp_fig2_mixing_distributions",
        "Excess mixing predictor distributions",
        "Distributions of the four cluster-level mixing predictors used as "
        "explanatory variables in the mixing-predictor count models: SIMD-quintile, "
        "age, sex, and joint age-sex-profile excess-mixing scores. Each is the "
        "observed-minus-expected pair-discordance score for that attribute within "
        "non-singleton clusters, on the percentage-point scale.",
    ),
    (
        "supp_fig3_observed_expected_matrices",
        "Observed-minus-expected pair-probability matrices",
        "Heatmaps of observed-minus-expected pair probabilities for SIMD-quintile "
        "pairs (left) and age-band pairs (right), averaged across all non-singleton "
        "clusters. Cells show mean excess pair probability in percentage points "
        "relative to the lineage- and calendar-window-matched expectation. "
        "Same-quintile SIMD pairs are positive throughout, peaking at quintile 1 × "
        "quintile 1 (most deprived, +0.3 pp). Same-age-band pairs are positive on "
        "the diagonal, peaking among young adults (20-24, +0.2 pp). All cells are "
        "annotated with the numerical excess value (pp).",
    ),
    (
        "supp_fig4_deprivation_domain_outcomes",
        "SIMD-domain deprivation effects on count outcomes",
        "Coefficient plot of the per-1-SD effect of each SIMD subdomain "
        "(income, employment, education, health, access, crime, housing) on the "
        "four count-model components: cluster-size hurdle (odds ratio), positive "
        "cluster size (ZTNB count ratio), geographic-spread hurdle (odds ratio), "
        "and positive geographic spread (ZTNB count ratio). Coefficients are "
        "estimated in domain-specific models adjusting for the surveillance and "
        "lineage-window covariates. Housing and crime show the strongest negative "
        "associations with positive cluster size and positive geographic spread; "
        "access deprivation is the only domain positively associated with the "
        "geographic-spread hurdle and positive geographic spread.",
    ),
    (
        "supp_fig5_deprivation_domain_mixing",
        "SIMD-domain deprivation effects on mixing outcomes",
        "Four-panel coefficient plot of the per-1-SD effect of each SIMD subdomain "
        "on (A) domain-quintile excess mixing, (B) age excess mixing, (C) sex "
        "excess mixing, and (D) joint age-sex profile excess mixing. Mixing "
        "outcomes are observed-minus-expected pair-discordance scores in "
        "percentage points; the x-axis is the adjusted percentage-point change "
        "per 1 SD higher domain deprivation. Education and crime deprivation are "
        "associated with greater domain-quintile mixing; access and housing with "
        "the reverse pattern. Across age, sex, and joint age-sex mixing, access "
        "deprivation behaves opposite to the other six domains.",
    ),
    (
        "supp_fig6_deprivation_domain_wave_mixing",
        "Wave-specific SIMD-domain deprivation effects on demographic mixing",
        "Three stacked heatmaps showing the per-1-SD effect of each SIMD "
        "subdomain on (A) age, (B) sex, and (C) joint age-sex excess mixing for "
        "each epidemic wave. Rows are SIMD subdomains (Overall, Income, "
        "Employment, Education, Health, Access, Crime, Housing); columns are "
        "wave groups (B.1.177, Alpha, Delta, BA.1, BA.2, BA.5). Cells are "
        "coloured by adjusted percentage-point change with shared colour scale "
        "across all three panels (±5 pp), and annotated with the numeric value. "
        "Age mixing is positive for most domain × wave cells in earlier waves, "
        "with access as a consistent negative outlier; sex mixing turns negative "
        "in BA.2 and BA.5 for most domains.",
    ),
    (
        "supp_fig7_mixing_domain_outcomes",
        "Domain-specific mixing-predictor effects on ZTNB cluster outcomes",
        "Two-panel heatmap of the per-1-SD effect of four cluster-level mixing "
        "predictors (domain-quintile, age, sex, age-sex profile excess mixing) on "
        "(A) the positive cluster-size ZTNB count ratio and (B) the positive "
        "geographic-spread ZTNB count ratio, fit separately in each of eight "
        "SIMD subdomains (rows). A single ratio-scale colour bar is shared across "
        "panels and capped at ratio 5 (the largest observed cell ratio in either "
        "panel). The domain-quintile column saturates the upper triangle in both "
        "panels (cell values 2.8-3.4), and all cells are annotated with the raw "
        "ratio. The hurdle component of geographic spread is omitted from this "
        "figure and reported in Supplementary Table 2. Note: between-row "
        "variation is very small because each row represents a separate "
        "per-subdomain model in which the age/sex/age-sex mixing predictors are "
        "identical observed variables; only the domain-specific deprivation and "
        "domain-excess-mixing pair differs across rows, and Scotland's SIMD "
        "subdomains are highly correlated so per-domain coefficients land in "
        "essentially the same place.",
    ),
    (
        "supp_fig8_deprivation_size_adjusted",
        "Size-adjusted positive geographic spread",
        "Coefficient plot comparing SIMD-deprivation and surveillance-covariate "
        "effects on the positive ZTNB geographic-spread model with and without "
        "additional adjustment for log cluster size. The SIMD point estimate "
        "flips direction once cluster size is conditioned on (count ratio "
        "0.851, 95% CI 0.792-0.915 unadjusted; 1.027, 95% CI 1.010-1.044 "
        "size-adjusted), showing that the unadjusted negative geographic-spread "
        "association is explained by deprivation's association with cluster size.",
    ),
    (
        "supp_fig9_deprivation_loglinear",
        "Deprivation log-linear vs hurdle/ZTNB count models",
        "Coefficient plot contrasting SIMD-deprivation effects on cluster size "
        "and geographic spread from a single-component log-linear (Poisson-style) "
        "model with the corresponding hurdle (odds ratio) and ZTNB (count ratio) "
        "components of the two-part main model. Log-linear estimates are "
        "substantially attenuated (cluster size geometric mean ratio 0.992, "
        "geographic spread 1.001) because they average over the structural mass "
        "at the count minimum, masking the within-component associations that "
        "the two-part model separates.",
    ),
    (
        "supp_fig10_mixing_loglinear",
        "Mixing-predictor log-linear vs hurdle/ZTNB count models",
        "Same contrast as Supplementary Figure 9 but for the four mixing "
        "predictors (SIMD-quintile, age, sex, age-sex profile excess mixing). "
        "Log-linear estimates attenuate the mixing-predictor effects on count "
        "outcomes, again because of averaging across the hurdle and "
        "positive-count components of the cluster-size / geographic-spread "
        "distributions.",
    ),
]


SUPPLEMENTARY_TABLE_DESCRIPTIONS: list[tuple[str, str, str, str]] = [
    (
        "supp_table_fig4_wave_mixing_hurdle_geographic_spread.csv",
        "Wave-specific mixing-predictor effects on the geographic-spread hurdle "
        "(companion to Figure 4)",
        "Window-clustered hurdle (binomial GLM with logit link) odds ratios for "
        "the four mixing predictors (SIMD-quintile, age, sex, age-sex profile "
        "excess mixing) by epidemic wave. This component is omitted from "
        "Figure 4 because the SIMD coefficient reaches an odds ratio of ~29,000 "
        "in the Alpha wave (95% CI 6,450-129,898), making a heatmap "
        "uninformative.",
        "wave_group",
    ),
    (
        "supp_table_fig7_domain_mixing_hurdle_geographic_spread.csv",
        "SIMD-domain mixing-predictor effects on the geographic-spread hurdle "
        "(companion to Supplementary Figure 7)",
        "Window-clustered hurdle odds ratios for the four mixing predictors "
        "(domain-quintile, age, sex, age-sex profile excess mixing) in each "
        "per-domain model. The crime and education rows have point estimates "
        "but no cluster-robust SE, CI or p-value because the window-clustered "
        "sandwich variance estimator failed numerically for those two hurdle "
        "fits (Hessian effectively singular under the heavy outcome imbalance, "
        "88% of clusters being multi-datazone). Point estimates remain valid.",
        "domain",
    ),
]


def _format_markdown_row(values: list[str]) -> str:
    return "| " + " | ".join(values) + " |"


def _format_ratio_with_ci(row: pd.Series) -> str:
    ratio = row.get("ratio")
    lo = row.get("ratio_ci_low")
    hi = row.get("ratio_ci_high")
    if pd.isna(ratio):
        return "—"
    if pd.isna(lo) or pd.isna(hi):
        return f"{ratio:.3g} (CI unavailable)"
    if abs(ratio) >= 100:
        return f"{ratio:.0f} ({lo:.0f}–{hi:.0f})"
    return f"{ratio:.3g} ({lo:.3g}–{hi:.3g})"


def _format_p_value(p: float) -> str:
    if pd.isna(p):
        return "—"
    if p == 0:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.1e}"
    return f"{p:.3g}"


def _hurdle_table_to_markdown(df: pd.DataFrame, group_col: str) -> list[str]:
    """Render a hurdle-results dataframe as a compact markdown table."""
    pred_labels = {
        **{t: MIXING_LABELS[t.replace("_excess_mixing_z", "")] for t in MIXING_PREDICTOR_TERMS},
        **{k: DOMAIN_MIXING_PREDICTOR_LABELS[k] for k in DOMAIN_MIXING_PREDICTOR_ORDER},
    }
    group_labels = {
        "domain": DOMAIN_LABELS,
    }
    has_predictor_col = "predictor" in df.columns
    header = [group_col.replace("_", " ").title(), "Mixing predictor", "Ratio (95% CI)", "p", "Notes"]
    lines = [_format_markdown_row(header)]
    lines.append(_format_markdown_row(["---"] * len(header)))
    for _, row in df.iterrows():
        if has_predictor_col and isinstance(row["predictor"], str):
            label = pred_labels.get(row["predictor"], row["predictor"])
        else:
            label = pred_labels.get(row["term"], row["term"])
        notes = row.get("notes") if "notes" in df.columns else ""
        if isinstance(notes, float) and pd.isna(notes):
            notes = ""
        group_val = row[group_col]
        if group_col in group_labels:
            group_val = group_labels[group_col].get(group_val, group_val)
        lines.append(
            _format_markdown_row(
                [
                    str(group_val),
                    str(label),
                    _format_ratio_with_ci(row),
                    _format_p_value(row.get("p_value", float("nan"))),
                    str(notes),
                ]
            )
        )
    return lines


def write_supplementary_files_document(out_dir: Path) -> None:
    """Write a single markdown file with captions for every supplementary
    figure and the supplementary tables (with the dropped hurdle results)."""
    lines: list[str] = []
    lines.append("# Part 1 Supplementary Files")
    lines.append("")
    lines.append(
        "Companion document to the Part 1 manuscript. Contains the figure "
        "captions for every supplementary figure produced by `make_figures.py`, "
        "and the supplementary tables for hurdle geographic-spread results that "
        "were dropped from Figure 4 and Supplementary Figure 7 (those are now "
        "ZTNB-only heatmaps). Figures are saved as PDF, PNG, and TIFF in this "
        "directory."
    )
    lines.append("")

    # Figures.
    lines.append("## Supplementary Figures")
    lines.append("")
    for idx, (stem, title, caption) in enumerate(SUPPLEMENTARY_FIGURE_CAPTIONS, start=1):
        figure_path = out_dir / f"{stem}.png"
        if not figure_path.exists():
            continue
        lines.append(f"### Supplementary Figure {idx}: {title}")
        lines.append("")
        lines.append(f"**File:** `{stem}` (PDF/PNG/TIFF)")
        lines.append("")
        lines.append(caption)
        lines.append("")

    # Tables.
    lines.append("## Supplementary Tables")
    lines.append("")
    for table_idx, (filename, title, blurb, group_col) in enumerate(
        SUPPLEMENTARY_TABLE_DESCRIPTIONS, start=1
    ):
        table_path = out_dir / filename
        if not table_path.exists():
            continue
        df = pd.read_csv(table_path)
        lines.append(f"### Supplementary Table {table_idx}: {title}")
        lines.append("")
        lines.append(f"**File:** `{filename}`")
        lines.append("")
        lines.append(blurb)
        lines.append("")
        lines.extend(_hurdle_table_to_markdown(df, group_col=group_col))
        lines.append("")
        # Note the columns omitted from the markdown view.
        lines.append(
            "*Full coefficient, standard error, z, n_observations, and "
            "n_events columns are in the companion CSV.*"
        )
        lines.append("")

    out_path = out_dir / "part1_supplementary_files.md"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Drawing primitives
# ---------------------------------------------------------------------------


def draw_ratio_panel(
    ax,
    df: pd.DataFrame,
    terms: list[str],
    colours: dict[str, str],
    *,
    title: str,
    show_ylabels: bool,
    xlim: tuple[float, float],
    xlabel: str | None = None,
    ticks: list[float] | None = None,
) -> None:
    from matplotlib.ticker import NullFormatter, NullLocator

    y_positions = np.arange(len(terms))[::-1]
    position = dict(zip(terms, y_positions))
    for term in terms:
        row = df[df["term"] == term]
        if row.empty:
            continue
        row = row.iloc[0]
        y = position[term]
        ax.plot(
            [row["ratio_ci_low"], row["ratio_ci_high"]],
            [y, y],
            color=colours[term],
            linewidth=1.1,
            solid_capstyle="round",
        )
        ax.scatter(
            row["ratio"],
            y,
            color=colours[term],
            edgecolor="white",
            linewidth=0.3,
            s=18,
            zorder=3,
        )

    ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    if ticks is None:
        if xlim[1] >= 3.9:
            ticks = [0.8, 1.0, 1.5, 2.0, 3.0, 4.0]
        else:
            ticks = [0.9, 1.0, 1.5, 2.0, 3.0]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{tick:g}" for tick in ticks])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_title(title, pad=4)
    ax.set_ylim(-0.7, len(terms) - 0.3)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([TERM_LABELS[t] for t in terms] if show_ylabels else [])
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)


def draw_difference_panel(
    ax,
    df: pd.DataFrame,
    terms: list[str],
    colours: dict[str, str],
    *,
    title: str,
    show_ylabels: bool,
    xlim: tuple[float, float],
    xlabel: str | None = None,
) -> None:
    y_positions = np.arange(len(terms))[::-1]
    position = dict(zip(terms, y_positions))
    for term in terms:
        row = df[df["term"] == term]
        if row.empty:
            continue
        row = row.iloc[0]
        y = position[term]
        ax.plot(
            [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
            [y, y],
            color=colours[term],
            linewidth=1.1,
            solid_capstyle="round",
        )
        ax.scatter(
            row["coefficient_percentage_points"],
            y,
            color=colours[term],
            edgecolor="white",
            linewidth=0.3,
            s=18,
            zorder=3,
        )

    ax.axvline(0.0, color="#666666", linestyle="--", linewidth=0.7)
    ax.set_xlim(*xlim)
    ax.set_title(title, pad=4)
    ax.set_ylim(-0.7, len(terms) - 0.3)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([TERM_LABELS[t] for t in terms] if show_ylabels else [])
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)


def domain_effect_rows(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["term"].astype(str).eq(df["domain"].astype(str) + "_deprivation_z")].copy()


def domain_mixing_predictor_key(row: pd.Series) -> str | None:
    term = str(row["term"])
    domain = str(row["domain"])
    if term == f"{domain}_domain_excess_mixing_z":
        return "domain_quintile"
    if term == "age_excess_mixing_z":
        return "age"
    if term == "sex_excess_mixing_z":
        return "sex"
    if term == "age_sex_excess_mixing_z":
        return "age_sex"
    return None


# ---------------------------------------------------------------------------
# MAIN FIGURE 1 — Line 1 (deprivation) overall
# ---------------------------------------------------------------------------


def plot_deprivation_overall(
    style,
    count_results: pd.DataFrame,
    mixing_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Single combined figure: deprivation effects on counts (top) and mixing (bottom).

    Top row (4 panels): cluster size hurdle OR, cluster size positive CR,
    geographic spread hurdle OR, geographic spread positive CR. Each panel
    shows deprivation + four surveillance covariates per 1 SD higher covariate.

    Bottom row (4 panels): SIMD-quintile, age-band, sex, joint
    SIMD-age-sex profile excess discordance (pp). Each panel shows deprivation
    + four surveillance covariates + log cluster size.
    """
    import matplotlib.pyplot as plt

    colours = term_colours(style)
    primary_terms = primary_terms_for_results(count_results)
    mixing_terms = mixing_terms_for_results(mixing_results)
    count_outcomes = COUNT_OUTCOMES
    count_components = COUNT_COMPONENTS
    mixing_outcomes = ["simd", "age", "sex", "profile"]

    fig, axes = style.new_figure(
        width="double",
        height_in=6.4,
        nrows=2,
        ncols=4,
        sharex=False,
        font_scale=0.82,
    )

    # Top row — count outcomes.
    count_xlim = (0.75, 4.0)
    count_ticks = [0.8, 1.0, 1.5, 2.0, 3.0, 4.0]
    panel_idx = 0
    for i_outcome, outcome in enumerate(count_outcomes):
        for j_component, component in enumerate(count_components):
            ax = axes[0, panel_idx]
            sub = count_results[
                (count_results["outcome"] == outcome)
                & (count_results["component"] == component)
                & (count_results["term"].isin(primary_terms))
            ]
            draw_ratio_panel(
                ax,
                sub,
                primary_terms,
                colours,
                title=f"{OUTCOME_LABELS[outcome]}\n{COMPONENT_LABELS[component]}",
                show_ylabels=(panel_idx == 0),
                xlim=count_xlim,
                xlabel=(
                    "Odds ratio"
                    if component == "hurdle_binary"
                    else "Count ratio"
                ),
                ticks=count_ticks,
            )
            panel_idx += 1

    # Bottom row — mixing outcomes.
    mix_ci_min = float(mixing_results["ci_low_percentage_points"].min())
    mix_ci_max = float(mixing_results["ci_high_percentage_points"].max())
    mix_limit = max(2.0, np.ceil(max(abs(mix_ci_min), abs(mix_ci_max)) * 1.05))
    mixing_xlim = (-mix_limit, mix_limit)

    for idx, outcome in enumerate(mixing_outcomes):
        ax = axes[1, idx]
        sub = mixing_results[mixing_results["outcome"] == outcome]
        draw_difference_panel(
            ax,
            sub,
            mixing_terms,
            colours,
            title=f"{MIXING_LABELS[outcome]} excess mixing",
            show_ylabels=(idx == 0),
            xlim=mixing_xlim,
            xlabel="pp per 1 SD",
        )

    style.add_panel_labels(axes.ravel(), x=-0.18, y=1.18, size=9)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.92, bottom=0.10, hspace=0.55, wspace=0.20)
    save_all(style, fig, out_dir / "fig1_deprivation_overall", "double", 6.4)
    plt.close("all")


# ---------------------------------------------------------------------------
# MAIN FIGURE 2 — Line 1 (deprivation) by epidemic wave
# ---------------------------------------------------------------------------


def plot_deprivation_wave_specific(
    style,
    wave_count_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Wave-specific deprivation effects on the four count-model components."""
    from matplotlib.ticker import NullFormatter, NullLocator
    import matplotlib.pyplot as plt

    colour = style.SIMD_DOMAIN_PALETTE["overall"]
    outcomes = COUNT_OUTCOMES
    components = COUNT_COMPONENTS
    data = wave_count_results[
        wave_count_results["term"].eq("deprivation_z")
        & wave_count_results["outcome"].isin(outcomes)
        & wave_count_results["component"].isin(components)
    ].copy()
    waves = [wave for wave in WAVE_ORDER if wave in set(data["wave_group"])]
    if data.empty or not waves:
        return

    xlims: dict[str, tuple[float, float]] = {}
    ticks: dict[str, list[float]] = {}
    for component in components:
        sub = data[data["component"] == component]
        ci_min = float(sub["ratio_ci_low"].min())
        ci_max = float(sub["ratio_ci_high"].max())
        lower = max(0.5, np.floor(ci_min * 10) / 10)
        upper = min(3.0, np.ceil(ci_max * 10) / 10)
        if component == "hurdle_binary":
            lower = min(lower, 0.8)
            upper = max(upper, 1.2)
            tick_candidates = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
        else:
            lower = min(lower, 0.6)
            upper = max(upper, 2.7)
            tick_candidates = [0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
        xlims[component] = (lower, upper)
        ticks[component] = [tick for tick in tick_candidates if lower <= tick <= upper]

    fig, axes = style.new_figure(
        width="double",
        height_in=4.8,
        nrows=len(outcomes),
        ncols=2,
        sharex=False,
        font_scale=0.85,
    )
    y_positions = np.arange(len(waves))[::-1]
    pos = dict(zip(waves, y_positions))

    for idx, outcome in enumerate(outcomes):
        for jdx, component in enumerate(components):
            ax = axes[idx, jdx]
            sub = data[(data["outcome"] == outcome) & (data["component"] == component)]
            for wave in waves:
                row = sub[sub["wave_group"] == wave]
                if row.empty:
                    continue
                row = row.iloc[0]
                y = pos[wave]
                ax.plot(
                    [row["ratio_ci_low"], row["ratio_ci_high"]],
                    [y, y],
                    color=colour,
                    linewidth=1.1,
                    solid_capstyle="round",
                )
                ax.scatter(
                    row["ratio"],
                    y,
                    color=colour,
                    edgecolor="white",
                    linewidth=0.3,
                    s=20,
                    zorder=3,
                )
            ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
            ax.set_xscale("log")
            ax.set_xlim(*xlims[component])
            ax.set_xticks(ticks[component])
            ax.set_xticklabels([f"{tick:g}" for tick in ticks[component]])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_title(f"{OUTCOME_LABELS[outcome]}: {COMPONENT_LABELS[component]}", pad=4)
            ax.set_yticks(y_positions)
            ax.set_yticklabels(waves if jdx == 0 else [])
            if idx == len(outcomes) - 1:
                ax.set_xlabel(
                    "Odds ratio per 1 SD higher SIMD deprivation"
                    if component == "hurdle_binary"
                    else "ZTNB count ratio per 1 SD higher SIMD deprivation"
                )
            ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes.ravel(), x=-0.08, y=1.15, size=9)
    fig.subplots_adjust(left=0.13, right=0.99, top=0.93, bottom=0.10, hspace=0.42, wspace=0.13)
    save_all(style, fig, out_dir / "fig2_deprivation_wave_specific", "double", 4.8)
    plt.close("all")


# ---------------------------------------------------------------------------
# MAIN FIGURE 3 — Line 2 (mixing predictor) overall
# ---------------------------------------------------------------------------


def plot_mixing_overall(
    style,
    mixing_predictor_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Mixing-predictor effects on the three estimable count-model components.

    Cluster-size hurdle is omitted because mixing scores are undefined for
    singletons and the hurdle is not estimable (see Methods).

    Panel layout (1 row × 3 cols):
        A — Cluster size, positive ZTNB count ratio
        B — Geographic spread, hurdle OR
        C — Geographic spread, positive ZTNB count ratio
    """
    from matplotlib.ticker import NullFormatter, NullLocator
    import matplotlib.pyplot as plt

    colours = term_colours(style)
    data = mixing_predictor_results[
        mixing_predictor_results["term"].isin(MIXING_PREDICTOR_TERMS)
    ].copy()
    if data.empty:
        return

    panels: list[tuple[str, str, str]] = [
        ("cluster_size", "positive_zero_truncated_count", "Cluster size: ZTNB count ratio"),
        ("geographic_dispersion", "hurdle_binary", "Geographic spread: Hurdle odds"),
        ("geographic_dispersion", "positive_zero_truncated_count", "Geographic spread: ZTNB count ratio"),
    ]

    fig, axes = style.new_figure(
        width="double",
        height_in=3.4,
        nrows=1,
        ncols=3,
        sharex=False,
        font_scale=0.85,
    )

    for idx, (outcome, component, title) in enumerate(panels):
        ax = axes[idx]
        sub = data[(data["outcome"] == outcome) & (data["component"] == component)]
        if sub.empty:
            ax.text(0.5, 0.5, "Not estimable", transform=ax.transAxes, ha="center", va="center", color="#666666")
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        ci_min = float(sub["ratio_ci_low"].min())
        ci_max = float(sub["ratio_ci_high"].max())
        lower = max(0.5, np.floor(ci_min * 10.0) / 10.0)
        upper = max(1.5, np.ceil(ci_max * 10.0) / 10.0)
        lower = min(lower, 0.8)
        if upper >= 10:
            tick_candidates = [0.8, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0]
        elif upper >= 5:
            tick_candidates = [0.8, 1.0, 1.5, 2.0, 3.0, 5.0]
        else:
            tick_candidates = [0.8, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0]
        ticks = [tick for tick in tick_candidates if lower <= tick <= upper]
        if not ticks:
            ticks = [1.0]

        y_positions = np.arange(len(MIXING_PREDICTOR_TERMS))[::-1]
        pos = dict(zip(MIXING_PREDICTOR_TERMS, y_positions))
        for term in MIXING_PREDICTOR_TERMS:
            row = sub[sub["term"] == term]
            if row.empty:
                continue
            row = row.iloc[0]
            y = pos[term]
            ax.plot(
                [row["ratio_ci_low"], row["ratio_ci_high"]],
                [y, y],
                color=colours[term],
                linewidth=1.1,
                solid_capstyle="round",
            )
            ax.scatter(
                row["ratio"],
                y,
                color=colours[term],
                edgecolor="white",
                linewidth=0.3,
                s=18,
                zorder=3,
            )

        ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_xscale("log")
        ax.set_xlim(lower, upper)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{tick:g}" for tick in ticks])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_title(title, pad=4)
        ax.set_ylim(-0.7, len(MIXING_PREDICTOR_TERMS) - 0.3)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(
            [TERM_LABELS[t] for t in MIXING_PREDICTOR_TERMS] if idx == 0 else []
        )
        ax.set_xlabel(
            "Odds ratio per 1 SD\nhigher excess mixing"
            if component == "hurdle_binary"
            else "Count ratio per 1 SD\nhigher excess mixing"
        )
        ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes, x=-0.10, y=1.15, size=9)
    fig.subplots_adjust(left=0.20, right=0.98, top=0.91, bottom=0.20, wspace=0.18)
    save_all(style, fig, out_dir / "fig3_mixing_overall", "double", 3.4)
    plt.close("all")


# ---------------------------------------------------------------------------
# MAIN FIGURE 4 — Line 2 (mixing predictor) by epidemic wave
# ---------------------------------------------------------------------------


def plot_mixing_wave_specific(
    style,
    wave_mixing_predictor_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Wave-specific mixing-predictor effects as heatmap (ZTNB only).

    Two panels (ZTNB cluster-size count ratio and ZTNB geographic-spread count
    ratio), each a heatmap with waves as rows and the four mixing predictors as
    columns. The hurdle geographic-spread component is reported as a
    supplementary table because its SIMD coefficient blows up to ratio ~30 000
    in the Alpha wave and is uninformative on a heatmap.
    """
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.pyplot as plt

    data = wave_mixing_predictor_results[
        wave_mixing_predictor_results["term"].isin(MIXING_PREDICTOR_TERMS)
    ].copy()
    waves = [wave for wave in WAVE_ORDER if wave in set(data["wave_group"])]
    if data.empty or not waves:
        return

    panels: list[tuple[str, str, str]] = [
        ("cluster_size", "positive_zero_truncated_count", "Cluster size"),
        ("geographic_dispersion", "positive_zero_truncated_count", "Geographic spread"),
    ]
    values = data[(data["outcome"].isin([panels[0][0], panels[1][0]]))
               & (data["component"].isin([panels[0][1], panels[1][1]]))]
    vmax = max(2, float(np.nanmax(np.abs(values["ratio"]))))
    vmin = min(0, float(np.nanmin(np.abs(values["ratio"]))))
    vcenter = 1
    norm = TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

    fig, axes = style.new_figure(
        width="double",
        height_in=4.0,
        nrows=1,
        ncols=2,
        font_scale=0.80,
    )
    fig.subplots_adjust(left=0.10, right=0.86, top=0.88, bottom=0.22, wspace=0.12)

    image = None
    for idx, (outcome, component, title) in enumerate(panels):
        ax = axes[idx]
        sub = data[(data["outcome"] == outcome) & (data["component"] == component)]
        if sub.empty:
            matrix = pd.DataFrame(np.nan, index=waves, columns=MIXING_PREDICTOR_TERMS)
        else:
            matrix = (
                sub.pivot_table(
                    index="wave_group",
                    columns="term",
                    values="ratio",
                    aggfunc="first",
                )
                .reindex(index=waves, columns=MIXING_PREDICTOR_TERMS)
            )
        image = ax.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", norm=norm, aspect="auto")
        ax.set_title(title, pad=6)
        ax.set_yticks(np.arange(len(waves)))
        ax.set_yticklabels(waves if idx == 0 else [])
        ax.set_xticks(np.arange(len(MIXING_PREDICTOR_TERMS)))
        ax.set_xticklabels(
            [MIXING_LABELS[t.replace("_excess_mixing_z", "")] for t in MIXING_PREDICTOR_TERMS],
            rotation=35,
            ha="right",
        )
        ax.tick_params(length=0)
        for y in np.arange(len(waves) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.6)
        for x in np.arange(len(MIXING_PREDICTOR_TERMS) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.6)

    assert image is not None
    cbar_ax = fig.add_axes([0.885, 0.24, 0.020, 0.58])
    cbar = fig.colorbar(image, cax=cbar_ax, extend="both")
    cbar.set_label("ZTNB count ratio per 1 SD higher excess mixing", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)
    style.add_panel_labels(axes, x=-0.10, y=1.10, size=9)
    save_all(style, fig, out_dir / "fig4_mixing_wave_specific", "double", 4.0)
    plt.close("all")

# ---------------------------------------------------------------------------
# SUPPLEMENTARY FIGURES
# ---------------------------------------------------------------------------

def binned_percent(values: pd.Series, bins: list[float], labels: list[str]) -> pd.DataFrame:
    cats = pd.cut(values, bins=bins, labels=labels, include_lowest=True, right=True)
    pct = cats.value_counts(sort=False, normalize=True).mul(100)
    out = pct.rename("percent").reset_index()
    return out.rename(columns={out.columns[0]: "bin"})


def plot_outcome_distributions(style, cluster_table: pd.DataFrame, out_dir: Path) -> None:
    """Cluster size, duration, and geographic spread among non-singletons."""
    import matplotlib.pyplot as plt

    grey = "#6f6f6f"
    non_singleton = cluster_table.loc[cluster_table["cluster_size"] > 1].copy()
    if non_singleton.empty:
        return

    count_specs = [
        (
            "cluster_size",
            "Cluster size",
            [-np.inf, 2.5, 3.5, 5.5, 10.5, 20.5, 50.5, np.inf],
            ["2", "3", "4-5", "6-10", "11-20", "21-50", ">50"],
        ),
        (
            "duration_days",
            "Duration (days)",
            [-np.inf, 0.5, 1.5, 2.5, 5.5, 10.5, 15.5, np.inf],
            ["0", "1", "2", "3-5", "6-10", "11-15", ">15"],
        ),
        (
            "cluster_n_datazones",
            "Distinct datazones",
            [-np.inf, 1.5, 2.5, 3.5, 5.5, 10.5, 20.5, 50.5, np.inf],
            ["1", "2", "3", "4-5", "6-10", "11-20", "21-50", ">50"],
        ),
    ]

    fig, axes = style.new_figure(
        width="double",
        height_in=2.8,
        nrows=1,
        ncols=3,
        font_scale=0.85,
    )
    for ax, (col, title, bins, labels) in zip(axes.ravel(), count_specs):
        data = binned_percent(non_singleton[col], bins, labels)
        ax.bar(data["bin"].astype(str), data["percent"], color=grey, width=0.78)
        ax.set_title(title, pad=4)
        ax.set_ylabel("Clusters (%)" if ax is axes.ravel()[0] else "")
        ax.set_ylim(0, max(30, data["percent"].max() * 1.15))
        ax.tick_params(axis="x", rotation=45)
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes.ravel(), x=-0.18, y=1.16, size=9)
    fig.subplots_adjust(left=0.09, right=0.99, top=0.90, bottom=0.20, wspace=0.30)
    save_all(style, fig, out_dir / "supp_fig1_outcome_distributions", "double", 2.8)
    plt.close("all")


def plot_mixing_distributions(style, cluster_table: pd.DataFrame, out_dir: Path) -> None:
    """Distributions of excess discordance for SIMD, age, sex, and joint profile mixing."""
    import matplotlib.pyplot as plt

    grey = "#6f6f6f"
    non_singleton = cluster_table.loc[cluster_table["cluster_size"] > 1].copy()
    if non_singleton.empty:
        return

    mixing_specs = [
        ("simd_excess_discordance", "SIMD quintile"),
        ("age_excess_discordance", "Age band"),
        ("sex_excess_discordance", "Sex"),
    ]
    if "profile_excess_discordance" in non_singleton.columns:
        mixing_specs.append(("profile_excess_discordance", "Joint profile"))

    mixing_bins = np.arange(-100, 101, 10)

    def histogram_percent(values: pd.Series) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        clean = values.dropna().to_numpy(dtype=float)
        counts, edges = np.histogram(clean, bins=mixing_bins)
        total = counts.sum()
        percent = counts / total * 100 if total else counts.astype(float)
        centres = (edges[:-1] + edges[1:]) / 2
        widths = np.diff(edges)
        return centres, widths, percent

    n_panels = len(mixing_specs)
    fig, axes = style.new_figure(
        width="double",
        height_in=2.8,
        nrows=1,
        ncols=n_panels,
        font_scale=0.85,
    )
    flat_axes = axes.ravel() if hasattr(axes, "ravel") else [axes]
    for ax, (col, title) in zip(flat_axes, mixing_specs):
        centres, widths, percent = histogram_percent(non_singleton[col] * 100)
        ax.bar(centres, percent, width=widths * 0.92, color=grey, align="center")
        ax.axvline(0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_title(title, pad=4)
        ax.set_xlabel("Excess mixing (pp)")
        ax.set_ylabel("Clusters (%)" if ax is flat_axes[0] else "")
        ax.set_xlim(mixing_bins[0], mixing_bins[-1])
        ax.set_xticks([-100, -50, 0, 50, 100])
        ax.set_ylim(0, max(30, percent.max() * 1.15))
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(flat_axes, x=-0.18, y=1.16, size=9)
    fig.subplots_adjust(left=0.09, right=0.99, top=0.90, bottom=0.20, wspace=0.32)
    save_all(style, fig, out_dir / "supp_fig2_mixing_distributions", "double", 2.8)
    plt.close("all")


def plot_observed_expected_matrices(
    style,
    matrices: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Observed-minus-expected pair probability matrices for SIMD and age."""
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.pyplot as plt

    def category_key(value: object) -> tuple[int, str]:
        text = str(value)
        if text == "75+":
            return (75, text)
        digits = "".join(ch if ch.isdigit() else " " for ch in text).split()
        return (int(digits[0]) if digits else 999, text)

    specs = [
        ("simd", "SIMD quintile"),
        ("age", "Age band"),
    ]
    fig, axes = style.new_figure(
        width="double",
        height_in=3.7,
        nrows=1,
        ncols=2,
        font_scale=0.8,
    )
    overall = matrices[matrices["wave_group"] == "Overall"].copy()
    vmax = max(0.1, float(np.nanmax(np.abs(overall["excess_percentage_points"]))))
    vmin = min(-0.1, float(np.nanmin(np.abs(overall["excess_percentage_points"]))))
    vcenter = 0
    norm = TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)
    image = None
    for idx, (variable, title) in enumerate(specs):
        ax = axes[idx]
        sub = overall[overall["variable"] == variable]
        row_order = sorted(sub["category_i"].astype(str).unique(), key=category_key)
        col_order = sorted(sub["category_j"].astype(str).unique(), key=category_key)
        matrix = (
            sub.assign(
                category_i=sub["category_i"].astype(str),
                category_j=sub["category_j"].astype(str),
            )
            .pivot_table(
                index="category_i",
                columns="category_j",
                values="excess_percentage_points",
                aggfunc="first",
            )
            .reindex(index=row_order, columns=col_order)
        )
        image = ax.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", norm=norm, aspect="equal")
        ax.set_title(f"{title}: observed - expected", pad=4)
        ax.set_xticks(np.arange(len(col_order)))
        ax.set_yticks(np.arange(len(row_order)))
        if idx == 1:
            ax.set_xticklabels(col_order, rotation=45, ha="right")
        else:
            ax.set_xticklabels(col_order)
        ax.set_yticklabels(row_order)
        ax.set_xlabel(title)
        ax.set_ylabel(title if idx == 0 else "")
        ax.tick_params(length=0)
        for y in np.arange(len(row_order) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.5)
        for x in np.arange(len(col_order) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.5)

    assert image is not None
    fig.subplots_adjust(left=0.08, right=0.84, top=0.87, bottom=0.21, wspace=0.18)
    cbar_ax = fig.add_axes([0.875, 0.28, 0.024, 0.50])
    cbar = fig.colorbar(image, cax=cbar_ax, extend="both")
    cbar.set_label("Observed - expected pair probability (pp)")
    style.add_panel_labels(axes, x=-0.1, y=1.12, size=9)
    save_all(style, fig, out_dir / "supp_fig3_observed_expected_matrices", "double", 3.7)
    plt.close("all")


def plot_deprivation_domain_outcomes(
    style,
    domain_outcomes: pd.DataFrame,
    out_dir: Path,
) -> None:
    """SIMD-domain deprivation effects on count outcomes."""
    from matplotlib.ticker import NullFormatter, NullLocator
    import matplotlib.pyplot as plt

    colours = style.SIMD_DOMAIN_PALETTE
    outcomes = COUNT_OUTCOMES
    components = COUNT_COMPONENTS
    data = domain_effect_rows(domain_outcomes)
    ci_min = float(data["ratio_ci_low"].min())
    ci_max = float(data["ratio_ci_high"].max())
    xlim = (max(0.5, ci_min * 0.95), min(1.5, ci_max * 1.05))
    if xlim[0] > 0.9:
        xlim = (0.9, xlim[1])
    if xlim[1] < 1.1:
        xlim = (xlim[0], 1.1)

    fig, axes = style.new_figure(
        width="double",
        height_in=4.8,
        nrows=len(outcomes),
        ncols=2,
        sharex=True,
        font_scale=0.85,
    )
    y_positions = np.arange(len(DOMAIN_ORDER))[::-1]
    pos = dict(zip(DOMAIN_ORDER, y_positions))
    for idx, outcome in enumerate(outcomes):
        for jdx, component in enumerate(components):
            ax = axes[idx, jdx]
            sub = data[(data["outcome"] == outcome) & (data["component"] == component)]
            for domain in DOMAIN_ORDER:
                row = sub[sub["domain"] == domain]
                if row.empty:
                    continue
                row = row.iloc[0]
                y = pos[domain]
                ax.plot(
                    [row["ratio_ci_low"], row["ratio_ci_high"]],
                    [y, y],
                    color=colours[domain],
                    linewidth=1.0,
                    solid_capstyle="round",
                )
                ax.scatter(
                    row["ratio"],
                    y,
                    color=colours[domain],
                    edgecolor="white",
                    linewidth=0.3,
                    s=18,
                    zorder=3,
                )
            ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
            ax.set_xscale("log")
            ax.set_xlim(*xlim)
            ticks = [tick for tick in [0.8, 0.9, 1.0, 1.1, 1.2] if xlim[0] <= tick <= xlim[1]]
            ax.set_xticks(ticks)
            ax.set_xticklabels([f"{tick:g}" for tick in ticks])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_title(f"{OUTCOME_LABELS[outcome]}: {COMPONENT_LABELS[component]}", pad=4)
            ax.set_yticks(y_positions)
            ax.set_yticklabels([DOMAIN_LABELS[d] for d in DOMAIN_ORDER] if jdx == 0 else [])
            ax.set_xlabel(
                (
                    "Odds ratio per 1 SD higher domain deprivation"
                    if component == "hurdle_binary"
                    else "ZTNB count ratio per 1 SD higher domain deprivation"
                )
                if idx == len(outcomes) - 1
                else ""
            )
            ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes.ravel(), x=-0.08, y=1.15, size=9)
    fig.subplots_adjust(left=0.18, right=0.99, top=0.93, bottom=0.09, hspace=0.42, wspace=0.12)
    save_all(style, fig, out_dir / "supp_fig4_deprivation_domain_outcomes", "double", 4.8)
    plt.close("all")


def plot_deprivation_domain_mixing(
    style,
    domain_mixing: pd.DataFrame,
    domain_demo: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Domain deprivation effects on SIMD-quintile, age, sex, joint mixing."""
    import matplotlib.pyplot as plt

    colours = style.SIMD_DOMAIN_PALETTE
    domain_quintile = domain_effect_rows(domain_mixing).copy()
    domain_quintile["panel"] = "domain_quintile"
    demo = domain_effect_rows(domain_demo).copy()
    demo["panel"] = demo["mixing"]
    data = pd.concat([domain_quintile, demo], ignore_index=True)
    panels = [
        ("domain_quintile", "Domain-quintile mixing"),
        ("age", "Age mixing"),
        ("sex", "Sex mixing"),
        ("age_sex", "Joint age-sex mixing"),
    ]
    ci_min = float(data["ci_low_percentage_points"].min())
    ci_max = float(data["ci_high_percentage_points"].max())
    limit = max(2.5, np.ceil(max(abs(ci_min), abs(ci_max)) * 2) / 2)
    xlim = (-limit, limit)

    fig, axes = style.new_figure(
        width="double",
        height_in=5.0,
        nrows=2,
        ncols=2,
        sharex=True,
        font_scale=0.85,
    )
    y_positions = np.arange(len(DOMAIN_ORDER))[::-1]
    pos = dict(zip(DOMAIN_ORDER, y_positions))

    for idx, (panel, title) in enumerate(panels):
        ax = axes.ravel()[idx]
        sub = data[data["panel"] == panel]
        for domain in DOMAIN_ORDER:
            row = sub[sub["domain"] == domain]
            if row.empty:
                continue
            row = row.iloc[0]
            y = pos[domain]
            ax.plot(
                [row["ci_low_percentage_points"], row["ci_high_percentage_points"]],
                [y, y],
                color=colours[domain],
                linewidth=1.1,
                solid_capstyle="round",
            )
            ax.scatter(
                row["coefficient_percentage_points"],
                y,
                color=colours[domain],
                edgecolor="white",
                linewidth=0.3,
                s=19,
                zorder=3,
            )
        ax.axvline(0.0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_xlim(*xlim)
        ax.set_title(title, pad=4)
        ax.set_yticks(y_positions)
        ax.set_yticklabels([DOMAIN_LABELS[d] for d in DOMAIN_ORDER] if idx % 2 == 0 else [])
        ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    style.add_panel_labels(axes.ravel(), x=-0.08, y=1.15, size=9)
    fig.supxlabel(
        "Change in excess mixing (pp per 1 SD higher domain deprivation)",
        y=0.04,
        fontsize=8,
    )
    fig.subplots_adjust(left=0.18, right=0.99, top=0.90, bottom=0.14, hspace=0.32, wspace=0.12)
    save_all(style, fig, out_dir / "supp_fig5_deprivation_domain_mixing", "double", 5.0)
    plt.close("all")


def plot_deprivation_domain_wave_mixing(
    style,
    wave_domain_demo: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Per-wave heatmap of domain deprivation effects on demographic mixing."""
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.pyplot as plt

    data = domain_effect_rows(wave_domain_demo)
    mixings = ["age", "sex", "age_sex"]
    waves = [wave for wave in WAVE_ORDER if wave in set(data["wave_group"])]
    fig, axes = style.new_figure(
        width="double",
        height_in=6.2,
        nrows=3,
        ncols=1,
        font_scale=0.82,
    )
    vmax = max(1, float(np.nanmax(np.abs(data["coefficient_percentage_points"]))))
    vmin = min(-1, float(np.nanmin(np.abs(data["coefficient_percentage_points"]))))
    vcenter = 0
    norm = TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

    image = None
    for idx, mixing in enumerate(mixings):
        ax = axes[idx]
        sub = data[data["mixing"] == mixing]
        matrix = (
            sub.pivot_table(
                index="domain",
                columns="wave_group",
                values="coefficient_percentage_points",
                aggfunc="first",
            )
            .reindex(index=DOMAIN_ORDER, columns=waves)
        )
        image = ax.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", norm=norm, aspect="auto")
        ax.set_title(f"{MIXING_LABELS[mixing]} mixing", pad=4)
        ax.set_yticks(np.arange(len(DOMAIN_ORDER)))
        ax.set_yticklabels([DOMAIN_LABELS[d] for d in DOMAIN_ORDER])
        ax.set_xticks(np.arange(len(waves)))
        ax.set_xticklabels(waves, ha="center")
        ax.tick_params(length=0)
        for y in np.arange(len(DOMAIN_ORDER) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.6)
        for x in np.arange(len(waves) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.6)

    assert image is not None
    fig.subplots_adjust(left=0.17, right=0.84, top=0.93, bottom=0.12, hspace=0.46)
    cbar_ax = fig.add_axes([0.875, 0.20, 0.022, 0.62])
    cbar = fig.colorbar(image, cax=cbar_ax, extend="both")
    cbar.set_label("pp per 1 SD higher domain deprivation")
    style.add_panel_labels(axes, x=-0.1, y=1.10, size=9)
    save_all(style, fig, out_dir / "supp_fig6_deprivation_domain_wave_mixing", "double", 6.2)
    plt.close("all")


def plot_mixing_domain_outcomes(
    style,
    domain_mixing_predictor_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Heatmap of domain-specific mixing predictors on count outcomes (ZTNB).

    Two panels (ZTNB cluster-size count ratio and ZTNB geographic-spread count
    ratio). The hurdle geographic-spread component is reported as a
    supplementary table; the domain-quintile log-ratios are 3-10x the
    magnitude of age/sex/age-sex log-ratios, so the colour scale is anchored
    to the non-domain-quintile predictors and domain-quintile cells saturate
    (the colour bar extends with triangles). Cell annotations carry the exact
    ratio.
    """
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.pyplot as plt

    panels = [
        ("cluster_size", "positive_zero_truncated_count", "Cluster size\nZTNB count ratio"),
        ("geographic_dispersion", "positive_zero_truncated_count", "Geographic spread\nZTNB count ratio"),
    ]

    data = domain_mixing_predictor_results.copy()
    data["predictor"] = data.apply(domain_mixing_predictor_key, axis=1)
    data = data[data["predictor"].isin(DOMAIN_MIXING_PREDICTOR_ORDER)].copy()
    if data.empty:
        return

    values = data[(data["outcome"].isin([panels[0][0], panels[1][0]]))
               & (data["component"].isin([panels[0][1], panels[1][1]]))]
    vmax = max(2, float(np.nanmax(np.abs(values["ratio"]))))
    vmin = min(0, float(np.nanmin(np.abs(values["ratio"]))))
    vcenter = 1
    norm = TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

    fig, axes = style.new_figure(
        width="double",
        height_in=4.4,
        nrows=1,
        ncols=2,
        font_scale=0.80,
    )
    image = None
    for idx, (outcome, component, title) in enumerate(panels):
        ax = axes[idx]
        sub = data[(data["outcome"] == outcome) & (data["component"] == component)]
        if sub.empty:
            matrix = pd.DataFrame(
                np.nan,
                index=DOMAIN_ORDER,
                columns=DOMAIN_MIXING_PREDICTOR_ORDER,
            )
        else:
            matrix = (
                sub.pivot_table(
                    index="domain",
                    columns="predictor",
                    values="ratio",
                    aggfunc="first",
                )
                .reindex(index=DOMAIN_ORDER, columns=DOMAIN_MIXING_PREDICTOR_ORDER)
            )
        image = ax.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", norm=norm, aspect="auto")
        ax.set_title(title, pad=6)
        ax.set_yticks(np.arange(len(DOMAIN_ORDER)))
        ax.set_yticklabels(
            [DOMAIN_LABELS[d] for d in DOMAIN_ORDER] if idx == 0 else []
        )
        ax.set_xticks(np.arange(len(DOMAIN_MIXING_PREDICTOR_ORDER)))
        ax.set_xticklabels(
            [DOMAIN_MIXING_PREDICTOR_LABELS[k] for k in DOMAIN_MIXING_PREDICTOR_ORDER],
            rotation=35,
            ha="right",
        )
        ax.tick_params(length=0)
        for y in np.arange(len(DOMAIN_ORDER) + 1) - 0.5:
            ax.axhline(y, color="white", linewidth=0.6)
        for x in np.arange(len(DOMAIN_MIXING_PREDICTOR_ORDER) + 1) - 0.5:
            ax.axvline(x, color="white", linewidth=0.6)
        if sub.empty:
            ax.text(
                (len(DOMAIN_MIXING_PREDICTOR_ORDER) - 1) / 2,
                (len(DOMAIN_ORDER) - 1) / 2,
                "Not estimable",
                ha="center",
                va="center",
                color="#666666",
                fontsize=8,
            )

    assert image is not None
    fig.subplots_adjust(left=0.13, right=0.86, top=0.86, bottom=0.22, wspace=0.10)
    cbar_ax = fig.add_axes([0.885, 0.24, 0.020, 0.58])
    cbar = fig.colorbar(image, cax=cbar_ax, extend="both")
    cbar.set_label("ZTNB count ratio per 1 SD higher excess mixing", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)
    style.add_panel_labels(axes, x=-0.10, y=1.12, size=9)
    save_all(style, fig, out_dir / "supp_fig7_mixing_domain_outcomes", "double", 4.4)
    plt.close("all")


def plot_deprivation_size_adjusted(
    style,
    count_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Size-adjusted positive geographic spread."""
    import matplotlib.pyplot as plt

    colours = term_colours(style)
    size_adjusted_terms = mixing_terms_for_results(count_results)
    outcome = SIZE_ADJUSTED_OUTCOMES[0]
    fig, ax = style.new_figure(
        width="onehalf",
        height_in=3.0,
        font_scale=0.85,
    )
    sub = count_results[
        (count_results["outcome"] == outcome)
        & (count_results["component"] == "positive_zero_truncated_count")
    ]
    draw_ratio_panel(
        ax,
        sub,
        size_adjusted_terms,
        colours,
        title=OUTCOME_LABELS[outcome],
        show_ylabels=True,
        xlim=(0.85, 3.2),
        xlabel="ZTNB count ratio per 1 SD higher covariate",
    )
    style.add_panel_labels([ax], x=-0.08, y=1.15, size=9)
    fig.subplots_adjust(left=0.32, right=0.98, top=0.86, bottom=0.19)
    save_all(style, fig, out_dir / "supp_fig8_deprivation_size_adjusted", "onehalf", 3.0)
    plt.close("all")


def _draw_loglinear_comparison(
    style,
    count_results: pd.DataFrame,
    loglinear_results: pd.DataFrame,
    terms: list[str],
    out_path: Path,
    *,
    xlabel: str,
) -> None:
    from matplotlib.ticker import NullFormatter, NullLocator
    import matplotlib.pyplot as plt

    outcomes = COUNT_OUTCOMES
    model_colours = {
        "Log-linear": "#666666",
        "Hurdle": "#4e79a7",
        "ZTNB positive": "#f28e2b",
    }
    markers = {"Log-linear": "o", "Hurdle": "s", "ZTNB positive": "^"}

    log_map = {"cluster_size": "cluster_size", "geographic_dispersion": "geographic_dispersion"}
    pieces = []
    log = loglinear_results[
        loglinear_results["model"].isin(log_map)
        & loglinear_results["term"].isin(terms)
    ].copy()
    log["outcome"] = log["model"].map(log_map).fillna(log["model"])
    log["model_type"] = "Log-linear"
    log = log.rename(
        columns={
            "geometric_mean_ratio": "ratio",
            "ci_low": "ratio_ci_low",
            "ci_high": "ratio_ci_high",
        }
    )
    pieces.append(log[["outcome", "term", "model_type", "ratio", "ratio_ci_low", "ratio_ci_high"]])

    hurdle = count_results[
        (count_results["outcome"].isin(outcomes))
        & (count_results["component"] == "hurdle_binary")
        & (count_results["term"].isin(terms))
    ].copy()
    hurdle["model_type"] = "Hurdle"
    pieces.append(hurdle[["outcome", "term", "model_type", "ratio", "ratio_ci_low", "ratio_ci_high"]])

    ztnb = count_results[
        (count_results["outcome"].isin(outcomes))
        & (count_results["component"] == "positive_zero_truncated_count")
        & (count_results["term"].isin(terms))
    ].copy()
    ztnb["model_type"] = "ZTNB positive"
    pieces.append(ztnb[["outcome", "term", "model_type", "ratio", "ratio_ci_low", "ratio_ci_high"]])
    comp = pd.concat(pieces, ignore_index=True)
    if comp.empty:
        return

    ci_min = float(comp["ratio_ci_low"].min())
    ci_max = float(comp["ratio_ci_high"].max())
    lower = min(0.6, max(0.25, np.floor(ci_min * 10.0) / 10.0))
    upper = max(1.5, min(30.0, np.ceil(ci_max * 10.0) / 10.0))
    tick_candidates = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0]
    ticks = [tick for tick in tick_candidates if lower <= tick <= upper]

    fig, axes = style.new_figure(
        width="double",
        height_in=4.0,
        nrows=1,
        ncols=len(outcomes),
        sharex=True,
        font_scale=0.85,
        layout="constrained",
    )
    offsets = {"Log-linear": -0.18, "Hurdle": 0.0, "ZTNB positive": 0.18}
    y_positions = np.arange(len(terms))[::-1]
    pos = dict(zip(terms, y_positions))
    for idx, outcome in enumerate(outcomes):
        ax = axes[idx]
        sub = comp[comp["outcome"] == outcome]
        for model_type in ["Log-linear", "Hurdle", "ZTNB positive"]:
            model_sub = sub[sub["model_type"] == model_type]
            for _, row in model_sub.iterrows():
                y = pos[row["term"]] + offsets[model_type]
                ax.plot(
                    [row["ratio_ci_low"], row["ratio_ci_high"]],
                    [y, y],
                    color=model_colours[model_type],
                    linewidth=0.9,
                    solid_capstyle="round",
                )
                ax.scatter(
                    row["ratio"],
                    y,
                    color=model_colours[model_type],
                    marker=markers[model_type],
                    s=17,
                    zorder=3,
                    label=model_type,
                )
        ax.axvline(1.0, color="#666666", linestyle="--", linewidth=0.7)
        ax.set_xscale("log")
        ax.set_xlim(lower, upper)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{tick:g}" for tick in ticks], ha="center")
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_title(OUTCOME_LABELS[outcome], pad=4)
        ax.set_ylim(-0.8, len(terms) - 0.2)
        ax.set_yticks(y_positions)
        ax.set_yticklabels([TERM_LABELS[t] for t in terms] if idx == 0 else [])
        ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    fig.supxlabel(xlabel, x=0.575, fontsize=8)
    handles, labels = axes[1].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    fig.legend(
        unique.values(),
        unique.keys(),
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.575, -0.065),
    )
    style.add_panel_labels(axes, x=-0.1, y=1.1, size=9)
    save_all(style, fig, out_path, "double", 4.0)
    plt.close("all")


def plot_deprivation_loglinear(
    style,
    count_results: pd.DataFrame,
    loglinear_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    _draw_loglinear_comparison(
        style,
        count_results,
        loglinear_results,
        PRIMARY_TERMS,
        out_dir / "supp_fig9_deprivation_loglinear",
        xlabel="Model-specific ratio per 1 SD higher covariate",
    )


def plot_mixing_loglinear(
    style,
    count_results: pd.DataFrame,
    loglinear_results: pd.DataFrame,
    out_dir: Path,
) -> None:
    _draw_loglinear_comparison(
        style,
        count_results,
        loglinear_results,
        MIXING_PREDICTOR_TERMS,
        out_dir / "supp_fig10_mixing_loglinear",
        xlabel="Model-specific ratio per 1 SD higher excess mixing",
    )


# ---------------------------------------------------------------------------
# Top-level run
# ---------------------------------------------------------------------------


def run(
    root: Path,
    tables_dir: Path | None = None,
    out_dir: Path | None = None,
    cache_dir: Path | None = None,
) -> None:
    setup_environment()
    style = load_style(root)

    main_dir = root / "part1"
    if tables_dir is None:
        tables_dir = main_dir / "tables"
    if out_dir is None:
        out_dir = main_dir / "manuscript" / "figures"
    if cache_dir is None:
        cache_dir = main_dir / "cache"

    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Required core tables
    # ------------------------------------------------------------------
    count_results = pd.read_csv(tables_dir / "hurdle_count_model_results.csv")
    mixing_results = pd.read_csv(tables_dir / "mixing_model_results.csv")
    cluster_table = pd.read_parquet(cache_dir / "cluster_table.parquet")

    # ------------------------------------------------------------------
    # Main figures
    # ------------------------------------------------------------------
    # fig1 — Line 1 overall (counts + mixing)
    plot_deprivation_overall(style, count_results, mixing_results, out_dir)

    # fig2 — Line 1 by wave
    wave_count_path = tables_dir / "wave_specific_hurdle_count_model_results.csv"
    if wave_count_path.exists():
        plot_deprivation_wave_specific(style, pd.read_csv(wave_count_path), out_dir)

    # fig3 — Line 2 overall (mixing predictors)
    mixing_predictor_path = tables_dir / "mixing_predictor_hurdle_count_model_results.csv"
    mixing_predictor_count_results = None
    if mixing_predictor_path.exists():
        mixing_predictor_count_results = pd.read_csv(mixing_predictor_path)
        plot_mixing_overall(style, mixing_predictor_count_results, out_dir)

    # fig4 — Line 2 by wave
    wave_mixing_predictor_path = (
        tables_dir / "wave_specific_mixing_predictor_hurdle_count_model_results.csv"
    )
    if wave_mixing_predictor_path.exists():
        wmp = pd.read_csv(wave_mixing_predictor_path)
        plot_mixing_wave_specific(style, wmp, out_dir)
        # Companion supplementary table: the hurdle geographic-spread component
        # dropped from the figure.
        write_hurdle_geographic_spread_table(
            wmp,
            out_dir / "supp_table_fig4_wave_mixing_hurdle_geographic_spread.csv",
            group_cols=["wave_group"],
        )

    # ------------------------------------------------------------------
    # Supplementary figures
    # ------------------------------------------------------------------
    # supp_fig1 — count outcome distributions
    plot_outcome_distributions(style, cluster_table, out_dir)

    # supp_fig2 — excess mixing distributions
    plot_mixing_distributions(style, cluster_table, out_dir)

    # supp_fig3 — observed-expected matrices
    obs_exp_path = tables_dir / "observed_expected_mixing_matrices.csv"
    if obs_exp_path.exists():
        plot_observed_expected_matrices(style, pd.read_csv(obs_exp_path), out_dir)

    # supp_fig4 — deprivation × domain on count outcomes
    domain_outcome_path = tables_dir / "simd_domain_hurdle_count_model_results.csv"
    if domain_outcome_path.exists():
        plot_deprivation_domain_outcomes(style, pd.read_csv(domain_outcome_path), out_dir)

    # supp_fig5 — deprivation × domain on mixing outcomes
    domain_mixing = None
    domain_demo = None
    domain_mixing_path = tables_dir / "simd_domain_quintile_mixing_model_results.csv"
    if domain_mixing_path.exists():
        domain_mixing = pd.read_csv(domain_mixing_path)
    domain_demo_path = tables_dir / "simd_domain_demographic_mixing_model_results.csv"
    if domain_demo_path.exists():
        domain_demo = pd.read_csv(domain_demo_path)
    if domain_mixing is not None and domain_demo is not None:
        plot_deprivation_domain_mixing(style, domain_mixing, domain_demo, out_dir)

    # supp_fig6 — wave × domain demographic mixing heatmap
    wave_domain_demo_path = (
        tables_dir / "wave_specific_domain_demographic_mixing_model_results.csv"
    )
    if wave_domain_demo_path.exists():
        plot_deprivation_domain_wave_mixing(style, pd.read_csv(wave_domain_demo_path), out_dir)

    # supp_fig7 — mixing × domain on count outcomes
    domain_mixing_predictor_path = (
        tables_dir / "simd_domain_mixing_predictor_hurdle_count_model_results.csv"
    )
    if domain_mixing_predictor_path.exists():
        dmp = pd.read_csv(domain_mixing_predictor_path)
        plot_mixing_domain_outcomes(style, dmp, out_dir)
        # Companion supplementary table: the hurdle geographic-spread component
        # dropped from the figure.
        write_hurdle_geographic_spread_table(
            dmp,
            out_dir / "supp_table_fig7_domain_mixing_hurdle_geographic_spread.csv",
            group_cols=["domain"],
        )

    # supp_fig8 — size-adjusted positive geographic spread
    plot_deprivation_size_adjusted(style, count_results, out_dir)

    # supp_fig9 — deprivation log-linear vs hurdle/ZTNB
    loglinear_path = tables_dir / "loglinear_count_model_results.csv"
    if loglinear_path.exists():
        plot_deprivation_loglinear(style, count_results, pd.read_csv(loglinear_path), out_dir)

    # supp_fig10 — mixing-predictor log-linear vs hurdle/ZTNB
    mixing_loglinear_path = tables_dir / "mixing_predictor_loglinear_count_model_results.csv"
    if mixing_loglinear_path.exists() and mixing_predictor_count_results is not None:
        plot_mixing_loglinear(
            style,
            mixing_predictor_count_results,
            pd.read_csv(mixing_loglinear_path),
            out_dir,
        )

    # Supplementary-files document: captions for every supp figure + the
    # supplementary tables for the hurdle geographic-spread results.
    write_supplementary_files_document(out_dir)

    print(f"Wrote manuscript figures to {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing the model result CSV tables to plot. "
            "Defaults to part1/tables (primary results). "
            "Pass the --tables-dir used for a sensitivity run to plot those results."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory to write figures into. Defaults to part1/manuscript/figures.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Directory containing main_cluster_table.parquet. Defaults to part1/cache.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    run(
        root,
        tables_dir=args.tables_dir.resolve() if args.tables_dir else None,
        out_dir=args.out_dir.resolve() if args.out_dir else None,
        cache_dir=args.cache_dir.resolve() if args.cache_dir else None,
    )


if __name__ == "__main__":
    main()
