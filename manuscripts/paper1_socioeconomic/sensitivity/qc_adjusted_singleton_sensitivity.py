"""QC-adjusted sensitivity analysis for singleton logistic models (Paper 1, Suppl. S1).

Pre-specified in Statistical Analysis / Supplementary Methods S1.

Why this exists
---------------
The primary singleton analyses (Fig 4) restrict the analytic frame to Nextclade
``good`` genomes to prevent mediocre- and bad-QC assemblies from being pushed
artefactually into singleton clusters (inflated pairwise TN93 distances → weak
EpiLink edges → isolation in the Leiden graph). The question this sensitivity
analysis answers is: once QC tier is explicitly adjusted for in the regression,
does the headline deprivation gradient in singleton status change materially?

Trigger condition (pre-specified)
----------------------------------
The analysis is activated (i.e. QC-adjusted estimates are moved into Table 2
and discussed in Limitations) if EITHER:

  (i)  the Q1-versus-Q5 singleton OR changes by more than 10 % (absolute) in
       any variant epoch when mediocre and bad genomes are reintroduced, OR
  (ii) the direction of the quintile gradient reverses in any epoch (i.e. a
       quintile OR that was < 1 in the primary model crosses above 1 in the
       QC-adjusted model, or vice versa).

Otherwise the sensitivity findings are reported only in Supplementary Table S1
with the note that no material shift was detected.

Analysis strategy
-----------------
A new cluster-level regression frame is built from the full dataset (all
Nextclade QC tiers — good, mediocre, bad — retained, one row per
(window_id, cluster_id) at the primary Leiden resolution). Because cluster
membership itself changes when bad/mediocre sequences are included (some singletons
gain new members; new artefactual singletons may appear), ``is_singleton`` is
re-derived from this expanded frame rather than imported from the primary frame.

Two continuous cluster-level QC covariates are added:

  qc_frac_mediocre  fraction of sequences in the cluster with Nextclade QC
                    status ``mediocre`` (range 0–1; 0 in a clean cluster)
  qc_frac_bad       fraction of sequences in the cluster with Nextclade QC
                    status ``bad`` (range 0–1)
  reference: all sequences in cluster have QC ``good`` (both fractions = 0)

Singleton logistic GLMs are refitted epoch-by-epoch with the same specification
as the primary models:

  is_singleton ~ Q1..Q4 (vs Q5) + log-prop-seq + cr(time, df=3)
                 + qc_frac_mediocre + qc_frac_bad

Primary estimates are loaded from ``tables/fig4_singleton_ors.csv``; if that
file does not exist the primary models are re-fit from scratch.

Outputs (sensitivity/output/)
------------------------------
  qc_adjusted_singleton_ors.csv
      Tidy table with columns
      [term, estimate, std_error, conf_low, conf_high, z, p_value,
       quintile, epoch, model]
      where ``model`` ∈ {"primary", "qc_adjusted"}.

  qc_adjusted_trigger_report.txt
      Plain-text epoch-by-epoch comparison: Q1 OR, absolute % change,
      gradient-reversal flag, and the final TRIGGERED / NOT TRIGGERED verdict.

  fig_qc_adjusted_singleton_comparison.{pdf,png}
      Five-panel forest plot. Each panel = one VOC epoch. Within each panel,
      Q1–Q4 ORs (vs Q5) are shown for both models:
        • filled circles  = primary (good-QC only)
        • open diamonds   = QC-adjusted (all tiers + QC covariates)
      Quintiles are colour-coded with the shared SIMD palette.
      A vertical dashed line marks OR = 1. When the trigger condition fires,
      the panel title is annotated with "⚠ TRIGGERED".
"""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path
from functools import lru_cache

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from manuscripts.common import data, stats, style
from manuscripts.paper1_socioeconomic.models import simd_models


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TRIGGER_THRESHOLD_PCT: float = 10.0  # % change in Q1 OR that triggers action

# QC tier labels as they appear in the nextclade_qc column.
ALL_QC_TIERS: tuple[str, ...] = ("good", "mediocre", "bad")


# ---------------------------------------------------------------------------
# 1.  Build the all-QC cluster regression frame
# ---------------------------------------------------------------------------

@lru_cache(maxsize=4)
def build_allqc_cluster_regression_frame(
    resolution: float = data.PRIMARY_RESOLUTION,
    min_size: int = 1,
) -> pd.DataFrame:
    """One-row-per-cluster frame that includes all Nextclade QC tiers.

    Mirrors ``simd_models.build_cluster_regression_frame`` but:

    * Loads with ``qc=("good", "mediocre", "bad")`` so no sequences are dropped.
    * Explicitly includes ``nextclade_qc`` in the loaded columns.
    * Aggregates two extra cluster-level QC fraction covariates:
      ``qc_frac_mediocre`` and ``qc_frac_bad``.

    Because cluster membership expands when bad/mediocre sequences are added,
    ``is_singleton`` is re-derived here from the new ``n_sequences`` count and
    should NOT be imported from the primary frame.
    """

    cols = [
        "window_id", "cluster_id", "resolution", "sequence_id",
        "wn_mid_date", "wn_prop_sequenced", "who_voc", "pango_lineage",
        "datazone", "dz_simd_rank", "dz_simd_quintile",
        "dz_simd_income_rank", "dz_simd_employment_rank", "dz_simd_education_rank",
        "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank",
        "dz_simd_housing_rank",
        "nextclade_qc",  # ← needed for QC covariate aggregation
    ]

    df = data.load_analysis_columns(
        cols,
        resolution=resolution,
        qc=ALL_QC_TIERS,   # load all three tiers
    )

    # Pre-compute binary QC tier flags at sequence level for easy aggregation.
    df["_is_mediocre"] = (df["nextclade_qc"] == "mediocre").astype(float)
    df["_is_bad"]      = (df["nextclade_qc"] == "bad").astype(float)

    grp = df.groupby(["window_id", "cluster_id"], observed=True)

    def _mode(s: pd.Series):
        m = s.mode()
        return m.iloc[0] if len(m) > 0 else np.nan

    out = grp.agg(
        n_sequences        = ("sequence_id",         "nunique"),
        wn_mid_date        = ("wn_mid_date",          "first"),
        wn_prop_sequenced  = ("wn_prop_sequenced",    "first"),
        who_voc            = ("who_voc",               _mode),
        pango_lineage      = ("pango_lineage",         _mode),
        simd_quintile_mode = ("dz_simd_quintile",      _mode),
        simd_rank_mean     = ("dz_simd_rank",         "mean"),
        # SIMD domain means (same as primary frame)
        dz_simd_income_rank     = ("dz_simd_income_rank",     "mean"),
        dz_simd_employment_rank = ("dz_simd_employment_rank", "mean"),
        dz_simd_education_rank  = ("dz_simd_education_rank",  "mean"),
        dz_simd_health_rank     = ("dz_simd_health_rank",     "mean"),
        dz_simd_access_rank     = ("dz_simd_access_rank",     "mean"),
        dz_simd_crime_rank      = ("dz_simd_crime_rank",      "mean"),
        dz_simd_housing_rank    = ("dz_simd_housing_rank",    "mean"),
        # ── QC fraction covariates ────────────────────────────────────
        qc_frac_mediocre   = ("_is_mediocre",        "mean"),
        qc_frac_bad        = ("_is_bad",             "mean"),
    ).reset_index()

    out["is_singleton"] = (out["n_sequences"] == 1).astype(int)

    if min_size > 1:
        out = out[out["n_sequences"] >= min_size].reset_index(drop=True)

    return out


# ---------------------------------------------------------------------------
# 2.  QC-adjusted singleton epoch model
# ---------------------------------------------------------------------------


def singleton_epoch_model_qc_adjusted(
    frame: pd.DataFrame,
    *,
    min_prop_seq: float = 1e-3,
    include_time_spline: bool = True,
    time_spline_df: int = simd_models.TIME_SPLINE_DF_EPOCH,
):
    """Logistic GLM with explicit QC-tier covariates.

    Specification (per epoch)::

        is_singleton ~ Q1..Q4 (vs Q5)
                     + log-prop-seq (standardised)
                     + cr(wn_mid_date, df=3)       [if include_time_spline]
                     + qc_frac_mediocre
                     + qc_frac_bad

    Parameters
    ----------
    frame : pd.DataFrame
        Already restricted to a single epoch.  Must contain the columns
        produced by ``build_allqc_cluster_regression_frame`` — in particular
        ``qc_frac_mediocre`` and ``qc_frac_bad``.
    min_prop_seq : float
        Clusters with ``wn_prop_sequenced`` below this threshold are dropped
        (same guard as the primary model).
    include_time_spline : bool
        If True, add a natural-cubic-spline on ``wn_mid_date`` to absorb
        within-epoch temporal drift (default True, same as primary).
    time_spline_df : int
        Degrees of freedom for the within-epoch spline (default 3).

    Returns
    -------
    (fit, df_model) or (None, df) if the epoch sub-frame is too sparse.
    """
    required = [
        "simd_quintile_mode", "is_singleton",
        "wn_prop_sequenced", "wn_mid_date",
        "qc_frac_mediocre", "qc_frac_bad",
    ]
    df = frame.dropna(subset=required).copy()
    df = df[df["wn_prop_sequenced"] >= min_prop_seq]

    if len(df) < 50 or df["simd_quintile_mode"].nunique() < 2:
        return None, df

    # ── Design matrix ────────────────────────────────────────────────────
    q_dummies = stats.one_hot(
        df["simd_quintile_mode"].astype(int), reference=5, prefix="q"
    )

    log_prop = np.log(df["wn_prop_sequenced"].astype(float))
    log_prop_z = (log_prop - log_prop.mean()) / log_prop.std()

    X = q_dummies.copy()
    X["log_prop_seq_z"] = log_prop_z.values

    if include_time_spline:
        spline = simd_models._time_spline(df["wn_mid_date"], df=time_spline_df)
        X = pd.concat([X, spline], axis=1)

    # QC fraction covariates — already on [0, 1]; no further standardisation
    # needed, but we reset_index to align with X.
    X["qc_frac_mediocre"] = df["qc_frac_mediocre"].values
    X["qc_frac_bad"]      = df["qc_frac_bad"].values

    df_model = pd.concat(
        [df[["is_singleton"]].reset_index(drop=True), X.reset_index(drop=True)],
        axis=1,
    )

    fit = stats.logit_singleton(df_model, predictors=X.columns.tolist())
    return fit, df_model


# ---------------------------------------------------------------------------
# 3.  Loop over epochs → tidy OR table (QC-adjusted)
# ---------------------------------------------------------------------------


def build_qc_adjusted_singleton_epoch_table(
    frame: pd.DataFrame,
    *,
    epochs: list[str] | None = None,
) -> pd.DataFrame:
    """Fit the QC-adjusted logistic model in each epoch; return tidy OR table.

    Columns mirror ``simd_models.build_singleton_epoch_table`` so the two
    tables can be row-bound for comparison:
    [term, estimate, std_error, conf_low, conf_high, z, p_value,
     quintile, epoch]
    """
    frame = frame.copy()
    frame["epoch"] = data.assign_epoch(frame["wn_mid_date"])
    if epochs is None:
        epochs = [lbl for lbl, *_ in data.VOC_EPOCHS]

    out: list[pd.DataFrame] = []
    for epoch in epochs:
        sub = frame[frame["epoch"] == epoch]
        fit, _ = singleton_epoch_model_qc_adjusted(sub)
        if fit is None:
            print(f"  [warn] QC-adjusted model did not converge for epoch={epoch!r} "
                  f"(n={len(sub):,}); skipping.")
            continue
        tidy = stats.tidy_glm(fit)
        tidy = tidy[tidy["term"].str.startswith("q_")].copy()
        tidy["quintile"] = tidy["term"].str.replace("q_", "").astype(int)
        tidy["epoch"]    = epoch
        out.append(tidy)

    return (
        pd.concat(out, ignore_index=True)
        if out
        else pd.DataFrame(columns=[
            "term", "estimate", "std_error", "conf_low", "conf_high",
            "z", "p_value", "quintile", "epoch",
        ])
    )


# ---------------------------------------------------------------------------
# 4.  Trigger-condition evaluation
# ---------------------------------------------------------------------------


def evaluate_trigger(
    primary: pd.DataFrame,
    sensitivity: pd.DataFrame,
    *,
    quintile: int = 1,
    threshold_pct: float = TRIGGER_THRESHOLD_PCT,
) -> pd.DataFrame:
    """Return a per-epoch comparison table and an overall triggered flag.

    Parameters
    ----------
    primary : pd.DataFrame
        Tidy OR table from the primary (good-QC) singleton models.
    sensitivity : pd.DataFrame
        Tidy OR table from the QC-adjusted (all-tier) singleton models.
    quintile : int
        Which quintile's OR to evaluate for the trigger condition
        (default: Q1, the most deprived).
    threshold_pct : float
        Percentage-change threshold for condition (i).

    Returns
    -------
    pd.DataFrame with columns:
        epoch, or_primary, or_qc_adjusted, abs_pct_change,
        gradient_reversal, trigger_i, trigger_ii, triggered
    """
    p1 = primary[primary["quintile"] == quintile].set_index("epoch")["estimate"]
    s1 = sensitivity[sensitivity["quintile"] == quintile].set_index("epoch")["estimate"]

    epochs = sorted(set(p1.index) & set(s1.index),
                    key=lambda e: list(p1.index).index(e))

    rows = []
    for epoch in epochs:
        or_p = float(p1[epoch])
        or_s = float(s1[epoch])
        abs_pct = abs(or_s - or_p) / or_p * 100
        # gradient-reversal: primary OR < 1 (less likely singleton) but
        # sensitivity OR ≥ 1, or vice versa.
        reversal = (or_p < 1.0) != (or_s < 1.0)
        t_i  = abs_pct > threshold_pct
        t_ii = reversal
        rows.append(dict(
            epoch            = epoch,
            or_primary       = or_p,
            or_qc_adjusted   = or_s,
            abs_pct_change   = abs_pct,
            gradient_reversal= reversal,
            trigger_i        = t_i,
            trigger_ii       = t_ii,
            triggered        = t_i or t_ii,
        ))
    return pd.DataFrame(rows)


def format_trigger_report(
    report: pd.DataFrame,
    *,
    quintile: int = 1,
    threshold_pct: float = TRIGGER_THRESHOLD_PCT,
) -> str:
    """Format the trigger-condition evaluation as a readable plain-text report."""
    any_triggered = report["triggered"].any()
    verdict = "TRIGGERED" if any_triggered else "NOT TRIGGERED"

    lines = [
        "=" * 72,
        "QC-ADJUSTED SINGLETON SENSITIVITY ANALYSIS — TRIGGER REPORT",
        "=" * 72,
        "",
        f"Evaluated quintile : Q{quintile} (most deprived) vs Q5 (reference)",
        f"Threshold           : >{threshold_pct:.0f}% absolute change in OR triggers condition (i)",
        f"Overall verdict     : {verdict}",
        "",
        "Condition (i)  : |OR_qc_adjusted − OR_primary| / OR_primary × 100 > {:.0f}%".format(threshold_pct),
        "Condition (ii) : gradient reversal (OR crosses 1.0 between models)",
        "",
        "{:<20s}  {:>9s}  {:>12s}  {:>10s}  {:>9s}  {:>9s}  {:>10s}".format(
            "Epoch", "OR_primary", "OR_qc_adj", "Δ%", "Cond.(i)", "Cond.(ii)", "Triggered"
        ),
        "-" * 90,
    ]
    for _, row in report.iterrows():
        lines.append(
            "{:<20s}  {:>9.4f}  {:>12.4f}  {:>10.1f}  {:>9s}  {:>9s}  {:>10s}".format(
                row["epoch"],
                row["or_primary"],
                row["or_qc_adjusted"],
                row["abs_pct_change"],
                "YES ⚠" if row["trigger_i"]   else "no",
                "YES ⚠" if row["trigger_ii"]  else "no",
                "TRIGGERED" if row["triggered"] else "—",
            )
        )
    lines += [
        "-" * 90,
        "",
    ]

    if any_triggered:
        triggered_epochs = report.loc[report["triggered"], "epoch"].tolist()
        lines += [
            "ACTION: Material shift detected in: " + ", ".join(triggered_epochs),
            "",
            textwrap.fill(
                "Per pre-specification (Statistical Analysis § QC-adjusted sensitivity "
                "analysis), the QC-adjusted estimates should be placed alongside the "
                "primary estimates in Table 2 and the discrepancy should be discussed in "
                "the Strengths and Limitations subsection with reference to the mechanism "
                "(differential mediocre/bad-QC composition by SIMD quintile).",
                width=72,
            ),
        ]
    else:
        lines += [
            "No material shift detected.",
            "",
            textwrap.fill(
                "Per pre-specification, the QC-adjusted estimates are reported in "
                "Supplementary Table S1 with the note that no material shift was detected. "
                "The headline results in Table 2 / Fig 4 remain as primary.",
                width=72,
            ),
        ]

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 5.  Comparison figure
# ---------------------------------------------------------------------------


def make_comparison_figure(
    combined: pd.DataFrame,
    trigger_report: pd.DataFrame,
) -> plt.Figure:
    """Five-panel forest plot comparing primary vs QC-adjusted singleton ORs.

    Layout
    ------
    One panel per VOC epoch (left → right chronological order).
    Within each panel, Q1–Q4 ORs (vs Q5) are plotted on the x-axis and
    quintile on the y-axis.

    Markers
    -------
    • Filled circles      : primary (good-QC only)
    • Open diamonds (◇)   : QC-adjusted (all tiers + QC covariates)

    Quintile colours follow :data:`style.SIMD_QUINTILE_PALETTE`.
    A dashed vertical line marks OR = 1 (reference).
    Panel titles are annotated "⚠ TRIGGERED" in red when the trigger
    condition fires for that epoch.
    """
    epochs = [lbl for lbl, *_ in data.VOC_EPOCHS if lbl in set(combined["epoch"])]
    n_epochs = len(epochs)
    quintiles_plotted = [1, 2, 3, 4]

    fig, axes = style.new_figure(
        width="double",
        height_in=4.0,
        nrows=1,
        ncols=n_epochs,
        sharey=True,
        gridspec_kw={"wspace": 0.12},
    )
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])

    # y-offsets so primary and QC-adjusted markers don't overlap
    Y_OFFSET = {"primary": -0.15, "qc_adjusted": +0.15}
    MARKER    = {"primary": "o",   "qc_adjusted": "D"}
    MFILL     = {"primary": True,  "qc_adjusted": False}
    LABEL     = {"primary": "Good-QC only (primary)",
                 "qc_adjusted": "All QC + covariates (sensitivity)"}

    legend_handles: dict[str, object] = {}

    triggered_epochs = (
        set(trigger_report.loc[trigger_report["triggered"], "epoch"])
        if trigger_report is not None
        else set()
    )

    for ax, epoch in zip(axes, epochs):
        sub = combined[combined["epoch"] == epoch]
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
                color  = style.SIMD_QUINTILE_PALETTE[int(q)]
                mfc    = color if MFILL[model] else "white"
                y_pos  = float(q) + Y_OFFSET[model]
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

        title_color = "#cc0000" if epoch in triggered_epochs else "black"
        suffix = " ⚠" if epoch in triggered_epochs else ""
        ax.set_title(epoch + suffix, color=title_color, fontsize=8.5)

    axes[0].set_ylabel("SIMD quintile (vs Q5)")
    fig.supxlabel("OR vs. Q5 (least deprived)", y=-0.02)

    # Legend below the figure using the first two handles
    handles = [legend_handles.get("primary"), legend_handles.get("qc_adjusted")]
    labels  = [LABEL["primary"], LABEL["qc_adjusted"]]
    handles_clean = [(h, l) for h, l in zip(handles, labels) if h is not None]
    if handles_clean:
        hh, ll = zip(*handles_clean)
        fig.legend(
            hh, ll,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.18),
            ncol=2,
            fontsize=7.5,
            frameon=False,
        )

    return fig


# ---------------------------------------------------------------------------
# 6.  Main entry point
# ---------------------------------------------------------------------------


def main(
    out_dir: Path | None = None,
    *,
    force_refit_primary: bool = False,
    save_png: bool = True,
    save_pdf: bool = True,
) -> dict[str, Path]:
    """Run the full QC-adjusted sensitivity analysis and write outputs.

    Parameters
    ----------
    out_dir : Path, optional
        Directory where outputs are written.  Defaults to
        ``manuscripts/paper1_socioeconomic/sensitivity/output/``.
    force_refit_primary : bool
        If True, re-fit the primary (good-QC) singleton models even if
        ``tables/fig4_singleton_ors.csv`` already exists.
    save_png, save_pdf : bool
        Format flags for the comparison figure.

    Returns
    -------
    dict mapping output name → Path of each written file.
    """
    style.set_theme()
    paths_cfg = data.Paths.from_config()
    paper1_root = paths_cfg.root / "manuscripts" / "paper1_socioeconomic"
    out_dir = Path(out_dir) if out_dir else paper1_root / "sensitivity" / "output"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Primary estimates ─────────────────────────────────────────────
    primary_table_path = paper1_root / "tables" / "fig4_singleton_ors.csv"
    if primary_table_path.exists() and not force_refit_primary:
        print(f"Loading primary ORs from {primary_table_path}")
        primary_tab = pd.read_csv(primary_table_path)
    else:
        print("Fitting primary (good-QC) singleton models …")
        good_frame = simd_models.build_cluster_regression_frame()
        primary_tab = simd_models.build_singleton_epoch_table(good_frame)
        primary_tab.to_csv(primary_table_path, index=False)
        print(f"  Wrote {primary_table_path}")

    # ── 2. QC-adjusted estimates ─────────────────────────────────────────
    print("Building all-QC cluster regression frame …")
    allqc_frame = build_allqc_cluster_regression_frame()
    print(f"  Frame shape: {allqc_frame.shape[0]:,} cluster-window rows "
          f"({allqc_frame['qc_frac_bad'].gt(0).sum():,} clusters contain ≥1 bad-QC sequence; "
          f"{allqc_frame['qc_frac_mediocre'].gt(0).sum():,} contain ≥1 mediocre-QC sequence)")

    print("Fitting QC-adjusted singleton models …")
    sens_tab = build_qc_adjusted_singleton_epoch_table(allqc_frame)

    # ── 3. Trigger evaluation ────────────────────────────────────────────
    print("Evaluating trigger condition …")
    trigger_df = evaluate_trigger(primary_tab, sens_tab)
    report_txt = format_trigger_report(trigger_df)
    print(report_txt)

    # ── 4. Combined OR table ─────────────────────────────────────────────
    primary_tab_out = primary_tab.copy()
    primary_tab_out["model"] = "primary"
    sens_tab_out = sens_tab.copy()
    sens_tab_out["model"] = "qc_adjusted"
    combined = pd.concat([primary_tab_out, sens_tab_out], ignore_index=True)

    # ── 5. Write outputs ─────────────────────────────────────────────────
    saved: dict[str, Path] = {}

    csv_path = out_dir / "qc_adjusted_singleton_ors.csv"
    combined.to_csv(csv_path, index=False)
    saved["combined_ors_csv"] = csv_path
    print(f"Wrote {csv_path}")

    report_path = out_dir / "qc_adjusted_trigger_report.txt"
    report_path.write_text(report_txt)
    saved["trigger_report"] = report_path
    print(f"Wrote {report_path}")

    fig = make_comparison_figure(combined, trigger_df)
    fig_stem = out_dir / "fig_qc_adjusted_singleton_comparison"
    fig_paths = style.save_figure(
        fig, fig_stem,
        width="double",
        save_pdf=save_pdf,
        save_png=save_png,
    )
    saved.update(fig_paths)
    print(f"Wrote figure: " + ", ".join(str(p) for p in fig_paths.values()))

    return saved


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--output", type=Path, default=None,
        help="Output directory (default: sensitivity/output/ relative to paper1 root).",
    )
    ap.add_argument(
        "--force-refit-primary", action="store_true",
        help="Re-fit primary good-QC models even if fig4_singleton_ors.csv exists.",
    )
    ap.add_argument("--no-png",  action="store_true", help="Skip PNG output.")
    ap.add_argument("--no-pdf",  action="store_true", help="Skip PDF output.")
    args = ap.parse_args()

    written = main(
        out_dir=args.output,
        force_refit_primary=args.force_refit_primary,
        save_png=not args.no_png,
        save_pdf=not args.no_pdf,
    )
    print("\nAll outputs:")
    for k, v in written.items():
        print(f"  {k}: {v}")
