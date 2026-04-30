import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

def main():
    df = pd.read_csv("regression_per_wave_results_adj_resolution.csv")

    # ── Shared config ─────────────────────────────────────────────────────────────
    WAVE_ORDER = [
        "WV1_B.1.177_C108360", "WV2_B.1.1.7_C574152", "WV3_AY.4_C983568",
        "WV4_BA.2_C479360",    "WV5_BA.2_C85080",      "WV6_BA.5.2_C23016",
    ]
    WAVE_FULL = ["WV1 – B.1.177", "WV2 – B.1.1.7", "WV3 – AY.4",
                 "WV4 – BA.2",    "WV5 – BA.2*",    "WV6 – BA.5.2"]
    WAVE_COLORS = ["#D4603A", "#C99A28", "#3E8B55", "#2E75B0", "#6B4DA0", "#B03070"]

    CAPTION = (
        "Adjusted for Leiden resolution, resolution², proportion sequenced, "
        "and their interaction. GEE with exchangeable correlation;\n"
        "reference categories: age 40–59, male, SIMD quintile 3."
    )

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.axisbelow": True,
        "axes.grid": True,
        "grid.color": "#E2E8F0",
        "grid.linewidth": 0.7,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })

    NULL_KW = dict(color="#94A3B8", linewidth=0.9, linestyle="--", zorder=1)

    def sig_marker(p):
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return ""

    def wave_row(df, term, wave):
        rows = df[(df["wave"] == wave) & (df["term"] == term)]
        return rows.iloc[0] if not rows.empty else None

    # Shared legend handles (used by all three figures)
    def make_legend_handles():
        sig_h   = mlines.Line2D([], [], marker="o", color="#475569",
                                markerfacecolor="#475569", markersize=7, linewidth=1.6,
                                label="Significant (p < 0.05)")
        nosig_h = mlines.Line2D([], [], marker="o", color="#475569",
                                markerfacecolor="white", markersize=7, linewidth=1.6,
                                markeredgewidth=1.5, label="Non-significant (p ≥ 0.05)")
        ref_h   = mlines.Line2D([], [], marker="D", color="#94A3B8",
                                markerfacecolor="white", markersize=6, linewidth=0,
                                markeredgewidth=1.2, label="Reference category")
        return [sig_h, nosig_h, ref_h]


    # ═══════════════════════════════════════════════════════════════════════════════
    # FIGURE 1 – SIMD deprivation gradient
    # ═══════════════════════════════════════════════════════════════════════════════
    SIMD_MAP = {
        1: "C(dz_simd_quintile, Treatment(3))[T.1]",
        2: "C(dz_simd_quintile, Treatment(3))[T.2]",
        3: None,
        4: "C(dz_simd_quintile, Treatment(3))[T.4]",
        5: "C(dz_simd_quintile, Treatment(3))[T.5]",
    }
    XLABELS_SIMD = ["Q1\n(Most deprived)", "Q2", "Q3\n(Ref.)", "Q4", "Q5\n(Least deprived)"]

    fig1, axes1 = plt.subplots(2, 3, figsize=(11, 6.5), sharey=True, sharex=True)
    axes1_flat = axes1.flatten()

    for wi, (wave, color) in enumerate(zip(WAVE_ORDER, WAVE_COLORS)):
        ax = axes1_flat[wi]

        coefs, lo, hi, pvals = [], [], [], []
        for q in [1, 2, 3, 4, 5]:
            term = SIMD_MAP[q]
            if term is None:
                coefs.append(0.0); lo.append(0.0); hi.append(0.0); pvals.append(1.0)
            else:
                row = wave_row(df, term, wave)
                if row is None:
                    coefs.append(np.nan); lo.append(np.nan); hi.append(np.nan); pvals.append(1.0)
                else:
                    coefs.append(row["coef"]); lo.append(row["ci_lo"])
                    hi.append(row["ci_hi"]); pvals.append(row["pvalue"])

        xs = np.arange(1, 6)
        coefs, lo, hi = np.array(coefs), np.array(lo), np.array(hi)

        ax.axhline(0, **NULL_KW)
        ax.fill_between(xs, lo, hi, alpha=0.17, color=color)
        ax.plot(xs, coefs, color=color, linewidth=2.0, zorder=4)

        for i, (q, p) in enumerate(zip([1,2,3,4,5], pvals)):
            if q == 3:
                ax.plot(xs[i], 0, marker="D", color="#94A3B8", markersize=7,
                        zorder=5, markeredgewidth=1.2, markeredgecolor="#64748B",
                        markerfacecolor="white")
            else:
                sig = p < 0.05
                ax.plot(xs[i], coefs[i], "o", color=color, markersize=7, zorder=5,
                        markerfacecolor=color if sig else "white",
                        markeredgewidth=1.5 if not sig else 0)

        ax.set_title(WAVE_FULL[wi], fontsize=10.5, fontweight="bold",
                     color="#1E293B", pad=5)
        ax.yaxis.set_major_locator(mticker.MultipleLocator(0.1))
        ax.tick_params(labelsize=9)

    # x-axis tick labels: bottom row only
    for ax in axes1[0]:
        plt.setp(ax.get_xticklabels(), visible=False)
    for ax in axes1[1]:
        ax.set_xticks(range(1, 6))
        ax.set_xticklabels(XLABELS_SIMD, fontsize=9)

    # y-axis labels: left column only
    for ax in [axes1[0,0], axes1[1,0]]:
        ax.set_ylabel("GEE coefficient (log-odds)\nvs SIMD Q3 reference", fontsize=10)

    fig1.legend(handles=make_legend_handles(), fontsize=9.5, framealpha=0.8,
                loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0),
                handletextpad=0.6, columnspacing=1.2)
    fig1.suptitle(
        "Socioeconomic deprivation (SIMD) and SARS-CoV-2 genomic clustering\nacross six epidemic waves",
        fontsize=13, fontweight="bold")
    fig1.text(0.5, -0.05, CAPTION, ha="center", fontsize=8.5,
              color="#64748B", style="italic")
    fig1.tight_layout(rect=[0, 0.07, 1, 1])
    fig1.savefig("figures/fig1_simd_v2.png", dpi=180, bbox_inches="tight")
    plt.close(fig1)
    print("fig1_simd.png done")


    # ═══════════════════════════════════════════════════════════════════════════════
    # FIGURE 2 – Age group effects (2 × 3 panel forest plot)
    # ═══════════════════════════════════════════════════════════════════════════════
    AGE_GROUPS = [
        ("C(age_group, Treatment('40\u201359'))[T.00-09]", "0–9"),
        ("C(age_group, Treatment('40\u201359'))[T.10-19]", "10–19"),
        ("C(age_group, Treatment('40\u201359'))[T.20-39]", "20–39"),
        ("__ref__",                                          "40–59\n(Ref.)"),
        ("C(age_group, Treatment('40\u201359'))[T.60-74]", "60–74"),
        ("C(age_group, Treatment('40\u201359'))[T.elderly]", "Elderly\n(75+)"),
    ]
    N_AGE = len(AGE_GROUPS)

    fig2, axes2 = plt.subplots(2, 3, figsize=(11, 6.5))
    axes2_flat = axes2.flatten()

    # compute shared xlim across all panels
    all_lo, all_hi = [], []
    for wave in WAVE_ORDER:
        for term, _ in AGE_GROUPS:
            if term == "__ref__": continue
            row = wave_row(df, term, wave)
            if row is not None:
                all_lo.append(row["ci_lo"]); all_hi.append(row["ci_hi"])
    XMIN = min(all_lo) - 0.04
    XMAX = max(all_hi) + 0.09   # small buffer for stars

    for wi, (wave, color) in enumerate(zip(WAVE_ORDER, WAVE_COLORS)):
        ax = axes2_flat[wi]
        ax.axvline(0, **NULL_KW)
        ax.set_xlim(XMIN, XMAX)

        for yi, (term, label) in enumerate(AGE_GROUPS):
            y = N_AGE - 1 - yi

            if term == "__ref__":
                ax.plot(0, y, marker="D", color="#94A3B8", markersize=7, zorder=5,
                        markeredgewidth=1.2, markeredgecolor="#64748B", markerfacecolor="white")
                continue

            row = wave_row(df, term, wave)
            if row is None: continue

            sig = row["pvalue"] < 0.05
            ax.plot([row["ci_lo"], row["ci_hi"]], [y, y],
                    color=color, linewidth=1.6, zorder=3, solid_capstyle="round")
            ax.plot(row["coef"], y, "o", color=color, markersize=8, zorder=5,
                    markerfacecolor=color if sig else "white",
                    markeredgewidth=1.8 if not sig else 0)

            star = sig_marker(row["pvalue"])
            if star:
                ax.text(row["ci_hi"] + 0.012, y, star, va="center",
                        fontsize=8.5, color=color, fontweight="bold")

        ax.set_yticks(range(N_AGE))
        ax.set_yticklabels([l for _, l in reversed(AGE_GROUPS)], fontsize=9)
        ax.set_title(WAVE_FULL[wi], fontsize=10.5, fontweight="bold",
                     color="#1E293B", pad=5)
        ax.set_ylim(-0.55, N_AGE - 0.45)
        ax.tick_params(labelsize=9)
        ax.grid(axis="y", linewidth=0)

        # x-label on bottom row only
        if wi >= 3:
            ax.set_xlabel("GEE coefficient (log-odds) vs 40–59", fontsize=9.5)
        else:
            plt.setp(ax.get_xticklabels(), visible=True)   # keep ticks, hide if crowded

    fig2.legend(handles=make_legend_handles(), fontsize=9.5, framealpha=0.8,
                loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0),
                handletextpad=0.6, columnspacing=1.2)
    fig2.suptitle(
        "Age group effects on SARS-CoV-2 genomic clustering\nacross six epidemic waves",
        fontsize=13, fontweight="bold")
    fig2.text(0.5, -0.05, CAPTION, ha="center", fontsize=8.5,
              color="#64748B", style="italic")
    fig2.tight_layout(rect=[0, 0.07, 1, 1])
    fig2.savefig("figures/fig2_age_v2.png", dpi=180, bbox_inches="tight")
    plt.close(fig2)
    print("fig2_age.png done")


    # ═══════════════════════════════════════════════════════════════════════════════
    # FIGURE 3 – Sex effect
    # ═══════════════════════════════════════════════════════════════════════════════
    fig3, ax3 = plt.subplots(figsize=(9, 4.2))

    # Compute tight data range first, then allocate annotation space via axes transform
    data_max = max(
        wave_row(df, "is_female", w)["ci_hi"]
        for w in WAVE_ORDER if wave_row(df, "is_female", w) is not None
    )
    data_min = min(
        wave_row(df, "is_female", w)["ci_lo"]
        for w in WAVE_ORDER if wave_row(df, "is_female", w) is not None
    )
    XPAD = (data_max - data_min) * 0.07
    XMIN3, XMAX3 = data_min - XPAD, data_max + XPAD

    ax3.axvline(0, **NULL_KW)
    ax3.axvspan(-0.015, 0.015, alpha=0.07, color="#94A3B8", zorder=0, label="_nolegend_")

    for wi, (wave, color) in enumerate(zip(WAVE_ORDER, WAVE_COLORS)):
        y = len(WAVE_ORDER) - 1 - wi
        row = wave_row(df, "is_female", wave)
        if row is None: continue

        sig = row["pvalue"] < 0.05
        # CI line + caps
        ax3.plot([row["ci_lo"], row["ci_hi"]], [y, y],
                 color=color, linewidth=2.2, zorder=3, solid_capstyle="butt")
        for x_cap in [row["ci_lo"], row["ci_hi"]]:
            ax3.plot([x_cap, x_cap], [y-0.13, y+0.13], color=color, lw=1.6)
        # Central point
        ax3.plot(row["coef"], y, "o", color=color, markersize=10, zorder=5,
                 markerfacecolor=color if sig else "white",
                 markeredgewidth=2.0 if not sig else 0)

        # Annotation using axes transform so it's always right-aligned outside the data range
        p_str = "p < 0.001" if row["pvalue"] < 0.001 else f"p = {row['pvalue']:.3f}"
        star  = sig_marker(row["pvalue"])
        ann = (
            f"{row['coef']:+.3f}  "
            f"[{row['ci_lo']:+.3f}, {row['ci_hi']:+.3f}]   "
            f"{p_str}"
            + (f"  {star}" if star else "")
        )
        ax3.annotate(
            ann, xy=(1.02, y),
            xycoords=("axes fraction", "data"),
            fontsize=8.5, va="center", color="#1E293B",
            fontweight="bold" if sig else "normal",
            annotation_clip=False,
        )

    ax3.set_yticks(range(len(WAVE_ORDER)))
    ax3.set_yticklabels(list(reversed(WAVE_FULL)), fontsize=10.5)
    ax3.set_xlabel("GEE coefficient (log-odds)  —  female vs male", fontsize=11)
    ax3.set_xlim(XMIN3, XMAX3)
    ax3.set_ylim(-0.65, len(WAVE_ORDER) - 0.35)
    ax3.tick_params(left=False, labelsize=10)
    ax3.spines["left"].set_visible(False)
    ax3.grid(axis="y", linewidth=0)
    ax3.grid(axis="x")

    fig3.legend(handles=make_legend_handles(), fontsize=9.5, framealpha=0.8,
                loc="lower center", ncol=3, bbox_to_anchor=(0.38, -0.01),
                handletextpad=0.6, columnspacing=1.2)
    ax3.set_title(
        "Sex effect on SARS-CoV-2 genomic clustering across six epidemic waves",
        fontsize=13, fontweight="bold", pad=10)
    fig3.text(0.38, -0.09, CAPTION, ha="center", fontsize=8.5,
              color="#64748B", style="italic")
    fig3.tight_layout(rect=[0, 0.06, 1, 1])
    fig3.savefig("figures/fig3_sex_v2.png", dpi=180, bbox_inches="tight")
    plt.close(fig3)
    print("fig3_sex.png done")

if __name__ == "__main__":
    main()