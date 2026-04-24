"""
plot_model_summary.py
---------------------
Read a Bambi / ArviZ posterior summary CSV and produce publication-quality
separate figures for a Negative Binomial regression model.

Usage
-----
    python plot_model_summary.py                      # uses embedded CSV string
    python plot_model_summary.py summary.csv          # read from file
    python plot_model_summary.py summary.csv ./figs/  # custom output directory
    python plot_model_summary.py summary.csv . png    # change extension

Output files (PDF by default; swap extension for PNG/SVG)
    fig1_main_effects.pdf
    fig2_wave_effects.pdf
    fig3_simd_quintile.pdf
    fig4_interaction_heatmap.pdf
    fig5_interaction_forest.pdf
"""

import sys
import io
import re
import textwrap
import pathlib

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── 0.  Global style ──────────────────────────────────────────────────────────

DARK_BG   = "#0f1117"
PANEL_BG  = "#161b27"
GRID_COL  = "#1e293b"
SPINE_COL = "#334155"
TEXT_COL  = "#cbd5e1"
MUTED_COL = "#64748b"

C_POS_SIG = "#22d3ee"    # cyan  – positive, HDI excludes 0
C_NEG_SIG = "#f87171"    # red   – negative, HDI excludes 0
C_POS_NS  = "#60a5fa"    # blue  – positive, HDI spans 0
C_NEG_NS  = "#fb923c"    # amber – negative, HDI spans 0

mpl.rcParams.update({
    "figure.facecolor":  DARK_BG,
    "axes.facecolor":    PANEL_BG,
    "axes.edgecolor":    SPINE_COL,
    "axes.labelcolor":   TEXT_COL,
    "axes.titlecolor":   TEXT_COL,
    "axes.grid":         True,
    "axes.axisbelow":    True,
    "grid.color":        GRID_COL,
    "grid.linewidth":    0.7,
    "xtick.color":       MUTED_COL,
    "ytick.color":       TEXT_COL,
    "xtick.labelsize":   8,
    "ytick.labelsize":   9.5,
    "font.family":       "sans-serif",
    "font.sans-serif":   ["DejaVu Sans"],
    "text.color":        TEXT_COL,
    "legend.facecolor":  PANEL_BG,
    "legend.edgecolor":  SPINE_COL,
    "legend.fontsize":   8.5,
    "savefig.facecolor": DARK_BG,
    "savefig.dpi":       180,
})

# ── 1.  Embedded CSV ──────────────────────────────────────────────────────────

EMBEDDED_CSV = textwrap.dedent("""
term,mean,sd,hdi_2.5%,hdi_97.5%,mcse_mean,mcse_sd,ess_bulk,ess_tail,r_hat
alpha,0.562,0.017,0.532,0.596,0.000,0.000,5954.0,2755.0,1.0
Intercept,0.732,0.352,-0.014,1.371,0.013,0.008,750.0,1244.0,1.0
median_age,-0.005,0.002,-0.009,-0.001,0.000,0.000,5232.0,2785.0,1.0
age_diversity,0.124,0.004,0.116,0.132,0.000,0.000,4704.0,3278.0,1.0
frac_female,-0.093,0.099,-0.284,0.099,0.001,0.002,5423.0,3079.0,1.0
wave[WV2_B.1.1.7_C71769],-0.827,0.377,-1.571,-0.082,0.013,0.008,802.0,1506.0,1.0
wave[WV3_AY.4_C122946],-0.905,0.359,-1.612,-0.203,0.013,0.008,733.0,1160.0,1.0
wave[WV4_BA.2_C59920],-0.968,0.361,-1.649,-0.224,0.013,0.007,791.0,1401.0,1.0
wave[WV5_BA.2_C10635],-0.586,0.386,-1.370,0.154,0.014,0.008,787.0,1199.0,1.0
wave[WV6_BA.5.2_C2877],-0.829,0.430,-1.698,-0.009,0.014,0.008,984.0,1734.0,1.0
"C(simd_quintile_mode, Treatment(3))[1]",-0.454,0.405,-1.240,0.357,0.014,0.008,805.0,1418.0,1.0
"C(simd_quintile_mode, Treatment(3))[2]",0.032,0.438,-0.784,0.944,0.015,0.008,908.0,1704.0,1.0
"C(simd_quintile_mode, Treatment(3))[4]",-0.694,0.494,-1.617,0.351,0.016,0.009,922.0,1550.0,1.0
"C(simd_quintile_mode, Treatment(3))[5]",-0.434,0.616,-1.607,0.783,0.020,0.012,1001.0,1644.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV2_B.1.1.7_C71769, 1]",0.368,0.455,-0.540,1.224,0.015,0.008,917.0,1746.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV2_B.1.1.7_C71769, 2]",0.514,0.488,-0.425,1.488,0.016,0.008,989.0,1698.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV2_B.1.1.7_C71769, 4]",0.597,0.553,-0.472,1.715,0.017,0.009,1009.0,1768.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV2_B.1.1.7_C71769, 5]",0.105,0.675,-1.265,1.374,0.020,0.012,1093.0,1930.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV3_AY.4_C122946, 1]",0.429,0.436,-0.392,1.306,0.015,0.008,887.0,1620.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV3_AY.4_C122946, 2]",-0.111,0.469,-1.083,0.752,0.015,0.008,960.0,1708.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV3_AY.4_C122946, 4]",1.145,0.526,0.096,2.181,0.017,0.009,988.0,1726.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV3_AY.4_C122946, 5]",0.934,0.647,-0.340,2.179,0.020,0.012,1010.0,1676.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV4_BA.2_C59920, 1]",1.427,0.433,0.557,2.250,0.014,0.008,899.0,1522.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV4_BA.2_C59920, 2]",0.205,0.470,-0.738,1.090,0.015,0.008,1044.0,1659.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV4_BA.2_C59920, 4]",1.038,0.527,0.057,2.148,0.016,0.009,1029.0,1739.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV4_BA.2_C59920, 5]",0.905,0.638,-0.358,2.102,0.020,0.012,1037.0,1862.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV5_BA.2_C10635, 1]",-0.259,0.481,-1.170,0.720,0.015,0.009,1030.0,1497.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV5_BA.2_C10635, 2]",-0.157,0.507,-1.177,0.787,0.016,0.009,1063.0,1609.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV5_BA.2_C10635, 4]",0.739,0.568,-0.367,1.847,0.018,0.010,1019.0,1917.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV5_BA.2_C10635, 5]",-0.040,0.713,-1.354,1.454,0.020,0.013,1238.0,1614.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV6_BA.5.2_C2877, 1]",-0.382,0.541,-1.353,0.752,0.016,0.009,1170.0,2098.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV6_BA.5.2_C2877, 2]",0.370,0.580,-0.767,1.500,0.016,0.008,1297.0,2079.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV6_BA.5.2_C2877, 4]",0.420,0.673,-0.924,1.730,0.018,0.010,1409.0,2303.0,1.0
"wave:C(simd_quintile_mode, Treatment(3))[WV6_BA.5.2_C2877, 5]",0.057,0.786,-1.530,1.583,0.022,0.013,1218.0,2128.0,1.0
simd_quintile_std,1.698,0.056,1.592,1.805,0.001,0.001,4529.0,3303.0,1.0
""").strip()

# ── 2.  Load & classify ───────────────────────────────────────────────────────

WAVE_LABELS = {
    "WV2_B.1.1.7_C71769": "WV2 · Alpha (B.1.1.7)",
    "WV3_AY.4_C122946":   "WV3 · Delta (AY.4)",
    "WV4_BA.2_C59920":    "WV4 · Omicron (BA.2)",
    "WV5_BA.2_C10635":    "WV5 · Omicron (BA.2)",
    "WV6_BA.5.2_C2877":   "WV6 · Omicron (BA.5.2)",
}
SIMD_LABELS = {
    "1": "Q1 (most deprived)",
    "2": "Q2",
    "4": "Q4",
    "5": "Q5 (least deprived)",
}


def load_summary(path=None) -> pd.DataFrame:
    src = open(path) if path else io.StringIO(EMBEDDED_CSV)
    df = pd.read_csv(src)
    df.columns = df.columns.str.strip()
    df["term"] = df["term"].str.strip()
    return df


def classify(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        t = r["term"]
        sig = not (r["hdi_2.5%"] <= 0 <= r["hdi_97.5%"])

        m = re.search(r"wave:C\(simd_quintile_mode.*?\)\[(.+?),\s*(\d)\]", t)
        if m:
            wave_key, q = m.group(1), m.group(2)
            label = (f"{WAVE_LABELS.get(wave_key, wave_key)} "
                     f"× {SIMD_LABELS.get(q, 'Q'+q)}")
            rows.append({**r, "group": "interaction", "label": label, "sig": sig})
            continue

        m = re.search(r"^wave\[(.+?)\]$", t)
        if m:
            label = WAVE_LABELS.get(m.group(1), m.group(1)) + " vs WV1"
            rows.append({**r, "group": "wave", "label": label, "sig": sig})
            continue

        m = re.search(r"C\(simd_quintile_mode.*?\)\[(\d)\]", t)
        if m:
            q = m.group(1)
            label = f"{SIMD_LABELS.get(q, 'Q'+q)} vs Q3"
            rows.append({**r, "group": "simd", "label": label, "sig": sig})
            continue

        clean = {
            "alpha":             "α  (dispersion)",
            "Intercept":         "Intercept",
            "median_age":        "Median age",
            "age_diversity":     "Age diversity",
            "frac_female":       "Fraction female",
            "simd_quintile_std": "SIMD quintile SD",
        }.get(t, t)
        rows.append({**r, "group": "main", "label": clean, "sig": sig})

    return pd.DataFrame(rows)


def point_color(row) -> str:
    if row["sig"]:
        return C_POS_SIG if row["mean"] > 0 else C_NEG_SIG
    return C_POS_NS if row["mean"] > 0 else C_NEG_NS


# ── 3.  Shared legend & subtitle ─────────────────────────────────────────────

LEGEND_HANDLES = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=C_POS_SIG,
           markersize=8, label="Positive — HDI excludes 0"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=C_NEG_SIG,
           markersize=8, label="Negative — HDI excludes 0"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor=PANEL_BG, markeredgecolor=C_POS_NS,
           markeredgewidth=1.8, markersize=8, label="Positive — HDI spans 0"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor=PANEL_BG, markeredgecolor=C_NEG_NS,
           markeredgewidth=1.8, markersize=8, label="Negative — HDI spans 0"),
    Line2D([0], [0], linestyle="--", color="#475569",
           linewidth=1.4, label="Null  (β = 0)"),
]

SUBTITLE = (
    "Negative Binomial model  ·  outcome: n_sequences_minus_one  "
    "·  reference: Wave 1 × SIMD Q3\n"
    "Points = posterior mean  ·  error bars = 95 % HDI"
)


# ── 4.  Forest-plot builder ───────────────────────────────────────────────────

def make_forest(
    sub: pd.DataFrame,
    title: str,
    *,
    sort: bool = True,
    show_ess: bool = True,
    fig_width: float = 9.5,
    row_height: float = 0.60,
    ax_left: float = 0.36,    # fraction: left margin for term labels
    ax_width: float = 0.50,   # fraction: width of the plotting area
) -> plt.Figure:
    """Return a standalone Figure with a single horizontal forest plot."""

    if sort:
        sub = sub.sort_values("mean").reset_index(drop=True)
    else:
        sub = sub.reset_index(drop=True)

    n = len(sub)
    fig_h = max(3.8, n * row_height + 2.4)

    fig = plt.figure(figsize=(fig_width, fig_h), facecolor=DARK_BG)
    # axes positioned manually so term labels sit in figure space to the left
    ax = fig.add_axes([ax_left, 0.13, ax_width, 0.70])
    ax.set_facecolor(PANEL_BG)

    xlo = sub["hdi_2.5%"].min()
    xhi = sub["hdi_97.5%"].max()
    span = xhi - xlo

    # alternating row bands
    for i in range(n):
        ax.axhspan(i - 0.5, i + 0.5,
                   color="#1a2030" if i % 2 == 0 else PANEL_BG, zorder=0)

    # zero reference line
    ax.axvline(0, color="#475569", linewidth=1.3,
               linestyle="--", zorder=1, alpha=0.85)

    for i, row in sub.iterrows():
        col = point_color(row)

        # CI bar
        ax.plot([row["hdi_2.5%"], row["hdi_97.5%"]], [i, i],
                color=col, linewidth=2.2, solid_capstyle="round",
                alpha=0.55, zorder=2)
        # end caps
        for xv in (row["hdi_2.5%"], row["hdi_97.5%"]):
            ax.plot([xv, xv], [i - 0.22, i + 0.22],
                    color=col, linewidth=2.0, zorder=3)
        # point estimate: solid if sig, hollow if not
        if row["sig"]:
            ax.plot(row["mean"], i, "o", markersize=7, color=col, zorder=4)
        else:
            ax.plot(row["mean"], i, "o", markersize=7,
                    markerfacecolor=PANEL_BG, markeredgecolor=col,
                    markeredgewidth=2.0, zorder=4)

        # mean annotation to the right of the plot area
        ax.annotate(
            f"{row['mean']:+.3f}",
            xy=(xhi + span * 0.05, i),
            xycoords="data",
            fontsize=8, color=col,
            fontweight="bold" if row["sig"] else "normal",
            va="center", ha="left",
            annotation_clip=False,
        )
        # ESS annotation further right
        if show_ess:
            ax.annotate(
                f"ESS {int(row['ess_bulk'])}",
                xy=(xhi + span * 0.28, i),
                xycoords="data",
                fontsize=7, color=MUTED_COL,
                va="center", ha="left",
                annotation_clip=False,
            )

    # column header labels at the top of annotation columns
    top = n - 0.05
    ax.annotate("Mean", xy=(xhi + span * 0.05, top), xycoords="data",
                fontsize=7.5, color=MUTED_COL, style="italic",
                va="bottom", ha="left", annotation_clip=False)
    if show_ess:
        ax.annotate("ESS (bulk)", xy=(xhi + span * 0.28, top), xycoords="data",
                    fontsize=7.5, color=MUTED_COL, style="italic",
                    va="bottom", ha="left", annotation_clip=False)

    # y-axis: term labels
    ax.set_yticks(range(n))
    ax.set_yticklabels(sub["label"].tolist(), fontsize=9.5)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlim(xlo - span * 0.08,
                xhi + span * (0.52 if show_ess else 0.28))
    ax.set_xlabel("Log-scale coefficient  (95 % HDI)", fontsize=9, labelpad=7)
    ax.tick_params(axis="y", length=0, pad=6)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(SPINE_COL)
    ax.spines["bottom"].set_color(SPINE_COL)
    ax.grid(axis="x", color=GRID_COL, linewidth=0.7, zorder=0)
    ax.grid(axis="y", visible=False)

    # figure-level title and subtitle
    fig.text(0.01, 0.97, title,
             fontsize=13, fontweight="bold", color=TEXT_COL, va="top")
    fig.text(0.01, 0.93, SUBTITLE,
             fontsize=7.5, color=MUTED_COL, va="top", linespacing=1.55)

    # legend at the bottom
    fig.legend(handles=LEGEND_HANDLES, loc="lower center",
               bbox_to_anchor=(0.5, -0.05), ncol=3,
               framealpha=0.55, fontsize=8)

    return fig


# ── 5.  Interaction heatmap ───────────────────────────────────────────────────

def make_heatmap(inter: pd.DataFrame) -> plt.Figure:
    """Return a standalone Figure with the Wave × SIMD posterior-mean heatmap."""

    WAVES = ["WV2", "WV3", "WV4", "WV5", "WV6"]
    QDEP  = ["Q1\n(most dep.)", "Q2", "Q4", "Q5\n(least dep.)"]
    Q_MAP = {"Q1": 0, "Q2": 1, "Q4": 2, "Q5": 3}
    W_MAP = {"WV2": 0, "WV3": 1, "WV4": 2, "WV5": 3, "WV6": 4}

    mat     = np.full((5, 4), np.nan)
    sig_mat = np.zeros((5, 4), dtype=bool)
    sd_mat  = np.full((5, 4), np.nan)

    for _, row in inter.iterrows():
        wm = re.search(r"(WV\d)", row["label"])
        qm = re.search(r"× (Q\d)", row["label"])
        if wm and qm:
            wi = W_MAP.get(wm.group(1))
            qi = Q_MAP.get(qm.group(1))
            if wi is not None and qi is not None:
                mat[wi, qi]     = row["mean"]
                sig_mat[wi, qi] = row["sig"]
                sd_mat[wi, qi]  = row["sd"]

    vmax = np.nanmax(np.abs(mat))
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "rwb",
        [(0.97, 0.35, 0.35), (0.12, 0.14, 0.20), (0.13, 0.83, 0.91)],
        N=256,
    )

    fig, ax = plt.subplots(figsize=(7.5, 5.5), facecolor=DARK_BG)
    ax.set_facecolor(PANEL_BG)
    fig.subplots_adjust(left=0.12, right=0.83, top=0.80, bottom=0.16)

    im = ax.imshow(mat, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(4));  ax.set_xticklabels(QDEP, fontsize=9.5)
    ax.set_yticks(range(5));  ax.set_yticklabels(WAVES, fontsize=9.5)
    ax.set_xlabel("SIMD quintile", fontsize=10, labelpad=8)
    ax.set_ylabel("Wave", fontsize=10, labelpad=8)
    ax.grid(False)

    for wi in range(5):
        for qi in range(4):
            v  = mat[wi, qi]
            sd = sd_mat[wi, qi]
            if np.isnan(v):
                continue
            bright = abs(v) > vmax * 0.35
            tc = "white" if bright else TEXT_COL
            fw = "bold" if sig_mat[wi, qi] else "normal"
            star = "*" if sig_mat[wi, qi] else ""
            ax.text(qi, wi - 0.13, f"{v:+.2f}{star}",
                    ha="center", va="center",
                    fontsize=11, color=tc, fontweight=fw)
            ax.text(qi, wi + 0.23, f"SD = {sd:.2f}",
                    ha="center", va="center",
                    fontsize=7.5, color=tc, alpha=0.75)

    for x in np.arange(-0.5, 4, 1):
        ax.axvline(x, color=SPINE_COL, linewidth=0.8)
    for y in np.arange(-0.5, 5, 1):
        ax.axhline(y, color=SPINE_COL, linewidth=0.8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cb = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.ax.tick_params(labelsize=8, colors=MUTED_COL)
    cb.outline.set_edgecolor(SPINE_COL)
    cb.set_label("Posterior mean  (log scale)", fontsize=8, color=MUTED_COL)

    fig.text(0.01, 0.97,
             "Wave × SIMD interaction  –  posterior mean",
             fontsize=13, fontweight="bold", color=TEXT_COL, va="top")
    fig.text(0.01, 0.92,
             SUBTITLE + "\n* = 95 % HDI excludes 0",
             fontsize=7.5, color=MUTED_COL, va="top", linespacing=1.55)

    return fig


# ── 6.  Save helper ───────────────────────────────────────────────────────────

def save(fig: plt.Figure, out_dir: pathlib.Path, stem: str, ext: str = "pdf"):
    out = out_dir / f"{stem}.{ext}"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved → {out}")


# ── 7.  Entry point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    csv_path = sys.argv[1] if len(sys.argv) > 1 else None
    out_dir  = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else pathlib.Path(".")
    ext      = sys.argv[3] if len(sys.argv) > 3 else "pdf"

    out_dir.mkdir(parents=True, exist_ok=True)

    df_raw = load_summary(csv_path)
    df     = classify(df_raw)

    main  = df[df["group"] == "main"]
    wave  = df[df["group"] == "wave"]
    simd  = df[df["group"] == "simd"]
    inter = df[df["group"] == "interaction"]

    print("Saving figures …")

    save(
        make_forest(main, "Main effects",
                    sort=True, show_ess=True,
                    fig_width=9.5, row_height=0.68),
        out_dir, "fig1_main_effects", ext,
    )
    save(
        make_forest(wave, "Wave effects  (vs WV1)",
                    sort=True, show_ess=False,
                    fig_width=9.0, row_height=0.62),
        out_dir, "fig2_wave_effects", ext,
    )
    save(
        make_forest(simd, "SIMD quintile effects  (vs Q3)",
                    sort=True, show_ess=False,
                    fig_width=9.0, row_height=0.62),
        out_dir, "fig3_simd_quintile", ext,
    )
    save(
        make_heatmap(inter),
        out_dir, "fig4_interaction_heatmap", ext,
    )
    save(
        make_forest(inter.sort_values("mean"), "Wave × SIMD interactions",
                    sort=False, show_ess=False,
                    fig_width=9.5, row_height=0.44),
        out_dir, "fig5_interaction_forest", ext,
    )

    print("Done.")