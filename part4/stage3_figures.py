"""
Part 4 – Stage 3: Figures for the Alpha emergence and F5 vs L2 counterfactual case study.

Fig 1: Alpha seeding chain — cluster size timeline with HB spread
Fig 2: S:N501Y explosive rise — observed frequency + logistic fits + key events
Fig 3: Counterfactual projections — four scenarios vs observed + hospital occupancy
Fig 4: Growth rate comparison — Alpha under F5 vs L2 vs B.1.177 decline
Fig 5: Lineage displacement — B.1.177 vs Alpha vs total cases, Oct 2020 – Apr 2021
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
from scipy.stats import linregress

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from utils.policy import POLICY_PERIODS_PD

OUT_T = Path(__file__).resolve().parent / "tables"
OUT_F = Path(__file__).resolve().parent / "figures"
OUT_F.mkdir(parents=True, exist_ok=True)

PERIOD_COLOURS = {
    "P3":"#aec6cf","T1":"#ffb347","F5":"#e8735a","L2":"#c23b22",
    "SL":"#e8d87b","L3":"#b4d8b4","L21":"#77dd77","L0":"#a8d8ea",
    "NN":"#cfcfc4","OM":"#b19cd9","FE":"#dfd3c3","PR":"#f0f0f0",
}
PERIOD_LABELS = dict(zip(POLICY_PERIODS_PD["period_code"], POLICY_PERIODS_PD["period_label"]))

POLICY_VLINES = [
    ("2020-10-02", "#ffb347", "T1", "Pre-tier"),
    ("2020-11-02", "#e8735a", "F5", "Five-tier"),
    ("2021-01-05", "#c23b22", "L2", "Lockdown 2"),
    ("2021-04-02", "#77dd77", "SL", "Stay local"),
]

HB_COLOURS = {
    "Greater Glasgow and Clyde": "#e74c3c",
    "Lothian":                   "#2980b9",
    "Lanarkshire":               "#27ae60",
    "Grampian":                  "#f39c12",
    "Tayside":                   "#8e44ad",
    "Ayrshire and Arran":        "#16a085",
    "Fife":                      "#d35400",
    "Highland":                  "#2c3e50",
    "Forth Valley":              "#95a5a6",
    "Borders":                   "#1abc9c",
    "Dumfries and Galloway":     "#e67e22",
    "Western Isles":             "#34495e",
    "Orkney":                    "#7f8c8d",
    "Shetland":                  "#bdc3c7",
}

def add_period_shading(ax, wpm, xmin="2020-10-01", xmax="2021-04-30", alpha=0.12):
    wpm_sub = wpm[(wpm["wn_mid_date"] >= pd.Timestamp(xmin)) &
                  (wpm["wn_mid_date"] <= pd.Timestamp(xmax))]
    for p, grp in wpm_sub.groupby("policy_period"):
        ax.axvspan(grp["wn_mid_date"].min(), grp["wn_mid_date"].max(),
                   color=PERIOD_COLOURS.get(str(p),"#ddd"), alpha=alpha, zorder=0)

def add_policy_vlines(ax, include=None):
    for date_str, col, code, label in POLICY_VLINES:
        if include and code not in include:
            continue
        ax.axvline(pd.Timestamp(date_str), color=col, lw=1.6, ls="--", alpha=0.85, zorder=4)

def fmt_date_axis(ax):
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right", fontsize=9)

# ── Load tables ───────────────────────────────────────────────────────────
print("Loading tables …")
wpm         = pd.read_csv(OUT_T / "part4_window_period_map.csv", parse_dates=["wn_mid_date"])
chain       = pd.read_csv(OUT_T / "part4_alpha_cluster_chain.csv", parse_dates=["first_seq_date","last_seq_date","wn_mid_date"])
overlaps    = pd.read_csv(OUT_T / "part4_alpha_chain_overlaps.csv")
traj        = pd.read_csv(OUT_T / "part4_mutation_trajectories.csv", parse_dates=["wn_mid_date"])
proj        = pd.read_csv(OUT_T / "part4_counterfactual_projections.csv", parse_dates=["date"])
hosp        = pd.read_csv(OUT_T / "part4_scotland_hospital.csv", parse_dates=["date"])
lin_comp    = pd.read_csv(OUT_T / "part4_lineage_composition.csv", parse_dates=["wn_mid_date"])
alpha_wkly  = pd.read_csv(OUT_T / "part4_alpha_clusters_weekly.csv", parse_dates=["wn_mid_date"])
b1177_wkly  = pd.read_csv(OUT_T / "part4_b1177_clusters_weekly.csv", parse_dates=["wn_mid_date"])
growth_p    = pd.read_csv(OUT_T / "part4_growth_params.csv", parse_dates=["t0"])
cf_proj     = pd.read_csv(OUT_T / "part4_counterfactual_projections.csv", parse_dates=["date"])

n501 = traj[(traj["mutation"]=="S:N501Y") & (traj["mut_type"]=="amino_acid")].sort_values("window_idx")
s222 = traj[(traj["mutation"]=="S:A222V") & (traj["mut_type"]=="amino_acid")].sort_values("window_idx")

# Growth params
gp = growth_p.set_index("label")
r_f5  = float(gp.loc["Alpha_F5","slope"])
i_f5  = float(gp.loc["Alpha_F5","intercept"])
t0_f5 = pd.Timestamp(gp.loc["Alpha_F5","t0"])
r_l2  = float(gp.loc["Alpha_L2","slope"])
i_l2  = float(gp.loc["Alpha_L2","intercept"])
t0_l2 = pd.Timestamp(gp.loc["Alpha_L2","t0"])
r_b1177 = float(gp.loc["B1177_L2_decline","slope"])

# 95 % CI on slopes (from positive-test weighted binomial GLM)
ci_lo_f5 = float(gp.loc["Alpha_F5","ci_lo"])
ci_hi_f5 = float(gp.loc["Alpha_F5","ci_hi"])
ci_lo_l2 = float(gp.loc["Alpha_L2","ci_lo"])
ci_hi_l2 = float(gp.loc["Alpha_L2","ci_hi"])

# Doubling-time CIs
d_f5     = float(gp.loc["Alpha_F5","doubling_days"])
d_lo_f5  = float(gp.loc["Alpha_F5","doubling_ci_lo"])
d_hi_f5  = float(gp.loc["Alpha_F5","doubling_ci_hi"])
d_l2     = float(gp.loc["Alpha_L2","doubling_days"])
d_lo_l2  = float(gp.loc["Alpha_L2","doubling_ci_lo"])
d_hi_l2  = float(gp.loc["Alpha_L2","doubling_ci_hi"])

pct_reduction = (1 - r_l2 / r_f5) * 100   # how many % slower L2 is than F5


# ════════════════════════════════════════════════════════════════════════════
# FIG 1: Alpha seeding chain
# ════════════════════════════════════════════════════════════════════════════
print("Generating fig1: Alpha seeding chain …")

early = chain[chain["first_seq_date"] <= "2021-01-10"].copy()

fig, axes = plt.subplots(2, 1, figsize=(15, 10))
fig.suptitle("Alpha (B.1.1.7) Emergence in Scotland: The Seeding Chain\n"
             "November 2020 – January 2021", fontsize=13, fontweight="bold")

# Panel A: cluster size bubbles by first_seq_date, coloured by primary HB
ax = axes[0]
add_period_shading(ax, wpm, "2020-10-20", "2021-01-15")
add_policy_vlines(ax, include=["F5","L2"])

for _, row in early.iterrows():
    col = HB_COLOURS.get(row["primary_hb"], "#aaaaaa")
    ax.scatter(row["first_seq_date"], row["size"],
               s=max(row["size"]*8, 30), color=col, alpha=0.75,
               edgecolors="white", lw=0.5, zorder=3)

# Annotate the index case and key clusters
annots = {
    "W016|B.1.1.7|R0.3|S0": ("Index case\nGlasgow, essential worker", "left"),
    "W017|B.1.1.7|R0.3|C1": ("First cluster\n(5 seqs, 2 HBs)", "right"),
    "W021|B.1.1.7|R0.3|C1": ("GGC chain peak\n(33 seqs)", "right"),
    "W022|B.1.1.7|R0.3|C6": ("Dec 8 explosion\n(67 seqs, 8 HBs)", "right"),
    "W025|B.1.1.7|R0.3|C25": ("National superspreading\n(135 seqs, 10 HBs)", "left"),
}
for cid, (label, side) in annots.items():
    row = early[early["cluster_id"]==cid]
    if row.empty:
        continue
    row = row.iloc[0]
    xoff = pd.Timedelta(days=3 if side=="right" else -3)
    ha = "left" if side=="right" else "right"
    ax.annotate(label, xy=(row["first_seq_date"], row["size"]),
                xytext=(row["first_seq_date"]+xoff, row["size"]+5),
                fontsize=7.5, ha=ha, va="bottom",
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.8),
                zorder=5)

# HB legend
hb_shown = early["primary_hb"].unique()
handles = [mpatches.Patch(color=HB_COLOURS.get(hb,"#aaa"), label=hb)
           for hb in sorted(hb_shown) if hb in HB_COLOURS]
ax.legend(handles=handles, fontsize=7.5, loc="upper left", ncol=2, title="Primary health board")
ax.set_ylabel("Cluster size (sequences)", fontsize=10)
ax.set_title("A  Alpha cluster sizes by first sequence date", fontsize=11, loc="left")
ax.set_yscale("symlog", linthresh=10)
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax.grid(axis="y", lw=0.5, alpha=0.3)
fmt_date_axis(ax)

# Panel B: health board spread per week (unique HBs carrying Alpha clusters)
ax = axes[1]
add_period_shading(ax, wpm, "2020-10-20", "2021-01-15")
add_policy_vlines(ax, include=["F5","L2"])

hb_spread = (alpha_wkly[alpha_wkly["wn_mid_date"] <= "2021-01-15"]
             .groupby("wn_mid_date")["n_hb"].max().reset_index())

# Stacked bar: n clusters per HB per week
import pyarrow.parquet as pq
clust_data = pq.read_table(
    str(REPO / "data/processed/scotland_clustering_analysis_dataset.parquet"),
    columns=["cluster_id","wn_mid_date","pango_lineage","dz_health_board","resolution","nextclade_qc"]
).to_pandas()
clust_data = clust_data[(clust_data["resolution"]==0.3) & (clust_data["nextclade_qc"]=="good")]
clust_data["wn_mid_date"] = pd.to_datetime(clust_data["wn_mid_date"])
alpha_hb = clust_data[
    (clust_data["pango_lineage"].str.startswith("B.1.1.7", na=False)) &
    (clust_data["wn_mid_date"] >= "2020-10-20") &
    (clust_data["wn_mid_date"] <= "2021-01-15")
]
# Per week: unique HBs with at least one Alpha cluster
hb_per_week = (alpha_hb.groupby(["wn_mid_date","dz_health_board"])["cluster_id"]
               .nunique().reset_index(name="n_clusters"))

# Plot cumulative unique HBs over time
wkly_hbs = hb_per_week.groupby("wn_mid_date")["dz_health_board"].nunique().reset_index(name="n_hbs_active")
ax.bar(wkly_hbs["wn_mid_date"], wkly_hbs["n_hbs_active"],
       width=6, color="#e74c3c", alpha=0.8, label="Health boards with Alpha clusters")
ax.set_ylabel("Health boards with active Alpha clusters", fontsize=10)
ax.set_title("B  Geographic spread of Alpha: health boards per week", fontsize=11, loc="left")
ax.set_yticks(range(0, 15, 2))
ax.legend(fontsize=9)
ax.grid(axis="y", lw=0.5, alpha=0.3)
fmt_date_axis(ax)

for date_str, col, code, label in POLICY_VLINES:
    if code in ["F5","L2"]:
        for a in axes:
            a.text(pd.Timestamp(date_str), a.get_ylim()[1]*0.97, f" {code}",
                   color=col, fontsize=8, fontweight="bold", va="top")

plt.tight_layout()
fig.savefig(OUT_F / "fig1_alpha_seeding_chain.png", dpi=180, bbox_inches="tight")
plt.close()
print("  Saved fig1")
del clust_data, alpha_hb


# ════════════════════════════════════════════════════════════════════════════
# FIG 2: S:N501Y explosive rise with logistic fits
# ════════════════════════════════════════════════════════════════════════════
print("Generating fig2: S:N501Y trajectory …")

fig, ax = plt.subplots(figsize=(15, 6))
add_period_shading(ax, wpm, "2020-10-01", "2021-04-15")
add_policy_vlines(ax, include=["T1","F5","L2","SL"])

# Observed data
n501_plot = n501[(n501["wn_mid_date"] >= "2020-10-01") & (n501["wn_mid_date"] <= "2021-04-15")]
ax.scatter(n501_plot["wn_mid_date"], n501_plot["frequency"],
           color="#c0392b", s=40, zorder=5, label="Observed S:N501Y frequency", marker="o")
ax.plot(n501_plot["wn_mid_date"], n501_plot["frequency"],
        color="#c0392b", lw=1.5, zorder=4, alpha=0.7)

# Logistic fit curves
fit_dates = pd.date_range("2020-11-03", "2021-03-10", freq="3D")
anchor_logit = np.log(0.0066 / 0.9934)
anchor_date  = pd.Timestamp("2020-11-03")

# F5 fit + 95% CI band
f5_curve_x = pd.date_range("2020-11-03","2021-01-10",freq="3D")
f5_curve_y    = np.array([1/(1+np.exp(-(r_f5    *(d-anchor_date).days + anchor_logit))) for d in f5_curve_x])
f5_curve_y_lo = np.array([1/(1+np.exp(-(ci_lo_f5*(d-anchor_date).days + anchor_logit))) for d in f5_curve_x])
f5_curve_y_hi = np.array([1/(1+np.exp(-(ci_hi_f5*(d-anchor_date).days + anchor_logit))) for d in f5_curve_x])
ax.fill_between(f5_curve_x, f5_curve_y_lo, f5_curve_y_hi, color="#e8735a", alpha=0.18, zorder=2)
ax.plot(f5_curve_x, f5_curve_y, color="#e8735a", lw=2, ls="--", zorder=3,
        label=f"Binomial GLM under F5 (doubling = {d_f5:.1f}d, 95% CI {d_lo_f5:.1f}–{d_hi_f5:.1f}d)")

# L2 fit + 95% CI band
t0_l2_anchor_logit = np.log(0.604 / 0.396)
l2_curve_x = pd.date_range("2021-01-05","2021-03-10",freq="3D")
l2_t0 = pd.Timestamp("2021-01-05")
l2_curve_y    = np.array([1/(1+np.exp(-(r_l2    *(d-l2_t0).days + t0_l2_anchor_logit))) for d in l2_curve_x])
l2_curve_y_lo = np.array([1/(1+np.exp(-(ci_lo_l2*(d-l2_t0).days + t0_l2_anchor_logit))) for d in l2_curve_x])
l2_curve_y_hi = np.array([1/(1+np.exp(-(ci_hi_l2*(d-l2_t0).days + t0_l2_anchor_logit))) for d in l2_curve_x])
ax.fill_between(l2_curve_x, l2_curve_y_lo, l2_curve_y_hi, color="#c23b22", alpha=0.18, zorder=2)
ax.plot(l2_curve_x, l2_curve_y, color="#c23b22", lw=2, ls="--", zorder=3,
        label=f"Binomial GLM under L2 (doubling = {d_l2:.1f}d, 95% CI {d_lo_l2:.1f}–{d_hi_l2:.1f}d)")

# Dec 8 explosion annotation
ax.annotate("Dec 8\n5.5× jump in 1 week\n(3.2% → 17.7%)",
            xy=(pd.Timestamp("2020-12-08"), 0.177),
            xytext=(pd.Timestamp("2020-12-15"), 0.32),
            fontsize=8.5, ha="left",
            arrowprops=dict(arrowstyle="->", color="black", lw=1.2),
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9))

# Index case annotation
ax.annotate("Index case: Glasgow City\n4 Nov · essential worker",
            xy=(pd.Timestamp("2020-11-04"), 0.007),
            xytext=(pd.Timestamp("2020-10-20"), 0.08),
            fontsize=8, ha="left",
            arrowprops=dict(arrowstyle="->", color="#e74c3c", lw=1.2),
            color="#e74c3c",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#e74c3c", alpha=0.9))

# L2 imposition marker
ax.axvline(pd.Timestamp("2021-01-05"), color="#c23b22", lw=2, ls="-", alpha=0.9, zorder=3)
ax.text(pd.Timestamp("2021-01-05"), 0.65, " L2 imposed\n 60% Alpha",
        color="#c23b22", fontsize=9, va="bottom", fontweight="bold")

ax.set_ylabel("Proportion of sequences carrying S:N501Y", fontsize=11)
ax.set_xlabel("Window mid-date", fontsize=11)
ax.set_title("S:N501Y (Alpha) Frequency in Scotland: Explosive Rise During F5\n"
             "From index case to majority within 9 weeks", fontsize=12, fontweight="bold")
ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax.set_ylim(-0.02, 1.05)
ax.legend(fontsize=9, loc="upper left")
ax.grid(axis="y", lw=0.5, alpha=0.3)
fmt_date_axis(ax)

# Period label strip at bottom
for p in ["T1","F5","L2","SL"]:
    row_p = POLICY_PERIODS_PD[POLICY_PERIODS_PD["period_code"]==p]
    if row_p.empty: continue
    mid = row_p["start_date"].iloc[0] + (row_p["end_date"].iloc[0]-row_p["start_date"].iloc[0])/2
    ax.text(mid, -0.018, f"{p}", ha="center", fontsize=8, fontweight="bold",
            color=PERIOD_COLOURS.get(p,"#444"), transform=ax.get_xaxis_transform())

plt.tight_layout()
fig.savefig(OUT_F / "fig2_n501y_explosive_rise.png", dpi=180, bbox_inches="tight")
plt.close()
print("  Saved fig2")


# ════════════════════════════════════════════════════════════════════════════
# FIG 3: Counterfactual projections + hospital occupancy
# ════════════════════════════════════════════════════════════════════════════
print("Generating fig3: Counterfactual …")

SCEN_STYLES = {
    "Actual (F5 → L2 on 5 Jan)": ("#c0392b", "-",  2.5, "Actual timeline"),
    "L2 from 2 Nov (immediate)": ("#2980b9", "--", 2.0, "L2 from 2 Nov"),
    "L2 from 2 Dec":             ("#8e44ad", "--", 1.8, "L2 from 2 Dec"),
    "L2 from 8 Dec (explosion)": ("#27ae60", "--", 1.8, "L2 from 8 Dec"),
}

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 11), sharex=True)
fig.suptitle("Counterfactual Analysis: Would Earlier L2 Have Changed Alpha's Trajectory?\n"
             "Logistic model under five-tier (F5) vs lockdown (L2) growth rates",
             fontsize=12, fontweight="bold")

add_period_shading(ax1, wpm, "2020-10-20", "2021-03-15")
add_policy_vlines(ax1, include=["F5","L2"])

# Observed S:N501Y
n501_cf = n501[(n501["wn_mid_date"] >= "2020-10-20") & (n501["wn_mid_date"] <= "2021-03-15")]
ax1.scatter(n501_cf["wn_mid_date"], n501_cf["frequency"],
            color="#c0392b", s=45, zorder=6, label="Observed S:N501Y", marker="o")

for scen, (col, ls, lw, lab) in SCEN_STYLES.items():
    sub = cf_proj[cf_proj["scenario"]==scen]
    ax1.plot(sub["date"], sub["frequency"], color=col, ls=ls, lw=lw,
             label=lab, zorder=4, alpha=0.9)

# Mark 50% threshold
ax1.axhline(0.50, color="gray", lw=1, ls=":", alpha=0.7)
ax1.text(pd.Timestamp("2020-10-22"), 0.51, "50% dominance", color="gray", fontsize=8)

# 5 Jan annotation
ax1.axvline(pd.Timestamp("2021-01-05"), color="#c23b22", lw=1.5, ls="-", alpha=0.7)

# Annotate key dates for each scenario
date_50 = {}
for scen in SCEN_STYLES:
    sub = cf_proj[(cf_proj["scenario"]==scen) & (cf_proj["frequency"]>=0.50)]
    if len(sub):
        date_50[scen] = sub["date"].min()
col_scen = list(SCEN_STYLES.keys())
for i, (scen, d50) in enumerate(date_50.items()):
    col = SCEN_STYLES[scen][0]
    ax1.axvline(d50, color=col, lw=0.8, ls=":", alpha=0.5, zorder=2)

ax1.set_ylabel("S:N501Y (Alpha) frequency", fontsize=11)
ax1.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax1.set_ylim(-0.03, 1.05)
ax1.legend(fontsize=9, loc="upper left")
ax1.grid(axis="y", lw=0.5, alpha=0.3)
ax1.set_title("A  Alpha frequency trajectory under each scenario", fontsize=11, loc="left")

# Panel B: hospital occupancy
add_period_shading(ax2, wpm, "2020-10-20", "2021-03-15")
add_policy_vlines(ax2, include=["F5","L2"])

hosp_plot = hosp[(hosp["date"] >= "2020-10-20") & (hosp["date"] <= "2021-03-15")]
hosp_7d = hosp_plot.set_index("date")["hb_hospital_occupancy"].rolling(7, center=True).mean().reset_index()
ax2.fill_between(hosp_plot["date"], hosp_plot["hb_hospital_occupancy"],
                 alpha=0.2, color="#c0392b", label="Daily hospital occupancy")
ax2.plot(hosp_7d["date"], hosp_7d["hb_hospital_occupancy"],
         color="#c0392b", lw=2.5, label="7-day rolling average", zorder=5)

ax2.axvline(pd.Timestamp("2021-01-21"), color="darkred", lw=1.2, ls=":", alpha=0.7)
ax2.text(pd.Timestamp("2021-01-21"), ax2.get_ylim()[1] if ax2.get_ylim()[1]>0 else 2200,
         " Peak: 2,049\n 21 Jan", color="darkred", fontsize=8.5, va="top")

# Add Alpha frequency (right axis) for context
ax2r = ax2.twinx()
ax2r.plot(n501_cf["wn_mid_date"], n501_cf["frequency"],
          color="#8e44ad", lw=1.5, ls="-.", alpha=0.7, label="Observed S:N501Y")
ax2r.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax2r.set_ylabel("Alpha frequency (S:N501Y)", fontsize=9, color="#8e44ad")
ax2r.tick_params(axis="y", colors="#8e44ad")
ax2r.set_ylim(-0.05, 1.15)

ax2.set_ylabel("Scotland hospital occupancy (beds)", fontsize=11)
ax2.set_title("B  Hospital occupancy and Alpha frequency (actual timeline)", fontsize=11, loc="left")
lines2, labels2 = ax2.get_legend_handles_labels()
lines2r, labels2r = ax2r.get_legend_handles_labels()
ax2.legend(lines2+lines2r, labels2+labels2r, fontsize=9, loc="upper left")
ax2.grid(axis="y", lw=0.5, alpha=0.3)
ax2.set_xlabel("Date", fontsize=11)
fmt_date_axis(ax2)

# Delta-frequency table inset in ax1
delay_text = "Days to 50% Alpha dominance:\n"
for scen, d50 in date_50.items():
    scen_short = scen.split("(")[0].strip()
    delay_text += f"  {scen_short}: {d50.strftime('%d %b')}\n"
ax1.text(0.99, 0.02, delay_text.strip(), transform=ax1.transAxes,
         fontsize=8, ha="right", va="bottom",
         bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.9))

plt.tight_layout()
fig.savefig(OUT_F / "fig3_counterfactual_projections.png", dpi=180, bbox_inches="tight")
plt.close()
print("  Saved fig3")


# ════════════════════════════════════════════════════════════════════════════
# FIG 4: Growth rate comparison panel
# ════════════════════════════════════════════════════════════════════════════
print("Generating fig4: Growth rate comparison …")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle(f"Alpha Growth Rate Under F5 vs L2:\n"
             f"Only {pct_reduction:.0f}% Slower Under Lockdown — Policy Difference was Modest",
             fontsize=12, fontweight="bold")

# Panel A: logit(freq) vs time, fitted lines
ax = axes[0]
anchor_logit = np.log(0.0066/0.9934)
anchor_date  = pd.Timestamp("2020-11-03")

# F5 observed
f5_obs = n501[(n501["wn_mid_date"] >= "2020-11-03") & (n501["wn_mid_date"] <= "2021-01-04") &
              (n501["frequency"] > 0.005) & (n501["frequency"] < 0.995)].copy()
f5_obs["logit"] = np.log(f5_obs["frequency"] / (1 - f5_obs["frequency"]))
f5_obs["days"]  = (f5_obs["wn_mid_date"] - anchor_date).dt.days

l2_obs = n501[(n501["wn_mid_date"] >= "2021-01-05") & (n501["wn_mid_date"] <= "2021-03-16") &
              (n501["frequency"] > 0.005) & (n501["frequency"] < 0.995)].copy()
l2_obs["logit"] = np.log(l2_obs["frequency"] / (1 - l2_obs["frequency"]))
l2_obs["days"]  = (l2_obs["wn_mid_date"] - anchor_date).dt.days

ax.scatter(f5_obs["days"], f5_obs["logit"], color="#e8735a", s=50, zorder=5, label="Alpha: F5 period")
ax.scatter(l2_obs["days"], l2_obs["logit"], color="#c23b22", s=50, marker="s", zorder=5, label="Alpha: L2 period")

# Fit lines
x_range = np.linspace(-5, 110, 200)
ax.plot(x_range, r_f5 * x_range + anchor_logit, color="#e8735a", lw=2, ls="--",
        label=f"F5 GLM: r={r_f5:.4f}/day\n(doubling {d_f5:.1f}d, CI {d_lo_f5:.1f}–{d_hi_f5:.1f}d)")
# L2 line anchored at Jan 5
l2_anchor_logit = np.log(0.604/0.396)
l2_anchor_days  = (pd.Timestamp("2021-01-05") - anchor_date).days
x_l2 = np.linspace(l2_anchor_days, l2_anchor_days+70, 100)
ax.plot(x_l2, r_l2 * (x_l2-l2_anchor_days) + l2_anchor_logit, color="#c23b22", lw=2, ls="--",
        label=f"L2 GLM: r={r_l2:.4f}/day\n(doubling {d_l2:.1f}d, CI {d_lo_l2:.1f}–{d_hi_l2:.1f}d)")

# L2 vs F5 diff annotation
ax.text(0.5, 0.05,
        f"L2 growth rate = {r_l2/r_f5:.0%} of F5 rate\n"
        f"L2 is only {pct_reduction:.0f}% slower",
        transform=ax.transAxes, ha="center", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", fc="#fff3cd", ec="orange", alpha=0.95))

# x-axis: label key dates
key_days = [(pd.Timestamp(d)-anchor_date).days for d in ["2020-11-03","2020-12-08","2021-01-05","2021-02-02","2021-03-01"]]
key_labels = ["3 Nov","8 Dec","5 Jan","2 Feb","1 Mar"]
ax.set_xticks(key_days)
ax.set_xticklabels(key_labels, fontsize=8, rotation=30)
ax.set_xlabel("Date", fontsize=10)
ax.set_ylabel("logit(S:N501Y frequency)", fontsize=10)
ax.set_title("A  Logistic growth: F5 vs L2 phase", fontsize=11, loc="left")
ax.legend(fontsize=8.5, loc="lower right")
ax.grid(lw=0.5, alpha=0.3)
ax.axvline((pd.Timestamp("2021-01-05")-anchor_date).days, color="#c23b22", lw=1.2, ls=":", alpha=0.7)

# Panel B: bar chart comparing rates + doubling times, with 95% CI error bars
ax = axes[1]
labels = ["Alpha\nunder F5", "Alpha\nunder L2", "B.1.177\ndecline L2"]
slopes = [r_f5, r_l2, abs(r_b1177)]
doub   = [d_f5, d_l2, np.log(2)/abs(r_b1177)]
colors = ["#e8735a", "#c23b22", "#2980b9"]
# CI error bars for Alpha models (none for B.1.177 OLS)
yerr_lo = [r_f5 - ci_lo_f5, r_l2 - ci_lo_l2, 0]
yerr_hi = [ci_hi_f5 - r_f5, ci_hi_l2 - r_l2, 0]
bars = ax.bar(labels, slopes, color=colors, alpha=0.85, edgecolor="white", lw=1.5,
              yerr=[yerr_lo, yerr_hi], capsize=6,
              error_kw={"lw": 1.8, "color": "#333333", "capthick": 1.8})
for bar, d_val, dlo, dhi in zip(bars, doub,
                                 [d_lo_f5, d_lo_l2, None],
                                 [d_hi_f5, d_hi_l2, None]):
    if dlo is not None:
        label_txt = f"{d_val:.1f}d\n({dlo:.1f}–{dhi:.1f}d)"
    else:
        label_txt = f"{d_val:.1f}d"
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(yerr_hi) + 0.0015,
            label_txt, ha="center", va="bottom", fontsize=9.5, fontweight="bold")

ax.set_ylabel("Growth/decline rate (per day)", fontsize=10)
ax.set_title("B  Growth rate comparison\n(bars = point estimate, error bars = 95% CI, labels = doubling/halving days)",
             fontsize=10, loc="left")
ax.axhline(0, color="black", lw=0.8)

# Add annotation for B.1.177
ax.text(2, slopes[2]+0.002, "(halving\n= decline)", ha="center", fontsize=8, color="#2980b9")

ax2b = ax.twinx()
ax2b.set_ylim(0, ax.get_ylim()[1])
ax2b.set_yticks([])
ax.grid(axis="y", lw=0.5, alpha=0.3)

# Highlight the small F5→L2 difference
ax.annotate("", xy=(1, slopes[1]), xytext=(0, slopes[0]),
            arrowprops=dict(arrowstyle="<->", color="black", lw=1.5))
ax.text(0.5, (slopes[0]+slopes[1])/2 + 0.001,
        f"−{pct_reduction:.0f}%", ha="center", fontsize=10, fontweight="bold")

plt.tight_layout()
fig.savefig(OUT_F / "fig4_growth_rate_comparison.png", dpi=180, bbox_inches="tight")
plt.close()
print("  Saved fig4")


# ════════════════════════════════════════════════════════════════════════════
# FIG 5: Lineage displacement — B.1.177 vs Alpha vs total cases
# ════════════════════════════════════════════════════════════════════════════
print("Generating fig5: Lineage displacement …")

lin_comp["wn_mid_date"] = pd.to_datetime(lin_comp["wn_mid_date"])
lin_sub = lin_comp[(lin_comp["wn_mid_date"] >= "2020-09-01") & (lin_comp["wn_mid_date"] <= "2021-04-15")].copy()
lin_sub["group"] = lin_sub["pango_lineage"].apply(
    lambda x: "B.1.1.7 (Alpha)" if str(x).startswith("B.1.1.7")
    else ("B.1.177 (+ sublin.)" if str(x).startswith("B.1.177") else "Other"))
lg = lin_sub.groupby(["wn_mid_date","group"])["n"].sum().reset_index()
pv = lg.pivot(index="wn_mid_date", columns="group", values="n").fillna(0)
pct = pv.div(pv.sum(axis=1), axis=0)

wpm_sub = wpm[(wpm["wn_mid_date"] >= "2020-09-01") & (wpm["wn_mid_date"] <= "2021-04-15")]

fig, axes = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
fig.suptitle("Lineage Displacement in Scotland: B.1.177 → Alpha (B.1.1.7)\n"
             "September 2020 – April 2021", fontsize=12, fontweight="bold")

# Panel A: absolute sequence counts (stacked)
ax = axes[0]
add_period_shading(ax, wpm, "2020-09-01", "2021-04-15")
add_policy_vlines(ax, include=["T1","F5","L2","SL"])

lin_cols = {"B.1.177 (+ sublin.)":"#f08080", "B.1.1.7 (Alpha)":"#6c3483", "Other":"#cccccc"}
bottom = np.zeros(len(pv))
for gname in ["Other","B.1.177 (+ sublin.)","B.1.1.7 (Alpha)"]:
    if gname in pv.columns:
        ax.bar(pv.index, pv[gname].values, bottom=bottom,
               color=lin_cols[gname], alpha=0.85, width=6, label=gname)
        bottom += pv[gname].values

ax.set_ylabel("Sequences per window", fontsize=11)
ax.set_title("A  Absolute sequence counts by lineage group", fontsize=11, loc="left")
ax.legend(fontsize=9, loc="upper left")
ax.grid(axis="y", lw=0.5, alpha=0.3)

# Panel B: proportion + hospital occupancy overlay
ax = axes[1]
add_period_shading(ax, wpm, "2020-09-01", "2021-04-15")
add_policy_vlines(ax, include=["T1","F5","L2","SL"])

bottom = np.zeros(len(pct))
for gname in ["Other","B.1.177 (+ sublin.)","B.1.1.7 (Alpha)"]:
    if gname in pct.columns:
        ax.bar(pct.index, pct[gname].values, bottom=bottom,
               color=lin_cols[gname], alpha=0.85, width=6)
        bottom += pct[gname].values

ax.set_ylabel("Proportion of sequences", fontsize=11)
ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax.set_title("B  Lineage proportions + hospital occupancy", fontsize=11, loc="left")

# Hospital overlay
axr = ax.twinx()
hosp_b = hosp[(hosp["date"] >= "2020-09-01") & (hosp["date"] <= "2021-04-15")]
hosp_7d = hosp_b.set_index("date")["hb_hospital_occupancy"].rolling(7, center=True).mean().reset_index()
axr.plot(hosp_7d["date"], hosp_7d["hb_hospital_occupancy"],
         color="black", lw=2.5, label="Hospital occupancy (7d avg)", zorder=5)
axr.set_ylabel("Hospital occupancy (beds)", fontsize=10)
axr.legend(fontsize=9, loc="upper right")
axr.grid(False)

# Period labels
for p in ["P3","T1","F5","L2","SL","L3"]:
    row_p = POLICY_PERIODS_PD[POLICY_PERIODS_PD["period_code"]==p]
    if row_p.empty: continue
    mid = pd.Timestamp(row_p["start_date"].iloc[0]) + (pd.Timestamp(row_p["end_date"].iloc[0]) - pd.Timestamp(row_p["start_date"].iloc[0]))/2
    if mid < pd.Timestamp("2020-09-01") or mid > pd.Timestamp("2021-04-15"): continue
    ax.text(mid, -0.06, p, ha="center", fontsize=8, fontweight="bold",
            color=PERIOD_COLOURS.get(p,"#444"), transform=ax.get_xaxis_transform())

ax.set_xlabel("Window mid-date", fontsize=11)
ax.grid(axis="y", lw=0.5, alpha=0.3)
fmt_date_axis(ax)

plt.tight_layout()
fig.savefig(OUT_F / "fig5_lineage_displacement.png", dpi=180, bbox_inches="tight")
plt.close()
print("  Saved fig5")

print("\nStage 3 complete. All figures saved to:", OUT_F)
