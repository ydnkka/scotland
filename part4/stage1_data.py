"""
Part 4 – Stage 1: Data extraction
Extracts the Alpha seeding chain, mutation frequency trajectory,
cluster evolution, and hospital occupancy data needed for analysis and figures.
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from utils.policy import attach_period_pandas

OUT = Path(__file__).resolve().parent / "tables"
OUT.mkdir(parents=True, exist_ok=True)

NEXTCLADE_TSV = REPO / "data/raw/cog_all_scotland_nextclade.tsv"
DATA_PARQUET  = REPO / "data/processed/scotland_clustering_analysis_dataset.parquet"
HB_PARQUET    = REPO / "data/processed/scotland_hb_daily_trends.parquet"

# ── 1. Load core cluster data ─────────────────────────────────────────────
print("Loading cluster data …")
cols = ["sequence_id","cluster_id","cluster_size","cluster_n_datazones",
        "cluster_start_date","cluster_end_date","collection_date",
        "pango_lineage","resolution","nextclade_qc",
        "age_band","age_midpoint","sex","is_female","is_vaccinated",
        "dz_health_board","dz_local_authority","dz_simd_quintile",
        "window_id","window_idx","wn_mid_date","test_reason","test_type"]
df = pq.read_table(str(DATA_PARQUET), columns=cols).to_pandas()
df = df[(df["resolution"]==0.3) & (df["nextclade_qc"]=="good")].copy()
df["collection_date"] = pd.to_datetime(df["collection_date"])
df["wn_mid_date"]     = pd.to_datetime(df["wn_mid_date"])
df = attach_period_pandas(df, "wn_mid_date")
print(f"  {len(df):,} rows")

# Save seq→window map
seq_map = df[["sequence_id","window_id","window_idx","wn_mid_date","policy_period"]].drop_duplicates("sequence_id")
seq_map.to_parquet(OUT / "seq_window_map.parquet", index=False)

# Window → period map
wpm = (df.groupby(["window_id","window_idx","wn_mid_date","policy_period",
                    "policy_period_label","policy_intensity"])
         .size().reset_index(name="n_seqs").sort_values("window_idx"))
wpm.to_csv(OUT / "part4_window_period_map.csv", index=False)

# ── 2. Alpha seeding chain ───────────────────────────────────────────────
print("Building Alpha seeding chain …")
alpha = df[df["pango_lineage"].str.startswith("B.1.1.7", na=False)]
alpha_clust = (alpha.groupby("cluster_id").agg(
    size              =("cluster_size","first"),
    n_datazones       =("cluster_n_datazones","first"),
    first_seq_date    =("collection_date","min"),
    last_seq_date     =("collection_date","max"),
    window_id         =("window_id","first"),
    window_idx        =("window_idx","first"),
    wn_mid_date       =("wn_mid_date","first"),
    primary_hb        =("dz_health_board", lambda x: x.mode()[0]),
    n_health_boards   =("dz_health_board","nunique"),
    n_local_auth      =("dz_local_authority","nunique"),
    local_authorities =("dz_local_authority", lambda x: "; ".join(sorted(x.unique()))),
    test_reason_top   =("test_reason", lambda x: x.value_counts(dropna=False).index[0]
                         if len(x)>0 else "unknown"),
    policy_period     =("policy_period","first"),
).reset_index().sort_values("first_seq_date"))

# Sequence overlap between consecutive Alpha clusters (chain evidence)
early_chain = alpha_clust[
    (alpha_clust["first_seq_date"] <= "2020-12-15") &
    (alpha_clust["size"] >= 1)
]["cluster_id"].tolist()

overlap_rows = []
cluster_seqs = {cid: set(alpha[alpha["cluster_id"]==cid]["sequence_id"]) for cid in early_chain}
for i, c1 in enumerate(early_chain):
    for c2 in early_chain[i+1:]:
        ov = cluster_seqs[c1] & cluster_seqs[c2]
        if ov:
            overlap_rows.append({
                "cluster_a": c1, "cluster_b": c2,
                "n_shared": len(ov),
                "pct_of_a": 100*len(ov)/max(len(cluster_seqs[c1]),1),
                "pct_of_b": 100*len(ov)/max(len(cluster_seqs[c2]),1),
            })

alpha_clust.to_csv(OUT / "part4_alpha_cluster_chain.csv", index=False)
pd.DataFrame(overlap_rows).to_csv(OUT / "part4_alpha_chain_overlaps.csv", index=False)
print(f"  {len(alpha_clust)} Alpha clusters, {len(overlap_rows)} pairwise overlaps in early chain")

# Early Alpha demographic/area summaries for narrative interpretation.
# These are cluster-level summaries over overlapping windows, not unique case
# counts. They are intended to describe the early surveillance signal.
def mode_or_na(values):
    clean = values.dropna()
    if clean.empty:
        return pd.NA
    return clean.astype(str).value_counts().index[0]


def vaccination_profile(values):
    clean = pd.to_numeric(values, errors="coerce").dropna()
    if clean.empty:
        return "unknown"
    prop = clean.mean()
    if prop <= 0:
        return "none vaccinated"
    if prop >= 1:
        return "all vaccinated"
    return "mixed vaccination"


def weighted_profile(sub, col, weight_col="cluster_size", max_items=3):
    valid = sub.dropna(subset=[col])
    if valid.empty:
        return ""
    totals = valid.groupby(col, dropna=False)[weight_col].sum().sort_values(ascending=False)
    denom = totals.sum()
    return "; ".join(
        f"{idx} ({100 * val / denom:.1f}%)"
        for idx, val in totals.head(max_items).items()
    )


def count_profile(sub, col, max_items=3):
    valid = sub.dropna(subset=[col])
    if valid.empty:
        return ""
    totals = valid[col].astype(str).value_counts()
    denom = totals.sum()
    return "; ".join(
        f"{idx} ({100 * val / denom:.1f}%)"
        for idx, val in totals.head(max_items).items()
    )


alpha_cluster_demo = (
    alpha.groupby("cluster_id")
    .agg(
        cluster_size=("sequence_id", "nunique"),
        window_id=("window_id", "first"),
        window_idx=("window_idx", "first"),
        first_seq_date=("collection_date", "min"),
        last_seq_date=("collection_date", "max"),
        mean_age_midpoint=("age_midpoint", "mean"),
        predominant_age_band=("age_band", mode_or_na),
        predominant_sex=("sex", mode_or_na),
        cluster_prop_female=("is_female", "mean"),
        predominant_simd_quintile=("dz_simd_quintile", mode_or_na),
        n_vaccination_known=("is_vaccinated", "count"),
        cluster_prop_vaccinated=("is_vaccinated", "mean"),
        vaccination_profile=("is_vaccinated", vaccination_profile),
    )
    .reset_index()
)

def simd_label(value):
    if pd.isna(value):
        return pd.NA
    try:
        q = int(float(value))
    except (TypeError, ValueError):
        return str(value)
    if q == 1:
        return "SIMD 1 most deprived"
    if q == 5:
        return "SIMD 5 least deprived"
    return f"SIMD {q}"


alpha_cluster_demo["predominant_simd_label"] = (
    alpha_cluster_demo["predominant_simd_quintile"].map(simd_label)
)

phase_specs = [
    ("Cryptic GGC chain", "W016-W021", [f"W{i:03d}" for i in range(16, 22)]),
    ("Multi-region expansion", "W022-W024", [f"W{i:03d}" for i in range(22, 25)]),
    ("F5/L2 bridge", "W025", ["W025"]),
]
phase_rows = []
for phase, windows_label, windows in phase_specs:
    sub = alpha_cluster_demo[alpha_cluster_demo["window_id"].isin(windows)].copy()
    seq_sub = alpha[alpha["window_id"].isin(windows)].copy()
    if sub.empty:
        continue
    unique_seq = seq_sub.drop_duplicates("sequence_id").copy()
    unique_seq["simd_label"] = unique_seq["dz_simd_quintile"].map(simd_label)
    unique_seq["vaccination_status"] = pd.Series(
        np.select(
            [
                unique_seq["is_vaccinated"].eq(1),
                unique_seq["is_vaccinated"].eq(0),
            ],
            ["vaccinated", "unvaccinated"],
            default="unknown",
        ),
        index=unique_seq.index,
    )
    age_valid = sub.dropna(subset=["mean_age_midpoint"])
    age_weight = age_valid["cluster_size"].sum()
    weighted_age = (
        np.average(age_valid["mean_age_midpoint"], weights=age_valid["cluster_size"])
        if age_weight else np.nan
    )
    unique_mean_age = unique_seq["age_midpoint"].mean()
    phase_rows.append({
        "phase": phase,
        "windows": windows_label,
        "n_clusters": len(sub),
        "summed_cluster_size": int(sub["cluster_size"].sum()),
        "unique_sequences": int(seq_sub["sequence_id"].nunique()),
        "overlap_duplicate_memberships": int(
            sub["cluster_size"].sum() - seq_sub["sequence_id"].nunique()
        ),
        "weighted_mean_age": round(float(weighted_age), 1) if pd.notna(weighted_age) else np.nan,
        "dominant_age_bands": weighted_profile(sub, "predominant_age_band"),
        "predominant_sex_profile": weighted_profile(sub, "predominant_sex", max_items=2),
        "simd_quintile_profile": weighted_profile(sub, "predominant_simd_label"),
        "vaccination_profile": weighted_profile(sub, "vaccination_profile"),
        "unique_mean_age": round(float(unique_mean_age), 1) if pd.notna(unique_mean_age) else np.nan,
        "unique_age_band_profile": count_profile(unique_seq, "age_band"),
        "unique_sex_profile": count_profile(unique_seq, "sex", max_items=2),
        "unique_simd_quintile_profile": count_profile(unique_seq, "simd_label"),
        "unique_vaccination_status_profile": count_profile(unique_seq, "vaccination_status"),
        "first_seq_date": sub["first_seq_date"].min(),
        "last_seq_date": sub["last_seq_date"].max(),
    })

pd.DataFrame(phase_rows).to_csv(
    OUT / "part4_alpha_phase_demographic_summary.csv", index=False
)

# Wider context: all Alpha cluster sizes by week (for growth plot)
alpha_weekly = (alpha.groupby(["wn_mid_date","cluster_id"])
                .agg(size=("cluster_size","first"), n_hb=("dz_health_board","nunique"))
                .reset_index())
alpha_weekly.to_csv(OUT / "part4_alpha_clusters_weekly.csv", index=False)

# ── 3. B.1.177 cluster sizes by week (for displacement context) ──────────
b1177 = df[df["pango_lineage"].str.startswith("B.1.177", na=False)]
b1177_weekly = (b1177.groupby(["wn_mid_date","cluster_id"])
                .agg(size=("cluster_size","first"))
                .reset_index())
b1177_weekly.to_csv(OUT / "part4_b1177_clusters_weekly.csv", index=False)

# ── 4. Lineage composition per window ────────────────────────────────────
lin_comp = (df.groupby(["wn_mid_date","pango_lineage"])["sequence_id"]
              .count().reset_index(name="n"))
lin_comp.to_csv(OUT / "part4_lineage_composition.csv", index=False)

# ── 5. Hospital occupancy (Scotland total) ───────────────────────────────
print("Loading hospital data …")
hb = pq.read_table(str(HB_PARQUET)).to_pandas()
hb["date"] = pd.to_datetime(hb["date"])
scot_hosp = (hb.groupby("date")[["hb_daily_positive","hb_hospital_occupancy"]]
               .sum().reset_index())
scot_hosp = scot_hosp[(scot_hosp["date"] >= "2020-09-01") & (scot_hosp["date"] <= "2021-04-30")]
scot_hosp.to_csv(OUT / "part4_scotland_hospital.csv", index=False)
print(f"  Hospital data: {len(scot_hosp)} days")

# ── 6. Save focal sequence IDs for stage 2 ───────────────────────────────
all_seqids = set(df["sequence_id"].unique())
with open(OUT / "all_seqids.txt", "w") as f:
    f.write("\n".join(all_seqids))
print(f"  All seq IDs: {len(all_seqids):,}")

print("\nStage 1 complete.")
