"""Supplementary Part 2 descriptive question tables.

This script answers secondary, thesis-introduction questions that are useful
for interpretation but are not part of the core Part 2 figure set:

1. Is vaccination-status completeness patterned by wave, age, SIMD, or
   sequencing fraction?
2. Do vaccination-profile groups differ in cluster structure and demographic
   mixing within waves?
3. Are booster coverage and dose recency gradients consistent across SIMD
   domains within waves?

Run from the repository root:

    conda run -n PhD python part2/supplementary_questions.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _bootstrap_repo_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in [here, *here.parents]:
        if (candidate / "config.yaml").exists():
            root_str = str(candidate)
            if root_str not in sys.path:
                sys.path.insert(0, root_str)
            return candidate
    raise FileNotFoundError("Could not locate config.yaml.")


ROOT = _bootstrap_repo_root()
TABLE_DIR = ROOT / "part2" / "tables"
CACHE_DIR = ROOT / "part2" / "cache"

from utils.data import load_analysis_columns_pandas  # noqa: E402
from cluster_characterisation import (  # noqa: E402
    DEFAULT_SIMD_RANK_MAX,
    PRIMARY_RESOLUTION,
    QC_DEFAULT,
    SIMD_DOMAIN_ORDER,
    SIMD_DOMAIN_SPECS,
    SIMD_LABEL_ORDER,
    WAVE_ORDER,
    add_sequence_simd_domain_quintiles,
    assign_wave,
    vaccination_age_group,
)


CASE_COLUMNS = [
    "sequence_id",
    "cluster_id",
    "window_idx",
    "wn_mid_date",
    "wn_prop_sequenced",
    "collection_date",
    "pango_lineage",
    "age_band",
    "sex",
    "is_vaccinated",
    "dz_simd_rank",
    "dz_simd_income_rank",
    "dz_simd_employment_rank",
    "dz_simd_education_rank",
    "dz_simd_health_rank",
    "dz_simd_access_rank",
    "dz_simd_crime_rank",
    "dz_simd_housing_rank",
]


def _simd_quintile_label(quintile: pd.Series) -> pd.Categorical:
    labels = {
        1: "1 most deprived",
        2: "2",
        3: "3",
        4: "4",
        5: "5 least deprived",
    }
    return pd.Categorical(
        quintile.map(labels),
        categories=SIMD_LABEL_ORDER,
        ordered=True,
    )


def load_case_rows() -> pd.DataFrame:
    """Load one row per unique sequence/case with vaccination completeness fields."""
    seq = load_analysis_columns_pandas(
        columns=CASE_COLUMNS,
        resolution=PRIMARY_RESOLUTION,
        qc=QC_DEFAULT,
    )
    seq["collection_date"] = pd.to_datetime(seq["collection_date"])
    seq["wn_mid_date"] = pd.to_datetime(seq["wn_mid_date"])
    seq["is_vaccinated"] = pd.to_numeric(seq["is_vaccinated"], errors="coerce")
    seq["wn_prop_sequenced"] = pd.to_numeric(seq["wn_prop_sequenced"], errors="coerce")

    cases = (
        seq.sort_values(["sequence_id", "window_idx", "collection_date"])
        .drop_duplicates("sequence_id")
        .copy()
    )
    cases["wave_group"] = pd.Categorical(
        cases["pango_lineage"].astype(str).map(assign_wave),
        categories=WAVE_ORDER,
        ordered=True,
    )
    cases["vaccination_age_group"] = cases["age_band"].map(vaccination_age_group)
    cases = add_sequence_simd_domain_quintiles(
        cases,
        simd_rank_max=DEFAULT_SIMD_RANK_MAX,
    )
    cases["simd_quintile_label"] = _simd_quintile_label(cases["simd_quintile"])
    cases["vaccination_known"] = cases["is_vaccinated"].notna()

    quartile_labels = [
        "Q1 lowest sequencing fraction",
        "Q2",
        "Q3",
        "Q4 highest sequencing fraction",
    ]
    cases["sequencing_fraction_quartile"] = pd.qcut(
        cases["wn_prop_sequenced"],
        q=4,
        labels=quartile_labels,
        duplicates="drop",
    )
    return cases


def _summarise_case_group(
    cases: pd.DataFrame,
    group_cols: list[str],
    stratum_type: str,
) -> pd.DataFrame:
    if group_cols:
        grouped = cases.groupby(group_cols, observed=True, dropna=False)
        out = grouped.agg(
            n_cases=("sequence_id", "nunique"),
            n_vaccination_known=("vaccination_known", "sum"),
            n_vaccinated=("is_vaccinated", lambda x: int((x == 1).sum())),
            mean_window_seq_fraction=("wn_prop_sequenced", "mean"),
            median_window_seq_fraction=("wn_prop_sequenced", "median"),
        ).reset_index()
    else:
        out = pd.DataFrame(
            [
                {
                    "n_cases": cases["sequence_id"].nunique(),
                    "n_vaccination_known": int(cases["vaccination_known"].sum()),
                    "n_vaccinated": int((cases["is_vaccinated"] == 1).sum()),
                    "mean_window_seq_fraction": cases["wn_prop_sequenced"].mean(),
                    "median_window_seq_fraction": cases["wn_prop_sequenced"].median(),
                }
            ]
        )

    out.insert(0, "stratum_type", stratum_type)
    out["pct_vaccination_known"] = 100.0 * out["n_vaccination_known"] / out["n_cases"]
    out["pct_vaccinated_among_known"] = np.divide(
        100.0 * out["n_vaccinated"],
        out["n_vaccination_known"],
        out=np.full(len(out), np.nan),
        where=out["n_vaccination_known"].to_numpy() > 0,
    )
    return out


def vaccination_missingness_summary(cases: pd.DataFrame) -> pd.DataFrame:
    """Summarise vaccination-status completeness across key strata."""
    specs = [
        ("overall", []),
        ("wave", ["wave_group"]),
        ("age_group", ["vaccination_age_group"]),
        ("sex", ["sex"]),
        ("overall_simd", ["simd_quintile_label"]),
        ("sequencing_fraction_quartile", ["sequencing_fraction_quartile"]),
        ("wave_by_simd", ["wave_group", "simd_quintile_label"]),
        ("wave_by_age_group", ["wave_group", "vaccination_age_group"]),
        ("wave_by_sequencing_fraction", ["wave_group", "sequencing_fraction_quartile"]),
    ]
    frames = [
        _summarise_case_group(cases, group_cols, stratum_type)
        for stratum_type, group_cols in specs
    ]
    return pd.concat(frames, ignore_index=True)


def load_joined_cluster_tables() -> pd.DataFrame:
    """Join Part 2 vaccination cluster table to Part 1 category/mixing table."""
    vacc_cols = [
        "cluster_id",
        "wave_group",
        "cluster_vaccination_profile",
        "vaccination_mixing_category",
        "cluster_prop_vaccinated",
        "cluster_prop_boosted_vaccinated_members",
        "mean_days_since_vaccination",
        "median_days_since_vaccination",
        "mean_dose_vaccinated_members",
    ]
    vacc_cols.extend(
        spec["quintile_col"]
        for spec in SIMD_DOMAIN_SPECS.values()
        if "quintile_col" in spec
    )

    cat_cols = [
        "cluster_id",
        "cluster_size",
        "cluster_n_datazones",
        "cluster_size_category",
        "geographic_dispersion_category",
        "simd_mixing_category",
        "age_mixing_category",
        "sex_mixing_category",
        "profile_mixing_category",
        "simd_excess_discordance",
        "age_excess_discordance",
        "sex_excess_discordance",
        "profile_excess_discordance",
    ]

    vacc = pd.read_parquet(CACHE_DIR / "vaccination_cluster_table.parquet", columns=vacc_cols)
    cats = pd.read_parquet(CACHE_DIR / "cluster_categories.parquet", columns=cat_cols)
    joined = cats.merge(vacc, on="cluster_id", how="left", validate="one_to_one")
    joined["wave_group"] = pd.Categorical(
        joined["wave_group"],
        categories=WAVE_ORDER,
        ordered=True,
    )
    return joined


def _pct_eq(values: pd.Series, label: str) -> float:
    values = values.astype("string")
    return 100.0 * values.eq(label).mean()


def _summarise_cluster_groups(clusters: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    grouped = clusters.groupby(group_cols, observed=True, dropna=False)
    return grouped.agg(
        n_clusters=("cluster_id", "size"),
        n_cases=("cluster_size", "sum"),
        median_cluster_size=("cluster_size", "median"),
        median_datazones=("cluster_n_datazones", "median"),
        pct_large_or_very_large=(
            "cluster_size_category",
            lambda x: _pct_eq(x, "large") + _pct_eq(x, "very large"),
        ),
        pct_large_or_very_large_dispersion=(
            "geographic_dispersion_category",
            lambda x: _pct_eq(x, "large dispersion") + _pct_eq(x, "very large dispersion"),
        ),
        mean_prop_vaccinated=("cluster_prop_vaccinated", "mean"),
        mean_prop_boosted_vaccinated_members=(
            "cluster_prop_boosted_vaccinated_members",
            "mean",
        ),
        median_days_since_vaccination=("median_days_since_vaccination", "median"),
        mean_simd_excess=("simd_excess_discordance", "mean"),
        mean_age_excess=("age_excess_discordance", "mean"),
        mean_sex_excess=("sex_excess_discordance", "mean"),
        mean_profile_excess=("profile_excess_discordance", "mean"),
        pct_simd_less_mix=("simd_mixing_category", lambda x: _pct_eq(x, "less mix")),
        pct_age_more_mix=("age_mixing_category", lambda x: _pct_eq(x, "more mix")),
        pct_sex_more_mix=("sex_mixing_category", lambda x: _pct_eq(x, "more mix")),
        pct_profile_more_mix=("profile_mixing_category", lambda x: _pct_eq(x, "more mix")),
    ).reset_index()


def simd_domain_vaccination_gradients(
    clusters: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute wave-specific Q5-Q1 vaccination gradients for each SIMD domain."""
    quintile_frames = []
    gradient_records = []

    for domain in SIMD_DOMAIN_ORDER:
        spec = SIMD_DOMAIN_SPECS[domain]
        qcol = spec["quintile_col"]
        label = spec["label"]
        work = clusters.dropna(subset=[qcol]).copy()
        work[qcol] = pd.to_numeric(work[qcol], errors="coerce")

        by_q = (
            work.groupby(["wave_group", qcol], observed=True)
            .agg(
                n_clusters=("cluster_id", "size"),
                mean_prop_vaccinated=("cluster_prop_vaccinated", "mean"),
                mean_prop_boosted_vaccinated_members=(
                    "cluster_prop_boosted_vaccinated_members",
                    "mean",
                ),
                mean_days_since_vaccination=("mean_days_since_vaccination", "mean"),
                median_days_since_vaccination=("median_days_since_vaccination", "median"),
            )
            .reset_index()
            .rename(columns={qcol: "simd_quintile"})
        )
        by_q.insert(0, "simd_domain", domain)
        by_q.insert(1, "simd_domain_label", label)
        quintile_frames.append(by_q)

        for wave, sub in by_q.groupby("wave_group", observed=True):
            q1 = sub[sub["simd_quintile"] == 1]
            q5 = sub[sub["simd_quintile"] == 5]
            if q1.empty or q5.empty:
                continue
            q1 = q1.iloc[0]
            q5 = q5.iloc[0]
            gradient_records.append(
                {
                    "wave_group": wave,
                    "simd_domain": domain,
                    "simd_domain_label": label,
                    "n_clusters_q1": int(q1["n_clusters"]),
                    "n_clusters_q5": int(q5["n_clusters"]),
                    "q5_minus_q1_prop_vaccinated_pp": 100.0
                    * (q5["mean_prop_vaccinated"] - q1["mean_prop_vaccinated"]),
                    "q5_minus_q1_booster_among_vaccinated_pp": 100.0
                    * (
                        q5["mean_prop_boosted_vaccinated_members"]
                        - q1["mean_prop_boosted_vaccinated_members"]
                    ),
                    "q5_minus_q1_days_since_vaccination": (
                        q5["mean_days_since_vaccination"]
                        - q1["mean_days_since_vaccination"]
                    ),
                }
            )

    return pd.concat(quintile_frames, ignore_index=True), pd.DataFrame(gradient_records)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    print("Part 2 supplementary questions")
    print("=" * 34)

    print("Loading case-level rows for vaccination completeness summaries...")
    cases = load_case_rows()
    missingness = vaccination_missingness_summary(cases)
    missingness_path = TABLE_DIR / "supp_vaccination_missingness_summary.csv"
    missingness.to_csv(missingness_path, index=False)
    print(f"  Saved {missingness_path.relative_to(ROOT)}")

    print("Joining vaccination profiles to cluster mixing categories...")
    clusters = load_joined_cluster_tables()

    profile_summary = _summarise_cluster_groups(
        clusters,
        ["wave_group", "cluster_vaccination_profile"],
    )
    profile_path = TABLE_DIR / "supp_vaccination_profile_cluster_mixing_summary.csv"
    profile_summary.to_csv(profile_path, index=False)
    print(f"  Saved {profile_path.relative_to(ROOT)}")

    mixing_summary = _summarise_cluster_groups(
        clusters,
        ["wave_group", "vaccination_mixing_category"],
    )
    mixing_path = TABLE_DIR / "supp_vaccination_mixing_demographic_summary.csv"
    mixing_summary.to_csv(mixing_path, index=False)
    print(f"  Saved {mixing_path.relative_to(ROOT)}")

    by_quintile, gradients = simd_domain_vaccination_gradients(clusters)
    by_quintile_path = TABLE_DIR / "supp_simd_domain_vaccination_by_quintile.csv"
    gradients_path = TABLE_DIR / "supp_simd_domain_vaccination_gradients.csv"
    by_quintile.to_csv(by_quintile_path, index=False)
    gradients.to_csv(gradients_path, index=False)
    print(f"  Saved {by_quintile_path.relative_to(ROOT)}")
    print(f"  Saved {gradients_path.relative_to(ROOT)}")

    print("\nQuick checks:")
    overall = missingness[missingness["stratum_type"] == "overall"].iloc[0]
    print(
        "  Vaccination status known for "
        f"{overall['pct_vaccination_known']:.1f}% of unique sequenced cases."
    )
    print(
        "  Joined cluster table covers "
        f"{len(clusters):,} non-singleton clusters with vaccination and mixing fields."
    )
    print("Part 2 supplementary tables complete.")


if __name__ == "__main__":
    main()
