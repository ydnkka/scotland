#!/usr/bin/env python3
"""
Consolidate per-group cluster parquets and merge with all metadata.

Steps:
1.  Concatenate all per-group long parquets from cluster_long_dir.
2.  Add singleton assignments (lineage groups of size 1, not processed by tn93).
3.  Merge with sequence metadata, SIMD, testing, vaccination.
4.  Compute window-level aggregates using theoretical window boundaries.
5.  Attach per-datazone cumulative vaccination coverage (merge_asof, backward).
6.  Attach per-datazone cumulative sequencing fraction at time of sampling.
7.  Join 7-day rolling test positivity per datazone.
8.  Compute cluster-level descriptors (size, geographic spread, temporal duration).
9.  Compute derived individual- and area-level variables.
10. Join health-board daily trends (hospital, ICU, reinfections) via merge_asof.
11. Anonymise patient_id.
12. Write final analysis dataset parquet.

Output columns (~80 total):

Window-level:
    window_idx, window_id, wn_start_date, wn_mid_date, wn_end_date,
    wn_no_sequences, wn_positive_tests, wn_prop_sequenced

Sequence / cluster identifiers:
    sequence_id, patient_id, resolution, cluster_id

Cluster descriptors:
    cluster_size, cluster_n_datazones, cluster_start_date, cluster_end_date,
    cluster_duration_days

Sample-level:
    collection_date, datazone, dz_xcoord, dz_ycoord,
    sex, is_female, age_band, age_group, age_midpoint,
    is_vaccinated, vacc_dose_number, vacc_date_prior,
    vacc_product_name, vacc_booster, days_since_vaccination,
    test_type, test_reason_raw, test_reason, s_gene_status,
    policy_period, policy_period_label, policy_era,
    pango_lineage, clade, who_voc, nextclade_qc

Datazone sociodemographic:
    dz_population, dz_working_age_population,
    dz_simd_rank, dz_simd_quintile, dz_simd_decile, dz_simd_vigintile,
    dz_simd_income_rank, dz_simd_employment_rank, dz_simd_education_rank,
    dz_simd_health_rank, dz_simd_access_rank, dz_simd_crime_rank,
    dz_simd_housing_rank,
    dz_urban_rural_class, dz_local_authority, dz_local_authority_code,
    dz_health_board, dz_health_board_code

Datazone daily testing (on collection_date):
    dz_total_tests, dz_positive_tests, dz_negative_tests,
    dz_pcr_positive_tests, dz_lfd_positive_tests, dz_care_home_tests,
    dz_test_positivity, dz_7d_test_positivity

Datazone vaccination (daily new + cumulative):
    dz_total_vaccinated, dz_cum_vaccinated, dz_cum_prop_vaccinated

Datazone cumulative surveillance:
    dz_cum_sequences, dz_cum_positive_tests, dz_cum_prop_sequenced,
    dz_cum_incidence_per_capita

Health-board daily trends (on or before collection_date):
    hb_daily_positive, hb_cumulative_positive,
    hb_hospital_admissions, hb_hospital_occupancy,
    hb_icu_admissions, hb_icu_occupancy_lt28d, hb_icu_occupancy_ge28d,
    hb_daily_reinfections, hb_reinfection_rate

Usage:
    python3 method/05_consolidate.py
    python3 method/05_consolidate.py --config config.yaml --root /path/to/repo
"""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import timedelta
from pathlib import Path

import pandas as pd
import yaml


def setup_logging(level: str = "INFO") -> None:
    """Configure timestamped console logging for consolidation."""
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def load_config(path: Path) -> dict:
    """Load the YAML pipeline configuration file."""
    with open(path) as f:
        return yaml.safe_load(f)


def filter_nextclade_qc(
    metadata: pd.DataFrame,
    required_status: str,
) -> pd.DataFrame:
    """Return only the configured Nextclade QC cohort."""
    if "nextclade_qc" not in metadata.columns:
        raise KeyError(
            "Sequence metadata lacks 'nextclade_qc'. Rebuild it with "
            "method/01_prep_metadata.py."
        )
    qc_status = str(required_status).strip()
    if not qc_status:
        raise ValueError("tn93.nextclade_qc must be a non-empty status.")
    qc_matches = (
        metadata["nextclade_qc"]
        .astype("string")
        .str.strip()
        .str.casefold()
        .eq(qc_status.casefold())
    )
    n_before_qc = len(metadata)
    filtered = metadata.loc[qc_matches].copy()
    logging.info(
        "Nextclade QC filter (%s): retained %d/%d metadata rows; dropped %d.",
        qc_status,
        len(filtered),
        n_before_qc,
        n_before_qc - len(filtered),
    )
    if filtered.empty:
        raise ValueError(
            f"No sequence metadata rows have Nextclade QC status {qc_status!r}."
        )
    return filtered


def build_windows(
    min_date: pd.Timestamp,
    max_date: pd.Timestamp,
    window_size: timedelta,
    step: timedelta,
) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Build half-open rolling windows from minimum to maximum date."""
    if max_date < min_date + window_size:
        return []
    starts = pd.date_range(start=min_date, end=max_date - window_size, freq=step)
    return [(s, s + window_size) for s in starts]


def load_cluster_parquets(cluster_long_dir: Path) -> pd.DataFrame:
    """Load and concatenate all per-group cluster assignment parquet files."""
    files = sorted(cluster_long_dir.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"No parquet files in {cluster_long_dir}")
    logging.info("Loading %d cluster parquets from %s", len(files), cluster_long_dir)
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


def build_singletons(
    metadata: pd.DataFrame,
    windows: list[tuple[pd.Timestamp, pd.Timestamp]],
    resolutions: list[float],
) -> pd.DataFrame:
    """Assign lineage-window groups of size one to singleton clusters."""
    rows = []
    for i, (start, end) in enumerate(windows, start=1):
        window_id = f"W{str(i).zfill(3)}"
        wdf = metadata[
            (metadata["collection_date"] >= start) & (metadata["collection_date"] < end)
        ]
        for lin, gdf in wdf.groupby("pango_lineage", sort=False):
            if gdf["sequence_id"].nunique() == 1:
                sid = gdf["sequence_id"].iloc[0]
                for res in resolutions:
                    rows.append(
                        {
                            "sequence_id": sid,
                            "window_id": window_id,
                            "resolution": res,
                            "cluster_id": f"{window_id}|{lin}|R{res}|S0",
                        }
                    )
    return pd.DataFrame(rows)


def build_window_info(
    metadata: pd.DataFrame,
    testing: pd.DataFrame,
    windows: list[tuple[pd.Timestamp, pd.Timestamp]],
) -> pd.DataFrame:
    """Compute per-window sequence counts and sequencing fractions."""
    rows = []
    for i, (start, end) in enumerate(windows, start=1):
        window_id = f"W{str(i).zfill(3)}"
        wdf = metadata[
            (metadata["collection_date"] >= start) & (metadata["collection_date"] < end)
        ]
        tdf = testing[
            (testing["collection_date"] >= start) & (testing["collection_date"] < end)
        ]
        pos_tests = tdf["dz_positive_tests"].sum()
        rows.append(
            {
                "window_id": window_id,
                "wn_no_sequences": wdf["sequence_id"].nunique(),
                "wn_positive_tests": pos_tests,
                "wn_prop_sequenced": wdf["sequence_id"].nunique() / pos_tests
                if pos_tests > 0
                else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def build_dz_cumulative_sequencing(
    metadata: pd.DataFrame,
    testing: pd.DataFrame,
) -> pd.DataFrame:
    """Per-(datazone, date): cumulative sequences, positive tests, and sequencing fraction.

    For each row in the final dataset the returned frame provides:

    dz_cum_sequences
        Total unique genomes collected from that datazone on or before this date.
    dz_cum_positive_tests
        Total positive PCR/LFD tests recorded in that datazone on or before this date.
    dz_cum_prop_sequenced
        dz_cum_sequences / dz_cum_positive_tests — the fraction of the local epidemic
        that has been genomically characterised by the time this sample was taken.

    The two sources (metadata and testing) may record disjoint date sets for a given
    datazone, so an outer join followed by per-datazone forward-fill is used to build
    a unified timeline.  Rows with zero cumulative positive tests yield NaN for the
    proportion.

    Note: cumulative sequences are counted from the QC-filtered metadata before any
    SIMD filtering, so the numerator reflects retained good-quality surveillance data.
    """
    seq_daily = (
        metadata.groupby(["datazone", "collection_date"])["sequence_id"]
        .nunique()
        .reset_index(name="_new_seqs")
        .sort_values(["datazone", "collection_date"])
    )
    seq_daily["dz_cum_sequences"] = seq_daily.groupby("datazone")["_new_seqs"].cumsum()
    seq_daily = seq_daily.drop(columns=["_new_seqs"])

    tests_daily = (
        testing[["datazone", "collection_date", "dz_positive_tests"]]
        .sort_values(["datazone", "collection_date"])
        .copy()
    )
    tests_daily["dz_cum_positive_tests"] = tests_daily.groupby("datazone")[
        "dz_positive_tests"
    ].cumsum()
    tests_daily = tests_daily.drop(columns=["dz_positive_tests"])

    # Outer-join to cover dates present in either source, then forward-fill so every
    # date inherits the most recent cumulative count from earlier in the timeline.
    combined = seq_daily.merge(
        tests_daily, on=["datazone", "collection_date"], how="outer"
    ).sort_values(["datazone", "collection_date"])
    combined["dz_cum_sequences"] = (
        combined.groupby("datazone")["dz_cum_sequences"].ffill().fillna(0)
    )
    combined["dz_cum_positive_tests"] = (
        combined.groupby("datazone")["dz_cum_positive_tests"].ffill().fillna(0)
    )
    combined["dz_cum_prop_sequenced"] = combined["dz_cum_sequences"] / combined[
        "dz_cum_positive_tests"
    ].replace(0.0, float("nan"))
    return combined[
        [
            "datazone",
            "collection_date",
            "dz_cum_sequences",
            "dz_cum_positive_tests",
            "dz_cum_prop_sequenced",
        ]
    ]


def build_dz_cumulative_vaccination(
    vaccination: pd.DataFrame,
    simd: pd.DataFrame,
) -> pd.DataFrame:
    """Running cumulative vaccination total and population coverage per (datazone, date).

    ``vaccination`` is the daily aggregate produced by ``prep_vaccination`` in
    01_prep_metadata.py: one row per (vaccination_date, datazone) where
    ``dz_total_vaccinated`` is the count of *unique patients* vaccinated in that
    datazone on that day.

    Important caveat: the daily aggregate counts unique patients per day, so a
    patient receiving dose 2 is counted again on the dose-2 date.  The cumulative
    sum therefore slightly over-estimates the number of distinct individuals ever
    vaccinated (it counts vaccination events, not people).  A precise count would
    require tracking the first-dose date per patient, which is available in the raw
    vaccination CSV but not in this aggregate.

    Returns a frame sorted by (datazone, vaccination_date) ready for
    ``pd.merge_asof(direction='backward')`` so each sequence inherits the
    cumulative coverage as of its collection_date.
    """
    pop = simd.set_index("datazone")["dz_population"]
    vacc = (
        vaccination[["datazone", "vaccination_date", "dz_total_vaccinated"]]
        .sort_values(["datazone", "vaccination_date"])
        .copy()
    )
    vacc["dz_cum_vaccinated"] = vacc.groupby("datazone")["dz_total_vaccinated"].cumsum()
    vacc["dz_cum_prop_vaccinated"] = vacc["dz_cum_vaccinated"] / vacc["datazone"].map(
        pop
    )
    return vacc[
        ["datazone", "vaccination_date", "dz_cum_vaccinated", "dz_cum_prop_vaccinated"]
    ]


def build_dz_rolling_positivity(testing: pd.DataFrame) -> pd.DataFrame:
    """7-day rolling test positivity per (datazone, collection_date).

    For each (datazone, date) row in the testing aggregate, computes the sum of
    positive and total tests over the preceding 7 rows (calendar gaps are not
    filled — this is a row-based window over whatever dates are present in the
    testing data).  A datazone-date pair that has fewer than 7 prior rows uses
    all available history (``min_periods=1``).  Division by zero total tests
    yields NaN.

    Returns columns: datazone, collection_date, dz_7d_test_positivity.
    Intended for a direct left-merge on (collection_date, datazone) into the
    main dataset, so no merge_asof is needed.
    """
    df = (
        testing[["datazone", "collection_date", "dz_positive_tests", "dz_total_tests"]]
        .sort_values(["datazone", "collection_date"])
        .copy()
    )
    grp = df.groupby("datazone")
    pos_7d = grp["dz_positive_tests"].transform(
        lambda x: x.rolling(7, min_periods=1).sum()
    )
    tot_7d = grp["dz_total_tests"].transform(
        lambda x: x.rolling(7, min_periods=1).sum()
    )
    df["dz_7d_test_positivity"] = pos_7d / tot_7d.replace(0.0, float("nan"))
    return df[["datazone", "collection_date", "dz_7d_test_positivity"]]


def build_cluster_descriptors(ds: pd.DataFrame) -> pd.DataFrame:
    """Compute per-cluster summary statistics.

    Each ``cluster_id`` already encodes the window, lineage, and resolution
    (e.g. ``W042|BA.2|R0.3|C017``), so grouping on it yields counts that are
    resolution-specific.

    cluster_size
        Number of unique sequences assigned to this cluster.
    cluster_n_datazones
        Number of distinct datazones represented in the cluster — a coarse
        geographic spread metric.
    cluster_start_date / cluster_end_date
        Earliest and latest collection dates of sequences in the cluster.
    cluster_duration_days
        cluster_end_date − cluster_start_date in whole days.  Zero for clusters
        where all sequences were collected on the same day (including singletons).

    The returned frame is merged back into ``ds`` on ``cluster_id``.
    """
    agg = (
        ds.groupby("cluster_id", sort=False)
        .agg(
            cluster_size=("sequence_id", "nunique"),
            cluster_n_datazones=("datazone", "nunique"),
            cluster_start_date=("collection_date", "min"),
            cluster_end_date=("collection_date", "max"),
        )
        .reset_index()
    )
    agg["cluster_duration_days"] = (
        agg["cluster_end_date"] - agg["cluster_start_date"]
    ).dt.days
    return agg[
        [
            "cluster_id",
            "cluster_size",
            "cluster_n_datazones",
            "cluster_start_date",
            "cluster_end_date",
            "cluster_duration_days",
        ]
    ]


def main() -> int:
    """Run consolidation and write the final clustering analysis dataset."""
    ap = argparse.ArgumentParser(
        description="Consolidate cluster parquets into analysis dataset."
    )
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()

    setup_logging(args.log_level)
    cfg = load_config(args.root / args.config)
    pipe = cfg["pipeline"]
    required_nextclade_qc = cfg["tn93"]["nextclade_qc"]
    proc = {k: args.root / v for k, v in cfg["data"]["processed"].items()}

    resolutions = [float(x) for x in pipe["leiden_resolutions"]]
    window_size = timedelta(weeks=pipe["window_size_weeks"])
    step = timedelta(weeks=pipe["step_weeks"])

    # ── Load all processed parquets ───────────────────────────────────────────
    metadata = pd.read_parquet(proc["metadata"])
    metadata = filter_nextclade_qc(metadata, required_nextclade_qc)
    metadata["collection_date"] = pd.to_datetime(metadata["collection_date"])

    testing = pd.read_parquet(proc["testing"])
    testing["collection_date"] = pd.to_datetime(testing["collection_date"])

    vaccination = pd.read_parquet(proc["vaccination"])
    vaccination["vaccination_date"] = pd.to_datetime(vaccination["vaccination_date"])

    simd = pd.read_parquet(proc["simd"])

    hb_trends = pd.read_parquet(proc["hb_trends"])
    hb_trends["date"] = pd.to_datetime(hb_trends["date"])

    windows = build_windows(
        metadata["collection_date"].min(),
        metadata["collection_date"].max(),
        window_size,
        step,
    )
    logging.info("%d time windows", len(windows))

    # Guard: build_singletons below adds assignments only for groups of size == 1,
    # which are skipped by 02_gen_tn93_commands.py when min_group_size >= 2.
    # If min_group_size < 2, single-sequence groups also pass through cluster_pairwise.py
    # (as fallback isolates), creating duplicate rows in the final dataset.
    min_group_size = pipe.get("min_group_size", 2)
    if min_group_size < 2:
        logging.warning(
            "min_group_size=%d: single-sequence lineage groups may appear in both "
            "cluster parquets (via fallback_isolates) and the singletons built here. "
            "Check for duplicate sequence_id \u00d7 window_id \u00d7 resolution rows.",
            min_group_size,
        )

    # ── Cluster assignments ───────────────────────────────────────────────────
    clustering = load_cluster_parquets(proc["cluster_long_dir"])
    allowed_sequence_ids = set(metadata["sequence_id"].astype(str))
    unexpected_cluster_ids = ~clustering["sequence_id"].astype(str).isin(
        allowed_sequence_ids
    )
    if unexpected_cluster_ids.any():
        n_unexpected = clustering.loc[unexpected_cluster_ids, "sequence_id"].nunique()
        raise ValueError(
            f"Cluster assignments contain {n_unexpected} sequence(s) outside the "
            f"configured Nextclade QC cohort {str(required_nextclade_qc)!r}. "
            "Remove stale group, pairwise, and cluster outputs, then rebuild from "
            "method/02_gen_tn93_commands.py."
        )
    singletons = build_singletons(metadata, windows, resolutions)
    logging.info(
        "Cluster rows: %d + %d singletons = %d total",
        len(clustering),
        len(singletons),
        len(clustering) + len(singletons),
    )
    clustering = pd.concat([clustering, singletons], ignore_index=True)

    window_info = build_window_info(metadata, testing, windows)

    # ── Merge metadata + SIMD ─────────────────────────────────────────────────
    # Inner join with SIMD; log any sequences silently dropped from datazones absent
    # in the SIMD table so the data loss is always visible in the run log.
    n_meta_before = len(metadata)
    full_meta = metadata.merge(simd, on="datazone", how="inner")
    n_dropped = n_meta_before - len(full_meta)
    if n_dropped:
        logging.warning(
            "SIMD inner join dropped %d/%d sequences (%.1f%%) from datazones "
            "absent in SIMD.",
            n_dropped,
            n_meta_before,
            100.0 * n_dropped / n_meta_before,
        )

    # Retain the per-patient most-recent-prior-vaccination date in the final dataset
    # (renamed from vaccination_date to avoid collision with the daily aggregate join
    # that follows). This preserves the individual-level lookback computed in step 01.
    full_meta = full_meta.rename(columns={"vaccination_date": "vacc_date_prior"})

    # ── Main wide join ────────────────────────────────────────────────────────
    ds = (
        clustering.merge(full_meta, on="sequence_id", how="inner")
        .merge(window_info, on="window_id", how="inner")
        .merge(testing, on=["collection_date", "datazone"], how="left")
        .merge(
            # Rename vaccination_date -> collection_date so the daily vaccination
            # aggregate joins on sample collection date, adding dz_total_vaccinated
            # and the mean/median age/dose columns for that specific day.
            vaccination.rename(columns={"vaccination_date": "collection_date"}),
            on=["collection_date", "datazone"],
            how="left",
        )
    )

    # Window index (integer) — extracted before the window-bounds merge.
    ds["window_idx"] = ds["window_id"].str.extract(r"(\d+)", expand=False).astype(int)

    # Use *theoretical* window boundaries rather than the empirical min/max of sample
    # dates within each window; the latter produces narrower apparent ranges when a
    # window has sparse data and is confusing to interpret.
    win_bounds = pd.DataFrame(
        [
            {
                "window_id": f"W{str(i).zfill(3)}",
                "wn_start_date": start,
                "wn_end_date": end,
                "wn_mid_date": start + (end - start) / 2,
            }
            for i, (start, end) in enumerate(windows, start=1)
        ]
    )
    ds = ds.merge(win_bounds, on="window_id", how="left")

    # ── Cumulative vaccination coverage (merge_asof, backward-looking) ────────
    # dz_cum_vaccinated: total vaccination events in this datazone on or before
    #   this sequence's collection_date.
    # dz_cum_prop_vaccinated: dz_cum_vaccinated / dz_population.
    # merge_asof requires both frames sorted by the join key (collection_date).
    vacc_cum = build_dz_cumulative_vaccination(vaccination, simd)
    vacc_cum_for_merge = vacc_cum.rename(
        columns={"vaccination_date": "collection_date"}
    ).sort_values("collection_date")
    ds = ds.sort_values("collection_date").reset_index(drop=True)
    ds = pd.merge_asof(
        ds,
        vacc_cum_for_merge[
            [
                "datazone",
                "collection_date",
                "dz_cum_vaccinated",
                "dz_cum_prop_vaccinated",
            ]
        ],
        on="collection_date",
        by="datazone",
        direction="backward",
    )

    # ── Cumulative sequencing fraction per datazone ───────────────────────────
    # dz_cum_prop_sequenced: fraction of all positive tests ever recorded in this
    # datazone that have a linked genome sequence, as of each sequence's collection
    # date. Provides a local measure of genomic surveillance intensity over time.
    cum_seq = build_dz_cumulative_sequencing(metadata, testing)
    cum_seq_for_merge = cum_seq.sort_values("collection_date")
    ds = pd.merge_asof(
        ds,  # already sorted by collection_date
        cum_seq_for_merge,
        on="collection_date",
        by="datazone",
        direction="backward",
    )

    # ── 7-day rolling test positivity per datazone ────────────────────────────
    # Joined via a standard left-merge on (collection_date, datazone) since the
    # rolling stats are pre-computed on the testing aggregate's exact date grid.
    roll_pos = build_dz_rolling_positivity(testing)
    ds = ds.merge(roll_pos, on=["collection_date", "datazone"], how="left")

    # ── Cluster descriptors ───────────────────────────────────────────────────
    # Computed now (after all joins have fixed the sequence x resolution rows)
    # and merged back on cluster_id.
    cluster_desc = build_cluster_descriptors(ds)
    ds = ds.merge(cluster_desc, on="cluster_id", how="left")

    # ── Derived individual-level variables ────────────────────────────────────
    # days_since_vaccination: number of days between the most-recent prior dose
    # and the sample collection date. NaN for unvaccinated cases (vacc_date_prior
    # is NaT when is_vaccinated == 0).
    ds["days_since_vaccination"] = (
        pd.to_datetime(ds["collection_date"]) - pd.to_datetime(ds["vacc_date_prior"])
    ).dt.days

    # ── Derived area-level variables ──────────────────────────────────────────
    # dz_test_positivity: same-day positivity rate in the datazone. Uses the
    # dz_positive_tests and dz_total_tests columns already joined from testing.
    ds["dz_test_positivity"] = ds["dz_positive_tests"] / ds["dz_total_tests"].replace(
        0.0, float("nan")
    )

    # dz_cum_incidence_per_capita: cumulative positive tests per head of
    # population, a measure of cumulative epidemic intensity in the datazone.
    ds["dz_cum_incidence_per_capita"] = (
        ds["dz_cum_positive_tests"] / ds["dz_population"]
    )

    # dz_population_density: residents per km², combining dz_population (SIMD)
    # and dz_area_km2 (shapefile StdAreaKm2 via prep_geography).  NaN where
    # the shapefile lacked an area for the datazone.
    if "dz_area_km2" in ds.columns:
        ds["dz_population_density"] = ds["dz_population"] / ds["dz_area_km2"]

    # ── Health-board daily trends join (merge_asof, backward-looking) ─────────
    # Joins the health-board daily surveillance aggregate (hospital admissions,
    # ICU occupancy, reinfection rate, etc.) to each sequence using the most
    # recent HB report on or before the sequence's collection date.
    #
    # Requires dz_health_board_code in ds (from SIMD HBcode column, populated by
    # prep_simd in 01_prep_metadata.py if present in the SIMD release). Skipped
    # with a warning if the column is absent.
    if "dz_health_board_code" in ds.columns:
        hb_cols = [c for c in hb_trends.columns if c not in ("date", "hb_code")]
        hb_for_merge = hb_trends.rename(
            columns={"date": "collection_date", "hb_code": "dz_health_board_code"}
        )[["collection_date", "dz_health_board_code"] + hb_cols].sort_values(
            "collection_date"
        )
        ds = ds.sort_values("collection_date").reset_index(drop=True)
        ds = pd.merge_asof(
            ds,
            hb_for_merge,
            on="collection_date",
            by="dz_health_board_code",
            direction="backward",
        )
        logging.info("Joined HB trends (%d columns).", len(hb_cols))
    else:
        logging.warning(
            "dz_health_board_code not found in dataset (SIMD release may lack HBcode); "
            "skipping health-board trends join."
        )

    # ── Anonymise patient_id ──────────────────────────────────────────────────
    uid = ds["patient_id"].unique()
    id_map = {old: f"P{i:06d}" for i, old in enumerate(uid, start=1)}
    ds["patient_id"] = ds["patient_id"].map(id_map)

    # ── Column order ──────────────────────────────────────────────────────────
    column_order = [
        # Window-level identifiers and summary stats
        "window_idx",
        "window_id",
        "wn_start_date",
        "wn_mid_date",
        "wn_end_date",
        "wn_no_sequences",
        "wn_positive_tests",
        "wn_prop_sequenced",
        # Sequence / cluster identifiers
        "sequence_id",
        "patient_id",
        "resolution",
        "cluster_id",
        # Cluster descriptors
        "cluster_size",
        "cluster_n_datazones",
        "cluster_start_date",
        "cluster_end_date",
        "cluster_duration_days",
        # Sample-level fields
        "collection_date",
        "datazone",
        "dz_xcoord",
        "dz_ycoord",
        "sex",
        "is_female",
        "age_band",
        "age_group",
        "age_midpoint",
        "is_vaccinated",
        "vacc_dose_number",
        "vacc_date_prior",
        "vacc_product_name",
        "vacc_booster",
        "days_since_vaccination",
        "test_type",
        "test_reason_raw",
        "test_reason",
        "s_gene_status",
        "policy_period",
        "policy_period_label",
        "policy_era",
        "is_reinfection",
        "pango_lineage",
        "clade",
        "who_voc",
        "nextclade_qc",
        # Datazone sociodemographic attributes
        "dz_population",
        "dz_working_age_population",
        "dz_area_km2",
        "dz_population_density",
        "dz_simd_rank",
        "dz_simd_quintile",
        "dz_simd_decile",
        "dz_simd_vigintile",
        "dz_simd_income_rank",
        "dz_simd_employment_rank",
        "dz_simd_education_rank",
        "dz_simd_health_rank",
        "dz_simd_access_rank",
        "dz_simd_crime_rank",
        "dz_simd_housing_rank",
        "dz_urban_rural_class",
        "dz_local_authority",
        "dz_local_authority_code",
        "dz_health_board",
        "dz_health_board_code",
        # Datazone daily testing counts (on collection_date)
        "dz_total_tests",
        "dz_positive_tests",
        "dz_negative_tests",
        "dz_pcr_positive_tests",
        "dz_lfd_positive_tests",
        "dz_care_home_tests",
        "dz_test_positivity",
        "dz_7d_test_positivity",
        # Datazone vaccination: daily new + cumulative coverage
        "dz_total_vaccinated",
        "dz_cum_vaccinated",
        "dz_cum_prop_vaccinated",
        # Datazone cumulative surveillance
        "dz_cum_sequences",
        "dz_cum_positive_tests",
        "dz_cum_prop_sequenced",
        "dz_cum_incidence_per_capita",
        # Health-board daily trends (most-recent report on or before collection_date)
        "hb_daily_positive",
        "hb_cumulative_positive",
        "hb_hospital_admissions",
        "hb_hospital_occupancy",
        "hb_icu_admissions",
        "hb_icu_occupancy_lt28d",
        "hb_icu_occupancy_ge28d",
        "hb_daily_reinfections",
        "hb_reinfection_rate",
    ]
    required_columns = {
        "age_group",
        "test_type",
        "test_reason_raw",
        "test_reason",
        "policy_period",
        "policy_period_label",
        "policy_era",
    }
    missing_required = sorted(required_columns - set(ds.columns))
    if missing_required:
        raise KeyError(
            "Dataset is missing expected derived metadata columns after consolidation: "
            f"{missing_required}. Rebuild processed metadata with method/01_prep_metadata.py."
        )

    final_qc_matches = (
        ds["nextclade_qc"]
        .astype("string")
        .str.strip()
        .str.casefold()
        .eq(str(required_nextclade_qc).strip().casefold())
    )
    if not final_qc_matches.all():
        raise AssertionError(
            "Final analysis rows fall outside the configured Nextclade QC cohort."
        )
    ds = ds[[c for c in column_order if c in ds.columns]].reset_index(drop=True)

    out_path: Path = proc["analysis_dataset"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_parquet(out_path, index=False, compression="zstd")
    
    logging.info(
        "Analysis dataset: %d rows \u00d7 %d cols \u2192 %s", len(ds), len(ds.columns), out_path
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
