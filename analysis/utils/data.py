"""Data-loading helpers for the clustering manuscripts — Polars edition.

Drop-in replacement for the pandas version. Public signatures are identical;
return types are ``pl.DataFrame`` / ``pl.Series`` throughout.

Key behavioural notes
---------------------
* ``assign_epoch`` returns a ``pl.Series`` with dtype ``pl.Enum(cats)`` —
  the Polars equivalent of pandas' ordered ``Categorical``.  Downstream
  code that relied on ``pd.Categorical`` ordering should use
  ``pl.Enum``-aware comparisons or cast to ``pl.Utf8`` first.
* Week truncation in ``_weekly_dominant_voc`` uses Monday-start ISO weeks
  (``dt.truncate("1w")``), a one-day shift from the pandas ``W-SUN``
  convention.  Epoch boundaries derived from the data are unaffected in
  practice.
* ``lru_cache`` arguments must remain hashable; the ``qc`` parameter stays
  a ``tuple[str, ...]`` for this reason.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
import ruptures as rpt
import polars as pl
import yaml

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

#: Primary Leiden resolution used in headline results.
PRIMARY_RESOLUTION: float = 0.3

WAVE = None


@dataclass(frozen=True)
class DOMAINS:
    overall:    str = "dz_simd_rank"
    income:     str = "dz_simd_income_rank"
    employment: str = "dz_simd_employment_rank"
    education:  str = "dz_simd_education_rank"
    health:     str = "dz_simd_health_rank"
    access:     str = "dz_simd_access_rank"
    crime:      str = "dz_simd_crime_rank"
    housing:    str = "dz_simd_housing_rank"


def _simd_domain() -> dict[str, str]:
    """Return a mapping of SIMD domain names to parquet column names."""
    return {name: col for name, col in asdict(DOMAINS()).items()}


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


def repo_root(start: Path | None = None) -> Path:
    """Walk up from *start* (default: this file) until ``config.yaml`` is found."""
    p = (start or Path(__file__)).resolve()
    for cand in [p, *p.parents]:
        if (cand / "config.yaml").exists():
            return cand
    raise FileNotFoundError("Could not locate config.yaml in any parent directory.")


@dataclass(frozen=True)
class Paths:
    root: Path
    analysis_dataset: Path

    @classmethod
    def from_config(cls, root: Path | None = None) -> "Paths":
        root = root or repo_root()
        with open(root / "config.yaml") as f:
            cfg = yaml.safe_load(f)
        proc = cfg["data"]["processed"]
        return cls(
            root=root,
            analysis_dataset=root / proc["analysis_dataset"],
        )


# ---------------------------------------------------------------------------
# Column-level loaders
# ---------------------------------------------------------------------------


def load_analysis_columns(
    columns: Iterable[str],
    *,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: Iterable[str] | None = ("good",),
    paths: Paths | None = None,
) -> pl.DataFrame:
    """Read a narrow slice of the master sequence-level parquet.

    Parameters
    ----------
    columns:
        Names of columns to read; ``resolution`` and ``nextclade_qc`` are
        added automatically when filtering is requested.
    resolution:
        If provided, rows are restricted to that Leiden resolution.
    qc:
        If provided, rows are restricted to these Nextclade QC statuses.
    paths:
        If provided, use these paths instead of resolving from config.
    """
    paths = paths or Paths.from_config()
    need = set(columns)
    if resolution is not None:
        need.add("resolution")
    if qc is not None:
        need.add("nextclade_qc")

    df = pl.read_parquet(paths.analysis_dataset, columns=sorted(need))

    if resolution is not None:
        df = df.filter(pl.col("resolution") == resolution)
    if qc is not None:
        df = df.filter(pl.col("nextclade_qc").is_in(list(qc)))

    return df

# ---------------------------------------------------------------------------
# Aggregation-level loaders
# ---------------------------------------------------------------------------

def _standardise(values):
    """Negate z-score so higher values = greater deprivation."""
    return -(values - values.mean()) / values.std()

@lru_cache(maxsize=4)
def load_cluster_features(
    resolution: float = PRIMARY_RESOLUTION,
    qc: tuple[str, ...] = ("good",),
) -> pl.DataFrame:
    """Return one row per (window_id, cluster_id) with size, date, lineage, and SIMD features."""

    cols = [
        "window_id", "window_idx", "cluster_id", "sequence_id", "collection_date",
        "wn_mid_date", "wn_prop_sequenced", "pango_lineage", "who_voc", "nextclade_qc",
        "dz_simd_rank", "dz_simd_quintile", "dz_simd_decile",
        "dz_simd_income_rank", "dz_simd_employment_rank", "dz_simd_education_rank",
        "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank",
        "dz_simd_housing_rank", "age_midpoint", "is_female", "is_vaccinated",
    ]

    df = load_analysis_columns(cols, resolution=resolution, qc=qc)

    df = df.with_columns(
        (pl.col("nextclade_qc") == "mediocre").cast(pl.Float64).alias("_is_mediocre"),
        (pl.col("nextclade_qc") == "bad").cast(pl.Float64).alias("_is_bad"),
    )

    # Build aggregation expressions ----------------------------------------
    agg_exprs: list[pl.Expr] = [
        pl.col("sequence_id").n_unique().alias("n_sequences"),
        pl.col("window_idx").first(),
        pl.col("wn_mid_date").first(),
        pl.col("wn_prop_sequenced").first(),
        pl.col("who_voc").first(),
        pl.col("pango_lineage").first(),
        pl.col("_is_mediocre").mean().alias("qc_frac_mediocre"),
        pl.col("_is_bad").mean().alias("qc_frac_bad"),
        # Demographics
        pl.col("age_midpoint").median().alias("median_age"),
        pl.col("age_midpoint").std().fill_null(0).alias("age_diversity"),
        pl.col("is_female").mean().alias("frac_female"),
        pl.col("is_vaccinated").mean().alias("frac_vaccinated"),
        # SIMD quintile / decile summary
        # mode() returns all modal values; sort().first() picks the smallest on ties
        pl.col("dz_simd_quintile").drop_nulls().mode().sort().first().alias("simd_quintile_mode"),
        pl.col("dz_simd_quintile").std().alias("simd_quintile_std"),
        pl.col("dz_simd_decile").drop_nulls().mode().sort().first().alias("simd_decile_mode"),
        pl.col("dz_simd_decile").std().fill_null(0).alias("simd_decile_std"),
    ]

    for dom, col in _simd_domain().items():
        agg_exprs.append(pl.col(col).mean().alias(f"simd_{dom}_mean"))

    out = df.group_by(["window_id", "cluster_id"]).agg(agg_exprs)

    out = out.with_columns(
        (pl.col("n_sequences") == 1).cast(pl.Int8).alias("is_singleton"),
    )
    for dom, col in _simd_domain().items():
        out = out.with_columns(_standardise(out[f"simd_{dom}_mean"]).alias(f"{dom}_zscore"))

    # Wave assignment --------------------------------------------------------
    _, wdf = detect_waves(df)
    global WAVE
    if WAVE is None:
        WAVE = get_voc_waves(wdf)
    out = out.with_columns(
        assign_wave(out["wn_mid_date"], wdf).alias("wave")
    )
    out = out.with_columns(pl.col("wave").fill_null("unknown").alias("wave"))
    out = out.with_columns(pl.col("who_voc").fill_null("non_voc").alias("who_voc"))

    # Pandas reports NaN std for singletons. Treat as 0 (no within-cluster mixing).
    out = out.with_columns(
        pl.when(pl.col("simd_quintile_std").is_null() & (pl.col("is_singleton") == 1))
        .then(pl.lit(0.0))
        .otherwise(pl.col("simd_quintile_std"))
        .alias("simd_quintile_std")
    )

    # Log sequencing proportion (natural log) --------------------------------
    out = out.with_columns(
        pl.col("wn_prop_sequenced").log().alias("log_seq_prop")
    )

    assert out["log_seq_prop"].is_finite().all(), (
        "Non-finite offset values — check wn_prop_sequenced > 0"
    )

    # Ordered integer columns (equivalent to pandas ordered Categorical).
    # Cast to pl.Enum if downstream code needs categorical ordering; kept as
    # Int32 here so arithmetic/statsmodels usage is unaffected.
    out = out.with_columns(
        pl.col("simd_quintile_mode").cast(pl.Int32),
        pl.col("simd_decile_mode").cast(pl.Int32),
    )

    out = out.with_columns((pl.col("n_sequences") - 1).alias("n_sequences_minus_one"))

    out = out.select([
        "window_id", "window_idx", "wn_mid_date", "wn_prop_sequenced", "log_seq_prop",
        "cluster_id", "n_sequences", "n_sequences_minus_one",
        "median_age", "age_diversity", "frac_female", "frac_vaccinated",
        "simd_decile_mode", "simd_decile_std", "simd_quintile_mode", "simd_quintile_std",
        *list(f"{dom}_zscore" for dom in _simd_domain().keys()),
        "pango_lineage", "who_voc", "qc_frac_mediocre", "qc_frac_bad", "wave",
    ])

    return out


@lru_cache(maxsize=4)
def load_individual_features(
    qc: tuple[str, ...] = ("good",),
    format: Literal["aggregate", "long"] = "aggregate",
) -> pl.DataFrame:
    """Return sequence-level features in aggregate or long format.

    Parameters
    ----------
    qc:
        QC filter values to pass to load_analysis_columns.
    format:
        - "aggregate" : one row per sequence_id, non_singleton_k/n summarised
                        across all (window_id, resolution) combinations.
        - "long"      : one row per (sequence_id, resolution, window_id), with binary
                        in_non_singleton for use in mixed-effects models.
    """
    _NON_SINGLETON_RE = r"\|C\d+$"
    cols = [
        "window_id", "window_idx", "cluster_id", "patient_id", "sequence_id", "resolution",
        "wn_mid_date", "wn_prop_sequenced", "pango_lineage", "who_voc", "nextclade_qc",
        "datazone", "dz_simd_rank", "dz_simd_quintile", "dz_simd_decile",
        "dz_simd_income_rank", "dz_simd_employment_rank", "dz_simd_education_rank",
        "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank",
        "dz_simd_housing_rank", "age_band", "age_midpoint", "is_female", "is_vaccinated", "collection_date",
    ]
    df = load_analysis_columns(cols, resolution=None, qc=qc)
    df = df.with_columns(
        pl.col("cluster_id").cast(pl.Utf8).str.contains(_NON_SINGLETON_RE)
        .cast(pl.UInt8).alias("_non_singleton")
    )

    # ── shared: stable per-sequence columns ─────────────────────────────────
    stable = [
        "patient_id", "collection_date",
        "age_band", "age_midpoint", "is_female", "is_vaccinated", "who_voc",
        "datazone", "dz_simd_quintile", "dz_simd_decile",
        *list(_simd_domain().values()), "pango_lineage",
    ]

    # ── wave detection (needed in both paths) ────────────────────────────────
    _, wdf = detect_waves(df)
    global WAVE
    if WAVE is None:
        WAVE = get_voc_waves(wdf)

    # ════════════════════════════════════════════════════════════════════════
    if format == "long":
        return _build_long(df, stable, wdf)
    else:
        return _build_aggregate(df, stable, wdf)


def _attach_stable_and_wave(
    out: pl.DataFrame,
    wdf: pl.DataFrame,
) -> pl.DataFrame:
    """Add z-scores, wave, and voc fill-null — shared by both paths."""
    for dom in _simd_domain().keys():
        col = "dz_simd_rank" if dom == "overall" else f"dz_simd_{dom}_rank"
        out = out.with_columns(_standardise(out[col]).alias(f"{dom}_zscore"))

    out = out.with_columns(
        assign_wave(out["collection_date"], wdf).alias("wave")
    )
    out = out.with_columns(
        pl.col("wave").fill_null("unknown"),
        pl.col("who_voc").fill_null("non_voc"),
    )
    return out


def _build_aggregate(
    df: pl.DataFrame,
    stable: list[str],
    wdf: pl.DataFrame,
) -> pl.DataFrame:
    """One row per sequence_id — original behaviour."""
    agg_exprs: list[pl.Expr] = (
        [pl.col(c).first() for c in stable]
        + [
            pl.col("_non_singleton").sum().cast(pl.Int64).alias("non_singleton_k"),
            pl.col("_non_singleton").count().alias("non_singleton_n"),
            pl.col("wn_prop_sequenced").mean(),
            pl.col("wn_prop_sequenced").log().mean().alias("log_seq_prop"),
            pl.col("window_id").n_unique().alias("n_windows"),
            pl.col("resolution").n_unique().alias("n_resolutions"),
        ]
    )
    out = df.group_by("sequence_id").agg(agg_exprs)
    out = _attach_stable_and_wave(out, wdf)

    assert out["log_seq_prop"].is_finite().all(), (
        "Non-finite offset values — check wn_prop_sequenced > 0"
    )
    return out.select([
        "patient_id", "sequence_id", "collection_date", "log_seq_prop",
        "age_band", "age_midpoint", "is_female", "is_vaccinated", "pango_lineage", "who_voc",
        "datazone", "dz_simd_quintile", "dz_simd_decile",
        *[f"{dom}_zscore" for dom in _simd_domain().keys()],
        "non_singleton_k", "non_singleton_n", "n_windows", "n_resolutions", "wave",
    ])


def _build_long(
    df: pl.DataFrame,
    stable: list[str],
    wdf: pl.DataFrame,
) -> pl.DataFrame:
    """One row per (sequence_id, resolution, window_id)."""
    agg_exprs: list[pl.Expr] = (
        [pl.col(c).first() for c in stable]
        + [
            pl.col("_non_singleton").max().alias("in_non_singleton"),
            pl.col("wn_mid_date").first().alias("wn_mid_date"),
            pl.col("window_idx").first().alias("window_idx"),
            pl.col("wn_prop_sequenced").mean(),
            pl.col("wn_prop_sequenced").log().mean().alias("log_seq_prop"),
        ]
    )
    out = df.group_by(["sequence_id", "resolution", "window_id"]).agg(agg_exprs)
    out = _attach_stable_and_wave(out, wdf)

    assert out["log_seq_prop"].is_finite().all(), (
        "Non-finite offset values — check wn_prop_sequenced > 0"
    )

    age_map = {
        '00-04': '00-09', '05-09': '00-09',
        '10-14': '10-19', '15-19': '10-19',
        '20-24': '20-39', '25-29': '20-39', '30-34': '20-39', '35-39': '20-39',
        '40-44': '40–59', '45-49': '40–59', '50-54': '40–59', '55-59': '40–59',
        '60-64': '60-74', '65-69': '60-74', '70-74': '60-74',
        '75+': 'elderly',
    }

    out = out.with_columns(
        pl.col("age_band").replace(age_map).alias("age_group")
    )
    return out.select([
        "patient_id", "sequence_id", "resolution", "window_id", "window_idx",
        "wn_mid_date", "collection_date", "log_seq_prop",
        "age_band", "age_group", "is_female", "age_midpoint", "is_vaccinated", "pango_lineage", "who_voc",
        "datazone", "dz_simd_quintile", "dz_simd_decile",
        *[f"{dom}_zscore" for dom in _simd_domain().keys()],
        "in_non_singleton", "wave",
    ])


# ---------------------------------------------------------------------------
# Derive epidemic waves
# ---------------------------------------------------------------------------
def detect_waves(
    sequences: pl.DataFrame,
    *,
    smooth_window: int = 7,
    min_peak_height: float = 100,
    min_wave_duration: int = 21,
    min_trough_drop_frac: float = 0.0,
    penalty: float = 13.0,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Detect epidemic waves from a sequences DataFrame.

    Parameters
    ----------
    sequences : pl.DataFrame
        One row per sequence. Must contain:
          - ``collection_date``  (Date)
          - ``sequence_id``      (any, used for counting)
          - ``pango_lineage``    (str)
    smooth_window : int
        Rolling mean window size for pre-smoothing.
    min_peak_height : float
        Minimum smoothed count for a wave's peak to be retained.
    min_wave_duration : int
        Minimum number of days for a wave to be retained.
    min_trough_drop_frac : float
        Trough must drop by at least this fraction below the smaller adjacent
        peak for two segments to be kept as separate waves.
    penalty : float
        PELT penalty (higher → fewer breakpoints → fewer waves).

    Returns
    -------
    labelled : pl.DataFrame
        Gap-filled daily counts with ``smoothed_count``, ``outbreak_id``,
        and ``phase`` columns.
    waves : pl.DataFrame
        One row per detected wave with summary metadata and ``dominant_lineage``
        / ``dominant_lineage_count`` columns.
    """

    # ------------------------------------------------------------------ #
    # 0.  Validate inputs                                                 #
    # ------------------------------------------------------------------ #
    required = {"collection_date", "sequence_id", "pango_lineage"}
    missing  = required - set(sequences.columns)
    if missing:
        raise ValueError(f"sequences is missing required columns: {missing}")

    # ------------------------------------------------------------------ #
    # 1.  Compute daily counts and gap-fill                               #
    # ------------------------------------------------------------------ #

    sequences = sequences.with_columns(
        pl.col("collection_date").cast(pl.Date)
    )
    daily_counts = (
        sequences
        .group_by("collection_date")
        .agg(pl.col("sequence_id").n_unique().alias("count"))
        .sort("collection_date")
    )

    all_dates = pl.DataFrame({
        "collection_date": pl.date_range(
            daily_counts["collection_date"].min(),
            daily_counts["collection_date"].max(),
            interval="1d",
            eager=True,
        )
    })

    df_full = (
        all_dates
        .join(daily_counts, on="collection_date", how="left")
        .with_columns(pl.col("count").fill_null(0))
    )

    # ------------------------------------------------------------------ #
    # 2.  Centred rolling mean                                            #
    # ------------------------------------------------------------------ #
    half_window = smooth_window // 2
    df_full = (
        df_full
        .with_columns(
            pl.col("count")
            .cast(pl.Float64)
            .rolling_mean(window_size=smooth_window, min_samples=1)
            .shift(-half_window)
            .alias("smoothed_count")
        )
        .with_columns(pl.col("smoothed_count").forward_fill().backward_fill())
    )

    dates    = df_full["collection_date"].to_list()
    smoothed = df_full["smoothed_count"].to_numpy()
    n        = len(smoothed)

    # ------------------------------------------------------------------ #
    # 3.  PELT change-point detection                                     #
    # ------------------------------------------------------------------ #
    algo     = rpt.Pelt(model="rbf", min_size=max(2, min_wave_duration // 2)).fit(smoothed.reshape(-1, 1))
    raw_bkps = algo.predict(pen=penalty)
    bkps     = [0] + raw_bkps[:-1] + [n]
    segments = list(zip(bkps[:-1], bkps[1:]))

    # ------------------------------------------------------------------ #
    # 4.  Locate peak within each segment                                 #
    # ------------------------------------------------------------------ #
    seg_peaks = [
        (start, end - 1, int(np.argmax(smoothed[start:end])) + start)
        for start, end in segments
    ]

    # ------------------------------------------------------------------ #
    # 5.  Merge segments whose trough is too shallow                      #
    # ------------------------------------------------------------------ #
    def _should_merge(left: tuple, right: tuple) -> bool:
        _, l_end, l_peak = left
        r_start, _, r_peak = right
        trough = smoothed[l_end: r_start + 1]
        if len(trough) == 0:
            return True
        smaller_peak = min(smoothed[l_peak], smoothed[r_peak])
        if smaller_peak <= 0:
            return True
        return (smaller_peak - trough.min()) / smaller_peak < min_trough_drop_frac

    merged = [seg_peaks[0]]
    for seg in seg_peaks[1:]:
        if _should_merge(merged[-1], seg):
            prev_start, _, prev_peak = merged[-1]
            _, seg_end, seg_peak     = seg
            new_peak   = prev_peak if smoothed[prev_peak] >= smoothed[seg_peak] else seg_peak
            merged[-1] = (prev_start, seg_end, new_peak)
        else:
            merged.append(seg)

    # ------------------------------------------------------------------ #
    # 6.  Filter by height and duration                                   #
    # ------------------------------------------------------------------ #
    cleaned = [
        (s, e, p) for s, e, p in merged
        if (e - s + 1) >= min_wave_duration
        and smoothed[p] >= min_peak_height
        and s <= p <= e
    ]

    # ------------------------------------------------------------------ #
    # 7.  Build output                                                    #
    # ------------------------------------------------------------------ #
    _empty_waves = pl.DataFrame({
        "outbreak_id":         pl.Series([], dtype=pl.Int64),
        "start_date":          pl.Series([], dtype=pl.Date),
        "peak_date":           pl.Series([], dtype=pl.Date),
        "end_date":            pl.Series([], dtype=pl.Date),
        "duration_days":       pl.Series([], dtype=pl.Int64),
        "peak_smoothed_count": pl.Series([], dtype=pl.Float64),
    })

    if not cleaned:
        return (
            df_full.with_columns([
                pl.lit(None).cast(pl.Int64).alias("outbreak_id"),
                pl.lit(None).cast(pl.String).alias("phase"),
            ]),
            _empty_waves,
        )

    outbreak_id = [None] * n
    phase       = [None] * n
    records     = []

    for num, (start_idx, end_idx, peak_idx) in enumerate(cleaned, start=1):
        for i in range(start_idx, end_idx + 1):
            outbreak_id[i] = num
            phase[i]       = "growth" if i <= peak_idx else "decline"

        records.append({
            "outbreak_id":         num,
            "start_date":          dates[start_idx],
            "peak_date":           dates[peak_idx],
            "end_date":            dates[end_idx],
            "duration_days":       end_idx - start_idx + 1,
            "peak_smoothed_count": float(smoothed[peak_idx]),
        })

    labelled = df_full.with_columns([
        pl.Series("outbreak_id", outbreak_id),
        pl.Series("phase",       phase),
    ])
    waves = pl.DataFrame(records)

    # ------------------------------------------------------------------ #
    # 8.  Dominant pango lineage per wave                                #
    # ------------------------------------------------------------------ #
    dominant = (
        sequences
        .select(["collection_date", "pango_lineage"])
        .join(
            labelled.select(["collection_date", "outbreak_id"]),
            on="collection_date",
            how="left",
        )
        .filter(pl.col("outbreak_id").is_not_null())
        .group_by(["outbreak_id", "pango_lineage"])
        .agg(pl.len().alias("n"))
        .sort(["outbreak_id", "n"], descending=[False, True])
        .group_by("outbreak_id", maintain_order=True)
        .first()
        .rename({"pango_lineage": "dominant_lineage", "n": "dominant_lineage_count"})
    )

    waves = waves.join(dominant, on="outbreak_id", how="left")

    return labelled, waves

def get_voc_waves(waves) -> list[tuple[str, pl.Date, pl.Date]]:
    wave_labels = []
    for i, row in enumerate(waves.iter_rows(named=True)):
        wave_labels.append((
            f"WV{row['outbreak_id']}_{row["dominant_lineage"]}_C{row['dominant_lineage_count']}",
            row["start_date"], row["end_date"],
        ))
    return wave_labels


def assign_wave(
    dates: pl.Series,
    waves: pl.DataFrame
) -> pl.Series:
    """Assign each row to a wave label based on its date.

    Returns a ``pl.Series`` with dtype ``pl.Enum(cats)`` — the ordered
    categorical equivalent of ``pd.Categorical(..., ordered=True)``.
    Null values indicate rows that fall outside all wave windows.
    """
    wv = get_voc_waves(waves)
    cats = [lbl for lbl, *_ in wv]
    wave_enum = pl.Enum(cats)

    # Build label column via successive when/otherwise passes
    label_expr: pl.Expr = pl.lit(None, dtype=pl.Utf8)
    for lbl, s, e in wv:
        label_expr = (
            pl.when(
                (pl.col("_date") >= s) & (pl.col("_date") <= e)
            )
            .then(pl.lit(lbl))
            .otherwise(label_expr)
        )

    result = (
        pl.DataFrame({"_date": dates.cast(pl.Date)})
        .with_columns(label_expr.alias("wave"))
        ["wave"]
        .cast(wave_enum)
    )
    return result


def __getattr__(name: str):
    """Lazy module-level attribute so ``data.VOC_EPOCHS`` uses derived values."""
    if name == "WAVES":
        if WAVE is None:
            df  = load_analysis_columns(
                ["collection_date", "sequence_id", "pango_lineage"],
                resolution=PRIMARY_RESOLUTION
            )
            _, wdf = detect_waves(df)
            return get_voc_waves(wdf)
        return WAVE
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")