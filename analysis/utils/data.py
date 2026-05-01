"""Data-loading helpers."""

from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import date
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Literal, Sequence

import polars as pl
import geopandas as gpd
import yaml

from .policy import attach_period

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

PRIMARY_RESOLUTION: float = 0.3

# Obtained from wave_dates.py
WAVE = [
    ("pre-Alpha (B.1.177)",       "2020-08-10",  "2021-02-08"),
    ("Alpha (B.1.1.7)",           "2020-11-30",  "2021-06-14"),
    ("Delta (AY.*/B.1.617.2)",    "2021-05-03",  "2021-12-27"),
    ("Omicron (BA.1)",            "2021-12-06",  "2022-03-28"),
    ("Omicron (BA.2)",            "2022-01-31",  "2022-07-04"),
    ("Omicron (BA.4)",            "2022-05-23",  "2022-10-24"),
    ("Omicron (BA.5 / BE.*)",     "2022-05-23",  "2022-12-19"),
    ("Omicron (BQ.1)",            "2022-10-03",  "2023-02-06"),
    ("Omicron (XBB)",             "2022-12-12",  "2023-02-06"),
]


@dataclass(frozen=True)
class DOMAINS:
    overall: str = "dz_simd_rank"
    income: str = "dz_simd_income_rank"
    employment: str = "dz_simd_employment_rank"
    education: str = "dz_simd_education_rank"
    health: str = "dz_simd_health_rank"
    access: str = "dz_simd_access_rank"
    crime: str = "dz_simd_crime_rank"
    housing: str = "dz_simd_housing_rank"


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
    geography: Path

    @classmethod
    def from_config(cls, root: Path = None) -> "Paths":
        root = root or repo_root()
        with open(root / "config.yaml") as f:
            cfg = yaml.safe_load(f)
        proc = cfg["data"]["processed"]
        return cls(
            root=root,
            analysis_dataset=root / proc["analysis_dataset"],
            geography=root / proc["geography"]
        )



# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _normalise_waves(
        waves: Sequence[tuple[str, date | str, date | str]]
) -> list[tuple[str, date, date]]:
    """Coerce wave tuples into date-backed tuples."""
    normalised: list[tuple[str, date, date]] = []
    for label, start, end in waves:
        start_date = date.fromisoformat(start) if isinstance(start, str) else start
        end_date = date.fromisoformat(end) if isinstance(end, str) else end
        if end_date < start_date:
            raise ValueError(f"Wave {label!r} ends before it starts.")
        normalised.append((label, start_date, end_date))
    return normalised


WAVES = _normalise_waves(WAVE)


def _assign_wave(dates: pl.Series) -> pl.Series:
    """Assign each row to a configured epidemic wave label."""
    if not WAVES:
        return pl.Series("wave", ["unknown"] * len(dates), dtype=pl.String)

    cats = [label for label, *_ in WAVES] + ["unknown"]
    wave_enum = pl.Enum(cats)

    label_expr: pl.Expr = pl.lit("unknown", dtype=pl.Utf8)
    for label, start, end in WAVES:
        label_expr = (
            pl.when((pl.col("_date") >= start) & (pl.col("_date") <= end))
            .then(pl.lit(label))
            .otherwise(label_expr)
        )

    return (
        pl.DataFrame({"_date": dates.cast(pl.Date)})
        .with_columns(label_expr.alias("wave"))
        .get_column("wave")
        .cast(wave_enum)
    )


def _with_wave(df: pl.DataFrame, date_col: str) -> pl.DataFrame:
    """Attach configured epidemic wave labels using the given date column."""
    return df.with_columns(_assign_wave(df[date_col]).alias("wave"))


def _with_policy(df: pl.DataFrame, date_col: str) -> pl.DataFrame:
    """Attach configured policy-period labels using the given date column."""
    return attach_period(df, date_col)


QCStatus = Literal["good", "mediocre", "bad"]
_VALID_QC: frozenset[str] = frozenset({"good", "mediocre", "bad"})


def _validate_qc(qc: Iterable[QCStatus] | QCStatus) -> None:
    if qc is not None:
        if isinstance(qc, str):
            qc: Iterable[QCStatus] = tuple((qc,))
        else:
            qc: Iterable[QCStatus] = tuple(qc)
        invalid = set(qc) - _VALID_QC
        if invalid:
            raise ValueError(
                f"Invalid QC status values: {invalid!r}. "
                f"Must be one or more of {sorted(_VALID_QC)}"
            )


def _standardise(values):
    """Negate z-score so higher values = greater deprivation."""
    return -(values - values.mean()) / values.std()


def _shannon_entropy(counts):
    """Shannon entropy (bits) from a series of counts."""
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts / total
    # mask zeros to avoid log(0)
    p = p.filter(p > 0)
    return -(p * p.log(base=2)).sum()

# ---------------------------------------------------------------------------
# Column-level loaders
# ---------------------------------------------------------------------------

def load_analysis_columns(
    columns: Iterable[str],
    *,
    all_cols: bool = False,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: Iterable[QCStatus] | QCStatus = "good"
) -> pl.DataFrame:
    """Read a narrow slice of the master sequence-level parquet.

    Parameters
    ----------
    columns:
        Names of columns to read; ``resolution`` and ``nextclade_qc`` are
        added automatically when filtering is requested.
        ``sequence_id``, ``collection_date``, ``wave``,
        ``policy_period``, ``policy_period_label``,  and ``policy_intensity``
        are also added automatically.
    all_cols:
        If True, ignore *columns* and read all columns.  This is not recommended
        for general use, but can be useful for exploratory analysis or when
        debugging.
    resolution:
        If provided, rows are restricted to that Leiden resolution.
    qc:
        If provided, rows are restricted to these Nextclade QC statuses.
        Accepted values: ``"good"``, ``"mediocre"``, ``"bad"``.
        Pass ``None`` to skip QC filtering entirely.

    Notes
    -----
    ``Available columns``:
    'window_idx', 'window_id', 'wn_start_date', 'wn_mid_date',
    'wn_end_date', 'wn_no_sequences', 'wn_positive_tests',
    'wn_prop_sequenced', 'sequence_id', 'patient_id', 'resolution',
    'cluster_id', 'cluster_size', 'cluster_n_datazones',
    'cluster_start_date', 'cluster_end_date', 'cluster_duration_days',
    'collection_date', 'datazone', 'dz_xcoord', 'dz_ycoord', 'sex',
    'is_female', 'age_band', 'age_midpoint', 'is_vaccinated',
    'vacc_dose_number', 'vacc_date_prior', 'vacc_product_name',
    'vacc_booster', 'days_since_vaccination', 's_gene_status',
    'is_reinfection', 'pango_lineage', 'clade', 'who_voc', 'nextclade_qc',
    'dz_population', 'dz_working_age_population', 'dz_area_km2',
    'dz_population_density', 'dz_simd_rank', 'dz_simd_quintile',
    'dz_simd_decile', 'dz_simd_vigintile', 'dz_simd_income_rank',
    'dz_simd_employment_rank', 'dz_simd_education_rank',
    'dz_simd_health_rank', 'dz_simd_access_rank', 'dz_simd_crime_rank',
    'dz_simd_housing_rank', 'dz_urban_rural_class', 'dz_local_authority',
    'dz_local_authority_code', 'dz_health_board', 'dz_health_board_code',
    'dz_total_tests', 'dz_positive_tests', 'dz_negative_tests',
    'dz_pcr_positive_tests', 'dz_lfd_positive_tests', 'dz_care_home_tests',
    'dz_test_positivity', 'dz_7d_test_positivity', 'dz_total_vaccinated',
    'dz_cum_vaccinated', 'dz_cum_prop_vaccinated', 'dz_cum_sequences',
    'dz_cum_positive_tests', 'dz_cum_prop_sequenced',
    'dz_cum_incidence_per_capita', 'hb_daily_positive',
    'hb_cumulative_positive', 'hb_hospital_admissions',
    'hb_hospital_occupancy', 'hb_icu_admissions', 'hb_icu_occupancy_lt28d',
    'hb_icu_occupancy_ge28d', 'hb_daily_reinfections',
    'hb_reinfection_rate'

    See ``data/processed/analysis_dataset_description.md`` for details.

    Raises
    ------
    ValueError
        If any value in ``qc`` is not one of the accepted QC statuses.
    """
    _validate_qc(qc)

    paths = Paths.from_config()

    if all_cols:
        return pl.read_parquet(paths.analysis_dataset)

    need = {"sequence_id", "collection_date"}
    need = need.union(columns)

    if resolution is not None:
        need.add("resolution")
    if qc is not None:
        need.add("nextclade_qc")

    df = pl.read_parquet(paths.analysis_dataset, columns=list(need))

    if resolution is not None:
        df = df.filter(pl.col("resolution") == resolution)

    if qc is not None:
        if isinstance(qc, str):
            qc: Iterable[QCStatus] = tuple((qc,))
        df = df.filter(pl.col("nextclade_qc").is_in(list(qc)))

    df = _with_wave(df, "collection_date")
    df = _with_policy(df, "collection_date")

    return df


def load_datazone_info(columns: Iterable[str]) -> gpd.GeoDataFrame:
    """Read a narrow slice of the datazone information parquet.

    Parameters
    ----------
    columns:
        Names of columns to read; ``datazone`` and ``geometry`` columns are
        added automatically.

    Notes
    -----
    ``Available columns``:
    'geometry', 'dz_xcoord', 'dz_ycoord', 'dz_area_km2', 'dz_population',
    'dz_working_age_population', 'dz_urban_rural_class',
    'dz_local_authority_code', 'dz_local_authority', 'dz_health_board_code',
    'dz_health_board', 'dz_simd_rank', 'dz_simd_quintile', 'dz_simd_decile',
    'dz_simd_vigintile', 'dz_simd_income_rank', 'dz_simd_employment_rank',
    'dz_simd_education_rank', 'dz_simd_health_rank', 'dz_simd_access_rank',
    'dz_simd_crime_rank', 'dz_simd_housing_rank'
    """
    paths = Paths.from_config()
    need = {"datazone", "geometry"}
    need = need.union(columns)
    return gpd.read_parquet(paths.geography, columns=list(need))


# ---------------------------------------------------------------------------
# Aggregation-level loaders
# ---------------------------------------------------------------------------


@lru_cache(maxsize=4)
def load_cluster_features(
        qc: Iterable[QCStatus] | QCStatus = "good",
        min_cluster_size: int = 5,
) -> pl.DataFrame:
    """Return one row per (window_id, cluster_id) with size, date, lineage, and SIMD features."""
    _validate_qc(qc)

    cols = [
        "window_id", "window_idx", "cluster_id", "sequence_id", "resolution",
        "wn_mid_date", "wn_prop_sequenced", "pango_lineage", "who_voc", "nextclade_qc",
        "dz_simd_rank", "dz_simd_quintile", "dz_simd_decile", "collection_date",
        "dz_simd_income_rank", "dz_simd_employment_rank", "dz_simd_education_rank",
        "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank",
        "dz_simd_housing_rank", "age_midpoint", "age_band", "is_female", "is_vaccinated",
    ]

    df = load_analysis_columns(cols, resolution=None, qc=qc)

    # Build aggregation expressions ----------------------------------------
    agg_exprs: list[pl.Expr] = [
        pl.col("sequence_id").n_unique().alias("n_sequences"),
        pl.col("resolution").first().alias("resolution"),
        pl.col("window_id").first(),
        pl.col("window_idx").first(),
        pl.col("wn_mid_date").first(),
        pl.col("wn_prop_sequenced").first(),
        pl.col("who_voc").first(),
        pl.col("pango_lineage").first(),
        # Demographics
        pl.col("age_midpoint").median().alias("median_age"),
        pl.col("age_midpoint").mean().alias("mean_age"),
        pl.col("is_female").mean().alias("frac_female"),
        pl.col("is_vaccinated").mean().fill_null(0).alias("frac_vaccinated"),
        # SIMD quintile / decile summary
        # mode() returns all modal values; sort().first() picks the smallest on ties
        pl.col("dz_simd_quintile").drop_nulls().mode().sort().first().alias("simd_quintile_mode"),
        pl.col("dz_simd_quintile").std().fill_null(0).alias("simd_quintile_std"),
        pl.col("dz_simd_decile").drop_nulls().mode().sort().first().alias("simd_decile_mode"),
        pl.col("dz_simd_decile").std().fill_null(0).alias("simd_decile_std"),
    ]

    for dom, col in _simd_domain().items():
        agg_exprs.append(pl.col(col).mean().alias(f"simd_{dom}_mean"))

    out = df.group_by("cluster_id").agg(agg_exprs)

    # Compute entropy per cluster from age_band counts --------------------------------------
    age_entropy = (
        df.group_by(["cluster_id", "age_band"])
        .agg(pl.len().alias("count"))
        .group_by("cluster_id")
        .agg(
            pl.col("count")
            .map_batches(_shannon_entropy, return_dtype=pl.Float64, returns_scalar=True)
            .alias("age_diversity")
        )
    )

    out = out.join(age_entropy, on="cluster_id", how="left")

    for dom, col in _simd_domain().items():
        out = out.with_columns(_standardise(out[f"simd_{dom}_mean"]).alias(f"{dom}_zscore"))

    out = _with_wave(out, "wn_mid_date")
    out = _with_policy(out, "wn_mid_date")
    out = out.with_columns(pl.col("who_voc").fill_null("other").alias("who_voc"))
    out = out.with_columns((pl.col("n_sequences") - 1).alias("n_sequences_minus_one"))

    out = out.select([
        "window_id", "window_idx", "wn_mid_date", "wn_prop_sequenced",
        "cluster_id", "n_sequences", "n_sequences_minus_one", "resolution",
        "median_age", "mean_age", "age_diversity", "frac_female", "frac_vaccinated",
        "simd_decile_mode", "simd_decile_std", "simd_quintile_mode", "simd_quintile_std",
        *list(f"{dom}_zscore" for dom in _simd_domain().keys()),
        "pango_lineage", "who_voc", "wave",
        "policy_period", "policy_period_label", "policy_intensity",
    ])

    out = out.filter(pl.col("n_sequences") >= min_cluster_size)

    return out


@lru_cache(maxsize=4)
def load_individual_features(
        qc: Iterable[QCStatus] | QCStatus = "good"
) -> pl.DataFrame:
    """Return sequence × resolution level features for GLMM.

    One row per (sequence_id, resolution). Clustering outcome is summarised
    across all windows a sequence appears in at that resolution, giving:

    - prop_clustered : proportion of windows where sequence was in a
                       non-singleton cluster (continuous 0–1)
    - ever_clustered : 1 if clustered in any window at this resolution,
                       0 if singleton in all windows (binary)

    Parameters
    ----------
    qc:
        QC filter values to pass to load_analysis_columns.
    """
    _validate_qc(qc)

    _NON_SINGLETON_RE = r"\|C\d+$"

    cols = [
        "window_id", "window_idx", "cluster_id", "sequence_id", "resolution",
        "wn_start_date", "wn_end_date", "wn_prop_sequenced", "pango_lineage", "who_voc", "nextclade_qc",
        "dz_simd_rank", "dz_simd_quintile", "dz_simd_decile",
        "dz_simd_income_rank", "dz_simd_employment_rank", "dz_simd_education_rank",
        "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank",
        "dz_simd_housing_rank", "age_band", "age_midpoint", "is_female",
        "is_vaccinated", "collection_date",
    ]

    df = load_analysis_columns(cols, resolution=None, qc=qc)

    df = df.with_columns(
        pl.col("cluster_id")
        .cast(pl.Utf8)
        .str.contains(_NON_SINGLETON_RE)
        .cast(pl.Int64)
        .alias("_non_singleton")
    )

    # --- stable per-sequence columns (take first — constant within sequence) -
    stable = [
        "collection_date", "age_band", "age_midpoint", "is_female",
        "is_vaccinated", "who_voc", "dz_simd_quintile",
        "dz_simd_decile", "pango_lineage", *list(_simd_domain().values())
    ]

    # --- aggregate to sequence × resolution ----------------------------------
    # prop_clustered: mean of binary outcome across windows
    # ever_clustered: max of binary outcome across windows
    agg_exprs: list[pl.Expr] = (
            [pl.col(c).first() for c in stable]
            + [
                pl.col("_non_singleton").mean().alias("prop_clustered"),
                pl.col("_non_singleton").max().cast(pl.Int64).alias("ever_clustered"),
                pl.col("wn_prop_sequenced").mean(),
                pl.col("window_id").n_unique().alias("n_windows"),
                pl.col("wn_start_date").min().alias("wn_start_date"),
                pl.col("wn_end_date").max().alias("wn_end_date"),
            ]
    )

    out = df.group_by(["sequence_id", "resolution"]).agg(agg_exprs)

    for dom in _simd_domain().keys():
        col = "dz_simd_rank" if dom == "overall" else f"dz_simd_{dom}_rank"
        out = out.with_columns(_standardise(out[col]).alias(f"{dom}_zscore"))

    out = _with_wave(out, "collection_date")
    out = _with_policy(out, "collection_date")
    out = out.with_columns(pl.col("who_voc").fill_null("other").alias("who_voc"))

    # --- age group mapping ---------------------------------------------------
    age_map = {
        "00-04": "00-09", "05-09": "00-09",
        "10-14": "10-19", "15-19": "10-19",
        "20-24": "20-39", "25-29": "20-39", "30-34": "20-39", "35-39": "20-39",
        "40-44": "40-59", "45-49": "40-59", "50-54": "40-59", "55-59": "40-59",
        "60-64": "60-74", "65-69": "60-74", "70-74": "60-74",
        "75+": "elderly",
    }
    out = out.with_columns(
        pl.col("age_band").replace(age_map).alias("age_group")
    )

    return out.select([
        "sequence_id", "resolution", "collection_date", "wn_end_date", "wn_start_date",
        "wn_prop_sequenced", "age_band", "age_group", "age_midpoint", "is_female", "is_vaccinated",
        "pango_lineage", "who_voc", "n_windows", "wave", "policy_period",
        "policy_period_label", "policy_intensity", "prop_clustered", "ever_clustered",
        "dz_simd_quintile", "dz_simd_decile", *list(f"{dom}_zscore" for dom in _simd_domain().keys())
    ])
