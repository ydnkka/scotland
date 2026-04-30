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

WAVE = [
    ("Alpha_n27129", "2020-11-04", "2021-08-11"),
    ("Delta_n99781", "2021-04-02", "2022-03-12"),
    ("Omicron_n144698", "2021-09-01", "2023-02-10")
]
    # Actual Omicron first sequence date is "2021-02-20",
    # move it forward because of long period of minimal transmission


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
    datazones: Path

    @classmethod
    def from_config(cls, root: Path = None) -> "Paths":
        root = root or repo_root()
        with open(root / "config.yaml") as f:
            cfg = yaml.safe_load(f)
        proc = cfg["data"]["processed"]
        return cls(
            root=root,
            analysis_dataset=root / proc["analysis_dataset"],
            datazones=root / "data/processed/datazones_information.parquet"
        )


# ---------------------------------------------------------------------------
# Column-level loaders
# ---------------------------------------------------------------------------


QCStatus = Literal["good", "mediocre", "bad"]
_VALID_QC: frozenset[str] = frozenset({"good", "mediocre", "bad"})
def _validate_qc(qc: Iterable[QCStatus] | QCStatus) -> None:
    if qc is not None:
        if isinstance(qc, str):
            qc :Iterable[QCStatus] = tuple((qc,))
        else:
            qc: Iterable[QCStatus] = tuple(qc)
        invalid = set(qc) - _VALID_QC
        if invalid:
            raise ValueError(
                f"Invalid QC status values: {invalid!r}. "
                f"Must be one or more of {sorted(_VALID_QC)}"
            )

def load_analysis_columns(
    columns: Iterable[str],
    *,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: Iterable[QCStatus] | QCStatus = "good"
) -> pl.DataFrame:
    """Read a narrow slice of the master sequence-level parquet.

    Parameters
    ----------
    columns:
        Names of columns to read; ``resolution`` and ``nextclade_qc`` are
        added automatically when filtering is requested.
        ``sequence_id``, ``collection_date`` are also added automatically.
    resolution:
        If provided, rows are restricted to that Leiden resolution.
    qc:
        If provided, rows are restricted to these Nextclade QC statuses.
        Accepted values: ``"good"``, ``"mediocre"``, ``"bad"``.
        Pass ``None`` to skip QC filtering entirely.

    Raises
    ------
    ValueError
        If any value in ``qc`` is not one of the accepted QC statuses.
    """
    _validate_qc(qc)

    paths = Paths.from_config()

    need = set(columns)
    if resolution is not None:
        need.add("resolution")
    if qc is not None:
        need.add("nextclade_qc")
        
    need = need.union(["sequence_id", "collection_date"])

    df = pl.read_parquet(paths.analysis_dataset, columns=sorted(need))

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
    """Read a narrow slice of the datazone information parquet."""
    paths = Paths.from_config()
    need = {"datazone", "geometry"}
    need = need.union(columns)
    return gpd.read_parquet(paths.datazones, columns=sorted(need))


# ---------------------------------------------------------------------------
# Aggregation-level loaders
# ---------------------------------------------------------------------------

def _standardise(values):
    """Negate z-score so higher values = greater deprivation."""
    return -(values - values.mean()) / values.std()

def _shannon_entropy(counts: pl.Series) -> float:
    """Shannon entropy (bits) from a series of counts."""
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts / total
    # mask zeros to avoid log(0)
    p = p.filter(p > 0)
    return -(p * p.log(base=2)).sum()

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
        "75+":   "elderly",
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


def assign_wave(dates: pl.Series) -> pl.Series:
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
    return df.with_columns(assign_wave(df[date_col]).alias("wave"))


def _with_policy(df: pl.DataFrame, date_col: str) -> pl.DataFrame:
    """Attach configured policy-period labels using the given date column."""
    return attach_period(df, date_col)
