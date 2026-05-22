"""Data-loading helpers."""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Literal, Sequence

import geopandas as gpd
import pandas as pd
import polars as pl
import yaml

from .policy import attach_period

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

PRIMARY_RESOLUTION: float = 0.3


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
    simd: Path
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
            simd=root / proc["simd"],
            geography=root / proc["geography"]
        )



# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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

# ---------------------------------------------------------------------------
# Column-level loaders
# ---------------------------------------------------------------------------

def load_analysis_columns(
    columns: Iterable[str] | None = None,
    all_cols: bool = False,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: Iterable[QCStatus] | QCStatus = "good",
    add_policy: bool = False,
) -> pl.DataFrame:
    """Read a narrow slice of the master sequence-level parquet.

    Parameters
    ----------
    columns:
        Names of columns to read; ``resolution`` and ``nextclade_qc`` are
        added automatically when filtering is requested.
        ``sequence_id``, ``collection_date``, and ``pango_lineage`` are also added automatically.
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
    add_policy:
        If True, attach policy period labels ``policy_period``,
        ``policy_period_label``,  and ``policy_intensity``
        using the configured policy periods.

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

    need = {"sequence_id", "collection_date", "pango_lineage"}
    if columns is not None:
        need = need.union(columns)

    if resolution is not None:
        need.add("resolution")
    if qc is not None:
        need.add("nextclade_qc")

    lf = pl.scan_parquet(paths.analysis_dataset)

    if resolution is not None:
        lf = lf.filter(pl.col("resolution") == resolution)

    if qc is not None:
        if isinstance(qc, str):
            qc: Iterable[QCStatus] = tuple((qc,))
        lf = lf.filter(pl.col("nextclade_qc").is_in(list(qc)))

    if not all_cols:
        lf = lf.select(list(need))

    df = lf.collect()

    if add_policy:
        df = attach_period(df, "collection_date")

    return df


def load_analysis_columns_pandas(
    columns: Iterable[str] | None = None,
    all_cols: bool = False,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: Iterable[QCStatus] | QCStatus = "good",
    add_policy: bool = False,
) -> pd.DataFrame:
    """Pandas wrapper around ``load_analysis_columns``."""
    return load_analysis_columns(
        columns=columns,
        all_cols=all_cols,
        resolution=resolution,
        qc=qc,
        add_policy=add_policy,
    ).to_pandas()


def load_simd_columns(
    columns: Iterable[str] | None = None,
    all_cols: bool = False,
) -> pl.DataFrame:
    """Read a narrow slice of the processed SIMD parquet."""
    paths = Paths.from_config()
    lf = pl.scan_parquet(paths.simd)

    if not all_cols:
        need = {"datazone"}
        if columns is not None:
            need = need.union(columns)
        lf = lf.select(list(need))

    return lf.collect()


def load_simd_columns_pandas(
    columns: Iterable[str] | None = None,
    all_cols: bool = False,
) -> pd.DataFrame:
    """Pandas wrapper around ``load_simd_columns``."""
    return load_simd_columns(columns=columns, all_cols=all_cols).to_pandas()


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
