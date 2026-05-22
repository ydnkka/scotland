from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import geopandas as gpd
import pandas as pd
import yaml

from .policy import attach_period


PRIMARY_RESOLUTION: float = 0.3

QCStatus = Literal["good", "mediocre", "bad"]
VALID_QC_STATUSES: set[str] = {"good", "mediocre", "bad"}


def repo_root(start: Path | None = None) -> Path:
    """Walk up from *start* until ``config.yaml`` is found."""
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
    def from_config(cls, root: Path | None = None) -> "Paths":
        root = root or repo_root()

        with open(root / "config.yaml") as f:
            cfg = yaml.safe_load(f)

        proc = cfg["data"]["processed"]

        return cls(
            root=root,
            analysis_dataset=root / proc["analysis_dataset"],
            geography=root / proc["geography"],
        )


def _normalise_qc(qc: QCStatus | Iterable[QCStatus] | None) -> list[str] | None:
    """Return QC values as a list, validating accepted statuses."""
    if qc is None:
        return None

    qc_values: list[str] = [str(qc)] if isinstance(qc, str) else [str(x) for x in qc]

    invalid = set(qc_values) - VALID_QC_STATUSES
    if invalid:
        raise ValueError(
            f"Invalid QC status value(s): {sorted(invalid)}. "
            f"Accepted values are: {sorted(VALID_QC_STATUSES)}."
        )

    return qc_values


def select_window_stride(
    windows: Iterable[int],
    stride: int = 1,
    *,
    offset: int = 0,
) -> list[int]:
    """Select rolling-window indices by position from sorted unique windows.

    Examples
    --------
    ``select_window_stride(range(1, 25), 2)`` returns windows
    ``[1, 3, 5, ...]``. This is intentionally position-based, equivalent to
    ``sorted_windows[offset::stride]``, so notebooks do not need to assume that
    retained windows can be described by a modulo rule.
    """
    if stride < 1:
        raise ValueError("stride must be a positive integer")
    if offset < 0:
        raise ValueError("offset must be non-negative")
    if offset >= stride:
        raise ValueError("offset must be less than stride")

    sorted_windows = sorted(pd.Series(list(windows)).dropna().astype(int).unique())
    return sorted_windows[offset::stride]


def apply_window_stride(
    df: pd.DataFrame,
    stride: int = 1,
    *,
    offset: int = 0,
    window_col: str = "window_idx",
    window_id_col: str = "window_id",
    renumber: bool = True,
) -> pd.DataFrame:
    """Filter to a rolling-window stride and optionally renumber retained windows."""
    if window_col not in df.columns:
        raise KeyError(f"{window_col!r} is required for window stride filtering")

    retained = select_window_stride(df[window_col], stride=stride, offset=offset)
    old_to_new = {old: new + 1 for new, old in enumerate(retained)}

    out = df.loc[df[window_col].isin(retained)].copy()
    if renumber:
        out[window_col] = out[window_col].map(old_to_new)
        if window_id_col in out.columns:
            out[window_id_col] = out[window_col].apply(lambda x: f"W{x:03d}")
    return out.reset_index(drop=True)


def load_analysis_columns(
    columns: Iterable[str] | None = None,
    all_cols: bool = False,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: QCStatus | Iterable[QCStatus] | None = "good",
    add_policy: bool = False,
    window_stride: int | None = None,
    window_offset: int = 0,
    renumber_windows: bool = True,
) -> pd.DataFrame:
    """Read a narrow slice of the master sequence-level parquet.

    Parameters
    ----------
    columns:
        Names of columns to read. ``sequence_id``, ``collection_date``, and
        ``pango_lineage`` are added automatically. ``resolution`` and
        ``nextclade_qc`` are also added automatically when filtering is used.
    all_cols:
        If True, ignore ``columns`` and read all columns.
    resolution:
        If provided, rows are restricted to that Leiden resolution.
    qc:
        If provided, rows are restricted to these Nextclade QC statuses.
        Accepted values are ``"good"``, ``"mediocre"``, and ``"bad"``.
        Pass ``None`` to skip QC filtering.
    add_policy:
        If True, attach policy period labels ``policy_period``,
        ``policy_period_label``, and ``policy_intensity`` using the configured
        policy periods.
    window_stride:
        If provided, retain sorted unique ``window_idx`` values by positional
        stride using ``windows[window_offset::window_stride]``. For example,
        ``window_stride=2`` keeps original windows 1, 3, 5, ... when the source
        windows are numbered consecutively from 1.
    window_offset:
        Positional offset into the sorted unique windows before striding.
    renumber_windows:
        If True with ``window_stride``, renumber retained ``window_idx`` values
        to 1..N and rebuild ``window_id`` where present.

    Notes
    -----
    ``Available columns``:
    "window_idx", "window_id", "wn_start_date", "wn_mid_date",
    "wn_end_date", "wn_no_sequences", "wn_positive_tests",
    "wn_prop_sequenced", "sequence_id", "patient_id", "resolution",
    "cluster_id", "cluster_size", "cluster_n_datazones",
    "cluster_start_date", "cluster_end_date", "cluster_duration_days",
    "collection_date", "datazone", "dz_xcoord", "dz_ycoord", "sex",
    "is_female", "age_band", "age_midpoint", "is_vaccinated",
    "vacc_dose_number", "vacc_date_prior", "vacc_product_name",
    "vacc_booster", "days_since_vaccination", "s_gene_status",
    "is_reinfection", "pango_lineage", "clade", "who_voc", "nextclade_qc",
    "dz_population", "dz_working_age_population", "dz_area_km2",
    "dz_population_density", "dz_simd_rank", "dz_simd_quintile",
    "dz_simd_decile", "dz_simd_vigintile", "dz_simd_income_rank",
    "dz_simd_employment_rank", "dz_simd_education_rank",
    "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank",
    "dz_simd_housing_rank", "dz_urban_rural_class", "dz_local_authority",
    "dz_local_authority_code", "dz_health_board", "dz_health_board_code",
    "dz_total_tests", "dz_positive_tests", "dz_negative_tests",
    "dz_pcr_positive_tests", "dz_lfd_positive_tests", "dz_care_home_tests",
    "dz_test_positivity", "dz_7d_test_positivity", "dz_total_vaccinated",
    "dz_cum_vaccinated", "dz_cum_prop_vaccinated", "dz_cum_sequences",
    "dz_cum_positive_tests", "dz_cum_prop_sequenced",
    "dz_cum_incidence_per_capita", "hb_daily_positive",
    "hb_cumulative_positive", "hb_hospital_admissions",
    "hb_hospital_occupancy", "hb_icu_admissions", "hb_icu_occupancy_lt28d",
    "hb_icu_occupancy_ge28d", "hb_daily_reinfections",
    "hb_reinfection_rate"

    Returns
    -------
    pandas.DataFrame

    Raises
    ------
    ValueError
        If any value in ``qc`` is not one of the accepted QC statuses.
    """
    paths = Paths.from_config()
    qc_values = _normalise_qc(qc)

    need = {"sequence_id", "collection_date", "pango_lineage"}
    requested = set(columns or [])

    if columns is not None:
        need.update(requested)

    if window_stride is not None:
        need.add("window_idx")
        if "window_id" in requested:
            need.add("window_id")

    if resolution is not None:
        need.add("resolution")

    if qc_values is not None:
        need.add("nextclade_qc")

    read_columns = None if all_cols else list(need)

    df = pd.read_parquet(paths.analysis_dataset, columns=read_columns)

    if resolution is not None:
        df = df.loc[df["resolution"] == resolution]

    if qc_values is not None:
        df = df.loc[df["nextclade_qc"].isin(qc_values)]

    if window_stride is not None:
        df = apply_window_stride(
            df,
            stride=window_stride,
            offset=window_offset,
            renumber=renumber_windows,
        )

    df = df.reset_index(drop=True)

    if add_policy:
        df = attach_period(df, "collection_date")

    return df


def load_datazone_info(columns: Iterable[str]) -> gpd.GeoDataFrame:
    """Read a narrow slice of the datazone information parquet.

    Parameters
    ----------
    columns:
        Names of columns to read. ``datazone`` and ``geometry`` are added
        automatically.

    Notes
    -----
    ``Available columns``:
    "geometry", "dz_xcoord", "dz_ycoord", "dz_area_km2", "dz_population",
    "dz_working_age_population", "dz_urban_rural_class",
    "dz_local_authority_code", "dz_local_authority", "dz_health_board_code",
    "dz_health_board", "dz_simd_rank", "dz_simd_quintile", "dz_simd_decile",
    "dz_simd_vigintile", "dz_simd_income_rank", "dz_simd_employment_rank",
    "dz_simd_education_rank", "dz_simd_health_rank", "dz_simd_access_rank",
    "dz_simd_crime_rank", "dz_simd_housing_rank"

    Returns
    -------
    geopandas.GeoDataFrame
    """
    paths = Paths.from_config()

    need = {"datazone", "geometry"}
    need.update(columns)

    return gpd.read_parquet(paths.geography, columns=list(need))
