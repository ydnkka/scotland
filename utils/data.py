"""Shared data-loading helpers for the Scotland clustering analysis.

The functions in this module resolve repository paths from ``config.yaml``,
load selected columns from the processed analysis dataset, and apply the
rolling-window stride used by the notebooks.
"""

from __future__ import annotations

from collections.abc import Iterable as IterableABC
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import geopandas as gpd
import numpy as np
import pandas as pd
import yaml

from .policy import attach_period


__all__ = [
    "QCStatus",
    "CLADES",
    "CLADE_PALETTE",
    "VALID_QC_STATUSES",
    "load_analysis_columns",
    "load_pairwise_edges",
    "load_datazone_info",
    "pango_lineages_for_clades",
]


PRIMARY_RESOLUTION: float = 0.3

QCStatus = Literal["good", "mediocre", "bad"]
VALID_QC_STATUSES: set[str] = {"good", "mediocre", "bad"}

SIMD_GROUP_SIZES: dict[str, int] = {
    "quintile": 5,
    "decile": 10,
    "vigintile": 20,
}
SIMD_BASE_GROUP_COLUMNS: dict[str, int] = {
    "dz_simd_quintile": 5,
    "dz_simd_decile": 10,
    "dz_simd_vigintile": 20,
}
SIMD_DOMAIN_RANK_COLUMNS: dict[str, str] = {
    "income": "dz_simd_income_rank",
    "employment": "dz_simd_employment_rank",
    "education": "dz_simd_education_rank",
    "health": "dz_simd_health_rank",
    "access": "dz_simd_access_rank",
    "crime": "dz_simd_crime_rank",
    "housing": "dz_simd_housing_rank",
}
SIMD_GROUP_RANK_COLUMNS: dict[str, str] = {
    **{col: "dz_simd_rank" for col in SIMD_BASE_GROUP_COLUMNS},
    **{
        f"dz_simd_{domain}_{group_name}": rank_col
        for domain, rank_col in SIMD_DOMAIN_RANK_COLUMNS.items()
        for group_name in SIMD_GROUP_SIZES
    },
}
SIMD_GROUP_COLUMNS: dict[str, int] = {
    **SIMD_BASE_GROUP_COLUMNS,
    **{
        f"dz_simd_{domain}_{group_name}": n_groups
        for domain in SIMD_DOMAIN_RANK_COLUMNS
        for group_name, n_groups in SIMD_GROUP_SIZES.items()
    },
}
SIMD_COMPUTED_GROUP_COLUMNS: set[str] = (
    set(SIMD_GROUP_COLUMNS) - set(SIMD_BASE_GROUP_COLUMNS)
)

CLADES: dict[str, str] = {
    "20B": "20B",
    "20A": "20A",
    "20E": "20E (EU1)",
    "20I": "20I (Alpha)",
    "21K": "21K (Omicron)",
    "21J": "21J (Delta)",
    "21I": "21I (Delta)",
    "21L": "21L (Omicron)",
    "22B": "22B (Omicron)",
    "22A": "22A (Omicron)",
    "22C": "22C (Omicron)",
    "22E": "22E (Omicron)",
}

CLADE_PALETTE: dict[str, str] = {
    "20A": "#4477AA",
    "20B": "#66CCEE",
    "20E (EU1)": "#EE6677",
    "20I (Alpha)": "#117733",
    "21I (Delta)": "#AA3377",
    "21J (Delta)": "#CCBB44",
    "21K (Omicron)": "#63A227",
    "21L (Omicron)": "#BBBBBB",
    "22A (Omicron)": "#777777",
    "22B (Omicron)": "#EE7733",
    "22C (Omicron)": "#882255",
    "22E (Omicron)": "#332288",
    "Other": "#DDDDDD",
}

TEST_REASON_MAP = {
    "symptomatic-citizen": "symptomatic_citizen",
    "I have coronavirus symptoms": "symptomatic_citizen",
    "I live~ work or study in a lockdown area with a coronavirus outbreak": "symptomatic_citizen",
    "symptomatic-essential-worker": "symptomatic_essential_worker",
    "Im an essential worker": "symptomatic_essential_worker",
    "scotland-wales-keyworker": "symptomatic_essential_worker",
    "wales-keyworker": "symptomatic_essential_worker",
    "test-for-contact-tracing": "contact_tracing",
    "test-for-contact-tracing-app": "contact_tracing",
    "test-for-contact-self-referral": "contact_tracing",
    "for-symptomatic-household-member": "contact_tracing",
    "Ive been in contact with a person who has tested positive for coronavirus and Ive been asked to take a test by a contact tracer (Northern Ireland and Scotland)": "contact_tracing",
    "Ive been in contact with a person who has tested positive for coronavirus and have since developed symptoms": "contact_tracing",
    "confirmatory-positive-test": "confirmatory",
    "confirmatory-other-reason": "confirmatory",
    "confirmatory-test-unclear": "confirmatory",
    "confirmatory-test-borders": "confirmatory",
    "told-to-order-repeat-test": "confirmatory",
    "self-isolation-support-grant": "isolation_scheme",
    "isolation-testing-home": "isolation_scheme",
    "isolation-testing-facility": "isolation_scheme",
    "gp-healthcare-request": "clinical",
    "antiviral-order": "clinical",
    "dental-patient-testing": "clinical",
    "I have been told to have a test before I go into hospital~ for example~ for surgery": "clinical",
    "zoe-symptom-study": "surveillance_research",
    "contact-testing-study": "surveillance_research",
    "events-research-programme": "surveillance_research",
    "serial-testing": "surveillance_research",
    "ntrg-member": "surveillance_research",
    "local-council-request": "local_outbreak",
    "attended-outbreak-venue": "local_outbreak",
    "community-testing": "local_outbreak",
    "scotland-university": "local_outbreak",
    "wales-university": "local_outbreak",
    "green-traveller": "travel",
    "other": "other",
    "Other": "other",
    "none": "other",
    "do-not-know": "other",
    "general-cta-referral": "other",
    "personal-assistant": "other",
    "Im a visiting professional": "other",
    "asymptomatic-home-order": "other",
}

POLICY_ERA_BY_PERIOD = {
    "E0": "early_restriction_easing",
    "L1": "early_restriction_easing",
    "P1": "early_restriction_easing",
    "P2": "early_restriction_easing",
    "P3": "early_restriction_easing",
    "T1": "autumn_winter_restrictions",
    "F5": "autumn_winter_restrictions",
    "L2": "autumn_winter_restrictions",
    "SL": "spring_summer_2021_easing",
    "L3": "spring_summer_2021_easing",
    "L21": "spring_summer_2021_easing",
    "L0": "spring_summer_2021_easing",
    "NN": "near_normal_delta",
    "OM": "omicron_response",
    "FE": "omicron_response",
    "PR": "post_restriction",
}

# Aligning 5-year census/health bands to Scottish infectious disease brackets
# "Because raw data was aggregated in 5-year intervals, age 15 was grouped into the 15-24 young-adult band rather than the 5-15 school-age band, closely approximating Public Health Scotland's youth surveillance frameworks."
AGE_GROUP_MAP = {
    "00-04": "00-04",  # Infants/Toddlers (RSV/Rotavirus)
    "05-09": "05-14",  # School-age (Primary school mixing)
    "10-14": "05-14",  # School-age (Secondary school mixing)
    
    # NOTE: 15-19 includes 15 (school) and 16-19 (young adult). 
    # In Scottish data, 15-19 is usually kept together or grouped with 20-24 
    # to capture the broader "Youth/Higher Education" transition.
    "15-19": "15-24",  
    "20-24": "15-24",  # Young Adults / University / High Social Mixing
    
    "25-29": "25-64",  # Working-age Adults
    "30-34": "25-64",
    "35-39": "25-64",
    "40-44": "25-64",
    "45-49": "25-64",
    "50-54": "25-64",
    "55-59": "25-64",
    "60-64": "25-64",
    
    "65-69": "65-74",  # "Young-Old" / Post-retirement / Elevated COVID Risk
    "70-74": "65-74",  
    
    "75+"  : "75+",    # Older Adults / Highest Clinical Vulnerability
}

def repo_root(start: Path | None = None) -> Path:
    """Walk up from *start* until ``config.yaml`` is found."""
    p = (start or Path(__file__)).resolve()

    for cand in [p, *p.parents]:
        if (cand / "config.yaml").exists():
            return cand

    raise FileNotFoundError("Could not locate config.yaml in any parent directory.")


@dataclass(frozen=True)
class Paths:
    """Resolved paths to commonly used processed data files."""

    root: Path
    analysis_dataset: Path
    pairwise_distances_dataset: Path
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
            pairwise_distances_dataset=(
                root
                / proc.get(
                    "pairwise_distances_dataset",
                    "data/processed/pairwise_distances_dataset",
                )
            ),
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


def _as_list(value: object | Iterable[object] | None) -> list[object] | None:
    """Return scalar or iterable input as a list, treating strings as scalars."""
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    if isinstance(value, IterableABC):
        return list(value)
    return [value]


def _normalise_str_values(
    values: object | Iterable[object] | None,
    *,
    name: str,
) -> list[str] | None:
    """Return string filter values with duplicates removed in input order."""
    raw_values = _as_list(values)
    if raw_values is None:
        return None

    out: list[str] = []
    for value in raw_values:
        if value is None:
            raise ValueError(f"{name} cannot contain None")
        out.append(str(value))

    return list(dict.fromkeys(out))


def _normalise_clades(
    clades: object | Iterable[object] | None,
) -> list[str] | None:
    """Return raw Nextclade clade labels, accepting display labels as aliases."""
    values = _normalise_str_values(clades, name="clades")
    if values is None:
        return None

    display_to_raw = {display: raw for raw, display in CLADES.items()}
    return [display_to_raw.get(value, value) for value in values]


def _normalise_window_ids(
    windows: object | Iterable[object] | None,
) -> list[str] | None:
    """Normalise window IDs to the processed pairwise format, e.g. ``W095``."""
    raw_windows = _as_list(windows)
    if raw_windows is None:
        return None

    out: list[str] = []
    for window in raw_windows:
        if window is None:
            raise ValueError("windows cannot contain None")

        if isinstance(window, (int, np.integer)):
            window_idx = int(window)
            if window_idx < 1:
                raise ValueError("window indices must be positive")
            out.append(f"W{window_idx:03d}")
            continue

        value = str(window).strip()
        if not value:
            raise ValueError("windows cannot contain empty strings")

        upper_value = value.upper()
        if upper_value.startswith("W") and upper_value[1:].isdigit():
            window_idx = int(upper_value[1:])
            if window_idx < 1:
                raise ValueError("window indices must be positive")
            out.append(f"W{window_idx:03d}")
        elif value.isdigit():
            window_idx = int(value)
            if window_idx < 1:
                raise ValueError("window indices must be positive")
            out.append(f"W{window_idx:03d}")
        else:
            out.append(value)

    return list(dict.fromkeys(out))


def _select_window_stride(
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


def _apply_window_stride(
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

    retained = _select_window_stride(df[window_col], stride=stride, offset=offset)
    old_to_new = {old: new + 1 for new, old in enumerate(retained)}

    out = df.loc[df[window_col].isin(retained)].copy()
    if renumber:
        out[window_col] = out[window_col].map(old_to_new)
        if window_id_col in out.columns:
            out[window_id_col] = out[window_col].apply(lambda x: f"W{x:03d}")
    return out.reset_index(drop=True)


def _requested_simd_group_columns(
    columns: Iterable[str] | None,
    *,
    all_cols: bool = False,
) -> set[str]:
    """Return requested SIMD grouping columns that support population weighting."""
    if all_cols:
        return set(SIMD_GROUP_COLUMNS)
    if columns is None:
        return set()
    return set(columns) & set(SIMD_GROUP_COLUMNS)


def _computed_simd_group_columns(
    simd_cols: Iterable[str],
    *,
    weighted: bool,
) -> set[str]:
    """Return group columns that must be computed rather than read as stored."""
    simd_cols = set(simd_cols)
    if weighted:
        return simd_cols
    return simd_cols & SIMD_COMPUTED_GROUP_COLUMNS


def _required_simd_group_source_columns(
    simd_cols: Iterable[str],
    *,
    weighted: bool,
) -> set[str]:
    """Return rank/population columns needed to compute SIMD groups."""
    simd_cols = set(simd_cols)
    if not simd_cols:
        return set()

    required = {SIMD_GROUP_RANK_COLUMNS[col] for col in simd_cols}
    if weighted:
        required.add("dz_population")
    return required


def _apply_simd_groups(
    df: pd.DataFrame,
    simd_cols: Iterable[str],
    *,
    weighted: bool,
    pop_col: str = "dz_population",
) -> pd.DataFrame:
    """Replace requested SIMD group columns with rank-derived groups.

    The ranking itself is unchanged: Data Zones are sorted by the matching
    Scottish Government SIMD rank. With ``weighted=True``, the sorted zones are
    split into equal population shares; otherwise, they are split into equal
    numbers of Data Zones.
    """
    simd_cols = list(simd_cols)
    if not simd_cols:
        return df

    unknown = set(simd_cols) - set(SIMD_GROUP_COLUMNS)
    if unknown:
        raise ValueError(f"Unsupported SIMD grouping column(s): {sorted(unknown)}")

    missing = _required_simd_group_source_columns(simd_cols, weighted=weighted) - set(
        df.columns
    )
    if missing:
        raise KeyError(
            f"SIMD grouping requires source columns: {sorted(missing)}"
        )

    out = df.copy()
    for rank_col in sorted({SIMD_GROUP_RANK_COLUMNS[col] for col in simd_cols}):
        rank_cols = [
            col for col in simd_cols if SIMD_GROUP_RANK_COLUMNS[col] == rank_col
        ]
        valid = out[rank_col].notna()
        if weighted:
            valid &= out[pop_col].notna()

        if not valid.all():
            missing_n = int((~valid).sum())
            raise ValueError(
                f"Cannot compute SIMD groups from {rank_col!r} with {missing_n} "
                "missing rank/population row(s)."
            )

        if weighted and (out[pop_col] < 0).any():
            raise ValueError(f"Negative values found in {pop_col}")

        ordered = out.sort_values(rank_col, ascending=True)
        if weighted:
            total_pop = ordered[pop_col].sum()
            if total_pop <= 0:
                raise ValueError("Total population must be greater than zero")
            cum_prop = ordered[pop_col].cumsum() / total_pop
        else:
            cum_prop = pd.Series(
                np.arange(1, len(ordered) + 1, dtype=float) / len(ordered),
                index=ordered.index,
            )

        for col in rank_cols:
            n_groups = SIMD_GROUP_COLUMNS[col]
            group_values = np.ceil(cum_prop * n_groups).astype(int).clip(1, n_groups)
            out.loc[ordered.index, col] = group_values

    return out


def _simd_group_lookup(
    paths: Paths,
    simd_cols: Iterable[str],
    *,
    weighted: bool,
) -> pd.DataFrame:
    """Build a datazone lookup with requested SIMD grouping columns."""
    simd_cols = list(simd_cols)
    source_cols = _required_simd_group_source_columns(simd_cols, weighted=weighted)
    lookup = pd.read_parquet(
        paths.geography,
        columns=list(source_cols),
    )
    lookup = _apply_simd_groups(lookup, simd_cols, weighted=weighted)
    return lookup.reset_index()[["datazone", *simd_cols]]


def _attach_simd_groups(
    df: pd.DataFrame,
    simd_cols: Iterable[str],
    *,
    paths: Paths,
    weighted: bool,
) -> pd.DataFrame:
    """Overwrite SIMD group columns in an analysis frame using a datazone lookup."""
    simd_cols = list(simd_cols)
    if not simd_cols:
        return df
    if "datazone" not in df.columns:
        raise KeyError("'datazone' is required to attach SIMD groups")

    out = df.copy()
    lookup = _simd_group_lookup(paths, simd_cols, weighted=weighted).set_index(
        "datazone"
    )

    for col in simd_cols:
        out[col] = out["datazone"].map(lookup[col])

    return out


def load_analysis_columns(
    columns: Iterable[str] | None = None,
    all_cols: bool = False,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: QCStatus | Iterable[QCStatus] | None = "good",
    add_policy: bool = False,
    window_stride: int | None = None,
    window_offset: int = 0,
    renumber_windows: bool = True,
    weighted_simd: bool = True,
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
    weighted_simd:
        If True, replace any requested SIMD grouping columns with
        population-weighted groups while retaining the same column names. This
        includes the overall columns ``dz_simd_quintile``, ``dz_simd_decile``,
        and ``dz_simd_vigintile``, plus computed domain columns such as
        ``dz_simd_income_quintile``. The weighting uses all datazones in the
        configured geography parquet, not just rows retained after
        sequence-level filters. If False, computed domain group columns are
        still derived from their rank columns, but without population weights.

    Notes
    -----
    ``Available columns``:
    "window_idx", "window_id", "wn_start_date", "wn_mid_date",
    "wn_end_date", "wn_no_sequences", "wn_positive_tests",
    "wn_prop_sequenced", "sequence_id", "patient_id", "resolution",
    "cluster_id", "cluster_size", "cluster_n_datazones",
    "cluster_start_date", "cluster_end_date", "cluster_duration_days",
    "collection_date", "datazone", "dz_xcoord", "dz_ycoord", "sex",
    "is_female", "age_band", "age_group", "age_midpoint", "is_vaccinated",
    "vacc_dose_number", "vacc_date_prior", "vacc_product_name",
    "vacc_booster", "days_since_vaccination", "test_reason",
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

    Computed on request:
    ``dz_simd_<domain>_<group>``, where domain is one of
    ``income``, ``employment``, ``education``, ``health``, ``access``,
    ``crime``, or ``housing`` and group is ``quintile``, ``decile``, or
    ``vigintile``.

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
    simd_cols = _requested_simd_group_columns(requested, all_cols=all_cols)
    computed_simd_cols = _computed_simd_group_columns(
        simd_cols,
        weighted=weighted_simd,
    )
    computed_age_cols = {"age_group"} & requested

    if columns is not None:
        need.update(requested - computed_simd_cols - computed_age_cols)

    if "age_group" in requested:
        need.add("age_band")

    if window_stride is not None:
        need.add("window_idx")
        if "window_id" in requested:
            need.add("window_id")

    if resolution is not None:
        need.add("resolution")

    if qc_values is not None:
        need.add("nextclade_qc")

    output_columns = set(need) | computed_simd_cols
    if computed_simd_cols:
        need.add("datazone")

    read_columns = None if all_cols else list(need)

    df = pd.read_parquet(paths.analysis_dataset, columns=read_columns)

    if all_cols or "age_group" in requested:
        df["age_group"] = df["age_band"].map(AGE_GROUP_MAP).fillna("unknown")
        if not all_cols and "age_band" not in requested:
            df = df.drop(columns="age_band")

    if resolution is not None:
        df = df.loc[df["resolution"] == resolution]

    if qc_values is not None:
        df = df.loc[df["nextclade_qc"].isin(qc_values)]

    if window_stride is not None:
        df = _apply_window_stride(
            df,
            stride=window_stride,
            offset=window_offset,
            renumber=renumber_windows,
        )

    df = df.reset_index(drop=True)

    if computed_simd_cols:
        df = _attach_simd_groups(
            df,
            computed_simd_cols,
            paths=paths,
            weighted=weighted_simd,
        )
        if not all_cols:
            drop_cols = {"datazone"} - output_columns
            df = df.drop(columns=[c for c in drop_cols if c in df.columns])

    if add_policy:
        df = attach_period(df, "collection_date")
        df["policy_era"] = (
            df["policy_period"].map(POLICY_ERA_BY_PERIOD).fillna(df["policy_period"])
        )

    if "test_reason" in df.columns:
        df["test_reason"] = df["test_reason"].replace(TEST_REASON_MAP).fillna("other")

    fill_values = {
        "is_vaccinated": 0.0,
        "days_since_vaccination": -1,
        "vacc_booster": 0,
        "who_voc": "Other/Non-VOC",
    }

    for col, value in fill_values.items():
        if col in df.columns:
            df[col] = df[col].fillna(value)

    return df


def pango_lineages_for_clades(
    clades: str | Iterable[str],
    *,
    windows: str | int | Iterable[str | int] | None = None,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: QCStatus | Iterable[QCStatus] | None = "good",
    paths: Paths | None = None,
) -> list[str]:
    """Resolve Nextclade clade labels to Pango lineages in the analysis dataset.

    Parameters
    ----------
    clades:
        One or more raw clade labels such as ``"21L"``. Display labels from
        ``CLADES`` such as ``"21L (Omicron)"`` are also accepted.
    windows:
        Optional rolling-window IDs or indices used to restrict the lookup.
        Integer-like inputs are normalised to the ``W001`` style used by the
        pairwise dataset.
    resolution:
        Leiden resolution filter applied to the analysis dataset before
        resolving lineages. Pass ``None`` to skip this filter.
    qc:
        Nextclade QC filter applied before resolving lineages. Pass ``None`` to
        skip this filter.
    paths:
        Optional resolved path bundle. Mostly useful for tests.

    Returns
    -------
    list[str]
        Sorted unique Pango lineage labels associated with the requested clades.
    """
    paths = paths or Paths.from_config()
    clade_values = _normalise_clades(clades)
    if clade_values is None:
        raise ValueError("clades must contain at least one value")

    window_ids = _normalise_window_ids(windows)
    qc_values = _normalise_qc(qc)

    need = {"clade", "pango_lineage"}
    if window_ids is not None:
        need.add("window_id")
    if resolution is not None:
        need.add("resolution")
    if qc_values is not None:
        need.add("nextclade_qc")

    df = pd.read_parquet(paths.analysis_dataset, columns=list(need))

    if resolution is not None:
        df = df.loc[df["resolution"] == resolution]
    if qc_values is not None:
        df = df.loc[df["nextclade_qc"].isin(qc_values)]
    if window_ids is not None:
        df = df.loc[df["window_id"].isin(window_ids)]

    lineages = sorted(
        df.loc[df["clade"].isin(clade_values), "pango_lineage"]
        .dropna()
        .astype(str)
        .unique()
    )

    if not lineages:
        detail = f"clades={clade_values}"
        if window_ids is not None:
            detail += f", windows={window_ids}"
        raise ValueError(f"No Pango lineages found for {detail}.")

    return lineages


def _resolve_pairwise_lineages(
    *,
    clades: str | Iterable[str] | None,
    pango_lineages: str | Iterable[str] | None,
    windows: list[str] | None,
    resolution: float | None,
    qc: QCStatus | Iterable[QCStatus] | None,
    paths: Paths,
) -> list[str] | None:
    """Resolve direct and clade-derived Pango lineage filters."""
    direct_lineages = _normalise_str_values(
        pango_lineages,
        name="pango_lineages",
    )

    clade_lineages = None
    if clades is not None:
        clade_lineages = pango_lineages_for_clades(
            clades,
            windows=windows,
            resolution=resolution,
            qc=qc,
            paths=paths,
        )

    if direct_lineages is None:
        return clade_lineages
    if clade_lineages is None:
        return direct_lineages

    retained = sorted(set(direct_lineages) & set(clade_lineages))
    if not retained:
        raise ValueError(
            "No Pango lineages remain after intersecting direct lineage and "
            "clade-derived lineage filters."
        )
    return retained


def _pairwise_filter_expression(
    *,
    window_ids: list[str] | None,
    pango_lineages: list[str] | None,
    compatibility_threshold: float | None,
    score_column: str,
):
    """Build a PyArrow dataset filter expression."""
    import pyarrow.compute as pc

    expr = None

    if window_ids is not None:
        window_expr = pc.field("window_id").isin(window_ids)
        expr = window_expr if expr is None else expr & window_expr

    if pango_lineages is not None:
        lineage_expr = pc.field("pango_lineage").isin(pango_lineages)
        expr = lineage_expr if expr is None else expr & lineage_expr

    if compatibility_threshold is not None:
        score_expr = pc.field(score_column) > float(compatibility_threshold)
        expr = score_expr if expr is None else expr & score_expr

    return expr


def _read_pairwise_dataset(
    path: Path,
    *,
    columns: list[str],
    window_ids: list[str] | None,
    pango_lineages: list[str] | None,
    compatibility_threshold: float | None,
    score_column: str,
) -> pd.DataFrame:
    """Read pairwise edges with PyArrow dataset projection and filter pushdown."""
    import pyarrow.dataset as ds

    dataset = ds.dataset(path, format="parquet")
    table = dataset.to_table(
        columns=columns,
        filter=_pairwise_filter_expression(
            window_ids=window_ids,
            pango_lineages=pango_lineages,
            compatibility_threshold=compatibility_threshold,
            score_column=score_column,
        ),
    )

    return table.to_pandas()


def load_pairwise_edges(
    columns: Iterable[str] | None = None,
    *,
    windows: str | int | Iterable[str | int] | None = None,
    clades: str | Iterable[str] | None = None,
    pango_lineages: str | Iterable[str] | None = None,
    compatibility_threshold: float | None = 0.001,
    score_column: str = "epilink_compatibility",
    pairwise_dataset: str | Path | None = None,
    clade_resolution: float | None = PRIMARY_RESOLUTION,
    clade_qc: QCStatus | Iterable[QCStatus] | None = "good",
) -> pd.DataFrame:
    """Load pairwise EpiLink edges using PyArrow pushdown filters.

    Parameters
    ----------
    columns:
        Pairwise columns to return. Defaults to ``window_id``,
        ``pango_lineage``, ``id1``, ``id2``, and ``epilink_compatibility``.
    windows:
        Optional window IDs or indices. Inputs such as ``95``, ``"95"``,
        ``"W95"``, and ``"W095"`` are all normalised to ``"W095"``.
    clades:
        Optional Nextclade clade labels. These are resolved to Pango lineages by
        looking up ``clade`` and ``pango_lineage`` in the individual analysis
        dataset, then applying the resulting lineage list to the pairwise
        dataset.
    pango_lineages:
        Optional direct Pango lineage filter. If both ``clades`` and
        ``pango_lineages`` are supplied, the filters are intersected.
    compatibility_threshold:
        Sparsification threshold for ``score_column``. Rows are retained where
        ``score_column > compatibility_threshold``. Pass ``None`` to skip score
        sparsification.
    score_column:
        Pairwise compatibility score column used for sparsification.
    pairwise_dataset:
        Optional override for the pairwise parquet dataset path.
    clade_resolution, clade_qc:
        Filters applied to the analysis dataset when resolving ``clades`` to
        Pango lineages. They do not filter the pairwise rows directly.

    Returns
    -------
    pandas.DataFrame
        Pairwise edge rows matching the requested filters.

    Notes
    -----
    Broad filters can still return very large edge sets. Prefer combining
    ``windows``, ``clades`` or ``pango_lineages``, and
    ``compatibility_threshold`` for interactive notebook work.
    """
    if compatibility_threshold is not None and compatibility_threshold < 0:
        raise ValueError("compatibility_threshold must be non-negative or None")

    paths = Paths.from_config()
    pairwise_path = (
        Path(pairwise_dataset) if pairwise_dataset else paths.pairwise_distances_dataset
    )
    if not pairwise_path.exists():
        raise FileNotFoundError(f"Pairwise dataset not found: {pairwise_path}")

    selected_columns = list(
        columns
        if columns is not None
        else ["window_id", "pango_lineage", "id1", "id2", score_column]
    )

    window_ids = _normalise_window_ids(windows)
    lineage_values = _resolve_pairwise_lineages(
        clades=clades,
        pango_lineages=pango_lineages,
        windows=window_ids,
        resolution=clade_resolution,
        qc=clade_qc,
        paths=paths,
    )

    return _read_pairwise_dataset(
        pairwise_path,
        columns=selected_columns,
        window_ids=window_ids,
        pango_lineages=lineage_values,
        compatibility_threshold=compatibility_threshold,
        score_column=score_column,
    )


def load_datazone_info(
    columns: Iterable[str],
    weighted_simd: bool = True,
) -> gpd.GeoDataFrame:
    """Read a narrow slice of the datazone information parquet.

    Parameters
    ----------
    columns:
        Names of columns to read. ``datazone`` and ``geometry`` are added
        automatically.
    weighted_simd:
        If True, replace any requested SIMD grouping columns with
        population-weighted groups while retaining the same column names. This
        includes the overall columns ``dz_simd_quintile``, ``dz_simd_decile``,
        and ``dz_simd_vigintile``, plus computed domain columns such as
        ``dz_simd_income_quintile``. If False, computed domain group columns
        are still derived from their rank columns, but without population
        weights.

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

    Computed on request:
    ``dz_simd_<domain>_<group>``, where domain is one of
    ``income``, ``employment``, ``education``, ``health``, ``access``,
    ``crime``, or ``housing`` and group is ``quintile``, ``decile``, or
    ``vigintile``.

    Returns
    -------
    geopandas.GeoDataFrame
    """
    paths = Paths.from_config()

    need = {"datazone", "geometry"}
    requested = set(columns)
    simd_cols = _requested_simd_group_columns(requested)
    computed_simd_cols = _computed_simd_group_columns(
        simd_cols,
        weighted=weighted_simd,
    )
    need.update(requested - computed_simd_cols)

    output_columns = set(need)
    if computed_simd_cols:
        output_columns |= computed_simd_cols
        need.update(
            _required_simd_group_source_columns(
                computed_simd_cols,
                weighted=weighted_simd,
            )
        )

    df = gpd.read_parquet(paths.geography, columns=list(need))

    if computed_simd_cols:
        df = _apply_simd_groups(df, computed_simd_cols, weighted=weighted_simd)
        drop_cols = (
            _required_simd_group_source_columns(
                computed_simd_cols,
                weighted=weighted_simd,
            )
            - output_columns
        )
        df = df.drop(columns=[c for c in drop_cols if c in df.columns])

    return df  # type: ignore[return-value]
