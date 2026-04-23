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
from datetime import date, timedelta
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import polars as pl
import yaml


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

#: Primary Leiden resolution used in headline results.
PRIMARY_RESOLUTION: float = 0.3

#: Hardcoded fallback VOC epochs (used when the parquet is unavailable).
VOC_EPOCHS_DEFAULT: list[tuple[str, str, str]] = [
    ("Pre-VOC",        "2020-07-01", "2020-11-30"),
    ("Alpha",          "2020-12-01", "2021-05-31"),
    ("Delta",          "2021-06-01", "2021-12-15"),
    ("Omicron BA.1",   "2021-12-16", "2022-02-28"),
    ("Omicron BA.2+",  "2022-03-01", "2023-02-28"),
]


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


@lru_cache(maxsize=4)
def load_cluster_features(
    min_size: int = 1,
    resolution: float = PRIMARY_RESOLUTION,
    qc: tuple[str, ...] = ("good",),
) -> pl.DataFrame:
    """Return one row per (window_id, cluster_id) with size, date, lineage, and SIMD features."""

    cols = [
        "window_id", "window_idx", "cluster_id", "sequence_id",
        "wn_mid_date", "wn_prop_sequenced", "who_voc", "pango_lineage", "nextclade_qc",
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
        pl.col("age_midpoint").std().alias("age_diversity"),
        pl.col("is_female").mean().alias("frac_female"),
        pl.col("is_vaccinated").mean().alias("frac_vaccinated"),
        # SIMD quintile / decile summary
        # mode() returns all modal values; sort().first() picks the smallest on ties
        pl.col("dz_simd_quintile").drop_nulls().mode().sort().first().alias("simd_quintile_mode"),
        pl.col("dz_simd_quintile").std().alias("simd_quintile_std"),
        pl.col("dz_simd_decile").drop_nulls().mode().sort().first().alias("simd_decile_mode"),
        pl.col("dz_simd_decile").std().alias("simd_decile_std"),
    ]

    for dom, col in _simd_domain().items():
        agg_exprs.append(pl.col(col).mean().alias(f"simd_{dom}_mean"))

    out = df.group_by(["window_id", "cluster_id"]).agg(agg_exprs)

    out = out.with_columns(
        (pl.col("n_sequences") == 1).cast(pl.Int8).alias("is_singleton"),
    )

    # Epoch assignment --------------------------------------------------------
    out = out.with_columns(
        assign_epoch(out["wn_mid_date"]).alias("epoch")
    )
    out = out.filter(pl.col("epoch").is_not_null())

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

    if min_size > 1:
        out = out.filter(pl.col("n_sequences") >= min_size)

    return out


@lru_cache(maxsize=4)
def load_individual_features(
    qc: tuple[str, ...] = ("good",),
) -> pl.DataFrame:
    """Return one row per sequence_id, averaging window-level metrics across
    all (window_id, resolution) combinations the sequence appears in.
    """
    _NON_SINGLETON_RE = r"C\d+$"
    cols = [
        "window_id", "window_idx", "cluster_id", "patient_id", "sequence_id", "resolution",
        "wn_mid_date", "wn_prop_sequenced", "who_voc", "pango_lineage", "nextclade_qc",
        "datazone", "dz_simd_rank", "dz_simd_quintile", "dz_simd_decile",
        "dz_simd_income_rank", "dz_simd_employment_rank", "dz_simd_education_rank",
        "dz_simd_health_rank", "dz_simd_access_rank", "dz_simd_crime_rank",
        "dz_simd_housing_rank", "age_band", "is_female", "is_vaccinated", "collection_date",
    ]
    df = load_analysis_columns(cols, resolution=None, qc=qc)
    df = df.with_columns(
        pl.col("cluster_id").cast(pl.Utf8).str.contains(_NON_SINGLETON_RE)
        .cast(pl.UInt8).alias("_non_singleton")
    )
    stable = [
        "patient_id", "collection_date",
        "age_band", "is_female", "is_vaccinated",
        "datazone", "dz_simd_quintile", "dz_simd_decile",
        *list(_simd_domain().values()), "pango_lineage",
    ]
    agg_exprs: list[pl.Expr] = (
        [pl.col(c).first() for c in stable]
        + [
            pl.col("_non_singleton").mean().alias("non_singleton_cluster_fraction"),
            pl.col("_non_singleton").sum().cast(pl.Int64).alias("non_singleton_k"),
            pl.col("_non_singleton").count().alias("non_singleton_n"),
            pl.col("wn_prop_sequenced").mean(),
            pl.col("wn_prop_sequenced").log().mean().alias("log_seq_prop"),
            pl.col("window_id").n_unique().alias("n_windows"),
        ]
    )
    out = df.group_by("sequence_id").agg(agg_exprs)

    # Epoch assignment --------------------------------------------------------
    out = out.with_columns(
        assign_epoch(out["collection_date"]).alias("epoch")
    )
    out = out.filter(pl.col("epoch").is_not_null())
    assert out["log_seq_prop"].is_finite().all(), (
        "Non-finite offset values — check wn_prop_sequenced > 0"
    )
    return out


# ---------------------------------------------------------------------------
# VOC epoch derivation from data
# ---------------------------------------------------------------------------


def _is_ba1(lineage: object) -> bool:
    """True for BA.1 and any BA.1.* sub-lineage (but not BA.10, BA.11, …)."""
    if not isinstance(lineage, str):
        return False
    return lineage == "BA.1" or lineage.startswith("BA.1.")


def _weekly_dominant_voc(
    dates: pl.Series,
    voc: pl.Series,
    *,
    dominance_threshold: float,
) -> pl.DataFrame:
    """For each ISO week, return the dominant WHO VOC and its share.

    Weeks whose dominant label's share is below *dominance_threshold*, or
    whose dominant label is missing / "None", are marked "Pre-VOC".
    """
    df = pl.DataFrame({
        "date": dates.cast(pl.Date),
        "voc":  voc.fill_null("None").cast(pl.Utf8),
    })

    # Truncate to Monday-start week (ISO).  One day off from pandas W-SUN but
    # immaterial for epoch boundary derivation.
    df = df.with_columns(
        pl.col("date").dt.truncate("1w").alias("week")
    )

    # Count rows per (week, voc), then compute each label's share within the week
    counts = (
        df.group_by(["week", "voc"])
        .agg(pl.len().alias("n"))
        .join(
            df.group_by("week").agg(pl.len().alias("total")),
            on="week",
        )
        .with_columns(
            (pl.col("n") / pl.col("total")).alias("share")
        )
    )

    # Dominant = the label with highest share in each week
    dominant = (
        counts
        .group_by("week")
        .agg(
            pl.col("voc").sort_by("share", descending=True).first().alias("dominant"),
            pl.col("share").max().alias("share"),
        )
        .sort("week")
    )

    # Weeks with no clear leader → "Pre-VOC"
    dominant = dominant.with_columns(
        pl.when(
            (pl.col("share") < dominance_threshold)
            | pl.col("dominant").is_in(["None", ""])
        )
        .then(pl.lit("Pre-VOC"))
        .otherwise(pl.col("dominant"))
        .alias("dominant")
    )

    return dominant


def _contiguous_runs(
    weekly: pl.DataFrame,
    min_weeks: int,
) -> list[tuple[str, date, date]]:
    """Group consecutive weeks with the same dominant label.

    Runs shorter than *min_weeks* are merged into the preceding run (if any)
    or the following run (for a short leading segment).
    """
    if weekly.is_empty():
        return []

    rows = list(weekly.sort("week").iter_rows(named=True))

    runs: list[list] = []
    cur_label: str = rows[0]["dominant"]
    cur_start: date = rows[0]["week"]
    cur_end:   date = cur_start

    for row in rows[1:]:
        if row["dominant"] == cur_label:
            cur_end = row["week"]
        else:
            runs.append([cur_label, cur_start, cur_end])
            cur_label = row["dominant"]
            cur_start = row["week"]
            cur_end   = row["week"]
    runs.append([cur_label, cur_start, cur_end])

    def _weeks(run: list) -> int:
        return int(((run[2] - run[1]).days // 7) + 1)

    # Absorb short runs into their longer neighbour
    cleaned: list[list] = []
    for rn in runs:
        if _weeks(rn) < min_weeks and cleaned:
            cleaned[-1][2] = rn[2]
        else:
            cleaned.append(rn)

    # A short leading run has no predecessor — fold it forward instead
    if len(cleaned) >= 2 and _weeks(cleaned[0]) < min_weeks:
        cleaned[1][1] = cleaned[0][1]
        cleaned = cleaned[1:]

    # Collapse adjacent runs that now share a label after merging
    collapsed: list[list] = []
    for rn in cleaned:
        if collapsed and collapsed[-1][0] == rn[0]:
            collapsed[-1][2] = rn[2]
        else:
            collapsed.append(rn)

    return [(lbl, s, e) for lbl, s, e in collapsed]


def _split_omicron_by_lineage(
    epochs: list[tuple[str, date, date]],
    df: pl.DataFrame,
) -> list[tuple[str, date, date]]:
    """Split the single Omicron epoch into BA.1 vs BA.2+ using pango_lineage.

    The split week is the first week in which BA.1's share falls below 50%.
    If no crossover exists, the epoch is kept as "Omicron" unchanged.
    """
    refined: list[tuple[str, date, date]] = []

    for lbl, s, e in epochs:
        if lbl != "Omicron":
            refined.append((lbl, s, e))
            continue

        sub = df.filter(
            (pl.col("collection_date") >= s) & (pl.col("collection_date") <= e)
        )
        if len(sub) < 100:
            refined.append((lbl, s, e))
            continue

        sub = sub.with_columns(
            pl.col("collection_date").cast(pl.Date).dt.truncate("1w").alias("_week"),
            pl.col("pango_lineage")
            .map_elements(_is_ba1, return_dtype=pl.Boolean)
            .alias("_ba1"),
        )

        ba1_share = (
            sub.group_by("_week")
            .agg(pl.col("_ba1").cast(pl.Float64).mean().alias("ba1_share"))
            .sort("_week")
        )

        below = ba1_share.filter(pl.col("ba1_share") < 0.5)
        if below.is_empty() or len(ba1_share) < 4:
            refined.append((lbl, s, e))
            continue

        split_week: date = below["_week"][0]
        refined.append(("Omicron BA.1", s, split_week - timedelta(days=1)))
        refined.append(("Omicron BA.2+", split_week, e))

    return refined


@lru_cache(maxsize=4)
def derive_voc_epochs_from_data(
    *,
    resolution: float = PRIMARY_RESOLUTION,
    dominance_threshold: float = 0.5,
    min_weeks: int = 3,
    split_omicron: bool = True,
) -> list[tuple[str, str, str]]:
    """Derive VOC epochs from the actual data.

    Algorithm
    ---------
    1. Aggregate QC-passing sequences by ISO week × ``who_voc``.
    2. In each week, find the dominant VOC.  Weeks where no label holds at
       least *dominance_threshold* of sequences are labelled "Pre-VOC".
    3. Group consecutive weeks sharing a dominant label into a single epoch.
       Runs shorter than *min_weeks* are merged into an adjacent epoch.
    4. If *split_omicron*, split the Omicron epoch into BA.1 and BA.2+
       using ``pango_lineage`` (crossover = first week BA.1 share < 50%).

    Returns
    -------
    list of ``(label, "YYYY-MM-DD", "YYYY-MM-DD")``.  Falls back to
    :data:`VOC_EPOCHS_DEFAULT` if the master parquet cannot be loaded.
    """
    try:
        cols = ["sequence_id", "collection_date", "who_voc"]
        if split_omicron:
            cols.append("pango_lineage")
        df = load_analysis_columns(cols, resolution=resolution)
    except Exception:  # pragma: no cover — path/config failures
        return list(VOC_EPOCHS_DEFAULT)

    df = df.filter(pl.col("collection_date").is_not_null()).unique("sequence_id")
    if df.is_empty():
        return list(VOC_EPOCHS_DEFAULT)

    weekly = _weekly_dominant_voc(
        df["collection_date"],
        df["who_voc"],
        dominance_threshold=dominance_threshold,
    )
    epochs = _contiguous_runs(weekly, min_weeks=min_weeks)
    if split_omicron:
        epochs = _split_omicron_by_lineage(epochs, df)
    if not epochs:
        return list(VOC_EPOCHS_DEFAULT)

    return [
        (lbl, s.strftime("%Y-%m-%d"), e.strftime("%Y-%m-%d"))
        for lbl, s, e in epochs
    ]


def get_voc_epochs(*, from_data: bool = True) -> list[tuple[str, str, str]]:
    """Return the authoritative list of VOC epochs.

    By default the epochs are derived from the dataset via
    :func:`derive_voc_epochs_from_data`; set ``from_data=False`` to force
    the hardcoded :data:`VOC_EPOCHS_DEFAULT`.
    """
    if not from_data:
        return list(VOC_EPOCHS_DEFAULT)
    return list(derive_voc_epochs_from_data())


def assign_epoch(dates: pl.Series) -> pl.Series:
    """Assign each row to a VOC epoch label based on its date.

    Returns a ``pl.Series`` with dtype ``pl.Enum(cats)`` — the ordered
    categorical equivalent of ``pd.Categorical(..., ordered=True)``.
    Null values indicate rows that fall outside all epoch windows.
    """
    epochs = get_voc_epochs()
    cats = [lbl for lbl, *_ in epochs]
    epoch_enum = pl.Enum(cats)

    # Build label column via successive when/otherwise passes
    label_expr: pl.Expr = pl.lit(None, dtype=pl.Utf8)
    for lbl, s, e in epochs:
        start_d = date.fromisoformat(s)
        end_d   = date.fromisoformat(e)
        label_expr = (
            pl.when(
                (pl.col("_date") >= start_d) & (pl.col("_date") <= end_d)
            )
            .then(pl.lit(lbl))
            .otherwise(label_expr)
        )

    result = (
        pl.DataFrame({"_date": dates.cast(pl.Date)})
        .with_columns(label_expr.alias("epoch"))
        ["epoch"]
        .cast(epoch_enum)
    )
    return result


def __getattr__(name: str):
    """Lazy module-level attribute so ``data.VOC_EPOCHS`` uses derived values."""
    if name == "VOC_EPOCHS":
        return get_voc_epochs()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")