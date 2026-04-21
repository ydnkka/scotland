"""Data-loading helpers for the clustering manuscripts.

The analysis dataset is ~7.4M sequence-rows; eager loading of every column is
wasteful and often impossible on a laptop. These helpers read only the columns
each figure needs, and preserve the usual filters (Nextclade QC, a chosen
primary Leiden resolution, valid WHO VOC labels) in one place so the papers
stay internally consistent.

All paths are resolved from the repository's top-level `config.yaml`, not
hard-coded, so moving the repo does not break the figure scripts.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

#: Primary Leiden resolution used in headline results. Override per figure if
#: the figure is an explicit sensitivity analysis over resolution.
PRIMARY_RESOLUTION: float = 0.3

#: Hardcoded WHO variant-of-concern epochs used as a fallback when the
#: master parquet is not available. All live code should use
#: :func:`get_voc_epochs` (or the :data:`VOC_EPOCHS` lazy attribute), which
#: prefers epochs derived from the actual `who_voc` × `collection_date`
#: distribution in the data — see :func:`derive_voc_epochs_from_data`.
VOC_EPOCHS_DEFAULT: list[tuple[str, str, str]] = [
    # (label, start_date, end_date)
    ("Pre-VOC",        "2020-07-01", "2020-11-30"),
    ("Alpha",          "2020-12-01", "2021-05-31"),
    ("Delta",          "2021-06-01", "2021-12-15"),
    ("Omicron BA.1",   "2021-12-16", "2022-02-28"),
    ("Omicron BA.2+",  "2022-03-01", "2023-02-28"),
]

# Canonical WHO VOCs we retain in scripts (others collapsed into 'Other').
KEPT_VOCS: tuple[str, ...] = ("Alpha", "Delta", "Omicron")


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


def repo_root(start: Path | None = None) -> Path:
    """Walk up from `start` (default: this file) until `config.yaml` is found."""
    p = (start or Path(__file__)).resolve()
    for cand in [p] + list(p.parents):
        if (cand / "config.yaml").exists():
            return cand
    raise FileNotFoundError("Could not locate config.yaml in any parent directory.")


@dataclass(frozen=True)
class Paths:
    root: Path
    analysis_dataset: Path
    @classmethod
    def from_config(cls, root: Path = None) -> "Paths":
        root = root or repo_root()
        with open(root / "config.yaml") as f:
            cfg = yaml.safe_load(f)
        proc = cfg["data"]["processed"]
        return cls(
            root=root,
            analysis_dataset=root / proc["analysis_dataset"]
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
) -> pd.DataFrame:
    """Read a narrow slice of the master sequence-level parquet.

    Parameters
    ----------
    columns : list of str
        Names of columns to read; `resolution` and `nextclade_qc` are added
        automatically if filtering is requested.
    resolution : float or None
        If provided, rows are restricted to that Leiden resolution.
    qc : iterable of str or None
        If provided, rows are restricted to these Nextclade QC statuses.
    paths : Paths or None
        If provided, use these paths instead of resolving from config.
    """
    paths: Paths = paths or Paths.from_config()
    need = set(columns)
    if resolution is not None:
        need.add("resolution")
    if qc is not None:
        need.add("nextclade_qc")
    df = pd.read_parquet(paths.analysis_dataset, columns=sorted(need))
    if resolution is not None:
        df = df[df["resolution"] == resolution]
    if qc is not None:
        df = df[df["nextclade_qc"].isin(list(qc))]
    return df.reset_index(drop=True)


@lru_cache(maxsize=4)
def load_cluster_demographic_features(
    min_size: int = 1,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: Iterable[str] | None = ("good",)
) -> pd.DataFrame:
    cols = [
        "window_idx", "window_id", "wn_mid_date",
        "sequence_id", "cluster_id", "pango_lineage",
        "age_midpoint", "is_female", "is_vaccinated", "vacc_dose_number",
    ]
    df = load_analysis_columns(cols, resolution=resolution, qc=qc)
    grp = df.groupby(["window_id", "cluster_id"], observed=True)

    out = grp.agg(
        n_sequences=("sequence_id", "nunique"),
        median_age=("age_midpoint", "median"),
        age_diversity=("age_midpoint", "std"),
        frac_female=("is_female", "mean"),
        frac_vaccinated=("is_vaccinated", "mean"),
        mean_vacc_dose=("vacc_dose_number", "mean"),
        wn_mid_date=("wn_mid_date", "first"),
        window_idx=("window_idx", "first"),
        pango_lineage=("pango_lineage", "first")
    ).reset_index()
    out["is_singleton"] = (out["n_sequences"] == 1).astype(int)
    out = out[out["n_sequences"] >= min_size]
    return out


@lru_cache(maxsize=4)
def load_cluster_simd_features(
    min_size: int = 1,
    resolution: float | None = PRIMARY_RESOLUTION,
    qc: Iterable[str] | None = ("good",)
) -> pd.DataFrame:
    """Return one row per (window_id, cluster_id) with its size, date, lineage, and SIMD features.
    """
    cols = [
        "window_idx", "window_id", "wn_mid_date",
        "sequence_id", "cluster_id", "pango_lineage",
        "datazone", "dz_simd_rank", "dz_simd_quintile", "dz_simd_decile",
    ]
    df = load_analysis_columns(cols, resolution=resolution, qc=qc)
    grp = df.groupby(["window_id", "cluster_id"], observed=True)
    def _mode(s: pd.Series):
        m = s.mode()
        return m.iloc[0] if len(m) > 0 else np.nan

    out = grp.agg(
        n_sequences=("sequence_id", "nunique"),
        n_datazones=("datazone", "nunique"),
        simd_rank_mean=("dz_simd_rank", "mean"),
        simd_rank_median=("dz_simd_rank", "median"),
        simd_quintile_mode=("dz_simd_quintile", _mode),
        simd_decile_mode=("dz_simd_decile", _mode),
        frac_deprived_q1=("dz_simd_quintile", lambda x: (x == 1).mean()),
        wn_mid_date=("wn_mid_date", "first"),
        window_idx=("window_idx", "first"),
        pango_lineage=("pango_lineage", lambda x: x.mode().iloc[0] if len(x) > 0 else None),
    ).reset_index()
    out["is_singleton"] = (out["n_sequences"] == 1).astype(int)
    out = out[out["n_sequences"] >= min_size]
    return out


# ---------------------------------------------------------------------------
# VOC epoch derivation from data
# ---------------------------------------------------------------------------


def _is_ba1(lineage) -> bool:
    """True for BA.1 and any BA.1.* sub-lineage (but not BA.10, BA.11, …)."""
    if not isinstance(lineage, str):
        return False
    return lineage == "BA.1" or lineage.startswith("BA.1.")


def _weekly_dominant_voc(
    dates: pd.Series,
    voc: pd.Series,
    *,
    dominance_threshold: float,
) -> pd.DataFrame:
    """For each ISO week, return the dominant WHO VOC and its share.

    Weeks whose dominant label's share is below `dominance_threshold`, or
    whose dominant label is missing / "None", are marked "Pre-VOC".
    """
    weeks = pd.to_datetime(dates).dt.to_period("W-SUN").dt.start_time
    labels = voc.fillna("None").astype(str).replace({"": "None"})
    tbl = pd.crosstab(weeks, labels)
    shares = tbl.div(tbl.sum(axis=1), axis=0)
    dominant = shares.idxmax(axis=1)
    share = shares.max(axis=1)
    out = pd.DataFrame(
        {"week": shares.index, "dominant": dominant.values, "share": share.values}
    )
    undetermined = (out["share"] < dominance_threshold) | out["dominant"].isin(["None"])
    out.loc[undetermined, "dominant"] = "Pre-VOC"
    return out.sort_values("week").reset_index(drop=True)


def _contiguous_runs(
    weekly: pd.DataFrame, min_weeks: int
) -> list[tuple[str, pd.Timestamp, pd.Timestamp]]:
    """Group consecutive weeks with the same dominant label.

    Runs shorter than `min_weeks` are merged into the preceding run (if any)
    or the following run (for a short leading segment). This keeps a one-week
    VOC blip from spawning a spurious epoch.
    """
    if weekly.empty:
        return []
    runs: list[list] = []  # (label, start, end)
    cur_label = weekly["dominant"].iloc[0]
    cur_start = weekly["week"].iloc[0]
    cur_end = cur_start
    for _, row in weekly.iloc[1:].iterrows():
        if row["dominant"] == cur_label:
            cur_end = row["week"]
        else:
            runs.append([cur_label, cur_start, cur_end])
            cur_label = row["dominant"]
            cur_start = row["week"]
            cur_end = row["week"]
    runs.append([cur_label, cur_start, cur_end])

    def _weeks(run) -> int:
        return int(((run[2] - run[1]).days // 7) + 1)

    # Absorb short runs into their longer neighbour.
    cleaned: list[list] = []
    for rn in runs:
        if _weeks(rn) < min_weeks and cleaned:
            cleaned[-1][2] = rn[2]
        else:
            cleaned.append(rn)
    # A short leading run now has no predecessor; fold it forward instead.
    if len(cleaned) >= 2 and _weeks(cleaned[0]) < min_weeks:
        cleaned[1][1] = cleaned[0][1]
        cleaned = cleaned[1:]
    # Collapse any adjacent runs that now share a label post-merge.
    collapsed: list[list] = []
    for rn in cleaned:
        if collapsed and collapsed[-1][0] == rn[0]:
            collapsed[-1][2] = rn[2]
        else:
            collapsed.append(rn)
    return [(lbl, s, e) for lbl, s, e in collapsed]


def _split_omicron_by_lineage(
    epochs: list[tuple[str, pd.Timestamp, pd.Timestamp]],
    df: pd.DataFrame,
) -> list[tuple[str, pd.Timestamp, pd.Timestamp]]:
    """Split the single Omicron epoch into BA.1 vs BA.2+ using pango_lineage.

    The split week is the first week in which BA.1's share falls below 50%.
    If no such crossover exists (or the Omicron run is too short), the epoch
    is kept as "Omicron" unchanged.
    """
    refined: list[tuple[str, pd.Timestamp, pd.Timestamp]] = []
    for lbl, s, e in epochs:
        if lbl != "Omicron":
            refined.append((lbl, s, e))
            continue
        sub = df[(df["collection_date"] >= s) & (df["collection_date"] <= e)]
        if len(sub) < 100:
            refined.append((lbl, s, e))
            continue
        week = pd.to_datetime(sub["collection_date"]).dt.to_period("W-SUN").dt.start_time
        ba1_share = (
            sub.assign(_week=week.values, _ba1=sub["pango_lineage"].map(_is_ba1))
            .groupby("_week")["_ba1"].mean().sort_index()
        )
        below = ba1_share[ba1_share < 0.5]
        if below.empty or len(ba1_share) < 4:
            refined.append((lbl, s, e))
            continue
        split_week = below.index[0]
        refined.append(("Omicron BA.1", s, split_week - pd.Timedelta(days=1)))
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
    1. Aggregate QC-passing sequences (one row per ``sequence_id``) by ISO
       week × ``who_voc``.
    2. In each week, find the dominant VOC. Weeks where no label holds at
       least ``dominance_threshold`` of the sequenced cases are labelled
       "Pre-VOC".
    3. Group consecutive weeks sharing a dominant label into a single epoch.
       Runs shorter than ``min_weeks`` are merged into an adjacent epoch.
    4. If ``split_omicron``, split the Omicron epoch into BA.1 and BA.2+
       using ``pango_lineage`` (crossover = first week BA.1 share < 50%).

    Returns
    -------
    list of ``(label, "YYYY-MM-DD", "YYYY-MM-DD")``. Falls back to
    :data:`VOC_EPOCHS_DEFAULT` if the master parquet cannot be loaded.
    """
    try:
        cols = ["sequence_id", "collection_date", "who_voc"]
        if split_omicron:
            cols.append("pango_lineage")
        df = load_analysis_columns(cols, resolution=resolution)
    except Exception:  # pragma: no cover — path/config failures
        return list(VOC_EPOCHS_DEFAULT)

    df = df.dropna(subset=["collection_date"]).drop_duplicates("sequence_id")
    if df.empty:
        return list(VOC_EPOCHS_DEFAULT)

    weekly = _weekly_dominant_voc(
        df["collection_date"], df["who_voc"],
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
    :func:`derive_voc_epochs_from_data`; set ``from_data=False`` to force the
    hardcoded :data:`VOC_EPOCHS_DEFAULT`. If data-driven derivation fails for
    any reason (missing parquet, unexpected schema), the default is returned.
    """
    if not from_data:
        return list(VOC_EPOCHS_DEFAULT)
    return list(derive_voc_epochs_from_data())


def __getattr__(name: str):
    """Lazy module-level attribute so `data.VOC_EPOCHS` uses derived values.

    Existing call sites (``from manuscripts.utils.data import VOC_EPOCHS``
    and ``data.VOC_EPOCHS``) are resolved through this hook on first access
    and cached thereafter by :func:`derive_voc_epochs_from_data`.
    """
    if name == "VOC_EPOCHS":
        return get_voc_epochs()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# ---------------------------------------------------------------------------
# Derived fields & helpers
# ---------------------------------------------------------------------------


def assign_epoch(dates: pd.Series) -> pd.Categorical:
    """Assign each row to a VOC epoch label based on its date."""
    epochs = get_voc_epochs()
    labels = pd.Series(np.nan, index=dates.index, dtype=object)
    for label, s, e in epochs:
        mask = (dates >= pd.Timestamp(s)) & (dates <= pd.Timestamp(e))
        labels.loc[mask] = label
    cats = [lbl for lbl, *_ in epochs]
    return pd.Categorical(labels, categories=cats, ordered=True)


def collapse_voc(voc: pd.Series) -> pd.Categorical:
    """Collapse rare VOC labels into 'Other' for plotting consistency."""
    keep = list(KEPT_VOCS)
    out = voc.where(voc.isin(keep), "Other")
    return pd.Categorical(out, categories=keep + ["Other"], ordered=True)


def canonical_cluster_sequence_counts(
    resolution: float = PRIMARY_RESOLUTION,
    qc: Iterable[str] | None = ("good",),
    *,
    paths: Paths | None = None,
) -> pd.DataFrame:
    """Return one row per (window_id, cluster_id) with its size, date, and lineage.

    Useful for cluster-level analyses where the sequence-level parquet is
    too wide. Cheaper than loading the full master parquet.
    """
    cols = ["window_id", "resolution", "cluster_id", "sequence_id",
            "wn_mid_date", "pango_lineage", "who_voc"]
    df = load_analysis_columns(cols, resolution=resolution, paths=paths)
    if qc is not None:
        df = df[df["nextclade_qc"].isin(list(qc))]
    grp = df.groupby(["window_id", "cluster_id"], observed=True)
    out = grp.agg(
        n_sequences=("sequence_id", "nunique"),
        wn_mid_date=("wn_mid_date", "first"),
        pango_lineage=("pango_lineage", "first"),
        who_voc=("who_voc", "first"),
    ).reset_index()
    out["is_singleton"] = (out["n_sequences"] == 1).astype(int)
    return out
