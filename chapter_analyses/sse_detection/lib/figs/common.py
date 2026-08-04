"""Shared plumbing for SSE detection figures."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

from ..sse.config import (
    BAYESIAN_OUTPUT_DIR,
    FIGURE_DIR,
    PROJECT_ROOT,
    SSE_OUTPUT_DIR,
    TABLE_DIR,
)

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import (  # noqa: F401
    CLADES,
    load_policy_calendar,
    policy_era_labels,
    policy_order,
)
from utils.style import *

DEFAULT_TABLE_DIR = SSE_OUTPUT_DIR
DEFAULT_RESULT_TABLE_DIR = TABLE_DIR

HIGH_PRIORITY = frozenset(
    {
        "high_priority_both_axes",
        "high_priority_burst",
        "high_priority_burden",
    }
)

_POLICY_DAILY = load_policy_calendar()
POLICY_ORDER = {
    "policy_era": policy_order("policy_era"),
    "period_code": policy_order("policy_period"),
    "policy_period": policy_order("policy_period"),
}
POLICY_LABELS = policy_era_labels()
POLICY_STRINGENCY = (
    _POLICY_DAILY.groupby("policy_era", sort=False, observed=True)["stringency_index"]
    .mean()
    .to_dict()
)

POLICY_STRINGENCY_CMAP = plt.get_cmap("RdYlGn_r")
POLICY_STRINGENCY_NORM = Normalize(vmin=1, vmax=100)
POLICY_COLORS = {
    period: POLICY_STRINGENCY_CMAP(POLICY_STRINGENCY_NORM(stringency))
    for period, stringency in POLICY_STRINGENCY.items()
}


@dataclass(frozen=True)
class Paths:
    table_dir: Path = DEFAULT_TABLE_DIR
    figure_dir: Path = FIGURE_DIR
    bayesian_result_dir: Path = BAYESIAN_OUTPUT_DIR
    result_table_dir: Path = DEFAULT_RESULT_TABLE_DIR


def read_table(paths: Paths, name: str) -> pd.DataFrame:
    searched = []
    for directory in dict.fromkeys((paths.table_dir, paths.result_table_dir)):
        for suffix, reader in (("parquet", pd.read_parquet), ("csv", pd.read_csv)):
            path = directory / f"{name}.{suffix}"
            searched.append(path)
            if path.exists():
                return reader(path)
    locations = ", ".join(str(path) for path in searched)
    raise FileNotFoundError(f"Missing table {name!r}; searched: {locations}")


def latex_table_path(paths: Paths, name: str) -> Path:
    """Return a LaTeX table-fragment path under the figure directory."""
    paths.figure_dir.mkdir(parents=True, exist_ok=True)
    return paths.figure_dir / f"{name}.tex"


def styled_save_figure(
    fig: Figure,
    paths: Paths,
    name: str,
    *,
    width: WIDTHS = "double",
) -> dict[str, Path]:
    paths.figure_dir.mkdir(parents=True, exist_ok=True)
    return save_figure(
        fig,
        paths.figure_dir / name,
        width=width,
        save_pdf=True,
        save_png=True,
    )


def add_policy_bands(ax: Axes, window_coverage: pd.DataFrame) -> None:
    if "policy_era" not in window_coverage.columns:
        return
    work = window_coverage[["wn_mid_date", "policy_era"]].dropna().copy()
    if work.empty:
        return
    work = work.sort_values("wn_mid_date")
    half_window = pd.Timedelta(days=7)
    for period, group in work.groupby("policy_era", sort=False):
        start = group["wn_mid_date"].min() - half_window
        end = group["wn_mid_date"].max() + half_window
        ax.axvspan(
            start,
            end,
            color=POLICY_COLORS.get(str(period), "#f0f0f0"),
            alpha=0.18,
            lw=0,
            zorder=0,
        )


def sort_by_policy(df: pd.DataFrame, column: str) -> pd.DataFrame:
    order = {policy: idx for idx, policy in enumerate(POLICY_ORDER.get(column, []))}
    out = df.copy()
    out["_policy_sort"] = out[column].astype(str).map(order).fillna(999)
    return out.sort_values(["_policy_sort", column]).drop(columns="_policy_sort")


def date_axis(ax: Axes) -> None:
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=DEFAULT_TABLE_DIR,
        help="Directory containing SSE output tables.",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=FIGURE_DIR,
        help="Directory for generated SSE figures.",
    )
    parser.add_argument(
        "--bayesian-result-dir",
        type=Path,
        default=BAYESIAN_OUTPUT_DIR,
        help="Directory containing fitted Bayesian model outputs.",
    )
    parser.add_argument(
        "--result-table-dir",
        type=Path,
        default=DEFAULT_RESULT_TABLE_DIR,
        help="Directory for generated result tables.",
    )


def paths_from_args(args: argparse.Namespace) -> Paths:
    return Paths(
        table_dir=args.table_dir,
        figure_dir=args.figure_dir,
        bayesian_result_dir=args.bayesian_result_dir,
        result_table_dir=args.result_table_dir,
    )


def wilson(k: pd.Series, n: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Calculate the 95% Wilson score interval for binomial proportions.

    Args:
        k: Series of success counts.
        n: Series of total trials.

    Returns:
        A tuple of (lower_bound, upper_bound) as pandas Series, clipped to [0, 1].
    """
    z = 1.959963984540054
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (centre - half).clip(lower=0), (centre + half).clip(upper=1)
