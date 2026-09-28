"""Shared plumbing for genomic-network artifacts."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import add_policy_bands as add_policy_bands  # noqa: F401
from utils import (
    load_policy_calendar,
    policy_era_labels,
    policy_order,
)
from utils import window_idx_from_id as window_idx_from_id
from utils.style import WIDTHS
from utils.style import add_panel_labels as add_panel_labels
from utils.style import new_figure as new_figure
from utils.style import save_figure

DEFAULT_TABLE_DIR = PROJECT_ROOT / "analyses/genomic_networks/results/tables"
FIGURE_DIR = PROJECT_ROOT / "results/figures"

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

ATTRIBUTE_ORDER = [
    "Sex",
    "Age band",
    "Age group",
    "SIMD quintile",
    "Urban/rural class",
    "Health board",
    "Local authority",
]

POLICY_STRINGENCY_CMAP = plt.get_cmap("RdYlGn_r")
POLICY_STRINGENCY_NORM = Normalize(vmin=0, vmax=100)
POLICY_COLORS = {
    period: POLICY_STRINGENCY_CMAP(POLICY_STRINGENCY_NORM(stringency))
    for period, stringency in POLICY_STRINGENCY.items()
}


@dataclass(frozen=True)
class Paths:
    table_dir: Path
    figure_dir: Path = FIGURE_DIR
    publication_table_dir: Path = PROJECT_ROOT / "results/tables"


def read_table(paths: Paths, name: str) -> pd.DataFrame:
    parquet_path = paths.table_dir / f"{name}.parquet"
    csv_path = paths.table_dir / f"{name}.csv"
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    if csv_path.exists():
        return pd.read_csv(csv_path)
    raise FileNotFoundError(f"Missing table: {name}")


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


def ordered_policy_values(values: pd.Series, column: str) -> list[object]:
    order = {policy: idx for idx, policy in enumerate(POLICY_ORDER.get(column, []))}
    unique = [value for value in pd.unique(values.dropna())]
    return sorted(unique, key=lambda value: (order.get(str(value), 999), str(value)))


def sort_by_policy(df: pd.DataFrame, column: str) -> pd.DataFrame:
    order = {policy: idx for idx, policy in enumerate(POLICY_ORDER.get(column, []))}
    out = df.copy()
    out["_policy_sort"] = out[column].astype(str).map(order).fillna(999)
    return out.sort_values(["_policy_sort", column]).drop(columns="_policy_sort")


def date_axis(ax: Axes) -> None:
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--figure-dir", type=Path, default=FIGURE_DIR)
    parser.add_argument(
        "--publication-table-dir",
        type=Path,
        default=PROJECT_ROOT / "results/tables",
        help="Directory for generated LaTeX tables.",
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=DEFAULT_TABLE_DIR,
        help="Directory containing analyses/genomic_networks result tables.",
    )


def paths_from_args(args: argparse.Namespace) -> Paths:
    return Paths(
        publication_table_dir=args.publication_table_dir,
        table_dir=args.table_dir,
        figure_dir=args.figure_dir,
    )


__all__ = [
    "Paths",
    "add_common_args",
    "paths_from_args",
    "read_table",
    "styled_save_figure",
    "sort_by_policy",
    "date_axis",
    "add_panel_labels",
    "new_figure",
    "add_policy_bands",
    "DEFAULT_TABLE_DIR",
    "FIGURE_DIR",
    "ATTRIBUTE_ORDER",
    "POLICY_ORDER",
    "POLICY_LABELS",
    "POLICY_COLORS",
    "POLICY_STRINGENCY_CMAP",
    "POLICY_STRINGENCY_NORM",
    "ordered_policy_values",
    "window_idx_from_id",
]
