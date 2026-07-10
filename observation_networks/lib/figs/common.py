"""Shared plumbing for Chapter 4 observation-network artifacts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TABLE_DIR = Path(
    "/Users/ydnkka/Desktop/PhD Project/projects/scotland/"
    "observation_networks/results/tables"
)
FIGURE_DIR = REPO_ROOT / "thesis/figures/observation_networks"
LATEX_TABLE_DIR = REPO_ROOT / "thesis/tables/observation_networks"

POLICY_ORDER = [
    "P2",
    "P3",
    "T1",
    "F5",
    "L2",
    "SL",
    "L3",
    "L21",
    "L0",
    "NN",
    "OM",
    "FE",
    "PR",
]

ATTRIBUTE_ORDER = [
    "Sex",
    "Age band",
    "Age group",
    "SIMD quintile",
    "Urban/rural class",
    "Health board",
    "Local authority",
]

# Match the Chapter 2 policy strip in scotland/surveillance/policy_sequences_over_time.py.
POLICY_STRINGENCY = {
    "P2": 76.45619047619047,
    "P3": 67.64571428571442,
    "T1": 64.80999999999996,
    "F5": 70.88968750000001,
    "L2": 85.53586206896566,
    "SL": 69.90375000000003,
    "L3": 58.330000000000005,
    "L21": 56.13063492063501,
    "L0": 52.77999999999999,
    "NN": 31.637857142857182,
    "OM": 34.83857142857142,
    "FE": 19.786666666666626,
    "PR": 8.184418604651126,
}
POLICY_STRINGENCY_CMAP = plt.get_cmap("RdYlGn_r")
POLICY_STRINGENCY_NORM = Normalize(vmin=1, vmax=100)
POLICY_COLORS = {
    period: POLICY_STRINGENCY_CMAP(POLICY_STRINGENCY_NORM(stringency))
    for period, stringency in POLICY_STRINGENCY.items()
}


@dataclass(frozen=True)
class Paths:
    table_dir: Path
    figure_dir: Path = FIGURE_DIR
    latex_table_dir: Path = LATEX_TABLE_DIR


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_table(paths: Paths, name: str) -> pd.DataFrame:
    parquet_path = paths.table_dir / f"{name}.parquet"
    csv_path = paths.table_dir / f"{name}.csv"
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    if csv_path.exists():
        return pd.read_csv(csv_path)
    raise FileNotFoundError(f"Missing table: {name}")


def save_figure(
    fig: mpl.figure.Figure,
    paths: Paths,
    name: str,
    *,
    tight: bool = True,
) -> None:
    paths.figure_dir.mkdir(parents=True, exist_ok=True)
    if tight:
        fig.tight_layout()
    fig.savefig(paths.figure_dir / f"{name}.png", bbox_inches="tight")
    fig.savefig(paths.figure_dir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )


def window_idx_from_id(values: pd.Series) -> pd.Series:
    extracted = values.astype(str).str.extract(r"(\d+)")[0]
    return pd.to_numeric(extracted, errors="coerce")


def ordered_policy_values(values: pd.Series) -> list[object]:
    order = {policy: idx for idx, policy in enumerate(POLICY_ORDER)}
    unique = [value for value in pd.unique(values.dropna())]
    return sorted(unique, key=lambda value: (order.get(str(value), 999), str(value)))


def sort_by_policy(df: pd.DataFrame, column: str = "policy_period") -> pd.DataFrame:
    order = {policy: idx for idx, policy in enumerate(POLICY_ORDER)}
    out = df.copy()
    out["_policy_sort"] = out[column].astype(str).map(order).fillna(999)
    return out.sort_values(["_policy_sort", column]).drop(columns="_policy_sort")


def add_policy_bands(ax: mpl.axes.Axes, window_coverage: pd.DataFrame) -> None:
    if "policy_period" not in window_coverage.columns:
        return
    work = window_coverage[["wn_mid_date", "policy_period"]].dropna().copy()
    if work.empty:
        return
    work = work.sort_values("wn_mid_date")
    half_window = pd.Timedelta(days=7)
    for period, group in work.groupby("policy_period", sort=False):
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


def date_axis(ax: mpl.axes.Axes) -> None:
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=DEFAULT_TABLE_DIR,
        help="Directory containing observation_networks result tables.",
    )


def paths_from_args(args: argparse.Namespace) -> Paths:
    return Paths(table_dir=args.table_dir)

