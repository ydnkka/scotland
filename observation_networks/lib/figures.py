"""Figure helpers for Chapter 4 observation/network outputs."""

from __future__ import annotations

from pathlib import Path
import sys

from matplotlib.ticker import PercentFormatter
import pandas as pd
import seaborn as sns

from .config import FIGURES_DIR, PROJECT_ROOT


if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import new_figure, save_figure, set_theme  # noqa: E402


def _date_col(df: pd.DataFrame) -> pd.Series:
    for col in ("wn_mid_date", "window_mid_date"):
        if col in df.columns:
            return pd.to_datetime(df[col], errors="coerce")
    if "window_idx" in df.columns:
        return df["window_idx"]
    raise KeyError("No window date or index column found for plotting")


def plot_window_coverage(
    window_coverage: pd.DataFrame,
    *,
    out_path: Path = FIGURES_DIR / "window_coverage",
) -> dict[str, Path]:
    """Plot rolling-window sequences, positives, and sequencing proportion."""
    set_theme()
    work = window_coverage.sort_values("window_idx").copy()
    x = _date_col(work)

    fig, ax = new_figure("double", height_in=3.1)
    ax.plot(x, work["wn_no_sequences"], color="#35618f", lw=1.3, label="Sequences")
    ax.plot(
        x,
        work["wn_positive_tests"],
        color="#8a8a8a",
        lw=1.0,
        alpha=0.75,
        label="Positive tests",
    )
    ax.set_ylabel("Rolling-window count")
    ax.set_xlabel("")
    ax.legend(loc="upper left")

    ax2 = ax.twinx()
    ax2.plot(
        x,
        work["wn_prop_sequenced"],
        color="#b0473c",
        lw=1.2,
        label="Sequenced proportion",
    )
    ax2.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax2.set_ylabel("Sequenced")
    ax2.legend(loc="upper right")
    fig.autofmt_xdate()
    return save_figure(fig, out_path, width="double", save_png=True)


def plot_clade_frequencies(
    clade_window_counts: pd.DataFrame,
    *,
    out_path: Path = FIGURES_DIR / "clade_window_frequencies",
) -> dict[str, Path]:
    """Plot stacked clade proportions over rolling windows."""
    work = clade_window_counts.sort_values(["window_idx", "clade"]).copy()
    if work.empty:
        raise ValueError("clade_window_counts is empty")

    pivot = work.pivot_table(
        index="window_idx",
        columns="clade",
        values="proportion",
        aggfunc="sum",
        fill_value=0.0,
    )
    fig, ax = new_figure("double", height_in=3.0)
    ax.stackplot(
        pivot.index, pivot.T.to_numpy(), labels=[str(c) for c in pivot.columns]
    )
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_xlabel("Window index")
    ax.set_ylabel("Sequence proportion")
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        title="Clade",
    )
    return save_figure(fig, out_path, width="double", save_png=True)


def plot_cluster_size_summary(
    cluster_window_summary: pd.DataFrame,
    *,
    out_path: Path = FIGURES_DIR / "cluster_size_summary",
) -> dict[str, Path]:
    """Plot median, 90th percentile, and maximum cluster size by window."""
    work = cluster_window_summary.sort_values("window_idx")
    fig, ax = new_figure("double", height_in=3.0)
    ax.plot(work["window_idx"], work["median_cluster_size"], label="Median", lw=1.2)
    ax.plot(
        work["window_idx"], work["p90_cluster_size"], label="90th percentile", lw=1.2
    )
    ax.plot(work["window_idx"], work["max_cluster_size"], label="Maximum", lw=1.0)
    ax.set_xlabel("Window index")
    ax.set_ylabel("Cluster size")
    ax.legend(loc="upper left")
    return save_figure(fig, out_path, width="double", save_png=True)


def plot_assortativity_over_time(
    assortativity_summary: pd.DataFrame,
    *,
    out_path: Path = FIGURES_DIR / "compatibility_assortativity",
) -> dict[str, Path]:
    """Plot attribute assortativity over windows."""
    if "window_id" not in assortativity_summary.columns:
        raise KeyError("assortativity_summary needs 'window_id'")
    work = assortativity_summary.copy()
    work["window_idx_plot"] = (
        work["window_id"].astype(str).str.extract(r"(\d+)").astype(float)
    )
    work = work.sort_values(["attribute", "window_idx_plot"])

    fig, ax = new_figure("double", height_in=3.1)
    sns.lineplot(
        data=work,
        x="window_idx_plot",
        y="assortativity",
        hue="attribute_label",
        ax=ax,
        lw=1.1,
        errorbar=None,
    )
    ax.axhline(0, color="#777777", lw=0.7, ls=":")
    ax.set_xlabel("Window index")
    ax.set_ylabel("Assortativity")
    ax.legend(title="Attribute", loc="center left", bbox_to_anchor=(1.01, 0.5))
    return save_figure(fig, out_path, width="double", save_png=True)


def plot_degree_assortativity_over_time(
    degree_assortativity: pd.DataFrame,
    *,
    out_path: Path = FIGURES_DIR / "compatibility_degree_assortativity",
) -> dict[str, Path]:
    """Plot degree/strength assortativity diagnostics over windows."""
    if "window_id" not in degree_assortativity.columns:
        raise KeyError("degree_assortativity needs 'window_id'")
    work = degree_assortativity.copy()
    work["window_idx_plot"] = (
        work["window_id"].astype(str).str.extract(r"(\d+)").astype(float)
    )
    metrics = {
        "degree_assortativity": "Degree",
        "weighted_degree_assortativity": "Degree, edge-weighted",
        "strength_assortativity": "Strength, edge-weighted",
    }
    available = [metric for metric in metrics if metric in work.columns]
    if not available:
        raise KeyError("degree_assortativity has no assortativity metric columns")

    long = work.melt(
        id_vars=["window_idx_plot"],
        value_vars=available,
        var_name="metric",
        value_name="assortativity",
    )
    long["metric_label"] = long["metric"].map(metrics)
    long = long.sort_values(["metric_label", "window_idx_plot"])

    fig, ax = new_figure("double", height_in=3.0)
    sns.lineplot(
        data=long,
        x="window_idx_plot",
        y="assortativity",
        hue="metric_label",
        ax=ax,
        lw=1.1,
        errorbar=None,
    )
    ax.axhline(0, color="#777777", lw=0.7, ls=":")
    ax.set_xlabel("Window index")
    ax.set_ylabel("Degree/strength assortativity")
    ax.legend(title="", loc="center left", bbox_to_anchor=(1.01, 0.5))
    return save_figure(fig, out_path, width="double", save_png=True)


def plot_transition_window_summary(
    transition_window_summary: pd.DataFrame,
    *,
    out_path: Path = FIGURES_DIR / "transition_graph_window_summary",
) -> dict[str, Path]:
    """Plot temporal transition graph node and edge counts by window."""
    work = transition_window_summary.sort_values("window_idx")
    fig, ax = new_figure("double", height_in=3.0)
    ax.plot(work["window_idx"], work["n_nodes"], label="Nodes", lw=1.2)
    ax.plot(work["window_idx"], work["n_out_edges"], label="Outgoing edges", lw=1.2)
    ax.set_xlabel("Retained window index")
    ax.set_ylabel("Count")
    ax.legend(loc="upper left")
    return save_figure(fig, out_path, width="double", save_png=True)
