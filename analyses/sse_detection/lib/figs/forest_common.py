"""Shared input and layout helpers for Bayesian forest plots."""

from __future__ import annotations

import pandas as pd

from .common import Paths

OUTLIER_LABELS = ["Orkney", "Western Isles"]
MODEL_OUTCOMES = (
    ("candidate", "Candidate status"),
    ("burst_score", "Burst score"),
    ("burden_score", "Burden score"),
)


def _consolidated_table(paths: Paths, name: str) -> pd.DataFrame:
    table_dir = paths.bayesian_result_dir / "consolidated_tables"
    searched = []
    for suffix, reader in (("parquet", pd.read_parquet), ("csv", pd.read_csv)):
        path = table_dir / f"{name}.{suffix}"
        searched.append(path)
        if path.exists():
            return reader(path)
    locations = ", ".join(str(path) for path in searched)
    raise FileNotFoundError(
        f"Missing consolidated table {name!r}; searched: {locations}"
    )


def _load_mixing_tables(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame]:
    return (
        _consolidated_table(paths, "mixing_logistic_consolidated_results"),
        _consolidated_table(paths, "mixing_linear_consolidated_results"),
    )


def _load_composition_tables(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame]:
    return (
        _consolidated_table(paths, "composition_logistic_consolidated_results"),
        _consolidated_table(paths, "composition_linear_consolidated_results"),
    )


def _finish_legend(fig, axes, *, ncol: int, y: float) -> None:
    handles, labels = axes[0].get_legend_handles_labels()
    for ax in axes:
        legend = ax.get_legend()
        if legend:
            legend.remove()
    fig.legend(
        handles=handles,
        labels=labels,
        loc="lower center",
        bbox_to_anchor=(0.5, y),
        ncol=ncol,
        frameon=False,
    )


def _label_rows(axes, labels: tuple[str, ...], *, x: float = -0.55) -> None:
    for ax, label in zip(axes, labels):
        ax.annotate(
            label,
            xy=(x, 0.5),
            xycoords="axes fraction",
            ha="center",
            va="center",
            rotation=90,
            fontweight="bold",
        )


ADJUSTED_TERM_TYPES = ["continuous_adjuster", "random_intercept"]
MIXING_SCALES = ["null_standardised", "observed"]
VALUE_COLUMNS = ["plot_estimate", "plot_hdi95_low", "plot_hdi95_high"]


def _load_mixing_table(paths: Paths) -> pd.DataFrame:
    return pd.concat(
        [
            _consolidated_table(paths, "mixing_logistic_consolidated_results"),
            _consolidated_table(paths, "mixing_linear_consolidated_results"),
        ],
        ignore_index=True,
    )


def _load_composition_table(paths: Paths) -> pd.DataFrame:
    return pd.concat(
        [
            _consolidated_table(paths, "composition_logistic_consolidated_results"),
            _consolidated_table(paths, "composition_linear_consolidated_results"),
        ],
        ignore_index=True,
    )


def _complete_plot_rows(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out[VALUE_COLUMNS] = out[VALUE_COLUMNS].apply(pd.to_numeric, errors="coerce")
    return out.dropna(subset=["plot_label", "plot_variable_order", *VALUE_COLUMNS])


def _mixing_label_order(df: pd.DataFrame) -> list[str]:
    outcomes = [outcome for outcome, _ in MODEL_OUTCOMES]
    work = df.loc[
        df["outcome"].isin(outcomes)
        & df["plot_term_type"].isin(ADJUSTED_TERM_TYPES)
        & df["plot_scale"].isin(MIXING_SCALES)
    ].copy()
    work = _complete_plot_rows(work)
    labels = (
        work[["plot_label", "plot_variable_order"]]
        .drop_duplicates()
        .assign(
            plot_variable_order=lambda x: pd.to_numeric(
                x["plot_variable_order"], errors="coerce"
            )
        )
        .dropna(subset=["plot_variable_order"])
        .sort_values(["plot_variable_order", "plot_label"])
        .drop_duplicates(subset=["plot_label"])
    )
    return labels["plot_label"].astype(str).tolist()


def _composition_label_order(df: pd.DataFrame) -> list[str]:
    outcomes = [outcome for outcome, _ in MODEL_OUTCOMES]
    work = df.loc[
        df["outcome"].isin(outcomes)
        & df["plot_term_type"].isin(ADJUSTED_TERM_TYPES)
        & ~df["plot_label"].isin(OUTLIER_LABELS)
    ].copy()
    work = _complete_plot_rows(work)
    labels = (
        work[["plot_label", "plot_panel", "plot_variable_order", "plot_level_order"]]
        .drop_duplicates()
        .assign(
            plot_variable_order=lambda x: pd.to_numeric(
                x["plot_variable_order"], errors="coerce"
            )
        )
        .dropna(subset=["plot_variable_order"])
        .sort_values(
            ["plot_panel", "plot_variable_order", "plot_level_order", "plot_label"]
        )
        .drop_duplicates(subset=["plot_label"])
    )
    return labels["plot_label"].astype(str).tolist()
