"""Forest plots for Bayesian regression summaries."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

PROJECT_ROOT = Path(__file__).resolve().parents[4]

if __package__ in {None, ""} and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def plot_mixing_forest(
    ax: Axes,
    df: pd.DataFrame,
    *,
    outcome: str | None = None,
    plot_scale: str | list | None = "null_standardised",
    term_type: str | list | None = "mixing_entropy",
    model_kinds: str | tuple[str, str] | None = None,
    offset: float = 0.16,
    exclude: list[str] | None = None,
):
    """
    Forest plot for consolidated mixing results.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    df : pandas.DataFrame
        One of the mixing consolidated tables.
    outcome : str or None
        Required when df contains multiple outcomes, e.g. mixing_linear.
    plot_scale : {"observed", "null_standardised", "null", None, list}
        If None, keeps all available scales and plots each scale/model combination.
        Examples: "mixing_entropy", "continuous_adjuster",
        "random_intercept", "random_effect_sd", "intercept", "residual_sd".
        If None, keeps all term types.
    model_kinds : {"primary", "expanded"}, tuple, or None
        If None, keeps whatever model kinds remain after filtering.
    offset : float
        Offset between model series.
    exclude : list of str or None
        List of plot_label values to exclude from the plot.
    """

    df = df.copy()

    if exclude is not None:
        df = df[~df["plot_label"].isin(exclude)]

    required = {
        "outcome",
        "model_set",
        "plot_model_kind",
        "plot_scale",
        "plot_term_type",
        "plot_label",
        "plot_variable_order",
        "plot_estimate",
        "plot_hdi95_low",
        "plot_hdi95_high",
        "plot_reference_value",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required consolidated-table columns: {missing}")

    if df.empty:
        raise ValueError("Input dataframe is empty.")

    def _as_list(value):
        if value is None:
            return None
        if isinstance(value, str):
            return [value]
        return list(value)

    def _normalise_scale(value):
        if value == "null":
            return "null_standardised"
        return value

    # Outcome filtering
    available_outcomes = sorted(df["outcome"].dropna().astype(str).unique())
    if outcome is None:
        if len(available_outcomes) > 1:
            raise ValueError(
                "Dataframe contains multiple outcomes. Specify outcome=... "
                f"Available outcomes: {available_outcomes}"
            )
    else:
        df = df.loc[df["outcome"].astype(str).eq(str(outcome))].copy()

    # Term filtering
    term_types = _as_list(term_type)
    if term_types is not None:
        df = df.loc[df["plot_term_type"].isin(term_types)].copy()

    # Scale filtering
    scales = _as_list(plot_scale)
    if scales is not None:
        scales = [_normalise_scale(scale) for scale in scales]
        df = df.loc[df["plot_scale"].isin(scales)].copy()

    # Model-kind filtering
    kinds = _as_list(model_kinds)
    if kinds is not None:
        df = df.loc[df["plot_model_kind"].isin(kinds)].copy()

    if df.empty:
        raise ValueError("No rows remain after filtering.")

    # Drop rows with missing plot values.
    value_cols = ["plot_estimate", "plot_hdi95_low", "plot_hdi95_high"]
    df[value_cols] = df[value_cols].apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=["plot_label", "plot_variable_order", *value_cols]).copy()

    if df.empty:
        raise ValueError("No rows remain with complete plotting values.")

    # One plotted series is one scale x model-kind combination.
    series_cols = ["plot_scale", "plot_model_kind"]
    series = (
        df[series_cols]
        .drop_duplicates()
        .sort_values(series_cols)
        .to_dict(orient="records")
    )

    if not series:
        raise ValueError("No model series available after filtering.")

    # Each series must have at most one row per plot label.
    for spec in series:
        mask = df["plot_scale"].eq(spec["plot_scale"]) & df["plot_model_kind"].eq(
            spec["plot_model_kind"]
        )
        sub = df.loc[mask]
        dup = sub["plot_label"].duplicated(keep=False)
        if dup.any():
            bad = sorted(sub.loc[dup, "plot_label"].astype(str).unique())
            raise ValueError(
                "Duplicate rows within a plotted series. "
                f"Series={spec}, duplicated plot_label values={bad}. "
                "Filter more specifically, usually by outcome, plot_scale, or term_type."
            )

    # plot_variable_order is a group order, not a unique y-position.
    tick_df = (
        df[["plot_label", "plot_variable_order"]]
        .drop_duplicates()
        .assign(
            plot_variable_order=lambda x: pd.to_numeric(
                x["plot_variable_order"], errors="coerce"
            )
        )
        .dropna(subset=["plot_variable_order"])
        .sort_values(["plot_variable_order", "plot_label"])
        .reset_index(drop=True)
    )

    if tick_df.empty:
        raise ValueError("No valid plot labels/order values remain.")

    # Put first sorted label at top.
    n_rows = len(tick_df)
    tick_df["_y"] = list(range(n_rows - 1, -1, -1))
    label_to_y = dict(zip(tick_df["plot_label"], tick_df["_y"]))

    n_series = len(series)
    if n_series == 1:
        shifts = [0.0]
    else:
        midpoint = (n_series - 1) / 2
        shifts = [(midpoint - i) * offset for i in range(n_series)]

    colors = [
        "royalblue",
        "darkorange",
        "seagreen",
        "firebrick",
        "mediumpurple",
        "sienna",
    ]
    markers = ["s", "o", "^", "D", "v", "P"]

    def _scale_label(scale):
        return {
            "observed": "Observed",
            "null_standardised": "Null-standardised",
        }.get(str(scale), str(scale))

    def _kind_label(kind):
        return str(kind).title()

    unique_scales = df["plot_scale"].dropna().astype(str).nunique()
    unique_kinds = df["plot_model_kind"].dropna().astype(str).nunique()

    def _series_label(spec):
        scale = spec["plot_scale"]
        kind = spec["plot_model_kind"]
        if unique_scales > 1 and unique_kinds > 1:
            return f"{_scale_label(scale)} {_kind_label(kind)}"
        if unique_scales > 1:
            return _scale_label(scale)
        if unique_kinds > 1:
            return _kind_label(kind)
        return f"{_scale_label(scale)} {_kind_label(kind)}"

    for i, spec in enumerate(series):
        mask = df["plot_scale"].eq(spec["plot_scale"]) & df["plot_model_kind"].eq(
            spec["plot_model_kind"]
        )
        sub = df.loc[mask].copy()
        if sub.empty:
            continue

        y = sub["plot_label"].map(label_to_y).astype(float) + shifts[i]
        x = sub["plot_estimate"]
        lower_err = sub["plot_estimate"] - sub["plot_hdi95_low"]
        upper_err = sub["plot_hdi95_high"] - sub["plot_estimate"]

        bad_interval = (lower_err < 0) | (upper_err < 0)
        if bad_interval.any():
            bad = sub.loc[bad_interval, "plot_label"].tolist()
            raise ValueError(f"Invalid intervals for labels: {bad}")

        ax.errorbar(
            x=x,
            y=y,
            xerr=[lower_err, upper_err],
            fmt=markers[i % len(markers)],
            color=colors[i % len(colors)],
            ecolor=colors[i % len(colors)],
            linestyle="none",
            capsize=3,
            ms=5,
            label=_series_label(spec),
        )

    reference_values = df["plot_reference_value"].dropna().unique()

    reference = float(reference_values[0]) if len(reference_values) else 0.0

    ax.axvline(reference, color="red", linestyle="--", alpha=0.7, linewidth=1)
    ax.set_yticks(tick_df["_y"])
    ax.set_yticklabels(tick_df["plot_label"])
    ax.set_ylim(-0.5 - offset, n_rows - 0.5 + offset)

    if reference == 1.0:
        # ax.set_xscale("log")
        # set_readable_or_ticks(ax)
        ax.set_xlabel("Odds ratio")
    else:
        ax.set_xlabel("Coefficient")
    ax.legend(loc="best")
    return ax


def plot_composition_forest(
    ax: Axes,
    df: pd.DataFrame,
    *,
    outcome: str | None = None,
    panel: str | list[str] | None = None,
    variables: str | list[str] | None = None,
    term_type: str | list[str] | None = "categorical_contrast",
    model_kinds: str | list[str] | None = None,
    label_col: str = "plot_label",
    offset: float = 0.16,
    exclude: list[str] | None = None,
):
    df = df.copy()

    if exclude is not None:
        df = df[~df["plot_label"].isin(exclude)]

    required = {
        "outcome",
        "model_set",
        "plot_model_kind",
        "plot_term_type",
        "plot_variable",
        "plot_panel",
        "plot_variable_order",
        "plot_reference_value",
        label_col,
        "mean",
        "hdi95_lb",
        "hdi95_ub",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if df.empty:
        raise ValueError("Input dataframe is empty.")

    def _as_list(value):
        if value is None:
            return None
        if isinstance(value, str):
            return [value]
        return list(value)

    available_outcomes = sorted(df["outcome"].dropna().astype(str).unique())
    if outcome is None:
        if len(available_outcomes) > 1:
            raise ValueError(
                "Dataframe contains multiple outcomes. Specify outcome=... "
                f"Available outcomes: {available_outcomes}"
            )
    else:
        df = df.loc[df["outcome"].astype(str).eq(str(outcome))].copy()

    term_types = _as_list(term_type)
    if term_types is not None:
        df = df.loc[df["plot_term_type"].isin(term_types)].copy()

    panels = _as_list(panel)
    if panels is not None:
        df = df.loc[df["plot_panel"].isin(panels)].copy()

    variables = _as_list(variables)
    if variables is not None:
        df = df.loc[df["plot_variable"].isin(variables)].copy()

    kinds = _as_list(model_kinds)
    if kinds is not None:
        df = df.loc[df["plot_model_kind"].isin(kinds)].copy()

    if df.empty:
        raise ValueError("No rows remain after filtering.")

    reference_values = df["plot_reference_value"].dropna().unique()
    reference = float(reference_values[0]) if len(reference_values) else 0.0
    logistic_scale = np.isclose(reference, 1.0)

    # Use log-scale posterior columns for logistic OR plotting to avoid
    # rounded OR_hdi columns where OR_mean can fall outside rounded HDI bounds.
    for col in ["mean", "hdi95_lb", "hdi95_ub"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if logistic_scale:
        df["_estimate"] = np.exp(df["mean"])
        df["_low"] = np.exp(df["hdi95_lb"])
        df["_high"] = np.exp(df["hdi95_ub"])
    else:
        df["_estimate"] = df["mean"]
        df["_low"] = df["hdi95_lb"]
        df["_high"] = df["hdi95_ub"]

    df = df.dropna(
        subset=[label_col, "plot_variable_order", "_estimate", "_low", "_high"]
    ).copy()

    if df.empty:
        raise ValueError("No rows remain with complete plotting values.")

    series = (
        df[["plot_model_kind"]]
        .drop_duplicates()
        .sort_values("plot_model_kind")
        .to_dict(orient="records")
    )

    for spec in series:
        sub = df.loc[df["plot_model_kind"].eq(spec["plot_model_kind"])]
        dup = sub[label_col].duplicated(keep=False)
        if dup.any():
            bad = sorted(sub.loc[dup, label_col].astype(str).unique())
            raise ValueError(
                f"Duplicate labels within {spec['plot_model_kind']}: {bad}. "
                "Filter more specifically, or use label_col='plot_contrast_label'."
            )

    sort_cols = ["plot_panel", "plot_variable_order"]
    if "plot_level_order" in df.columns:
        sort_cols.append("plot_level_order")
    sort_cols.append(label_col)

    tick_cols = [label_col, "plot_panel", "plot_variable_order"]
    if "plot_level_order" in df.columns:
        tick_cols.append("plot_level_order")

    tick_df = (
        df[tick_cols]
        .drop_duplicates()
        .assign(
            plot_variable_order=lambda x: pd.to_numeric(
                x["plot_variable_order"], errors="coerce"
            )
        )
        .dropna(subset=["plot_variable_order"])
        .sort_values(sort_cols)
        .reset_index(drop=True)
    )

    n_rows = len(tick_df)
    tick_df["_y"] = list(range(n_rows - 1, -1, -1))
    label_to_y = dict(zip(tick_df[label_col], tick_df["_y"]))

    n_series = len(series)
    shifts = (
        [0.0]
        if n_series == 1
        else [((n_series - 1) / 2 - i) * offset for i in range(n_series)]
    )

    colors = ["royalblue", "darkorange", "seagreen", "firebrick"]
    markers = ["s", "o", "^", "D"]

    for i, spec in enumerate(series):
        sub = df.loc[df["plot_model_kind"].eq(spec["plot_model_kind"])].copy()
        if sub.empty:
            continue

        y = sub[label_col].map(label_to_y).astype(float) + shifts[i]
        lower_err = sub["_estimate"] - sub["_low"]
        upper_err = sub["_high"] - sub["_estimate"]

        # Tiny negatives can still happen from floating point only.
        lower_err = lower_err.clip(lower=0)
        upper_err = upper_err.clip(lower=0)

        ax.errorbar(
            x=sub["_estimate"],
            y=y,
            xerr=[lower_err, upper_err],
            fmt=markers[i % len(markers)],
            color=colors[i % len(colors)],
            ecolor=colors[i % len(colors)],
            linestyle="none",
            capsize=3,
            ms=5,
            label=str(spec["plot_model_kind"]).title(),
        )

    ax.axvline(reference, color="red", linestyle="--", alpha=0.7, linewidth=1)

    ax.set_yticks(tick_df["_y"])
    ax.set_yticklabels(tick_df[label_col])
    ax.set_ylim(-0.5 - offset, n_rows - 0.5 + offset)

    if logistic_scale:
        # ax.set_xscale("log")
        # set_readable_or_ticks(ax)
        ax.set_xlabel("Odds ratio")
    else:
        ax.set_xlabel("Coefficient")

    ax.legend(loc="best")
    return ax
