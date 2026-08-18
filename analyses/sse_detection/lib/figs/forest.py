"""Forest plots for Bayesian regression summaries."""

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

PROJECT_ROOT = Path(__file__).resolve().parents[4]

if __package__ in {None, ""} and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


MIXING_FEATURE_LABELS: dict[str, str] = {
    "sex_entropy": "Sex",
    "age_entropy": "Age",
    "simd_entropy": "SIMD",
    "urban_rural_entropy": "Urban/rural",
    "health_board_entropy": "Health board",
}
DEFAULT_MIXING_FEATURE_ORDER: tuple[str, ...] = (
    "sex_entropy",
    "age_entropy",
    "simd_entropy",
    "urban_rural_entropy",
    "health_board_entropy",
)
RANDOM_EFFECT_PREFIX = "1|"


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
    label_order: Sequence[str] | None = None,
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
    if label_order is None:
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
    else:
        labels = list(dict.fromkeys(str(label) for label in label_order))
        extras = [
            label
            for label in sorted(df["plot_label"].dropna().astype(str).unique())
            if label not in labels
        ]
        tick_df = pd.DataFrame({"plot_label": [*labels, *extras]})

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
    label_order: Sequence[str] | None = None,
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
        "plot_estimate",
        "plot_hdi95_low",
        "plot_hdi95_high",
        "plot_reference_value",
        label_col,
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

    value_cols = ["plot_estimate", "plot_hdi95_low", "plot_hdi95_high"]
    for col in value_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(
        subset=[label_col, "plot_variable_order", *value_cols]
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

    if label_order is None:
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
    else:
        labels = list(dict.fromkeys(str(label) for label in label_order))
        extras = [
            label
            for label in sorted(df[label_col].dropna().astype(str).unique())
            if label not in labels
        ]
        tick_df = pd.DataFrame({label_col: [*labels, *extras]})

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
        lower_err = sub["plot_estimate"] - sub["plot_hdi95_low"]
        upper_err = sub["plot_hdi95_high"] - sub["plot_estimate"]

        # Tiny negatives can still happen from floating point only.
        lower_err = lower_err.clip(lower=0)
        upper_err = upper_err.clip(lower=0)

        ax.errorbar(
            x=sub["plot_estimate"],
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


def _collect_mixing_forest_rows(
    model_config: Any | str | Path,
    *,
    family: str,
    outcome: str,
    feature_order: Sequence[str],
) -> tuple[pd.DataFrame, list[str]]:
    """Collect legacy mixing fixed-effect rows from per-model summary files."""
    model_sets = (
        ("observed", "primary", "observed_primary"),
        ("observed", "expanded", "observed_expanded"),
        ("null", "primary", "null_primary"),
        ("null", "expanded", "null_expanded"),
    )
    rows: list[pd.DataFrame] = []
    missing: list[str] = []

    for scale, model_kind, model_set in model_sets:
        summary_path = _summary_path(
            model_config,
            domain="mixing",
            model_set=model_set,
            family=family,
            outcome=outcome,
        )
        if summary_path is None or not summary_path.exists():
            missing.append(model_set)
            continue
        summary = pd.read_csv(summary_path)
        fixed = _tidy_mixing_summary(
            summary,
            family=family,
            scale=scale,
            model_kind=model_kind,
            model_set=model_set,
            feature_order=feature_order,
            summary_path=summary_path,
        )
        rows.append(fixed)

    if not rows:
        return pd.DataFrame(), missing

    out = pd.concat(rows, ignore_index=True)
    out["row_id"] = out["feature"]
    out["feature"] = pd.Categorical(
        out["feature"],
        categories=list(feature_order),
        ordered=True,
    )
    out["model_kind"] = _model_kind_categorical(out["model_kind"].tolist())
    out["scale"] = pd.Categorical(out["scale"], categories=["observed", "null"])
    return out.sort_values(["scale", "feature", "model_kind"]), missing


def _tidy_mixing_summary(
    summary: pd.DataFrame,
    *,
    family: str,
    scale: str,
    model_kind: str,
    model_set: str,
    feature_order: Sequence[str],
    summary_path: Path,
) -> pd.DataFrame:
    _require_summary_columns(summary, family=family, path=summary_path)
    out = _fixed_effect_rows(summary)
    out["feature"] = out["parameter"].map(_mixing_feature_stem)
    out = out.loc[out["feature"].isin(feature_order)].copy()
    out["scale"] = scale
    out["model_kind"] = model_kind
    out["model_set"] = model_set
    out["summary_path"] = str(summary_path)
    return _select_plot_columns(out, family=family)


def _summary_path(
    model_config: Any | str | Path,
    *,
    domain: str,
    model_set: str,
    family: str,
    outcome: str,
) -> Path | None:
    if isinstance(model_config, (str, Path)):
        result_dir = Path(model_config)
        parts = [result_dir, Path(domain)]
        if family == "linear":
            parts.append(Path(_slug(outcome)))
        parts.append(Path(_slug(model_set)))
        return Path(*parts) / "summary.csv"

    try:
        frame = model_config.select(
            domain=domain,
            outcome=outcome,
            model_set=model_set,
        )
    except (AttributeError, KeyError, ValueError):
        return None
    return Path(frame.output_dir) / "summary.csv"


def _fixed_effect_rows(summary: pd.DataFrame) -> pd.DataFrame:
    return summary.loc[
        ~summary["parameter"].eq("Intercept")
        & ~summary["parameter"].str.startswith(RANDOM_EFFECT_PREFIX, na=False)
    ].copy()


def _select_plot_columns(
    out: pd.DataFrame,
    *,
    family: str,
) -> pd.DataFrame:
    keep = [
        "scale",
        "model_kind",
        "model_set",
        "feature",
        "parameter",
        "mean",
        "hdi95_lb",
        "hdi95_ub",
        "OR_mean",
        "OR_hdi95_lb",
        "OR_hdi95_ub",
        "P(OR > 1 | data)",
        "P(OR < 1 | data)",
        "P(beta > 0 | data)",
        "P(beta < 0 | data)",
        "summary_path",
    ]
    if family == "linear":
        keep = [col for col in keep if not col.startswith("OR_")]
    return out.loc[:, [col for col in keep if col in out.columns]]


def _require_summary_columns(summary: pd.DataFrame, *, family: str, path: Path) -> None:
    required = {"parameter", "mean", "hdi95_lb", "hdi95_ub"}
    if family == "logistic":
        required |= {"OR_mean", "OR_hdi95_lb", "OR_hdi95_ub"}
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")


def _mixing_feature_stem(parameter: str) -> str:
    for suffix in ("_obs_x10", "_z"):
        if str(parameter).endswith(suffix):
            return str(parameter)[: -len(suffix)]
    return str(parameter)


def _model_kind_categorical(values: list[object]) -> pd.Categorical:
    return pd.Categorical(values, categories=["primary", "expanded"], ordered=True)


def _slug(value: str) -> str:
    return str(value).replace("/", "_").replace(" ", "_")
