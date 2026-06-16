"""Forest plots for Bayesian regression summaries.

The functions in this module read the ``summary.csv`` files written by the
Bayesian regression notebooks and turn fixed-effect summaries into
publication-style forest plots.  Two result shapes are supported:

* cluster-level mixing models, with observed and null-standardised entropy panels;
* sequence-level composition models, with cleaned categorical contrasts.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence
import re
import warnings
import sys

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]

if __package__ in {None, ""} and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import style  # noqa: E402


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

COMPOSITION_VARIABLE_LABELS: dict[str, str] = {
    "sex": "Sex",
    "age_band": "Age band",
    "dz_simd_quintile": "SIMD",
    "urban_rural_class": "Urban/rural",
    "dz_urban_rural_class": "Urban/rural",
    "health_board": "Health board",
    "dz_health_board": "Health board",
}

COMPOSITION_REFERENCE_LABELS: dict[str, str] = {
    "sex": "Male",
    "age_band": "20-24",
    "dz_simd_quintile": "1",
    "urban_rural_class": "Large Urban Areas",
    "dz_urban_rural_class": "Large Urban Areas",
    "health_board": "Greater Glasgow and Clyde",
    "dz_health_board": "Greater Glasgow and Clyde",
}

DEFAULT_COMPOSITION_VARIABLE_ORDER: tuple[str, ...] = (
    "sex",
    "age_band",
    "dz_simd_quintile",
    "dz_urban_rural_class",
    "dz_health_board",
)

COMPOSITION_PANEL_GROUPS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    (
        "sociodemographic",
        "Sociodemographic",
        ("sex", "age_band", "dz_simd_quintile"),
    ),
    (
        "geographic",
        "Geographic",
        (
            "urban_rural_class",
            "dz_urban_rural_class",
            "health_board",
            "dz_health_board",
        ),
    ),
)

DEFAULT_MODEL_LABELS: dict[str, str] = {
    "primary": "Primary",
    "expanded": "Expanded",
}

DEFAULT_COLORS: dict[str, str] = {
    "primary": "#4C78A8",
    "expanded": "#E6A23C",
}

RANDOM_EFFECT_PREFIX = "1|"
_CATEGORICAL_RE = re.compile(r"^C\((?P<inside>.+?)\)\[(?P<level>.+)\]$")
_BARE_CATEGORICAL_RE = re.compile(
    r"^(?P<variable>[A-Za-z_][A-Za-z0-9_]*)\[(?P<level>.+)\]$"
)
_TREATMENT_VARIABLE_RE = re.compile(
    r"^(?P<variable>[^,]+),\s*Treatment\(reference=(?P<reference>.+)\)$"
)


@dataclass(frozen=True)
class ForestPlotResult:
    """Container returned by forest plotting helpers."""

    fig: Figure
    axes: np.ndarray
    plot_data: pd.DataFrame
    missing_model_sets: tuple[str, ...]


def plot_mixing_primary_expanded_forest(
    model_config: Any | str | Path,
    *,
    family: str = "logistic",
    outcome: str = "candidate",
    feature_order: Sequence[str] = DEFAULT_MIXING_FEATURE_ORDER,
    feature_labels: Mapping[str, str] = MIXING_FEATURE_LABELS,
    model_labels: Mapping[str, str] = DEFAULT_MODEL_LABELS,
    width: style.WIDTHS = "double",
    height_in: float = 3.4,
    context: style.CONTEXTS = "paper",
    font_scale: float = 1.0,
    colors: Mapping[str, str] | None = None,
    point_size: float = 5.0,
    interval_lw: float = 1.4,
    dodge: float = 0.13,
    xlim: tuple[float, float] | None = None,
    title: str | None = None,
) -> ForestPlotResult:
    """Plot primary and expanded cluster-level mixing estimates.

    Logistic models are shown as odds ratios on a log scale. Linear models are
    shown as posterior coefficients on the native outcome scale.
    """
    family = _normalise_family(family)
    rows, missing = _collect_mixing_forest_rows(
        model_config,
        family=family,
        outcome=outcome,
        feature_order=feature_order,
    )
    if rows.empty:
        raise FileNotFoundError("No fitted mixing summary.csv files were found.")

    colors = dict(colors or DEFAULT_COLORS)
    fig, axes = style.new_figure(
        nrows=1,
        ncols=2,
        sharey=True,
        width=width,
        height_in=height_in,
        layout="constrained",
        context=context,
        font_scale=font_scale,
    )

    panel_specs = (
        (
            "observed",
            "Observed entropy",
            "OR per 0.1-bit increase"
            if family == "logistic"
            else "Coefficient per 0.1-bit increase",
        ),
        (
            "null",
            "Null-standardised entropy",
            "OR per 1 SD increase"
            if family == "logistic"
            else "Coefficient per 1 SD increase",
        ),
    )
    y_lookup = _y_lookup(feature_order)
    x_cols = _estimate_columns(family)

    for ax, (scale, panel_title, x_label) in zip(axes, panel_specs):
        panel = rows.loc[rows["scale"].eq(scale)].copy()
        _draw_paired_forest_panel(
            ax,
            panel,
            y_lookup=y_lookup,
            row_order=feature_order,
            row_labels=feature_labels,
            model_labels=model_labels,
            colors=colors,
            x_cols=x_cols,
            reference=1.0 if family == "logistic" else 0.0,
            point_size=point_size,
            interval_lw=interval_lw,
            dodge=dodge,
        )
        if family == "logistic":
            ax.set_xscale("log")
            _set_readable_or_ticks(ax)
        ax.set_xlabel(x_label)
        ax.set_title(panel_title)
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)
        ax.tick_params(axis="y", length=0)
        if xlim is not None:
            ax.set_xlim(*xlim)

    axes[1].tick_params(labelleft=False)
    _add_model_legend(fig, axes[0], colors, model_labels, point_size)
    _finish_forest_figure(fig, axes, title=title)
    _warn_missing(missing, "mixing")
    plt.close(fig)
    return ForestPlotResult(fig, np.asarray(axes), rows, tuple(missing))


def plot_logistic_mixing_primary_expanded_forest(
    model_config: Any | str | Path,
    **kwargs: Any,
) -> ForestPlotResult:
    """Backward-compatible wrapper for logistic mixing forests."""
    return plot_mixing_primary_expanded_forest(
        model_config,
        family="logistic",
        outcome=kwargs.pop("outcome", "candidate"),
        **kwargs,
    )


def plot_composition_primary_expanded_forest(
    model_config: Any | str | Path,
    *,
    family: str = "logistic",
    outcome: str = "candidate",
    variable_order: Sequence[str] = DEFAULT_COMPOSITION_VARIABLE_ORDER,
    variable_labels: Mapping[str, str] = COMPOSITION_VARIABLE_LABELS,
    reference_labels: Mapping[str, str] = COMPOSITION_REFERENCE_LABELS,
    model_labels: Mapping[str, str] = DEFAULT_MODEL_LABELS,
    width: style.WIDTHS = "double",
    height_per_contrast: float = 0.33,
    min_height_in: float = 3.4,
    context: style.CONTEXTS = "paper",
    font_scale: float = 1.0,
    colors: Mapping[str, str] | None = None,
    point_size: float = 5.0,
    interval_lw: float = 1.4,
    dodge: float = 0.13,
    xlim: tuple[float, float] | None = None,
    title: str | None = None,
) -> ForestPlotResult:
    """Plot primary and expanded sequence-level composition contrasts.

    Categorical Bambi/formulae parameter names are parsed into labels such as
    ``Age band: 30-34 vs 20-24`` before plotting.
    """
    family = _normalise_family(family)
    rows, missing = _collect_composition_forest_rows(
        model_config,
        family=family,
        outcome=outcome,
        variable_order=variable_order,
        variable_labels=variable_labels,
        reference_labels=reference_labels,
    )
    if rows.empty:
        raise FileNotFoundError("No fitted composition summary.csv files were found.")

    rows = _assign_composition_panels(rows)
    panel_rows = [
        (
            panel_id,
            panel_title,
            rows.loc[rows["panel"].eq(panel_id)].copy(),
        )
        for panel_id, panel_title, _ in COMPOSITION_PANEL_GROUPS
    ]
    panel_rows = [spec for spec in panel_rows if not spec[2].empty]
    if not panel_rows:
        raise FileNotFoundError("No fitted composition summaries matched the panels.")

    max_panel_rows = max(
        panel["contrast_id"].drop_duplicates().size for _, _, panel in panel_rows
    )
    height_in = max(min_height_in, height_per_contrast * max_panel_rows + 1.2)
    colors = dict(colors or DEFAULT_COLORS)
    fig, axes = style.new_figure(
        width=width,
        height_in=height_in,
        nrows=1,
        ncols=len(panel_rows),
        layout="constrained",
        context=context,
        font_scale=font_scale,
        gridspec_kw={"wspace": 0.18},
    )
    axes = np.atleast_1d(axes)

    x_cols = _estimate_columns(family)
    reference_value = 1.0 if family == "logistic" else 0.0
    x_label = (
        "Odds ratio vs reference category"
        if family == "logistic"
        else "Coefficient vs reference category"
    )
    for ax, (_, panel_title, panel) in zip(axes, panel_rows):
        row_order = panel["contrast_id"].drop_duplicates().tolist()
        row_labels = dict(zip(panel["contrast_id"], panel["display_label"]))
        _draw_paired_forest_panel(
            ax,
            panel,
            y_lookup=_y_lookup(row_order),
            row_order=row_order,
            row_labels=row_labels,
            model_labels=model_labels,
            colors=colors,
            x_cols=x_cols,
            reference=reference_value,
            point_size=point_size,
            interval_lw=interval_lw,
            dodge=dodge,
        )
        if xlim is None:
            _set_panel_xlim(ax, panel, x_cols=x_cols, reference=reference_value)
        else:
            ax.set_xlim(*xlim)
        if family == "logistic":
            ax.set_xscale("log")
            _set_readable_or_ticks(ax)
        ax.set_xlabel(x_label)
        ax.set_title(panel_title)
        ax.grid(axis="x", color="#E6E6E6", lw=0.6)
        ax.tick_params(axis="y", length=0, pad=2)
        _add_composition_group_guides(ax, panel, row_order)

    _add_model_legend(fig, axes[0], colors, model_labels, point_size)
    _finish_forest_figure(fig, axes, title=title)
    _warn_missing(missing, "composition")
    plt.close(fig)
    return ForestPlotResult(fig, axes, rows, tuple(missing))


def plot_logistic_composition_primary_expanded_forest(
    model_config: Any | str | Path,
    **kwargs: Any,
) -> ForestPlotResult:
    """Convenience wrapper for logistic sequence-composition forests."""
    return plot_composition_primary_expanded_forest(
        model_config,
        family="logistic",
        outcome=kwargs.pop("outcome", "candidate"),
        **kwargs,
    )


def _draw_paired_forest_panel(
    ax: Axes,
    panel: pd.DataFrame,
    *,
    y_lookup: Mapping[str, int],
    row_order: Sequence[str],
    row_labels: Mapping[str, str],
    model_labels: Mapping[str, str],
    colors: Mapping[str, str],
    x_cols: tuple[str, str, str],
    reference: float,
    point_size: float,
    interval_lw: float,
    dodge: float,
) -> None:
    ax.set_ylim(-0.6, len(row_order) - 0.4)
    ax.set_yticks(
        [y_lookup[row] for row in row_order],
        [row_labels.get(row, row) for row in row_order],
    )
    ax.axvline(reference, color="#555555", lw=0.8, ls="--", zorder=0)

    if panel.empty:
        ax.text(
            0.5,
            0.5,
            "No fitted summaries",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="#666666",
        )
        return

    estimate_col, lower_col, upper_col = x_cols
    offsets = {"primary": dodge, "expanded": -dodge}
    for model_kind in ("primary", "expanded"):
        data = panel.loc[panel["model_kind"].eq(model_kind)]
        if data.empty:
            continue
        y = data["row_id"].map(y_lookup).astype(float) + offsets[model_kind]
        color = colors[model_kind]
        ax.hlines(
            y=y,
            xmin=data[lower_col],
            xmax=data[upper_col],
            color=color,
            lw=interval_lw,
            alpha=0.95,
            zorder=2,
        )
        ax.plot(
            data[estimate_col],
            y,
            "o",
            color=color,
            markersize=point_size,
            label=model_labels.get(model_kind, model_kind.title()),
            zorder=3,
        )


def _assign_composition_panels(rows: pd.DataFrame) -> pd.DataFrame:
    out = rows.copy()
    panel_lookup = {
        variable: panel_id
        for panel_id, _, variables in COMPOSITION_PANEL_GROUPS
        for variable in variables
    }
    out["panel"] = out["variable"].map(panel_lookup).fillna("other")
    out["display_label"] = out.apply(
        lambda row: _composition_display_label(row["variable"], row["level"]),
        axis=1,
    )
    out["reference_display_label"] = out.apply(
        lambda row: _composition_display_label(row["variable"], row["reference"]),
        axis=1,
    )
    return out


def _add_composition_group_guides(
    ax: Axes,
    panel: pd.DataFrame,
    row_order: Sequence[str],
) -> None:
    if not row_order:
        return

    y_lookup = _y_lookup(row_order)
    guide_x = -0.45
    groups = panel.drop_duplicates("contrast_id").groupby("variable", sort=False)
    for _, group in groups:
        y_values = [y_lookup[row] for row in group["contrast_id"]]
        ymin, ymax = min(y_values), max(y_values)
        reference = str(group["reference_display_label"].iloc[0])
        ax.text(
            guide_x,
            (ymin + ymax) / 2,
            f"Ref: {reference}",
            transform=ax.get_yaxis_transform(),
            rotation=90,
            ha="center",
            va="center",
            color="#555555",
            fontsize="small",
            linespacing=0.9,
            clip_on=False,
        )
        if ymin > 0:
            ax.axhline(
                ymin - 0.5,
                color="#EEEEEE",
                lw=0.7,
                zorder=0,
            )


def _set_panel_xlim(
    ax: Axes,
    panel: pd.DataFrame,
    *,
    x_cols: tuple[str, str, str],
    reference: float,
) -> None:
    _, lower_col, upper_col = x_cols
    bounds = pd.concat(
        [
            panel[lower_col],
            panel[upper_col],
            pd.Series([reference], dtype=float),
        ],
        ignore_index=True,
    )
    bounds = pd.to_numeric(bounds, errors="coerce")
    bounds = bounds[np.isfinite(bounds)]
    if bounds.empty:
        return

    xmin = float(bounds.min())
    xmax = float(bounds.max())
    if reference > 0 and xmin > 0:
        log_min, log_max = np.log([xmin, xmax])
        pad = max((log_max - log_min) * 0.08, 0.04)
        ax.set_xlim(float(np.exp(log_min - pad)), float(np.exp(log_max + pad)))
        return

    span = xmax - xmin
    pad = max(span * 0.08, 0.05)
    ax.set_xlim(xmin - pad, xmax + pad)


def _collect_mixing_forest_rows(
    model_config: Any | str | Path,
    *,
    family: str,
    outcome: str,
    feature_order: Sequence[str],
) -> tuple[pd.DataFrame, list[str]]:
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
        out["feature"], categories=list(feature_order), ordered=True
    )
    out["model_kind"] = _model_kind_categorical(out["model_kind"])
    out["scale"] = pd.Categorical(out["scale"], categories=["observed", "null"])
    return out.sort_values(["scale", "feature", "model_kind"]), missing


def _collect_composition_forest_rows(
    model_config: Any | str | Path,
    *,
    family: str,
    outcome: str,
    variable_order: Sequence[str],
    variable_labels: Mapping[str, str],
    reference_labels: Mapping[str, str],
) -> tuple[pd.DataFrame, list[str]]:
    rows: list[pd.DataFrame] = []
    missing: list[str] = []

    for model_kind, model_set in (("primary", "primary"), ("expanded", "expanded")):
        summary_path = _summary_path(
            model_config,
            domain="composition",
            model_set=model_set,
            family=family,
            outcome=outcome,
        )
        if summary_path is None or not summary_path.exists():
            missing.append(model_set)
            continue
        summary = pd.read_csv(summary_path)
        rows.append(
            _tidy_composition_summary(
                summary,
                family=family,
                model_kind=model_kind,
                model_set=model_set,
                variable_order=variable_order,
                variable_labels=variable_labels,
                reference_labels=reference_labels,
                summary_path=summary_path,
            )
        )

    if not rows:
        return pd.DataFrame(), missing

    out = pd.concat(rows, ignore_index=True)
    order_map = {variable: i for i, variable in enumerate(variable_order)}
    out["_variable_order"] = out["variable"].map(order_map).fillna(len(order_map))
    out["_level_sort"] = out["level"].map(_level_sort_key)
    out = out.sort_values(
        ["_variable_order", "variable", "_level_sort", "level", "model_kind"]
    ).drop(columns=["_variable_order", "_level_sort"])
    out["contrast_id"] = out["variable"] + "=" + out["level"].astype(str)
    out["row_id"] = out["contrast_id"]
    out["model_kind"] = _model_kind_categorical(out["model_kind"])
    return out, missing


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


def _tidy_composition_summary(
    summary: pd.DataFrame,
    *,
    family: str,
    model_kind: str,
    model_set: str,
    variable_order: Sequence[str],
    variable_labels: Mapping[str, str],
    reference_labels: Mapping[str, str],
    summary_path: Path,
) -> pd.DataFrame:
    _require_summary_columns(summary, family=family, path=summary_path)
    out = _fixed_effect_rows(summary)
    parsed = out["parameter"].map(parse_categorical_parameter).apply(pd.Series)
    out = pd.concat([out, parsed], axis=1)
    out = out.loc[out["variable"].notna()].copy()
    if variable_order:
        out = out.loc[out["variable"].isin(variable_order)].copy()
    out["reference"] = out.apply(
        lambda row: (
            row["reference"]
            if pd.notna(row["reference"]) and row["reference"] != ""
            else reference_labels.get(row["variable"], "reference")
        ),
        axis=1,
    )
    out["variable_label"] = out["variable"].map(variable_labels).fillna(out["variable"])
    out["contrast_label"] = out.apply(
        lambda row: f"{row['variable_label']}: {row['level']} vs {row['reference']}",
        axis=1,
    )
    out["model_kind"] = model_kind
    out["model_set"] = model_set
    out["summary_path"] = str(summary_path)
    return _select_plot_columns(
        out,
        family=family,
        extra_cols=[
            "variable",
            "variable_label",
            "level",
            "reference",
            "contrast_label",
        ],
    )


def parse_categorical_parameter(parameter: str) -> dict[str, str | None]:
    """Parse Bambi/formulae categorical parameter names.

    Supported examples include::

        C(age_band, Treatment(reference='20-24'))[30-34]
        C(dz_simd_quintile, Treatment(reference=1))[2]
        age_band[30-34]
        age_band[T.30-34]

    Unparseable parameters return ``None`` fields so callers can filter them.
    """
    parameter = str(parameter)
    match = _CATEGORICAL_RE.match(parameter)
    if match:
        inside = match.group("inside")
        level = _clean_level(match.group("level"))
        treatment = _TREATMENT_VARIABLE_RE.match(inside)
        if treatment:
            variable = treatment.group("variable").strip()
            reference = _clean_reference(treatment.group("reference"))
        else:
            variable = inside.strip()
            reference = None
        return {"variable": variable, "level": level, "reference": reference}

    match = _BARE_CATEGORICAL_RE.match(parameter)
    if match:
        return {
            "variable": match.group("variable").strip(),
            "level": _clean_level(match.group("level")),
            "reference": None,
        }

    return {"variable": None, "level": None, "reference": None}


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
    extra_cols: Sequence[str] = (),
) -> pd.DataFrame:
    keep = [
        *extra_cols,
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


def _estimate_columns(family: str) -> tuple[str, str, str]:
    if family == "logistic":
        return "OR_mean", "OR_hdi95_lb", "OR_hdi95_ub"
    return "mean", "hdi95_lb", "hdi95_ub"


def _mixing_feature_stem(parameter: str) -> str:
    for suffix in ("_obs_x10", "_z"):
        if str(parameter).endswith(suffix):
            return str(parameter)[: -len(suffix)]
    return str(parameter)


def _add_model_legend(
    fig: Figure,
    ax: Axes,
    colors: Mapping[str, str],
    model_labels: Mapping[str, str],
    point_size: float,
) -> None:
    handles = [
        ax.plot(
            [],
            [],
            marker="o",
            linestyle="",
            color=colors[key],
            label=model_labels.get(key, key.title()),
            markersize=point_size,
        )[0]
        for key in ("primary", "expanded")
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.08),
        frameon=False,
    )


def _finish_forest_figure(
    fig: Figure,
    axes: np.ndarray,
    *,
    title: str | None,
) -> None:
    if title:
        fig.suptitle(title)
    for ax in axes.ravel():
        ax.spines["left"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)


def _set_readable_or_ticks(ax: Axes) -> None:
    """Use sparse, readable odds-ratio ticks on a log axis."""
    xmin, xmax = ax.get_xlim()
    max_ticks = 5
    if xmax <= 1.35:
        preferred = np.array([0.67, 0.85, 1.0, 1.1, 1.2])
    elif xmax <= 3.2:
        preferred = np.array([0.67, 0.85, 1.0, 1.5, 2.5])
    else:
        preferred = np.array([0.67, 1.0, 1.5, 2.5, 5.0])
    ticks = preferred[(preferred >= xmin) & (preferred <= xmax)]
    if len(ticks) >= 2:
        ax.xaxis.set_major_locator(FixedLocator(list(ticks)))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        return

    candidates = np.array(
        [
            0.25,
            0.33,
            0.5,
            0.67,
            0.75,
            0.85,
            0.9,
            1.0,
            1.1,
            1.2,
            1.5,
            2.0,
            2.5,
            3.0,
            4.0,
            5.0,
        ]
    )
    ticks = candidates[(candidates >= xmin) & (candidates <= xmax)]
    if 1.0 >= xmin and 1.0 <= xmax:
        ticks = np.sort(np.unique(np.append(ticks, 1.0)))
    if len(ticks) > max_ticks:
        if 1.0 in ticks:
            others = ticks[ticks != 1.0]
            keep = np.linspace(0, len(others) - 1, max_ticks - 1).round().astype(int)
            ticks = np.sort(np.append(others[keep], 1.0))
        else:
            keep = np.linspace(0, len(ticks) - 1, max_ticks).round().astype(int)
            ticks = ticks[keep]
    if len(ticks) < 2:
        return
    ax.xaxis.set_major_locator(FixedLocator(list(ticks)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())


def _warn_missing(missing: Sequence[str], domain: str) -> None:
    if missing:
        warnings.warn(
            f"Skipped missing {domain} summaries: " + ", ".join(missing),
            stacklevel=2,
        )


def _y_lookup(row_order: Sequence[str]) -> dict[str, int]:
    return {row: i for i, row in enumerate(reversed(row_order))}


def _model_kind_categorical(values: pd.Series) -> pd.Categorical:
    return pd.Categorical(values, categories=["primary", "expanded"], ordered=True)


def _normalise_family(family: str) -> str:
    family = str(family).lower()
    if family not in {"logistic", "linear"}:
        raise ValueError("family must be 'logistic' or 'linear'.")
    return family


def _clean_reference(value: str) -> str:
    value = str(value).strip()
    if (value.startswith("'") and value.endswith("'")) or (
        value.startswith('"') and value.endswith('"')
    ):
        value = value[1:-1]
    return _clean_level(value)


def _clean_level(value: str) -> str:
    value = str(value).strip()
    if value.startswith("T."):
        value = value[2:]
    if (value.startswith("'") and value.endswith("'")) or (
        value.startswith('"') and value.endswith('"')
    ):
        value = value[1:-1]
    return value


def _composition_display_label(variable: str, value: Any) -> str:
    text = str(value)
    if variable == "dz_simd_quintile":
        return f"SIMD Q{text}"
    if variable in {"urban_rural_class", "dz_urban_rural_class"}:
        return _shorten_urban_rural_label(text)
    if variable in {"health_board", "dz_health_board"}:
        return _shorten_health_board_label(text)
    return text


def _shorten_urban_rural_label(value: str) -> str:
    replacements = {
        "Large Urban Areas": "Large urban",
        "Other Urban Areas": "Other urban",
        "Accessible Small Towns": "Accessible town",
        "Remote Small Towns": "Remote town",
        "Very Remote Small Towns": "Very remote town",
        "Accessible Rural": "Accessible rural",
        "Remote Rural": "Remote rural",
        "Very Remote Rural": "Very remote rural",
    }
    return replacements.get(value, value.replace(" Areas", "").replace(" Small", ""))


def _shorten_health_board_label(value: str) -> str:
    text = value.removeprefix("NHS ").removeprefix("Nhs ")
    return text.replace(" & ", " and ")


def _level_sort_key(value: str) -> tuple[int, float | str]:
    text = str(value)
    try:
        return (0, float(text))
    except ValueError:
        return (1, text)


def _slug(value: str) -> str:
    return str(value).replace("/", "_").replace(" ", "_")
