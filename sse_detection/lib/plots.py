from __future__ import annotations

import re
from typing import Any, Mapping, Sequence

from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from utils.style import (
    set_theme,
    WIDTHS,
    CONTEXTS,
    new_figure,
    add_panel_labels,
    lighten
)

from .palettes import (
    DYNAMIC_ORDER,
    ROLE_ORDER,
    ROLE_PALETTE,
    WAVE_GROUPS,
    WAVE_GROUP_PALETTE,
)


# ---------------------------------------------------------------------------
# Cluster sizes
# ---------------------------------------------------------------------------

def plot_cluster_size_distribution(
    df: pd.DataFrame,
    size_col: str = "cluster_size",
    *,
    by: str | None = "pango_lineage",
    min_size: int = 1,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 3.5,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    complementary: bool = True,
) -> Figure:
    """Two-panel plot of cluster-size distributions.

    Left panel:
        Complementary ECDF / CCDF of size_col on log-log axes.

    Right panel:
        Violin plot of log10(size_col), stratified by group.

    If ``by="pango_lineage"``, lineages matching WAVE_GROUPS are shown as
    individual wave groups and all remaining lineages are pooled as "Other".
    """

    df = df.loc[
        (df[size_col] >= min_size)
        & (df[size_col] > 0)
    ].copy()

    def _assign_wave_group(lineage: Any) -> str:
        if pd.isna(lineage):
            return "Other"

        lineage = str(lineage)

        for group_name, matcher in WAVE_GROUPS.items():
            if matcher(lineage):
                return group_name

        return "Other"

    high_contrast_palette = [
        "#000000",
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#999999",
        "#332288",
        "#88CCEE",
        "#44AA99",
    ]

    # ------------------------------------------------------------------
    # Build plotting group
    # ------------------------------------------------------------------
    if by is None or by not in df.columns:
        df["_plot_group"] = "All clusters"
        group_order = ["All clusters"]
        palette = {"All clusters": "#000000"}

    elif by == "pango_lineage":
        df["_plot_group"] = df[by].apply(_assign_wave_group)

        observed_groups = set(df["_plot_group"])

        group_order = [
            group_name
            for group_name in WAVE_GROUPS
            if group_name in observed_groups
        ]

        if "Other" in observed_groups:
            group_order.append("Other")

        palette = {
            group: WAVE_GROUP_PALETTE.get(group, "#8C8C8C")
            for group in group_order
        }

    else:
        df["_plot_group"] = df[by].astype(str)

        group_order = sorted(df["_plot_group"].dropna().unique())

        palette = {
            group: high_contrast_palette[i % len(high_contrast_palette)]
            for i, group in enumerate(group_order)
        }

    df["_log10_cluster_size"] = np.log10(df[size_col])

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = new_figure(
        nrows = 1,
        ncols = 2,
        gridspec_kw={"width_ratios": [1.15, 1.0]},
        layout="constrained",
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale
    )

    ax_ecdf = axes[0]
    ax_violin = axes[1]

    # ------------------------------------------------------------------
    # Left: complementary ECDF / CCDF
    # ------------------------------------------------------------------
    sns.ecdfplot(
        data=df,
        x=size_col,
        hue="_plot_group",
        hue_order=group_order,
        palette=palette,
        complementary=complementary,
        stat="proportion",
        linewidth=1.5,
        ax=ax_ecdf,
    )

    ax_ecdf.set_xscale("log")
    ax_ecdf.set_yscale("log")
    ax_ecdf.set_xlabel("Cluster size")
    ax_ecdf.set_ylabel("P(X ≥ cluster size)")

    leg = ax_ecdf.get_legend()
    if leg is not None:
        leg.set_title("")
        leg.set_frame_on(False)

    # ------------------------------------------------------------------
    # Right: violin plot on log10-transformed cluster size
    # ------------------------------------------------------------------
    sns.violinplot(
        data=df,
        x="_plot_group",
        y="_log10_cluster_size",
        hue="_plot_group",
        order=group_order,
        palette=palette,
        cut=0,
        inner="quartile",
        linewidth=0.8,
        ax=ax_violin,
    )

    ax_violin.set_xlabel("")
    ax_violin.set_ylabel("Cluster size")
    ax_violin.tick_params(axis="x", rotation=45)


    # Convert log10 tick labels back to original cluster sizes
    smin = np.floor(df["_log10_cluster_size"].min())
    smax = np.ceil(df["_log10_cluster_size"].max())

    log_ticks = np.arange(smin, smax + 1)
    size_ticks = 10 ** log_ticks

    # Violin plot: axis is already log10(cluster_size)
    ax_violin.set_yticks(log_ticks)
    ax_violin.set_yticklabels([f"{int(t):g}" for t in size_ticks])

    # ECDF plot: axis is raw cluster_size, displayed on log scale
    ax_ecdf.set_xticks(size_ticks)
    ax_ecdf.set_xticklabels([f"{int(t):g}" for t in size_ticks])

    add_panel_labels([ax_ecdf, ax_violin])
    plt.close(fig)

    return fig


# ---------------------------------------------------------------------------
# Layer 1
# ---------------------------------------------------------------------------


def plot_role_dynamic_heatmap(
    candidates: pd.DataFrame,
    *,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 5.0,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Counts of ``sse_role`` x ``sse_onward_dynamic`` for SSE candidates.

    Cells are coloured on a log scale and annotated with raw counts.
    """
    role_order = [r for r in ROLE_ORDER if r in candidates["sse_role"].dropna().unique()]
    dyn_order = [d for d in DYNAMIC_ORDER if d in candidates["sse_onward_dynamic"].dropna().unique()]
    heat = (
        pd.crosstab(candidates["sse_role"], candidates["sse_onward_dynamic"])
        .reindex(index=role_order, columns=dyn_order, fill_value=0)
    ).T

    fig, ax = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale
    )
    sns.heatmap(
        np.log10(heat + 1),
        annot=heat,
        fmt="d",
        cmap="YlGnBu",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "log10(n + 1)"},
        ax=ax,
    )
    ax.set_ylabel("Onward dynamic")
    ax.set_xlabel("Node role")
    ax.tick_params(axis="x", rotation=35)
    ax.tick_params(axis="y", rotation=0)

    plt.close(fig)
    return fig


def plot_candidate_rate_over_time(
    node_stats: pd.DataFrame,
    *,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 5.0,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Per-window candidate rate, with role composition stacked underneath.

    Top panel: % of nodes flagged ``sse_candidate`` per window, with the
    raw candidate count as a light bar in the background.
    Bottom panel: stacked counts of candidates by ``sse_role`` per window.
    """
    if "window_idx" not in node_stats.columns:
        raise KeyError("node_stats needs 'window_idx'")

    summary = (
        node_stats.groupby(["window_idx", "wn_mid_date"], as_index=False)
        .agg(
            n_nodes=("cluster_id", "nunique"),
            n_candidates=("sse_candidate", "sum"),
        )
    )
    summary["candidate_share"] = summary["n_candidates"] / summary["n_nodes"]

    role_counts = (
        node_stats.loc[node_stats["sse_candidate"]]
        .groupby(["wn_mid_date", "sse_role"], as_index=False)
        .size()
        .rename(columns={"size": "n"})
    )
    role_pivot = (
        role_counts.pivot(index="wn_mid_date", columns="sse_role", values="n")
        .fillna(0)
        .sort_index()
    )
    role_pivot = role_pivot.reindex(
        columns=[r for r in ROLE_ORDER if r in role_pivot.columns],
        fill_value=0,
    )

    fig, axes = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        nrows=2,
        sharex=True,
        context=context,
        font_scale=font_scale,
        layout="constrained",
    )
    ax = axes[0]
    ax.bar(
        summary["wn_mid_date"],
        summary["n_candidates"],
        width=5,
        color=lighten("#C75C2C", 0.2),
        alpha=0.85,
        label="candidate nodes",
    )
    ax2 = ax.twinx()
    ax2.plot(
        summary["wn_mid_date"],
        100 * summary["candidate_share"],
        color="#14151F",
        linewidth=1.8,
        label="candidate share",
    )
    ax.set_ylabel("candidates (n)")
    ax2.set_ylabel("candidate share (%)")
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, loc="best", ncol=1, frameon=False)

    ax = axes[1]
    if len(role_pivot.columns):
        colors = [ROLE_PALETTE.get(r, "#8C8C8C") for r in role_pivot.columns]
        ax.stackplot(
            role_pivot.index,
            [role_pivot[c].to_numpy() for c in role_pivot.columns],
            labels=list(role_pivot.columns),
            colors=colors,
            alpha=0.86,
        )
        ax.legend(loc="upper left", ncol=3, frameon=False)
    ax.set_ylabel("candidates by role")
    ax.set_xlabel("window midpoint")
    fig.autofmt_xdate()

    add_panel_labels(list(axes))

    plt.close(fig)
    return fig


def plot_core_metric_space(
    node_stats: pd.DataFrame,
    *,
    height_in: float = 5.0,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    min_size: int = 1,
) -> Figure:
    """Scatter of core amplification vs onward dissemination, faceted by SSE candidate."""
    set_theme(
        context=context,
        font_scale=font_scale,
    )

    plot_df = node_stats.loc[node_stats["cluster_size"].ge(min_size)].copy()

    clip_hi = plot_df["cluster_size"].quantile(0.995)
    plot_df["marker_size"] = np.sqrt(
        plot_df["cluster_size"].clip(lower=1, upper=clip_hi)
    )

    g = sns.relplot(
        data=plot_df,
        x="core_amplification_score",
        y="onward_dissemination_score",
        size="marker_size",
        hue="sse_role",
        palette=ROLE_PALETTE,
        col="sse_candidate",
        kind="scatter",
        alpha=0.3,
        sizes=(12, 220),
        height=height_in,
        facet_kws={"sharex": True, "sharey": True},
        legend=False,
    )

    g.set_axis_labels("Core amplification score", "Onward dissemination score")
    g.set_titles("SSE candidate: {col_name}")

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            label=role.replace("_", " ").capitalize(),
            markerfacecolor=ROLE_PALETTE[role],
            markeredgecolor=ROLE_PALETTE[role],
            markersize=7,
            alpha=1,
        )
        for role in plot_df["sse_role"].dropna().unique()
        if role in ROLE_PALETTE
    ]

    g.figure.legend(
        handles=handles,
        title="SSE role",
        loc="center right",
        bbox_to_anchor=(0.65, 0.5),
        frameon=False,
    )

    add_panel_labels(list(g.axes.flat))

    plt.close(g.figure)
    return g.figure


def plot_composite_distributions(
    node_stats: pd.DataFrame,
    *,
    columns: Sequence[str] = (
            "cluster_size",
            "core_amplification_score",
            "onward_dissemination_score",
            "mixing_score",
    ),
    nrows: int = 2,
    ncols: int = 2,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 2.5,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    min_size: int = 1,
) -> Figure:
    """Overlaid KDE of each composite score for candidates vs background.

    Visualises whether each component genuinely separates the two groups
    rather than just shifting the mean.
    """
    if "sse_candidate" not in node_stats.columns:
        raise KeyError("node_stats needs 'sse_candidate'")

    columns = [s for s in columns if s in node_stats.columns]
    if not columns:
        raise ValueError("None of the requested score columns are present.")

    if "cluster_size" in columns:
        node_stats = node_stats.copy()
        node_stats["cluster_size"] = np.log(node_stats["cluster_size"])

    fig, axes = new_figure(
        nrows=nrows,
        ncols=ncols,
        layout="constrained",
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
    )

    axes = axes.flatten()

    for ax, col in zip(axes, columns):
        for label, color, mask in [
            ("background", "#8C8C8C", (
                ~node_stats["sse_candidate"] & node_stats["cluster_size"].gt(min_size))
                ),
            ("candidate", "#C75C2C", node_stats["sse_candidate"]),
        ]:
            values = node_stats[mask][col].dropna().to_numpy()
            if len(values) < 5:
                continue
            sns.kdeplot(
                values, ax=ax, fill=True, color=color, alpha=0.35,
                linewidth=1.2, label=label, common_norm=False,
            )
        ax.set_xlabel(col.replace("_", " "))
        if col == "cluster_size":
            ax.set_xlabel("log(cluster size)")
        if col == "mixing_score":
            ax.set_xlabel("Observed socio-geodemographic entropy")
        ax.set_ylabel("density")
    axes[0].legend(loc="best", frameon=False)
    plt.close(fig)
    return fig


def plot_socio_demo_breakdown(
    node_stats: pd.DataFrame,
    col: str = "top_simd_quintiles",
    score: str = "simd_entropy_obs",
    xlabels: tuple[str, str] = ("Deprivation mixing", "SIMD quintile (1 = most deprived)"),
    *,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float = 4.0,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
    min_size: int = 1,
) -> tuple[Figure, pd.DataFrame]:
    """
    Plot the distribution of mixing scores for candidate vs background nodes,
    and the distribution of class-label frequencies for candidate vs background nodes.
    """

    fig, axes = new_figure(
        ncols=2,
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
        layout="constrained",
    )

    ax = axes[0]

    for label, color, mask in [
        (
            "background",
            "#8C8C8C",
            (
                ~node_stats["sse_candidate"]
                & node_stats["cluster_size"].gt(min_size)
            ),
        ),
        (
            "candidate",
            "#C75C2C",
            node_stats["sse_candidate"],
        ),
    ]:
        values = node_stats[mask][score].dropna().to_numpy()

        if len(values) < 5:
            continue

        sns.kdeplot(
            values,
            ax=ax,
            fill=True,
            color=color,
            alpha=0.35,
            linewidth=1.2,
            label=label,
            common_norm=False,
        )

    ax.set_xlabel(xlabels[0])
    ax.set_ylabel("density")
    ax.legend(loc="best", frameon=False)

    ax = axes[1]

    def _parse_freq_counts(s: str):
        """
        Parse strings like:
        'class label one (3); class label two (10); class label three (1)'

        Returns:
        [('class label one', 3), ('class label two', 10), ...]
        """
        if not isinstance(s, str) or not s.strip():
            return []

        out = []

        for part in s.split(";"):
            part = part.strip()

            if not part:
                continue

            match = re.match(r"^(.*?)\s*\((\d+)\)\s*$", part)

            if match is None:
                continue

            label = match.group(1).strip()
            count = int(match.group(2))

            if label and count > 0:
                out.append((label, count))

        return out

    records = []

    for _, row in node_stats.iterrows():
        candidate = row["sse_candidate"]

        for q, n in _parse_freq_counts(row[col]):
            records.append({
                "q": q,
                "candidate": candidate,
                "n": n,
            })

    if records:
        share = (
            pd.DataFrame(records)
            .groupby(["q", "candidate"], as_index=False)["n"]
            .sum()
        )

        denom = share.groupby("candidate")["n"].transform("sum")
        share["frac"] = share["n"] / denom

        order = sorted(share["q"].unique())

        bar_h = 0.4
        y = np.arange(len(order))

        for off, candidate_value, lbl, color in [
            (-bar_h / 2, False, "background", "#8C8C8C"),
            ( bar_h / 2, True,  "candidate",  "#C75C2C"),
        ]:
            sub = (
                share.loc[share["candidate"] == candidate_value]
                .set_index("q")
                .reindex(order, fill_value=0)
            )

            ax.barh(
                y + off,
                sub["frac"].to_numpy(),
                height=bar_h,
                color=color,
                label=lbl,
            )

        ax.set_yticks(y)
        ax.set_yticklabels(order)

    ax.set_ylabel(xlabels[1])
    ax.set_xlabel("Fraction of nodes")

    add_panel_labels(axes)

    plt.close(fig)
    return fig, pd.DataFrame(records)


def plot_socio_demo_candidate_background_diff(
    results_df: pd.DataFrame,
    *,
    q_col: str = "q",
    diff_col: str = "diff_candidate_minus_background",
    p_col: str | None = "p_adj_bh",
    xlabel: str = "Candidate − background fraction",
    ylabel: str = "Age band",
    title: str | None = None,
    order: str | list[str] = "age",
    ax: Axes | None = None,
    annotate: bool = True,
    as_percent: bool = True,
    sig_alpha: float = 0.05,
) -> Axes | Figure:
    """
    Plot signed percentage-point differences between candidate and background distributions.

    Positive values mean over-represented among candidates.
    Negative values mean under-represented among candidates.
    """

    df = results_df.copy()

    if q_col not in df or diff_col not in df:
        raise ValueError(f"`results_df` must contain `{q_col}` and `{diff_col}`.")

    df = df.dropna(subset=[q_col, diff_col])

    def _age_sort_key(label):
        label = str(label)

        if label.endswith("+"):
            return int(label.replace("+", ""))

        match = re.match(r"^(\d+)-(\d+)$", label)
        if match:
            return int(match.group(1))

        return label

    if order == "age":
        df = df.sort_values(q_col, key=lambda s: s.map(_age_sort_key))
    elif order == "effect":
        df = df.sort_values(diff_col)
    elif order == "abs_effect":
        df = df.sort_values(diff_col, key=lambda s: s.abs())
    elif isinstance(order, list):
        df[q_col] = pd.Categorical(df[q_col], categories=order, ordered=True)
        df = df.sort_values(q_col)
    else:
        raise ValueError("`order` must be 'age', 'effect', 'abs_effect', or a list.")

    plot_values = df[diff_col].to_numpy()

    if as_percent:
        plot_values = plot_values * 100
        if xlabel:
            xlabel = xlabel + " percentage points"

    y = np.arange(len(df))

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, max(3, 0.35 * len(df))))
    else:
        fig = ax.figure

    colors = np.where(plot_values >= 0, "#C75C2C", "#8C8C8C")

    ax.barh(
        y,
        plot_values,
        color=colors,
        height=0.75,
    )

    ax.axvline(0, color="black", linewidth=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels(df[q_col].astype(str))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title is not None:
        ax.set_title(title)

    if annotate:
        x_pad = max(abs(plot_values).max() * 0.03, 0.05)

        for i, value in enumerate(plot_values):
            ha = "left" if value >= 0 else "right"
            x = value + x_pad if value >= 0 else value - x_pad

            label = f"{value:+.1f}"

            if p_col is not None and p_col in df.columns:
                p = df.iloc[i][p_col]
                if pd.notna(p) and p < sig_alpha:
                    label += "*"

            ax.text(
                x,
                i,
                label,
                va="center",
                ha=ha,
                fontsize=8,
            )

    max_abs = max(abs(plot_values).max(), 0.01)
    ax.set_xlim(-max_abs * 1.2, max_abs * 1.2)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.close()

    return ax if ax is not None else fig


# ---------------------------------------------------------------------------
# Regression output figures and manuscript tables
# ---------------------------------------------------------------------------


_COMPOSITION_ORDER = [
    "sex",
    "age_band",
    "simd_quintile",
    "urban_rural_class",
    "health_board",
]

_MIXING_ORDER = [
    "sex_entropy_z",
    "age_entropy_z",
    "simd_entropy_z",
    "urban_rural_entropy_z",
    "health_board_entropy_z",
    "sex_entropy_obs",
    "age_entropy_obs",
    "simd_entropy_obs",
    "urban_rural_entropy_obs",
    "health_board_entropy_obs",
]

_PRETTY_LABELS = {
    "sex": "Sex",
    "age_band": "Age band",
    "simd_quintile": "SIMD quintile",
    "urban_rural_class": "Urban/rural class",
    "health_board": "Health board",
    "dz_simd_quintile": "SIMD quintile",
    "dz_urban_rural_class": "Urban/rural class",
    "dz_health_board": "Health board",
    "sex_entropy_obs": "Sex entropy",
    "age_entropy_obs": "Age entropy",
    "simd_entropy_obs": "SIMD entropy",
    "urban_rural_entropy_obs": "Urban/rural entropy",
    "health_board_entropy_obs": "Health-board entropy",
    "sex_entropy_z": "Sex entropy z-score",
    "age_entropy_z": "Age entropy z-score",
    "simd_entropy_z": "SIMD entropy z-score",
    "urban_rural_entropy_z": "Urban/rural entropy z-score",
    "health_board_entropy_z": "Health-board entropy z-score",
    "primary": "Primary",
    "expanded": "Expanded",
    "single": "Single predictor",
    "joint": "Joint",
    "composition": "Composition",
    "node_mixing": "Node mixing",
}


def _pretty_text(value: Any, label_map: Mapping[str, str] | None = None) -> str:
    """Human-readable label with conservative project-specific replacements."""
    if pd.isna(value):
        return ""
    text = str(value)
    if label_map and text in label_map:
        return label_map[text]
    if text in _PRETTY_LABELS:
        return _PRETTY_LABELS[text]
    return text.replace("_", " ").strip().capitalize()


def _term_level(term: Any) -> str:
    """Extract the displayed contrast level from a simple patsy term."""
    if pd.isna(term):
        return ""
    term = str(term)
    match = re.search(r"\[T\.(.*)\]$", term)
    if match:
        return match.group(1)
    return _pretty_text(term)


def _filter_regression_table(
    table: pd.DataFrame,
    *,
    domain: str | Sequence[str] | None = None,
    model_set: str | Sequence[str] | None = None,
    predictor_set: str | Sequence[str] | None = None,
    predictors: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Filter a regression output table without mutating the caller's data."""
    df = table.copy()
    string_cols = df.select_dtypes(include=["object", "string"]).columns
    for col in string_cols:
        df[col] = df[col].map(lambda x: x.strip() if isinstance(x, str) else x)

    for col, value in [
        ("domain", domain),
        ("model_set", model_set),
        ("predictor_set", predictor_set),
    ]:
        if value is None or col not in df.columns:
            continue
        values = [value] if isinstance(value, str) else list(value)
        df = df.loc[df[col].isin(values)]

    if predictors is not None and "predictor" in df.columns:
        df = df.loc[df["predictor"].isin(predictors)]

    return df.copy()


def _with_regression_display_labels(
    df: pd.DataFrame,
    *,
    label_map: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Add stable display labels used by both figures and tables."""
    out = df.copy()

    if "label" in out.columns:
        out["_predictor_label"] = out["label"].map(
            lambda x: _pretty_text(x, label_map)
        )
    elif "predictor" in out.columns:
        out["_predictor_label"] = out["predictor"].map(
            lambda x: _pretty_text(x, label_map)
        )
    else:
        out["_predictor_label"] = ""

    if "predictor" in out.columns:
        predictor_fallback = out["predictor"].map(lambda x: _pretty_text(x, label_map))
        all_mixing = out["_predictor_label"].str.lower().eq("all mixing predictors")
        out.loc[all_mixing, "_predictor_label"] = predictor_fallback.loc[all_mixing]
        if "term" in out.columns:
            term_fallback = out["term"].map(lambda x: _pretty_text(x, label_map))
            out.loc[all_mixing, "_predictor_label"] = term_fallback.loc[all_mixing]

    if "term" in out.columns:
        term_label = out["term"].map(_term_level)
    else:
        term_label = pd.Series("", index=out.index)
    out["_term_label"] = term_label

    if "reference" in out.columns:
        out["_contrast_label"] = np.where(
            out["_term_label"].eq(out["_predictor_label"]) | out["_term_label"].eq(""),
            out["_predictor_label"],
            out["_term_label"] + " vs " + out["reference"].astype(str),
        )
    else:
        out["_contrast_label"] = np.where(
            out["_term_label"].eq(out["_predictor_label"]) | out["_term_label"].eq(""),
            out["_predictor_label"],
            out["_term_label"],
        )

    if {"model_set", "predictor_set"}.issubset(out.columns):
        out["_model_label"] = (
            out["model_set"].map(lambda x: _pretty_text(x, label_map))
            + "\n"
            + out["predictor_set"].map(lambda x: _pretty_text(x, label_map))
        )
    elif "model_set" in out.columns:
        out["_model_label"] = out["model_set"].map(lambda x: _pretty_text(x, label_map))
    elif "predictor_set" in out.columns:
        out["_model_label"] = out["predictor_set"].map(
            lambda x: _pretty_text(x, label_map)
        )
    else:
        out["_model_label"] = ""

    return out


def _predictor_sort_key(values: pd.Series) -> pd.Series:
    order = {name: i for i, name in enumerate(_COMPOSITION_ORDER + _MIXING_ORDER)}
    return values.map(lambda x: order.get(str(x), len(order)))


def _regression_sort_source(df: pd.DataFrame) -> pd.Series:
    """Predictor-like series for stable display ordering."""
    if "predictor" not in df.columns:
        return pd.Series("", index=df.index)
    source = df["predictor"].astype(str).copy()
    if "term" in df.columns:
        use_term = source.eq("all_mixing") | source.eq("all_mixing_predictors")
        source.loc[use_term] = df.loc[use_term, "term"].astype(str)
    return source


def _format_p_value(
    value: Any,
    *,
    threshold: float = 0.001,
    digits: int = 3,
) -> str:
    if pd.isna(value):
        return ""
    value = float(value)
    if value < threshold:
        return f"<{threshold:.{digits}f}"
    return f"{value:.{digits}f}"


def _format_number(value: Any, *, digits: int = 2) -> str:
    if pd.isna(value):
        return ""
    return f"{float(value):.{digits}f}"


def _format_int(value: Any) -> str:
    if pd.isna(value):
        return ""
    return f"{int(round(float(value))):,}"


def _format_or_ci(row: pd.Series, *, digits: int = 2) -> str:
    if row[["odds_ratio", "or_low", "or_high"]].isna().any():
        return ""
    return (
        f"{row['odds_ratio']:.{digits}f} "
        f"({row['or_low']:.{digits}f}-{row['or_high']:.{digits}f})"
    )


def plot_regression_wald_heatmap(
    wald_df: pd.DataFrame,
    *,
    domain: str | Sequence[str] | None = None,
    model_set: str | Sequence[str] | None = None,
    predictor_set: str | Sequence[str] | None = None,
    predictors: Sequence[str] | None = None,
    p_col: str = "p_adj_bh",
    row_col: str | None = None,
    label_map: Mapping[str, str] | None = None,
    cap_neg_log10_p: float = 20.0,
    annotate_p: bool = True,
    title: str | None = None,
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float | None = None,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Heatmap of omnibus Wald evidence across regression specifications.

    The colour scale is ``-log10(p)`` using the BH-adjusted column when
    available. Exact zero p-values are capped for display and annotated as
    ``<1e-k`` rather than being dropped.
    """
    df = _filter_regression_table(
        wald_df,
        domain=domain,
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
    )
    if df.empty:
        raise ValueError("No Wald rows remain after filtering.")

    if p_col not in df.columns:
        fallback = "P>chi2"
        if fallback not in df.columns:
            raise ValueError(f"`wald_df` must contain `{p_col}` or `{fallback}`.")
        p_col = fallback

    df = _with_regression_display_labels(df, label_map=label_map)
    if row_col is None:
        row_col = "term" if "term" in df.columns else "predictor"
    if row_col not in df.columns:
        raise ValueError(f"`wald_df` does not contain row column `{row_col}`.")

    df["_row_id"] = df[row_col].astype(str)
    row_labels = (
        df.assign(_sort_source=_regression_sort_source(df))
        .sort_values("_sort_source", key=_predictor_sort_key)
        .drop_duplicates("_row_id")
        .set_index("_row_id")["_predictor_label"]
    )
    missing_label = row_labels.str.lower().eq("all mixing predictors")
    if missing_label.any() and "term" in df.columns:
        replacements = (
            df.drop_duplicates("_row_id").set_index("_row_id")["term"].map(_pretty_text)
        )
        row_labels.loc[missing_label] = replacements.loc[missing_label]

    df["_p_for_plot"] = pd.to_numeric(df[p_col], errors="coerce")
    positive = df.loc[df["_p_for_plot"].gt(0), "_p_for_plot"]
    min_positive = positive.min() if not positive.empty else 10 ** (-cap_neg_log10_p)
    zero_floor = min(min_positive, 10 ** (-cap_neg_log10_p))
    df["_neg_log10_p"] = -np.log10(df["_p_for_plot"].replace(0, zero_floor))
    df["_neg_log10_p"] = df["_neg_log10_p"].clip(upper=cap_neg_log10_p)

    pivot = df.pivot_table(
        index="_row_id",
        columns="_model_label",
        values="_neg_log10_p",
        aggfunc="first",
    )
    row_order = [row for row in row_labels.index if row in pivot.index]
    pivot = pivot.reindex(row_order)

    if {"model_set", "predictor_set"}.issubset(df.columns):
        col_order_df = (
            df[["model_set", "predictor_set", "_model_label"]]
            .drop_duplicates()
            .assign(
                _model_set_order=lambda x: x["model_set"].map(
                    {"primary": 0, "expanded": 1}
                ).fillna(99),
                _predictor_set_order=lambda x: x["predictor_set"].map(
                    {"single": 0, "joint": 1}
                ).fillna(99),
            )
            .sort_values(["_model_set_order", "_predictor_set_order", "model_set", "predictor_set"])
        )
        pivot = pivot.reindex(columns=col_order_df["_model_label"])

    annot = None
    if annotate_p:
        p_pivot = df.pivot_table(
            index="_row_id",
            columns="_model_label",
            values=p_col,
            aggfunc="first",
        ).reindex(index=pivot.index, columns=pivot.columns)
        annot = p_pivot.apply(
            lambda col: col.map(
                lambda x: "" if pd.isna(x) else (
                    f"<1e-{int(cap_neg_log10_p)}"
                    if float(x) == 0 else _format_p_value(x)
                )
            )
        )

    if height_in is None:
        height_in = max(2.5, 0.35 * len(pivot.index) + 1.2)

    fig, ax = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
    )
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="YlOrRd",
        vmin=0,
        vmax=cap_neg_log10_p,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "-log10(FDR-adjusted p)"},
        annot=annot,
        fmt="" if annotate_p else ".2f",
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels([row_labels.get(row, row) for row in pivot.index], rotation=0)
    ax.tick_params(axis="x", rotation=0)
    if title is not None:
        ax.set_title(title)

    plt.close(fig)
    return fig


def plot_regression_odds_ratio_forest(
    odds_df: pd.DataFrame,
    *,
    domain: str | Sequence[str] | None = None,
    model_set: str | Sequence[str] | None = "primary",
    predictor_set: str | Sequence[str] | None = "single",
    predictors: Sequence[str] | None = None,
    terms: Sequence[str] | None = None,
    p_col: str = "p_value",
    label_map: Mapping[str, str] | None = None,
    sort_by: str = "table",
    max_rows: int | None = 40,
    sig_alpha: float = 0.05,
    title: str | None = None,
    xlabel: str = "Odds ratio for candidate-node membership",
    width: WIDTHS = "double",
    width_in: float | None = None,
    height_in: float | None = None,
    context: CONTEXTS = "paper",
    font_scale: float = 1.0,
) -> Figure:
    """Forest plot of odds ratios and 95% confidence intervals."""
    required = {"odds_ratio", "or_low", "or_high", "term"}
    missing = required - set(odds_df.columns)
    if missing:
        raise ValueError(f"`odds_df` is missing required columns: {sorted(missing)}")

    df = _filter_regression_table(
        odds_df,
        domain=domain,
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
    )
    if terms is not None:
        df = df.loc[df["term"].isin(terms)].copy()
    if df.empty:
        raise ValueError("No odds-ratio rows remain after filtering.")

    df = _with_regression_display_labels(df, label_map=label_map)
    for col in ["odds_ratio", "or_low", "or_high"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["odds_ratio", "or_low", "or_high"]
    )
    df = df.loc[df[["odds_ratio", "or_low", "or_high"]].gt(0).all(axis=1)]
    if df.empty:
        raise ValueError("No finite positive odds ratios remain after filtering.")

    if sort_by == "table":
        sort_cols = []
        if "model_set" in df.columns:
            df["_model_set_order"] = df["model_set"].map({"primary": 0, "expanded": 1}).fillna(99)
            sort_cols.append("_model_set_order")
        if "predictor_set" in df.columns:
            df["_predictor_set_order"] = df["predictor_set"].map({"single": 0, "joint": 1}).fillna(99)
            sort_cols.append("_predictor_set_order")
        if "predictor" in df.columns:
            df["_predictor_order"] = _predictor_sort_key(_regression_sort_source(df))
            sort_cols.append("_predictor_order")
        sort_cols.append("term")
        df = df.sort_values(sort_cols)
    elif sort_by == "odds_ratio":
        df = df.sort_values("odds_ratio")
    elif sort_by == "abs_log_or":
        df = df.assign(_abs_log_or=np.log(df["odds_ratio"]).abs()).sort_values("_abs_log_or")
    elif sort_by == "p_value":
        if p_col not in df.columns:
            raise ValueError(f"`odds_df` does not contain `{p_col}`.")
        df = df.sort_values(p_col)
    else:
        raise ValueError("`sort_by` must be 'table', 'odds_ratio', 'abs_log_or', or 'p_value'.")

    if max_rows is not None and len(df) > max_rows:
        df = df.tail(max_rows) if sort_by == "abs_log_or" else df.head(max_rows)

    df = df.reset_index(drop=True)
    y = np.arange(len(df))
    colors = np.full(len(df), "#6C6F73", dtype=object)
    if p_col in df.columns:
        colors = np.where(pd.to_numeric(df[p_col], errors="coerce") < sig_alpha, "#C75C2C", "#6C6F73")

    if height_in is None:
        height_in = max(3.0, 0.28 * len(df) + 1.4)

    fig, ax = new_figure(
        width=width,
        width_in=width_in,
        height_in=height_in,
        context=context,
        font_scale=font_scale,
    )

    xerr = np.vstack([
        df["odds_ratio"].to_numpy() - df["or_low"].to_numpy(),
        df["or_high"].to_numpy() - df["odds_ratio"].to_numpy(),
    ])
    ax.errorbar(
        df["odds_ratio"],
        y,
        xerr=xerr,
        fmt="none",
        ecolor="#5B5F66",
        elinewidth=1.0,
        capsize=2,
        zorder=1,
    )
    ax.scatter(df["odds_ratio"], y, c=colors, s=28, zorder=2)

    ax.axvline(1, color="black", linewidth=0.8)
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.set_yticks(y)
    ax.set_yticklabels(df["_contrast_label"])
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("")
    if title is not None:
        ax.set_title(title)

    finite_values = df[["or_low", "or_high"]].to_numpy().ravel()
    finite_values = finite_values[np.isfinite(finite_values) & (finite_values > 0)]
    if finite_values.size:
        lo, hi = finite_values.min(), finite_values.max()
        pad = np.exp(0.08 * (np.log(hi) - np.log(lo) if hi > lo else 1))
        ax.set_xlim(lo / pad, hi * pad)

    ax.grid(axis="x", color="#D9DDE3", linewidth=0.5)
    plt.close(fig)
    return fig


def make_regression_wald_table(
    wald_df: pd.DataFrame,
    *,
    domain: str | Sequence[str] | None = None,
    model_set: str | Sequence[str] | None = None,
    predictor_set: str | Sequence[str] | None = None,
    predictors: Sequence[str] | None = None,
    p_col: str = "P>chi2",
    p_adj_col: str = "p_adj_bh",
    label_map: Mapping[str, str] | None = None,
    digits: int = 2,
) -> pd.DataFrame:
    """Return a manuscript-facing omnibus Wald table.

    The table keeps model identifiers and sample sizes explicit while adding
    formatted chi-square, raw p, and BH-adjusted p columns for display.
    """
    df = _filter_regression_table(
        wald_df,
        domain=domain,
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
    )
    if df.empty:
        return pd.DataFrame()

    df = _with_regression_display_labels(df, label_map=label_map)
    df["_predictor_order"] = _predictor_sort_key(_regression_sort_source(df))
    df = df.sort_values(
        [c for c in ["domain", "model_set", "predictor_set"] if c in df.columns]
        + (["_predictor_order"] if "predictor" in df.columns else []),
    )

    out = pd.DataFrame(index=df.index)
    if "domain" in df.columns:
        out["Domain"] = df["domain"].map(lambda x: _pretty_text(x, label_map))
    if "model_set" in df.columns:
        out["Model set"] = df["model_set"].map(lambda x: _pretty_text(x, label_map))
    if "predictor_set" in df.columns:
        out["Specification"] = df["predictor_set"].map(lambda x: _pretty_text(x, label_map))
    out["Predictor"] = df["_predictor_label"]
    if "reference" in df.columns:
        out["Reference/scale"] = df["reference"].astype(str)
    if "df" in df.columns:
        out["df"] = df["df"].map(_format_int)
    if "chi2" in df.columns:
        out["Wald chi-square"] = df["chi2"].map(lambda x: _format_number(x, digits=digits))
    if p_col in df.columns:
        out["P value"] = df[p_col].map(_format_p_value)
    if p_adj_col in df.columns:
        out["BH-adjusted P value"] = df[p_adj_col].map(_format_p_value)
    for source, target in [
        ("n_model_rows", "Model rows"),
        ("n_sequences", "Sequences"),
        ("n_nodes", "Nodes"),
    ]:
        if source in df.columns:
            out[target] = df[source].map(_format_int)

    return out.reset_index(drop=True)


def make_regression_odds_ratio_table(
    odds_df: pd.DataFrame,
    *,
    domain: str | Sequence[str] | None = None,
    model_set: str | Sequence[str] | None = None,
    predictor_set: str | Sequence[str] | None = None,
    predictors: Sequence[str] | None = None,
    terms: Sequence[str] | None = None,
    p_col: str = "p_value",
    label_map: Mapping[str, str] | None = None,
    digits: int = 2,
) -> pd.DataFrame:
    """Return a manuscript-facing coefficient odds-ratio table."""
    df = _filter_regression_table(
        odds_df,
        domain=domain,
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
    )
    if terms is not None and "term" in df.columns:
        df = df.loc[df["term"].isin(terms)].copy()
    if df.empty:
        return pd.DataFrame()

    df = _with_regression_display_labels(df, label_map=label_map)
    if "predictor" in df.columns:
        df["_predictor_order"] = _predictor_sort_key(_regression_sort_source(df))
        df = df.sort_values(
            [c for c in ["domain", "model_set", "predictor_set", "_predictor_order", "term"] if c in df.columns]
        )

    out = pd.DataFrame(index=df.index)
    if "domain" in df.columns:
        out["Domain"] = df["domain"].map(lambda x: _pretty_text(x, label_map))
    if "model_set" in df.columns:
        out["Model set"] = df["model_set"].map(lambda x: _pretty_text(x, label_map))
    if "predictor_set" in df.columns:
        out["Specification"] = df["predictor_set"].map(lambda x: _pretty_text(x, label_map))
    out["Predictor"] = df["_predictor_label"]
    out["Contrast"] = df["_contrast_label"]
    if "reference" in df.columns:
        out["Reference/scale"] = df["reference"].astype(str)
    out["Odds ratio (95% CI)"] = df.apply(lambda row: _format_or_ci(row, digits=digits), axis=1)
    if p_col in df.columns:
        out["P value"] = df[p_col].map(_format_p_value)
    for source, target in [
        ("n_model_rows", "Model rows"),
        ("n_sequences", "Sequences"),
        ("n_nodes", "Nodes"),
    ]:
        if source in df.columns:
            out[target] = df[source].map(_format_int)

    return out.reset_index(drop=True)


def make_regression_fit_table(
    fit_df: pd.DataFrame,
    *,
    domain: str | Sequence[str] | None = None,
    model_set: str | Sequence[str] | None = None,
    predictor_set: str | Sequence[str] | None = None,
    predictors: Sequence[str] | None = None,
    label_map: Mapping[str, str] | None = None,
    digits: int = 3,
) -> pd.DataFrame:
    """Return a manuscript-facing model fit comparison table."""
    df = _filter_regression_table(
        fit_df,
        domain=domain,
        model_set=model_set,
        predictor_set=predictor_set,
        predictors=predictors,
    )
    if df.empty:
        return pd.DataFrame()

    df = _with_regression_display_labels(df, label_map=label_map)
    if "predictor" in df.columns:
        df["_predictor_order"] = _predictor_sort_key(_regression_sort_source(df))
        df = df.sort_values(
            [c for c in ["domain", "model_set", "predictor_set", "_predictor_order"] if c in df.columns]
        )

    out = pd.DataFrame(index=df.index)
    if "domain" in df.columns:
        out["Domain"] = df["domain"].map(lambda x: _pretty_text(x, label_map))
    if "model_set" in df.columns:
        out["Model set"] = df["model_set"].map(lambda x: _pretty_text(x, label_map))
    if "predictor_set" in df.columns:
        out["Specification"] = df["predictor_set"].map(lambda x: _pretty_text(x, label_map))
    out["Predictor"] = df["_predictor_label"]
    if "r2_mcfadden" in df.columns:
        out["McFadden pseudo-R2"] = df["r2_mcfadden"].map(lambda x: _format_number(x, digits=digits))
    for source, target in [
        ("aic", "AIC"),
        ("bic_llf", "BIC"),
        ("log_likelihood", "Log likelihood"),
    ]:
        if source in df.columns:
            out[target] = df[source].map(lambda x: _format_number(x, digits=1))
    for source, target in [
        ("n_model_rows", "Model rows"),
        ("n_sequences", "Sequences"),
        ("n_nodes", "Nodes"),
    ]:
        if source in df.columns:
            out[target] = df[source].map(_format_int)
    if "converged" in df.columns:
        out["Converged"] = df["converged"].map(lambda x: "" if pd.isna(x) else str(bool(x)))

    return out.reset_index(drop=True)
