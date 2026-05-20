from __future__ import annotations

import re

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from utils import style

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
    size_col="cluster_size",
    *,
    by ="pango_lineage",
    min_size=1,
    width="double",
    width_in =None,
    height_in =3.5,
    context="paper",
    font_scale=1,
    complementary = True,
) -> plt.Figure:
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

    def _assign_wave_group(lineage: object) -> str:
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
    fig, axes = style.new_figure(
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
    ax_ecdf =axes[0]
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

    style.add_panel_labels([ax_ecdf, ax_violin])
    plt.close(fig)

    return fig


# ---------------------------------------------------------------------------
# Layer 1
# ---------------------------------------------------------------------------


def plot_role_dynamic_heatmap(
    candidates: pd.DataFrame,
    *,
    width="double",
    width_in=None,
    height_in=5,
    context="paper",
    font_scale= 1,
) -> plt.Figure:
    """Counts of ``sse_role`` x ``sse_onward_dynamic`` for SSE candidates.

    Cells are coloured on a log scale and annotated with raw counts.
    """
    role_order = [r for r in ROLE_ORDER if r in candidates["sse_role"].dropna().unique()]
    dyn_order = [d for d in DYNAMIC_ORDER if d in candidates["sse_onward_dynamic"].dropna().unique()]
    heat = (
        pd.crosstab(candidates["sse_role"], candidates["sse_onward_dynamic"])
        .reindex(index=role_order, columns=dyn_order, fill_value=0)
    ).T

    fig, ax = style.new_figure(
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
    width="double",
    width_in=None,
    height_in=5,
    context="paper",
    font_scale=1,
) -> plt.Figure:
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

    fig, axes = style.new_figure(
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
        color=style.lighten("#C75C2C", 0.2),
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

    style.add_panel_labels(list(axes))

    plt.close(fig)
    return fig


def plot_core_metric_space(
    node_stats: pd.DataFrame,
    *,
    height_in=3.5,
    context="paper",
    font_scale=1,
    min_size=1
) -> plt.Figure:
    """Scatter of core amplification vs onward dissemination, faceted by SSE candidate."""
    style.set_theme(
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
    )

    g.set_axis_labels("Core amplification score", "Onward dissemination score")
    g.set_titles("SSE candidate: {col_name}")

    if g._legend is not None:
        g._legend.remove()

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

    style.add_panel_labels(g.axes.flat)

    plt.close(g.figure)
    return g.figure


def plot_composite_distributions(
    node_stats: pd.DataFrame,
    *,
    columns=(
            "cluster_size",
            "core_amplification_score",
            "onward_dissemination_score",
    ),
    nrows=1,
    ncols=3,
    width="double",
    width_in=None,
    height_in=2.5,
    context="paper",
    font_scale=1,
    min_size=1
) -> plt.Figure:
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

    fig, axes = style.new_figure(
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
            ("background", "#8C8C8C", (~node_stats["sse_candidate"] &
                                       node_stats["cluster_size"].gt(min_size))),
            ("candidate", "#C75C2C", node_stats["sse_candidate"]),
        ]:
            values = node_stats.loc[mask, col].dropna().to_numpy()
            if len(values) < 5:
                continue
            sns.kdeplot(
                values, ax=ax, fill=True, color=color, alpha=0.35,
                linewidth=1.2, label=label, common_norm=False,
            )
        ax.set_xlabel(col.replace("_", " "))
        if col == "cluster_size":
            ax.set_xlabel("log(cluster size)")
        ax.set_ylabel("density")
    axes[0].legend(loc="best", frameon=False)
    plt.close(fig)
    return fig


def plot_socio_demo_breakdown(
    node_stats: pd.DataFrame,
    col="top_simd_quintiles",
    score="simd_entropy_obs",
    xlabels=("Deprivation mixing", "SIMD quintile (1 = most deprived)"),
    *,
    width="double",
    width_in=None,
    height_in=4,
    context="paper",
    font_scale=1,
    min_size=1
) -> tuple[plt.Figure, pd.DataFrame]:
    """
    Plot the distribution of mixing scores for candidate vs background nodes,
    and the distribution of class-label frequencies for candidate vs background nodes.
    """

    fig, axes = style.new_figure(
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
        values = node_stats.loc[mask, score].dropna().to_numpy()

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

    style.add_panel_labels(axes)

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
    order: str = "age",
    ax=None,
    annotate: bool = True,
    as_percent: bool = True,
    sig_alpha: float = 0.05,
) -> plt.Axes | plt.Figure:
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
