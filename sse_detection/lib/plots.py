"""Figure functions for the SSE-detection output notebook.

Each function returns a matplotlib ``Figure`` so the notebook can decide
how to display and save. All functions use ``utils.style.new_figure`` for
size/font/spine consistency with the rest of the Scotland clustering work.

Figure list (matches the four-question structure of the notebook):

Overview
    * :func:`plot_sequence_volume_timeline`
    * :func:`plot_cluster_size_ccdf`
Layer 1 (node-level SSE signatures)
    * :func:`plot_role_dynamic_heatmap`
    * :func:`plot_candidate_rate_over_time`
    * :func:`plot_metric_space_scatter`
    * :func:`plot_composite_score_distributions`
Layer 2 (meta-cluster weekly growth)
    * :func:`plot_meta_cluster_trajectories`
    * :func:`plot_norm_change_histogram`
    * :func:`plot_threshold_sensitivity`
Layer 1 x Layer 2
    * :func:`plot_layer_concordance`
Spatial / demographic
    * :func:`plot_simd_breakdown`
Methods / robustness
    * :func:`plot_null_comparison`
"""

from __future__ import annotations

from typing import Iterable, Sequence

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from utils import style

from .palettes import (
    DYNAMIC_ORDER,
    ROLE_ORDER,
    ROLE_PALETTE,
)


# ---------------------------------------------------------------------------
# Overview
# ---------------------------------------------------------------------------


def plot_sequence_volume_timeline(
    node_stats: pd.DataFrame,
    *,
    by: str = "who_voc_plot",
    title: str = "Cluster volume over time",
    width: str = "slide",
    height_in: float = 4.2,
) -> plt.Figure:
    """Stacked area of cluster volume per window, coloured by ``by``.

    "Volume" is summed ``cluster_size`` per (window, group), which gives a
    sense of how many sequences belong to each lineage in each window.
    """
    if by not in node_stats.columns:
        raise KeyError(f"{by!r} not in node_stats columns")
    if "wn_mid_date" not in node_stats.columns:
        raise KeyError("node_stats needs a 'wn_mid_date' column")

    g = (
        node_stats.groupby(["wn_mid_date", by], as_index=False)["cluster_size"]
        .sum()
        .rename(columns={"cluster_size": "n_sequences"})
    )
    pivot = (
        g.pivot(index="wn_mid_date", columns=by, values="n_sequences")
        .fillna(0)
        .sort_index()
    )
    voc_order = [c for c in style.WHO_VOC_PALETTE if c in pivot.columns]
    extra = [c for c in pivot.columns if c not in voc_order]
    ordered = voc_order + sorted(extra)
    colors = [style.WHO_VOC_PALETTE.get(c, "#8C8C8C") for c in ordered]

    fig, ax = style.new_figure(width=width, height_in=height_in, context="talk", font_scale=0.8)
    ax.stackplot(
        pivot.index,
        [pivot[c].to_numpy() for c in ordered],
        labels=ordered,
        colors=colors,
        alpha=0.88,
    )
    ax.set_ylabel("clustered sequences per window")
    ax.set_xlabel("window midpoint")
    ax.set_title(title)
    ax.legend(loc="upper left", ncol=min(5, len(ordered)), frameon=False)
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    fig.autofmt_xdate()
    fig.tight_layout()
    return fig


def plot_cluster_size_ccdf(
    node_stats: pd.DataFrame,
    *,
    by: str | None = "who_voc_plot",
    min_size: int = 1,
    title: str = "Cluster-size complementary CDF",
    width: str = "slide",
    height_in: float = 4.6,
) -> plt.Figure:
    """Log-log CCDF of ``cluster_size``, optionally split by ``by``.

    The tail is the relevant region for superspreading; this plot makes it
    visible without any thresholding.
    """
    df = node_stats.loc[node_stats["cluster_size"] >= min_size].copy()

    fig, ax = style.new_figure(width=width, height_in=height_in, context="talk", font_scale=0.85)

    def _ccdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if len(values) == 0:
            return np.array([]), np.array([])
        x = np.sort(values)
        n = len(x)
        # P(X >= x) using rank
        ccdf = 1.0 - (np.arange(1, n + 1) - 1) / n
        return x, ccdf

    if by is None or by not in df.columns:
        x, y = _ccdf(df["cluster_size"].to_numpy())
        ax.plot(x, y, color="#3A6EA5", lw=1.6)
    else:
        groups = sorted(df[by].dropna().unique())
        # Use WHO VOC ordering when applicable.
        if by.startswith("who_voc"):
            voc_order = [c for c in style.WHO_VOC_PALETTE if c in groups]
            extra = [c for c in groups if c not in voc_order]
            groups = voc_order + sorted(extra)
        for g in groups:
            sub = df.loc[df[by] == g, "cluster_size"].to_numpy()
            if len(sub) == 0:
                continue
            x, y = _ccdf(sub)
            color = (
                style.WHO_VOC_PALETTE.get(g, "#8C8C8C")
                if by.startswith("who_voc")
                else None
            )
            ax.plot(x, y, lw=1.4, alpha=0.9, label=g, color=color)
        ax.legend(loc="lower left", ncol=2, frameon=False)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("cluster size")
    ax.set_ylabel("P(X >= cluster size)")
    ax.set_title(title)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Layer 1
# ---------------------------------------------------------------------------


def plot_role_dynamic_heatmap(
    candidates: pd.DataFrame,
    *,
    title: str = "SSE-like candidates by role and onward dynamic",
    width: str = "slide",
    height_in: float = 5.6,
) -> plt.Figure:
    """Counts of ``sse_role`` x ``sse_onward_dynamic`` for SSE candidates.

    Cells are coloured on a log scale and annotated with raw counts.
    """
    role_order = [r for r in ROLE_ORDER if r in candidates["sse_role"].dropna().unique()]
    dyn_order = [d for d in DYNAMIC_ORDER if d in candidates["sse_onward_dynamic"].dropna().unique()]
    heat = (
        pd.crosstab(candidates["sse_role"], candidates["sse_onward_dynamic"])
        .reindex(index=role_order, columns=dyn_order, fill_value=0)
    )

    fig, ax = style.new_figure(width=width, height_in=height_in, context="talk", font_scale=0.72)
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
    ax.set_xlabel("onward dynamic")
    ax.set_ylabel("node role")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    return fig


def plot_candidate_rate_over_time(
    node_stats: pd.DataFrame,
    *,
    title: str = "SSE-candidate rate per sliding window",
    width: str = "slide",
    height_in: float = 5.0,
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
        width=width, height_in=height_in, nrows=2, sharex=True,
        context="talk", font_scale=0.78,
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
    ax.set_title(title)
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, loc="upper left", ncol=2, frameon=False)

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
    fig.tight_layout()
    return fig


def plot_metric_space_scatter(
    candidates: pd.DataFrame,
    *,
    x: str = "core_amplification_score",
    y: str = "onward_dissemination_score",
    color_by: str = "sse_role",
    size_col: str = "cluster_size",
    sample_cap: int = 8000,
    random_state: int = 42,
    title: str | None = None,
    width: str = "slide",
    height_in: float = 5.4,
) -> plt.Figure:
    """Scatter of ``x`` vs ``y`` for SSE candidates, coloured by ``color_by``.

    Point size encodes ``cluster_size``. The plot is downsampled to
    ``sample_cap`` rows by default so it remains legible.
    """
    df = candidates.dropna(subset=[x, y]).copy()
    if len(df) > sample_cap:
        df = df.sample(sample_cap, random_state=random_state)
    if size_col in df.columns:
        clip_hi = df[size_col].quantile(0.99)
        df["_size"] = np.sqrt(df[size_col].clip(lower=1, upper=clip_hi)) * 9.0
    else:
        df["_size"] = 18.0

    fig, ax = style.new_figure(width=width, height_in=height_in, context="talk", font_scale=0.82)
    if color_by in df.columns:
        for group, sub in df.groupby(color_by):
            ax.scatter(
                sub[x], sub[y],
                s=sub["_size"], alpha=0.28,
                color=ROLE_PALETTE.get(group, "#8C8C8C") if color_by == "sse_role" else None,
                label=group, edgecolor="none",
            )
        ax.legend(loc="upper left", ncol=2, markerscale=1.6, frameon=False)
    else:
        ax.scatter(df[x], df[y], s=df["_size"], alpha=0.35, color="#3A6EA5", edgecolor="none")
    ax.set_xlabel(x.replace("_", " "))
    ax.set_ylabel(y.replace("_", " "))
    ax.set_title(title or f"{x} vs {y} (SSE candidates)")
    if df[x].between(0, 1.01).all():
        ax.set_xlim(0, 1.02)
    if df[y].between(0, 1.01).all():
        ax.set_ylim(0, 1.02)
    fig.tight_layout()
    return fig


def plot_composite_score_distributions(
    node_stats: pd.DataFrame,
    *,
    scores: Sequence[str] = (
        "core_amplification_score",
        "onward_dissemination_score",
        "mixing_score",
    ),
    title: str = "Composite scores: candidates vs background",
    width: str = "slide",
    height_in: float = 3.6,
) -> plt.Figure:
    """Overlaid KDE of each composite score for candidates vs background.

    Visualises whether each component genuinely separates the two groups
    rather than just shifting the mean.
    """
    if "sse_candidate" not in node_stats.columns:
        raise KeyError("node_stats needs 'sse_candidate'")

    scores = [s for s in scores if s in node_stats.columns]
    if not scores:
        raise ValueError("None of the requested score columns are present.")

    fig, axes = style.new_figure(
        width=width, height_in=height_in, ncols=len(scores),
        context="talk", font_scale=0.78,
    )
    if len(scores) == 1:
        axes = [axes]

    for ax, score in zip(axes, scores):
        for label, color, mask in [
            ("background", "#8C8C8C", ~node_stats["sse_candidate"]),
            ("candidate", "#C75C2C", node_stats["sse_candidate"]),
        ]:
            values = node_stats.loc[mask, score].dropna().to_numpy()
            if len(values) < 5:
                continue
            sns.kdeplot(
                values, ax=ax, fill=True, color=color, alpha=0.35,
                linewidth=1.2, label=label, common_norm=False, clip=(0, 1),
            )
        ax.set_xlabel(score.replace("_", " "))
        ax.set_ylabel("density")
        ax.set_xlim(0, 1)
        ax.legend(loc="upper left", frameon=False)
    axes[0].figure.suptitle(title, y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Layer 2
# ---------------------------------------------------------------------------


def plot_meta_cluster_trajectories(
    weekly_growth: pd.DataFrame,
    meta_summary: pd.DataFrame | None = None,
    *,
    top_n: int = 12,
    rank_by: str = "n_sequences",
    title: str = "Meta-cluster cumulative size with SSE weeks marked",
    width: str = "slide",
    height_in: float = 5.4,
) -> plt.Figure:
    """Cumulative-size trajectories for the largest meta-clusters.

    Each panel shows one meta-cluster: x-axis is calendar week, y-axis is
    cumulative size (linear), and weeks flagged ``is_sse`` are highlighted.
    Weeks with no newly observed sequences are omitted so inactive flat
    stretches do not dominate the panel range.

    ``weekly_growth`` is expected to be the full weekly table (not the
    SSE-only filter) so the full curve is visible.
    """
    required = {"meta_cluster_id", "week", "new_sequences", "cc_size", "is_sse"}
    missing = sorted(required - set(weekly_growth.columns))
    if missing:
        cols = ", ".join(repr(c) for c in missing)
        raise KeyError(f"weekly_growth must include required column(s): {cols}")

    active_weekly = weekly_growth.loc[weekly_growth["new_sequences"] > 0].copy()
    if active_weekly.empty:
        raise ValueError("No weeks with new_sequences > 0 available to plot.")

    if meta_summary is None:
        # fall back to rank by max cc_size in the weekly table
        rank = (
            active_weekly.groupby("meta_cluster_id")["cc_size"]
            .max()
            .sort_values(ascending=False)
            .head(top_n)
            .index
        )
    else:
        rank = (
            meta_summary.sort_values(rank_by, ascending=False)
            .head(top_n)["meta_cluster_id"]
            .to_numpy()
        )
    active_meta_ids = set(active_weekly["meta_cluster_id"].dropna().unique())
    rank = [meta_id for meta_id in rank if meta_id in active_meta_ids]
    if not rank:
        raise ValueError("No meta-clusters with new_sequences > 0 available to plot.")

    sub = active_weekly.loc[active_weekly["meta_cluster_id"].isin(rank)].copy()
    sub["week"] = pd.to_datetime(sub["week"])

    ncols = min(3, len(rank))
    nrows = int(np.ceil(len(rank) / ncols))
    fig, axes = style.new_figure(
        width=width, height_in=height_in * nrows / 3,
        nrows=nrows, ncols=ncols, sharex=False,
        context="talk", font_scale=0.7,
    )
    axes = np.atleast_2d(axes).flatten()

    for ax, meta_id in zip(axes, rank):
        m = sub.loc[sub["meta_cluster_id"] == meta_id].sort_values("week")
        ax.plot(m["week"], m["cc_size"], color="#3A6EA5", linewidth=1.5)
        sse = m.loc[m["is_sse"]]
        if not sse.empty:
            ax.scatter(
                sse["week"], sse["cc_size"],
                color="#C75C2C", s=24, zorder=3, label="SSE week",
            )
        ax.set_title(meta_id, fontsize=9)
        locator = mdates.AutoDateLocator(minticks=3, maxticks=5)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        ax.tick_params(axis="x", labelrotation=35)
        ax.set_ylabel("cum. size")

    for ax in axes[len(rank):]:
        ax.set_visible(False)

    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    return fig


def plot_norm_change_histogram(
    weekly_growth: pd.DataFrame,
    *,
    threshold: float = 9.0,
    bins: int = 60,
    title: str = "Distribution of weekly normalised change",
    width: str = "slide",
    height_in: float = 3.8,
) -> plt.Figure:
    """Histogram (log y) of ``norm_change`` with the SSE threshold marked.

    ``weekly_growth`` should be the full table, not the SSE-only filter.
    """
    vals = weekly_growth["norm_change"].dropna().to_numpy()
    vals = vals[np.isfinite(vals)]

    fig, ax = style.new_figure(width=width, height_in=height_in, context="talk", font_scale=0.82)
    ax.hist(vals, bins=bins, color="#3A6EA5", edgecolor="white", alpha=0.85)
    ax.axvline(threshold, color="#C75C2C", linestyle="--", linewidth=1.4, label=f"SSE threshold = {threshold:g}")
    ax.set_yscale("log")
    ax.set_xlabel("normalised weekly change")
    ax.set_ylabel("weeks (log)")
    ax.set_title(title)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    return fig


def plot_threshold_sensitivity(
    weekly_growth: pd.DataFrame,
    *,
    thresholds: Iterable[float] = (3, 5, 7, 9, 11, 13, 15, 20),
    title: str = "SSE-week count vs threshold",
    width: str = "onehalf",
    height_in: float = 3.4,
) -> plt.Figure:
    """Sweep the ``flag_sse`` threshold and count the resulting SSE weeks.

    Lets you see how brittle the choice of 9 is.
    """
    vals = weekly_growth["norm_change"].dropna().to_numpy()
    thr = np.asarray(list(thresholds), dtype=float)
    counts = np.array([(vals > t).sum() for t in thr])
    affected_meta = []
    for t in thr:
        affected_meta.append(
            weekly_growth.loc[weekly_growth["norm_change"] > t, "meta_cluster_id"].nunique()
        )

    fig, ax = style.new_figure(width=width, height_in=height_in, context="talk", font_scale=0.82)
    ax.plot(thr, counts, "-o", color="#3A6EA5", label="SSE weeks")
    ax.plot(thr, affected_meta, "-s", color="#C75C2C", label="unique meta-clusters")
    ax.axvline(9.0, color="#14151F", linestyle="--", linewidth=1.0, alpha=0.6)
    ax.set_xlabel("threshold")
    ax.set_ylabel("count")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Layer 1 x Layer 2 concordance
# ---------------------------------------------------------------------------


def plot_layer_concordance(
    node_stats: pd.DataFrame,
    weekly_growth: pd.DataFrame,
    *,
    log_axes: bool = True,
    title: str = "Layer-1 vs Layer-2 SSE evidence per meta-cluster",
    width: str = "slide",
    height_in: float = 5.0,
) -> plt.Figure:
    """For each meta-cluster, plot Layer-1 candidate count vs Layer-2 SSE-week count.

    Meta-clusters on the diagonal are agreed; off-diagonal cases are the
    interpretively interesting failure modes.
    """
    if "meta_cluster_id" not in node_stats.columns:
        raise KeyError("node_stats needs 'meta_cluster_id'")

    layer1 = (
        node_stats.groupby("meta_cluster_id")["sse_candidate"]
        .sum()
        .rename("n_layer1_candidates")
    )
    layer2 = (
        weekly_growth.loc[weekly_growth["is_sse"]]
        .groupby("meta_cluster_id")
        .size()
        .rename("n_layer2_sse_weeks")
    )
    joined = pd.concat([layer1, layer2], axis=1).fillna(0)
    joined["agree"] = (joined["n_layer1_candidates"] > 0) & (joined["n_layer2_sse_weeks"] > 0)
    joined["layer1_only"] = (joined["n_layer1_candidates"] > 0) & (joined["n_layer2_sse_weeks"] == 0)
    joined["layer2_only"] = (joined["n_layer1_candidates"] == 0) & (joined["n_layer2_sse_weeks"] > 0)

    fig, ax = style.new_figure(width=width, height_in=height_in, context="talk", font_scale=0.82)
    masks = [
        ("neither", (joined["n_layer1_candidates"] == 0) & (joined["n_layer2_sse_weeks"] == 0), "#dddddd"),
        ("layer 1 only", joined["layer1_only"], "#3A6EA5"),
        ("layer 2 only", joined["layer2_only"], "#7C8A43"),
        ("both", joined["agree"], "#C75C2C"),
    ]
    for label, mask, color in masks:
        sub = joined.loc[mask]
        if sub.empty:
            continue
        jitter_x = np.random.default_rng(0).uniform(-0.18, 0.18, len(sub))
        jitter_y = np.random.default_rng(1).uniform(-0.18, 0.18, len(sub))
        ax.scatter(
            sub["n_layer1_candidates"] + jitter_x,
            sub["n_layer2_sse_weeks"] + jitter_y,
            s=12, alpha=0.55, color=color, edgecolor="none", label=label,
        )

    if log_axes:
        ax.set_xscale("symlog", linthresh=1)
        ax.set_yscale("symlog", linthresh=1)

    ax.set_xlabel("Layer-1 SSE candidate nodes per meta-cluster")
    ax.set_ylabel("Layer-2 SSE weeks per meta-cluster")
    ax.set_title(title)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Spatial / demographic
# ---------------------------------------------------------------------------


def plot_simd_breakdown(
    node_stats: pd.DataFrame,
    *,
    quintile_col: str = "top_simd_quintiles",
    title: str = "SIMD-entropy z-score: candidates vs background",
    width: str = "slide",
    height_in: float = 4.4,
) -> plt.Figure:
    """Compare SIMD entropy z-score distributions for candidates vs background.

    The detector's ``simd_entropy_z`` column captures whether a cluster is
    more or less SIMD-diverse than expected for its size and window. This
    figure asks whether SSE candidates concentrate in deprivation-mixed
    clusters or not.
    """
    score = "simd_entropy_z"
    if score not in node_stats.columns:
        raise KeyError(f"node_stats needs '{score}'")

    fig, axes = style.new_figure(
        width=width, height_in=height_in, ncols=2,
        context="talk", font_scale=0.78,
    )

    ax = axes[0]
    for label, color, mask in [
        ("background", "#8C8C8C", ~node_stats["sse_candidate"]),
        ("candidate", "#C75C2C", node_stats["sse_candidate"]),
    ]:
        values = node_stats.loc[mask, score].dropna().to_numpy()
        if len(values) < 5:
            continue
        sns.kdeplot(
            values, ax=ax, fill=True, color=color, alpha=0.35,
            linewidth=1.2, label=label, common_norm=False,
        )
    ax.axvline(0, color="#14151F", linestyle="--", linewidth=1.0, alpha=0.5)
    ax.axvline(1.96, color="#C75C2C", linestyle=":", linewidth=1.0, alpha=0.7)
    ax.axvline(-1.96, color="#C75C2C", linestyle=":", linewidth=1.0, alpha=0.7)
    ax.set_xlabel("SIMD entropy z-score")
    ax.set_ylabel("density")
    ax.set_title("Entropy z-score")
    ax.legend(loc="upper left", frameon=False)

    ax = axes[1]
    if quintile_col in node_stats.columns:
        def _modal_quintile(s: str) -> str | float:
            if not isinstance(s, str) or not s:
                return np.nan
            return s.split(";")[0].split(" ")[0]

        mq = node_stats[quintile_col].map(_modal_quintile)
        df = pd.DataFrame({"q": mq, "candidate": node_stats["sse_candidate"]})
        df = df.dropna(subset=["q"])
        share = (
            df.groupby(["q", "candidate"])
            .size()
            .rename("n")
            .reset_index()
        )
        denom = share.groupby("candidate")["n"].transform("sum")
        share["frac"] = share["n"] / denom
        order = sorted(share["q"].unique())
        bar_w = 0.4
        x = np.arange(len(order))
        for off, (lbl, color) in zip(
            [-bar_w / 2, bar_w / 2],
            [("background", "#8C8C8C"), ("candidate", "#C75C2C")],
        ):
            sub = share.loc[share["candidate"] == (lbl == "candidate")].set_index("q").reindex(order, fill_value=0)
            ax.bar(x + off, sub["frac"].to_numpy(), bar_w, color=color, label=lbl)
        ax.set_xticks(x)
        ax.set_xticklabels(order)
        ax.set_xlabel("modal SIMD quintile (1 = most deprived)")
        ax.set_ylabel("fraction of nodes")
        ax.set_title("Modal quintile composition")
        ax.legend(loc="upper right", frameon=False)
    else:
        ax.set_visible(False)

    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Methods / robustness
# ---------------------------------------------------------------------------


def plot_null_comparison(
    node_stats: pd.DataFrame,
    *,
    obs_col: str = "downstream_entropy_norm",
    null_mean_col: str = "downstream_entropy_norm_null_mean",
    z_col: str = "downstream_entropy_norm_z",
    title: str | None = None,
    width: str = "slide",
    height_in: float = 4.6,
) -> plt.Figure:
    """Observed metric vs null mean, coloured by significance.

    Confirms the null distribution is doing what it claims and shows where
    the observed metric departs from it.
    """
    cols = [obs_col, null_mean_col]
    df = node_stats[[c for c in cols if c in node_stats.columns]].copy()
    if df.shape[1] < 2:
        raise KeyError(
            f"Missing required columns. Need {obs_col!r} and {null_mean_col!r}."
        )
    df = df.dropna()
    if z_col in node_stats.columns:
        df["_z"] = node_stats.loc[df.index, z_col]
    else:
        df["_z"] = np.nan

    fig, ax = style.new_figure(width=width, height_in=height_in, context="talk", font_scale=0.82)
    lo = min(df[obs_col].min(), df[null_mean_col].min())
    hi = max(df[obs_col].max(), df[null_mean_col].max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=1.0, alpha=0.5)
    if df["_z"].notna().any():
        sc = ax.scatter(
            df[null_mean_col], df[obs_col],
            c=df["_z"], cmap="coolwarm", vmin=-3, vmax=3,
            s=10, alpha=0.55, edgecolor="none",
        )
        fig.colorbar(sc, ax=ax, label="z-score")
    else:
        ax.scatter(
            df[null_mean_col], df[obs_col],
            s=10, alpha=0.55, color="#3A6EA5", edgecolor="none",
        )
    ax.set_xlabel(f"null mean ({null_mean_col})")
    ax.set_ylabel(f"observed ({obs_col})")
    ax.set_title(title or f"Observed vs null: {obs_col}")
    fig.tight_layout()
    return fig
