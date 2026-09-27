"""Chapter 6 descriptive figures and tables; uses saved fits and detector labels."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import FuncFormatter

from utils import load_analysis_columns
from utils.latex_tables import latex_escape

from ..characterisation import (
    ATTRIBUTES,
    ENTROPY_ATTRIBUTES,
    SIZE_BANDS,
    build_summaries,
    validate_summaries,
)
from ..sse.config import TRANSITION_WINDOW_STRIDE
from ..sse.io import write_table
from .common import (
    Paths,
    add_common_args,
    new_figure,
    paths_from_args,
    read_table,
    styled_save_figure,
)

BUILDER_NAME = "ch6_size_profiles"
COLORS = {"CSC": "#0072B2", "Background": "#D55E00"}
SCENARIO_ORDER = ("p_0.010", "p_0.025", "primary", "p_0.100", "size_10", "size_20")
SCENARIO_LABELS = ("p≤0.01", "p≤0.025", "Primary", "p≤0.10", "Size≥10", "Size≥20")
CATEGORY_ORDERS = {
    "sex": ["Female", "Male"],
    "age": ["00-04", "05-14", "15-24", "25-64", "65-74", "75+"],
    "simd": ["1", "2", "3", "4", "5"],
    "urban_rural": [
        "Large Urban Areas",
        "Other Urban Areas",
        "Accessible Small Towns",
        "Remote Small Towns",
        "Accessible Rural",
        "Remote Rural",
    ],
}
SHORT_LABELS = {
    "Large Urban Areas": "Large urban",
    "Other Urban Areas": "Other urban",
    "Accessible Small Towns": "Accessible towns",
    "Remote Small Towns": "Remote towns",
    "Accessible Rural": "Accessible rural",
    "Remote Rural": "Remote rural",
    "Greater Glasgow and Clyde": "Glasgow and Clyde",
    "Dumfries and Galloway": "Dumfries and Galloway",
}


def prepare(paths: Paths) -> dict[str, pd.DataFrame]:
    nodes = read_table(paths, "cluster_table")
    columns = [
        "cluster_id",
        "window_id",
        "sequence_id",
        *[column for _, column, _ in ATTRIBUTES],
    ]
    records = load_analysis_columns(columns, window_stride=TRANSITION_WINDOW_STRIDE)
    tables = build_summaries(nodes, records)
    validate_summaries(tables, scotland_baseline=True)
    for key, table in tables.items():
        write_table(table, paths.result_table_dir, f"tab_ch6_{key}")
    return tables


def load_saved_summaries(paths: Paths) -> dict[str, pd.DataFrame]:
    """Load existing aggregates for presentation-only rendering."""
    names = (
        "size_summary", "size_cumulative", "composition", "composition_contrasts", "entropy"
    )
    tables = {name: read_table(paths, f"tab_ch6_{name}") for name in names}
    validate_summaries(tables, scotland_baseline=True)
    return tables


def _finish(fig, paths, name):
    outputs = styled_save_figure(fig, paths, name)
    plt.close(fig)
    return outputs


def plot_cumulative(tables, paths):
    fig, axes = new_figure(width="double", nrows=2, ncols=2, height_in=5.4)
    columns = (
        "cumulative_clusters",
        "cluster_percent",
        "cumulative_records",
        "record_percent",
    )
    labels = (
        "Cumulative clusters",
        "Cumulative clusters (%)",
        "Cumulative sequence records",
        "Cumulative sequence records (%)",
    )
    for i, (ax, col, label) in enumerate(zip(axes.flat, columns, labels)):
        for group, block in tables["size_cumulative"].groupby("group"):
            ax.step(
                np.r_[5, block.cluster_size],
                np.r_[0, block[col]],
                where="post",
                color=COLORS[group],
                label=group,
            )
        ax.set_xscale("log")
        ax.set_xticks(
            [6, 10, 20, 50, 100, 500, 1700],
            ["6", "10", "20", "50", "100", "500", "1,700"],
        )
        ax.set_xlabel("Cluster size (sequences)")
        ax.set_ylabel(label)
        ax.set_title(chr(65 + i), loc="left", fontweight="bold")
        ax.grid(axis="y", alpha=0.18)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
    axes[0, 0].legend()
    fig.tight_layout()
    return _finish(fig, paths, "fig_ch6_size_cumulative")


def _heatmap(fig, ax, matrix, limit, title, *, annotate=False):
    cmap = plt.get_cmap("RdBu").copy()
    cmap.set_bad("#eeeeee")
    mesh = ax.imshow(
        matrix.to_numpy(float),
        cmap=cmap,
        norm=TwoSlopeNorm(0, -limit, limit),
        aspect="auto",
    )
    ax.set_yticks(
        range(len(matrix)),
        [SHORT_LABELS.get(str(x), str(x)) for x in matrix.index],
        fontsize=7,
    )
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, fontsize=8)
    ax.set_title(title, loc="left", fontsize=9, fontweight="bold")
    ax.tick_params(length=0)
    if annotate:
        for row in range(len(matrix)):
            for col in range(len(matrix.columns)):
                value = matrix.iloc[row, col]
                ax.text(
                    col,
                    row,
                    f"{value:+.1f}" if pd.notna(value) else "–",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if abs(value) > 0.65 * limit else "#222222",
                )
    return mesh


def plot_composition(tables, paths, *, weighting="records", mode="size"):
    data = tables["composition_contrasts"]
    data = data.loc[data.weighting.eq(weighting)].copy()
    if mode == "size":
        data = data.loc[data.scenario.eq("primary") & data.size_band.ne("All")]
        xcol, order, labels = "size_band", SIZE_BANDS, SIZE_BANDS
        name = (
            "fig_ch6_composition_by_size"
            if weighting == "records"
            else "fig_app_ch6_composition_by_size"
        )
    elif mode == "route":
        data = data.loc[data.scenario.eq("route")]
        xcol, order, labels = (
            "candidate_group",
            ["Burst only", "Burden only", "Both axes"],
            ["Burst\nonly", "Burden\nonly", "Both\naxes"],
        )
        name = "fig_app_ch6_route_composition"
    else:
        data = data.loc[data.scenario.isin(SCENARIO_ORDER) & data.size_band.eq("All")]
        xcol, order, labels = "scenario", SCENARIO_ORDER, SCENARIO_LABELS
        name = "fig_app_ch6_sensitivity_composition"
    fig = plt.figure(figsize=(7.2, 7.0), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 1.8])
    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1:]),
    ]
    limit = max(1, np.ceil(data.difference_pp.abs().max()))
    for i, ((key, _, label), ax) in enumerate(zip(ATTRIBUTES, axes)):
        block = data.loc[data.attribute.eq(key)]
        matrix = block.pivot(
            index="category", columns=xcol, values="difference_pp"
        ).reindex(columns=order)
        matrix = matrix.reindex(CATEGORY_ORDERS.get(key, sorted(matrix.index)))
        matrix.columns = labels
        if key == "simd":
            matrix.index = ["Q" + x for x in matrix.index]
        mesh = _heatmap(
            fig, ax, matrix, limit, f"{chr(65 + i)}  {label}", annotate=True
        )
        if mode == "sensitivity":
            ax.tick_params(axis="x", labelrotation=45)
    fig.colorbar(
        mesh,
        ax=axes,
        location="bottom",
        shrink=0.7,
        pad=0.025,
        aspect=40,
        label="Candidate minus background (percentage points)",
    )
    return _finish(fig, paths, name)


def plot_entropy(tables, paths, *, mode="size"):
    data = tables["entropy"]
    if mode == "size":
        data = data.loc[data.scenario.eq("primary") & data.size_band.ne("All")]
        xcol, order, labels = "size_band", SIZE_BANDS, SIZE_BANDS
        name = "fig_ch6_entropy_by_size"
    else:
        data = data.loc[data.scenario.isin(SCENARIO_ORDER) & data.size_band.eq("All")]
        xcol, order, labels = "scenario", SCENARIO_ORDER, SCENARIO_LABELS
        name = "fig_app_ch6_sensitivity_entropy"
    fig, axes = new_figure(width="double", nrows=6, ncols=2, height_in=9.2)
    for row, (key, _, label) in enumerate(ENTROPY_ATTRIBUTES):
        for col, scale in enumerate(("obs", "z")):
            ax = axes[row, col]
            for offset, group in [(-0.12, "Background"), (0.12, "CSC")]:
                block = (
                    data.loc[
                        data.attribute.eq(key)
                        & data.scale.eq(scale)
                        & data.group.eq(group)
                    ]
                    .set_index(xcol)
                    .reindex(order)
                )
                x = np.arange(len(order)) + offset
                ax.errorbar(
                    x,
                    block["median"],
                    yerr=np.vstack(
                        [block["median"] - block.q25, block.q75 - block["median"]]
                    ),
                    fmt="o",
                    markersize=3,
                    capsize=2,
                    lw=1,
                    color=COLORS[group],
                    label=group,
                )
            ax.set_title(
                label + (" · observed" if scale == "obs" else " · null-standardised"),
                fontsize=8,
                loc="left",
            )
            ax.set_xticks(
                range(len(order)), labels if row == 5 else [""] * len(order), fontsize=7
            )
            if mode != "size" and row == 5:
                ax.tick_params(axis="x", labelrotation=40)
            if scale == "obs":
                ax.set_ylim(-0.03, 1.03)
            else:
                ax.axhline(0, color="#777777", ls=":", lw=0.7)
            ax.grid(axis="y", alpha=0.18)
            ax.tick_params(axis="y", labelsize=9)
    if mode == "size":
        handles, legend_labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(
            handles, legend_labels, loc="lower center",
            bbox_to_anchor=(0.5, 0.005), fontsize=9, ncol=2, frameon=False,
        )
        fig.tight_layout(h_pad=0.8, rect=(0, 0.04, 1, 1))
    else:
        axes[0, 0].legend(loc="lower left", fontsize=8, ncol=2)
        fig.tight_layout(h_pad=0.8)
    return _finish(fig, paths, name)


def plot_route_entropy(tables, paths):
    data = tables["entropy"].loc[lambda t: t.scenario.eq("route")]
    groups = [
        "Background",
        "Burden-eligible background",
        "Burst only",
        "Burden only",
        "Both axes",
    ]
    counts = tables["size_summary"].loc[lambda t: t.scenario.eq("route")].set_index("group").n_clusters
    short = [f"{g.replace('Burden-eligible background', 'Burden background')} ({counts[g]:,})" for g in groups]
    colors = ["#888888", "#444444", "#D55E00", "#0072B2", "#009E73"]
    fig, axes = new_figure(width="double", nrows=6, ncols=2, height_in=9.4)
    for row, (key, _, label) in enumerate(ENTROPY_ATTRIBUTES):
        for col, scale in enumerate(("obs", "z")):
            ax = axes[row, col]
            block = (
                data.loc[data.attribute.eq(key) & data.scale.eq(scale)]
                .set_index("group")
                .reindex(groups)
            )
            for y, (group, color) in enumerate(zip(groups, colors)):
                b = block.loc[group]
                ax.errorbar(
                    b["median"],
                    y,
                    xerr=[[b["median"] - b.q25], [b.q75 - b["median"]]],
                    fmt="o",
                    markersize=3,
                    capsize=2,
                    color=color,
                )
            ax.set_yticks(range(5), short if col == 0 else [""] * 5, fontsize=8)
            ax.invert_yaxis()
            ax.set_title(
                label + (" · observed" if scale == "obs" else " · null-standardised"),
                fontsize=8,
                loc="left",
            )
            ax.tick_params(axis="x", labelsize=9)
            if scale == "obs":
                ax.set_xlim(-0.03, 1.03)
            else:
                ax.axvline(0, color="#777777", ls=":", lw=0.7)
    fig.tight_layout(h_pad=0.8)
    return _finish(fig, paths, "fig_app_ch6_route_entropy")


def build_figures(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    outputs = [
        plot_cumulative(tables, paths),
        plot_composition(tables, paths),
        plot_entropy(tables, paths),
        plot_composition(tables, paths, weighting="clusters"),
        plot_composition(tables, paths, mode="route"),
        plot_route_entropy(tables, paths),
        plot_composition(tables, paths, mode="sensitivity"),
        plot_entropy(tables, paths, mode="sensitivity"),
    ]
    return {"outputs": outputs}


def _table(
    path: Path,
    caption: str,
    label: str,
    headers: list[str],
    rows: list[list],
    *,
    long=False,
):
    path.parent.mkdir(parents=True, exist_ok=True)
    spec = (
        "@{}"
        + "l" * min(2, len(headers))
        + "r" * (len(headers) - min(2, len(headers)))
        + "@{}"
    )
    header = " & ".join(latex_escape(x) for x in headers) + r" \\"
    body = [" & ".join(latex_escape(x) for x in row) + r" \\" for row in rows]
    if long:
        lines = [
            r"\begingroup\small",
            r"\setlength{\tabcolsep}{4pt}",
            r"\begin{longtable}{" + spec + "}",
            r"\caption{" + caption + r"}\label{" + label + r"}\\",
            r"\toprule",
            header,
            r"\midrule\endfirsthead",
            r"\multicolumn{"
            + str(len(headers))
            + r"}{l}{\tablename\ \thetable\ continued}\\",
            r"\toprule",
            header,
            r"\midrule\endhead",
            r"\bottomrule\endfoot",
            *body,
            r"\end{longtable}\endgroup",
        ]
    else:
        lines = [
            r"\begin{table}[htbp]\centering",
            r"\caption{" + caption + r"}\label{" + label + "}",
            r"\begin{thesistablebody}{" + spec + "}",
            r"\toprule",
            header,
            r"\midrule",
            *body,
            r"\bottomrule\end{thesistablebody}",
            r"\end{table}",
        ]
    path.write_text("\n".join(lines) + "\n")


def write_size_composition_table(directory: Path, comp: pd.DataFrame) -> Path:
    """Lay out saved percentages by category, with three columns per size band."""
    path = directory / "tab_app_ch6_size_composition.tex"
    directory.mkdir(parents=True, exist_ok=True)
    band_labels = [band.replace("-", "--").replace("100+", r"$\geq100$") for band in SIZE_BANDS]
    header = [
        "Category & " + " & ".join(r"\multicolumn{3}{c}{" + band + "}" for band in band_labels) + r" \\",
        "".join(r"\cmidrule(lr){" + f"{2 + 3*i}-{4 + 3*i}" + "}" for i in range(5)),
        " & " + " & ".join([r"CSC & BG & $\Delta$"] * 5) + r" \\",
        r"\midrule",
    ]
    lines = [
        r"\begin{landscape}\begingroup",
        r"\renewcommand{\thesistablesetup}{\small\setstretch{1}\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.0}}",
        r"\begin{longtable}{@{}p{0.22\linewidth}*{15}{r}@{}}",
        r"\caption[Composition within size bands]{Composition within size bands. CSC and BG give the percentage of sequence records in each category for CSC and background clusters; $\Delta$ is CSC minus background in percentage points. Denominators are the record counts in Table~\ref{tab:ch6_size_overview}; all five attributes are complete.}\label{tab:app_ch6_size_composition}\\",
        r"\toprule", *header, r"\endfirsthead",
        r"\multicolumn{16}{l}{\tablename\ \thetable\ continued}\\",
        r"\toprule", *header, r"\endhead",
        r"\midrule\multicolumn{16}{r}{\small\itshape Continued on next page}\\\endfoot",
        r"\bottomrule\endlastfoot",
    ]
    for attr_index, (key, _, label) in enumerate(ATTRIBUTES):
        b = comp.loc[
            comp.scenario.eq("primary") & comp.weighting.eq("records")
            & comp.attribute.eq(key) & comp.size_band.ne("All")
        ]
        categories = CATEGORY_ORDERS.get(key, sorted(b.category.unique()))
        if key == "health_board":
            lines.append(r"\newpage")
        elif attr_index:
            lines.append(r"\addlinespace[0.4em]")
        lines.append(r"\multicolumn{16}{l}{\textbf{" + label + r"}}\\*")
        for cat in categories:
            display = "Q" + cat if key == "simd" else SHORT_LABELS.get(cat, cat)
            values = [display]
            for band in SIZE_BANDS:
                r = b.loc[b.category.eq(cat) & b.size_band.eq(band)].iloc[0]
                values.extend([f"{r.candidate_percent:.1f}", f"{r.background_percent:.1f}", f"{r.difference_pp:+.1f}"])
            lines.append(" & ".join(values) + r" \\")
    lines.extend([r"\end{longtable}\endgroup\end{landscape}"])
    path.write_text("\n".join(lines) + "\n")
    return path


def build_tables(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    # The central table builder passes its publication table directory as figure_dir.
    directory = paths.figure_dir
    size = tables["size_summary"]
    rows = []
    for band in ("All", *SIZE_BANDS):
        for group in ("CSC", "Background"):
            r = size.loc[
                size.scenario.eq("primary")
                & size.size_band.eq(band)
                & size.group.eq(group)
            ].iloc[0]
            rows.append(
                [
                    band,
                    group,
                    f"{r.n_clusters:,}",
                    f"{r.n_records:,}",
                    f"{r.size_median:g} ({r.size_q25:g}--{r.size_q75:g})",
                ]
            )
    _table(
        directory / "tab_ch6_size_overview.tex",
        "Cluster sizes and sequence records in the primary comparison. Records count sequence appearances in retained windows, not unique people. The interquartile range (IQR) describes the middle half of cluster sizes.",
        "tab:ch6_size_overview",
        ["Size band", "Group", "Clusters", "Records", "Median (IQR)"],
        rows,
    )
    rows = []
    for scenario, label in zip(SCENARIO_ORDER, SCENARIO_LABELS):
        for group in ("CSC", "Background"):
            r = size.loc[
                size.scenario.eq(scenario)
                & size.size_band.eq("All")
                & size.group.eq(group)
            ].iloc[0]
            rows.append(
                [
                    label.replace("≤", " <= ").replace("≥", " >= "),
                    group,
                    f"{r.n_clusters:,}",
                    f"{r.n_records:,}",
                    f"{r.size_median:g}",
                ]
            )
    _table(
        directory / "tab_ch6_sensitivity_sizes.tex",
        "Descriptive sensitivity comparisons using saved scores. Primary uses $p\\leq0.05$ and size $\\geq6$. Alternative p cutoffs retain the same testing floor; size restrictions retain primary labels without recalibration.",
        "tab:ch6_sensitivity_sizes",
        ["Comparison", "Group", "Clusters", "Records", "Median size"],
        rows,
    )
    rows = []
    for r in size.loc[size.scenario.eq("route")].itertuples():
        rows.append(
            [
                r.group.replace("Burden-eligible background", "Burden background"),
                f"{r.n_clusters:,}",
                f"{r.n_records:,}",
                f"{r.size_median:g} ({r.size_q25:g}--{r.size_q75:g})",
            ]
        )
    _table(
        directory / "tab_ch6_route_sizes.tex",
        "Cluster sizes by screening route. Burden background is the subset of background clusters with positive observed downstream burden; it is the comparator for burden-only and both-axis candidates. Rows therefore do not form a disjoint partition.",
        "tab:ch6_route_sizes",
        ["Group", "Clusters", "Records", "Median (IQR)"],
        rows,
    )
    write_size_composition_table(directory, tables["composition_contrasts"])
    return {"outputs": list(directory.glob("*ch6*size*.tex"))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--tables-only", action="store_true")
    parser.add_argument(
        "--from-saved-tables", action="store_true",
        help="Render existing tab_ch6_* aggregates without rebuilding or writing them.",
    )
    args = parser.parse_args()
    paths = paths_from_args(args)
    tables = load_saved_summaries(paths) if args.from_saved_tables else None
    if args.tables_only:
        paths = Paths(
            table_dir=paths.table_dir,
            figure_dir=paths.figure_dir.parent / "tables",
            result_table_dir=paths.result_table_dir,
        )
        build_tables(paths, tables=tables)
    else:
        build_figures(paths, tables=tables)


if __name__ == "__main__":
    main()
