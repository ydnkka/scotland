"""Observed-vs-expected socioeconomic and demographic mixing matrices."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from .cluster_outcome_models import analysis_dataset_path, repo_root
    from .wave_specific_domain_demographic_mixing import assign_wave, WAVE_LABELS, WAVE_ORDER
except ImportError:
    from cluster_outcome_models import analysis_dataset_path, repo_root
    from wave_specific_domain_demographic_mixing import assign_wave, WAVE_LABELS, WAVE_ORDER


QC_DEFAULT = "good"

VARIABLES = {
    "simd": {
        "column": "dz_simd_quintile",
        "label": "SIMD quintile",
        "levels": [1, 2, 3, 4, 5],
    },
    "age": {
        "column": "age_band",
        "label": "Age band",
        "levels": [
            "00-04",
            "05-09",
            "10-14",
            "15-19",
            "20-24",
            "25-29",
            "30-34",
            "35-39",
            "40-44",
            "45-49",
            "50-54",
            "55-59",
            "60-64",
            "65-69",
            "70-74",
            "75+",
        ],
    },
}

SEQUENCE_COLUMNS = [
    "cluster_id",
    "sequence_id",
    "resolution",
    "window_id",
    "pango_lineage",
    "nextclade_qc",
    "dz_simd_quintile",
    "age_band",
]


def read_sequence_rows(root: Path, qc: str | None) -> pd.DataFrame:
    filters = None if qc is None else [("nextclade_qc", "==", qc)]
    df = pd.read_parquet(
        analysis_dataset_path(root),
        columns=SEQUENCE_COLUMNS,
        filters=filters,
        engine="pyarrow",
    )
    for col in ["cluster_id", "sequence_id", "window_id", "pango_lineage", "nextclade_qc", "age_band"]:
        df[col] = df[col].astype("category")
    df["resolution_label"] = df["resolution"].map(lambda x: f"{x:.1f}")
    df["wave_group"] = df["pango_lineage"].astype(str).map(assign_wave).astype("category")
    return df[df["wave_group"].isin(WAVE_ORDER)].copy()


def observed_ordered_pairs(
    cluster_counts: pd.DataFrame,
    levels: list,
) -> pd.DataFrame:
    wide = (
        cluster_counts
        .pivot_table(
            index=["cluster_id", "wave_group", "window_id", "pango_lineage", "resolution_label"],
            columns="category",
            values="n",
            aggfunc="sum",
            fill_value=0,
            observed=True,
        )
        .reindex(columns=levels, fill_value=0)
        .astype(np.int32)
        .reset_index()
    )
    wide["n_valid"] = wide[levels].sum(axis=1)
    wide = wide[wide["n_valid"] >= 2].copy()

    rows = []
    for left in levels:
        for right in levels:
            values = wide[left].astype(np.int64) * wide[right].astype(np.int64)
            if left == right:
                values = wide[left].astype(np.int64) * (wide[left].astype(np.int64) - 1)
            by_wave = values.groupby(wide["wave_group"], observed=True).sum()
            for wave, n_pairs in by_wave.items():
                rows.append(
                    {
                        "wave_group": wave,
                        "category_i": left,
                        "category_j": right,
                        "observed_pairs": float(n_pairs),
                    }
                )
    return pd.DataFrame(rows), wide


def expected_ordered_pairs(
    cluster_wide: pd.DataFrame,
    stratum_counts: pd.DataFrame,
    levels: list,
) -> pd.DataFrame:
    cluster_wide = cluster_wide.copy()
    cluster_wide["ordered_pairs"] = cluster_wide["n_valid"] * (cluster_wide["n_valid"] - 1)
    stratum_cols = ["wave_group", "window_id", "pango_lineage", "resolution_label"]
    stratum_pair_totals = (
        cluster_wide.groupby(stratum_cols, observed=True)["ordered_pairs"]
        .sum()
        .rename("cluster_ordered_pairs")
        .reset_index()
    )

    stratum_wide = (
        stratum_counts
        .pivot_table(
            index=stratum_cols,
            columns="category",
            values="n",
            aggfunc="sum",
            fill_value=0,
            observed=True,
        )
        .reindex(columns=levels, fill_value=0)
        .astype(np.int64)
        .reset_index()
    )
    stratum_wide["stratum_n"] = stratum_wide[levels].sum(axis=1)
    stratum_wide = stratum_wide.merge(stratum_pair_totals, on=stratum_cols, how="inner")
    denom = stratum_wide["stratum_n"] * (stratum_wide["stratum_n"] - 1)

    rows = []
    for left in levels:
        for right in levels:
            numerator = stratum_wide[left].astype(np.float64) * stratum_wide[right].astype(np.float64)
            if left == right:
                numerator = stratum_wide[left].astype(np.float64) * (
                    stratum_wide[left].astype(np.float64) - 1
                )
            expected = stratum_wide["cluster_ordered_pairs"] * numerator / denom
            expected = expected.replace([np.inf, -np.inf], np.nan).fillna(0)
            by_wave = expected.groupby(stratum_wide["wave_group"], observed=True).sum()
            for wave, n_pairs in by_wave.items():
                rows.append(
                    {
                        "wave_group": wave,
                        "category_i": left,
                        "category_j": right,
                        "expected_pairs": float(n_pairs),
                    }
                )
    return pd.DataFrame(rows)


def build_matrix_for_variable(seq: pd.DataFrame, variable: str) -> pd.DataFrame:
    spec = VARIABLES[variable]
    levels = spec["levels"]
    work = seq.dropna(subset=[spec["column"]]).copy()
    work["category"] = work[spec["column"]]

    cluster_counts = (
        work.groupby(
            ["cluster_id", "wave_group", "window_id", "pango_lineage", "resolution_label", "category"],
            observed=True,
        )
        .size()
        .rename("n")
        .reset_index()
    )
    observed, cluster_wide = observed_ordered_pairs(cluster_counts, levels)

    stratum_counts = (
        work.groupby(
            ["wave_group", "window_id", "pango_lineage", "resolution_label", "category"],
            observed=True,
        )
        .size()
        .rename("n")
        .reset_index()
    )
    expected = expected_ordered_pairs(cluster_wide, stratum_counts, levels)

    matrix = observed.merge(
        expected,
        on=["wave_group", "category_i", "category_j"],
        how="outer",
    ).fillna({"observed_pairs": 0, "expected_pairs": 0})
    matrix["variable"] = variable
    matrix["variable_label"] = spec["label"]

    overall = (
        matrix.groupby(["variable", "variable_label", "category_i", "category_j"], observed=True)[
            ["observed_pairs", "expected_pairs"]
        ]
        .sum()
        .reset_index()
    )
    overall["wave_group"] = "Overall"
    matrix = pd.concat([matrix, overall], ignore_index=True)

    totals = (
        matrix.groupby(["variable", "wave_group"], observed=True)[["observed_pairs", "expected_pairs"]]
        .sum()
        .rename(columns={"observed_pairs": "total_observed_pairs", "expected_pairs": "total_expected_pairs"})
        .reset_index()
    )
    matrix = matrix.merge(totals, on=["variable", "wave_group"], how="left")
    matrix["observed_probability"] = matrix["observed_pairs"] / matrix["total_observed_pairs"]
    matrix["expected_probability"] = matrix["expected_pairs"] / matrix["total_expected_pairs"]
    matrix["excess_probability"] = matrix["observed_probability"] - matrix["expected_probability"]
    matrix["excess_percentage_points"] = matrix["excess_probability"] * 100
    matrix["observed_expected_ratio"] = matrix["observed_probability"] / matrix["expected_probability"]
    matrix["wave_label"] = matrix["wave_group"].map(lambda w: "Overall" if w == "Overall" else WAVE_LABELS.get(w, w))
    return matrix[
        [
            "variable",
            "variable_label",
            "wave_group",
            "wave_label",
            "category_i",
            "category_j",
            "observed_pairs",
            "expected_pairs",
            "observed_probability",
            "expected_probability",
            "excess_probability",
            "excess_percentage_points",
            "observed_expected_ratio",
        ]
    ]


def plot_heatmap(matrix: pd.DataFrame, variable: str, wave: str, out_base: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    levels = VARIABLES[variable]["levels"]
    sub = matrix[(matrix["variable"] == variable) & (matrix["wave_group"] == wave)]
    sub = sub.assign(
        category_i=sub["category_i"].astype(str),
        category_j=sub["category_j"].astype(str),
    )
    level_keys = [str(level) for level in levels]
    values = (
        sub.pivot(index="category_i", columns="category_j", values="excess_percentage_points")
        .reindex(index=level_keys, columns=level_keys)
    )
    vmax = max(0.5, float(np.nanmax(np.abs(values.to_numpy()))))

    width = 4.2 if variable == "simd" else 6.6
    fig, ax = plt.subplots(figsize=(width, width))
    im = ax.imshow(values.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(levels)))
    ax.set_yticks(range(len(levels)))
    ax.set_xticklabels(levels, rotation=45, ha="right")
    ax.set_yticklabels(levels)
    ax.set_xlabel(VARIABLES[variable]["label"])
    ax.set_ylabel(VARIABLES[variable]["label"])
    ax.set_title(f"{VARIABLES[variable]['label']}: observed - expected ({wave})")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Excess pair probability (pp)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_simd_by_wave(matrix: pd.DataFrame, out_base: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/scotland-mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/scotland-xdg-cache")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    levels = VARIABLES["simd"]["levels"]
    level_keys = [str(level) for level in levels]
    waves = [w for w in WAVE_ORDER if w in set(matrix["wave_group"])]
    vmax = max(
        0.5,
        float(
            np.nanmax(
                np.abs(
                    matrix[
                        (matrix["variable"] == "simd")
                        & (matrix["wave_group"].isin(waves))
                    ]["excess_percentage_points"]
                )
            )
        ),
    )

    fig, axes = plt.subplots(3, 3, figsize=(8, 8), sharex=True, sharey=True)
    for ax, wave in zip(axes.flat, waves):
        sub = matrix[(matrix["variable"] == "simd") & (matrix["wave_group"] == wave)]
        sub = sub.assign(
            category_i=sub["category_i"].astype(str),
            category_j=sub["category_j"].astype(str),
        )
        values = (
            sub.pivot(index="category_i", columns="category_j", values="excess_percentage_points")
            .reindex(index=level_keys, columns=level_keys)
        )
        im = ax.imshow(values.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(WAVE_LABELS.get(wave, wave))
        ax.set_xticks(range(len(levels)))
        ax.set_yticks(range(len(levels)))
        ax.set_xticklabels(levels)
        ax.set_yticklabels(levels)

    for ax in axes.flat[len(waves):]:
        ax.axis("off")

    cbar = fig.colorbar(im, ax=axes, shrink=0.7)
    cbar.set_label("Excess pair probability (pp)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.supxlabel("SIMD quintile")
    fig.supylabel("SIMD quintile")
    fig.subplots_adjust(left=0.1, right=0.86, bottom=0.08, top=0.94, wspace=0.12, hspace=0.28)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def run(root: Path, qc: str | None) -> None:
    tables_dir = root / "part1" / "tables"
    figures_dir = root / "part1" / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("Reading sequence rows", flush=True)
    seq = read_sequence_rows(root, qc)

    outputs = []
    for variable in VARIABLES:
        print(f"Building observed/expected matrix for {variable}", flush=True)
        outputs.append(build_matrix_for_variable(seq, variable))

    matrix = pd.concat(outputs, ignore_index=True)
    matrix.to_csv(tables_dir / "observed_expected_mixing_matrices.csv", index=False)

    plot_heatmap(matrix, "simd", "Overall", figures_dir / "observed_expected_simd_matrix_overall")
    plot_heatmap(matrix, "age", "Overall", figures_dir / "observed_expected_age_matrix_overall")
    plot_simd_by_wave(matrix, figures_dir / "observed_expected_simd_matrix_by_wave")

    print(f"Wrote {tables_dir / 'observed_expected_mixing_matrices.csv'}", flush=True)
    print(f"Wrote {figures_dir / 'observed_expected_simd_matrix_overall.png'}", flush=True)
    print(f"Wrote {figures_dir / 'observed_expected_age_matrix_overall.png'}", flush=True)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument("--qc", default=QC_DEFAULT)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    qc = None if str(args.qc).lower() == "none" else str(args.qc)
    run(args.root.resolve(), qc)


if __name__ == "__main__":
    main()
