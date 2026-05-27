"""Render clade-stratified association significance figures.

The figures use the clade-sensitivity Wald tables and show where candidate-node
composition and internal-mixing predictors are significant within clade groups.
"""

from __future__ import annotations

import sys
import re
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils import CLADES, new_figure, save_figure, add_panel_labels  # noqa: E402


RESULTS_DIR = PROJECT_ROOT / "sse_detection" / "results" / "sensitivity_clade"
FIG_DIR = PROJECT_ROOT / "sse_detection" / "figures"
ANALYSIS_DATASET = (
    PROJECT_ROOT / "data" / "processed" / "scotland_clustering_analysis_dataset.parquet"
)

ALPHA = 0.05
SCORE_CAP = 6.0  # Cap -log10(q) to avoid extreme values dominating the color scale in the heatmaps.

COMPOSITION_ROWS = [
    ("sex", "Sex"),
    ("age_band", "Age band"),
    ("simd_quintile", "SIMD quintile"),
    ("urban_rural_class", "Urban/rural class"),
    ("health_board", "Health board"),
]

MIXING_ROWS = [
    ("sex", "Sex entropy"),
    ("age", "Age entropy"),
    ("simd", "SIMD entropy"),
    ("urban_rural", "Urban/rural entropy"),
    ("health_board", "Health-board entropy"),
]

MODEL_PANELS = [
    ("primary_single", "Primary single-predictor", "primary", "single"),
    ("primary_joint", "Primary joint", "primary", "joint"),
    ("expanded_single", "Expanded single-predictor", "expanded", "single"),
    ("expanded_joint", "Expanded joint", "expanded", "joint"),
]

COMPOSITION_EFFECT_SPECS = [
    (
        "age_band",
        "Age-band composition odds ratios by clade group",
        "fig17_clade_composition_or_age_band",
        [
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
    ),
    (
        "sex",
        "Sex composition odds ratios by clade group",
        "fig18_clade_composition_or_sex",
        ["Male", "Female"],
    ),
    (
        "simd_quintile",
        "SIMD-quintile composition odds ratios by clade group",
        "fig19_clade_composition_or_simd_quintile",
        ["1", "2", "3", "4", "5"],
    ),
    (
        "urban_rural_class",
        "Urban-rural composition odds ratios by clade group",
        "fig20_clade_composition_or_urban_rural_class",
        [
            "Large Urban Areas",
            "Other Urban Areas",
            "Accessible Small Towns",
            "Remote Small Towns",
            "Accessible Rural",
            "Remote Rural",
        ],
    ),
    (
        "health_board",
        "Health-board composition odds ratios by clade group",
        "fig21_clade_composition_or_health_board",
        [
            "Greater Glasgow and Clyde",
            "Ayrshire and Arran",
            "Borders",
            "Dumfries and Galloway",
            "Fife",
            "Forth Valley",
            "Grampian",
            "Highland",
            "Lanarkshire",
            "Lothian",
            "Orkney",
            "Shetland",
            "Tayside",
            "Western Isles",
        ],
    ),
]

COMPOSITION_REFERENCE_LEVELS = {
    "sex": "Male",
    "age_band": "30-34",
    "simd_quintile": "1",
    "urban_rural_class": "Large Urban Areas",
    "health_board": "Greater Glasgow and Clyde",
}

COMPOSITION_TERM_TOKENS = {
    "sex": "C(sex,",
    "age_band": "C(age_band,",
    "simd_quintile": "C(dz_simd_quintile,",
    "urban_rural_class": "C(dz_urban_rural_class,",
    "health_board": "C(dz_health_board,",
}


def _clean_strings(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(col).strip() for col in out.columns]
    for col in out.select_dtypes(include=["object", "string"]).columns:
        present = out[col].notna()
        out.loc[present, col] = out.loc[present, col].astype(str).str.strip()
    return out


def _canonical_mixing_predictor(value: object) -> str:
    return (
        str(value)
        .replace("_entropy_z", "")
        .replace("_entropy_obs_x10", "")
        .replace("_entropy_obs", "")
    )


def _read_wald(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return _clean_strings(pd.read_csv(path, skipinitialspace=True))


def _clade_summary() -> pd.DataFrame:
    data = pd.read_parquet(
        ANALYSIS_DATASET,
        columns=["sequence_id", "collection_date", "clade"],
    ).drop_duplicates("sequence_id")
    data["collection_date"] = pd.to_datetime(data["collection_date"])
    data["clade_group"] = data["clade"].map(CLADES).fillna("Other")
    summary = (
        data.groupby("clade_group", dropna=False)
        .agg(
            n_sequences=("sequence_id", "nunique"),
            start=("collection_date", "min"),
            end=("collection_date", "max"),
            median=("collection_date", "median"),
        )
        .reset_index()
    )
    summary["_other"] = summary["clade_group"].eq("Other")
    summary = summary.sort_values(["_other", "median", "clade_group"]).drop(
        columns="_other"
    )
    return summary


def _clade_labels(summary: pd.DataFrame, clades: Iterable[str]) -> list[str]:
    summary = summary.copy()
    summary["median"] = pd.to_datetime(summary["median"])
    summary["median"] = summary["median"].apply(
        lambda dt: dt.strftime("%b %Y") if pd.notna(dt) else "N/A"
    )
    lookup = summary.set_index("clade_group")
    labels = []
    for clade in clades:
        if clade in lookup.index and pd.notna(lookup.loc[clade, "median"]):
            median = lookup.loc[clade, "median"]
            labels.append(f"{median}::{clade}")
        else:
            labels.append(str(clade))
    return labels

def _pivot_q(
    wald: pd.DataFrame,
    *,
    domain: str,
    model_set: str,
    predictor_set: str,
    rows: list[tuple[str, str]],
    clade_order: list[str],
) -> pd.DataFrame:
    data = wald.loc[
        wald["domain"].eq(domain)
        & wald["model_set"].eq(model_set)
        & wald["predictor_set"].eq(predictor_set)
    ].copy()

    row_keys = [key for key, _ in rows]
    if domain == "node_mixing":
        source_col = "term" if predictor_set == "joint" else "predictor"
        data["_row_key"] = data[source_col].map(_canonical_mixing_predictor)
    else:
        data["_row_key"] = data["predictor"]

    data = data.loc[data["_row_key"].isin(row_keys)].copy()
    matrix = data.pivot_table(
        index="_row_key",
        columns="clade_group",
        values="p_adj_bh",
        aggfunc="min",
    )
    return matrix.reindex(index=row_keys, columns=clade_order)


def _score(q_values: pd.DataFrame) -> np.ndarray:
    q = q_values.astype(float).to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        scores = -np.log10(np.clip(q, 1e-300, 1.0))
    return np.clip(scores, 0, SCORE_CAP)


def _draw_heatmap_grid(
    panels: list[tuple[str, pd.DataFrame]],
    *,
    row_labels: list[str],
    clade_labels: list[str],
    cmap: str,
    output_stem: str,
) -> None:
    fig, axes = new_figure(
        width="double",
        height_in=4,
        nrows=2,
        ncols=2,
        layout="constrained",
        sharex=True,
        sharey=True,
    )

    norm = Normalize(0, SCORE_CAP)
    im = None
    for panel_idx, (ax, (panel_title, q_values)) in enumerate(
        zip(axes.ravel(), panels)
    ):
        scores = _score(q_values)
        im = ax.imshow(scores, cmap=cmap, norm=norm, aspect="auto")
        # ax.set_title(panel_title)
        ax.set_xticks(np.arange(len(clade_labels)))
        ax.set_xticklabels(clade_labels, rotation=90, ha="center")
        ax.set_yticks(np.arange(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=10)
        ax.set_xlim(-0.5, len(clade_labels) - 0.5)
        ax.set_ylim(len(row_labels) - 0.5, -0.5)
        ax.set_xticks(np.arange(-0.5, len(clade_labels), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.2)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.tick_params(axis="both", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if panel_idx < 2:
            ax.tick_params(labelbottom=False)

        sig = q_values.astype(float).to_numpy() < ALPHA
        for row, col in np.argwhere(sig):
            ax.scatter(
                col,
                row,
                s=42,
                marker="*",
                facecolor="black",
                edgecolor="white",
                linewidth=0.7,
                zorder=3,
            )

    if im is not None:
        cbar = fig.colorbar(im, ax=axes.ravel(), shrink=0.82, pad=0.02)
        cbar.set_label("Evidence strength (-log10 p-value)")
        threshold = -np.log10(ALPHA)
        cbar.ax.axhline(threshold, color="#000000", linewidth=2.0)

    add_panel_labels(axes.ravel())
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    save_figure(
        fig,
        FIG_DIR / f"{output_stem}",
        width="double",
        save_png=True,
        save_pdf=True,
    )
    plt.close(fig)


def _summary_table(
    panels: list[tuple[str, pd.DataFrame]],
    *,
    rows: list[tuple[str, str]],
    clade_summary: pd.DataFrame,
    output_name: str,
) -> None:
    labels = dict(rows)
    pieces = []
    for predictor_set, values in panels:
        long = (
            values.rename_axis("predictor")
            .reset_index()
            .melt(id_vars="predictor", var_name="clade_group", value_name="q_value")
        )
        long["predictor_set"] = predictor_set
        pieces.append(long)
    out = pd.concat(pieces, ignore_index=True)
    out["label"] = out["predictor"].map(labels)
    out["significant"] = out["q_value"].astype(float) < ALPHA
    out = out.merge(clade_summary, on="clade_group", how="left")
    out.to_csv(FIG_DIR / output_name, index=False)


def _term_level(term: object) -> str:
    match = re.search(r"\[T\.(.*)\]$", str(term))
    return match.group(1) if match else str(term)


def _read_odds(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    out = _clean_strings(pd.read_csv(path, skipinitialspace=True))
    numeric_cols = ["odds_ratio", "or_low", "or_high", "p_value"]
    for col in numeric_cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def _composition_or_matrix(
    odds: pd.DataFrame,
    *,
    predictor: str,
    model_set: str,
    predictor_set: str,
    levels: list[str],
    clade_order: list[str],
) -> tuple[pd.DataFrame, str]:
    data = odds.loc[
        odds["domain"].eq("composition")
        & odds["model_set"].eq(model_set)
        & odds["predictor_set"].eq(predictor_set)
    ].copy()
    if predictor_set == "joint":
        token = COMPOSITION_TERM_TOKENS[predictor]
        data = data.loc[data["term"].astype(str).str.contains(token, regex=False)]
    else:
        data = data.loc[data["predictor"].eq(predictor)]
    if data.empty:
        return pd.DataFrame(index=levels, columns=clade_order, dtype=float), ""

    reference = COMPOSITION_REFERENCE_LEVELS[predictor]
    data["level"] = data["term"].map(_term_level)
    matrix = data.pivot_table(
        index="level",
        columns="clade_group",
        values="odds_ratio",
        aggfunc="first",
    ).reindex(index=levels, columns=clade_order)
    if reference in matrix.index:
        matrix.loc[reference, clade_order] = 1.0
    return matrix, reference


def _or_scores(or_values: pd.DataFrame, *, log2_cap: float = 2.0) -> np.ndarray:
    values = or_values.astype(float).to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        scores = np.log2(values)
    return np.clip(scores, -log2_cap, log2_cap)


def _draw_or_heatmap_grid(
    panels: list[tuple[str, pd.DataFrame]],
    *,
    row_labels: list[str],
    clade_labels: list[str],
    title: str,
    output_stem: str,
    log2_cap: float = 2.0,
) -> None:
    n_rows = len(row_labels)
    height_in = min(9.5, max(4.2, 3.4 + 0.28 * n_rows))
    fig, axes = new_figure(
        width="double",
        height_in=height_in,
        nrows=2,
        ncols=2,
        layout="constrained",
        sharex=True,
        sharey=True,
    )

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#E5E7EB")
    norm = Normalize(-log2_cap, log2_cap)
    im = None
    for panel_idx, (ax, (panel_title, or_values)) in enumerate(
        zip(axes.ravel(), panels)
    ):
        scores = np.ma.masked_invalid(_or_scores(or_values, log2_cap=log2_cap))
        im = ax.imshow(scores, cmap=cmap, norm=norm, aspect="auto")
        ax.set_title(panel_title)
        ax.set_xticks(np.arange(len(clade_labels)))
        ax.set_xticklabels(clade_labels, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(np.arange(len(row_labels)))
        ax.set_yticklabels(row_labels)
        ax.set_xlim(-0.5, len(clade_labels) - 0.5)
        ax.set_ylim(n_rows - 0.5, -0.5)
        ax.set_xticks(np.arange(-0.5, len(clade_labels), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.0)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.tick_params(axis="both", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if panel_idx < 2:
            ax.tick_params(labelbottom=False)

    fig.suptitle(title)
    if im is not None:
        cbar = fig.colorbar(im, ax=axes.ravel(), shrink=0.82, pad=0.02)
        ticks = [-log2_cap, -1.0, 0.0, 1.0, log2_cap]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f"{2 ** tick:g}" for tick in ticks])
        cbar.set_label("Odds ratio vs reference")
        cbar.ax.axhline(0, color="#111827", linewidth=1.2)

    add_panel_labels(axes.ravel())
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    save_figure(
        fig,
        FIG_DIR / output_stem,
        width="double",
        height_in=height_in,
        save_png=True,
        save_pdf=True,
    )
    plt.close(fig)

def _composition_or_summary_table(
    panels: list[tuple[str, pd.DataFrame]],
    *,
    predictor: str,
    clade_summary: pd.DataFrame,
    output_name: str,
) -> None:
    pieces = []
    for panel_name, values in panels:
        long = (
            values.rename_axis("level")
            .reset_index()
            .melt(id_vars="level", var_name="clade_group", value_name="odds_ratio")
        )
        long["model"] = panel_name
        long["predictor"] = predictor
        pieces.append(long)
    out = pd.concat(pieces, ignore_index=True)
    out = out.merge(clade_summary, on="clade_group", how="left")
    out.to_csv(FIG_DIR / output_name, index=False)


def main() -> None:
    clade_summary = _clade_summary()
    clade_order = clade_summary["clade_group"].tolist()
    labels = _clade_labels(clade_summary, clade_order)

    composition = _read_wald(RESULTS_DIR / "composition_wald.csv")
    comp_primary_single = _pivot_q(
        composition,
        domain="composition",
        model_set="primary",
        predictor_set="single",
        rows=COMPOSITION_ROWS,
        clade_order=clade_order,
    )
    comp_primary_joint = _pivot_q(
        composition,
        domain="composition",
        model_set="primary",
        predictor_set="joint",
        rows=COMPOSITION_ROWS,
        clade_order=clade_order,
    )
    comp_expanded_single = _pivot_q(
        composition,
        domain="composition",
        model_set="expanded",
        predictor_set="single",
        rows=COMPOSITION_ROWS,
        clade_order=clade_order,
    )
    comp_expanded_joint = _pivot_q(
        composition,
        domain="composition",
        model_set="expanded",
        predictor_set="joint",
        rows=COMPOSITION_ROWS,
        clade_order=clade_order,
    )
    comp_panels = [
        ("primary_single", comp_primary_single),
        ("primary_joint", comp_primary_joint),
        ("expanded_single", comp_expanded_single),
        ("expanded_joint", comp_expanded_joint),
    ]
    _draw_heatmap_grid(
        [
            ("Primary single-predictor", comp_primary_single),
            ("Primary joint", comp_primary_joint),
            ("Expanded single-predictor", comp_expanded_single),
            ("Expanded joint", comp_expanded_joint),
        ],
        row_labels=[label for _, label in COMPOSITION_ROWS],
        clade_labels=labels,
        cmap="YlGnBu",
        output_stem="fig15_clade_composition_significance",
    )
    _summary_table(
        comp_panels,
        rows=COMPOSITION_ROWS,
        clade_summary=clade_summary,
        output_name="fig15_clade_composition_significance_table.csv",
    )

    mixing = _read_wald(RESULTS_DIR / "mixing_wald.csv")
    mix_primary_single = _pivot_q(
        mixing,
        domain="node_mixing",
        model_set="primary",
        predictor_set="single",
        rows=MIXING_ROWS,
        clade_order=clade_order,
    )
    mix_primary_joint = _pivot_q(
        mixing,
        domain="node_mixing",
        model_set="primary",
        predictor_set="joint",
        rows=MIXING_ROWS,
        clade_order=clade_order,
    )
    mix_expanded_single = _pivot_q(
        mixing,
        domain="node_mixing",
        model_set="expanded",
        predictor_set="single",
        rows=MIXING_ROWS,
        clade_order=clade_order,
    )
    mix_expanded_joint = _pivot_q(
        mixing,
        domain="node_mixing",
        model_set="expanded",
        predictor_set="joint",
        rows=MIXING_ROWS,
        clade_order=clade_order,
    )
    mix_panels = [
        ("primary_single", mix_primary_single),
        ("primary_joint", mix_primary_joint),
        ("expanded_single", mix_expanded_single),
        ("expanded_joint", mix_expanded_joint),
    ]
    _draw_heatmap_grid(
        [
            ("Primary single-predictor", mix_primary_single),
            ("Primary joint", mix_primary_joint),
            ("Expanded single-predictor", mix_expanded_single),
            ("Expanded joint", mix_expanded_joint),
        ],
        row_labels=[label for _, label in MIXING_ROWS],
        clade_labels=labels,
        cmap="YlOrRd",
        output_stem="fig16_clade_mixing_significance",
    )
    _summary_table(
        mix_panels,
        rows=MIXING_ROWS,
        clade_summary=clade_summary,
        output_name="fig16_clade_mixing_significance_table.csv",
    )

    odds = _read_odds(RESULTS_DIR / "composition_odds_ratios.csv")
    effect_labels = clade_order
    effect_outputs = []
    for predictor, title, output_stem, levels in COMPOSITION_EFFECT_SPECS:
        effect_panels = []
        reference = ""
        for panel_key, panel_title, model_set, predictor_set in MODEL_PANELS:
            matrix, panel_reference = _composition_or_matrix(
                odds,
                predictor=predictor,
                model_set=model_set,
                predictor_set=predictor_set,
                levels=levels,
                clade_order=clade_order,
            )
            reference = reference or panel_reference
            effect_panels.append((panel_title, matrix))

        row_labels = [
            f"{level} (ref)" if reference and level == reference else level
            for level in levels
        ]
        _draw_or_heatmap_grid(
            effect_panels,
            row_labels=row_labels,
            clade_labels=effect_labels,
            title=title,
            output_stem=output_stem,
        )
        _composition_or_summary_table(
            [(key, matrix) for (key, _, _, _), (_, matrix) in zip(MODEL_PANELS, effect_panels)],
            predictor=predictor,
            clade_summary=clade_summary,
            output_name=f"{output_stem}_table.csv",
        )
        effect_outputs.append(FIG_DIR / f"{output_stem}.png")

    print(FIG_DIR / "fig15_clade_composition_significance.png")
    print(FIG_DIR / "fig16_clade_mixing_significance.png")
    for path in effect_outputs:
        print(path)


if __name__ == "__main__":
    main()
