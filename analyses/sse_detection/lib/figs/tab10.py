"""Build composition within cluster-size bands table."""

from __future__ import annotations

import argparse

from ..characterisation import ATTRIBUTES, SIZE_BANDS
from .common import Paths, add_common_args, paths_from_args
from .size_profiles import CATEGORY_ORDERS, SHORT_LABELS, load_saved_summaries, prepare

TABLE_NAME = "tab_size_composition"


def build(paths: Paths, *, tables=None):
    if tables is None:
        tables = prepare(paths)
    comp = tables["composition_contrasts"]
    directory = paths.publication_table_dir
    """Lay out saved percentages by category, with three columns per size band."""
    path = directory / f"{TABLE_NAME}.tex"
    directory.mkdir(parents=True, exist_ok=True)
    band_labels = [
        band.replace("-", "--").replace("100+", r"$\geq100$") for band in SIZE_BANDS
    ]
    header = [
        "Category & "
        + " & ".join(r"\multicolumn{3}{c}{" + band + "}" for band in band_labels)
        + r" \\",
        "".join(r"\cmidrule(lr){" + f"{2 + 3*i}-{4 + 3*i}" + "}" for i in range(5)),
        " & " + " & ".join([r"CSC & BG & $\Delta$"] * 5) + r" \\",
        r"\midrule",
    ]
    lines = [
        r"\begin{landscape}\begingroup",
        r"\renewcommand{\thesistablesetup}{\small\setstretch{1}\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.0}}",
        r"\begin{longtable}{@{}p{0.22\linewidth}*{15}{r}@{}}",
        r"\caption[Composition within size bands]{Composition within size bands. CSC and BG give the percentage of sequence records in each category for CSC and background clusters; $\Delta$ is CSC minus background in percentage points. Denominators are the record counts in Table~\ref{tab:size_overview}; all five attributes are complete.}\label{tab:size_composition}\\",
        r"\toprule",
        *header,
        r"\endfirsthead",
        r"\multicolumn{16}{l}{\tablename\ \thetable\ continued}\\",
        r"\toprule",
        *header,
        r"\endhead",
        r"\midrule\multicolumn{16}{r}{\small\itshape Continued on next page}\\\endfoot",
        r"\bottomrule\endlastfoot",
    ]
    for attr_index, (key, _, label) in enumerate(ATTRIBUTES):
        b = comp.loc[
            comp.scenario.eq("primary")
            & comp.weighting.eq("records")
            & comp.attribute.eq(key)
            & comp.size_band.ne("All")
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
                values.extend(
                    [
                        f"{r.candidate_percent:.1f}",
                        f"{r.background_percent:.1f}",
                        f"{r.difference_pp:+.1f}",
                    ]
                )
            lines.append(" & ".join(values) + r" \\")
    lines.extend([r"\end{longtable}\endgroup\end{landscape}"])
    path.write_text("\n".join(lines) + "\n")
    return {"tex": path}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument(
        "--from-saved-tables",
        action="store_true",
        help="Render saved size aggregates without rebuilding or writing them.",
    )
    args = parser.parse_args()
    paths = paths_from_args(args)
    build(paths, tables=load_saved_summaries(paths) if args.from_saved_tables else None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
