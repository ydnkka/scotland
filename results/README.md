# Publication Results

Project-level publication builders live here. The individual figure and table
implementations remain inside the analysis packages under `analyses/*/lib/figs/`;
this package provides the central registry and command-line entry points.

```bash
python -m results.make_figures --list
python -m results.make_tables --list
python -m results.make_figures --skip-missing
python -m results.make_tables --skip-missing
```

Restrict a run by analysis domain or build a specific artifact by fully
qualified name:

```bash
python -m results.make_figures --domain surveillance
python -m results.make_tables --domain genomic_networks
python -m results.make_figures genomic_networks:fig09_compatibility_topology
python -m results.make_tables sse_detection:tab01_bayesian_model_specifications
```

Final PNG/PDF figures are written to `results/figures/`. Final LaTeX table
fragments are written to `results/tables/`. CSV and parquet tables stay under
each analysis package's own `results/tables/` directory.

Each analysis library uses `figNN.py` for one figure and `tabNN.py` for one table,
with a `build(...)` entry point and a `FIGURE_NAME` or `TABLE_NAME` constant.
The numbers are stable script identifiers, independent of manuscript numbering;
gaps can remain when figures are retired. Output stems and registry keys start
with their producer's identifier: `fig05.py` in genomic networks writes
`fig05_cluster_landscape.pdf/.png`, while its `tab03.py` writes
`tab03_cluster_period_summary.tex`. Descriptive suffixes distinguish scripts
with the same number in different domains. Filenames do not encode manuscript
chapter or supplement placement. The registry validates this correspondence.
Shared data preparation, formatting, and panel drawing live in named helper
modules alongside the scripts.

| Library | Figure scripts | Table scripts | Shared helpers |
| --- | --- | --- | --- |
| Surveillance | `fig01`, `fig02` | `tab01`–`tab08` | `sequence_data`, `policy_data` |
| Genomic networks | `fig01`–`fig10` | `tab01`–`tab07` | `common`, `table_common`, `assortativity_analysis`, `assortativity_panels` |
| SSE detection | `fig01`–`fig11`, `fig13`–`fig15`, `fig20`–`fig29` | `tab01`–`tab10` | `common`, `forest`, `forest_common`, `composition_profiles`, `size_profiles`, `graph_styles`, `table_common` |

SSE cluster composition uses only the demographic and health-board overviews.
The individual sex, age, SIMD, urban/rural, and health-board figures are retired.
The former bundled size-profile and composition-profile builder names have
been replaced with individual names shown by `--list`.

```bash
python -m results.make_figures sse_detection:fig20_cluster_composition_overview_demographic
python -m results.make_figures sse_detection:fig22_size_cumulative
python -m results.make_tables sse_detection:tab07_size_overview --table-dir tmp/publication-tables
python -m analyses.sse_detection.lib.figs.fig22 --from-saved-tables --figure-dir tmp/figures
python -m analyses.sse_detection.lib.figs.tab07 --from-saved-tables --publication-table-dir tmp/publication-tables
```

Surveillance exports also carry their script number. Bayesian table names use
`tab03_bayesian_fixed_effects_focal`, `tab05_bayesian_fixed_effects_full`, and
`tab06_bayesian_random_effect_sds` without main/appendix placement.
Surveillance CSV/parquet outputs and other analysis data companions still go to
the analysis-local table directories; `--table-dir` on the central command
selects the publication LaTeX destination. Individual table modules accept
`--publication-table-dir` for that destination (surveillance CSV/parquet modules
instead use `--table-dir`).

## Editing table captions

All publication table captions live in [`table_captions.yaml`](table_captions.yaml),
keyed by their exact LaTeX label (for example, `tab:policy_denominators`). Each
entry contains:

- `caption`: concise text printed with the table.
- `short_caption`: plain-text title used in the list of tables.
- `columns`: descriptions of the displayed columns or grouped subcolumns.
- `notes`: optional definitions and interpretation details.
- `source`: the table script to rebuild after editing.
- `caption_is_latex`: optional flag for captions containing LaTeX math or references.

Column descriptions and notes are editorial metadata and are not automatically
appended to the caption. Keep essential interpretation in `caption`; use the
other fields for longer explanations. Column headings and layout are set in the
table scripts.

Keep each caption on one line using YAML single quotes, which preserve LaTeX
backslashes. Double any apostrophe inside the quoted text (`''`). Plain-text
captions are escaped automatically, including `%`, `&`, and `<`. With
`caption_is_latex: true`, write LaTeX directly, including `\%` for percentages
and `Table~\ref{tab:label}` for references. Short captions are always plain text.

Table scripts load the entry using the shared helper:

```python
from utils import load_table_metadata, table_caption, write_latex_table

write_latex_table(
    output_path,
    **table_caption("tab:policy_denominators"),
    columns=columns,
    rows=rows,
)
column_descriptions = load_table_metadata("tab:policy_denominators")["columns"]
```

Data-dependent captions use `{{name}}` placeholders. For example, the pairwise
distance table supplies `n_summaries` and `n_pairs` as formatted strings to
`table_caption(...)`. These placeholders leave ordinary LaTeX braces intact;
missing labels or required values raise an error. The YAML path is resolved
relative to the helper, so loading does not depend on the working directory.

After editing, regenerate a table with its usual command, for example:

```bash
python -m results.make_tables genomic_networks:tab02_policy_denominators
```

## Migrating existing exports and syncing thesis assets

```bash
python -m results.migrate_asset_names --dry-run
python -m results.migrate_asset_names
```

Migration moves the existing PDF/PNG and LaTeX exports to their numbered names
without regenerating scientific results. It validates all old/new name collisions
before making changes; byte-identical duplicates can be consolidated, while
different existing exports are preserved with an explicit conflict error.
Re-running migration is safe when no old names remain.

`results/asset_manifest.json` records every publication stem, source module,
script identifier, format, and legacy name. Central figure/table builds refresh
it after successful builds. Surveillance CSV/parquet builders are analysis-local
and are excluded from this publication asset inventory.

From the thesis directory, synchronize and remove the old managed export names:

```bash
python assets/scotland/sync_scotland_assets.py --clean
```

The thesis sync validates all manifest-declared source exports before copying,
records their checksums and producers, and leaves unrelated authored files alone.
