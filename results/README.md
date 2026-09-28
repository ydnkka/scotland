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
python -m results.make_figures genomic_networks:fig_compatibility_topology
python -m results.make_tables sse_detection:tab_bayesian_model_specifications
```

Final PNG/PDF figures are written to `results/figures/`. Final LaTeX table
fragments are written to `results/tables/`. CSV and parquet tables stay under
each analysis package's own `results/tables/` directory.

Each analysis library uses `figNN.py` for one figure and `tabNN.py` for one table,
with a `build(...)` entry point and a `FIGURE_NAME` or `TABLE_NAME` constant.
The numbers are stable script identifiers, independent of manuscript numbering;
gaps can remain when figures are retired. Output stems and registry keys are
semantic `fig_...`/`tab_...` names without chapter or supplement prefixes.
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
python -m results.make_figures sse_detection:fig_cluster_composition_overview_demographic
python -m results.make_figures sse_detection:fig_size_cumulative
python -m results.make_tables sse_detection:tab_size_overview --table-dir tmp/publication-tables
python -m analyses.sse_detection.lib.figs.fig22 --from-saved-tables --figure-dir tmp/figures
python -m analyses.sse_detection.lib.figs.tab07 --from-saved-tables --publication-table-dir tmp/publication-tables
```

Surveillance builder names and output stems now also carry `fig_` or `tab_`.
Bayesian table names use `tab_bayesian_fixed_effects_focal`,
`tab_bayesian_fixed_effects_full`, and `tab_bayesian_random_effect_sds` without
main/appendix placement. Existing generated files are not renamed automatically.
Surveillance CSV/parquet outputs and other analysis data companions still go to
the analysis-local table directories; `--table-dir` on the central command
selects the publication LaTeX destination. Individual table modules accept
`--publication-table-dir` for that destination (surveillance CSV/parquet modules
instead use `--table-dir`).
