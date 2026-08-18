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
