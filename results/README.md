# Publication Results

Project-level publication builders live here.

```bash
python -m results.make_figures --list
python -m results.make_tables --list
python -m results.make_figures --skip-missing
python -m results.make_tables --skip-missing
```

Build a specific artifact by fully qualified name:

```bash
python -m results.make_figures genomic_networks:fig_compatibility_topology
python -m results.make_tables sse_detection:tab_bayesian_model_specifications
```

Final figures are written to `results/figures/`; final LaTeX table fragments
are written to `results/tables/`. CSV and parquet tables stay under each
analysis package's own `results/tables/` directory.
