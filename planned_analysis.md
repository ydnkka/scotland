Similar to `method/03_process_group.py`, I want to have a script that process one (window, lineage) group: tn93 CSV → EpiLink compatibility weights → weighted Assortativity Analysis.

For this, main variables of interests are age group, sex, geographies (data zone, council areas, and health boards), and adjusted SIMD index (as in `sse_detection`).

Will need a script that process one group and one that produce a list of commands like `04_gen_cluster_commands.py` to be run in parallel with method/parallel_run.sh, and one that consolidate the results to a long table.
