#!/usr/bin/env python3
"""
Generate the parallel command file for running cluster_pairwise.py.

Pairs each <stem>.parquet in pairwise_distances_dataset with the expected
<stem>.ids file in group_fasta_dir and emits one python3 command per group into
a command file for parallel_run.sh.

Usage:
    python3 method/04_gen_cluster_commands.py
    python3 method/04_gen_cluster_commands.py --config config.yaml --root /path/to/repo
"""

from __future__ import annotations

import argparse
import re
import shlex
import sys
from pathlib import Path

import yaml


def load_config(path: Path) -> dict:
    """Load the YAML pipeline configuration file."""
    with open(path) as f:
        return yaml.safe_load(f)


def main() -> int:
    """Generate the command file for per-group Leiden clustering jobs."""
    ap = argparse.ArgumentParser(
        description="Generate cluster processing command file."
    )
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument(
        "--include", type=str, default=None, help="Regex to filter group names"
    )
    ap.add_argument(
        "--exclude", type=str, default=None, help="Regex to exclude group names"
    )
    args = ap.parse_args()

    cfg = load_config(args.root / args.config)
    pipe = cfg["pipeline"]
    proc = {k: args.root / v for k, v in cfg["data"]["processed"].items()}

    pairwise_dir: Path = proc.get(
        "pairwise_distances_dataset",
        args.root / "data/processed/pairwise_distances_dataset",
    )
    group_fasta_dir: Path = proc["group_fasta_dir"]
    cluster_long_dir: Path = proc["cluster_long_dir"]
    script = args.root / "method" / "cluster_pairwise.py"
    if not pairwise_dir.exists():
        raise SystemExit(f"Pairwise dataset directory not found: {pairwise_dir}")

    inc_re = re.compile(args.include) if args.include else None
    exc_re = re.compile(args.exclude) if args.exclude else None

    resolutions = ",".join(str(r) for r in pipe["leiden_resolutions"])

    lines: list[str] = [f"mkdir -p {shlex.quote(str(cluster_long_dir))}"]
    n, missing_ids = 0, 0

    for pairwise_file in sorted(pairwise_dir.glob("*.parquet")):
        stem = pairwise_file.stem
        if inc_re and not inc_re.search(stem):
            continue
        if exc_re and exc_re.search(stem):
            continue
        ids_path = group_fasta_dir / f"{stem}.ids"
        cmd = (
            f"python3 {shlex.quote(str(script))}"
            f" --pairwise-file {shlex.quote(str(pairwise_file))}"
            f" --seq-ids {shlex.quote(str(ids_path))}"
            f" --out-long-dir {shlex.quote(str(cluster_long_dir))}"
            f" --resolutions {shlex.quote(resolutions)}"
            f" --seed {pipe['seed']}"
            f" --sparsification {pipe['sparsification']}"
        )
        if not ids_path.exists():
            missing_ids += 1
        lines.append(cmd)
        n += 1

    if n == 0:
        raise SystemExit("No pairwise parquet files matched the requested filters.")

    cmd_file: Path = proc["cluster_commands_file"]
    cmd_file.parent.mkdir(parents=True, exist_ok=True)
    cmd_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote {n} commands to {cmd_file}", file=sys.stderr)
    if missing_ids:
        print(
            f"{missing_ids} expected .ids files were not found while generating; "
            "commands still include --seq-ids and cluster_pairwise.py will infer "
            "nodes from pairwise endpoints only if the file is also missing at run time",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
