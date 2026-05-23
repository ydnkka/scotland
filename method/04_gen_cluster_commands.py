#!/usr/bin/env python3
"""
Generate the parallel command file for running 03_process_group.py on all tn93 groups.

Pairs each <stem>.csv in tn93_results_dir with its <stem>.ids in group_fasta_dir and
emits one python3 command per group into a command file for parallel_run.sh.

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

    tn93_results_dir: Path = proc["tn93_results_dir"]
    group_fasta_dir: Path = proc["group_fasta_dir"]
    cluster_long_dir: Path = proc["cluster_long_dir"]
    metadata: Path = proc["metadata"]
    script = args.root / "method" / "03_process_group.py"

    inc_re = re.compile(args.include) if args.include else None
    exc_re = re.compile(args.exclude) if args.exclude else None

    resolutions = ",".join(str(r) for r in pipe["leiden_resolutions"])

    lines: list[str] = [f"mkdir -p {shlex.quote(str(cluster_long_dir))}"]
    n, missing = 0, 0

    for csv in sorted(tn93_results_dir.glob("*.csv")):
        stem = csv.stem
        if inc_re and not inc_re.search(stem):
            continue
        if exc_re and exc_re.search(stem):
            continue
        ids_path = group_fasta_dir / f"{stem}.ids"
        if not ids_path.exists():
            print(f"# WARNING: missing .ids for {stem}", file=sys.stderr)
            missing += 1
            continue
        cmd = (
            f"python3 {shlex.quote(str(script))}"
            f" --tn93-csv {shlex.quote(str(csv))}"
            f" --seq-ids {shlex.quote(str(ids_path))}"
            f" --metadata {shlex.quote(str(metadata))}"
            f" --out-long-dir {shlex.quote(str(cluster_long_dir))}"
            f" --alignment-length {pipe['alignment_length']}"
            f" --resolutions {shlex.quote(resolutions)}"
            f" --seed {pipe['seed']}"
            f" --sparsification {pipe['sparsification']}"
        )
        lines.append(cmd)
        n += 1

    cmd_file: Path = proc["cluster_commands_file"]
    cmd_file.parent.mkdir(parents=True, exist_ok=True)
    cmd_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote {n} commands to {cmd_file}", file=sys.stderr)
    if missing:
        print(f"Skipped {missing} groups with missing .ids files", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
