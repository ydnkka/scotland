#!/usr/bin/env python3
"""
Generate per-(window, lineage) group .ids files and a tn93 command file.

For each 3-week window (1-week stride) and each pango_lineage with >= min_group_size
sequences, emits one compound shell command:

    samtools faidx -o <group>.fasta <all_seqs>.fasta.gz -r <group>.ids && \\
    tn93 -q -t <thr> -a <ambig> -o <group>.csv <group>.fasta

Usage:
  python3 clustering_scripts/02_gen_tn93_commands.py
  python3 clustering_scripts/02_gen_tn93_commands.py --config config.yaml --root /path/to/repo
"""

from __future__ import annotations

import argparse
import logging
import re
import shlex
import shutil
import sys
from datetime import timedelta
from pathlib import Path

import pandas as pd
import yaml


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def slugify(text: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "-", text.strip())
    return re.sub(r"-{2,}", "-", s).strip("-")


def build_windows(
    min_date: pd.Timestamp, max_date: pd.Timestamp,
    window_size: timedelta, step: timedelta,
) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    if max_date < min_date + window_size:
        return []
    starts = pd.date_range(start=min_date, end=max_date - window_size, freq=step)
    return [(s, s + window_size) for s in starts]


def load_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)

    if "sequence_id" in df.columns:
        df = df.drop_duplicates(subset=["sequence_id"], keep="first").set_index("sequence_id", drop=True)

    df["collection_date"] = pd.to_datetime(df["collection_date"], errors="coerce")
    n_before = len(df)
    df = df.dropna(subset=["collection_date", "pango_lineage"])
    if len(df) < n_before:
        logging.warning("Dropped %d rows with missing collection_date or pango_lineage.", n_before - len(df))
    return df


def write_ids_and_command(
    window_id: str, lineage: str, seq_ids: list[str],
    group_fasta_dir: Path, tn93_results_dir: Path,
    sequences_file: Path, tn93_threshold: float, tn93_ambiguous: str, tn93_quiet: bool,
) -> str:
    slug = slugify(lineage)
    name = f"{window_id}_{slug}_{len(seq_ids)}"

    ids_path = group_fasta_dir / f"{name}.ids"
    fasta_path = group_fasta_dir / f"{name}.fasta"
    csv_path = tn93_results_dir / f"{name}.csv"

    ids_path.write_text("\n".join(seq_ids) + "\n", encoding="utf-8")

    faidx_cmd = shlex.join(["samtools", "faidx", "-o", str(fasta_path), str(sequences_file), "-r", str(ids_path)])
    tn93_args = ["tn93"]
    if tn93_quiet:
        tn93_args.append("-q")
    tn93_args += ["-t", str(tn93_threshold), "-a", tn93_ambiguous, "-o", str(csv_path), str(fasta_path)]
    tn93_cmd = shlex.join(tn93_args)

    return f"{faidx_cmd} && {tn93_cmd}"


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate tn93 command file for windowed lineage groups.")
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()

    setup_logging(args.log_level)

    for exe in ("samtools", "tn93"):
        if shutil.which(exe) is None:
            logging.warning("%s not found on PATH — generated commands may fail.", exe)

    cfg = load_config(args.root / args.config)
    pipe = cfg["pipeline"]
    raw = {k: args.root / v for k, v in cfg["data"]["raw"].items()}
    proc = {k: args.root / v for k, v in cfg["data"]["processed"].items()}

    metadata = load_metadata(proc["metadata"])

    window_size = timedelta(weeks=pipe["window_size_weeks"])
    step = timedelta(weeks=pipe["step_weeks"])
    windows = build_windows(metadata["collection_date"].min(), metadata["collection_date"].max(), window_size, step)
    logging.info("Generated %d windows.", len(windows))

    group_fasta_dir: Path = proc["group_fasta_dir"]
    tn93_results_dir: Path = proc["tn93_results_dir"]
    group_fasta_dir.mkdir(parents=True, exist_ok=True)
    tn93_results_dir.mkdir(parents=True, exist_ok=True)

    commands: list[str] = [
        f"mkdir -p {shlex.quote(str(group_fasta_dir))} {shlex.quote(str(tn93_results_dir))}"
    ]
    n_groups = 0

    for i, (start, end) in enumerate(windows, start=1):
        wdf = metadata[(metadata["collection_date"] >= start) & (metadata["collection_date"] < end)]
        if wdf.empty:
            continue
        window_id = f"W{str(i).zfill(3)}"
        for lin, gdf in wdf.groupby("pango_lineage", sort=False):
            seq_ids = list(dict.fromkeys(gdf.index.tolist()))
            if len(seq_ids) < pipe["min_group_size"]:
                continue
            cmd = write_ids_and_command(
                window_id=window_id,
                lineage=str(lin),
                seq_ids=seq_ids,
                group_fasta_dir=group_fasta_dir,
                tn93_results_dir=tn93_results_dir,
                sequences_file=raw["sequences"],
                tn93_threshold=cfg["tn93"]["threshold"],
                tn93_ambiguous=cfg["tn93"]["ambiguous"],
                tn93_quiet=cfg["tn93"]["quiet"],
            )
            commands.append(cmd)
            n_groups += 1

    cmd_file: Path = proc["tn93_commands_file"]
    cmd_file.parent.mkdir(parents=True, exist_ok=True)
    cmd_file.write_text("\n".join(commands) + "\n", encoding="utf-8")

    logging.info("Wrote %d group commands to %s", n_groups, cmd_file)
    if n_groups > 1:
        logging.info("Example:\n  %s", commands[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
