"""Rename existing publication exports to match numbered producer scripts."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from .lib.config import RESULTS_DIR
from .lib.registry import publication_assets, write_asset_manifest


def migration_plan(results_dir: Path, assets: list[dict]) -> list[tuple[Path, Path]]:
    """Validate all collisions before moving any existing export."""
    moves = []
    for asset in assets:
        directory = results_dir / asset["directory"]
        for suffix in asset["formats"]:
            old = directory / f"{asset['legacy_stem']}.{suffix}"
            new = directory / f"{asset['stem']}.{suffix}"
            if not old.is_file():
                continue
            if new.exists():
                if not new.is_file() or hashlib.sha256(old.read_bytes()).digest() != hashlib.sha256(new.read_bytes()).digest():
                    raise FileExistsError(f"Conflicting old/new exports: {old} and {new}")
            moves.append((old, new))
    return moves


def migrate(results_dir: Path, *, dry_run: bool = False) -> int:
    moves = migration_plan(results_dir, publication_assets())
    for old, new in moves:
        print(f"{'Would rename' if dry_run else 'Renamed'} {old.name} -> {new.name}")
        if not dry_run:
            if new.exists():
                old.unlink()  # The preflight established byte-identical duplicate exports.
            else:
                old.rename(new)
    if not dry_run:
        manifest = write_asset_manifest(results_dir / "asset_manifest.json")
        print(f"Recorded publication producer mapping: {manifest}")
    return len(moves)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    count = migrate(args.results_dir.expanduser().resolve(), dry_run=args.dry_run)
    print(f"{count} exports {'would be renamed' if args.dry_run else 'migrated'}")


if __name__ == "__main__":
    main()
