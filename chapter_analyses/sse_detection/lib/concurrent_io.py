"""Small helpers for process-safe file writes and lock files."""

from __future__ import annotations

import contextlib
import os
from pathlib import Path
import socket
import time
from typing import Any, Iterator
from uuid import uuid4

try:
    import fcntl
except ImportError:  # pragma: no cover - fcntl is available on macOS/Linux.
    fcntl = None  # type: ignore[assignment]


class LockAlreadyHeldError(RuntimeError):
    """Raised when a fail-fast lock file is already held."""


def atomic_write_csv(frame: Any, path: Path | str, **kwargs: Any) -> None:
    """Write a CSV via a temporary file, then atomically replace the target."""
    with atomic_write_path(path) as tmp_path:
        frame.to_csv(tmp_path, **kwargs)


def atomic_write_parquet(frame: Any, path: Path | str, **kwargs: Any) -> None:
    """Write a parquet file via a temporary file, then atomically replace it."""
    with atomic_write_path(path) as tmp_path:
        frame.to_parquet(tmp_path, **kwargs)


def atomic_write_netcdf(data: Any, path: Path | str, **kwargs: Any) -> None:
    """Write a NetCDF file via a temporary file, then atomically replace it."""
    with atomic_write_path(path) as tmp_path:
        data.to_netcdf(tmp_path, **kwargs)


@contextlib.contextmanager
def atomic_write_path(path: Path | str) -> Iterator[Path]:
    """Yield a same-directory temporary path and replace the target on success."""
    final_path = Path(path)
    final_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = _temporary_sibling(final_path)
    try:
        yield tmp_path
        os.replace(tmp_path, final_path)
    finally:
        with contextlib.suppress(FileNotFoundError):
            tmp_path.unlink()


@contextlib.contextmanager
def exclusive_file_lock(lock_path: Path | str) -> Iterator[Path]:
    """Block on an advisory lock file until exclusive access is available."""
    path = Path(lock_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fcntl is None:
        with exclusive_create_lock(path):
            yield path
        return

    with path.open("a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        lock_handle.seek(0)
        lock_handle.truncate()
        lock_handle.write(_lock_metadata())
        lock_handle.flush()
        try:
            yield path
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


@contextlib.contextmanager
def exclusive_create_lock(
    lock_path: Path | str,
    *,
    details: str | None = None,
) -> Iterator[Path]:
    """Acquire a fail-fast lock by atomically creating a lock file."""
    path = Path(lock_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    token = uuid4().hex
    metadata = _lock_metadata(token=token, details=details)
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
    try:
        fd = os.open(path, flags, 0o644)
    except FileExistsError as exc:
        existing = _read_existing_lock(path)
        message = (
            f"Lock file already exists: {path}. Another process may be running "
            "the same model. Remove this file only if you are sure it is stale."
        )
        if existing:
            message = f"{message}\n\nExisting lock details:\n{existing}"
        raise LockAlreadyHeldError(message) from exc

    with os.fdopen(fd, "w", encoding="utf-8") as lock_handle:
        lock_handle.write(metadata)
        lock_handle.flush()

    try:
        yield path
    finally:
        _unlink_owned_lock(path, token)


def _temporary_sibling(path: Path) -> Path:
    token = f"{os.getpid()}.{time.time_ns()}.{uuid4().hex[:8]}"
    return path.with_name(f".{path.name}.{token}.tmp{path.suffix}")


def _lock_metadata(*, token: str | None = None, details: str | None = None) -> str:
    lines = [
        f"pid={os.getpid()}",
        f"host={socket.gethostname()}",
        f"started_at={time.strftime('%Y-%m-%dT%H:%M:%S%z')}",
    ]
    if token is not None:
        lines.append(f"token={token}")
    if details:
        lines.append("")
        lines.append(details)
    return "\n".join(lines) + "\n"


def _read_existing_lock(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8").strip()
    except OSError:
        return ""


def _unlink_owned_lock(path: Path, token: str) -> None:
    try:
        contents = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return
    except OSError:
        return
    if f"token={token}" in contents:
        with contextlib.suppress(FileNotFoundError):
            path.unlink()
