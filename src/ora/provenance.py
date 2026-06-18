"""Command and output provenance helpers."""

from __future__ import annotations

import glob
import hashlib
from pathlib import Path
from typing import Any

import pandas as pd


def command_manifest_table(command_manifest: dict[str, Any]) -> pd.DataFrame:
    """Flatten configured pipeline commands into a table."""

    rows = []
    for stage, spec in command_manifest.get("commands", {}).items():
        rows.append(
            {
                "stage": stage,
                "description": spec.get("description", ""),
                "command": spec.get("command", ""),
                "inputs": ",".join(str(item) for item in spec.get("inputs", [])),
                "outputs": ",".join(str(item) for item in spec.get("outputs", [])),
                "deferred": bool(spec.get("deferred", False)),
            }
        )
    return pd.DataFrame(rows)


def output_provenance_table(
    command_manifest: dict[str, Any],
    *,
    base_dir: str | Path = ".",
    checksum_max_bytes: int = 100 * 1024 * 1024,
) -> pd.DataFrame:
    """Summarize existence, size, mtime, and optional checksum for manifest outputs."""

    root = Path(base_dir)
    rows = []
    for stage, spec in command_manifest.get("commands", {}).items():
        if spec.get("deferred", False):
            continue
        for pattern in spec.get("outputs", []):
            matches = _resolve_matches(root, str(pattern))
            if not matches:
                rows.append(_missing_row(stage, str(pattern)))
                continue
            for path in matches:
                rows.append(_existing_row(stage, root, path, checksum_max_bytes=checksum_max_bytes))
    return pd.DataFrame(rows)


def _resolve_matches(root: Path, pattern: str) -> list[Path]:
    candidate = Path(pattern)
    search = str(candidate if candidate.is_absolute() else root / candidate)
    if any(token in pattern for token in ["*", "?", "["]):
        return [Path(match) for match in sorted(glob.glob(search))]
    path = candidate if candidate.is_absolute() else root / candidate
    return [path] if path.exists() else []


def _missing_row(stage: str, path: str) -> dict[str, object]:
    return {
        "stage": stage,
        "path": path,
        "exists": False,
        "size_bytes": 0,
        "mtime_utc": "",
        "sha256": "",
        "checksum_status": "missing",
    }


def _existing_row(stage: str, root: Path, path: Path, *, checksum_max_bytes: int) -> dict[str, object]:
    stat = path.stat()
    size_bytes = _directory_size(path) if path.is_dir() else int(stat.st_size)
    checksum_status = _checksum_status(path, size_bytes, checksum_max_bytes)
    return {
        "stage": stage,
        "path": str(path.relative_to(root) if path.is_relative_to(root) else path),
        "exists": True,
        "size_bytes": size_bytes,
        "mtime_utc": pd.Timestamp(stat.st_mtime, unit="s", tz="UTC").isoformat(),
        "sha256": _sha256(path) if checksum_status == "sha256" else "",
        "checksum_status": checksum_status,
    }


def _checksum_status(path: Path, size_bytes: int, checksum_max_bytes: int) -> str:
    if path.is_dir():
        return "directory"
    return "sha256" if size_bytes <= int(checksum_max_bytes) else "skipped_large_file"


def _directory_size(path: Path) -> int:
    return int(sum(item.stat().st_size for item in path.rglob("*") if item.is_file()))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
