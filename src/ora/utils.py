"""Small file and table helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def ensure_parent(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def write_json(data: Any, path: str | Path) -> Path:
    output = ensure_parent(path)
    output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return output


def require_columns(columns: list[str], required: list[str], context: str) -> None:
    missing = [col for col in required if col not in columns]
    if missing:
        joined = ", ".join(missing)
        raise KeyError(f"Missing required columns for {context}: {joined}")


def normalize_token(value: object) -> str:
    if value is None:
        return ""
    try:
        if value != value:
            return ""
    except Exception:
        pass
    text = str(value).strip().lower()
    for old, new in [
        ("_", " "),
        ("-", " "),
        ("\u2019", "'"),
        ("\u2018", "'"),
        ("\u201c", '"'),
        ("\u201d", '"'),
    ]:
        text = text.replace(old, new)
    return " ".join(text.split())

