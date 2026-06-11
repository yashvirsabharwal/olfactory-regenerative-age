"""Configuration loading helpers.

The project declares PyYAML for normal use, but this module includes a tiny
fallback parser so the core tests can run in minimal local runtimes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML config file.

    PyYAML is preferred. If unavailable, a conservative parser handles the simple
    nested dict/list/scalar YAML used by this repository's default configs.
    """

    config_path = Path(path)
    text = config_path.read_text(encoding="utf-8")
    try:
        import yaml  # type: ignore

        data = yaml.safe_load(text)
        return data or {}
    except ModuleNotFoundError:
        return _load_simple_yaml(text)


def project_path(path: str | Path, base_dir: str | Path | None = None) -> Path:
    """Resolve a path relative to a base directory or the current working directory."""

    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return Path(base_dir or ".").resolve() / candidate


def _load_simple_yaml(text: str) -> dict[str, Any]:
    lines = []
    for raw_line in text.splitlines():
        stripped = raw_line.split("#", 1)[0].rstrip()
        if stripped:
            indent = len(stripped) - len(stripped.lstrip(" "))
            lines.append((indent, stripped.lstrip(" ")))
    root: dict[str, Any] = {}
    stack: list[tuple[int, Any]] = [(-1, root)]
    pending_key: list[tuple[int, dict[str, Any], str]] = []

    for idx, (indent, content) in enumerate(lines):
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]

        if content.startswith("- "):
            value_text = content[2:].strip()
            if not isinstance(parent, list):
                if not pending_key:
                    raise ValueError(f"List item without list parent near: {content}")
                _, pending_parent, key = pending_key.pop()
                new_list: list[Any] = []
                pending_parent[key] = new_list
                stack.append((indent - 1, new_list))
                parent = new_list
            if ":" in value_text and not _looks_quoted(value_text):
                key, value = value_text.split(":", 1)
                item: dict[str, Any] = {key.strip(): _parse_scalar(value.strip())}
                parent.append(item)
                stack.append((indent, item))
            else:
                parent.append(_parse_scalar(value_text))
            continue

        if ":" not in content:
            raise ValueError(f"Unsupported YAML line: {content}")

        key, value = content.split(":", 1)
        key = key.strip()
        value = value.strip()
        if not isinstance(parent, dict):
            raise ValueError(f"Mapping item inside non-mapping near: {content}")

        if value == "":
            next_is_list = idx + 1 < len(lines) and lines[idx + 1][0] > indent and lines[idx + 1][1].startswith("- ")
            child: Any = [] if next_is_list else {}
            parent[key] = child
            stack.append((indent, child))
            if next_is_list:
                pending_key.append((indent, parent, key))
        else:
            parent[key] = _parse_scalar(value)

    return root


def _parse_scalar(text: str) -> Any:
    if text == "":
        return None
    if text in {"null", "Null", "NULL", "~"}:
        return None
    if text in {"true", "True", "TRUE"}:
        return True
    if text in {"false", "False", "FALSE"}:
        return False
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        if not inner:
            return []
        return [_parse_scalar(part.strip()) for part in inner.split(",")]
    if _looks_quoted(text):
        return text[1:-1]
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


def _looks_quoted(text: str) -> bool:
    return len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}

