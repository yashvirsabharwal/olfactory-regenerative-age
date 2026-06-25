"""Release artifact manifest helpers."""

from __future__ import annotations

import glob
import hashlib
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_RELEASE_ARTIFACTS: list[dict[str, object]] = [
    {"path": "README.md", "category": "git_tracked", "required_for_review": True, "notes": "Repository overview and rerun instructions."},
    {"path": "pyproject.toml", "category": "git_tracked", "required_for_review": True, "notes": "Base Python package dependencies."},
    {"path": "Makefile", "category": "git_tracked", "required_for_review": True, "notes": "Primary workflow entry points."},
    {"path": "configs/gateway.yaml", "category": "git_tracked", "required_for_review": True, "notes": "Dataset and output configuration."},
    {"path": "configs/models.yaml", "category": "git_tracked", "required_for_review": True, "notes": "ORA model configuration."},
    {"path": "configs/command_manifest.yaml", "category": "git_tracked", "required_for_review": True, "notes": "Configured command-to-output provenance."},
    {"path": "configs/gene_sets.yaml", "category": "git_tracked", "required_for_review": True, "notes": "Curated module definitions."},
    {"path": "docs/claim_ledger.md", "category": "git_tracked", "required_for_review": True, "notes": "Claim boundaries and audit checklist."},
    {"path": "docs/changelog.md", "category": "git_tracked", "required_for_review": True, "notes": "Milestone implementation log."},
    {"path": "docs/external_label_request_log.md", "category": "git_tracked", "required_for_review": True, "notes": "External label request draft and outcome log."},
    {"path": "docs/project_status_and_remaining_work.md", "category": "git_tracked", "required_for_review": True, "notes": "Consolidated research status and remaining-work tracker."},
    {"path": "manuscript/main.tex", "category": "git_tracked", "required_for_review": True, "notes": "Manuscript source."},
    {"path": "manuscript/references.bib", "category": "git_tracked", "required_for_review": True, "notes": "Bibliography source."},
    {"path": "data/raw/gateway.h5ad", "category": "source_data", "required_for_review": True, "notes": "Primary Gateway H5AD source; redistribute according to source terms."},
    {"path": "data/external/GSE184117_RAW.tar", "category": "external_archive", "required_for_review": True, "notes": "External raw 10x validation archive."},
    {"path": "data/external/GSE184117_series_matrix.txt.gz", "category": "external_archive", "required_for_review": True, "notes": "External GEO metadata."},
    {"path": "data/processed/gateway_hvg3003_4m.h5ad", "category": "large_artifact", "required_for_review": True, "notes": "Reduced full-4M scVI substrate."},
    {"path": "data/processed/gateway_scvi_full_4m_reduced.h5ad", "category": "large_artifact", "required_for_review": True, "notes": "Full-4M reduced scVI embedding artifact."},
    {"path": "results/models/gateway_scvi_full_4m_reduced/model.pt", "category": "large_artifact", "required_for_review": True, "notes": "Full-4M reduced scVI model weights."},
    {"path": "data/processed/gateway_scvi_stratified_250k_seed23.h5ad", "category": "large_artifact", "required_for_review": True, "notes": "250k seed23 scVI sensitivity artifact."},
    {"path": "results/models/gateway_scvi_stratified_250k_seed23/model.pt", "category": "large_artifact", "required_for_review": True, "notes": "250k seed23 scVI model weights."},
    {"path": "data/processed/gateway_scvi_lineage_basal_neural_100k.h5ad", "category": "large_artifact", "required_for_review": True, "notes": "Lineage-focused scVI sensitivity artifact."},
]

EXTERNAL_ARCHIVE_CHECKSUM_MAX_BYTES = 1024 * 1024 * 1024


def build_release_manifest(
    command_manifest: dict[str, Any],
    *,
    base_dir: str | Path = ".",
    checksum_max_bytes: int = 100 * 1024 * 1024,
    extra_artifacts: list[dict[str, object]] | None = None,
) -> pd.DataFrame:
    """Build a publication artifact manifest from configured outputs and extra artifacts."""

    root = Path(base_dir)
    extras = list(DEFAULT_RELEASE_ARTIFACTS)
    if extra_artifacts:
        extras.extend(extra_artifacts)
    overrides = {str(spec["path"]): spec for spec in extras}
    rows: list[dict[str, object]] = []
    seen: set[str] = set()

    for stage, spec in command_manifest.get("commands", {}).items():
        command = str(spec.get("command", ""))
        inputs = ",".join(str(item) for item in spec.get("inputs", []))
        deferred = bool(spec.get("deferred", False))
        for pattern in spec.get("outputs", []):
            matches = _resolve_matches(root, str(pattern))
            if not matches:
                override = overrides.get(str(pattern), {})
                rows.append(
                    _artifact_row(
                        root,
                        Path(str(pattern)),
                        category=str(override.get("category", "locally_generated")),
                        artifact_group=stage,
                        generating_command=command,
                        source_input=inputs,
                        required_for_review=bool(override.get("required_for_review", not deferred)),
                        deferred=deferred,
                        checksum_max_bytes=checksum_max_bytes,
                        archive_uri=str(override.get("archive_uri", "")),
                        notes=str(override.get("notes", "")),
                    )
                )
                seen.add(str(pattern))
                continue
            for path in matches:
                rel = str(path.relative_to(root) if path.is_relative_to(root) else path)
                override = overrides.get(rel, {})
                rows.append(
                    _artifact_row(
                        root,
                        path,
                        category=str(override.get("category", "locally_generated")),
                        artifact_group=stage,
                        generating_command=command,
                        source_input=inputs,
                        required_for_review=bool(override.get("required_for_review", not deferred)),
                        deferred=deferred,
                        checksum_max_bytes=checksum_max_bytes,
                        archive_uri=str(override.get("archive_uri", "")),
                        notes=str(override.get("notes", "")),
                    )
                )
                seen.add(rel)

    for spec in extras:
        path = str(spec["path"])
        if path in seen:
            continue
        rows.append(
            _artifact_row(
                root,
                Path(path),
                category=str(spec.get("category", "git_tracked")),
                artifact_group=str(spec.get("artifact_group", "release_extra")),
                generating_command=str(spec.get("generating_command", "tracked or externally sourced")),
                source_input=str(spec.get("source_input", "")),
                required_for_review=bool(spec.get("required_for_review", True)),
                deferred=bool(spec.get("deferred", False)),
                checksum_max_bytes=checksum_max_bytes,
                archive_uri=str(spec.get("archive_uri", "")),
                notes=str(spec.get("notes", "")),
            )
        )
        seen.add(path)

    manifest = pd.DataFrame(rows)
    if manifest.empty:
        return manifest
    category_order = {"git_tracked": 0, "source_data": 1, "external_archive": 2, "large_artifact": 3, "locally_generated": 4}
    manifest["_category_order"] = manifest["category"].map(category_order).fillna(9)
    manifest = manifest.sort_values(["_category_order", "category", "artifact_group", "path"]).drop(columns="_category_order")
    return manifest.reset_index(drop=True)


def render_release_manifest_markdown(manifest: pd.DataFrame) -> str:
    """Render a compact Markdown index for the release artifact manifest."""

    lines = [
        "# Release Artifact Manifest",
        "",
        "This file is generated from the command manifest plus publication-critical source, external, and large artifacts.",
        "",
    ]
    if manifest.empty:
        return "\n".join([*lines, "No artifacts found.", ""])

    counts = manifest.groupby(["category", "artifact_status"], dropna=False).size().reset_index(name="artifacts")
    lines.extend(["## Category Summary", "", "| Category | Status | Artifacts |", "| --- | --- | ---: |"])
    for _, row in counts.iterrows():
        lines.append(f"| {row['category']} | {row['artifact_status']} | {int(row['artifacts'])} |")

    required_missing = manifest[
        manifest["required_for_review"].astype(bool)
        & ~manifest["artifact_status"].astype(str).isin(["present", "archived"])
    ]
    lines.extend(["", "## Required Items Needing Attention", ""])
    if required_missing.empty:
        lines.append("All required-for-review artifacts are present or archived.")
    else:
        lines.extend(["| Category | Path | Status | Notes |", "| --- | --- | --- | --- |"])
        for _, row in required_missing.head(50).iterrows():
            lines.append(f"| {row['category']} | `{row['path']}` | {row['artifact_status']} | {row['notes']} |")
        if required_missing.shape[0] > 50:
            lines.append(f"| ... | ... | ... | {required_missing.shape[0] - 50} additional rows omitted. |")

    uri_missing = manifest[
        manifest["required_for_review"].astype(bool)
        & manifest["category"].astype(str).isin(["source_data", "external_archive", "large_artifact", "locally_generated"])
        & manifest["archive_uri"].fillna("").eq("")
    ]
    lines.extend(
        [
            "",
            "## Archive URI Status",
            "",
            f"- Required artifacts without stable archive URI: {uri_missing.shape[0]}",
            "- Fill `archive_uri` during M2.4 once a durable repository or reviewer-access location is assigned.",
            "",
        ]
    )
    return "\n".join(lines)


def _artifact_row(
    root: Path,
    path: Path,
    *,
    category: str,
    artifact_group: str,
    generating_command: str,
    source_input: str,
    required_for_review: bool,
    deferred: bool,
    checksum_max_bytes: int,
    archive_uri: str,
    notes: str,
) -> dict[str, object]:
    resolved = path if path.is_absolute() else root / path
    display_path = str(resolved.relative_to(root) if resolved.exists() and resolved.is_relative_to(root) else path)
    exists = resolved.exists()
    status = "deferred" if deferred else "present" if exists else "missing"
    if archive_uri:
        status = "archived"
    size_bytes = _path_size(resolved) if exists else 0
    checksum_limit = _checksum_max_bytes_for_category(category, checksum_max_bytes)
    checksum_status = _checksum_status(resolved, size_bytes, checksum_limit) if exists else "missing"
    return {
        "category": category,
        "artifact_group": artifact_group,
        "path": display_path,
        "artifact_status": status,
        "required_for_review": bool(required_for_review),
        "size_bytes": int(size_bytes),
        "sha256": _sha256(resolved) if checksum_status == "sha256" else "",
        "checksum_status": checksum_status,
        "generating_command": generating_command,
        "source_input": source_input,
        "archive_uri": archive_uri,
        "notes": notes,
    }


def _resolve_matches(root: Path, pattern: str) -> list[Path]:
    candidate = Path(pattern)
    search = str(candidate if candidate.is_absolute() else root / candidate)
    if any(token in pattern for token in ["*", "?", "["]):
        return [Path(match) for match in sorted(glob.glob(search))]
    path = candidate if candidate.is_absolute() else root / candidate
    return [path] if path.exists() else []


def _path_size(path: Path) -> int:
    if path.is_dir():
        return int(sum(item.stat().st_size for item in path.rglob("*") if item.is_file()))
    return int(path.stat().st_size)


def _checksum_status(path: Path, size_bytes: int, checksum_max_bytes: int) -> str:
    if path.is_dir():
        return "directory"
    return "sha256" if size_bytes <= int(checksum_max_bytes) else "skipped_large_file"


def _checksum_max_bytes_for_category(category: str, checksum_max_bytes: int) -> int:
    if category == "external_archive":
        return max(int(checksum_max_bytes), EXTERNAL_ARCHIVE_CHECKSUM_MAX_BYTES)
    return int(checksum_max_bytes)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
