"""Reviewer-access archive planning helpers."""

from __future__ import annotations

import pandas as pd


DESTINATION_BY_CATEGORY = {
    "git_tracked": "GitHub release plus optional Zenodo code DOI",
    "source_data": "Source repository citation; do not redistribute unless license permits",
    "external_archive": "GEO/SRA accession citation; optionally mirror checksum only",
    "large_artifact": "Institutional storage, OSF, Figshare, or Zenodo large-file record",
    "locally_generated": "Supplementary artifact archive or Zenodo/OSF project component",
}


def build_archive_review_package(release_manifest: pd.DataFrame) -> pd.DataFrame:
    """Convert the release manifest into reviewer-access action rows."""

    required = release_manifest[release_manifest["required_for_review"].astype(bool)].copy()
    if required.empty:
        return pd.DataFrame()
    rows = []
    for _, row in required.iterrows():
        category = str(row.get("category", ""))
        status = str(row.get("artifact_status", ""))
        raw_archive_uri = row.get("archive_uri", "")
        archive_uri = "" if pd.isna(raw_archive_uri) else str(raw_archive_uri)
        checksum_status = str(row.get("checksum_status", "") or "")
        rows.append(
            {
                "path": row.get("path", ""),
                "category": category,
                "artifact_status": status,
                "size_bytes": int(row.get("size_bytes", 0) or 0),
                "checksum_status": checksum_status,
                "sha256": row.get("sha256", ""),
                "current_archive_uri": archive_uri,
                "proposed_destination": DESTINATION_BY_CATEGORY.get(category, "Manual review"),
                "redistribution_review": _redistribution_review(category),
                "reviewer_access_action": _reviewer_access_action(category, status, archive_uri, checksum_status),
                "blocking_issue": _blocking_issue(category, status, archive_uri, checksum_status),
                "notes": row.get("notes", ""),
            }
        )
    output = pd.DataFrame(rows)
    dedupe_columns = [
        "path",
        "category",
        "artifact_status",
        "checksum_status",
        "current_archive_uri",
        "reviewer_access_action",
        "blocking_issue",
    ]
    output = output.drop_duplicates(subset=dedupe_columns, keep="first")
    return output.sort_values(["blocking_issue", "category", "path"], ascending=[False, True, True]).reset_index(drop=True)


def render_archive_review_markdown(package: pd.DataFrame) -> str:
    """Render a compact Markdown summary of archive actions."""

    lines = [
        "# Archive Reviewer-Access Plan",
        "",
        "Generated from `results/reports/release_artifact_manifest.tsv`.",
        "",
    ]
    if package.empty:
        return "\n".join([*lines, "No required artifacts found.", ""])
    blockers = package[package["blocking_issue"].astype(str).ne("")]
    lines.extend(
        [
            "## Summary",
            "",
            f"- Required artifacts: {package.shape[0]}",
            f"- Blocking archive/access issues: {blockers.shape[0]}",
            "",
            "## Blocking Items",
            "",
        ]
    )
    if blockers.empty:
        lines.append("No blocking archive/access issues remain.")
    else:
        lines.extend(["| Category | Path | Issue | Action |", "| --- | --- | --- | --- |"])
        for _, row in blockers.iterrows():
            lines.append(
                "| "
                f"{row['category']} | `{row['path']}` | {row['blocking_issue']} | "
                f"{row['reviewer_access_action']} |"
            )

    lines.extend(["", "## Destination Counts", "", "| Destination | Artifacts |", "| --- | ---: |"])
    counts = package.groupby("proposed_destination", dropna=False).size().reset_index(name="artifacts")
    for _, row in counts.iterrows():
        lines.append(f"| {row['proposed_destination']} | {int(row['artifacts'])} |")
    lines.append("")
    return "\n".join(lines)


def _redistribution_review(category: str) -> str:
    if category == "source_data":
        return "cite_source_not_redistribute_by_default"
    if category == "external_archive":
        return "public_accession_citation_preferred"
    if category == "large_artifact":
        return "license_and_size_review_required"
    if category == "locally_generated":
        return "archive_if_small_or_needed_for_review"
    return "safe_to_include_if_repository_policy_allows"


def _reviewer_access_action(category: str, status: str, archive_uri: str, checksum_status: str) -> str:
    if archive_uri:
        return "Verify URI resolves and cite in data availability."
    if status == "deferred":
        return "Recover or stage artifact, compute checksum if feasible, then assign reviewer-access URI."
    if status == "missing":
        return "Regenerate or document as intentionally unavailable before submission."
    if category in {"source_data", "external_archive"}:
        return "Cite public source accession/DOI and record checksum/source URL."
    if category == "large_artifact":
        return "Upload to approved large-artifact repository or institutional reviewer folder."
    if checksum_status == "skipped_large_file":
        return "Compute checksum using high-byte limit or external checksum job before archive."
    return "Include in repository release or supplementary artifact archive."


def _blocking_issue(category: str, status: str, archive_uri: str, checksum_status: str) -> str:
    if archive_uri:
        return ""
    if status in {"missing", "deferred"}:
        return f"required artifact is {status}"
    if category == "large_artifact":
        return "stable archive URI missing"
    if category in {"source_data", "external_archive"} and checksum_status == "skipped_large_file":
        return "large source checksum skipped locally"
    return ""
