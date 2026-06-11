"""Download helpers for Gateway source data."""

from __future__ import annotations

from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from .utils import ensure_parent


def download_direct_url(url: str, output_path: str | Path, chunk_size: int = 1024 * 1024) -> Path:
    """Download a URL to a local path using only the Python standard library."""

    output = ensure_parent(output_path)
    request = Request(url, headers={"User-Agent": "olfactory-regenerative-age/0.1"})
    with urlopen(request) as response, output.open("wb") as handle:  # noqa: S310 - user-supplied research URL
        while True:
            chunk = response.read(chunk_size)
            if not chunk:
                break
            handle.write(chunk)
    return output


def download_cellxgene_dataset(dataset_id: str, output_path: str | Path) -> Path:
    """Download a CELLxGENE source H5AD by dataset ID.

    This uses the public `cellxgene-census` helper when installed. The exact
    collection-to-dataset resolution step changes across CELLxGENE surfaces, so
    callers should pass a concrete dataset ID or direct H5AD URL.
    """

    output = ensure_parent(output_path)
    try:
        import cellxgene_census  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "cellxgene-census is required for dataset-id downloads. "
            "Install optional dependencies with `pip install -e '.[full]'`, "
            "or pass --url with a direct H5AD download URL."
        ) from exc

    download = getattr(cellxgene_census, "download_source_h5ad", None)
    if download is None:
        raise RuntimeError("Installed cellxgene-census does not expose download_source_h5ad().")

    try:
        download(dataset_id, to_path=str(output))
    except TypeError:
        # Older/newer releases have varied keyword spellings. Keep the fallback
        # narrow and explicit rather than hiding unrelated download errors.
        download(dataset_id=dataset_id, to_path=str(output))
    return output


def infer_download_mode(url: str | None, dataset_id: str | None) -> str:
    if url:
        parsed = urlparse(url)
        if parsed.scheme not in {"http", "https"}:
            raise ValueError("--url must be an http(s) URL.")
        return "url"
    if dataset_id:
        return "dataset_id"
    return "missing"

