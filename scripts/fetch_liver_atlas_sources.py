"""Download and verify the pinned public liver-atlas source files.

The manifest is the source of truth for URLs, byte lengths and SHA-256
digests.  Downloads stream to a sibling ``.part`` file and are moved into
place only after the complete digest and length checks pass.  ``--offline``
performs cache verification only and never contacts the network.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "config/liver-atlas-sources.json"
CHUNK_SIZE = 1024 * 1024


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload.get("sources"), list) or not payload["sources"]:
        raise ValueError("Source manifest has no sources")
    directory = payload.get("source_directory")
    if not isinstance(directory, str) or not directory:
        raise ValueError("Source manifest has no source_directory")
    return payload


def verify_cached(path: Path, source: Mapping[str, Any]) -> dict[str, Any]:
    """Verify an existing final cache entry against pinned bytes and digest."""

    if not path.exists():
        raise FileNotFoundError(path)
    expected_bytes = int(source["bytes"])
    actual_bytes = path.stat().st_size
    if actual_bytes != expected_bytes:
        raise ValueError(f"Cached {path.name} has {actual_bytes} bytes; expected {expected_bytes}")
    expected_sha = str(source["sha256"])
    actual_sha = sha256_file(path)
    if actual_sha != expected_sha:
        raise ValueError(f"Cached {path.name} has SHA-256 {actual_sha}; expected {expected_sha}")
    return {"id": source["id"], "file": str(path), "status": "verified", "bytes": actual_bytes, "sha256": actual_sha}


def _download(source: Mapping[str, Any], destination: Path, timeout: float = 60.0) -> dict[str, Any]:
    expected_bytes = int(source["bytes"])
    max_bytes = int(source.get("max_bytes", expected_bytes))
    if expected_bytes < 0 or max_bytes < expected_bytes:
        raise ValueError(f"Invalid byte limits for {source['id']}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        return verify_cached(destination, source)

    part = destination.with_name(destination.name + ".part")
    request = Request(str(source["url"]), headers={
        "User-Agent": "PSC-public-data-research/1.0",
        "Accept-Encoding": "identity",
    })
    digest = hashlib.sha256()
    received = 0
    try:
        with urlopen(request, timeout=timeout) as response, part.open("wb") as stream:
            declared = response.headers.get("Content-Length")
            if declared is not None and int(declared) > max_bytes:
                raise ValueError(f"Server Content-Length exceeds manifest limit for {destination.name}")
            while True:
                chunk = response.read(CHUNK_SIZE)
                if not chunk:
                    break
                received += len(chunk)
                if received > max_bytes:
                    raise ValueError(f"Downloaded {destination.name} exceeds manifest max_bytes")
                stream.write(chunk)
                digest.update(chunk)
            stream.flush()
            os.fsync(stream.fileno())
    except (HTTPError, URLError, OSError, ValueError):
        # Keep the partial file for diagnosis/retry.  It is never treated as
        # a valid source because only the final atomic rename creates the
        # cache entry used by analysis.
        raise
    if received != expected_bytes:
        raise ValueError(f"Downloaded {destination.name} has {received} bytes; expected {expected_bytes}")
    actual_sha = digest.hexdigest()
    if actual_sha != str(source["sha256"]):
        raise ValueError(f"Downloaded {destination.name} has SHA-256 {actual_sha}; expected {source['sha256']}")
    os.replace(part, destination)
    return verify_cached(destination, source)


def fetch_sources(
    manifest_path: Path = DEFAULT_MANIFEST,
    *,
    root: Path = ROOT,
    source_ids: Sequence[str] | None = None,
    offline: bool = False,
    timeout: float = 60.0,
) -> dict[str, Any]:
    manifest = load_manifest(manifest_path)
    wanted = set(source_ids or [])
    sources = [source for source in manifest["sources"] if not wanted or str(source["id"]) in wanted]
    if wanted and len(sources) != len(wanted):
        missing = sorted(wanted - {str(source["id"]) for source in sources})
        raise ValueError("Unknown source ID(s): " + ", ".join(missing))
    directory = root / str(manifest["source_directory"])
    results = []
    for source in sources:
        destination = directory / str(source["file"])
        if offline:
            result = verify_cached(destination, source)
            result["mode"] = "offline"
        else:
            result = _download(source, destination, timeout)
            result["mode"] = "download_or_cache"
        results.append(result)
    return {"manifest": str(manifest_path), "source_directory": str(directory), "offline": offline,
            "sources": results, "count": len(results)}


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--offline", action="store_true", help="Verify local final files without network access")
    parser.add_argument("--source-id", action="append", default=[], help="Fetch only this manifest source ID; repeatable")
    parser.add_argument("--timeout", type=float, default=60.0)
    args = parser.parse_args(argv)
    result = fetch_sources(args.manifest, source_ids=args.source_id, offline=args.offline, timeout=args.timeout)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError, HTTPError, URLError, OSError) as exc:
        print(f"liver-atlas source fetch failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
