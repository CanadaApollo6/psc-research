#!/usr/bin/env python3
"""Restore pinned public PSC blood inputs; never silently accept changed sources."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true", help="No network: require and hash all cached files")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "config/psc-blood-sources.json").read_text())
    cross = manifest["reused_gene_crosswalk"]
    if sha256(ROOT / cross["path"]) != cross["sha256"]:
        raise ValueError("The versioned complete Ensembl relationship export differs from its pin")
    for source in manifest["source_files"]:
        dest = (ROOT / source["path"]).resolve()
        if not dest.is_relative_to((ROOT / "data/raw").resolve()):
            raise ValueError("Downloads must remain under ignored data/raw/")
        if dest.exists():
            if dest.stat().st_size != source["bytes"] or sha256(dest) != source["sha256"]:
                raise ValueError(f"Existing input differs from its pin: {source['path']}")
            print("verified", source["path"], flush=True)
            continue
        if args.check_only:
            raise FileNotFoundError(source["path"])
        if not source["url"].startswith("https://"):
            raise ValueError("Expected a pinned public HTTPS URL")
        dest.parent.mkdir(parents=True, exist_ok=True)
        partial = dest.with_name(dest.name + ".partial")
        request = urllib.request.Request(source["url"], headers={"User-Agent":"psc-research-reproduction/1.0", "Accept-Encoding":"identity"})
        received = 0
        restored_at = datetime.now(timezone.utc).isoformat()
        try:
            with urllib.request.urlopen(request, timeout=90) as response, partial.open("wb") as output:
                while block := response.read(1024 * 1024):
                    received += len(block)
                    if received > source["bytes"]:
                        raise ValueError("Source exceeds its frozen byte count")
                    output.write(block)
            if received != source["bytes"] or sha256(partial) != source["sha256"]:
                raise ValueError(f"Upstream source changed: {source['url']}; retain original pins and investigate")
            partial.replace(dest)
        finally:
            partial.unlink(missing_ok=True)
        receipt = {"url": source["url"], "retrieved_at_utc": restored_at, "path": source["path"], "bytes": received, "sha256": sha256(dest), "status": "restored_exact_pinned_body"}
        receipt_path = ROOT / "work/research-path/restored-input-receipts.jsonl"
        receipt_path.parent.mkdir(parents=True, exist_ok=True)
        with receipt_path.open("a") as stream:
            stream.write(json.dumps(receipt, sort_keys=True) + "\n")
        print("restored", source["path"], flush=True)


if __name__ == "__main__":
    main()
