"""Fetch bounded, explicitly listed public inputs without replacing source pins."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import requests

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "config/ubash3a-hypothesis-test-plan.json"
RAW = ROOT / "data/raw/ubash3a-hypothesis-test"
MANIFEST = ROOT / "config/ubash3a-hypothesis-test-sources.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(spec, existing, offline, entrez_client, biostudies_client):
    target = RAW / spec["file"]
    if target.parent != RAW or not spec["source_id"]:
        raise ValueError("Source path must be a direct child of the target cache")
    old = existing.get(spec["source_id"])
    if target.exists():
        if old and sha(target) == old["sha256"] and target.stat().st_size == old["bytes"]:
            return old
        raise ValueError("Existing file is unpinned or changed: " + str(target))
    if offline:
        raise FileNotFoundError(target)
    temporary = target.with_suffix(target.suffix + ".part")
    if temporary.exists():
        raise ValueError("An earlier incomplete response needs inspection: " + str(temporary))
    record = {**spec, "retrieval_started_at_utc": datetime.now(timezone.utc).isoformat()}
    if spec.get("service") in {"entrez", "biostudies"}:
        client = entrez_client if spec["service"] == "entrez" else biostudies_client
        if not client:
            raise ValueError("Supply the corresponding installed research-skill client")
        payload = {**spec["request"], "save_raw": True, "raw_output_path": str(temporary), "max_items": 10, "timeout_sec": 45}
        env = {k: v for k, v in os.environ.items() if k not in {"NCBI_API_KEY", "NCBI_EUTILS_API_KEY", "NCBI_EMAIL"}}
        result = subprocess.run([sys.executable, str(client)], input=json.dumps(payload), text=True, capture_output=True, env=env, timeout=60)
        parsed = json.loads(result.stdout)
        if not parsed.get("ok"):
            raise RuntimeError(json.dumps(parsed))
        record["representation"] = "NCBI research-skill output cache; JSON is parsed and pretty-printed by the client, XML/text retains decoded response text." if spec["service"] == "entrez" else "BioStudies research-skill client saved decoded UTF-8 response text without JSON reformatting."
        record["client_sha256"] = sha(client)
        record["client_result"] = {k: v for k, v in parsed.items() if k not in {"raw_output_path", "text_head"}}
        time.sleep(.4)
    else:
        maximum = spec.get("maximum_bytes", 10_000_000)
        with requests.get(spec["url"], stream=True, timeout=(15, 60), headers={"User-Agent": "PSC-public-research/UBASH3A-hypothesis-test"}) as response:
            response.raise_for_status()
            record.update(status_code=response.status_code, final_url=response.url, content_type=response.headers.get("Content-Type"), declared_content_length=response.headers.get("Content-Length"))
            if response.headers.get("Content-Length") and int(response.headers["Content-Length"]) > maximum:
                raise ValueError("Declared file size exceeds fixed download bound")
            total = 0
            with temporary.open("wb") as stream:
                for block in response.iter_content(1_048_576):
                    total += len(block)
                    if total > maximum:
                        raise ValueError("Downloaded body exceeded fixed size bound")
                    stream.write(block)
        record["representation"] = "HTTP response body saved by requests iter_content; transfer decoding may be applied, content is not reformatted."
    if not temporary.is_file() or not temporary.stat().st_size:
        raise ValueError("Empty response")
    record.update(local_cache_path=str(target.relative_to(ROOT)), bytes=temporary.stat().st_size, sha256=sha(temporary), retrieved_at_utc=datetime.now(timezone.utc).isoformat())
    if spec.get("md5") and hashlib.md5(temporary.read_bytes()).hexdigest() != spec["md5"]:
        raise ValueError("Publisher MD5 mismatch")
    if old and (old["sha256"], old["bytes"]) != (record["sha256"], record["bytes"]):
        raise ValueError("Re-downloaded source differs from historical pin; response left in .part for inspection")
    temporary.replace(target)
    return old or record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("requests_file", type=Path)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--workers", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument("--entrez-client", type=Path)
    parser.add_argument("--biostudies-client", type=Path)
    args = parser.parse_args()
    specs = json.loads(args.requests_file.read_text())["requests"]
    if args.workers != 1 and any(s.get("service") == "entrez" for s in specs):
        raise ValueError("Entrez requests must be serial and rate-limited")
    if len({s["source_id"] for s in specs}) != len(specs):
        raise ValueError("Duplicate source IDs in request batch")
    manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {"manifest_id": "ubash3a-hypothesis-test-sources-v1", "created_at_utc": datetime.now(timezone.utc).isoformat(), "plan_sha256": sha(PLAN), "sources": [], "failed_attempts": []}
    if manifest["plan_sha256"] != sha(PLAN):
        raise ValueError("Plan changed")
    RAW.mkdir(parents=True, exist_ok=True)
    existing = {s["source_id"]: s for s in manifest["sources"]}
    failures = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = [(s, pool.submit(fetch, s, existing, args.offline, args.entrez_client, args.biostudies_client)) for s in specs]
        for spec, job in jobs:
            try:
                row = job.result()
                existing[spec["source_id"]] = row
                print(json.dumps({"source_id": spec["source_id"], "status": "verified", "bytes": row["bytes"]}), flush=True)
            except Exception as exc:
                failures += 1
                failure = {"source_id": spec["source_id"], "request": spec, "checked_at_utc": datetime.now(timezone.utc).isoformat(), "error": str(exc)}
                manifest["failed_attempts"].append(failure)
                print(json.dumps(failure), flush=True)
            manifest["sources"] = sorted(existing.values(), key=lambda r: r["source_id"])
            if not args.offline:
                MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
