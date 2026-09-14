"""Fetch only the three public files in the fixed initial TRAILS read audit."""

import argparse
import concurrent.futures
import hashlib
import json
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "config/ubash3a-raw-read-plan.json"


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(source):
    path = ROOT / source["local_path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    receipt = path.with_name(path.name + ".receipt.json")
    if path.exists():
        if not receipt.exists():
            raise RuntimeError(f"Existing file without receipt: {path}")
        previous = json.loads(receipt.read_text())
        assert path.stat().st_size == source["expected_bytes"]
        assert sha(path) == previous["sha256"]
        assert previous["url"] == source["url"]
        return previous
    partial = path.with_name(path.name + ".part")
    offset = partial.stat().st_size if partial.exists() else 0
    headers = {"Accept-Encoding": "identity"}
    if offset:
        headers.update(Range=f"bytes={offset}-", **{"If-Range": source["etag"]})
    started = datetime.now(timezone.utc).isoformat()
    with requests.get(source["url"], headers=headers, stream=True, timeout=(30, 120)) as response:
        response.raise_for_status()
        assert response.headers.get("ETag") == source["etag"], "Source ETag changed"
        if offset and response.status_code == 206:
            assert response.headers["Content-Range"].startswith(f"bytes {offset}-")
            mode = "ab"
        else:
            assert response.status_code == 200
            mode, offset = "wb", 0
        last_report, beginning = time.monotonic(), time.monotonic()
        total = offset
        with partial.open(mode) as stream:
            for chunk in response.iter_content(chunk_size=4 * 1024 * 1024):
                stream.write(chunk)
                total += len(chunk)
                if total > source["expected_bytes"]:
                    raise RuntimeError("Download exceeded fixed file size")
                if time.monotonic() - last_report >= 25:
                    print(json.dumps({"source": source["id"], "downloaded_bytes": total,
                                      "expected_bytes": source["expected_bytes"],
                                      "MB_per_second": round((total-offset)/1e6/(time.monotonic()-beginning), 2)}), flush=True)
                    last_report = time.monotonic()
        assert total == source["expected_bytes"], f"Incomplete download: {total}"
        result = {**source, "started_at_utc": started,
                  "completed_at_utc": datetime.now(timezone.utc).isoformat(),
                  "sha256": sha(partial), "size_bytes": total,
                  "integrity_status": "byte_count_and_transport_metadata_verified_decompression_pending",
                  "last_modified": response.headers.get("Last-Modified")}
    partial.replace(path)
    receipt.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"source": source["id"], "status": "download_complete", "size_bytes": total}), flush=True)
    return result


def fetch_reference_ranges(source):
    """Use eight bounded HTTP ranges of the identical pinned public reference."""
    path = ROOT / source["local_path"]
    if path.exists():
        return fetch(source)
    prefix = path.with_name(path.name + ".part")
    offset = prefix.stat().st_size if prefix.exists() else 0
    folder = path.with_name(path.name + ".ranges")
    folder.mkdir(parents=True, exist_ok=True)
    range_plan = folder / "plan.json"
    if range_plan.exists():
        plan = json.loads(range_plan.read_text())
        assert plan["prefix_bytes"] == offset
        assert plan["etag"] == source["etag"]
    else:
        step = (source["expected_bytes"] - offset + 7) // 8
        plan = {"prefix_bytes": offset, "etag": source["etag"],
                "ranges": [[a, min(a + step, source["expected_bytes"]) - 1]
                           for a in range(offset, source["expected_bytes"], step)]}
        range_plan.write_text(json.dumps(plan, indent=2) + "\n")
    started = datetime.now(timezone.utc).isoformat()

    def segment(item):
        i, (left, right) = item
        part = folder / f"segment-{i:02d}"
        previous = part.stat().st_size if part.exists() else 0
        expected = right - left + 1
        assert previous <= expected
        if previous < expected:
            headers = {"Accept-Encoding": "identity", "Range": f"bytes={left+previous}-{right}",
                       "If-Range": source["etag"]}
            with requests.get(source["url"], headers=headers, stream=True, timeout=(30, 120)) as response:
                response.raise_for_status()
                assert response.status_code == 206
                assert response.headers["Content-Range"] == f"bytes {left+previous}-{right}/{source['expected_bytes']}"
                assert response.headers.get("ETag") == source["etag"]
                with part.open("ab") as output:
                    for chunk in response.iter_content(1024 * 1024):
                        output.write(chunk)
                        if output.tell() > expected:
                            raise RuntimeError("Range exceeded planned size")
        assert part.stat().st_size == expected
        print(json.dumps({"source": source["id"], "completed_range": i, "bytes": expected}), flush=True)
        return part

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        parts = list(pool.map(segment, enumerate(plan["ranges"])))
    assembled = path.with_name(path.name + ".assembled.part")
    with assembled.open("wb") as output:
        for part in ([prefix] if offset else []) + parts:
            with part.open("rb") as stream:
                shutil.copyfileobj(stream, output, 8 * 1024 * 1024)
    assert assembled.stat().st_size == source["expected_bytes"]
    result = {**source, "started_at_utc": started, "completed_at_utc": datetime.now(timezone.utc).isoformat(),
              "sha256": sha(assembled), "size_bytes": assembled.stat().st_size,
              "download_method": "existing_prefix_then_eight_HTTP_ranges_with_exact_Content_Range_and_ETag_checks",
              "integrity_status": "byte_count_and_transport_metadata_verified_decompression_pending"}
    assembled.replace(path)
    path.with_name(path.name + ".receipt.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"source": source["id"], "status": "download_complete", "size_bytes": result["size_bytes"]}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-ranges", action="store_true")
    args = parser.parse_args()
    plan = json.loads(PLAN.read_text())
    assert plan["frozen_before_raw_target_inspection"]
    def dispatch(source):
        if args.reference_ranges and source["id"] == "GRCh38_primary":
            return fetch_reference_ranges(source)
        return fetch(source)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        receipts = list(pool.map(dispatch, plan["downloads"]))
    manifest = {"plan_sha256": sha(PLAN), "fetch_script_sha256": sha(Path(__file__)),
                "sources": receipts, "target_reads_analyzed": False}
    (ROOT / "config/ubash3a-raw-read-downloads.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
