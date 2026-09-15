#!/usr/bin/env python3
"""Retrieve exactly the frozen complete experiment, with transfer receipts."""
import concurrent.futures
import datetime
import hashlib
import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
SELECTION = ROOT / "config/ubash3a-proteomics-spectrum-selection.json"
RECEIPTS = ROOT / "work/ubash3a-proteomics-spectra/download-receipts"


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def fetch(record):
    path = ROOT / record["local_path"]
    receipt_path = RECEIPTS / (record["file_name"] + ".json")
    if path.exists() and receipt_path.exists():
        prior = json.loads(receipt_path.read_text())
        if prior.get("complete") and path.stat().st_size == record["expected_bytes"] and digest(path) == prior["sha256"]:
            return prior
        raise ValueError(f"Existing raw/receipt mismatch: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    receipt = {**record, "started_at_utc": now(), "complete": False,
               "publisher_checksum": None, "publisher_checksum_not_supplied": True}
    total, checksum = 0, hashlib.sha256()
    try:
        with requests.get(record["download_url"], stream=True, timeout=(30, 120)) as response:
            response.raise_for_status()
            assert response.status_code == 200
            receipt.update(status=response.status_code, response_url=response.url, headers=dict(response.headers))
            assert int(response.headers.get("Content-Length", record["expected_bytes"])) == record["expected_bytes"]
            with partial.open("wb") as handle:
                for block in response.iter_content(1024 * 1024):
                    total += len(block)
                    if total > record["expected_bytes"]:
                        raise ValueError("Transfer exceeded source byte count")
                    handle.write(block)
                    checksum.update(block)
        assert total == record["expected_bytes"], (record["file_name"], total)
        partial.replace(path)
        receipt.update(complete=True, downloaded_bytes=total, sha256=checksum.hexdigest(), completed_at_utc=now())
    except Exception as exc:
        receipt.update(error=str(exc), downloaded_bytes=total, ended_at_utc=now())
        receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
        raise
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"file": record["file_name"], "bytes": total, "complete": True}), flush=True)
    return receipt


if __name__ == "__main__":
    selection = json.loads(SELECTION.read_text())
    # Reference qualification is required before this larger transfer stage.
    specificity = json.loads((ROOT / "data/derived/ubash3a-proteomics-specificity/specificity-summary.json").read_text())
    verification = json.loads((ROOT / "reports/ubash3a-proteomics-specificity-verification.json").read_text())
    assert specificity["reference_sets"] == 5 and specificity["all_parent_checks_pass"]
    assert verification["all_checks_pass"]
    assert len(selection["files"]) == 14 and selection["download_bytes"] == 2387832397
    assert sum(row["expected_bytes"] for row in selection["files"]) == selection["download_bytes"]
    RECEIPTS.mkdir(parents=True, exist_ok=True)
    stage = {"started_at_utc": now(), "selection_sha256": digest(SELECTION),
             "specificity_verification_sha256": digest(ROOT / "reports/ubash3a-proteomics-specificity-verification.json"),
             "before_public_spectrum_scores": True, "download_budget_bytes": selection["download_bytes"],
             "local_storage_budget_gb": selection["local_storage_budget_gb"]}
    (ROOT / "config/ubash3a-proteomics-spectrum-download-start.json").write_text(json.dumps(stage, indent=2) + "\n")
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        receipts = list(pool.map(fetch, selection["files"]))
    result = {"complete": True, "completed_at_utc": now(), "files": receipts,
              "total_bytes": sum(row["downloaded_bytes"] for row in receipts)}
    (ROOT / "config/ubash3a-proteomics-spectrum-downloads.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"complete": True, "files": len(receipts), "total_bytes": result["total_bytes"]}), flush=True)
