#!/usr/bin/env python3
"""Consolidate retrieval provenance and verify every completed local input."""
import json
from pathlib import Path

from analyze_ubash3a_proteomics_specificity import sha256

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/ubash3a-proteomics-specificity"
SAFE_HEADERS = {"content-type", "content-length", "content-encoding", "content-range", "etag", "last-modified", "date",
                "x-total-results", "x-uniprot-release", "x-uniprot-release-date", "link"}


def relative(path):
    path = Path(path)
    return str(path.relative_to(ROOT)) if path.is_absolute() else str(path)


def compact_receipt(path):
    data = json.loads(path.read_text())
    result = data.get("result", {})
    local_path = data.get("local_path") or result.get("raw_output_path")
    artifact = None
    if local_path:
        local_path = relative(local_path)
        source = ROOT / local_path
        if source.exists() and data.get("sha256"):
            assert sha256(source) == data["sha256"], source
            artifact = {"local_path": local_path, "sha256": data["sha256"], "size_bytes": source.stat().st_size}
    headers = {key: value for key, value in data.get("headers", {}).items() if key.lower() in SAFE_HEADERS}
    responses = []
    for response in data.get("response_headers_captured", []):
        responses.append({"url": response["url"], "status": response["status"],
                          "headers": {key: value for key, value in response.get("headers", {}).items() if key.lower() in SAFE_HEADERS}})
    request = data.get("request")
    if request:
        request = {key: value for key, value in request.items() if key in {
            "base_url", "path", "params", "response_format", "timeout_sec", "save_raw", "max_items"}}
    return {"receipt_id": path.name, "original_receipt_path": relative(path),
            "original_receipt_sha256": sha256(path), "source_url": data.get("url"),
            "response_url": data.get("response_url"), "request": request,
            "started_at_utc": data.get("started_at_utc"),
            "completed_or_retrieved_at_utc": data.get("completed_at_utc", data.get("finished_at_utc", data.get("ended_at_utc", data.get("retrieved_at_utc")))),
            "status": data.get("status", result.get("status_code")), "headers": headers,
            "service_response_headers": responses, "artifact": artifact,
            "error": data.get("error", result.get("error")),
            "snapshot_format": data.get("snapshot_format"),
            "client_path": data.get("client_path"), "client_sha256": data.get("client_sha256"),
            "publisher_ftp_mdtm": data.get("publisher_mdtm"), "publisher_bytes": data.get("publisher_size_bytes")}


def main():
    receipts = [compact_receipt(path) for path in sorted(WORK.glob("*.receipt.json"))]
    failures = []
    for path in sorted(WORK.glob("*.initial-network-failure.json")):
        failures.append({"local_path": str(path.relative_to(ROOT)), "sha256": sha256(path),
                         "record": json.loads(path.read_text())})
    lock_path = ROOT / "config/ubash3a-proteomics-reference-lock.json"
    lock = json.loads(lock_path.read_text())
    for source in lock["sources"]:
        assert sha256(ROOT / source["local_path"]) == source["sha256"]
    downloads_path = ROOT / "config/ubash3a-proteomics-spectrum-downloads.json"
    downloads = json.loads(downloads_path.read_text())
    assert downloads["complete"] and len(downloads["files"]) == 14
    for file in downloads["files"]:
        path = ROOT / file["local_path"]
        assert path.stat().st_size == file["expected_bytes"] == file["downloaded_bytes"]
        assert sha256(path) == file["sha256"]
    reused = json.loads((ROOT / "config/ubash3a-proteomics-specificity-plan.json").read_text())["input_pins"]
    for path, expected in reused.items():
        assert sha256(ROOT / path) == expected
    report = {"study": "PXD000376", "qualification_receipts": receipts, "initial_network_failures": failures,
              "uniprot_complete_pagination_record": "config/ubash3a-proteomics-reference-lock.json",
              "reference_lock_sha256": sha256(lock_path), "reference_sources": lock["sources"],
              "prior_inputs_preserved": reused, "spectrum_download_manifest": str(downloads_path.relative_to(ROOT)),
              "spectrum_download_manifest_sha256": sha256(downloads_path),
              "verified_raw_files": len(downloads["files"]), "verified_raw_bytes": sum(file["downloaded_bytes"] for file in downloads["files"]),
              "all_recorded_input_hash_checks_pass": True,
              "snapshot_limitations": ["Service JSON was parsed and reserialized by the unchanged skill clients",
                  "UniProt FASTA was decoded as text and concatenated from every returned page with release and total checks",
                  "The original all-at-once UniProt stream failed and was not substituted for a complete export",
                  "Deposited pepXML files were range-read for complete search settings; their old peptide identifications were not used as new protein evidence",
                  "No publisher raw-spectrum checksum was supplied; source byte totals, HTTP metadata and computed SHA-256 are preserved",
                  "Large input caches and vendor binaries remain in ignored local storage; Git contains their provenance and small outputs",
                  "Live API endpoints can change releases; exact offline replay depends on the pinned local snapshots"]}
    path = ROOT / "config/ubash3a-proteomics-sources.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"source_receipts": len(receipts), "verified_reference_records": lock["total_records"],
                      "verified_raw_files": report["verified_raw_files"], "verified_raw_bytes": report["verified_raw_bytes"]}))


if __name__ == "__main__":
    main()
