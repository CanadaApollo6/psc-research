#!/usr/bin/env python3
"""Run the frozen Comet search over every qualified file in the selected set."""
import argparse
import datetime
import json
from pathlib import Path
import subprocess
import time

from analyze_ubash3a_proteomics_specificity import sha256

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/ubash3a-proteomics-spectra"


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def main(progressive=False):
    conversion_path = ROOT / "reports/ubash3a-proteomics-conversion.json"
    if progressive:
        order_path = ROOT / "config/ubash3a-proteomics-execution-order.json"
        assert json.loads(order_path.read_text())["changed_files_or_search_settings"] is False
        first = json.loads((ROOT / "reports/ubash3a-proteomics-first-conversion.json").read_text())
        assert first["all_conversion_and_content_checks_pass"]
        selection = json.loads((ROOT / "config/ubash3a-proteomics-spectrum-selection.json").read_text())
        conversion = {"files": [{"mzml_path": "work/ubash3a-proteomics-spectra/converted/" + Path(row["file_name"]).stem + ".mzML",
                                  "raw_file": row["file_name"]} for row in selection["files"]]}
    else:
        conversion = json.loads(conversion_path.read_text())
        assert conversion["complete_experiment"] and conversion["all_conversion_and_content_checks_pass"]
    database_path = ROOT / "config/ubash3a-proteomics-search-database.json"
    database = json.loads(database_path.read_text())
    qualification = ROOT / "reports/ubash3a-proteomics-comet-qualification.json"
    assert json.loads(qualification.read_text())["all_checks_pass"]
    plan_path = ROOT / "config/ubash3a-proteomics-search-plan.json"
    plan = json.loads(plan_path.read_text())
    parameters = ROOT / plan["parameter_file"]
    assert sha256(parameters) == plan["parameter_sha256"] == database["parameter_sha256"]
    assert sha256(ROOT / database["database_path"]) == database["database_sha256"]
    assert len(conversion["files"]) == 14
    tool = json.loads((ROOT / "config/ubash3a-proteomics-tool-install.json").read_text())["sources"][1]
    binary = ROOT / tool["installed_binary"]
    assert sha256(binary) == tool["binary_sha256"]
    started = now()
    start_path = ROOT / "config/ubash3a-proteomics-search-start.json"
    if not start_path.exists():
        start_path.write_text(json.dumps({"started_at_utc": started,
            "database_record_sha256": sha256(database_path), "database_sha256": database["database_sha256"],
            "parameter_sha256": sha256(parameters), "conversion_record_sha256": sha256(conversion_path) if conversion_path.exists() else None,
            "execution_order": "each complete file is independently qualified and pinned before scoring; whole-scope interpretation waits for all fourteen" if progressive else "all fourteen conversions completed before scoring",
            "execution_order_amendment_sha256": sha256(order_path) if progressive else None,
            "synthetic_qualification_sha256": sha256(qualification), "comet_binary_sha256": sha256(binary),
            "files": [row["mzml_path"] for row in conversion["files"]], "file_count": 14,
            "total_converted_MS2_spectra": conversion.get("total_MS2_spectra"),
            "candidate_scores_not_inspected_before_start": True}, indent=2) + "\n")
    else:
        prior = json.loads(start_path.read_text())
        assert prior["database_sha256"] == database["database_sha256"] and prior["parameter_sha256"] == sha256(parameters)
    log_dir = WORK / "search-logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for index, file in enumerate(conversion["files"]):
        if progressive:
            conversion_receipt = WORK / "conversion-logs" / (file["raw_file"] + ".qualification.json")
            while not conversion_receipt.exists():
                failure_path = WORK / "conversion-pipeline-failure.json"
                if failure_path.exists():
                    raise RuntimeError(f"Conversion pipeline failed; inspect {failure_path}")
                time.sleep(5)
            file = json.loads(conversion_receipt.read_text())
            assert file["binary_and_scan_checks_pass"] and file["all_vendor_MS2_scan_identifiers_equal_mzml"]
        mzml = ROOT / file["mzml_path"]
        assert sha256(mzml) == file["mzml_sha256"]
        output_paths = [mzml.with_suffix(suffix) for suffix in (".txt", ".pep.xml", ".pin")]
        receipt_path = log_dir / (mzml.stem + ".json")
        if receipt_path.exists():
            previous = json.loads(receipt_path.read_text())
            if previous.get("complete"):
                assert all(sha256(ROOT / path) == value for path, value in previous["output_sha256"].items())
                results.append(previous)
                continue
            raise RuntimeError(f"Preserved failed search requires review: {receipt_path}")
        args = [str(binary), "-P" + str(parameters), str(mzml)]
        record = {"file_index": index + 1, "input_mzml": file["mzml_path"], "input_sha256": file["mzml_sha256"],
                  "started_at_utc": now(), "args": args, "complete": False}
        if progressive:
            record["conversion_receipt_sha256_before_search"] = sha256(conversion_receipt)
        stdout_path = log_dir / (mzml.stem + ".stdout.txt")
        stderr_path = log_dir / (mzml.stem + ".stderr.txt")
        print(f"Searching frozen file {index + 1}/14: {mzml.name}", flush=True)
        with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
            result = subprocess.run(args, cwd=ROOT, stdout=stdout, stderr=stderr)
        record.update(ended_at_utc=now(), exit_code=result.returncode,
                      log_sha256={str(path.relative_to(ROOT)): sha256(path) for path in (stdout_path, stderr_path)})
        if result.returncode or not all(path.exists() and path.stat().st_size for path in output_paths):
            receipt_path.write_text(json.dumps(record, indent=2) + "\n")
            raise RuntimeError(f"Search failed: {mzml.name}; inspect preserved logs")
        record.update(complete=True, output_sha256={str(path.relative_to(ROOT)): sha256(path) for path in output_paths})
        receipt_path.write_text(json.dumps(record, indent=2) + "\n")
        results.append(record)
        print(f"Completed frozen file {index + 1}/14: {mzml.name}", flush=True)
    if progressive:
        while not conversion_path.exists():
            if (WORK / "conversion-pipeline-failure.json").exists():
                raise RuntimeError("Final conversion aggregation failed")
            time.sleep(5)
        complete_conversion = json.loads(conversion_path.read_text())
        assert complete_conversion["complete_experiment"] and len(complete_conversion["files"]) == 14
    report = {"complete_experiment": True, "finished_at_utc": now(), "files": results,
              "file_count": len(results), "search_start_sha256": sha256(start_path),
              "complete_conversion_record_sha256": sha256(conversion_path),
              "interpretation": "Search execution complete; candidate assessment and error-control analysis are separate"}
    (ROOT / "reports/ubash3a-proteomics-search-execution.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"all_14_searches_complete": True}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--progressive", action="store_true")
    args = parser.parse_args()
    main(args.progressive)
