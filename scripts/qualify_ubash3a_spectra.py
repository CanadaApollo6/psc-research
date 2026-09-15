#!/usr/bin/env python3
"""Convert the frozen raw files and independently inspect complete mzML content."""
import argparse
import base64
from collections import Counter
import datetime
import json
import math
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET
import zlib

import numpy as np

from analyze_ubash3a_proteomics_specificity import sha256, write_csv

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/ubash3a-proteomics-spectra"
NS = {"m": "http://psi.hupo.org/ms/mzml"}


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def terms(element):
    return {term.attrib["accession"]: term.attrib for term in element.findall("m:cvParam", NS)}


def decode_array(element, expected_length):
    cv = terms(element)
    if "MS:1000523" in cv:
        dtype = np.dtype("<f8")
    elif "MS:1000521" in cv:
        dtype = np.dtype("<f4")
    else:
        raise ValueError("Unsupported or missing binary precision")
    text = element.find("m:binary", NS).text or ""
    data = base64.b64decode("".join(text.split()), validate=True)
    if "MS:1000574" in cv:
        data = zlib.decompress(data)
    elif "MS:1000576" not in cv:
        raise ValueError("Unsupported or missing binary compression")
    array = np.frombuffer(data, dtype=dtype)
    if len(array) != expected_length or not np.all(np.isfinite(array)):
        raise ValueError("Binary array length or finite-value check failed")
    return array


def inspect_mzml(path):
    ids, rows = set(), []
    listed_count = None
    total_points, nonempty, centroid = 0, 0, 0
    charges, activation, filters = Counter(), Counter(), Counter()
    selected_examples = {}
    last_rt = -math.inf
    for event, element in ET.iterparse(path, events=("start", "end")):
        local_tag = element.tag.rsplit("}", 1)[-1]
        if event == "start" and local_tag == "spectrumList":
            listed_count = int(element.attrib["count"])
        if event != "end" or local_tag != "spectrum":
            continue
        identifier = element.attrib["id"]
        assert identifier not in ids, identifier
        ids.add(identifier)
        cv = terms(element)
        assert cv["MS:1000511"]["value"] == "2", identifier
        centroid += "MS:1000127" in cv
        scans = re.findall(r"(?:^|\s)scan=(\d+)(?:$|\s)", identifier)
        assert len(scans) == 1
        scan_id = int(scans[0])
        scan_element = element.find("m:scanList/m:scan", NS)
        scan_cv = terms(scan_element)
        rt_term = scan_cv["MS:1000016"]
        rt = float(rt_term["value"])
        if rt_term.get("unitAccession") == "UO:0000031":
            rt *= 60
        else:
            assert rt_term.get("unitAccession") == "UO:0000010", rt_term
        assert math.isfinite(rt) and rt >= last_rt
        last_rt = rt
        filter_value = scan_cv.get("MS:1000512", {}).get("value", "")
        filters[filter_value.split(" ", 1)[0] if filter_value else "missing"] += 1
        ions = element.findall("m:precursorList/m:precursor/m:selectedIonList/m:selectedIon", NS)
        assert len(ions) == 1, (identifier, len(ions))
        ion_cv = terms(ions[0])
        precursor = float(ion_cv["MS:1000744"]["value"])
        assert math.isfinite(precursor) and precursor > 0
        charge = int(ion_cv["MS:1000041"]["value"]) if "MS:1000041" in ion_cv else 0
        assert charge >= 0
        charges[charge] += 1
        activation_element = element.find("m:precursorList/m:precursor/m:activation", NS)
        activation_cv = terms(activation_element)
        activation_names = [term["name"] for accession, term in activation_cv.items() if accession != "MS:1000045"]
        activation_name = " | ".join(sorted(activation_names))
        activation[activation_name or "missing"] += 1
        arrays = {}
        expected = int(element.attrib["defaultArrayLength"])
        for array_element in element.findall("m:binaryDataArrayList/m:binaryDataArray", NS):
            array_cv = terms(array_element)
            if "MS:1000514" in array_cv:
                arrays["mz"] = decode_array(array_element, expected)
            elif "MS:1000515" in array_cv:
                arrays["intensity"] = decode_array(array_element, expected)
        assert set(arrays) == {"mz", "intensity"}
        mz, intensity = arrays["mz"], arrays["intensity"]
        assert np.all(mz > 0) and np.all(intensity >= 0)
        assert np.all(np.diff(mz) >= 0)
        nonempty += bool(expected)
        total_points += expected
        rows.append({"file": path.stem, "scan": scan_id, "spectrum_id": identifier,
                     "retention_seconds": rt, "precursor_mz": precursor, "charge": charge,
                     "peak_count": expected, "total_intensity": float(intensity.sum()),
                     "activation": activation_name, "instrument_filter": filter_value})
        # Keep the first three spectra as a predeclared format cross-check set.
        if len(selected_examples) < 3:
            selected_examples[scan_id] = {"mz": mz.tolist(), "intensity": intensity.tolist(),
                                          "precursor_mz": precursor, "retention_seconds": rt, "charge": charge}
        element.clear()
    assert listed_count is not None and listed_count == len(ids) == len(rows), (listed_count, len(ids))
    assert len({row["scan"] for row in rows}) == len(rows)
    assert centroid == len(rows), (centroid, len(rows))
    assert len(rows) > 0
    return {"spectra": len(rows), "declared_spectrum_count": listed_count,
            "nonempty_spectra": nonempty, "centroid_spectra": centroid,
            "total_fragment_points": total_points, "charges": dict(sorted(charges.items())),
            "activation": dict(activation), "filter_analyzers": dict(filters),
            "first_scan": rows[0]["scan"], "last_scan": rows[-1]["scan"],
            "binary_and_scan_checks_pass": True}, rows, selected_examples


def execute_one(record, binary):
    raw = ROOT / record["local_path"]
    assert raw.stat().st_size == record["downloaded_bytes"] and sha256(raw) == record["sha256"]
    output_dir = WORK / "converted"
    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir = WORK / "conversion-logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    mzml = output_dir / (raw.stem + ".mzML")
    metadata_path = output_dir / (raw.stem + "-metadata.json")
    started = now()
    # The publisher's development-build metadata writer repeats a complete
    # CountScanOrder pass inside every scan iteration. Use the separately pinned
    # single-pass metadata reader; the spectrum conversion itself is unchanged.
    args = [str(binary), "-i", str(raw), "-o", str(output_dir), "-f", "1", "-m", "2", "-L", "2"]
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    stdout_path, stderr_path = log_dir / (raw.name + ".stdout.txt"), log_dir / (raw.name + ".stderr.txt")
    stdout_path.write_text(result.stdout)
    stderr_path.write_text(result.stderr)
    if result.returncode:
        failure = {"file": raw.name, "args": args, "exit_code": result.returncode, "started_at_utc": started, "ended_at_utc": now()}
        (log_dir / (raw.name + ".failure.json")).write_text(json.dumps(failure, indent=2) + "\n")
        raise RuntimeError(f"Conversion failed: {raw.name}; inspect saved logs")
    repair = json.loads((ROOT / "config/ubash3a-proteomics-metadata-repair.json").read_text())
    validator = ROOT / repair["validator_dll"]
    assert sha256(validator) == repair["validator_sha256"]
    for dependency, expected in repair["runtime_dependency_sha256"].items():
        assert sha256(ROOT / dependency) == expected, dependency
    metadata_args = ["dotnet", str(validator), str(raw), str(metadata_path)]
    metadata_result = subprocess.run(metadata_args, cwd=ROOT, capture_output=True, text=True)
    (log_dir / (raw.name + ".metadata.stdout.txt")).write_text(metadata_result.stdout)
    (log_dir / (raw.name + ".metadata.stderr.txt")).write_text(metadata_result.stderr)
    if metadata_result.returncode:
        raise RuntimeError(f"Independent metadata loop failed for {raw.name}: {metadata_result.stderr}")
    metadata = json.loads(metadata_path.read_text())
    vendor_count = int(metadata["countsByMSOrder"]["2"])
    qc, rows, examples = inspect_mzml(mzml)
    assert vendor_count == qc["spectra"], (raw.name, vendor_count, qc["spectra"])
    assert metadata["ms2ScanNumbers"] == [row["scan"] for row in rows]
    expected_analyzer = set(qc["filter_analyzers"]) == {"FTMS"}
    assert expected_analyzer, (raw.name, qc["filter_analyzers"])
    assert set(qc["activation"]) == {"beam-type collision-induced dissociation"}, qc["activation"]
    write_csv(log_dir / (raw.name + ".scan-inventory.csv"), rows)
    (log_dir / (raw.name + ".first-three-spectra.json")).write_text(json.dumps(examples) + "\n")
    report = {"raw_file": raw.name, "raw_sha256": record["sha256"], "started_at_utc": started,
              "completed_at_utc": now(), "args": args, "exit_code": result.returncode,
              "mzml_path": str(mzml.relative_to(ROOT)), "mzml_sha256": sha256(mzml), "mzml_bytes": mzml.stat().st_size,
              "metadata_path": str(metadata_path.relative_to(ROOT)), "metadata_sha256": sha256(metadata_path),
              "vendor_MS2_count": vendor_count, "metadata_args": metadata_args,
              "metadata_method": metadata["metadataMethod"], "instrument_model": metadata["instrumentModel"],
              "vendor_count_equals_independent_mzml_count": True,
              "all_vendor_MS2_scan_identifiers_equal_mzml": True,
              "stdout_sha256": sha256(stdout_path), "stderr_sha256": sha256(stderr_path),
              "warning_lines": [line for line in (result.stdout + result.stderr).splitlines() if re.search(r"\bWARN(?:ING)?\b", line, re.I)],
              "completion_status_lines": [line for line in (result.stdout + result.stderr).splitlines() if "Processing completed" in line],
              **qc}
    (log_dir / (raw.name + ".qualification.json")).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"file": raw.name, "converted_and_checked": True, "spectra": qc["spectra"]}), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--first-only", action="store_true", help="Qualify the first frozen file before completing the set")
    args = parser.parse_args()
    tools = json.loads((ROOT / "config/ubash3a-proteomics-tool-install.json").read_text())
    tool = tools["sources"][0]
    binary = ROOT / tool["installed_binary"]
    assert sha256(binary) == tool["binary_sha256"]
    download_path = ROOT / "config/ubash3a-proteomics-spectrum-downloads.json"
    if args.first_only and not download_path.exists():
        selection = json.loads((ROOT / "config/ubash3a-proteomics-spectrum-selection.json").read_text())
        first_name = selection["files"][0]["file_name"]
        first_receipt = json.loads((WORK / "download-receipts" / (first_name + ".json")).read_text())
        assert first_receipt["complete"]
        chosen = [first_receipt]
    else:
        downloads = json.loads(download_path.read_text())
        assert downloads["complete"] and len(downloads["files"]) == 14
        chosen = downloads["files"][:1] if args.first_only else downloads["files"]
    results = []
    for file in chosen:
        existing = WORK / "conversion-logs" / (file["file_name"] + ".qualification.json")
        if existing.exists():
            prior = json.loads(existing.read_text())
            assert prior["raw_sha256"] == file["sha256"] and sha256(ROOT / prior["mzml_path"]) == prior["mzml_sha256"]
            results.append(prior)
        else:
            results.append(execute_one(file, binary))
    output = {"complete_experiment": not args.first_only, "files": results, "tool_binary_sha256": sha256(binary),
              "completed_at_utc": now(), "total_MS2_spectra": sum(row["spectra"] for row in results),
              "all_conversion_and_content_checks_pass": True,
              "qualification_limitation": "Independent mzML decoding and counts; both converter and metadata counts use the vendor reader, not an independent instrument measurement",
              "metadata_export_repair": "A separate single-pass vendor-metadata loop avoids the development writer's repeated whole-file scan loop and mislabeled MS2 field. All MS2 scan IDs and counts match independently parsed mzML; publisher code and repair are pinned."}
    name = "ubash3a-proteomics-first-conversion.json" if args.first_only else "ubash3a-proteomics-conversion.json"
    (ROOT / "reports" / name).write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"files": len(results), "total_MS2": output["total_MS2_spectra"], "complete_experiment": output["complete_experiment"]}), flush=True)
