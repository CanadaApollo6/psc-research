#!/usr/bin/env python3
"""Account for physical MS2 scan IDs and exported queries without absence claims."""
import csv
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from analyze_ubash3a_proteomics_specificity import sha256, write_csv

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/ubash3a-proteomics-spectra"
OUT = ROOT / "data/derived/ubash3a-proteomics-spectra"


def main():
    execution = json.loads((ROOT / "reports/ubash3a-proteomics-search-execution.json").read_text())
    assert execution["complete_experiment"] and execution["file_count"] == 14
    files, unassigned = [], []
    for item in execution["files"]:
        stem = Path(item["input_mzml"]).stem
        xml_name = next(name for name in item["output_sha256"] if name.endswith(".pep.xml"))
        xml_path = ROOT / xml_name
        assert sha256(xml_path) == item["output_sha256"][xml_name]
        stdout_name = next(name for name in item["log_sha256"] if name.endswith(".stdout.txt"))
        stdout_path = ROOT / stdout_name
        assert sha256(stdout_path) == item["log_sha256"][stdout_name]
        batches = [int(value) for value in re.findall(r"Load spectra:\s*(\d+)", stdout_path.read_text())]
        inventory_path = WORK / "conversion-logs" / (stem + ".raw.scan-inventory.csv")
        inventory = list(csv.DictReader(inventory_path.open()))
        ids = {int(row["scan"]) for row in inventory}
        assert len(ids) == len(inventory)
        queries, hit_queries = set(), set()
        nonidentical_scan_labels = 0
        maximum_mz_difference = maximum_retention_difference = 0.0
        for _, query in ET.iterparse(xml_path, events=("end",)):
            if query.tag.rsplit("}", 1)[-1] != "spectrum_query":
                continue
            comet_scan = int(query.attrib["start_scan"])
            matched = re.search(r"(?:^|\s)scan=([0-9]+)(?=\s|$)", query.attrib["spectrumNativeID"])
            assert matched is not None
            original_scan = int(matched.group(1))
            key = original_scan, int(query.attrib["assumed_charge"])
            assert query.attrib["start_scan"] == query.attrib["end_scan"]
            assert key not in queries and key[0] in ids
            # In this complete MS2-only conversion, the exported scan is the
            # one-based mzML spectrum position. Check rather than assume it.
            original = inventory[comet_scan - 1]
            assert int(original["scan"]) == original_scan
            assert original["spectrum_id"] == query.attrib["spectrumNativeID"]
            mz_implied = float(query.attrib["precursor_neutral_mass"]) / key[1] + 1.007276466621
            mz_difference = abs(mz_implied - float(original["precursor_mz"]))
            rt_difference = abs(float(query.attrib["retention_time_sec"]) - float(original["retention_seconds"]))
            assert mz_difference <= 0.001, (stem, comet_scan, mz_difference)
            assert rt_difference <= 0.051, (stem, comet_scan, rt_difference)
            assert int(original["charge"]) in (0, key[1])
            maximum_mz_difference = max(maximum_mz_difference, mz_difference)
            maximum_retention_difference = max(maximum_retention_difference, rt_difference)
            nonidentical_scan_labels += comet_scan != original_scan
            queries.add(key)
            if any(child.tag.rsplit("}", 1)[-1] == "search_hit" for child in query.iter()):
                hit_queries.add(key)
            query.clear()
        exported_scans = {scan for scan, _ in queries}
        assigned_scans = {scan for scan, _ in hit_queries}
        assert len(queries) <= sum(batches)
        for row in inventory:
            scan = int(row["scan"])
            if scan not in assigned_scans:
                unassigned.append({**row, "has_exported_query": scan in exported_scans,
                                   "interpretation": "No exported peptide assignment; filtering or no reported competitor is unresolved, not evidence of protein absence"})
        files.append({"file": stem, "original_converted_MS2": len(ids),
                      "engine_loaded_queries": sum(batches), "loaded_batches_json": json.dumps(batches),
                      "exported_XML_queries": len(queries), "exported_XML_queries_with_hits": len(hit_queries),
                      "exported_XML_queries_without_hits": len(queries - hit_queries),
                      "loaded_queries_not_exported": sum(batches) - len(queries),
                      "physical_MS2_with_exported_query": len(exported_scans),
                      "physical_MS2_with_reported_assignment": len(assigned_scans),
                      "physical_MS2_without_reported_assignment": len(ids - assigned_scans),
                      "all_exported_scan_ids_in_original_inventory": True,
                      "all_exported_scan_ordinals_match_original_mzML_positions": True,
                      "queries_with_different_exported_and_instrument_scan_labels": nonidentical_scan_labels,
                      "maximum_precursor_mz_difference": maximum_mz_difference,
                      "maximum_retention_time_difference_seconds": maximum_retention_difference,
                      "scan_inventory_sha256": sha256(inventory_path),
                      "XML_sha256": sha256(xml_path)})
    small_path = OUT / "scan-coverage-audit.csv"
    missing_path = WORK / "analysis/scans-without-exported-assignment.csv"
    write_csv(small_path, files)
    fields = list(unassigned[0]) if unassigned else ["file", "scan", "has_exported_query", "interpretation"]
    write_csv(missing_path, unassigned, fields)
    summary = {"all_checks_pass": True, "files": len(files),
               "original_converted_MS2": sum(row["original_converted_MS2"] for row in files),
               "engine_loaded_queries": sum(row["engine_loaded_queries"] for row in files),
               "exported_XML_queries": sum(row["exported_XML_queries"] for row in files),
               "loaded_queries_not_exported": sum(row["loaded_queries_not_exported"] for row in files),
               "physical_MS2_with_reported_assignment": sum(row["physical_MS2_with_reported_assignment"] for row in files),
               "physical_MS2_without_reported_assignment": len(unassigned),
               "scan_identifiers_accounted_for_without_assignment": True,
               "all_exported_scan_ordinals_native_IDs_and_precursor_retention_fields_checked": True,
               "maximum_precursor_mz_difference": max(row["maximum_precursor_mz_difference"] for row in files),
               "maximum_retention_time_difference_seconds": max(row["maximum_retention_time_difference_seconds"] for row in files),
               "field_precision_tolerances": {"precursor_mz": 0.001, "retention_seconds": 0.051},
               "absence_claim": False,
               "limitation": "Search output cannot uniquely assign every missing result to preprocessing or lack of a reported competitor; loaded queries and physical MS2 scans are different units.",
               "output_hashes": {str(path.relative_to(ROOT)): sha256(path) for path in (small_path, missing_path)}}
    assert summary["original_converted_MS2"] == summary["physical_MS2_with_reported_assignment"] + len(unassigned)
    (ROOT / "reports/ubash3a-proteomics-scan-coverage.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
