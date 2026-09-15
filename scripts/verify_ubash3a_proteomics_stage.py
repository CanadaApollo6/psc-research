#!/usr/bin/env python3
"""Final integrity, scope, regression and preservation checks for direction 3."""
import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urlparse

from analyze_ubash3a_proteomics_specificity import sha256

ROOT = Path(__file__).resolve().parents[1]
BASE = "c7bccc3c84cf8559b92759cfcec063d3f194cecf"
NAVIGATION = {"README.md", "docs/research-plan.md", "docs/data-sources.md",
              "docs/evidence-log.md", "docs/computational-followup-roadmap.md"}
REPORT = "reports/ubash3a-proteomics-stage-verification.json"


def read(name):
    return json.loads((ROOT / name).read_text())


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def task_file(name):
    return (name in NAVIGATION or name in {"requirements-ubash3a-proteomics.txt", "scripts/review_ubash3a_candidate_fragments.py"}
            or name.startswith(("config/ubash3a-proteomics-", "docs/ubash3a-proteomics-",
                                "reports/ubash3a-proteomics-", "data/derived/ubash3a-proteomics-",
                                "scripts/ProteomicsRawMetadata/", "scripts/ubash3a_proteomics_sources/"))
            or (name.startswith(("scripts/", "tests/")) and name.endswith(".py")
                and ("ubash3a_proteomics" in name or "ubash3a_spectra" in name or "ubash3a_comet" in name)))


def verify():
    baseline = {}
    for item in git("ls-tree", "-rz", BASE).split(b"\0"):
        if not item:
            continue
        metadata, name = item.split(b"\t", 1)
        mode, kind, oid = metadata.split()
        assert kind == b"blob"
        baseline[name.decode()] = oid.decode()
    changed = []
    for name, expected in baseline.items():
        data = (ROOT / name).read_bytes()
        actual = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if actual != expected:
            changed.append(name)
    assert set(changed) == NAVIGATION, changed
    assert len(baseline) == 775

    candidates = set(git("diff", "--name-only", "-z", BASE).decode().split("\0"))
    candidates.update(git("ls-files", "--others", "--exclude-standard", "-z").decode().split("\0"))
    candidates.discard("")
    candidates.add(REPORT)
    assert all(task_file(name) for name in candidates), sorted(name for name in candidates if not task_file(name))
    assert not any(name.startswith(("work/", "data/raw/", "data/private/")) for name in candidates)

    lock = read("config/ubash3a-proteomics-reference-lock.json")
    for source in lock["sources"]:
        assert sha256(ROOT / source["local_path"]) == source["sha256"]
    pagination = lock["uniprot_full_page_receipt"]
    for page in pagination["pages"]:
        assert sha256(ROOT / page["local_path"]) == page["sha256"]
    assert len(pagination["pages"]) == 466
    downloads = read("config/ubash3a-proteomics-spectrum-downloads.json")
    assert downloads["complete"] and len(downloads["files"]) == 14
    for file in downloads["files"]:
        path = ROOT / file["local_path"]
        assert path.stat().st_size == file["downloaded_bytes"] == file["expected_bytes"]
        assert sha256(path) == file["sha256"]
    execution = read("reports/ubash3a-proteomics-search-execution.json")
    assert execution["complete_experiment"] and execution["file_count"] == 14
    for file in execution["files"]:
        assert file["complete"] and file["exit_code"] == 0
        assert sha256(ROOT / file["input_mzml"]) == file["input_sha256"]
        for field in ("output_sha256", "log_sha256"):
            for name, expected in file[field].items():
                assert sha256(ROOT / name) == expected

    start = read("config/ubash3a-proteomics-analysis-start.json")
    assert start["all_14_searches_complete_before_any_new_candidate_score_inspection"]
    assert sha256(ROOT / "reports/ubash3a-proteomics-search-execution.json") == start["complete_search_execution_sha256"]
    assert sha256(ROOT / "config/ubash3a-proteomics-sources.json") == start["source_manifest_sha256"]
    assert sha256(ROOT / "config/ubash3a-proteomics-comet.params") == start["parameter_sha256"]
    format_repair = read("config/ubash3a-proteomics-verifier-format-repair.json")
    for name, expected in format_repair["primary_outputs_before_repair"].items():
        assert sha256(ROOT / name) == expected
    for source in format_repair["format_source_receipts"]:
        assert sha256(ROOT / source["local_path"]) == source["sha256"]

    required_reports = ["reports/ubash3a-proteomics-specificity-verification.json",
                        "reports/ubash3a-proteomics-comet-qualification.json",
                        "reports/ubash3a-proteomics-search-verification.json",
                        "reports/ubash3a-proteomics-scan-coverage.json",
                        "reports/ubash3a-proteomics-replay-verification.json"]
    for name in required_reports:
        assert read(name)["all_checks_pass"], name
    conversion = read("reports/ubash3a-proteomics-conversion.json")
    assert conversion["complete_experiment"] and conversion["all_conversion_and_content_checks_pass"]
    assert sum(len(row["warning_lines"]) for row in conversion["files"]) == 0
    summary = read("data/derived/ubash3a-proteomics-spectra/search-summary.json")
    coverage = read("reports/ubash3a-proteomics-scan-coverage.json")
    independent = read("reports/ubash3a-proteomics-search-verification.json")
    fragments = read("reports/ubash3a-proteomics-fragment-review.json")
    assert summary["converted_MS2"] == coverage["original_converted_MS2"] == conversion["total_MS2_spectra"] == 66766
    assert summary["Comet_loaded_spectrum_queries"] == coverage["engine_loaded_queries"] == coverage["exported_XML_queries"] == 55741
    assert summary["spectrum_winners"] == coverage["physical_MS2_with_reported_assignment"] == independent["search"]["independent_spectrum_winners"] == 54278
    assert summary["all_reported_assignments"] == independent["search"]["TXT_XML_assignment_rows"] == 269986
    assert summary["candidate_reported_assignments"] == fragments["candidate_assignment_count"] == independent["search"]["candidate_assignments_verified"] == 0
    assert summary["spectrum_q01_target_wins"] == independent["search"]["independent_spectrum_q01_targets"] == 11611
    assert summary["ordinary_reference_UBASH3A_spectrum_q01_wins"] == 0
    for field in (summary["large_output_hashes"], fragments["files"], coverage["output_hashes"]):
        for name, expected in field.items():
            assert sha256(ROOT / name) == expected

    test_status = read("work/ubash3a-proteomics-spectra/post-parser-repair-test-suite-status.json")
    assert test_status["exit_code"] == 0
    assert sha256(ROOT / test_status["log"]) == test_status["log_sha256"]
    text = (ROOT / test_status["log"]).read_text()
    count = re.search(r"Ran (\d+) tests", text)
    assert count and int(count.group(1)) == 238 and text.rstrip().endswith("OK")

    checked_links = 0
    for name in sorted(candidates):
        if not name.endswith(".md"):
            continue
        for target in re.findall(r"\]\(([^)]+)\)", (ROOT / name).read_text()):
            target = target.strip("<>")
            if urlparse(target).scheme or target.startswith("#"):
                continue
            target = unquote(target.split("#", 1)[0])
            if not target:
                continue
            resolved = (ROOT / name).parent / target
            assert resolved.exists() or resolved.resolve() == ROOT / REPORT, (name, target)
            checked_links += 1
    subprocess.run(["git", "diff", "--check"], cwd=ROOT, check=True)
    files = {name: {"sha256": sha256(ROOT / name), "size_bytes": (ROOT / name).stat().st_size}
             for name in sorted(candidates) if name != REPORT}
    assert all(value["size_bytes"] < 2_000_000 for value in files.values())
    report = {"all_checks_pass": True, "complete_bounded_roadmap_direction": 3,
              "verified_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "baseline_commit": BASE, "baseline_tracked_files": len(baseline),
              "prior_files_preserved_byte_for_byte": len(baseline) - len(changed),
              "intentional_navigation_changes": sorted(changed),
              "source_reference_records": lock["total_records"], "UniProt_pages_rechecked": len(pagination["pages"]),
              "complete_raw_files_rechecked": len(downloads["files"]),
              "complete_raw_bytes_rechecked": sum(file["downloaded_bytes"] for file in downloads["files"]),
              "all_14_search_outputs_and_input_hashes_rechecked": True,
              "original_analysis_start_source_manifest_and_parameters_preserved": True,
              "primary_output_bytes_unchanged_by_verifier_format_repair": True,
              "candidate_assignments": 0, "ordinary_UBASH3A_spectrum_q01_assignments": 0,
              "biological_absence_claim": False, "new_laboratory_or_personal_genome_work": False,
              "software_tests": {"passed": 238, **test_status},
              "verified_report_hashes": {name: sha256(ROOT / name) for name in required_reports},
              "relative_markdown_links_checked": checked_links,
              "delivery_file_count_including_this_record": len(candidates),
              "delivery_files": files, "self_excluded_from_file_hashes": REPORT,
              "completion_note": "Both failed analysis/verifier attempts remain preserved as historical logs; completed independent verification, scan coverage and replay records define the final status. No spectrum search or biological threshold was changed to obtain the result."}
    (ROOT / REPORT).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key not in {"delivery_files", "verified_report_hashes"}}, indent=2))


if __name__ == "__main__":
    verify()
