#!/usr/bin/env python3
"""Verify deterministic primary tables from the complete pinned local cache.

Small outputs are written to a separate directory. The primary spectrum
analyzer deterministically rebuilds its two large local tables in place;
both are hashed before and after. This does not rerun the search engine.
"""
import argparse
import contextlib
import datetime
import io
import json
from pathlib import Path

import analyze_ubash3a_proteomics_specificity as specificity
import analyze_ubash3a_proteomics_spectra as spectra

ROOT = Path(__file__).resolve().parents[1]


def compare_tree(original, replay):
    checks = []
    for path in sorted(replay.rglob("*")):
        if not path.is_file():
            continue
        source = original / path.relative_to(replay)
        assert source.is_file(), source
        expected, actual = specificity.sha256(source), specificity.sha256(path)
        assert expected == actual, (source, path)
        checks.append({"original": str(source.relative_to(ROOT)),
                       "replay": str(path.relative_to(ROOT)), "sha256": expected})
    return checks


def main(output):
    execution_path = ROOT / "reports/ubash3a-proteomics-search-execution.json"
    execution = json.loads(execution_path.read_text())
    assert execution["complete_experiment"] and execution["file_count"] == 14
    output = output.resolve()
    assert output.is_relative_to(ROOT / "work"), "Replay directory must stay in ignored local work storage"
    output.mkdir(parents=True, exist_ok=True)
    assert not any(output.iterdir()), "Preserve an existing replay; choose a new output directory"
    large = [ROOT / "work/ubash3a-proteomics-spectra/analysis" / name for name in
             ("all-spectrum-winners.csv", "all-peptide-winners.csv")]
    before = {str(path.relative_to(ROOT)): specificity.sha256(path) for path in large}
    sequence_output = output / "specificity"
    sequence_output.mkdir()
    metadata = specificity.audit_metadata(sequence_output)
    (sequence_output / "metadata-summary.json").write_text(json.dumps(metadata, indent=2) + "\n")
    specificity.analyze(sequence_output)
    sequence_checks = compare_tree(ROOT / "data/derived/ubash3a-proteomics-specificity", sequence_output)
    assert len(sequence_checks) == 19
    original_output = spectra.OUT
    try:
        spectra.OUT = output / "spectra"
        with contextlib.redirect_stdout(io.StringIO()):
            spectra.main()
    finally:
        spectra.OUT = original_output
    spectrum_checks = compare_tree(original_output, output / "spectra")
    assert len(spectrum_checks) == 5
    after = {str(path.relative_to(ROOT)): specificity.sha256(path) for path in large}
    assert before == after
    result = {"all_checks_pass": True,
              "completed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "specificity_primary_files_replayed_exactly": len(sequence_checks),
              "spectrum_primary_files_replayed_exactly": len(spectrum_checks),
              "small_file_checks": sequence_checks + spectrum_checks,
              "large_tables_rebuilt_with_identical_bytes": after,
              "scope": "Deterministic table reconstruction from existing pinned references and complete search outputs; no repeated spectrum search, fragment review, reference-lock creation or biological replication"}
    report = ROOT / "reports/ubash3a-proteomics-replay-verification.json"
    report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "small_file_checks"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=ROOT / "work/ubash3a-proteomics-spectra/deterministic-replay")
    main(parser.parse_args().output)
