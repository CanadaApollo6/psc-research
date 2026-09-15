#!/usr/bin/env python3
"""Qualify whole files as frozen downloads complete; aggregate only all fourteen."""
import datetime
import json
from pathlib import Path
import subprocess
import time

from analyze_ubash3a_proteomics_specificity import sha256
from qualify_ubash3a_spectra import execute_one

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/ubash3a-proteomics-spectra"


def main():
    selection = json.loads((ROOT / "config/ubash3a-proteomics-spectrum-selection.json").read_text())
    assert len(selection["files"]) == 14
    tool = json.loads((ROOT / "config/ubash3a-proteomics-tool-install.json").read_text())["sources"][0]
    binary = ROOT / tool["installed_binary"]
    assert sha256(binary) == tool["binary_sha256"]
    for index, item in enumerate(selection["files"]):
        receipt_path = WORK / "download-receipts" / (item["file_name"] + ".json")
        while not receipt_path.exists():
            time.sleep(5)
        receipt = json.loads(receipt_path.read_text())
        assert receipt["complete"] and receipt["downloaded_bytes"] == item["expected_bytes"]
        existing = WORK / "conversion-logs" / (item["file_name"] + ".qualification.json")
        if existing.exists():
            result = json.loads(existing.read_text())
            assert result["raw_sha256"] == receipt["sha256"] and sha256(ROOT / result["mzml_path"]) == result["mzml_sha256"]
        else:
            result = execute_one(receipt, binary)
        print(f"Qualified complete file {index + 1}/14: {item['file_name']}", flush=True)
    while not (ROOT / "config/ubash3a-proteomics-spectrum-downloads.json").exists():
        time.sleep(2)
    subprocess.run([str(ROOT / ".venv/bin/python"), str(ROOT / "scripts/qualify_ubash3a_spectra.py")], cwd=ROOT, check=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        (WORK / "conversion-pipeline-failure.json").write_text(json.dumps({
            "failed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "error": repr(exc), "all_fourteen_complete": False,
        }, indent=2) + "\n")
        raise
