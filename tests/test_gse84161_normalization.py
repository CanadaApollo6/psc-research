"""Guard and source-identity tests for the GSE84161 array runner."""

from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = ROOT / "work/gse84161-normalization/runtime/bin/Rscript"
SCRIPT = ROOT / "scripts/normalize_gse84161.R"
RUNTIME_PLAN = ROOT / "config/gse84161-array-runtime.json"


class TestGSE84161NormalizationGuard(unittest.TestCase):
    def test_summary_csv_and_entrez_medians_on_synthetic_inputs(self) -> None:
        result = subprocess.run(
            [str(RSCRIPT), str(ROOT / "tests/test_gse84161_normalization.R")],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("regression checks passed", result.stdout)

    def run_runner(self, mode: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(RSCRIPT), str(SCRIPT), "--plan", str(RUNTIME_PLAN), mode],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_check_only_preserves_all_thirty_metadata_ordered_inputs(self) -> None:
        result = self.run_runner("--check-only")
        self.assertEqual(result.returncode, 0, result.stderr)
        checked = json.loads(result.stdout)
        self.assertEqual(checked["samples"]["count"], 30)
        self.assertEqual(checked["samples"]["accessions"][0], "GSM2228386")
        self.assertEqual(checked["samples"]["accessions"][-1], "GSM2228415")
        self.assertEqual(checked["source"]["archive_member_count"], 30)
        self.assertTrue(checked["source"]["extracted_all_present"])
        self.assertFalse(checked["expression_values_read"])

    def test_execute_requires_frozen_input_gate(self) -> None:
        result = self.run_runner("--execute")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('normalization_run_status="frozen"', result.stderr)


if __name__ == "__main__":
    unittest.main()
