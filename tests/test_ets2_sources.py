"""Focused tests for the deterministic ETS2 source audit."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import audit_ets2_sources as audit  # noqa: E402


class TestETS2Sources(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = audit.build_audit(write_outputs=True)

    def test_both_manifests_and_public_md5_verify(self) -> None:
        manifests = self.report["source_manifests"]
        self.assertEqual(
            [item["path"] for item in manifests],
            ["config/ets2-benchmark-sources.json", "config/ets2-program-sources.json"],
        )
        for manifest in manifests:
            self.assertTrue(all(source["verified"] for source in manifest["sources"]))
        self.assertTrue(self.report["publisher_zip_md5"]["verified"])
        self.assertEqual(
            self.report["publisher_zip_md5"]["actual"],
            "9d6b7ddcf377cb194bc72732d807ae80",
        )

    def test_master_sets_reproduce_strict_fdr05_membership(self) -> None:
        expected = {
            "ETS2g1_UP": 668,
            "ETS2g1_DN": 928,
            "CHR21_DN": 141,
            "CHR21_UP": 215,
            "ETS2g2_UP": 716,
            "ETS2g2_DN": 875,
        }
        reproduction = self.report["gene_sets"]["fdr05_reproduction"]
        for name, count in expected.items():
            self.assertEqual(reproduction[name]["master_count"], count)
            self.assertEqual(reproduction[name]["reproduced_fdr05_count"], count)
            self.assertTrue(reproduction[name]["exact_member_set"])
        self.assertEqual(
            self.report["gene_sets"]["fdr_rule"],
            "strict adjusted P value < 0.05; direction determined by source logFC",
        )
        self.assertEqual(
            self.report["gene_sets"]["source_logfc_t_sign_agreement"]["s1_g1"],
            {"agree": 12144, "disagree": 0},
        )

    def test_workbook_literal_none_duplicates_and_footer_are_preserved(self) -> None:
        workbooks = self.report["workbooks"]
        self.assertEqual(workbooks["s1"]["valid_ensembl_rows"], 12144)
        self.assertEqual(workbooks["s1"]["worksheet_body_rows_after_header"], 12149)
        self.assertEqual(workbooks["s1"]["literal_none_rows"], 481)
        self.assertEqual(workbooks["s1"]["duplicate_non_none_symbols"]["Y_RNA"], 3)
        self.assertEqual(workbooks["s2"]["valid_ensembl_rows"], 10196)
        self.assertEqual(workbooks["s2"]["worksheet_body_rows_after_header"], 10201)
        self.assertEqual(workbooks["s2"]["literal_none_rows"], 199)
        self.assertEqual(workbooks["s2"]["duplicate_non_none_symbols"]["U2"], 10)

        notes_path = ROOT / "data/derived/ets2-source-notes.csv"
        with notes_path.open(encoding="utf-8", newline="") as stream:
            notes = list(csv.DictReader(stream))
        self.assertEqual(len(notes), 8)
        self.assertIn("reverse complement of ETS2", notes[-1]["value_original"])

    def test_rank_provenance_uses_one_to_one_identity(self) -> None:
        by_id = {item["rank_dataset_id"]: item for item in self.report["rank_datasets"]}
        ko = by_id["ko_g1"]["source_comparison"]
        self.assertEqual(ko["primary_unambiguous_match_count"], 11658)
        self.assertEqual(ko["primary_numeric_match_within_tolerance_count"], 11658)
        self.assertEqual(ko["primary_sign_mismatch_count"], 0)

        oe250 = by_id["oe_250"]["source_comparison"]
        oe500 = by_id["oe_500"]["source_comparison"]
        for comparison, mismatches in ((oe250, 36), (oe500, 42)):
            self.assertEqual(comparison["primary_unambiguous_match_count"], 9908)
            self.assertEqual(comparison["primary_numeric_match_within_tolerance_count"], 0)
            self.assertEqual(comparison["primary_sign_mismatch_count"], mismatches)
            self.assertEqual(comparison["ambiguous_shared_symbol_count"], 5)
            self.assertEqual(comparison["rank_only_symbol_count"], 66)

        self.assertEqual(
            by_id["meki_100"]["source_comparison"]["primary_unambiguous_match_count"],
            None,
        )
        self.assertEqual(
            by_id["meki_500"]["source_comparison_status"],
            "publisher_differential_table_not_pinned",
        )

    def test_rank_output_retains_all_five_vectors_and_order(self) -> None:
        path = ROOT / "data/derived/ets2-source-ranks.csv.gz"
        with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 11660 + 9979 + 9979 + 12095 + 12095)
        by_id = {}
        for row in rows:
            by_id.setdefault(row["rank_dataset_id"], []).append(row)
        for dataset_rows in by_id.values():
            self.assertEqual(
                [int(row["rank_position"]) for row in dataset_rows],
                list(range(1, len(dataset_rows) + 1)),
            )
            self.assertEqual(
                len({row["gene_symbol_original"] for row in dataset_rows}),
                len(dataset_rows),
            )
        self.assertEqual(by_id["ko_g1"][0]["gene_symbol_original"], "FSTL3")
        self.assertEqual(by_id["ko_g1"][-1]["gene_symbol_original"], "ETS2")

    def test_design_counts_and_literal_g4_factor_mismatch(self) -> None:
        designs = self.report["designs"]
        self.assertEqual(designs["crispr_ets2"]["all_names_count"], 38)
        self.assertEqual(designs["crispr_ets2"]["selected_model_sample_count"], 27)
        self.assertEqual(designs["crispr_ets2"]["complete_g1_ntc_pair_donors"], [2, 3, 5, 6, 7, 8, 9, 10])
        self.assertEqual(designs["oe_ets2"]["sample_count"], 32)
        self.assertEqual(designs["meki"]["sample_count"], 9)

        with (ROOT / "data/derived/ets2-source-design.csv").open(
            encoding="utf-8", newline=""
        ) as stream:
            rows = list(csv.DictReader(stream))
        g4 = next(row for row in rows if row["sample_name_literal"] == "1_g4")
        self.assertEqual(g4["target_as_coded"], "g2")
        self.assertEqual(g4["included_in_model"], "true")
        self.assertFalse(designs["crispr_ets2"]["raw_count_matrix_present_in_same_archive_folder"])
        self.assertFalse(designs["oe_ets2"]["raw_count_matrix_present_in_same_archive_folder"])
        self.assertFalse(designs["meki"]["raw_count_matrix_present_in_same_archive_folder"])

    def test_derived_outputs_are_hash_recorded_and_repeatable(self) -> None:
        paths = {
            "gene_sets": ROOT / "data/derived/ets2-source-gene-sets.csv",
            "differential": ROOT / "data/derived/ets2-source-differential.csv.gz",
            "ranks": ROOT / "data/derived/ets2-source-ranks.csv.gz",
            "design": ROOT / "data/derived/ets2-source-design.csv",
            "notes": ROOT / "data/derived/ets2-source-notes.csv",
        }
        for name, path in paths.items():
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                self.report["derived_output_sha256"][name],
            )
        second = audit.build_audit(write_outputs=True)
        self.assertEqual(second["derived_output_sha256"], self.report["derived_output_sha256"])
        report_path = ROOT / "reports/ets2-source-audit.json"
        on_disk = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertEqual(on_disk["derived_output_sha256"], self.report["derived_output_sha256"])


if __name__ == "__main__":
    unittest.main()
