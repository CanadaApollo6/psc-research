"""Focused tests for the bounded MEK-inhibitor reconstruction."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import reconstruct_ets2_inhibitor as reconstruction  # noqa: E402


class TestETS2InhibitorReconstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = reconstruction.build_reconstruction(write_outputs=True)

    def test_rank_vectors_are_complete_finite_and_strictly_ordered(self) -> None:
        summaries = {
            row["rank_dataset_id"]: row
            for row in self.report["rank_vectors"]["summaries"]
        }
        self.assertEqual(set(summaries), {"meki_100", "meki_500"})
        for row in summaries.values():
            self.assertEqual(row["rank_row_count"], 12095)
            self.assertEqual(row["finite_statistic_count"], 12095)
            self.assertEqual(row["unique_gene_symbol_count"], 12095)
            self.assertEqual(row["strict_descending_order"], "true")
            self.assertEqual(row["adjacent_tie_count"], 0)
            self.assertGreater(row["positive_count"], 0)
            self.assertGreater(row["negative_count"], 0)

    def test_cross_rank_identity_is_reported_without_treating_it_as_validation(self) -> None:
        overlap = self.report["rank_vectors"]["identity_overlaps"]
        meki_self = next(
            row for row in overlap
            if row["left_dataset_id"] == "meki_100"
            and row["right_dataset_id"] == "meki_500"
        )
        self.assertEqual(meki_self["shared_symbol_count"], 12095)
        meki_ko = next(
            row for row in overlap
            if row["left_dataset_id"] == "meki_100"
            and row["right_dataset_id"] == "ko_g1"
        )
        self.assertEqual(meki_ko["shared_symbol_count"], 11041)

    def test_signed_es_matches_archived_output_and_fig5_nes(self) -> None:
        compatibility = self.report["fgsea_signed_es_compatibility"]
        self.assertTrue(compatibility["all_nine_es_within_1e-9"])
        self.assertTrue(compatibility["all_nine_fig5_and_archive_nes_within_1e-9"])
        self.assertTrue(compatibility["all_nine_signed_es_agree_with_reported_fig5_nes"])
        self.assertTrue(compatibility["all_nine_sign_flipped_rank_es_signs_invert"])
        self.assertEqual(compatibility["max_abs_computed_es_minus_archive_es"], "0.000000000000163")
        self.assertEqual(compatibility["max_abs_fig5_minus_archive_nes"], "0.0000000000000002")

    def test_code_reconstruction_contains_explicit_steps_and_unresolved_export(self) -> None:
        steps = {
            row["step_id"]: row for row in self.report["reconstruction_steps"]
        }
        expected = {
            "input_count_matrix", "sample_names", "donor_and_drug_annotation",
            "expression_filter", "library_normalization", "donor_aware_design",
            "contrasts", "voom", "lm_fit", "contrasts_fit", "e_bayes",
            "top_table_statistic", "symbol_mapping", "symbol_collapse",
            "rank_serialization", "fgsea_parameters",
        }
        self.assertEqual(set(steps), expected)
        self.assertIn("topTable", steps["top_table_statistic"]["evidence_text_original"])
        self.assertEqual(steps["symbol_mapping"]["confidence"], "low_to_medium")
        self.assertIn("cannot be certified", steps["rank_serialization"]["limitation"])

    def test_source_fig5_rows_and_output_hashes_are_preserved(self) -> None:
        fig5 = self.report["source_data_fig5"]
        self.assertEqual(fig5["sheet"], "5j")
        self.assertEqual(fig5["worksheet_dimension"], "A1:C19")
        self.assertEqual(len(fig5["notes"]), 5)
        self.assertIn("treated versus vehicle", fig5["notes"][2]["value_original"])

        report_path = ROOT / "reports/ets2-inhibitor-reconstruction.json"
        on_disk = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertEqual(on_disk["derived_output_sha256"], self.report["derived_output_sha256"])
        for name, digest in self.report["derived_output_sha256"].items():
            path = ROOT / "data/derived" / f"ets2-inhibitor-reconstruction-{name}.csv"
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)

    def test_pathway_table_keeps_all_nine_original_sets(self) -> None:
        path = ROOT / "data/derived/ets2-inhibitor-reconstruction-pathways.csv"
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 9)
        self.assertTrue(all(row["computed_es_within_1e-9_of_archive_es"] == "true" for row in rows))
        self.assertTrue(all(row["computed_es_and_reported_nes_sign_agree"] == "true" for row in rows))
        self.assertEqual(
            {row["pathway_name_reported_original"] for row in rows},
            {
                "Macrophage activation", "Myeloid leukocyte activation",
                "Myeloid leukocyte migration", "Phagocytosis", "IL-1 production",
                "IL-6 production", "IL-8 production", "ROS production",
                "TNFSF cytokine production",
            },
        )


if __name__ == "__main__":
    unittest.main()
