"""Focused structural tests for the GSE84161 metadata audit."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import audit_gse84161 as audit  # noqa: E402


class TestGSE84161Metadata(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = audit.build_audit(write_outputs=True)
        cls.monocytes = cls.report["gse84161"]

    def test_changed_source_is_rejected_before_parsing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            changed = Path(directory) / "changed.soft"
            changed.write_text("^SERIES = GSE84161\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Source bytes changed for GSE84161"):
                audit.build_audit(write_outputs=False, gse84161_path=changed)

    def test_exact_six_arm_five_group_accounting_and_controls(self) -> None:
        groups = self.monocytes["replicate_groups"]
        self.assertEqual(self.monocytes["sample_count"], 30)
        self.assertEqual(self.monocytes["explicit_series_donor_count"], 5)
        self.assertEqual(
            [group["replicate_number"] for group in groups],
            [1, 2, 3, 4, 5],
        )
        self.assertEqual([group["arm_count"] for group in groups], [6] * 5)
        self.assertEqual(self.monocytes["planned_control_pair_count"], 20)
        self.assertEqual(self.monocytes["mapping_errors"], [])
        self.assertEqual(
            [arm["condition_key"] for arm in groups[0]["arms"]],
            audit.ARM_ORDER,
        )
        pairs = self.monocytes["planned_control_pairs"]
        by_treated_accession = {
            pair["treated_sample_accession"]: pair for pair in pairs
        }
        self.assertEqual(
            by_treated_accession["GSM2228387"]["control_sample_accession"],
            "GSM2228386",
        )
        self.assertEqual(
            by_treated_accession["GSM2228390"]["control_sample_accession"],
            "GSM2228389",
        )

    def test_nonsequential_superseries_block_and_duplicate_label_guard(self) -> None:
        comparison = self.report["series_comparison"]
        self.assertEqual(comparison["common_monocyte_sample_count"], 30)
        self.assertTrue(comparison["common_sample_ids_unique"])
        self.assertTrue(comparison["monocyte_metadata_exact_match"])
        self.assertEqual(comparison["superseries_numeric_gsm_step_values"], [1, 104])
        self.assertTrue(comparison["superseries_contains_nonsequential_gsm_block"])
        self.assertEqual(
            comparison["superseries_only_neutrophil_context"]["sample_count"],
            20,
        )

        parsed = audit.parse_soft(audit.GSE84161_SOFT)
        metadata = [audit._sample_metadata(sample) for sample in parsed["samples"]]
        duplicate = deepcopy(metadata[0])
        duplicate["sample_accession"] = "synthetic_duplicate_label"
        _, _, errors = audit._build_control_map(metadata + [duplicate])
        self.assertTrue(any("duplicate arm" in error for error in errors))

    def test_processed_subset_and_full_gpl_mapping_are_distinguished(self) -> None:
        processed = self.monocytes["processed_table"]
        platform = self.monocytes["platform_gpl570"]
        processed_annotation = self.monocytes["processed_annotation_via_current_gpl570"]
        full_annotation = self.monocytes["full_gpl570_annotation"]

        self.assertEqual(processed["processed_table_header"], ["ID_REF", "VALUE"])
        self.assertEqual(processed["idref_namespace"], "Affymetrix probe-set ID_REF")
        self.assertEqual(processed["value_label"], "RMA normalized signal intensity")
        self.assertEqual(processed["processed_row_count"], 19944)
        self.assertEqual(processed["processed_unique_id_count"], 19944)
        self.assertEqual(processed["processed_duplicate_id_count"], 0)
        self.assertTrue(processed["all_samples_same_ordered_ids"])
        self.assertEqual(platform["table_row_count"], 54675)
        self.assertEqual(platform["table_unique_id_count"], 54675)
        self.assertEqual(platform["table_duplicate_id_count"], 0)
        self.assertEqual(platform["table_last_ids"], ["AFFX-TrpnX-3_at", "AFFX-TrpnX-5_at", "AFFX-TrpnX-M_at"])
        self.assertIn("ENTREZ_GENE_ID", platform["annotation_columns"])
        self.assertIn("RefSeq Transcript ID", platform["annotation_columns"])
        self.assertEqual(processed_annotation["probe_rows_considered"], 19944)
        self.assertEqual(processed_annotation["blank_or_unresolved_entrez_fields"], 267)
        self.assertEqual(processed_annotation["single_entrez_fields"], 19263)
        self.assertEqual(processed_annotation["multiple_entrez_fields"], 414)
        self.assertEqual(full_annotation["probe_rows_considered"], 54675)
        self.assertEqual(full_annotation["blank_or_unresolved_entrez_fields"], 10541)
        self.assertEqual(full_annotation["single_entrez_fields"], 41834)
        self.assertEqual(full_annotation["multiple_entrez_fields"], 2300)
        self.assertEqual(full_annotation["unique_entrez_assignments"], 22600)
        self.assertEqual(full_annotation["unique_genes_from_single_id_probes"], 20486)
        self.assertEqual(processed_annotation["unique_genes_from_single_id_probes"], 19228)

    def test_entrez_tokens_are_deduplicated_within_one_probe(self) -> None:
        platform = {
            "annotation_by_id": {
                "probe_a": {"ENTREZ_GENE_ID": "1 /// 1 /// 2"},
                "probe_b": {"ENTREZ_GENE_ID": "2"},
            }
        }
        stats = audit._annotation_stats(["probe_a", "probe_b"], platform)
        self.assertEqual(stats["single_entrez_fields"], 1)
        self.assertEqual(stats["multiple_entrez_fields"], 1)
        self.assertEqual(stats["total_entrez_assignments"], 3)
        self.assertEqual(stats["unique_entrez_assignments"], 2)
        self.assertEqual(stats["entrez_ids_assigned_to_more_than_one_probe"], 1)
        self.assertEqual(stats["probe_assignments_in_multi_probe_entrez_groups"], 2)

    def test_missing_control_is_reported_without_filling_metadata(self) -> None:
        parsed = audit.parse_soft(audit.GSE84161_SOFT)
        metadata = [audit._sample_metadata(sample) for sample in parsed["samples"]]
        reduced = [
            item
            for item in metadata
            if item["sample_accession"] != "GSM2228386"
        ]
        groups, pairs, errors = audit._build_control_map(reduced)
        self.assertEqual(groups[0]["arm_count"], 5)
        self.assertEqual(len(pairs), 18)
        self.assertTrue(any("missing control untreated_unstim" in error for error in errors))
        self.assertTrue(all(item["status"] == "planned_unscored" for item in pairs))

    def test_processed_identifier_metadata_is_invariant_to_value_text(self) -> None:
        template = """^SERIES = TEST\n^SAMPLE = GSMTEST\n!Sample_data_row_count = 2\n#VALUE = RMA normalized signal intensity\n!sample_table_begin\nID_REF\tVALUE\nprobe_a\t{first}\nprobe_b\t{second}\n!sample_table_end\n"""
        with tempfile.TemporaryDirectory(dir=audit.OUT_WORK) as directory:
            first_path = Path(directory) / "first.soft"
            second_path = Path(directory) / "second.soft"
            first_path.write_text(template.format(first="1.0", second="2.0"), encoding="utf-8")
            second_path.write_text(template.format(first="9001.25", second="-77"), encoding="utf-8")
            first = audit.parse_soft(first_path)["samples"]
            second = audit.parse_soft(second_path)["samples"]
        first_table = first[0]["table"]
        second_table = second[0]["table"]
        self.assertEqual(first_table["header"], second_table["header"])
        self.assertEqual(first_table["row_count"], second_table["row_count"])
        self.assertEqual(first_table["feature_ids"], second_table["feature_ids"])
        self.assertEqual(first_table["id_digest"], second_table["id_digest"])
        self.assertEqual(
            audit._processed_consistency(first),
            audit._processed_consistency(second),
        )

    def test_truncated_soft_table_fails_instead_of_being_accepted(self) -> None:
        with tempfile.TemporaryDirectory(dir=audit.OUT_WORK) as directory:
            path = Path(directory) / "truncated.soft"
            path.write_text(
                "^SAMPLE = GSMTEST\n!sample_table_begin\nID_REF\tVALUE\nprobe_a\t1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "unterminated sample table"):
                audit.parse_soft(path)

    def test_source_contradictions_and_donor_ambiguity_are_retained(self) -> None:
        ambiguity_ids = {
            item["id"] for item in self.report["contradictions_and_ambiguities"]
        }
        self.assertIn("pretreatment_duration", ambiguity_ids)
        self.assertIn("donor_vs_replicate", ambiguity_ids)
        self.assertIn("variance_scope", ambiguity_ids)
        self.assertIn("heatmap_residual_quantity", ambiguity_ids)
        self.assertIn("historical_gpl_annotation", ambiguity_ids)
        timing = next(
            item
            for item in self.report["contradictions_and_ambiguities"]
            if item["id"] == "pretreatment_duration"
        )
        self.assertEqual(timing["kind"], "differing_specificity")
        self.assertEqual(
            self.monocytes["timing_evidence"],
            {
                "series_design_pretreatment": "1 hr",
                "sample_protocol_pretreatment": "at least 30 min",
            },
        )
        self.assertFalse(self.monocytes["sample_level_donor_field_present"])
        self.assertEqual(
            self.monocytes["donor_pairing_status"],
            "series-level five-donor statement; sample-level replicate number 1-5 only; donor identity not auditable from these records",
        )

    def test_neutrophil_inventory_is_separate_from_monocyte_evidence(self) -> None:
        neutrophils = self.report["gse84162"]["superseries_only_neutrophil_context"]
        self.assertEqual(neutrophils["accession"], "GSE84153")
        self.assertEqual(neutrophils["cell_context"], "purified human neutrophils")
        self.assertEqual(neutrophils["sample_count"], 20)
        self.assertEqual(neutrophils["explicit_donor_tags"], [
            "Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5",
        ])
        self.assertEqual(set(neutrophils["samples_per_donor"].values()), {4})
        self.assertFalse(neutrophils["substitution_for_monocytes"])
        self.assertIn("Entrez Gene IDs", neutrophils["processed_file_format"])

    def test_outputs_are_complete_csvs_and_metadata_only(self) -> None:
        sample_path = ROOT / "data/derived/gse84161-samples.csv"
        processing_path = ROOT / "data/derived/gse84161-processing-evidence.csv"
        report_path = ROOT / "reports/gse84161-metadata-audit.json"
        with sample_path.open(encoding="utf-8", newline="") as stream:
            samples = list(csv.DictReader(stream))
        with processing_path.open(encoding="utf-8", newline="") as stream:
            evidence = list(csv.DictReader(stream))
        self.assertEqual(len(samples), 30)
        self.assertEqual(len({row["sample_accession"] for row in samples}), 30)
        self.assertTrue(all(row["cel_available_route"] == "true" for row in samples))
        self.assertEqual(len(evidence), 11)
        self.assertTrue(all(set(row) == set(evidence[0]) for row in evidence))
        on_disk = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertFalse(on_disk["scope"]["expression_values_read"])
        self.assertFalse(on_disk["scope"]["program_scores_computed"])
        self.assertFalse(on_disk["scope"]["target_values_or_differences_read"])


if __name__ == "__main__":
    unittest.main()
