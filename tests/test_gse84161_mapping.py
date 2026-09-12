"""Focused structural tests for the GSE84161 identifier mapping audit."""

from __future__ import annotations

import csv
import gzip
import io
import json
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import map_gse84161 as mapping  # noqa: E402


class TestGSE84161Mapping(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = mapping.build_audit(write_outputs=True)

    def test_full_probe_and_crosswalk_universes_are_pinned(self) -> None:
        gpl = self.report["gpl570"]
        self.assertEqual(gpl["table_row_count"], 54675)
        self.assertEqual(gpl["table_unique_probe_count"], 54675)
        self.assertEqual(gpl["eligible_probe_rows"], 41834)
        self.assertEqual(gpl["eligible_unique_entrez_genes"], 20486)
        self.assertEqual(
            gpl["probe_status_counts"],
            {
                "ambiguous_multiple_entrez": 2300,
                "single_unique_entrez": 41834,
                "unmapped_blank_or_unresolved": 10541,
            },
        )

        crosswalk = self.report["ensembl116_crosswalk"]
        self.assertEqual(crosswalk["response_relationship_rows"], 91748)
        self.assertEqual(crosswalk["unique_ensembl_gene_ids"], 86411)
        self.assertEqual(crosswalk["count_response_value"], 86411)
        self.assertTrue(crosswalk["terminal_success_stamp"] == "[success]")
        self.assertEqual(crosswalk["ensembl_ids_without_entrez"], 54919)
        self.assertEqual(crosswalk["reciprocal_one_to_one_ensembl_genes"], 25210)
        self.assertEqual(crosswalk["duplicate_relationship_rows"], 0)

        array_universe = self.report["array_entrez_gene_universe"]
        self.assertEqual(array_universe["array_entrez_gene_count"], 20486)
        self.assertTrue(array_universe["all_genes_have_eligible_probe"])
        self.assertEqual(
            array_universe["status_counts"],
            {
                "ambiguous_ensembl_to_entrez": 176,
                "ambiguous_entrez_to_ensembl": 1258,
                "no_ensembl_relation": 1601,
                "reciprocal_one_to_one": 17451,
            },
        )

    def test_source_namespaces_denominators_and_stage_statuses(self) -> None:
        definitions = self.report["source_programs"]["definitions"]
        self.assertEqual(
            [item["source_namespace"] for item in definitions],
            ["Ensembl_gene_id"] * 4 + ["HGNC_symbol"] * 4,
        )
        self.assertEqual(
            [item["intended_member_count_after_ets2_exclusion"] for item in definitions],
            [927, 668, 875, 141, 815, 139, 436, 588],
        )

        coverage = {
            row["program_id"]: row for row in self.report["program_coverage"]
        }
        self.assertEqual(coverage["ets2_g1_dn"]["array_mapped_unique_entrez_gene_count"], 801)
        self.assertEqual(coverage["ets2_g1_dn"]["mapping_fraction"], "0.864077669903")
        self.assertEqual(coverage["ets2_g1_dn"]["coverage_gate_pass"], "true")
        self.assertEqual(coverage["ets2_g1_up"]["coverage_gate_pass"], "false")

        primary_crosswalk = json.loads(coverage["ets2_g1_dn"]["crosswalk_status_counts"])
        primary_array = json.loads(coverage["ets2_g1_dn"]["array_status_counts"])
        self.assertEqual(sum(primary_crosswalk.values()), 927)
        self.assertEqual(sum(primary_array.values()), 927)
        self.assertEqual(primary_crosswalk["ambiguous_ensembl_entrez_relation"], 74)
        self.assertEqual(primary_crosswalk["missing_ensembl_id_in_snapshot"], 2)
        self.assertEqual(primary_crosswalk["missing_entrez_relation"], 23)
        self.assertEqual(primary_array["mapped_crosswalk_no_eligible_array_probe"], 27)
        self.assertEqual(primary_array["not_array_mapped"], 99)

        derived = coverage[mapping.DERIVED_PROGRAM_ID]
        self.assertEqual(derived["parent_mapping_denominator"], 927)
        self.assertEqual(derived["parent_array_mapped_unique_entrez_gene_count"], 801)
        self.assertEqual(derived["parent_mapping_fraction"], "0.864077669903")
        self.assertEqual(derived["mapping_fraction"], "0.864077669903")
        self.assertEqual(derived["gate_evaluated_mapping_fraction"], "0.864077669903")
        self.assertEqual(derived["retained_mapped_unique_entrez_gene_count"], 652)
        self.assertEqual(derived["retained_mapping_fraction"], "0.70334412082")
        self.assertEqual(derived["mapped_comparator_overlap_unique_entrez_gene_count"], 149)
        self.assertEqual(derived["array_mapped_source_member_count"], 652)
        self.assertEqual(derived["coverage_gate_pass"], "true")

        derived_status = json.loads(derived["member_status_counts"])
        self.assertEqual(derived_status["mapped_comparator_overlap"], 149)
        self.assertEqual(derived_status["mapped_nonoverlap"], 652)

    def test_probe_entrez_tokens_are_deduplicated_and_stage_reasons_are_distinct(self) -> None:
        values, stats = mapping._probe_mapping(
            [
                {"ID": "probe_a", "ENTREZ_GENE_ID": "1 /// 1"},
                {"ID": "probe_b", "ENTREZ_GENE_ID": "1 /// 2"},
                {"ID": "probe_c", "ENTREZ_GENE_ID": "---"},
                {"ID": "probe_d", "ENTREZ_GENE_ID": "3 /// invalid"},
            ]
        )
        self.assertEqual(stats["status_counts"]["single_unique_entrez"], 1)
        self.assertEqual(stats["status_counts"]["ambiguous_multiple_entrez"], 1)
        self.assertEqual(stats["status_counts"]["unmapped_blank_or_unresolved"], 1)
        self.assertEqual(stats["status_counts"]["unmapped_invalid_entrez_token"], 1)
        self.assertEqual(values[0]["entrez_unique_tokens"], "1")
        self.assertEqual(values[0]["mapping_status"], "single_unique_entrez")

        maps = {
            "ensembl_to_entrez": {"ENSG00000000001": {"1"}},
            "entrez_to_ensembl": {"1": {"ENSG00000000001"}},
            "ensembl_universe": {
                "ENSG00000000001",
                "ENSG00000000003",
            },
            "symbol_to_ensembl": {
                "AMBIG": {"ENSG00000000001", "ENSG00000000002"},
                "ONE": {"ENSG00000000001"},
            },
            "reciprocal_ensembl_to_entrez": {"ENSG00000000001": "1"},
        }
        ambiguous = mapping._program_row_mapping(
            program_id="synthetic",
            source_member_index=1,
            source_row_number=2,
            source_value="AMBIG",
            source_namespace="HGNC_symbol",
            duplicate_source_member=False,
            crosswalk_maps=maps,
            probes_by_entrez={"1": ["probe_a"]},
        )
        self.assertEqual(ambiguous["crosswalk_mapping_status"], "ambiguous_hgnc_symbol")
        self.assertEqual(ambiguous["array_mapping_status"], "not_array_mapped")

        no_probe = mapping._program_row_mapping(
            program_id="synthetic",
            source_member_index=2,
            source_row_number=3,
            source_value="ONE",
            source_namespace="HGNC_symbol",
            duplicate_source_member=False,
            crosswalk_maps=maps,
            probes_by_entrez={},
        )
        self.assertEqual(no_probe["crosswalk_mapping_status"], "reciprocal_one_to_one")
        self.assertEqual(
            no_probe["array_mapping_status"],
            "mapped_crosswalk_no_eligible_array_probe",
        )

        missing_snapshot = mapping._program_row_mapping(
            program_id="synthetic",
            source_member_index=3,
            source_row_number=4,
            source_value="ENSG00000000002",
            source_namespace="Ensembl_gene_id",
            duplicate_source_member=False,
            crosswalk_maps=maps,
            probes_by_entrez={},
        )
        self.assertEqual(
            missing_snapshot["crosswalk_mapping_status"],
            "missing_ensembl_id_in_snapshot",
        )

        missing_entrez = mapping._program_row_mapping(
            program_id="synthetic",
            source_member_index=4,
            source_row_number=5,
            source_value="ENSG00000000003",
            source_namespace="Ensembl_gene_id",
            duplicate_source_member=False,
            crosswalk_maps=maps,
            probes_by_entrez={},
        )
        self.assertEqual(
            missing_entrez["crosswalk_mapping_status"],
            "missing_entrez_relation",
        )

    def test_outputs_have_stable_paths_widths_and_identifier_only_scope(self) -> None:
        expected = {
            "data/derived/gse84161-probe-map.csv.gz": (54675, 10, gzip.open),
            "data/derived/gse84161-gene-map.csv.gz": (91748, 16, gzip.open),
            "data/derived/gse84161-entrez-gene-universe.csv.gz": (20486, 8, gzip.open),
            "data/derived/gse84161-program-mapping.csv.gz": (5518, 22, gzip.open),
            "data/derived/gse84161-program-coverage.csv": (9, 29, open),
        }
        for relative, (row_count, width, opener) in expected.items():
            path = ROOT / relative
            with opener(path, "rt", encoding="utf-8", newline="") as stream:
                rows = list(csv.reader(stream))
            self.assertEqual(len(rows) - 1, row_count, relative)
            self.assertEqual({len(row) for row in rows}, {width}, relative)
            if relative.endswith("entrez-gene-universe.csv.gz"):
                entrez_ids = [row[0] for row in rows[1:]]
                self.assertEqual(entrez_ids, sorted(entrez_ids, key=lambda value: int(value)))
                self.assertEqual(len(entrez_ids), len(set(entrez_ids)))
                self.assertTrue(all(int(row[1]) > 0 for row in rows[1:]))
                self.assertEqual(
                    {row[6] for row in rows[1:]},
                    {
                        "reciprocal_one_to_one",
                        "ambiguous_entrez_to_ensembl",
                        "ambiguous_ensembl_to_entrez",
                        "no_ensembl_relation",
                    },
                )

        self.assertEqual(
            self.report["output_paths"],
            {
                "probe_map": "data/derived/gse84161-probe-map.csv.gz",
                "gene_map": "data/derived/gse84161-gene-map.csv.gz",
                "entrez_gene_universe": "data/derived/gse84161-entrez-gene-universe.csv.gz",
                "program_mapping": "data/derived/gse84161-program-mapping.csv.gz",
                "program_coverage": "data/derived/gse84161-program-coverage.csv",
                "report": "reports/gse84161-mapping-audit.json",
            },
        )
        self.assertFalse(self.report["scope"]["expression_values_read"])
        self.assertFalse(self.report["scope"]["array_intensities_read"])
        self.assertFalse(self.report["scope"]["program_scores_computed"])
        self.assertFalse(self.report["scope"]["target_values_or_differences_read"])
        self.assertFalse(self.report["scope"]["donor_or_treatment_effects_computed"])

    def test_cli_summary_is_json_and_reports_no_expression_read(self) -> None:
        # Exercise the command-line summary formatter without rebuilding the
        # large identifier tables a second time in this focused test suite.
        output = io.StringIO()
        original_build_audit = mapping.build_audit
        mapping.build_audit = lambda *, write_outputs=True: self.report
        try:
            with redirect_stdout(output):
                mapping.main()
        finally:
            mapping.build_audit = original_build_audit
        summary = json.loads(output.getvalue())
        self.assertEqual(summary["gpl570_probes"], 54675)
        self.assertEqual(summary["eligible_probe_rows"], 41834)
        self.assertEqual(summary["eligible_entrez_genes"], 20486)
        self.assertEqual(summary["crosswalk_relationship_rows"], 91748)
        self.assertEqual(summary["crosswalk_unique_ensembl"], 86411)
        self.assertEqual(summary["programs"], 9)
        self.assertFalse(summary["expression_values_read"])


if __name__ == "__main__":
    unittest.main()
