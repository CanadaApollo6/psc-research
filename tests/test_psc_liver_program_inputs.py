"""Identifier-only tests. Synthetic annotations, never expression values."""
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("prepare_psc_liver_program_inputs", ROOT/"scripts/prepare_psc_liver_program_inputs.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def ensg(number):
    return f"ENSG{number:011d}"


def relation(number, symbol, entrez=""):
    return {"ensembl_gene_id": ensg(number), "hgnc_symbol": symbol, "entrez_gene_id": entrez}


class IdentifierTests(unittest.TestCase):
    def test_case_sensitive_exact_symbol_only(self):
        reference = MODULE.build_reference([relation(1, "ABC", "1")])
        rows = MODULE.map_features(["ABC", "abc", "OLDABC", " ABC"], reference)
        self.assertEqual([row[MODULE.INTERFACES[0]] for row in rows], ["eligible", "no_exact_reference_symbol", "no_exact_reference_symbol", "no_exact_reference_symbol"])
        self.assertEqual([row["label_as_submitted"] for row in rows], ["ABC", "abc", "OLDABC", " ABC"])

    def test_full_reference_symbol_ambiguity_not_arbitrarily_resolved(self):
        reference = MODULE.build_reference([relation(1, "A", "1"), relation(2, "A", "1")])
        row = MODULE.map_features(["A"], reference)[0]
        self.assertEqual(row[MODULE.INTERFACES[0]], "symbol_has_multiple_ensembl")
        self.assertEqual(row["candidate_ensembl_ids"], [ensg(1), ensg(2)])
        normalized, candidates, status = MODULE.resolve_source_member(ensg(1), "Ensembl_gene_id", reference)
        self.assertEqual(status, "symbol_has_multiple_ensembl")

    def test_entrez_ambiguity_uses_genes_absent_from_input_features(self):
        reference = MODULE.build_reference([relation(1, "A", "7"), relation(2, "B", "7")])
        row = MODULE.map_features(["A"], reference)[0]
        self.assertEqual(row[MODULE.INTERFACES[0]], "eligible")
        self.assertEqual(row[MODULE.INTERFACES[1]], "entrez_has_multiple_ensembl")
        self.assertEqual(len(reference["entrez_ens"]["7"]), 2)

    def test_multiple_entrez_relations_remain_distinct(self):
        reference = MODULE.build_reference([relation(1, "A", "1"), relation(1, "A", "2")])
        row = MODULE.map_features(["A"], reference)[0]
        self.assertEqual(row[MODULE.INTERFACES[0]], "eligible")
        self.assertEqual(row[MODULE.INTERFACES[1]], "multiple_entrez_for_ensembl")
        self.assertEqual(row["candidate_entrez_ids"], ["1", "2"])

    def test_missing_entrez_does_not_mean_missing_measured_gene(self):
        reference = MODULE.build_reference([relation(1, "NONCODING")])
        row = MODULE.map_features(["NONCODING"], reference)[0]
        self.assertEqual(row[MODULE.INTERFACES[0]], "eligible")
        self.assertEqual(row[MODULE.INTERFACES[1]], "no_entrez_relation")
        self.assertEqual(row["source_row_index_0based"], 0)

    def test_reverse_symbol_and_source_feature_collision_not_collapsed(self):
        reference = MODULE.build_reference([relation(1, "A", "1"), relation(1, "B", "1")])
        rows = MODULE.map_features(["A", "B"], reference)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["source_feature_collision_for_ensembl"] for row in rows))
        self.assertTrue(all(row[MODULE.INTERFACES[0]] == "ensembl_has_multiple_symbols" for row in rows))

    def test_duplicate_source_features_rejected(self):
        reference = MODULE.build_reference([relation(1, "A", "1")])
        with self.assertRaises(ValueError):
            MODULE.map_features(["A", "A"], reference)

    def test_duplicate_reference_relationships_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.build_reference([relation(1, "A", "1"), relation(1, "A", "1")])

    def test_non_numeric_entrez_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.build_reference([relation(1, "A", "guess")])

    def test_numeric_version_removal_only(self):
        self.assertEqual(MODULE.stable_ensembl(ensg(1)+".12"), ensg(1))
        for bad in [ensg(1)+"_PAR_Y", ensg(1)+".x", " "+ensg(1), ensg(1).lower()]:
            with self.assertRaises(ValueError):
                MODULE.stable_ensembl(bad)

    def test_absent_and_present_without_symbol_differ(self):
        reference = MODULE.build_reference([relation(1, "", "1")])
        self.assertEqual(MODULE.resolve_source_member(ensg(1), "Ensembl_gene_id", reference)[2], "ensembl_present_without_hgnc_symbol")
        self.assertEqual(MODULE.resolve_source_member(ensg(2), "Ensembl_gene_id", reference)[2], "ensembl_id_absent_from_reference_snapshot")

    def test_ETS2_exclusion_preserved_for_both_namespaces(self):
        reference = MODULE.build_reference([relation(1, "A", "1")])
        for value, namespace in [(MODULE.ETS2, "Ensembl_gene_id"), ("ETS2", "HGNC_symbol")]:
            self.assertEqual(MODULE.resolve_source_member(value, namespace, reference)[2], "excluded_ETS2")

    def fixture_programs(self):
        programs = {"ets2_g1_dn": {"source_namespace": "Ensembl_gene_id", "source_values": [ensg(1), ensg(2), ensg(99), MODULE.ETS2], "intended_members": 3}}
        for program in MODULE.COMPARATORS:
            programs[program] = {"source_namespace": "HGNC_symbol", "source_values": ["A"], "intended_members": 1}
        return programs

    def test_missing_member_keeps_none_not_zero_and_index_zero_maps(self):
        reference = MODULE.build_reference([relation(1, "A", "1"), relation(2, "B", "2")])
        features = MODULE.map_features(["A", "B"], reference)
        rows, inventory, _ = MODULE.map_program_members(self.fixture_programs(), {"minimum_mapped_fraction": .5, "minimum_mapped_genes": 1}, reference, features)
        primary = [row for row in rows if row["program_id"] == "ets2_g1_dn"]
        self.assertEqual(primary[0][MODULE.INTERFACES[0]]["source_feature_index_0based"], 0)
        self.assertIsNone(primary[2][MODULE.INTERFACES[0]]["source_feature_index_0based"])
        coverage = next(row for row in inventory if row["program_id"] == "ets2_g1_dn" and row["interface"] == MODULE.INTERFACES[0])
        self.assertEqual(coverage["intended_source_members_after_ETS2_exclusion"], 3)
        self.assertEqual(coverage["mapped_unique_ensembl_ids"], 2)
        self.assertEqual(coverage["not_mapped_intended_members"], 1)

    def test_nonoverlap_gate_keeps_parent_denominator(self):
        reference = MODULE.build_reference([relation(1, "A", "1"), relation(2, "B", "2")])
        features = MODULE.map_features(["A", "B"], reference)
        _, inventory, sets = MODULE.map_program_members(self.fixture_programs(), {"minimum_mapped_fraction": .5, "minimum_mapped_genes": 1}, reference, features)
        sensitivity = next(row for row in inventory if row["program_id"] == "ets2_g1_dn_without_comparators" and row["interface"] == MODULE.INTERFACES[0])
        self.assertEqual(sensitivity["overlap_removed"], 1)
        self.assertEqual(sensitivity["mapped_unique_ensembl_ids"], 1)
        self.assertAlmostEqual(sensitivity["retained_fraction_of_original_intended_parent"], 1/3)
        self.assertAlmostEqual(sensitivity["parent_mapping_fraction_used_for_gate"], 2/3)
        self.assertTrue(sensitivity["mapping_gate_pass"])

    def test_background_excludes_all_original_mapped_programs_and_target(self):
        reference = MODULE.build_reference([relation(1, "A", "1"), relation(2, "B", "2"), {"ensembl_gene_id": MODULE.ETS2, "hgnc_symbol": "ETS2", "entrez_gene_id": "2114"}])
        features = MODULE.map_features(["A", "B", "ETS2", "UNKNOWN"], reference)
        sets = {interface: {"ets2_g1_dn": {ensg(1)}, "ets2_g1_dn_without_comparators": {ensg(1)}} for interface in MODULE.INTERFACES}
        pools, rows = MODULE.candidate_pools(features, sets)
        for interface in MODULE.INTERFACES:
            self.assertEqual(pools[interface]["pool_status_counts"], {"excluded_original_program_union": 1, "eligible_background_candidate": 1, "excluded_ETS2": 1, "identity_not_qualified": 1})
            self.assertEqual(pools[interface]["released_feature_axis_rows"], 4)
            self.assertFalse(pools[interface]["abundance_bins_or_draws_created"])
            self.assertTrue(pools[interface]["counts_normalization_denominator_not_filtered_by_this_pool"])
        self.assertEqual(len(rows), 8)

    def test_duplicate_program_members_cannot_collapse_silently(self):
        reference = MODULE.build_reference([relation(1, "A", "1"), relation(2, "B", "2")])
        features = MODULE.map_features(["A", "B"], reference)
        programs = self.fixture_programs()
        programs["ets2_g1_dn"]["source_values"] = [ensg(1), ensg(1)+".1"]
        programs["ets2_g1_dn"]["intended_members"] = 2
        with self.assertRaises(ValueError):
            MODULE.map_program_members(programs, {"minimum_mapped_fraction": .5, "minimum_mapped_genes": 1}, reference, features)

    def test_legacy_mapped_status_or_effect_fields_do_not_select_members(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"members.csv"
            path.write_text("source_value,mapping_status,legacy_effect\nA,mapped,do_not_use\nB,unmapped,do_not_use\n")
            selected = MODULE.selected_identifier_rows(path, ["source_value"])
            self.assertEqual(selected, [{"source_value": "A"}, {"source_value": "B"}])

    def test_module_has_no_effect_or_network_actions(self):
        text = (ROOT/"scripts/prepare_psc_liver_program_inputs.py").read_text()
        for option in ["--score", "--match", "--normalize", "--effect", "--model", "--acquire"]:
            self.assertNotIn('add_argument("'+option, text)
        self.assertNotIn("import requests", text)
        self.assertNotIn("import numpy", text)


    def test_failed_program_gate_still_excludes_mapped_gene_from_background(self):
        reference = MODULE.build_reference([relation(index, symbol, str(index)) for index, symbol in enumerate(["A", "B", "C", "D"], 1)])
        features = MODULE.map_features(["A", "B", "C", "D"], reference)
        programs = self.fixture_programs()
        programs["ets2_g1_dn"]["source_values"] = [ensg(1), ensg(2), MODULE.ETS2]
        programs["ets2_g1_dn"]["intended_members"] = 2
        programs["ets2_g1_up"] = {"source_namespace": "Ensembl_gene_id", "source_values": [ensg(3), ensg(99)], "intended_members": 2}
        for program, gene in [("ets2_g2_dn", 1), ("chr21_dn", 2)]:
            programs[program] = {"source_namespace": "Ensembl_gene_id", "source_values": [ensg(gene)], "intended_members": 1}
        self.assertEqual(len(programs), 8)
        _, inventory, sets = MODULE.map_program_members(programs, {"minimum_mapped_fraction": .8, "minimum_mapped_genes": 1}, reference, features)
        pools, rows = MODULE.candidate_pools(features, sets)
        for interface in MODULE.INTERFACES:
            companion = next(row for row in inventory if row["program_id"] == "ets2_g1_up" and row["interface"] == interface)
            self.assertFalse(companion["mapping_gate_pass"])
            self.assertEqual(pools[interface]["eight_original_mapped_program_union_genes"], 3)
            self.assertEqual(pools[interface]["background_candidate_genes"], 1)
            self.assertEqual(next(row for row in rows if row["label_as_submitted"] == "C" and row["interface"] == interface)["pool_status"], "excluded_original_program_union")

    def test_pin_change_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"input.txt"
            path.write_text("original")
            pin = MODULE.pin_record(path, Path(directory))
            MODULE.check_pin(Path(directory), pin)
            path.write_text("modified")
            with self.assertRaises(ValueError):
                MODULE.check_pin(Path(directory), pin)


class PreparedArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT/"data/derived/psc-liver-program-inputs.json").read_text())
        cls.summary = json.loads((ROOT/"work/psc-liver-program-inputs/preparation-summary.json").read_text())

    def test_public_read_only_verifier(self):
        result = MODULE.verify_preparation()
        self.assertEqual(result["status"], "all_preparation_and_closed_material_pins_verified")
        self.assertEqual(result["ignored_preparation_artifacts_checked"], 8)
        self.assertFalse(result["expression_analysis_performed"])

    def test_schema_and_no_inference_scope(self):
        self.assertEqual(self.data["schema_version"], 1)
        for flag in ["normalization_performed", "expression_matching_performed", "scores_or_effects_computed", "aliases_guessed", "ambiguous_IDs_collapsed", "missing_genes_zero_filled", "atlas_mapped_program_subsets_reused", "organoid_derived_genes_used", "historical_counting_GTF_reconstructed"]:
            self.assertFalse(self.data["scope"][flag])
        self.assertTrue(self.data["scope"]["root_method_and_analysis_freeze_required"])
        self.assertEqual(self.data["scope"]["network_requests"], 0)

    def test_complete_feature_axis_and_both_interfaces(self):
        self.assertEqual(self.data["source_feature_axis"]["released_source_rows"], 41096)
        self.assertEqual(self.data["source_feature_axis"]["literal_Ensembl_ID_rows"], 0)
        self.assertTrue(self.data["source_feature_axis"]["all_original_source_rows_preserved"])
        self.assertTrue(self.data["source_feature_axis"]["ID_eligibility_not_a_count_normalization_universe"])
        expected_eligible = [37778, 22562]
        for interface, eligible in zip(MODULE.INTERFACES, expected_eligible):
            counts = self.data["feature_identifier_status_counts"][interface]
            self.assertEqual(sum(counts.values()), 41096)
            self.assertEqual(counts["eligible"], eligible)

    def test_original_denominators_and_full_export_provenance(self):
        sources = self.data["source_provenance"]
        self.assertEqual([row["intended_source_denominator"] for row in sources["fixed_source_programs"]], [927, 668, 875, 141, 815, 139, 436, 588])
        self.assertEqual(sources["complete_original_source_memberships_including_ETS2"], 4590)
        self.assertEqual(sources["complete_intended_memberships_after_ETS2_exclusion"], 4589)
        self.assertFalse(sources["original_ETS2_archive"]["reread"])
        self.assertEqual(self.data["whole_human_reference"]["relationship_rows"], 91748)
        self.assertEqual(self.data["whole_human_reference"]["unique_ensembl_ids"], 86411)

    def test_exact_reconstruction_is_not_fresh_retrieval(self):
        reconstruction = self.data["source_provenance"]["reference_raw_response_reconstruction"]
        self.assertTrue(reconstruction["exact_pin_match"])
        self.assertFalse(reconstruction["new_network_retrieval"])
        self.assertFalse(reconstruction["direct_historical_raw_cache_reread"])
        self.assertFalse(reconstruction["original_raw_path_present"])
        self.assertEqual(reconstruction["candidate_bytes"], 2236970)
        self.assertEqual(reconstruction["candidate_sha256"], "a60cf0e0ac5782e7c39269a577d1ee24be08d26c100324fe2c38b77a7f81f167")

    def test_program_coverage_and_nonoverlap_keep_original_parent_gate(self):
        coverage = {(row["interface"], row["program_id"]): row for row in self.data["program_coverage"]}
        for interface, mapped, retained, removed in [(MODULE.INTERFACES[0], 848, 695, 153), (MODULE.INTERFACES[1], 825, 673, 152)]:
            primary = coverage[(interface, "ets2_g1_dn")]
            sensitivity = coverage[(interface, "ets2_g1_dn_without_comparators")]
            self.assertEqual(primary["mapped_unique_ensembl_ids"], mapped)
            self.assertEqual(primary["intended_source_members_after_ETS2_exclusion"], 927)
            self.assertEqual(sensitivity["mapped_unique_ensembl_ids"], retained)
            self.assertEqual(sensitivity["overlap_removed"], removed)
            self.assertEqual(sensitivity["parent_mapping_fraction_used_for_gate"], mapped/927)
            self.assertTrue(sensitivity["mapping_gate_pass"])
            self.assertLess(sensitivity["retained_fraction_of_original_intended_parent"], .8)
            self.assertFalse(sensitivity["effects_authorized"])

    def test_pool_sizes_and_no_matching_or_count_filter(self):
        for interface, eligible, union, background in [(MODULE.INTERFACES[0], 37778, 2944, 34833), (MODULE.INTERFACES[1], 22562, 2831, 19730)]:
            pool = self.data["candidate_matching_reference_pools"][interface]
            self.assertEqual(pool["ID_qualified_feature_rows"], eligible)
            self.assertEqual(pool["eight_original_mapped_program_union_genes"], union)
            self.assertEqual(pool["background_candidate_genes"], background)
            self.assertEqual(sum(pool["pool_status_counts"].values()), 41096)
            self.assertFalse(pool["matching_pool_sufficiency_observed"])
            self.assertFalse(pool["abundance_bins_or_draws_created"])
            self.assertTrue(pool["counts_normalization_denominator_not_filtered_by_this_pool"])

    def test_missing_indices_are_none_not_synthetic_zero(self):
        members = json.loads((ROOT/"work/psc-liver-program-inputs/program-member-identifier-audit.json").read_text())
        missing = []
        for member in members:
            for interface in MODULE.INTERFACES:
                resolution = member[interface]
                if resolution["status"] != "mapped":
                    missing.append(resolution)
                    self.assertIsNone(resolution["source_feature_index_0based"])
                    self.assertIsNone(resolution["ensembl_gene_id"])
        self.assertEqual(len(members), 4590)
        self.assertEqual(len(missing), 891)

    def test_clinical_and_historical_annotation_limits_retained(self):
        limits = self.data["inference_boundaries"]
        self.assertEqual({key: limits["sample_design_retained_from_closed_qualification"][key] for key in ["PSC", "ASC", "AIH"]}, {"PSC": 17, "ASC": 17, "AIH": 30})
        self.assertTrue(limits["current_reference_not_historical_GTF_or_biological_locus_reconstruction"])
        self.assertIn("six paired PSC patients", limits["GSE159676"])
        self.assertIn("PARK", limits["NoPSC_atlas"])
        self.assertTrue(limits["root_identity_policy_choice_required_before_matching_or_effects"])

    def test_source_and_ignored_artifact_pins(self):
        MODULE.verify_sources()
        for pin in self.data["ignored_artifact_pins"]:
            MODULE.check_pin(ROOT, pin)
        preservation = MODULE.verify_closed_materials()
        self.assertEqual(preservation["closed_liver_versioned_files"], 4)
        self.assertEqual(preservation["closed_organoid_versioned_files"], 4)
        self.assertEqual(preservation["closed_organoid_source_and_artifact_pins"], 67)
        self.assertTrue(preservation["all_checked_pins_unchanged"])


if __name__ == "__main__":
    unittest.main()
