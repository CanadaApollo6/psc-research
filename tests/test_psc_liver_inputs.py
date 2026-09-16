"""Synthetic qualification tests. No disease contrasts or target/program values."""
import csv
import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("qualify_psc_liver_inputs", ROOT / "scripts/qualify_psc_liver_inputs.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class LiverParserTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "input.txt"

    def write(self, text):
        self.path.write_text(text)
        return self.path

    def series(self, rows, columns=("GSM1", "GSM2"), end=True):
        header = "\t".join('"' + value + '"' for value in columns)
        text = '!Sample_geo_accession\t' + header + '\n!series_matrix_table_begin\n"ID_REF"\t' + header + "\n"
        text += "\n".join(rows) + "\n"
        if end:
            text += "!series_matrix_table_end\n"
        return self.write(text)

    def test_geo_keeps_original_diagnosis_and_skips_numeric_preview(self):
        path = self.write("^SAMPLE = GSM1\n!Sample_title = A\n!Sample_characteristics_ch1 = sample group: ASC\n#VALUE = Quantile normalized\n!sample_table_begin\nID_REF\tVALUE\np1\t999999\n!sample_table_end\n")
        record = MODULE.geo_records(path)[0]
        self.assertEqual(MODULE.characteristics(record), {"sample group": "ASC"})
        self.assertEqual(record["table"], [])
        self.assertEqual(record["column_descriptions"], {"VALUE": "Quantile normalized"})
        self.assertEqual(MODULE.geo_records(path, include_tables=True)[0]["table"][0]["VALUE"], "999999")

    def test_repeated_characteristic_not_silently_overwritten(self):
        record = {"fields": {"Sample_characteristics_ch1": ["x: a", "x: b"]}}
        with self.assertRaises(ValueError):
            MODULE.characteristics(record)

    def test_unclosed_soft_table_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.geo_records(self.write("^SAMPLE = GSM1\n!sample_table_begin\nID_REF\tVALUE\n"))

    def test_exact_join_reorders_but_never_guesses(self):
        samples = [{"accession": "GSM1", "fields": {"Sample_title": ["A"]}}, {"accession": "GSM2", "fields": {"Sample_title": ["B"]}}]
        self.assertEqual(MODULE.exact_sample_join(["B", "A"], samples, "Sample_title")[0]["accession"], "GSM2")
        with self.assertRaises(ValueError):
            MODULE.exact_sample_join(["a", "B"], samples, "Sample_title")
        with self.assertRaises(ValueError):
            MODULE.exact_sample_join(["A", "A"], samples, "Sample_title")
        with self.assertRaises(ValueError):
            MODULE.exact_sample_join(["A"], samples, "Sample_title")

    def test_duplicate_metadata_keys_rejected(self):
        samples = [{"accession": "GSM1", "fields": {"Sample_title": ["A"]}}, {"accession": "GSM2", "fields": {"Sample_title": ["A"]}}]
        with self.assertRaises(ValueError):
            MODULE.exact_sample_join(["A"], samples, "Sample_title")

    def test_complete_series_and_global_outlier_counts(self):
        table = MODULE.numeric_table(self.series(['"p1"\t5.2\t1963485', '"p2"\t6.3\t7.4']), series_matrix=True)
        qc, _, _ = MODULE.numeric_qc(table, "processed_array")
        self.assertEqual((qc["features"], qc["columns"], qc["cells"]), (2, 2, 4))
        self.assertEqual(qc["absolute_value_threshold_audit"]["100"], {"cells": 1, "features": 1, "columns": 1})
        self.assertFalse(qc["expression_scale_qualified"])
        self.assertEqual(table["anomalies"][0]["literal_as_submitted"], "1963485")

    def test_missing_series_end_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.numeric_table(self.series(["p1\t1\t2"], end=False), series_matrix=True)

    def test_metadata_order_mismatch_rejected(self):
        path = self.series(["p1\t1\t2"])
        path.write_text(path.read_text().replace('!Sample_geo_accession\t"GSM1"\t"GSM2"', '!Sample_geo_accession\t"GSM2"\t"GSM1"'))
        with self.assertRaises(ValueError):
            MODULE.numeric_table(path, series_matrix=True)

    def test_ragged_matrix_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.numeric_table(self.write("\tA\tB\np1\t1\n"))

    def test_duplicate_columns_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.numeric_table(self.write("\tA\tA\np1\t1\t2\n"))

    def test_rna_complete_universe_zero_rows_and_scientific_integer(self):
        table = MODULE.numeric_table(self.write("\tA\tB\ng1\t0\t0\ng2\t1e3\t2.0\ng3\t5\t3\n"))
        qc, _, _ = MODULE.numeric_qc(table, "submitted_raw_gene_counts")
        self.assertTrue(qc["count_compatible"])
        self.assertEqual(qc["all_zero_feature_rows"], 1)
        self.assertEqual(qc["library_total_counts_min_median_max"], [5.0, 505.0, 1005.0])
        self.assertEqual(qc["literal_audit"]["exponent_literals"], 1)

    def test_count_fractions_negatives_nonfinite_do_not_qualify(self):
        table = MODULE.numeric_table(self.write("\tA\tB\ng1\t-1\t0.5\ng2\tNA\tinf\ng3\tword\t1\n"))
        qc, _, _ = MODULE.numeric_qc(table, "submitted_raw_gene_counts")
        self.assertFalse(qc["count_compatible"])
        self.assertEqual(qc["nonfinite_cells"], 3)
        self.assertEqual(qc["negative_cells"], 1)
        self.assertEqual(qc["literal_audit"]["invalid_literals"], 1)

    def test_exact_count_literal_cannot_be_hidden_by_float_rounding(self):
        table = MODULE.numeric_table(self.write("\tA\ng1\t9007199254740993\ng2\t1e-400\n"))
        qc, _, _ = MODULE.numeric_qc(table, "submitted_raw_gene_counts")
        self.assertFalse(qc["count_compatible"])
        self.assertEqual(qc["literal_audit"]["exact_literal_abs_gt_float64_integer_limit"], 1)
        self.assertEqual(qc["literal_audit"]["exact_literal_noninteger_values"], 1)

    def test_duplicate_feature_labels_preserved_not_summed(self):
        table = MODULE.numeric_table(self.write("\tA\tB\ng1\t1\t1\ng1\t2\t2\n"))
        qc, _, groups = MODULE.numeric_qc(table, "submitted_raw_gene_counts")
        self.assertEqual(qc["features"], 2)
        self.assertEqual(qc["duplicate_feature_label_groups"], 1)
        self.assertEqual(groups, [["A", "B"]])

    def test_annotation_multi_mapping_is_not_arbitrarily_collapsed(self):
        annotation = {"metadata": {"!Annotation_date": "historic"}, "columns": ["ID", "Gene ID", "Gene symbol"],
                      "rows": [{"ID": "p1", "Gene ID": "1", "Gene symbol": "OLD1"},
                               {"ID": "p2", "Gene ID": "1 /// 2", "Gene symbol": "OLD1 /// OLD2"},
                               {"ID": "p3", "Gene ID": "", "Gene symbol": ""},
                               {"ID": "p4", "Gene ID": "1", "Gene symbol": "OLD1"}]}
        qc, join = MODULE.annotation_join(["p1", "p2", "p3", "p4", "missing"], annotation)
        self.assertEqual(qc["join_status_counts"], {"one_entrez_id": 2, "multiple_entrez_ids": 1, "no_entrez_id": 1, "absent_from_annotation": 1})
        self.assertEqual(qc["entrez_ids_with_multiple_measurement_features"], 1)
        self.assertEqual(join[1]["historical_entrez_ids"], ["1", "2"])
        self.assertFalse(qc["current_gene_identifier_reconciliation_performed"])

    def test_annotation_completeness_and_duplicate_ids(self):
        valid = "^Annotation\n!Annotation_date = 2016\n!platform_table_begin\nID\tGene ID\tGene symbol\np1\t1\tOLD\n!platform_table_end\n"
        self.assertEqual(len(MODULE.annotation_table(self.write(valid))["rows"]), 1)
        with self.assertRaises(ValueError):
            MODULE.annotation_table(self.write(valid.replace("!platform_table_end\n", "")))
        with self.assertRaises(ValueError):
            MODULE.annotation_table(self.write(valid.replace("!platform_table_end", "p1\t2\tOTHER\n!platform_table_end")))

    def test_gzip_crc_and_truncation_not_ignored(self):
        compressed = self.path.with_suffix(".gz")
        compressed.write_bytes(gzip.compress(b"\tA\ng1\t1\n")[:-4])
        with self.assertRaises((EOFError, OSError)):
            MODULE.numeric_table(compressed)

    def test_no_target_list_or_effect_cli(self):
        text = (ROOT / "scripts/qualify_psc_liver_inputs.py").read_text()
        for prohibited in ["--contrast", "--score", "--target", "--model"]:
            self.assertNotIn('add_argument("' + prohibited, text)


    def test_bounded_reader_never_overreads_its_cap(self):
        import io
        class Raw:
            def __init__(self):
                self.stream, self.requests = io.BytesIO(b"abcdef"), []
            def read(self, amount, decode_content=False):
                self.requests.append(amount)
                return self.stream.read(amount)
        class Response:
            raw = Raw()
        response = Response()
        size, digest, error = MODULE.copy_bounded_body(response, self.path, 4)
        self.assertEqual((size, self.path.read_bytes(), error), (4, b"abcd", None))
        self.assertEqual(response.raw.stream.tell(), 4)
        self.assertEqual(response.raw.requests, [4])

    def test_bounded_reader_preserves_partial_bytes_after_error(self):
        class Raw:
            calls = 0
            def read(self, amount, decode_content=False):
                self.calls += 1
                if self.calls == 1:
                    return b"abc"
                raise OSError("synthetic transfer failure")
        class Response:
            raw = Raw()
        size, _, error = MODULE.copy_bounded_body(Response(), self.path, 100)
        self.assertEqual(size, 3)
        self.assertEqual(self.path.read_bytes(), b"abc")
        self.assertIn("synthetic transfer failure", error)

    def test_quick_preview_does_not_define_feature_universe(self):
        table = MODULE.numeric_table(self.series(["p1\t5.2\t1963485", "p2\t6.3\t7.4"]), series_matrix=True)
        samples = [{"accession": "GSM2", "table": [{"ID_REF": "p1", "VALUE": "1963485"}]}]
        qc, errors = MODULE.quick_table_crosscheck(samples, table)
        self.assertEqual(qc["cached_quick_cells_checked"], 1)
        self.assertEqual(qc["quick_cells_abs_gt_100_also_in_full_matrix"], 1)
        self.assertEqual(errors, [])
        self.assertEqual(len(table["features"]), 2)


@unittest.skipUnless((ROOT / "data/derived/psc-liver-input-qualification.json").exists(), "Qualification snapshot not present")
class LiverQualificationSnapshotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads((ROOT / "data/derived/psc-liver-input-qualification.json").read_text())

    def test_complete_matrix_dimensions_and_original_groups(self):
        rna = self.result["datasets"]["GSE303271"]
        adult = self.result["datasets"]["GSE159676"]
        self.assertEqual((rna["numeric_qc"]["features"], rna["numeric_qc"]["columns"]), (41096, 64))
        self.assertEqual(rna["metadata_structure"]["original_condition_counts"], {"PSC": 17, "ASC": 17, "AIH": 30})
        self.assertEqual((adult["numeric_qc"]["features"], adult["numeric_qc"]["columns"]), (17046, 33))
        self.assertEqual(adult["metadata_structure"]["original_condition_counts"]["Liver tissue healthy"], 6)

    def test_adult_integrity_failure_and_pairing_preserved(self):
        adult = self.result["datasets"]["GSE159676"]
        self.assertFalse(adult["numeric_qc"]["expression_scale_qualified"])
        self.assertEqual(adult["numeric_qc"]["absolute_value_threshold_audit"]["100"], {"cells": 4829, "columns": 33, "features": 3912})
        self.assertEqual(adult["metadata_structure"]["psc_complete_ab_blocks"], 6)
        self.assertEqual(adult["source_methods"]["independence"]["reported_PSC_patients"], 6)
        self.assertEqual(adult["metadata_structure"]["value_column_declarations"], ["Quantile normalized"])
        self.assertFalse(adult["eligible_for_automatic_expression_analysis"])

    def test_historical_annotation_keeps_multimapping_and_missingness(self):
        annotation = self.result["datasets"]["GSE159676"]["feature_annotation"]
        self.assertEqual(annotation["annotation_feature_rows"], 33297)
        self.assertEqual(annotation["join_status_counts"], {"one_entrez_id": 16806, "multiple_entrez_ids": 115, "no_entrez_id": 125})
        self.assertEqual(annotation["unique_entrez_ids_among_one_id_features"], 16020)
        self.assertFalse(annotation["current_gene_identifier_reconciliation_performed"])

    def test_atlas_is_metadata_only_and_overlap_is_not_replication(self):
        atlas = self.result["NoPSC_atlas"]
        self.assertFalse(atlas["count_export_qualified"])
        self.assertFalse(atlas["donor_map_qualified"])
        self.assertFalse(atlas["full_feature_download_qualified"])
        self.assertFalse(atlas["bulk_download_performed"])
        self.assertEqual(atlas["privacy_and_scope"]["named_gene_expression_payloads_requested"], 0)
        self.assertEqual(atlas["source_design"]["nuclear_reported_donors"]["total"], 16)
        self.assertEqual(atlas["source_design"]["spatial_reported_donors"]["total"], 30)

    def test_no_effects_corrections_exclusions_or_covariate_fabrication(self):
        scope = self.result["scope"]
        for key in ["primary_gene_contrasts_computed", "target_or_program_scores_computed", "models_fit", "source_values_repaired", "features_or_samples_excluded", "genes_aggregated", "new_gene_alias_corrections"]:
            self.assertIs(scope[key], False)
        for dataset in self.result["datasets"].values():
            self.assertFalse(dataset["metadata_structure"]["group_level_covariates_assigned_to_individuals"])
            self.assertFalse(dataset["source_methods"]["clinical_covariates"]["exact_sample_linked_age_sex_fibrosis_stage_batch_recovered"])

    def test_budget_with_conservative_atlas_charge(self):
        budget = self.result["provenance"]["acquisition_budget"]
        self.assertEqual(budget["new_body_bytes_read_and_retained_total"], 23735752)
        self.assertEqual(budget["conservative_atlas_inclusive_budget_charge"], 26376304)
        self.assertLess(budget["conservative_atlas_inclusive_budget_charge"], budget["overall_limit_bytes"])
        self.assertLessEqual(budget["atlas_conservative_full_advertised_prefix_GET_charge"], budget["atlas_limit_bytes"])
        self.assertTrue(budget["within_all_allocations"])

    def test_public_result_has_no_sample_join_or_expression_values(self):
        text = json.dumps(self.result, sort_keys=True)
        self.assertNotIn('"matrix_column_id_as_submitted"', text)
        self.assertNotIn('"literal_as_submitted"', text)
        self.assertNotIn('"patient_age"', text)
        self.assertNotIn('"per_sample"', text)
        self.assertTrue(all(row["path"].startswith("work/psc-liver-qualification/") for row in self.result["ignored_artifact_pins"]))

    def test_independent_source_audit_passes(self):
        self.assertEqual(self.result["validation"]["independent_audit_check_status_counts"], {"PASS": 82})
        self.assertTrue(self.result["closed_organoid_preservation"]["closed_versioned_organoid_pins_unchanged"])
        self.assertEqual(self.result["closed_organoid_preservation"]["prior_organoid_source_and_ignored_artifact_pins_checked"], 67)


if __name__ == "__main__":
    unittest.main()
