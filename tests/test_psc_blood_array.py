"""Synthetic, offline checks for the PSC adult-array replication arm.

No test reads the real group expression data. Source-shaped synthetic matrices
exercise the standalone runner, including the full 275-adult/95-child selection.
"""
from __future__ import annotations

import copy
import csv
import gzip
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts import analyze_psc_blood_array as array

ANNOTATION = {"table_format": "geo_annotation", "probe_column": "ID",
              "entrez_column": "Gene ID", "symbol_column": "Gene symbol"}


def matrix_text(probes, values, records):
    stream = io.StringIO()
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(["!Series_geo_accession", array.STUDY])
    writer.writerow(["!Series_sample_id", " ".join(row["sample_accession"] for row in records)])
    writer.writerow(["!Sample_geo_accession", *[row["sample_accession"] for row in records]])
    for key in ("title", "source_name_ch1", "platform_id", "data_processing"):
        writer.writerow(["!Sample_" + key, *[row[key][0] for row in records]])
    for index in range(3):
        writer.writerow(["!Sample_characteristics_ch1", *[row["characteristics_ch1"][index] for row in records]])
    writer.writerow(["!Sample_data_row_count", *([str(len(probes))] * len(records))])
    writer.writerow(["!series_matrix_table_begin"])
    writer.writerow(["ID_REF", *[row["sample_accession"] for row in records]])
    for probe, row in zip(probes, values):
        writer.writerow([probe, *[format(value, ".17g") for value in row]])
    writer.writerow(["!series_matrix_table_end"])
    return stream.getvalue()


def synthetic_records(full=True):
    records = []
    reverse = {value: key for key, value in array.CONDITION_MAP.items()}
    for age, counts in (("adult", array.ADULT_COUNTS), ("child", array.CHILD_COUNTS)):
        for group, n in counts.items():
            for _ in range(n if full else 2):
                accession = f"GSM{1000000 + len(records)}"
                records.append({"sample_accession": accession, "title": ["Synthetic " + accession],
                                "source_name_ch1": [group + ", " + age], "platform_id": [array.PLATFORM],
                                "characteristics_ch1": ["condition: " + reverse[group], "age group: " + age,
                                                         "tissue: whole blood"],
                                "data_processing": ["Quantile normalised using lumi library in R 2.15."]})
    return records


def synthetic_fixture(root):
    records = synthetic_records()
    n = len(records)
    rng = np.random.default_rng(45723)
    values = rng.normal(6, 0.3, (7, n))
    is_psc = np.array([row["characteristics_ch1"][0] == "condition: primary sclerosing cholangitis" for row in records])
    is_child = np.array([row["characteristics_ch1"][1] == "age group: child" for row in records])
    values[0, is_psc] += 1.5
    values[1, is_psc] += 1.5
    values[2, is_psc] -= 1.0
    values[:3, is_child] += 4.0  # Pediatric effects must never enter an adult comparison.
    values[5, :] = 6.0
    linear = np.exp2(values)
    probes = ["p1", "p2", "p3", "p_multi", "p_no", "p_const", "p_missing"]
    matrix_path = root / "matrix.txt.gz"
    matrix_path.write_bytes(gzip.compress(matrix_text(probes, linear, records).encode(), mtime=0))
    annotation_path = root / "annotation.txt.gz"
    annotation_rows = [
        ["p1", "101", "SOURCE_A"], ["p2", "101", "SOURCE_ALIAS_A"],
        ["p3", "102", "SOURCE_B"], ["p_multi", "103 /// 104", "AMBIG"],
        ["p_no", "", "DO_NOT_MAP_BY_SYMBOL"], ["p_const", "105", "CONSTANT"],
        ["platform_only", "900", "NOT_MEASURED"],
    ]
    annot = ("^Annotation\n!Annotation_platform = GPL10558\n"
             "!Annotation_platform_title = platform title, not its accession\n"
             "!Annotation_platform_organism = Homo sapiens\n"
             "!Annotation_date = Aug 09 2016\n!platform_table_begin\nID\tGene ID\tGene symbol\n"
             + "".join("\t".join(row) + "\n" for row in annotation_rows) + "!platform_table_end\n")
    annotation_path.write_bytes(gzip.compress(annot.encode(), mtime=0))
    metadata_path = root / "metadata.json"
    # Source order intentionally differs: accession identity, not row position, is required.
    metadata_path.write_text(json.dumps({array.STUDY: records[::-1]}))
    def pin(path):
        return {"path": path.name, "sha256": array.sha256_file(path), "bytes": path.stat().st_size,
                "url": "https://example.org/synthetic/" + path.name, "retrieved_at_utc": "2026-09-15T00:00:00Z"}
    plan = {"schema_version": 1, "study_accession": array.STUDY, "platform_accession": array.PLATFORM,
            "status": "frozen", "allow_effects": True, "frozen_at_utc": "2026-09-15T00:00:00Z",
            "code_sha256": array.sha256_file(Path(array.__file__)),
            "expected_total_samples": 370, "expected_adult_counts": array.ADULT_COUNTS,
            "expected_child_counts": array.CHILD_COUNTS, "expected_probe_count": len(probes),
            "expected_annotation_rows": len(annotation_rows),
            "source_quantity": "source_normalized_linear_intensity", "methods": array.METHODS,
            "annotation": ANNOTATION,
            "sources": {"series_matrix": pin(matrix_path), "platform_annotation": pin(annotation_path),
                        "sample_metadata": pin(metadata_path)}}
    plan_path = root / "plan.json"
    plan_path.write_text(json.dumps(plan))
    return plan_path, plan, records, linear, probes


class TestSourceParsing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.plan_path, self.plan, self.records, self.linear, self.probes = synthetic_fixture(self.root)

    def test_complete_matrix_and_permuted_metadata_join(self):
        matrix = array.load_series_matrix(self.root / "matrix.txt.gz")
        metadata = array.load_sample_metadata(self.root / "metadata.json")
        rows, audit = array.join_samples(matrix, metadata, 370, array.ADULT_COUNTS, array.CHILD_COUNTS)
        np.testing.assert_array_equal(matrix.values, self.linear)
        self.assertEqual(rows[0]["metadata_row_number_1based"], 370)
        self.assertEqual(audit["adult_samples_retained"], 275)
        self.assertEqual(audit["pediatric_samples_excluded"], 95)
        self.assertEqual(audit["adults_with_exact_age"], 0)
        self.assertEqual(audit["adults_with_sex"], 0)

    def test_duplicate_matrix_feature_rejected(self):
        path = self.root / "bad.txt"
        probes = self.probes.copy()
        probes[-1] = probes[0]
        path.write_text(matrix_text(probes, self.linear, self.records))
        with self.assertRaisesRegex(ValueError, "Duplicate matrix probe"):
            array.load_series_matrix(path)

    def test_duplicate_sample_and_matrix_metadata_order_rejected(self):
        text = gzip.decompress((self.root / "matrix.txt.gz").read_bytes()).decode()
        path = self.root / "bad.txt"
        old = self.records[0]["sample_accession"]
        other = self.records[1]["sample_accession"]
        path.write_text(text.replace(old, other))
        with self.assertRaisesRegex(ValueError, "Duplicate matrix sample"):
            array.load_series_matrix(path)
        # Change only the metadata accession row; the matrix order remains intact.
        lines = text.splitlines()
        row_index = next(i for i, line in enumerate(lines) if line.startswith("!Sample_geo_accession"))
        row = lines[row_index].split("\t")
        row[1], row[2] = row[2], row[1]
        lines[row_index] = "\t".join(row)
        path.write_text("\n".join(lines) + "\n")
        with self.assertRaisesRegex(ValueError, "order differ"):
            array.load_series_matrix(path)

    def test_truncated_wrong_width_and_transposed_matrices_rejected(self):
        text = gzip.decompress((self.root / "matrix.txt.gz").read_bytes()).decode()
        path = self.root / "bad.txt"
        path.write_text(text.replace("!series_matrix_table_end\n", ""))
        with self.assertRaisesRegex(ValueError, "truncated"):
            array.load_series_matrix(path)
        path.write_text(text.replace("p1\t", "p1\t1\t", 1))
        with self.assertRaisesRegex(ValueError, "incorrect width"):
            array.load_series_matrix(path)
        path.write_text(text.replace("ID_REF\t", "SAMPLE_REF\t", 1))
        with self.assertRaisesRegex(ValueError, "ID_REF"):
            array.load_series_matrix(path)

    def test_corrupt_gzip_is_not_accepted(self):
        path = self.root / "truncated.gz"
        path.write_bytes((self.root / "matrix.txt.gz").read_bytes()[:-8])
        with self.assertRaises((EOFError, OSError)):
            array.load_series_matrix(path)

    def test_source_title_condition_and_unknown_age_mismatch_rejected(self):
        matrix = array.load_series_matrix(self.root / "matrix.txt.gz")
        for key, value in (("title", ["altered"]), ("characteristics_ch1", ["condition: PSC", "age group: adult", "tissue: whole blood"])):
            records = copy.deepcopy(self.records)
            records[0][key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "Exact metadata join failed"):
                array.join_samples(matrix, records, 370, array.ADULT_COUNTS, array.CHILD_COUNTS)
        records = copy.deepcopy(self.records)
        records[0]["characteristics_ch1"][1] = "age group: not known"
        path = self.root / "age.txt"
        path.write_text(matrix_text(self.probes, self.linear, records))
        other = array.load_series_matrix(path)
        with self.assertRaisesRegex(ValueError, "Unknown age"):
            array.join_samples(other, records, 370, array.ADULT_COUNTS, array.CHILD_COUNTS)

    def test_missing_or_duplicate_independent_sample_rejected(self):
        matrix = array.load_series_matrix(self.root / "matrix.txt.gz")
        for records in (self.records[:-1], self.records[:-1] + [self.records[0]]):
            with self.assertRaises(ValueError):
                array.join_samples(matrix, records, 370, array.ADULT_COUNTS, array.CHILD_COUNTS)

    def test_exact_duplicate_titles_rejected_without_guessing_donors(self):
        records = copy.deepcopy(self.records)
        records[1]["title"] = records[0]["title"]
        path = self.root / "titles.txt"
        path.write_text(matrix_text(self.probes, self.linear, records))
        matrix = array.load_series_matrix(path)
        with self.assertRaisesRegex(ValueError, "duplicate exact sample titles"):
            array.join_samples(matrix, records, 370, array.ADULT_COUNTS, array.CHILD_COUNTS)

    def test_quick_soft_reads_metadata_not_preview_expression(self):
        path = self.root / "quick.txt"
        path.write_text("^SAMPLE = GSM1\n!Sample_geo_accession = GSM1\n!Sample_series_id = GSE119600\n"
                        "!Sample_title = exact title\n!sample_table_begin\nID_REF\tVALUE\n"
                        "some_gene\tnot even parsed as a number\n!sample_table_end\n")
        rows = array.load_sample_metadata(path)
        self.assertEqual(rows[0]["title"], ["exact title"])
        self.assertNotIn("some_gene", rows[0])

    def test_annotation_identity_and_raw_fields(self):
        rows, audit = array.load_platform_annotation(self.root / "annotation.txt.gz", ANNOTATION)
        self.assertEqual(audit["platform_accessions_as_published"], ["GPL10558"])
        self.assertEqual(rows[3]["Gene ID"], "103 /// 104")
        self.assertIn("!Annotation_date = Aug 09 2016", audit["metadata_as_published"])

    def test_annotation_duplicate_platform_mismatch_and_truncation(self):
        text = gzip.decompress((self.root / "annotation.txt.gz").read_bytes()).decode()
        path = self.root / "annot.txt"
        for bad in (text.replace("p2\t", "p1\t"), text.replace("= GPL10558", "= GPL570"),
                    text.replace("!platform_table_end\n", "")):
            path.write_text(bad)
            with self.assertRaises(ValueError):
                array.load_platform_annotation(path, ANNOTATION)


class TestProbeMappingAndScale(unittest.TestCase):
    def test_entrez_exact_single_unique_id_only(self):
        for raw in ("123", " 123 ", "123 /// 123", "123;123"):
            with self.subTest(raw=raw):
                row = array.parse_entrez_field(raw)
                self.assertEqual(row["entrez_gene_id"], "123")
                self.assertEqual(row["entrez_gene_id_raw"], raw)
        for raw in ("", "---", "NA", "123 /// 456", "123;456", "123|456", "123,456", "123 /// ---",
                    "123 ///", "ENSG00000000123", "0123", "123.0", "-123", "0", "SYMBOL"):
            with self.subTest(raw=raw):
                self.assertEqual(array.parse_entrez_field(raw)["entrez_gene_id"], "")

    def test_union_join_retains_all_unmapped_and_platform_only_rows(self):
        probes = ["p2", "p1", "multi", "blank", "missing"]
        annotation = [{"ID": "p1", "Gene ID": "10", "Gene symbol": "SYM"},
                      {"ID": "p2", "Gene ID": "10", "Gene symbol": "OLD_LABEL"},
                      {"ID": "multi", "Gene ID": "10 /// 20", "Gene symbol": "SYM"},
                      {"ID": "blank", "Gene ID": "", "Gene symbol": "SYM"},
                      {"ID": "platform_only", "Gene ID": "20", "Gene symbol": "OTHER"}]
        rows, genes, audit = array.build_feature_join(probes, annotation, ANNOTATION)
        self.assertEqual(len(rows), 6)
        self.assertEqual([row["entrez_gene_id"] for row in genes], ["10"])
        self.assertEqual(genes[0]["eligible_probe_count"], 2)
        self.assertEqual(json.loads(genes[0]["gene_symbols_as_published_json"]), ["OLD_LABEL", "SYM"])
        self.assertEqual(audit["join_status_counts"], {"matrix_and_platform": 4, "matrix_only": 1, "platform_only": 1})

    def test_median_of_transformed_probes_not_log_of_linear_median(self):
        probes = ["p2", "p1", "unmapped"]
        annotation = [{"ID": "p1", "Gene ID": "10", "Gene symbol": "S"},
                      {"ID": "p2", "Gene ID": "10", "Gene symbol": "S"},
                      {"ID": "unmapped", "Gene ID": "", "Gene symbol": "S"}]
        rows, _, _ = array.build_feature_join(probes, annotation, ANNOTATION)
        linear = np.array([[2., 16.], [8., 4.], [1024., 1024.]])
        genes, values = array.aggregate_probes(array.log2_no_shift(linear), probes, rows)
        self.assertEqual(genes, ["10"])
        np.testing.assert_array_equal(values, [[2., 3.]])
        self.assertNotEqual(values[0, 0], np.log2(np.median(linear[:2, 0])))

    def test_no_effect_or_expression_based_probe_choice(self):
        annotation = [{"ID": f"p{i}", "Gene ID": "10", "Gene symbol": "S"} for i in range(3)]
        probes = [row["ID"] for row in annotation]
        rows, _, _ = array.build_feature_join(probes, annotation, ANNOTATION)
        logs = np.array([[100., -100.], [3., 3.], [4., 4.]])
        _, values = array.aggregate_probes(logs, probes, rows)
        np.testing.assert_array_equal(values, [[4., 3.]])

    def test_no_offsets_imputation_or_nonpositive_exclusions(self):
        for value in (0., -1., np.nan, np.inf, -np.inf):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "offsets are forbidden"):
                array.log2_no_shift(np.array([[value, 2.]]))
        values = np.array([[np.nextafter(0., 1.), 1., 16.]])
        np.testing.assert_array_equal(array.log2_no_shift(values), np.log2(values))
        audit = array.linear_value_audit(np.array([[0., -1., np.nan, 2.]]))
        self.assertEqual((audit["zero_cells"], audit["negative_cells"], audit["nonfinite_cells"]), (1, 1, 1))
        self.assertFalse(audit["log2_without_shift_valid"])


class TestStatistics(unittest.TestCase):
    def test_welch_matches_scipy_unequal_n_and_variance(self):
        rng = np.random.default_rng(931)
        case = rng.normal(2, 2.1, (30, 7))
        reference = rng.normal(0, 0.8, (30, 19))
        observed = array.pairwise_statistics(case, reference)
        expected = stats.ttest_ind(case, reference, equal_var=False, axis=1)
        np.testing.assert_allclose(observed["welch_t"], expected.statistic, rtol=3e-14)
        np.testing.assert_allclose(observed["welch_p"], expected.pvalue, rtol=3e-14)
        np.testing.assert_allclose(observed["welch_df"], expected.df, rtol=3e-14)
        np.testing.assert_allclose(observed["log2_intensity_difference"], case.mean(axis=1) - reference.mean(axis=1))

    def test_hc3_matches_independent_ols_sandwich_and_normal_inference(self):
        rng = np.random.default_rng(101)
        case, reference = rng.normal(2, 1.2, (8, 5)), rng.normal(0, 0.5, (8, 9))
        observed = array.pairwise_statistics(case, reference)
        x = np.column_stack([np.ones(14), [1.] * 5 + [0.] * 9])
        bread = np.linalg.inv(x.T @ x)
        leverage = np.sum(x * (x @ bread), axis=1)
        for i in range(8):
            y = np.r_[case[i], reference[i]]
            beta = bread @ x.T @ y
            residual = y - x @ beta
            meat = x.T @ (((residual / (1 - leverage)) ** 2)[:, None] * x)
            covariance = bread @ meat @ bread
            se = np.sqrt(covariance[1, 1])
            z = beta[1] / se
            self.assertAlmostEqual(observed["hc3_se"][i], se, places=13)
            self.assertAlmostEqual(observed["hc3_z"][i], z, places=13)
            self.assertAlmostEqual(observed["hc3_p"][i], 2 * stats.norm.sf(abs(z)), places=13)
        self.assertTrue(np.all(observed["hc3_se"] > observed["welch_se"]))

    def test_one_constant_group_remains_estimable_both_constant_are_missing(self):
        observed = array.pairwise_statistics(np.array([[2., 2., 2.], [3., 3., 3.]]),
                                            np.array([[1., 2., 3., 4.], [1., 1., 1., 1.]]))
        self.assertTrue(np.isfinite(observed["welch_p"][0]))
        self.assertTrue(np.isnan(observed["welch_p"][1]))
        self.assertTrue(np.isnan(observed["welch_bh_q"][1]))
        self.assertEqual(observed["welch_status"][1], "not_estimable_zero_or_invalid_variance")

    def test_bh_complete_family_missing_ties_and_bounds(self):
        np.testing.assert_allclose(array.bh_adjust(np.array([0.01, np.nan, 0.04])), [0.03, np.nan, 0.06], equal_nan=True)
        np.testing.assert_allclose(array.bh_adjust(np.array([0.01, 0.01, 1., 0.])), [0.04 / 3, 0.04 / 3, 1., 0.])
        self.assertTrue(np.all(np.isnan(array.bh_adjust(np.array([np.nan, np.nan])))))
        for p in (np.array([]), np.array([1.1]), np.array([-0.1]), np.array([np.inf])):
            with self.assertRaises(ValueError):
                array.bh_adjust(p)

    def test_tests_reject_missing_data_and_singleton_group(self):
        for case in (np.array([[1.]]), np.array([[1., np.nan]])):
            with self.assertRaises(ValueError):
                array.pairwise_statistics(case, np.array([[1., 2., 3.]]))


class TestStandaloneRun(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.plan_path, self.plan, self.records, self.linear, self.probes = synthetic_fixture(self.root)

    def run_fixture(self, label="run", validate_only=False):
        return array.run_analysis(self.plan_path, self.root / (label + "-public"),
                                  self.root / "work" / label, validate_only=validate_only, root=self.root)

    def test_full_synthetic_effects_complete_family_and_private_sample_export(self):
        audit = self.run_fixture()
        self.assertTrue(audit["group_effects_computed"])
        self.assertEqual(audit["effect_rows"], 3 * 4)
        with gzip.open(self.root / "run-public/array-gene-effects.csv.gz", "rt") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual({row["contrast"] for row in rows}, {a + "_vs_" + b for a, b in array.CONTRASTS})
        primary = {row["entrez_id"]: row for row in rows if row["contrast"] == "PSC_vs_Control"}
        self.assertGreater(float(primary["101"]["estimate"]), 1.)
        self.assertLess(float(primary["102"]["estimate"]), -0.5)
        self.assertEqual(primary["105"]["pvalue"], "")
        self.assertTrue(all(row["multiple_testing_family_size"] == "3" for row in rows))
        self.assertEqual(primary["101"]["n_case"], "45")
        self.assertEqual(primary["101"]["n_reference"], "47")
        for row in rows:
            self.assertEqual(row["pvalue"], row["welch_p"])
            self.assertEqual(row["padj"], row["welch_bh_q"])
        for path in (self.root / "run-public").iterdir():
            text = gzip.decompress(path.read_bytes()).decode() if path.suffix == ".gz" else path.read_text()
            self.assertNotIn("GSM1000000", text)
            self.assertNotIn("Synthetic GSM", text)
        gene_values = np.load(self.root / "work/run/array-adult-entrez-log2.npy", allow_pickle=False)
        self.assertEqual(gene_values.shape, (3, 275))
        self.assertTrue((self.root / "work/run/array-sample-join.csv").exists())

    def test_draft_validation_never_calls_effects_or_exports_expression(self):
        self.plan.update(status="draft", allow_effects=False)
        self.plan.pop("code_sha256")
        self.plan_path.write_text(json.dumps(self.plan))
        with patch.object(array, "analyse_contrasts", side_effect=AssertionError("No effects allowed")):
            audit = self.run_fixture(validate_only=True)
        self.assertEqual(audit["stage"], "validated_without_effects")
        self.assertFalse((self.root / "run-public/array-gene-effects.csv.gz").exists())
        self.assertFalse((self.root / "work/run/array-adult-entrez-log2.npy").exists())
        with self.assertRaisesRegex(ValueError, "Effects are blocked"):
            self.run_fixture(label="blocked")

    def test_strict_source_plan_and_code_hashes(self):
        with self.assertRaisesRegex(ValueError, "Outer plan SHA"):
            array.run_analysis(self.plan_path, self.root / "public", self.root / "work/run",
                               root=self.root, expected_plan_sha256="0" * 64)
        self.plan["code_sha256"] = "0" * 64
        self.plan_path.write_text(json.dumps(self.plan))
        with self.assertRaisesRegex(ValueError, "implementation SHA"):
            self.run_fixture()
        self.plan["code_sha256"] = array.sha256_file(Path(array.__file__))
        self.plan["sources"]["series_matrix"]["sha256"] = "0" * 64
        self.plan_path.write_text(json.dumps(self.plan))
        with self.assertRaisesRegex(ValueError, "Source SHA"):
            self.run_fixture()

    def test_nonpositive_unmapped_pediatric_cell_blocks_whole_run(self):
        # Even this cell, outside all analysis genes/adults, cannot trigger an offset or hidden exclusion.
        self.linear[4, -1] = 0.
        path = self.root / "matrix.txt.gz"
        path.write_bytes(gzip.compress(matrix_text(self.probes, self.linear, self.records).encode(), mtime=0))
        self.plan["sources"]["series_matrix"].update(sha256=array.sha256_file(path), bytes=path.stat().st_size)
        self.plan_path.write_text(json.dumps(self.plan))
        with patch.object(array, "analyse_contrasts", side_effect=AssertionError("No effects allowed")):
            audit = self.run_fixture()
        self.assertEqual(audit["stage"], "blocked_nonpositive_or_nonfinite_source_intensity")
        self.assertFalse(audit["group_effects_computed"])
        self.assertEqual(audit["numeric_audit"]["zero_cells"], 1)

    def test_offline_replay_outputs_are_byte_identical(self):
        self.run_fixture(label="a")
        self.run_fixture(label="b")
        for left in (self.root / "a-public").iterdir():
            right = self.root / "b-public" / left.name
            self.assertEqual(left.read_bytes(), right.read_bytes(), left.name)
        for left in (self.root / "work/a").iterdir():
            self.assertEqual(left.read_bytes(), (self.root / "work/b" / left.name).read_bytes(), left.name)

    def test_patient_work_directory_must_be_ignored_and_runs_not_overwritten(self):
        with self.assertRaisesRegex(ValueError, "ignored work"):
            array.run_analysis(self.plan_path, self.root / "public", self.root / "patients", root=self.root)
        self.run_fixture()
        with self.assertRaisesRegex(ValueError, "existing run directory"):
            self.run_fixture()

    def test_outer_array_section_and_frozen_timestamp_gate(self):
        self.plan_path.write_text(json.dumps({"array_analysis": self.plan, "other_arm": "unused"}))
        audit = self.run_fixture()
        self.assertTrue(audit["group_effects_computed"])
        for invalid in (None, "2026-09-15T00:00:00", "2026-09-15T01:00:00+01:00"):
            plan = copy.deepcopy(self.plan)
            plan["frozen_at_utc"] = invalid
            with self.assertRaisesRegex(ValueError, "timestamp"):
                array.validate_plan(plan, validate_only=False, code_path=Path(array.__file__))


if __name__ == "__main__":
    unittest.main()
