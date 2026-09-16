"""Synthetic offline tests for source-only organoid qualification.

No gene-specific study effects, donor contrasts, or program scores are read.
"""
from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts import qualify_psc_organoid_inputs as q


class SparseQualificationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.path = self.root / "matrix.mtx.gz"

    def tearDown(self):
        self.tmp.cleanup()

    def matrix(self, entries, dimensions="3 4 5", banner=None):
        text = (banner or "%%MatrixMarket matrix coordinate integer general") + "\n" + dimensions + "\n" + entries
        self.path.write_bytes(gzip.compress(text.encode(), mtime=0))
        return self.path

    def scan(self, entries="1 1 2\n3 1 1\n2 2 0\n2 3 4\n1 4 3\n", **kwargs):
        self.matrix(entries, **{key: kwargs.pop(key) for key in list(kwargs) if key in {"dimensions", "banner"}})
        return q.scan_matrix(self.path, [0, 0, 1, 1], 2, expected_shape=(3, 4), **kwargs)

    def test_all_entries_axis_and_sums(self):
        audit, counts, totals, detected, feature_cells = self.scan(chunk_lines=2)
        np.testing.assert_array_equal(counts, [[2, 3], [0, 4], [1, 0]])
        np.testing.assert_array_equal(totals, [3, 0, 4, 3])
        np.testing.assert_array_equal(detected, [2, 0, 1, 1])
        np.testing.assert_array_equal(feature_cells, [2, 1, 1])
        self.assertEqual(audit["explicit_zero_entries"], 1)
        self.assertEqual(audit["full_released_matrix_sum"], 10)
        self.assertTrue(audit["gzip_crc_and_eof_verified"])
        self.assertEqual(audit["zero_total_cells"], 1)

    def test_arbitrary_coordinate_order_is_supported(self):
        audit, counts, *_ = self.scan("2 3 4\n1 1 2\n1 4 3\n3 1 1\n2 2 0\n", chunk_lines=2)
        self.assertFalse(audit["strictly_column_major"])
        np.testing.assert_array_equal(counts, [[2, 3], [0, 4], [1, 0]])

    def test_chunks_do_not_change_outputs(self):
        left = self.scan(chunk_lines=1)
        right = self.scan(chunk_lines=20)
        self.assertEqual(left[0], right[0])
        for a, b in zip(left[1:], right[1:], strict=True):
            np.testing.assert_array_equal(a, b)

    def test_last_line_without_newline(self):
        self.assertEqual(self.scan("1 1 2\n3 1 1\n2 2 0\n2 3 4\n1 4 3")[0]["observed_entries"], 5)

    def test_duplicate_in_chunk_fails(self):
        with self.assertRaisesRegex(ValueError, "Duplicate sparse"):
            self.scan("1 1 1\n1 1 2\n2 2 1\n2 3 1\n1 4 1\n", chunk_lines=10)

    def test_duplicate_across_chunks_fails(self):
        with self.assertRaisesRegex(ValueError, "Duplicate sparse"):
            self.scan("1 1 1\n2 2 1\n1 1 2\n2 3 1\n1 4 1\n", chunk_lines=2)

    def test_multiple_coordinates_in_same_byte_do_not_collide(self):
        _, counts, *_ = self.scan("1 1 1\n2 1 2\n3 1 3\n1 2 4\n2 2 5\n", chunk_lines=2)
        np.testing.assert_array_equal(counts[:, 0], [5, 7, 3])

    def test_invalid_quantities_fail(self):
        for value in ["nan", "NaN", "inf", "-inf", "1.5", "1e0", "1e309", "-1", "9999999999999999999"]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.scan(f"1 1 {value}\n", dimensions="3 4 1")

    def test_integer_overflow_guard(self):
        with self.assertRaisesRegex(ValueError, "overflow"):
            self.scan("1 1 999999999999999999\n", dimensions="3 4 100")

    def test_bounds_fail(self):
        for row, col in [(0, 1), (4, 1), (1, 0), (1, 5), (-1, 2)]:
            with self.subTest(row=row, col=col), self.assertRaisesRegex(ValueError, "bounds"):
                self.scan(f"{row} {col} 1\n", dimensions="3 4 1")

    def test_wrong_stored_count_fails(self):
        with self.assertRaisesRegex(ValueError, "declared"):
            self.scan("1 1 1\n", dimensions="3 4 2")
        with self.assertRaisesRegex(ValueError, "declared"):
            self.scan("1 1 1\n2 1 1\n", dimensions="3 4 1")

    def test_axis_dimension_mismatch_fails(self):
        with self.assertRaisesRegex(ValueError, "axes"):
            self.scan("1 1 1\n", dimensions="3 3 1")

    def test_real_and_symmetric_formats_not_silently_converted(self):
        for banner in ["%%MatrixMarket matrix coordinate real general", "%%MatrixMarket matrix coordinate integer symmetric"]:
            with self.subTest(banner=banner), self.assertRaises(ValueError):
                self.scan("1 1 1\n", dimensions="3 4 1", banner=banner)

    def test_malformed_extra_field_fails(self):
        with self.assertRaises(ValueError):
            self.scan("1 1 1 1\n", dimensions="3 4 1")

    def test_truncated_gzip_detected(self):
        self.matrix("1 1 1\n", dimensions="3 4 1")
        self.path.write_bytes(self.path.read_bytes()[:-5])
        with self.assertRaises((EOFError, OSError)):
            q.scan_matrix(self.path, [0, 0, 1, 1], 2)

    def test_bad_gzip_crc_detected(self):
        self.matrix("1 1 1\n", dimensions="3 4 1")
        body = bytearray(self.path.read_bytes())
        body[-8] ^= 1
        self.path.write_bytes(body)
        with self.assertRaises(OSError):
            q.scan_matrix(self.path, [0, 0, 1, 1], 2)

    def test_gzip_output_is_deterministic(self):
        first = self.root / "one.tsv.gz"
        second = self.root / "two.tsv.gz"
        q.write_tsv(first, ["id", "value"], [["A", 2]])
        q.write_tsv(second, ["id", "value"], [["A", 2]])
        self.assertEqual(first.read_bytes(), second.read_bytes())
        self.assertEqual(q.read_gzip_rows(first), [["id", "value"], ["A", "2"]])


class SourceJoinTests(unittest.TestCase):
    def data(self):
        header = ["", "orig.ident", "Sample", "Treatment", "condition", "conditionTreatment",
                  "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT"]
        rows, family = [], []
        for disease, diagnosis in [("PSC", "primary sclerosing cholangitis (PSC)"), ("nonPSC", "non-PSC")]:
            donor = disease + "_1"
            condition = "PSC" if disease == "PSC" else "Non.PSC"
            for arm, treatment in [("CNT", "vehicle"), ("IL17", "IL17A")]:
                title = donor + "_" + arm
                rows.append([title + "_barcode", title, donor, arm, condition, condition + "_" + arm, "100", "5", "10", "2"])
                family.append({"accession": "GSM" + str(len(family)), "fields": {
                    "title": [title], "characteristics_ch1": ["patient diagnosis: " + diagnosis,
                    "passage: passage 3", "treatment: " + treatment]}})
        return [header, *rows], [row[0] for row in rows], family

    def test_exact_join_not_order(self):
        rows, barcodes, family = self.data()
        header, aligned, libraries, audit = q.qualify_metadata([rows[0], *rows[:0:-1]], barcodes, family)
        self.assertFalse(audit["metadata_order_matches_barcode_axis"])
        self.assertEqual([row[""] for row in aligned], barcodes)
        self.assertEqual(audit["donors"], 2)
        self.assertEqual(audit["libraries"], 4)
        self.assertEqual(audit["complete_donor_treatment_pairs"], 2)
        self.assertEqual(audit["source_value_corrections"], [])
        self.assertEqual(libraries[-1]["source_condition"], "Non.PSC")
        self.assertEqual(libraries[-1]["disease"], "nonPSC")

    def test_duplicate_metadata_barcode_fails(self):
        rows, barcodes, family = self.data()
        rows[2][0] = rows[1][0]
        with self.assertRaisesRegex(ValueError, "Duplicate metadata barcode"):
            q.qualify_metadata(rows, barcodes, family)

    def test_missing_barcode_fails(self):
        rows, barcodes, family = self.data()
        barcodes[-1] = "not_present"
        with self.assertRaisesRegex(ValueError, "sets differ"):
            q.qualify_metadata(rows, barcodes, family)

    def test_source_treatment_not_inferred_from_title(self):
        rows, barcodes, family = self.data()
        rows[1][3] = "IL17"
        with self.assertRaisesRegex(ValueError, "disagrees"):
            q.qualify_metadata(rows, barcodes, family)

    def test_source_donor_not_inferred_from_title(self):
        rows, barcodes, family = self.data()
        rows[1][2] = "PSC_2"
        with self.assertRaisesRegex(ValueError, "disagrees"):
            q.qualify_metadata(rows, barcodes, family)

    def test_incomplete_pair_fails(self):
        rows, barcodes, family = self.data()
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            q.qualify_metadata(rows[:-1], barcodes[:-1], family[:-1])

    def test_conflicting_GSM_accession_field_fails(self):
        rows, barcodes, family = self.data()
        family[0]["fields"]["geo_accession"] = ["GSM-other"]
        with self.assertRaisesRegex(ValueError, "accession field"):
            q.qualify_metadata(rows, barcodes, family)

    def test_duplicate_GSM_accession_fails(self):
        rows, barcodes, family = self.data()
        family[1]["accession"] = family[0]["accession"]
        with self.assertRaisesRegex(ValueError, "GSM accession"):
            q.qualify_metadata(rows, barcodes, family)

    def test_extra_conflicting_GSM_characteristic_fails(self):
        for value in ["treatment: IL17A", "patient diagnosis: non-PSC", "passage: passage 5"]:
            rows, barcodes, family = self.data()
            family[0]["fields"]["characteristics_ch1"].append(value)
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "characteristics"):
                q.qualify_metadata(rows, barcodes, family)

    def test_raw_matching_margins_alone_do_not_authorize_inference(self):
        rows = [{"nCount_RNA": "10", "nCount_SCT": "8", "nFeature_RNA": "4", "nFeature_SCT": "3"}]
        check = q.compare_metadata_counts(rows, np.asarray([10]), np.asarray([4]))
        self.assertTrue(check["interpretation"]["raw_RNA_margins_match"])
        self.assertFalse(check["interpretation"]["raw_count_inference_eligible"])

    def test_GSM_diagnosis_disagreement_fails(self):
        rows, barcodes, family = self.data()
        family[0]["fields"]["characteristics_ch1"][0] = "patient diagnosis: non-PSC"
        with self.assertRaisesRegex(ValueError, "characteristics"):
            q.qualify_metadata(rows, barcodes, family)

    def test_feature_identifiers_not_silently_collapsed(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            q.unique_index(["A", "A"], "feature")

    def test_SCT_integer_counts_are_not_raw_RNA(self):
        rows, barcodes, family = self.data()
        _, aligned, _, _ = q.qualify_metadata(rows, barcodes, family)
        check = q.compare_metadata_counts(aligned, np.full(4, 10), np.full(4, 2))
        self.assertEqual(check["nCount_SCT"]["matching_cells"], 4)
        self.assertEqual(check["nCount_RNA"]["matching_cells"], 0)
        self.assertTrue(check["interpretation"]["SCT_margins_match"])
        self.assertFalse(check["interpretation"]["raw_count_inference_eligible"])

    def test_ambiguous_assay_not_claimed_raw(self):
        rows = [{"nCount_RNA": "10", "nCount_SCT": "10", "nFeature_RNA": "2", "nFeature_SCT": "2"}]
        check = q.compare_metadata_counts(rows, np.asarray([10]), np.asarray([2]))
        self.assertFalse(check["interpretation"]["raw_count_inference_eligible"])
        self.assertEqual(check["interpretation"]["status"], "assay_quantity_unresolved")

    def test_scientific_notation_is_exact_and_source_unchanged(self):
        rows = [{"x": "1e+05"}, {"x": "1.0"}]
        np.testing.assert_array_equal(q.exact_nonnegative_integers(rows, "x"), [100000, 1])
        self.assertEqual(rows[0]["x"], "1e+05")

    def test_invalid_numeric_metadata_not_rounded(self):
        for value in ["1.5", "NaN", "Infinity", "-1", "1e30", "missing"]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                q.exact_nonnegative_integers([{"x": value}], "x")


class ArchiveAndOriginalBarcodeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def archive(self, name="safe.tsv.gz", kind=None, link=None):
        import io
        import tarfile
        path = self.root / "test.tar"
        with tarfile.open(path, "w") as archive:
            member = tarfile.TarInfo(name)
            if kind is not None:
                member.type = kind
                member.linkname = link or "outside"
            else:
                member.size = 3
            archive.addfile(member, io.BytesIO(b"abc") if member.isfile() else None)
        return path, {name: member.size}

    def test_regular_allowlisted_member_only(self):
        path, expected = self.archive()
        ledger = q.safe_extract_processed_archive(path, expected, self.root / "out")
        self.assertEqual(ledger[0]["bytes"], 3)
        self.assertEqual((self.root / "out/safe.tsv.gz").read_bytes(), b"abc")

    def test_path_traversal_and_absolute_path_rejected(self):
        for name in ["../outside", "/tmp/outside", "nested/file", "..\\outside"]:
            with self.subTest(name=name):
                path, expected = self.archive(name)
                with self.assertRaisesRegex(ValueError, "Unsafe"):
                    q.safe_extract_processed_archive(path, expected, self.root / "out")

    def test_symlink_hardlink_directory_rejected(self):
        import tarfile
        for kind in [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.DIRTYPE]:
            with self.subTest(kind=kind):
                path, expected = self.archive(kind=kind)
                with self.assertRaisesRegex(ValueError, "Unsafe"):
                    q.safe_extract_processed_archive(path, expected, self.root / "out")

    def test_unlisted_member_rejected(self):
        path, _ = self.archive()
        with self.assertRaisesRegex(ValueError, "allowlist"):
            q.safe_extract_processed_archive(path, {"other": 3}, self.root / "out")

    def test_member_size_mismatch_rejected(self):
        path, _ = self.archive()
        with self.assertRaisesRegex(ValueError, "size"):
            q.safe_extract_processed_archive(path, {"safe.tsv.gz": 2}, self.root / "out")

    def test_preexisting_destination_symlink_rejected(self):
        path, expected = self.archive()
        out = self.root / "out"
        out.mkdir()
        (out / "safe.tsv.gz").symlink_to(self.root / "outside")
        with self.assertRaisesRegex(ValueError, "destination"):
            q.safe_extract_processed_archive(path, expected, out)

    def test_exact_original_library_join_preserves_barcodes(self):
        rows = [{"": "AAA-1_7", "orig.ident": "PSC_1_CNT"}, {"": "BBB-1_7", "orig.ident": "PSC_1_CNT"}]
        indices, mapping = q.join_original_barcodes(rows, ["BBB-1", "CCC-1", "AAA-1"])
        self.assertEqual(indices, [2, 0])
        self.assertEqual(mapping[0]["source_merged_barcode"], "AAA-1_7")
        self.assertEqual(mapping[0]["source_original_barcode"], "AAA-1")
        self.assertEqual(mapping[0]["removed_merge_suffix"], "_7")

    def test_missing_original_source_cell_rejected(self):
        with self.assertRaisesRegex(ValueError, "absent"):
            q.join_original_barcodes([{"": "AAA-1_7", "orig.ident": "L"}], ["BBB-1"])

    def test_ambiguous_suffix_and_duplicate_cell_rejected(self):
        for cells in [["AAA-1", "BBB-1_7"], ["AAA-1_x", "BBB-1_7"],
                      ["AAA-1_7", "BBB-1_8"], ["AAA-1_7", "AAA-1_7"]]:
            with self.subTest(cells=cells), self.assertRaises(ValueError):
                q.join_original_barcodes([{ "": cell, "orig.ident": "L"} for cell in cells], ["AAA-1", "BBB-1"])


class RootScopeDecisionTests(unittest.TestCase):
    """Test the final binding gate with tiny mocked upstream-qualified inputs."""
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.work = Path(self.tmp.name)
        files = ["raw-rna-pseudobulk.tsv.gz", "raw-rna-features.tsv", "raw-rna-libraries.tsv", "raw-rna-cell-source-joins.tsv.gz"]
        for name in files:
            (self.work / name).write_text("synthetic upstream-validated input\n")
        (self.work / "raw-rna-features.tsv").write_text("source_feature_id\tsource_gene_label\tsource_feature_type\nE1\tA\tGene Expression\nE2\tA\tGene Expression\n")
        self.summary = {
            "original_RNA": {
                "work_artifacts": [{"file": name, "bytes": (self.work / name).stat().st_size, "sha256": q.sha256(self.work / name)} for name in files],
                "compressed_members_checked": 48, "all_member_gzip_crc_and_eof_verified": True,
                "retained_cells": 46343, "barcode_join": {"all_released_cells_found_once": True},
                "full_feature_count": 36601, "unique_feature_IDs": 36601, "Ensembl_gene_ID_rows": 36601,
                "all_feature_axes_exactly_identical": True,
                "full_sparse_QC": {key: 0 for key in ["duplicate_coordinates", "negative_entries", "nonfinite_entries", "noninteger_entries", "out_of_bounds_coordinates"]},
                "retained_count_total": 101, "retained_positive_entries": 52,
                "quantity_validation": {"interpretation": {"raw_RNA_margins_match": False},
                                        "nCount_RNA": {"source_sum": 100}, "nFeature_RNA": {"source_sum": 50}}},
            "joins": {"complete_donor_treatment_pairs": 8, "libraries": 16},
            "quantity_validation": {"interpretation": {"SCT_margins_match": True, "raw_RNA_margins_match": False}},
            "source_manifest": {"sha256": "synthetic-pin"},
            "validation": {"independent_verification": {"status": "passed", "total_mismatches": 0, "pseudobulk_count_values": 36601 * 16},
                           "offline_replay": {"status": "passed", "count_and_numeric_QC_mismatches": 0}}}

    def tearDown(self):
        self.tmp.cleanup()

    def test_new_raw_universe_not_historical_reconstruction(self):
        result = q.apply_root_scope_decision(self.summary, self.work)
        self.assertTrue(result["raw_RNA_count_input_qualified"])
        self.assertFalse(result["original_RNA"]["historical_RNA_metadata_exactly_reproduced"])
        self.assertFalse(result["ready_for_raw_count_inference"])
        decision = result["eligibility_decision"]
        self.assertFalse(decision["effects_authorized_by_this_qualification"])
        self.assertTrue(decision["root_final_analysis_freeze_still_required"])
        self.assertEqual(decision["counts_above_historical_RNA_metadata"], 1)
        self.assertEqual(decision["positive_entries_above_historical_RNA_metadata"], 2)
        self.assertEqual(len(decision["input_artifact_hashes"]), 4)
        self.assertEqual(result["original_RNA"]["duplicate_source_gene_label_groups"], 1)
        self.assertEqual(result["original_RNA"]["feature_rows_with_duplicated_source_gene_labels"], 2)

    def test_changed_bound_count_file_rejects(self):
        (self.work / "raw-rna-pseudobulk.tsv.gz").write_text("changed\n")
        with self.assertRaisesRegex(ValueError, "changed"):
            q.apply_root_scope_decision(self.summary, self.work)

    def test_independent_verification_is_required(self):
        del self.summary["validation"]["independent_verification"]
        result = q.apply_root_scope_decision(self.summary, self.work)
        self.assertFalse(result["raw_RNA_count_input_qualified"])
        self.assertEqual(result["eligibility_decision"]["decision"], "pending_qualification_checks")

    def test_sparse_QC_failure_blocks_eligibility(self):
        self.summary["original_RNA"]["full_sparse_QC"]["duplicate_coordinates"] = 1
        result = q.apply_root_scope_decision(self.summary, self.work)
        self.assertFalse(result["raw_RNA_count_input_qualified"])


if __name__ == "__main__":
    unittest.main()
