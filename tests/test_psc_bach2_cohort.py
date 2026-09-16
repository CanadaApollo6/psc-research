"""Synthetic cohort-interface tests; no network or real expression effects."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import tarfile
import tempfile
import threading
import unittest
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

from scripts import qualify_psc_bach2_cohort as q

NAME = "sample01_filtered_feature_bc_matrix.tar.gz"
URL = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-14013/" + NAME
BARE1 = "A" * 16 + "-1"
BARE2 = "C" * 16 + "-1"


class FakeResponse:
    def __init__(self, body=b"abc", status=200, headers=None, url=URL):
        self.body = io.BytesIO(body)
        self.status = status
        self.headers = {"Content-Length": str(len(body))} if headers is None else headers
        self.url = url
        self.read_calls = 0
    def read(self, size=-1):
        self.read_calls += 1
        return self.body.read(size)
    def close(self):
        pass
    def geturl(self):
        return self.url


class FakeOpener:
    def __init__(self, responses):
        self.responses = iter(responses)
        self.calls = []
    def open(self, request, timeout):
        self.calls.append(request.full_url)
        return next(self.responses)


def synthetic_sources(samples=("sample01", "sample02")):
    header = ["Source Name", "Characteristics[individual]", "Characteristics[sex]",
              "Characteristics[disease]", "Characteristics[organism part]", "Assay Name",
              "Characteristics[genotype]", "Factor Value[genotype]",
              "Derived Array Data File", "Derived Array Data File", "Derived Array Data File"]
    rows = []
    filenames = {"scrna_meta.tsv"}
    for index, sample in enumerate(samples, 1):
        filtered = sample + "_filtered_feature_bc_matrix.tar.gz"
        raw = sample + "_raw_feature_bc_matrix.tar.gz"
        rows.append(["Source " + str(index), str(index), "male", "PSC", "blood", "Assay " + str(index),
                     "SNP", "SNP", raw, filtered, "scrna_meta.tsv"])
        filenames.update({raw, filtered})
    stream = io.StringIO()
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerows([header, *rows])
    study = json.dumps({"files": [{"type": "file", "path": filename} for filename in sorted(filenames)]})
    expected = {sample + "_filtered_feature_bc_matrix.tar.gz" for sample in samples}
    return stream.getvalue(), study, expected


def interfaces():
    return {sample: {"genotype": "SNP", "sex": "male"} for sample in ("sample01", "sample02")}


def cell_text(rows=None):
    if rows is None:
        rows = [["sample01", "sample01_" + BARE1, "SNP", "male"],
                ["sample02", "sample02_" + BARE1, "SNP", "male"]]
    return "sample_name\tbarcode\tcondition\tsex\n" + "".join("\t".join(row) + "\n" for row in rows)


def archive_fixture(path: Path, sample="sample01", fractional=False):
    prefix = sample + "_filtered_feature_bc_matrix"
    matrix = ("%%MatrixMarket matrix coordinate " + ("real" if fractional else "integer") + " general\n"
              '%metadata_json: {"software_version": "cellranger-6.1.1", "format_version": 2}\n'
              "3 2 3\n1 1 " + ("9007199254740993.0000000000001" if fractional else "9007199254740993") + "\n3 1 7\n2 2 5\n")
    payloads = {"features.tsv.gz": b"g1\tBACH2\tGene Expression\ng2\tname\tGene Expression\na1\tADT_CD3\tAntibody Capture\n",
                "barcodes.tsv.gz": (BARE1 + "\n" + BARE2 + "\n").encode(),
                "matrix.mtx.gz": matrix.encode()}
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode="w") as archive:
        directory = tarfile.TarInfo(prefix + "/"); directory.type = tarfile.DIRTYPE
        archive.addfile(directory)
        for name, payload in payloads.items():
            value = gzip.compress(payload)
            member = tarfile.TarInfo(prefix + "/" + name); member.size = len(value)
            archive.addfile(member, io.BytesIO(value))
    path.write_bytes(gzip.compress(raw.getvalue()))
    return {"complete": True, "body_bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "url": URL, "response_url": URL, "started_at_utc": "test", "completed_at_utc": "test"}


class AcquisitionTests(unittest.TestCase):
    def test_frozen_exact_seven_and_total(self):
        self.assertEqual(len(q.EXPECTED), 7)
        self.assertEqual(sum(q.EXPECTED.values()), 298421791)
        self.assertNotIn("sample05_filtered_feature_bc_matrix.tar.gz", q.EXPECTED)

    def test_only_exact_ordinary_public_routes(self):
        q.check_url(URL, NAME)
        q.check_url("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/013/E-MTAB-14013/Files/" + NAME, NAME)
        for invalid in [URL.replace("https:", "http:"), URL + "?download=1", URL + "#x",
                        URL.replace("www.ebi.ac.uk", "other.example"), URL.replace("www.ebi.ac.uk", "www.ebi.ac.uk:9443"),
                        URL.replace("/biostudies/files/E-MTAB-14013/", "/private/"), URL.replace("www.ebi.ac.uk", "user@www.ebi.ac.uk")]:
            with self.subTest(url=invalid), self.assertRaises(ValueError):
                q.check_url(invalid, NAME)

    def test_sample05_and_raw_reads_cannot_be_requested(self):
        for name in ["sample05_filtered_feature_bc_matrix.tar.gz", "sample01_raw_feature_bc_matrix.tar.gz", "sample01.fastq.gz"]:
            with self.assertRaises(ValueError):
                q.check_url("https://www.ebi.ac.uk/biostudies/files/E-MTAB-14013/" + name, name)

    def test_exact_body_and_cache_no_redownload(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(q, "EXPECTED", {NAME: 3}):
            opener = FakeOpener([FakeResponse()])
            item = {"filename": NAME, "expected_bytes": 3, "url": URL}
            with patch.object(q, "build_opener", return_value=opener):
                first = q.acquire_one(item, Path(temp), threading.Event())
                second = q.acquire_one(item, Path(temp), threading.Event())
            self.assertEqual(first, second)
            self.assertTrue(first["complete"])
            self.assertEqual(first["body_bytes"], 3)
            self.assertEqual(len(opener.calls), 1)

    def test_short_body_sets_stop_and_cannot_retry(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(q, "EXPECTED", {NAME: 4}):
            opener = FakeOpener([FakeResponse(b"abc", headers={"Content-Length": "4"})])
            item = {"filename": NAME, "expected_bytes": 4, "url": URL}
            stop = threading.Event()
            with patch.object(q, "build_opener", return_value=opener):
                result = q.acquire_one(item, Path(temp), stop)
                with self.assertRaises(ValueError):
                    q.acquire_one(item, Path(temp), threading.Event())
            self.assertFalse(result["complete"])
            self.assertTrue(stop.is_set())
            self.assertEqual(result["body_bytes"], 3)
            self.assertEqual(len(opener.calls), 1)

    def test_oversized_missing_or_encoded_header_body_never_read(self):
        for headers in [{"Content-Length": "4"}, {}, {"Content-Length": "3", "Content-Encoding": "gzip"}]:
            with self.subTest(headers=headers), tempfile.TemporaryDirectory() as temp, patch.object(q, "EXPECTED", {NAME: 3}):
                response = FakeResponse(headers=headers)
                with patch.object(q, "build_opener", return_value=FakeOpener([response])):
                    result = q.acquire_one({"filename": NAME, "expected_bytes": 3, "url": URL}, Path(temp), threading.Event())
                self.assertFalse(result["complete"])
                self.assertEqual(response.read_calls, 0)
                self.assertEqual(result["body_bytes"], 0)

    def test_foreign_redirect_not_followed_and_error_body_not_read(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(q, "EXPECTED", {NAME: 3}):
            response = FakeResponse(status=302, headers={"Location": "https://other.example/" + NAME})
            opener = FakeOpener([response])
            with patch.object(q, "build_opener", return_value=opener):
                result = q.acquire_one({"filename": NAME, "expected_bytes": 3, "url": URL}, Path(temp), threading.Event())
            self.assertFalse(result["complete"])
            self.assertEqual(len(opener.calls), 1)
            self.assertEqual(response.read_calls, 0)

    def test_existing_unverified_source_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(q, "EXPECTED", {NAME: 3}):
            p = Path(temp) / NAME; p.write_bytes(b"keep")
            with self.assertRaises(ValueError):
                q.acquire_one({"filename": NAME, "expected_bytes": 3, "url": URL}, Path(temp), threading.Event())
            self.assertEqual(p.read_bytes(), b"keep")

    def test_global_stop_prevents_first_request(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(q, "EXPECTED", {NAME: 3}):
            stop = threading.Event(); stop.set()
            opener = FakeOpener([])
            with patch.object(q, "build_opener", return_value=opener):
                result = q.acquire_one({"filename": NAME, "expected_bytes": 3, "url": URL}, Path(temp), stop)
            self.assertFalse(result["complete"])
            self.assertFalse(opener.calls)


class InterfaceTests(unittest.TestCase):
    def test_repeated_sdrf_file_columns_preserved(self):
        sdrf, study, expected = synthetic_sources()
        result = q.source_interfaces(sdrf, study, expected)
        self.assertEqual(set(result), {"sample01", "sample02"})
        self.assertEqual(len(result["sample01"]["files"]), 3)
        self.assertEqual(result["sample02"]["individual"], "2")

    def test_sdrf_inventoried_file_and_genotype_factor_mismatch_rejected(self):
        sdrf, study, expected = synthetic_sources()
        with self.assertRaises(ValueError):
            q.source_interfaces(sdrf, study, {"different_filtered_feature_bc_matrix.tar.gz"})
        with self.assertRaises(ValueError):
            q.source_interfaces(sdrf.replace("SNP\tSNP", "SNP\tnoSNP", 1), study, expected)

    def test_duplicate_accession_local_donor_rejected(self):
        sdrf, study, expected = synthetic_sources()
        with self.assertRaises(ValueError):
            q.source_interfaces(sdrf.replace("Source 2\t2", "Source 2\t1"), study, expected)

    def test_bare_barcode_can_recur_between_samples(self):
        retained, summary = q.fixed_cohort_metadata(cell_text(), interfaces())
        self.assertEqual(retained["sample01"], retained["sample02"])
        self.assertEqual(summary["unique_sample_barcode_pairs"], 2)
        self.assertEqual(summary["bare_barcode_strings_repeated_across_samples"], 1)
        self.assertEqual(summary["RNA_SCT_margin_columns_present"], [])

    def test_duplicate_within_sample_barcode_rejected(self):
        text = cell_text() + "sample01\tsample01_" + BARE1 + "\tSNP\tmale\n"
        with self.assertRaises(ValueError):
            q.fixed_cohort_metadata(text, interfaces())

    def test_wrong_sample_prefix_and_suffix_repair_rejected(self):
        for text in [cell_text().replace("sample01_", "sample02_", 1), cell_text().replace(BARE1, "AAAA-1", 1)]:
            with self.assertRaises(ValueError):
                q.fixed_cohort_metadata(text, interfaces())

    def test_unknown_sample_or_source_condition_mismatch_rejected(self):
        for text in [cell_text().replace("sample01", "sample03"), cell_text().replace("SNP", "noSNP", 1), cell_text().replace("male", "female", 1)]:
            with self.assertRaises(ValueError):
                q.fixed_cohort_metadata(text, interfaces())

    def test_empty_donor_and_duplicate_metadata_header_rejected(self):
        for text in [cell_text().split("sample02\t")[0], cell_text().replace("condition\tsex", "condition\tcondition")]:
            with self.assertRaises(ValueError):
                q.fixed_cohort_metadata(text, interfaces())

    def test_membership_records_missing_and_extra_without_new_retention(self):
        result = q.cell_membership([BARE1, BARE2], {BARE1, "G" * 16 + "-1"})
        self.assertEqual(result["source_cells"], 2)
        self.assertEqual(result["retained_cells_matched"], 1)
        self.assertEqual(result["retained_cells_missing"], 1)
        self.assertEqual(result["source_cells_outside_fixed_retained_metadata"], 1)


class AxisAndAssayTests(unittest.TestCase):
    def setUp(self):
        self.axis = [["g1", "same", "Gene Expression"], ["g2", "same", "Gene Expression"]]

    def test_full_ordered_axis_preserves_duplicate_names(self):
        result = q.compare_rna_axes({"sample01": self.axis, "sample02": self.axis.copy()})
        self.assertTrue(result["all_eight_ordered_RNA_axes_identical"])
        self.assertEqual(result["RNA_feature_universe_rows_if_identical"], 2)
        self.assertFalse(result["assembly_performed"])

    def test_reordered_axis_stops_assembly_even_same_ids(self):
        result = q.compare_rna_axes({"sample01": self.axis, "sample02": list(reversed(self.axis))})
        self.assertFalse(result["all_eight_ordered_RNA_axes_identical"])
        self.assertTrue(result["comparisons"]["sample02"]["same_ID_set_as_reference"])
        self.assertIn("STOP", result["assembly_status"])

    def test_missing_feature_not_intersected_unioned_or_imputed(self):
        result = q.compare_rna_axes({"sample01": self.axis, "sample02": self.axis[:1]})
        self.assertFalse(result["all_eight_ordered_RNA_axes_identical"])
        self.assertEqual(result["comparisons"]["sample02"]["reference_IDs_absent_from_this_sample"], 1)
        self.assertIsNone(result["RNA_feature_universe_rows_if_identical"])

    def test_changed_source_name_not_reconciled(self):
        changed = [["g1", "new-name", "Gene Expression"], self.axis[1]]
        result = q.compare_rna_axes({"sample01": self.axis, "sample02": changed})
        self.assertEqual(result["comparisons"]["sample02"]["shared_IDs_with_name_or_type_difference"], 1)
        self.assertFalse(result["all_eight_ordered_RNA_axes_identical"])

    def test_duplicate_or_absent_rna_id_axis_rejected(self):
        for axis in [[], [self.axis[0], self.axis[0]]]:
            with self.assertRaises(ValueError):
                q.compare_rna_axes({"sample01": axis})

    def test_canonical_axis_hash_is_field_boundary_safe(self):
        self.assertNotEqual(q.axis_sha256([["a\tb", "c", "Gene Expression"]]),
                            q.axis_sha256([["a", "b\tc", "Gene Expression"]]))

    def test_all_and_retained_totals_stay_type_separate_and_exact(self):
        total = q.technical_totals({"Gene Expression": [9007199254740993, 5], "Antibody Capture": [7, 9]}, [BARE1, BARE2], {BARE1})
        self.assertEqual(total["Gene Expression"]["all_source_cells_total_units_exact"], "9007199254740998")
        self.assertEqual(total["Gene Expression"]["fixed_retained_cells_total_units_exact"], "9007199254740993")
        self.assertEqual(total["Antibody Capture"]["all_source_cells_total_units_exact"], "16")
        self.assertNotIn("Combined", total)

    def test_fraction_not_forced_into_integer_and_margin_shape_checked(self):
        result = q.technical_totals({"Gene Expression": [Fraction(1, 3)]}, [BARE1], {BARE1})
        self.assertEqual(result["Gene Expression"]["fixed_retained_cells_total_units_exact"], "1/3")
        with self.assertRaises(ValueError):
            q.technical_totals({"Gene Expression": [1]}, [BARE1, BARE2], {BARE1})

    def test_rna_bach2_does_not_create_adt_measurement(self):
        result = q.antibody_annotations({"sample01": [["g1", "BACH2", "Gene Expression"], ["ADT_CD3", "ADT_CD3", "Antibody Capture"]]})
        self.assertFalse(result["explicit_BACH2_named_Antibody_Capture_feature_present"])
        self.assertEqual(result["distinct_original_ID_name_type_triples"], 1)
        self.assertFalse(result["tag_units_are_absolute_protein_molecules"])

    def test_explicit_bach2_adt_is_only_source_label_evidence(self):
        result = q.antibody_annotations({"sample01": [["ADT_BACH2", "ADT_BACH2", "Antibody Capture"]]})
        self.assertTrue(result["explicit_BACH2_named_Antibody_Capture_feature_present"])
        self.assertFalse(result["antibody_reagent_identity_specificity_independently_reconciled"])
        self.assertFalse(result["RNA_to_protein_translation_effect_tested"])

    def test_opaque_and_similar_labels_not_inferred_as_bach2(self):
        rows = [["TAG001", "unknown", "Antibody Capture"], ["ADT_BACH2B", "ADT_BACH2B", "Antibody Capture"]]
        result = q.antibody_annotations({"sample01": rows, "sample02": rows})
        self.assertFalse(result["explicit_BACH2_named_Antibody_Capture_feature_present"])
        self.assertFalse(result["all_source_names_use_ADT_prefix"])
        self.assertEqual(result["targets_as_published"][0]["present_in_samples"], ["sample01", "sample02"])


class ReusedParserIntegrationTests(unittest.TestCase):
    def test_sample_specific_directory_exact_values_and_cache_replay(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp); archive = directory / NAME
            receipt = archive_fixture(archive)
            result, features, barcodes = q.qualify_one("sample01", archive, receipt, {BARE1}, directory / "work")
            replay, _, _ = q.qualify_one("sample01", archive, receipt, {BARE1}, directory / "work")
            self.assertEqual(result, replay)
            self.assertEqual(len(features), 3)
            self.assertEqual(len(barcodes), 2)
            self.assertEqual(result["technical_totals_by_feature_type"]["Gene Expression"]["all_source_cells_total_units_exact"], "9007199254740998")
            self.assertTrue(result["released_RNA_source_count_format_supported"])
            self.assertFalse(result["original_UMI_count_likelihood_gate_certified"])

    def test_unexpected_internal_sample_or_source_hash_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp); archive = directory / NAME
            receipt = archive_fixture(archive)
            with self.assertRaises(ValueError):
                q.qualify_one("sample02", archive, receipt, {BARE1}, directory / "wrong-sample")
            with self.assertRaises(ValueError):
                q.qualify_one("sample01", archive, dict(receipt, sha256="0" * 64), {BARE1}, directory / "bad-pin")

    def test_fractional_source_not_promoted_to_count_interface(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp); archive = directory / NAME
            receipt = archive_fixture(archive, fractional=True)
            result, _, _ = q.qualify_one("sample01", archive, receipt, {BARE1}, directory / "work")
            self.assertEqual(result["matrix"]["fractional_entries"], 1)
            self.assertFalse(result["released_RNA_source_count_format_supported"])


if __name__ == "__main__":
    unittest.main()
