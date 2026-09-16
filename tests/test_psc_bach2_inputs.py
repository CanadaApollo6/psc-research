"""Synthetic tests for source-only BACH2 qualification, with no network or counts."""
from __future__ import annotations

import copy
import csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import qualify_psc_bach2_inputs as q


def tsv(header, rows):
    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return out.getvalue()


def sdrf_fixture():
    header = ["Source Name", "Characteristics[individual]", "Characteristics[sex]",
              "Characteristics[disease]", "Characteristics[organism part]", "Assay Name",
              "Characteristics[genotype]", "Factor Value[genotype]",
              "Derived Array Data File", "Derived Array Data File", "Derived Array Data File"]
    rows = []
    inventory = {}
    for source, donor, sample, genotype in [("Arbitrary A", "x7", "batchB", "SNP"),
                                           ("Arbitrary B", "y2", "batchA", "noSNP")]:
        raw = sample + "_raw_feature_bc_matrix.tar.gz"
        filtered = sample + "_filtered_feature_bc_matrix.tar.gz"
        rows.append([source, donor, "male", "primary sclerosing cholangitis", "blood", source,
                     genotype, genotype, raw, filtered, "scrna_meta.tsv"])
        inventory.update({filename: {"type": "file", "path": filename, "size": 100}
                          for filename in (raw, filtered, "scrna_meta.tsv")})
    return header, rows, inventory


def cell_fixture():
    _, records = q.parse_sdrf(tsv(*sdrf_fixture()[:2]))
    samples = q.sample_interfaces(records, sdrf_fixture()[2])
    header = ["sample_name", "sex", "condition", "seurat_clusters", "paper_clusters",
              "paper_clusters_short", "UMAP_1", "UMAP_2", "barcode"]
    rows = []
    for sample, info in samples.items():
        for i in range(55):
            barcode = "".join("ACGT"[(i >> (2 * j)) & 3] for j in range(16)) + "-1"
            state = "C0" if i < 50 else "C15"
            label = "C0: CD4+ TN RTE" if i < 50 else "C15: MAIT"
            seurat = "0" if i < 50 else "16"
            rows.append([sample, "male", info["genotype"], seurat, label, state,
                         "not parsed", "not parsed", sample + "_" + barcode])
    return header, rows, samples


class ParsingTests(unittest.TestCase):
    def test_repeated_sdrf_columns_preserved(self):
        header, rows, inventory = sdrf_fixture()
        parsed_header, records = q.parse_sdrf(tsv(header, rows))
        self.assertEqual(parsed_header.count("Derived Array Data File"), 3)
        self.assertEqual(len(records[0]["files"]), 3)
        samples = q.sample_interfaces(records, inventory)
        self.assertEqual(samples["batchB"]["individual"], "x7")
        self.assertEqual(samples["batchA"]["genotype"], "noSNP")

    def test_row_order_does_not_define_join(self):
        header, rows, samples = cell_fixture()
        original, _ = q.summarize_cells(tsv(header, rows), samples)
        reversed_result, _ = q.summarize_cells(tsv(header, rows[::-1]), dict(reversed(list(samples.items()))))
        self.assertEqual(original, reversed_result)

    def test_all_state_counts_and_barcode_collisions(self):
        header, rows, samples = cell_fixture()
        result, private = q.summarize_cells(tsv(header, rows), samples)
        self.assertEqual(result["retained_cells"], 110)
        self.assertEqual(result["barcode_join"]["bare_barcodes_shared_across_samples"], 55)
        self.assertTrue(result["barcode_join"]["sample_plus_stripped_barcode_unique"])
        self.assertFalse(result["barcode_join"]["matrix_barcode_membership_verified"])
        self.assertEqual(result["states"][0]["donors_at_cell_count_threshold"]["50"], {"SNP": 1, "noSNP": 1})
        self.assertEqual(result["states"][1]["seurat_clusters"], ["16"])
        self.assertEqual(result["states"][1]["donors_at_cell_count_threshold"]["20"], {"SNP": 0, "noSNP": 0})
        self.assertIn("sample_to_source", private)
        self.assertNotIn("sample_to_source", result)
        self.assertFalse(result["umap_values_parsed_or_interpreted"])

    def test_bad_tsv_shape_rejected(self):
        with self.assertRaises(ValueError):
            q.read_tsv("a\tb\n1\n")

    def test_duplicated_donor_rejected(self):
        header, rows, _ = sdrf_fixture()
        rows[1][1] = rows[0][1]
        with self.assertRaisesRegex(ValueError, "Nonunique individual"):
            q.parse_sdrf(tsv(header, rows))

    def test_conflicting_genotype_factor_rejected(self):
        header, rows, _ = sdrf_fixture()
        rows[0][7] = "noSNP"
        with self.assertRaises(ValueError):
            q.parse_sdrf(tsv(header, rows))

    def test_undeclared_file_rejected(self):
        header, rows, inventory = sdrf_fixture()
        _, records = q.parse_sdrf(tsv(header, rows))
        inventory.pop(rows[0][-2])
        with self.assertRaises(ValueError):
            q.sample_interfaces(records, inventory)

    def test_conflicting_cell_annotation_rejected(self):
        header, rows, samples = cell_fixture()
        rows[0][2] = "noSNP"
        with self.assertRaises(ValueError):
            q.summarize_cells(tsv(header, rows), samples)

    def test_unknown_sample_rejected(self):
        header, rows, samples = cell_fixture()
        rows[0][0] = "expression_guessed_donor"
        with self.assertRaises(ValueError):
            q.summarize_cells(tsv(header, rows), samples)

    def test_duplicate_barcode_rejected(self):
        header, rows, samples = cell_fixture()
        rows.append(rows[0])
        with self.assertRaises(ValueError):
            q.summarize_cells(tsv(header, rows), samples)

    def test_missing_prefix_rejected(self):
        header, rows, samples = cell_fixture()
        rows[0][-1] = rows[0][-1].split("_", 1)[1]
        with self.assertRaises(ValueError):
            q.summarize_cells(tsv(header, rows), samples)

    def test_crosswalk_ambiguity_rejected(self):
        header, rows, samples = cell_fixture()
        rows[0][3] = "16"
        with self.assertRaisesRegex(ValueError, "not one-to-one"):
            q.summarize_cells(tsv(header, rows), samples)

    def test_absent_sample_rejected(self):
        header, rows, samples = cell_fixture()
        rows = [r for r in rows if r[0] == "batchA"]
        with self.assertRaises(ValueError):
            q.summarize_cells(tsv(header, rows), samples)

    def test_missing_required_state_rejected(self):
        header, rows, samples = cell_fixture()
        rows[0][4] = ""
        with self.assertRaises(ValueError):
            q.summarize_cells(tsv(header, rows), samples)

    def test_liver_no_genotype_or_cross_accession_inference(self):
        records = [{"individual": "1", "sex": "male", "age": "18", "age_unit": "year",
                    "genotype": None, "files": ["count_table.txt"]},
                   {"individual": "2", "sex": "female", "age": "59", "age_unit": "year",
                    "genotype": None, "files": ["count_table.txt"]}]
        result = q.summarize_liver(records)
        self.assertEqual(result["genotype_nonmissing"], 0)
        self.assertEqual(result["age_min"], 18)
        self.assertFalse(result["cross_accession_participant_join_established"])
        self.assertFalse(result["genetic_replication_eligible"])

    def test_duplicate_inventory_path_rejected(self):
        f = {"path": "x", "type": "file", "size": 1}
        with self.assertRaises(ValueError):
            q.file_inventory({"files": [[f, f]]})


class ClinicalTests(unittest.TestCase):
    def fixture(self):
        lines = ["Table S2: Clinical parameters of people with PSC, related to Figure 2A-J."]
        for i in range(8):
            group = "homozygous" if i < 4 else "non-carrier"
            ibd = "PSC-assoc. colitis" if i == 0 else "/"
            therapy = "Vedolizumab" if i == 1 else "/"
            values = [f"patient #{i+1}", group, "cis", str(1960+i), ibd, "/", therapy,
                      "5,6", "6,0", "negative", "negative", "negative"]
            lines.append("  ".join(values))
        return "\n".join(lines)

    def test_clinical_slashes_and_source_decimal_strings_preserved(self):
        summary, rows = q.parse_clinical_s2(self.fixture())
        group = summary["groups_as_published"]["homozygous"]
        self.assertEqual(group["IBD status"], {"/": 3, "PSC-assoc. colitis": 1})
        self.assertEqual(group["immunosuppr. therapy"], {"/": 3, "Vedolizumab": 1})
        self.assertEqual(group["transient elastography [kPa]"]["minimum_as_published"], "5,6")
        self.assertEqual(rows[0]["patient #"], "patient #1")
        self.assertNotIn("individual", rows[0])
        self.assertIn("not merged", summary["clinical_patient_to_accession_individual_join"])
        self.assertTrue(summary["birth_year_is_not_age_at_sampling"])

    def test_clinical_incomplete_group_rejected(self):
        text = self.fixture().rsplit("\n", 1)[0]
        with self.assertRaises(ValueError):
            q.parse_clinical_s2(text)

    def test_html_supplement_is_excluded_despite_http200(self):
        entry = {"complete": True, "method": "GET", "kind": "source_supplement"}
        self.assertEqual(q.response_evidence_status(entry, b"<html>challenge</html>"),
                         "not_requested_supplement_format_excluded")
        self.assertEqual(q.response_evidence_status(entry, b"%PDF-1.4\nstuff\n%%EOF\n"),
                         "complete_pdf_source")


class FakeResponse:
    def __init__(self, data=b"", status=200, headers=None, url="https://www.ebi.ac.uk/test"):
        self.data = data
        self.status = status
        self.headers = headers or {}
        self.url = url
        self.read_calls = 0
        self.closed = False

    def geturl(self):
        return self.url

    def read(self, size):
        self.read_calls += 1
        part, self.data = self.data[:size], self.data[size:]
        return part

    def close(self):
        self.closed = True


class FakeOpener:
    def __init__(self, responses):
        self.responses = list(responses)
        self.requests = []

    def open(self, request, timeout):
        self.requests.append(request)
        return self.responses.pop(0)


class AcquisitionTests(unittest.TestCase):
    def get_item(self, **kwargs):
        result = {"name": "metadata", "method": "GET", "kind": "sample_metadata",
                  "url": "https://www.ebi.ac.uk/test.tsv", "max_bytes": 100}
        result.update(kwargs)
        return result

    def run_fake(self, batch, responses):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            opener = FakeOpener(responses)
            with patch.object(q, "build_opener", return_value=opener), patch("builtins.print"):
                q.acquire(batch, directory)
            return q.read_ledger(directory / "request-ledger.jsonl"), opener

    def test_count_get_even_mislabeled_is_rejected(self):
        with self.assertRaises(ValueError):
            q.validate_request(self.get_item(url="https://ftp.ebi.ac.uk/sample_filtered_feature_bc_matrix.tar.gz"))
        with self.assertRaises(ValueError):
            q.validate_request(self.get_item(kind="count_inventory"))

    def test_unapproved_or_credentialed_endpoint_rejected(self):
        for url in ("https://private.example/x", "http://www.ebi.ac.uk/x", "https://user:password@www.ebi.ac.uk/x"):
            with self.assertRaises(ValueError):
                q.validate_url(url)

    def test_head_never_reads_body(self):
        response = FakeResponse(b"MUST NOT READ", headers={"Content-Length": "325315169"})
        records, _ = self.run_fake([self.get_item(method="HEAD", kind="count_inventory")], [response])
        self.assertEqual(response.read_calls, 0)
        self.assertEqual(records[0]["body_bytes"], 0)
        self.assertTrue(records[0]["complete"])

    def test_safe_public_redirect(self):
        first = FakeResponse(status=302, headers={"Location": "https://ftp.ebi.ac.uk/test.tsv"})
        second = FakeResponse(b"abc", headers={"Content-Length": "3"})
        records, opener = self.run_fake([self.get_item(expected_bytes=3)], [first, second])
        self.assertEqual(len(opener.requests), 2)
        self.assertEqual(first.read_calls, 0)
        self.assertTrue(records[0]["complete"])
        self.assertEqual(records[0]["sha256"], q.sha256(b"abc"))

    def test_redirect_cannot_fetch_count_body(self):
        first = FakeResponse(status=302, headers={"Location": "https://ftp.ebi.ac.uk/count_table.txt"})
        records, opener = self.run_fake([self.get_item()], [first])
        self.assertEqual(len(opener.requests), 1)
        self.assertFalse(records[0]["complete"])
        self.assertEqual(records[0]["body_bytes"], 0)

    def test_declared_over_limit_not_read(self):
        response = FakeResponse(b"secret", headers={"Content-Length": "101"})
        records, _ = self.run_fake([self.get_item()], [response])
        self.assertEqual(response.read_calls, 0)
        self.assertFalse(records[0]["complete"])

    def test_unknown_length_capped_and_preserved_incomplete(self):
        response = FakeResponse(b"a" * 200)
        records, _ = self.run_fake([self.get_item()], [response])
        self.assertEqual(records[0]["body_bytes"], 100)
        self.assertFalse(records[0]["complete"])

    def test_length_disagreement_keeps_failed_bytes(self):
        response = FakeResponse(b"abc", headers={"Content-Length": "3"})
        records, _ = self.run_fake([self.get_item(expected_bytes=4)], [response])
        self.assertEqual(records[0]["body_bytes"], 3)
        self.assertFalse(records[0]["complete"])

    def test_budget_is_cumulative(self):
        first = FakeResponse(b"a" * 9, headers={"Content-Length": "9"})
        second = FakeResponse(b"bc", headers={"Content-Length": "2"})
        with patch.object(q, "BUDGET", 10):
            records, _ = self.run_fake([self.get_item(name="a", max_bytes=9), self.get_item(name="b", max_bytes=9)], [first, second])
        self.assertEqual(sum(x["body_bytes"] for x in records), 9)
        self.assertEqual(second.read_calls, 0)
        self.assertFalse(records[-1]["complete"])

    def test_non200_metadata_body_not_read(self):
        response = FakeResponse(b"error", status=404)
        records, _ = self.run_fake([self.get_item()], [response])
        self.assertEqual(response.read_calls, 0)
        self.assertFalse(records[0]["complete"])

    def test_corrupted_cache_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            (directory / "x.body").write_bytes(b"bad")
            entry = {"name": "x", "complete": True, "method": "GET", "body_file": "x.body",
                     "body_bytes": 3, "sha256": q.sha256(b"yes")}
            with self.assertRaises(ValueError):
                q.verified_source("x", [entry], directory)


if __name__ == "__main__":
    unittest.main()
