"""Synthetic source/sparse/security checks; no real expression effects or network."""
from __future__ import annotations

import gzip
import io
import json
import sys
import tarfile
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import qualify_psc_bach2_counts as q


def make_tar(path, members):
    with tarfile.open(path, "w") as out:
        for name, data, kind in members:
            info = tarfile.TarInfo(name)
            info.type = kind
            if kind == tarfile.REGTYPE:
                info.size = len(data)
                out.addfile(info, io.BytesIO(data))
            else:
                if kind == tarfile.SYMTYPE:
                    info.linkname = "outside"
                out.addfile(info)


def fixture_members():
    values = {"features.tsv.gz": b"ENSG000001\tSame\tGene Expression\nENSG000002\tSame\tGene Expression\nAB1\tCD4\tAntibody Capture\n",
              "barcodes.tsv.gz": b"AAAAAAAAAAAAAAAA-1\nCCCCCCCCCCCCCCCC-1\n",
              "matrix.mtx.gz": b'%%MatrixMarket matrix coordinate integer general\n%metadata_json: {"software_version": "cellranger-6.1.1"}\n3 2 4\n1 1 9007199254740993\n3 1 4\n2 2 5\n3 2 6\n'}
    return [("sample05_filtered_feature_bc_matrix", b"", tarfile.DIRTYPE)] + [
        ("sample05_filtered_feature_bc_matrix/" + name, gzip.compress(data), tarfile.REGTYPE)
        for name, data in values.items()]


class ArchiveTests(unittest.TestCase):
    def test_safe_paths_and_canonicalization(self):
        self.assertEqual(q.canonical_member_name("./x/y"), "x/y")
        for bad in ("../x", "/x", "C:/x", "x\\y", "x/../y", "", "."):
            with self.assertRaises(ValueError):
                q.canonical_member_name(bad)

    def test_alias_collision_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "a.tar"
            make_tar(p, [("x", b"a", tarfile.REGTYPE), ("./x", b"b", tarfile.REGTYPE)])
            with self.assertRaises(ValueError):
                q.inventory_tar(p)

    def test_symlink_and_fifo_rejected(self):
        for kind in (tarfile.SYMTYPE, tarfile.FIFOTYPE, tarfile.LNKTYPE):
            with tempfile.TemporaryDirectory() as temp:
                p = Path(temp) / "a.tar"
                make_tar(p, [("x", b"", kind)])
                with self.assertRaises(ValueError):
                    q.inventory_tar(p)

    def test_path_ancestor_file_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "a.tar"
            make_tar(p, [("x", b"a", tarfile.REGTYPE), ("x/y", b"b", tarfile.REGTYPE)])
            with self.assertRaises(ValueError):
                q.inventory_tar(p)

    def test_member_count_and_size_caps(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "a.tar"
            make_tar(p, [("x", b"ab", tarfile.REGTYPE)])
            with patch.dict(q.CAPS, {"tar_regular_member_bytes": 1}):
                with self.assertRaises(ValueError):
                    q.inventory_tar(p)
            with patch.dict(q.CAPS, {"tar_members": 0}):
                with self.assertRaises(ValueError):
                    q.inventory_tar(p)

    def test_tar_trailer_required_and_nonzero_appendage_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "a.tar"
            make_tar(p, [("x", b"a", tarfile.REGTYPE)])
            original = p.read_bytes()
            p.write_bytes(original[:1024])
            with self.assertRaises(ValueError):
                q.inventory_tar(p)
            p.write_bytes(original + b"x" * 512)
            with self.assertRaises(ValueError):
                q.inventory_tar(p)

    def test_tar_header_checksum_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "a.tar"
            make_tar(p, [("x", b"a", tarfile.REGTYPE)])
            data = bytearray(p.read_bytes()); data[0] ^= 1; p.write_bytes(data)
            with self.assertRaises(tarfile.ReadError):
                q.inventory_tar(p)

    def test_outer_gzip_crc_and_short_eof(self):
        for mode in ("crc", "short"):
            with tempfile.TemporaryDirectory() as temp:
                p = Path(temp) / "a.gz"
                data = bytearray(gzip.compress(b"hello"))
                if mode == "crc":
                    data[-8] ^= 1
                else:
                    data = data[:-3]
                p.write_bytes(data)
                with self.assertRaises((gzip.BadGzipFile, EOFError)):
                    q.decompress_outer(p, Path(temp) / "work")

    def test_outer_uncompressed_cap(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "a.gz"; p.write_bytes(gzip.compress(b"a" * 101))
            with patch.dict(q.CAPS, {"outer_uncompressed_bytes": 100}):
                with self.assertRaises(ValueError):
                    q.decompress_outer(p, Path(temp) / "work")

    def test_complete_nested_integrity_and_schema(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp); p = directory / "a.tar"
            make_tar(p, fixture_members())
            inventory = q.inventory_tar(p)
            self.assertEqual(len(inventory), 4)
            prepared = q.prepare_members(p, inventory, directory / "work")
            self.assertEqual(set(prepared), {"features", "barcodes", "matrix"})
            self.assertTrue(all(x["nested_gzip_CRC_ISIZE_and_EOF_passed"] for x in prepared.values()))
            self.assertEqual(prepared, q.prepare_members(p, inventory, directory / "work"))
            features, fs = q.parse_features(Path(prepared["features"]["decoded_path"]))
            barcodes = q.parse_barcodes(Path(prepared["barcodes"]["decoded_path"]))
            matrix, sums = q.parse_matrix(Path(prepared["matrix"]["decoded_path"]), features, barcodes)
            self.assertEqual(fs["duplicate_feature_name_rows"], 1)
            self.assertEqual(matrix["sum_exact_by_feature_type"]["Gene Expression"], "9007199254740998")
            self.assertEqual(sums["Antibody Capture"], [4, 6])

    def test_corrupt_nested_crc_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "a.tar"
            members = fixture_members()
            name, data, kind = members[-1]
            bad = bytearray(data); bad[-8] ^= 1
            members[-1] = (name, bytes(bad), kind)
            make_tar(p, members)
            with self.assertRaises(gzip.BadGzipFile):
                q.prepare_members(p, q.inventory_tar(p), Path(temp) / "work")

    def test_bounded_copy_never_overwrites(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "member"
            p.write_bytes(b"preserve")
            with self.assertRaises(FileExistsError):
                q.bounded_copy(io.BytesIO(b"other"), p, 10)
            self.assertEqual(p.read_bytes(), b"preserve")


class NumericAndJoinTests(unittest.TestCase):
    def matrix(self, text):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "matrix"; p.write_text(text)
            return q.parse_matrix(p, [["g1", "a", "Gene Expression"], ["g2", "b", "Antibody Capture"]], ["barcode"])

    def test_exact_large_integer_and_fraction(self):
        self.assertEqual(q.exact_value("9007199254740993", "integer"), 9007199254740993)
        value = q.exact_value("9007199254740993.0000000000000000001", "real")
        self.assertEqual(value, Fraction(9007199254740993) + Fraction(1, 10**19))
        self.assertNotIsInstance(value, int)
        self.assertEqual(q.exact_value("1e2", "real"), 100)

    def test_invalid_numeric_tokens(self):
        for token in ("NaN", "inf", "-1", "1/2", "1e9999"):
            with self.assertRaises(ValueError):
                q.exact_value(token, "real")
        with self.assertRaises(ValueError):
            q.exact_value("1.0", "integer")

    def test_fractional_values_not_rounded_or_forced_counts(self):
        result, sums = self.matrix("%%MatrixMarket matrix coordinate real general\n2 1 2\n1 1 9007199254740993.0000000000000000001\n2 1 0.1\n")
        self.assertEqual(result["fractional_entries"], 2)
        self.assertEqual(sums["Antibody Capture"], [Fraction(1, 10)])

    def test_unsorted_unique_coordinates_accepted(self):
        result, _ = self.matrix("%%MatrixMarket matrix coordinate integer general\n2 1 2\n2 1 1\n1 1 2\n")
        self.assertFalse(result["strict_column_major_order"])
        self.assertEqual(result["duplicate_coordinates"], 0)

    def test_duplicate_coordinates_rejected_not_summed(self):
        with self.assertRaises(ValueError):
            self.matrix("%%MatrixMarket matrix coordinate integer general\n2 1 2\n1 1 1\n1 1 2\n")

    def test_out_of_bounds_and_truncated_coordinates(self):
        for body in ("3 1 1\n", "1 1 1\n"):
            with self.assertRaises(ValueError):
                self.matrix("%%MatrixMarket matrix coordinate integer general\n2 1 2\n" + body)

    def test_axis_dimensions_and_extra_entries_rejected(self):
        for text in ("2 2 0\n", "2 1 0\n1 1 1\n"):
            with self.assertRaises(ValueError):
                self.matrix("%%MatrixMarket matrix coordinate integer general\n" + text)

    def test_duplicate_feature_id_rejected_but_names_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "features"
            p.write_text("g1\tx\tGene Expression\ng1\ty\tGene Expression\n")
            with self.assertRaises(ValueError):
                q.parse_features(p)

    def test_fixed_metadata_subset_and_no_source_margins(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "metadata"
            p.write_text("sample_name\tbarcode\tpaper_clusters\nsample05\tsample05_AAAAAAAAAAAAAAAA-1\tC0\nsample01\tsample01_AAAAAAAAAAAAAAAA-1\tC1\n")
            result, retained = q.fixed_metadata_join(["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1"], p)
            self.assertEqual(result["matched_retained_barcodes"], 1)
            self.assertEqual(result["matrix_barcodes_outside_fixed_retained_metadata"], 1)
            self.assertEqual(result["RNA_SCT_count_margin_columns_present"], [])
            self.assertEqual(retained, {"AAAAAAAAAAAAAAAA-1"})

    def test_wrong_sample_barcode_prefix_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "metadata"
            p.write_text("sample_name\tbarcode\nsample05\tsample01_AAAAAAAAAAAAAAAA-1\n")
            with self.assertRaises(ValueError):
                q.fixed_metadata_join(["AAAAAAAAAAAAAAAA-1"], p)


class FakeResponse:
    def __init__(self, data=b"", status=200, headers=None):
        self.data = data; self.status = status; self.headers = headers or {}; self.read_calls = 0
    def geturl(self):
        return q.URL
    def close(self):
        pass
    def read(self, n):
        self.read_calls += 1
        result, self.data = self.data[:n], self.data[n:]
        return result


class FakeOpener:
    def __init__(self, responses):
        self.responses = list(responses); self.requests = []
    def open(self, request, timeout):
        self.requests.append(request)
        return self.responses.pop(0)


class AcquisitionTests(unittest.TestCase):
    def test_only_one_exact_archive_allowed(self):
        q.check_url(q.URL)
        for bad in (q.URL.replace("sample05", "sample01"), q.URL.replace("filtered", "raw"),
                    q.URL.replace("https:", "http:"), q.URL.replace("www.ebi.ac.uk", "private.example"), q.URL + "?token=x"):
            with self.assertRaises(ValueError):
                q.check_url(bad)

    def test_exact_body_then_offline_cached_no_second_fetch(self):
        with tempfile.TemporaryDirectory() as temp:
            response = FakeResponse(b"abcde", headers={"Content-Length": "5"})
            opener = FakeOpener([response])
            with patch.object(q, "build_opener", return_value=opener), patch.object(q, "EXPECTED_BYTES", 5):
                result = q.acquire(Path(temp))
                repeated = q.acquire(Path(temp))
            self.assertEqual(result, repeated)
            self.assertEqual(len(opener.requests), 1)
            self.assertEqual(result["body_bytes"], 5)

    def test_bad_length_blocks_body(self):
        with tempfile.TemporaryDirectory() as temp:
            response = FakeResponse(b"abcde", headers={"Content-Length": "6"})
            with patch.object(q, "build_opener", return_value=FakeOpener([response])), patch.object(q, "EXPECTED_BYTES", 5):
                with self.assertRaises(ValueError):
                    q.acquire(Path(temp))
            self.assertEqual(response.read_calls, 0)

    def test_unauthorized_redirect_does_not_fetch_second_archive(self):
        with tempfile.TemporaryDirectory() as temp:
            response = FakeResponse(status=302, headers={"Location": q.URL.replace("sample05", "sample01")})
            opener = FakeOpener([response])
            with patch.object(q, "build_opener", return_value=opener):
                with self.assertRaises(ValueError):
                    q.acquire(Path(temp))
            self.assertEqual(len(opener.requests), 1)
            self.assertEqual(response.read_calls, 0)

    def test_short_body_receipt_retained_and_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            response = FakeResponse(b"abc", headers={"Content-Length": "5"})
            opener = FakeOpener([response])
            with patch.object(q, "build_opener", return_value=opener), patch.object(q, "EXPECTED_BYTES", 5):
                with self.assertRaises(ValueError):
                    q.acquire(Path(temp))
                with self.assertRaises(ValueError):
                    q.acquire(Path(temp))
            receipt = json.loads((Path(temp) / "acquisition-receipt.json").read_text())
            self.assertEqual(receipt["body_bytes"], 3)
            self.assertFalse(receipt["complete"])
            self.assertEqual(len(opener.requests), 1)


class AdditionalSafetyTests(unittest.TestCase):
    def test_cached_member_source_lineage_not_just_local_hash(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp); source = directory / "source.tar"
            members = fixture_members(); make_tar(source, members)
            q.prepare_members(source, q.inventory_tar(source), directory / "work")
            name, data, kind = members[1]
            members[1] = (name, gzip.compress(b"CHANGED\n"), kind)
            make_tar(source, members)
            with self.assertRaises(ValueError):
                q.prepare_members(source, q.inventory_tar(source), directory / "work")

    def test_duplicate_exporter_comments_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "matrix"
            p.write_text('%%MatrixMarket matrix coordinate integer general\n%metadata_json: {"v":1}\n%metadata_json: {"v":2}\n1 1 1\n1 1 1\n')
            with self.assertRaises(ValueError):
                q.parse_matrix(p, [["g", "gene", "Gene Expression"]], ["barcode"])

    def test_definition_metadata_cap_is_cumulative(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            first = FakeResponse(b"a" * 800000, headers={"Content-Length": "800000"})
            second = FakeResponse(b"b" * 300000, headers={"Content-Length": "300000"})
            opener = FakeOpener([first, second])
            with patch.object(q, "build_opener", return_value=opener):
                r1 = q.acquire_definition(directory)
                r2 = q.acquire_definition(directory, name="count-definition-current")
            self.assertEqual(r1["body_bytes"] + r2["body_bytes"], 800000)
            self.assertEqual(second.read_calls, 0)
            self.assertFalse(r2["complete_response"])


if __name__ == "__main__":
    unittest.main()
