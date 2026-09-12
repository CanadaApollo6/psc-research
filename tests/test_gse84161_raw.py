"""Focused synthetic tests for the metadata-only GSE84161 raw audit."""

from __future__ import annotations

import gzip
import io
import struct
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import audit_gse84161_raw as audit  # noqa: E402


def _cel_bytes(
    *,
    magic: int = audit.EXPECTED_CEL_MAGIC,
    version: int = audit.EXPECTED_CEL_VERSION,
    chip: str = audit.EXPECTED_CHIP_HEADER,
) -> bytes:
    header = (
        f"Cols={audit.EXPECTED_CEL_COLUMNS}\n"
        f"Rows={audit.EXPECTED_CEL_ROWS}\n"
        f"TotalX={audit.EXPECTED_CEL_COLUMNS}\n"
        f"TotalY={audit.EXPECTED_CEL_ROWS}\n"
        f"DatHeader={chip}.1sq\n"
    ).encode("latin-1")
    fixed = struct.pack(
        "<6I",
        magic,
        version,
        audit.EXPECTED_CEL_COLUMNS,
        audit.EXPECTED_CEL_ROWS,
        audit.EXPECTED_CEL_CELLS,
        len(header),
    )
    return fixed + header + b"\0" * 32


def _gzip_bytes(payload: bytes) -> bytes:
    output = io.BytesIO()
    with gzip.GzipFile(fileobj=output, mode="wb") as stream:
        stream.write(payload)
    return output.getvalue()


def _write_tar(path: Path, members: dict[str, bytes]) -> None:
    with tarfile.open(path, mode="w") as archive:
        for name, payload in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))


class TestGSE84161Raw(unittest.TestCase):
    def _audit(self, archive_path: Path, names: list[str]) -> dict[str, object]:
        _, digest = audit.sha256_file(archive_path)
        return audit.audit_archive(
            archive_path=archive_path,
            expected_archive_sha256=digest,
            expected_members=names,
            member_order=names,
        )

    def test_valid_members_reach_gzip_eof_and_expected_cel_header(self) -> None:
        names = [
            "GSMTEST1_1_HGU133P.CEL.gz",
            "GSMTEST2_2_HGU133P.CEL.gz",
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "valid.tar"
            _write_tar(path, {name: _gzip_bytes(_cel_bytes()) for name in names})
            result = self._audit(path, names)

        self.assertEqual(result["archive_member_count"], 2)
        self.assertTrue(result["all_gzip_eof_crc_checks_pass"])
        self.assertTrue(result["all_CEL_v4_headers_match_expected_array"])
        self.assertEqual(result["dimensions"], [[1164, 1164, 1354896]])
        self.assertEqual(len(result["members"]), 2)
        self.assertEqual(result["members"][0]["cel_magic"], "64")
        self.assertEqual(result["members"][0]["chip_header_match"], "HG-U133_Plus_2")

    def test_archive_hash_is_checked_before_member_read(self) -> None:
        name = "GSMTEST1_1_HGU133P.CEL.gz"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "hash.tar"
            _write_tar(path, {name: _gzip_bytes(_cel_bytes())})
            with self.assertRaises(audit.ArchiveHashMismatch):
                audit.audit_archive(path, "0" * 64, [name])

    def test_corrupted_gzip_crc_is_rejected(self) -> None:
        name = "GSMTEST1_1_HGU133P.CEL.gz"
        corrupted = bytearray(_gzip_bytes(_cel_bytes()))
        corrupted[-1] ^= 0x01
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "crc.tar"
            _write_tar(path, {name: bytes(corrupted)})
            with self.assertRaises(audit.GzipIntegrityError):
                self._audit(path, [name])

    def test_wrong_cel_header_is_rejected(self) -> None:
        name = "GSMTEST1_1_HGU133P.CEL.gz"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "header.tar"
            _write_tar(path, {name: _gzip_bytes(_cel_bytes(version=3))})
            with self.assertRaises(audit.CelHeaderError):
                self._audit(path, [name])

    def test_missing_expected_member_is_rejected(self) -> None:
        names = [
            "GSMTEST1_1_HGU133P.CEL.gz",
            "GSMTEST2_2_HGU133P.CEL.gz",
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "missing.tar"
            _write_tar(path, {names[0]: _gzip_bytes(_cel_bytes())})
            with self.assertRaises(audit.MemberSetError):
                self._audit(path, names)


if __name__ == "__main__":
    unittest.main()
