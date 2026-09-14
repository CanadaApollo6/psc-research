"""Integrity failures must preserve historical source caches."""

import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("read_cache", Path(__file__).resolve().parents[1] / "scripts/fetch_ubash3a_upf1_reads.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ReadCacheIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.cache = self.root / "raw"
        self.cache.mkdir()
        self.sample = {"sample_id": "sample", "bam_url": "https://example.invalid/sample.bam"}
        self.file = self.cache / "sample.UBASH3A.bam"
        self.file.write_bytes(b"original")
        self.pin = {"bytes": 8, "sha256": hashlib.sha256(b"original").hexdigest()}
        self.old = {"local_files": {"raw/sample.UBASH3A.bam": self.pin},
                    "source_metadata_before": {"bytes": 100, "etag": "old", "last_modified": "then"}}

    def tearDown(self):
        self.temp.cleanup()

    def test_modified_cache_is_rejected_without_network_or_overwrite(self):
        self.file.write_bytes(b"modified")
        with patch.object(m, "ROOT", self.root), patch.object(m, "metadata") as net:
            with self.assertRaisesRegex(ValueError, "Changed locus cache"):
                m.fetch(self.sample, {}, [1, 10], old=self.old, cache=self.cache)
            net.assert_not_called()
        self.assertEqual(self.file.read_bytes(), b"modified")

    def test_partial_cache_offline_is_missing_not_a_new_zero_result(self):
        self.old["local_files"]["raw/missing.header.sam"] = self.pin
        with patch.object(m, "ROOT", self.root), patch.object(m, "metadata") as net:
            with self.assertRaises(FileNotFoundError):
                m.fetch(self.sample, {}, [1, 10], old=self.old, offline=True, cache=self.cache)
            net.assert_not_called()
        self.assertEqual(self.file.read_bytes(), b"original")

    def test_changed_remote_identity_stops_partial_cache_restoration(self):
        self.old["local_files"]["raw/missing.header.sam"] = self.pin
        index = self.cache / "source.bai"
        index.write_bytes(b"index")
        source = {"local_cache_path": "raw/source.bai", "sha256": hashlib.sha256(b"index").hexdigest()}
        with patch.object(m, "ROOT", self.root), patch.object(m, "metadata", return_value={"bytes": 100, "etag": "changed", "last_modified": "then"}), patch.object(m.requests, "get") as get:
            with self.assertRaisesRegex(ValueError, "Remote BAM identity differs"):
                m.fetch(self.sample, source, [1, 10], old=self.old, cache=self.cache)
            get.assert_not_called()
        self.assertEqual(self.file.read_bytes(), b"original")
        self.assertFalse((self.cache / "missing.header.sam").exists())


if __name__ == "__main__":
    unittest.main()
