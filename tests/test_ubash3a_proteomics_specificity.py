import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("proteomics_specificity", Path(__file__).resolve().parents[1] / "scripts/analyze_ubash3a_proteomics_specificity.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ProteomicsSpecificityTests(unittest.TestCase):
    def test_overlapping_hits_are_not_lost(self):
        self.assertEqual(list(MODULE.match_positions("AAAAA", "AAA")), [0, 1, 2])

    def test_il_equivalence_is_separate_from_literal_match(self):
        self.assertEqual(list(MODULE.match_positions("KLLHLESGIRK", "IHLESGIR")), [])
        self.assertEqual(list(MODULE.match_positions("KLLHLESGIRK", "IHLESGIR", True)), [2])

    def test_proline_blocks_tryptic_cleavage(self):
        self.assertFalse(MODULE.boundary("AKPA", 2))
        self.assertTrue(MODULE.boundary("AKAA", 2))

    def test_c_terminal_extension_and_nested_candidate(self):
        short = MODULE.match_record("KIHLESGIRKP", "IHLESGIR", 1)
        long = MODULE.match_record("KIHLESGIRKP", "IHLESGIRKP", 1)
        self.assertEqual(short["following_residue"], "K")
        self.assertEqual(short["internal_tryptic_cleavages"], 0)
        self.assertEqual(long["internal_tryptic_cleavages"], 1)
        self.assertTrue(short["fits_frozen_candidate_digest"])
        self.assertTrue(long["fits_frozen_candidate_digest"])

    def test_nontyptic_match_is_preserved_but_not_qualified(self):
        record = MODULE.match_record("AAIHLESGIRAA", "IHLESGIR", 2)
        self.assertFalse(record["n_tryptic"])
        self.assertFalse(record["fits_frozen_candidate_digest"])
        self.assertEqual(record["start_1based"], 3)

    def test_full_parent_sequence_is_not_truncated(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.fa"
            path.write_text(">one preserved header\r\nACDX\r\nIHLESGIR\r\n>two\nABC*\n")
            self.assertEqual(list(MODULE.fasta(path)), [("one preserved header", "ACDXIHLESGIR"), ("two", "ABC*")])
            path.write_text("ACD\n>broken\nACD\n")
            with self.assertRaises(ValueError):
                list(MODULE.fasta(path))

    def test_uppercase_raw_and_escaped_settings(self):
        prefix = "<search_summary>Search name: Example\nFile Name(s):\n C:\\data\\CD4.RAW\n D:\\data\\a.raw\nEnzyme: A &amp; B\n</search_summary>"
        group, names, text = MODULE.extract_settings(prefix)
        self.assertEqual(group, "Example")
        self.assertEqual(names, ["CD4.RAW", "a.raw"])
        self.assertIn("A & B", text)
        with self.assertRaises(ValueError):
            MODULE.extract_settings(prefix[:-20])


if __name__ == "__main__":
    unittest.main()
