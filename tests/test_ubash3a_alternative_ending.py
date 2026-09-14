"""Independent small examples for the prospective peptide search dictionary."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_ubash3a_alternative_ending import digest_peptides


class CandidatePeptideTests(unittest.TestCase):
    def test_proline_prevents_lysine_cleavage(self):
        self.assertEqual(list(digest_peptides("AKPRAK", (0,), 1, 30)),
                         [(0,4,0,"AKPR"),(4,6,0,"AK")])

    def test_one_missed_cleavage_preserves_a_terminal_fragment(self):
        self.assertEqual(list(digest_peptides("AKAA", (0,1), 1, 30)),
                         [(0,2,0,"AK"),(0,4,1,"AKAA"),(2,4,0,"AA")])

    def test_length_bounds_are_inclusive(self):
        self.assertEqual(list(digest_peptides("AAAAAAK", (0,), 7, 7)), [(0,7,0,"AAAAAAK")])
        self.assertEqual(list(digest_peptides("AAAAAAK", (0,), 8, 20)), [])


if __name__ == "__main__":
    unittest.main()
