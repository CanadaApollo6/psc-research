"""Synthetic CIGAR and quality truths for the fixed long-read locus test."""

import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location("upf1_reads", Path(__file__).resolve().parents[1] / "scripts/analyze_ubash3a_upf1_reads.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ReadStructureTests(unittest.TestCase):
    def test_splice_boundaries_are_one_based_exonic_endpoints(self):
        chain, anchors, blocks = m.structure(100, [(0, 30), (3, 70), (0, 25)])
        self.assertEqual(chain, [(130, 201)])
        self.assertEqual(anchors, [(30, 25)])
        self.assertEqual(blocks, [(101, 130), (201, 225)])

    def test_deletion_is_not_a_splice_and_deleted_bases_are_not_aligned(self):
        chain, _, blocks = m.structure(100, [(0, 30), (2, 70), (0, 25)])
        self.assertEqual(chain, [])
        self.assertFalse(m.covered(blocks, 131, 133))
        self.assertTrue(m.covered(blocks, 201, 203))

    def test_softclip_and_insertion_do_not_advance_reference(self):
        chain, anchors, _ = m.structure(100, [(4, 10), (0, 20), (1, 5), (0, 10), (3, 70), (0, 25)])
        self.assertEqual(chain, [(130, 201)])
        self.assertEqual(anchors, [(10, 25)])

    def test_contiguous_match_and_mismatch_operations_share_an_anchor(self):
        chain, anchors, _ = m.structure(100, [(7, 10), (8, 2), (0, 18), (3, 70), (7, 25)])
        self.assertEqual(chain, [(130, 201)])
        self.assertEqual(anchors, [(30, 25)])

    def test_multiple_splices_preserve_full_chain(self):
        chain, anchors, _ = m.structure(100, [(0, 30), (3, 70), (0, 25), (3, 75), (0, 40)])
        self.assertEqual(chain, [(130, 201), (225, 301)])
        self.assertEqual(anchors, [(30, 25), (25, 40)])

    def test_quality_threshold_and_flags(self):
        self.assertEqual(m.exclusions(0, 20), [])
        self.assertEqual(m.exclusions(16, 20), [])  # Alignment direction is recorded, not inferred as RNA strand.
        self.assertEqual(m.exclusions(0, 19), ["MAPQ_below_20"])
        self.assertEqual(m.exclusions(256 | 2048, 60), ["secondary", "supplementary"])
        self.assertEqual(m.exclusions(512 | 1024, 60), ["QC_failed", "duplicate_flag"])

    def test_codon_coverage_checks_all_three_bases(self):
        self.assertTrue(m.covered([(100, 101), (102, 103)], 100, 102))
        self.assertFalse(m.covered([(100, 101), (103, 104)], 100, 102))

    def test_unsupported_backward_cigar_fails(self):
        with self.assertRaises(ValueError):
            m.structure(100, [(9, 10)])


if __name__ == "__main__":
    unittest.main()
