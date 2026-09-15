"""Adversarial full-allele and reference-orientation checks for PRKD2."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))

from audit_prkd2_missing_alleles import candidate_edits
from coverage_repair_common import Reference
from prkd2_expanded_common import key_and_sign, normalize_edit, qtl_identity


class Prkd2AlleleTests(unittest.TestCase):
    def setUp(self):
        self.reference = Reference('chr19', 1, 'GCAAACGTACGT')

    def row(self, a, b, position=2, source_row=1):
        return {'source_row': source_row, 'SNP': f'row{source_row}', 'pos': position,
                'allele_0': a, 'allele_1': b}

    def test_source_audit_includes_snvs(self):
        candidates, counts = candidate_edits([self.row('C', 'T')], self.reference)
        self.assertEqual(candidates[(2, 'C', 'T')], {(1, 'row1')})
        self.assertEqual(counts['SNV_rows'], 1)
        self.assertEqual(counts['orientations_considered'], 2)

    def test_source_risk_allele_is_not_assumed_alt(self):
        candidates, _ = candidate_edits([self.row('T', 'C')], self.reference)
        self.assertIn((2, 'C', 'T'), candidates)

    def test_source_audit_keeps_both_possible_indel_orientations(self):
        candidates, counts = candidate_edits([self.row('C', 'CA')], self.reference)
        self.assertEqual(set(candidates), {(2, 'C', 'CA'), (2, 'CA', 'C')})
        self.assertEqual(counts['reference_valid_orientations'], 2)

    def test_source_audit_rejects_both_reference_mismatches(self):
        candidates, counts = candidate_edits([self.row('A', 'T')], self.reference)
        self.assertFalse(candidates)
        self.assertEqual(counts['reference_invalid_orientations'], 2)

    def test_source_audit_retains_multiple_source_candidates(self):
        candidates, _ = candidate_edits([self.row('C', 'T'), self.row('T', 'C', source_row=2)], self.reference)
        self.assertEqual(len(candidates[(2, 'C', 'T')]), 2)

    def test_equivalent_repeat_insertions_left_align(self):
        self.assertEqual(normalize_edit(self.reference, 4, 'A', 'AA'), (2, 'C', 'CA'))

    def test_equivalent_repeat_deletions_left_align(self):
        self.assertEqual(normalize_edit(self.reference, 3, 'AA', 'A'), (2, 'CA', 'C'))

    def test_qtl_beta_coding_tracks_lexical_effect_allele(self):
        positive = qtl_identity('chr19_2_C_T', self.reference)
        negative = qtl_identity('chr19_2_C_A', self.reference)
        self.assertEqual(positive['canonical_beta_sign'], 1)
        self.assertEqual(negative['canonical_beta_sign'], -1)
        self.assertEqual(negative['snp'], 'chr19:2:A:C')

    def test_non_snv_key_preserves_ref_alt_orientation(self):
        self.assertEqual(key_and_sign(2, 'C', 'CA'), ('chr19:2:C:CA', 1))
        self.assertEqual(key_and_sign(2, 'CA', 'C'), ('chr19:2:CA:C', 1))

    def test_qtl_reference_mismatch_is_rejected(self):
        with self.assertRaises(ValueError):
            qtl_identity('chr19_2_A_T', self.reference)

    def test_qtl_wrong_chromosome_is_rejected(self):
        with self.assertRaises(ValueError):
            qtl_identity('chr18_2_C_T', self.reference)

    def test_qtl_ambiguous_snv_is_rejected(self):
        with self.assertRaises(ValueError):
            qtl_identity('chr19_2_C_N', self.reference)


if __name__ == '__main__':
    unittest.main()
