"""Meaningful identity and missing-data checks for the exploratory follow-up."""

import json
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from analyze_coverage_repair_alternative import eligible_qtl
from coloc_common import ChainMap
from coverage_repair_common import (Reference, component_weights, coverage_parts,
                                    exact_qtl_status, normalize, reciprocal_edit,
                                    unordered_edit_candidates)


class DeletionIdentityTests(unittest.TestCase):
    def test_shifted_repeat_deletions_have_one_exact_identity(self):
        reference = Reference('chr1', 1, 'GGTCACACACAGG')
        left = normalize(reference, 4, 'CAC', 'C')
        right = normalize(reference, 8, 'CAC', 'C')
        self.assertEqual(left, (3, 'TCA', 'T'))
        self.assertEqual(left, right)
        self.assertEqual(reference.haplotype(4, 'CAC', 'C'), reference.haplotype(*left))

    def test_different_repeat_lengths_are_not_the_same_variant(self):
        ref = Reference('chr1', 1, 'GGTCACACACAGG')
        self.assertNotEqual(normalize(ref, 4, 'CAC', 'C'), normalize(ref, 4, 'CACAC', 'C'))

    def test_full_padded_representation_matches_minimal_alleles(self):
        ref = Reference('chr1', 1, 'GGTCACACACAGG')
        self.assertEqual(normalize(ref, 3, 'TCACACACAG', 'TCACACAG'), normalize(ref, 4, 'CAC', 'C'))

    def test_ref_mismatch_and_unbounded_repeat_are_rejected(self):
        ref = Reference('chr1', 100, 'AAAAC')
        with self.assertRaisesRegex(ValueError, 'REF'):
            normalize(ref, 101, 'AG', 'A')
        with self.assertRaisesRegex(ValueError, 'boundary'):
            normalize(ref, 102, 'AA', 'A')

    def test_unordered_gwas_alleles_do_not_silently_define_ref(self):
        ref = Reference('chr1', 1, 'GGTCACACACAGG')
        candidates = unordered_edit_candidates(ref, 3, 'T', 'TCA')
        self.assertEqual(len(candidates), 2)
        self.assertIn((3, 'TCA', 'T'), candidates)
        self.assertIn((3, 'T', 'TCA'), candidates)

    def test_full_allele_mapping_rejects_a_gap_even_with_matching_anchor(self):
        forward = ChainMap('chain 1 chr1 20 + 0 20 chr1 21 + 0 21 1\n5 0 1\n15\n')
        reverse = ChainMap('chain 1 chr1 21 + 0 21 chr1 20 + 0 20 1\n5 1 0\n15\n')
        a = Reference('chr1', 1, 'ACGTACGTACGTACGTACGT')
        b = Reference('chr1', 1, 'ACGTATCGTACGTACGTACGT')
        with self.assertRaisesRegex(ValueError, 'Noncontiguous'):
            reciprocal_edit(a, b, forward, reverse, 4, 'TAC', 'T', flank=1)

    def test_negative_strand_uses_reverse_complement_of_full_allele(self):
        chain = ChainMap('chain 1 chr1 12 + 0 12 chr1 12 - 0 12 1\n12\n')
        a = Reference('chr1', 1, 'AAACGTCCGTTA')
        b = Reference('chr1', 1, 'TAACGGACGTTT')
        mapped, bases, strand = reciprocal_edit(a, b, chain, chain, 4, 'CGT', 'C', flank=2)
        self.assertEqual(mapped, (7, 'ACG', 'G'))
        self.assertEqual((bases, strand), (7, '-'))


class MissingnessAndCoverageTests(unittest.TestCase):
    def test_wrong_alt_at_multiallelic_position_does_not_fill_missing_effect(self):
        row = dict(position='100', ref='A', alt='G', gene_id='GENE1')
        status, _, exact, measured = exact_qtl_status([row], 100, 'A', 'C', 'GENE1')
        self.assertEqual(status, 'position_measured_only_with_different_alleles')
        self.assertEqual((exact, measured), ([], []))

    def test_other_gene_measurement_is_distinct_from_absent_position(self):
        row = dict(position='100', ref='A', alt='C', gene_id='GENE2.1')
        self.assertEqual(exact_qtl_status([row], 100, 'A', 'C', 'GENE1')[0], 'exact_allele_measured_only_for_other_genes')
        self.assertEqual(exact_qtl_status([], 100, 'A', 'C', 'GENE1')[0], 'position_absent_across_all_genes')

    def test_coverage_uses_all_source_evidence_and_rejects_double_counting(self):
        weights = component_weights(np.log([.50, .34, .16]))
        self.assertAlmostEqual(coverage_parts(weights, [1, 0, 0], [0, 1, 0])[2], .84)
        with self.assertRaisesRegex(ValueError, 'already counted'):
            coverage_parts(weights, [1, 0, 0], [1, 0, 0])
        with self.assertRaisesRegex(ValueError, 'full normalized'):
            coverage_parts(weights[:2], [1, 0], [0, 1])
        with self.assertRaises(ValueError): component_weights([1, np.nan])

    def test_alternative_snp_filter_keeps_missing_r2_but_rejects_bad_reported_r2(self):
        plan = dict(gene_id='GENE', start_grch38_0based=0, end_grch38_exclusive=200)
        row = dict(gene_id='GENE.1', molecular_trait_id='GENE', chromosome='2', position='100', ref='A', alt='C',
                   variant='chr2_100_A_C', beta='.1', se='.05', pvalue='.0455', maf='.1', r2='NA', an='1000', ma_samples='90')
        self.assertAlmostEqual(eligible_qtl(row, set(), plan)['beta_canonical'], .1)
        with self.assertRaisesRegex(ValueError, 'R2'):
            eligible_qtl(dict(row, r2='.2'), set(), plan)
        with self.assertRaisesRegex(ValueError, 'MAF'):
            eligible_qtl(dict(row, maf='.005'), set(), plan)
        with self.assertRaisesRegex(ValueError, 'palindromic'):
            eligible_qtl(dict(row, alt='T', variant='chr2_100_A_T'), set(), plan)

    def test_real_recovered_records_do_not_turn_incomplete_components_into_passes(self):
        root = Path(__file__).resolve().parents[1]
        import pandas as pd
        coverage = pd.read_csv(root / 'data/derived/coverage-repair-pfkfb3-signal-coverage.csv')
        strongest = coverage[coverage.component.eq(1)]
        self.assertEqual(len(strongest), 2)
        self.assertFalse(strongest.qtl_90pct_gate.any())
        np.testing.assert_allclose(strongest.recoverable_signal_coverage, [.52233561719003, .8855374318799414], atol=1e-12)
        identities = json.loads((root / 'data/derived/coverage-repair-exact-identities.json').read_text())
        self.assertEqual([r['recoverable'] for r in identities['pfkfb3']], [False, True, True])
        self.assertTrue(all(r['bcftools_matches'] for r in identities['independent_normalization']))


if __name__ == '__main__': unittest.main()
