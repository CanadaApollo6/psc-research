import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from match_comparison_variants import choose_pair, dosage_r2, maf, parse_vcf
from prepare_matching_covariates import mutation_class, static_pair
from summarize_matched_comparison import LFC, matched_contrast, score_variant_table


class MatchedComparisonTests(unittest.TestCase):
    def test_missing_and_constant_genotypes_do_not_become_zero_ld(self):
        self.assertEqual(dosage_r2([0, 0, 0], [0, 1, 2], 2), (None, 3))
        self.assertEqual(dosage_r2([0, 1, np.nan], [0, 1, 2], 3), (None, 2))
        self.assertAlmostEqual(dosage_r2([0, 1, 2], [2, 1, 0], 3)[0], 1)
        with self.assertRaises(ValueError):
            dosage_r2([0, 1], [1, 2, 0])
        self.assertEqual(maf([0, 1, np.nan, 2]), (.5, 3))
        self.assertEqual(maf([np.nan]), (None, 0))

    def test_vcf_is_aligned_by_coordinates_alleles_and_sample_identity(self):
        head = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\tA\tC\n'
        body = '2\t101\t.\tG\tT\t.\tPASS\t.\tGT:DP\t1|1:8\t0|1:9\t./.:0\n'
        targets = [{'rsid': 'rs1', 'chromosome': '2', 'position_grch37_1based': 101, 'reference_bases': 'G', 'alternate_bases': 'T'}]
        result, ids = parse_vcf(head + body, ['A', 'B', 'C'], targets)
        np.testing.assert_equal(result['rs1'], [1, 2, np.nan])
        self.assertEqual(ids['rs1'], '.')
        for bad in [body.replace('\t101\t', '\t100\t'), body.replace('\tG\tT\t', '\tG\tA\t'), body.replace('1|1:8', '1|2:8')]:
            with self.assertRaises(ValueError):
                parse_vcf(head + bad, ['A', 'B', 'C'], targets)
        with self.assertRaises(ValueError):
            parse_vcf(head + body, ['A', 'missing'], targets)

    def test_pair_selection_skips_mutually_correlated_top_comparators(self):
        ordered = [{'rsid': 'rs1'}, {'rsid': 'rs2'}, {'rsid': 'rs3'}]
        def r2(a, b):
            return (.9 if {a, b} == {'rs1', 'rs2'} else .01), 503
        selected = choose_pair(ordered, r2, .1)
        self.assertEqual([r['rsid'] for r in selected[:2]], ['rs1', 'rs3'])
        self.assertIsNone(choose_pair(ordered, lambda a, b: (None, 300), .1))

    def test_static_matching_keeps_boundary_and_mutation_requirements(self):
        anchor = {'mutation_class': 'transition', 'noncoding_context': 'intron', 'position_grch38_1based': 300000, 'nearest_tss_distance_bp': 100}
        limits = {'maximum_distance_from_anchor_bp': 250000}
        self.assertTrue(static_pair(anchor, {**anchor, 'position_grch38_1based': 550000}, limits))
        self.assertFalse(static_pair(anchor, {**anchor, 'position_grch38_1based': 550001}, limits))
        self.assertFalse(static_pair(anchor, {**anchor, 'mutation_class': 'transversion'}, limits))
        self.assertEqual(mutation_class('T', 'C'), 'transition')
        self.assertEqual(mutation_class('C', 'A'), 'transversion')

    def test_gene_score_uses_absolute_before_median_and_keeps_same_universe(self):
        genes = [{'gene_id': 'ENSG1', 'gene_name': 'gene1'}, {'gene_id': 'ENSG2', 'gene_name': 'gene2'}]
        rows = [{'gene_id': g+'.1', 'track_name': t, 'raw_score': value, 'variant_scorer': LFC, 'ontology_curie': 'CL:test'}
                for g, t, value in [('ENSG1','track1',-.4),('ENSG1','track2',.2),('ENSG2','track1',.1),('ENSG2','track2',.1)]]
        frame = pd.DataFrame(rows)
        ranked, _ = score_variant_table(frame, genes, ['track1','track2'], 'CL:test')
        self.assertEqual(ranked.iloc[0].gene_id, 'ENSG1')
        self.assertAlmostEqual(ranked.iloc[0].gene_effect, .3)
        self.assertAlmostEqual(ranked.iloc[0].signed_median, -.1)
        for malformed in [frame.iloc[:-1], pd.concat([frame, frame.iloc[[0]]]), frame.assign(raw_score=[0,0,0,np.nan])]:
            with self.assertRaises(ValueError):
                score_variant_table(malformed, genes, ['track1','track2'], 'CL:test')

    def test_contrast_keeps_missing_data_unavailable(self):
        self.assertAlmostEqual(matched_contrast(.4, .1, .3), .2)
        self.assertIsNone(matched_contrast(.4, None, .3))
        self.assertIsNone(matched_contrast(.4, np.nan, .3))


if __name__ == '__main__':
    unittest.main()
