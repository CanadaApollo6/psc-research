import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_ubash3a_strand import motifs, regtools_strand
from summarize_ubash3a_junctions import (UnestimableError, aggregate_donors, hc3_linear,
                                       parse_genotypes, parse_junction)


class PublicJunctionTests(unittest.TestCase):
    def test_sparse_response_must_be_complete_before_missing_sample_means_zero(self):
        target = {'chromosome': 'chr21', 'start_0based': 100, 'end_0based_exclusive': 200,
                  'biological_strand': '+'}
        header = 'DataSource:Type\tsnaptron_id\tchromosome\tstart\tend\tlength\tstrand\tleft_motif\tright_motif\tsamples\tsamples_count\tcoverage_sum\n'
        row = 'srav3h:I\t1\tchr21\t101\t200\t100\t+\tGT\tAG\t,7:3,8:2\t2\t5\n'
        counts, _ = parse_junction(header + row, target)
        self.assertEqual(counts, {'7': 3, '8': 2})
        for invalid in [header, 'Traceback: server error', header + row.replace('7:3,8:2', '7:3,7:2'),
                        header + row.replace('\t2\t5\n', '\t2\t6\n'),
                        header + row.replace('\t101\t', '\t100\t')]:
            with self.assertRaises((ValueError, KeyError)):
                parse_junction(invalid, target)

    def test_vcf_missing_rsid_does_not_replace_coordinate_and_allele_validation(self):
        head = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n'
        row = '21\t43855067\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1/1\t./.\n'
        self.assertEqual(parse_genotypes(head + row), {'A': 0, 'B': 1, 'C': 2, 'D': None})
        for invalid in [row.replace('43855067', '42434957'), row.replace('\tA\tC\t', '\tA\tG\t'),
                        row.replace('1/1', '1/2')]:
            with self.assertRaises(ValueError):
                parse_genotypes(head + invalid)

    def test_repeated_runs_are_combined_before_donor_ratio(self):
        a = {'run': 'run1', 'center': 'lab', 'canonical': 7, 'plus29': 1}
        b = {'run': 'run2', 'center': 'lab', 'canonical': 18, 'plus29': 2}
        donors = aggregate_donors({'donor1': [a, b], 'donor2': [{**a, 'run': 'run3'}]},
                                 {'donor1': {'pop': 'CEU'}}, {'donor1': 0})
        self.assertEqual(len(donors), 2)
        self.assertEqual((donors[0]['rna_runs'], donors[0]['canonical'], donors[0]['plus29']), (2, 25, 3))
        self.assertEqual(donors[0]['pair_total'], 28)
        self.assertEqual(donors[0]['dosage_c'], 0)
        self.assertIsNone(donors[1]['dosage_c'])

    def test_hc3_matches_simple_regression_and_rejects_singleton_population(self):
        x = np.column_stack([np.ones(5), np.arange(5)])
        beta, cov = hc3_linear(x, [1, 2, 1, 4, 5])
        np.testing.assert_allclose(beta, [.6, 1])
        self.assertAlmostEqual(cov[1, 1], .08 + .02 * (.4 / .7) ** 2)
        x = [[1, 0, 0], [1, 1, 0], [1, 2, 0], [1, 0, 0], [1, 1, 1]]
        with self.assertRaises(UnestimableError):
            hc3_linear(x, [0, .1, .2, .1, .5])

    def test_sequence_orientation_and_paired_end_library_flags(self):
        self.assertEqual(motifs('CCCCGTCCAGTT', 0, 4, 10), ('GT', 'AG', 'CT', 'AC'))
        with self.assertRaises(ValueError):
            motifs('CCCCGTCCAGTT', 0, 4, 100)
        for flag in (83, 163):
            self.assertEqual(regtools_strand(flag, 1), '+')
            self.assertEqual(regtools_strand(flag, 2), '-')
        for flag in (99, 147):
            self.assertEqual(regtools_strand(flag, 1), '-')
            self.assertEqual(regtools_strand(flag, 2), '+')


if __name__ == '__main__':
    unittest.main()
