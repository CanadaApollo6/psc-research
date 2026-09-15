"""Failure cases that could otherwise turn missing evidence into a positive result."""

import json
import sys
import unittest
import zlib
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from analyze_prkd2_recovery import overlap, parse_blueprint, parse_dice, target_status
from fetch_prkd2_recovery_catalogue import exact_target
from fetch_prkd2_recovery_ibdverse import scan_member, variant_position


class OriginalSourceParsing(unittest.TestCase):
    def test_dice_does_not_invent_se_from_rounded_statistic(self):
        row = parse_dice(['chr19\t47203356\trs112445263\tC\tA\t.\tPASS\tGene=ENSG00000105287.8;GeneSymbol=PRKD2;Pvalue=0.224;Beta=0.22;Statistic=1.225;FDR=1\r\n'])[0]
        self.assertEqual(row['source_se'], '')
        self.assertEqual(row['source_pvalue'], '0.224')
        self.assertEqual(row['source_statistic'], '1.225')

    def test_original_dice_rejects_another_gene(self):
        with self.assertRaises(ValueError):
            parse_dice(['chr19\t1\trs1\tC\tA\t.\tPASS\tGene=ENSG000001052870.8;Pvalue=1;Beta=0;Statistic=0;FDR=1'])

    def test_blueprint_keeps_complete_indel_and_nonsignificant_row(self):
        row = parse_blueprint(['19:47185096_C_CT rs1 ENSG00000105287.8 8.1e-1 -0.004 1.0 0.9 0.14 0.12'])[0]
        self.assertEqual((row['ref'], row['alt']), ('C', 'CT'))
        self.assertEqual(row['source_pvalue'], '8.1e-1')
        self.assertEqual(row['source_se'], '0.12')

    def test_missing_blueprint_standard_error_is_not_silently_dropped(self):
        with self.assertRaises(ValueError):
            parse_blueprint(['19:47185096_C_CT rs1 ENSG00000105287.8 8.1e-1 -0.004 1.0 0.9 0.14'])


class ExactIdentityAndMissingness(unittest.TestCase):
    def setUp(self):
        self.target = {'chromosome': '19', 'position': 46700099, 'ref': 'C', 'alt': 'A'}
        self.row = {'gene_id': 'ENSG00000105287', **self.target, 'beta': '.2', 'se': '.1', 'pvalue': '.05'}

    def test_access_failure_and_absence_are_distinct(self):
        self.assertIsNone(target_status(None, self.target)['matching_rows'])
        self.assertEqual(target_status(pd.DataFrame(), self.target)['matching_rows'], 0)

    def test_same_coordinate_other_chromosome_is_not_match(self):
        row = {**self.row, 'chromosome': '18'}
        self.assertEqual(target_status(pd.DataFrame([row]), self.target)['matching_rows'], 0)
        self.assertFalse(exact_target(row, 'ENSG00000105287'))

    def test_same_coordinate_other_allele_or_gene_is_not_match(self):
        for change in [{'alt': 'G'}, {'gene_id': 'ENSG00000105288'}]:
            self.assertEqual(target_status(pd.DataFrame([{**self.row, **change}]), self.target)['matching_rows'], 0)

    def test_distinct_probe_effects_are_preserved(self):
        result = target_status(pd.DataFrame([self.row, {**self.row, 'beta': '-.1'}]), self.target)
        self.assertEqual(result['status'], 'multiple-distinct-measurements')
        self.assertEqual(result['matching_rows'], 2)

    def test_control_coordinate_is_verified_against_prior_raw_identity(self):
        root = Path(__file__).resolve().parents[1]
        plan = json.loads((root / 'config/prkd2-recovery-followup-plan.json').read_text())
        control = next(r for r in plan['targets'] if r['rsid'] == 'rs313839')
        self.assertEqual((control['position'], control['ref'], control['alt']), (46718300, 'C', 'G'))
        self.assertNotIn(46702450, [r['position'] for r in plan['targets']])


class FullEvidenceCoverage(unittest.TestCase):
    def setUp(self):
        self.gwas = pd.DataFrame({'snp': ['a'], 'source_row': [0]})
        self.full = pd.DataFrame({'source_row': [0, 1], 'N': [1000, 1000], 'MAF': [.2, .2],
                                  'case_fraction': [.5, .5], 'pvalue': [.05, .05]})
        self.rna = pd.DataFrame({'snp': ['a', 'unmatched'], 'in_fixed_region': [True, True],
                                'source_beta': ['.1', '2'], 'source_se': ['.1', '.1']})

    def test_dominant_omitted_rna_weight_is_kept_in_denominator(self):
        result = overlap(self.rna, self.gwas, self.full)
        self.assertLess(result['RNA_marginal_weight_available'], 1e-30)
        self.assertAlmostEqual(result['full_PSC_marginal_weight_available'], .5)
        self.assertFalse(result['multi_signal_model_qualified'])

    def test_missing_se_blocks_rna_coverage_instead_of_dropping_row(self):
        self.rna.loc[1, 'source_se'] = ''
        result = overlap(self.rna, self.gwas, self.full)
        self.assertIsNone(result['RNA_marginal_weight_available'])
        self.assertIsNone(result['coverage_gate'])

    def test_outside_window_rna_cannot_inflate_shared_evidence(self):
        self.rna.loc[1, 'snp'] = 'a'
        self.rna.loc[1, 'in_fixed_region'] = False
        result = overlap(self.rna, self.gwas, self.full)
        self.assertLess(result['RNA_marginal_weight_available'], 1e-30)

    def test_exact_duplicate_source_row_is_not_additional_evidence(self):
        self.rna['source_row'] = [0, 1]
        once = overlap(self.rna, self.gwas, self.full)
        duplicate = self.rna.iloc[[1]].copy()
        duplicate['source_row'] = 2
        twice = overlap(pd.concat([self.rna, duplicate]), self.gwas, self.full)
        self.assertEqual(twice['exact_duplicate_source_rows'], 1)
        self.assertEqual(twice['distinct_source_association_rows'], 2)
        self.assertAlmostEqual(once['RNA_marginal_weight_available'], twice['RNA_marginal_weight_available'])


class CompleteAuthorArchive(unittest.TestCase):
    def setUp(self):
        self.plain = (b'phenotype_id\tvariant_id\tpval_nominal\tslope\tslope_se\n'
                      b'ENSG00000105287.1\tchr19:46700099:C:A\t0.8\t0.1\t0.2\n'
                      b'ENSG000001052870\tchr19:46700099:C:A\t0.01\t0.1\t0.2\n'
                      b'ENSG00000105287.1\tchr19:46718300:C:G\t0.4\t0.1\t0.2\n')
        enc = zlib.compressobj(wbits=-15)
        self.compressed = enc.compress(self.plain) + enc.flush()

    def scan(self, data, crc=None):
        return scan_member([data[:10], data[10:]], 'ENSG00000105287', {46700099},
                           zlib.crc32(self.plain) if crc is None else crc, len(self.plain))

    def test_complete_scan_preserves_weak_associations_and_distinguishes_gene(self):
        result, _, genes, points = self.scan(self.compressed)
        self.assertEqual(result['source_rows'], 3)
        self.assertEqual(len(genes), 2)
        self.assertEqual(len(points), 2)
        self.assertIn(b'\t0.8\t', genes[0])

    def test_incomplete_or_bad_crc_member_cannot_prove_absence(self):
        for data, crc in [(self.compressed[:-4], None), (self.compressed, 0)]:
            with self.assertRaises(ValueError):
                self.scan(data, crc)

    def test_positions_preserve_chromosome_and_require_explicit_coordinate(self):
        self.assertEqual(variant_position('chr19_46700099_C_A'), 46700099)
        self.assertIsNone(variant_position('chr18_46700099_C_A'))
        self.assertIsNone(variant_position('rs112445263'))


if __name__ == '__main__':
    unittest.main()
