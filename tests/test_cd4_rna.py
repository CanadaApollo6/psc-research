import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from summarize_cd4_rna import direction_agreement, exact_allele_rows, gene_id, measured_numbers, missing_pair_reason, unique_gene_measurement


class Cd4RnaTests(unittest.TestCase):
    def setUp(self):
        self.gene = 'ENSG00000000001'
        self.row = {'variant':'chr2_100_G_C', 'chromosome':'2', 'position':'100', 'ref':'G', 'alt':'C',
                    'rsid':'rs1', 'gene_id':self.gene, 'molecular_trait_id':self.gene,
                    'beta':'.5', 'se':'.1', 'pvalue':'.00001', 'an':'176', 'ac':'30', 'r2':'NA', 'median_tpm':'NA'}

    def test_exact_allele_join_does_not_guess_palindromic_orientation_or_rsid(self):
        swapped = {**self.row, 'variant':'chr2_100_C_G', 'ref':'C', 'alt':'G'}
        self.assertEqual(exact_allele_rows([swapped], 'chr2', 100, 'G', 'C'), [])
        self.assertEqual(exact_allele_rows([{**self.row,'rsid':'rs99'}], '2', 100, 'G', 'C')[0]['rsid'], 'rs99')
        with self.assertRaises(ValueError): exact_allele_rows([{**self.row,'position':'101'}], '2',100,'G','C')

    def test_alias_duplicates_collapse_but_conflicting_measurements_do_not(self):
        row, count, aliases = unique_gene_measurement([self.row,{**self.row,'rsid':'rs2'}],self.gene)
        self.assertEqual(count,2); self.assertEqual(aliases,'rs1;rs2'); self.assertEqual(row['beta'],'.5')
        with self.assertRaises(ValueError): unique_gene_measurement([self.row,{**self.row,'beta':'-.5'}],self.gene)

    def test_gene_version_is_removed_without_gene_symbol_or_transcript_substitution(self):
        self.assertEqual(gene_id(self.gene+'.12'), self.gene)
        for wrong in ['IL2RA','ENST00000000001',self.gene+'-PAR-Y']:
            with self.assertRaises(ValueError): gene_id(wrong)

    def test_missing_measurements_keep_distinct_unavailable_states(self):
        self.assertEqual(missing_pair_reason([], None, None),'exact_variant_absent')
        self.assertEqual(missing_pair_reason([self.row],None,None),'gene_annotation_absent')
        self.assertEqual(missing_pair_reason([self.row],{},False),'outside_source_cis_window')
        self.assertEqual(missing_pair_reason([self.row],{},True),'gene_variant_pair_absent')

    def test_uncertainty_and_allele_counts_cannot_be_invalid_or_fabricated(self):
        result = measured_numbers(self.row,88)
        self.assertEqual(result['sample_size_from_an'],88)
        self.assertIsNone(result['imputation_r2']); self.assertIsNone(result['median_tpm'])
        self.assertAlmostEqual(result['wald_95_lower'],.304)
        for bad in [{'beta':'nan'},{'se':'0'},{'pvalue':'-1'},{'an':'175'},{'an':'178'},{'ac':'177'}]:
            with self.assertRaises(ValueError): measured_numbers({**self.row,**bad},88)
        self.assertTrue(measured_numbers({**self.row,'pvalue':'0'},88)['pvalue_reported_as_zero'])

    def test_mixed_model_tracks_and_zero_observed_beta_do_not_create_direction_agreement(self):
        self.assertIsNone(direction_agreement('mixed_or_zero',.4))
        self.assertIsNone(direction_agreement('increase',0))
        self.assertTrue(direction_agreement('increase',.4))
        self.assertFalse(direction_agreement('increase',-.4))


if __name__ == '__main__':
    unittest.main()
