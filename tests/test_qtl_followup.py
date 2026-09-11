import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from fetch_qtl_followup import canonical_subset
from summarize_qtl_followup import beta_for_allele, junction_coordinates, select_gene_row, strand_audit


class QtlAlignmentTests(unittest.TestCase):
    def test_palindromic_pair_needs_explicit_effect_and_risk_alleles(self):
        self.assertAlmostEqual(beta_for_allele(.91, 'G', 'C', 'G', 'C'), -.91)
        self.assertAlmostEqual(beta_for_allele(.91, 'G', 'C', 'G', 'G'), .91)
        with self.assertRaises(ValueError):
            beta_for_allele(.91, 'A', 'C', 'G', 'C')
        with self.assertRaises(ValueError):
            beta_for_allele(float('nan'), 'G', 'C', 'G', 'C')

    def test_rsid_or_position_alone_cannot_select_other_alternate_allele(self):
        row = {'variant': 'chr21_42434957_A_C', 'gene_id': 'GENE', 'molecular_trait_id': 'GENE',
               'chromosome': '21', 'position': '42434957', 'ref': 'A', 'alt': 'C'}
        other = {**row, 'variant': 'chr21_42434957_A_G', 'alt': 'G'}
        self.assertEqual(select_gene_row([other, row], row['variant'], 'GENE'), row)
        with self.assertRaises(ValueError):
            select_gene_row([{**row, 'position': '42434956'}], row['variant'], 'GENE')
        with self.assertRaises(ValueError):
            select_gene_row([row, row], row['variant'], 'GENE')

    def test_splice_conversion_preserves_conflicting_strand(self):
        canonical = junction_coordinates('21:42434954:42437488:clu_31404_-')
        alternate = junction_coordinates('21:42434983:42437488:clu_31404_-')
        self.assertEqual(alternate, ('21', 42434983, 42437487, '-'))
        self.assertEqual(alternate[1]-canonical[1], 29)
        self.assertEqual(alternate[2], canonical[2])

    def test_strand_audit_does_not_count_join_duplicates_or_multigene_rows(self):
        r = {'phenotype_id': '21:1:4:clu_1_-', 'strand': '1', 'gene_count': '1'}
        other = {'phenotype_id': '21:5:8:clu_2_+', 'strand': '1', 'gene_count': '1'}
        result = strand_audit([r, r, other, {**other, 'phenotype_id': '21:9:12:clu_3_-', 'gene_count': '2'}])
        self.assertEqual(result['unique_single_gene_phenotypes'], 2)
        self.assertEqual(result['discordant_fraction'], .5)

    def test_empty_filtered_result_remains_header_only(self):
        self.assertEqual(canonical_subset('variant\tbeta', []), b'variant\tbeta\n')
        self.assertNotIn(b'0.0', canonical_subset('variant\tbeta', []))


if __name__ == '__main__':
    unittest.main()
