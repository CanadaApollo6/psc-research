import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from coloc_common import ChainMap, reciprocal_base
from fetch_coloc_regions import read_gwas_range
from prepare_coloc_inputs import canonical_key, eligible_gwas
from prepare_coloc_platform_inputs import platform_gwas


class ColocInputTests(unittest.TestCase):
    def setUp(self):
        self.forward = ChainMap('chain 100 chr1 1000 + 100 130 chr1 1000 + 200 231 1\n10 5 6\n15\n')
        self.reverse = ChainMap('chain 100 chr1 1000 + 200 231 chr1 1000 + 100 130 1\n10 6 5\n15\n')

    def test_chain_gaps_boundaries_and_reciprocity_are_explicit(self):
        self.assertEqual(reciprocal_base(self.forward,self.reverse,'chr1',109),('chr1',209,'+'))
        self.assertIsNone(reciprocal_base(self.forward,self.reverse,'chr1',110))
        self.assertEqual(reciprocal_base(self.forward,self.reverse,'chr1',115),('chr1',216,'+'))
        self.assertIsNone(reciprocal_base(self.forward,self.reverse,'chr1',130))

    def test_reverse_chain_and_ambiguous_mappings_are_not_guessed(self):
        forward=ChainMap('chain 100 chr1 1000 + 100 110 chr1 1000 - 200 210 1\n10\n')
        reverse=ChainMap('chain 100 chr1 1000 + 790 800 chr1 1000 - 890 900 1\n10\n')
        self.assertEqual(reciprocal_base(forward,reverse,'chr1',100),('chr1',799,'-'))
        ambiguous=ChainMap('chain 100 chr1 1000 + 100 110 chr1 1000 + 200 210 1\n10\n\nchain 90 chr1 1000 + 100 110 chr1 1000 + 300 310 2\n10\n')
        self.assertIsNone(reciprocal_base(ambiguous,self.reverse,'chr1',100))

    def test_canonical_identity_preserves_distinct_alleles(self):
        self.assertEqual(canonical_key('chr1',201,'T','A'),canonical_key('chr1',201,'A','T'))
        self.assertNotEqual(canonical_key('chr1',201,'A','T'),canonical_key('chr1',201,'A','C'))
        with self.assertRaises(ValueError):canonical_key('chr1',201,'A','AT')

    def test_gwas_risk_orientation_and_palindromic_frequency_checks(self):
        row={'chr':'1','pos':'101','SNP':'rsTest','allele_0':'A','allele_1':'T','freq_1_controls':'0.2','p':'0.01','or':'1.2','platform':'both'}
        region={'gene_name':'test','start_grch38_0based':190,'end_grch38_exclusive':250}
        ref={'chr1:201:A:T':{'effect_frequency':0.2}}
        result=eligible_gwas(row,region,self.forward,self.reverse,ref)
        self.assertEqual(result['status'],'eligible'); self.assertEqual(result['canonical_sign'],1)
        flipped={**row,'allele_0':'T','allele_1':'A','freq_1_controls':'0.8'}
        result=eligible_gwas(flipped,region,self.forward,self.reverse,ref)
        self.assertEqual(result['status'],'eligible'); self.assertEqual(result['canonical_sign'],-1)
        self.assertEqual(eligible_gwas(row,region,self.forward,self.reverse,{})['status'],'palindromic_reference_frequency_unavailable')
        self.assertEqual(eligible_gwas({**row,'freq_1_controls':'0.5'},region,self.forward,self.reverse,ref)['status'],'palindromic_frequency_disagreement')
        self.assertEqual(eligible_gwas({**row,'platform':'omni'},region,self.forward,self.reverse,ref)['status'],'platform_sample_size_unavailable')

    def test_only_complete_gwas_lines_can_enter_a_regional_subset(self):
        line=b'1 rsTest 101 A T 0.2 0.2 0.2 0.99 1.2 0.1 0.01 both\n'
        data=b'partial row\n'+line+b'1 next partial'
        result=read_gwas_range(data,1000)
        self.assertEqual(len(result),1)
        self.assertEqual(result[0][0],(1,101))
        self.assertEqual(result[0][1],1012)

    def test_platform_amendment_uses_documented_counts_and_preserves_source(self):
        raw={'chr':'1','pos':'101','SNP':'rsTest','allele_0':'A','allele_1':'T','freq_1_controls':'0.2','p':'0.01','or':'1.2','platform':'omni'}
        region={'gene_name':'test','start_grch38_0based':190,'end_grch38_exclusive':250}
        reference={'chr1:201:A:T':{'effect_frequency':0.2}}
        counts={'omni':{'N':11386,'cases':2170,'controls':9216}}
        result=platform_gwas(raw,region,self.forward,self.reverse,reference,counts)
        self.assertEqual(result['status'],'eligible')
        self.assertEqual((result['N'],result['cases'],result['controls']),(11386,2170,9216))
        self.assertEqual((result['platform'],raw['platform']),('omni','omni'))
        result=platform_gwas({**raw,'platform':'unknown'},region,self.forward,self.reverse,reference,counts)
        self.assertEqual(result['status'],'platform_sample_size_unavailable')


if __name__=='__main__':
    unittest.main()
