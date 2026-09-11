"""Guard allele orientation, ambiguous mappings, and interval boundaries."""

import importlib.util
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('normalize_benchmark', ROOT / 'scripts/normalize_benchmark.py')
normalizer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(normalizer)


class BenchmarkIntegrityTests(unittest.TestCase):
    def test_multiallelic_variant_uses_explicit_published_alternate(self):
        reference, alleles = normalizer.validate_alleles({'allele_string': 'A/C/G'}, 'C', 'A', 'A')
        self.assertEqual(reference, 'A')
        self.assertEqual(alleles, ['A', 'C', 'G'])
        with self.assertRaisesRegex(ValueError, 'Selected alternate'):
            normalizer.validate_alleles({'allele_string': 'A/C/G'}, 'T', 'A', 'A')

    def test_model_reference_mismatch_fails(self):
        with self.assertRaisesRegex(ValueError, 'Reference-base mismatch'):
            normalizer.validate_alleles({'allele_string': 'A/C'}, 'C', 'A', 'G')

    def test_reference_risk_allele_reverses_expression_sign(self):
        evidence = [{'gwas_risk_allele': 'A', 'functional_risk_allele_as_published': 'A', 'functional_beta': '-0.5'}]
        risk, expected, _ = normalizer.interpret_direction('A', 'C', evidence)
        self.assertEqual((risk, expected), ('A', 1))

    def test_incompatible_risk_label_is_not_guessed(self):
        evidence = [{'gwas_risk_allele': 'A', 'functional_risk_allele_as_published': 'A', 'functional_beta': '-0.5'}]
        risk, expected, note = normalizer.interpret_direction('C', 'G', evidence)
        self.assertEqual(risk, 'A')
        self.assertIsNone(expected)
        self.assertIn('excluded', note)

    def test_missing_risk_direction_stays_missing(self):
        risk, expected, note = normalizer.interpret_direction('G', 'A', [])
        self.assertIsNone(risk)
        self.assertIsNone(expected)

    def test_multiple_primary_mappings_are_rejected(self):
        mapping = {'assembly_name':'GRCh38','seq_region_name':'21','coord_system':'chromosome','strand':1,'start':100,'end':100}
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            normalizer.primary_mapping({'mappings':[mapping, dict(mapping)]}, 'GRCh38', '21')

    def test_reverse_strand_tss_and_half_open_interval(self):
        def gene(identifier, start, end, strand):
            return {'gene_id':identifier, 'seq_region_name':'21','start':start,'end':end,'strand':strand,'biotype':'protein_coding'}
        records = [gene('left', 96, 96, 1), gene('right', 105, 105, 1), gene('outside', 106, 106, 1), gene('reverse', 95, 100, -1)]
        start0, end0, genes = normalizer.candidate_genes(records, '21', 101, 10)
        self.assertEqual((start0, end0), (95, 105))
        self.assertEqual({g['gene_id'] for g in genes}, {'left','right','reverse'})
        self.assertEqual(genes[0]['gene_id'], 'reverse')
        self.assertEqual(genes[0]['tss_1based'], 100)


if __name__ == '__main__':
    unittest.main()
