"""Checks for scientific data corruption risks in source normalization."""

import hashlib
import importlib.util
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


catalog = module('build_catalog')
fetcher = module('fetch_sources')


class DataIntegrityTests(unittest.TestCase):
    def test_secondary_signal_keeps_locus_columns(self):
        body = ET.fromstring('<tbody><tr><td rowspan="2">chr3</td><td rowspan="2">GENE</td><td>1</td><td>rs10</td></tr><tr><td>2</td><td>rs20</td></tr><tr><td>chr4</td><td>NEXT</td><td>1</td><td>rs30</td></tr></tbody>')
        rows = list(catalog.expanded_rows(body))
        self.assertEqual(rows[1], ['chr3', 'GENE', '2', 'rs20'])
        self.assertEqual(rows[2], ['chr4', 'NEXT', '1', 'rs30'])

    def test_gene_context_does_not_leak_between_loci(self):
        import openpyxl
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'source.xlsx'
            book = openpyxl.Workbook()
            sheet = book.active
            sheet.title = 'Supplementary_Data_2'
            sheet.append(['PSC GWAS Locus'])
            sheet.append(['Chr', 'Lead GWAS SNP', 'eGene', 'Molecular QTL Type', 'Tissue', 'Tissue state', 'PP4', 'Risk allele', 'Beta', 'p-value', 'Risk allele', 'Beta', 'p-value'])
            sheet.append([1, 'rs10', 'GENE1', 'eQTL', 'T cell', 'resting', 0.9, 'A', 0.2, 1e-9, 'A', -0.5, 1e-7])
            sheet.append([None, None, 'GENE1', 'eQTL', 'B cell', 'resting', 0.8, None, None, None, 'A', -0.4, 1e-6])
            sheet.append([2, 'rs20', 'GENE2', 'eQTL', 'Monocyte', 'resting', 0.95, 'G', 0.4, 1e-8, 'G', 0.5, 1e-9])
            book.save(path)
            book.close()
            records = catalog.parse_colocalisation(path)
        self.assertEqual(records[1]['lead_gwas_snp'], 'rs10')
        self.assertEqual(records[1]['gwas_risk_allele'], 'A')
        self.assertEqual(records[2]['lead_gwas_snp'], 'rs20')
        self.assertEqual(records[2]['gwas_risk_allele'], 'G')

    def test_repeated_libraries_do_not_become_distinct_donor_labels(self):
        text = '^SAMPLE = GSM1\n!Sample_title = PSC001_SC_5pr\n!Sample_description = 3pr scRNA-seq\n^SAMPLE = GSM2\n!Sample_title = PSC001_SN_3pr\n!Sample_description = 3pr scRNA-seq\n'
        records = catalog.parse_geo(text, 'GSEtest')
        self.assertEqual(len(records), 2)
        self.assertEqual(len({row['donor_label_from_title'] for row in records}), 1)
        self.assertIn('5pr', records[0]['flags'])
        self.assertIn('nuclei', records[1]['flags'])

    def test_modified_input_is_rejected(self):
        source = {'id': 'fixture', 'max_bytes': 100, 'sha256': hashlib.sha256(b'original').hexdigest()}
        with self.assertRaisesRegex(ValueError, 'Checksum changed'):
            fetcher.validate(b'changed', source)


if __name__ == '__main__':
    unittest.main()
