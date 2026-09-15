"""Guard against incomplete sources, false donor counts and inferred model facts."""
import json
import sys
import tempfile
import unittest
import zlib
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from analyze_ibdverse_model_audit import sample_counts
from fetch_ibdverse_model_members import decode_member, extract_gene
from read_ibdverse_model_metadata import MetadataReader, read_column

ROOT = Path(__file__).resolve().parents[1]


class CompleteSourceRequirements(unittest.TestCase):
    def setUp(self):
        self.plain = (b'phenotype_id\tvariant_id\tpval_nominal\tqval\tqvals\n'
                      b'ENSG00000105287.8\tchr19:1:C:A\t0.8\t0.900000\t\n'
                      b'ENSG000001052870\tchr19:2:C:A\t0.0001\t0.02\t\n')
        encoder = zlib.compressobj(wbits=-15)
        self.compressed = encoder.compress(self.plain) + encoder.flush()
        self.member = {'uncompressed_bytes': len(self.plain), 'crc32': zlib.crc32(self.plain)}

    def test_exact_gene_matching_preserves_weak_rows_numeric_strings_and_blank_fields(self):
        selected, total, kept = extract_gene(decode_member(self.compressed, self.member), 'ENSG00000105287')
        self.assertEqual((total, kept), (2, 1))
        self.assertIn(b'0.8\t0.900000\t\n', selected)
        self.assertNotIn(b'ENSG000001052870', selected)

    def test_truncated_or_appended_compression_cannot_establish_gene_absence(self):
        for corrupt in [self.compressed[:-3], self.compressed + b'additional unscanned member']:
            with self.subTest(corrupt=corrupt), self.assertRaises(ValueError):
                decode_member(corrupt, self.member)

    def test_bad_crc_and_wrong_size_reject_complete_looking_member(self):
        for change in [{'crc32': 0}, {'uncompressed_bytes': len(self.plain) + 1}]:
            with self.subTest(change=change), self.assertRaises(ValueError):
                decode_member(self.compressed, self.member | change)

    def test_duplicate_header_or_truncated_row_is_not_silently_accepted(self):
        for data in [b'phenotype_id\tvariant_id\tphenotype_id\n', b'phenotype_id\tvariant_id\nENSG00000105287\n']:
            with self.subTest(data=data), self.assertRaises(ValueError):
                extract_gene(data, 'ENSG00000105287')


class MetadataMissingness(unittest.TestCase):
    def category(self, codes):
        source = h5py.File('synthetic-metadata', 'w', driver='core', backing_store=False)
        group = source.create_group('category')
        group.create_dataset('categories', data=np.array([b'donor_A', b'donor_B']))
        group.create_dataset('codes', data=np.array(codes))
        return source, group

    def test_missing_category_is_not_last_donor(self):
        source, group = self.category([0, -1, 1])
        with source:
            self.assertEqual(list(read_column(group)), ['donor_A', None, 'donor_B'])

    def test_invalid_category_cannot_become_a_valid_individual(self):
        for invalid in [-2, 2]:
            source, group = self.category([invalid])
            with source, self.assertRaises(ValueError):
                read_column(group)

    def test_metadata_reader_forbids_unbounded_reads_and_network_fallback_offline(self):
        plan = json.loads((ROOT / 'config/ibdverse-model-audit-h5ad-plan.json').read_text())
        reader = MetadataReader(plan, offline=True)
        with self.assertRaises(ValueError):
            reader.read()
        with self.assertRaises(ValueError):
            reader.seek(-1)
        # A deterministic synthetic gap, not a scientific HDF5 data request.
        reader.record = {'blocks': []}
        with patch('read_ibdverse_model_metadata.fetch', side_effect=AssertionError('Network forbidden')):
            with self.assertRaises(ValueError):
                reader.block(0)

    def run_sample_case(self, mapping):
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory)
            prefix = 'E-MTAB-16986-restored-'
            (raw / (prefix + 'sample_metadata.tsv.data')).write_text('sanger_sample_id\tbiopsy_type\nsample_A\tblood\n')
            (raw / (prefix + 'individual_metadata.tsv.data')).write_text('individual_id\tdisease_status\ndonor_A\tCD\n')
            (raw / (prefix + 'sample_individual_mapping.tsv.data')).write_text(mapping)
            with patch('analyze_ibdverse_model_audit.RAW', raw):
                return sample_counts()

    def test_duplicate_mapping_does_not_inflate_donor_count(self):
        with self.assertRaises(ValueError):
            self.run_sample_case('sanger_sample_id\tindividual_id\nsample_A\tdonor_A\nsample_A\tdonor_A\n')

    def test_missing_sample_mapping_does_not_disappear_from_denominator(self):
        with self.assertRaises(ValueError):
            self.run_sample_case('sanger_sample_id\tindividual_id\nsample_B\tdonor_A\n')


class ScientificQualification(unittest.TestCase):
    def test_eligibility_and_conditional_pcs_do_not_qualify_a_shared_cause_posterior(self):
        audit = json.loads((ROOT / 'reports/ibdverse-model-audit.json').read_text())
        self.assertFalse(audit['ready_for_new_shared_cause_posterior'])
        self.assertEqual(audit['new_H4_posteriors'], 0)
        for row in audit['selected_source_summary']:
            self.assertGreater(row['RNA_eligible_public_donor_labels'], 0)
            self.assertGreater(row['expression_PCs_if_RNA_count_equals_model_N_and_five_genotype_PCs_only'], 0)
            self.assertFalse(row['actual_executed_model_N_observed'])
            self.assertFalse(row['actual_selected_PCs_observed'])


if __name__ == '__main__':
    unittest.main()
