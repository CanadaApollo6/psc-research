"""Failure guards for missingness, strand identity and BED coordinate conversion."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from analyze_prkd2_mechanism_peaks import overlaps
from audit_prkd2_mechanism_training import accessions, specimen_ids
from summarize_prkd2_mechanism import one_value, complete_rank


class MechanismGuards(unittest.TestCase):
    def test_missing_gene_does_not_shrink_rank_universe(self):
        self.assertEqual(complete_rank([{'gene_id': 'a', 'v': 1.0}, {'gene_id': 'b', 'v': None}], 'v'), {})

    def test_gene_rank_ties_use_id(self):
        result = complete_rank([{'gene_id': 'b', 'v': .1}, {'gene_id': 'a', 'v': .1}], 'v')
        self.assertEqual(result, {'a': 1, 'b': 2})

    def test_duplicate_score_identity_is_rejected(self):
        frame = pd.DataFrame([{'track': 'RNA', 'value': 1}, {'track': 'RNA', 'value': 2}])
        with self.assertRaises(ValueError): one_value(frame, track='RNA')

    def test_missing_score_is_not_zero(self):
        self.assertIsNone(one_value(pd.DataFrame([{'track': 'RNA'}]), track='ATAC'))

    def test_strand_specific_score_identity(self):
        frame = pd.DataFrame([{'name': 'CAGE', 'strand': '+', 'value': 1}, {'name': 'CAGE', 'strand': '-', 'value': -1}])
        self.assertEqual(one_value(frame, name='CAGE', strand='-')['value'], -1)

    def test_bed_exclusive_end_and_inclusive_flank(self):
        # Base 101 is [100,101); BED [90,100) does not include that base.
        self.assertFalse(overlaps(90, 100, 100, 101))
        self.assertTrue(overlaps(100, 101, 100, 101))
        # Base 1001 ±250 maps to [750,1251), including the flank endpoints.
        self.assertTrue(overlaps(750, 751, 750, 1251))
        self.assertTrue(overlaps(1250, 1251, 750, 1251))
        self.assertFalse(overlaps(1251, 1252, 750, 1251))

    def test_replicates_share_donor_and_missing_donor_is_unresolved(self):
        exp = {'replicates': [
            {'library': {'biosample': {'accession': 'ENCBS188BKX', 'donor': '/human-donors/ENCDO265AAA/'}}},
            {'library': {'biosample': {'accession': 'ENCBS865RXK', 'donor': {'accession': 'ENCDO265AAA'}}}},
            {'library': {'biosample': {'accession': 'ENCBS609ENC'}}}]}
        biosamples, donors, unresolved = specimen_ids(exp)
        self.assertEqual(len(biosamples), 3)
        self.assertEqual(donors, ['ENCDO265AAA'])
        self.assertEqual(unresolved, 1)

    def test_accession_lists_preserve_membership_without_substring_hits(self):
        self.assertEqual(accessions('ENCSR012PII,ENCSR400VWA,ENCSR012PII,', 'ENCSR'), ['ENCSR012PII', 'ENCSR400VWA'])


if __name__ == '__main__': unittest.main()
