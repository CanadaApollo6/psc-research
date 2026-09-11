import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from splice_context import donor_contexts, junction_in_scope, position_index


class SpliceContextTests(unittest.TestCase):
    def test_transcript_numbering_and_base_conversion(self):
        exon = {'id': 'shared-exon', 'start': 100, 'end': 122}
        following = {'id': 'next-exon', 'start': 160, 'end': 180}
        gene = {'strand': 1, 'assembly_name': 'GRCh38', 'Transcript': [
            {'id': 'short', 'version': 1, 'biotype': 'protein_coding', 'Exon': [exon, following]},
            {'id': 'long', 'version': 1, 'biotype': 'protein_coding', 'is_canonical': 1,
             'Exon': [{'id': 'extra', 'start': 10, 'end': 20}, exon, following]},
        ]}
        rows = donor_contexts(gene, 125)
        self.assertEqual([r['upstream_exon_number'] for r in rows], [1, 2])
        self.assertTrue(all(r['junction_start_0based'] == 122 and r['intron_offset_1based'] == 3 for r in rows))
        self.assertTrue(all(r['model_donor_position_0based'] == 121 for r in rows))
        self.assertEqual(position_index(121, 100, 200), 21)
        with self.assertRaises(ValueError):
            position_index(200, 100, 200)
        with self.assertRaises(ValueError):
            position_index(122, 100, 200, 128)

    def test_junction_scope_preserves_direct_and_skipping_events(self):
        self.assertTrue(junction_in_scope(122, 159, '+', 122, 159))
        self.assertTrue(junction_in_scope(122, 200, '+', 122, 159))
        self.assertTrue(junction_in_scope(100, 159, '+', 122, 159))
        self.assertTrue(junction_in_scope(100, 200, '+', 122, 159))
        self.assertFalse(junction_in_scope(125, 150, '+', 122, 159))
        self.assertFalse(junction_in_scope(122, 159, '-', 122, 159))


if __name__ == '__main__':
    unittest.main()
