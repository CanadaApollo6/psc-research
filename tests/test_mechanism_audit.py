import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from fetch_gwas_audit_rows import complete_rows

try:
    import pandas as pd
    from summarize_splice_followup import aligned_delta
except ImportError:
    pd = None


class OriginalGwasRowsTests(unittest.TestCase):
    def test_range_fragments_cannot_become_association_rows(self):
        row = b'19 rs313839 47221557 G C 0.84 0.87 0.84 0.999 1.322 0.051 2.11865164971472e-08 omni'
        rows = complete_rows(b'partial leading fragment\n' + row + b'\ntruncated trailing row')
        self.assertEqual(rows, [((19, 47221557), 'rs313839', row.decode())])
        with self.assertRaises(ValueError):
            complete_rows(b'fragment\n19 wrong_schema\nfragment')


@unittest.skipIf(pd is None, 'Install requirements-alphagenome.lock for model-summary checks')
class FeatureAlignmentTests(unittest.TestCase):
    def test_reordered_and_absent_junctions_do_not_create_false_effects(self):
        reference = pd.DataFrame({'junction': ['a', 'b', 'd'], 'reference': [10., 20., 40.]})
        alternate = pd.DataFrame({'junction': ['b', 'a', 'c'], 'alternate': [21., 8., 30.]})
        rows = aligned_delta(reference, alternate, ['junction']).set_index('junction')
        self.assertEqual(rows.loc['a', 'delta_alt_minus_ref'], -2.)
        self.assertEqual(rows.loc['b', 'delta_alt_minus_ref'], 1.)
        self.assertTrue(pd.isna(rows.loc['c', 'delta_alt_minus_ref']))
        self.assertTrue(pd.isna(rows.loc['d', 'delta_alt_minus_ref']))
        with self.assertRaises(ValueError):
            aligned_delta(pd.concat([reference, reference]), alternate, ['junction'])


if __name__ == '__main__':
    unittest.main()
