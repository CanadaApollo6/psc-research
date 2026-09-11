"""Check ranking versus direction and the common-gene baseline."""

import importlib.util
import unittest
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    pd = None

ROOT = Path(__file__).resolve().parents[1]
if pd is not None:
    spec = importlib.util.spec_from_file_location('summarize_pilot', ROOT / 'scripts/summarize_pilot.py')
    summarizer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(summarizer)


@unittest.skipIf(pd is None, 'Install requirements-alphagenome.lock for model-summary checks')
class PilotSummaryTests(unittest.TestCase):
    def test_absolute_ranking_does_not_change_signed_effect(self):
        scores = pd.DataFrame([
            {'gene_id': gene, 'raw_score': score, 'ontology_curie': 'cell', 'variant_scorer': summarizer.LFC}
            for gene, score in [('A', -0.4), ('A', -0.2), ('B', 0.8), ('B', -0.1), ('outside', 10.0)]
        ])
        candidates = pd.DataFrame([
            {'gene_id': gene, 'gene_name': gene, 'biotype': 'protein_coding', 'distance_to_variant': distance}
            for gene, distance in [('A', 10), ('B', 20), ('unscored', 1)]
        ])
        ranking, baseline = summarizer.rank_genes(scores, candidates, 'cell')
        self.assertEqual(ranking.gene_id.tolist(), ['B', 'A'])
        self.assertAlmostEqual(ranking.iloc[0].median_absolute_lfc, 0.45)
        self.assertAlmostEqual(ranking.iloc[1].median_signed_lfc, -0.3)
        self.assertEqual(baseline.gene_id, 'A')
        self.assertNotIn('unscored', ranking.gene_id.tolist())


if __name__ == '__main__':
    unittest.main()
