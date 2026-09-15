import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from summarize_psc_blood_replication import adjust_pvalues, conjunction, reciprocal_crosswalk, build_conjunctions


class TestBloodReplication(unittest.TestCase):
    def test_adjustment_matches_independent_library(self):
        p = np.array([0.001, 0.7, 0.04, 0.001, 1.0, 0.03])
        for method, external in [("BH", "fdr_bh"), ("BY", "fdr_by")]:
            np.testing.assert_allclose(adjust_pvalues(p, method), multipletests(p, method=external)[1], atol=1e-15)

    def test_conjunction_retains_missing_and_opposed(self):
        p = [[0.001, 0.02], [0.001, 0.001], [0.01, np.nan], [0.03, 0.01], [0.0, 0.01]]
        b = [[1, 2], [1, -1], [1, 1], [-1, -2], [0, 1]]
        result, status = conjunction(p, b)
        np.testing.assert_allclose(result, [0.02, 1, 1, 0.03, 1])
        self.assertEqual(status[2], "unavailable_constituent")
        self.assertEqual(status[1], "discordant_or_zero_direction")

    def test_more_constituents_cannot_strengthen_uncorrected_evidence(self):
        p = np.array([[0.001, 0.002, 0.3], [0.1, 0.2, 0.05]])
        b = np.ones_like(p)
        self.assertTrue(np.all(conjunction(p, b)[0] >= conjunction(p[:, :2], b[:, :2])[0]))

    def test_family_keeps_unavailable_tests(self):
        np.testing.assert_allclose(adjust_pvalues([0.01, np.nan]), [0.02, 1])

    def test_bad_probability_rejected(self):
        for p in [[-0.1], [1.1], [np.inf]]:
            with self.assertRaises(ValueError):
                adjust_pvalues(p)
        with self.assertRaises(ValueError):
            conjunction([[0.1, 1.5]], [[1, 1]])

    def test_global_identity_not_expression_selected(self):
        f = pd.DataFrame([("ENSG1", "1", "A"), ("ENSG2", "2", "B"), ("ENSG3", "2", "B"), ("ENSG4", "3", "C"), ("ENSG4", "4", "C"), ("ENSG5", "", "D")], columns=["ensembl_gene_id", "entrez_gene_id", "hgnc_symbol"])
        mapped = reciprocal_crosswalk(f).set_index("ensembl_gene_id")
        self.assertEqual(mapped.loc["ENSG1", "entrez_gene_id"], "1")
        self.assertEqual(mapped.loc["ENSG2", "identity_status"], "ambiguous_entrez_to_ensembl")
        self.assertEqual(mapped.loc["ENSG4", "identity_status"], "ambiguous_ensembl_to_entrez")
        self.assertEqual(mapped.loc["ENSG5", "identity_status"], "no_entrez_relation")

    def test_repeated_source_relation_does_not_change_mapping(self):
        f = pd.DataFrame([("ENSG1", "1", "A")]*2, columns=["ensembl_gene_id", "entrez_gene_id", "hgnc_symbol"])
        self.assertEqual(len(reciprocal_crosswalk(f)), 1)

    def test_constituent_p_not_overwritten(self):
        f = pd.DataFrame({"p1": [0.01], "p2": [np.nan], "b1": [1.0], "b2": [1.0]})
        out = build_conjunctions(f, {"rep": [("p1", "b1"), ("p2", "b2")]})
        self.assertTrue(np.isnan(out.loc[0, "p2"]))
        self.assertEqual(out.loc[0, "rep_pvalue"], 1)

    def test_assembly_does_not_drop_missing_pvalues(self):
        from summarize_psc_blood_replication import assemble_results
        rna = pd.DataFrame({"gene_id": ["ENSG1", "ENSG2", "ENSG3"], "eligible": [True, True, False], "pvalue": [0.001, np.nan, 0.9], "effect": [1.0, 1.0, 0.1]})
        cross = pd.DataFrame([(f"ENSG{i}", str(i), f"G{i}") for i in range(1, 4)], columns=["ensembl_gene_id", "entrez_gene_id", "hgnc_symbol"])
        contrast_names = {"control": "PSC_vs_Control", "uc": "PSC_vs_UC", "pbc": "PSC_vs_PBC", "cd": "PSC_vs_CD"}
        array = pd.DataFrame([{"entrez_id": str(i), "contrast": c, "pvalue": 0.01, "estimate": 0.2} for i in range(1, 4) for c in contrast_names.values()])
        schema = {"rna_gene_id": "gene_id", "rna_filter": "eligible", "rna_pvalue": "pvalue", "rna_effect": "effect", "array_gene_id": "entrez_id", "array_contrast": "contrast", "array_pvalue": "pvalue", "array_effect": "estimate", "array_contrasts": contrast_names}
        result, coverage = assemble_results(rna, array, cross, schema)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(coverage), 3)
        self.assertEqual(result.loc[1, "replication_status"], "unavailable_constituent")
        self.assertEqual(result.loc[0, "replication_qvalue_bh"], 0.02)
        self.assertTrue(np.isnan(result.loc[1, "rna_pvalue"]))
        with self.assertRaises(ValueError):
            assemble_results(rna, array.iloc[:-1], cross, schema)
