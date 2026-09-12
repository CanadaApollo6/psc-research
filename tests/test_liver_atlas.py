"""Synthetic checks for the donor-aware liver-atlas expression core."""

from __future__ import annotations

import sys
import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_liver_atlas import (  # noqa: E402
    GENES,
    aggregate_adata,
    build_contrasts,
    combined_exclusion_counts,
    prepare_metadata,
    resolve_target_features,
    summarize_expression,
    validate_frozen_thresholds,
    validate_counts_matrix,
    verify_frozen_source,
)


def make_adata(rows, counts, feature_names=("ETS2", "PRKD2", "PFKFB3", "BCL2L11")):
    obs = pd.DataFrame(rows, index=[f"cell-{i}" for i in range(len(rows))])
    counts = np.asarray(counts, dtype=np.int64)
    raw_var = pd.DataFrame(
        {"feature_name": list(feature_names)},
        index=[f"ENSG{i:011d}" for i in range(len(feature_names))],
    )
    # X intentionally has normalized, non-integer values.  The analysis must
    # use the separate integer raw object for both target UMIs and denominators.
    normalized = counts.astype(float) / np.maximum(counts.sum(axis=1, keepdims=True), 1.0)
    var = pd.DataFrame(index=[f"v{i}" for i in range(len(feature_names))])
    result = ad.AnnData(X=normalized, obs=obs, var=var)
    result.raw = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs.copy(), var=raw_var)
    return result


def row(donor, disease, assay="10x 3' v3", library="L1", author="hepatocyte", cell_type="hepatocyte", suspension="nucleus"):
    return {
        "donor_id": donor,
        "author_cell_type": author,
        "cell_type": cell_type,
        "disease": disease,
        "assay": assay,
        "suspension_type": suspension,
        "library_uuid": library,
    }


class LiverAtlasAggregationTests(unittest.TestCase):
    def test_frozen_thresholds_cannot_drift(self):
        plan = {"primary_min_cells_per_donor_cell_type": 20,
                "sensitivity_min_cells_per_donor_cell_type": 50,
                "minimum_donors_per_cohort_for_contrast": 3}
        self.assertEqual(validate_frozen_thresholds(plan), {"primary_min_cells": 20, "sensitivity_min_cells": 50,
                                                             "minimum_donors_for_contrast": 3})
        for field, value in [("primary_min_cells_per_donor_cell_type", 19),
                             ("sensitivity_min_cells_per_donor_cell_type", 49),
                             ("minimum_donors_per_cohort_for_contrast", 2)]:
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, "thresholds changed"):
                    validate_frozen_thresholds({**plan, field: value})

    def test_frozen_source_hash_catches_same_size_override_with_different_basename(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "renamed-input.h5ad"
            path.write_bytes(b"abd")
            source = {"bytes": 3, "sha256": "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"}
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                verify_frozen_source(path, source)
            path.write_bytes(b"abc")
            self.assertEqual(verify_frozen_source(path, source)["bytes"], 3)

    def test_exclusion_counts_add_across_datasets(self):
        results = [SimpleNamespace(qc={"excluded_by_reason": {"cell_type_unknown": 4, "author_cell_type_doublet": 2}}),
                   SimpleNamespace(qc={"excluded_by_reason": {"cell_type_unknown": 7, "author_cell_type_doublet": 1}})]
        self.assertEqual(combined_exclusion_counts(results), {"author_cell_type_doublet": 3, "cell_type_unknown": 11})

    def test_repeated_library_and_donor_rows_are_combined_with_all_gene_denominator(self):
        rows = [
            row("D1", "primary sclerosing cholangitis", library="L1"),
            row("D1", "primary sclerosing cholangitis", library="L1"),
            row("D1", "primary sclerosing cholangitis", library="L2"),
            row("D1", "primary sclerosing cholangitis", assay="10x 5' v2", library="L3"),
        ]
        data = make_adata(rows, [[2, 0, 0, 3], [0, 3, 0, 0], [1, 0, 0, 0], [4, 0, 0, 0]])
        result = aggregate_adata(data, "sn", chunk_size=2)
        self.assertEqual(len(result.pseudobulk), 2)
        same_assay = result.pseudobulk[result.pseudobulk.assay.eq("10x 3' v3")].iloc[0]
        self.assertEqual(same_assay.cell_count, 3)
        self.assertEqual(same_assay.library_count, 2)
        self.assertEqual(same_assay.all_feature_umis, 9)
        self.assertEqual(same_assay.ETS2_umis, 3)
        self.assertEqual(same_assay.ETS2_detected_cells, 2)
        self.assertAlmostEqual(same_assay.ETS2_cpm, 333333.3333333333)
        self.assertEqual(result.qc["count_source"], "raw.X")
        self.assertEqual(result.qc["retained_cells"], 4)

    def test_assay_strata_do_not_merge(self):
        rows = [
            row("D1", "normal", assay="10x 3' v3"),
            row("D1", "normal", assay="10x 5' v2"),
        ]
        result = aggregate_adata(make_adata(rows, [[1, 0, 0, 0], [5, 0, 0, 0]]), "sn")
        self.assertEqual(set(result.pseudobulk.assay), {"10x 3' v3", "10x 5' v2"})
        self.assertEqual(set(result.pseudobulk.cell_count), {1})

    def test_missing_and_zero_genes_are_distinguished(self):
        rows = [row("D1", "normal"), row("D1", "normal")]
        data = make_adata(rows, [[0, 0, 4], [0, 0, 2]], feature_names=("ETS2", "PRKD2", "PFKFB3"))
        result = aggregate_adata(data, "sn")
        group = result.pseudobulk.iloc[0]
        self.assertEqual(group["PRKD2_status"], "present")
        self.assertEqual(group["PRKD2_umis"], 0)
        self.assertEqual(group["PRKD2_detection_fraction"], 0.0)
        self.assertEqual(group["UBASH3A_status"], "missing_feature_name")
        self.assertTrue(pd.isna(group["UBASH3A_cpm"]))
        summary = summarize_expression(result.pseudobulk, result.cell_accounting, thresholds=(1,), gene_audit=result.gene_audit)
        zero = summary[(summary.gene.eq("PRKD2")) & summary.cohort.eq("Control")].iloc[0]
        missing = summary[(summary.gene.eq("UBASH3A")) & summary.cohort.eq("Control")].iloc[0]
        self.assertEqual(zero.status, "evaluable_descriptive")
        self.assertEqual(zero.donor_weighted_median_log2_cpm1, 0.0)
        self.assertEqual(missing.status, "missing_feature_name")

    def test_metadata_exclusions_remain_in_accounting_and_leave_count_denominator(self):
        rows = [
            row("D1", "normal"),
            row("D1", "normal", author="Doublet hepatocyte"),
            row("D1", "normal", cell_type="unknown"),
            row("D1", "normal", author="hybrid--label"),
        ]
        result = aggregate_adata(make_adata(rows, [[1, 0, 0, 0]] * 4), "sn")
        self.assertEqual(result.qc["retained_cells"], 1)
        self.assertEqual(result.qc["excluded_cells"], 3)
        self.assertEqual(int(result.cell_accounting.cell_count.sum()), 4)
        self.assertEqual(int(result.cell_accounting.excluded_cell_count.sum()), 3)

    def test_raw_counts_are_required_and_normalized_x_is_not_used(self):
        rows = [row("D1", "normal"), row("D1", "normal")]
        data = make_adata(rows, [[3, 0, 0, 0], [0, 0, 0, 0]])
        result = aggregate_adata(data, "sn")
        self.assertEqual(result.pseudobulk.iloc[0].all_feature_umis, 3)
        data.raw = None
        with self.assertRaisesRegex(ValueError, "raw.*missing"):
            aggregate_adata(data, "sn")

    def test_minimum_donor_gate_and_sn_processing_flag(self):
        rows = []
        counts = []
        for donor in ("P1", "P2", "P3"):
            rows.extend([row(donor, "primary sclerosing cholangitis")] * 2)
            counts.extend([[2, 0, 0, 0], [0, 1, 0, 0]])
        for donor in ("C1", "C2"):
            rows.extend([row(donor, "normal")] * 2)
            counts.extend([[1, 0, 0, 0], [0, 1, 0, 0]])
        result = aggregate_adata(make_adata(rows, counts), "sn")
        contrasts = build_contrasts(result.pseudobulk, result.cell_accounting, genes=("ETS2",), thresholds=(2,))
        contrast = contrasts.iloc[0]
        self.assertEqual(contrast.status, "insufficient_donors")
        self.assertEqual(contrast.psc_eligible_donors, 3)
        self.assertEqual(contrast.control_eligible_donors, 2)

    def test_cross_assay_contrast_is_blocked_instead_of_fallback(self):
        rows = []
        counts = []
        for donor in ("P1", "P2", "P3"):
            rows.append(row(donor, "primary sclerosing cholangitis", assay="10x 5' v2", suspension="cell"))
            counts.append([2, 0, 0, 0])
        for donor in ("C1", "C2", "C3"):
            rows.append(row(donor, "normal", assay="10x 5' v1", suspension="cell"))
            counts.append([1, 0, 0, 0])
        result = aggregate_adata(make_adata(rows, counts), "sc")
        contrasts = build_contrasts(result.pseudobulk, result.cell_accounting, genes=("ETS2",), thresholds=(1,))
        self.assertEqual(len(contrasts), 1)
        self.assertEqual(contrasts.iloc[0].status, "blocked_cross_chemistry")
        self.assertFalse(bool(contrasts.iloc[0].shared_assay))

    def test_count_qc_rejects_fractional_negative_nonfinite_and_bad_cell_ids(self):
        with self.assertRaisesRegex(ValueError, "non-integer"):
            validate_counts_matrix(sparse.csr_matrix([[1.5]]))
        with self.assertRaisesRegex(ValueError, "negative"):
            validate_counts_matrix(sparse.csr_matrix([[-1.0]]))
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            validate_counts_matrix(sparse.csr_matrix([[np.nan]]))
        data = make_adata([row("D1", "normal"), row("D1", "normal")], [[1, 0, 0, 0], [1, 0, 0, 0]])
        data.obs_names = ["same", "same"]
        with self.assertRaisesRegex(ValueError, "duplicate cell IDs"):
            aggregate_adata(data, "sn")

    def test_feature_map_requires_unique_exact_feature_name_matches(self):
        raw_var = pd.DataFrame({"feature_name": ["ETS2", "ETS2"]}, index=["a", "b"])
        audit = resolve_target_features(raw_var, 2, genes=("ETS2", "PRKD2"))
        self.assertEqual(audit.loc[audit.gene.eq("ETS2"), "status"].item(), "ambiguous_feature_name")
        self.assertEqual(audit.loc[audit.gene.eq("PRKD2"), "status"].item(), "missing_feature_name")


if __name__ == "__main__":
    unittest.main()
