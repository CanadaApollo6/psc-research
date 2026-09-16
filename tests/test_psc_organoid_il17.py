"""Synthetic-only organoid analysis tests. Never read real counts or target effects."""
from __future__ import annotations

import copy
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, norm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts import analyze_psc_organoid_il17 as m


def metadata_fixture() -> pd.DataFrame:
    rows = []
    for donor, disease in m.DONOR_DISEASE.items():
        for treatment in ["CNT", "IL17"]:
            rows.append({"library_id": donor + "_" + treatment, "donor_id": donor,
                "disease": disease, "treatment": treatment, "n_cells": 200,
                "original_optional_column": "NA"})
    rows[-1]["n_cells"] = 46343 - 15 * 200
    return pd.DataFrame(rows).set_index("library_id", drop=False)


def synthetic_counts(n_features: int = 100) -> pd.DataFrame:
    metadata = metadata_fixture()
    rng = np.random.default_rng(7613)
    baseline = rng.uniform(500, 4000, n_features)
    donor_effects = dict(zip(sorted(m.DONOR_DISEASE), rng.normal(0, 0.3, 8)))
    means = np.asarray([baseline * np.exp(donor_effects[row.donor_id]) for row in metadata.itertuples()])
    psc_treat = (metadata.disease.eq("PSC") & metadata.treatment.eq("IL17")).to_numpy()
    treat = metadata.treatment.eq("IL17").to_numpy()
    means[psc_treat, 0] *= 8
    means[treat, 1] *= 4
    means[psc_treat, 1] /= 4
    values = rng.negative_binomial(35, 35 / (35 + means))
    return pd.DataFrame(values, index=metadata.index, columns=[f"synthetic_{i}" for i in range(n_features)])


def plan_fixture() -> dict:
    pin = {"path": "not-read", "sha256": "a" * 64}
    return {"status": "frozen", "frozen_at_utc": "2026-01-01T00:00:00Z", "organoid": {
        "series": "GSE239283", "population": "all_published_retained_organoid_cells",
        "donor_disease": dict(m.DONOR_DISEASE), "expected_library_count": 16,
        "expected_cell_count": 46343, "expected_feature_count": 100,
        "feature_namespace": "synthetic_original_IDs", "feature_universe_description": "synthetic_complete_universe",
        "filter": dict(m.FILTER_SETTINGS), "model": copy.deepcopy(m.MODEL_SETTINGS),
        "sensitivity": dict(m.SENSITIVITY_SETTINGS), "n_cpus": 1, "software": m.current_software(),
        "counts": dict(pin), "features": dict(pin), "libraries": dict(pin), "count_quantity_qualification": dict(pin),
        "source_artifacts": [{**pin, "url": "https://example.org/synthetic", "retrieved_at_utc": "2026-01-01T00:00:00Z"}],
        "code_artifacts": [dict(pin)], "ready_for_inference": True, "raw_count_gate_passed": True,
        "count_quantity": "raw_RNA_counts"}}


def disk_fixture(directory: Path) -> Path:
    meta = metadata_fixture(); counts = synthetic_counts()
    counts.iloc[:, -1] = 0
    counts.iloc[:, -2] = 1
    ids = counts.columns.tolist()
    frame = counts.T.rename_axis("source_feature_id").reset_index()
    count_path = directory / "counts.tsv.gz"
    frame.to_csv(count_path, sep="\t", index=False)
    feature_path = directory / "features.tsv"
    pd.DataFrame({"source_feature_id": ids[::-1], "literal_source_symbol": ["NA"] * len(ids),
        "unknown_original_annotation": ["unchanged"] * len(ids)}).to_csv(feature_path, sep="\t", index=False)
    meta["total_released_counts"] = counts.sum(axis=1)
    sample_path = directory / "libraries.tsv"
    meta.iloc[::-1].to_csv(sample_path, sep="\t", index=False)
    qualification_path = directory / "quantity-qualification.json"
    qualification_path.write_text('{"synthetic_only": true, "raw_count_gate_passed": true}')
    plan = plan_fixture(); arm = plan["organoid"]
    for name, path in [("counts", count_path), ("features", feature_path), ("libraries", sample_path), ("count_quantity_qualification", qualification_path)]:
        arm[name] = {"path": str(path), "sha256": m.sha256(path), "bytes": path.stat().st_size}
    arm["source_artifacts"] = [{**arm["counts"], "url": "https://example.org/synthetic_counts",
        "retrieved_at_utc": "2026-01-01T00:00:00Z"}]
    script = Path(m.__file__).resolve()
    arm["code_artifacts"] = [{"path": str(script), "sha256": m.sha256(script)}]
    plan_path = directory / "plan.json"; plan_path.write_text(json.dumps(plan))
    return plan_path


class TestSourceIdentity(unittest.TestCase):
    def test_literal_NA_annotation_full_feature_identity_and_shuffled_join(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp); plan = json.loads(disk_fixture(base).read_text())["organoid"]
            data = m.read_counts(Path(plan["counts"]["path"]), Path(plan["features"]["path"]))
            self.assertEqual(data.features.source_feature_id.tolist(), data.counts.columns.tolist())
            self.assertTrue(data.features.literal_source_symbol.eq("NA").all())
            self.assertTrue(data.features.unknown_original_annotation.eq("unchanged").all())
            metadata = m.read_table(Path(plan["libraries"]["path"]))
            joined = m.align_samples(data.counts, metadata, m.DONOR_DISEASE, 46343)
            self.assertEqual(joined.index.tolist(), data.counts.index.tolist())
            self.assertTrue(joined.original_optional_column.eq("NA").all())
            self.assertEqual(data.counts.iloc[:, -1].sum(), 0)

    def test_duplicate_header_and_wrong_row_width_fail(self):
        for contents in ["source_feature_id\tA\tA\ng\t1\t2\n", "source_feature_id\tA\ng\t1\t2\n", "source_feature_id\tA\ng\n"]:
            with self.subTest(contents=contents), tempfile.TemporaryDirectory() as temp:
                path = Path(temp) / "bad.tsv"; path.write_text(contents)
                with self.assertRaises(m.AnalysisError): m.read_table(path)

    def test_duplicate_feature_ids_and_nonexact_features_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            count = Path(temp) / "c.tsv"; feature = Path(temp) / "f.tsv"
            count.write_text("source_feature_id\tA\tB\ng1\t1\t2\ng2\t3\t4\n")
            for text in ["source_feature_id\ng1\ng1\n", "source_feature_id\ng1\ng3\n"]:
                feature.write_text(text)
                with self.assertRaises(m.AnalysisError): m.read_counts(count, feature)
            feature.write_text("source_feature_id\ng1\ng2\n")
            count.write_text("source_feature_id\tA\tB\ng1\t1\t2\ng1\t3\t4\n")
            with self.assertRaises(m.AnalysisError): m.read_counts(count, feature)

    def test_source_column_swap_and_missing_or_extra_library_fail(self):
        counts = synthetic_counts(); metadata = metadata_fixture()
        for bad in [counts.rename(index={counts.index[0]: "unrecognized_source_column"}), counts.iloc[1:],
                    pd.concat([counts, counts.iloc[[0]]])]:
            with self.subTest(), self.assertRaises(m.AnalysisError): m.align_samples(bad, metadata)
        bad = pd.concat([metadata, metadata.iloc[[0]]])
        with self.assertRaises(m.AnalysisError): m.align_samples(counts, bad)

    def test_unpaired_or_inconsistent_disease_and_missing_metadata_fail(self):
        for field, value in [("treatment", "IL17"), ("disease", "nonPSC"), ("donor_id", ""), ("disease", np.nan), ("n_cells", 0), ("n_cells", 2.5)]:
            metadata = metadata_fixture().copy()
            metadata[field] = metadata[field].astype(object)
            metadata.iloc[0, metadata.columns.get_loc(field)] = value
            with self.subTest(field=field, value=value), self.assertRaises(m.AnalysisError): m.validate_pairs(metadata)
        with self.assertRaises(m.AnalysisError): m.validate_pairs(metadata_fixture().iloc[:-1])

    def test_frozen_donor_map_and_released_cell_count_enforced(self):
        metadata = metadata_fixture(); counts = synthetic_counts()
        with self.assertRaises(m.AnalysisError): m.align_samples(counts, metadata, {"wrong": "PSC"})
        with self.assertRaises(m.AnalysisError): m.align_samples(counts, metadata, m.DONOR_DISEASE, 46342)

    def test_invalid_counts_never_rounded_imputed_or_counted_as_library(self):
        for value in [-1, 0.5, np.nan, np.inf, "NA", True, 2**54]:
            array = np.asarray([[1, 2], [3, value]], dtype=object)
            with self.subTest(value=value), self.assertRaises(m.AnalysisError): m.validate_counts(array)
        for array in [np.zeros((2, 3)), np.ones(4), np.empty((0, 2)), np.ones((2, 2), dtype=bool), [[2**53-1, 2], [1, 1]]]:
            with self.subTest(array=array), self.assertRaises(m.AnalysisError): m.validate_counts(array)
        np.testing.assert_array_equal(m.validate_counts([[0, 1], [2, 3]]), [[0, 1], [2, 3]])


    def test_fractional_or_underflow_text_counts_never_round_into_integer_data(self):
        for token in ["9.9999999999999999", "10.0000000000000001", "-1e-400", "1e-400", "1.0", "1e3", "9007199254740992", " 10", "+10"]:
            with self.subTest(token=token), self.assertRaises(m.AnalysisError):
                m.validate_counts([["1", token], ["2", "3"]])
        np.testing.assert_array_equal(m.validate_counts([["00010", "0"], ["12", "34"]]), [[10, 0], [12, 34]])


class TestDesignAndMath(unittest.TestCase):
    def test_rank_interaction_orientation_no_redundant_disease_and_leverage(self):
        metadata = metadata_fixture(); design, qc = m.build_design(metadata)
        self.assertEqual((qc["rank"], qc["n_columns"], qc["residual_df"]), (10, 10, 6))
        self.assertEqual(qc["max_identical_full_design_rows"], 1)
        self.assertEqual(qc["cooks_eligible_libraries_at_3_identical_rows"], 0)
        self.assertEqual(qc["donors_by_disease"], {"PSC": 4, "nonPSC": 4})
        self.assertAlmostEqual(qc["leverage_max"], 0.625)
        self.assertNotIn("disease_PSC", design.columns)
        beta = np.zeros(10); beta[-2:] = [2, 3]
        y = design.to_numpy() @ beta
        changes, disease = m.paired_changes(metadata, pd.DataFrame({"g": y}, index=metadata.index))
        self.assertAlmostEqual(changes.loc[disease.eq("PSC"), "g"].mean() - changes.loc[disease.eq("nonPSC"), "g"].mean(), 3)
        with self.assertRaises(m.AnalysisError): m.check_design(design.assign(disease_PSC=metadata.disease.eq("PSC").astype(float)))

    def test_cell_number_does_not_change_design_or_donor_sample_size(self):
        meta = metadata_fixture(); design, _ = m.build_design(meta)
        meta.n_cells *= 10_000
        new, qc = m.build_design(meta)
        pd.testing.assert_frame_equal(design, new)
        self.assertEqual(qc["n_donors"], 8); self.assertEqual(qc["n_libraries"], 16)

    def test_count_filter_boundary_and_fixed_universe(self):
        counts = pd.DataFrame(np.ones((16, 3), dtype=int), columns=["pass", "fail", "zero"])
        counts.loc[:3, "pass"] = 10; counts.loc[:2, "fail"] = 10; counts.zero = 0
        filtering = m.expression_filter(counts)
        self.assertEqual(filtering.eligible.tolist(), [True, False, False])
        self.assertEqual(filtering.libraries_count_ge_threshold.tolist(), [4, 3, 0])
        self.assertTrue(filtering.loc["zero", "all_zero_released"])

    def test_ratio_math_even_basis_is_median_logs_not_median_ratios(self):
        counts = pd.DataFrame([[1, 9, 0], [9, 1, 8], [4, 4, 2]])
        factors, basis = m.ratio_size_factors(counts)
        ratios = counts.iloc[:, :2].to_numpy() / np.prod(counts.iloc[:, :2], axis=0).to_numpy()**(1/3)
        np.testing.assert_allclose(factors, np.sqrt(ratios[:, 0] * ratios[:, 1]))
        self.assertEqual(basis, 2)
        self.assertNotAlmostEqual(factors[0], np.median(ratios[0]))
        with self.assertRaises(m.AnalysisError): m.ratio_size_factors(pd.DataFrame([[0, 1], [1, 0]]))

    def test_log_cpm_uses_complete_unfiltered_released_totals(self):
        counts = pd.DataFrame([[10, 90, 900], [20, 180, 1800]])
        result = m.log2_cpm(counts)
        np.testing.assert_allclose(result[:, 0], np.log2(1+10_000))
        self.assertNotAlmostEqual(result[0, 0], m.log2_cpm(counts.iloc[:, :1])[0, 0])

    def test_fixed_BH_BY_match_statsmodels_with_NA_correction_only(self):
        p = np.asarray([0.001, np.nan, 0.04, 0, 1, 0.02, np.nan, 0.04])
        for method, reference in [("bh", "fdr_bh"), ("by", "fdr_by")]:
            np.testing.assert_allclose(m.fixed_adjust(p, method), multipletests(np.where(np.isnan(p), 1, p), method=reference)[1], rtol=1e-14)
        self.assertTrue(np.isnan(p[[1, 6]]).all())
        self.assertEqual(m.fixed_adjust(p)[1], 1)
        self.assertGreater(m.fixed_adjust(p)[0], multipletests(p[np.isfinite(p)], method="fdr_bh")[1][0])
        for values in [[-0.1], [1.1], [np.inf], [[0.1]]]:
            with self.assertRaises(m.AnalysisError): m.fixed_adjust(values)

    def test_welch_matches_scipy_including_unequal_variances_and_df(self):
        index = [f"d{i}" for i in range(8)]
        disease = pd.Series(["PSC"]*4 + ["nonPSC"]*4, index=index)
        values = np.array([[1, 7], [2, 8], [5, 9], [8, 10], [-1, 1], [0, 20], [0.5, -10], [1, 40]])
        changes = pd.DataFrame(values, index=index, columns=["g1", "g2"])
        result = m.welch_changes(changes, disease)
        expected = ttest_ind(values[:4], values[4:], equal_var=False)
        np.testing.assert_allclose(result.t, expected.statistic)
        np.testing.assert_allclose(result.pvalue, expected.pvalue)
        np.testing.assert_allclose(result.df, expected.df)
        np.testing.assert_allclose(result.ci95_lower, expected.confidence_interval().low)
        self.assertTrue(result.n_PSC_donors.eq(4).all())

    def test_zero_variance_and_NA_cases_never_leak_into_nominal_tests(self):
        index = [f"d{i}" for i in range(8)]; disease = pd.Series(["PSC"]*4 + ["nonPSC"]*4, index=index)
        changes = pd.DataFrame({"g": [1]*4+[0]*4}, index=index)
        result = m.welch_changes(changes, disease)
        self.assertTrue(result.pvalue.isna().all()); self.assertTrue(result.q_bh_sensitivity.eq(1).all())
        changes.iloc[0, 0] = np.nan
        for method in [m.welch_changes, m.leave_one_donor_out, m.allocation_sensitivity]:
            with self.assertRaises(m.AnalysisError): method(changes, disease)
        with self.assertRaises(m.AnalysisError): m.welch_changes(changes.fillna(0), disease.iloc[::-1])

    def test_pairing_never_uses_row_position_or_cells(self):
        metadata = metadata_fixture().iloc[::-1]
        y = pd.DataFrame({"g": np.where(metadata.treatment.eq("IL17"), 3, 0)}, index=metadata.index)
        changes, disease = m.paired_changes(metadata, y)
        self.assertEqual(len(changes), 8); self.assertTrue(changes.g.eq(3).all())
        self.assertEqual(disease.value_counts().to_dict(), {"PSC": 4, "nonPSC": 4})
        with self.assertRaises(m.AnalysisError): m.paired_changes(metadata, y.iloc[::-1])

    def test_leaveout_removes_one_whole_donor_no_filter_and_preserves_features(self):
        index = [f"d{i}" for i in range(8)]; disease = pd.Series(["PSC"]*4 + ["nonPSC"]*4, index=index)
        changes = pd.DataFrame({"g1": [2, 3, 4, 5, 0, 1, 0, 1], "allzero": np.zeros(8)}, index=index)
        summary, all_values = m.leave_one_donor_out(changes, disease)
        self.assertEqual(all_values.shape, (8, 2)); self.assertEqual(summary.index.tolist(), changes.columns.tolist())
        for donor in index:
            kept = changes.drop(index=donor); labels = disease.drop(index=donor)
            expected = kept[labels.eq("PSC")].mean() - kept[labels.eq("nonPSC")].mean()
            np.testing.assert_allclose(all_values.loc[donor], expected)
        self.assertEqual(summary.loc["g1", "n_same_nonzero_sign"], 8)
        self.assertEqual(summary.loc["allzero", "n_same_nonzero_sign"], 0)

    def test_all70_allocations_complements_ties_resolution_and_zero(self):
        index = [f"d{i}" for i in range(8)]; disease = pd.Series(["PSC"]*4 + ["nonPSC"]*4, index=index)
        changes = pd.DataFrame({"separated": [1]*4 + [0]*4, "zero": np.zeros(8), "tiny": [1e-14]*4+[0]*4}, index=index)
        summary, effects, assignments = m.allocation_sensitivity(changes, disease)
        self.assertEqual(len(effects), 70); self.assertEqual(assignments.is_observed.sum(), 1)
        np.testing.assert_allclose(effects.to_numpy(), -effects.to_numpy()[::-1])
        self.assertEqual(summary.loc["separated", "n_as_or_more_extreme"], 2)
        self.assertEqual(summary.loc["separated", "descriptive_tail_fraction"], 2/70)
        self.assertEqual(summary.loc["zero", "descriptive_tail_fraction"], 1)
        self.assertEqual(summary.loc["tiny", "descriptive_tail_fraction"], 1)

    def test_contrast_covariance_retains_nonzero_cross_term(self):
        design, _ = m.build_design(metadata_fixture())
        x = design.to_numpy(); mu = np.arange(16) * 20 + 500
        covariance = m.nb_covariance(design, mu, 0.12)
        contrasts = m.contrast_vectors(design)
        t = contrasts["nonPSC_treatment"]; i = contrasts["interaction"]; p = contrasts["PSC_treatment"]
        actual = p @ covariance @ p
        self.assertAlmostEqual(actual, t @ covariance @ t + i @ covariance @ i + 2 * t @ covariance @ i)
        self.assertGreater(abs(2 * t @ covariance @ i), 0.01)
        weights = mu/(1+mu*0.12); information = x.T @ np.diag(weights) @ x
        inverse = np.linalg.solve(information+np.eye(10)*1e-6, np.eye(10))
        np.testing.assert_allclose(covariance, inverse @ information @ inverse)


class TestNegativeBinomial(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metadata = metadata_fixture(); cls.counts = synthetic_counts()
        cls.design, _ = m.build_design(cls.metadata)
        cls.fit = m.fit_negative_binomial(cls.counts, cls.metadata, cls.design)

    def test_synthetic_positive_PSC_interaction_and_negative_second_feature(self):
        table = self.fit.results["interaction"]
        self.assertGreater(table.log2FoldChange.iloc[0], 2)
        self.assertLess(table.log2FoldChange.iloc[1], -1)
        self.assertEqual(table.index.tolist(), self.counts.columns.tolist())
        np.testing.assert_allclose(table.q_bh, m.fixed_adjust(table.pvalue))
        np.testing.assert_allclose(table.q_by, m.fixed_adjust(table.pvalue, "by"))
        np.testing.assert_allclose(table.pvalue, 2*norm.sf(np.abs(table.stat)))

    def test_two_within_contrasts_are_joint_model_not_independent_and_correct_covariance(self):
        p = self.fit.results["PSC_treatment"]; c = self.fit.results["nonPSC_treatment"]; i = self.fit.results["interaction"]
        np.testing.assert_allclose(p.log2FoldChange, c.log2FoldChange+i.log2FoldChange)
        np.testing.assert_allclose(p.lfcSE**2, c.lfcSE**2+i.lfcSE**2 + 2*self.fit.feature_diagnostics.treatment_interaction_covariance_log2, rtol=1e-8)

    def test_no_cooks_exclusion_no_replacement_and_full_flags_preserved(self):
        fit = self.fit
        self.assertFalse(fit.library_diagnostics.cooks_filter_eligible.any())
        self.assertFalse(fit.feature_diagnostics.cooks_automatic_exclusion.any())
        self.assertEqual(fit.cooks.shape, self.counts.shape)
        self.assertEqual(fit.coefficients.shape, (100, 10))
        self.assertIn(fit.diagnostics["dispersion_fit_type"], ["parametric", "mean"])
        for flag in ["genewise_converged", "MAP_converged", "LFC_converged"]:
            self.assertIn(flag, fit.feature_diagnostics)

    def test_misalignment_or_changed_design_fails_before_fit(self):
        with self.assertRaises(m.AnalysisError): m.fit_negative_binomial(self.counts.iloc[::-1], self.metadata, self.design)
        bad = self.design.copy(); bad.iloc[0, -1] = 1
        with self.assertRaises(m.AnalysisError): m.fit_negative_binomial(self.counts, self.metadata, bad)


class TestCountBoundaryEstimability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metadata = metadata_fixture(); cls.counts = synthetic_counts()
        cls.counts.loc[cls.metadata.disease.eq("PSC"), "synthetic_0"] = 0
        cls.counts.loc[cls.metadata.disease.eq("PSC") & cls.metadata.treatment.eq("IL17"), "synthetic_1"] = 0
        cls.counts.loc[cls.metadata.donor_id.eq("PSC_4"), "synthetic_2"] = 0
        cls.counts.loc[cls.metadata.disease.eq("nonPSC") & cls.metadata.treatment.eq("IL17"), "synthetic_3"] = 0
        cls.design, _ = m.build_design(cls.metadata)
        cls.support = m.count_arm_support(cls.counts, cls.metadata)
        cls.fit = m.fit_negative_binomial(cls.counts, cls.metadata, cls.design)

    def test_all_zero_disease_group_is_unidentified_not_a_finite_interaction(self):
        self.assertTrue(m.expression_filter(self.counts).loc["synthetic_0", "eligible"])
        self.assertEqual(self.support.loc["synthetic_0", "count_sum_PSC_CNT"], 0)
        self.assertEqual(self.support.loc["synthetic_0", "count_sum_PSC_IL17"], 0)
        self.assertEqual(self.support.loc["synthetic_0", "n_zero_total_donor_pairs_PSC"], 4)
        for label in ["interaction", "PSC_treatment"]:
            row = self.fit.results[label].loc["synthetic_0"]
            self.assertFalse(row.count_supported_contrast)
            self.assertEqual(row.zero_required_count_arms, "PSC_CNT|PSC_IL17")
            for field in ["log2FoldChange", "lfcSE", "stat", "pvalue", "ci95_lower_log2fc", "ci95_upper_log2fc"]:
                self.assertTrue(pd.isna(row[field]))
            self.assertEqual(row.q_bh, 1); self.assertEqual(row.q_by, 1)
            # The native regularized package reports finite values despite the
            # nonexistent disease-group response information in the source data.
            self.assertTrue(np.isfinite(row.pydeseq2_unmasked_pvalue))
        self.assertTrue(self.fit.results["nonPSC_treatment"].loc["synthetic_0", "count_supported_contrast"])
        self.assertTrue(np.isfinite(self.fit.results["nonPSC_treatment"].loc["synthetic_0", "pvalue"]))

    def test_one_zero_treatment_arm_is_boundary_and_other_group_can_be_tested(self):
        for gene, failed, valid, arm in [("synthetic_1", "PSC_treatment", "nonPSC_treatment", "PSC_IL17"),
                                         ("synthetic_3", "nonPSC_treatment", "PSC_treatment", "nonPSC_IL17")]:
            self.assertTrue(m.expression_filter(self.counts).loc[gene, "eligible"])
            for label in ["interaction", failed]:
                row = self.fit.results[label].loc[gene]
                self.assertFalse(row.count_supported_contrast)
                self.assertEqual(row.zero_required_count_arms, arm)
                self.assertTrue(np.isnan(row.pvalue)); self.assertEqual(row.q_bh, 1)
                self.assertTrue(np.isfinite(row.pydeseq2_unmasked_pvalue))
            row = self.fit.results[valid].loc[gene]
            self.assertTrue(row.count_supported_contrast)
            self.assertEqual(row.pvalue, row.pydeseq2_unmasked_pvalue)

    def test_zero_total_individual_pair_is_flagged_not_silently_deleted(self):
        self.assertEqual(self.support.loc["synthetic_2", "n_zero_total_donor_pairs_PSC"], 1)
        for table in self.fit.results.values():
            self.assertTrue(table.loc["synthetic_2", "count_supported_contrast"])
        self.assertEqual(self.fit.library_diagnostics.shape[0], 16)
        self.assertEqual(self.fit.cooks.shape[0], 16)

    def test_family_and_native_covariance_are_preserved_after_contrast_specific_mask(self):
        for table in self.fit.results.values():
            self.assertEqual(len(table), len(self.counts.columns))
            np.testing.assert_allclose(table.q_bh, multipletests(table.pvalue.fillna(1), method="fdr_bh")[1])
            np.testing.assert_allclose(table.q_by, multipletests(table.pvalue.fillna(1), method="fdr_by")[1])
            self.assertIn("pydeseq2_unmasked_finite_family_padj", table)
        p = self.fit.results["PSC_treatment"]; c = self.fit.results["nonPSC_treatment"]; i = self.fit.results["interaction"]
        np.testing.assert_allclose(p.pydeseq2_unmasked_log2FoldChange,
            c.pydeseq2_unmasked_log2FoldChange + i.pydeseq2_unmasked_log2FoldChange)
        self.assertEqual(self.fit.diagnostics["count_boundary_unavailable"], {"interaction": 3, "nonPSC_treatment": 1, "PSC_treatment": 2})


class TestFrozenExecution(unittest.TestCase):
    def test_software_raw_quantity_freeze_and_method_gates(self):
        plan = plan_fixture(); self.assertIs(m.validate_plan(plan), plan["organoid"])
        changes = [lambda p: p.update(status="draft"), lambda p: p.update(frozen_at_utc="2999-01-01T00:00:00Z"),
            lambda p: p.update(frozen_at_utc="2026-01-01"), lambda p: p["organoid"]["software"].update(numpy="wrong"),
            lambda p: p["organoid"]["model"].update(cooks_filter=True), lambda p: p["organoid"]["filter"].update(min_libraries=3),
            lambda p: p["organoid"].update(raw_count_gate_passed=False), lambda p: p["organoid"].update(count_quantity="SCT_corrected_counts"),
            lambda p: p["organoid"].update(ready_for_inference=False), lambda p: p["organoid"].update(code_artifacts=[]),
            lambda p: p["organoid"].update(source_artifacts=[])]
        for change in changes:
            test = copy.deepcopy(plan); change(test)
            with self.subTest(change=change), self.assertRaises(m.AnalysisError): m.validate_plan(test)
        plan["organoid"].update(raw_count_gate_passed=False, count_quantity="unqualified")
        m.validate_plan(plan, for_inference=False)

    def test_plan_allows_qualified_other_raw_universe_and_preserves_unknown_fields(self):
        plan = plan_fixture(); plan["organoid"]["expected_feature_count"] = 36601
        plan["unknown_original"] = {"field": "NA"}; before = copy.deepcopy(plan)
        m.validate_plan(plan)
        self.assertEqual(plan, before)

    def test_pin_changes_lengths_missing_and_partial_sources_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "file"; path.write_text("source")
            pin = {"path": str(path), "sha256": m.sha256(path), "bytes": 6}
            self.assertEqual(m.check_pin(pin), path)
            for bad in [{"path": str(path)}, {**pin, "sha256": "0"*64}, {**pin, "bytes": 7}]:
                with self.assertRaises(m.AnalysisError): m.check_pin(bad)
            prefix = path.with_suffix(".prefix"); prefix.write_bytes(path.read_bytes())
            with self.assertRaises(m.AnalysisError): m.check_pin({**pin, "path": str(prefix)})

    def test_validate_only_never_calls_effect_functions_and_leaves_all_source_bytes(self):
        with tempfile.TemporaryDirectory() as temp, tempfile.TemporaryDirectory(dir=REPO/"work") as work:
            base = Path(temp); source = base/"source"; source.mkdir(); output = base/"output"
            plan_path = disk_fixture(source); document = json.loads(plan_path.read_text())
            document["organoid"].update(raw_count_gate_passed=False, ready_for_inference=False, count_quantity="SCT_or_unknown")
            plan_path.write_text(json.dumps(document))
            original = {path: m.sha256(path) for path in source.iterdir()}
            with patch.object(m, "fit_negative_binomial", side_effect=AssertionError("No real fit")), patch.object(m, "paired_changes", side_effect=AssertionError("No effects")), patch.object(m, "log2_cpm", side_effect=AssertionError("No effects")):
                qc = m.run_analysis(plan_path, output, Path(work)/"private", validate_only=True)
            self.assertFalse(qc["effects_computed"]); self.assertEqual(qc["eligible_features"], 98)
            self.assertFalse((output/"organoid-primary.csv.gz").exists())
            self.assertTrue((output/"organoid-feature-universe.csv.gz").exists())
            for path, digest in original.items(): self.assertEqual(m.sha256(path), digest)

    def test_changed_software_sources_code_and_library_totals_fail_before_effects(self):
        with tempfile.TemporaryDirectory() as temp, tempfile.TemporaryDirectory(dir=REPO/"work") as work:
            base = Path(temp); source = base/"source"; source.mkdir(); plan_path = disk_fixture(source)
            original = json.loads(plan_path.read_text())
            for kind in ["counts", "features", "libraries", "count_quantity_qualification"]:
                bad = copy.deepcopy(original); bad["organoid"][kind]["sha256"] = "0"*64; plan_path.write_text(json.dumps(bad))
                with self.subTest(kind=kind), self.assertRaises(m.AnalysisError): m.run_analysis(plan_path, base/"output", Path(work)/"private", validate_only=True)
            bad = copy.deepcopy(original); bad["organoid"]["code_artifacts"][0]["sha256"] = "0"*64; plan_path.write_text(json.dumps(bad))
            with self.assertRaises(m.AnalysisError): m.run_analysis(plan_path, base/"output", Path(work)/"private", validate_only=True)
            bad = copy.deepcopy(original); sample_path = Path(bad["organoid"]["libraries"]["path"])
            table = m.read_table(sample_path); table.loc[0, "total_released_counts"] = "1"; table.to_csv(sample_path, sep="\t", index=False)
            bad["organoid"]["libraries"] = {"path": str(sample_path), "sha256": m.sha256(sample_path)}; plan_path.write_text(json.dumps(bad))
            with self.assertRaisesRegex(m.AnalysisError, "library sums"): m.run_analysis(plan_path, base/"output", Path(work)/"private", validate_only=True)

    def test_occupied_outputs_never_keep_stale_effects_as_a_new_success(self):
        with tempfile.TemporaryDirectory() as temp, tempfile.TemporaryDirectory(dir=REPO/"work") as work:
            base = Path(temp); source = base/"source"; source.mkdir(); plan_path = disk_fixture(source)
            output = base/"output"; output.mkdir()
            stale = output/"organoid-qc.json"; stale.write_text('{"status":"complete","effects_computed":true}')
            digest = m.sha256(stale)
            with self.assertRaisesRegex(m.AnalysisError, "must both be new"):
                m.run_analysis(plan_path, output, Path(work)/"private", validate_only=True)
            self.assertEqual(m.sha256(stale), digest)
            self.assertFalse((output/"organoid-validation.json").exists())
            self.assertFalse((Path(work)/"private").exists())

    def test_existing_destination_aliases_cannot_mutate_sources_or_publish_donors(self):
        with tempfile.TemporaryDirectory() as temp, tempfile.TemporaryDirectory(dir=REPO/"work") as work:
            base = Path(temp); source = base/"source"; source.mkdir(); plan_path = disk_fixture(source)
            count_path = source/"counts.tsv.gz"; original = m.sha256(count_path)
            output = base/"output"; output.mkdir()
            (output/"organoid-feature-universe.csv.gz").symlink_to(count_path)
            with self.assertRaises(m.AnalysisError): m.run_analysis(plan_path, output, Path(work)/"private", validate_only=True)
            self.assertEqual(m.sha256(count_path), original)
            private = Path(work)/"private2"; private.mkdir()
            leaked = base/"leaked-donors.csv"; (private/"selected-libraries.csv").symlink_to(leaked)
            with self.assertRaises(m.AnalysisError): m.run_analysis(plan_path, base/"output2", private, validate_only=True)
            self.assertFalse(leaked.exists())

    def test_exclusive_writers_reject_symlinks_hardlinks_and_existing_files(self):
        import os
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp); target = base/"source"; target.write_text("unchanged")
            digest = m.sha256(target)
            aliases = [base/"symbolic.csv", base/"hard.csv", base/"ordinary.csv"]
            aliases[0].symlink_to(target); os.link(target, aliases[1]); aliases[2].write_text("old")
            for alias in aliases:
                with self.subTest(alias=alias), self.assertRaises(FileExistsError): m._write_csv(pd.DataFrame({"g": [1]}), alias)
                with self.subTest(alias=alias), self.assertRaises(FileExistsError): m._write_json(alias, {"x": 1})
            zipped = base/"symbolic.csv.gz"; zipped.symlink_to(target)
            with self.assertRaises(FileExistsError): m._write_csv(pd.DataFrame({"g": [1]}), zipped)
            self.assertEqual(m.sha256(target), digest)

    def test_failed_fit_leaves_no_completion_artifact_and_cannot_mix_on_rerun(self):
        with tempfile.TemporaryDirectory() as temp, tempfile.TemporaryDirectory(dir=REPO/"work") as work:
            base = Path(temp); source = base/"source"; source.mkdir(); plan_path = disk_fixture(source)
            output = base/"output"; private = Path(work)/"private"
            with patch.object(m, "fit_negative_binomial", side_effect=RuntimeError("synthetic_fit_failure")):
                with self.assertRaisesRegex(RuntimeError, "synthetic_fit_failure"): m.run_analysis(plan_path, output, private)
            for name in ["organoid-primary.csv.gz", "organoid-qc.json", "organoid-output-manifest.json"]:
                self.assertFalse((output/name).exists())
            before = {path: m.sha256(path) for path in output.iterdir()}
            with self.assertRaises(m.AnalysisError): m.run_analysis(plan_path, output, private, validate_only=True)
            for path, digest in before.items(): self.assertEqual(m.sha256(path), digest)

    def test_synthetic_full_pipeline_coverage_privacy_and_exact_replay(self):
        with tempfile.TemporaryDirectory() as temp, tempfile.TemporaryDirectory(dir=REPO/"work") as work:
            base = Path(temp); source = base/"source"; source.mkdir(); plan_path = disk_fixture(source)
            qc = m.run_analysis(plan_path, base/"output", Path(work)/"private", m.sha256(plan_path))
            primary = pd.read_csv(base/"output/organoid-primary.csv.gz", keep_default_na=False)
            self.assertEqual(len(primary), 100); self.assertEqual(qc["eligible_features"], 98)
            self.assertEqual(primary.iloc[-1].status, "excluded_fixed_count_filter")
            self.assertEqual(primary.iloc[-1].pvalue, ""); self.assertEqual(primary.iloc[-1].literal_source_symbol, "NA")
            self.assertFalse(set(metadata_fixture().index) & set(primary.columns))
            within = pd.read_csv(base/"output/organoid-within-group.csv.gz")
            self.assertEqual(len(within), 200)
            for suffix in ["welch", "leave-one-donor-out", "allocations"]:
                self.assertEqual(len(pd.read_csv(base/f"output/organoid-{suffix}.csv.gz")), 100)
            self.assertEqual(len(pd.read_csv(Path(work)/"private/all-70-allocation-effects.csv.gz")), 70)
            self.assertEqual(len(pd.read_csv(Path(work)/"private/donor-paired-log2cpm-changes.csv.gz")), 8)
            manifest = json.loads((base/"output/organoid-output-manifest.json").read_text())
            for name, pin in manifest["outputs"].items(): self.assertEqual(m.sha256(base/"output"/name), pin["sha256"])
            # A second isolated synthetic fit also tests deterministic serialization.
            m.run_analysis(plan_path, base/"replay", Path(work)/"replay", m.sha256(plan_path))
            for name, pin in manifest["outputs"].items(): self.assertEqual(m.sha256(base/"replay"/name), pin["sha256"])


if __name__ == "__main__":
    unittest.main()
