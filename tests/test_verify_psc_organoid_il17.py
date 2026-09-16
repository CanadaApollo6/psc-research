"""Synthetic tests for the independent verifier. No real matrices/effects opened."""
from __future__ import annotations

import ast
import copy
import gzip
import importlib.metadata
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import optimize, stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts import verify_psc_organoid_il17 as v


def metadata() -> pd.DataFrame:
    rows = [{"library_id": donor + "_" + treatment, "donor_id": donor, "disease": group,
             "treatment": treatment, "n_cells": "200", "total_released_counts": "1"}
            for donor, group in v.DONORS.items() for treatment in ["CNT", "IL17"]]
    rows[-1]["n_cells"] = str(46343 - 15*200)
    return pd.DataFrame(rows).set_index("library_id", drop=False)


def synthetic_counts(n: int = 112) -> np.ndarray:
    rng = np.random.default_rng(785262)
    mu = rng.uniform(200, 4000, n)[:, None] * np.exp(rng.normal(0, .15, 16))
    mu[0, np.arange(1, 8, 2)] *= 4
    if n > 1:
        mu[1, np.arange(9, 16, 2)] *= 3
    return rng.negative_binomial(20, 20/(20+mu))


def current_versions() -> dict:
    return {name: ".".join(map(str, sys.version_info[:3])) if name == "python" else importlib.metadata.version(name)
            for name in v.SOFTWARE}


def disk_fixture(directory: Path) -> Path:
    values = synthetic_counts()
    values[2, 8:] = 0                 # nonPSC has no response information
    values[3, np.arange(0, 8, 2)] = 0  # PSC response has an infinite-MLE boundary
    values[4, :2] = 0                 # zero-total individual donor pair, flagged
    values[-1] = 0; values[-2] = 1
    meta = metadata(); ids = [f"synthetic_{i:03d}" for i in range(len(values))]
    count_path = directory / "counts.tsv.gz"
    pd.DataFrame(values, index=ids, columns=meta.library_id).rename_axis("source_feature_id").reset_index().to_csv(count_path, sep="\t", index=False)
    feature_path = directory / "features.tsv"
    pd.DataFrame({"source_feature_id": ids[::-1], "literal_symbol": ["NA"]*len(ids),
                  "source_version": ["unchanged.7"]*len(ids)}).to_csv(feature_path, sep="\t", index=False)
    meta["total_released_counts"] = values.sum(axis=0).astype(str)
    libraries_path = directory / "libraries.tsv"
    meta.iloc[::-1].to_csv(libraries_path, sep="\t", index=False)
    qual_path = directory / "quantity.json"
    qual_path.write_text('{"synthetic_only":true,"raw_count_gate_passed":true}')
    pin = lambda p: {"path": str(p), "sha256": v.digest(p), "bytes": p.stat().st_size}
    arm = {"series": "GSE239283", "population": "all_published_retained_organoid_cells", "donor_disease": dict(v.DONORS),
           "expected_library_count": 16, "expected_cell_count": 46343, "expected_feature_count": len(values),
           "feature_namespace": "synthetic_ID", "feature_universe_description": "synthetic_complete_universe",
           "n_cpus": 1, "ready_for_inference": True, "raw_count_gate_passed": True, "count_quantity": "raw_RNA_counts",
           "filter": dict(v.FILTER), "model": dict(v.MODEL), "sensitivity": dict(v.SENSITIVITY),
           "software": current_versions(), "counts": pin(count_path), "features": pin(feature_path),
           "libraries": pin(libraries_path), "count_quantity_qualification": pin(qual_path),
           "source_artifacts": [{**pin(count_path), "url": "https://example.org/synthetic-only",
                                 "retrieved_at_utc": "2026-01-01T00:00:00Z"}],
           "code_artifacts": [pin(REPO / "scripts/analyze_psc_organoid_il17.py"), pin(Path(v.__file__).resolve())]}
    plan = {"status": "frozen", "frozen_at_utc": "2026-01-02T00:00:00Z", "organoid": arm}
    path = directory / "plan.json"; path.write_text(json.dumps(plan))
    return path


class TestIndependentMath(unittest.TestCase):
    def test_no_analysis_or_pydeseq_helpers_imported(self):
        tree = ast.parse(Path(v.__file__).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [x.name for x in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            else:
                continue
            self.assertFalse(any("analyze_psc_organoid" in name or "pydeseq2" in name for name in names))

    def test_literal_sources_ragged_duplicate_axes_and_exact_integer_guards(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.tsv"
            path.write_text("ID\tlabel\ng1\tNA\ng2\t\n")
            frame = v.read_table(path, "\t")
            self.assertEqual(frame.label.tolist(), ["NA", ""])
            for content in ["ID\tID\ng1\t2\n", "ID\tlabel\ng1\n", "ID\tlabel\ng1\t2\t3\n"]:
                path.write_text(content)
                with self.assertRaises(v.VerificationError): v.read_table(path, "\t")
        for token in ["1.0", "1e3", "-1e-400", "1e-400", " 1", "+1", "NA", "nan", "inf", "9007199254740992"]:
            with self.subTest(token=token), self.assertRaises(v.VerificationError):
                v.integer_tokens(np.array([[token]]))
        np.testing.assert_array_equal(v.integer_tokens(np.array([["00010", "0"]])), [[10, 0]])

    def test_design_rank_leverage_alias_direction_and_units(self):
        m = metadata(); x, names, donors, leverage = v.make_design(m)
        self.assertEqual(x.shape, (16, 10)); self.assertEqual(np.linalg.matrix_rank(x), 10)
        np.testing.assert_allclose(leverage, .625, atol=1e-13)
        self.assertEqual(len(donors), 8)
        extra = np.column_stack([x, (m.disease == "PSC").astype(float)])
        self.assertEqual(np.linalg.matrix_rank(extra), 10)
        beta = np.zeros(10); beta[-2:] = [-2, 3]
        fitted = x @ beta
        for donor in donors:
            mask = m.donor_id == donor
            self.assertEqual(fitted[mask][1] - fitted[mask][0], 1 if donor.startswith("PSC") else -2)
        changed = m.copy(); changed.n_cells = "500000"
        np.testing.assert_array_equal(v.make_design(changed)[0], x)
        for column, value in [("disease", "nonPSC"), ("treatment", "IL17"), ("donor_id", "unknown")]:
            bad = m.copy(); bad.iloc[0, bad.columns.get_loc(column)] = value
            with self.assertRaises(v.VerificationError): v.make_design(bad)

    def test_filter_boundary_and_even_log_median_normalization(self):
        y = np.full((6, 16), 20)
        y[0] = np.arange(1, 17) * 10
        y[1] = np.arange(16, 0, -1) * 25
        y[2] = 0; y[2, :4] = 10
        y[3] = 0; y[3, :3] = 10
        y[4] = 1; y[5] = 0
        eligible, info, sf, basis = v.filter_and_ratio(y)
        np.testing.assert_array_equal(eligible, [True, True, True, False, False, False])
        self.assertEqual(basis, 2)
        ratio = y[:2] / np.exp(np.log(y[:2]).mean(axis=1))[:, None]
        np.testing.assert_allclose(sf, np.sqrt(ratio[0] * ratio[1]), rtol=1e-14)
        self.assertNotAlmostEqual(sf[0], np.median(ratio[:, 0]))
        self.assertEqual(info["libraries_count_ge_threshold"][2], 4)
        self.assertEqual(info["libraries_count_ge_threshold"][3], 3)
        y[:2, 0] = 0
        with self.assertRaises(v.VerificationError): v.filter_and_ratio(y)

    def test_whole_family_BH_BY_preserves_missing_and_zero(self):
        p = np.array([0, .002, np.nan, .11, .02, 1, .002, np.nan])
        original = p.copy()
        for by, method in [(False, "fdr_bh"), (True, "fdr_by")]:
            np.testing.assert_allclose(v.adjust(p, by), multipletests(np.nan_to_num(p, nan=1), method=method)[1], rtol=1e-14)
            self.assertEqual(v.adjust(p, by)[2], 1)
        np.testing.assert_array_equal(p, original)
        self.assertGreater(v.adjust(p)[1], multipletests(p[np.isfinite(p)], method="fdr_bh")[1][1])
        self.assertEqual(len(v.adjust([])), 0)
        for bad in [[np.inf], [-.1], [1.1], [[.1]]]:
            with self.assertRaises(v.VerificationError): v.adjust(bad)

    def test_logCPM_uses_all_raw_features_and_exact_donor_pairing(self):
        m = metadata()
        y = synthetic_counts(12)
        selected = np.array([True]*4 + [False]*8)
        _, _, donors, _ = v.make_design(m)
        actual = v.changes_from_counts(y, m, donors, selected)
        expression = np.log2(y * (1e6 / y.sum(axis=0)) + 1)
        expected = expression[:4, 1::2] - expression[:4, ::2]
        np.testing.assert_allclose(actual, expected, atol=3e-14)
        wrong = v.changes_from_counts(y[:4], m, donors, np.ones(4, bool))
        self.assertGreater(np.max(np.abs(wrong - actual)), 1e-3)
        order = np.arange(16)[::-1]
        np.testing.assert_allclose(v.changes_from_counts(y[:, order], m.iloc[order], donors, selected), actual)

    def test_welch_actual_df_and_student_intervals(self):
        y = np.array([[1, 2, 5, 8, -1, 0, .5, 1], [7, 8, 9, 10, 1, 20, -10, 40],
                      [1, 1, 1, 1, 0, 1, 2, 3], [1, 2, 3, 4, 0, 1, 2, 3]])
        result = v.donor_sensitivity(y, np.array([1]*4 + [0]*4, bool))["welch"]
        reference = stats.ttest_ind(y[:, :4], y[:, 4:], axis=1, equal_var=False)
        np.testing.assert_allclose(result["t"], reference.statistic)
        np.testing.assert_allclose(result["pvalue"], reference.pvalue)
        np.testing.assert_allclose(result["df"], reference.df)
        np.testing.assert_allclose(result["ci95_lower"], reference.confidence_interval().low)
        np.testing.assert_allclose(result["ci95_upper"], reference.confidence_interval().high)
        self.assertEqual(result["df"][2], 3); self.assertEqual(result["df"][3], 6)
        normal_low = result["effect_log2_cpm"] - stats.norm.ppf(.975) * result["se"]
        self.assertTrue(np.all(result["ci95_lower"] < normal_low))

    def test_degenerate_welch_missingness_and_no_case_deletion(self):
        y = np.array([[1]*4 + [0]*4, [0]*8], dtype=float)
        result = v.donor_sensitivity(y, np.array([1]*4 + [0]*4, bool))
        self.assertTrue(np.isnan(result["welch"]["pvalue"]).all())
        np.testing.assert_array_equal(result["welch"]["q_bh_sensitivity"], 1)
        self.assertEqual(result["allocation_summary"]["descriptive_tail_fraction"][1], 1)
        y[0, 0] = np.nan
        with self.assertRaises(v.VerificationError): v.donor_sensitivity(y, np.array([1]*4 + [0]*4, bool))

    def test_all_leaveouts_whole_donors_and_allocation_ties_complements(self):
        y = np.array([[10, 11, 12, 13, 0, 1, 2, 3], [0]*8, [10, 0, 0, 0, 1, 1, 1, 1]], float)
        groups = np.array([1]*4 + [0]*4, bool)
        result = v.donor_sensitivity(y, groups)
        self.assertEqual(result["leaveout"].shape, (8, 3))
        self.assertLess(result["leaveout"][0, 2], 0)
        for i in range(8):
            keep = np.arange(8) != i
            expected = y[:, groups & keep].mean(axis=1) - y[:, ~groups & keep].mean(axis=1)
            np.testing.assert_allclose(result["leaveout"][i], expected)
        tuples = result["assignment_indices"]
        self.assertEqual(len(tuples), 70); self.assertEqual(len(set(tuples)), 70)
        self.assertEqual(result["allocation_summary"]["n_as_or_more_extreme"][0], 2)
        self.assertEqual(result["allocation_summary"]["descriptive_tail_fraction"][0], 2/70)
        self.assertEqual(result["allocation_summary"]["descriptive_tail_fraction"][1], 1)
        for i, chosen in enumerate(tuples):
            complement = tuple(j for j in range(8) if j not in chosen)
            np.testing.assert_allclose(result["allocations"][i], -result["allocations"][tuples.index(complement)])
        shifted = v.donor_sensitivity(y + 1e-13, groups)
        np.testing.assert_array_equal(shifted["allocation_summary"]["n_as_or_more_extreme"], result["allocation_summary"]["n_as_or_more_extreme"])
        self.assertNotIn("q_bh", result["allocation_summary"])

    def test_unstudentized_allocations_not_welch_statistics(self):
        y = np.array([[0, 2, 3, 100, 10, 11, 12, 13]], float)
        result = v.donor_sensitivity(y, np.array([1]*4 + [0]*4, bool))
        for i, chosen in enumerate(result["assignment_indices"]):
            other = [j for j in range(8) if j not in chosen]
            self.assertAlmostEqual(result["allocations"][i, 0], y[0, list(chosen)].mean() - y[0, other].mean())
        self.assertNotAlmostEqual(result["allocations"][0, 0], result["welch"]["t"][0])

    def test_covariance_three_contrasts_and_full_cross_term(self):
        x = v.make_design(metadata())[0]
        beta = np.zeros((1, 10)); beta[0, 0] = np.log(200); beta[0, -2:] = [-.2, .7]
        factors = np.exp(np.linspace(-.3, .3, 16)); alpha = np.array([.1])
        means = factors * np.exp(x @ beta[0])
        information = x.T @ np.diag(means / (1 + .1*means)) @ x
        h = np.linalg.inv(information + 1e-6*np.eye(10))
        expected = h @ information @ h
        np.testing.assert_allclose(v.conditional_covariance(x, means, .1), expected, rtol=1e-12, atol=1e-13)
        reconstructed = v.nb_reconstruction(x, factors, means[None], beta, alpha)
        tables = reconstructed["contrasts"]
        self.assertAlmostEqual(tables["PSC_treatment"]["log2FoldChange"][0], .5/np.log(2))
        self.assertAlmostEqual(tables["interaction"]["log2FoldChange"][0], .7/np.log(2))
        expected_var = expected[-2, -2] + expected[-1, -1] + 2*expected[-2, -1]
        self.assertAlmostEqual(tables["PSC_treatment"]["lfcSE"][0]**2, expected_var/np.log(2)**2)
        wrong = tables["nonPSC_treatment"]["lfcSE"][0]**2 + tables["interaction"]["lfcSE"][0]**2
        self.assertNotAlmostEqual(wrong, tables["PSC_treatment"]["lfcSE"][0]**2)
        self.assertAlmostEqual(tables["interaction"]["pvalue"][0], 2*stats.norm.sf(abs(tables["interaction"]["stat"][0])))

    def test_independent_NB_objective_gradient_hessian_and_optimization(self):
        x = v.make_design(metadata())[0]; y = synthetic_counts(1)[0].astype(float)
        factors = np.exp(np.linspace(-.4, .4, 16)); alpha = .06
        beta = np.linalg.lstsq(x, np.log(y/factors), rcond=None)[0]
        loss, gradient, hessian = v.nb_objective(beta, y, factors, x, alpha)
        numeric_gradient = optimize._numdiff.approx_derivative(lambda b: np.array([v.nb_objective(b, y, factors, x, alpha)[0]]), beta).ravel()
        numeric_hessian = optimize._numdiff.approx_derivative(lambda b: v.nb_objective(b, y, factors, x, alpha)[1], beta)
        np.testing.assert_allclose(gradient, numeric_gradient, atol=2e-6, rtol=2e-6)
        np.testing.assert_allclose(hessian, numeric_hessian, atol=3e-6, rtol=2e-6)
        self.assertGreater(np.linalg.eigvalsh(hessian).min(), 0)
        fitted = v.optimize_conditional(y, factors, x, alpha)
        self.assertLessEqual(fitted["objective"], loss + 1e-8)
        self.assertLess(fitted["newton_decrement"], 1e-5)
        beta2 = fitted["beta"] + np.linspace(-.1, .1, 10)
        mu1, mu2 = factors*np.exp(x@fitted["beta"]), factors*np.exp(x@beta2)
        nll_difference = -stats.nbinom.logpmf(y, 1/alpha, 1/(1+alpha*mu2)).sum() + stats.nbinom.logpmf(y, 1/alpha, 1/(1+alpha*mu1)).sum()
        ridge_difference = 1e-6*(np.dot(beta2, beta2)-np.dot(fitted["beta"], fitted["beta"]))/2
        self.assertAlmostEqual(v.nb_objective(beta2, y, factors, x, alpha)[0] - fitted["objective"], nll_difference+ridge_difference, places=8)

    def test_zero_group_and_arm_mask_only_affected_contrasts_keep_fixed_family(self):
        counts = synthetic_counts(6)
        counts[0, 8:] = 0                    # entirely absent nonPSC group
        counts[1, np.arange(0, 8, 2)] = 0     # PSC CNT arm absent
        counts[2, :2] = 0                    # individual nuisance boundary only
        counts[3, :8] = 0                    # entirely absent PSC group
        counts[4, np.arange(9, 16, 2)] = 0    # nonPSC IL17 arm absent
        sums, masks = v.count_boundary_support(counts, metadata())
        np.testing.assert_array_equal(masks["interaction"]["supported"], [False, False, True, False, False, True])
        np.testing.assert_array_equal(masks["PSC_treatment"]["supported"], [True, False, True, False, True, True])
        np.testing.assert_array_equal(masks["nonPSC_treatment"]["supported"], [False, True, True, True, False, True])
        self.assertEqual(sums["n_zero_total_donor_pairs_PSC"][2], 1)
        self.assertEqual(masks["interaction"]["zero_required_count_arms"][0], "nonPSC_CNT|nonPSC_IL17")
        native = {label: {"log2FoldChange": np.ones(6), "lfcSE": np.ones(6), "stat": np.ones(6),
                          "pvalue": np.full(6, .001), "pydeseq2_finite_family_padj": np.full(6, .001)} for label in masks}
        adjusted = v.apply_boundary_support(native, masks)
        for label in adjusted:
            supported = masks[label]["supported"]; table = adjusted[label]
            self.assertTrue(np.isnan(table["pvalue"][~supported]).all())
            self.assertTrue(np.isnan(table["ci95_lower_log2fc"][~supported]).all())
            np.testing.assert_array_equal(table["q_bh"][~supported], 1)
            np.testing.assert_array_equal(table["pydeseq2_unmasked_pvalue"], .001)
            np.testing.assert_allclose(table["q_bh"], multipletests(np.where(supported, .001, 1), method="fdr_bh")[1])
            np.testing.assert_array_equal(native[label]["pvalue"], .001)

    def test_conditional_panel_uses_ID_hash_and_inputs_not_effects_or_row_order(self):
        y = synthetic_counts(50); ids = [f"ID{i}" for i in range(50)]; eligible = np.ones(50, bool)
        y[0, 0] = 0; eligible[1] = False
        first = v.refit_indices(y, ids, eligible)
        self.assertEqual(len(first), 16); self.assertNotIn(0, first); self.assertNotIn(1, first)
        reverse = np.arange(50)[::-1]
        second = v.refit_indices(y[reverse], [ids[i] for i in reverse], eligible[reverse])
        self.assertEqual([ids[i] for i in first], [ids[reverse[i]] for i in second])

    def test_nonfinite_native_fit_keeps_rows_and_correction_family(self):
        x = v.make_design(metadata())[0]
        beta = np.zeros((3, 10)); beta[:, 0] = 5; beta[1, 0] = np.nan
        tables = v.nb_reconstruction(x, np.ones(16), np.full((3, 16), 150), beta, np.array([.1, .1, np.nan]))["contrasts"]
        for table in tables.values():
            self.assertEqual(len(table["pvalue"]), 3)
            self.assertTrue(np.isnan(table["pvalue"][[1, 2]]).all())
            np.testing.assert_array_equal(table["q_bh"][[1, 2]], 1)
            self.assertTrue(np.isnan(table["pydeseq2_finite_family_padj"][[1, 2]]).all())

    def test_numeric_audit_rejects_missing_P_drop_and_tiny_P_corruption(self):
        audit = v.Audit()
        with self.assertRaises(v.VerificationError): audit.close("lost_NA", [0], [np.nan])
        with self.assertRaises(v.VerificationError): audit.p("tiny_p", [1e-50], [1e-60])
        audit.p("exact_zero", [0, np.nan], [0, np.nan])


class TestFrozenNativeInterface(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # The analysis module is used only to PRODUCE a synthetic interoperability
        # fixture. None of its mathematical helpers supply expected values.
        from scripts import analyze_psc_organoid_il17 as analysis
        cls.storage = tempfile.TemporaryDirectory(prefix="organoid-verifier-synthetic-")
        cls.work_storage = tempfile.TemporaryDirectory(prefix="verifier-synthetic-", dir=REPO / "work")
        cls.directory = Path(cls.storage.name); cls.work = Path(cls.work_storage.name) / "native"
        cls.plan = disk_fixture(cls.directory)
        cls.public = cls.directory / "public"
        analysis.run_analysis(cls.plan, cls.public, cls.work, v.digest(cls.plan))

    @classmethod
    def tearDownClass(cls):
        cls.storage.cleanup(); cls.work_storage.cleanup()

    def test_full_independent_reconstruction_and_conditional_fits(self):
        result = v.verify_run(self.plan, v.digest(self.plan), self.public, self.work)
        self.assertEqual(result["status"], "independent_numerical_verification_passed")
        self.assertFalse(result["independent_analysis_helpers_imported"])
        self.assertFalse(result["full_dispersion_model_refitted"])
        self.assertEqual(result["counts"]["features"], 112)
        self.assertEqual(result["counts"]["eligible"], 110)
        conditional = result["conditional_fixed_dispersion_NB"]["results"]
        self.assertEqual(len(conditional), 16)
        self.assertTrue(all(row["status"] == "conditional_optimization_agrees" for row in conditional))
        self.assertIn(result["dispersion_fit_type"], {"parametric", "mean"})
        self.assertGreater(len(result["checks"]), 250)
        self.assertFalse((self.public / "verification.json").exists())

    def test_native_finite_boundary_values_are_preserved_but_not_tested(self):
        primary = v.read_table(self.public / "organoid-primary.csv.gz").set_index("source_feature_id")
        within = v.read_table(self.public / "organoid-within-group.csv.gz")
        for feature in ["synthetic_002", "synthetic_003"]:
            self.assertNotEqual(primary.loc[feature, "pydeseq2_unmasked_pvalue"], "")
            self.assertEqual(primary.loc[feature, "pvalue"], "")
            self.assertEqual(float(primary.loc[feature, "q_bh"]), 1)
            self.assertEqual(primary.loc[feature, "status"], "unavailable_boundary_or_unidentified_group_log_response")
        for feature, unaffected in [("synthetic_002", "PSC_treatment"), ("synthetic_003", "nonPSC_treatment")]:
            row = within.loc[(within.source_feature_id == feature) & (within.contrast == unaffected)].iloc[0]
            self.assertEqual(row.count_supported_contrast, "True")
            self.assertNotEqual(row.pvalue, "")
        self.assertEqual(primary.loc["synthetic_004", "count_supported_contrast"], "True")
        self.assertEqual(float(primary.loc["synthetic_004", "n_zero_total_donor_pairs_PSC"]), 1)

    def test_structure_only_does_not_read_or_compute_effects(self):
        with patch.object(v, "nb_reconstruction", side_effect=AssertionError("effect")), patch.object(v, "changes_from_counts", side_effect=AssertionError("effect")):
            result = v.verify_run(self.plan, v.digest(self.plan), self.public, self.work, structure_only=True)
        self.assertEqual(result["status"], "structural_verification_passed")
        self.assertNotIn("conditional_fixed_dispersion_NB", result)

    def test_plan_source_code_software_and_freeze_guards(self):
        plan = json.loads(self.plan.read_text())
        changes = [lambda p: p.update(status="draft"), lambda p: p.update(frozen_at_utc="2099-01-01T00:00:00Z"),
                   lambda p: p["organoid"].update(raw_count_gate_passed=False),
                   lambda p: p["organoid"].update(count_quantity="SCT_corrected"),
                   lambda p: p["organoid"]["filter"].update(min_libraries=3),
                   lambda p: p["organoid"]["model"].update(cooks_filter=True),
                   lambda p: p["organoid"]["sensitivity"].update(allocation="Welch_studentized"),
                   lambda p: p["organoid"]["software"].update(pydeseq2="0.0"),
                   lambda p: p["organoid"]["counts"].update(sha256="0"*64),
                   lambda p: p["organoid"]["code_artifacts"].pop(),
                   lambda p: p["organoid"]["source_artifacts"][0].update(retrieved_at_utc="2026-01-03T00:00:00Z")]
        for change in changes:
            bad = copy.deepcopy(plan); change(bad)
            with self.subTest(change=change), self.assertRaises(v.VerificationError): v.validate_frozen_plan(bad)
        with self.assertRaises(v.VerificationError): v.verify_run(self.plan, "0"*64, self.public, self.work)

    def test_builtin_mean_fallback_permitted_other_fallback_rejected(self):
        qc_path = self.public / "organoid-qc.json"
        table_path = self.public / "organoid-primary.csv.gz"
        manifest_path = self.public / "organoid-output-manifest.json"
        old = {p: p.read_bytes() for p in [qc_path, table_path, manifest_path]}
        try:
            for trend in ["mean", "unapproved_loess"]:
                qc = json.loads(old[qc_path])
                qc["negative_binomial"]["dispersion_fit_type"] = trend
                qc["negative_binomial"]["dispersion_fallback_used"] = trend == "mean"
                qc_path.write_text(json.dumps(qc))
                table = v.read_table(table_path)
                table.loc[table.eligible == "True", "dispersion_fit_type"] = trend
                table.to_csv(table_path, index=False)
                manifest = json.loads(old[manifest_path])
                for path in [qc_path, table_path]:
                    manifest["outputs"][path.name] = {"sha256": v.digest(path), "bytes": path.stat().st_size}
                manifest_path.write_text(json.dumps(manifest))
                if trend == "mean":
                    # This checks acceptance of the allowed contract, not a second
                    # independent dispersion fit or evidence that relabeling is valid.
                    self.assertEqual(v.verify_run(self.plan, v.digest(self.plan), self.public, self.work)["dispersion_fit_type"], "mean")
                else:
                    with self.assertRaisesRegex(v.VerificationError, "Unapproved dispersion fallback"):
                        v.verify_run(self.plan, v.digest(self.plan), self.public, self.work)
        finally:
            for path, content in old.items(): path.write_bytes(content)

    def test_receipt_overwrite_and_protected_destination_refused_before_verification(self):
        arguments = ["verify", "--plan", str(self.plan), "--plan-sha256", v.digest(self.plan),
                     "--output-dir", str(self.public), "--work-dir", str(self.work)]
        for destination in [self.plan, self.public / "new-receipt.json", self.work / "new-receipt.json"]:
            with self.subTest(destination=destination), patch.object(sys, "argv", arguments + ["--verification-json", str(destination)]), patch.object(v, "verify_run") as run:
                with self.assertRaises(v.VerificationError): v.main()
                run.assert_not_called()
            if destination != self.plan: self.assertFalse(destination.exists())

    def test_corrupt_native_change_detected_without_an_output_manifest(self):
        path = self.work / "donor-paired-log2cpm-changes.csv.gz"
        original = path.read_bytes()
        try:
            frame = pd.read_csv(path, keep_default_na=False)
            frame.iloc[0, 1] = float(frame.iloc[0, 1]) + .5
            frame.to_csv(path, index=False)
            with self.assertRaisesRegex(v.VerificationError, "donor_changes:all_values"):
                v.verify_run(self.plan, v.digest(self.plan), self.public, self.work)
        finally:
            path.write_bytes(original)

    def test_corrupt_public_hash_detected(self):
        path = self.public / "organoid-primary.csv.gz"
        original = path.read_bytes()
        try:
            path.write_bytes(original + b"corruption")
            with self.assertRaisesRegex(v.VerificationError, "hash mismatch"):
                v.verify_run(self.plan, v.digest(self.plan), self.public, self.work)
        finally:
            path.write_bytes(original)


if __name__ == "__main__":
    unittest.main()
