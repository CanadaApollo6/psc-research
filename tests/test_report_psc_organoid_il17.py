"""Tests for post-fit aggregation only; no primary model modifications."""
import copy
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from scripts import report_psc_organoid_il17 as m


def fixture():
    rows = []
    for i in range(4):
        eligible = i < 3
        row = dict(zip(m.ANNOTATION, [f"G{i+1}", f"label{i+1}", "Gene Expression", str(eligible),
                                     "4" if eligible else "0", "16" if eligible else "0",
                                     str(not eligible), "100" if eligible else "0"]))
        row.update(log2FoldChange=["2", "-1", "", ""][i], pvalue=[".001", ".02", "", ""][i],
                   q_bh=[".01", ".04", "1", ""][i], q_by=[".03", ".08", "1", ""][i],
                   lfcSE=".4", stat="2", ci95_lower_log2fc="-2", ci95_upper_log2fc="3",
                   status="ok" if i < 2 else "unavailable", count_supported_contrast=str(i < 2),
                   zero_required_count_arms="nonPSC_CNT" if i == 2 else "",
                   genewise_converged="False" if i == 1 else "True", MAP_converged="False" if i == 1 else "True",
                   LFC_converged="True", n_libraries_cooks_above_diagnostic_cutoff="0", n_nonfinite_cooks="0",
                   cooks_automatic_exclusion="False", n_zero_total_donor_pairs_PSC="0",
                   n_zero_total_donor_pairs_nonPSC="1" if i == 2 else "0",
                   pydeseq2_unmasked_log2FoldChange="2", pydeseq2_unmasked_lfcSE=".4",
                   pydeseq2_unmasked_stat="5", pydeseq2_unmasked_pvalue=".001",
                   pydeseq2_unmasked_finite_family_padj=".01")
        rows.append(row)
    p = pd.DataFrame(rows, dtype=str)
    within = []
    for group in ("nonPSC", "PSC"):
        part = p.copy()
        part["contrast"] = group + "_treatment"
        if group == "PSC":
            part.loc[2, ["log2FoldChange", "pvalue", "q_bh", "q_by", "count_supported_contrast", "status"]] = [".1", ".4", ".8", ".9", "True", "ok"]
        within.append(part)
    w = p[m.ANNOTATION].copy()
    values = {"mean_change_PSC": ["2", ".5", ".2", ""], "mean_change_nonPSC": ["1", "1", "0", ""],
              "effect_log2_cpm": ["1", "-.5", ".2", ""], "pvalue": [".1", ".2", ".8", ""],
              "q_bh_sensitivity": [".8", ".8", ".8", ""], "q_by_sensitivity": ["1", "1", "1", ""]}
    for key, val in values.items():
        w[key] = val
    for key, val in (("se", ".5"), ("df", "4"), ("ci95_lower", "-1"), ("ci95_upper", "2")):
        w[key] = val
    leave = p[m.ANNOTATION].copy()
    for key, val in {"full_effect_log2_cpm": "1", "minimum_leaveout_effect": ".1", "maximum_leaveout_effect": "2",
                     "n_estimable_omissions": "8", "n_same_nonzero_sign": "8", "n_opposite_sign": "0", "n_zero_effect": "0"}.items():
        leave[key] = val
    leave.loc[1, ["n_same_nonzero_sign", "n_opposite_sign"]] = ["6", "2"]
    allocations = p[m.ANNOTATION].copy()
    allocations["n_allocations"] = "70"
    allocations["n_as_or_more_extreme"] = ["2", "4", "70", ""]
    allocations["descriptive_tail_fraction"] = [str(2/70), str(4/70), "1", ""]
    tables = {"primary": p, "within": pd.concat(within, ignore_index=True), "welch": w,
              "leaveout": leave, "allocations": allocations, "universe": p[m.ANNOTATION].copy()}
    qc = {"released_features": 4, "eligible_features": 3, "ratio_basis_positive_in_all_16_libraries": 2,
          "donors": 8, "libraries": 16, "released_cells": 100, "plan_sha256": "test", "fit_warning_count": 0,
          "negative_binomial": {"dispersion_fit_type": "parametric"}, "effects_computed": True}
    verify = {"status": "independent_numerical_verification_passed", "plan_sha256": "test", "checks": ["a"],
              "full_dispersion_model_refitted": False,
              "conditional_fixed_dispersion_NB": {"results": [
                  {"status": "conditional_optimization_agrees", "optimizer_success": False,
                   "maximum_contrast_difference_in_native_SE": .001, "twice_objective_improvement": .0001, "newton_decrement": 1e-8},
                  {"status": "unavailable_outside_scope"}]}}
    return tables, qc, verify


class TestOrganoidAggregateReport(unittest.TestCase):
    def setUp(self):
        self.tables, self.qc, self.verification = fixture()

    def test_complete_family_and_source_rows(self):
        summary, exports = m.summarize(self.tables, self.qc, self.verification)
        self.assertEqual(summary["feature_coverage"]["released"], 4)
        self.assertEqual(summary["interaction"]["family_size"], 3)
        self.assertEqual(summary["interaction"]["finite_pvalues"], 2)
        self.assertEqual(len(exports["bh-candidates-with-diagnostics.csv.gz"]), 2)

    def test_bh_by_signs_are_interactions_not_induction(self):
        summary, _ = m.summarize(self.tables, self.qc, self.verification)
        self.assertEqual(summary["interaction"]["BH"], {"count": 2, "positive": 1, "negative": 1})
        self.assertEqual(summary["interaction"]["BY"]["count"], 1)

    def test_flags_union_not_double_counted(self):
        summary, exports = m.summarize(self.tables, self.qc, self.verification)
        self.assertEqual(summary["diagnostics"]["any_native_nonconvergence"], 1)
        self.assertEqual(summary["diagnostics"]["genewise_nonconverged"], 1)
        self.assertEqual(summary["diagnostics"]["MAP_nonconverged"], 1)
        self.assertEqual(len(exports["native-flagged-features.csv.gz"]), 1)

    def test_masks_preserve_other_group_and_raw_values(self):
        summary, exports = m.summarize(self.tables, self.qc, self.verification)
        self.assertEqual(summary["diagnostics"]["boundary_masks_by_contrast"], {"interaction": 1, "nonPSC_treatment": 1})
        mask = exports["masked-contrasts.csv.gz"]
        self.assertEqual(set(mask["q_bh"]), {"1"})
        self.assertEqual(set(mask["pydeseq2_unmasked_pvalue"]), {".001"})
        self.assertEqual(set(mask["pvalue"]), {""})

    def test_sensitivities_do_not_replace_frozen_family(self):
        summary, _ = m.summarize(self.tables, self.qc, self.verification)
        self.assertEqual(summary["BH_candidate_sensitivities"]["all_eight_leaveout_signs_match_full_paired_change"], 1)
        self.assertEqual(summary["BH_candidate_sensitivities"]["same_NB_and_Welch_nonzero_direction"], 2)
        self.assertEqual(summary["welch_sensitivity"]["BH"]["count"], 0)
        self.assertEqual(summary["interaction"]["BH"]["count"], 2)

    def test_conditional_failure_status_is_visible(self):
        summary, _ = m.summarize(self.tables, self.qc, self.verification)
        self.assertEqual(summary["verification"]["conditional_optimizer_success_false"], 1)
        self.assertEqual(summary["verification"]["conditional_panel_status"]["unavailable_outside_scope"], 1)
        self.assertFalse(summary["verification"]["full_dispersion_model_refitted"])

    def test_duplicate_id_rejected(self):
        self.tables["primary"].loc[1, m.ID] = "G1"
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            m.summarize(self.tables, self.qc, self.verification)

    def test_annotation_changes_rejected(self):
        self.tables["welch"].loc[0, "source_gene_label"] = "changed"
        with self.assertRaisesRegex(ValueError, "Annotation"):
            m.summarize(self.tables, self.qc, self.verification)

    def test_order_and_missing_feature_rejected(self):
        self.tables["allocations"] = self.tables["allocations"].iloc[::-1]
        with self.assertRaisesRegex(ValueError, "order"):
            m.summarize(self.tables, self.qc, self.verification)

    def test_malformed_boolean_not_truthy(self):
        with self.assertRaisesRegex(ValueError, "boolean"):
            m.truth(pd.Series(["False", "unexpected"]))
        self.assertFalse(m.truth(pd.Series(["False"])).iloc[0])

    def test_bad_numeric_not_silently_coerced(self):
        with self.assertRaises(ValueError):
            m.number(pd.Series(["broken"]))

    def test_missing_p_requires_both_correction_values_one(self):
        self.tables["primary"].loc[2, "q_by"] = ".9"
        with self.assertRaisesRegex(ValueError, "correction-only"):
            m.summarize(self.tables, self.qc, self.verification)

    def test_empty_candidate_family_retained(self):
        self.tables["primary"]["q_bh"] = ["1", "1", "1", ""]
        self.tables["primary"]["q_by"] = ["1", "1", "1", ""]
        summary, exports = m.summarize(self.tables, self.qc, self.verification)
        self.assertEqual(summary["interaction"]["BH"]["count"], 0)
        self.assertIsNone(summary["BH_candidate_sensitivities"]["allocation_fraction_min"])
        self.assertTrue(exports["bh-candidates-with-diagnostics.csv.gz"].empty)

    def test_new_destination_and_symlink_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            inputs = base / "inputs"
            inputs.mkdir()
            with self.assertRaisesRegex(ValueError, "exists"):
                m.destinations(inputs, base / "figure", inputs)
            alias = base / "alias"
            alias.symlink_to(inputs, target_is_directory=True)
            with self.assertRaisesRegex(ValueError, "Protected"):
                m.destinations(base / "out", alias / "figure", inputs)
            m.destinations(base / "new", base / "fig", inputs)

    def test_lossless_deterministic_gzip_and_exclusive_write(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "a.gz", Path(tmp) / "b.gz"
            m.csv_gz(self.tables["primary"], a)
            m.csv_gz(self.tables["primary"], b)
            self.assertEqual(a.read_bytes(), b.read_bytes())
            roundtrip = pd.read_csv(a, compression="gzip", dtype=str, keep_default_na=False)
            pd.testing.assert_frame_equal(roundtrip, self.tables["primary"])
            with gzip.open(a, "rt") as f:
                self.assertIn("pydeseq2_unmasked_pvalue", f.readline())
            with self.assertRaises(FileExistsError):
                m.csv_gz(self.tables["primary"], a)

    def test_verified_manifest_and_file_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            for key, filename in m.TABLES.items():
                m.csv_gz(self.tables[key], base / filename)
            (base / "organoid-qc.json").write_text(json.dumps(self.qc))
            manifest = {"plan_sha256": "test", "source_values_mutated": False,
                        "outputs": {p.name: {"sha256": m.digest(p), "bytes": p.stat().st_size}
                                    for p in base.iterdir()}}
            (base / "organoid-output-manifest.json").write_text(json.dumps(manifest))
            receipt = copy.deepcopy(self.verification)
            receipt["public_manifest_sha256"] = m.digest(base / "organoid-output-manifest.json")
            verify = base / "verification.json"
            verify.write_text(json.dumps(receipt))
            tables, qc, _ = m.read_verified_inputs(base, verify)
            self.assertEqual(len(tables["primary"]), 4)
            self.assertEqual(qc["eligible_features"], 3)
            with (base / m.TABLES["primary"]).open("ab") as f:
                f.write(b"bad")
            with self.assertRaisesRegex(ValueError, "Changed input"):
                m.read_verified_inputs(base, verify)

    def test_summary_is_finite_json_and_source_unchanged(self):
        before = {key: frame.copy(deep=True) for key, frame in self.tables.items()}
        summary, _ = m.summarize(self.tables, self.qc, self.verification)
        json.dumps(summary, allow_nan=False)
        for key in before:
            pd.testing.assert_frame_equal(before[key], self.tables[key])


if __name__ == "__main__":
    unittest.main()
