"""Generated synthetic math/guard fixtures only. No real expression values loaded."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("psc_liver_program_method",ROOT/"scripts/analyze_psc_liver_programs.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
CORE, STRICT = MODULE.INTERFACES


def component(name="all",source=12,indices=None):
    return {"name":name,"source_members":source,"indices":list(range(10)) if indices is None else list(indices)}


def program(program_id="SYN_A",interface=CORE,parts=None,kind="unsigned_equal_gene",parent_gate=None):
    parts = [component()] if parts is None else parts
    return {"program_id":program_id,"interface":interface,"components":parts,"score_kind":kind,
            "coverage":MODULE.coverage_status(parts,kind,parent_gate)}


def synthetic_fixture(replicates=67):
    units = [{"unit_id":f"SYNTHETIC_SAMPLE_{i:02d}","group":group} for i,group in enumerate(
        ["PSC"]*3+["ASC"]*3+["AIH"]*4)]
    rng = np.random.Generator(np.random.PCG64(7113))
    counts = rng.integers(1,1000,size=(24,len(units))).astype(np.float64)
    counts[0,:] = 0  # A measured zero gene is retained, unlike a missing member.
    specs = []
    for interface in MODULE.INTERFACES:
        specs.extend([
            program("SYN_A",interface,parts=[component(indices=range(10 if interface==CORE else 9))]),
            program("SYN_B",interface,parts=[component(source=12,indices=range(10,20))]),
            program("SYN_SIGNED",interface,parts=[component("positive",5,range(4)),component("negative",10,range(10,18))],kind="balanced_signed"),
            program("SYN_UNAVAILABLE",interface,parts=[component(source=20,indices=[])])])
    method = MODULE.draft_plan_method_without_runtime()
    method["zero_variance_policy"] = "welch_if_one_positive"
    method["bootstrap"].update(seed=90210,replicates=replicates)
    return counts,specs,units,method


class CoverageTests(unittest.TestCase):
    def test_exact_eighty_percent_and_boundary_below(self):
        self.assertEqual(MODULE.coverage_status([component(source=15,indices=range(12))],"unsigned_equal_gene")["score_status"],"available")
        self.assertEqual(MODULE.coverage_status([component(source=15,indices=range(11))],"unsigned_equal_gene")["score_status"],"unavailable")

    def test_ten_member_boundary(self):
        for n,expected in ((9,"unavailable"),(10,"available")):
            self.assertEqual(MODULE.coverage_status([component(source=n,indices=range(n))],"unsigned_equal_gene")["score_status"],expected)

    def test_duplicate_and_negative_indices_fail(self):
        for indices in ([0]*10,range(-1,9),[True]+list(range(1,10))):
            with self.assertRaises(MODULE.GuardError):
                MODULE.coverage_status([component(indices=indices)],"unsigned_equal_gene")

    def test_denominator_cannot_shrink_or_zero(self):
        for denominator in (0,True,9):
            with self.assertRaises(MODULE.GuardError):
                MODULE.coverage_status([component(source=denominator)],"unsigned_equal_gene")

    def test_parent_gate_not_new_nonoverlap_fraction(self):
        gate = MODULE.coverage_status([component(source=100,indices=range(10))],"unsigned_equal_gene",True)
        self.assertEqual(gate["score_status"],"available")
        self.assertEqual(gate["gate_basis"],"parent_gate_plus_10_retained")
        self.assertEqual(MODULE.coverage_status([component(source=100)],"unsigned_equal_gene",False)["score_status"],"unavailable")

    def test_nonoverlap_still_has_ten_retained_gate(self):
        self.assertEqual(MODULE.coverage_status([component(source=100,indices=range(9))],"unsigned_equal_gene",True)["score_status"],"unavailable")

    def test_signed_each_component_not_only_union(self):
        parts = [component("positive",10,range(9)),component("negative",10,range(20,27))]
        self.assertEqual(MODULE.coverage_status(parts,"balanced_signed")["score_status"],"unavailable")

    def test_signed_combined_ten_not_ten_per_component(self):
        parts = [component("positive",5,range(5)),component("negative",5,range(5,10))]
        self.assertEqual(MODULE.coverage_status(parts,"balanced_signed")["score_status"],"available")

    def test_missing_signed_component_never_zero_filled(self):
        parts = [component("positive",10,range(10)),component("negative",5,[])]
        self.assertIn("empty_required_component:negative",MODULE.coverage_status(parts,"balanced_signed")["reasons"])
        with self.assertRaises(MODULE.GuardError):
            MODULE.coverage_status(parts[:1],"balanced_signed")

    def test_signed_cross_component_overlap_rejected(self):
        with self.assertRaises(MODULE.GuardError):
            MODULE.coverage_status([component("positive",10,range(10)),component("negative",10,range(10))],"balanced_signed")


class NormalizationScoringTests(unittest.TestCase):
    def test_hand_computed_full_universe_cpm(self):
        counts = np.array([[1,2],[3,6],[6,2]])
        actual,totals = MODULE.normalize_counts(counts,(3,2))
        np.testing.assert_array_equal(totals,[10,10])
        np.testing.assert_allclose(actual,[[math.log2(100001),math.log2(200001)],
            [math.log2(300001),math.log2(600001)],[math.log2(600001),math.log2(200001)]],rtol=1e-14)

    def test_nonprogram_rows_change_program_score(self):
        counts = np.ones((12,3))
        before,_ = MODULE.normalize_counts(counts,(12,3))
        counts[11,:] = 1000
        after,_ = MODULE.normalize_counts(counts,(12,3))
        self.assertTrue((MODULE.score_program(after,program()) < MODULE.score_program(before,program())).all())

    def test_observed_zero_gene_and_zero_program_retained(self):
        counts = np.zeros((12,3)); counts[11,:]=10
        logged,_ = MODULE.normalize_counts(counts,(12,3))
        score = MODULE.score_program(logged,program())
        np.testing.assert_array_equal(score,np.zeros(3))
        self.assertEqual(program()["coverage"]["mapped_unique_total"],10)

    def test_zero_library_fails_not_dropped(self):
        with self.assertRaises(MODULE.DataIntegrityError):
            MODULE.normalize_counts(np.zeros((2,3)),(2,3))

    def test_negative_nonfinite_noninteger_missing_counts_fail(self):
        for value in (-1,np.nan,np.inf,0.5,None,2**54):
            counts = [[1,value],[1,1]]
            with self.subTest(value=value),self.assertRaises(MODULE.DataIntegrityError):
                MODULE.normalize_counts(counts,(2,2))

    def test_exact_integer_total_rejects_lost_overflow_unit(self):
        with self.assertRaises(MODULE.DataIntegrityError):
            MODULE.normalize_counts(np.array([[2**53],[1]],dtype=np.float64),(2,1))
        with self.assertRaises(MODULE.DataIntegrityError):
            MODULE.normalize_counts(np.full((2048,1),2**53,dtype=np.uint64),(2048,1))
        _,totals=MODULE.normalize_counts(np.array([[2**53],[0]],dtype=np.uint64),(2,1))
        self.assertEqual(totals[0],2**53)

    def test_integer_limit_checked_before_lossy_float_conversion(self):
        for dtype in (np.int64,np.uint64):
            with self.assertRaises(MODULE.DataIntegrityError):
                MODULE.normalize_counts(np.array([[2**53+1],[0]],dtype=dtype),(2,1))

    def test_boolean_complex_and_string_recodings_rejected(self):
        for counts in (np.ones((2,2),dtype=bool),np.ones((2,2),dtype=complex),np.full((2,2),"1")):
            with self.assertRaises(MODULE.DataIntegrityError):MODULE.normalize_counts(counts,(2,2))

    def test_shape_guard(self):
        with self.assertRaises(MODULE.DataIntegrityError):
            MODULE.normalize_counts(np.ones((2,2)),(3,2))

    def test_sample_permutation_keeps_corresponding_values(self):
        counts = np.arange(1,37).reshape(12,3)
        a,_ = MODULE.normalize_counts(counts,(12,3))
        b,_ = MODULE.normalize_counts(counts[:,[2,0,1]],(12,3))
        np.testing.assert_array_equal(a[:,[2,0,1]],b)

    def test_mean_log_not_log_of_mean(self):
        logged = np.arange(36,dtype=float).reshape(12,3)
        np.testing.assert_array_equal(MODULE.score_program(logged,program()),np.mean(logged[:10],axis=0))

    def test_balanced_signed_not_global_equal_gene(self):
        logged = np.zeros((20,3)); logged[:4,:]=10; logged[10:18,:]=2
        spec = program(parts=[component("positive",5,range(4)),component("negative",10,range(10,18))],kind="balanced_signed")
        np.testing.assert_array_equal(MODULE.score_program(logged,spec),[4,4,4])
        self.assertNotEqual(4,(4*10-8*2)/12)

    def test_coverage_failure_returns_none_not_zero(self):
        self.assertIsNone(MODULE.score_program(np.zeros((12,3)),program(parts=[component(indices=range(9))])))

    def test_nonfinite_score_fails_not_nanmean(self):
        logged = np.ones((12,3)); logged[0,0]=np.nan
        with self.assertRaises(MODULE.DataIntegrityError):
            MODULE.score_program(logged,program())


class WelchTests(unittest.TestCase):
    def test_welch_matches_scipy_and_fractional_df(self):
        a,b = [1,2,8,11],[0,1,2,3,4,6]
        result = MODULE.welch_summary(a,b,"welch_if_one_positive")
        oracle = stats.ttest_ind(a,b,equal_var=False)
        self.assertAlmostEqual(result["p_raw"],float(oracle.pvalue),14)
        self.assertAlmostEqual(result["degrees_of_freedom"],float(oracle.df),14)
        self.assertNotEqual(result["degrees_of_freedom"],len(a)+len(b)-2)
        a_term,b_term = np.var(a,ddof=1)/len(a),np.var(b,ddof=1)/len(b)
        expected_df=(a_term+b_term)**2/(a_term*a_term/(len(a)-1)+b_term*b_term/(len(b)-1))
        expected_width=stats.t.ppf(.975,expected_df)*math.sqrt(a_term+b_term)
        self.assertAlmostEqual(result["ci_high"]-result["effect"],expected_width,13)

    def test_one_zero_variance_mathematical_welch(self):
        result=MODULE.welch_summary([0,0,0],[1,2,3],"welch_if_one_positive")
        self.assertEqual(result["degrees_of_freedom"],2)
        self.assertAlmostEqual(result["p_raw"],1-math.sqrt(6/7),14)
        self.assertEqual(result["zero_variance_groups"],["a"])
        self.assertEqual(result["welch_status"],"available")

    def test_one_zero_policy_conservative_is_explicit(self):
        result=MODULE.welch_summary([0,0,0],[1,2,3],"unavailable_if_either_zero")
        self.assertEqual(result["effect"],-2)
        self.assertIsNone(result["p_raw"])
        self.assertEqual(result["welch_status"],"unavailable_zero_variance_policy")

    def test_both_zero_equal_and_different_descriptive_only(self):
        for b in ([1,1,1],[7,7,7]):
            result=MODULE.welch_summary([1,1,1],b,"welch_if_one_positive")
            self.assertEqual(result["effect"],1-b[0])
            for key in ("p_raw","degrees_of_freedom","standard_error","ci_low","ci_high"):
                self.assertIsNone(result[key])

    def test_missing_and_single_observation_status(self):
        zero=MODULE.welch_summary([], [1,2],"welch_if_one_positive")
        one=MODULE.welch_summary([4],[1,2],"welch_if_one_positive")
        self.assertIsNone(zero["effect"])
        self.assertEqual(one["effect"],2.5)
        self.assertIsNone(one["p_raw"])
        self.assertEqual(one["welch_status"],"unavailable_n_below_2")

    def test_nonfinite_is_data_failure_not_missing_p(self):
        for bad in (np.nan,np.inf):
            with self.assertRaises(MODULE.DataIntegrityError):
                MODULE.welch_summary([1,bad],[2,3],"welch_if_one_positive")

    def test_tiny_positive_variance_is_not_thresholded(self):
        result=MODULE.welch_summary(np.array([1,2,4])*1e-150,np.array([2,4,9,10])*1e-150,"welch_if_one_positive")
        self.assertEqual(result["welch_status"],"available")
        self.assertGreater(result["variance_a"],0)
        self.assertAlmostEqual(result["p_raw"],float(stats.ttest_ind([1,2,4],[2,4,9,10],equal_var=False).pvalue),13)

    def test_contrast_reversal_changes_only_direction(self):
        a,b=[1,3,8,12],[0,1,2]
        x=MODULE.welch_summary(a,b,"welch_if_one_positive")
        y=MODULE.welch_summary(b,a,"welch_if_one_positive")
        self.assertAlmostEqual(x["p_raw"],y["p_raw"],15)
        self.assertAlmostEqual(x["degrees_of_freedom"],y["degrees_of_freedom"],15)
        self.assertEqual(x["effect"],-y["effect"])
        self.assertAlmostEqual(x["ci_low"],-y["ci_high"],14)

    def test_exact_float_constants_do_not_invent_variance(self):
        result=MODULE.welch_summary(np.full(30,2.1),np.full(17,4.2),"welch_if_one_positive")
        self.assertEqual(result["variance_a"],0)
        self.assertEqual(result["variance_b"],0)
        self.assertEqual(result["mean_a"],2.1)
        self.assertIsNone(result["p_raw"])

    def test_nonconstant_variance_underflow_retains_effect_not_false_zero(self):
        result=MODULE.welch_summary(np.array([1,2,4])*1e-200,[1,2,4],"welch_if_one_positive")
        self.assertEqual(result["effect_status"],"available")
        self.assertTrue(math.isfinite(result["effect"]))
        self.assertEqual(result["welch_status"],"unavailable_numeric_variance_underflow")
        self.assertNotIn("a",result["zero_variance_groups"])
        for key in ("variance_a","p_raw","degrees_of_freedom","ci_low","ci_high"):
            self.assertIsNone(result[key])

    def test_variance_overflow_is_unavailable_inference_not_invalid_input(self):
        result=MODULE.welch_summary([-1e200,1e200],[-2e200,2e200],"welch_if_one_positive")
        self.assertEqual(result["effect"],0)
        self.assertEqual(result["welch_status"],"unavailable_numeric_variance")
        self.assertIsNone(result["variance_a"])
        for key in ("p_raw","degrees_of_freedom","ci_low","ci_high"):
            self.assertIsNone(result[key])
        json.dumps(result,allow_nan=False)

    def test_tail_underflow_never_becomes_observed_zero(self):
        result=MODULE.welch_summary(np.arange(30.)*1e-14,np.full(17,2.1),"welch_if_one_positive")
        self.assertEqual(result["welch_status"],"unavailable_numeric_tail_underflow")
        self.assertTrue(result["p_underflow_to_zero"])
        self.assertTrue(math.isfinite(result["effect"]))
        for key in ("p_raw","degrees_of_freedom","ci_low","ci_high"):
            self.assertIsNone(result[key])

    def test_nonfinite_test_arithmetic_keeps_finite_descriptive_effect(self):
        for function,value in (("sf",np.nan),("ppf",np.inf)):
            with patch.object(MODULE.stats.t,function,return_value=value):
                result=MODULE.welch_summary([0,1,3],[2,4,5],"welch_if_one_positive")
            self.assertAlmostEqual(result["effect"],-7/3)
            self.assertEqual(result["welch_status"],"unavailable_numeric_inference")
            for key in ("p_raw","degrees_of_freedom","ci_low","ci_high"):
                self.assertIsNone(result[key])

    def test_unknown_policy_rejected(self):
        with self.assertRaises(MODULE.GuardError):
            MODULE.welch_summary([1,2],[3,4],"guess")


class CorrectionTests(unittest.TestCase):
    def test_known_bh_and_by_with_missing(self):
        values=[.01,.04,.03,None]
        np.testing.assert_allclose(MODULE.adjust_p_values(values,"BH"),[.04,4*.04/3,4*.04/3,1])
        expected=np.minimum(1,np.array([.04,4*.04/3,4*.04/3,1])*sum(1/i for i in range(1,5)))
        np.testing.assert_allclose(MODULE.adjust_p_values(values,"BY"),expected)

    def test_missing_not_dropped_from_denominator(self):
        self.assertNotEqual(MODULE.adjust_p_values([.01,None],"BH")[0],MODULE.adjust_p_values([.01],"BH")[0])
        self.assertEqual(MODULE.adjust_p_values([None,None],"BH"),[1,1])

    def test_ties_zero_one_preserve_order(self):
        values=[1,0,.2,.2,.05]
        adjusted=MODULE.adjust_p_values(values,"BH")
        self.assertEqual(adjusted[0],1); self.assertEqual(adjusted[1],0)
        self.assertEqual(adjusted[2],adjusted[3])
        self.assertTrue(all(0<=v<=1 for v in adjusted))

    def test_invalid_p_not_coerced_to_missing(self):
        for value in (np.nan,np.inf,-.01,1.01,True):
            with self.subTest(value=value),self.assertRaises(MODULE.DataIntegrityError):
                MODULE.adjust_p_values([value],"BH")

    def test_full_family_raw_null_and_strict_no_q(self):
        rows=[{"program_id":"P","contrast_id":c[0],"interface":interface,"p_raw":p}
              for interface in MODULE.INTERFACES for c,p in zip(MODULE.CONTRASTS,[.01,None,.05])]
        receipt=MODULE.correction_family(rows,["P"])
        self.assertEqual(receipt["family_size"],3)
        self.assertIsNone(rows[1]["p_raw"]); self.assertIsNone(rows[1]["q_bh"])
        self.assertEqual(rows[1]["p_for_correction"],1)
        self.assertTrue(rows[1]["p_was_substituted"])
        for row in rows[3:]:
            self.assertIsNone(row["q_bh"]); self.assertIsNone(row["p_for_correction"])
            self.assertFalse(row["discovery_eligible"])

    def test_incomplete_or_duplicated_family_rejected(self):
        rows=[{"program_id":"P","contrast_id":MODULE.CONTRASTS[0][0],"interface":CORE,"p_raw":.01}]
        for bad in (rows,rows*3):
            with self.assertRaises(MODULE.GuardError):
                MODULE.correction_family(bad,["P"])


class ResamplingTests(unittest.TestCase):
    def setUp(self):
        self.units=[{"unit_id":"SYN"+str(i),"group":group} for i,group in enumerate(["PSC"]*3+["ASC"]*2+["AIH"]*4)]

    def test_joint_bootstrap_reproducible_and_group_stratified(self):
        a=MODULE.stratified_draws(self.units,100,991)
        b=MODULE.stratified_draws(self.units,100,991)
        for group in MODULE.GROUPS:
            np.testing.assert_array_equal(a[group],b[group])
            self.assertEqual(a[group].shape,(100,sum(u["group"]==group for u in self.units)))
            self.assertTrue(all(self.units[i]["group"]==group for i in a[group].ravel()))

    def test_joint_contrasts_preserve_linear_relation(self):
        draws=MODULE.stratified_draws(self.units,51,23)
        d=MODULE.bootstrap_differences(np.arange(9.),draws)
        np.testing.assert_allclose(d["PSC_minus_AIH"],d["PSC_minus_ASC"]+d["ASC_minus_AIH"],rtol=0,atol=2e-15)

    def test_constant_draws_retained_no_bootstrap_p(self):
        d=MODULE.bootstrap_differences(np.ones(9),MODULE.stratified_draws(self.units,99,3))
        result=MODULE.bootstrap_summary(d["PSC_minus_AIH"],99)
        self.assertEqual(result["valid_replicates"],99)
        self.assertEqual(result["status"],"available_degenerate")
        self.assertEqual(result["ci_low"],0)
        self.assertFalse(result["p_value_computed"])

    def test_constant_decimal_bootstrap_has_no_rounding_effect(self):
        units=[{"unit_id":"SYN"+str(i),"group":group} for i,group in enumerate(["PSC"]*17+["ASC"]*17+["AIH"]*30)]
        values=MODULE.bootstrap_differences(np.full(64,2.1),MODULE.stratified_draws(units,19,100))
        for vector in values.values():np.testing.assert_array_equal(vector,np.zeros(19))

    def test_bootstrap_no_bad_draw_deletion(self):
        with self.assertRaises(MODULE.DataIntegrityError):
            MODULE.bootstrap_summary([1,np.nan,2],3)
        with self.assertRaises(MODULE.DataIntegrityError):
            MODULE.bootstrap_summary([1,2],3)
        self.assertEqual(MODULE.bootstrap_summary(None,100)["status"],"unavailable_score")

    def test_every_omission_including_uninvolved_units(self):
        rows=MODULE.omission_summaries(np.arange(9.),self.units,MODULE.CONTRASTS[0],"welch_if_one_positive")
        self.assertEqual(len(rows),9)
        outside=[r for r in rows if not r["in_contrast"]]
        self.assertEqual(len(outside),2)
        self.assertTrue(all(r["status"]=="unchanged_outside_contrast" and r["effect_delta"]==0 for r in outside))
        self.assertEqual(rows[0]["statistics"]["n_a"],2)
        self.assertEqual(rows[0]["statistics"]["n_b"],4)

    def test_public_omission_denominators_exclude_every_outside_unit(self):
        units=[{"unit_id":"SYN_OMIT"+str(i),"group":g} for i,g in enumerate(["PSC"]*17+["ASC"]*17+["AIH"]*30)]
        for contrast,expected in zip(MODULE.CONTRASTS,(47,34,47)):
            rows=MODULE.omission_summaries(np.arange(64.),units,contrast,"welch_if_one_positive")
            public=MODULE.summarize_omissions(rows)
            self.assertEqual(len(rows),64)
            self.assertEqual(public["relevant_unit_denominator"],expected)
            self.assertEqual(public["defined_effect_rows"],expected)
            self.assertNotIn("all_unit_rows",public)
            self.assertNotIn("unchanged_outside_contrast_rows",public)
            expected_only=MODULE.summarize_omissions([row for row in rows if row["in_contrast"]])
            self.assertEqual(public,expected_only)
            for row in rows:
                if not row["in_contrast"]:
                    row["statistics"]["effect"]=1e100
                    row["effect_delta"]=1e100
                    row["strict_sign_reversal"]=True
            self.assertEqual(MODULE.summarize_omissions(rows),expected_only)

    def test_unavailable_omission_effects_do_not_imply_zero_reversals(self):
        rows=MODULE.omission_summaries(None,self.units,MODULE.CONTRASTS[0],"welch_if_one_positive")
        public=MODULE.summarize_omissions(rows)
        self.assertEqual(public["relevant_unit_denominator"],7)
        self.assertEqual(public["defined_effect_rows"],0)
        self.assertIsNone(public["strict_sign_reversal_rows"])

    def test_root_bootstrap_seed_order_and_size_on_generated_units(self):
        method=MODULE.draft_plan_method_without_runtime()["bootstrap"]
        self.assertEqual(method["replicates"],10000)
        self.assertEqual(method["seed"],2026091602)
        self.assertEqual(method["group_order"],["PSC","ASC","AIH"])
        draws=MODULE.stratified_draws(self.units,method["replicates"],method["seed"])
        rng=np.random.Generator(np.random.PCG64(2026091602))
        for group in ("PSC","ASC","AIH"):
            indices=[i for i,u in enumerate(self.units) if u["group"]==group]
            expected=rng.choice(indices,size=(10000,len(indices)),replace=True)
            np.testing.assert_array_equal(draws[group],expected)

    def test_omission_coverage_unavailable_not_reselected(self):
        rows=MODULE.omission_summaries(None,self.units,MODULE.CONTRASTS[0],"welch_if_one_positive")
        self.assertEqual(len(rows),len(self.units))
        self.assertTrue(all(r["statistics"] is None for r in rows))

    def test_omission_small_n_preserves_defined_effect(self):
        units=[{"unit_id":"SYN"+str(i),"group":g} for i,g in enumerate(["PSC","PSC","ASC","AIH","AIH"])]
        row=MODULE.omission_summaries(np.arange(5.),units,MODULE.CONTRASTS[0],"welch_if_one_positive")[0]
        self.assertIsNotNone(row["statistics"]["effect"])
        self.assertIsNone(row["statistics"]["p_raw"])


class SyntheticPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.counts,cls.specs,cls.units,cls.method=synthetic_fixture()
        cls.result=MODULE.analyze_counts(cls.counts,cls.specs,cls.units,cls.method,cls.counts.shape)

    def test_complete_status_inventory_and_family(self):
        public=self.result["public"]
        self.assertEqual(len(public["endpoint_rows"]),4*3*2)
        self.assertEqual(public["core_family_size"],4*3)
        self.assertEqual(len(self.result["omission_rows"]),10*4*3*2)
        self.assertEqual(public["unavailable_core_p_values"],3)

    def test_core_available_strict_unavailable_same_original_denominator(self):
        rows=[r for r in self.result["public"]["endpoint_rows"] if r["program_id"]=="SYN_A"]
        self.assertTrue(all(r["coverage"]["components"][0]["original_source_members"]==12 for r in rows))
        self.assertTrue(all(r["effect"] is not None for r in rows if r["interface"]==CORE))
        self.assertTrue(all(r["effect"] is None and r["q_bh"] is None for r in rows if r["interface"]==STRICT))

    def test_public_has_no_unit_ids_or_per_unit_values(self):
        text=json.dumps(self.result["public"],allow_nan=False)
        self.assertNotIn("SYNTHETIC_SAMPLE_",text)
        self.assertNotIn('"unit_id"',text)
        self.assertNotIn('"values"',text)
        self.assertNotIn('"indices"',text)

    def test_every_available_program_three_contrast_identity(self):
        for interface in MODULE.INTERFACES:
            for program_id in ("SYN_A","SYN_B","SYN_SIGNED"):
                rows=[r for r in self.result["public"]["endpoint_rows"] if r["interface"]==interface and r["program_id"]==program_id]
                if rows[0]["effect"] is not None:
                    self.assertAlmostEqual(rows[0]["effect"],rows[1]["effect"]+rows[2]["effect"],13)

    def test_no_score_recomputation_during_omissions(self):
        counts,specs,units,method=synthetic_fixture(5)
        with patch.object(MODULE,"score_program",wraps=MODULE.score_program) as scorer:
            MODULE.analyze_counts(counts,specs,units,method,counts.shape)
        self.assertEqual(scorer.call_count,len(specs))

    def test_corrupt_scores_fail_whole_computation(self):
        counts=self.counts.copy();counts[0,0]=np.nan
        with self.assertRaises(MODULE.DataIntegrityError):
            MODULE.analyze_counts(counts,self.specs,self.units,self.method,counts.shape)

    def test_inference_unavailable_keeps_complete_core_correction_family(self):
        counts,specs,units,method=synthetic_fixture(5)
        with patch.object(MODULE.stats.t,"sf",return_value=np.nan):
            result=MODULE.analyze_counts(counts,specs,units,method,counts.shape)
        rows=result["public"]["endpoint_rows"]
        core=[row for row in rows if row["interface"]==CORE]
        self.assertEqual(len(rows),24)
        self.assertEqual(result["correction"]["family_size"],12)
        self.assertTrue(all(row["p_raw"] is None and row["q_bh"] is None and row["q_by"] is None for row in rows))
        self.assertTrue(all(row["p_for_correction"]==1 for row in core))
        self.assertTrue(any(row["effect"] is not None for row in core))
        self.assertTrue(all(row["degrees_of_freedom"] is None and row["ci_low"] is None for row in rows))

    def test_full_universe_never_replaced_by_program_rows(self):
        counts=self.counts[:20,:]
        with self.assertRaises(MODULE.DataIntegrityError):
            MODULE.analyze_counts(counts,self.specs,self.units,self.method,self.counts.shape)


class MetadataMappingTests(unittest.TestCase):
    def fixture(self):
        rows=[{"program_id":"SYN_EXT","source_member_index":str(i+1),"source_value":f"SYN_G{i}",
            "source_identifier_namespace":"MSigDB_geneSymbols","published_direction":"UNSIGNED",
            "membership_weight":"1"} for i in range(12)]
        features=[{"label_as_submitted":f"SYN_G{i}","source_row_index_0based":i,
            CORE:"eligible",STRICT:"eligible" if i<9 else "no_entrez_relation"} for i in range(10)]
        endpoint={"program_id":"SYN_EXT","score_kind":"unsigned_equal_gene",
            "original_component_members":{"all":12},"direction_literals":{"all":"UNSIGNED"},
            "source_identifier_namespace":"MSigDB_geneSymbols","target_gene_exclusions":[]}
        return rows,features,endpoint

    def test_exact_mapping_no_alias_or_missing_index_zero(self):
        rows,features,endpoint=self.fixture()
        specs,audit=MODULE.map_external_rows(rows,endpoint,features)
        self.assertEqual(specs[(CORE,"SYN_EXT")]["coverage"]["score_status"],"available")
        self.assertEqual(specs[(STRICT,"SYN_EXT")]["coverage"]["score_status"],"unavailable")
        self.assertEqual(audit[0]["feature_index"],0)
        self.assertIsNone(audit[10]["feature_index"])

    def test_repeated_source_ids_cannot_raise_weights(self):
        rows,features,endpoint=self.fixture()
        rows[0]["membership_weight"]="2"
        with self.assertRaises(MODULE.GuardError):
            MODULE.map_external_rows(rows,endpoint,features)

    def test_target_exclusion_not_automatically_extended(self):
        rows,features,endpoint=self.fixture();endpoint["target_gene_exclusions"]=["ETS2"]
        with self.assertRaises(MODULE.GuardError):
            MODULE.map_external_rows(rows,endpoint,features)

    def test_duplicate_symbols_not_collapsed(self):
        rows,features,endpoint=self.fixture();rows[1]["source_value"]=rows[0]["source_value"]
        with self.assertRaises(MODULE.GuardError):
            MODULE.map_external_rows(rows,endpoint,features)

    def test_unknown_source_polarity_not_guessed(self):
        rows,features,endpoint=self.fixture();rows[0]["published_direction"]="UP"
        with self.assertRaises(MODULE.GuardError):
            MODULE.map_external_rows(rows,endpoint,features)

    def test_source_case_not_folded_or_alias_mapped(self):
        rows,features,endpoint=self.fixture();rows[0]["source_value"]="syn_g0"
        specs,audit=MODULE.map_external_rows(rows,endpoint,features)
        self.assertIsNone(audit[0]["feature_index"])
        self.assertEqual(specs[(CORE,"SYN_EXT")]["coverage"]["score_status"],"unavailable")


class GuardTests(unittest.TestCase):
    def setUp(self):
        base=ROOT/MODULE.WORK_REL/"synthetic-tests"
        base.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=base)
        self.root=Path(self.temp.name)
        for rel in (MODULE.SCRIPT_REL,MODULE.TEST_REL):
            path=self.root/rel;path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text("# Generated synthetic implementation-pin fixture\n")
        self.plan=MODULE.draft_plan(self.root)
        self.plan.update(state="frozen",panel_approved=True,expected_compiled_programs_sha256=MODULE.canonical_sha([]))
        self.plan["root_authorization"]={"approver_role":"root","reference":"SYNTHETIC_AUTH_NOT_REAL",
            "approved_at":"SYNTHETIC_TIMESTAMP","real_execution_authorized":True}
        self.plan["decisions"]={key:{"status":"approved","reference":"SYNTHETIC_DECISION_NOT_REAL"} for key in MODULE.DECISION_NAMES}
        self.plan["unit"]["assumption_approved"]=True
        self.plan["method"]["zero_variance_policy"]="welch_if_one_positive"
        self.plan["method"]["bootstrap"]["seed"]=MODULE.BOOTSTRAP_SEED
        review=self.root/"work/synthetic-independent-review.json"
        review.parent.mkdir(parents=True,exist_ok=True)
        MODULE.dump_json(review,{"status":"pass_synthetic_and_guard_review",
            "implementation_sha256":{p["path"]:p["sha256"] for p in self.plan["implementation_pins"]},
            "real_execution_reviewed":False})
        self.review_path=review
        self.plan["independent_review"]={"status":"accepted_synthetic_and_guard_review","reference":"SYNTHETIC_REVIEW",
            "artifact_pin":{"path":str(review.relative_to(self.root)),"sha256":MODULE.file_sha(review)}}
        self.plan_path=self.root/"config/psc-liver-program-plan.json"
        self.plan_path.parent.mkdir(parents=True,exist_ok=True)
        self.output=self.root/MODULE.RUNS_REL/"synthetic-run"
        self.write_plan()

    def tearDown(self):
        self.temp.cleanup()

    def write_plan(self):
        MODULE.dump_json(self.plan_path,self.plan)
        self.sha=MODULE.file_sha(self.plan_path)

    def verify_synthetic_pin(self,root,pin):
        if pin["path"] in (MODULE.RAW_MATRIX_PIN["path"],MODULE.MATRIX_PIN["path"]):
            return self.root/"SYNTHETIC_MATRIX_NEVER_OPENED.npy"
        return MODULE.verify_pin(root,pin)

    def invoke_guards(self,plan=None,expected_sha=None,ack=None):
        if plan is not None:
            self.plan=plan;self.write_plan()
        with patch.object(MODULE,"load_count_matrix",side_effect=AssertionError("REAL NUMERIC LOADER MUST NOT RUN")) as loader:
            result=MODULE.authorize_execution(self.plan_path,expected_sha or self.sha,
                MODULE.EXECUTION_ACK if ack is None else ack,self.output,self.root)
            self.assertEqual(loader.call_count,0)
            return result

    def test_unapproved_plan_refused_before_metadata_or_values(self):
        self.plan["state"]="draft_not_frozen";self.write_plan()
        with patch.object(MODULE,"source_only_preflight",side_effect=AssertionError("must not reach metadata")):
            with self.assertRaises(MODULE.GuardError):self.invoke_guards()
        self.assertFalse(self.output.exists())

    def test_missing_explicit_execution_action(self):
        with self.assertRaises(MODULE.GuardError):self.invoke_guards(ack="")
        self.assertFalse(self.output.exists())

    def test_changed_frozen_plan_sha_refused(self):
        with self.assertRaises(MODULE.GuardError):self.invoke_guards(expected_sha="0"*64)

    def test_missing_root_authorization_refused(self):
        self.plan["root_authorization"]["real_execution_authorized"]=False;self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_unresolved_scientific_choice_refused(self):
        self.plan["decisions"]["zero_variance"]["status"]="proposed_unresolved";self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_unresolved_seed_refused(self):
        self.plan["method"]["bootstrap"]["seed"]=None;self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_other_seeds_and_zero_policies_refused_in_real_plan(self):
        for changed in ("seed","zero","group_order","omissions"):
            plan=deepcopy(self.plan)
            if changed=="seed":plan["method"]["bootstrap"]["seed"]=303271
            elif changed=="zero":plan["method"]["zero_variance_policy"]="unavailable_if_either_zero"
            elif changed=="group_order":plan["method"]["bootstrap"]["group_order"]=["AIH","ASC","PSC"]
            else:plan["method"]["omissions"]["public_stability_units"]="all_units"
            with self.assertRaises(MODULE.GuardError):MODULE.validate_method(plan)

    def test_root_method_reference_required(self):
        self.plan["root_method_decision_reference"]="UNREVIEWED_ALTERNATIVE";self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_hold_source_refused_not_zero_endpoint(self):
        self.plan["panel"].append({"program_id":"SYN_HELD","status":"HOLD","kind":"external_symbol_csv",
            "source_approval_reference":None});self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_missing_fixed_endpoint_refused(self):
        self.plan["panel"].pop(0);self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_source_denominator_cannot_be_amended_silently(self):
        self.plan["panel"][0]["original_source_members"]=800;self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_known_repeats_or_unapproved_independence_blocks(self):
        for field,value in (("assumption_approved",False),("known_repeat_status","known_repeats_unresolved")):
            modified=deepcopy(self.plan);modified["unit"][field]=value
            with self.assertRaises(MODULE.GuardError):self.invoke_guards(plan=modified)

    def test_runtime_or_implementation_pin_change_blocks(self):
        self.plan["software"]["numpy"]="SYNTHETIC_OTHER_VERSION";self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()
        self.plan["software"]=MODULE.software_versions()
        self.plan["implementation_pins"][0]["sha256"]="0"*64;self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_stale_independent_review_blocks(self):
        record=MODULE.read_json(self.review_path);record["implementation_sha256"]={}
        MODULE.dump_json(self.review_path,record)
        self.plan["independent_review"]["artifact_pin"]["sha256"]=MODULE.file_sha(self.review_path);self.write_plan()
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_reused_output_even_empty_is_refused(self):
        self.output.mkdir(parents=True)
        with self.assertRaises(MODULE.GuardError):self.invoke_guards()

    def test_protected_or_symlink_output_refused(self):
        for output in (self.root/"data/derived",self.root/"reports/new",self.root/MODULE.RUNS_REL):
            with self.assertRaises(MODULE.GuardError):MODULE.fresh_output_path(output,self.root)
        allowed=self.root/MODULE.RUNS_REL;allowed.mkdir(parents=True)
        target=self.root/"outside";target.mkdir()
        (allowed/"link").symlink_to(target,target_is_directory=True)
        with self.assertRaises(MODULE.GuardError):MODULE.fresh_output_path(allowed/"link/new",self.root)

    def test_inputs_cannot_escape_repository(self):
        with self.assertRaises(MODULE.GuardError):MODULE.relative_input(self.root,"../outside")
        with self.assertRaises(MODULE.GuardError):MODULE.relative_input(self.root,"/outside")

    def test_registered_metadata_symlink_aliases_refused_before_hash(self):
        payload=self.root/MODULE.MATRIX_PIN["path"]
        payload.parent.mkdir(parents=True,exist_ok=True)
        payload.write_bytes(b"SYNTHETIC_NOT_REAL_EXPRESSION")
        alias=self.root/"config/SYN_METADATA.json";alias.symlink_to(payload)
        with patch.object(MODULE,"file_sha",side_effect=AssertionError("payload alias hashed")):
            with self.assertRaises(MODULE.GuardError):
                MODULE.verify_pin(self.root,{"path":"config/SYN_METADATA.json","sha256":"0"*64})
        directory=self.root/"SYN_LINKED_INPUT_DIR";directory.symlink_to(payload.parent,target_is_directory=True)
        with patch.object(MODULE,"file_sha",side_effect=AssertionError("payload ancestry alias hashed")):
            with self.assertRaises(MODULE.GuardError):
                MODULE.verify_pin(self.root,{"path":"SYN_LINKED_INPUT_DIR/"+payload.name,"sha256":"0"*64})

    def test_registered_metadata_hardlink_payload_alias_refused_before_hash(self):
        payload=self.root/MODULE.MATRIX_PIN["path"]
        payload.parent.mkdir(parents=True,exist_ok=True)
        payload.write_bytes(b"SYNTHETIC_NOT_REAL_EXPRESSION")
        alias=self.root/"config/SYN_METADATA.json";alias.hardlink_to(payload)
        with patch.object(MODULE,"file_sha",side_effect=AssertionError("payload hardlink hashed")):
            with self.assertRaises(MODULE.GuardError):
                MODULE.verify_pin(self.root,{"path":"config/SYN_METADATA.json","sha256":"0"*64})

    def test_plan_hardlinks_to_either_payload_refused_before_json_read(self):
        for pin in (MODULE.MATRIX_PIN,MODULE.RAW_MATRIX_PIN):
            payload=self.root/pin["path"]
            payload.parent.mkdir(parents=True,exist_ok=True)
            payload.write_bytes(b"SYNTHETIC_NOT_REAL_EXPRESSION")
            self.plan_path.unlink()
            self.plan_path.hardlink_to(payload)
            with patch.object(MODULE,"read_json",side_effect=AssertionError("plan alias body opened")):
                self.assertEqual(MODULE.main(["--preflight","--plan",str(self.plan_path)],self.root),2)

    def test_proposal_implementation_hardlinks_refused_before_hash(self):
        original=MODULE.file_sha
        for rel in (MODULE.SCRIPT_REL,MODULE.TEST_REL):
            payload=self.root/MODULE.MATRIX_PIN["path"]
            payload.parent.mkdir(parents=True,exist_ok=True)
            payload.write_bytes(b"SYNTHETIC_NOT_REAL_EXPRESSION")
            implementation=self.root/rel
            implementation.unlink();implementation.hardlink_to(payload)
            def guarded_hash(path):
                if Path(path).samefile(payload):raise AssertionError("implementation alias body hashed")
                return original(path)
            with patch.object(MODULE,"file_sha",side_effect=guarded_hash):
                self.assertEqual(MODULE.main(["--write-proposal"],self.root),2)
            implementation.unlink();implementation.write_text("# Generated synthetic file\n")

    def test_changed_metadata_pin_refused_before_loader(self):
        with patch.object(MODULE,"source_only_preflight",side_effect=MODULE.GuardError("changed metadata pin")):
            with self.assertRaises(MODULE.GuardError):self.invoke_guards()
        self.assertFalse(self.output.exists())

    def test_changed_compiled_memberships_refused_before_payload_hash(self):
        with patch.object(MODULE,"source_only_preflight",return_value=({},[{"changed":True}],[],[])):
            with self.assertRaises(MODULE.GuardError):self.invoke_guards()
        self.assertFalse(self.output.exists())

    def test_changed_payload_pin_before_numeric_loader_or_output_creation(self):
        original=MODULE.verify_pin
        def verify(root,pin):
            if pin["path"]==MODULE.RAW_MATRIX_PIN["path"]:
                raise MODULE.GuardError("changed raw payload pin")
            return original(root,pin)
        with patch.object(MODULE,"source_only_preflight",return_value=({},[],[],[])),patch.object(MODULE,"verify_pin",side_effect=verify):
            with self.assertRaises(MODULE.GuardError):self.invoke_guards()
        self.assertFalse(self.output.exists())

    def test_complete_guard_fixture_never_loads_values_or_creates_output(self):
        original=MODULE.verify_pin
        calls=[]
        def verify(root,pin):
            calls.append(pin["path"])
            if pin["path"] in (MODULE.RAW_MATRIX_PIN["path"],MODULE.MATRIX_PIN["path"]):
                return self.root/"SYNTHETIC_MATRIX_NEVER_OPENED.npy"
            return original(root,pin)
        with patch.object(MODULE,"source_only_preflight",return_value=({},[],[],[])),patch.object(MODULE,"verify_pin",side_effect=verify):
            result=self.invoke_guards()
        self.assertEqual(result[4],self.output)
        self.assertEqual(calls[-2:],[MODULE.RAW_MATRIX_PIN["path"],MODULE.MATRIX_PIN["path"]])
        self.assertFalse(self.output.exists())

    def test_planning_output_rejects_symlink_ancestors_before_writes(self):
        fake=self.root/"SYNTHETIC_PLANNING_ROOT"
        (fake/"reports").mkdir(parents=True)
        (fake/"work").symlink_to(fake/"reports",target_is_directory=True)
        with patch.object(MODULE,"draft_plan",side_effect=AssertionError("planner reached through protected symlink")):
            self.assertEqual(MODULE.main(["--write-proposal"],fake),2)
        self.assertEqual(list((fake/"reports").iterdir()),[])

    def test_plan_argument_cannot_open_matrix_or_arbitrary_json(self):
        for path in (self.root/MODULE.MATRIX_PIN["path"],self.root/MODULE.RAW_MATRIX_PIN["path"],
                     self.root/"work/old-column-qc.json"):
            with patch.object(MODULE,"read_json",side_effect=AssertionError("payload read before path guard")):
                self.assertEqual(MODULE.main(["--preflight","--plan",str(path)],self.root),2)
        link=self.root/"config/psc-liver-program-freeze.json"
        link.symlink_to(self.plan_path)
        with self.assertRaises(MODULE.GuardError):MODULE.checked_plan_path(link,self.root)

    def test_proposal_records_root_choices_without_execution_authorization(self):
        proposal=MODULE.draft_plan(self.root)
        self.assertEqual(proposal["state"],"draft_not_frozen")
        self.assertTrue(proposal["panel_approved"])  # prospective root panel decision, not execution
        self.assertEqual(len(proposal["panel"]),10)
        self.assertFalse(proposal["root_authorization"]["real_execution_authorized"])
        self.assertIsNone(proposal["expected_compiled_programs_sha256"])
        self.assertEqual(proposal["method"]["bootstrap"]["seed"],2026091602)
        self.assertEqual(proposal["method"]["zero_variance_policy"],"welch_if_one_positive")
        self.assertTrue(proposal["unit"]["assumption_approved"])
        self.assertEqual(proposal["root_method_decision_reference"],MODULE.ROOT_METHOD_DECISION)
        self.assertTrue(all(d["status"]=="approved" for d in proposal["decisions"].values()))

    def test_noop_cli_and_incomplete_execution_cannot_load(self):
        with patch.object(MODULE,"load_count_matrix",side_effect=AssertionError("must not load")) as loader:
            with self.assertRaises(SystemExit):MODULE.main([],self.root)
            self.assertEqual(MODULE.main(["--execute"],self.root),2)
        self.assertEqual(loader.call_count,0)

    def test_json_missing_values_are_null_and_nonfinite_or_duplicate_keys_rejected(self):
        target=self.root/"test.json"
        MODULE.dump_json(target,{"p_raw":None})
        self.assertEqual(MODULE.read_json(target),{"p_raw":None})
        for text in ('{"p":NaN}','{"p":1,"p":2}'):
            target.write_text(text)
            with self.assertRaises(MODULE.GuardError):MODULE.read_json(target)
        with self.assertRaises(ValueError):MODULE.dump_json(target,{"p":np.nan})

    def test_failed_execution_receipt_cannot_be_complete(self):
        approved=(self.plan,[],[],[],self.output,self.root/"SYN_NO_VALUES.npy",{})
        with patch.object(MODULE,"authorize_execution",return_value=approved),patch.object(MODULE,"load_count_matrix",side_effect=MODULE.DataIntegrityError("synthetic refused load; no real values")):
            with self.assertRaises(MODULE.DataIntegrityError):
                MODULE.execute(self.plan_path,self.sha,MODULE.EXECUTION_ACK,self.output,self.root)
        failed=MODULE.read_json(self.output/"FAILED.json")
        self.assertEqual(failed["status"],"failed_not_complete")
        self.assertFalse(failed["access_audit"]["real_expression_values_loaded"])
        self.assertFalse((self.output/"completion-receipt.json").exists())
        self.assertEqual(len(failed["endpoint_status_rows"]),len(self.plan["panel"])*3*2)


class SourceOnlyPreflightTests(unittest.TestCase):
    def test_source_only_current_metadata_never_calls_math_or_real_loader(self):
        with patch.object(MODULE,"load_count_matrix",side_effect=AssertionError("real loader reached")) as loader, \
             patch.object(MODULE,"normalize_counts",side_effect=AssertionError("normalization reached")), \
             patch.object(MODULE,"analyze_counts",side_effect=AssertionError("analysis reached")):
            original_file_sha=MODULE.file_sha
            def safe_hash(path):
                if str(path).endswith(("GSE303271-submitted-float64.npy","GSE303271_raw_counts.txt.gz")):
                    raise AssertionError("Source-only preflight opened expression payload bytes")
                return original_file_sha(path)
            with patch.object(MODULE,"file_sha",side_effect=safe_hash):
                receipt,specs,units,audit=MODULE.source_only_preflight(root=ROOT)
        self.assertEqual(loader.call_count,0)
        self.assertFalse(receipt["real_expression_payload_opened"])
        self.assertFalse(receipt["real_analysis_executed"])
        self.assertFalse(receipt["root_freeze_written"])
        self.assertEqual(receipt["prospectively_approved_core_family_size"],30)
        self.assertEqual(receipt["supplied_plan_state"],"draft_not_frozen")
        self.assertFalse(receipt["execution_authorization_checked"])
        self.assertNotIn("frozen_discovery_family_size",receipt)
        self.assertEqual(receipt["external_target_gene_membership"][MODULE.ECM_PROGRAM],
            {"literal_symbol":"ETS2","present_in_original_source":False,"excluded":False})
        self.assertEqual(receipt["groups"],{"PSC":17,"ASC":17,"AIH":30})
        self.assertEqual(len(specs),20)
        self.assertEqual(len(receipt["endpoint_status_rows"]),60)
        self.assertEqual(len(audit),642)
        self.assertEqual(receipt["external_membership_rows_read"],321)
        self.assertTrue(receipt["does_not_fulfill_held_IL17_response_transfer"])
        self.assertEqual(len(units),64)
        self.assertTrue(all(row["p_raw"] is None and row["p_for_correction"] is None for row in receipt["endpoint_status_rows"]))




class ExternalSourceGuardTests(unittest.TestCase):
    def setUp(self):
        base=ROOT/MODULE.WORK_REL/"synthetic-tests"
        base.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=base)
        self.root=Path(self.temp.name)
        self.endpoint=deepcopy(MODULE.EXTERNAL_SOURCE_REGISTRY[MODULE.ECM_PROGRAM])
        self.endpoint.update(kind="external_symbol_csv",status="qualified",source_approval_reference="SYNTHETIC_TEST_ONLY")
        # Source-schema scalars only, no actual gene or effect values.
        self.handoff={"schema_version":1,"program_id":MODULE.ECM_PROGRAM,
            "config_path":self.endpoint["source_manifest_pin"]["path"],
            "final_config_sha256":self.endpoint["source_manifest_pin"]["sha256"],
            "membership_path":self.endpoint["membership_pin"]["path"],
            "membership_sha256":self.endpoint["membership_pin"]["sha256"],
            "membership_bytes":125834,"membership_rows":321,"complete_original_denominator":321,
            "original_source_identifier_count":340,"source_direction":"UNSIGNED",
            "source_version":self.endpoint["source_version"],"scoring_or_endpoint_authorized":False,
            "source_qualification_status":"FINAL_partial_ECM_qualified_epithelial_IL17A_response_HOLD"}
        self.handoff_path=self.root/"SYNTHETIC_SCHEMA_ONLY.json"
        MODULE.dump_json(self.handoff_path,self.handoff)
        self.manifest_path=self.root/"SYNTHETIC_QUALIFIED_MANIFEST.json"
        self.manifest={"schema_version":1,"qualification_status":"partial_ECM_qualified_epithelial_IL17A_response_HOLD",
            "source_pass_closed":True,"membership_output":{**self.endpoint["membership_pin"],"rows":321},
            "programs":[{"program_id":MODULE.ECM_PROGRAM,"status":"qualified_complete_official_pathway_membership_source_only",
                "official_program_name":"REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION","msigdb_systematic_id":"M610",
                "exact_source_id":"R-HSA-1474244","msigdb_version":"2026.1.Hs","reactome_release":"95",
                "collection":"C2:CP:REACTOME","organism":"Homo sapiens","source_identifier_namespace":"Human_Ensembl_Gene_ID",
                "member_namespace":"MSigDB_geneSymbols","expected_member_count":321,"expected_source_identifier_count":340,
                "direction":"UNSIGNED","membership_weight":"1","published_response_directions_available":False,
                "scoring_authorized":False,"source_keys":{"html":"ecm_html","json":"ecm_json"}}],
            "sources":{key:{**pin,"complete":True,"status":200} for key,pin in MODULE.EXTERNAL_RAW_SOURCE_PINS.items()}}
        MODULE.dump_json(self.manifest_path,self.manifest)

    def tearDown(self):self.temp.cleanup()

    def mock_pin(self,root,pin):
        if pin["path"]==self.endpoint["schema_handoff_pin"]["path"]:return self.handoff_path
        if pin["path"]==self.endpoint["source_manifest_pin"]["path"]:return self.manifest_path
        return self.root/"SYNTHETIC_SOURCE_NO_REAL_ROWS.csv"

    def test_qualified_handoff_binds_manifest_membership_and_semantics(self):
        with patch.object(MODULE,"verify_pin",side_effect=self.mock_pin) as verify:
            result=MODULE.validate_external_source(self.endpoint,self.root)
        self.assertEqual(result,self.root/"SYNTHETIC_SOURCE_NO_REAL_ROWS.csv")
        self.assertEqual(verify.call_count,5)

    def test_count_or_arbitrary_paths_refused_before_byte_access(self):
        for pin_name in ("source_manifest_pin","membership_pin"):
            for path in (MODULE.MATRIX_PIN["path"],MODULE.RAW_MATRIX_PIN["path"],"work/arbitrary.csv"):
                endpoint=deepcopy(self.endpoint);endpoint[pin_name]["path"]=path
                with patch.object(MODULE,"verify_pin",side_effect=AssertionError("must reject before byte access")):
                    with self.assertRaises(MODULE.GuardError):MODULE.validate_external_source(endpoint,self.root)

    def test_wrong_source_version_namespace_denominator_or_polarity_refused(self):
        changes={"source_version":"OTHER","source_identifier_namespace":"published_Entrez",
                 "original_component_members":{"all":320},"direction_literals":{"all":"UP"},
                 "score_kind":"balanced_signed","target_gene_exclusions":["ETS2"]}
        for key,value in changes.items():
            endpoint=deepcopy(self.endpoint);endpoint[key]=value
            with patch.object(MODULE,"verify_pin",side_effect=AssertionError("source opened before descriptor bound")):
                with self.assertRaises(MODULE.GuardError):MODULE.validate_external_source(endpoint,self.root)

    def test_unregistered_or_held_signature_not_an_empty_program(self):
        endpoint=deepcopy(self.endpoint);endpoint["program_id"]="epithelial_il17a_alone_response"
        with patch.object(MODULE,"verify_pin",side_effect=AssertionError("unqualified input read")):
            with self.assertRaises(MODULE.GuardError):MODULE.validate_external_source(endpoint,self.root)

    def test_handoff_mismatched_manifest_or_csv_pin_refused(self):
        for key in ("membership_sha256","final_config_sha256","source_qualification_status"):
            changed=dict(self.handoff);changed[key]="changed"
            MODULE.dump_json(self.handoff_path,changed)
            with patch.object(MODULE,"verify_pin",side_effect=self.mock_pin):
                with self.assertRaises(MODULE.GuardError):MODULE.validate_external_source(self.endpoint,self.root)

    def test_manifest_semantics_are_checked_not_only_handoff(self):
        for key,value in (("direction","UP"),("expected_member_count",320),("member_namespace","Entrez")):
            changed=deepcopy(self.manifest);changed["programs"][0][key]=value
            MODULE.dump_json(self.manifest_path,changed)
            with patch.object(MODULE,"verify_pin",side_effect=self.mock_pin):
                with self.assertRaises(MODULE.GuardError):MODULE.validate_external_source(self.endpoint,self.root)

    def test_manifest_csv_binding_and_complete_raw_sources(self):
        for modification in ("csv","raw"):
            changed=deepcopy(self.manifest)
            if modification=="csv":changed["membership_output"]["path"]=MODULE.MATRIX_PIN["path"]
            else:changed["sources"]["ecm_html"]["complete"]=False
            MODULE.dump_json(self.manifest_path,changed)
            with patch.object(MODULE,"verify_pin",side_effect=self.mock_pin):
                with self.assertRaises(MODULE.GuardError):MODULE.validate_external_source(self.endpoint,self.root)

    def test_json_dictionary_key_order_does_not_define_signed_polarity(self):
        rows,features,endpoint=MetadataMappingTests().fixture()
        endpoint.update(score_kind="balanced_signed",original_component_members={"positive":6,"negative":6},
                        direction_literals={"positive":"UP","negative":"DOWN"})
        for i,row in enumerate(rows):row["published_direction"]="UP" if i<6 else "DOWN"
        endpoint=json.loads(json.dumps(endpoint,sort_keys=True))
        specs,audit=MODULE.map_external_rows(rows,endpoint,features)
        self.assertEqual([c["name"] for c in specs[(CORE,"SYN_EXT")]["components"]],["positive","negative"])
        self.assertEqual(len(audit),24)


if __name__ == "__main__":
    unittest.main()
