"""Synthetic tests of the independent verifier; no real RNA values are opened."""
from __future__ import annotations
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock

PATH = Path(__file__).resolve().parents[1]/"scripts/verify_psc_liver_programs.py"
spec = importlib.util.spec_from_file_location("independent_liver_program_verifier", PATH)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


class TestIndependentNormalization(unittest.TestCase):
    def test_full_feature_universe_not_program_subset(self):
        result = v.logcpm_reference([[10, 20], [90, 180]])
        self.assertEqual(result["library_totals"], [100, 200])
        expected = math.log2(100001)
        self.assertAlmostEqual(result["logcpm"][0][0], expected)
        self.assertAlmostEqual(result["logcpm"][0][1], expected)
        self.assertNotEqual(expected, v.logcpm_reference([[10, 20]])["logcpm"][0][0])

    def test_zero_feature_is_measured_zero_and_not_removed(self):
        result = v.logcpm_reference([[0, 0], [10, 20]])
        self.assertEqual(result["feature_rows"], 2)
        self.assertEqual(result["logcpm"][0], [0, 0])
        self.assertEqual(v.score_reference(result["logcpm"], [0]), [0, 0])

    def test_tiny_expression_remains_positive(self):
        result = v.logcpm_reference([[1], [2**52]])
        expected = math.log1p(1e6/(2**52+1))/math.log(2)
        self.assertGreater(result["logcpm"][0][0], 0)
        self.assertEqual(result["logcpm"][0][0], expected)

    def test_invalid_counts_never_recode_or_drop(self):
        for bad in (-1, 0.5, float("nan"), float("inf"), True, 2**54):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                v.logcpm_reference([[bad, 2], [2, 3]])

    def test_integer_literal_above_float_safety_cannot_round_down(self):
        with self.assertRaises(ValueError):v.logcpm_reference([[2**53+1]])

    def test_library_total_above_float_safety_cannot_round_down(self):
        with self.assertRaises(ValueError):v.logcpm_reference([[2**53],[1]])
        result=v.logcpm_reference([[2**53-1],[1]])
        self.assertEqual(result["library_totals"],[float(2**53)])

    def test_zero_library_rejected(self):
        with self.assertRaisesRegex(ValueError, "library"):
            v.logcpm_reference([[0, 1], [0, 2]])

    def test_empty_ragged_axes_rejected(self):
        for matrix in ([], [[]], [[1, 2], [3]]):
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                v.logcpm_reference(matrix)

    def test_exact_80_percent_gate_and_minimum(self):
        self.assertTrue(v.mapping_gate(20, 16))
        self.assertFalse(v.mapping_gate(20, 15))
        self.assertFalse(v.mapping_gate(10, 9))
        self.assertTrue(v.mapping_gate(927, 848))

    def test_bad_denominators_rejected(self):
        for denominator, mapped in ((0, 0), (10, 11), (10, -1), (True, 1)):
            with self.subTest(denominator=denominator), self.assertRaises(ValueError):
                v.mapping_gate(denominator, mapped)

    def test_unsigned_means_are_not_source_direction_inverted(self):
        matrix = [[2, 4], [6, 8], [10, 12]]
        self.assertEqual(v.score_reference(matrix, [0, 1]), [4, 6])
        self.assertEqual(v.score_reference(matrix, [1, 0]), [4, 6])

    def test_signed_math_explicit_equal_component_not_equal_all_gene(self):
        matrix = [[2, 4], [6, 8], [10, 12]]
        self.assertEqual(v.score_reference(matrix, [0, 1], [2]), [-3, -3])
        self.assertNotEqual(v.score_reference(matrix, [0, 1], [2]), [-2/3, 0])

    def test_absent_conflicting_repeated_components_rejected(self):
        for positive, negative in (([], None), ([0], []), ([0, 0], None),
                                    ([0], [0]), ([4], None), ([-1], None), ([True], None)):
            with self.subTest(positive=positive, negative=negative), self.assertRaises(ValueError):
                v.score_reference([[2, 3], [5, 6]], positive, negative)


class TestIndependentWelch(unittest.TestCase):
    def test_unequal_sizes_variances_hand_formula(self):
        r = v.welch_reference([1,2,3], [3,5,7,9])
        self.assertEqual(r["effect"], -4)
        self.assertAlmostEqual(r["se"], math.sqrt(2))
        self.assertAlmostEqual(r["df"], 216/53)
        self.assertAlmostEqual(r["t"], -math.sqrt(8))
        self.assertLess(r["ci_low"], r["effect"])
        self.assertGreater(r["ci_high"], r["effect"])
        self.assertTrue(0 < r["p_raw"] < 1)

    def test_one_zero_arm_uses_other_degrees_of_freedom(self):
        r = v.welch_reference([2,2], [0,2])
        self.assertEqual(r["status"], "ok")
        self.assertAlmostEqual(r["df"], 1)
        self.assertAlmostEqual(r["t"], 1)
        # Cauchy(df1) has a closed-form upper tail.
        self.assertAlmostEqual(r["p_raw"], 1-2*math.atan(1)/math.pi)
        self.assertAlmostEqual(r["ci_high"], 1+math.tan(math.pi*0.475), places=8)

    def test_both_zero_variances_equal_or_different_stay_missing(self):
        for right in ([1,1], [2,2]):
            r = v.welch_reference([1,1], right)
            self.assertEqual(r["status"], "zero_standard_error")
            for field in ("p_raw", "df", "t", "ci_low", "ci_high"):
                self.assertIsNone(r[field])
            self.assertIsNotNone(r["effect"])

    def test_constant_nonbinary_float_does_not_invent_variance(self):
        r = v.welch_reference([2.1]*30, [2.1]*17)
        self.assertEqual(r["sd_left"], 0)
        self.assertEqual(r["sd_right"], 0)
        self.assertIsNone(r["p_raw"])

    def test_tiny_positive_variance_does_not_trigger_epsilon_cutoff(self):
        r = v.welch_reference([x*1e-200 for x in [1,2,3]],
                              [x*1e-200 for x in [3,5,7,9]])
        ordinary = v.welch_reference([1,2,3], [3,5,7,9])
        self.assertEqual(r["status"], "ok")
        self.assertAlmostEqual(r["df"], ordinary["df"])
        self.assertAlmostEqual(r["p_raw"], ordinary["p_raw"])
        self.assertTrue(math.isclose(r["se"], ordinary["se"]*1e-200, rel_tol=1e-14))

    def test_symmetry_preserves_p_and_reverses_effect_and_interval(self):
        a, b = [1,2,3], [3,5,7,9]
        x, y = v.welch_reference(a,b), v.welch_reference(b,a)
        self.assertAlmostEqual(x["p_raw"], y["p_raw"])
        self.assertAlmostEqual(x["effect"], -y["effect"])
        self.assertAlmostEqual(x["ci_low"], -y["ci_high"])

    def test_one_unit_is_not_an_estimated_population_variance(self):
        r = v.welch_reference([1], [2,3])
        self.assertEqual(r["status"], "insufficient_units")
        self.assertIsNone(r["p_raw"])
        self.assertEqual(r["effect"], -1.5)

    def test_nonfinite_scores_rejected(self):
        for bad in (float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                v.welch_reference([1,bad], [2,3])


class TestIndependentCorrection(unittest.TestCase):
    def test_missing_only_in_correction_full_denominator(self):
        raw = [0.01, None, 0.04]
        result = v.adjust_reference(raw)
        self.assertEqual(raw, [0.01, None, 0.04])
        self.assertEqual(result["family_size"], 3)
        self.assertEqual(result["p_for_correction"], [0.01, 1, 0.04])
        for a,b in zip(result["q_bh"], [0.03, 1, 0.06]): self.assertAlmostEqual(a,b)
        for a,b in zip(result["q_by"], [0.055, 1, 0.11]): self.assertAlmostEqual(a,b)

    def test_ties_zeros_ones_and_order(self):
        p = [1, 0, 0.01, 0.01]
        q = v.adjust_reference(p)
        self.assertEqual(q["q_bh"][0:2], [1,0])
        self.assertEqual(q["q_bh"][2], q["q_bh"][3])
        self.assertAlmostEqual(q["q_bh"][2], 0.04/3)

    def test_invalid_p_is_not_missing(self):
        for p in ([], [float("nan")], [-0.01], [1.01], [True]):
            with self.subTest(p=p), self.assertRaises(ValueError): v.adjust_reference(p)

    @staticmethod
    def family(interface=v.CORE):
        rows=[]
        raw=[0.01, 0.2, 0.3, None, None, None]
        correction=v.adjust_reference(raw)
        for i,(program, (a,b)) in enumerate((p,c) for p in ("A","HOLD") for c in v.CONTRASTS):
            rows.append({"interface":interface, "program_id":program,"left":a,"right":b,
                         "status":"ok" if raw[i] is not None else "unavailable_source", "p_raw":raw[i],
                         "q_bh":None if raw[i] is None or interface==v.STRICT else correction["q_bh"][i],
                         "q_by":None if raw[i] is None or interface==v.STRICT else correction["q_by"][i]})
        return rows

    def test_full_status_grid_including_source_hold(self):
        self.assertEqual(v.assert_family(self.family(), ["A","HOLD"]), 6)

    def test_missing_extra_duplicate_and_reordered_rows_rejected(self):
        rows=self.family()
        for bad in (rows[:-1], rows+[rows[0]], rows[1:]+rows[:1], rows[:-1]+[rows[0]]):
            with self.subTest(length=len(bad)), self.assertRaises(ValueError):
                v.assert_family(bad,["A","HOLD"])

    def test_selected_family_q_is_rejected(self):
        rows=self.family()
        rows[0]["q_bh"]=0.03
        with self.assertRaises(ValueError): v.assert_family(rows,["A","HOLD"])

    def test_missing_p_not_replaced_in_public_table(self):
        rows=self.family(); rows[-1]["p_raw"]=1
        with self.assertRaises(ValueError): v.assert_family(rows,["A","HOLD"])

    def test_strict_has_no_discovery_q(self):
        rows=self.family(v.STRICT)
        v.assert_family(rows,["A","HOLD"],v.STRICT)
        rows[0]["q_bh"]=0.1
        with self.assertRaises(ValueError): v.assert_family(rows,["A","HOLD"],v.STRICT)


class TestIndependentResampling(unittest.TestCase):
    def setUp(self):
        self.groups=["PSC","AIH","ASC","PSC","ASC","AIH","AIH"]
        self.scores=[1,5,3,2,4,6,7]
        self.draws=v.stratified_indices(self.groups,31,716)

    def test_strata_sizes_and_deterministic_replay(self):
        self.assertEqual(self.draws,v.stratified_indices(self.groups,31,716))
        for group,draws in self.draws.items():
            self.assertEqual(len(draws),31)
            for row in draws:
                self.assertEqual(len(row),self.groups.count(group))
                self.assertTrue(all(self.groups[i]==group for i in row))

    def test_linear_percentile_hand_interpolation(self):
        self.assertAlmostEqual(v.percentile_linear([0,10,20,30],0.025),0.75)
        self.assertEqual(v.percentile_linear([0,10,20,30],0.975),29.25)
        self.assertEqual(v.percentile_linear([3],0.025),3)

    def test_joint_draws_preserve_affine_program_relationship(self):
        a=v.bootstrap_reference(self.scores,self.groups,self.draws,"PSC","AIH")
        b=v.bootstrap_reference([2*x+9 for x in self.scores],self.groups,self.draws,"PSC","AIH")
        for x,y in zip(a["draws"],b["draws"]): self.assertAlmostEqual(y,2*x)
        self.assertAlmostEqual(b["ci_low"],2*a["ci_low"])

    def test_constant_draws_are_retained_not_filtered(self):
        r=v.bootstrap_reference([0]*7,self.groups,self.draws,"PSC","AIH")
        self.assertEqual(r["draws"],[0]*31)
        self.assertEqual(r["ci_low"],0)
        self.assertEqual(r["ci_high"],0)

    def test_cross_stratum_or_wrong_width_is_rejected(self):
        draws=copy.deepcopy(self.draws); draws["PSC"][0][0]=1
        with self.assertRaises(ValueError):
            v.bootstrap_reference(self.scores,self.groups,draws,"PSC","AIH")
        draws=copy.deepcopy(self.draws); draws["PSC"][0].pop()
        with self.assertRaises(ValueError):
            v.bootstrap_reference(self.scores,self.groups,draws,"PSC","AIH")

    def test_omits_every_unit_including_unchanged_asc(self):
        ids=[f"S{i}" for i in range(7)]
        rows=v.omissions_reference(self.scores,self.groups,ids,"PSC","AIH")
        self.assertEqual([r["unit_id"] for r in rows],ids)
        baseline=1.5-6
        for row in rows:
            if row["omitted_group"]=="ASC":
                self.assertEqual(row["status"],"not_in_contrast")
                self.assertEqual(row["effect"],baseline)
        self.assertEqual(rows[0]["effect"],2-6)

    def test_repeated_unit_id_fails_without_silent_pseudoreplication(self):
        with self.assertRaises(ValueError):
            v.omissions_reference(self.scores,self.groups,["same"]*7,"PSC","AIH")


class TestIndependentGuards(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.work=self.root/"work/psc-liver-program-plan"
        self.work.mkdir(parents=True)
        self.source=self.root/"source.txt"; self.source.write_text("source")
        self.pin={"bytes":6,"sha256":hashlib.sha256(b"source").hexdigest()}

    def tearDown(self): self.temp.cleanup()

    def test_source_pin_read_only_and_mutation_detected(self):
        v.check_pin(self.root,"source.txt",self.pin)
        self.assertEqual(self.source.read_text(),"source")
        self.source.write_text("other!")
        with self.assertRaises(ValueError): v.check_pin(self.root,"source.txt",self.pin)

    def test_fresh_work_child_only(self):
        output=self.work/"new"
        v.assert_fresh_private(self.root,output,[self.source])
        self.assertFalse(output.exists())
        output.mkdir()
        with self.assertRaises(ValueError): v.assert_fresh_private(self.root,output,[self.source])
        with self.assertRaises(ValueError):
            v.assert_fresh_private(self.root,self.root/"data/derived/public",[self.source])

    def test_symlink_input_and_output_parent_rejected(self):
        link=self.root/"link.txt"; link.symlink_to(self.source)
        with self.assertRaises(ValueError): v.check_pin(self.root,"link.txt",self.pin)
        alias=self.work/"alias"; alias.symlink_to(self.root, target_is_directory=True)
        with self.assertRaises(ValueError): v.assert_fresh_private(self.root,alias/"new",[self.source])

    def test_path_escape_and_source_output_overlap_rejected(self):
        for path in ("../source.txt",str(self.source)):
            with self.assertRaises(ValueError): v.relative_file(self.root,path)
        with self.assertRaises(ValueError):
            v.assert_fresh_private(self.root,self.work/"new",[self.work/"new/counts"])

    @staticmethod
    def approval_fixture():
        plan={"schema_version":1,"frozen":True,"real_expression_authorized":True,
              "normalization_feature_rows":41096,"primary_interface":v.CORE,"sensitivity_interface":v.STRICT,
              "contrasts":[list(c) for c in v.CONTRASTS],
              "programs":[{"program_id":"A","source_status":"qualified","source_sha256":"f"*64},
                          {"program_id":"HOLD","source_status":"HOLD","members":None,
                           "effect_status":"unavailable_source"}]}
        data=json.dumps(plan).encode()
        approval={"root_approved":True,"real_expression_authorized":True,
                  "plan_sha256":hashlib.sha256(data).hexdigest()}
        return plan,approval,data

    def test_frozen_byte_bound_approval_oracle(self):
        plan,approval,data=self.approval_fixture()
        v.assert_approval(plan,approval,data)
        for key in ("root_approved","real_expression_authorized"):
            bad=copy.deepcopy(approval); bad[key]=False
            with self.assertRaises(ValueError): v.assert_approval(plan,bad,data)
        with self.assertRaises(ValueError): v.assert_approval(plan,approval,data+b" ")

    def test_hold_not_empty_or_fallback_and_universe_not_subset(self):
        plan,approval,data=self.approval_fixture()
        plan["programs"][1]["members"]=[]
        with self.assertRaises(ValueError): v.assert_approval(plan,approval,data)
        plan,approval,data=self.approval_fixture(); plan["normalization_feature_rows"]=37778
        with self.assertRaises(ValueError): v.assert_approval(plan,approval,data)

    def test_source_only_rejects_expanded_whitelist_before_opening_counts(self):
        scope=self.root/"scope.json"
        scope.write_text(json.dumps({"source_summary_pins":{"real-counts.npy":self.pin}}))
        with self.assertRaisesRegex(ValueError,"whitelist"):
            v.source_only(self.root,scope)

    def test_cli_has_no_real_counts_or_execution_mode(self):
        with self.assertRaises(SystemExit): v.main(["--execute"])
        with self.assertRaises(SystemExit): v.main(["--matrix",str(self.source)])

    def test_builtin_synthetic_never_opens_source_values(self):
        with mock.patch.object(Path,"read_bytes",side_effect=AssertionError("file open")), \
             mock.patch.object(Path,"read_text",side_effect=AssertionError("file open")):
            result=v.builtin_synthetic()
        self.assertTrue(result["no_real_count_or_score_rows_opened"])
        self.assertEqual(result["full_family_size_including_unavailable"],6)

    def test_verifier_does_not_import_author(self):
        import ast
        tree=ast.parse(PATH.read_text())
        names=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import): names.extend(x.name for x in node.names)
            if isinstance(node,ast.ImportFrom): names.append(node.module or "")
        self.assertFalse(any("analyze_psc_liver_programs" in name for name in names))


class TestSyntheticNumericalBundle(unittest.TestCase):
    @staticmethod
    def fixture():
        counts=[[1,2,4,3,5,8],[9,8,6,7,5,2],[90]*6]
        groups=["PSC","PSC","ASC","ASC","AIH","AIH"]
        ids=[f"S{i}" for i in range(6)]
        normal=v.logcpm_reference(counts)
        indices=v.stratified_indices(groups,7,716)
        programs=[]; observed=[]; inference=[]
        for tier in (v.CORE,v.STRICT):
            for name,eligible in (("A",True),("MISSING",False)):
                programs.append({"interface":tier,"program_id":name,"coverage_available":eligible,
                                 "positive":[0],"negative":None})
                scores=v.score_reference(normal["logcpm"],[0]) if eligible else None
                rows=[]
                for left,right in v.CONTRASTS:
                    stat=(v.welch_reference([scores[i] for i,g in enumerate(groups) if g==left],
                                           [scores[i] for i,g in enumerate(groups) if g==right])
                          if eligible else None)
                    rows.append({"left":left,"right":right,"welch":stat,
                                 "bootstrap":v.bootstrap_reference(scores,groups,indices,left,right) if eligible else None,
                                 "omissions":v.omissions_reference(scores,groups,ids,left,right) if eligible else None})
                    if eligible:
                        rows[-1]["omission_summary"]=v.omission_summary_reference(rows[-1]["omissions"],stat["effect"])
                    else:
                        n_in=groups.count(left)+groups.count(right)
                        rows[-1]["omission_summary"]={"in_contrast_rows":n_in,"stability_denominator":n_in,
                            "defined_in_contrast_effect_rows":0,"effect_min":None,"effect_max":None,
                            "max_absolute_effect_change":None,"strict_sign_reversal_rows":None}
                    inference.append({"interface":tier,"program_id":name,"left":left,"right":right,
                                      "status":stat["status"] if eligible else "unavailable_score",
                                      "p_raw":stat["p_raw"] if eligible else None})
                observed.append({"interface":tier,"program_id":name,"scores":scores,"contrasts":rows})
        adjusted=v.adjust_reference([r["p_raw"] for r in inference if r["interface"]==v.CORE])
        for i,row in enumerate(inference):
            missing=row["interface"]==v.STRICT or row["p_raw"] is None
            row.update(q_bh=None if missing else adjusted["q_bh"][i],
                       q_by=None if missing else adjusted["q_by"][i])
        return {"schema_version":1,"synthetic_only":True,"counts":counts,"groups":groups,
                "unit_ids":ids,"program_order":["A","MISSING"],"programs":programs,"replicates":7,"seed":716,
                "observed":{"normalization":normal,"bootstrap_indices":indices,"programs":observed,
                            "inference":inference,"full_core_correction":adjusted}}

    def test_complete_synthetic_bundle(self):
        result=v.verify_numeric_bundle(self.fixture())
        self.assertEqual(result["full_status_rows"],12)
        self.assertEqual(result["full_Core_family_size"],6)
        self.assertEqual(result["scalar_checks"],684)
        self.assertIs(result["real_expression_reviewed"],False)

    def test_without_synthetic_declaration_refused(self):
        bundle=self.fixture(); bundle["synthetic_only"]=False
        with self.assertRaises(ValueError): v.verify_numeric_bundle(bundle)

    def test_tampered_normalization_score_and_interval_detected(self):
        for mutation in ("normalization","score","interval"):
            bundle=self.fixture()
            if mutation=="normalization":bundle["observed"]["normalization"]["logcpm"][0][0]+=0.001
            elif mutation=="score":bundle["observed"]["programs"][0]["scores"][0]+=0.001
            else:bundle["observed"]["programs"][0]["contrasts"][0]["welch"]["ci_high"]+=0.001
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):v.verify_numeric_bundle(bundle)

    def test_joint_bootstrap_and_omission_errors_detected(self):
        bundle=self.fixture();bundle["observed"]["bootstrap_indices"]["PSC"][0][0]=4
        with self.assertRaises(ValueError):v.verify_numeric_bundle(bundle)
        bundle=self.fixture()
        bundle["observed"]["programs"][0]["contrasts"][0]["omissions"][0]["effect"]+=0.1
        with self.assertRaises(ValueError):v.verify_numeric_bundle(bundle)

    def test_strict_or_unavailable_row_loss_detected(self):
        bundle=self.fixture();bundle["observed"]["programs"].pop()
        with self.assertRaises(ValueError):v.verify_numeric_bundle(bundle)
        bundle=self.fixture();bundle["programs"].pop()
        with self.assertRaises(ValueError):v.verify_numeric_bundle(bundle)

    def test_full_adjustment_vectors_compared(self):
        bundle=self.fixture();bundle["observed"]["full_core_correction"]["p_for_correction"][-1]=0.5
        with self.assertRaises(ValueError):v.verify_numeric_bundle(bundle)


class TestIndependentECMIdentity(unittest.TestCase):
    @staticmethod
    def source():
        return [{"source_value":f"G{i}","source_member_index":str(i+1),
                 "program_id":"reactome_ecm_organization","source_identifier_namespace":"MSigDB_geneSymbols",
                 "membership_weight":"1","published_direction":"UNSIGNED"} for i in range(12)]
    @staticmethod
    def features():
        return [{"label_as_submitted":f"G{i}","source_row_index_0based":i,
                 v.CORE:"eligible",v.STRICT:"eligible" if i<10 else "no_entrez_relation"}
                for i in range(12)]
    def test_identity_sensitivity_never_changes_source_denominator(self):
        r=v.ecm_mapping_reference(self.features(),self.source(),12)
        self.assertEqual(r[v.CORE]["mapped_members"],12)
        self.assertEqual(r[v.STRICT]["mapped_members"],10)
        self.assertEqual(r[v.STRICT]["original_source_members"],12)
        self.assertTrue(r[v.STRICT]["coverage_available"])
    def test_alias_case_miss_not_zero_or_entrez_bridge(self):
        rows=self.source();rows[0]["source_value"]="g0";rows[0]["published_ncbi_gene_id"]="999"
        r=v.ecm_mapping_reference(self.features(),rows,12)
        self.assertIsNone(r[v.CORE]["audit"][0]["feature_index"])
        self.assertEqual(r[v.CORE]["audit"][0]["mapping_status"],"no_exact_submitted_symbol")
    def test_repeated_source_symbol_or_row_index_rejected(self):
        source=self.source();source[-1]["source_value"]="G0"
        with self.assertRaises(ValueError):v.ecm_mapping_reference(self.features(),source,12)
        feature=self.features();feature[-1]["source_row_index_0based"]=0
        with self.assertRaises(ValueError):v.ecm_mapping_reference(feature,self.source(),12)
    def test_source_ids_do_not_multiply_members_or_invert_direction(self):
        source=self.source();source[0]["original_source_ids"]="ENS1;ENS2;ENS3"
        r=v.ecm_mapping_reference(self.features(),source,12)
        self.assertEqual(r[v.CORE]["mapped_members"],12)
        source[0]["published_direction"]="DOWN"
        with self.assertRaises(ValueError):v.ecm_mapping_reference(self.features(),source,12)


class TestRootFixedConventions(unittest.TestCase):
    def test_nonfinite_inference_retains_finite_effect_but_not_p(self):
        row=v.welch_reference([1,1], [1e-310,2e-310])
        self.assertEqual(row["effect"],1)
        self.assertEqual(row["status"],"unavailable_numerical_inference")
        self.assertIsNone(row["p_raw"])
        self.assertIsNone(row["ci_low"])
        self.assertIsNone(row["df"])
        self.assertEqual(v.adjust_reference([row["p_raw"]])["p_for_correction"],[1])
    def test_finite_statistic_tail_underflow_is_missing_not_observed_zero(self):
        row=v.welch_reference([i*1e-150 for i in range(17)],[1]*30)
        self.assertEqual(row["status"],"unavailable_numerical_inference")
        self.assertEqual(row["effect"],-1)
        self.assertIsNone(row["p_raw"])
        self.assertIsNone(row["df"])

    def test_public_omission_denominators_are_47_34_47_not_64(self):
        groups=["PSC"]*17+["ASC"]*17+["AIH"]*30
        scores=list(range(64));ids=[f"S{i}" for i in range(64)]
        for (left,right),expected in zip(v.CONTRASTS,[47,34,47]):
            rows=v.omissions_reference(scores,groups,ids,left,right)
            a=[scores[i] for i,g in enumerate(groups) if g==left]
            b=[scores[i] for i,g in enumerate(groups) if g==right]
            summary=v.omission_summary_reference(rows,sum(a)/len(a)-sum(b)/len(b))
            self.assertEqual(len(rows),64)
            self.assertNotIn("all_unit_rows",summary)
            self.assertEqual(summary["stability_denominator"],expected)
            self.assertEqual(summary["defined_in_contrast_effect_rows"],expected)
    def test_unchanged_third_diagnosis_rows_do_not_affect_stability(self):
        rows=[{"status":"omitted","effect":2},{"status":"omitted","effect":3},
              {"status":"not_in_contrast","effect":-1e20}]
        summary=v.omission_summary_reference(rows,1)
        self.assertEqual(summary["effect_min"],2)
        self.assertEqual(summary["strict_sign_reversal_rows"],0)
        self.assertEqual(summary["stability_denominator"],2)
    def test_sign_reversal_avoids_product_underflow(self):
        summary=v.omission_summary_reference([{"status":"omitted","effect":-1e-200}],1e-200)
        self.assertEqual(summary["strict_sign_reversal_rows"],1)


if __name__=="__main__": unittest.main()
