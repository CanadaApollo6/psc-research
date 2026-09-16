"""Synthetic flag-capture and authentication tests; never rerun real fits here."""
from __future__ import annotations
import copy
import io
import json
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import Mock,patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts import audit_macromap_independent_optimizer as d


def synthetic_arrays():
    return {"train_x_raw":np.zeros((27,5)),"train_class_index":np.tile(np.arange(9,dtype=int),3),
        "train_line":np.repeat(["private_line_A","private_line_B","private_line_C"],9),
        "scaler_mean":np.zeros(5),"scaler_scale":np.ones(5),"scaler_constant":np.ones(5,dtype=bool),
        "predictor_names":np.asarray(["inflammation","interferon_gamma","oxidative_stress","apoptotic_signaling",d.VARIANTS[0]]),
        "baseline_coefficients":np.zeros((5,9)),"extended_coefficients":np.zeros((6,9))}


def fake_solver(x,labels,lines):
    return {"coefficients":np.zeros((x.shape[1]+1,9)),"objective":float(np.log(9)),
        "gradient_inf":5e-10,"optimizer_success":False,"iterations":3}


def synthetic_bundle(root):
    directory=root/"work/macromap-design/synthetic-complete";directory.mkdir(parents=True)
    verifier=root/d.VERIFIER_RELATIVE;verifier.parent.mkdir(parents=True)
    verifier.write_bytes((ROOT/d.VERIFIER_RELATIVE).read_bytes())
    manifest=[];pins={}
    for scope in d.SCOPES:
        for protocol in d.PROTOCOLS:
            for time in d.TIMES:
                for fold in range(5):
                    for variant in d.VARIANTS:
                        relative=f"{scope}/{protocol}-t{time}-fold{fold}/{variant}/fit-arrays.npz"
                        arrays=synthetic_arrays();arrays["predictor_names"]=np.asarray(["inflammation","interferon_gamma","oxidative_stress","apoptotic_signaling",variant])
                        buffer=io.BytesIO();np.savez_compressed(buffer,**arrays);body=buffer.getvalue()
                        path=directory/relative;path.parent.mkdir(parents=True);path.write_bytes(body)
                        sha=d.digest_bytes(body);pins[relative]=sha
                        manifest.append({"identity_scope":scope,"protocol":protocol,"time":time,"fold":fold,
                            "variant":variant,"status":"available","fit_npz":relative,"fit_npz_sha256":sha})
    manifest_path=directory/"fit-manifest.json";manifest_path.write_text(json.dumps(manifest))
    pins["fit-manifest.json"]=d.digest_bytes(manifest_path.read_bytes())
    plan={"analysis_id":"macromap-fixed-program-context-v1","status":"frozen","execute_real_expression":True,
        "runtime":d.runtime_versions(),"thread_environment":{},"output_directory":str(directory.relative_to(root)),
        "code_and_protocol_artifacts":[{"path":d.VERIFIER_RELATIVE,"sha256":d.FROZEN_VERIFIER_SHA256}],
        "script_sha256":"a"*64,"runner_script_sha256":"b"*64,"prepared_summary_sha256":"c"*64}
    plan_path=root/"synthetic-plan.json";plan_path.write_text(json.dumps(plan));plan_sha=d.digest_bytes(plan_path.read_bytes())
    complete={"stage":"execution_complete_all_prespecified_statuses_retained","plan_sha256":plan_sha,
        **{key:plan[key] for key in ["runtime","script_sha256","runner_script_sha256","prepared_summary_sha256"]},
        "artifact_sha256":pins}
    complete_path=directory/"execution-complete.json";complete_path.write_text(json.dumps(complete))
    return {"root":root,"directory":directory,"plan_path":plan_path,"plan_sha":plan_sha,
        "complete_path":complete_path,"complete_sha":d.digest_bytes(complete_path.read_bytes()),"manifest":manifest}


def run_bundle(bundle):
    return d.run_diagnostic(bundle["plan_path"],bundle["plan_sha"],bundle["complete_sha"],d.FROZEN_VERIFIER_SHA256,root=bundle["root"])


class OptimizerCaptureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.frozen=d.load_frozen_verifier(ROOT/d.VERIFIER_RELATIVE,d.FROZEN_VERIFIER_SHA256)

    def test_false_flags_preserved_even_when_frozen_numeric_criteria_pass(self):
        solver=Mock(side_effect=fake_solver)
        module=types.SimpleNamespace(weighted_scaler=self.frozen.weighted_scaler,independent_softmax_fit=solver)
        identity={"identity_scope":d.SCOPES[0],"protocol":d.PROTOCOLS[0],"time":6,"fold":0,"variant":d.VARIANTS[0]}
        rows=d.capture_cell(synthetic_arrays(),identity,module)
        self.assertEqual(solver.call_count,2)
        for row in rows:
            self.assertIs(row["optimizer_success"],False);self.assertTrue(row["frozen_numeric_criteria_pass"])
            self.assertTrue(row["coefficient_agreement"]);self.assertEqual(row["iterations"],3)
        self.assertEqual([call.args[0].shape[1] for call in solver.call_args_list],[4,5])
        self.assertTrue(all(not call.kwargs for call in solver.call_args_list))
        self.assertNotIn("private_line",json.dumps(rows))

    def test_constant_names_verified_but_unit_scale_is_not_always_constant(self):
        arrays=synthetic_arrays()
        arrays["train_x_raw"]=np.zeros((36,5));arrays["train_x_raw"][:,0]=np.repeat([-1.,1.,-1.,1.],9)
        arrays["train_class_index"]=np.tile(np.arange(9,dtype=int),4)
        arrays["train_line"]=np.repeat(["private_line_A","private_line_B","private_line_C","private_line_D"],9)
        arrays["scaler_constant"][0]=False
        identity=dict(zip(d.IDENTITY_FIELDS,[d.SCOPES[0],d.PROTOCOLS[0],6,0,d.VARIANTS[0]]))
        rows=d.capture_cell(arrays,identity,self.frozen)
        self.assertEqual([row["constant_training_feature_count"] for row in rows],[3,4])
        self.assertTrue(all("inflammation" not in row["constant_training_feature_names"] for row in rows))
        self.assertIn(d.VARIANTS[0],rows[1]["constant_training_feature_names"])
        arrays["scaler_constant"][0]=True
        with self.assertRaises(AssertionError):d.capture_cell(arrays,identity,self.frozen)

    def test_same_frozen_real_helper_on_synthetic_uniform_intercept_optimum(self):
        identity=dict(zip(d.IDENTITY_FIELDS,[d.SCOPES[0],d.PROTOCOLS[0],6,0,d.VARIANTS[0]]))
        rows=d.capture_cell(synthetic_arrays(),identity,self.frozen)
        self.assertEqual(len(rows),2)
        for row in rows:
            self.assertIs(row["optimizer_success"],True)
            self.assertEqual(row["iterations"],0)
            self.assertLess(row["gradient_inf"],1e-12)
            self.assertAlmostEqual(row["objective"],np.log(9),places=12)
            self.assertTrue(row["coefficient_agreement"])

    def test_scaler_tampering_rejected_before_any_independent_fit(self):
        arrays=synthetic_arrays();arrays["scaler_mean"][0]=1
        solver=Mock(side_effect=fake_solver)
        module=types.SimpleNamespace(weighted_scaler=self.frozen.weighted_scaler,independent_softmax_fit=solver)
        with self.assertRaises(AssertionError):d.capture_cell(arrays,{},module)
        solver.assert_not_called()

    def test_coefficient_disagreement_is_not_hidden_by_true_success_flag(self):
        arrays=synthetic_arrays();arrays["baseline_coefficients"][0,0]=.01
        def successful(*args):return {**fake_solver(*args),"optimizer_success":True}
        module=types.SimpleNamespace(weighted_scaler=self.frozen.weighted_scaler,independent_softmax_fit=successful)
        identity=dict(zip(d.IDENTITY_FIELDS,[d.SCOPES[0],d.PROTOCOLS[0],6,0,d.VARIANTS[0]]))
        rows=d.capture_cell(arrays,identity,module)
        self.assertTrue(rows[0]["optimizer_success"]);self.assertFalse(rows[0]["coefficient_agreement"])
        self.assertFalse(rows[0]["frozen_numeric_criteria_pass"])
        self.assertTrue(rows[1]["frozen_numeric_criteria_pass"])

    def test_frozen_helper_failure_never_triggers_retry_or_new_optimizer(self):
        failure=Mock(side_effect=ValueError("synthetic helper failure"))
        module=types.SimpleNamespace(weighted_scaler=self.frozen.weighted_scaler,independent_softmax_fit=failure)
        identity=dict(zip(d.IDENTITY_FIELDS,[d.SCOPES[0],d.PROTOCOLS[0],6,0,d.VARIANTS[0]]))
        with self.assertRaises(ValueError):d.capture_cell(synthetic_arrays(),identity,module)
        self.assertEqual(failure.call_count,1)


class AuthenticationTests(unittest.TestCase):
    def test_complete80_cell_flag_capture_with_original_false_flags_and_no_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            bundle=synthetic_bundle(Path(directory))
            frozen=d.load_frozen_verifier(ROOT/d.VERIFIER_RELATIVE,d.FROZEN_VERIFIER_SHA256)
            solver=Mock(side_effect=fake_solver)
            module=types.SimpleNamespace(weighted_scaler=frozen.weighted_scaler,independent_softmax_fit=solver)
            with patch.object(d,"load_frozen_verifier",return_value=module):result=run_bundle(bundle)
            self.assertEqual(solver.call_count,160);self.assertEqual(len(result["rows"]),160)
            self.assertEqual(result["summary"]["optimizer_success_false"],160)
            self.assertEqual(result["summary"]["frozen_numeric_criteria_pass"],160)
            self.assertEqual(result["existing_numerical_verification_status"],"unchanged")
            self.assertNotIn("private_line",json.dumps(result))
            self.assertNotIn("coefficients",result["rows"][0])

    def test_npz_mutation_rejected_before_loading_frozen_solver(self):
        with tempfile.TemporaryDirectory() as directory:
            bundle=synthetic_bundle(Path(directory))
            path=bundle["directory"]/bundle["manifest"][-1]["fit_npz"];path.write_bytes(path.read_bytes()+b"mutation")
            with patch.object(d,"load_frozen_verifier") as loader:
                with self.assertRaises(d.DiagnosticError):run_bundle(bundle)
                loader.assert_not_called()

    def test_manifest_mutation_and_completion_pin_mutation_fail_before_solver(self):
        with tempfile.TemporaryDirectory() as directory:
            bundle=synthetic_bundle(Path(directory))
            with patch.object(d,"load_frozen_verifier") as loader:
                changed={**bundle,"complete_sha":"0"*64}
                with self.assertRaises(d.DiagnosticError):run_bundle(changed)
                path=bundle["directory"]/"fit-manifest.json";path.write_text(json.dumps(bundle["manifest"][:-1]))
                with self.assertRaises(d.DiagnosticError):run_bundle(bundle)
                loader.assert_not_called()

    def test_missing_duplicate_failed_or_unknown_cells_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest=synthetic_bundle(Path(directory))["manifest"]
            wrong=[manifest[:-1],manifest+[manifest[0]],manifest[1:]+[manifest[1]]]
            failed=copy.deepcopy(manifest);failed[0]["status"]="optimizer_failed";wrong.append(failed)
            unknown=copy.deepcopy(manifest);unknown[0]["variant"]="posthoc_new_program";wrong.append(unknown)
            for value in wrong:
                with self.assertRaises(d.DiagnosticError):d.check_manifest(value)

    def test_verifier_authentication_and_alias_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);path=root/"frozen.py";path.write_bytes((ROOT/d.VERIFIER_RELATIVE).read_bytes())
            path.write_bytes(path.read_bytes()+b"\n# mutation")
            with self.assertRaises(d.DiagnosticError):d.load_frozen_verifier(path,d.FROZEN_VERIFIER_SHA256)
            with self.assertRaises(d.DiagnosticError):d.load_frozen_verifier(path,d.digest_bytes(path.read_bytes()))
            alias=root/"alias";alias.symlink_to(path)
            for value in [alias,root/"alias/../output",root/"../output"]:
                with self.assertRaises(d.DiagnosticError):d.safe_path(value)

    def test_runtime_and_thread_environment_gate(self):
        with tempfile.TemporaryDirectory() as directory:
            bundle=synthetic_bundle(Path(directory));path=bundle["plan_path"]
            original=json.loads(path.read_text())
            for key,value in [("runtime",{"python":"wrong"}),("thread_environment",{"DIAGNOSTIC_REQUIRED_TEST_THREAD":"1"})]:
                plan={**original,key:value};path.write_text(json.dumps(plan));changed={**bundle,"plan_sha":d.digest_bytes(path.read_bytes())}
                with patch.object(d,"load_frozen_verifier") as loader:
                    with self.assertRaises(d.DiagnosticError):run_bundle(changed)
                    loader.assert_not_called()

    def test_CLI_never_overwrites_an_existing_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);reports=root/"reports";reports.mkdir();output=reports/"existing.json";output.write_text('original')
            argv=['audit','--root',str(root),'--plan',str(root/'none'),'--plan-sha256','a'*64,
                '--execution-complete-sha256','b'*64,'--verifier-sha256',d.FROZEN_VERIFIER_SHA256,'--output',str(output)]
            with patch.object(sys,'argv',argv),patch.object(d,'run_diagnostic') as run:
                with self.assertRaises(d.DiagnosticError):d.main()
                run.assert_not_called()
            self.assertEqual(output.read_text(),'original')


if __name__=="__main__":unittest.main()
