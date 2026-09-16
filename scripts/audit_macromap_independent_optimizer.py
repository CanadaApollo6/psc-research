#!/usr/bin/env python3
"""Capture raw independent-optimizer flags without changing frozen verification.

This new post-fit diagnostic reruns only the SAME hash-pinned independent
softmax helper and scaler. It does not fit primary models, change thresholds,
choose optimizers, recompute endpoints, or replace the frozen verification.
Reports contain aggregate fit identities, never training line IDs or arrays.
Run under the thread environment specified in the frozen plan.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import os
import platform
import re
import types
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy

ROOT=Path(__file__).resolve().parents[1]
VERIFIER_RELATIVE="scripts/verify_macromap_program_context.py"
FROZEN_VERIFIER_SHA256="31b69fca6bcbe39e9204d3816203b39c5cb96762ef2e56ec9562dede77c71a7a"
SCOPES=("mapped_primary","all_source_lines_sensitivity")
PROTOCOLS=("Mod_SmartSeq2","NEB")
VARIANTS=("ets2_g1_dn","ets2_g1_dn_without_comparators")
TIMES=(6,24)
IDENTITY_FIELDS=("identity_scope","protocol","time","fold","variant")


class DiagnosticError(ValueError):
    pass


def require(condition: bool,message: str) -> None:
    if not condition:raise DiagnosticError(message)


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def safe_path(path: Path) -> Path:
    require(".." not in path.parts,"Literal parent traversal is forbidden")
    absolute=Path(os.path.abspath(path))
    require(not any(part.is_symlink() for part in [absolute,*absolute.parents]),"Aliased input/output path")
    return absolute


def contained(directory: Path,relative: str) -> Path:
    name=Path(relative)
    require(not name.is_absolute() and bool(name.parts) and ".." not in name.parts,"Invalid artifact path")
    path=safe_path(directory/name)
    require(path.is_relative_to(safe_path(directory)),"Artifact escapes completed bundle")
    return path


def read_pinned(path: Path,expected: str) -> bytes:
    require(isinstance(expected,str) and bool(re.fullmatch(r"[a-f0-9]{64}",expected)),"Exact SHA-256 required")
    body=safe_path(path).read_bytes()
    require(digest_bytes(body)==expected,"Pinned artifact changed")
    return body


def load_frozen_verifier(path: Path,expected: str):
    require(expected==FROZEN_VERIFIER_SHA256,"This diagnostic supports only the frozen verifier")
    body=read_pinned(path,expected)
    # Compile the very bytes that were authenticated, rather than reopening a
    # potentially changed pathname. This executes trusted, externally pinned
    # project code, not code from NPZ/JSON inputs. No main() call occurs.
    module=types.ModuleType("_frozen_macromap_optimizer_diagnostic")
    module.__file__=str(safe_path(path))
    exec(compile(body,module.__file__,"exec"),module.__dict__)
    return module


def runtime_versions() -> dict:
    return {"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__}


def check_manifest(manifest: list[dict]) -> None:
    expected={(s,p,t,f,g) for s in SCOPES for p in PROTOCOLS for t in TIMES for f in range(5) for g in VARIANTS}
    keys=[tuple(row[key] for key in IDENTITY_FIELDS) for row in manifest]
    require(len(keys)==80 and len(set(keys))==80 and set(keys)==expected,"Require all80 unique frozen fit cells")
    require(all(row["status"]=="available" for row in manifest),"Only the80 available primary fit cells may be audited")


def capture_cell(arrays: dict,identity: dict,frozen) -> list[dict]:
    """Use the exact scaling/refit recipe from frozen verify_fit_arrays()."""
    x=np.asarray(arrays["train_x_raw"],float)
    labels=arrays["train_class_index"]
    lines=arrays["train_line"].tolist()
    require(x.shape==(len(labels),5) and len(lines)==len(labels),"Unexpected training dimensions")
    require(np.issubdtype(labels.dtype,np.integer) and set(labels)==set(range(9)),"Invalid frozen class axis")
    mean,scale,_=frozen.weighted_scaler(x,lines)
    # Identical saved-scaler tolerances and independent transformation.
    np.testing.assert_allclose(arrays["scaler_mean"],mean,rtol=1e-10,atol=1e-12)
    np.testing.assert_allclose(arrays["scaler_scale"],scale,rtol=1e-10,atol=1e-12)
    # weighted_scaler returns a unit fallback but not its raw constant mask.
    # Reconstruct its exact equal-line population-variance calculation. A
    # returned scale of one alone does NOT imply a constant training feature.
    groups={line:np.flatnonzero(np.asarray(lines)==line) for line in dict.fromkeys(lines)}
    variance=sum(((x[group]-mean)**2).mean(axis=0) for group in groups.values())/len(groups)
    constant=np.sqrt(variance)<=1e-12
    saved_constant=np.asarray(arrays["scaler_constant"])
    require(saved_constant.shape==(5,) and saved_constant.dtype==np.dtype(bool),"Missing saved constant-feature flags")
    np.testing.assert_array_equal(saved_constant,constant)
    feature_names=["inflammation","interferon_gamma","oxidative_stress","apoptotic_signaling",identity["variant"]]
    np.testing.assert_array_equal(arrays["predictor_names"],feature_names)
    standardized=(x-mean)/scale
    rows=[]
    for model,p in [("baseline",4),("extended",5)]:
        # No options, alternate solver, retry or tolerance override. The frozen
        # function returns False success flags when its own numeric gate passes.
        fit=frozen.independent_softmax_fit(standardized[:,:p],labels,lines)
        coefficients=np.asarray(fit["coefficients"],float)
        saved=np.asarray(arrays[model+"_coefficients"],float)
        require(coefficients.shape==saved.shape==(p+1,9) and np.isfinite(saved).all(),"Invalid coefficient shape/value")
        require(isinstance(fit["optimizer_success"],(bool,np.bool_)),"Missing original solver-success flag")
        objective=float(fit["objective"]);gradient=float(fit["gradient_inf"])
        finite=bool(np.isfinite(objective) and np.isfinite(gradient) and np.isfinite(coefficients).all())
        require(finite,"Frozen helper returned nonfinite diagnostics")
        agreement=bool(np.allclose(saved,coefficients,rtol=2e-5,atol=2e-5,equal_nan=False))
        # This repeats, rather than tightens or relaxes, the existing gate.
        numeric_pass=bool(gradient<=1e-7 and agreement)
        raw_success=bool(fit["optimizer_success"])
        rows.append({**{key:identity[key] for key in IDENTITY_FIELDS},"model":model,
            "optimizer_success":raw_success,"iterations":int(fit["iterations"]),
            "gradient_inf":gradient,"objective":objective,
            "coefficient_max_absolute_error":float(np.max(np.abs(saved-coefficients))),
            "coefficient_agreement":agreement,"frozen_numeric_criteria_pass":numeric_pass,
            "solver_flag_status":"reported_success" if raw_success else "reported_unsuccessful",
            "training_rows":len(x),"training_feature_count":p,
            "constant_training_feature_count":int(np.sum(constant[:p])),
            "constant_training_feature_names":[feature_names[j] for j in range(p) if constant[j]],
            "saved_constant_flags_verified":True})
    return rows


def run_diagnostic(plan_path: Path,plan_sha256: str,completion_sha256: str,verifier_sha256: str,*,root: Path=ROOT) -> dict:
    root=safe_path(root)
    plan=json.loads(read_pinned(plan_path,plan_sha256))
    require(plan.get("status")=="frozen" and plan.get("execute_real_expression") is True and
        plan.get("analysis_id")=="macromap-fixed-program-context-v1","Not the separately authorized frozen analysis")
    require(plan["runtime"]==runtime_versions(),"Numerical runtime differs from freeze")
    require(all(os.environ.get(key)==value for key,value in plan["thread_environment"].items()),"Use the frozen thread/hash-seed environment before Python starts")
    code_pin=next((row for row in plan["code_and_protocol_artifacts"] if row["path"]==VERIFIER_RELATIVE),None)
    require(code_pin is not None and code_pin["sha256"]==verifier_sha256==FROZEN_VERIFIER_SHA256,"Verifier pin differs from root freeze")
    directory=contained(root,plan["output_directory"])
    require(directory.is_relative_to(root/"work/macromap-design"),"Expected the frozen private primary bundle")
    receipt=json.loads(read_pinned(directory/"execution-complete.json",completion_sha256))
    require(receipt["stage"]=="execution_complete_all_prespecified_statuses_retained" and receipt["plan_sha256"]==plan_sha256,"Foreign/incomplete primary receipt")
    require(receipt["runtime"]==plan["runtime"] and receipt["script_sha256"]==plan["script_sha256"] and
        receipt["runner_script_sha256"]==plan["runner_script_sha256"] and
        receipt["prepared_summary_sha256"]==plan["prepared_summary_sha256"],"Primary receipt provenance differs from freeze")
    pins=receipt["artifact_sha256"]
    manifest_sha=pins["fit-manifest.json"]
    manifest=json.loads(read_pinned(directory/"fit-manifest.json",manifest_sha))
    check_manifest(manifest)
    authenticated=[]
    # Authenticate every NPZ before invoking any numerical helper. The exact
    # authenticated bytes are then loaded with pickle disabled.
    for spec in sorted(manifest,key=lambda row:tuple(row[key] for key in IDENTITY_FIELDS)):
        expected=f"{spec['identity_scope']}/{spec['protocol']}-t{spec['time']}-fold{spec['fold']}/{spec['variant']}/fit-arrays.npz"
        require(spec["fit_npz"]==expected and spec["fit_npz_sha256"]==pins[expected],"NPZ identity/manifest/completion pin differs")
        authenticated.append((spec,read_pinned(contained(directory,expected),pins[expected])))
    frozen=load_frozen_verifier(root/VERIFIER_RELATIVE,verifier_sha256)
    rows=[]
    for spec,body in authenticated:
        with np.load(io.BytesIO(body),allow_pickle=False) as state:arrays=dict(state)
        rows.extend(capture_cell(arrays,spec,frozen))
    require(len(rows)==160,"Incomplete independent optimizer capture")
    numeric_pass=sum(row["frozen_numeric_criteria_pass"] for row in rows)
    success_counts=Counter(row["optimizer_success"] for row in rows)
    return {"schema":"macromap-independent-optimizer-flag-capture-v1",
        "purpose":"post_fit_flag_capture_not_replacement_verification",
        "status":"captured_frozen_numeric_criteria_pass" if numeric_pass==160 else "captured_numeric_discrepancy",
        "existing_numerical_verification_status":"unchanged",
        "created_at_utc":datetime.now(timezone.utc).isoformat(),
        "provenance":{"plan_sha256":plan_sha256,"primary_completion_sha256":completion_sha256,
            "fit_manifest_sha256":manifest_sha,"frozen_verifier_sha256":verifier_sha256,
            "diagnostic_wrapper_sha256":digest_bytes(Path(__file__).read_bytes()),
            "runtime":runtime_versions(),"thread_environment":dict(plan["thread_environment"]),
            "fit_npz_pins":[{**{key:spec[key] for key in IDENTITY_FIELDS},"sha256":spec["fit_npz_sha256"]} for spec,_ in authenticated]},
        "metric_definitions":{"optimizer_success":"Raw success flag from the same frozen independent helper rerun, not the primary optimizer.",
            "gradient_inf":"Returned independent-fit gradient norm, not the saved primary-model gradient.",
            "objective":"Returned independent weighted mean cross-entropy plus the frozen ridge penalty.",
            "coefficient_max_absolute_error":"Maximum difference from saved primary coefficients, including intercepts."},
        "unchanged_criteria":{"independent_gradient_limit":1e-7,"coefficient_rtol":2e-5,"coefficient_atol":2e-5,
            "scaler_rtol":1e-10,"scaler_atol":1e-12,"constant_training_sd_threshold":1e-12,
            "method":"frozen independent_softmax_fit defaults",
            "raw_solver_success_is_separate_from_numeric_acceptance":True},
        "summary":{"fit_cells":80,"models":160,"optimizer_success_true":success_counts[True],
            "optimizer_success_false":success_counts[False],"frozen_numeric_criteria_pass":numeric_pass,
            "models_with_constant_training_features":sum(row["constant_training_feature_count"]>0 for row in rows),
            "constant_training_feature_name_counts":dict(Counter(name for row in rows for name in row["constant_training_feature_names"])),
            "all_saved_constant_flags_verified":all(row["saved_constant_flags_verified"] for row in rows),
            "max_gradient_inf":max(row["gradient_inf"] for row in rows),
            "max_coefficient_absolute_error":max(row["coefficient_max_absolute_error"] for row in rows)},
        "rows":rows,
        "limitations":["These are flags from rerunning the identical frozen independent helper, not stored original SciPy result objects.",
            "False optimizer-success flags are preserved even when frozen gradient/coefficient criteria pass.",
            "No primary refit, alternate optimizer, threshold change, endpoint recomputation or interpretation rescue is performed.",
            "The frozen helper does not return the SciPy message/status code; none is invented."]}


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=ROOT)
    parser.add_argument("--plan",type=Path,required=True)
    parser.add_argument("--plan-sha256",required=True)
    parser.add_argument("--execution-complete-sha256",required=True)
    parser.add_argument("--verifier-sha256",required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    output=safe_path(args.output)
    require(not output.exists(),"Never overwrite an existing diagnostic")
    require(output.is_relative_to(safe_path(args.root)/"reports"),"Public aggregate diagnostic belongs under reports")
    result=run_diagnostic(args.plan,args.plan_sha256,args.execution_complete_sha256,args.verifier_sha256,root=args.root)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("x",encoding="utf-8") as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write("\n")
    print(json.dumps({"status":result["status"],**result["summary"]},sort_keys=True))
    if result["summary"]["frozen_numeric_criteria_pass"]!=160:raise SystemExit(2)


if __name__=="__main__":main()
