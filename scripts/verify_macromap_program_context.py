#!/usr/bin/env python3
"""Independent MacroMap checks. No imports from primary scoring/model helpers.

Default CLI checks source-only prepared identities and folds. Reading fit-value
artifacts requires a separately frozen, externally hash-pinned execution plan.
This is stimulus context in healthy-derived macrophages, not PSC validation.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[1]
PROTOCOLS = ("Mod_SmartSeq2", "NEB")
TIMES = (6, 24)
STIMULI = ("CIL", "IFNB", "IFNG", "IL4", "LIL10", "MBP", "P3C", "R848", "sLPS")
COMPARATORS = ("inflammation", "interferon_gamma", "oxidative_stress", "apoptotic_signaling")
ORIGINAL_PROGRAMS = ("ets2_g1_dn", "ets2_g1_up", "ets2_g2_dn", "chr21_dn") + COMPARATORS
NONOVERLAP = "ets2_g1_dn_without_comparators"


class VerificationError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def pinned(path: Path, expected: str) -> Path:
    require(isinstance(expected, str) and bool(re.fullmatch(r"[a-f0-9]{64}", expected)), "Missing exact SHA-256")
    require(path.is_file() and digest(path) == expected, "Pinned artifact changed or missing")
    return path


def read_rows(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = csv.reader(stream, delimiter=delimiter, strict=True)
        header = next(rows, [])
        require(bool(header) and len(set(header)) == len(header), "Missing or duplicate physical CSV header")
        result = []
        for row in rows:
            require(len(row) == len(header), "Ragged CSV row")
            result.append(dict(zip(header, row)))
    return result


def exact_keys(rows: list[dict]) -> Counter:
    return Counter((r["protocol"], r["run"], r["line"], int(r["time"]), r["stimulus"]) for r in rows)


def reconstruct_folds(samples: list[dict], seed: int = 20260916, k: int = 5) -> dict:
    line_groups = defaultdict(set); runs = defaultdict(set)
    for row in samples:
        line_groups[row["line"]].add((row["protocol"], row["run"]))
        runs[row["protocol"]].add(row["run"])
    require(all(len(group) == 1 for group in line_groups.values()), "Line crosses run/protocol")
    assignments = {}
    for protocol in sorted(runs):
        require(len(runs[protocol]) >= k, "Insufficient runs for frozen folds")
        ordered = sorted(runs[protocol], key=lambda run: (hashlib.sha256(f"{seed}:{protocol}:{run}".encode()).hexdigest(), run))
        for ordinal, run in enumerate(ordered):
            assignments[protocol, run] = ordinal % k
    return assignments


def reconstruct_pairs(samples: list[dict], include_nomatch: bool = False) -> list[dict]:
    by_line = defaultdict(dict)
    for row in samples:
        require(row["condition"] not in by_line[row["line"]], "Duplicate source culture")
        by_line[row["line"]][row["condition"]] = row
    result = []
    for line, conditions in sorted(by_line.items()):
        info = next(iter(conditions.values()))
        for stimulus in STIMULI + ("PIC",):
            for time in TIMES:
                treated = conditions.get(f"{stimulus}_{time}"); control = conditions.get(f"Ctrl_{time}")
                status = ("excluded_NOMATCH_identity" if info["HipsciID"] == "NOMATCH" and not include_nomatch
                    else "unavailable_mock_identity" if stimulus == "PIC"
                    else "missing_stimulus_and_control" if treated is None and control is None
                    else "missing_stimulus" if treated is None
                    else "missing_time_matched_control" if control is None else "available")
                result.append({"protocol": info["protocol"], "run": info["run"], "line": line,
                    "time": time, "stimulus": stimulus, "status": status,
                    "control_id": control["sample_id"] if control else "",
                    "treatment_id": treated["sample_id"] if treated else "",
                    "control_index": int(control["sample_index"]) if control else "",
                    "treatment_index": int(treated["sample_index"]) if treated else ""})
    return result


def reconstruct_prediction_cohort(pairs: list[dict]) -> list[dict]:
    valid = [row for row in pairs if row["status"] == "available"]
    times = defaultdict(set)
    for row in valid: times[row["line"]].add(int(row["time"]))
    return [row for row in valid if times[row["line"]] == set(TIMES)]


def verify_preparation(directory: Path, expected_sha256: str, root: Path = ROOT) -> dict:
    """Metadata/annotation only. Counts are hashed as opaque bytes, never parsed."""
    summary = json.loads(pinned(directory/"source-only-readiness.json", expected_sha256).read_text())
    require(summary["ready_for_inference"] is False and summary["real_program_scores_computed"] is False
        and summary["real_effects_or_predictions_computed"] is False, "Not a source-only preparation")
    for name, sha in summary["artifact_sha256"].items():
        require(Path(name).name == name and name not in {".", ".."}, "Unsafe prepared artifact name")
        pinned(directory/name, sha)
    for name, pin in summary["source_pins"].items():
        path = root/name
        require(path.resolve().is_relative_to(root.resolve()), "Source pin escapes repository")
        pinned(path, pin["sha256"])
        require(path.stat().st_size == pin["bytes"], "Source byte count differs")
    samples = read_rows(directory/"sample-map.csv")
    metadata = read_rows(root/"data/raw/macromap-qualification/worker-methods/sample-metadata.csv")
    original = {r["SampleID"]+"_"+r["RunID"]: r for r in metadata}
    header = (root/"data/raw/macromap-qualification/MacroMap_raw_expression.header.txt").read_text().rstrip("\r\n").split("\t")[6:]
    require(len(original) == len(metadata) == len(samples) and set(original) == set(header), "Source sample identities differ")
    require([r["sample_id"] for r in samples] == header, "Prepared sample order differs from source header")
    grouped = defaultdict(list)
    for index, row in enumerate(samples):
        source = original[row["sample_id"]]
        require(all(row[key] == value for key, value in source.items()), "Original metadata value changed")
        suffix = "_" + source["Stimulus_Hours"]
        require(source["SampleID"].endswith(suffix), "Invalid source condition suffix")
        line = source["SampleID"][:-len(suffix)]
        require(row["line"] == line and row["condition"] == source["Stimulus_Hours"] and
            row["run"] == source["RunID"] and row["protocol"] == source["Library_prep"] and
            int(row["sample_index"]) == index and row["identity_mapped"] == str(source["HipsciID"] != "NOMATCH"), "Derived sample identity changed")
        grouped[line].append(row)
    for rows in grouped.values():
        require(all(len({row[field] for row in rows}) == 1 for field in metadata[0] if field not in {"SampleID", "Stimulus_Hours"}), "Within-line metadata changes")
    mapped_identities = defaultdict(set)
    for row in samples:
        if row["HipsciID"] != "NOMATCH": mapped_identities[row["HipsciID"]].add(row["line"])
    require(all(len(value) == 1 for value in mapped_identities.values()), "HipSci identity has multiple source lines")
    expected_pairs = reconstruct_pairs(samples)
    actual_pairs = read_rows(directory/"pair-registry.csv")
    require(exact_keys(expected_pairs) == exact_keys(actual_pairs), "Prepared pair registry differs")
    actual_by_key = {(r["line"], int(r["time"]), r["stimulus"]): r for r in actual_pairs}
    for row in expected_pairs:
        other = actual_by_key[row["line"], row["time"], row["stimulus"]]
        require(all(str(value) == other[key] for key, value in row.items()), "Pair status/index/source label differs")
    predictive = reconstruct_prediction_cohort(expected_pairs)
    require(exact_keys(predictive) == exact_keys(read_rows(directory/"predictive-pairs.csv")), "Prediction cohort differs")
    folds = reconstruct_folds(samples)
    saved_folds = read_rows(directory/"run-folds.csv")
    require(len(saved_folds) == len(folds) and {(r["protocol"], r["run"]): int(r["fold"]) for r in saved_folds} == folds, "Whole-run fold assignment differs")
    all_gates = []
    for protocol in PROTOCOLS:
        for time in TIMES:
            for fold in range(5):
                train = [r for r in predictive if r["protocol"] == protocol and r["time"] == time and folds[protocol,r["run"]] != fold]
                test = [r for r in predictive if r["protocol"] == protocol and r["time"] == time and folds[protocol,r["run"]] == fold]
                require(not ({r["line"] for r in train} & {r["line"] for r in test}), "Fold line leakage")
                require(bool(test), "Empty test fold")
                require(all(len({r["line"] for r in train if r["stimulus"] == s}) >= 3 and len({r["run"] for r in train if r["stimulus"] == s}) >= 3 for s in STIMULI), "Insufficient training class lines/runs")
                all_gates.append(True)
    features = read_rows(directory/"features.csv")
    require([int(r["feature_index"]) for r in features] == list(range(len(features))), "Feature indices differ")
    feature_ids = [r["feature_id"] for r in features]
    require(len(set(feature_ids)) == len(feature_ids), "Duplicate original feature ID")
    universe = read_rows(root/"data/raw/macromap-qualification/raw-gene-universe.tsv", "\t")
    require([r["Geneid"] for r in universe] == feature_ids, "Prepared feature order differs from released universe")
    annotation = {}
    for line in (root/"data/raw/macromap-qualification/worker-methods/github_gencode_v27_annotation.body").read_text().splitlines():
        fields = line.split()
        require(len(fields) == 6 and fields[4] not in annotation, "Invalid original annotation")
        annotation[fields[4]] = fields[5]
    require(set(feature_ids) <= set(annotation), "Source feature lacks original annotation")
    require([annotation[value] for value in feature_ids] == [r["feature_name"] for r in features], "Prepared gene label differs from original annotation")
    absent = set(annotation)-set(feature_ids)
    require(len(absent) == 45 and all(value.endswith("_PAR_Y") for value in absent), "Reference-only feature boundary differs")
    stable = [re.sub(r"\.\d+$", "", value) for value in feature_ids]
    require(len(set(stable)) == len(stable), "Stable-ID collision")
    require(stable == [r["stable_feature_id"] for r in features], "Stable-ID normalization differs")
    membership = read_rows(root/"data/derived/ets2-program-membership.csv")
    prepared_map = read_rows(directory/"program-mapping.csv")
    by_id = {value: i for i,value in enumerate(stable)}; by_name = defaultdict(list)
    for i,row in enumerate(features): by_name[row["feature_name"]].append(i)
    mapped = {}; source_counts = {}
    for program in ORIGINAL_PROGRAMS:
        copies = [[r["source_value"] for r in sorted((r for r in membership if r["dataset_key"] == modality and r["program_id"] == program), key=lambda r:int(r["source_index"]))] for modality in ["sc", "sn"]]
        require(copies[0] == copies[1] and len(copies[0]) == len(set(copies[0])), "Historical source membership copies differ")
        selected = set(); intended = 0
        for value in copies[0]:
            match = re.fullmatch(r"(ENSG[0-9]+)(?:\.[0-9]+)?", value)
            key = match.group(1) if match else value
            if key in {"ENSG00000157557", "ETS2"}: continue
            intended += 1
            indices = ([by_id[key]] if key in by_id else []) if match else by_name.get(key, [])
            if len(indices) == 1: selected.add(indices[0])
        mapped[program] = selected; source_counts[program] = intended
        saved = {int(r["feature_index"]) for r in prepared_map if r["program_id"] == program and r["mapping_status"] == "mapped"}
        require(saved == selected, "Independently mapped program members differ")
    comparator_union = set().union(*(mapped[p] for p in COMPARATORS))
    no_overlap = mapped["ets2_g1_dn"] - comparator_union
    saved = {int(r["feature_index"]) for r in prepared_map if r["program_id"] == NONOVERLAP and r["mapping_status"] == "mapped_nonoverlap"}
    require(saved == no_overlap, "Nonoverlap sensitivity membership differs")
    inventory = {r["program_id"]:r for r in read_rows(directory/"mapping-inventory.csv")}
    for program in ORIGINAL_PROGRAMS:
        gate = len(mapped[program]) >= 10 and len(mapped[program])/source_counts[program] >= .8
        require(inventory[program]["mapping_gate_pass"] == str(gate) and int(inventory[program]["intended_unique_source_genes"]) == source_counts[program] and int(inventory[program]["mapped_unique_genes"]) == len(mapped[program]), "Program gate/count differs")
    require(inventory[NONOVERLAP]["mapping_gate_pass"] == str(len(no_overlap)>=10 and inventory["ets2_g1_dn"]["mapping_gate_pass"] == "True"), "Parent nonoverlap gate differs")
    return {"status":"source_only_verified", "preparation_sha256":expected_sha256, "samples":len(samples), "features":len(features),
        "source_lines":len(grouped), "mapped_lines":len(mapped_identities), "predictive_lines":len({r["line"] for r in predictive}),
        "runs":len(folds), "fold_gates":len(all_gates), "programs":len(mapped)+1,
        "overlap_removed":len(mapped["ets2_g1_dn"])-len(no_overlap), "effects_read":False, "real_scores_or_fits_computed":False}


def training_reference(log_values: np.ndarray, pairs: list[dict], folds: dict, protocol: str, time: int, fold: int) -> tuple[np.ndarray,list[int]]:
    chosen = {}
    runs = set()
    for row in pairs:
        if row["status"] != "available" or row["protocol"] != protocol or int(row["time"]) != time or folds[protocol,row["run"]] == fold: continue
        index = int(row["control_index"])
        require(row["line"] not in chosen or chosen[row["line"]] == index, "Line has multiple reference controls")
        chosen[row["line"]] = index; runs.add(row["run"])
    require(len(chosen) >= 3 and len(runs) >= 3, "Insufficient reference controls/runs")
    indices = [chosen[line] for line in sorted(chosen)]
    require(len(set(indices)) == len(indices), "Reference control reused across lines")
    require(min(indices)>=0 and max(indices)<log_values.shape[1], "Bad control index")
    values = np.asarray(log_values[:,indices], dtype=float)
    require(np.isfinite(values).all(), "Invalid training controls")
    # Explicit per-line accumulation instead of a primary scoring helper.
    return sum(values[:,i] for i in range(len(indices)))/len(indices), indices


def independent_matches(reference: np.ndarray, ids: list[str], target: list[int], excluded: list[int], *, program: str, protocol: str, time: int, fold: int, bins: int = 20, draws: int = 1000) -> dict:
    require(len(reference)==len(ids) and len(set(ids))==len(ids) and np.isfinite(reference).all(), "Bad matching universe")
    target = sorted(set(target)); exclusions = set(excluded)
    require(set(target)<=exclusions and target and bins>0 and draws>0, "Bad target/matching settings")
    order = sorted(range(len(ids)), key=lambda i:(float(reference[i]), re.sub(r"\.\d+$", "", ids[i]), ids[i]))
    groups = [list(map(int,g)) for g in np.array_split(order,bins)]
    counts = [len(set(group)&set(target)) for group in groups]
    pools = [[i for i in group if i not in exclusions] for group in groups]
    underflow = [i for i in range(bins) if counts[i]>len(pools[i])]
    text = f"20260912:macromap:protocol={protocol}:time={time}:fold={fold}:{program}"
    seed = int(hashlib.sha256(text.encode()).hexdigest()[:16],16)
    feature_bin = np.empty(len(ids),dtype=int)
    for i,group in enumerate(groups): feature_bin[group]=i
    metadata = {"feature_bin":feature_bin,"target_bin_counts":{i:n for i,n in enumerate(counts) if n},
        "pool_bin_counts":{i:len(pool) for i,pool in enumerate(pools)},"seed":seed,"underflow_bins":underflow}
    if underflow: return {**metadata,"available":False,"reason":"insufficient_control_pool"}
    rng = np.random.Generator(np.random.PCG64(seed)); chosen = []
    for _ in range(draws):
        row=[]
        for need,pool in zip(counts,pools):
            if need: row.extend(sorted(map(int,rng.choice(pool,size=need,replace=False))))
        chosen.append(row)
    selection = np.asarray(chosen,dtype='<i8')
    background = np.zeros(len(ids))
    for row in chosen:
        for gene in row: background[gene] += 1/(draws*len(target))
    return {**metadata,"available":True,"reason":"available","selections":selection,"background_weights":background,
        "selection_sha256":hashlib.sha256(selection.tobytes()).hexdigest()}


def weighted_scaler(x: np.ndarray, lines: list[str]) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    x=np.asarray(x,dtype=float)
    require(x.ndim==2 and len(x)==len(lines) and len(lines)>0 and np.isfinite(x).all(), "Bad training predictors")
    groups=defaultdict(list)
    for i,line in enumerate(lines): groups[line].append(i)
    weights=np.asarray([1/(len(groups)*len(groups[line])) for line in lines])
    mean=sum(x[group].mean(axis=0) for group in groups.values())/len(groups)
    variance=sum(((x[group]-mean)**2).mean(axis=0) for group in groups.values())/len(groups)
    sd=np.sqrt(variance)
    return mean,np.where(sd<=1e-12,1.,sd),weights


def logits_logp(coefficients: np.ndarray, x: np.ndarray) -> np.ndarray:
    logits=np.column_stack([np.ones(len(x)),x]) @ coefficients
    return logits-logsumexp(logits,axis=1)[:,None]


def softmax_calculus(coefficients: np.ndarray, x: np.ndarray, labels: np.ndarray, weights: np.ndarray, penalty: float=.1) -> tuple[float,np.ndarray,np.ndarray]:
    """Independent rowwise likelihood, analytic score and full Hessian."""
    x=np.asarray(x,float); coefficients=np.asarray(coefficients,float); labels=np.asarray(labels)
    n,p=x.shape; k=coefficients.shape[1]
    require(coefficients.shape==(p+1,k) and labels.shape==(n,) and weights.shape==(n,), "Softmax shape mismatch")
    require(np.issubdtype(labels.dtype,np.integer) and set(labels)<=set(range(k)) and np.isfinite(x).all() and np.isfinite(coefficients).all(), "Bad softmax inputs")
    require(np.isfinite(weights).all() and np.all(weights>0) and np.isclose(sum(weights),1) and penalty>0, "Mean-loss weights/penalty invalid")
    gradient=np.zeros_like(coefficients); hessian=np.zeros(((p+1)*k,(p+1)*k)); objective=0.
    for vector,y,w in zip(x,labels,weights):
        design=np.r_[1.,vector]; eta=design @ coefficients
        normalizer=logsumexp(eta); probability=np.exp(eta-normalizer)
        objective+=w*(normalizer-eta[y])
        error=probability.copy(); error[y]-=1
        gradient+=w*np.outer(design,error)
        class_cov=np.diag(probability)-np.outer(probability,probability)
        hessian+=w*np.kron(np.outer(design,design),class_cov)
    objective+=penalty/2*np.square(coefficients[1:]).sum()
    gradient[1:]+=penalty*coefficients[1:]
    hessian[k:,k:]+=penalty*np.eye(p*k)
    return float(objective),gradient,hessian


def independent_softmax_fit(x: np.ndarray, labels: np.ndarray, lines: list[str], k: int=9, penalty: float=.1) -> dict:
    """Trust-exact fit in an orthonormal sum-zero class basis, not primary LBFGS."""
    x=np.asarray(x,float); labels=np.asarray(labels,dtype=int)
    require(set(labels)==set(range(k)), "Every frozen class is required")
    _,_,weights=weighted_scaler(x,lines)
    basis=np.zeros((k,k-1))
    for j in range(k-1):
        basis[:j+1,j]=1/np.sqrt((j+1)*(j+2)); basis[j+1,j]=-(j+1)/np.sqrt((j+1)*(j+2))
    projection=np.kron(np.eye(x.shape[1]+1),basis)
    shape=(x.shape[1]+1,k-1)
    def terms(theta):
        coef=theta.reshape(shape) @ basis.T
        value,gradient,hessian=softmax_calculus(coef,x,labels,weights,penalty)
        return value,gradient@basis,projection.T@hessian@projection
    def fun(theta):
        value,gradient,_=terms(theta);return value,gradient.ravel()
    fit=minimize(fun,np.zeros(np.prod(shape)),jac=True,hess=lambda theta:terms(theta)[2],method='trust-exact',options={'gtol':1e-9,'maxiter':300})
    coefficients=fit.x.reshape(shape)@basis.T
    value,gradient,_=softmax_calculus(coefficients,x,labels,weights,penalty)
    require(np.isfinite(value) and np.max(np.abs(gradient))<=1e-7, "Independent optimum not qualified")
    return {'coefficients':coefficients,'objective':value,'gradient_inf':float(np.max(np.abs(gradient))),'optimizer_success':bool(fit.success),'iterations':int(fit.nit)}


def loss_records(records: list[dict], expected: list[dict]) -> tuple[dict,dict]:
    require(exact_keys(records)==exact_keys(expected), "Missing or extra expected OOF prediction")
    require(all(value==1 for value in exact_keys(records).values()), "Duplicate prediction identity")
    grouped=defaultdict(list); line_locations=defaultdict(set)
    for row in records:
        require(row['stimulus'] in STIMULI and int(row['time']) in TIMES and row['protocol'] in PROTOCOLS,'Unknown prediction stratum')
        require(all(np.isfinite(row[key]) and row[key]>=0 for key in ['baseline_loss','extended_loss']),'Unavailable prediction cannot be omitted')
        grouped[row['protocol'],row['run'],row['line'],int(row['time'])].append(row)
        line_locations[row['line']].add((row['protocol'],row['run']))
    require(all(len(locations)==1 for locations in line_locations.values()),'Predicted line crosses runs/protocols')
    line_loss={key:np.mean([[r['baseline_loss'],r['extended_loss']] for r in rows],axis=0) for key,rows in grouped.items()}
    populations={p:{key[2] for key in line_loss if key[0]==p} for p in PROTOCOLS}
    for p in PROTOCOLS:
        require(populations[p] and all({key[2] for key in line_loss if key[0]==p and key[3]==t}==populations[p] for t in TIMES),'Missing line time')
        for t in TIMES:require({r['stimulus'] for r in records if r['protocol']==p and int(r['time'])==t}==set(STIMULI),'Missing evaluation class')
    return line_loss,populations


def aggregate_losses(records: list[dict], expected: list[dict]) -> dict:
    line_loss,populations=loss_records(records,expected)
    cells={}
    for p in PROTOCOLS:
        for t in TIMES:
            losses=np.asarray([v for key,v in line_loss.items() if key[0]==p and key[3]==t])
            base,extended=losses.mean(axis=0)
            cells[p,t]={'baseline_loss':float(base),'extended_loss':float(extended),'gain':float(base-extended)}
    total=sum(len(lines) for lines in populations.values())
    per_line=[]
    for p in PROTOCOLS:
        for line in sorted(populations[p]):
            losses=np.asarray([v for key,v in line_loss.items() if key[0]==p and key[2]==line])
            per_line.append((losses[:,0]-losses[:,1]).mean())
    return {'primary_gain':float(np.mean(per_line)),'cells':cells,'protocol_weights':{p:len(populations[p])/total for p in PROTOCOLS},'n_lines':total}


def conditional_run_bootstrap(records: list[dict], expected: list[dict], draws: int=5000, seed: int=2026091601) -> dict:
    line_loss,populations=loss_records(records,expected)
    point=aggregate_losses(records,expected)
    runs={p:sorted({key[1] for key in line_loss if key[0]==p}) for p in PROTOCOLS}
    require(min(map(len,runs.values()))>=3 and draws>0,'Insufficient bootstrap runs/draws')
    rng=np.random.Generator(np.random.PCG64(seed)); values=np.empty(draws); cells={key:np.empty(draws) for key in point['cells']}
    # Expand selected runs and all their lines explicitly, rather than using the
    # primary implementation's sufficient-sum/count matrix calculation.
    for b in range(draws):
        total=0.
        for p in PROTOCOLS:
            selected=[runs[p][i] for i in rng.integers(len(runs[p]),size=len(runs[p]))]
            time_values=[]
            for t in TIMES:
                expanded=[float(loss[0]-loss[1]) for run in selected for key,loss in line_loss.items() if key[0]==p and key[1]==run and key[3]==t]
                mean=float(np.mean(expanded)); cells[p,t][b]=mean;time_values.append(mean)
            total+=point['protocol_weights'][p]*np.mean(time_values)
        values[b]=total
    return {'point':point,'interval':np.quantile(values,[.025,.975],method='linear').tolist(),
        'cell_intervals':{key:np.quantile(value,[.025,.975],method='linear').tolist() for key,value in cells.items()},
        'draws':draws,'seed':seed,'conditional_fixed_predictions_only':True}


def verify_fit_arrays(arrays: dict[str,np.ndarray], refit: bool=True) -> dict:
    """Future read-only audit of one saved fold's numerical model state.

    This cannot prove that saved training predictors came from the right matched
    references. Source/design/matching manifests must be checked separately.
    """
    required={'train_x','test_x','train_y','test_y','train_lines','test_lines','train_runs','test_runs','scaler_mean','scaler_scale','baseline_coefficients','extended_coefficients','baseline_logp','extended_logp'}
    require(required<=arrays.keys(),'Incomplete numerical fit-audit state')
    x=np.asarray(arrays['train_x'],float); test=np.asarray(arrays['test_x'],float)
    require(x.shape==(len(arrays['train_y']),5) and test.shape==(len(arrays['test_y']),5),'Wrong predictor columns')
    require(not set(arrays['train_lines'])&set(arrays['test_lines']) and not set(arrays['train_runs'])&set(arrays['test_runs']),'Held-out line/run leakage')
    require(set(arrays['train_y'])==set(range(9)), 'Training class absence')
    for label in range(9):
        members=np.flatnonzero(arrays['train_y']==label)
        require(len(set(arrays['train_lines'][members]))>=3 and len(set(arrays['train_runs'][members]))>=3,'Insufficient training class support')
    mean,scale,weights=weighted_scaler(x,arrays['train_lines'].tolist())
    np.testing.assert_allclose(arrays['scaler_mean'],mean,rtol=1e-10,atol=1e-12)
    np.testing.assert_allclose(arrays['scaler_scale'],scale,rtol=1e-10,atol=1e-12)
    tr=(x-mean)/scale; te=(test-mean)/scale; fit_checks={}
    for model,p in [('baseline',4),('extended',5)]:
        coef=arrays[model+'_coefficients']; logp=logits_logp(coef,te[:,:p])
        require(coef.shape==(p+1,9) and np.isfinite(coef).all() and np.isfinite(logp).all(),'Invalid fitted coefficients/logP')
        np.testing.assert_allclose(coef.mean(axis=1),0,atol=1e-10)
        np.testing.assert_allclose(logp,arrays[model+'_logp'],rtol=1e-9,atol=1e-10)
        objective,gradient,_=softmax_calculus(coef,tr[:,:p],arrays['train_y'],weights)
        require(np.max(np.abs(gradient))<=1e-6,'Saved fit is not an accepted stationary point')
        comparison=None
        if refit:
            own=independent_softmax_fit(tr[:,:p],arrays['train_y'],arrays['train_lines'].tolist())
            comparison=float(np.max(np.abs(coef-own['coefficients'])))
            np.testing.assert_allclose(coef,own['coefficients'],rtol=2e-5,atol=2e-5)
        fit_checks[model]={'objective':objective,'gradient_inf':float(np.max(np.abs(gradient))),'independent_coefficient_max_error':comparison}
    return {'status':'fit_state_verified','n_train_rows':len(x),'n_test_rows':len(test),'models':fit_checks}


def axis_digest(values: list[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def canonical_pair(row: dict) -> str:
    return json.dumps([row["protocol"],row["run"],row["line"],int(row["time"]),row["stimulus"]],separators=(",",":"))


def deterministic_subset(values: list[str], n: int) -> list[str]:
    require(len(values)==len(set(values)) and n>0,"Nonunique audit sampling frame")
    return sorted(values,key=lambda value:(hashlib.sha256(("macromap-independent-v1:"+value).encode()).hexdigest(),value))[:n]


def verify_normalized_subset(source: Path, normalized: np.ndarray, feature_ids: list[str], sample_ids: list[str], n: int=8) -> dict:
    """Read every gene for eight hash-selected samples, not eight genes.

    Pure read-only function. Real-source invocation is allowed only through the
    root-frozen CLI path. Synthetic fixtures may call this function directly.
    No source value is rewritten and unselected count tokens are not validated.
    """
    chosen=deterministic_subset(sample_ids,n); indices=[sample_ids.index(value) for value in chosen]
    require(normalized.shape==(len(feature_ids),len(sample_ids)),"Cache axes disagree")
    counts=np.empty((len(feature_ids),len(indices)),dtype=float);totals=[0]*len(indices);observed=0
    number=re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
    with gzip.open(source,"rt",encoding="utf-8") as stream:
        header=stream.readline().rstrip("\r\n").split("\t")
        require(header==["Chr","Geneid","Start","End","Strand","Length"]+sample_ids,"Raw source axis mismatch")
        for row,line in enumerate(stream):
            fields=line.rstrip("\r\n").split("\t")
            require(row<len(feature_ids) and len(fields)==len(sample_ids)+6 and fields[1]==feature_ids[row],"Raw feature/count-field axis mismatch")
            for j,index in enumerate(indices):
                token=fields[index+6].strip()
                require(bool(number.fullmatch(token)),"Invalid selected raw count syntax")
                try: value=Decimal(token)
                except InvalidOperation as exc: raise VerificationError("Invalid selected raw count") from exc
                require(value.is_finite() and 0<=value<2**53 and value==value.to_integral_value(),"Selected count not an exact nonnegative integer")
                integer=int(value);counts[row,j]=integer;totals[j]+=integer
            observed+=1
    require(observed==len(feature_ids) and all(0<value<2**53 for value in totals),"Raw rows/totals unavailable")
    expected=np.log2(1+(counts/np.asarray(totals,dtype=float))*1_000_000)
    actual=np.asarray(normalized[:,indices])
    require(np.isfinite(actual).all(),"Selected normalized values unavailable")
    np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=2e-12)
    return {"status":"normalized_source_subset_verified","sample_count":len(indices),"gene_count":len(feature_ids),
        "selection_sha256":axis_digest(chosen),"max_absolute_error":float(np.max(np.abs(actual-expected))),
        "scope":"all_released_genes_in_eight_hash_selected_samples_not_every_count_token"}


def verify_context_arrays(arrays: dict, audit: dict, normalized: np.ndarray, design: dict, *,
                          scope: str, protocol: str, time: int, fold: int, draws: int=1000, bins: int=20) -> tuple[dict,dict]:
    """Independent controls -> bins/draws -> sample means -> paired changes.

    Audits four training and four held-out predictive keys at fixed fold zero in
    both protocols/times. All nine programs are retained, including failures.
    Other folds' numerical state is verified separately against saved weights.
    """
    programs=list(ORIGINAL_PROGRAMS)+( [NONOVERLAP] )
    pairs=reconstruct_pairs(design["samples"],include_nomatch=scope!="mapped_primary")
    predictive=reconstruct_prediction_cohort(pairs)
    relevant=[r for r in predictive if r["protocol"]==protocol and int(r["time"])==time]
    full=[r for r in pairs if r["status"]=="available" and r["protocol"]==protocol and int(r["time"])==time]
    expected_indices=sorted({int(r[key]) for r in full for key in ["treatment_index","control_index"]})
    np.testing.assert_array_equal(arrays["selected_sample_indices"],expected_indices)
    np.testing.assert_array_equal(arrays["selected_sample_ids"],[design["sample_ids"][i] for i in expected_indices])
    np.testing.assert_array_equal(arrays["program_ids"],programs)
    np.testing.assert_array_equal(arrays["pool_exclusion_indices"],design["pool_exclusion_indices"])
    for key,values in [("feature_axis_sha256",design["feature_ids"]),("sample_axis_sha256",design["sample_ids"])]:
        require(str(arrays[key].item())==axis_digest(values),"Stratum axis hash differs")
    targets=[design["program_indices"][program] for program in programs]
    np.testing.assert_array_equal(arrays["target_gene_indices"],[i for group in targets for i in group])
    np.testing.assert_array_equal(arrays["target_gene_offsets"],np.cumsum([0]+[len(group) for group in targets]))
    reference,controls=training_reference(normalized,pairs,design["folds"],protocol,time,fold)
    np.testing.assert_allclose(arrays["reference"],reference,rtol=2e-13,atol=2e-13)
    np.testing.assert_array_equal(arrays["training_control_indices"],controls)
    np.testing.assert_array_equal(arrays["training_control_ids"],[design["sample_ids"][i] for i in controls])
    require(audit["scope"]==scope and audit["protocol"]==protocol and int(audit["time"])==time and int(audit["heldout_fold"])==fold,"Stratum identity differs")
    selected=[]
    for heldout in [False,True]:
        group=[r for r in relevant if (design["folds"][protocol,r["run"]]==fold)==heldout]
        by_key={canonical_pair(row):row for row in group}
        selected.extend(by_key[key] for key in deterministic_subset(list(by_key),4))
    score_indices=sorted({int(r[key]) for r in selected for key in ["control_index","treatment_index"]})
    positions={value:i for i,value in enumerate(expected_indices)}
    small_positions={value:i for i,value in enumerate(score_indices)}
    values=np.asarray(normalized[:,score_indices],float)
    paired={canonical_pair(row):np.full(len(programs),np.nan) for row in selected};failed=[]
    for j,program in enumerate(programs):
        target=design["program_indices"][program]
        own=independent_matches(reference,design["feature_ids"],target,design["pool_exclusion_indices"],program=program,
            protocol=protocol,time=time,fold=fold,bins=bins,draws=draws)
        source=audit["matching"][program]
        require(source["available"]==own["available"] and source["seed"]==own["seed"] and source["reason"]==own["reason"],"Matching availability/seed/reason differs")
        for key in ["target_bin_counts","pool_bin_counts"]:
            require({int(k):int(v) for k,v in source[key].items()}==own[key],"Matching target/pool counts differ")
        require(source["underflow_bins"]==own["underflow_bins"],"Matching underflow bins differ")
        np.testing.assert_array_equal(arrays["feature_bin"],own["feature_bin"])
        weights=np.zeros(len(reference));weights[target]=1/len(target)
        np.testing.assert_allclose(arrays["target_weights"][j],weights,rtol=0,atol=1e-15)
        raw=values[target].mean(axis=0);detected=(values[target]>0).mean(axis=0)
        at=[positions[index] for index in score_indices]
        np.testing.assert_allclose(arrays["sample_raw_scores"][at,j],raw,rtol=2e-12,atol=2e-12)
        np.testing.assert_allclose(arrays["sample_detected_fractions"][at,j],detected,rtol=0,atol=2e-12)
        status="mapping_failed" if not design["mapping_gates"][program] else own["reason"]
        require(audit["program_status"][program]==status,"Program eligibility status differs")
        require(bool(arrays["matching_available"][j])==(status=="available"),"Available matching mask differs")
        if own["available"]:
            require(source["selection_sha256"]==own["selection_sha256"],"Realized matching draw hash differs")
        if status!="available":
            require(np.isnan(arrays["background_weights"][j]).all() and np.isnan(arrays["sample_adjusted_scores"][at,j]).all(),"Unavailable matched score/weights replaced by a value")
            failed.append(program);continue
        np.testing.assert_allclose(arrays["background_weights"][j],own["background_weights"],rtol=1e-12,atol=1e-15)
        # Independent explicit per-draw matched means on the bounded sample panel.
        background=np.mean([values[indices].mean(axis=0) for indices in own["selections"]],axis=0)
        adjusted=raw-background
        np.testing.assert_allclose(arrays["sample_background_scores"][at,j],background,rtol=2e-12,atol=2e-12)
        np.testing.assert_allclose(arrays["sample_adjusted_scores"][at,j],adjusted,rtol=2e-12,atol=2e-12)
        for row in selected:
            paired[canonical_pair(row)][j]=adjusted[small_positions[row["treatment_index"]]]-adjusted[small_positions[row["control_index"]]]
    return {"status":"bounded_context_verified","identity_scope":scope,"protocol":protocol,"time":time,"fold":fold,
        "programs":len(programs),"paired_keys":len(selected),"paired_key_selection_sha256":axis_digest(sorted(paired)),
        "failed_programs":failed},paired


def pair_rows_from_arrays(arrays: dict, prefix: str) -> list[dict]:
    keys=["protocol","run","line","time","stimulus","control_id","treatment_id","control_index","treatment_index"]
    size=len(arrays[prefix+"_line"])
    require(all(arrays[prefix+"_"+key].shape==(size,) for key in keys),"Fit identity axis shape mismatch")
    return [{key:(int(arrays[prefix+"_"+key][i]) if key in {"time","control_index","treatment_index"} else str(arrays[prefix+"_"+key][i]))
             for key in keys} for i in range(size)]


def verify_fit_bundle(arrays: dict, audit: dict, expected_train: list[dict], expected_test: list[dict], *,
                      variant: str, paired_subset: dict | None=None, refit: bool=True) -> tuple[dict,list[dict]]:
    """Checks the actual v1 artifact schema and each expected source-linked row."""
    require(audit["schema"]=="macromap-private-audit-v1","Unknown fit audit schema")
    programs=list(COMPARATORS)+[variant];order=list(ORIGINAL_PROGRAMS)+[NONOVERLAP]
    np.testing.assert_array_equal(arrays["feature_programs"],programs)
    np.testing.assert_array_equal(arrays["class_order"],STIMULI)
    pairs={}
    for prefix,expected in [("train",expected_train),("test",expected_test)]:
        actual=pair_rows_from_arrays(arrays,prefix);pairs[prefix]=actual
        require([canonical_pair(row) for row in actual]==[canonical_pair(row) for row in expected],"Fit keys differ from expected metadata/folds")
        for row,source in zip(actual,expected):require(all(str(row[key])==str(source[key]) for key in row),"Fit source sample identity/index differs")
        labels=np.asarray([STIMULI.index(row["stimulus"]) for row in actual],dtype=int)
        np.testing.assert_array_equal(arrays[prefix+"_class_index"],labels)
        if prefix+"_pair_ids" in arrays:np.testing.assert_array_equal(arrays[prefix+"_pair_ids"],[canonical_pair(row) for row in actual])
        if paired_subset:
            for i,row in enumerate(actual):
                if canonical_pair(row) in paired_subset:
                    expected_x=paired_subset[canonical_pair(row)][[order.index(program) for program in programs]]
                    np.testing.assert_allclose(arrays[prefix+"_x_raw"][i],expected_x,rtol=2e-10,atol=2e-12,equal_nan=True)
    if audit["status"]!="available":
        require("baseline_log_probabilities" not in arrays and "extended_log_probabilities" not in arrays,"Failed fit has substituted predictions")
        return {"status":"expected_failed_fit_retained","fit_status":audit["status"]},[]
    adapted={"train_x":arrays["train_x_raw"],"test_x":arrays["test_x_raw"],
        "train_y":arrays["train_class_index"],"test_y":arrays["test_class_index"],
        "train_lines":arrays["train_line"],"test_lines":arrays["test_line"],
        "train_runs":arrays["train_run"],"test_runs":arrays["test_run"],
        **{key:arrays[key] for key in ["scaler_mean","scaler_scale","baseline_coefficients","extended_coefficients"]},
        "baseline_logp":arrays["baseline_log_probabilities"],"extended_logp":arrays["extended_log_probabilities"]}
    result=verify_fit_arrays(adapted,refit=refit)
    mean,scale,weights=weighted_scaler(adapted["train_x"],adapted["train_lines"].tolist())
    np.testing.assert_allclose(arrays["train_line_weights"],weights,rtol=1e-12,atol=1e-15)
    for prefix in ["train","test"]:np.testing.assert_allclose(arrays[prefix+"_x_standardized"],(adapted[prefix+"_x"]-mean)/scale,rtol=1e-12,atol=1e-12)
    losses={}
    for model in ["baseline","extended"]:
        loss=-adapted[model+"_logp"][np.arange(len(expected_test)),adapted["test_y"]]
        np.testing.assert_allclose(arrays[model+"_loss"],loss,rtol=1e-12,atol=1e-12);losses[model]=loss
        require(audit["fits"][model]["status"]=="available","Unqualified optimizer success")
        np.testing.assert_allclose(audit["fits"][model]["objective"],result["models"][model]["objective"],rtol=1e-9,atol=1e-10)
    np.testing.assert_allclose(arrays["gain"],losses["baseline"]-losses["extended"],rtol=1e-12,atol=1e-12)
    records=[{**row,"baseline_loss":float(losses["baseline"][i]),"extended_loss":float(losses["extended"][i]),
        "gain":float(losses["baseline"][i]-losses["extended"][i])} for i,row in enumerate(expected_test)]
    return result,records


def load_independent_design(directory: Path) -> dict:
    samples=read_rows(directory/"sample-map.csv")
    for row in samples:row["sample_index"]=int(row["sample_index"])
    features=read_rows(directory/"features.csv");members=read_rows(directory/"program-mapping.csv")
    programs=list(ORIGINAL_PROGRAMS)+[NONOVERLAP]
    indices={program:sorted({int(row["feature_index"]) for row in members if row["program_id"]==program
        and row["mapping_status"] in {"mapped","mapped_nonoverlap"}}) for program in programs}
    excluded=set().union(*(set(indices[program]) for program in ORIGINAL_PROGRAMS))
    excluded|={i for i,row in enumerate(features) if row["stable_feature_id"]=="ENSG00000157557" or row["feature_name"]=="ETS2"}
    return {"samples":samples,"sample_ids":[row["sample_id"] for row in samples],"feature_ids":[row["feature_id"] for row in features],
        "program_indices":indices,"pool_exclusion_indices":sorted(excluded),"folds":reconstruct_folds(samples),
        "mapping_gates":{row["program_id"]:row["mapping_gate_pass"]=="True" for row in read_rows(directory/"mapping-inventory.csv")}}


def contained_artifact(directory: Path, relative: str) -> Path:
    part=Path(relative)
    require(not part.is_absolute() and ".." not in part.parts and bool(part.parts),"Unsafe audit artifact path")
    path=directory/part
    require(path.resolve().is_relative_to(directory.resolve()),"Audit artifact escaped its bundle")
    require(not any(parent.is_symlink() for parent in [path,*path.parents]),"Audit artifact is aliased")
    return path


def reference_response(values: list[float], lines: list[str], runs: list[str], *, interval: bool=False,
                       draws: int=5000) -> dict:
    x=np.asarray(values,float)
    require(len(x)==len(lines)==len(runs) and len(lines)==len(set(lines)),"Repeated/misaligned response line")
    if not np.isfinite(x).all():return {"status":"missing_response","median":None}
    if len(x)<3:return {"status":"insufficient_paired_lines","median":None}
    ordered=sorted(set(runs));remaining=[[i for i,value in enumerate(runs) if value!=run] for run in ordered]
    loo=[float(np.median(np.r_[x[:i],x[i+1:]])) for i in range(len(x))] if len(x)>=4 else []
    lor=[float(np.median(x[group])) for group in remaining if len(group)>=3]
    result={"status":"available_descriptive","lines":len(x),"runs":len(ordered),"median":float(np.median(x)),
        "mean":float(np.mean(x)),"minimum":float(min(x)),"maximum":float(max(x)),
        "positive":int(sum(x>0)),"negative":int(sum(x<0)),"zero":int(sum(x==0)),
        "leave_one_line_out_median_range":[min(loo),max(loo)] if loo else None,
        "leave_one_run_out_median_range":[min(lor),max(lor)] if lor else None}
    if interval and len(ordered)>=3:
        rng=np.random.Generator(np.random.PCG64(2026091601));means=[];medians=[]
        groups=[[i for i,value in enumerate(runs) if value==run] for run in ordered]
        for _ in range(draws):
            indices=[i for draw in rng.integers(len(ordered),size=len(ordered)) for i in groups[draw]]
            means.append(float(np.mean(x[indices])));medians.append(float(np.median(x[indices])))
        result.update(mean_interval=np.quantile(means,[.025,.975],method="linear").tolist(),
                      median_interval=np.quantile(medians,[.025,.975],method="linear").tolist())
    return result


def verify_response_panel(panel: list[dict], rows: list[dict], *, interval_draws: int=5000) -> dict:
    """All points/statuses; CI audit fixed to CIL6h in both protocols, all programs.

    This interval subset is a computation check, not an endpoint filter. Every
    response cell remains in the original full panel. No pooled response CI.
    """
    expected={(p,t,s,g,q) for p in (*PROTOCOLS,"pooled_descriptive") for t in TIMES
        for s in STIMULI+("PIC",) for g in (*ORIGINAL_PROGRAMS,NONOVERLAP)
        for q in ([None] if s=="PIC" else ["raw_change","background_change","adjusted_change"])}
    actual=[(r["protocol"],int(r["time"]),r["stimulus"],r["program_id"],r.get("quantity")) for r in panel]
    require(len(actual)==len(set(actual)) and set(actual)==expected,"Incomplete response panel or extra endpoint")
    checked_intervals=0
    for item in panel:
        p,t,s,g,q=item["protocol"],int(item["time"]),item["stimulus"],item["program_id"],item.get("quantity")
        if s=="PIC":require(item["status"]=="unavailable_mock_identity","PIC received a qualified effect");continue
        selected=[r for r in rows if int(r["time"])==t and r["stimulus"]==s and r["program_id"]==g
            and (p=="pooled_descriptive" or r["protocol"]==p)]
        values=[float(r[q]) for r in selected];lines=[r["line"] for r in selected]
        runs=[r["protocol"]+":"+r["run"] for r in selected]
        ci=(p in PROTOCOLS and t==6 and s=="CIL")
        own=reference_response(values,lines,runs,interval=ci,draws=interval_draws)
        require(item["status"]==own["status"],"Response eligibility status differs")
        require(item["input_statuses"]==dict(Counter(r["status"] for r in selected)),"Response source statuses changed")
        require(item["expected_paired_lines"]==len(selected) and item["expected_paired_runs"]==len(set(runs)),"Response expected counts differ")
        for key,value in own.items():
            if key=="status":continue
            if value is None:require(item.get(key) is None,"Unavailable summary replaced")
            else:np.testing.assert_allclose(item[key],value,rtol=2e-11,atol=2e-12)
        if "mean_interval" in own:checked_intervals+=1
        if p=="pooled_descriptive":
            require(item.get("median_interval") is None and item.get("mean_interval") is None,"Unplanned pooled response CI")
    return {"complete_cells":len(panel),"all_points_statuses_checked":True,"protocol_CIL6h_interval_cells_checked":checked_intervals,
        "pooled_response_intervals":False,"other_response_interval_values":"not_numerically_recomputed"}


def verify_engine_bundle(directory: Path, normalized: np.ndarray, design: dict, *, matching_draws: int=1000,
                         bootstrap_draws: int=5000, refit: bool=True) -> dict:
    """Read-only engine bundle check; caller authenticates freeze and all files.

    Pure-input entry also supports explicitly synthetic tests with fewer draws.
    The real CLI never changes the fixed1000/5000 settings.
    """
    manifest=json.loads((directory/"fit-manifest.json").read_text())
    scopes=["mapped_primary","all_source_lines_sensitivity"];variants=["ets2_g1_dn",NONOVERLAP]
    expected={(s,p,t,f,g) for s in scopes for p in PROTOCOLS for t in TIMES for f in range(5) for g in variants}
    keys=[(r["identity_scope"],r["protocol"],int(r["time"]),int(r["fold"]),r["variant"]) for r in manifest]
    require(len(keys)==len(set(keys)) and set(keys)==expected,"Expected fit cells lost or added")
    by_key=dict(zip(keys,manifest));fit_checks=[];contexts=[];endpoints=[];response_checks=[]
    program_order=list(ORIGINAL_PROGRAMS)+[NONOVERLAP]
    for scope in scopes:
        all_pairs=reconstruct_pairs(design["samples"],include_nomatch=scope!="mapped_primary")
        predictive=reconstruct_prediction_cohort(all_pairs)
        saved_pairs=read_rows(directory/scope/"expected-prediction-pairs.csv")
        require(exact_keys(saved_pairs)==exact_keys(predictive),"Saved predictive cohort changed")
        losses={g:[] for g in variants};expected_responses=[]
        for p in PROTOCOLS:
            for t in TIMES:
                eligible=[r for r in all_pairs if r["status"]=="available" and r["protocol"]==p and r["time"]==t]
                model_pairs=[r for r in predictive if r["protocol"]==p and r["time"]==t]
                for f in range(5):
                    folder=directory/scope/f"{p}-t{t}-fold{f}"
                    with np.load(folder/"score-reference-arrays.npz",allow_pickle=False) as state:arrays=dict(state)
                    audit=json.loads((folder/"score-reference-audit.json").read_text())
                    indices=sorted({r[key] for r in eligible for key in ["treatment_index","control_index"]})
                    np.testing.assert_array_equal(arrays["selected_sample_indices"],indices)
                    np.testing.assert_array_equal(arrays["selected_sample_ids"],[design["sample_ids"][i] for i in indices])
                    np.testing.assert_array_equal(arrays["program_ids"],program_order)
                    controls={r["line"]:r["control_index"] for r in eligible if design["folds"][p,r["run"]]!=f}
                    np.testing.assert_array_equal(arrays["training_control_indices"],[controls[line] for line in sorted(controls)])
                    paired=None
                    if f==0:
                        checked,paired=verify_context_arrays(arrays,audit,normalized,design,scope=scope,protocol=p,time=t,fold=f,draws=matching_draws)
                        contexts.append(checked)
                    positions={index:i for i,index in enumerate(indices)}
                    for row in eligible:
                        if design["folds"][p,row["run"]]!=f:continue
                        tr,ctrl=positions[row["treatment_index"]],positions[row["control_index"]]
                        for j,g in enumerate(program_order):
                            expected_responses.append({**row,"program_id":g,"heldout_fold":f,"status":audit["program_status"][g],
                                **{q:float(arrays[key][tr,j]-arrays[key][ctrl,j]) for q,key in [
                                    ("raw_change","sample_raw_scores"),("background_change","sample_background_scores"),("adjusted_change","sample_adjusted_scores")]}})
                    train=[r for r in model_pairs if design["folds"][p,r["run"]]!=f]
                    test=[r for r in model_pairs if design["folds"][p,r["run"]]==f]
                    for g in variants:
                        spec=by_key[scope,p,t,f,g]
                        expected_npz=str((folder/g/"fit-arrays.npz").relative_to(directory))
                        require(spec["fit_npz"]==expected_npz,"Fit artifact path/context mismatch")
                        npz=pinned(contained_artifact(directory,spec["fit_npz"]),spec["fit_npz_sha256"])
                        with np.load(npz,allow_pickle=False) as state:fit=dict(state)
                        fit_audit=json.loads(contained_artifact(directory,spec["fit_audit"]).read_text())
                        require(all(fit_audit[k]==spec[k] for k in ["identity_scope","protocol","time","fold","variant","status"]),"Fit manifest/audit identity mismatch")
                        # All folds: require saved X to equal paired changes from that
                        # same outer-fold score state, never the global OOF table.
                        columns=[program_order.index(name) for name in (*COMPARATORS,g)]
                        for prefix,rows in [("train",train),("test",test)]:
                            matrix=np.asarray([arrays["sample_adjusted_scores"][positions[r["treatment_index"]],columns]-
                                arrays["sample_adjusted_scores"][positions[r["control_index"]],columns] for r in rows]).reshape(len(rows),5)
                            np.testing.assert_allclose(fit[prefix+"_x_raw"],matrix,rtol=1e-13,atol=1e-13,equal_nan=True)
                        checked,records=verify_fit_bundle(fit,fit_audit,train,test,variant=g,paired_subset=paired,refit=refit)
                        fit_checks.append({"identity_scope":scope,"protocol":p,"time":t,"fold":f,"variant":g,**checked})
                        losses[g].extend(records)
        response_rows=read_rows(directory/scope/"line-responses.csv")
        to_key=lambda row:(canonical_pair(row),row["program_id"])
        observed={to_key(row):row for row in response_rows}
        require(len(observed)==len(response_rows)==len(expected_responses) and set(observed)=={to_key(r) for r in expected_responses},"Lost/extra response row")
        for row in expected_responses:
            other=observed[to_key(row)]
            require(other["status"]==row["status"] and int(other["heldout_fold"])==row["heldout_fold"],"Response status/fold changed")
            for q in ["raw_change","background_change","adjusted_change"]:np.testing.assert_allclose(float(other[q]),row[q],rtol=1e-12,atol=1e-12,equal_nan=True)
        panel=json.loads((directory/scope/"response-panel.json").read_text())
        response_checks.append({"identity_scope":scope,**verify_response_panel(panel,response_rows,interval_draws=bootstrap_draws)})
        for g in variants:
            endpoint=json.loads((directory/scope/g/"incremental-endpoint.json").read_text())
            if exact_keys(losses[g])!=exact_keys(predictive):
                require(endpoint["status"]!="available","Failed expected fit was silently omitted")
                endpoints.append({"identity_scope":scope,"variant":g,"status":"incomplete_expected_predictions_retained"});continue
            own=conditional_run_bootstrap(losses[g],predictive,draws=bootstrap_draws)
            require(endpoint["status"]=="available" and endpoint["protocol_weights"]==own["point"]["protocol_weights"],"Primary availability/protocol weights differ")
            np.testing.assert_allclose(endpoint["primary_gain"],own["point"]["primary_gain"],rtol=1e-11,atol=1e-12)
            np.testing.assert_allclose(endpoint["interval"],own["interval"],rtol=1e-11,atol=1e-12)
            require(endpoint["replicates"]==bootstrap_draws and endpoint["seed"]==2026091601 and
                endpoint["interval_status"]=="fixed_prediction_conditional_only","Predictive interval convention differs")
            for cell,point in own["point"]["cells"].items():
                for key,value in point.items():np.testing.assert_allclose(endpoint["cells"][str(cell)][key],value,rtol=1e-11,atol=1e-12)
                np.testing.assert_allclose(endpoint["cell_intervals"][str(cell)],own["cell_intervals"][cell],rtol=1e-11,atol=1e-12)
            csv_rows=read_rows(directory/scope/g/"heldout-loss-records.csv")
            require(exact_keys(csv_rows)==exact_keys(losses[g]),"Held-out loss export keys differ")
            saved={canonical_pair(row):row for row in csv_rows}
            for row in losses[g]:
                for key in ["baseline_loss","extended_loss","gain"]:np.testing.assert_allclose(float(saved[canonical_pair(row)][key]),row[key],rtol=1e-12,atol=1e-12)
            influence=json.loads((directory/scope/g/"run-influence.json").read_text())
            line_losses,_=loss_records(losses[g],predictive)
            run_keys={(p,r["run"]) for p in PROTOCOLS for r in predictive if r["protocol"]==p}
            require(len(influence)==len(run_keys) and {(r["protocol"],r["omitted_run"]) for r in influence}==run_keys,"Missing leave-run influence cell")
            for omitted in influence:
                total=0.
                for p in PROTOCOLS:
                    time_means=[]
                    for t in TIMES:
                        values=[value[0]-value[1] for key,value in line_losses.items() if key[0]==p and key[3]==t
                            and (key[0],key[1])!=(omitted["protocol"],omitted["omitted_run"])]
                        require(bool(values),"Empty leave-run protocol/time cell")
                        time_means.append(float(np.mean(values)))
                    total+=own["point"]["protocol_weights"][p]*float(np.mean(time_means))
                require(omitted["status"]=="available","Unavailable leave-run influence cell")
                np.testing.assert_allclose(omitted["gain"],total,rtol=1e-11,atol=1e-12)
            endpoints.append({"identity_scope":scope,"variant":g,"status":"primary_point_cells_conditional_intervals_and_run_influence_verified"})
    return {"status":"independently_verified_with_declared_source_subset","fit_cells":len(fit_checks),"fit_checks":fit_checks,
        "source_score_contexts":contexts,"endpoints":endpoints,"responses":response_checks,
        "scope_limits":["Independent raw normalization uses eight samples, not every source count.",
            "Independent reference/bin/draw/score reconstruction uses fold0 only; all folds have identity/X/model/loss checks.",
            "Response CI values independently replayed for CIL6h only, all programs/quantities in both protocols.",
            "Bootstrap conditions on fixed OOF scoring/models; it does not refit the complete procedure."]}


def verify_execution(directory: Path, plan: dict, plan_sha256: str, complete_sha256: str, preparation: Path, *, root: Path=ROOT) -> dict:
    """Authenticated read-only real-output path. No import of primary methods."""
    import platform
    import scipy
    require(plan.get("status")=="frozen" and plan.get("execute_real_expression") is True and plan.get("analysis_id")=="macromap-fixed-program-context-v1","Execution is not explicitly authorized")
    require(directory.resolve()==(root/plan["output_directory"]).resolve() and preparation.resolve()==(root/plan["prepared_directory"]).resolve(),"Plan input/output location differs")
    require(directory.resolve().is_relative_to((root/"work/macromap-design").resolve()),"Execution not in private work area")
    pinned(root/"scripts/analyze_macromap_program_context.py",plan["script_sha256"])
    pinned(root/"scripts/run_macromap_program_context.py",plan["runner_script_sha256"])
    runtime={"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__}
    require(plan["runtime"]==runtime,"Verification runtime differs from frozen numerical runtime")
    complete=json.loads(pinned(directory/"execution-complete.json",complete_sha256).read_text())
    require(complete["stage"]=="execution_complete_all_prespecified_statuses_retained" and complete["plan_sha256"]==plan_sha256,"Incomplete/foreign execution receipt")
    for key in ["script_sha256","runner_script_sha256","runtime","prepared_summary_sha256","method_contract"]:
        require(complete[key]==plan[key],"Execution/code/software/preparation pin differs")
    start=datetime.fromisoformat(complete["started_at_utc"]);finish=datetime.fromisoformat(complete["completed_at_utc"])
    freeze=datetime.fromisoformat(plan["frozen_at_utc"])
    require(freeze<=start<=finish<=datetime.now(timezone.utc),"Execution was outside its freeze chronology")
    require(not (directory/"execution-failed.json").exists(),"Conflicting failed/completed receipts")
    files={str(path.relative_to(directory)) for path in directory.rglob("*") if path.is_file() and path.name!="execution-complete.json"}
    require(files==set(complete["artifact_sha256"]),"Complete artifact inventory differs from directory")
    for name,sha in complete["artifact_sha256"].items():pinned(contained_artifact(directory,name),sha)
    require({"execution-start.json","fit-manifest.json","released-log2cpm.f64","released-log2cpm.f64.json","released-log2cpm.f64.lock"}<=files,"Missing complete artifact inventory")
    summary=json.loads(pinned(preparation/"source-only-readiness.json",plan["prepared_summary_sha256"]).read_text())
    require(complete["source_pins"]==summary["source_pins"],"Execution original-source pins differ")
    design=load_independent_design(preparation);cache=complete["cache"]
    require(cache["complete"] is True and cache["dtype"]=="<f8" and cache["shape"]==[len(design["feature_ids"]),len(design["sample_ids"])],"Invalid normalized cache state")
    require(json.loads((directory/"released-log2cpm.f64.json").read_text())==cache and cache["all_library_totals_positive"] is True,"Cache sidecar differs from completed receipt")
    require((directory/"released-log2cpm.f64").stat().st_size==8*int(np.prod(cache["shape"])),"Cache byte extent differs from axes")
    require(cache["cache_sha256"]==complete["artifact_sha256"]["released-log2cpm.f64"] and
        cache["feature_id_sha256"]==axis_digest(design["feature_ids"]) and cache["sample_id_sha256"]==axis_digest(design["sample_ids"]),"Normalized cache axis/hash differs")
    raw_relative="data/raw/macromap-qualification/MacroMap_raw_expression.txt.gz"
    require(cache["source_sha256"]==summary["source_pins"][raw_relative]["sha256"],"Normalized source not frozen raw matrix")
    normalized=np.memmap(directory/"released-log2cpm.f64",dtype="<f8",mode="r",shape=tuple(cache["shape"]))
    try:
        source_check=verify_normalized_subset(root/raw_relative,normalized,design["feature_ids"],design["sample_ids"])
        result=verify_engine_bundle(directory,normalized,design)
    finally:normalized._mmap.close()
    return {**result,"normalization":source_check,"effects_read":True,"verification_is_read_only":True,
        "execution_receipt_sha256":complete_sha256,"plan_sha256":plan_sha256}


def _frozen_plan(path: Path, expected_sha256: str) -> dict:
    plan=json.loads(pinned(path,expected_sha256).read_text())
    require(plan.get('status')=='frozen','Real-output verification requires root freeze')
    stamp=datetime.fromisoformat(plan['frozen_at_utc'].replace('Z','+00:00'))
    require(stamp.tzinfo is not None and stamp<=datetime.now(timezone.utc),'Invalid freeze timestamp')
    return plan


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preparation',type=Path,required=True)
    parser.add_argument('--preparation-sha256',required=True)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--execution-directory',type=Path,help='Optional completed private v1 bundle; post-freeze only')
    parser.add_argument('--execution-complete-sha256')
    parser.add_argument('--verifier-sha256',help='Externally supplied frozen independent verifier pin')
    parser.add_argument('--plan',type=Path)
    parser.add_argument('--plan-sha256')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    plan=None
    if args.execution_directory:
        require(args.plan is not None and args.plan_sha256 is not None,'Missing root freeze pin')
        plan=_frozen_plan(args.plan,args.plan_sha256)
        require(plan.get('execute_real_expression') is True,'No root permission for real numerical outputs')
        pinned(Path(__file__),args.verifier_sha256)
        require(args.preparation_sha256==plan['prepared_summary_sha256'],'Preparation differs from freeze')
    else:
        require(not any([args.plan,args.plan_sha256,args.execution_complete_sha256,args.verifier_sha256]),'Numerical audit flags require execution-directory')
    require('..' not in args.output.parts and not any(path.is_symlink() for path in [args.output,*args.output.parents]),'Verification output path is aliased')
    require(not args.output.exists(),'Do not overwrite verification output')
    protected=[args.preparation,args.root/'data/raw']+([args.execution_directory] if args.execution_directory else [])
    require(not any(args.output.resolve().is_relative_to(path.resolve()) for path in protected),'Keep verification report outside immutable input bundles')
    source=verify_preparation(args.preparation,args.preparation_sha256,args.root)
    result=source
    if plan is not None:
        result=verify_execution(args.execution_directory,plan,args.plan_sha256,args.execution_complete_sha256,args.preparation,root=args.root)
        result['source_preparation']=source
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
    print(json.dumps({'status':result['status'],'effects_read':result['effects_read']}))


if __name__=='__main__':main()
