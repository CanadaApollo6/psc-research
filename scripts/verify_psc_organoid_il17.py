#!/usr/bin/env python3
"""Independent numerical verifier for the frozen GSE239283 organoid analysis.

No analysis-module imports, fitting of dispersions, downloads or input mutation.
All calculations use the pinned native inputs and exported fit arrays. This does
not establish exchangeability, NB small-sample calibration or clinical efficacy.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import re
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

REPO = Path(__file__).resolve().parents[1]
DONORS = {"PSC_4": "PSC", "PSC_6": "PSC", "PSC_8": "PSC", "PSC_9": "PSC",
          "nonPSC_2": "nonPSC", "nonPSC_3": "nonPSC", "nonPSC_4": "nonPSC", "nonPSC_5": "nonPSC"}
RIDGE = 1e-6
REFIT_N = 16
REFIT_SELECTION = "first16_by_sha256_utf8(PSC-organoid-conditional-NB-v1: concatenated_with_exact_source_ID)_among_eligible_genes_with_all16_counts>=10"
FILTER = {"min_count": 10, "min_libraries": 4}
MODEL = {
    "design": "Intercept+donor+treatment_IL17+disease_PSC:treatment_IL17",
    "primary_contrast": "PSC_response_minus_nonPSC_response",
    "size_factors_fit_type": "ratio", "fit_type": "parametric",
    "dispersion_fallback": "builtin_mean_if_parametric_fails", "refit_cooks": False,
    "cooks_filter": False, "independent_filter": False, "lfc_shrinkage": False,
    "alpha": 0.05, "multiplicity": "BH_all_fixed_eligible_features_NA_as_1_correction_only",
    "dependence_sensitivity": "BY_same_family",
    "count_boundary_rule": "group_response_needs_both_arms_positive_interaction_needs_all4_arms_positive",
    "zero_total_donor_pairs": "flag_only_no_new_feature_filter",
}
SENSITIVITY = {
    "transformation": "log2(CPM+1)", "library_size": "all_unfiltered_released_features",
    "unit": "one_IL17_minus_CNT_change_per_donor", "test": "Welch_t_with_Satterthwaite_df",
    "leaveout": "one_donor_out_effect_and_sign_only_fixed_features",
    "allocation": "all_70_4v4_labels_absolute_mean_change_difference_descriptive_only",
    "allocation_tie_atol": 1e-12, "allocation_tie_rtol": 1e-12,
}
SOFTWARE = {"python", "numpy", "pandas", "scipy", "pydeseq2", "anndata", "formulaic", "joblib"}
PUBLIC_FILES = ["organoid-feature-universe.csv.gz", "organoid-primary.csv.gz",
                "organoid-within-group.csv.gz", "organoid-welch.csv.gz",
                "organoid-leave-one-donor-out.csv.gz", "organoid-allocations.csv.gz", "organoid-qc.json"]
PRIVATE_FILES = ["execution-plan.json", "selected-libraries.csv", "design-matrix.csv",
                 "library-normalization.csv", "NB-natural-log-coefficients.csv.gz",
                 "NB-all-native-feature-estimates-and-flags.csv.gz", "all-cooks-distances.csv.gz",
                 "donor-paired-log2cpm-changes.csv.gz", "all-donor-leaveout-effects.csv.gz",
                 "all-70-allocation-effects.csv.gz", "all-70-allocation-labels.csv", "fit-warnings.json"]


class VerificationError(ValueError):
    """A source, freeze, output or independent arithmetic check failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def utc(value: Any) -> datetime:
    try:
        stamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (ValueError, TypeError, AttributeError) as error:
        raise VerificationError("Missing/invalid UTC timestamp") from error
    require(stamp.tzinfo is not None and stamp.utcoffset().total_seconds() == 0,
            "Timestamp must explicitly use UTC")
    require(stamp <= datetime.now(timezone.utc), "Future timestamp")
    return stamp


def check_pin(spec: dict, repo: Path = REPO) -> Path:
    require(isinstance(spec, dict) and bool(spec.get("path")), "Missing artifact path")
    require(bool(re.fullmatch(r"[0-9a-f]{64}", str(spec.get("sha256", "")))), "Invalid SHA-256 pin")
    p = Path(spec["path"])
    p = (p if p.is_absolute() else repo / p).resolve()
    require(p.is_file() and not p.name.endswith(".prefix"), "Missing complete pinned artifact")
    require(digest(p) == spec["sha256"], f"Artifact hash mismatch: {p.name}")
    require("bytes" not in spec or p.stat().st_size == spec["bytes"], f"Artifact length mismatch: {p.name}")
    return p


def exact_names(values: Any, label: str) -> list[str]:
    names = list(values)
    require(bool(names) and all(isinstance(x, str) and x and x == x.strip() for x in names), f"Empty/padded {label}")
    require(len(set(names)) == len(names), f"Duplicate {label}")
    return names


def read_table(path: Path, separator: str = ",") -> pd.DataFrame:
    require(separator in {",", "\t"}, "Unsupported separator")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter=separator, strict=True)
        columns = exact_names(next(reader, []), "table header")
        rows = list(reader)
    require(bool(rows) and all(len(row) == len(columns) for row in rows), "Empty/ragged table")
    return pd.DataFrame(rows, columns=columns, dtype=str)


def indexed(path: Path, id_column: str, expected: list[str], separator: str = ",",
            permit_reorder: bool = False) -> pd.DataFrame:
    frame = read_table(path, separator)
    require(id_column in frame, f"Missing identifier column {id_column}")
    names = exact_names(frame[id_column], id_column)
    require(set(names) == set(expected) if permit_reorder else names == expected,
            f"Non-exact axis: {path.name}")
    return frame.set_index(id_column, drop=False).loc[expected]


def numeric(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    require(set(columns).issubset(frame.columns), "Missing numeric output columns")
    try:
        return frame[columns].replace("", np.nan).to_numpy(dtype=float)
    except (ValueError, TypeError) as error:
        raise VerificationError("Non-numeric output token (only blank is numerical NA)") from error


def integer_tokens(values: Any, positive: bool = False) -> np.ndarray:
    array = np.asarray(values)
    require(array.ndim in {1, 2} and array.size > 0, "Empty/incompatible integer array")
    output = np.empty(array.shape, dtype=np.int64)
    for position, token in np.ndenumerate(array):
        require(isinstance(token, str) and bool(re.fullmatch(r"[0-9]+", token)), "Non-exact integer token")
        token = token.lstrip("0") or "0"
        require(len(token) <= 16 and int(token) <= 2**53 - 1, "Integer beyond exact validation range")
        value = int(token)
        require(not positive or value > 0, "Expected positive integer")
        output[position] = value
    return output


@dataclass
class Audit:
    checks: list[dict] = field(default_factory=list)

    def equal(self, name: str, actual: Any, expected: Any) -> None:
        a, b = np.asarray(actual), np.asarray(expected)
        require(a.shape == b.shape and np.array_equal(a, b), f"Exact mismatch: {name}")
        self.checks.append({"check": name, "values": int(a.size), "kind": "exact"})

    def close(self, name: str, actual: Any, expected: Any, rtol: float = 2e-8,
              atol: float = 2e-10) -> None:
        a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        require(a.shape == b.shape, f"Shape mismatch: {name}")
        ok = np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
        require(bool(ok.all()), f"Numerical/missingness mismatch: {name}; {int((~ok).sum())} values")
        finite = np.isfinite(a) & np.isfinite(b)
        self.checks.append({"check": name, "values": int(a.size), "finite_values": int(finite.sum()),
                            "kind": "numeric", "rtol": rtol, "atol": atol,
                            "maximum_absolute_difference": float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else None})

    def p(self, name: str, actual: Any, expected: Any) -> None:
        # No ordinary absolute tolerance that would forgive wrong tiny P/q values.
        self.close(name, actual, expected, rtol=2e-7, atol=float(np.nextafter(0., 1.)))


def adjust(pvalues: Any, by: bool = False) -> np.ndarray:
    """Step-up correction implemented as a backward running scalar minimum."""
    p = np.asarray(pvalues, dtype=float)
    require(p.ndim == 1 and not np.isinf(p).any() and not ((p < 0) | (p > 1)).any(), "Invalid P family")
    m = len(p)
    ranked = sorted(range(m), key=lambda i: (1. if np.isnan(p[i]) else p[i], i))
    result = np.empty(m)
    multiplier = m * (sum(1. / j for j in range(1, m+1)) if by else 1.)
    running = 1.
    for k in range(m, 0, -1):
        index = ranked[k-1]
        value = 1. if np.isnan(p[index]) else p[index]
        running = min(running, multiplier * value / k)
        result[index] = running
    return result


def filter_and_ratio(counts: np.ndarray) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, int]:
    """Gene-major counts; normalization basis comes from eligible genes only."""
    require(counts.ndim == 2 and counts.shape[1] == 16 and np.isfinite(counts).all()
            and (counts >= 0).all() and (counts == np.floor(counts)).all(), "Invalid complete count matrix")
    high = np.count_nonzero(counts >= 10, axis=1)
    eligible = high >= 4
    basis = eligible & np.all(counts > 0, axis=1)
    require(eligible.any() and basis.any(), "Empty eligible/all-positive normalization basis; no fallback")
    log_counts = np.log(counts[basis].astype(float))
    log_ratios = log_counts - np.sum(log_counts, axis=1)[:, None] / 16
    ordered = np.sort(log_ratios, axis=0)
    n = len(ordered)
    center = ordered[n//2] if n % 2 else (ordered[n//2-1] + ordered[n//2]) / 2
    factors = np.exp(center)
    diagnostics = {"libraries_count_ge_threshold": high,
                   "libraries_nonzero": np.count_nonzero(counts, axis=1),
                   "all_zero_released": np.all(counts == 0, axis=1),
                   "mean_raw_count": counts.mean(axis=1)}
    return eligible, diagnostics, factors, n


def make_design(libraries: pd.DataFrame) -> tuple[np.ndarray, list[str], list[str], np.ndarray]:
    required = {"library_id", "donor_id", "disease", "treatment", "n_cells", "total_released_counts"}
    require(required.issubset(libraries), "Incomplete library metadata")
    exact_names(libraries.library_id, "library IDs")
    donors = sorted(set(libraries.donor_id))
    require(set(donors) == set(DONORS) and len(libraries) == 16, "Not the fixed eight-donor panel")
    for donor in donors:
        rows = libraries.loc[libraries.donor_id == donor]
        require(len(rows) == 2 and sorted(rows.treatment) == ["CNT", "IL17"]
                and set(rows.disease) == {DONORS[donor]}, "Incomplete/inconsistent donor pair")
    integer_tokens(libraries.n_cells.to_numpy(), positive=True)
    column_names = ["Intercept"] + ["donor_" + d for d in donors[1:]] + ["treatment_IL17", "disease_PSC:treatment_IL17"]
    x = np.zeros((16, 10)); x[:, 0] = 1
    for i, row in enumerate(libraries.itertuples(index=False)):
        if row.donor_id != donors[0]:
            x[i, donors.index(row.donor_id)] = 1
        x[i, -2] = row.treatment == "IL17"
        x[i, -1] = row.treatment == "IL17" and row.disease == "PSC"
    require(np.linalg.matrix_rank(x) == 10, "Incorrect rank")
    leverage = np.diag(x @ np.linalg.solve(x.T @ x, x.T))
    require(np.allclose(leverage, 0.625, atol=1e-12, rtol=0), "Unexpected balanced-design leverage")
    return x, column_names, donors, leverage


def changes_from_counts(counts: np.ndarray, libraries: pd.DataFrame,
                        donors: list[str], eligible: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=0)
    require((totals > 0).all(), "Zero-total raw library")
    # Compute full-universe denominators before subsetting any genes.
    expression = np.log1p(1_000_000 * counts[eligible].astype(float) / totals) / np.log(2)
    pairs = []
    for donor in donors:
        idx = np.flatnonzero(libraries.donor_id.to_numpy() == donor)
        treatment = libraries.treatment.to_numpy()[idx]
        ctrl, stim = idx[treatment == "CNT"][0], idx[treatment == "IL17"][0]
        pairs.append(expression[:, stim] - expression[:, ctrl])
    return np.column_stack(pairs)


def donor_sensitivity(changes: np.ndarray, psc: np.ndarray) -> dict[str, Any]:
    y = np.asarray(changes, dtype=float); psc = np.asarray(psc, dtype=bool)
    require(y.ndim == 2 and y.shape[1] == 8 and psc.shape == (8,) and psc.sum() == 4
            and np.isfinite(y).all(), "Invalid eight complete donor changes; no case deletion")
    a, b = y[:, psc], y[:, ~psc]
    ma, mb = a.mean(axis=1), b.mean(axis=1)
    effect = ma - mb
    va = np.sum((a - ma[:, None])**2, axis=1) / 12
    vb = np.sum((b - mb[:, None])**2, axis=1) / 12
    se = np.sqrt(va + vb)
    usable = ((va > 0) | (vb > 0)) & (se > 32 * np.finfo(float).eps * np.maximum(1, np.abs(effect)))
    df = np.full(len(y), np.nan); t = df.copy()
    w = np.divide(va, va + vb, out=np.zeros_like(va), where=usable)
    df[usable] = 3 / (w[usable]**2 + (1-w[usable])**2)
    t[usable] = effect[usable] / se[usable]
    probability = 2 * special.stdtr(df, -np.abs(t))
    margin = special.stdtrit(df, 0.975) * se
    welch = {"mean_change_PSC": ma, "mean_change_nonPSC": mb, "effect_log2_cpm": effect,
             "se": se, "t": t, "df": df, "pvalue": probability,
             "q_bh_sensitivity": adjust(probability), "q_by_sensitivity": adjust(probability, by=True),
             "ci95_lower": effect-margin, "ci95_upper": effect+margin,
             "n_PSC_donors": np.full(len(y), 4), "n_nonPSC_donors": np.full(len(y), 4)}
    welch_status = np.where(usable, "tested_sensitivity", "zero_or_numerically_degenerate_variance_unavailable")
    leaveout = np.empty((8, len(y)))
    for omitted in range(8):
        keep = np.arange(8) != omitted
        leaveout[omitted] = y[:, psc & keep].mean(axis=1) - y[:, ~psc & keep].mean(axis=1)
    full_sign = np.sign(effect)
    leave_summary = {"full_effect_log2_cpm": effect, "minimum_leaveout_effect": leaveout.min(axis=0),
                     "maximum_leaveout_effect": leaveout.max(axis=0), "n_estimable_omissions": np.full(len(y), 8),
                     "n_same_nonzero_sign": ((np.sign(leaveout) == full_sign) & (leaveout != 0) & (effect != 0)).sum(axis=0),
                     "n_opposite_sign": ((np.sign(leaveout) == -full_sign) & (leaveout != 0) & (effect != 0)).sum(axis=0),
                     "n_zero_effect": (leaveout == 0).sum(axis=0)}
    # Independently enumerate bit masks, then use the documented lexicographic order.
    tuples = sorted(tuple(j for j in range(8) if mask & (1 << j)) for mask in range(256) if mask.bit_count() == 4)
    allocations = np.empty((70, len(y)))
    for i, chosen in enumerate(tuples):
        mask = np.isin(np.arange(8), chosen)
        allocations[i] = y[:, mask].mean(axis=1) - y[:, ~mask].mean(axis=1)
    tolerance = 1e-12 * (1 + np.abs(effect))
    extreme = np.sum(np.abs(allocations) >= np.maximum(0, np.abs(effect) - tolerance), axis=0)
    allocation_summary = {"observed_mean_change_difference": effect, "n_allocations": np.full(len(y), 70),
                          "n_as_or_more_extreme": extreme, "descriptive_tail_fraction": extreme / 70}
    require((extreme >= 2).all() and (extreme % 2 == 0).all(), "Complement/tail invariant failed")
    return {"welch": welch, "welch_status": welch_status, "leaveout": leaveout,
            "leaveout_summary": leave_summary, "allocations": allocations,
            "allocation_summary": allocation_summary, "assignment_indices": tuples}


def conditional_covariance(x: np.ndarray, means: np.ndarray, dispersion: float) -> np.ndarray:
    require(np.isfinite(means).all() and (means > 0).all() and np.isfinite(dispersion)
            and dispersion > 0, "Invalid fitted means/dispersion")
    w = 1 / (1 / means + dispersion)
    gram = x.T @ (w[:, None] * x)
    penalized = gram + RIDGE * np.eye(x.shape[1])
    left = np.linalg.solve(penalized, gram)
    return np.linalg.solve(penalized, left.T).T


def nb_reconstruction(x: np.ndarray, factors: np.ndarray, counts: np.ndarray,
                      beta: np.ndarray, dispersions: np.ndarray) -> dict[str, Any]:
    n = len(counts)
    require(beta.shape == (n, x.shape[1]) and dispersions.shape == (n,), "Native NB array shape mismatch")
    vectors = {"interaction": np.eye(x.shape[1])[-1], "nonPSC_treatment": np.eye(x.shape[1])[-2],
               "PSC_treatment": np.eye(x.shape[1])[-2:].sum(axis=0)}
    result = {key: {field: np.full(n, np.nan) for field in ["log2FoldChange", "lfcSE", "stat", "pvalue"]}
              for key in vectors}
    covariance_ti = np.full(n, np.nan)
    means = np.exp(beta @ x.T) * factors
    for g in range(n):
        if not (np.isfinite(beta[g]).all() and np.isfinite(dispersions[g]) and dispersions[g] > 0
                and np.isfinite(means[g]).all() and (means[g] > 0).all()):
            continue
        covariance = conditional_covariance(x, means[g], float(dispersions[g]))
        covariance_ti[g] = covariance[-2, -1] / np.log(2)**2
        for label, vector in vectors.items():
            effect = vector @ beta[g]
            se = np.sqrt(vector @ covariance @ vector)
            result[label]["log2FoldChange"][g] = effect / np.log(2)
            result[label]["lfcSE"][g] = se / np.log(2)
            result[label]["stat"][g] = effect / se
            result[label]["pvalue"][g] = 2 * special.ndtr(-abs(effect / se))
    for table in result.values():
        table["baseMean"] = (counts / factors).mean(axis=1)
        table["q_bh"] = adjust(table["pvalue"]); table["q_by"] = adjust(table["pvalue"], by=True)
        finite = np.isfinite(table["pvalue"])
        finite_q = np.full(n, np.nan); finite_q[finite] = adjust(table["pvalue"][finite])
        table["pydeseq2_finite_family_padj"] = finite_q
        half_width = special.ndtri(.975) * table["lfcSE"]
        table["ci95_lower_log2fc"] = table["log2FoldChange"] - half_width
        table["ci95_upper_log2fc"] = table["log2FoldChange"] + half_width
    return {"contrasts": result, "covariance_treatment_interaction": covariance_ti, "means": means}


def count_boundary_support(counts: np.ndarray, libraries: pd.DataFrame) -> tuple[dict, dict]:
    """Gene-major raw-count masks; no fitted values enter this decision."""
    sums, masks = {}, {}
    for group in ["nonPSC", "PSC"]:
        members = libraries.disease.to_numpy() == group
        for treatment in ["CNT", "IL17"]:
            selection = members & (libraries.treatment.to_numpy() == treatment)
            sums["count_sum_" + group + "_" + treatment] = counts[:, selection].sum(axis=1)
        pair_totals = np.column_stack([counts[:, libraries.donor_id.to_numpy() == donor].sum(axis=1)
                                       for donor in sorted(DONORS) if DONORS[donor] == group])
        sums["n_zero_total_donor_pairs_" + group] = np.count_nonzero(pair_totals == 0, axis=1)
    for label, groups in [("interaction", ["nonPSC", "PSC"]), ("nonPSC_treatment", ["nonPSC"]), ("PSC_treatment", ["PSC"])]:
        arms = [group + "_" + treatment for group in groups for treatment in ["CNT", "IL17"]]
        absent = np.column_stack([sums["count_sum_" + arm] == 0 for arm in arms])
        masks[label] = {"supported": ~absent.any(axis=1),
                        "zero_required_count_arms": ["|".join(arm for arm, zero in zip(arms, row) if zero) for row in absent]}
    return sums, masks


def apply_boundary_support(contrasts: dict, masks: dict) -> dict:
    """Preserve reconstructed native diagnostics, mask ordinary Wald inference."""
    result = {}
    for label, original in contrasts.items():
        table = {key: value.copy() for key, value in original.items()}
        for key in ["log2FoldChange", "lfcSE", "stat", "pvalue"]:
            table["pydeseq2_unmasked_" + key] = table[key].copy()
            table[key][~masks[label]["supported"]] = np.nan
        table["pydeseq2_unmasked_finite_family_padj"] = table.pop("pydeseq2_finite_family_padj")
        for by, key in [(False, "q_bh"), (True, "q_by")]:
            table[key] = adjust(table["pvalue"], by=by)
        margin = special.ndtri(.975) * table["lfcSE"]
        table["ci95_lower_log2fc"] = table["log2FoldChange"] - margin
        table["ci95_upper_log2fc"] = table["log2FoldChange"] + margin
        result[label] = table
    return result


def refit_indices(counts: np.ndarray, ids: list[str], eligible: np.ndarray) -> list[int]:
    candidates = np.flatnonzero(eligible & np.all(counts >= 10, axis=1)).tolist()
    return sorted(candidates, key=lambda i: (hashlib.sha256(("PSC-organoid-conditional-NB-v1:" + ids[i]).encode()).hexdigest(), ids[i]))[:REFIT_N]


def nb_objective(beta: np.ndarray, y: np.ndarray, factors: np.ndarray, x: np.ndarray,
                 alpha: float) -> tuple[float, np.ndarray, np.ndarray]:
    """Conditional NB NLL up to beta-independent constants, with the pinned ridge."""
    r = 1 / alpha
    eta = np.log(factors) + x @ beta - np.log(r)
    p = special.expit(eta)
    loss = float(np.sum(r * np.logaddexp(0, eta) + y * np.logaddexp(0, -eta)) + RIDGE * np.dot(beta, beta) / 2)
    score = x.T @ ((y+r) * p-y) + RIDGE * beta
    observed_weights = (y+r) * p * special.expit(-eta)
    hessian = x.T @ (observed_weights[:, None] * x) + RIDGE * np.eye(x.shape[1])
    return loss, score, hessian


def optimize_conditional(y: np.ndarray, factors: np.ndarray, x: np.ndarray, alpha: float) -> dict:
    require(y.shape == factors.shape == (len(x),) and np.isfinite(y).all() and (y >= 0).all()
            and np.isfinite(factors).all() and (factors > 0).all() and np.isfinite(alpha) and alpha > 0,
            "Invalid conditional NB inputs")
    start = np.linalg.lstsq(x, np.log((y + .5) / factors), rcond=None)[0]
    objective = lambda b: nb_objective(b, y, factors, x, alpha)
    fitted = optimize.minimize(lambda b: objective(b)[0], start, jac=lambda b: objective(b)[1],
                              hess=lambda b: objective(b)[2], method="trust-exact",
                              options={"gtol": 1e-8, "maxiter": 1000})
    loss, score, hessian = objective(fitted.x)
    # The optimizer can flag finite precision at an otherwise stationary solution.
    decrement = float(np.sqrt(max(0., score @ np.linalg.solve(hessian, score))))
    require(np.isfinite(fitted.x).all() and decrement <= 1e-5, "Independent conditional optimizer not stationary")
    return {"beta": fitted.x, "objective": loss, "newton_decrement": decrement,
            "optimizer_success": bool(fitted.success), "iterations": int(fitted.nit)}


def validate_frozen_plan(plan: dict, repo: Path = REPO) -> tuple[dict, dict[str, Path], list[Path]]:
    require(plan.get("status") == "frozen", "Real verification requires a frozen plan")
    frozen = utc(plan.get("frozen_at_utc"))
    arm = plan.get("organoid", {})
    require(arm.get("series") == "GSE239283" and arm.get("population") == "all_published_retained_organoid_cells", "Wrong population")
    require(arm.get("donor_disease") == DONORS and arm.get("expected_library_count") == 16
            and arm.get("expected_cell_count") == 46343, "Wrong biological-unit freeze")
    require(type(arm.get("expected_feature_count")) is int and arm["expected_feature_count"] > 0, "Missing feature-universe size")
    require(bool(arm.get("feature_namespace")) and bool(arm.get("feature_universe_description")), "Missing feature provenance")
    require(type(arm.get("n_cpus")) is int and arm["n_cpus"] > 0, "Invalid n_cpus")
    require(arm.get("ready_for_inference") is True and arm.get("raw_count_gate_passed") is True
            and arm.get("count_quantity") == "raw_RNA_counts", "Unqualified raw-RNA gate")
    for name, expected in [("filter", FILTER), ("model", MODEL), ("sensitivity", SENSITIVITY)]:
        actual = arm.get(name, {})
        require(all(k in actual and type(actual[k]) is type(v) and actual[k] == v for k, v in expected.items()),
                f"Wrong frozen {name} settings")
    software = arm.get("software", {})
    require(SOFTWARE.issubset(software) and software["pydeseq2"] == "0.5.4", "Unpinned software")
    for name, version in software.items():
        actual = ".".join(map(str, sys.version_info[:3])) if name == "python" else importlib.metadata.version(name)
        require(actual == version, f"Software mismatch: {name}")
    paths = {key: check_pin(arm.get(key, {}), repo) for key in ["counts", "features", "libraries", "count_quantity_qualification"]}
    others = []
    for kind in ["source_artifacts", "code_artifacts"]:
        require(isinstance(arm.get(kind), list) and bool(arm[kind]), f"No {kind}")
        these = []
        for spec in arm[kind]:
            if kind == "source_artifacts":
                require(bool(spec.get("url")) and utc(spec.get("retrieved_at_utc")) <= frozen, "Source acquired after freeze/no URL")
            these.append(check_pin(spec, repo))
        require(len(set(these)) == len(these), "Duplicate artifact pins")
        others.extend(these)
    code_paths = [check_pin(spec, repo) for spec in arm["code_artifacts"]]
    require(Path(__file__).resolve() in code_paths, "Verifier itself must be frozen before effects")
    require((repo / "scripts/analyze_psc_organoid_il17.py").resolve() in code_paths, "Primary script not pinned")
    return arm, paths, others


def public_axis(frame: pd.DataFrame, features: pd.DataFrame, eligible: np.ndarray, diagnostics: dict,
                audit: Audit, label: str, id_column: str) -> None:
    require(id_column in frame, "Public feature identifier missing")
    audit.equal(label + ":all_source_IDs", frame[id_column], features[id_column])
    for column in features.columns:
        require(column in frame, "Source annotation dropped")
        audit.equal(label + ":annotation:" + column, frame[column], features[column])
    audit.equal(label + ":eligibility", frame.eligible, np.where(eligible, "True", "False"))
    for column, expected in diagnostics.items():
        if expected.dtype.kind == "b":
            audit.equal(label + ":" + column, frame[column], np.where(expected, "True", "False"))
        else:
            audit.close(label + ":" + column, numeric(frame, [column]).ravel(), expected)
    if "status" in frame:
        audit.equal(label + ":excluded_status", frame.loc[~eligible, "status"], ["excluded_fixed_count_filter"] * int((~eligible).sum()))
    require(not ({"donor_id", "omitted_donor_id", "PSC_label_donors_json"} & set(frame.columns)), "Private donor axis in public table")


def check_numeric_fields(frame: pd.DataFrame, eligible: np.ndarray, expected: dict,
                         audit: Audit, label: str) -> None:
    for column, values in expected.items():
        actual = numeric(frame, [column]).ravel()
        compare = audit.p if column in {"pvalue", "q_bh", "q_by", "q_bh_sensitivity", "q_by_sensitivity"} or column.endswith("_pvalue") or column.endswith("_padj") else audit.close
        compare(label + ":" + column, actual[eligible], values)
        require(np.isnan(actual[~eligible]).all(), f"Excluded gene has numerical result: {label}:{column}")


def verify_run(plan_path: Path, plan_sha256: str, output_dir: Path, work_dir: Path,
               structure_only: bool = False, repo: Path = REPO) -> dict:
    """Read-only verifier. The caller alone decides whether to save the receipt."""
    plan_path = plan_path.resolve(); output_dir = output_dir.resolve(); work_dir = work_dir.resolve()
    require(bool(re.fullmatch(r"[0-9a-f]{64}", str(plan_sha256))) and digest(plan_path) == plan_sha256, "Plan hash mismatch/missing")
    plan = json.loads(plan_path.read_text()); arm, inputs, other_pins = validate_frozen_plan(plan, repo)
    require(work_dir.is_relative_to(repo / "work"), "Native donor outputs must stay in repository work/")
    require(not output_dir.is_relative_to(repo / "work") and not output_dir.is_relative_to(repo / "data/raw"), "Public/native input paths overlap")
    require(not work_dir.is_relative_to(output_dir) and not output_dir.is_relative_to(work_dir), "Work/public paths overlap")
    require(all(not p.is_relative_to(output_dir) and not p.is_relative_to(work_dir) for p in list(inputs.values()) + other_pins + [plan_path]), "Run directories contain frozen inputs")
    audit = Audit()
    id_column = arm.get("feature_id_column", "source_feature_id")
    raw = read_table(inputs["counts"], arm["counts"].get("separator", "\t"))
    require(raw.columns[0] == id_column and raw.shape[1] == 17, "Counts are not feature x 16 library")
    ids = exact_names(raw[id_column], "source feature IDs"); libraries_ids = exact_names(raw.columns[1:], "raw libraries")
    require(len(ids) == arm["expected_feature_count"], "Incorrect released feature count")
    counts = integer_tokens(raw.iloc[:, 1:].to_numpy())
    totals = np.asarray([sum(map(int, counts[:, j])) for j in range(16)])
    require((totals > 0).all() and (totals <= 2**53-1).all(), "Invalid exact library total")
    features = indexed(inputs["features"], id_column, ids, arm["features"].get("separator", "\t"), permit_reorder=True)
    libraries = indexed(inputs["libraries"], "library_id", libraries_ids, arm["libraries"].get("separator", "\t"), permit_reorder=True)
    x, columns, donors, leverage = make_design(libraries)
    audit.equal("source:library_totals", integer_tokens(libraries.total_released_counts.to_numpy(), positive=True), totals)
    audit.equal("source:retained_cells", integer_tokens(libraries.n_cells.to_numpy(), positive=True).sum(), 46343)
    eligible, diagnostics, factors, basis_size = filter_and_ratio(counts)
    eligible_ids = [g for g, yes in zip(ids, eligible) if yes]
    # This selection deliberately occurs before any fitted effects are read.
    selected_for_refit = refit_indices(counts, ids, eligible)
    count_support, support_masks = count_boundary_support(counts[eligible], libraries)
    feature_public = read_table(output_dir / "organoid-feature-universe.csv.gz")
    public_axis(feature_public, features, eligible, diagnostics, audit, "universe", id_column)
    audit.equal("frozen_work_plan", digest(work_dir / "execution-plan.json"), plan_sha256)
    design_file = indexed(work_dir / "design-matrix.csv", "library_id", libraries_ids)
    audit.equal("design:column_names", list(design_file.columns), ["library_id"] + columns)
    audit.equal("design:entries", numeric(design_file, columns), x)
    selected_file = indexed(work_dir / "selected-libraries.csv", "library_id", libraries_ids)
    for column in libraries.columns:
        if column == "n_cells":
            audit.equal("selected:" + column, integer_tokens(selected_file[column].to_numpy()), integer_tokens(libraries[column].to_numpy()))
        else:
            audit.equal("selected:" + column, selected_file[column], libraries[column])
    summary = {"schema_version": 1, "status": "structural_verification_passed", "plan_sha256": plan_sha256,
               "verifier_sha256": digest(Path(__file__).resolve()), "verification_time_utc": datetime.now(timezone.utc).isoformat(),
               "independent_analysis_helpers_imported": False, "full_dispersion_model_refitted": False,
               "counts": {"features": len(ids), "eligible": len(eligible_ids), "ratio_basis": basis_size,
                          "donors": 8, "libraries": 16, "retained_cells": 46343, "unfiltered_RNA_total": int(sum(totals))},
               "design": {"rank": 10, "residual_df": 6, "leverage": .625, "identical_row_replicates": 1},
               "conditional_refit_selection": {"rule": REFIT_SELECTION, "source_feature_ids": [ids[i] for i in selected_for_refit]},
               "limits": ["raw_assay_provenance_relies_on_reviewed_pinned_qualification_not_integer_values",
                          "no_recovery_of_historical_RNA_metadata_feature_universe",
                          "conditional_Wald_reconstruction_not_independent_dispersion_estimation",
                          "no_small_sample_calibration_or_observational_exchangeability_proof",
                          "eight_donors_not_cells_libraries_genes_or_omissions_are_replicates",
                          "not_a_causal_or_clinical_validation"]}
    if structure_only:
        summary["checks"] = audit.checks
        return summary
    manifest = json.loads((output_dir / "organoid-output-manifest.json").read_text())
    audit.equal("manifest:plan", manifest["plan_sha256"], plan_sha256)
    audit.equal("manifest:public_set", sorted(manifest["outputs"]), sorted(PUBLIC_FILES))
    audit.equal("manifest:privacy", manifest["donor_level_outputs"], "ignored_work_directory_only")
    audit.equal("manifest:no_mutation", manifest["source_values_mutated"], False)
    for name, pin in manifest["outputs"].items():
        check_pin({**pin, "path": str(output_dir / name)}, repo)
    qc = json.loads((output_dir / "organoid-qc.json").read_text())
    for key, expected in [("status", "complete"), ("effects_computed", True), ("plan_sha256", plan_sha256),
                          ("donors", 8), ("libraries", 16), ("released_cells", 46343),
                          ("released_features", len(ids)), ("eligible_features", len(eligible_ids)),
                          ("ratio_basis_positive_in_all_16_libraries", basis_size), ("raw_count_gate_passed", True),
                          ("count_quantity", "raw_RNA_counts"), ("donor_values_in_public_outputs", False)]:
        audit.equal("qc:" + key, qc[key], expected)
    primary_script = (repo / "scripts/analyze_psc_organoid_il17.py").resolve()
    audit.equal("qc:script", qc["script_sha256"], digest(primary_script))
    audit.equal("manifest:script", manifest["script_sha256"], qc["script_sha256"])
    for qkey, akey in [("software", "software"), ("source_pins", "source_artifacts"), ("code_pins", "code_artifacts"),
                       ("filter", "filter"), ("model", "model"), ("sensitivity", "sensitivity")]:
        require(qc[qkey] == arm[akey], "QC no longer matches frozen " + qkey)
    require(qc["derived_input_pins"] == {key: arm[key] for key in inputs}, "QC input pins differ")
    for key, expected in [("rank", 10), ("residual_df", 6), ("n_columns", 10), ("n_libraries", 16),
                          ("n_donors", 8), ("max_identical_full_design_rows", 1), ("cooks_eligible_libraries_at_3_identical_rows", 0)]:
        audit.equal("qc:design:" + key, qc["design"][key], expected)
    audit.close("qc:leverage", [qc["design"]["leverage_min"], qc["design"]["leverage_max"]], [.625, .625])
    normalization = indexed(work_dir / "library-normalization.csv", "library_id", libraries_ids)
    for column, values in {"size_factor": factors, "design_leverage": leverage,
                           "fitted_feature_library_total": counts[eligible].sum(axis=0),
                           "unfiltered_released_library_total": totals}.items():
        audit.close("normalization:" + column, numeric(normalization, [column]).ravel(), values, rtol=1e-12, atol=1e-12)
    audit.equal("normalization:no_cook_exclusion", normalization.cooks_filter_eligible, ["False"]*16)
    native = indexed(work_dir / "NB-all-native-feature-estimates-and-flags.csv.gz", "source_feature_id", eligible_ids)
    beta_frame = indexed(work_dir / "NB-natural-log-coefficients.csv.gz", "source_feature_id", eligible_ids)
    audit.equal("NB:coefficient_names", list(beta_frame.columns), ["source_feature_id"] + columns)
    beta = numeric(beta_frame, columns); alpha = numeric(native, ["dispersions"]).ravel()
    nb = nb_reconstruction(x, factors, counts[eligible], beta, alpha)
    nb["contrasts"] = apply_boundary_support(nb["contrasts"], support_masks)
    public = {}
    for key, name in [("interaction", "organoid-primary.csv.gz"), ("welch", "organoid-welch.csv.gz"),
                      ("leaveout", "organoid-leave-one-donor-out.csv.gz"), ("allocations", "organoid-allocations.csv.gz")]:
        public[key] = read_table(output_dir / name)
        public_axis(public[key], features, eligible, diagnostics, audit, key, id_column)
    within = read_table(output_dir / "organoid-within-group.csv.gz")
    audit.equal("within:contrast_order", within.contrast, ["nonPSC_treatment"]*len(ids) + ["PSC_treatment"]*len(ids))
    for label in ["nonPSC_treatment", "PSC_treatment"]:
        public[label] = within.loc[within.contrast == label].copy()
        public_axis(public[label], features, eligible, diagnostics, audit, label, id_column)
    for label, expected in nb["contrasts"].items():
        check_numeric_fields(public[label], eligible, expected, audit, label)
        status = np.where(np.isfinite(expected["pvalue"]), "tested_NB_Wald_small_n_calibration_limited", "fit_pvalue_unavailable")
        status = np.where(support_masks[label]["supported"], status, "unavailable_boundary_or_unidentified_group_log_response")
        audit.equal(label + ":test_status", public[label].loc[eligible, "status"], status)
        audit.equal(label + ":count_support", public[label].loc[eligible, "count_supported_contrast"], np.where(support_masks[label]["supported"], "True", "False"))
        audit.equal(label + ":zero_required_arms", public[label].loc[eligible, "zero_required_count_arms"], support_masks[label]["zero_required_count_arms"])
    primary = public["interaction"]
    check_numeric_fields(primary, eligible, count_support, audit, "count_arm_support")
    check_numeric_fields(primary, eligible, {"treatment_interaction_covariance_log2": nb["covariance_treatment_interaction"]}, audit, "NB")
    for qtype in ["bh", "by"]:
        selected = np.zeros(len(ids), dtype=bool)
        selected[eligible] = np.isfinite(nb["contrasts"]["interaction"]["pvalue"]) & (nb["contrasts"]["interaction"]["q_" + qtype] <= .05)
        audit.equal("NB:significance:" + qtype, primary["interaction_q_" + qtype + "_le_0_05"], np.where(selected, "True", "False"))
    fitted_qc = qc["negative_binomial"]
    require(fitted_qc["dispersion_fit_type"] in {"parametric", "mean"}, "Unapproved dispersion fallback")
    audit.equal("dispersion:fallback_flag", fitted_qc["dispersion_fallback_used"], fitted_qc["dispersion_fit_type"] == "mean")
    for key, expected in [("normalization_basis_features", basis_size), ("fixed_family_size", len(eligible_ids)),
                          ("automatic_cooks_exclusions", 0), ("cooks_eligible_libraries", 0)]:
        audit.equal("NB:qc:" + key, fitted_qc[key], expected)
    audit.equal("NB:dispersion_type", primary.loc[eligible, "dispersion_fit_type"], [fitted_qc["dispersion_fit_type"]]*len(eligible_ids))
    audit.equal("NB:no_automatic_exclusion", primary.loc[eligible, "cooks_automatic_exclusion"], ["False"]*len(eligible_ids))
    for column in ["genewise_dispersions", "fitted_dispersions", "dispersions", "_genewise_converged", "_MAP_converged", "_LFC_converged", "_outlier_genes", "non_zero"]:
        if column in native:
            audit.equal("native_flag:" + column, primary.loc[eligible, column.lstrip("_")], native[column])
    for column in ["genewise_converged", "MAP_converged", "LFC_converged"]:
        if column in primary:
            audit.equal("NB:nonconvergence:" + column, fitted_qc[column + "_false"], int((primary.loc[eligible, column] == "False").sum()))
    cooks_file = indexed(work_dir / "all-cooks-distances.csv.gz", "library_id", libraries_ids)
    audit.equal("Cook:feature_axis", list(cooks_file.columns), ["library_id"] + eligible_ids)
    cooks = numeric(cooks_file, eligible_ids)
    cutoff = stats.f.ppf(.99, 10, 6)
    audit.close("Cook:diagnostic_cutoff", fitted_qc["cooks_cutoff_diagnostic_only"], cutoff)
    check_numeric_fields(primary, eligible, {"max_cooks_distance": cooks.max(axis=0),
                         "n_libraries_cooks_above_diagnostic_cutoff": (cooks > cutoff).sum(axis=0),
                         "n_nonfinite_cooks": (~np.isfinite(cooks)).sum(axis=0)}, audit, "Cook")
    for label, values in nb["contrasts"].items():
        audit.equal("NB:finite_family:" + label, fitted_qc["finite_pvalues"][label], int(np.isfinite(values["pvalue"]).sum()))
        audit.equal("NB:native_finite_family:" + label, fitted_qc["native_unmasked_finite_pvalues"][label], int(np.isfinite(values["pydeseq2_unmasked_pvalue"]).sum()))
        audit.equal("NB:boundary_unavailable:" + label, fitted_qc["count_boundary_unavailable"][label], int((~support_masks[label]["supported"]).sum()))
    changes = changes_from_counts(counts, libraries, donors, eligible)
    psc = np.asarray([DONORS[d] == "PSC" for d in donors])
    sensitivity = donor_sensitivity(changes, psc)
    changes_file = indexed(work_dir / "donor-paired-log2cpm-changes.csv.gz", "donor_id", donors)
    audit.equal("donor_changes:axis", list(changes_file.columns), ["donor_id"] + eligible_ids + ["disease"])
    audit.equal("donor_changes:disease", changes_file.disease, [DONORS[d] for d in donors])
    audit.close("donor_changes:all_values", numeric(changes_file, eligible_ids), changes.T)
    check_numeric_fields(public["welch"], eligible, sensitivity["welch"], audit, "Welch")
    audit.equal("Welch:status", public["welch"].loc[eligible, "status"], sensitivity["welch_status"])
    audit.equal("Welch:finite_count", qc["welch_finite_pvalues"], int(np.isfinite(sensitivity["welch"]["pvalue"]).sum()))
    check_numeric_fields(public["leaveout"], eligible, sensitivity["leaveout_summary"], audit, "leaveout")
    audit.equal("leaveout:status", public["leaveout"].loc[eligible, "status"], ["diagnostic_only_shared_donors_not_replications_or_intervals"]*len(eligible_ids))
    leave_file = indexed(work_dir / "all-donor-leaveout-effects.csv.gz", "omitted_donor_id", donors)
    audit.equal("leaveout:feature_axis", list(leave_file.columns), ["omitted_donor_id"] + eligible_ids)
    audit.close("leaveout:all_values", numeric(leave_file, eligible_ids), sensitivity["leaveout"])
    allocation_ids = list(map(str, range(70)))
    allocation_file = indexed(work_dir / "all-70-allocation-effects.csv.gz", "allocation_id", allocation_ids)
    audit.equal("allocation:feature_axis", list(allocation_file.columns), ["allocation_id"] + eligible_ids)
    audit.close("allocation:all_values", numeric(allocation_file, eligible_ids), sensitivity["allocations"])
    check_numeric_fields(public["allocations"], eligible, sensitivity["allocation_summary"], audit, "allocation")
    audit.equal("allocation:status", public["allocations"].loc[eligible, "status"], ["observational_exchangeability_sensitivity_not_discovery_P"]*len(eligible_ids))
    labels = indexed(work_dir / "all-70-allocation-labels.csv", "allocation_id", allocation_ids)
    for index, chosen in enumerate(sensitivity["assignment_indices"]):
        audit.equal("allocation:labels:" + str(index), json.loads(labels.iloc[index].PSC_label_donors_json), [donors[j] for j in chosen])
        observed = tuple(np.flatnonzero(psc)) == chosen
        audit.equal("allocation:observed:" + str(index), labels.iloc[index].is_observed, str(observed))
    audit.equal("allocation:qc_count", qc["allocation_count"], 70)
    audit.close("allocation:resolution", qc["allocation_minimum_possible_two_sided_fraction"], 2/70)
    warnings = json.loads((work_dir / "fit-warnings.json").read_text())
    audit.equal("NB:warning_count", qc["fit_warning_count"], len(warnings))
    conditional = []
    lookup = {source_index: eligible_index for eligible_index, source_index in enumerate(np.flatnonzero(eligible))}
    for source_index in selected_for_refit:
        i = lookup[source_index]
        row = {"source_feature_id": ids[source_index]}
        if not (np.isfinite(beta[i]).all() and np.isfinite(alpha[i]) and alpha[i] > 0):
            row["status"] = "unavailable_native_fit_no_conditional_comparison"
        elif nb["means"][i].min() <= .5 or np.abs(beta[i]).max() >= 30:
            row["status"] = "outside_unclamped_interior_conditional_scope"
        else:
            refit = optimize_conditional(counts[source_index].astype(float), factors, x, float(alpha[i]))
            objective, _, _ = nb_objective(beta[i], counts[source_index], factors, x, float(alpha[i]))
            native_gap = float(objective - refit["objective"])
            standardized = max(abs((beta[i, -1] - refit["beta"][-1]) / (nb["contrasts"]["interaction"]["lfcSE"][i] * np.log(2))),
                               abs((beta[i, -2] - refit["beta"][-2]) / (nb["contrasts"]["nonPSC_treatment"]["lfcSE"][i] * np.log(2))),
                               abs((beta[i, -2:].sum() - refit["beta"][-2:].sum()) / (nb["contrasts"]["PSC_treatment"]["lfcSE"][i] * np.log(2))))
            require(native_gap >= -1e-6 and 2*native_gap <= 1e-3 and standardized <= .02,
                    "Conditional NB stationary-fit disagreement on prespecified dense gene")
            row.update({"status": "conditional_optimization_agrees", "twice_objective_improvement": 2*native_gap,
                        "maximum_contrast_difference_in_native_SE": float(standardized),
                        "maximum_natural_log_coefficient_difference": float(np.max(np.abs(beta[i] - refit["beta"]))),
                        "optimizer_success": refit["optimizer_success"], "newton_decrement": refit["newton_decrement"]})
        conditional.append(row)
    # Catch mutation during verification as well as mismatched input-at-start pins.
    for spec in [arm[key] for key in inputs] + arm["source_artifacts"] + arm["code_artifacts"]:
        check_pin(spec, repo)
    audit.equal("plan:unchanged_after_verification", digest(plan_path), plan_sha256)
    receipt_inputs = {name: {"sha256": digest(work_dir / name), "bytes": (work_dir / name).stat().st_size} for name in PRIVATE_FILES}
    summary.update({"status": "independent_numerical_verification_passed", "checks": audit.checks,
                    "dispersion_fit_type": fitted_qc["dispersion_fit_type"],
                    "conditional_fixed_dispersion_NB": {"selection_before_effects": True, "ridge": RIDGE,
                        "criteria": {"maximum_twice_objective_improvement": .001, "maximum_contrast_difference_in_native_SE": .02},
                        "results": conditional, "dispersion_estimation_repeated": False},
                    "public_manifest_sha256": digest(output_dir / "organoid-output-manifest.json"),
                    "private_array_receipts": receipt_inputs})
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True, help="Existing primary public output directory")
    parser.add_argument("--work-dir", type=Path, required=True, help="Existing native private output directory")
    parser.add_argument("--structure-only", action="store_true", help="No effects read or computed; still requires a frozen plan")
    parser.add_argument("--verification-json", type=Path, help="Optional new receipt, only when requested; never overwrite")
    args = parser.parse_args()
    if args.verification_json:
        destination = args.verification_json.resolve()
        require(not args.verification_json.exists() and not args.verification_json.is_symlink(), "Receipt destination already exists")
        protected = [args.output_dir.resolve(), args.work_dir.resolve(), REPO / "data/raw", REPO / "scripts", REPO / "tests", REPO / "config", REPO / "docs"]
        require(not any(destination.is_relative_to(p) for p in protected) and destination != args.plan.resolve(), "Receipt destination overlaps protected inputs/results")
    report = verify_run(args.plan, args.plan_sha256, args.output_dir, args.work_dir, args.structure_only)
    if args.verification_json:
        with args.verification_json.open("x", encoding="utf-8") as stream:
            json.dump(report, stream, sort_keys=True, indent=2, allow_nan=False)
            stream.write("\n")
    print(json.dumps({"status": report["status"], "checks": len(report["checks"]), "counts": report["counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
