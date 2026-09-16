#!/usr/bin/env python3
"""Frozen GSE239283 donor-paired IL-17A analysis. No downloads or source mutation.

The experimental units are eight donors, not 46,343 retained organoid cells.
Primary: NB log2 ratio-of-ratios, PSC response minus procedure-control response.
Sensitivity: equal-donor mean difference of paired log2(CPM+1) changes.
Library/donor-level values, all Cook distances and fits stay under ignored work/.
Public tables retain every released feature, including untested/unavailable rows.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import itertools
import json
import re
import sys
import warnings
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution, norm, t as t_distribution

REPO = Path(__file__).resolve().parents[1]
DONOR_DISEASE = {**{f"PSC_{n}": "PSC" for n in (4, 6, 8, 9)},
                 **{f"nonPSC_{n}": "nonPSC" for n in (2, 3, 4, 5)}}
FILTER_SETTINGS = {"min_count": 10, "min_libraries": 4}
MODEL_SETTINGS = {
    "design": "Intercept+donor+treatment_IL17+disease_PSC:treatment_IL17",
    "primary_contrast": "PSC_response_minus_nonPSC_response",
    "size_factors_fit_type": "ratio", "fit_type": "parametric",
    "dispersion_fallback": "builtin_mean_if_parametric_fails",
    "refit_cooks": False, "cooks_filter": False, "independent_filter": False,
    "lfc_shrinkage": False, "alpha": 0.05,
    "multiplicity": "BH_all_fixed_eligible_features_NA_as_1_correction_only",
    "dependence_sensitivity": "BY_same_family",
    "count_boundary_rule": "group_response_needs_both_arms_positive_interaction_needs_all4_arms_positive",
    "zero_total_donor_pairs": "flag_only_no_new_feature_filter",
}
SENSITIVITY_SETTINGS = {
    "transformation": "log2(CPM+1)", "library_size": "all_unfiltered_released_features",
    "unit": "one_IL17_minus_CNT_change_per_donor", "test": "Welch_t_with_Satterthwaite_df",
    "leaveout": "one_donor_out_effect_and_sign_only_fixed_features",
    "allocation": "all_70_4v4_labels_absolute_mean_change_difference_descriptive_only",
    "allocation_tie_atol": 1e-12, "allocation_tie_rtol": 1e-12,
}
REQUIRED_SOFTWARE = {"python", "numpy", "pandas", "scipy", "pydeseq2", "anndata", "formulaic", "joblib"}


class AnalysisError(ValueError):
    """A frozen source, design or analysis contract was violated."""


@dataclass
class CountData:
    counts: pd.DataFrame  # libraries x unchanged source-feature IDs
    features: pd.DataFrame  # original columns retained, indexed by exact source ID


@dataclass
class NBFit:
    results: dict[str, pd.DataFrame]
    feature_diagnostics: pd.DataFrame
    library_diagnostics: pd.DataFrame
    cooks: pd.DataFrame
    coefficients: pd.DataFrame
    feature_fit_state: pd.DataFrame  # every native per-feature estimate/flag, unrenamed
    diagnostics: dict[str, Any]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_pin(spec: dict[str, Any], base: Path = REPO) -> Path:
    expected = spec.get("sha256", "")
    if not isinstance(expected, str) or not re.fullmatch(r"[a-f0-9]{64}", expected):
        raise AnalysisError("Every artifact needs an exact lowercase SHA-256 pin")
    path = Path(spec.get("path", ""))
    path = path if path.is_absolute() else base / path
    if not path.is_file() or path.name.endswith(".prefix"):
        raise AnalysisError(f"Pinned complete artifact is missing: {path}")
    if sha256(path) != expected:
        raise AnalysisError(f"SHA-256 mismatch: {path}")
    if "bytes" in spec and path.stat().st_size != spec["bytes"]:
        raise AnalysisError(f"Byte-length mismatch: {path}")
    return path.resolve()


def _unique(values: Any, label: str) -> list[str]:
    result = list(values)
    if any(not isinstance(v, str) or not v or v != v.strip() for v in result):
        raise AnalysisError(f"Missing, non-string or padded {label}")
    if len(set(result)) != len(result):
        raise AnalysisError(f"Duplicate {label}; resolve identity before joining")
    return result


def read_table(path: Path, sep: str = "\t") -> pd.DataFrame:
    """Literal strings including NA; reject duplicate physical headers/row widths."""
    if sep not in (",", "\t"):
        raise AnalysisError("Pinned table separator must be comma or tab")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter=sep, strict=True)
        header = next(reader, [])
        _unique(header, "physical table columns")
        if not header:
            raise AnalysisError("Empty table")
        rows = []
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise AnalysisError(f"Incorrect table row width at line {line}")
            rows.append(row)
    if not rows:
        raise AnalysisError("Table has no rows")
    return pd.DataFrame(rows, columns=header)


def validate_counts(values: Any) -> np.ndarray:
    """Validate textual integers before conversion; never round a source token."""
    array = np.asarray(values)
    if array.ndim != 2 or not all(array.shape) or array.dtype.kind == "b":
        raise AnalysisError("Counts must be a nonempty numeric libraries-by-features matrix")
    maximum = 2**53 - 1
    if array.dtype.kind in {"O", "U", "S"}:
        exact = np.empty(array.shape, dtype=np.int64)
        for index, value in np.ndenumerate(array):
            if isinstance(value, str):
                # Source pseudobulk TSVs must use integer tokens, not decimals or
                # exponent notation that could round or underflow through float64.
                if not re.fullmatch(r"[0-9]+", value):
                    raise AnalysisError("Text count tokens must be exact nonnegative integers; never round/impute")
                digits = value.lstrip("0") or "0"
                if len(digits) > 16 or (len(digits) == 16 and digits > str(maximum)):
                    raise AnalysisError("Counts exceed the exact integer validation range")
                number = int(digits)
            elif isinstance(value, (bool, np.bool_)):
                raise AnalysisError("Boolean values are not counts")
            elif isinstance(value, (int, np.integer)):
                number = int(value)
            elif isinstance(value, (float, np.floating)):
                if not np.isfinite(value) or value < 0 or value != np.floor(value):
                    raise AnalysisError("Counts must be finite nonnegative integers; never round/impute")
                number = int(value)
            else:
                raise AnalysisError("Invalid count value")
            if not 0 <= number <= maximum:
                raise AnalysisError("Count outside the finite nonnegative exact-integer range")
            exact[index] = number
        numbers = exact.astype(np.float64)
    else:
        if array.dtype.kind not in {"i", "u", "f"}:
            raise AnalysisError("Invalid count value type")
        numbers = array.astype(np.float64)
        if not np.isfinite(numbers).all() or (numbers < 0).any() or not np.equal(numbers, np.floor(numbers)).all():
            raise AnalysisError("Counts must be finite nonnegative integers; never round/impute")
        if (numbers > maximum).any():
            raise AnalysisError("Counts exceed the exact integer validation range")
        exact = numbers.astype(np.int64)
    totals = numbers.sum(axis=1)
    if (totals > maximum).any():
        raise AnalysisError("Counts exceed the exact integer validation range")
    if (totals <= 0).any():
        raise AnalysisError("A library has zero total counts")
    return exact


def read_counts(count_path: Path, feature_path: Path, count_sep: str = "\t",
                feature_sep: str = "\t", id_column: str = "source_feature_id") -> CountData:
    frame = read_table(count_path, count_sep)
    features = read_table(feature_path, feature_sep)
    if frame.columns[0] != id_column or len(frame.columns) < 2 or id_column not in features:
        raise AnalysisError("Expected explicit source feature ID and library count columns")
    ids = _unique(frame[id_column], "count feature IDs")
    _unique(features[id_column], "annotation feature IDs")
    if set(ids) != set(features[id_column]):
        raise AnalysisError("Non-exact source feature join")
    features = features.set_index(id_column, drop=False).loc[ids].copy()
    features.index.name = "feature_index"
    values = validate_counts(frame.iloc[:, 1:].to_numpy().T)
    return CountData(pd.DataFrame(values, index=frame.columns[1:], columns=ids), features)


def validate_pairs(metadata: pd.DataFrame, expected_donors: dict[str, str] | None = None,
                   expected_per_group: int | None = 4) -> None:
    required = {"library_id", "donor_id", "disease", "treatment", "n_cells"}
    if not required.issubset(metadata.columns):
        raise AnalysisError("Missing library_id/donor_id/disease/treatment/n_cells metadata")
    _unique(metadata.index, "metadata library index")
    _unique(metadata.library_id, "metadata library IDs")
    if list(metadata.index) != list(metadata.library_id):
        raise AnalysisError("Metadata index must equal library_id")
    for field in ["donor_id", "disease", "treatment"]:
        if metadata[field].isna().any() or any(not isinstance(v, str) or not v or v != v.strip() for v in metadata[field]):
            raise AnalysisError(f"Missing/invalid {field}")
    if set(metadata.disease) != {"PSC", "nonPSC"} or set(metadata.treatment) != {"CNT", "IL17"}:
        raise AnalysisError("Exactly PSC/nonPSC and CNT/IL17 are required")
    # The same exact-token guard applies to source cell counts; each must be >0.
    validate_counts(metadata[["n_cells"]].to_numpy())
    diseases = {}
    for donor, group in metadata.groupby("donor_id", sort=True):
        if len(group) != 2 or set(group.treatment) != {"CNT", "IL17"} or group.disease.nunique() != 1:
            raise AnalysisError("Each donor must have exactly one CNT and one IL17 library and one disease")
        diseases[donor] = group.disease.iloc[0]
    if expected_donors is not None and diseases != expected_donors:
        raise AnalysisError("Source donor/disease assignments differ from the frozen eight donors")
    group_counts = pd.Series(diseases).value_counts().to_dict()
    if expected_per_group is not None and group_counts != {"PSC": expected_per_group, "nonPSC": expected_per_group}:
        raise AnalysisError("Wrong independent donor counts (cells and libraries are not donors)")


def align_samples(counts: pd.DataFrame, metadata: pd.DataFrame,
                  expected_donors: dict[str, str] | None = None,
                  expected_cells: int | None = None) -> pd.DataFrame:
    _unique(counts.index, "count library IDs")
    _unique(counts.columns, "count feature IDs")
    if "library_id" not in metadata:
        raise AnalysisError("No explicit library ID column")
    _unique(metadata.library_id, "metadata library IDs")
    if set(counts.index) != set(metadata.library_id):
        raise AnalysisError("Non-exact library join; source/count columns are not interchangeable")
    selected = metadata.set_index("library_id", drop=False).loc[counts.index].copy()
    validate_pairs(selected, expected_donors)
    selected["n_cells"] = pd.to_numeric(selected.n_cells, errors="raise").astype(np.int64)
    if expected_cells is not None and int(selected.n_cells.sum()) != expected_cells:
        raise AnalysisError("Released source cell count differs from the frozen count")
    return selected


def build_design(metadata: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    validate_pairs(metadata)
    donors = sorted(metadata.donor_id.unique())
    columns = {"Intercept": np.ones(len(metadata))}
    for donor in donors[1:]:
        columns[f"donor_{donor}"] = metadata.donor_id.eq(donor).to_numpy(dtype=float)
    columns["treatment_IL17"] = metadata.treatment.eq("IL17").to_numpy(dtype=float)
    columns["disease_PSC:treatment_IL17"] = (metadata.disease.eq("PSC") & metadata.treatment.eq("IL17")).to_numpy(dtype=float)
    design = pd.DataFrame(columns, index=metadata.index)
    diagnostic = check_design(design)
    diagnostic.update({"columns": list(design.columns), "n_donors": len(donors),
        "donors_by_disease": metadata.drop_duplicates("donor_id").disease.value_counts().to_dict(),
        "libraries_by_disease_and_treatment": metadata.groupby(["disease", "treatment"]).size().unstack().to_dict(orient="index"),
        "no_disease_main_effect": True, "experimental_unit": "donor"})
    return design, diagnostic


def check_design(design: pd.DataFrame | np.ndarray) -> dict[str, Any]:
    x = np.asarray(design, dtype=float)
    if x.ndim != 2 or not all(x.shape) or not np.isfinite(x).all():
        raise AnalysisError("Invalid model design")
    n, p = x.shape
    rank = int(np.linalg.matrix_rank(x))
    if rank != p or n <= p:
        raise AnalysisError("Rank-deficient or saturated model design")
    q, _ = np.linalg.qr(x, mode="reduced")
    leverage = np.sum(q*q, axis=1)
    if (leverage >= 1 - 1e-10).any():
        raise AnalysisError("Unit-leverage sample")
    _, replicate_counts = np.unique(x, axis=0, return_counts=True)
    return {"n_libraries": n, "n_columns": p, "rank": rank, "residual_df": n-p,
        "condition_number": float(np.linalg.cond(x)), "leverage_min": float(leverage.min()),
        "leverage_max": float(leverage.max()), "max_identical_full_design_rows": int(replicate_counts.max()),
        "cooks_eligible_libraries_at_3_identical_rows": int(replicate_counts[replicate_counts >= 3].sum())}


def contrast_vectors(design: pd.DataFrame) -> dict[str, np.ndarray]:
    treatment = np.zeros(design.shape[1]); interaction = treatment.copy()
    treatment[design.columns.get_loc("treatment_IL17")] = 1
    interaction[design.columns.get_loc("disease_PSC:treatment_IL17")] = 1
    return {"interaction": interaction, "nonPSC_treatment": treatment, "PSC_treatment": treatment + interaction}


def expression_filter(counts: pd.DataFrame, min_count: int = 10, min_libraries: int = 4) -> pd.DataFrame:
    if min_count < 1 or not 1 <= min_libraries <= len(counts):
        raise AnalysisError("Invalid count filter")
    values = validate_counts(counts.to_numpy())
    above = (values >= min_count).sum(axis=0)
    return pd.DataFrame({"eligible": above >= min_libraries,
        "libraries_count_ge_threshold": above, "libraries_nonzero": (values > 0).sum(axis=0),
        "all_zero_released": (values == 0).all(axis=0), "mean_raw_count": values.mean(axis=0)}, index=counts.columns)


def fixed_adjust(pvalues: Any, method: str = "bh") -> np.ndarray:
    """Whole supplied fixed family. Only correction inputs replace unavailable P by 1."""
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or np.isinf(p).any() or ((p < 0) | (p > 1)).any():
        raise AnalysisError("Invalid P vector")
    if method not in {"bh", "by"}:
        raise AnalysisError("Only BH and BY are supported")
    if not len(p):
        return np.empty(0)
    correction = np.where(np.isnan(p), 1.0, p)
    order = np.argsort(correction, kind="stable")
    m = len(p)
    factor = np.sum(1 / np.arange(1, m + 1)) if method == "by" else 1.0
    q = np.minimum(1, np.minimum.accumulate((correction[order] * m * factor / np.arange(1, m + 1))[::-1])[::-1])
    result = np.empty_like(q); result[order] = q
    return result


def ratio_size_factors(counts: pd.DataFrame) -> tuple[np.ndarray, int]:
    values = validate_counts(counts.to_numpy())
    usable = (values > 0).all(axis=0)
    if not usable.any():
        raise AnalysisError("No all-positive eligible normalization basis; iterative fallback prohibited")
    logs = np.log(values[:, usable].astype(float))
    factors = np.exp(np.median(logs - logs.mean(axis=0), axis=1))
    if not np.isfinite(factors).all() or (factors <= 0).any():
        raise AnalysisError("Invalid ratio size factors")
    return factors, int(usable.sum())


def log2_cpm(counts: pd.DataFrame) -> np.ndarray:
    """Call on the complete UNFILTERED RELEASED universe, before eligibility subset."""
    values = validate_counts(counts.to_numpy())
    return np.log2(1 + values.astype(float) / values.sum(axis=1)[:, None] * 1_000_000)


def paired_changes(metadata: pd.DataFrame, outcomes: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    validate_pairs(metadata)
    if not metadata.index.equals(outcomes.index):
        raise AnalysisError("Paired outcome library order is not exactly aligned")
    _unique(outcomes.columns, "paired feature IDs")
    if not np.isfinite(outcomes.to_numpy(dtype=float)).all():
        raise AnalysisError("Nonfinite paired outcomes; no implicit case deletion")
    changes = []; diseases = []; donors = sorted(metadata.donor_id.unique())
    for donor in donors:
        rows = metadata.loc[metadata.donor_id.eq(donor)]
        by_treatment = rows.set_index("treatment").library_id
        changes.append(outcomes.loc[by_treatment["IL17"]].to_numpy() - outcomes.loc[by_treatment["CNT"]].to_numpy())
        diseases.append(rows.disease.iloc[0])
    index = pd.Index(donors, name="donor_id")
    return pd.DataFrame(changes, index=index, columns=outcomes.columns), pd.Series(diseases, index=index, name="disease")


def _change_groups(changes: pd.DataFrame, disease: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    _unique(changes.index, "change donor IDs"); _unique(changes.columns, "change feature IDs")
    if not changes.index.equals(disease.index) or set(disease) != {"PSC", "nonPSC"}:
        raise AnalysisError("Non-exact donor-group alignment")
    y = changes.to_numpy(dtype=float)
    if not np.isfinite(y).all():
        raise AnalysisError("Nonfinite donor changes; do not leak missing cases")
    a = y[disease.eq("PSC").to_numpy()]; b = y[disease.eq("nonPSC").to_numpy()]
    if min(len(a), len(b)) < 2:
        raise AnalysisError("At least two independent donors per group are needed")
    return a, b


def welch_changes(changes: pd.DataFrame, disease: pd.Series) -> pd.DataFrame:
    a, b = _change_groups(changes, disease)
    effect = a.mean(axis=0) - b.mean(axis=0)
    va = a.var(axis=0, ddof=1) / len(a); vb = b.var(axis=0, ddof=1) / len(b)
    variance = va + vb; se = np.sqrt(variance)
    denominator = va**2 / (len(a)-1) + vb**2 / (len(b)-1)
    valid = (denominator > 0) & (se > np.finfo(float).eps * 32 * np.maximum(1, np.abs(effect)))
    df = np.divide(variance**2, denominator, out=np.full_like(variance, np.nan), where=valid)
    stat = np.divide(effect, se, out=np.full_like(effect, np.nan), where=valid)
    p = 2 * t_distribution.sf(np.abs(stat), df)
    margin = t_distribution.ppf(0.975, df) * se
    return pd.DataFrame({"mean_change_PSC": a.mean(axis=0), "mean_change_nonPSC": b.mean(axis=0),
        "effect_log2_cpm": effect, "se": se, "t": stat, "df": df, "pvalue": p,
        "q_bh_sensitivity": fixed_adjust(p), "q_by_sensitivity": fixed_adjust(p, "by"),
        "ci95_lower": effect-margin, "ci95_upper": effect+margin,
        "n_PSC_donors": len(a), "n_nonPSC_donors": len(b),
        "status": np.where(valid, "tested_sensitivity", "zero_or_numerically_degenerate_variance_unavailable")}, index=changes.columns)


def leave_one_donor_out(changes: pd.DataFrame, disease: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    a, b = _change_groups(changes, disease)
    if (len(a), len(b)) != (4, 4):
        raise AnalysisError("Prespecified leaveouts require the full four-by-four donor panel")
    full = a.mean(axis=0) - b.mean(axis=0)
    estimates = []
    for donor in changes.index:
        keep = changes.index != donor
        pa, pb = _change_groups(changes.loc[keep], disease.loc[keep])
        estimates.append(pa.mean(axis=0) - pb.mean(axis=0))
    stack = np.asarray(estimates)
    table = pd.DataFrame(stack, index=changes.index, columns=changes.columns)
    table.index.name = "omitted_donor_id"
    summary = pd.DataFrame({"full_effect_log2_cpm": full, "minimum_leaveout_effect": stack.min(axis=0),
        "maximum_leaveout_effect": stack.max(axis=0), "n_estimable_omissions": len(changes),
        "n_same_nonzero_sign": ((np.sign(stack) == np.sign(full)) & (stack != 0) & (full != 0)).sum(axis=0),
        "n_opposite_sign": ((np.sign(stack) == -np.sign(full)) & (stack != 0) & (full != 0)).sum(axis=0),
        "n_zero_effect": (stack == 0).sum(axis=0),
        "status": "diagnostic_only_shared_donors_not_replications_or_intervals"}, index=changes.columns)
    return summary, table


def allocation_sensitivity(changes: pd.DataFrame, disease: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    a, b = _change_groups(changes, disease)
    if (len(a), len(b)) != (4, 4):
        raise AnalysisError("Exact allocation sensitivity requires four-by-four donors")
    y = changes.to_numpy(dtype=float)
    full = a.mean(axis=0) - b.mean(axis=0)
    assignments = []; differences = []
    for allocation, chosen in enumerate(itertools.combinations(range(8), 4)):
        mask = np.zeros(8, dtype=bool); mask[list(chosen)] = True
        differences.append(y[mask].mean(axis=0) - y[~mask].mean(axis=0))
        assignments.append({"allocation_id": allocation,
            "PSC_label_donors_json": json.dumps(changes.index[mask].tolist()),
            "is_observed": np.array_equal(mask, disease.eq("PSC").to_numpy())})
    stack = np.asarray(differences)
    tolerance = SENSITIVITY_SETTINGS["allocation_tie_atol"] + SENSITIVITY_SETTINGS["allocation_tie_rtol"] * np.abs(full)
    exceed = (np.abs(stack) >= np.maximum(0, np.abs(full)-tolerance)).sum(axis=0)
    result = pd.DataFrame({"observed_mean_change_difference": full, "n_allocations": 70,
        "n_as_or_more_extreme": exceed, "descriptive_tail_fraction": exceed / 70,
        "status": "observational_exchangeability_sensitivity_not_discovery_P"}, index=changes.columns)
    return result, pd.DataFrame(stack, index=pd.Index(range(70), name="allocation_id"), columns=changes.columns), pd.DataFrame(assignments)


def nb_covariance(design: Any, fitted_mean: Any, dispersion: float) -> np.ndarray:
    """Pinned PyDESeq2 0.5.4 Wald covariance H M H, natural-log coefficients.

    For a contrast c, var(c'beta)=c'Cov c. In particular PSC treatment requires
    BOTH treatment/interaction variances AND twice their covariance.
    """
    x = np.asarray(design, dtype=float); mu = np.asarray(fitted_mean, dtype=float)
    if mu.shape != (len(x),) or not np.isfinite(mu).all() or (mu <= 0).any() or not np.isfinite(dispersion) or dispersion <= 0:
        raise AnalysisError("Invalid NB fitted means or dispersion")
    weights = mu / (1 + mu * dispersion)
    information = (x.T * weights) @ x
    inverse = np.linalg.inv(information + np.eye(x.shape[1]) * 1e-6)
    return inverse @ information @ inverse


def count_arm_support(counts: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Raw-count support, not expression-selected filtering or a convergence test.

    A wholly zero arm has no finite unconstrained group log-response MLE. If both
    group arms are zero its response is unidentified. An individual zero-total
    donor pair is flagged, but does not by itself make the common group response
    unidentified when both group arms have support from the remaining donors.
    """
    if not counts.index.equals(metadata.index):
        raise AnalysisError("Count-support sample alignment mismatch")
    validate_pairs(metadata)
    values = validate_counts(counts.to_numpy())
    support = pd.DataFrame(index=counts.columns)
    for disease in ["nonPSC", "PSC"]:
        zero_pairs = np.zeros(len(counts.columns), dtype=int)
        for treatment in ["CNT", "IL17"]:
            selected = (metadata.disease.eq(disease) & metadata.treatment.eq(treatment)).to_numpy()
            support[f"count_sum_{disease}_{treatment}"] = values[selected].sum(axis=0)
        for donor in metadata.loc[metadata.disease.eq(disease), "donor_id"].unique():
            selected = metadata.donor_id.eq(donor).to_numpy()
            zero_pairs += values[selected].sum(axis=0) == 0
        support[f"n_zero_total_donor_pairs_{disease}"] = zero_pairs
    return support


def mask_count_boundary(table: pd.DataFrame, support: pd.DataFrame, contrast: str) -> pd.DataFrame:
    """Keep native values, mask only nonidentified/boundary log-response inference."""
    requirements = {"nonPSC_treatment": ["nonPSC_CNT", "nonPSC_IL17"],
        "PSC_treatment": ["PSC_CNT", "PSC_IL17"],
        "interaction": ["nonPSC_CNT", "nonPSC_IL17", "PSC_CNT", "PSC_IL17"]}
    if contrast not in requirements or not table.index.equals(support.index):
        raise AnalysisError("Boundary-mask contrast or feature alignment mismatch")
    result = table.copy()
    arms = requirements[contrast]
    positive = support[[f"count_sum_{arm}" for arm in arms]].to_numpy() > 0
    passed = positive.all(axis=1)
    # These native finite ridge/clamped results may be diagnostic but are not
    # evidence that the unconstrained response was identified by observed counts.
    for field in ["log2FoldChange", "lfcSE", "stat", "pvalue"]:
        result[f"pydeseq2_unmasked_{field}"] = result[field]
        result.loc[~passed, field] = np.nan
    result = result.rename(columns={"pydeseq2_finite_family_padj": "pydeseq2_unmasked_finite_family_padj"})
    result["count_supported_contrast"] = passed
    result["zero_required_count_arms"] = ["|".join(arm for arm, present in zip(arms, row) if not present) for row in positive]
    result.loc[~passed, "status"] = "unavailable_boundary_or_unidentified_group_log_response"
    result["q_bh"] = fixed_adjust(result.pvalue)
    result["q_by"] = fixed_adjust(result.pvalue, "by")
    result["ci95_lower_log2fc"] = result.log2FoldChange - norm.ppf(0.975) * result.lfcSE
    result["ci95_upper_log2fc"] = result.log2FoldChange + norm.ppf(0.975) * result.lfcSE
    return result


def fit_negative_binomial(counts: pd.DataFrame, metadata: pd.DataFrame,
                          design: pd.DataFrame, n_cpus: int = 1) -> NBFit:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.utils import n_or_more_replicates

    if importlib.metadata.version("pydeseq2") != "0.5.4":
        raise AnalysisError("This implementation requires PyDESeq2 0.5.4")
    if not counts.index.equals(metadata.index) or not counts.index.equals(design.index):
        raise AnalysisError("NB libraries must be exactly aligned")
    expected_design, _ = build_design(metadata)
    if not design.equals(expected_design):
        raise AnalysisError("NB design differs from the paired prespecified design")
    validate_counts(counts.to_numpy()); _unique(counts.columns, "NB feature IDs")
    expected_factors, basis = ratio_size_factors(counts)
    support = count_arm_support(counts, metadata)
    if not isinstance(n_cpus, int) or isinstance(n_cpus, bool) or n_cpus < 1:
        raise AnalysisError("Positive integer n_cpus required")
    cooks_eligible = n_or_more_replicates(design, 3).to_numpy(dtype=bool)
    if cooks_eligible.any():
        raise AnalysisError("Unexpected identical full-design row replicates")
    results = {}
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        dds = DeseqDataSet(counts=counts, metadata=metadata, design=design,
            fit_type="parametric", size_factors_fit_type="ratio", refit_cooks=False,
            n_cpus=n_cpus, quiet=True)
        if not np.array_equal(dds.obsm["design_matrix"].to_numpy(), design.to_numpy()):
            raise AnalysisError("PyDESeq2 altered the explicit design")
        dds.fit_size_factors(fit_type="ratio")
        factors = dds.obs["size_factors"].to_numpy(dtype=float)
        if not np.allclose(factors, expected_factors, rtol=1e-12, atol=1e-12):
            raise AnalysisError("Size factors differ from the pinned ratio convention")
        dds.fit_genewise_dispersions(); dds.fit_dispersion_trend()
        if dds.uns.get("disp_function_type") not in {"parametric", "mean"}:
            raise AnalysisError("Unapproved dispersion fallback")
        dds.fit_dispersion_prior(); dds.fit_MAP_dispersions(); dds.fit_LFC(); dds.calculate_cooks()
        for label, contrast in contrast_vectors(design).items():
            stats = DeseqStats(dds, contrast=contrast, cooks_filter=False,
                independent_filter=False, alpha=0.05, n_cpus=n_cpus, quiet=True)
            stats.summary()
            table = stats.results_df.copy()
            if not table.index.equals(counts.columns):
                raise AnalysisError("PyDESeq2 silently removed/reordered features")
            table = table.rename(columns={"padj": "pydeseq2_finite_family_padj"})
            table["q_bh"] = fixed_adjust(table.pvalue); table["q_by"] = fixed_adjust(table.pvalue, "by")
            table["ci95_lower_log2fc"] = table.log2FoldChange - norm.ppf(0.975) * table.lfcSE
            table["ci95_upper_log2fc"] = table.log2FoldChange + norm.ppf(0.975) * table.lfcSE
            table["status"] = np.where(table.pvalue.notna(), "tested_NB_Wald_small_n_calibration_limited", "fit_pvalue_unavailable")
            results[label] = table
    if not np.array_equal(dds.X, counts.to_numpy()):
        raise AnalysisError("Source counts changed during NB fitting")
    cook_values = np.asarray(dds.layers["cooks"])
    cutoff = float(f_distribution.ppf(0.99, design.shape[1], len(design)-design.shape[1]))
    feature_qc = pd.DataFrame({"max_cooks_distance": np.max(cook_values, axis=0),
        "n_libraries_cooks_above_diagnostic_cutoff": (cook_values > cutoff).sum(axis=0),
        "n_nonfinite_cooks": (~np.isfinite(cook_values)).sum(axis=0),
        "cooks_automatic_exclusion": False, "dispersion_fit_type": dds.uns["disp_function_type"]}, index=counts.columns)
    for column in ["genewise_dispersions", "fitted_dispersions", "dispersions", "_genewise_converged", "_MAP_converged", "_LFC_converged", "_outlier_genes", "non_zero"]:
        if column in dds.var:
            feature_qc[column.lstrip("_")] = dds.var[column].to_numpy()
    coefficients = dds.varm["LFC"].copy()
    if not isinstance(coefficients, pd.DataFrame):
        coefficients = pd.DataFrame(coefficients, index=counts.columns, columns=design.columns)
    # Validate library contrasts against the complete covariance, not sum-of-SEs.
    means = np.exp(design.to_numpy() @ coefficients.to_numpy().T) * factors[:, None]
    ti = design.columns.get_loc("treatment_IL17"); ii = design.columns.get_loc("disease_PSC:treatment_IL17")
    covariance = np.full(len(counts.columns), np.nan)
    for gene in range(len(counts.columns)):
        if not np.isfinite(means[:, gene]).all() or not (means[:, gene] > 0).all():
            continue
        dispersion = float(dds.var["dispersions"].iloc[gene])
        if not np.isfinite(dispersion) or dispersion <= 0:
            continue
        cov = nb_covariance(design, means[:, gene], dispersion) / np.log(2)**2
        covariance[gene] = cov[ti, ii]
        for label, contrast in contrast_vectors(design).items():
            se = float(np.sqrt(contrast @ cov @ contrast))
            published_se = float(results[label].lfcSE.iloc[gene])
            if np.isfinite(published_se) and not np.isclose(se, published_se, rtol=1e-8, atol=1e-10):
                raise AnalysisError("NB contrast covariance check failed")
    feature_qc["treatment_interaction_covariance_log2"] = covariance
    if not np.allclose(results["PSC_treatment"].log2FoldChange,
                       results["nonPSC_treatment"].log2FoldChange + results["interaction"].log2FoldChange,
                       rtol=1e-10, atol=1e-12, equal_nan=True):
        raise AnalysisError("Within-PSC contrast is not treatment plus interaction")
    # Apply source-count support only after validating unmodified native contrasts.
    # All eligible genes stay in every correction family; unavailable P maps to 1
    # only for adjustment. Group-supported contrasts are retained independently.
    feature_qc = feature_qc.join(support, validate="one_to_one")
    results = {label: mask_count_boundary(table, support, label) for label, table in results.items()}
    q, _ = np.linalg.qr(design.to_numpy(), mode="reduced")
    library_qc = pd.DataFrame({"library_id": counts.index, "size_factor": factors,
        "fitted_feature_library_total": counts.sum(axis=1).to_numpy(),
        "design_leverage": (q*q).sum(axis=1), "cooks_filter_eligible": cooks_eligible})
    diagnostic = {"normalization_basis_features": basis, "normalization": "exp(median(log_ratios))_no_fallback",
        "dispersion_fit_type": dds.uns["disp_function_type"], "dispersion_fallback_used": dds.uns["disp_function_type"] == "mean",
        "cooks_cutoff_diagnostic_only": cutoff, "automatic_cooks_exclusions": 0,
        "cooks_eligible_libraries": int(cooks_eligible.sum()), "fixed_family_size": len(counts.columns),
        "finite_pvalues": {label: int(table.pvalue.notna().sum()) for label, table in results.items()},
        "count_boundary_unavailable": {label: int((~table.count_supported_contrast).sum()) for label, table in results.items()},
        "native_unmasked_finite_pvalues": {label: int(table.pydeseq2_unmasked_pvalue.notna().sum()) for label, table in results.items()},
        "warnings": [{"category": w.category.__name__, "message": str(w.message)} for w in captured]}
    for column in ["genewise_converged", "MAP_converged", "LFC_converged"]:
        if column in feature_qc:
            diagnostic[column + "_false"] = int(feature_qc[column].eq(False).sum())
    return NBFit(results, feature_qc, library_qc,
        pd.DataFrame(cook_values, index=counts.index, columns=counts.columns), coefficients,
        dds.var.copy(), diagnostic)


def current_software() -> dict[str, str]:
    return {package: (".".join(map(str, sys.version_info[:3])) if package == "python"
        else importlib.metadata.version(package)) for package in sorted(REQUIRED_SOFTWARE)}


def _timestamp(value: Any, label: str) -> datetime:
    try:
        stamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (TypeError, ValueError, AttributeError) as error:
        raise AnalysisError(f"Invalid {label} timestamp") from error
    if stamp.tzinfo is None or stamp.utcoffset().total_seconds() != 0 or stamp > datetime.now(timezone.utc):
        raise AnalysisError(f"{label} must be a UTC timestamp not in the future")
    return stamp


def validate_plan(plan: dict[str, Any], for_inference: bool = True) -> dict[str, Any]:
    if plan.get("status") != "frozen":
        raise AnalysisError("An explicitly frozen plan is required")
    _timestamp(plan.get("frozen_at_utc"), "freeze")
    arm = plan.get("organoid", {})
    for label, expected in [("filter", FILTER_SETTINGS), ("model", MODEL_SETTINGS), ("sensitivity", SENSITIVITY_SETTINGS)]:
        actual = arm.get(label, {})
        if any(key not in actual or actual[key] != value for key, value in expected.items()):
            raise AnalysisError(f"Unsupported or unfrozen {label} setting")
    if arm.get("series") != "GSE239283" or arm.get("population") != "all_published_retained_organoid_cells":
        raise AnalysisError("This analysis requires GSE239283 all published retained organoid cells")
    if arm.get("donor_disease") != DONOR_DISEASE or arm.get("expected_library_count") != 16 or arm.get("expected_cell_count") != 46343:
        raise AnalysisError("Frozen eight-donor/16-library/46343-cell panel differs")
    feature_count = arm.get("expected_feature_count")
    if not isinstance(feature_count, int) or isinstance(feature_count, bool) or feature_count <= 0:
        raise AnalysisError("A qualified complete feature-universe size must be frozen")
    if not arm.get("feature_namespace") or not arm.get("feature_universe_description"):
        raise AnalysisError("Freeze the source namespace and actual released feature universe")
    if not isinstance(arm.get("n_cpus"), int) or isinstance(arm.get("n_cpus"), bool) or arm["n_cpus"] < 1:
        raise AnalysisError("Explicit positive n_cpus required")
    software = arm.get("software", {})
    if not REQUIRED_SOFTWARE.issubset(software) or software.get("pydeseq2") != "0.5.4":
        raise AnalysisError("Missing software pins, or PyDESeq2 is not version 0.5.4")
    for package, expected in software.items():
        actual = ".".join(map(str, sys.version_info[:3])) if package == "python" else importlib.metadata.version(package)
        if actual != expected:
            raise AnalysisError(f"Software pin mismatch: {package}")
    for kind in ["counts", "features", "libraries", "count_quantity_qualification"]:
        if not isinstance(arm.get(kind), dict) or not arm[kind].get("path") or not arm[kind].get("sha256"):
            raise AnalysisError(f"Missing pinned {kind} artifact")
    sources = arm.get("source_artifacts")
    if not isinstance(sources, list) or not sources:
        raise AnalysisError("Original source artifacts must be pinned separately from pseudobulk")
    for source in sources:
        if not source.get("url") or not source.get("sha256") or not source.get("path"):
            raise AnalysisError("Source URL/path/hash required")
        _timestamp(source.get("retrieved_at_utc"), "source retrieval")
    code = arm.get("code_artifacts")
    if not isinstance(code, list) or not code:
        raise AnalysisError("Pinned code artifacts required")
    if for_inference and (arm.get("ready_for_inference") is not True or arm.get("raw_count_gate_passed") is not True or arm.get("count_quantity") != "raw_RNA_counts"):
        raise AnalysisError("Inference blocked: integer values alone do not establish raw RNA counts; qualify the assay first")
    return arm


def _write_json(path: Path, document: Any) -> None:
    # Exclusive creation rejects existing symlink/hardlink destinations too.
    with path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(document, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    # Empty CSV fields mark missing numerical outputs, but not literal source IDs.
    if path.suffix == ".gz":
        import io
        with path.open("xb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                frame.to_csv(text, index=False, float_format="%.17g", na_rep="", lineterminator="\n")
    else:
        with path.open("x", encoding="utf-8", newline="") as text:
            frame.to_csv(text, index=False, float_format="%.17g", na_rep="", lineterminator="\n")


def _feature_join(features: pd.DataFrame, results: pd.DataFrame) -> pd.DataFrame:
    if not results.index.is_unique or not set(results.index).issubset(features.index):
        raise AnalysisError("Effect rows differ from the frozen feature universe")
    if set(features.columns) & set(results.columns):
        raise AnalysisError("Derived columns collide with original feature annotations")
    return features.join(results, how="left", validate="one_to_one")


def run_analysis(plan_path: Path, output_dir: Path, work_dir: Path,
                 plan_sha256: str | None = None, validate_only: bool = False) -> dict[str, Any]:
    plan_path = plan_path.resolve()
    plan_digest = sha256(plan_path)
    if plan_sha256 is not None and plan_digest != plan_sha256:
        raise AnalysisError("Execution plan SHA-256 mismatch")
    plan = json.loads(plan_path.read_text())
    arm = validate_plan(plan, for_inference=not validate_only)
    paths = {name: check_pin(arm[name]) for name in ["counts", "features", "libraries", "count_quantity_qualification"]}
    source_paths = [check_pin(spec) for spec in arm["source_artifacts"]]
    code_paths = [check_pin(spec) for spec in arm["code_artifacts"]]
    if Path(__file__).resolve() not in code_paths:
        raise AnalysisError("The executing analysis script must itself be pinned")
    if len(code_paths) != len(set(code_paths)) or len(source_paths) != len(set(source_paths)):
        raise AnalysisError("Repeated artifact pins")
    output_dir = output_dir.resolve(); work_dir = work_dir.resolve()
    if not work_dir.is_relative_to(REPO / "work"):
        raise AnalysisError("Library/donor work must stay under ignored repository work/")
    if output_dir.is_relative_to(REPO / "data/raw") or output_dir.is_relative_to(REPO / "work"):
        raise AnalysisError("Public outputs must be separate from raw/ and ignored work/")
    if output_dir == work_dir or output_dir.is_relative_to(work_dir) or work_dir.is_relative_to(output_dir):
        raise AnalysisError("Public output and donor-level work directories must be disjoint")
    pinned = set(paths.values()) | set(source_paths) | set(code_paths) | {plan_path}
    if any(path.is_relative_to(output_dir) or path.is_relative_to(work_dir) for path in pinned):
        raise AnalysisError("Output/work directories must not contain pinned input artifacts")
    if output_dir.exists() or work_dir.exists():
        raise AnalysisError("Output and work directories must both be new; never mix runs or reuse completed artifacts")
    data = read_counts(paths["counts"], paths["features"], arm["counts"].get("separator", "\t"),
        arm["features"].get("separator", "\t"), arm.get("feature_id_column", "source_feature_id"))
    metadata = read_table(paths["libraries"], arm["libraries"].get("separator", "\t"))
    selected = align_samples(data.counts, metadata, DONOR_DISEASE, arm["expected_cell_count"])
    if len(data.features) != arm["expected_feature_count"] or len(data.counts) != arm["expected_library_count"]:
        raise AnalysisError("Counts differ from the frozen complete dimensions")
    totals = data.counts.sum(axis=1).to_numpy()
    if "total_released_counts" not in selected:
        raise AnalysisError("Pseudobulk must preserve its independently recorded library totals")
    source_totals = validate_counts(selected[["total_released_counts"]].to_numpy()).ravel()
    if not np.array_equal(totals, source_totals):
        raise AnalysisError("Pseudobulk library sums disagree with the source qualification")
    design, design_qc = build_design(selected)
    filtering = expression_filter(data.counts, **FILTER_SETTINGS)
    features = _feature_join(data.features, filtering)
    reserved_result_columns = {
        "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "pydeseq2_finite_family_padj", "q_bh", "q_by",
        "ci95_lower_log2fc", "ci95_upper_log2fc", "status", "contrast", "dispersion_fit_type",
        "genewise_dispersions", "fitted_dispersions", "dispersions", "genewise_converged", "MAP_converged", "LFC_converged",
        "outlier_genes", "non_zero", "max_cooks_distance", "n_libraries_cooks_above_diagnostic_cutoff", "n_nonfinite_cooks",
        "cooks_automatic_exclusion", "treatment_interaction_covariance_log2", "interaction_q_bh_le_0_05", "interaction_q_by_le_0_05",
        "mean_change_PSC", "mean_change_nonPSC", "effect_log2_cpm", "se", "t", "df", "q_bh_sensitivity", "q_by_sensitivity",
        "ci95_lower", "ci95_upper", "n_PSC_donors", "n_nonPSC_donors", "full_effect_log2_cpm", "minimum_leaveout_effect",
        "maximum_leaveout_effect", "n_estimable_omissions", "n_same_nonzero_sign", "n_opposite_sign", "n_zero_effect",
        "observed_mean_change_difference", "n_allocations", "n_as_or_more_extreme", "descriptive_tail_fraction",
        "pydeseq2_unmasked_log2FoldChange", "pydeseq2_unmasked_lfcSE", "pydeseq2_unmasked_stat", "pydeseq2_unmasked_pvalue",
        "pydeseq2_unmasked_finite_family_padj", "count_supported_contrast", "zero_required_count_arms",
        "count_sum_nonPSC_CNT", "count_sum_nonPSC_IL17", "count_sum_PSC_CNT", "count_sum_PSC_IL17",
        "n_zero_total_donor_pairs_nonPSC", "n_zero_total_donor_pairs_PSC"}
    if set(data.features.columns) & reserved_result_columns:
        raise AnalysisError("Original annotation columns collide with derived result names; resolve namespaces before inference")
    eligible_ids = features.index[features.eligible]
    if not len(eligible_ids):
        raise AnalysisError("No features pass the fixed count filter")
    eligible_counts = data.counts.loc[:, eligible_ids]
    _, ratio_basis = ratio_size_factors(eligible_counts)
    # Exclusive directory/file creation prevents stale completed outputs and
    # preexisting destination aliases from mutating inputs or publishing donors.
    output_dir.mkdir(parents=True, exist_ok=False); work_dir.mkdir(parents=True, exist_ok=False)
    with (work_dir / "execution-plan.json").open("xb") as copy_stream:
        copy_stream.write(plan_path.read_bytes())
    _write_csv(selected.reset_index(drop=True), work_dir / "selected-libraries.csv")
    _write_csv(design.rename_axis("library_id").reset_index(), work_dir / "design-matrix.csv")
    _write_csv(features.reset_index(drop=True), output_dir / "organoid-feature-universe.csv.gz")
    qc = {"analysis": "GSE239283_donor_paired_PSC_by_IL17_interaction", "plan_sha256": plan_digest,
        "script_sha256": sha256(Path(__file__).resolve()), "software": arm["software"],
        "count_quantity": arm.get("count_quantity"), "raw_count_gate_passed": arm.get("raw_count_gate_passed") is True,
        "feature_namespace": arm["feature_namespace"], "feature_universe_description": arm["feature_universe_description"],
        "released_features": len(features), "eligible_features": len(eligible_ids), "donors": 8,
        "libraries": len(selected), "released_cells": int(selected.n_cells.sum()), "design": design_qc,
        "ratio_basis_positive_in_all_16_libraries": ratio_basis,
        "released_library_total_min": int(totals.min()), "released_library_total_max": int(totals.max()),
        "cells_per_library_min": int(selected.n_cells.min()), "cells_per_library_max": int(selected.n_cells.max()),
        "source_pins": arm["source_artifacts"], "derived_input_pins": {key: arm[key] for key in paths},
        "code_pins": arm["code_artifacts"], "filter": arm["filter"], "model": arm["model"],
        "sensitivity": arm["sensitivity"], "effects_computed": False,
        "interpretation_limits": [
            "eight_donors_not_cells_or_libraries_are_independent_units",
            "four_vs_four_small_n_NB_Wald_calibration_limited",
            "all_zero_required_count_arm_makes_group_log_response_boundary_or_unidentified_not_a_finite_Wald_test",
            "zero_total_individual_donor_pairs_flagged_without_new_filter_or_donor_exclusion",
            "PSC_status_not_randomized_and_nonPSC_are_procedure_controls_not_healthy_volunteers",
            "donor_blocking_does_not_remove_age_by_treatment_or_clinical_response_confounding",
            "published_IL17_responses_known_not_untouched_validation",
            "all_released_cells_no_treatment_selected_clusters",
            "released_feature_universe_not_a_claim_about_complete_historical_raw_libraries",
            "BH_requires_dependence_conditions_BY_sensitivity_uses_same_family",
            "within_group_significance_or_DEG_count_differences_do_not_test_interaction",
            "Cook_distances_diagnostic_only_no_automatic_exclusion_or_count_replacement",
            "Welch_and_NB_are_different_expression_scales_not_independent_replications",
            "leaveout_ranges_not_confidence_intervals_or_independent_replications",
            "allocation_fraction_not_randomized_disease_causal_P_or_discovery_criterion",
            "ex_vivo_response_not_treatment_efficacy_or_a_PSC_cause"]}
    if validate_only:
        qc["status"] = "structural_validation_only_no_effects"
        _write_json(output_dir / "organoid-validation.json", qc)
        return qc
    with (work_dir / "pydeseq2-fit.log").open("x") as log, redirect_stdout(log), redirect_stderr(log):
        fit = fit_negative_binomial(eligible_counts, selected, design, arm["n_cpus"])
    fit.library_diagnostics["unfiltered_released_library_total"] = totals
    _write_csv(fit.library_diagnostics, work_dir / "library-normalization.csv")
    _write_csv(fit.cooks.rename_axis("library_id").reset_index(), work_dir / "all-cooks-distances.csv.gz")
    _write_csv(fit.coefficients.rename_axis("source_feature_id").reset_index(), work_dir / "NB-natural-log-coefficients.csv.gz")
    _write_csv(fit.feature_fit_state.rename_axis("source_feature_id").reset_index(), work_dir / "NB-all-native-feature-estimates-and-flags.csv.gz")
    _write_json(work_dir / "fit-warnings.json", fit.diagnostics["warnings"])
    # Full released totals are used before selecting fixed eligible outcomes.
    complete_logcpm = pd.DataFrame(log2_cpm(data.counts), index=data.counts.index, columns=data.counts.columns)
    changes, disease = paired_changes(selected, complete_logcpm.loc[:, eligible_ids])
    _write_csv(changes.join(disease).reset_index(), work_dir / "donor-paired-log2cpm-changes.csv.gz")
    welch = welch_changes(changes, disease)
    leaveout, leaveout_values = leave_one_donor_out(changes, disease)
    allocations, allocation_values, allocation_labels = allocation_sensitivity(changes, disease)
    _write_csv(leaveout_values.reset_index(), work_dir / "all-donor-leaveout-effects.csv.gz")
    _write_csv(allocation_values.reset_index(), work_dir / "all-70-allocation-effects.csv.gz")
    _write_csv(allocation_labels, work_dir / "all-70-allocation-labels.csv")
    primary = _feature_join(features, fit.results["interaction"].join(fit.feature_diagnostics, validate="one_to_one"))
    primary.loc[~primary.eligible, "status"] = "excluded_fixed_count_filter"
    primary["interaction_q_bh_le_0_05"] = primary.eligible & primary.pvalue.notna() & primary.q_bh.le(0.05)
    primary["interaction_q_by_le_0_05"] = primary.eligible & primary.pvalue.notna() & primary.q_by.le(0.05)
    _write_csv(primary.reset_index(drop=True), output_dir / "organoid-primary.csv.gz")
    within = []
    for contrast in ["nonPSC_treatment", "PSC_treatment"]:
        table = _feature_join(features, fit.results[contrast])
        table["contrast"] = contrast; table.loc[~table.eligible, "status"] = "excluded_fixed_count_filter"
        within.append(table)
    _write_csv(pd.concat(within).reset_index(drop=True), output_dir / "organoid-within-group.csv.gz")
    for name, table in [("welch", welch), ("leave-one-donor-out", leaveout), ("allocations", allocations)]:
        full = _feature_join(features, table)
        full.loc[~full.eligible, "status"] = "excluded_fixed_count_filter"
        _write_csv(full.reset_index(drop=True), output_dir / f"organoid-{name}.csv.gz")
    qc.update({"effects_computed": True, "status": "complete", "negative_binomial": {key: value for key, value in fit.diagnostics.items() if key != "warnings"},
        "fit_warning_count": len(fit.diagnostics["warnings"]), "welch_finite_pvalues": int(welch.pvalue.notna().sum()),
        "allocation_count": 70, "allocation_minimum_possible_two_sided_fraction": 2/70,
        "donor_values_in_public_outputs": False})
    _write_json(output_dir / "organoid-qc.json", qc)
    names = ["organoid-feature-universe.csv.gz", "organoid-primary.csv.gz", "organoid-within-group.csv.gz",
        "organoid-welch.csv.gz", "organoid-leave-one-donor-out.csv.gz", "organoid-allocations.csv.gz", "organoid-qc.json"]
    manifest = {"plan_sha256": plan_digest, "script_sha256": qc["script_sha256"],
        "outputs": {name: {"sha256": sha256(output_dir / name), "bytes": (output_dir / name).stat().st_size} for name in names},
        "donor_level_outputs": "ignored_work_directory_only", "source_values_mutated": False}
    _write_json(output_dir / "organoid-output-manifest.json", manifest)
    return qc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--plan-sha256", help="Optional external plan pin")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--validate-only", action="store_true", help="Structural/QC/count-filter checks only; no effects")
    args = parser.parse_args()
    work_existed = args.work_dir.exists() or args.work_dir.is_symlink()
    try:
        qc = run_analysis(args.plan, args.output_dir, args.work_dir, args.plan_sha256, args.validate_only)
    except Exception as error:
        # Never create or modify a preexisting/rejected work directory, including
        # its aliases. An inference failure in this run's fresh work area may log.
        if not work_existed and args.work_dir.is_dir() and args.work_dir.resolve().is_relative_to(REPO / "work"):
            try:
                _write_json(args.work_dir / "organoid-failure.json", {"status": "not_completed", "error_type": type(error).__name__, "error": str(error)})
            except FileExistsError:
                pass
        raise
    print(json.dumps({key: qc[key] for key in ["status", "effects_computed", "donors", "libraries", "released_features", "eligible_features"]}))


if __name__ == "__main__":
    main()
