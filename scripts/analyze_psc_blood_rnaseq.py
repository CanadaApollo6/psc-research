#!/usr/bin/env python3
"""Frozen, genome-wide GSE177044 Norwegian PSC blood RNA-seq analysis.

This script does not download data. It accepts only checksum-pinned inputs from
an explicitly frozen plan. Sample-level joins, design matrices and fit logs go
under ignored work/. Public outputs contain only per-gene or aggregate results.
Source gene labels are not silently converted to validated gene symbols.
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
import warnings
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution, norm

REPO = Path(__file__).resolve().parents[1]
PRIMARY_PLATES = ("PS449", "PS450", "PS451", "PS452", "PS453")
PRIMARY_COUNTS = {"Control": 77, "PSC": 43, "PSCUC": 177}
MODEL_SETTINGS = {
    "contrast": "PSC_vs_Control",
    "covariates": ["age_centered", "sex", "plate"],
    "size_factors_fit_type": "ratio",
    "fit_type": "parametric",
    "dispersion_fallback": "mean_if_parametric_fails",
    "refit_cooks": False,
    "independent_filter": False,
    "cooks_filter": True,
    "alpha": 0.05,
    "bh_family": "all_count_eligible_genes_NA_as_1_for_adjustment_only",
    "lfc_shrinkage": False,
}
FILTER_SETTINGS = {"min_count": 10, "min_samples": 30}
SENSITIVITY_SETTINGS = {
    "transformation": "log2(1+CPM)",
    "library_size": "all_source_features",
    "covariance": "HC3",
    "reference_distribution": "normal",
    "leave_one_plate_out": "OLS_effect_sign_only",
    "subgroup": "joint_three_level_diagnosis_OLS_HC3",
}


class AnalysisError(ValueError):
    """An input or design does not satisfy the prespecified analysis contract."""


@dataclass
class CountData:
    counts: pd.DataFrame  # samples x source genes, in the original file order
    features: pd.DataFrame  # indexed by unchanged source gene ID


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_pin(spec: dict[str, Any], base: Path = REPO) -> Path:
    """Validate a complete local file, not an HTTP range/header prefix."""
    expected = spec.get("sha256", "")
    if not isinstance(expected, str) or not re.fullmatch(r"[a-f0-9]{64}", expected):
        raise AnalysisError("Every input requires a lowercase SHA-256 pin")
    path = Path(spec["path"])
    path = path if path.is_absolute() else base / path
    if not path.is_file():
        raise AnalysisError(f"Pinned input is missing: {path}")
    if path.name.endswith(".prefix"):
        raise AnalysisError("A partial source prefix is not an analysis input")
    actual = sha256(path)
    if actual != expected:
        raise AnalysisError(f"SHA-256 mismatch for {path}: {actual}")
    if "bytes" in spec and path.stat().st_size != spec["bytes"]:
        raise AnalysisError(f"Byte length mismatch for {path}")
    return path.resolve()


def _nonempty_unique(values: Any, label: str) -> list[str]:
    values = list(values)
    if any(not isinstance(v, str) or not v or v != v.strip() for v in values):
        raise AnalysisError(f"{label} contains missing, non-string or padded IDs")
    if len(set(values)) != len(values):
        raise AnalysisError(f"{label} contains duplicate IDs")
    return values


def versionless_ensembl(value: str) -> str:
    """Strip an Ensembl gene version only; do not infer an ID from a name."""
    match = re.fullmatch(r"(ENSG[0-9]{11})(?:\.[0-9]+)?", value)
    return match.group(1) if match else ""


def validate_counts(values: Any) -> np.ndarray:
    """Require finite nonnegative integer counts and nonzero sample libraries."""
    array = np.asarray(values)
    if array.ndim != 2 or not all(array.shape):
        raise AnalysisError("Counts must be a nonempty samples-by-features matrix")
    if array.dtype.kind == "b":
        raise AnalysisError("Boolean values are not raw read counts")
    try:
        numbers = np.asarray(array, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise AnalysisError("Nonnumeric or missing count value") from error
    if not np.isfinite(numbers).all():
        raise AnalysisError("Nonfinite or missing count value")
    if (numbers < 0).any() or not np.equal(numbers, np.floor(numbers)).all():
        raise AnalysisError("Counts must be nonnegative integers; do not round")
    if (numbers > 2**53 - 1).any() or (numbers.sum(axis=1) > 2**53 - 1).any():
        raise AnalysisError("Counts exceed the exact-integer validation range")
    if (numbers.sum(axis=1) == 0).any():
        raise AnalysisError("At least one sample has an all-zero library")
    return numbers.astype(np.int64)


def read_counts(path: Path) -> CountData:
    """Read every feature/count; reject duplicate physical CSV header names."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        header = next(csv.reader(stream), [])
    if header[:2] != ["Geneid", "gene_name"] or len(header) < 3:
        raise AnalysisError("Expected Geneid,gene_name followed by sample columns")
    _nonempty_unique(header, "Physical count-file columns")
    frame = pd.read_csv(path, dtype={"Geneid": str, "gene_name": str}, keep_default_na=False)
    if list(frame.columns) != header or not isinstance(frame.index, pd.RangeIndex):
        raise AnalysisError("CSV header or row width changed during parsing")
    genes = _nonempty_unique(frame["Geneid"], "Source gene IDs")
    if not genes:
        raise AnalysisError("Count file has no source features")
    counts = validate_counts(frame.iloc[:, 2:].to_numpy().T)
    features = pd.DataFrame({
        "source_gene_id": genes,
        "source_gene_name": frame["gene_name"].tolist(),
        "ensembl_gene_id_versionless": [versionless_ensembl(g) for g in genes],
        "source_feature_row": np.arange(1, len(genes) + 1),
    }, index=pd.Index(genes, name="source_gene_id_index"))
    bases = features["ensembl_gene_id_versionless"]
    features["ensembl_id_status"] = np.where(bases.eq(""), "not_an_ensembl_gene_id",
        np.where(bases.duplicated(keep=False), "ambiguous_versionless_id", "unique_versionless_id"))
    features["source_name_status"] = np.where(features.source_gene_name.eq(""), "missing_as_published", "as_published_unreconciled")
    return CountData(pd.DataFrame(counts, index=header[2:], columns=genes), features)


def parse_metadata(document: dict[str, Any], series: str = "GSE177044") -> pd.DataFrame:
    """Parse the archived GEO fields without using contact country as cohort."""
    records = document.get(series)
    if not isinstance(records, list) or not records:
        raise AnalysisError(f"Missing metadata series {series}")
    rows: list[dict[str, Any]] = []
    for record in records:
        titles = record.get("title", [])
        if not isinstance(titles, list) or len(titles) != 1:
            raise AnalysisError("Each GEO sample must have exactly one title")
        chars: dict[str, str] = {}
        for field in record.get("characteristics_ch1", []):
            if ":" not in field:
                raise AnalysisError("Unparseable GEO characteristic")
            key, value = field.split(":", 1)
            key = key.strip().casefold()
            if not key or key in chars:
                raise AnalysisError("Repeated or empty GEO characteristic key")
            chars[key] = value.strip()
        plate = chars.get("sequencing plate", chars.get("sequencing plate number (used for batch correction)", ""))
        rows.append({
            "sample_id": titles[0],
            "sample_accession": record.get("sample_accession", ""),
            "source_disease": chars.get("disease", ""),
            "age_source": chars.get("age", ""),
            "sex_source": chars.get("sex", ""),
            "plate": plate,
            "source_characteristics_json": json.dumps(record.get("characteristics_ch1", []), ensure_ascii=False),
        })
    _nonempty_unique([r["sample_id"] for r in rows], "GEO sample titles")
    _nonempty_unique([r["sample_accession"] for r in rows], "GEO sample accessions")
    return pd.DataFrame(rows).set_index("sample_id", drop=False)


def read_metadata(path: Path, series: str = "GSE177044") -> pd.DataFrame:
    """Accept the pinned original GEO quick SOFT or its complete JSON extract."""
    if path.suffix == ".json":
        return parse_metadata(json.loads(path.read_text()), series)
    records: list[dict[str, Any]] = []
    record: dict[str, Any] | None = None
    for line in path.read_text().splitlines():
        if line.startswith("^SAMPLE = "):
            record = {"sample_accession": line.split(" = ", 1)[1]}
            records.append(record)
        elif record is not None and line.startswith("!Sample_") and " = " in line:
            key, value = line.split(" = ", 1)
            key = key.removeprefix("!Sample_")
            if key == "geo_accession" and value != record["sample_accession"]:
                raise AnalysisError("GEO accession fields disagree")
            record.setdefault(key, []).append(value)
    if any(r.get("series_id") != [series] for r in records):
        raise AnalysisError("Unexpected or missing sample series assignment in GEO SOFT")
    return parse_metadata({series: records}, series)


def align_primary_samples(counts: pd.DataFrame, metadata: pd.DataFrame,
                          plates: tuple[str, ...] = PRIMARY_PLATES,
                          expected_counts: dict[str, int] | None = None,
                          expected_metadata_samples: int | None = None,
                          ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select by frozen plates; require the entire primary file to join exactly.

    Every count-file sample must be primary and every selected metadata sample
    must be in the count file. Excluded GEO records remain in the work audit.
    """
    _nonempty_unique(counts.index, "Count sample IDs")
    _nonempty_unique(counts.columns, "Count feature IDs")
    _nonempty_unique(metadata.index, "Metadata sample IDs")
    if expected_metadata_samples is not None and len(metadata) != expected_metadata_samples:
        raise AnalysisError("Unexpected number of source GEO metadata records")
    audit = metadata.copy()
    audit["selected"] = audit.plate.isin(plates)
    audit["selection_status"] = np.where(audit.selected, "primary_Norwegian_plate", "excluded_nonprimary_plate")
    audit["in_count_file"] = audit.index.isin(counts.index)
    selected = audit.loc[audit.selected].copy()
    if set(counts.index) != set(selected.index):
        absent = sorted(set(selected.index) - set(counts.index))
        extra = sorted(set(counts.index) - set(selected.index))
        raise AnalysisError(f"Non-exact primary sample join: missing={absent[:8]}, extra={extra[:8]}")
    # This is label alignment, not an assumption that GEO/file orders agree.
    selected = selected.loc[counts.index].copy()
    if set(selected.plate) != set(plates):
        raise AnalysisError("Not every prespecified sequencing plate is represented")
    if not selected.source_disease.isin(["Control", "PSC", "PSCUC"]).all():
        raise AnalysisError("Unexpected primary disease label; UC must not be pooled")
    actual_counts = selected.source_disease.value_counts().to_dict()
    if expected_counts is not None and actual_counts != expected_counts:
        raise AnalysisError(f"Unexpected source disease counts: {actual_counts}")
    selected["disease"] = selected.source_disease.replace({"PSCUC": "PSC"})
    try:
        selected["age"] = pd.to_numeric(selected.age_source, errors="raise")
    except (TypeError, ValueError) as error:
        raise AnalysisError("Missing or invalid primary age") from error
    if not np.isfinite(selected.age.to_numpy(dtype=float)).all() or not selected.age.between(0, 120).all():
        raise AnalysisError("Missing or implausible primary age")
    selected["sex"] = selected.sex_source.str.casefold()
    if not selected.sex.isin(["female", "male"]).all():
        raise AnalysisError("Missing or unrecognized primary sex; no imputation")
    return selected, audit


def check_design(design: pd.DataFrame | np.ndarray) -> dict[str, Any]:
    x = np.asarray(design, dtype=float)
    if x.ndim != 2 or not np.isfinite(x).all():
        raise AnalysisError("Nonfinite or nonmatrix design")
    n, p = x.shape
    rank = int(np.linalg.matrix_rank(x))
    if rank != p:
        raise AnalysisError(f"Rank-deficient design: rank={rank}, columns={p}")
    if n <= p:
        raise AnalysisError("No residual degrees of freedom")
    q, _ = np.linalg.qr(x, mode="reduced")
    leverage = np.sum(q * q, axis=1)
    if np.any(1.0 - leverage <= 1e-10):
        raise AnalysisError("HC3 undefined: unit-leverage sample")
    return {"n_samples": n, "n_columns": p, "rank": rank,
            "residual_df": n - p, "condition_number": float(np.linalg.cond(x)),
            "max_leverage": float(leverage.max())}


def build_design(metadata: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Explicit reference levels make the positive effect PSC minus Control."""
    if set(metadata.disease) != {"Control", "PSC"}:
        raise AnalysisError("Both PSC and Control are required")
    if set(metadata.sex) != {"female", "male"}:
        raise AnalysisError("Both sex levels are required by the frozen model")
    age = metadata.age.to_numpy(dtype=float)
    if not np.isfinite(age).all():
        raise AnalysisError("Age is missing or nonfinite")
    plates = sorted(set(metadata.plate))
    if any(not p for p in plates):
        raise AnalysisError("Plate is missing")
    values: dict[str, Any] = {"Intercept": np.ones(len(metadata)),
        "disease_PSC": metadata.disease.eq("PSC").to_numpy(dtype=float),
        "age_centered": age - age.mean(),
        "sex_male": metadata.sex.eq("male").to_numpy(dtype=float)}
    for plate in plates[1:]:
        values[f"plate_{plate}"] = metadata.plate.eq(plate).to_numpy(dtype=float)
    design = pd.DataFrame(values, index=metadata.index)
    diagnostic = check_design(design)
    diagnostic.update({"columns": list(design.columns), "disease_reference": "Control",
        "sex_reference": "female", "plate_reference": plates[0], "age_center": float(age.mean()),
        "plate_by_disease": pd.crosstab(metadata.plate, metadata.disease).to_dict(orient="index"),
        "sex_by_disease": pd.crosstab(metadata.sex, metadata.disease).to_dict(orient="index")})
    # Report the remaining disease variation after adjusting all nuisance terms.
    d = design["disease_PSC"].to_numpy()
    nuisance = design.drop(columns="disease_PSC").to_numpy()
    residual = d - nuisance @ np.linalg.lstsq(nuisance, d, rcond=None)[0]
    diagnostic["disease_residual_sum_squares"] = float(residual @ residual)
    diagnostic["disease_nuisance_r_squared"] = float(1 - residual @ residual / np.sum((d - d.mean()) ** 2))
    return design, diagnostic


def expression_filter(counts: pd.DataFrame, min_count: int = 10,
                      min_samples: int = 10) -> pd.DataFrame:
    if min_count < 1 or min_samples < 1:
        raise AnalysisError("Invalid expression filter")
    values = validate_counts(counts.to_numpy())
    n_above = (values >= min_count).sum(axis=0)
    return pd.DataFrame({"eligible": n_above >= min_samples,
        "samples_count_ge_threshold": n_above,
        "nonzero_samples": (values > 0).sum(axis=0),
        "all_zero_selected": (values == 0).all(axis=0),
        "mean_raw_count_selected": values.mean(axis=0)}, index=counts.columns)


def adjust_pvalues(pvalues: Any, method: str = "bh", family_size: int | None = None) -> np.ndarray:
    """Keep NA positions; optional fixed family size also counts missing tests."""
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or np.isinf(p).any() or ((p < 0) | (p > 1)).any():
        raise AnalysisError("P-values must be one-dimensional and in [0,1] or NA")
    valid = np.flatnonzero(np.isfinite(p))
    n = len(valid)
    m = n if family_size is None else family_size
    if not isinstance(m, (int, np.integer)) or m < n:
        raise AnalysisError("Multiplicity family is smaller than the observed tests")
    result = np.full(p.shape, np.nan)
    if method not in {"bh", "bonferroni"}:
        raise AnalysisError("Unknown multiplicity method")
    if not n:
        return result
    if method == "bonferroni":
        result[valid] = np.minimum(1, m * p[valid])
    else:
        order = valid[np.argsort(p[valid], kind="stable")]
        adjusted = p[order] * m / np.arange(1, n + 1)
        result[order] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    return result


def log2_cpm(counts: pd.DataFrame) -> np.ndarray:
    """Use every original source feature in each selected sample's library."""
    values = validate_counts(counts.to_numpy())
    library_sizes = values.sum(axis=1, dtype=np.int64)
    return np.log2(1 + values.astype(float) / library_sizes[:, None] * 1_000_000)


def fit_ols_hc3(design: pd.DataFrame | np.ndarray, outcomes: Any,
                coefficient: int = 1, chunk_size: int = 2048) -> pd.DataFrame:
    """Vectorized OLS/HC3 with two-sided asymptotic-normal inference.

    Variance for coefficient j is sum_i [(X'X)^-1 X']_ji^2 *
    residual_i^2/(1-h_ii)^2. No homoskedastic variance or empirical Bayes step.
    """
    check_design(design)
    x = np.asarray(design, dtype=float)
    y = np.asarray(outcomes, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    if y.ndim != 2 or len(y) != len(x) or not np.isfinite(y).all():
        raise AnalysisError("Nonfinite outcomes or inconsistent sample alignment")
    if not 0 <= coefficient < x.shape[1] or chunk_size < 1:
        raise AnalysisError("Invalid OLS coefficient or chunk size")
    q, r = np.linalg.qr(x, mode="reduced")
    projection = np.linalg.solve(r, q.T)
    leverage = np.sum(q * q, axis=1)
    weights = projection[coefficient] ** 2 / (1 - leverage) ** 2
    effect = np.empty(y.shape[1])
    se = np.empty(y.shape[1])
    for start in range(0, y.shape[1], chunk_size):
        stop = min(start + chunk_size, y.shape[1])
        block = y[:, start:stop]
        coefficients = projection @ block
        residual = block - x @ coefficients
        effect[start:stop] = coefficients[coefficient]
        se[start:stop] = np.sqrt(weights @ (residual * residual))
    # A deterministic/perfect-fit response is not given an invented p-value.
    defined = se > np.finfo(float).eps * np.maximum(1, np.abs(effect)) * 32
    statistic = np.divide(effect, se, out=np.full_like(effect, np.nan), where=defined)
    pvalues = 2 * norm.sf(np.abs(statistic))
    quantile = norm.ppf(0.975)
    return pd.DataFrame({"effect_log2_cpm": effect, "hc3_se": se,
        "normal_z": statistic, "pvalue": pvalues, "padj_bh": adjust_pvalues(np.where(np.isfinite(pvalues), pvalues, 1.0)),
        "ci95_lower": np.where(defined, effect - quantile * se, np.nan),
        "ci95_upper": np.where(defined, effect + quantile * se, np.nan),
        "status": np.where(defined, "tested", "zero_residual_variance_not_tested")})


def leave_one_plate_out(metadata: pd.DataFrame, outcomes: Any, gene_ids: Any,
                        full_effects: Any) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Diagnostic OLS signs only; keep the full-model feature universe fixed."""
    y = np.asarray(outcomes, dtype=float)
    genes = _nonempty_unique(gene_ids, "Leave-one-plate-out gene IDs")
    full = np.asarray(full_effects, dtype=float)
    if y.shape != (len(metadata), len(genes)) or full.shape != (len(genes),):
        raise AnalysisError("Leave-one-plate-out alignment mismatch")
    result = pd.DataFrame({"source_gene_id": genes, "full_ols_effect": full})
    estimates: list[np.ndarray] = []
    diagnostics: dict[str, Any] = {}
    for plate in sorted(set(metadata.plate)):
        keep = metadata.plate.ne(plate).to_numpy()
        try:
            design, diagnostic = build_design(metadata.loc[keep])
            effects = np.linalg.lstsq(design.to_numpy(), y[keep], rcond=None)[0][1]
            diagnostics[plate] = {"status": "estimable", **diagnostic}
        except AnalysisError as error:
            effects = np.full(len(genes), np.nan)
            diagnostics[plate] = {"status": "not_estimable", "reason": str(error)}
        result[f"effect_without_{plate}"] = effects
        estimates.append(effects)
    stack = np.stack(estimates)
    result["n_estimable_omissions"] = np.isfinite(stack).sum(axis=0)
    result["n_same_nonzero_sign"] = ((np.sign(stack) == np.sign(full)) & (stack != 0) & (full != 0)).sum(axis=0)
    result["n_opposite_sign"] = ((np.sign(stack) == -np.sign(full)) & (stack != 0) & (full != 0)).sum(axis=0)
    result["interpretation"] = "diagnostic_only_not_independent_replications"
    return result, diagnostics


def fit_subgroup_sensitivity(metadata: pd.DataFrame, outcomes: Any,
                             gene_ids: Any) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Joint Control/PSC-alone/PSCUC model, sharing controls and covariates."""
    genes = _nonempty_unique(gene_ids, "Subgroup gene IDs")
    y = np.asarray(outcomes, dtype=float)
    if y.shape != (len(metadata), len(genes)):
        raise AnalysisError("Subgroup sample/feature alignment mismatch")
    if set(metadata.source_disease) != {"Control", "PSC", "PSCUC"}:
        raise AnalysisError("Three-level subgroup model requires every source diagnosis")
    base, _ = build_design(metadata)
    design = base.drop(columns="disease_PSC").copy()
    design.insert(1, "diagnosis_PSC_alone", metadata.source_disease.eq("PSC").astype(float))
    design.insert(2, "diagnosis_PSCUC", metadata.source_disease.eq("PSCUC").astype(float))
    diagnostic = check_design(design)
    tables = []
    for coefficient, label in [(1, "psc_alone"), (2, "pscuc")]:
        table = fit_ols_hc3(design, y, coefficient=coefficient).add_prefix(f"ols_{label}_")
        table.index = genes
        tables.append(table)
    diagnostic.update({"columns": list(design.columns),
        "source_disease_counts": metadata.source_disease.value_counts().to_dict(),
        "shared_control_and_nuisance_coefficients": True,
        "interpretation": "sensitivity_shared_controls_not_independent_replications"})
    return pd.concat(tables, axis=1), diagnostic


def ratio_size_factors(counts: pd.DataFrame) -> tuple[np.ndarray, int]:
    """Independently reproduce PyDESeq2's median-log-ratio convention.

    Taking exp(median(log ratios)) is the pinned library's convention, including
    the even-feature case. Do not silently switch to poscounts or iterative.
    """
    values = validate_counts(counts.to_numpy())
    usable = (values > 0).all(axis=0)
    if not usable.any():
        raise AnalysisError("Median-of-ratios unavailable: every eligible gene has a zero")
    log_counts = np.log(values[:, usable].astype(float))
    factors = np.exp(np.median(log_counts - log_counts.mean(axis=0), axis=1))
    if not np.isfinite(factors).all() or not (factors > 0).all():
        raise AnalysisError("Invalid median-of-ratios normalization factors")
    return factors, int(usable.sum())


def fit_negative_binomial(counts: pd.DataFrame, metadata: pd.DataFrame,
                          design: pd.DataFrame, n_cpus: int = 1,
                          ) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """PyDESeq2 with an explicitly allowed parametric-to-mean trend fallback.

    No count replacement, LFC shrinkage or independent filtering. NA Wald tests
    stay NA. Only the fixed-family BH correction input replaces them with one.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.utils import n_or_more_replicates

    if not counts.index.equals(metadata.index) or not counts.index.equals(design.index):
        raise AnalysisError("NB input sample order is not exactly aligned")
    check_design(design)
    validate_counts(counts.to_numpy())
    expected_factors, normalization_genes = ratio_size_factors(counts)
    if n_cpus < 1:
        raise AnalysisError("n_cpus must be positive")
    dds = DeseqDataSet(counts=counts, metadata=metadata, design=design,
        fit_type="parametric", size_factors_fit_type="ratio", refit_cooks=False,
        n_cpus=n_cpus, quiet=True)
    if not np.array_equal(dds.obsm["design_matrix"].to_numpy(), design.to_numpy()):
        raise AnalysisError("PyDESeq2 design differs from the explicit aligned design")
    captured: list[warnings.WarningMessage]
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        dds.fit_size_factors(fit_type="ratio")
        factors = dds.obs["size_factors"].to_numpy(dtype=float)
        if not np.allclose(factors, expected_factors, rtol=1e-12, atol=1e-12):
            raise AnalysisError("PyDESeq2 size factors differ from pinned ratio method")
        dds.fit_genewise_dispersions()
        dds.fit_dispersion_trend()
        if dds.uns.get("disp_function_type") not in {"parametric", "mean"}:
            raise AnalysisError("Unexpected dispersion fit type")
        dds.fit_dispersion_prior()
        dds.fit_MAP_dispersions()
        dds.fit_LFC()
        dds.calculate_cooks()
        cooks_outlier = np.asarray(dds.cooks_outlier(), dtype=bool)
        contrast = np.zeros(design.shape[1])
        contrast[design.columns.get_loc("disease_PSC")] = 1.0
        stats = DeseqStats(dds, contrast=contrast, alpha=0.05,
            cooks_filter=True, independent_filter=False, n_cpus=n_cpus, quiet=True)
        stats.summary()
    if not np.array_equal(dds.X, counts.to_numpy()):
        raise AnalysisError("Source counts changed during the model fit")
    results = stats.results_df.copy()
    if not results.index.equals(counts.columns):
        raise AnalysisError("NB feature order changed or features were silently removed")
    qvalues = adjust_pvalues(results.pvalue.to_numpy(dtype=float))
    if not np.allclose(results.padj, qvalues, rtol=1e-10, atol=1e-14, equal_nan=True):
        raise AnalysisError("PyDESeq2 BH adjustment does not reproduce")
    results["pydeseq2_padj"] = results["padj"]
    results["padj"] = adjust_pvalues(np.where(np.isfinite(results.pvalue), results.pvalue, 1.0))
    results["dispersion_fit_type"] = dds.uns["disp_function_type"]
    results["cooks_outlier"] = cooks_outlier
    cooks_sample_mask = n_or_more_replicates(design, 3).to_numpy(dtype=bool)
    results["max_cooks_distance_all_samples"] = np.max(dds.layers["cooks"], axis=0)
    results["max_cooks_distance_filter_eligible_samples"] = (
        np.max(dds.layers["cooks"][cooks_sample_mask], axis=0)
        if cooks_sample_mask.any() else np.full(len(results), np.nan))
    for field in ["genewise_dispersions", "fitted_dispersions", "dispersions",
                  "_genewise_converged", "_MAP_converged", "_LFC_converged"]:
        if field in dds.var:
            results[field.lstrip("_")] = dds.var[field].to_numpy()
    finite = np.isfinite(results.pvalue.to_numpy(dtype=float))
    results["model_status"] = np.where(cooks_outlier, "cooks_filtered_pvalue_unavailable",
        np.where(finite, "tested", "fit_pvalue_unavailable"))
    if dds.uns["disp_function_type"] == "mean":
        results["model_status"] = results.model_status + "_prespecified_mean_dispersion_fallback"
    results["ci95_lower_log2fc"] = results.log2FoldChange - norm.ppf(0.975) * results.lfcSE
    results["ci95_upper_log2fc"] = results.log2FoldChange + norm.ppf(0.975) * results.lfcSE
    sf_table = pd.DataFrame({"sample_id": counts.index, "size_factor": factors,
        "library_size_all_fitted_features": counts.sum(axis=1).to_numpy(),
        "cooks_filter_eligible_sample": cooks_sample_mask})
    diagnostics = {
        "method": "PyDESeq2_NB_Wald_PSC_minus_Control",
        "dispersion_trend": dds.uns["disp_function_type"],
        "normalization_basis_genes_positive_in_every_sample": normalization_genes,
        "normalization_convention": "exp(median(log(count/geometric_mean)))_over_eligible_positive_everywhere_genes",
        "size_factor_min": float(factors.min()), "size_factor_max": float(factors.max()),
        "fitted_features": len(results), "finite_pvalues": int(finite.sum()),
        "cooks_filtered_features": int(cooks_outlier.sum()),
        "cooks_filter_eligible_samples": int(cooks_sample_mask.sum()),
        "cooks_filter_sample_rule": "at_least_3_identical_full_design_rows_including_age",
        "cooks_distance_cutoff": float(f_distribution.ppf(0.99, design.shape[1], len(design) - design.shape[1])),
        "bh_family_size": len(results),
        "pydeseq2_original_bh_family_size": int(finite.sum()),
        "bh_missing_test_policy": "source_pvalue_NA_preserved_correction_input_one",
        "dispersion_fallback_used": dds.uns["disp_function_type"] == "mean",
        "warnings": [{"category": w.category.__name__, "message": str(w.message)} for w in captured],
    }
    for column in ["LFC_converged", "genewise_converged", "MAP_converged"]:
        if column in results:
            diagnostics[f"{column}_false"] = int(results[column].eq(False).sum())
    return results, sf_table, diagnostics


def validate_plan(plan: dict[str, Any]) -> dict[str, Any]:
    if plan.get("status") != "frozen":
        raise AnalysisError("Plan status must be frozen before any expression is read")
    stamp = plan.get("frozen_at_utc", "")
    try:
        frozen = datetime.fromisoformat(stamp.replace("Z", "+00:00"))
    except (ValueError, TypeError, AttributeError) as error:
        raise AnalysisError("Missing or invalid plan freeze timestamp") from error
    if frozen.tzinfo is None or frozen > datetime.now(timezone.utc):
        raise AnalysisError("Plan freeze timestamp must have a time zone and be in the past")
    arm = plan.get("rnaseq", {})
    selection = arm.get("selection", {})
    if selection.get("plates") != list(PRIMARY_PLATES):
        raise AnalysisError("Primary analysis must use exactly plates PS449 through PS453")
    if selection.get("source_disease_counts") != PRIMARY_COUNTS or selection.get("metadata_samples") != 1035:
        raise AnalysisError("The frozen primary source sample counts do not match GSE177044")
    if arm.get("metadata_series") != "GSE177044":
        raise AnalysisError("Only the prespecified GSE177044 series is supported")
    if arm.get("model") != MODEL_SETTINGS:
        raise AnalysisError("Model options differ from the supported frozen NB model")
    if arm.get("sensitivity") != SENSITIVITY_SETTINGS:
        raise AnalysisError("Sensitivity options differ from the supported frozen method")
    filtering = arm.get("filter", {})
    if filtering != FILTER_SETTINGS:
        raise AnalysisError("Count filter must be frozen at 10 counts in 30 selected samples")
    if arm.get("expected_sample_count") != 297:
        raise AnalysisError("The primary source must contain all 297 selected samples")
    if arm.get("expected_feature_count") != 63677 or arm.get("annotation_columns") != ["Geneid", "gene_name"]:
        raise AnalysisError("Expected all 63,677 source features with original annotation columns")
    if not isinstance(arm.get("counts"), list) or len(arm["counts"]) != 1:
        raise AnalysisError("Exactly one primary PSC count file must be pinned")
    pins = arm.get("software", {})
    required = {"python", "numpy", "pandas", "scipy", "pydeseq2"}
    if not required.issubset(pins):
        raise AnalysisError("Required Python/numpy/pandas/scipy/PyDESeq2 version pins are missing")
    for package, expected in pins.items():
        actual = ".".join(map(str, sys.version_info[:3])) if package == "python" else importlib.metadata.version(package)
        if actual != expected:
            raise AnalysisError(f"Software pin mismatch: {package}={actual}, expected {expected}")
    if not isinstance(arm.get("n_cpus"), int) or arm["n_cpus"] < 1:
        raise AnalysisError("An explicit positive n_cpus is required")
    for kind, source in [("counts", arm["counts"][0]), ("metadata", arm.get("metadata", {}))]:
        if not source.get("url") or not source.get("retrieved_at_utc"):
            raise AnalysisError(f"Source URL and retrieval timestamp required for {kind}")
    labels = arm.get("known_target_source_labels", [])
    _nonempty_unique(labels, "Known benchmark source labels")
    if not {"IFIT1", "G0S2"}.issubset(labels):
        raise AnalysisError("The source's IFIT1 and G0S2 benchmarks must be labeled known")
    return arm


def _write_json(path: Path, document: Any) -> None:
    path.write_text(json.dumps(document, sort_keys=True, indent=2, ensure_ascii=False,
        allow_nan=False) + "\n")


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.17g", na_rep="NA", lineterminator="\n")


def run_analysis(plan_path: Path, output_dir: Path, work_dir: Path,
                 plan_sha256: str | None = None) -> dict[str, Any]:
    """Execute after validating the freeze and every input/software pin."""
    plan_path = plan_path.resolve()
    plan_digest = sha256(plan_path)
    if plan_sha256 is not None and plan_digest != plan_sha256:
        raise AnalysisError("Execution plan SHA-256 mismatch")
    plan = json.loads(plan_path.read_text())
    arm = validate_plan(plan)
    counts_path = check_pin(arm["counts"][0])
    metadata_path = check_pin(arm["metadata"])
    if "script" in arm:
        if check_pin(arm["script"]) != Path(__file__).resolve():
            raise AnalysisError("Pinned analysis script is not the executing script")
    output_dir = output_dir.resolve()
    work_dir = work_dir.resolve()
    if not work_dir.is_relative_to(REPO / "work"):
        raise AnalysisError("Sample-level work directory must be under ignored repository work/")
    if output_dir.is_relative_to(REPO / "data/raw"):
        raise AnalysisError("Outputs must not overwrite original raw inputs")
    if output_dir == work_dir or output_dir.is_relative_to(work_dir):
        raise AnalysisError("Public outputs and sample-level work must use distinct directories")
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "execution-plan.json").write_bytes(plan_path.read_bytes())
    metadata = read_metadata(metadata_path, arm["metadata_series"])
    data = read_counts(counts_path)
    if len(data.features) != arm["expected_feature_count"]:
        raise AnalysisError("Count file does not contain the complete expected feature universe")
    selected, sample_audit = align_primary_samples(data.counts, metadata,
        tuple(arm["selection"]["plates"]), arm["selection"]["source_disease_counts"],
        arm["selection"]["metadata_samples"])
    design, design_qc = build_design(selected)
    _write_csv(sample_audit.reset_index(drop=True), work_dir / "sample-selection-audit.csv")
    _write_csv(selected.reset_index(drop=True), work_dir / "selected-samples.csv")
    _write_csv(design.rename_axis("sample_id").reset_index(), work_dir / "design-matrix.csv")
    filter_table = expression_filter(data.counts, **arm["filter"])
    eligible_ids = filter_table.index[filter_table.eligible]
    if not len(eligible_ids):
        raise AnalysisError("No features pass the prespecified count filter")
    features = data.features.join(filter_table, validate="one_to_one")
    features["gene_id"] = features.ensembl_gene_id_versionless
    features["known_benchmark_source_label"] = features.source_gene_name.isin(arm["known_target_source_labels"])
    features["benchmark_label_status"] = np.where(features.known_benchmark_source_label,
        "known_not_discovery_literal_source_label_only", "not_a_prespecified_benchmark_source_label")
    # Count fitting sees all and only the fixed eligible source features.
    with (work_dir / "pydeseq2-fit.log").open("w") as log, redirect_stdout(log), redirect_stderr(log):
        nb_results, size_factors, nb_qc = fit_negative_binomial(
            data.counts.loc[:, eligible_ids], selected, design, arm["n_cpus"])
    size_factors["library_size_all_source_features"] = data.counts.sum(axis=1).to_numpy()
    _write_csv(size_factors, work_dir / "sample-normalization.csv")
    # Full-source libraries are calculated before subsetting the eligible genes.
    full_logcpm = log2_cpm(data.counts)
    eligibility = filter_table.eligible.to_numpy()
    outcomes = full_logcpm[:, eligibility]
    ols = fit_ols_hc3(design, outcomes)
    ols.index = eligible_ids
    ols = ols.add_prefix("ols_")
    leaveout, leaveout_qc = leave_one_plate_out(selected, outcomes, eligible_ids,
        ols.ols_effect_log2_cpm.to_numpy())
    subgroup, subgroup_qc = fit_subgroup_sensitivity(selected, outcomes, eligible_ids)
    primary = features.join(nb_results, validate="one_to_one").join(ols, validate="one_to_one").join(subgroup, validate="one_to_one")
    primary.loc[~primary.eligible, "model_status"] = "excluded_prespecified_count_filter"
    for status_column in ["ols_status", "ols_psc_alone_status", "ols_pscuc_status"]:
        primary.loc[~primary.eligible, status_column] = "excluded_prespecified_count_filter"
    # Membership is fixed before effects. NA and filtered rows remain explicit.
    primary["bh_q_le_0_05"] = primary.padj.le(0.05) & primary.padj.notna()
    full_leaveout = features[["source_gene_id", "gene_id", "eligible"]].merge(
        leaveout, on="source_gene_id", how="left", validate="one_to_one", sort=False)
    _write_csv(features.reset_index(drop=True), output_dir / "rnaseq-feature-universe.csv")
    _write_csv(primary.reset_index(drop=True), output_dir / "rnaseq-primary.csv")
    _write_csv(full_leaveout, output_dir / "rnaseq-leave-one-plate-out.csv")
    qc = {
        "analysis": "GSE177044_Norwegian_PSC_vs_Control",
        "plan_sha256": plan_digest,
        "script_sha256": sha256(Path(__file__).resolve()),
        "sources": {kind: arm[kind] for kind in ["counts", "metadata"]},
        "software": arm["software"],
        "source_metadata_samples": len(metadata), "selected_samples": len(selected),
        "excluded_metadata_samples": int((~sample_audit.selected).sum()),
        "count_file_samples": len(data.counts), "complete_primary_sample_join": True,
        "source_features": len(features), "eligible_features": int(features.eligible.sum()),
        "all_zero_features": int(features.all_zero_selected.sum()),
        "missing_source_names": int(features.source_gene_name.eq("").sum()),
        "versionless_id_status_counts": features.ensembl_id_status.value_counts().to_dict(),
        "source_disease_counts": selected.source_disease.value_counts().to_dict(),
        "raw_count_validity": "all_values_finite_nonnegative_integers_no_imputation_no_replacement",
        "filter": arm["filter"], "model": arm["model"], "sensitivity": arm["sensitivity"],
        "design": design_qc, "negative_binomial": nb_qc,
        "ols_finite_pvalues": int(ols.ols_pvalue.notna().sum()),
        "leave_one_plate_out": leaveout_qc,
        "subgroup_sensitivity": subgroup_qc,
        "interpretation_limits": [
            "observational_total_blood_association_not_causal_or_cell_intrinsic",
            "PSC_vs_UC_country_plate_confounding_precludes_a_pooled_contrast",
            "source_participant_labels_not_independently_verified_donor_identities",
            "cell_composition_treatment_and_disease_severity_not_adjusted",
            "known_source_targets_and_prior_candidates_are_not_new_discoveries",
            "missing_statistical_result_is_not_a_biological_null",
            "PyDESeq2_Cooks_filter_only_checks_samples_in_at_least_3_identical_full_design_rows_including_age",
            "leave_one_plate_out_fits_share_samples_and_are_not_independent_replications",
            "source_gene_name_is_as_published_unreconciled_not_an_identifier_map",
        ],
    }
    _write_json(output_dir / "rnaseq-qc.json", qc)
    output_names = ["rnaseq-feature-universe.csv", "rnaseq-primary.csv",
                    "rnaseq-leave-one-plate-out.csv", "rnaseq-qc.json"]
    manifest = {"plan_sha256": plan_digest, "script_sha256": qc["script_sha256"],
        "outputs": {name: {"sha256": sha256(output_dir / name),
                           "bytes": (output_dir / name).stat().st_size} for name in output_names},
        "sample_values_in_public_outputs": False,
        "sample_work_files": ["sample-selection-audit.csv", "selected-samples.csv",
                              "design-matrix.csv", "sample-normalization.csv"],
    }
    _write_json(output_dir / "rnaseq-output-manifest.json", manifest)
    return qc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--plan-sha256", help="Optional external SHA-256 pin for the frozen plan")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        qc = run_analysis(args.plan, args.output_dir, args.work_dir, args.plan_sha256)
    except Exception as error:
        # Do not leak per-sample exception text into a public aggregate report.
        if args.work_dir.resolve().is_relative_to(REPO / "work"):
            args.work_dir.mkdir(parents=True, exist_ok=True)
            _write_json(args.work_dir / "rnaseq-failure.json", {
                "status": "not_completed", "error_type": type(error).__name__, "error": str(error)})
        raise
    print(json.dumps({"status": "complete", "selected_samples": qc["selected_samples"],
                      "source_features": qc["source_features"], "eligible_features": qc["eligible_features"]}))


if __name__ == "__main__":
    main()
