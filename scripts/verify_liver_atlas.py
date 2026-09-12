"""Verify the reproducible liver-atlas expression tables and optional inputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from analyze_liver_atlas import (
    GENES,
    MIN_DONORS_FOR_CONTRAST,
    PRIMARY_MIN_CELLS,
    ROOT,
    SENSITIVITY_MIN_CELLS,
    aggregate_adata,
    load_frozen_plan,
    sha256_file,
)


PUBLIC_OUTPUTS = (
    "data/derived/liver-atlas-expression-cell-accounting.csv",
    "data/derived/liver-atlas-expression-gene-audit.csv",
    "data/derived/liver-atlas-expression-summary.csv",
    "data/derived/liver-atlas-expression-contrasts.csv",
)
DONOR_CACHE = "work/liver-atlas-expression-donor-pseudobulk.csv"


def _fail(message: str) -> None:
    raise ValueError(message)


def _assert_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        _fail(f"{label} is missing columns: {', '.join(missing)}")


def _assert_integer_series(frame: pd.DataFrame, column: str, label: str) -> None:
    values = frame[column].dropna().to_numpy(dtype=float)
    if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all() or (values < 0).any():
        _fail(f"{label}.{column} contains invalid nonnegative integer values")


def _verify_donor_cache(frame: pd.DataFrame, plan: dict[str, Any]) -> None:
    required = ["dataset_key", "dataset_id", "modality", "assay", "donor_id", "cohort", "author_cell_type",
                "cell_count", "all_feature_umis", "positive_all_feature_umi_denominator"]
    _assert_columns(frame, required, "donor pseudobulk")
    keys = ["dataset_key", "dataset_id", "modality", "suspension_type", "assay", "donor_id", "cohort", "author_cell_type"]
    if frame.duplicated(keys).any():
        _fail("Donor pseudobulk has duplicate donor/assay/cell-type groups")
    _assert_integer_series(frame, "cell_count", "donor pseudobulk")
    _assert_integer_series(frame, "all_feature_umis", "donor pseudobulk")
    for gene in GENES:
        for suffix in ("umis", "detected_cells"):
            _assert_integer_series(frame, f"{gene}_{suffix}", "donor pseudobulk")
        if (frame[f"{gene}_detected_cells"] > frame.cell_count).any():
            _fail(f"{gene} detected cells exceed group cell count")
        present = frame[f"{gene}_status"].eq("present")
        if (frame.loc[~present, f"{gene}_umis"] != 0).any() or frame.loc[~present, f"{gene}_cpm"].notna().any():
            _fail(f"Missing/ambiguous {gene} was represented as measured expression")
        measured = frame.loc[present & frame.positive_all_feature_umi_denominator]
        expected = measured[f"{gene}_umis"] / measured.all_feature_umis * 1_000_000
        if len(measured) and not np.allclose(measured[f"{gene}_cpm"], expected, rtol=0, atol=1e-9):
            _fail(f"{gene} CPM does not use the all-feature denominator")
        finite_log = measured[f"{gene}_log2_cpm1"].dropna().to_numpy(dtype=float)
        if len(finite_log) and not np.isfinite(finite_log).all():
            _fail(f"{gene} log2(CPM+1) contains nonfinite values")


def _verify_accounting(frame: pd.DataFrame, plan: dict[str, Any]) -> dict[str, int]:
    required = ["dataset_key", "modality", "suspension_type", "assay", "cohort", "author_cell_type",
                "ontology_cell_type", "exclusion_reason", "cell_count", "retained_cell_count", "excluded_cell_count"]
    _assert_columns(frame, required, "cell accounting")
    _assert_integer_series(frame, "cell_count", "cell accounting")
    _assert_integer_series(frame, "retained_cell_count", "cell accounting")
    _assert_integer_series(frame, "excluded_cell_count", "cell accounting")
    if not np.array_equal(frame.cell_count.to_numpy(), (frame.retained_cell_count + frame.excluded_cell_count).to_numpy()):
        _fail("Cell accounting retained plus excluded counts do not equal total counts")
    observed = {str(key): int(value) for key, value in frame.groupby("dataset_key").cell_count.sum().items()}
    expected = {str(spec["key"]): int(spec["expected_observations"]) for spec in plan["datasets"]}
    if observed != expected:
        _fail(f"Cell accounting does not reconcile to frozen observations: {observed} vs {expected}")
    return observed


def _verify_summary(frame: pd.DataFrame) -> None:
    required = ["dataset_key", "modality", "suspension_type", "assay", "author_cell_type", "cohort", "gene",
                "cell_threshold", "status", "eligible_donors", "observed_donors",
                "donor_weighted_median_log2_cpm1", "donor_weighted_median_detection_fraction",
                "observed_detection_fraction"]
    _assert_columns(frame, required, "expression summary")
    if not frame.gene.isin(GENES).all():
        _fail("Expression summary contains an unexpected gene")
    if not set(frame.cell_threshold.dropna().astype(int)).issubset({PRIMARY_MIN_CELLS, SENSITIVITY_MIN_CELLS}):
        _fail("Expression summary contains a non-frozen cell threshold")
    keys = ["dataset_key", "modality", "suspension_type", "assay", "author_cell_type", "cohort", "gene", "cell_threshold"]
    if frame.duplicated(keys).any():
        _fail("Expression summary has duplicate stratum/gene/threshold rows")
    for _, row in frame.iterrows():
        if row["status"] in {"missing_feature_name", "ambiguous_feature_name"}:
            if pd.notna(row["donor_weighted_median_log2_cpm1"]):
                _fail("Missing/ambiguous gene has a numeric expression summary")
        if row["eligible_donors"] < 0 or row["observed_donors"] < row["eligible_donors"]:
            _fail("Expression summary donor counts are inconsistent")


def _verify_contrasts(frame: pd.DataFrame) -> dict[str, int]:
    required = ["dataset_key", "modality", "suspension_type", "author_cell_type", "assay", "gene", "cell_threshold",
                "status", "shared_assay", "psc_eligible_donors", "control_eligible_donors"]
    _assert_columns(frame, required, "expression contrasts")
    if not frame.gene.isin(GENES).all():
        _fail("Expression contrasts contains an unexpected gene")
    keys = ["dataset_key", "modality", "suspension_type", "author_cell_type", "assay", "gene", "cell_threshold"]
    if frame.duplicated(keys).any():
        _fail("Expression contrasts has duplicate rows")
    evaluable = frame.status.astype(str).str.startswith("evaluable")
    if (evaluable & ~frame.shared_assay.astype(bool)).any():
        _fail("A contrast passed while assay matching was false")
    if (evaluable & (frame.psc_eligible_donors < MIN_DONORS_FOR_CONTRAST)).any() or (evaluable & (frame.control_eligible_donors < MIN_DONORS_FOR_CONTRAST)).any():
        _fail("A contrast passed without the minimum donor gate")
    if (evaluable & frame.modality.eq("sc")).any():
        _fail("An SC contrast was reported as biologically evaluable")
    sn_evaluable = evaluable & frame.modality.eq("sn")
    if sn_evaluable.any() and not frame.loc[sn_evaluable, "status"].astype(str).str.contains("residual_processing", regex=False).all():
        _fail("SN contrast lacks the residual processing-confounding flag")
    return {str(key): int(value) for key, value in frame.status.value_counts().items()}


def verify(
    root: Path = ROOT,
    *,
    require_inputs: bool = False,
    require_donor_cache: bool = False,
    chunk_size: int = 8192,
) -> dict[str, Any]:
    plan_path = root / "config/liver-atlas-plan.json"
    plan, manifest = load_frozen_plan(plan_path)
    audit_path = root / "reports/liver-atlas-expression-audit.json"
    if not audit_path.exists():
        _fail(f"Missing audit report: {audit_path}")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("plan_sha256") != sha256_file(plan_path):
        _fail("Expression audit plan hash does not match frozen plan")
    manifest_path = root / plan["source_manifest"]
    if audit.get("source_manifest_sha256") != sha256_file(manifest_path):
        _fail("Expression audit manifest hash does not match frozen manifest")
    for relative, field in (
        ("scripts/analyze_liver_atlas.py", "analysis_script_sha256"),
        ("requirements-liver-atlas.txt", "requirements_sha256"),
    ):
        if audit.get(field) != sha256_file(root / relative):
            _fail(f"Expression analysis provenance changed: {relative}")
    expected_outputs = dict(audit.get("outputs_sha256", {}))
    checked_outputs: list[str] = []
    missing_cache: list[str] = []
    for relative in PUBLIC_OUTPUTS + (DONOR_CACHE,):
        path = root / relative
        expected_hash = expected_outputs.get(relative)
        if not path.exists():
            if relative == DONOR_CACHE and not require_donor_cache:
                missing_cache.append(relative)
                continue
            _fail(f"Missing expression output: {relative}")
        if expected_hash is None or sha256_file(path) != expected_hash:
            _fail(f"Expression output hash changed: {relative}")
        checked_outputs.append(relative)

        if relative == DONOR_CACHE:
            _verify_donor_cache(pd.read_csv(path), plan)

    accounting = pd.read_csv(root / PUBLIC_OUTPUTS[0])
    observed_cells = _verify_accounting(accounting, plan)
    gene_audit = pd.read_csv(root / PUBLIC_OUTPUTS[1])
    _assert_columns(gene_audit, ["dataset_key", "gene", "status", "match_count"], "gene audit")
    if set(gene_audit.gene) != set(GENES) or gene_audit.duplicated(["dataset_key", "gene"]).any():
        _fail("Gene audit is incomplete or duplicated")
    summary = pd.read_csv(root / PUBLIC_OUTPUTS[2])
    _verify_summary(summary)
    contrasts = pd.read_csv(root / PUBLIC_OUTPUTS[3])
    contrast_counts = _verify_contrasts(contrasts)
    if audit.get("cell_accounting", {}).get("ledger_cells") != sum(observed_cells.values()):
        _fail("Audit ledger total does not match public accounting")
    if audit.get("contrast_status_counts") != contrast_counts:
        _fail("Audit contrast status counts do not match the public contrast table")

    input_qc: list[dict[str, Any]] = []
    if require_inputs:
        try:
            import anndata as ad
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("anndata is required for input verification") from exc
        source_by_file = {str(source["file"]): source for source in manifest["sources"]}
        for spec in plan["datasets"]:
            source = source_by_file[str(spec["file"])]
            path = root / manifest["source_directory"] / str(spec["file"])
            if not path.exists():
                _fail(f"Missing frozen atlas input: {path}")
            if path.stat().st_size != int(source["bytes"]) or sha256_file(path) != str(source["sha256"]):
                _fail(f"Frozen atlas input hash or length changed: {path}")
            data = ad.read_h5ad(path, backed="r")
            try:
                result = aggregate_adata(data, str(spec["key"]), expected_observations=int(spec["expected_observations"]), chunk_size=chunk_size)
                input_qc.append({"dataset_key": str(spec["key"]), **result.qc})
            finally:
                manager = getattr(data, "file", None)
                if manager is not None and hasattr(manager, "close"):
                    manager.close()
    return {
        "status": "all_expression_outputs_verified",
        "checked_outputs": checked_outputs,
        "ignored_donor_cache_missing": missing_cache,
        "cell_counts": observed_cells,
        "contrast_status_counts": contrast_counts,
        "input_qc": input_qc,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-inputs", action="store_true", help="Recheck both full raw matrices in chunks")
    parser.add_argument("--require-donor-cache", action="store_true", help="Require ignored work pseudobulk cache")
    parser.add_argument("--chunk-size", type=int, default=8192)
    args = parser.parse_args(argv)
    print(json.dumps(verify(require_inputs=args.require_inputs, require_donor_cache=args.require_donor_cache,
                            chunk_size=args.chunk_size), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
