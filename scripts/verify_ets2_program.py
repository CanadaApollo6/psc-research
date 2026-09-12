"""Verify ETS2 program aggregate outputs and optionally pinned full inputs.

The default mode is deliberately portable: it checks the audit, public CSV
hashes and aggregate-table contracts without opening the raw H5AD matrices or
the ignored donor/random-set cache.  ``--require-cache`` adds the ignored
cache hashes/schema.  ``--require-inputs`` is a separate, expensive check of
the frozen plan, source ZIP and both raw atlas matrices using the same
streaming aggregation validator as the analysis.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

try:  # Script execution from ``scripts/``.
    from analyze_ets2_program import (
        AUTHOR_CELL_TYPES,
        COMPARATOR_PROGRAM_IDS,
        FROZEN_PLAN_SHA256,
        PROGRAM_IDS,
        ROOT,
        _source_zip_from_manifest,
        aggregate_all_features,
        load_frozen_plan,
        sha256_file,
        verify_pinned_source,
    )
except ModuleNotFoundError:  # Package-style invocation from the repository root.
    from scripts.analyze_ets2_program import (
        AUTHOR_CELL_TYPES,
        COMPARATOR_PROGRAM_IDS,
        FROZEN_PLAN_SHA256,
        PROGRAM_IDS,
        ROOT,
        _source_zip_from_manifest,
        aggregate_all_features,
        load_frozen_plan,
        sha256_file,
        verify_pinned_source,
    )


PUBLIC_OUTPUT_BASENAMES: tuple[str, ...] = (
    "ets2-program-mapping.csv",
    "ets2-program-membership.csv",
    "ets2-program-matching-qc.csv",
    "ets2-program-context-summary.csv",
    "ets2-program-contrasts.csv",
    "ets2-program-correlations.csv",
)


def _fail(message: str) -> None:
    raise ValueError(message)


def _assert_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        _fail(f"{label} is missing columns: {', '.join(missing)}")


def _bool_values(values: pd.Series, label: str) -> np.ndarray:
    """Parse CSV booleans without treating the string ``False`` as truthy."""

    normalized = values.astype(str).str.strip().str.casefold()
    if not normalized.isin({"true", "false"}).all():
        _fail(f"{label} contains values other than true/false")
    return normalized.eq("true").to_numpy(dtype=bool)


def _recorded_paths(hashes: Mapping[str, Any], basenames: Sequence[str], label: str) -> dict[str, str]:
    """Resolve one recorded hash entry per required public basename."""

    candidates: dict[str, list[str]] = {basename: [] for basename in basenames}
    for name, digest in hashes.items():
        basename = Path(str(name)).name
        if basename in candidates:
            candidates[basename].append(str(name))
        if not isinstance(digest, str) or len(digest) != 64:
            _fail(f"{label} has an invalid SHA-256 entry for {name}")
    result: dict[str, str] = {}
    for basename, names in candidates.items():
        if len(names) != 1:
            _fail(f"{label} must contain exactly one entry for {basename}")
        result[basename] = names[0]
    return result


def _resolve_path(root: Path, recorded_name: str) -> Path:
    path = Path(recorded_name)
    return path if path.is_absolute() else root / path


def _verify_csv_hashes(root: Path, audit: Mapping[str, Any], *, require_cache: bool) -> dict[str, Any]:
    output_hashes = audit.get("outputs_sha256")
    if not isinstance(output_hashes, Mapping):
        _fail("ETS2 audit is missing outputs_sha256")
    public_entries = _recorded_paths(output_hashes, PUBLIC_OUTPUT_BASENAMES, "Public output hashes")
    checked_public: list[str] = []
    for basename, recorded_name in public_entries.items():
        path = _resolve_path(root, recorded_name)
        if not path.exists():
            _fail(f"Missing public ETS2 output: {path}")
        if sha256_file(path) != str(output_hashes[recorded_name]):
            _fail(f"Public ETS2 output hash changed: {path}")
        checked_public.append(str(recorded_name))

    ignored_hashes = audit.get("ignored_outputs_sha256", {})
    if not isinstance(ignored_hashes, Mapping):
        _fail("ETS2 audit is missing ignored_outputs_sha256")
    checked_cache: list[str] = []
    missing_cache: list[str] = []
    for recorded_name, expected_hash in sorted(ignored_hashes.items()):
        path = _resolve_path(root, str(recorded_name))
        if not path.exists():
            if require_cache:
                _fail(f"Missing required ignored ETS2 cache: {path}")
            missing_cache.append(str(recorded_name))
            continue
        if sha256_file(path) != str(expected_hash):
            _fail(f"Ignored ETS2 cache hash changed: {path}")
        checked_cache.append(str(recorded_name))
    return {"checked_public": checked_public, "checked_cache": checked_cache, "missing_cache": missing_cache}


def _verify_mapping(path: Path, datasets: Sequence[str]) -> dict[str, int]:
    frame = pd.read_csv(path)
    _assert_columns(
        frame,
        ["dataset_key", "program_id", "intended_unique_source_genes", "mapped_unique_genes", "mapping_fraction", "mapping_gate_pass"],
        "ETS2 mapping table",
    )
    if set(frame.dataset_key.astype(str)) != set(datasets) or set(frame.program_id.astype(str)) != set(PROGRAM_IDS):
        _fail("ETS2 mapping table does not cover both datasets and all nine programs")
    if frame.duplicated(["dataset_key", "program_id"]).any():
        _fail("ETS2 mapping table has duplicate dataset/program rows")
    for column in ("intended_unique_source_genes", "mapped_unique_genes"):
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all() or (values < 0).any():
            _fail(f"ETS2 mapping {column} contains invalid counts")
    fractions = pd.to_numeric(frame.mapping_fraction, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(fractions).all() or (fractions < 0).any() or (fractions > 1).any():
        _fail("ETS2 mapping fractions are outside [0, 1]")
    if (frame.mapped_unique_genes > frame.intended_unique_source_genes).any():
        _fail("ETS2 mapped genes exceed intended genes")
    return {str(key): int(value) for key, value in pd.Series(_bool_values(frame.mapping_gate_pass, "ETS2 mapping gate")).value_counts().items()}


def _verify_membership(path: Path, datasets: Sequence[str]) -> int:
    frame = pd.read_csv(path)
    _assert_columns(frame, ["dataset_key", "program_id", "source_value", "mapping_status"], "ETS2 membership table")
    if not set(frame.dataset_key.astype(str)).issubset(set(datasets)):
        _fail("ETS2 membership contains an unexpected dataset")
    if not set(frame.program_id.astype(str)).issubset(set(PROGRAM_IDS)):
        _fail("ETS2 membership contains an unexpected program")
    if frame.duplicated(["dataset_key", "program_id", "source_value"]).any():
        _fail("ETS2 membership has duplicate source-value rows")
    if frame.mapping_status.astype(str).str.contains("alias", case=False, regex=False).any():
        _fail("ETS2 membership contains an alias-based mapping status")
    return int(len(frame))


def _verify_matching(path: Path, datasets: Sequence[str], audit: Mapping[str, Any]) -> dict[str, int]:
    frame = pd.read_csv(path)
    _assert_columns(
        frame,
        ["dataset_key", "program_id", "bin", "bin_feature_count", "target_count", "control_pool_count", "required_controls_sufficient", "underflow", "matching_available", "control_pool_reference_mean"],
        "ETS2 matching QC table",
    )
    if set(frame.dataset_key.astype(str)) != set(datasets) or set(frame.program_id.astype(str)) != set(PROGRAM_IDS):
        _fail("ETS2 matching QC does not cover both datasets and all nine programs")
    if frame.duplicated(["dataset_key", "program_id", "bin"]).any():
        _fail("ETS2 matching QC has duplicate dataset/program/bin rows")
    for (dataset_key, program_id), group in frame.groupby(["dataset_key", "program_id"], sort=False):
        if set(pd.to_numeric(group.bin, errors="coerce").astype(int)) != set(range(20)):
            _fail(f"ETS2 matching QC does not contain exactly 20 bins for {dataset_key}/{program_id}")
        if (group.target_count < 0).any() or (group.control_pool_count < 0).any() or (group.bin_feature_count < 0).any():
            _fail(f"ETS2 matching QC has negative counts for {dataset_key}/{program_id}")
        required_sufficient = _bool_values(group.required_controls_sufficient, "ETS2 matching required_controls_sufficient")
        underflow = _bool_values(group.underflow, "ETS2 matching underflow")
        matching_available = _bool_values(group.matching_available, "ETS2 matching matching_available")
        if ((group.control_pool_count >= group.target_count).to_numpy(dtype=bool) != required_sufficient).any():
            _fail(f"ETS2 matching control sufficiency is inconsistent for {dataset_key}/{program_id}")
        if (underflow != (group.control_pool_count < group.target_count).to_numpy(dtype=bool)).any():
            _fail(f"ETS2 matching underflow is inconsistent for {dataset_key}/{program_id}")
        if matching_available.any() and underflow.any():
            _fail(f"ETS2 matching is available despite a control-pool underflow for {dataset_key}/{program_id}")
    expected_rows = int(audit.get("matching_qc_rows", -1))
    if expected_rows != len(frame):
        _fail("Audit matching_qc_rows does not match the public matching table")
    return {f"{dataset}/{program}": int(len(group)) for (dataset, program), group in frame.groupby(["dataset_key", "program_id"])}


def _verify_contexts(path: Path, datasets: Sequence[str], audit: Mapping[str, Any]) -> int:
    frame = pd.read_csv(path)
    _assert_columns(
        frame,
        [
            "dataset_key", "dataset_id", "modality", "suspension_type", "assay", "author_cell_type",
            "program_id", "threshold", "status", "raw_summary_status", "observed_donors", "eligible_donors",
            "score_eligible_donors", "mapped_program_genes", "median_program_score", "median_program_mean_log2_cpm1",
        ],
        "ETS2 context summary",
    )
    if set(frame.dataset_key.astype(str)) != set(datasets):
        _fail("ETS2 context summary does not cover both datasets")
    if set(frame.author_cell_type.astype(str)) != set(AUTHOR_CELL_TYPES):
        _fail("ETS2 context summary is missing a frozen author population")
    if set(frame.program_id.astype(str)) != set(PROGRAM_IDS):
        _fail("ETS2 context summary is missing a frozen program")
    if set(pd.to_numeric(frame.threshold, errors="coerce").dropna().astype(int)) != {20, 50}:
        _fail("ETS2 context summary contains a non-frozen threshold")
    keys = ["dataset_key", "dataset_id", "modality", "suspension_type", "assay", "author_cell_type", "program_id", "threshold"]
    if frame.duplicated(keys).any():
        _fail("ETS2 context summary has duplicate context/program/threshold rows")
    context_keys = ["dataset_key", "dataset_id", "modality", "suspension_type", "assay"]
    for context_key, group in frame.groupby(context_keys, sort=False):
        if set(group.author_cell_type.astype(str)) != set(AUTHOR_CELL_TYPES):
            _fail(f"ETS2 context summary lacks a frozen author row for {context_key}")
        for author, author_group in group.groupby("author_cell_type", sort=False):
            if set(author_group.program_id.astype(str)) != set(PROGRAM_IDS) or set(author_group.threshold.astype(int)) != {20, 50}:
                _fail(f"ETS2 context summary is not the full program/threshold Cartesian product for {context_key}/{author}")
    for column in ("observed_donors", "eligible_donors", "score_eligible_donors"):
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all() or (values < 0).any():
            _fail(f"ETS2 context summary {column} contains invalid counts")
    if (frame.score_eligible_donors > frame.eligible_donors).any() or (frame.eligible_donors > frame.observed_donors).any():
        _fail("ETS2 context donor gates are inconsistent")
    if int(audit.get("context_rows", -1)) != len(frame):
        _fail("Audit context_rows does not match the public context summary")
    return int(len(frame))


def _verify_contrasts(path: Path, audit: Mapping[str, Any]) -> dict[str, int]:
    frame = pd.read_csv(path)
    _assert_columns(
        frame,
        ["dataset_key", "program_id", "threshold", "left_population", "right_population", "status", "paired_donors", "left_observed_donors", "right_observed_donors", "left_score_eligible_donors", "right_score_eligible_donors"],
        "ETS2 paired contrast table",
    )
    if set(frame.program_id.astype(str)) - set(PROGRAM_IDS):
        _fail("ETS2 paired contrasts contain an unexpected program")
    if set(frame.left_population.astype(str)) - {"ActMac", "LAM-like", "Monocyte"} or set(frame.right_population.astype(str)) - {"Kupffer"}:
        _fail("ETS2 paired contrast contains an unexpected population")
    if frame.duplicated(["dataset_key", "modality", "suspension_type", "assay", "program_id", "threshold", "left_population"]).any():
        _fail("ETS2 paired contrasts have duplicate rows")
    numeric = ["paired_donors", "left_observed_donors", "right_observed_donors", "left_score_eligible_donors", "right_score_eligible_donors"]
    for column in numeric:
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all() or (values < 0).any():
            _fail(f"ETS2 contrast {column} contains invalid counts")
    if (frame.paired_donors > frame.left_score_eligible_donors).any() or (frame.paired_donors > frame.right_score_eligible_donors).any():
        _fail("ETS2 paired donors exceed score-eligible donors")
    counts = {str(key): int(value) for key, value in frame.status.value_counts().items()}
    if audit.get("contrast_status_counts") != counts:
        _fail("Audit contrast status counts do not match the public contrast table")
    return counts


def _verify_correlations(path: Path, audit: Mapping[str, Any]) -> dict[str, int]:
    frame = pd.read_csv(path)
    _assert_columns(frame, ["dataset_key", "comparison_id", "threshold", "author_cell_type", "status", "donors", "spearman_rho"], "ETS2 correlation table")
    # Correlation rows are emitted by comparison ID; the primary program is
    # implicit in the table contract and no population p-values are allowed.
    if "comparison_id" not in frame.columns or "comparison" not in frame.columns:
        _fail("ETS2 correlation table is missing comparison identifiers")
    if set(frame.author_cell_type.astype(str)) - set(AUTHOR_CELL_TYPES):
        _fail("ETS2 correlations contain an unexpected author population")
    if frame.duplicated(["dataset_key", "modality", "suspension_type", "author_cell_type", "threshold", "comparison_id"]).any():
        _fail("ETS2 correlation table has duplicate context/comparison rows")
    for column in frame.columns:
        if "pvalue" in column.casefold() or "p_value" in column.casefold():
            _fail("ETS2 correlation output must not contain population p-values")
    counts = {str(key): int(value) for key, value in frame.status.value_counts().items()}
    if audit.get("correlation_status_counts") != counts:
        _fail("Audit correlation status counts do not match the public correlation table")
    return counts


def _verify_cache_schema(root: Path, audit: Mapping[str, Any]) -> list[str]:
    ignored = audit.get("ignored_outputs_sha256", {})
    checked: list[str] = []
    for name in sorted(str(value) for value in ignored if str(value).endswith(".npz")):
        path = _resolve_path(root, name)
        with np.load(path, allow_pickle=False) as cache:
            required = {"counts", "feature_index", "feature_id", "donor_id", "group_id"}
            if not required.issubset(set(cache.files)):
                _fail(f"ETS2 all-feature cache is missing schema fields: {path}")
            counts = np.asarray(cache["counts"])
            feature_index = np.asarray(cache["feature_index"])
            feature_id = np.asarray(cache["feature_id"])
            donor_id = np.asarray(cache["donor_id"])
            group_id = np.asarray(cache["group_id"])
            if counts.ndim != 2 or counts.shape[1] != len(feature_index) or len(feature_index) != len(feature_id):
                _fail(f"ETS2 all-feature cache dimensions are inconsistent: {path}")
            if counts.shape[0] != len(donor_id) or len(donor_id) != len(group_id):
                _fail(f"ETS2 all-feature cache donor dimensions are inconsistent: {path}")
            if counts.dtype.kind not in "biuf" or (counts < 0).any() or not np.isfinite(counts).all():
                _fail(f"ETS2 all-feature cache contains invalid counts: {path}")
        checked.append(name)
    return checked


def _verify_inputs(root: Path, plan_path: Path, audit: Mapping[str, Any], chunk_size: int) -> list[dict[str, Any]]:
    plan, _program_manifest, _ = load_frozen_plan(plan_path, root=root)
    for relative, field in (
        ("scripts/analyze_ets2_program.py", "analysis_script_sha256"),
        ("scripts/analyze_liver_atlas.py", "atlas_helper_script_sha256"),
        ("requirements-ets2-program.txt", "requirements_sha256"),
    ):
        if audit.get(field) != sha256_file(root / relative):
            _fail(f"ETS2 provenance hash changed: {relative}")
    benchmark_manifest = json.loads((root / "config/ets2-benchmark-sources.json").read_text(encoding="utf-8"))
    source_zip = _source_zip_from_manifest(root, benchmark_manifest)
    recorded_source = audit.get("source_zip", {})
    if recorded_source.get("sha256") != sha256_file(source_zip) or int(recorded_source.get("bytes", -1)) != source_zip.stat().st_size:
        _fail("ETS2 audit source ZIP fingerprint does not match the pinned source")
    atlas_manifest = json.loads((root / "config/liver-atlas-sources.json").read_text(encoding="utf-8"))
    source_by_file = {str(item["file"]): item for item in atlas_manifest["sources"]}
    input_qc: list[dict[str, Any]] = []
    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - requirements pin anndata
        raise RuntimeError("anndata is required for --require-inputs") from exc
    for spec in plan["datasets"]:
        key = str(spec["key"])
        filename = str(spec["file"])
        source = source_by_file[filename]
        path = root / str(atlas_manifest["source_directory"]) / filename
        verify_pinned_source(path, source)
        data = ad.read_h5ad(path, backed="r")
        try:
            aggregate = aggregate_all_features(
                data,
                key,
                dataset_id=str(spec["dataset_id"]),
                version_id=str(spec.get("version_id", "")),
                source_file=str(path.relative_to(root)),
                expected_observations=int(spec["expected_observations"]),
                author_cell_types=plan["author_cell_types"],
                cohort=str(plan["cohort"]),
                thresholds=plan["cell_thresholds"],
                chunk_size=chunk_size,
            )
            input_qc.append({**aggregate.qc, "feature_ids": int(len(aggregate.feature_map))})
        finally:
            manager = getattr(data, "file", None)
            if manager is not None and hasattr(manager, "close"):
                manager.close()
    audit_by_key = {
        str(item.get("dataset_key")): item
        for item in audit.get("datasets", [])
        if isinstance(item, Mapping) and item.get("dataset_key")
    }
    for result in input_qc:
        expected = audit_by_key.get(str(result["dataset_key"]))
        if expected is None:
            _fail(f"ETS2 audit has no dataset QC row for {result['dataset_key']}")
        for field in (
            "raw_rows", "raw_features", "stored_nonzero_values", "total_cells",
            "retained_cells_all_metadata", "excluded_cells_all_metadata", "selected_psc_cells",
            "selected_donor_groups", "selected_donors",
        ):
            if int(expected.get(field, -1)) != int(result.get(field, -2)):
                _fail(f"ETS2 audit dataset QC changed for {result['dataset_key']}/{field}")
    return input_qc


def verify(
    root: Path = ROOT,
    *,
    audit_path: Path | None = None,
    plan_path: Path | None = None,
    require_cache: bool = False,
    require_inputs: bool = False,
    chunk_size: int = 8192,
) -> dict[str, Any]:
    """Verify public ETS2 outputs, with opt-in ignored-cache/full-input checks."""

    root = Path(root)
    audit_path = audit_path or (root / "reports/ets2-program-audit.json")
    if not audit_path.exists():
        _fail(f"Missing ETS2 audit report: {audit_path}")
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        _fail(f"ETS2 audit is not valid JSON: {exc}")
    if audit.get("analysis_id") != "psc-ets2-program-v1":
        _fail("Unexpected ETS2 audit analysis_id")
    if audit.get("plan_sha256") != FROZEN_PLAN_SHA256:
        _fail("ETS2 audit does not reference the frozen plan SHA-256")
    if not isinstance(audit.get("limitations"), str) or not audit["limitations"].strip():
        _fail("ETS2 audit limitations are missing or malformed")
    file_result = _verify_csv_hashes(root, audit, require_cache=require_cache)
    outputs = audit.get("outputs_sha256", {})
    public_paths = _recorded_paths(outputs, PUBLIC_OUTPUT_BASENAMES, "Public output hashes")
    dataset_keys = [str(item["dataset_key"]) for item in audit.get("datasets", []) if isinstance(item, Mapping) and "dataset_key" in item]
    dataset_keys = dataset_keys or ["sc", "sn"]
    if set(dataset_keys) != {"sc", "sn"}:
        _fail("ETS2 audit dataset list is not exactly sc and sn")
    mapping_counts = _verify_mapping(_resolve_path(root, public_paths["ets2-program-mapping.csv"]), dataset_keys)
    membership_rows = _verify_membership(_resolve_path(root, public_paths["ets2-program-membership.csv"]), dataset_keys)
    matching_counts = _verify_matching(_resolve_path(root, public_paths["ets2-program-matching-qc.csv"]), dataset_keys, audit)
    context_rows = _verify_contexts(_resolve_path(root, public_paths["ets2-program-context-summary.csv"]), dataset_keys, audit)
    contrast_counts = _verify_contrasts(_resolve_path(root, public_paths["ets2-program-contrasts.csv"]), audit)
    correlation_counts = _verify_correlations(_resolve_path(root, public_paths["ets2-program-correlations.csv"]), audit)
    input_qc: list[dict[str, Any]] = []
    if require_inputs:
        input_qc = _verify_inputs(root, plan_path or (root / "config/ets2-program-plan.json"), audit, chunk_size)
    if require_cache:
        cache_schema = _verify_cache_schema(root, audit)
    else:
        cache_schema = []
    return {
        "status": "all_ets2_program_outputs_verified",
        "mode": "full_inputs" if require_inputs else "audit_only",
        "checked_public_outputs": file_result["checked_public"],
        "checked_cache_files": file_result["checked_cache"],
        "missing_cache_files": file_result["missing_cache"],
        "cache_schema_checked": cache_schema,
        "mapping_gate_counts": mapping_counts,
        "membership_rows": membership_rows,
        "matching_groups": matching_counts,
        "context_rows": context_rows,
        "contrast_status_counts": contrast_counts,
        "correlation_status_counts": correlation_counts,
        "input_qc": input_qc,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--audit-path", type=Path, default=None)
    parser.add_argument("--plan", type=Path, default=None)
    parser.add_argument("--require-cache", action="store_true", help="Require and schema-check ignored count/cache files")
    parser.add_argument("--require-inputs", action="store_true", help="Recheck the pinned ZIP and both full raw H5AD inputs")
    parser.add_argument("--chunk-size", type=int, default=8192)
    args = parser.parse_args(argv)
    result = verify(
        args.root,
        audit_path=args.audit_path,
        plan_path=args.plan,
        require_cache=args.require_cache,
        require_inputs=args.require_inputs,
        chunk_size=args.chunk_size,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
