"""Stream a donor-aware descriptive analysis of the PSC liver atlas.

The analysis deliberately reads ``AnnData.raw.X`` in row chunks.  It never
uses the normalized ``AnnData.X`` matrix as counts and never materializes the
full expression matrix.  The public entry point writes small, deterministic
CSV tables; the functions below are also usable with synthetic AnnData-like
objects in tests.

This is a descriptive expression analysis.  Cells are retained only when
the frozen metadata rules permit a donor-level aggregation.  Repeated
libraries from the same donor are combined within dataset, assay and author
cell type.  Cell-level tests and causal or treatment interpretations are
outside this module.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
GENES: tuple[str, ...] = ("ETS2", "PRKD2", "UBASH3A", "PFKFB3", "BCL2L11")
REQUIRED_OBS_FIELDS: tuple[str, ...] = (
    "donor_id",
    "author_cell_type",
    "cell_type",
    "disease",
    "assay",
    "suspension_type",
    "library_uuid",
)
COHORTS: tuple[str, ...] = ("PSC", "Control", "PBC", "Other/unknown")
PRIMARY_MIN_CELLS = 20
SENSITIVITY_MIN_CELLS = 50
MIN_DONORS_FOR_CONTRAST = 3


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest for ``path``."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _clean(value: Any) -> str:
    """Convert a scalar metadata value to a stable, blank-aware string."""

    if value is None:
        return ""
    try:
        missing = pd.isna(value)
        if missing is pd.NA or (isinstance(missing, (bool, np.bool_)) and bool(missing)):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _is_missing(value: Any) -> bool:
    text = _clean(value)
    return not text or text.casefold() in {"nan", "none", "na", "n/a", "null"}


def canonical_cohort(value: Any) -> str:
    """Map the frozen disease labels to PSC, Control, PBC or Other/unknown."""

    text = _clean(value).casefold()
    compact = re.sub(r"[^a-z0-9]+", " ", text).strip()
    if "primary sclerosing cholangitis" in compact or compact in {"psc", "psc liver"}:
        return "PSC"
    if "primary biliary cholangitis" in compact or compact in {"pbc", "pbc liver"}:
        return "PBC"
    if compact == "normal" or compact.startswith("normal "):
        return "Control"
    if compact in {"control", "healthy", "healthy control", "reference", "uninjured"}:
        return "Control"
    return "Other/unknown"


def infer_modality(dataset_key: str, suspension_type: Any, assay: Any) -> str:
    """Infer the sc/sn reporting stratum from explicit metadata.

    The dataset key is only a final fallback for a missing/unclassified
    suspension label.  Known contradictory cell/nucleus labels are rejected
    by ``prepare_metadata`` so sc and sn cannot be merged silently.
    """

    suspension = _clean(suspension_type).casefold()
    assay_text = _clean(assay).casefold()
    if any(token in suspension for token in ("nucleus", "nuclei", "snrna", "single-nucleus", "single nucleus")):
        return "sn"
    if any(token in suspension for token in ("cell", "scrna", "single-cell", "single cell")):
        return "sc"
    if any(token in assay_text for token in ("nucleus", "nuclei", "snrna", "single-nucleus", "single nucleus")):
        return "sn"
    if any(token in assay_text for token in ("cell", "scrna", "single-cell", "single cell")):
        return "sc"
    if dataset_key in {"sc", "sn"}:
        return dataset_key
    return "unknown"


def _exclusion_reasons(row: Mapping[str, Any]) -> str:
    """Return all fixed metadata exclusion reasons for one cell."""

    reasons: list[str] = []
    required = ("donor_id", "author_cell_type", "assay", "disease")
    if any(_is_missing(row.get(field)) for field in required):
        reasons.append("missing_required_metadata")

    ontology = _clean(row.get("cell_type"))
    if not ontology or ontology.casefold() == "unknown":
        reasons.append("cell_type_unknown")

    author = _clean(row.get("author_cell_type"))
    if "doublet" in author.casefold():
        reasons.append("author_cell_type_doublet")
    if "--" in author:
        reasons.append("author_cell_type_double_hyphen")
    return ";".join(reasons)


def prepare_metadata(obs: pd.DataFrame, dataset_key: str, dataset_id: str | None = None) -> pd.DataFrame:
    """Normalize the frozen observation fields while preserving source labels.

    ``obs`` stays in its original row order.  Missing and excluded cells are
    retained in the returned ledger, but are marked ``retained=False`` and
    are never sent to count aggregation.
    """

    if not isinstance(obs, pd.DataFrame):
        obs = pd.DataFrame(obs)
    missing = [field for field in REQUIRED_OBS_FIELDS if field not in obs.columns]
    if missing:
        raise ValueError("AnnData obs is missing frozen metadata fields: " + ", ".join(missing))

    # Work with strings so categorical and nullable columns have identical
    # grouping behavior across pandas/anndata versions.
    result = pd.DataFrame(index=obs.index.copy())
    for field in REQUIRED_OBS_FIELDS:
        result[field] = obs[field].map(_clean)

    result["dataset_key"] = str(dataset_key)
    result["dataset_id"] = str(dataset_id or dataset_key)
    result["cohort"] = result["disease"].map(canonical_cohort)
    result["modality"] = [
        infer_modality(dataset_key, suspension, assay)
        for suspension, assay in zip(result["suspension_type"], result["assay"])
    ]
    expected_modality = dataset_key if dataset_key in {"sc", "sn"} else None
    contradictory = sorted(set(result.loc[result["modality"].isin({"sc", "sn"}), "modality"]))
    if expected_modality and contradictory and contradictory != [expected_modality]:
        raise ValueError(
            f"Dataset {dataset_key!r} contains contradictory sc/sn suspension labels: {contradictory}"
        )
    result["exclusion_reason"] = [
        _exclusion_reasons(row) for row in result[list(REQUIRED_OBS_FIELDS)].to_dict("records")
    ]
    result["retained"] = result["exclusion_reason"].eq("")
    result["cell_id"] = [str(value) for value in obs.index]
    result["suspension_type"] = result["suspension_type"].replace("", "__MISSING__")
    result["assay"] = result["assay"].replace("", "__MISSING__")
    result["author_cell_type"] = result["author_cell_type"].replace("", "__MISSING__")
    result["donor_id"] = result["donor_id"].replace("", "__MISSING__")
    result["ontology_cell_type"] = result["cell_type"].replace("", "__MISSING__")
    return result.reset_index(drop=True)


def _matrix_format(matrix: Any) -> str:
    if hasattr(matrix, "getformat"):
        try:
            return str(matrix.getformat())
        except Exception:  # pragma: no cover - unusual backed implementations
            pass
    value = getattr(matrix, "format", None)
    if value is not None:
        return str(value)
    try:
        from scipy import sparse

        if sparse.isspmatrix_csr(matrix):
            return "csr"
        if sparse.issparse(matrix):
            return str(matrix.getformat())
    except ImportError:  # pragma: no cover - scipy is a locked dependency
        pass
    return "dense"


def _is_sparse(matrix: Any) -> bool:
    try:
        from scipy import sparse

        return bool(sparse.issparse(matrix)) or _matrix_format(matrix) == "csr"
    except ImportError:  # pragma: no cover
        return _matrix_format(matrix) == "csr"


def _chunk(matrix: Any, start: int, stop: int) -> Any:
    value = matrix[start:stop]
    if _is_sparse(value):
        return value
    return np.asarray(value)


def _values_for_validation(value: Any) -> np.ndarray:
    if _is_sparse(value):
        return np.asarray(value.data)
    return np.asarray(value)


def _validate_numeric_values(values: np.ndarray, label: str) -> tuple[float | None, float | None, int]:
    """Validate one matrix chunk and return minimum, maximum and nonzeros."""

    array = np.asarray(values)
    if array.dtype.kind not in "biuf":
        raise ValueError(f"{label} contains unsupported dtype {array.dtype}")
    if array.size == 0:
        return None, None, 0
    if array.dtype.kind in "fc" and not np.isfinite(array).all():
        raise ValueError(f"{label} contains nonfinite values")
    if (array < 0).any():
        raise ValueError(f"{label} contains negative counts")
    if array.dtype.kind in "fc":
        if not np.equal(array, np.floor(array)).all():
            raise ValueError(f"{label} contains non-integer counts")
        if (array > np.iinfo(np.int64).max).any():
            raise ValueError(f"{label} contains counts outside int64 range")
    return float(np.min(array)), float(np.max(array)), int(np.count_nonzero(array))


def iter_validated_count_chunks(
    matrix: Any, chunk_size: int = 4096, label: str = "raw.X"
) -> Iterable[tuple[int, int, Any, dict[str, Any]]]:
    """Yield validated row chunks without densifying the complete matrix."""

    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    if not hasattr(matrix, "shape") or len(matrix.shape) != 2:
        raise ValueError(f"{label} is not a two-dimensional matrix")
    n_rows, _ = matrix.shape
    if _matrix_format(matrix) not in {"csr", "dense"}:
        raise ValueError(f"{label} must be CSR or dense, found {_matrix_format(matrix)!r}")
    for start in range(0, int(n_rows), chunk_size):
        stop = min(start + chunk_size, int(n_rows))
        value = _chunk(matrix, start, stop)
        if _is_sparse(value):
            # Backed CSR datasets return a CSR matrix per slice.  Check its
            # structure before using values; this catches malformed indices
            # without converting the slice to dense.
            if _matrix_format(value) != "csr":
                raise ValueError(f"{label} slice is not CSR")
            if hasattr(value, "check_format"):
                try:
                    value.check_format(full_check=True)
                except Exception as exc:
                    raise ValueError(f"{label} has invalid CSR indices") from exc
            if len(value.indptr) != (stop - start) + 1:
                raise ValueError(f"{label} has invalid CSR indptr length")
            if len(value.data) != int(value.indptr[-1]):
                raise ValueError(f"{label} has invalid CSR indptr/data length")
            if len(value.indices) and (
                int(np.min(value.indices)) < 0 or int(np.max(value.indices)) >= int(matrix.shape[1])
            ):
                raise ValueError(f"{label} has an out-of-range CSR index")
            # Unsorted CSR indices are valid storage.  Aggregation follows
            # the matrix multiplication semantics and therefore correctly
            # sums repeated entries in a row; structural/index checks above
            # are the relevant integrity checks here.
        values = _values_for_validation(value)
        minimum, maximum, nonzero = _validate_numeric_values(values, label)
        yield start, stop, value, {"min": minimum, "max": maximum, "nonzero": nonzero}


def validate_counts_matrix(matrix: Any, chunk_size: int = 4096, label: str = "raw.X") -> dict[str, Any]:
    """Validate all stored count values and CSR structure in row chunks."""

    n_rows, n_cols = getattr(matrix, "shape", (None, None))
    if n_rows is None or n_cols is None:
        raise ValueError(f"{label} has no matrix shape")
    minimums: list[float] = []
    maximums: list[float] = []
    nonzero = 0
    chunks = 0
    for _, _, _, stats in iter_validated_count_chunks(matrix, chunk_size, label):
        chunks += 1
        nonzero += stats["nonzero"]
        if stats["min"] is not None:
            minimums.append(stats["min"])
            maximums.append(stats["max"])
    return {
        "rows": int(n_rows),
        "features": int(n_cols),
        "chunks": chunks,
        "format": _matrix_format(matrix),
        "dtype": str(getattr(matrix, "dtype", "unknown")),
        "minimum_stored_value": min(minimums) if minimums else 0.0,
        "maximum_stored_value": max(maximums) if maximums else 0.0,
        "stored_nonzero_values": nonzero,
    }


def resolve_raw_counts(adata: Any) -> tuple[Any, pd.DataFrame, str]:
    """Return the frozen raw count matrix and its feature table.

    The plan explicitly requires ``raw/X`` and ``raw/var``.  Layers and
    normalized ``X`` are intentionally not accepted as a fallback.
    """

    raw = getattr(adata, "raw", None)
    if raw is None:
        raise ValueError("AnnData.raw is missing; normalized X cannot be used as counts")
    matrix = getattr(raw, "X", None)
    var = getattr(raw, "var", None)
    if matrix is None or var is None:
        raise ValueError("AnnData.raw must contain both X and var")
    if matrix is getattr(adata, "X", None):
        raise ValueError("AnnData.raw.X is the same object as normalized AnnData.X")
    if _matrix_format(matrix) not in {"csr", "dense"}:
        raise ValueError(f"AnnData.raw.X must be CSR or dense, found {_matrix_format(matrix)!r}")
    return matrix, pd.DataFrame(var).copy(), "raw.X"


def _feature_map(raw_var: pd.DataFrame, n_features: int) -> pd.DataFrame:
    """Build the exact ``feature_name`` map required by the frozen plan."""

    if len(raw_var) != n_features:
        raise ValueError("raw/var and raw/X feature dimensions disagree")
    if "feature_name" not in raw_var.columns:
        raise ValueError("raw/var is missing the frozen feature_name column")
    feature_ids = pd.Index([_clean(value) for value in raw_var.index])
    if feature_ids.has_duplicates or any(not value for value in feature_ids):
        raise ValueError("raw/var feature identifiers are missing or duplicated")
    names = [_clean(value) for value in raw_var["feature_name"]]
    result = pd.DataFrame({"feature_index": np.arange(n_features, dtype=np.int64),
                           "feature_id": feature_ids.astype(str),
                           "feature_name": names})
    return result


def resolve_target_features(raw_var: pd.DataFrame, n_features: int, genes: Sequence[str] = GENES) -> pd.DataFrame:
    """Return one audit row per target gene, including missing/ambiguous genes."""

    features = _feature_map(raw_var, n_features)
    records: list[dict[str, Any]] = []
    for gene in genes:
        matches = features.index[features["feature_name"].eq(gene)].tolist()
        if len(matches) == 1:
            match = features.iloc[matches[0]]
            records.append({"gene": gene, "status": "present", "feature_index": int(match.feature_index),
                            "feature_id": match.feature_id, "feature_name": match.feature_name,
                            "match_count": 1})
        elif len(matches) == 0:
            records.append({"gene": gene, "status": "missing_feature_name", "feature_index": None,
                            "feature_id": None, "feature_name": gene, "match_count": 0})
        else:
            records.append({"gene": gene, "status": "ambiguous_feature_name", "feature_index": None,
                            "feature_id": ";".join(features.iloc[matches].feature_id),
                            "feature_name": gene, "match_count": len(matches)})
    return pd.DataFrame(records)


def _stable_unique(values: Iterable[Any]) -> list[str]:
    return sorted({_clean(value) for value in values if not _is_missing(value)})


def _group_columns() -> list[str]:
    return ["dataset_key", "dataset_id", "modality", "suspension_type", "assay", "donor_id", "cohort", "author_cell_type"]


def _stratum_columns() -> list[str]:
    return ["dataset_key", "dataset_id", "modality", "suspension_type", "assay", "author_cell_type"]


@dataclass
class DatasetResult:
    dataset_key: str
    pseudobulk: pd.DataFrame
    cell_accounting: pd.DataFrame
    gene_audit: pd.DataFrame
    qc: dict[str, Any]


def aggregate_adata(
    adata: Any,
    dataset_key: str,
    *,
    dataset_id: str | None = None,
    dataset_version_id: str | None = None,
    source_file: str | None = None,
    expected_observations: int | None = None,
    chunk_size: int = 4096,
    genes: Sequence[str] = GENES,
) -> DatasetResult:
    """Aggregate raw counts by donor, assay and original author cell type."""

    n_obs = int(getattr(adata, "n_obs", getattr(adata, "shape", (0,))[0]))
    n_vars = int(getattr(adata, "n_vars", getattr(adata, "shape", (0, 0))[1]))
    if expected_observations is not None and n_obs != int(expected_observations):
        raise ValueError(f"{dataset_key} has {n_obs} observations; expected {expected_observations}")
    obs = pd.DataFrame(getattr(adata, "obs")).copy()
    if len(obs) != n_obs:
        raise ValueError("AnnData obs row count differs from n_obs")
    obs_names = pd.Index([str(value) for value in getattr(adata, "obs_names", obs.index)])
    if obs_names.has_duplicates:
        raise ValueError(f"{dataset_key} has duplicate cell IDs")
    obs.index = obs_names
    metadata = prepare_metadata(obs, dataset_key, dataset_id)
    raw_matrix, raw_var, count_source = resolve_raw_counts(adata)
    if raw_matrix.shape[0] != n_obs:
        raise ValueError("raw/X and obs dimensions disagree")
    if raw_matrix.shape[1] != n_vars and getattr(getattr(adata, "raw", None), "n_vars", raw_matrix.shape[1]) != raw_matrix.shape[1]:
        raise ValueError("raw/X feature dimension is inconsistent")
    normalized_x = getattr(adata, "X", None)
    normalized_shape = list(getattr(normalized_x, "shape", ())) if normalized_x is not None else []
    if normalized_x is None or not normalized_shape or normalized_shape[0] != n_obs:
        raise ValueError("normalized AnnData.X is missing or has the wrong observation dimension")
    if len(normalized_shape) != 2 or normalized_shape[1] != int(raw_matrix.shape[1]):
        raise ValueError("raw/X feature dimensions disagree")

    feature_audit = resolve_target_features(raw_var, int(raw_matrix.shape[1]), genes)
    target_indices = [
        int(value) if pd.notna(value) and status == "present" else None
        for value, status in zip(feature_audit["feature_index"], feature_audit["status"])
    ]

    group_columns = _group_columns()
    retained = metadata[metadata["retained"]].copy()
    if len(retained):
        # A missing donor is excluded by prepare_metadata; this assertion
        # prevents a future rule change from pooling anonymous cells.
        if (retained["donor_id"] == "__MISSING__").any():
            raise ValueError("Retained cells contain a missing donor ID")
        keys = [tuple(row[column] for column in group_columns) for row in retained.to_dict("records")]
    else:
        keys = []
    unique_keys = sorted(set(keys), key=lambda value: tuple(str(item) for item in value))
    code_by_key = {key: code for code, key in enumerate(unique_keys)}
    retained_codes = np.asarray([code_by_key[key] for key in keys], dtype=np.int64)
    retained_row_indices = retained.index.to_numpy(dtype=np.int64)

    n_groups = len(unique_keys)
    target_umis = np.zeros((n_groups, len(genes)), dtype=np.int64)
    detected_cells = np.zeros((n_groups, len(genes)), dtype=np.int64)
    all_feature_umis = np.zeros(n_groups, dtype=np.int64)
    cell_counts = np.zeros(n_groups, dtype=np.int64)
    qc = validate_counts_matrix_header(raw_matrix, n_obs, chunk_size, count_source)
    validation_chunks = 0
    minimums: list[float] = []
    maximums: list[float] = []
    stored_nonzero = 0
    row_cursor = 0
    for start, stop, chunk, stats in iter_validated_count_chunks(raw_matrix, chunk_size, count_source):
        validation_chunks += 1
        stored_nonzero += stats["nonzero"]
        if stats["min"] is not None:
            minimums.append(stats["min"])
            maximums.append(stats["max"])
        row_codes = np.full(stop - start, -1, dtype=np.int64)
        if len(retained_row_indices):
            left = np.searchsorted(retained_row_indices, start, side="left")
            right = np.searchsorted(retained_row_indices, stop, side="left")
            if right > left:
                row_codes[retained_row_indices[left:right] - start] = retained_codes[left:right]
        if not np.any(row_codes >= 0):
            continue
        if _is_sparse(chunk):
            row_totals = np.asarray(chunk.sum(axis=1)).ravel()
        else:
            row_totals = np.asarray(chunk).sum(axis=1)
        row_totals = np.asarray(row_totals, dtype=np.int64)
        valid_rows = row_codes >= 0
        codes = row_codes[valid_rows]
        np.add.at(cell_counts, codes, 1)
        np.add.at(all_feature_umis, codes, row_totals[valid_rows])
        values_by_gene: list[np.ndarray] = []
        for index in target_indices:
            if index is None:
                values_by_gene.append(np.zeros(stop - start, dtype=np.int64))
            elif _is_sparse(chunk):
                values_by_gene.append(np.asarray(chunk[:, index].toarray()).ravel().astype(np.int64, copy=False))
            else:
                values_by_gene.append(np.asarray(chunk)[:, index].astype(np.int64, copy=False))
        for gene_index, values in enumerate(values_by_gene):
            values = values[valid_rows]
            np.add.at(target_umis[:, gene_index], codes, values)
            np.add.at(detected_cells[:, gene_index], codes, values > 0)
        row_cursor = stop

    if row_cursor < n_obs and validation_chunks == 0 and n_obs:
        raise ValueError("raw/X yielded no chunks")
    qc.update({
        "rows": n_obs,
        "features": int(raw_matrix.shape[1]),
        "chunks": validation_chunks,
        "format": _matrix_format(raw_matrix),
        "dtype": str(getattr(raw_matrix, "dtype", "unknown")),
        "minimum_stored_value": min(minimums) if minimums else 0.0,
        "maximum_stored_value": max(maximums) if maximums else 0.0,
        "stored_nonzero_values": stored_nonzero,
        "count_source": count_source,
        "normalized_matrix": "AnnData.X (validated as a separate source; never aggregated)",
        "normalized_matrix_shape": normalized_shape,
    })

    # Inventory includes excluded cells and is the ledger used to make
    # missingness visible in summary tables.
    accounting = _cell_accounting(metadata)
    records: list[dict[str, Any]] = []
    gene_status = dict(zip(feature_audit["gene"], feature_audit["status"]))
    gene_ids = dict(zip(feature_audit["gene"], feature_audit["feature_id"]))
    for group_index, key in enumerate(unique_keys):
        values = dict(zip(group_columns, key))
        total = int(all_feature_umis[group_index])
        cells = int(cell_counts[group_index])
        row = {
            **values,
            "dataset_version_id": dataset_version_id or "",
            "source_file": source_file or "",
            "group_id": "|".join(str(value) for value in key),
            "cell_count": cells,
            "all_feature_umis": total,
            "positive_all_feature_umi_denominator": total > 0,
            "minimum_cells_primary": PRIMARY_MIN_CELLS,
            "minimum_cells_sensitivity": SENSITIVITY_MIN_CELLS,
            "passes_primary_cell_gate": cells >= PRIMARY_MIN_CELLS,
            "passes_sensitivity_cell_gate": cells >= SENSITIVITY_MIN_CELLS,
            "valid_primary_group": cells >= PRIMARY_MIN_CELLS and total > 0,
            "valid_sensitivity_group": cells >= SENSITIVITY_MIN_CELLS and total > 0,
        }
        group_meta = retained[retained_codes == group_index]
        row["library_ids"] = ";".join(_stable_unique(group_meta["library_uuid"]))
        row["library_count"] = len(_stable_unique(group_meta["library_uuid"]))
        row["ontology_cell_types"] = ";".join(_stable_unique(group_meta["cell_type"]))
        row["disease_labels"] = ";".join(_stable_unique(group_meta["disease"]))
        for gene_index, gene in enumerate(genes):
            present = gene_status[gene] == "present"
            umi = int(target_umis[group_index, gene_index])
            detected = int(detected_cells[group_index, gene_index])
            cpm = float(umi / total * 1_000_000) if present and total > 0 else math.nan
            row[f"{gene}_status"] = gene_status[gene]
            row[f"{gene}_feature_id"] = gene_ids[gene] or ""
            row[f"{gene}_feature_present"] = present
            row[f"{gene}_umis"] = umi if present else 0
            row[f"{gene}_detected_cells"] = detected if present else 0
            row[f"{gene}_detection_fraction"] = detected / cells if present and cells else (0.0 if present else math.nan)
            row[f"{gene}_cpm"] = cpm
            row[f"{gene}_log2_cpm1"] = float(np.log2(cpm + 1.0)) if np.isfinite(cpm) else math.nan
        records.append(row)
    pseudobulk = pd.DataFrame(records)
    if not pseudobulk.empty:
        pseudobulk = pseudobulk.sort_values(group_columns).reset_index(drop=True)
    qc["total_cells"] = n_obs
    qc["retained_cells"] = int(metadata["retained"].sum())
    qc["excluded_cells"] = int((~metadata["retained"]).sum())
    qc["retained_donor_groups"] = len(pseudobulk)
    qc["excluded_by_reason"] = {
        str(key): int(value) for key, value in metadata.loc[~metadata["retained"], "exclusion_reason"].value_counts().items()
    }
    return DatasetResult(dataset_key, pseudobulk, accounting, feature_audit.assign(dataset_key=dataset_key), qc)


def validate_counts_matrix_header(matrix: Any, expected_rows: int, chunk_size: int, label: str) -> dict[str, Any]:
    """Validate dimensions before the streaming pass begins."""

    if int(matrix.shape[0]) != int(expected_rows):
        raise ValueError(f"{label} and obs dimensions disagree")
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    return {}


def _cell_accounting(metadata: pd.DataFrame) -> pd.DataFrame:
    columns = ["dataset_key", "dataset_id", "modality", "suspension_type", "assay", "cohort", "author_cell_type", "ontology_cell_type", "exclusion_reason"]
    records: list[dict[str, Any]] = []
    for key, group in metadata.groupby(columns, dropna=False, sort=True):
        values = dict(zip(columns, key if isinstance(key, tuple) else (key,)))
        donor_values = _stable_unique(group["donor_id"])
        library_values = _stable_unique(group["library_uuid"])
        records.append({**values,
                        "cell_count": int(len(group)),
                        "retained_cell_count": int(group["retained"].sum()),
                        "excluded_cell_count": int((~group["retained"]).sum()),
                        "donor_count_observed": len(donor_values),
                        "library_count_observed": len(library_values),
                        "donor_ids_observed": ";".join(donor_values),
                        "library_ids_observed": ";".join(library_values)})
    if not records:
        return pd.DataFrame(columns=columns + ["cell_count", "retained_cell_count", "excluded_cell_count", "donor_count_observed", "library_count_observed", "donor_ids_observed", "library_ids_observed"])
    return pd.DataFrame(records).sort_values(columns).reset_index(drop=True)


def _gene_long_values(row: Mapping[str, Any], gene: str) -> tuple[float, float, int, int, str]:
    status = str(row[f"{gene}_status"])
    if status != "present":
        return math.nan, math.nan, 0, 0, status
    return (float(row[f"{gene}_cpm"]) if pd.notna(row[f"{gene}_cpm"]) else math.nan,
            float(row[f"{gene}_log2_cpm1"]) if pd.notna(row[f"{gene}_log2_cpm1"])
            else math.nan,
            int(row[f"{gene}_umis"]), int(row[f"{gene}_detected_cells"]), status)


def _summary_row(
    groups: pd.DataFrame,
    all_groups: pd.DataFrame,
    base: Mapping[str, Any],
    gene: str,
    threshold: int,
    gene_status_override: str | None = None,
) -> dict[str, Any]:
    present = bool(len(groups) and groups[f"{gene}_status"].eq("present").all())
    all_present_status = (all_groups[f"{gene}_status"].iloc[0] if len(all_groups)
                          else gene_status_override or "missing_feature_name")
    candidate = groups[groups[f"{gene}_status"].eq("present")].copy() if len(groups) else groups
    logs = candidate[f"{gene}_log2_cpm1"].astype(float).to_numpy() if len(candidate) else np.asarray([])
    cpms = candidate[f"{gene}_cpm"].astype(float).to_numpy() if len(candidate) else np.asarray([])
    detections = candidate[f"{gene}_detection_fraction"].astype(float).to_numpy() if len(candidate) else np.asarray([])
    umis = candidate[f"{gene}_umis"].astype(np.int64).to_numpy() if len(candidate) else np.asarray([], dtype=np.int64)
    if all_present_status != "present":
        status = all_present_status
    elif not len(groups):
        status = "no_eligible_donors"
    else:
        status = "evaluable_descriptive"
    row = {**base, "gene": gene, "cell_threshold": threshold,
           "status": status, "gene_status": all_present_status,
           "gene_feature_present": all_present_status == "present",
           "eligible_donors": int(candidate["donor_id"].nunique()) if len(candidate) else 0,
           "observed_donors": int(all_groups["donor_id"].nunique()) if len(all_groups) else 0,
           "retained_cells": int(candidate["cell_count"].sum()) if len(candidate) else 0,
           "observed_retained_cells": int(all_groups["cell_count"].sum()) if len(all_groups) else 0,
           "donors_with_gene_umis_gt0": int((umis > 0).sum()) if len(umis) else 0,
           "donor_weighted_median_log2_cpm1": float(np.median(logs)) if len(logs) else math.nan,
           "donor_weighted_median_cpm": float(np.median(cpms)) if len(cpms) else math.nan,
           "donor_weighted_median_detection_fraction": float(np.median(detections)) if len(detections) else math.nan,
           "observed_detection_fraction": float(candidate[f"{gene}_detected_cells"].sum() / candidate["cell_count"].sum()) if len(candidate) and candidate["cell_count"].sum() else math.nan,
           "minimum_log2_cpm1": float(np.min(logs)) if len(logs) else math.nan,
           "maximum_log2_cpm1": float(np.max(logs)) if len(logs) else math.nan,
           "all_feature_umis": int(candidate["all_feature_umis"].sum()) if len(candidate) else 0,
           "positive_denominator_donors": int(candidate["positive_all_feature_umi_denominator"].sum()) if len(candidate) else 0}
    return row


def summarize_expression(
    pseudobulk: pd.DataFrame,
    cell_accounting: pd.DataFrame | None = None,
    *,
    genes: Sequence[str] = GENES,
    thresholds: Sequence[int] = (PRIMARY_MIN_CELLS, SENSITIVITY_MIN_CELLS),
    gene_audit: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build condition-level donor summaries for every observed stratum."""

    if cell_accounting is None:
        cell_accounting = pd.DataFrame()
    strata_columns = _stratum_columns()
    if len(cell_accounting):
        inventory = cell_accounting[strata_columns + ["cohort"]].drop_duplicates()
    elif len(pseudobulk):
        inventory = pseudobulk[strata_columns + ["cohort"]].drop_duplicates()
    else:
        return pd.DataFrame()
    gene_status_by_gene: dict[str, str] = {}
    if gene_audit is not None and len(gene_audit):
        gene_status_by_gene = dict(zip(gene_audit["gene"], gene_audit["status"]))
    elif len(pseudobulk):
        for gene in genes:
            column = f"{gene}_status"
            if column in pseudobulk:
                gene_status_by_gene[gene] = str(pseudobulk[column].iloc[0])
    records: list[dict[str, Any]] = []
    for stratum, inventory_group in inventory.groupby(strata_columns, sort=True, dropna=False):
        base = dict(zip(strata_columns, stratum if isinstance(stratum, tuple) else (stratum,)))
        all_for_stratum = pseudobulk
        for column, value in base.items():
            all_for_stratum = all_for_stratum[all_for_stratum[column].eq(value)]
        for cohort in sorted(inventory_group["cohort"].unique(), key=lambda value: COHORTS.index(value) if value in COHORTS else len(COHORTS)):
            all_groups = all_for_stratum[all_for_stratum["cohort"].eq(cohort)]
            for threshold in thresholds:
                eligible = all_groups[(all_groups["cell_count"] >= threshold) & all_groups["positive_all_feature_umi_denominator"]]
                for gene in genes:
                    records.append(_summary_row(eligible, all_groups, {**base, "cohort": cohort}, gene, int(threshold),
                                                gene_status_by_gene.get(gene)))
    result = pd.DataFrame(records)
    if len(result):
        result = result.sort_values(strata_columns + ["cohort", "cell_threshold", "gene"]).reset_index(drop=True)
    return result


def _filter_groups(pseudobulk: pd.DataFrame, base: Mapping[str, Any], assay: str, cohort: str, threshold: int) -> pd.DataFrame:
    subset = pseudobulk
    for column, value in base.items():
        subset = subset[subset[column].eq(value)]
    return subset[(subset["assay"].eq(assay)) & subset["cohort"].eq(cohort) &
                  (subset["cell_count"] >= threshold) & subset["positive_all_feature_umi_denominator"]]


def _contrast_values(groups: pd.DataFrame, gene: str) -> np.ndarray:
    values = groups.loc[groups[f"{gene}_status"].eq("present"), f"{gene}_log2_cpm1"].astype(float).to_numpy()
    return values[np.isfinite(values)]


def _leave_one_out_differences(psc: pd.DataFrame, control: pd.DataFrame, gene: str,
                               min_donors: int = MIN_DONORS_FOR_CONTRAST) -> tuple[float | None, float | None, int]:
    values: list[float] = []
    if len(psc) - 1 >= min_donors and len(control) >= min_donors:
        for index in range(len(psc)):
            left = _contrast_values(psc.drop(psc.index[index]), gene)
            right = _contrast_values(control, gene)
            if len(left) >= min_donors and len(right) >= min_donors:
                values.append(float(np.median(left) - np.median(right)))
    if len(control) - 1 >= min_donors and len(psc) >= min_donors:
        for index in range(len(control)):
            left = _contrast_values(psc, gene)
            right = _contrast_values(control.drop(control.index[index]), gene)
            if len(left) >= min_donors and len(right) >= min_donors:
                values.append(float(np.median(left) - np.median(right)))
    return (min(values) if values else None, max(values) if values else None, len(values))


def build_contrasts(
    pseudobulk: pd.DataFrame,
    cell_accounting: pd.DataFrame | None = None,
    *,
    genes: Sequence[str] = GENES,
    thresholds: Sequence[int] = (PRIMARY_MIN_CELLS, SENSITIVITY_MIN_CELLS),
    min_donors: int = MIN_DONORS_FOR_CONTRAST,
) -> pd.DataFrame:
    """Build assay-matched descriptive PSC-Control contrasts.

    When PSC and Control cells occur only in different assay labels, one
    blocked row per gene and reporting stratum is retained.  This makes the
    SC chemistry confounding visible instead of allowing a cross-chemistry
    fallback.  SN rows that pass numerical gates retain an explicit residual
    processing-confounding flag supplied by the metadata audit.
    """

    if cell_accounting is None:
        cell_accounting = pd.DataFrame()
    strata_columns = ["dataset_key", "dataset_id", "modality", "suspension_type", "author_cell_type"]
    if len(cell_accounting):
        inventory = cell_accounting[strata_columns + ["assay", "cohort"]].drop_duplicates()
    elif len(pseudobulk):
        inventory = pseudobulk[strata_columns + ["assay", "cohort"]].drop_duplicates()
    else:
        return pd.DataFrame()
    gene_status_by_gene: dict[str, str] = {}
    if len(pseudobulk):
        for gene in genes:
            column = f"{gene}_status"
            if column in pseudobulk:
                gene_status_by_gene[gene] = str(pseudobulk[column].iloc[0])
    records: list[dict[str, Any]] = []
    for stratum, inv in inventory.groupby(strata_columns, sort=True, dropna=False):
        base = dict(zip(strata_columns, stratum if isinstance(stratum, tuple) else (stratum,)))
        psc_assays = sorted(inv.loc[inv.cohort.eq("PSC"), "assay"].unique())
        control_assays = sorted(inv.loc[inv.cohort.eq("Control"), "assay"].unique())
        shared_assays = sorted(set(psc_assays) & set(control_assays))
        assays = shared_assays if shared_assays else ["__NO_SHARED_ASSAY__"]
        for assay in assays:
            for threshold in thresholds:
                for gene in genes:
                    row: dict[str, Any] = {**base, "assay": assay, "gene": gene, "cell_threshold": int(threshold),
                                           "comparison": "PSC_minus_Control", "psc_assays_observed": ";".join(psc_assays),
                                           "control_assays_observed": ";".join(control_assays),
                                           "shared_assay": assay != "__NO_SHARED_ASSAY__",
                                           "min_donors_required": int(min_donors),
                                           "psc_median_log2_cpm1": math.nan, "control_median_log2_cpm1": math.nan,
                                           "psc_minus_control_median_difference": math.nan,
                                           "psc_median_detection_fraction": math.nan, "control_median_detection_fraction": math.nan,
                                           "psc_observed_detection_fraction": math.nan, "control_observed_detection_fraction": math.nan,
                                           "psc_eligible_donors": 0, "control_eligible_donors": 0,
                                           "leave_one_donor_out_min_difference": math.nan,
                                           "leave_one_donor_out_max_difference": math.nan,
                                           "leave_one_donor_out_evaluable_count": 0,
                                           "biological_interpretation": "unavailable"}
                    if not psc_assays or not control_assays:
                        row.update(status="missing_cohort", reason="PSC or Control assay stratum is absent")
                        records.append(row)
                        continue
                    if not shared_assays:
                        row.update(status="blocked_cross_chemistry", reason="PSC and Control have no exact shared assay; cross-chemistry fallback is prohibited",
                                    biological_interpretation="unavailable_due_to_assay_chemistry_confounding")
                        records.append(row)
                        continue
                    psc = _filter_groups(pseudobulk, base, assay, "PSC", int(threshold))
                    control = _filter_groups(pseudobulk, base, assay, "Control", int(threshold))
                    row["psc_eligible_donors"] = int(psc["donor_id"].nunique())
                    row["control_eligible_donors"] = int(control["donor_id"].nunique())
                    psc_values = _contrast_values(psc, gene)
                    control_values = _contrast_values(control, gene)
                    status_gene = "present"
                    if len(pseudobulk):
                        candidates = pseudobulk
                        for column, value in base.items():
                            candidates = candidates[candidates[column].eq(value)]
                        candidates = candidates[candidates.assay.eq(assay)]
                        if len(candidates) and candidates[f"{gene}_status"].iloc[0] != "present":
                            status_gene = str(candidates[f"{gene}_status"].iloc[0])
                        elif not len(candidates):
                            status_gene = gene_status_by_gene.get(gene, "missing_feature_name")
                    if status_gene != "present":
                        row.update(status="gene_missing", reason=status_gene)
                        records.append(row)
                        continue
                    if len(psc_values) < min_donors or len(control_values) < min_donors:
                        row.update(status="insufficient_donors", reason="At least three eligible donors per cohort are required")
                        records.append(row)
                        continue
                    difference = float(np.median(psc_values) - np.median(control_values))
                    loo_min, loo_max, loo_count = _leave_one_out_differences(psc, control, gene, min_donors)
                    row.update(status="evaluable_descriptive_residual_processing_confounded" if base["modality"] == "sn" else "evaluable_descriptive",
                               biological_interpretation="residual_processing_confounding_flagged_for_sn" if base["modality"] == "sn" else "descriptive_only",
                               psc_median_log2_cpm1=float(np.median(psc_values)),
                               control_median_log2_cpm1=float(np.median(control_values)),
                               psc_minus_control_median_difference=difference,
                               psc_median_detection_fraction=float(np.median(psc[f"{gene}_detection_fraction"])),
                               control_median_detection_fraction=float(np.median(control[f"{gene}_detection_fraction"])),
                               psc_observed_detection_fraction=float(psc[f"{gene}_detected_cells"].sum() / psc["cell_count"].sum()),
                               control_observed_detection_fraction=float(control[f"{gene}_detected_cells"].sum() / control["cell_count"].sum()),
                               leave_one_donor_out_min_difference=loo_min if loo_min is not None else math.nan,
                               leave_one_donor_out_max_difference=loo_max if loo_max is not None else math.nan,
                               leave_one_donor_out_evaluable_count=loo_count)
                    records.append(row)
    result = pd.DataFrame(records)
    if len(result):
        result = result.sort_values(strata_columns + ["assay", "cell_threshold", "gene"]).reset_index(drop=True)
    return result


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_frozen_plan(plan_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = _read_json(plan_path)
    manifest_path = (plan_path.parent.parent / plan["source_manifest"]).resolve() if not Path(plan["source_manifest"]).is_absolute() else Path(plan["source_manifest"])
    if not manifest_path.exists():
        raise ValueError(f"Missing liver-atlas source manifest: {manifest_path}")
    manifest = _read_json(manifest_path)
    expected_manifest_hash = plan.get("source_manifest_sha256")
    actual_manifest_hash = sha256_file(manifest_path)
    if expected_manifest_hash and actual_manifest_hash != expected_manifest_hash:
        raise ValueError("Frozen liver-atlas source manifest changed")
    if tuple(plan.get("genes", [])) != tuple(GENES):
        raise ValueError("Frozen liver-atlas gene list differs from the five fixed candidates")
    return plan, manifest


def validate_frozen_thresholds(plan: Mapping[str, Any]) -> dict[str, int]:
    """Require the run to use the exact cell and donor gates in the plan."""

    frozen = {
        "primary_min_cells": int(plan.get("primary_min_cells_per_donor_cell_type", -1)),
        "sensitivity_min_cells": int(plan.get("sensitivity_min_cells_per_donor_cell_type", -1)),
        "minimum_donors_for_contrast": int(plan.get("minimum_donors_per_cohort_for_contrast", -1)),
    }
    expected = {"primary_min_cells": PRIMARY_MIN_CELLS,
                "sensitivity_min_cells": SENSITIVITY_MIN_CELLS,
                "minimum_donors_for_contrast": MIN_DONORS_FOR_CONTRAST}
    if frozen != expected:
        raise ValueError(f"Frozen liver-atlas thresholds changed: {frozen}")
    return frozen


def verify_frozen_source(path: Path, source: Mapping[str, Any]) -> dict[str, Any]:
    """Verify a source against its manifest record, regardless of local name."""

    expected_bytes = int(source["bytes"])
    if path.stat().st_size != expected_bytes:
        raise ValueError(f"Source byte length changed for {path}")
    expected_sha = str(source["sha256"])
    actual_sha = sha256_file(path)
    if actual_sha != expected_sha:
        raise ValueError(f"Source SHA-256 changed for {path}")
    return {"bytes": expected_bytes, "sha256": actual_sha}


def _dataset_path(root: Path, manifest: Mapping[str, Any], spec: Mapping[str, Any]) -> Path:
    source_directory = Path(str(manifest.get("source_directory", "data/raw/liver-atlas")))
    return root / source_directory / str(spec["file"])


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def combined_exclusion_counts(results: Sequence[DatasetResult]) -> dict[str, int]:
    """Add exclusion ledgers across datasets without overwriting shared keys."""

    counts: Counter[str] = Counter()
    for result in results:
        counts.update(result.qc.get("excluded_by_reason", {}))
    return dict(sorted((str(key), int(value)) for key, value in counts.items()))


def run_analysis(
    *,
    root: Path = ROOT,
    plan_path: Path | None = None,
    dataset_paths: Mapping[str, Path] | None = None,
    output_directory: Path | None = None,
    audit_path: Path | None = None,
    chunk_size: int = 4096,
) -> dict[str, Any]:
    """Run both frozen atlas datasets and write derived tables."""

    plan_path = plan_path or (root / "config/liver-atlas-plan.json")
    plan, manifest = load_frozen_plan(plan_path)
    frozen_thresholds = validate_frozen_thresholds(plan)
    output_directory = output_directory or (root / "data/derived")
    audit_path = audit_path or (root / "reports/liver-atlas-expression-audit.json")
    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - dependency is pinned for the real run
        raise RuntimeError("anndata is required for the atlas analysis") from exc

    results: list[DatasetResult] = []
    source_records = {str(item["file"]): item for item in manifest.get("sources", [])}
    for spec in plan["datasets"]:
        dataset_key = str(spec["key"])
        path = Path(dataset_paths[dataset_key]) if dataset_paths and dataset_key in dataset_paths else _dataset_path(root, manifest, spec)
        if not path.exists():
            raise FileNotFoundError(path)
        manifest_source = source_records.get(str(spec["file"]))
        if manifest_source is None:
            raise ValueError(f"Frozen source manifest has no record for {spec['file']}")
        verify_frozen_source(path, manifest_source)
        adata = ad.read_h5ad(path, backed="r")
        try:
            results.append(aggregate_adata(adata, dataset_key,
                                           dataset_id=str(spec.get("dataset_id", dataset_key)),
                                           dataset_version_id=str(spec.get("version_id", "")),
                                           source_file=str(path.relative_to(root)) if path.is_relative_to(root) else str(path),
                                           expected_observations=int(spec.get("expected_observations", 0)) or None,
                                           chunk_size=chunk_size))
        finally:
            file_manager = getattr(adata, "file", None)
            if file_manager is not None and hasattr(file_manager, "close"):
                file_manager.close()

    pseudobulk = pd.concat([result.pseudobulk for result in results], ignore_index=True)
    accounting = pd.concat([result.cell_accounting for result in results], ignore_index=True)
    gene_audit = pd.concat([result.gene_audit for result in results], ignore_index=True)
    summaries = pd.concat([summarize_expression(result.pseudobulk, result.cell_accounting,
                                                 gene_audit=result.gene_audit) for result in results], ignore_index=True)
    contrasts = pd.concat([build_contrasts(result.pseudobulk, result.cell_accounting) for result in results], ignore_index=True)
    expected_cells_by_dataset = {result.dataset_key: int(result.qc["total_cells"]) for result in results}
    observed_cells_by_dataset = {str(key): int(value) for key, value in accounting.groupby("dataset_key")["cell_count"].sum().items()}
    if observed_cells_by_dataset != expected_cells_by_dataset:
        raise ValueError(f"Cell-accounting ledger does not reconcile to AnnData observations: {observed_cells_by_dataset} vs {expected_cells_by_dataset}")
    public_accounting = accounting.drop(columns=["donor_ids_observed", "library_ids_observed"], errors="ignore")
    outputs = {
        "work/liver-atlas-expression-donor-pseudobulk.csv": pseudobulk,
        "data/derived/liver-atlas-expression-cell-accounting.csv": public_accounting,
        "data/derived/liver-atlas-expression-gene-audit.csv": gene_audit,
        "data/derived/liver-atlas-expression-summary.csv": summaries,
        "data/derived/liver-atlas-expression-contrasts.csv": contrasts,
    }
    output_paths: dict[str, str] = {}
    for relative, frame in outputs.items():
        path = root / relative
        if output_directory != root / "data/derived":
            path = output_directory / Path(relative).name
        _write_csv(frame, path)
        output_paths[str(path.relative_to(root)) if path.is_relative_to(root) else str(path)] = sha256_file(path)

    exclusion_counts = combined_exclusion_counts(results)
    donors_by_dataset = {
        result.dataset_key: set(result.pseudobulk["donor_id"].astype(str))
        for result in results
    }
    donor_overlap = {
        "dataset_donor_counts": {key: len(value) for key, value in sorted(donors_by_dataset.items())},
        "shared_donor_counts": {
            f"{left}_vs_{right}": len(donors_by_dataset[left] & donors_by_dataset[right])
            for left in sorted(donors_by_dataset)
            for right in sorted(donors_by_dataset)
            if left < right
        },
    }
    dataset_specs = {str(spec["key"]): spec for spec in plan["datasets"]}
    frozen_source_by_dataset = {
        key: {
            "source_file": str(spec.get("file", "")),
            "source_bytes": int(source_records[str(spec["file"])]["bytes"]),
            "source_sha256": str(source_records[str(spec["file"])]["sha256"]),
        }
        for key, spec in dataset_specs.items()
    }
    audit = {
        "analysis_id": plan.get("analysis_id", "psc-liver-atlas-five-gene-descriptive-v1"),
        "plan_sha256": sha256_file(plan_path),
        "source_manifest_sha256": sha256_file((plan_path.parent.parent / plan["source_manifest"]).resolve()),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "requirements_sha256": sha256_file(root / "requirements-liver-atlas.txt") if (root / "requirements-liver-atlas.txt").exists() else None,
        "genes": list(GENES),
        "count_source": "raw.X",
        "normalized_matrix_policy": "AnnData.X was checked as a separate matrix and never used as UMI counts",
        "thresholds": {"primary_min_cells": PRIMARY_MIN_CELLS, "sensitivity_min_cells": SENSITIVITY_MIN_CELLS,
                       "minimum_donors_for_contrast": MIN_DONORS_FOR_CONTRAST},
        "datasets": [{"dataset_key": result.dataset_key, **frozen_source_by_dataset[result.dataset_key], **result.qc} for result in results],
        "cell_accounting": {
            "total_cells": int(sum(result.qc["total_cells"] for result in results)),
            "retained_cells": int(sum(result.qc["retained_cells"] for result in results)),
            "excluded_cells": int(sum(result.qc["excluded_cells"] for result in results)),
            "ledger_cells": int(accounting["cell_count"].sum()),
            "excluded_by_reason": exclusion_counts,
        },
        "donor_overlap": donor_overlap,
        "contrast_status_counts": {str(key): int(value) for key, value in contrasts["status"].value_counts().items()} if len(contrasts) else {},
        "outputs_sha256": output_paths,
        "limitations": "Donor-level descriptive expression only; no cell-level p-values, causal interpretation or treatment inference. SC cross-chemistry PSC-Control contrasts are blocked. SN contrasts retain residual processing-confounding flags.",
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True, allow_nan=False, default=_json_default) + "\n", encoding="utf-8")
    return audit


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=ROOT / "config/liver-atlas-plan.json")
    parser.add_argument("--dataset", action="append", metavar="KEY=PATH",
                        help="Override one frozen dataset path; may be repeated")
    parser.add_argument("--output-directory", type=Path, default=None)
    parser.add_argument("--audit-path", type=Path, default=None)
    parser.add_argument("--chunk-size", type=int, default=4096)
    args = parser.parse_args(argv)
    overrides: dict[str, Path] = {}
    for item in args.dataset or []:
        if "=" not in item:
            parser.error("--dataset must be KEY=PATH")
        key, value = item.split("=", 1)
        if not key or not value:
            parser.error("--dataset must be KEY=PATH")
        overrides[key] = Path(value)
    audit = run_analysis(plan_path=args.plan, dataset_paths=overrides or None,
                         output_directory=args.output_directory, audit_path=args.audit_path,
                         chunk_size=args.chunk_size)
    print(json.dumps(audit, indent=2, sort_keys=True, default=_json_default))


if __name__ == "__main__":
    main()
