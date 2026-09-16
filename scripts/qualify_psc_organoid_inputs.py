#!/usr/bin/env python3
"""Acquire and qualify GSE239283 released inputs without testing any gene effect.

Native environment: .venv/bin/python scripts/qualify_psc_organoid_inputs.py --acquire
All donor joins, cell QC, and count intermediates remain in ignored work/.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/psc-organoid-il17"
WORK = ROOT / "work/psc-organoid-inputs"
OUTPUT = ROOT / "data/derived/psc-organoid-input-qualification.json"
REPORT = ROOT / "reports/psc-organoid-input-qualification.md"
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239283/suppl/"
PREFIX = "GSE239283_extrahepatic_organoids_"
ACQUISITION_LIMIT_BYTES = 200_000_000
EXPECTED_PUBLIC_SOURCE_SHA256 = {'A900-methods.docx': 'b839464e46e3e74ce8bd03168c11b6736271e0df0ff6bb3e0bcde0e6ae39d69a', 'GSE239283-family.soft.gz': 'c403f41250c2320bc70f258037081004ece5c3770064252fc3b1390b9c9f1d40', 'GSE239283-series.txt': '9d4272909db6033ceac85315ebfaa68bcd4db96b4e5047bf521eb474a7f8cfc2', 'GSE239283-suppl-listing.html': '6459f4f0c2ac100dbf17b0ec46758806767841e882f3d4808f502237d6a51fa2', 'GSE239283_RAW.tar': 'eaac0234ea17cc610464385b81453a7322c610f410a8456482514e0b3288551b', 'GSE239283_extrahepatic_organoids_barcodes.tsv.gz': '5691998ebcd984ff169bab1340dfeadb9de2de66d0a1d6ec772f8c220afc5ce9', 'GSE239283_extrahepatic_organoids_genes.tsv.gz': 'b3f1b4c141d90c7d6445484fc312db7c2c4efc39f5f712ca3cc40ca10ae53f96', 'GSE239283_extrahepatic_organoids_matrix.mtx.gz': '48ff435b1e6ac6374932fe8a38953b925b6f6c1fbd9875296b772ef983b38b9c', 'GSE239283_extrahepatic_organoids_metadata.tsv.gz': '319739ce31587080c69f3892304281b701daf1bae7089b6974b1182312bf0ace', 'PMC11150034.xml': '6903bd9b6f603375d934c1208aa11a6d834f675da9d9a05a8270fa1219f228d6', 'filelist.txt': '0ba40aef7dc8c3ffba10a493da2febaeb73281083f447feb0ec44ea5e6166efb'}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def acquire(raw=RAW):
    """Copy verified discovery bodies, then GET only the three missing full inputs."""
    import requests

    raw = Path(raw)
    raw.mkdir(parents=True, exist_ok=True)
    manifest_path = raw / "source-manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {"schema_version": 1, "acquisition_limit_bytes": ACQUISITION_LIMIT_BYTES,
                    "raw_reads_acquired": False, "sources": []}
    sources = manifest["sources"]
    discovery = ROOT / "data/raw/public-data-discovery-options"
    receipt_path = discovery / "source-receipts.json"
    receipts = json.loads(receipt_path.read_text())
    reuse = {"GSE239283-series.txt": "GSE239283-series.txt",
             "GSE239283-family.soft.gz": "GSE239283-family.soft.gz",
             "GSE239283-organoid-metadata.tsv.gz": PREFIX + "metadata.tsv.gz",
             "PMC11150034.xml": "PMC11150034.xml",
             "GSE239283-suppl-listing.html": "GSE239283-suppl-listing.html"}
    existing = {row["file"]: row for row in sources}
    for old_name, new_name in reuse.items():
        if new_name in existing:
            continue
        matching = [row for row in receipts if row.get("file") == old_name and row.get("status") == 200]
        if len(matching) != 1:
            raise ValueError(f"Need one original receipt: {old_name}")
        receipt = matching[0]
        original = discovery / old_name
        if original.stat().st_size != receipt["bytes"] or sha256(original) != receipt["sha256"]:
            raise ValueError(f"Original cached body does not match receipt: {old_name}")
        shutil.copyfile(original, raw / new_name)
        sources.append({"file": new_name, "url": receipt["url"],
                        "retrieved_at_utc": receipt["retrieved_at_utc"],
                        "copied_at_utc": utc_now(), "bytes": receipt["bytes"],
                        "sha256": receipt["sha256"], "acquisition": "verified_cache_copy",
                        "copied_from": str(original.relative_to(ROOT)),
                        "original_receipts_file": str(receipt_path.relative_to(ROOT)),
                        "original_receipts_sha256": sha256(receipt_path),
                        "original_receipt": receipt})
        dump_json(manifest_path, manifest)
    for suffix in ["genes.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz"]:
        name = PREFIX + suffix
        if any(row["file"] == name for row in sources):
            continue
        url = BASE_URL + name
        start = utc_now()
        tmp = raw / (name + ".partial")
        maximum = ACQUISITION_LIMIT_BYTES - sum(row["bytes"] for row in sources)
        with requests.get(url, stream=True, timeout=(30, 120),
                          headers={"Accept-Encoding": "identity"}) as response:
            response.raise_for_status()
            if response.status_code != 200 or "Content-Range" in response.headers:
                raise ValueError("Complete non-range response required")
            declared = response.headers.get("Content-Length")
            if declared is not None and int(declared) > maximum:
                raise ValueError("Acquisition would exceed authorized byte limit")
            count = 0
            digest = hashlib.sha256()
            with tmp.open("wb") as stream:
                for block in response.iter_content(1024 * 1024):
                    count += len(block)
                    if count > maximum:
                        raise ValueError("Acquisition exceeded authorized byte limit")
                    stream.write(block)
                    digest.update(block)
            if declared is not None and int(declared) != count:
                raise ValueError("Response body length differs from Content-Length")
            if suffix == "matrix.mtx.gz" and count != 142_381_971:
                raise ValueError("Matrix byte size differs from pinned discovery header")
            tmp.replace(raw / name)
            sources.append({"file": name, "url": url, "response_url": response.url,
                            "started_at_utc": start, "retrieved_at_utc": utc_now(),
                            "status": response.status_code, "bytes": count,
                            "sha256": digest.hexdigest(), "acquisition": "complete_GET",
                            "headers": {key: value for key, value in response.headers.items()
                                        if key.lower() != "set-cookie"}})
        dump_json(manifest_path, manifest)
        print(json.dumps({"acquired": name, "bytes": count, "sha256": digest.hexdigest()}), flush=True)
    for receipt in sources:
        path = raw / receipt["file"]
        if path.stat().st_size != receipt["bytes"] or sha256(path) != receipt["sha256"]:
            raise ValueError(f"Pinned source mismatch: {path}")
    manifest["all_source_bytes"] = sum(row["bytes"] for row in sources)
    manifest["new_download_bytes"] = sum(row["bytes"] for row in sources if row["acquisition"] == "complete_GET")
    manifest["complete"] = True
    dump_json(manifest_path, manifest)
    return manifest




def fetch_source_file(name, url, manifest, raw=RAW, expected_bytes=None):
    import requests
    raw = Path(raw)
    old = next((row for row in manifest["sources"] if row["file"] == name), None)
    if old is not None:
        if (raw / name).stat().st_size != old["bytes"] or sha256(raw / name) != old["sha256"]:
            raise ValueError("Existing source differs from pinned receipt")
        return old
    remaining = manifest["acquisition_limit_bytes"] - sum(row["bytes"] for row in manifest["sources"])
    start = utc_now()
    partial = raw / (name + ".partial")
    with requests.get(url, stream=True, timeout=(30, 120), headers={"Accept-Encoding": "identity"}) as response:
        response.raise_for_status()
        if response.status_code != 200 or "Content-Range" in response.headers:
            raise ValueError("Complete non-range response required")
        length = response.headers.get("Content-Length")
        if length is not None and int(length) > remaining:
            raise ValueError("Request would exceed authorized cumulative processed-input limit")
        count = 0
        digest = hashlib.sha256()
        with partial.open("wb") as stream:
            for block in response.iter_content(1024 * 1024):
                count += len(block)
                if count > remaining:
                    raise ValueError("Request exceeded authorized cumulative processed-input limit")
                stream.write(block)
                digest.update(block)
        if length is not None and int(length) != count:
            raise ValueError("Source body length mismatch")
        if expected_bytes is not None and count != expected_bytes:
            raise ValueError("Source differs from expected byte size")
        partial.replace(raw / name)
        receipt = {"file": name, "url": url, "response_url": response.url,
                   "started_at_utc": start, "retrieved_at_utc": utc_now(), "status": response.status_code,
                   "bytes": count, "sha256": digest.hexdigest(), "acquisition": "complete_GET",
                   "headers": {key: value for key, value in response.headers.items() if key.lower() != "set-cookie"}}
    if receipt["sha256"] != EXPECTED_PUBLIC_SOURCE_SHA256[name]:
        raise ValueError("Downloaded source differs from qualified public source pin")
    manifest["sources"].append(receipt)
    manifest["all_source_bytes"] = sum(row["bytes"] for row in manifest["sources"])
    manifest["new_download_bytes"] = sum(row["bytes"] for row in manifest["sources"] if row["acquisition"] == "complete_GET")
    dump_json(raw / "source-manifest.json", manifest)
    print(json.dumps({"acquired": name, "bytes": count, "sha256": digest.hexdigest()}), flush=True)
    return receipt


def acquire_original(raw=RAW, inventory_only=False):
    """Acquire original processed 10x exports only, never FASTQ/raw reads."""
    raw = Path(raw)
    manifest = json.loads((raw / "source-manifest.json").read_text())
    manifest["acquisition_limit_bytes"] = 1_000_000_000
    manifest["authorization"] = {
        "initial_limit_bytes": ACQUISITION_LIMIT_BYTES,
        "expanded_limit_bytes": 1_000_000_000,
        "reason": "Merged matrix has SCT rather than raw RNA assay margins; root authorized original processed exports after inventory/format qualification.",
        "message_id": "agentmsg_4d716f46-058b-4e85-a240-59c4dc1c5868",
        "scope": "16 original processed RNA count matrix/feature/barcode exports or series TAR; no FASTQ/raw-read processing"}
    dump_json(raw / "source-manifest.json", manifest)
    fetch_source_file("A900-methods.docx", "https://links.lww.com/HC9/A900", manifest, raw,
                      expected_bytes=3_961_515)
    fetch_source_file("filelist.txt", BASE_URL + "filelist.txt", manifest, raw)
    filelist = (raw / "filelist.txt").read_text()
    print(filelist, flush=True)
    if inventory_only:
        return manifest
    names = []
    for line in filelist.splitlines():
        parts = line.split("\t")
        if len(parts) < 2 or not parts[0].startswith("File"):
            continue
        names.append(parts[1])
    expected = set()
    for record in parse_family(raw / "GSE239283-family.soft.gz"):
        for key, values in record["fields"].items():
            if key.startswith("supplementary_file_"):
                expected.update(value.rsplit("/", 1)[-1] for value in values)
    if len(expected) != 48 or len(names) != 48 or len(set(names)) != 48 or set(names) != expected:
        raise ValueError("Archive inventory does not exactly match 48 original GSM processed exports")
    if any(not re.fullmatch(r"GSM[0-9]+_(?:nonPSC|PSC)_[0-9]+_(?:CNT|IL17)\.(?:barcodes\.tsv|features\.tsv|matrix\.mtx)\.gz", name) for name in names):
        raise ValueError("Archive inventory includes a non-approved processed member")
    manifest["original_processed_inventory"] = {"file_count": len(names), "exact_GSM_supplement_links": True,
                                                 "only_processed_matrix_features_barcodes": True}
    dump_json(raw / "source-manifest.json", manifest)
    archive_rows = [line.split("\t") for line in filelist.splitlines() if line.startswith("Archive\t")]
    if len(archive_rows) != 1 or archive_rows[0][1] != "GSE239283_RAW.tar":
        raise ValueError("Expected exactly the source series archive")
    fetch_source_file("GSE239283_RAW.tar", BASE_URL + "GSE239283_RAW.tar", manifest, raw,
                      expected_bytes=int(archive_rows[0][3]))
    return manifest


def read_gzip_rows(path, delimiter="\t"):
    """Read the complete gzip member, thereby validating CRC and EOF."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.reader(stream, delimiter=delimiter))


def parse_family(path):
    """Preserve multivalued original GSM fields, without contact fields in output."""
    records = []
    current = None
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                current = {"accession": line.split(" = ", 1)[1], "fields": {}}
                records.append(current)
            elif line.startswith("^"):
                current = None
            elif current is not None and line.startswith("!Sample_") and " = " in line:
                key, value = line[len("!Sample_"):].split(" = ", 1)
                current["fields"].setdefault(key, []).append(value)
    return records


def unique_index(values, label):
    if any(not value for value in values):
        raise ValueError(f"Empty {label}")
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate {label}")
    return {value: index for index, value in enumerate(values)}


def qualify_metadata(metadata_rows, barcodes, family):
    """Exact barcode and original identity joins; never infer missing identity."""
    header, *rows = metadata_rows
    if len(header) != len(set(header)):
        raise ValueError("Duplicate metadata column")
    if not header or header[0] != "":
        raise ValueError("Expected original blank barcode column header")
    if any(len(row) != len(header) for row in rows):
        raise ValueError("Malformed metadata row")
    required = {"orig.ident", "Sample", "Treatment", "condition", "conditionTreatment",
                "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT"}
    if not required.issubset(header):
        raise ValueError("Missing original metadata field")
    metadata = [dict(zip(header, row, strict=True)) for row in rows]
    meta_index = unique_index([row[""] for row in metadata], "metadata barcode")
    unique_index(barcodes, "released barcode")
    if set(meta_index) != set(barcodes):
        raise ValueError("Metadata and released barcode sets differ")
    aligned = [metadata[meta_index[barcode]] for barcode in barcodes]
    unique_index([record["accession"] for record in family], "GSM accession")
    titles = []
    gsm = {}
    for record in family:
        fields = record["fields"]
        if "geo_accession" in fields and fields["geo_accession"] != [record["accession"]]:
            raise ValueError("GSM accession field differs from record identifier")
        if len(fields.get("title", [])) != 1:
            raise ValueError("Missing or multiple GSM title")
        title = fields["title"][0]
        titles.append(title)
        gsm[title] = record
    unique_index(titles, "GSM library title")
    if set(titles) != {row["orig.ident"] for row in aligned}:
        raise ValueError("GSM titles and metadata library IDs differ")
    libraries = []
    donor_arms = defaultdict(set)
    for title in titles:
        record = gsm[title]
        fields = record["fields"]
        cells = [row for row in aligned if row["orig.ident"] == title]
        identity = {key: {row[key] for row in cells}
                    for key in ["Sample", "Treatment", "condition", "conditionTreatment"]}
        if any(len(value) != 1 for value in identity.values()):
            raise ValueError("Inconsistent original library identity")
        original = {key: next(iter(value)) for key, value in identity.items()}
        donor = original["Sample"]
        treatment = original["Treatment"]
        if treatment not in {"CNT", "IL17"} or title != donor + "_" + treatment:
            raise ValueError("Library title disagrees with original donor/treatment")
        if original["condition"] not in {"PSC", "Non.PSC"}:
            raise ValueError("Unknown source condition")
        disease = {"PSC": "PSC", "Non.PSC": "nonPSC"}[original["condition"]]
        if not donor.startswith(disease + "_"):
            raise ValueError("Original condition disagrees with donor label")
        if original["conditionTreatment"] != original["condition"] + "_" + treatment:
            raise ValueError("Original conditionTreatment mismatch")
        characteristics = fields.get("characteristics_ch1", [])
        characteristic_values = defaultdict(list)
        for value in characteristics:
            if ": " in value:
                key, original_value = value.split(": ", 1)
                characteristic_values[key].append(original_value)
        expected_diagnosis = ("primary sclerosing cholangitis (PSC)" if disease == "PSC" else "non-PSC")
        expected_treatment = "vehicle" if treatment == "CNT" else "IL17A"
        for key, expected in {"patient diagnosis": expected_diagnosis,
                              "treatment": expected_treatment, "passage": "passage 3"}.items():
            if characteristic_values[key] != [expected]:
                raise ValueError("GSM characteristics disagree with original identity")
        donor_arms[donor].add(treatment)
        libraries.append({"library_id": title, "donor_id": donor, "disease": disease,
                          "treatment": treatment, "n_cells": len(cells),
                          "source_accession": record["accession"], "source_title": title,
                          "source_condition": original["condition"],
                          "source_Sample": donor, "source_Treatment": treatment,
                          "source_conditionTreatment": original["conditionTreatment"]})
    if any(arms != {"CNT", "IL17"} for arms in donor_arms.values()):
        raise ValueError("Incomplete donor-treatment pair")
    if len(libraries) != 2 * len(donor_arms):
        raise ValueError("More than one source library per donor-treatment")
    donor_diseases = defaultdict(set)
    for row in libraries:
        donor_diseases[row["donor_id"]].add(row["disease"])
    if any(len(values) != 1 for values in donor_diseases.values()):
        raise ValueError("Donor spans disease groups")
    return header, aligned, libraries, {
        "metadata_rows": len(aligned), "barcodes": len(barcodes),
        "exact_barcode_set_match": True,
        "metadata_order_matches_barcode_axis": [row[""] for row in metadata] == barcodes,
        "exact_GSM_library_title_set_match": True,
        "original_identity_columns": ["", "orig.ident", "Sample", "Treatment", "condition", "conditionTreatment"],
        "source_metadata_columns": header,
        "donors": len(donor_arms), "libraries": len(libraries),
        "complete_donor_treatment_pairs": len(donor_arms),
        "donors_per_disease": dict(sorted(Counter(next(iter(values)) for values in donor_diseases.values()).items())),
        "libraries_per_donor": 2,
        "identity_mapping": {"condition": {"PSC": "PSC", "Non.PSC": "nonPSC"},
                             "Treatment": {"CNT": "vehicle", "IL17": "IL-17A"},
                             "donor_id": "unchanged original Sample", "library_id": "unchanged orig.ident = exact GSM title"},
        "source_value_corrections": [],
    }


INTEGER_LINE = re.compile(rb"(?m)^[+-]?[0-9]{1,18}[ \t]+[+-]?[0-9]{1,18}[ \t]+[+-]?[0-9]{1,18}[ \t]*\r?$")


def scan_matrix(path, cell_library, n_libraries, expected_shape=None, chunk_lines=131072):
    """Check all sparse entries and aggregate only raw stored integers.

    A packed one-bit coordinate ledger detects duplicates, including across
    chunk boundaries and arbitrary coordinate order. It is not a dense count
    matrix. The only dense count matrix is features x source libraries.
    No gene is selected, tested, normalized, contrasted, or ranked.
    """
    import itertools
    import numpy as np

    cell_library = np.asarray(cell_library, dtype=np.int64)
    if chunk_lines <= 0:
        raise ValueError("chunk_lines must be positive")
    if n_libraries <= 0 or np.any(cell_library < 0) or np.any(cell_library >= n_libraries):
        raise ValueError("Invalid cell-library index")
    with gzip.open(path, "rb") as stream:
        banner = stream.readline().decode("ascii").rstrip("\r\n")
        if banner != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError("Expected integer general coordinate MatrixMarket")
        line = stream.readline()
        while line.startswith(b"%"):
            line = stream.readline()
        dims = line.split()
        if len(dims) != 3 or any(not item.isdigit() for item in dims):
            raise ValueError("Invalid MatrixMarket dimensions")
        n_features, n_cells, n_entries = map(int, dims)
        if n_features <= 0 or n_cells <= 0 or n_entries < 0:
            raise ValueError("Invalid MatrixMarket shape")
        if expected_shape is not None and (n_features, n_cells) != tuple(expected_shape):
            raise ValueError("Matrix dimensions disagree with feature/barcode axes")
        if len(cell_library) != n_cells:
            raise ValueError("Cell-library mapping dimension mismatch")
        n_coordinates = n_features * n_cells
        seen = np.zeros((n_coordinates + 7) // 8, dtype=np.uint8)
        library_sums = np.zeros((n_features, n_libraries), dtype=np.int64)
        cell_sums = np.zeros(n_cells, dtype=np.int64)
        cell_features = np.zeros(n_cells, dtype=np.int64)
        feature_cells = np.zeros(n_features, dtype=np.int64)
        stored = explicit_zeros = total = 0
        minimum = maximum = None
        previous_key = -1
        strictly_column_major = True
        for lines in iter(lambda: list(itertools.islice(stream, chunk_lines)), []):
            payload = b"".join(lines)
            residue = INTEGER_LINE.sub(b"", payload)
            if residue.strip(b"\r\n"):
                raise ValueError("Invalid entry: require finite integer tokens within 18 digits")
            values = np.fromstring(payload.decode("ascii"), sep=" ", dtype=np.int64)
            if values.size != 3 * len(lines):
                raise ValueError("Malformed sparse entry width")
            triples = values.reshape(-1, 3)
            rows, cols, counts = triples.T
            if np.any(rows < 1) or np.any(rows > n_features) or np.any(cols < 1) or np.any(cols > n_cells):
                raise ValueError("Sparse coordinate outside declared bounds")
            if np.any(counts < 0):
                raise ValueError("Negative sparse count")
            # This conservative bound makes every int64 aggregation overflow-safe.
            if counts.max(initial=0) > np.iinfo(np.int64).max // max(n_entries, 1):
                raise ValueError("Count exceeds overflow-safe aggregation bound")
            rows = rows - 1
            cols = cols - 1
            keys = cols * n_features + rows
            if keys[0] <= previous_key or np.any(keys[1:] <= keys[:-1]):
                strictly_column_major = False
            previous_key = int(keys[-1])
            unique_keys = np.unique(keys)
            if unique_keys.size != keys.size:
                raise ValueError("Duplicate sparse coordinate within chunk")
            byte_index = unique_keys // 8
            masks = np.left_shift(np.uint8(1), (unique_keys % 8).astype(np.uint8))
            if np.any(np.bitwise_and(seen[byte_index], masks)):
                raise ValueError("Duplicate sparse coordinate across chunks")
            np.bitwise_or.at(seen, byte_index, masks)
            stored += len(triples)
            if stored > n_entries:
                raise ValueError("More sparse entries than declared")
            explicit_zeros += int(np.count_nonzero(counts == 0))
            chunk_min, chunk_max = int(counts.min()), int(counts.max())
            minimum = chunk_min if minimum is None else min(minimum, chunk_min)
            maximum = chunk_max if maximum is None else max(maximum, chunk_max)
            total += int(counts.sum())
            np.add.at(cell_sums, cols, counts)
            positive = counts > 0
            np.add.at(cell_features, cols[positive], 1)
            np.add.at(feature_cells, rows[positive], 1)
            np.add.at(library_sums, (rows, cell_library[cols]), counts)
        # Reading through the final iteration validates gzip footer/CRC/EOF.
        if stored != n_entries:
            raise ValueError("Sparse entry count differs from declared count")
    if int(library_sums.sum()) != total or int(cell_sums.sum()) != total:
        raise ValueError("Aggregated sums disagree with full sparse total")
    audit = {"banner": banner, "features": n_features, "cells": n_cells,
             "declared_entries": n_entries, "observed_entries": stored,
             "positive_entries": stored - explicit_zeros, "explicit_zero_entries": explicit_zeros,
             "duplicate_coordinates": 0, "out_of_bounds_coordinates": 0,
             "nonfinite_entries": 0, "negative_entries": 0, "noninteger_entries": 0,
             "minimum_stored_value": minimum, "maximum_stored_value": maximum,
             "full_released_matrix_sum": total, "gzip_crc_and_eof_verified": True,
             "strictly_column_major": strictly_column_major,
             "coordinate_index_base": 1, "coordinate_ledger_bytes": len(seen),
             "all_cells_and_features_retained": True,
             "zero_total_cells": int(np.count_nonzero(cell_sums == 0)),
             "zero_total_features": int(np.count_nonzero(feature_cells == 0)),
             "library_sum_conservation": True, "cell_sum_conservation": True}
    return audit, library_sums, cell_sums, cell_features, feature_cells


def numeric_summary(values):
    import numpy as np
    values = np.asarray(values)
    quantiles = np.quantile(values, [0, .25, .5, .75, 1])
    return dict(zip(["minimum", "q25", "median", "q75", "maximum"], map(float, quantiles), strict=True))


def exact_nonnegative_integers(rows, column):
    from decimal import Decimal, InvalidOperation
    import numpy as np
    integers = []
    for row in rows:
        try:
            value = Decimal(row[column])
        except InvalidOperation as error:
            raise ValueError(f"Nonnumeric original metadata value in {column}") from error
        if (not value.is_finite() or value < 0 or value != value.to_integral_value()
                or value > np.iinfo(np.int64).max):
            raise ValueError(f"Noninteger original metadata value in {column}")
        integers.append(int(value))
    return np.asarray(integers, dtype=np.int64)


def compare_metadata_counts(rows, cell_sums, cell_features):
    """Test measurement provenance using every cell, never a selected gene."""
    import numpy as np
    result = {}
    for source_column, observed in [("nCount_RNA", cell_sums), ("nFeature_RNA", cell_features),
                                    ("nCount_SCT", cell_sums), ("nFeature_SCT", cell_features)]:
        source = exact_nonnegative_integers(rows, source_column)
        difference = observed - source
        result[source_column] = {"matching_cells": int(np.count_nonzero(difference == 0)),
                                 "differing_cells": int(np.count_nonzero(difference != 0)),
                                 "source_sum": sum(map(int, source)),
                                 "matrix_sum": sum(map(int, observed)),
                                 "minimum_matrix_minus_metadata": int(difference.min()),
                                 "maximum_matrix_minus_metadata": int(difference.max()),
                                 "all_cells_equal": bool(np.all(difference == 0))}
    raw_matches = result["nCount_RNA"]["all_cells_equal"] and result["nFeature_RNA"]["all_cells_equal"]
    sct_matches = result["nCount_SCT"]["all_cells_equal"] and result["nFeature_SCT"]["all_cells_equal"]
    result["interpretation"] = {
        "raw_RNA_margins_match": raw_matches, "SCT_margins_match": sct_matches,
        "raw_count_inference_eligible": False,
        "margin_only_assessment": True,
        "status": ("SCT_assay_margins_match_not_raw_RNA" if sct_matches and not raw_matches else
                   "RNA_assay_margins_match_export_provenance_pending" if raw_matches and not sct_matches else
                   "assay_quantity_unresolved"),
        "caution": "Matching cell sums and nonzero-feature counts identify an assay-consistent quantity, not the exact historical export command. Integer values alone do not prove unnormalized RNA counts."
    }
    return result


def write_tsv(path, header, rows):
    """Stable gzip header and source-order serialization; patient joins stay ignored."""
    import io
    path = Path(path)
    if path.suffix == ".gz":
        with path.open("wb") as binary:
            with gzip.GzipFile(fileobj=binary, mode="wb", filename="", mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                    writer.writerow(header)
                    writer.writerows(rows)
    else:
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)


def qualify(raw=RAW, work=WORK, output=OUTPUT, chunk_lines=131072):
    import numpy as np

    raw, work, output = Path(raw), Path(work), Path(output)
    work.mkdir(parents=True, exist_ok=True)
    manifest_path = raw / "source-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if not manifest.get("complete"):
        raise ValueError("Complete acquisition manifest required")
    for receipt in manifest["sources"]:
        path = raw / receipt["file"]
        if (path.stat().st_size != receipt["bytes"] or sha256(path) != receipt["sha256"]
                or receipt["sha256"] != EXPECTED_PUBLIC_SOURCE_SHA256[receipt["file"]]):
            raise ValueError(f"Pinned source mismatch: {path}")
    genes = read_gzip_rows(raw / (PREFIX + "genes.tsv.gz"))
    if not genes or any(len(row) != 2 for row in genes):
        raise ValueError("Expected two original gene-list fields")
    features = [row[0] for row in genes]
    unique_index(features, "source feature ID")
    barcodes_rows = read_gzip_rows(raw / (PREFIX + "barcodes.tsv.gz"))
    if any(len(row) != 1 for row in barcodes_rows):
        raise ValueError("Expected one-field barcode axis")
    barcodes = [row[0] for row in barcodes_rows]
    metadata = read_gzip_rows(raw / (PREFIX + "metadata.tsv.gz"), delimiter=",")
    family = parse_family(raw / "GSE239283-family.soft.gz")
    metadata_header, aligned, libraries, join_audit = qualify_metadata(metadata, barcodes, family)
    if join_audit["donors_per_disease"] != {"PSC": 4, "nonPSC": 4} or join_audit["libraries"] != 16:
        raise ValueError("Source donor/library design differs from discovery")
    library_index = {row["library_id"]: index for index, row in enumerate(libraries)}
    cell_library = np.asarray([library_index[row["orig.ident"]] for row in aligned], dtype=np.int64)
    audit, counts, cell_sums, cell_features, feature_cells = scan_matrix(
        raw / (PREFIX + "matrix.mtx.gz"), cell_library, len(libraries),
        expected_shape=(len(features), len(barcodes)), chunk_lines=chunk_lines)
    if (audit["features"], audit["cells"], audit["observed_entries"]) != (21562, 46343, 46282529):
        raise ValueError("Complete matrix differs from pinned discovery dimensions")
    quantity = compare_metadata_counts(aligned, cell_sums, cell_features)
    for index, library in enumerate(libraries):
        chosen = cell_library == index
        library["total_released_counts"] = int(counts[:, index].sum())
        library["features_positive_in_library"] = int(np.count_nonzero(counts[:, index]))
        library["cell_released_count_summary"] = numeric_summary(cell_sums[chosen])
        library["cell_released_feature_summary"] = numeric_summary(cell_features[chosen])
        library_rows = [row for row, flag in zip(aligned, chosen, strict=True) if flag]
        library["source_nCount_RNA_total"] = sum(map(int, exact_nonnegative_integers(library_rows, "nCount_RNA")))
        library["source_nCount_SCT_total"] = sum(map(int, exact_nonnegative_integers(library_rows, "nCount_SCT")))
    conditions = []
    for disease in sorted({row["disease"] for row in libraries}):
        for treatment in ["CNT", "IL17"]:
            chosen = [row for row in libraries if row["disease"] == disease and row["treatment"] == treatment]
            conditions.append({"disease": disease, "treatment": treatment,
                               "donors": len(chosen), "libraries": len(chosen),
                               "cells": sum(row["n_cells"] for row in chosen),
                               "total_released_counts": sum(row["total_released_counts"] for row in chosen)})
    feature_audit = {"rows": len(genes), "source_fields_per_row": 2,
                     "first_column_unique": True,
                     "second_column_unique": len({row[1] for row in genes}) == len(genes),
                     "both_columns_identical_rows": sum(left == right for left, right in genes),
                     "first_column_Ensembl_gene_ID_rows": sum(bool(re.fullmatch(r"ENSG[0-9]+(?:\.[0-9]+)?", feature)) for feature in features),
                     "namespace": "source gene-symbol-style labels in both columns; not stable Ensembl IDs",
                     "symbol_alias_resolution_performed": False,
                     "gene_annotation_release": None,
                     "gene_annotation_release_status": "not declared in inspected GEO fields",
                     "all_released_features_preserved": True,
                     "new_feature_exclusions": 0}
    # These are global count-availability summaries, not gene-specific results or filters.
    positive_libraries = np.count_nonzero(counts, axis=1)
    eligibility = {"features_by_number_of_positive_libraries": {str(int(k)): int(v) for k, v in zip(*np.unique(positive_libraries, return_counts=True), strict=True)},
                   "features_with_positive_total": int(np.count_nonzero(feature_cells)),
                   "features_with_zero_total": int(np.count_nonzero(feature_cells == 0)),
                   "count_thresholds_applied": False, "cell_exclusions_applied": 0,
                   "gene_exclusions_applied": 0, "donor_exclusions_applied": 0}
    write_tsv(work / "features.tsv", ["source_feature_id", "source_gene_field_1", "source_gene_field_2"],
              ([row[0], *row] for row in genes))
    write_tsv(work / "released-matrix-library-sums.tsv.gz", ["source_feature_id", *library_index],
              ([feature, *map(int, values)] for feature, values in zip(features, counts, strict=True)))
    library_columns = ["library_id", "donor_id", "disease", "treatment", "n_cells", "total_released_counts",
                       "features_positive_in_library", "source_accession", "source_title", "source_condition",
                       "source_Sample", "source_Treatment", "source_conditionTreatment",
                       "source_nCount_RNA_total", "source_nCount_SCT_total"]
    write_tsv(work / "libraries.tsv", library_columns, ([row[key] for key in library_columns] for row in libraries))
    write_tsv(work / "cell-source-joins.tsv.gz",
              ["source_barcode", *metadata_header[1:], "matrix_total", "matrix_positive_features"],
              ([row[""], *[row[key] for key in metadata_header[1:]], int(total), int(n_features)]
               for row, total, n_features in zip(aligned, cell_sums, cell_features, strict=True)))
    dump_json(work / "original-GSM-fields.json", [
        {"accession": row["accession"], "fields": {key: value for key, value in row["fields"].items() if not key.startswith("contact_")}}
        for row in family])
    artifact_paths = ["features.tsv", "released-matrix-library-sums.tsv.gz", "libraries.tsv",
                      "cell-source-joins.tsv.gz", "original-GSM-fields.json"]
    artifacts = [{"file": name, "bytes": (work / name).stat().st_size,
                  "sha256": sha256(work / name), "git_policy": "ignored work only"} for name in artifact_paths]
    summary = {"schema_version": 1, "study": "GSE239283", "stage": "source qualification only",
               "source_manifest": {"path": str(manifest_path.relative_to(ROOT)), "sha256": sha256(manifest_path)},
               "sources": [{key: value for key, value in row.items() if key in {"file", "url", "retrieved_at_utc", "bytes", "sha256", "acquisition", "copied_from"}}
                           for row in manifest["sources"]],
               "acquisition": {key: manifest[key] for key in ["acquisition_limit_bytes", "all_source_bytes", "new_download_bytes", "raw_reads_acquired"]},
               "matrix": audit, "feature_axis": feature_audit, "joins": join_audit,
               "quantity_validation": quantity, "libraries": libraries, "conditions": conditions,
               "global_QC": {"cell_count_summary": numeric_summary(cell_sums),
                             "cell_feature_summary": numeric_summary(cell_features)},
               "count_availability": eligibility,
               "work_artifacts": artifacts,
               "software": {"python": platform.python_version(), "numpy": np.__version__,
                            "script": "scripts/qualify_psc_organoid_inputs.py", "script_sha256": sha256(__file__)},
               "analysis_performed": {"full_sparse_QC": True, "all_feature_library_sums": True,
                                      "gene_wise_group_or_treatment_tests": False, "target_expression_inspection": False,
                                      "program_scores": False, "effect_contrasts": False, "model_fitting": False},
               "ready_for_raw_count_inference": False,
               "root_freeze_required": True,
               "qualification_status": "structurally_valid_assay_quantity_must_be_resolved"}
    dump_json(output, summary)
    print(json.dumps({"output": str(output), "matrix": audit, "quantity_validation": quantity,
                      "joins": {key: join_audit[key] for key in ["donors", "libraries", "complete_donor_treatment_pairs", "donors_per_disease"]}}, sort_keys=True), flush=True)
    return summary


def safe_extract_processed_archive(archive_path, expected_members, destination):
    """No tar path is trusted; only an exact basename allowlist of regular files."""
    import tarfile
    from pathlib import PurePosixPath

    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "r:") as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        if len(names) != len(set(names)) or set(names) != set(expected_members):
            raise ValueError("Tar member set differs from exact processed-file allowlist")
        for member in members:
            path = PurePosixPath(member.name)
            if (not member.isfile() or member.issym() or member.islnk() or member.issparse()
                    or path.is_absolute() or len(path.parts) != 1 or ".." in path.parts
                    or "\\" in member.name or member.name in {"", ".", ".."}):
                raise ValueError("Unsafe tar member type or path")
            if member.size != expected_members[member.name]:
                raise ValueError("Tar member size differs from publisher inventory")
        ledger = []
        for member in members:
            output = destination / member.name
            if output.is_symlink() or not output.resolve().is_relative_to(destination.resolve()):
                raise ValueError("Unsafe extraction destination")
            source = archive.extractfile(member)
            if source is None:
                raise ValueError("Missing regular tar member body")
            digest = hashlib.sha256()
            count = 0
            import tempfile
            with source, tempfile.NamedTemporaryFile(mode="wb", dir=destination, prefix=".extract-", delete=False) as stream:
                temporary = Path(stream.name)
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    count += len(block)
                    stream.write(block)
                    digest.update(block)
            if count != member.size:
                temporary.unlink()
                raise ValueError("Incomplete extracted member")
            if output.exists() and output.stat().st_size == count and sha256(output) == digest.hexdigest():
                temporary.unlink()  # Do not rewrite an identical file during offline replay.
            else:
                temporary.replace(output)
            ledger.append({"member": member.name, "bytes": count, "sha256": digest.hexdigest()})
    return ledger


def join_original_barcodes(merged_rows, original_barcodes):
    """Verify an explicit suffix-removal map, preserving both original strings.

    This is a validated technical barcode join within an already exact source
    library; it does not guess donor identities from barcode suffix numbers.
    """
    original_index = unique_index(original_barcodes, "original library barcode")
    indices = []
    mappings = []
    suffixes = set()
    for row in merged_rows:
        merged = row[""]
        if "_" not in merged:
            raise ValueError("Merged barcode lacks a distinct technical merge suffix")
        original, suffix = merged.rsplit("_", 1)
        if not suffix.isdigit() or not suffix or int(suffix) <= 0:
            raise ValueError("Unqualified merged barcode suffix")
        if original not in original_index:
            raise ValueError("Released cell is absent from its exact source library")
        indices.append(original_index[original])
        suffixes.add(suffix)
        mappings.append({"source_merged_barcode": merged, "source_original_barcode": original,
                         "removed_merge_suffix": "_" + suffix, "source_library": row["orig.ident"]})
    if len(set(indices)) != len(indices):
        raise ValueError("Multiple released cells map to one original source barcode")
    if len(suffixes) != 1:
        raise ValueError("A source library spans multiple technical merge suffixes")
    return indices, mappings


def qualify_original(raw=RAW, work=WORK, output=OUTPUT, chunk_lines=131072):
    import numpy as np

    raw, work, output = Path(raw), Path(work), Path(output)
    summary = json.loads(output.read_text())
    manifest_path = raw / "source-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for receipt in manifest["sources"]:
        file = raw / receipt["file"]
        if (file.stat().st_size != receipt["bytes"] or sha256(file) != receipt["sha256"]
                or receipt["sha256"] != EXPECTED_PUBLIC_SOURCE_SHA256[receipt["file"]]):
            raise ValueError(f"Pinned source mismatch: {file}")
    archive_path = raw / "GSE239283_RAW.tar"
    inventory_rows = [line.split("\t") for line in (raw / "filelist.txt").read_text().splitlines()
                      if line.startswith("File\t")]
    expected_members = {row[1]: int(row[3]) for row in inventory_rows}
    if len(expected_members) != 48:
        raise ValueError("Expected 48 original processed exports")
    member_ledger = safe_extract_processed_archive(archive_path, expected_members, raw / "original")
    archive_ledger = {"archive_file": archive_path.name, "archive_bytes": archive_path.stat().st_size,
                      "archive_sha256": sha256(archive_path), "member_count": len(member_ledger),
                      "all_members_regular_basename_files": True,
                      "exact_publisher_inventory_match": True, "members": member_ledger}
    dump_json(raw / "original-archive-members.json", archive_ledger)
    merged_barcodes = [row[0] for row in read_gzip_rows(raw / (PREFIX + "barcodes.tsv.gz"))]
    metadata = read_gzip_rows(raw / (PREFIX + "metadata.tsv.gz"), delimiter=",")
    family = parse_family(raw / "GSE239283-family.soft.gz")
    _, aligned, libraries, _ = qualify_metadata(metadata, merged_barcodes, family)
    retained_features = None
    raw_sums = []
    audits = []
    joined_rows = []
    all_features_detected = None
    total_original_cells = total_original_entries = total_original_counts = 0
    selected_count_array = np.zeros(len(aligned), dtype=np.int64)
    selected_feature_array = np.zeros(len(aligned), dtype=np.int64)
    merged_index = {row[""]: index for index, row in enumerate(aligned)}
    for library in libraries:
        title, accession = library["library_id"], library["source_accession"]
        stem = accession + "_" + title
        feature_path = raw / "original" / (stem + ".features.tsv.gz")
        features = read_gzip_rows(feature_path)
        if any(len(row) != 3 for row in features):
            raise ValueError("Expected original 10x feature ID, label and type")
        unique_index([row[0] for row in features], "original feature ID")
        if any(row[2] != "Gene Expression" for row in features):
            raise ValueError("Original exports contain a non-RNA feature type")
        if retained_features is None:
            retained_features = features
            all_features_detected = np.zeros(len(features), dtype=np.int64)
        elif features != retained_features:
            raise ValueError("Original feature axes differ; explicit reconciliation needed")
        barcode_rows = read_gzip_rows(raw / "original" / (stem + ".barcodes.tsv.gz"))
        if any(len(row) != 1 for row in barcode_rows):
            raise ValueError("Expected one original barcode per row")
        original_barcodes = [row[0] for row in barcode_rows]
        unique_index(original_barcodes, "original barcode axis")
        library_rows = [row for row in aligned if row["orig.ident"] == title]
        selected_indices, mapping = join_original_barcodes(library_rows, original_barcodes)
        selection = np.ones(len(original_barcodes), dtype=np.int64)
        selection[selected_indices] = 0
        audit, sums, total_per_cell, features_per_cell, feature_cells = scan_matrix(
            raw / "original" / (stem + ".matrix.mtx.gz"), selection, 2,
            expected_shape=(len(features), len(original_barcodes)), chunk_lines=chunk_lines)
        total_original_cells += len(original_barcodes)
        total_original_entries += audit["observed_entries"]
        total_original_counts += audit["full_released_matrix_sum"]
        selected_totals = total_per_cell[selected_indices]
        selected_detected = features_per_cell[selected_indices]
        comparisons = compare_metadata_counts(library_rows, selected_totals, selected_detected)
        all_features_detected += sums[:, 0] > 0
        raw_sums.append(sums[:, 0])
        for mapping_row, metadata_row, total, detected in zip(mapping, library_rows, selected_totals, selected_detected, strict=True):
            merged_idx = merged_index[mapping_row["source_merged_barcode"]]
            selected_count_array[merged_idx] = total
            selected_feature_array[merged_idx] = detected
            joined_rows.append({**mapping_row, "source_accession": accession,
                                "source_Sample": metadata_row["Sample"], "source_Treatment": metadata_row["Treatment"],
                                "original_matrix_total": int(total), "original_matrix_positive_features": int(detected),
                                "source_nCount_RNA": metadata_row["nCount_RNA"], "source_nFeature_RNA": metadata_row["nFeature_RNA"]})
        library["total_released_counts"] = int(sums[:, 0].sum())
        library["features_positive_in_library"] = int(np.count_nonzero(sums[:, 0]))
        library["source_original_cells"] = len(original_barcodes)
        library["source_cells_not_in_released_subset"] = len(original_barcodes) - len(selected_indices)
        library["source_full_library_total"] = audit["full_released_matrix_sum"]
        library["source_unselected_cells_total"] = int(sums[:, 1].sum())
        library["source_nCount_RNA_total"] = sum(map(int, exact_nonnegative_integers(library_rows, "nCount_RNA")))
        library["cell_raw_count_summary"] = numeric_summary(selected_totals)
        library["cell_raw_feature_summary"] = numeric_summary(selected_detected)
        library["technical_merge_suffix"] = mapping[0]["removed_merge_suffix"]
        audits.append({"library_id": title, "source_accession": accession, "matrix": audit,
                       "retained_cells": len(selected_indices), "unretained_source_cells": len(original_barcodes) - len(selected_indices),
                       "retained_cell_quantity_validation": comparisons,
                       "exact_original_barcode_join": True,
                       "feature_axis_identical_to_first_library": True})
        print(json.dumps({"original_library_completed": title, "full_source_cells": len(original_barcodes),
                          "retained_cells": len(selected_indices), "RNA_margins_match": comparisons["interpretation"]["raw_RNA_margins_match"]}), flush=True)
    counts = np.column_stack(raw_sums)
    retained_count_total = sum(row["total_released_counts"] for row in libraries)
    selected_cell_count_total = sum(map(int, selected_count_array))
    feature_ids = [row[0] for row in retained_features]
    quantity = compare_metadata_counts(aligned, selected_count_array, selected_feature_array)
    if retained_count_total != selected_cell_count_total:
        raise ValueError("Raw RNA selected-cell sum conservation failure")
    raw_library_cols = ["library_id", "donor_id", "disease", "treatment", "n_cells", "total_released_counts",
                        "features_positive_in_library", "source_accession", "source_title", "source_condition",
                        "source_Sample", "source_Treatment", "source_conditionTreatment", "source_original_cells",
                        "source_cells_not_in_released_subset", "source_full_library_total", "source_unselected_cells_total",
                        "source_nCount_RNA_total", "technical_merge_suffix"]
    write_tsv(work / "raw-rna-pseudobulk.tsv.gz", ["source_feature_id", *[row["library_id"] for row in libraries]],
              ([feature, *map(int, values)] for feature, values in zip(feature_ids, counts, strict=True)))
    write_tsv(work / "raw-rna-features.tsv", ["source_feature_id", "source_gene_label", "source_feature_type"], retained_features)
    write_tsv(work / "raw-rna-libraries.tsv", raw_library_cols, ([row[key] for key in raw_library_cols] for row in libraries))
    join_columns = list(joined_rows[0])
    write_tsv(work / "raw-rna-cell-source-joins.tsv.gz", join_columns,
              ([row[key] for key in join_columns] for row in joined_rows))
    dump_json(work / "raw-rna-library-audits.json", audits)
    raw_files = ["raw-rna-pseudobulk.tsv.gz", "raw-rna-features.tsv", "raw-rna-libraries.tsv",
                 "raw-rna-cell-source-joins.tsv.gz", "raw-rna-library-audits.json"]
    raw_artifacts = [{"file": file, "bytes": (work / file).stat().st_size,
                      "sha256": sha256(work / file), "git_policy": "ignored work only"} for file in raw_files]
    raw_conditions = []
    for disease in ["PSC", "nonPSC"]:
        for treatment in ["CNT", "IL17"]:
            chosen = [row for row in libraries if row["disease"] == disease and row["treatment"] == treatment]
            raw_conditions.append({"disease": disease, "treatment": treatment, "donors": len(chosen), "libraries": len(chosen),
                                   "cells": sum(row["n_cells"] for row in chosen),
                                   "total_released_counts": sum(row["total_released_counts"] for row in chosen)})
    ready = bool(quantity["interpretation"]["raw_RNA_margins_match"] and not quantity["interpretation"]["SCT_margins_match"])
    original = {"acquisition": "Complete original per-GSM processed 10x exports from GEO series TAR; no raw reads",
                "archive_member_ledger": {"file": str((raw / "original-archive-members.json").relative_to(ROOT)),
                                          "sha256": sha256(raw / "original-archive-members.json")},
                "libraries_checked": len(libraries), "compressed_members_checked": len(member_ledger),
                "all_member_gzip_crc_and_eof_verified": True,
                "all_feature_axes_exactly_identical": True,
                "full_feature_count": len(feature_ids),
                "feature_namespace": "Original Ensembl gene IDs; separate unmodified gene-symbol labels and Gene Expression type",
                "Ensembl_gene_ID_rows": sum(bool(re.fullmatch(r"ENSG[0-9]+(?:\.[0-9]+)?", feature)) for feature in feature_ids),
                "unique_feature_IDs": len(set(feature_ids)),
                "unique_source_gene_labels": len({row[1] for row in retained_features}),
                "annotation_release": None, "genome_build_source_declaration": "hg38",
                "full_original_cells": total_original_cells, "full_original_stored_entries": total_original_entries,
                "full_original_count_total": total_original_counts,
                "retained_cells": len(aligned), "original_cells_not_in_released_metadata": total_original_cells - len(aligned),
                "retained_feature_count": len(feature_ids), "retained_count_total": retained_count_total,
                "retained_positive_entries": sum(map(int, selected_feature_array)),
                "zero_total_features_in_retained_cells": int(np.count_nonzero(all_features_detected == 0)),
                "positive_total_features_in_retained_cells": int(np.count_nonzero(all_features_detected)),
                "features_by_number_of_positive_libraries": {str(int(k)): int(v) for k, v in zip(*np.unique(all_features_detected, return_counts=True), strict=True)},
                "selected_cell_global_count_summary": numeric_summary(selected_count_array),
                "selected_cell_global_feature_summary": numeric_summary(selected_feature_array),
                "quantity_validation": quantity, "libraries": libraries, "conditions": raw_conditions,
                "full_sparse_QC": {key: sum(item["matrix"][key] for item in audits) for key in
                                   ["duplicate_coordinates", "explicit_zero_entries", "out_of_bounds_coordinates", "nonfinite_entries", "negative_entries", "noninteger_entries", "zero_total_cells"]},
                "barcode_join": {"status": "exact within source library", "rule": "Remove one verified library-specific terminal _integer merge suffix; retain both original strings",
                                 "all_released_cells_found_once": True, "source_library_not_inferred_from_barcode": True},
                "selection": "All and only released merged-metadata cells; complete original feature axis; no newly chosen cell/gene/donor exclusions",
                "raw_RNA_count_input_qualified": ready,
                "work_artifacts": raw_artifacts}
    summary["original_RNA"] = original
    summary["source_manifest"] = {"path": str(manifest_path.relative_to(ROOT)), "sha256": sha256(manifest_path)}
    summary["sources"] = [{key: value for key, value in row.items() if key in {"file", "url", "request_url", "retrieved_at_utc", "bytes", "sha256", "acquisition", "copied_from"}}
                          for row in manifest["sources"]]
    summary["acquisition"] = {key: manifest[key] for key in ["acquisition_limit_bytes", "all_source_bytes", "new_download_bytes", "raw_reads_acquired", "authorization"]}
    summary["ready_for_raw_count_inference"] = False
    summary["raw_RNA_count_input_qualified"] = ready
    summary["qualification_status"] = "original_RNA_qualified_root_analysis_freeze_pending" if ready else "original_RNA_margins_unresolved"
    summary["software"]["script_sha256"] = sha256(__file__)
    dump_json(output, summary)
    return summary



def qualify_methods(raw=RAW):
    """Record bounded source declarations; never infer the eight donor covariates."""
    import xml.etree.ElementTree as ET
    import zipfile

    raw = Path(raw)
    with zipfile.ZipFile(raw / "A900-methods.docx") as archive:
        document = ET.fromstring(archive.read("word/document.xml"))
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    paragraphs = document.findall("w:body/w:p", ns)
    paragraph = "".join(paragraphs[15].itertext())
    declarations = [
        "Cell Ranger 5.0.0", "filtering cells with less than 200 nCounts or greater than 40% mitochondrial reads",
        "All samples were then normalized by SCTransform (v2) and integrated using Harmony (v0.1.1)",
        "Code for all analyses is available upon request."]
    if any(value not in paragraph for value in declarations):
        raise ValueError("Pinned Methods paragraph does not support recorded declarations")
    return {
        "article": {"pmcid": "PMC11150034", "doi": "10.1097/HC9.0000000000000454"},
        "source_methods_supplement": {"file": "A900-methods.docx", "url": "https://links.lww.com/HC9/A900",
                                      "sha256": sha256(raw / "A900-methods.docx"),
                                      "archive_member": "word/document.xml", "xpath": "/w:document/w:body/w:p[16]",
                                      "verified_method_fragments": declarations},
        "passage": 3, "intervention": {"source_treatment": "IL17", "molecule": "human recombinant IL-17A",
                                       "dose_ng_per_mL": 100, "duration_hours": 24,
                                       "control_source_treatment": "CNT", "control_description": "vehicle",
                                       "vehicle_composition": None},
        "genome_assembly_source_declaration": "hg38", "exact_reference_bundle": None,
        "gene_annotation_release": None, "allele_convention": "not applicable; no variant analysis",
        "pipeline_source_declarations": {"Cell_Ranger": "5.0.0", "Seurat": "v4", "SCTransform": "v2 (literal source wording)",
                                         "Harmony": "v0.1.1", "executed_commands_or_container": None},
        "source_cell_filter": {"exclude_if_nCounts_less_than": 200, "exclude_if_mitochondrial_percent_greater_than": 40,
                               "quantity_is_nCounts_not_nFeatures": True, "other_filters": None,
                               "historical_filter_reconstruction_performed": False},
        "historical_feature_prefilters": None, "merged_export_assay_slot_command": None,
        "original_count_export_command": None, "study_code": "Available upon request in inspected sources; not requested",
        "clinical_covariates": {"status": "no independently established sequencing-label to clinical-row crosswalk",
                                "sequenced_donor_keys": 8, "larger_clinical_table_PSC_rows": 10,
                                "larger_clinical_table_nonPSC_rows": 7,
                                "donors_with_independently_mapped_age_sex_source_site_cirrhosis": 0,
                                "unassigned_fields": ["age", "sex", "IBD", "cirrhosis", "procedure indication", "anatomical sampling site", "bile versus brush"],
                                "numeric_suffix_to_table_row_mapping_used": False,
                                "scope": "GEO metadata, article Methods/Table 1, and directly linked supplemental Methods; no claim about all uninspected material"},
        "control_group": "Non-PSC clinically indicated ERC procedure controls, not healthy volunteers",
        "sampling": {"source_context": "Extrahepatic bile duct organoids from brush and/or bile samples",
                     "article_declares_one_bile_derived_organoid_in_scRNA_or_NanoString_set": True,
                     "bile_donor_key": None,
                     "source_tension": "Text says all controls were brush-derived, while larger Table 1 non-PSC patient 6 says Common bile duct, Bile: cystic duct; neither is assigned to a sequencing key"},
        "limitations": ["Eight donors, not sixteen independent donor pairs; each donor contributes one CNT and one IL17 library",
                        "Number of culture wells/clones pooled within a source library is not established",
                        "Treatment-assignment or sequencing-batch randomization is not established in bounded sources",
                        "No invented age/sex/site/cirrhosis covariates or exclusions",
                        "Disease-versus-control differences remain observational despite the ex vivo treatment intervention"]}


def render_report(summary, report=REPORT):
    """Aggregate source/QC report; no cell-level joins or gene expression values."""
    matrix, raw = summary["matrix"], summary.get("original_RNA")
    if raw is None:
        raise ValueError("Complete original input qualification is required for final report")
    joins, quantity = summary["joins"], summary["quantity_validation"]
    original_quantity = raw["quantity_validation"]
    sections = [
        "# GSE239283 PSC organoid input qualification",
        "",
        "## Status and scope",
        "",
        "Source acquisition, numeric checks and donor–library joins are complete. No gene-wise test, target inspection, program score, treatment contrast or model fit was performed. An analysis freeze remains the root agent’s responsibility.",
        "",
        f"The design contains **{joins['donors']} donor keys** ({joins['donors_per_disease']['PSC']} PSC and {joins['donors_per_disease']['nonPSC']} non-PSC), **{joins['complete_donor_treatment_pairs']} complete donor-treatment pairs**, and **{joins['libraries']} source libraries**. Sixteen libraries are not sixteen donor pairs. Each donor contributes CNT and IL17 at passage 3. GEO and supplemental Methods specify 100 ng/mL recombinant human IL-17A for 24 hours versus vehicle; vehicle composition is not stated.",
        "",
        "## Critical assay distinction",
        "",
        f"The merged integer matrix has {matrix['features']:,} features × {matrix['cells']:,} cells and {matrix['observed_entries']:,} entries. Every cell’s total and positive-feature count exactly matches **nCount_SCT and nFeature_SCT**. Only {quantity['nCount_RNA']['matching_cells']:,} cell totals match nCount_RNA. Its total is {matrix['full_released_matrix_sum']:,}; the metadata RNA total is {quantity['nCount_RNA']['source_sum']:,}. The merged matrix must **not** be used as unnormalized RNA input to a count model. Integer values do not establish raw RNA provenance.",
        "",
        "Methods state SCTransform (v2) and Harmony (v0.1.1). The margin match supports an SCT-assay-consistent quantity, not recovery of the exact historical assay/slot export command. The complete merged input and its separately labeled library sums are preserved.",
        "",
        "## Original processed RNA exports",
        "",
        f"Under the expanded cumulative 1-GB processed-input authorization, the complete {raw['archive_member_ledger']['file'].split('/')[-1]} ledger accounts for 48 approved GEO TAR members: one matrix, feature list and barcode list per source library. The archive contains processed 10x count exports, not FASTQ reads. All members are regular basename-only files. Extraction rejects extra/duplicate members, path traversal, absolute paths, symlinks, hardlinks and unexpected sizes.",
        "",
        f"All 16 matrices have the same **{raw['full_feature_count']:,} unique Ensembl gene IDs**. Their {raw['unique_source_gene_labels']:,} distinct gene labels are retained separately; duplicate symbols are not collapsed. The original exports contain {raw['full_original_cells']:,} cells and {raw['full_original_stored_entries']:,} stored entries. Every released cell maps exactly once within its stated source library after a separately recorded technical merge-suffix removal. The analysis-preparation subset contains all **{raw['retained_cells']:,} released cells**. The other {raw['original_cells_not_in_released_metadata']:,} source cells remain outside that published subset, without a newly invented exclusion rule.",
        "",
        f"The all-feature retained-cell count total is **{raw['retained_count_total']:,}**, versus **{original_quantity['nCount_RNA']['source_sum']:,}** in historical nCount_RNA metadata: a difference of **{raw['retained_count_total'] - original_quantity['nCount_RNA']['source_sum']:,} counts** and **{raw['retained_positive_entries'] - original_quantity['nFeature_RNA']['source_sum']:,} positive entries**. Totals and positive-feature counts both match exactly in {original_quantity['nCount_RNA']['matching_cells']:,} cells and differ in {original_quantity['nCount_RNA']['differing_cells']:,}. Original-minus-metadata ranges are {original_quantity['nCount_RNA']['minimum_matrix_minus_metadata']} to {original_quantity['nCount_RNA']['maximum_matrix_minus_metadata']} counts and {original_quantity['nFeature_RNA']['minimum_matrix_minus_metadata']} to {original_quantity['nFeature_RNA']['maximum_matrix_minus_metadata']} features per cell. No source value was changed. No guessed historical gene filter was applied to force agreement. This is an unresolved historical feature-universe/export difference, not evidence that the integer SCT export was raw RNA.",
        "",
        f"Before effects, the root agent accepted the **complete original RNA feature universe for the same published retained cells**. This is a new released-raw-universe estimand, not exact reconstruction of the historical filtered RNA matrix. The small metadata discrepancy does not prove its cause. All {raw['full_feature_count']:,} unique source IDs are retained; {raw['duplicate_source_gene_label_groups']} duplicated symbol labels cover {raw['feature_rows_with_duplicated_source_gene_labels']} separate feature rows and are not collapsed. The JSON eligibility decision binds the exact count, feature, library and cell-join hashes and independent checks. A separate final root analysis freeze is still required; this qualification does not authorize effects.",
        "",
        "## Global and library QC",
        "",
        f"The merged matrix and all original matrices pass gzip CRC/EOF, dimensions, finite nonnegative integer entries, coordinate bounds and duplicate checks. Neither contains explicit stored zeros. Every feature and released cell is retained. The original-RNA subset has {raw['positive_total_features_in_retained_cells']:,} positive-total features and {raw['zero_total_features_in_retained_cells']:,} zero-total features; the merged SCT matrix has {matrix['zero_total_features']:,} zero-total features. These are availability counts, not gene filters or biological nulls.",
        "",
        "| Source library | Released cells | Original cells | Full-feature RNA total in released cells |",
        "|---|---:|---:|---:|"]
    for row in raw["libraries"]:
        sections.append(f"| {row['library_id']} | {row['n_cells']:,} | {row['source_original_cells']:,} | {row['total_released_counts']:,} |")
    sections += ["", "| Condition | Treatment | Donors | Cells | Full-feature RNA total |", "|---|---|---:|---:|---:|"]
    for row in raw["conditions"]:
        sections.append(f"| {row['disease']} | {row['treatment']} | {row['donors']} | {row['cells']:,} | {row['total_released_counts']:,} |")
    sections += [
        "", "Full released-matrix library totals are computed before any future feature filter. Cell barcodes, source joins and feature-by-library expression remain in ignored `work/psc-organoid-inputs/`.",
        "", "## Identity, annotation and release limits", "",
        "Non-PSC donors are clinically indicated ERC procedure controls, not healthy volunteers. The larger clinical Table 1 covers 10 PSC and 7 non-PSC patients. The inspected GEO, article Methods/Table 1 and supplemental Methods do not independently crosswalk its rows to the eight sequencing labels. Age, sex, IBD, cirrhosis, exact source site and brush/bile identity therefore remain unassigned. Numeric donor suffixes were not treated as table-row identifiers. One bile-derived organoid is described without its sequencing key. The text/table brush-versus-bile tension is preserved, not resolved by assumption.",
        "",
        "GEO declares hg38 and Cell Ranger 5.0.0. Methods declare Seurat v4, SCTransform v2 and Harmony v0.1.1. The exact genome reference bundle, annotation release, chemistry subversion, feature prefilters and executed export commands are not established. Source labels remain unchanged; the original Ensembl IDs are not replaced by current aliases. The merged gene list contains two identical symbol-style columns rather than stable Ensembl IDs.",
        "",
        "The supplemental cell-QC rule is **nCounts < 200 or mitochondrial reads > 40%**, not nFeatures < 200. It is not assumed to explain all absent cells or the historical RNA feature difference. No new cell, gene or donor exclusions were applied. Code is offered on author request in the inspected sources; no request or contact was made.",
        "",
        "## Reproduction and files", "",
        "```text", ".venv/bin/python scripts/qualify_psc_organoid_inputs.py --acquire",
        ".venv/bin/python scripts/qualify_psc_organoid_inputs.py --acquire-original",
        ".venv/bin/python scripts/qualify_psc_organoid_inputs.py --qualify --qualify-original --finalize",
        ".venv/bin/python -m unittest tests.test_psc_organoid_inputs -v", "```", "",
        "- Aggregate machine-readable result: `data/derived/psc-organoid-input-qualification.json`.",
        "- Source URLs, original retrieval times, byte counts and SHA-256: `data/raw/psc-organoid-il17/source-manifest.json` and the public JSON’s `sources` field.",
        "- Safe member inventory/hashes: `data/raw/psc-organoid-il17/original-archive-members.json`.",
        "- Raw-source matrix: `work/psc-organoid-inputs/raw-rna-pseudobulk.tsv.gz`, with `source_feature_id` then 16 exact source-library columns.",
        "- Axis and library/QC tables: `raw-rna-features.tsv` and `raw-rna-libraries.tsv` in the same ignored directory.",
        "- Cell crosswalk and independent verification stay ignored. The public report contains only study/library aggregates.",
        "",
        f"Acquired/copied source bodies total {summary['acquisition']['all_source_bytes']:,} bytes under the {summary['acquisition']['acquisition_limit_bytes']:,}-byte limit. No raw reads, controlled-access data, researcher contact, model prediction or clinical recommendation is involved.",
        "",
        "The sample size remains four donors per group. More cells do not create more independent donors. This stage establishes input provenance and numeric integrity, not an IL-17 disease mechanism or treatment benefit.", ""]
    validation = summary.get("validation", {})
    if validation:
        tests = validation["synthetic_tests"]
        repo = validation["repository_tests"]
        independent = validation["independent_verification"]
        replay = validation["offline_replay"]
        sections += ["## Validation", "",
                     f"All {tests['tests_run']} organoid synthetic tests pass. The independent SciPy sparse implementation reproduces all {independent['pseudobulk_count_values']:,} feature-by-library values and both per-cell RNA margins, with {independent['total_mismatches']} mismatches across {independent['total_comparisons']:,} comparisons. A second full pass uses different chunk boundaries; {replay['byte_identical_data_tables']} count/axis/join tables reproduce byte for byte.",
                     "",
                     f"The required full repository suite ran {repo['tests_run']} tests and reported {repo['errors']} errors from unavailable pre-existing ETS2/GSE84161 source caches or its separate Rscript runtime. These files were not downloaded, installed or altered. The organoid tests themselves pass; a claim that the entire repository suite passed would be incorrect.", ""]
    Path(report).write_text("\n".join(sections))



def apply_root_scope_decision(summary, work=WORK):
    """Qualify the full original count universe, not a historical filtered assay.

    Root accepted this scope before effects in the named message. This function
    never authorizes inference: the separate final root analysis freeze is still
    required. Independent checks and exact artifact pins must pass first.
    """
    work = Path(work)
    raw = summary["original_RNA"]
    required = ["raw-rna-pseudobulk.tsv.gz", "raw-rna-features.tsv", "raw-rna-libraries.tsv", "raw-rna-cell-source-joins.tsv.gz"]
    pinned = {row["file"]: row for row in raw["work_artifacts"]}
    bindings = []
    for name in required:
        receipt = pinned[name]
        file = work / name
        if file.stat().st_size != receipt["bytes"] or sha256(file) != receipt["sha256"]:
            raise ValueError("Root raw-universe decision cannot bind changed count/axis/join artifacts")
        bindings.append({"path": "work/psc-organoid-inputs/" + name,
                         "bytes": receipt["bytes"], "sha256": receipt["sha256"]})
    verification = summary.get("validation", {})
    independent = verification.get("independent_verification", {})
    replay = verification.get("offline_replay", {})
    checks = {
        "all_original_processed_members_verified": raw["compressed_members_checked"] == 48 and raw["all_member_gzip_crc_and_eof_verified"],
        "all_released_cells_joined_once": raw["retained_cells"] == 46343 and raw["barcode_join"]["all_released_cells_found_once"],
        "complete_original_unique_Ensembl_axis": raw["full_feature_count"] == raw["unique_feature_IDs"] == raw["Ensembl_gene_ID_rows"] == 36601 and raw["all_feature_axes_exactly_identical"],
        "all_original_numeric_checks_passed": all(raw["full_sparse_QC"][key] == 0 for key in ["duplicate_coordinates", "negative_entries", "nonfinite_entries", "noninteger_entries", "out_of_bounds_coordinates"]),
        "eight_complete_source_donor_pairs": summary["joins"]["complete_donor_treatment_pairs"] == 8 and summary["joins"]["libraries"] == 16,
        "independent_full_pseudobulk_verification_passed": independent.get("status") == "passed" and independent.get("total_mismatches") == 0 and independent.get("pseudobulk_count_values") == 36601 * 16,
        "different_chunk_replay_passed": replay.get("status") == "passed" and replay.get("count_and_numeric_QC_mismatches") == 0,
        "merged_SCT_is_not_the_raw_count_input": summary["quantity_validation"]["interpretation"]["SCT_margins_match"] and not summary["quantity_validation"]["interpretation"]["raw_RNA_margins_match"],
    }
    eligible = all(checks.values())
    raw["raw_RNA_count_input_qualified"] = eligible
    raw["source_count_input_verified"] = eligible
    raw["historical_RNA_metadata_exactly_reproduced"] = raw["quantity_validation"]["interpretation"]["raw_RNA_margins_match"]
    with (work / "raw-rna-features.tsv").open(newline="") as feature_stream:
        feature_labels = Counter(row["source_gene_label"] for row in csv.DictReader(feature_stream, delimiter="\t"))
    raw["duplicate_source_gene_label_groups"] = sum(value > 1 for value in feature_labels.values())
    raw["feature_rows_with_duplicated_source_gene_labels"] = sum(value for value in feature_labels.values() if value > 1)
    summary["raw_RNA_count_input_qualified"] = eligible
    summary["ready_for_raw_count_inference"] = False
    summary["qualification_status"] = "original_full_RNA_universe_qualified_root_final_freeze_pending" if eligible else "root_raw_universe_conditions_incomplete"
    summary["eligibility_decision"] = {
        "decision": "eligible_original_full_count_universe" if eligible else "pending_qualification_checks",
        "root_scope_authorization_message_id": "agentmsg_fd99b14b-415e-4065-8759-5a3e6cc709e9",
        "accepted_before_gene_effects": True,
        "scope": "Complete original Cell Ranger processed RNA feature universe for exactly the published retained cells; new explicit released-raw-universe estimand",
        "not_claimed": "Exact historical RNA-matrix/filter reconstruction",
        "quantity_provenance": "GEO per-GSM processed count-matrix/feature/barcode exports; GEO says Cell Ranger 5.0.0 performed demultiplexing, barcode processing and gene counting; Gene Expression feature types, integer sparse checks and original-cell joins verified",
        "historical_RNA_margin_match_required_for_this_new_scope": False,
        "historical_discrepancy_cause_established": False,
        "counts_above_historical_RNA_metadata": raw["retained_count_total"] - raw["quantity_validation"]["nCount_RNA"]["source_sum"],
        "positive_entries_above_historical_RNA_metadata": raw["retained_positive_entries"] - raw["quantity_validation"]["nFeature_RNA"]["source_sum"],
        "source_manifest_sha256": summary["source_manifest"]["sha256"],
        "input_artifact_hashes": bindings, "checks": checks,
        "merged_SCT_raw_count_model_eligible": False,
        "symbol_collapse_performed": False, "historical_feature_filter_guessed": False,
        "root_final_analysis_freeze_still_required": True,
        "effects_authorized_by_this_qualification": False}
    return summary


def finalize(raw=RAW, work=WORK, output=OUTPUT, report=REPORT):
    summary = json.loads(Path(output).read_text())
    summary["methods_qualification"] = qualify_methods(raw)
    summary["software"]["script_sha256"] = sha256(__file__)
    summary["software"]["tests"] = "tests/test_psc_organoid_inputs.py"
    summary["software"]["tests_sha256"] = sha256(ROOT / "tests/test_psc_organoid_inputs.py")
    summary["remaining_provenance_limits"] = [
        "Historical nCount_RNA/nFeature_RNA feature universe not reproduced; complete original feature axis retained without guessing a filter",
        "Exact reference/annotation release and executed export commands unavailable in inspected sources",
        "Exact age/sex/source-site/cirrhosis crosswalk for the eight sequenced donors not independently established",
        "Source library captures a donor-treatment culture; independent clone/well pooling and randomization are not fully specified",
        "Root prospective analysis freeze required before any gene effects"]
    verification_path = WORK / "validation-summary.json"
    if verification_path.exists():
        summary["validation"] = json.loads(verification_path.read_text())
    summary = apply_root_scope_decision(summary, work)
    dump_json(output, summary)
    render_report(summary, report)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", action="store_true")
    parser.add_argument("--qualify", action="store_true")
    parser.add_argument("--qualify-original", action="store_true")
    parser.add_argument("--finalize", action="store_true")
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--acquire-original", action="store_true")
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--chunk-lines", type=int, default=131072)
    parser.add_argument("--work-dir", type=Path, default=WORK)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.acquire:
        acquire()
    if args.acquire_original:
        acquire_original(inventory_only=args.inventory_only)
    if args.qualify:
        qualify(work=args.work_dir, output=args.output, chunk_lines=args.chunk_lines)
    if args.qualify_original:
        qualify_original(work=args.work_dir, output=args.output, chunk_lines=args.chunk_lines)
    if args.finalize:
        finalize(work=args.work_dir, output=args.output, report=args.report)
    if not any([args.acquire, args.qualify, args.acquire_original, args.qualify_original, args.finalize]):
        parser.error("Use acquisition and/or qualification flags")


if __name__ == "__main__":
    main()
