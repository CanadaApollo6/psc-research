#!/usr/bin/env python3
"""Bounded, source-only BACH2 input qualification. No expression analysis.

Run in the repository environment. `acquire` follows public redirects manually,
never reads HEAD bodies, and refuses GET for non-metadata/annotation resources.
`qualify` (offline) checks joins and writes only aggregate public results.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urljoin, urlsplit
from urllib.request import Request, HTTPRedirectHandler, build_opener

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/psc-bach2-qualification"
WORK = ROOT / "work/psc-bach2-qualification"
BUDGET = 15_000_000
ALLOWED_HOSTS = {"www.ebi.ac.uk", "ftp.ebi.ac.uk", "europepmc.org", "pmc.ncbi.nlm.nih.gov", "www.ncbi.nlm.nih.gov", "cdn.ncbi.nlm.nih.gov", "ars.els-cdn.com"}
GET_KINDS = {"sample_metadata", "processed_cell_annotation", "source_supplement", "source_metadata"}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def json_write(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def read_ledger(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [json.loads(x) for x in path.read_text().splitlines() if x.strip()]


def validate_url(url: str) -> None:
    parsed = urlsplit(url)
    if parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS:
        raise ValueError("Only allowlisted ordinary public HTTPS endpoints are permitted")
    if parsed.username or parsed.password:
        raise ValueError("Credentials in URLs are prohibited")


def validate_request(item: dict) -> None:
    validate_url(item["url"])
    if item["method"] not in {"GET", "HEAD"}:
        raise ValueError("Only GET and HEAD are permitted")
    if not re.fullmatch(r"[a-zA-Z0-9_.-]+", item["name"]):
        raise ValueError("Unsafe cache name")
    if item["method"] == "GET":
        if item["kind"] not in GET_KINDS:
            raise ValueError("Count or read bodies are not authorized")
        if any(s in urlsplit(item["url"]).path.lower() for s in ("feature_bc_matrix", "count_table", ".fastq", ".bam")):
            raise ValueError("Measurement bodies are not authorized")
        if not 0 < item.get("max_bytes", 0) <= BUDGET:
            raise ValueError("GET must have a positive bounded body allowance")


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def acquire(batch: list[dict], raw_dir: Path = RAW) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    ledger_path = raw_dir / "request-ledger.jsonl"
    ledger = read_ledger(ledger_path)
    opener = build_opener(NoRedirect())
    for item in batch:
        validate_request(item)
        if any(x["name"] == item["name"] and x.get("complete") for x in ledger):
            print(json.dumps({"cached": item["name"]}), flush=True)
            continue
        remaining = BUDGET - sum(x["body_bytes"] for x in ledger)
        if item["method"] == "GET" and remaining <= 0:
            raise RuntimeError("Cumulative body byte budget exhausted")
        entry = dict(item, started_at_utc=utc_now(), complete=False,
                     error=None, requests=[], body_bytes=0)
        response = None
        parts = []
        headers = {}
        url = item["url"]
        try:
            for _ in range(6):
                validate_request(dict(item, url=url))
                req = Request(url, method=item["method"], headers={
                    "User-Agent": "PSC-public-input-qualification/1.0",
                    "Accept-Encoding": "identity",
                })
                try:
                    response = opener.open(req, timeout=25)
                except HTTPError as exc:
                    response = exc
                headers = dict(response.headers.items())
                hop = {"url": url, "method": item["method"],
                       "status": response.status, "headers": headers}
                entry["requests"].append(hop)
                if response.status in {301, 302, 303, 307, 308}:
                    location = response.headers.get("Location")
                    response.close()
                    response = None
                    if not location:
                        raise ValueError("Redirect lacks Location")
                    url = urljoin(url, location)
                    continue
                entry["status"] = response.status
                entry["response_url"] = response.geturl()
                if item["method"] == "HEAD":
                    entry["complete"] = response.status == 200
                    break
                if response.status != 200:
                    raise ValueError("Non-200 metadata body not read")
                content_length = response.headers.get("Content-Length")
                limit = min(item["max_bytes"], remaining)
                if content_length and int(content_length) > limit:
                    raise ValueError("Content-Length exceeds body allowance; no body read")
                if response.headers.get("Content-Encoding", "identity") != "identity":
                    raise ValueError("Unexpected transport encoding; no body read")
                nread = 0
                eof = False
                while nread < limit:
                    part = response.read(min(65_536, limit - nread))
                    if not part:
                        eof = True
                        break
                    parts.append(part)
                    nread += len(part)
                # At the cap without an exact declared length, completion is unknown.
                entry["complete"] = eof or (content_length is not None and nread == int(content_length))
                if not entry["complete"]:
                    raise ValueError("Body allowance reached before confirmed completion")
                expected = item.get("expected_bytes")
                if expected is not None and nread != expected:
                    entry["complete"] = False
                    raise ValueError("Body length differs from published inventory")
                break
            else:
                raise ValueError("Redirect limit exceeded")
        except Exception as exc:
            entry["error"] = type(exc).__name__ + ": " + str(exc)
        finally:
            if response is not None:
                response.close()
        data = b"".join(parts)
        ordinal = len(ledger) + 1
        body_name = f"{ordinal:02d}-{item['name']}.body"
        (raw_dir / body_name).write_bytes(data)
        entry.update(body_file=body_name, body_bytes=len(data), sha256=sha256(data),
                     completed_at_utc=utc_now())
        with ledger_path.open("a") as out:
            out.write(json.dumps(entry, sort_keys=True) + "\n")
        ledger.append(entry)
        print(json.dumps(entry, sort_keys=True), flush=True)
    print(json.dumps({"total_new_body_bytes": sum(x["body_bytes"] for x in ledger), "budget_bytes": BUDGET}), flush=True)



def read_tsv(text: str) -> tuple[list[str], list[list[str]]]:
    """Keep repeated SDRF headers; reject malformed rows rather than dropping fields."""
    if "\x00" in text:
        raise ValueError("NUL byte in text input")
    reader = csv.reader(io.StringIO(text), delimiter="\t")
    try:
        header = next(reader)
    except StopIteration:
        raise ValueError("Empty TSV")
    if not header or not all(header):
        raise ValueError("Empty header field")
    rows = list(reader)
    if any(len(row) != len(header) for row in rows):
        raise ValueError("TSV field count mismatch")
    return header, rows


def field_values(header: list[str], row: list[str], name: str) -> list[str]:
    return [row[i] for i, key in enumerate(header) if key == name]


def single_field(header: list[str], row: list[str], name: str,
                 required: bool = True) -> str | None:
    values = field_values(header, row, name)
    if not values and not required:
        return None
    if len(values) != 1 or not values[0]:
        raise ValueError(f"Expected one nonempty {name}")
    return values[0]


def parse_sdrf(text: str) -> tuple[list[str], list[dict]]:
    header, rows = read_tsv(text)
    records = []
    for row in rows:
        record = {
            "source_name": single_field(header, row, "Source Name"),
            "individual": single_field(header, row, "Characteristics[individual]"),
            "sex": single_field(header, row, "Characteristics[sex]"),
            "disease": single_field(header, row, "Characteristics[disease]"),
            "organism_part": single_field(header, row, "Characteristics[organism part]"),
            "assay_name": single_field(header, row, "Assay Name"),
            "genotype": single_field(header, row, "Characteristics[genotype]", False),
            "genotype_factor": single_field(header, row, "Factor Value[genotype]", False),
            "age": single_field(header, row, "Characteristics[age]", False),
            "age_unit": single_field(header, row, "Unit[time unit]", False),
            "files": field_values(header, row, "Derived Array Data File"),
        }
        if record["genotype"] != record["genotype_factor"]:
            raise ValueError("Genotype characteristic/factor disagreement")
        if not record["files"] or not all(record["files"]):
            raise ValueError("Missing linked processed files")
        records.append(record)
    for key in ("source_name", "individual", "assay_name"):
        values = [r[key] for r in records]
        if len(values) != len(set(values)):
            raise ValueError(f"Nonunique {key}: unit of analysis is unresolved")
    return header, records


def file_inventory(study: dict) -> dict[str, dict]:
    """Traverse original BioStudies file nodes, rejecting ambiguous repeated paths."""
    found = {}
    def visit(node):
        if isinstance(node, dict):
            if node.get("type") == "file" and "path" in node:
                path = node["path"]
                if path in found:
                    raise ValueError("Duplicate file inventory path")
                found[path] = node
            for value in node.values():
                visit(value)
        elif isinstance(node, list):
            for value in node:
                visit(value)
    visit(study)
    return found


def sample_interfaces(records: list[dict], inventory: dict[str, dict]) -> dict[str, dict]:
    by_sample = {}
    for record in records:
        for filename in record["files"]:
            if filename not in inventory:
                raise ValueError("SDRF file absent from source inventory")
        filtered = [x for x in record["files"] if x.endswith("_filtered_feature_bc_matrix.tar.gz")]
        raw = [x for x in record["files"] if x.endswith("_raw_feature_bc_matrix.tar.gz")]
        if len(filtered) != 1 or len(raw) != 1 or "scrna_meta.tsv" not in record["files"]:
            raise ValueError("One raw, one filtered, and the cell annotation file required")
        sample = filtered[0].removesuffix("_filtered_feature_bc_matrix.tar.gz")
        if raw[0] != sample + "_raw_feature_bc_matrix.tar.gz":
            raise ValueError("Raw/filtered file prefix disagreement")
        if sample in by_sample:
            raise ValueError("Ambiguous archive/sample prefix")
        by_sample[sample] = dict(record, metadata_sample_name=sample,
                                 filtered_archive=filtered[0], raw_archive=raw[0])
    return by_sample


def summarize_cells(text: str, by_sample: dict[str, dict]) -> tuple[dict, dict]:
    header, rows = read_tsv(text)
    if len(header) != len(set(header)):
        raise ValueError("Repeated cell metadata header")
    required = {"sample_name", "sex", "condition", "seurat_clusters",
                "paper_clusters", "paper_clusters_short", "barcode"}
    if not required.issubset(header):
        raise ValueError("Cell annotation fields missing")
    lookup = {name: i for i, name in enumerate(header)}
    by_donor_state = Counter()
    donor_n = Counter()
    observed_samples = set()
    full_barcodes = set()
    stripped_keys = set()
    bare_barcode_counts = Counter()
    state_labels = defaultdict(set)
    seurat_to_paper = defaultdict(set)
    paper_to_seurat = defaultdict(set)
    missingness = Counter()
    conditions = Counter()
    sexes = Counter()
    for row in rows:
        r = {name: row[lookup[name]] for name in required}
        missingness.update(name for name, value in r.items() if not value)
        if not all(r.values()):
            raise ValueError("Missing required cell annotation")
        sample = r["sample_name"]
        if sample not in by_sample:
            raise ValueError("Cell sample has no explicitly linked SDRF archive")
        info = by_sample[sample]
        if r["condition"] != info["genotype"] or r["sex"] != info["sex"]:
            raise ValueError("Cell annotation/SDRF genotype or sex disagreement")
        barcode = r["barcode"]
        prefix = sample + "_"
        if not barcode.startswith(prefix):
            raise ValueError("Barcode lacks exact sample prefix")
        bare = barcode[len(prefix):]
        if not re.fullmatch(r"[ACGT]{16}-[0-9]+", bare):
            raise ValueError("Unrecognized source barcode form")
        if barcode in full_barcodes or (sample, bare) in stripped_keys:
            raise ValueError("Duplicate within-sample barcode")
        full_barcodes.add(barcode)
        stripped_keys.add((sample, bare))
        bare_barcode_counts[bare] += 1
        observed_samples.add(sample)
        donor_n[sample] += 1
        conditions[r["condition"]] += 1
        sexes[r["sex"]] += 1
        state = r["paper_clusters_short"]
        label = r["paper_clusters"]
        if not label.startswith(state + ": "):
            raise ValueError("Short/full paper cluster labels disagree")
        state_labels[state].add(label)
        seurat_to_paper[r["seurat_clusters"]].add(state)
        paper_to_seurat[state].add(r["seurat_clusters"])
        by_donor_state[sample, state] += 1
    if observed_samples != set(by_sample):
        raise ValueError("SDRF samples absent from retained cell metadata")
    if any(len(v) != 1 for v in state_labels.values()):
        raise ValueError("One short state has multiple labels")
    if any(len(v) != 1 for v in seurat_to_paper.values()) or any(len(v) != 1 for v in paper_to_seurat.values()):
        raise ValueError("Seurat/paper crosswalk is not one-to-one")
    genotype_donors = Counter(r["genotype"] for r in by_sample.values())
    states = []
    for state in sorted(state_labels, key=lambda x: (int(x[1:]), x)):
        counts = [by_donor_state[sample, state] for sample in by_sample]
        eligible = {}
        for threshold in (20, 50):
            eligible[str(threshold)] = {
                genotype: sum(by_donor_state[sample, state] >= threshold
                              for sample, info in by_sample.items() if info["genotype"] == genotype)
                for genotype in sorted(genotype_donors)
            }
        states.append({
            "paper_clusters_short": state,
            "paper_clusters": next(iter(state_labels[state])),
            "seurat_clusters": sorted(paper_to_seurat[state]),
            "retained_cells": sum(counts),
            "donors_with_cells": sum(n > 0 for n in counts),
            "min_cells_per_donor": min(counts), "max_cells_per_donor": max(counts),
            "donors_at_cell_count_threshold": eligible,
        })
    summary = {
        "metadata_columns": header,
        "retained_cells": len(rows),
        "retained_samples": len(observed_samples),
        "retained_donors": len({r["individual"] for r in by_sample.values()}),
        "donors_by_source_genotype": dict(sorted(genotype_donors.items())),
        "retained_cells_by_source_genotype": dict(sorted(conditions.items())),
        "sex_counts_cells": dict(sorted(sexes.items())),
        "missing_required_fields": {key: missingness[key] for key in sorted(required)},
        "barcode_join": {
            "rule": "derive sample_name from explicitly linked filtered archive basename; remove only sample_name + '_' from metadata barcode",
            "full_barcodes_unique": len(full_barcodes) == len(rows),
            "sample_plus_stripped_barcode_unique": len(stripped_keys) == len(rows),
            "distinct_bare_barcodes": len(bare_barcode_counts),
            "bare_barcodes_shared_across_samples": sum(n > 1 for n in bare_barcode_counts.values()),
            "matrix_barcode_membership_verified": False,
            "matrix_retention_denominator_verified": False,
        },
        "states": states,
        "state_annotations_include_naive_activation_subclusters": False,
        "umap_values_parsed_or_interpreted": False,
        "counts_or_gene_features_inspected": False,
    }
    private = {
        "sample_to_source": by_sample,
        "donor_state_retained_cells": [
            {"sample_name": sample, "paper_clusters_short": state,
             "retained_cells": by_donor_state[sample, state]}
            for sample in sorted(by_sample) for state in sorted(state_labels)
        ],
    }
    return summary, private


def summarize_liver(records: list[dict]) -> dict:
    ages = [int(x["age"]) for x in records if x["age"] is not None]
    return {
        "source_samples": len(records),
        "distinct_accession_local_individual_labels": len({r["individual"] for r in records}),
        "sex_counts": dict(sorted(Counter(r["sex"] for r in records).items())),
        "age_nonmissing": len(ages),
        "age_min": min(ages) if ages else None,
        "age_max": max(ages) if ages else None,
        "age_units": sorted({r["age_unit"] for r in records if r["age_unit"] is not None}),
        "genotype_nonmissing": sum(r["genotype"] is not None for r in records),
        "processed_files": sorted({f for r in records for f in r["files"]}),
        "count_header_and_sample_columns_inspected": False,
        "cross_accession_participant_join_established": False,
        "genetic_replication_eligible": False,
        "scope": "metadata and count-file HEAD only",
    }



def parse_clinical_s2(text: str) -> tuple[dict, list[dict]]:
    """Only Table S2, not any source expression-figure text, is interpreted."""
    pages = text.split("\f")
    matches = [page for page in pages if "Table S2: Clinical parameters" in page]
    if len(matches) != 1 or "related to Figure 2A-J" not in matches[0]:
        raise ValueError("Exactly one clinical Table S2 / Figure 2A-J is required")
    fields = ["patient #", "SNP-status", "gender", "year of birth", "IBD status",
              "other autoimmune disease", "immunosuppr. therapy",
              "transient elastography [kPa]", "MELD-Score", "ANAs", "SMA", "SLA/LP"]
    rows = []
    for line in matches[0].splitlines():
        if not re.match(r"\s*patient #[0-9]+\s", line):
            continue
        values = re.split(r"\s{2,}", line.strip())
        if len(values) != len(fields):
            raise ValueError("Unexpected Table S2 text extraction fields")
        rows.append(dict(zip(fields, values)))
    if len(rows) != 8 or len({r["patient #"] for r in rows}) != 8:
        raise ValueError("Table S2 must retain all eight unique source rows")
    groups = Counter(r["SNP-status"] for r in rows)
    if groups != {"homozygous": 4, "non-carrier": 4}:
        raise ValueError("Clinical source genotype group labels/counts differ")
    summary = {}
    for group in sorted(groups):
        selected = [r for r in rows if r["SNP-status"] == group]
        info = {"source_rows": len(selected)}
        for field in ("gender", "IBD status", "other autoimmune disease", "immunosuppr. therapy"):
            info[field] = dict(sorted(Counter(r[field] for r in selected).items()))
        for field in ("year of birth", "transient elastography [kPa]", "MELD-Score"):
            values = sorted((r[field] for r in selected), key=lambda x: float(x.replace(",", ".")))
            info[field] = {"minimum_as_published": values[0], "maximum_as_published": values[-1]}
        summary[group] = info
    return {
        "source": "Document S1, Table S2, PDF page 8, related to Figure 2A-J",
        "groups_as_published": summary,
        "slash_meaning": "Literal '/' retained; no legend in Table S2 defines it as absence versus unavailable/not applicable.",
        "number_conversion": "Decimal comma changed to point for ordering ranges only. Source range endpoint strings remain unchanged.",
        "birth_year_is_not_age_at_sampling": True,
        "clinical_patient_to_accession_individual_join": "Numbering and genotype group pattern agree, but explicit cross-resource patient/sample key is not supplied. Candidate only; not merged.",
        "clinical_group_imbalance_known": True,
        "remaining_unresolved": ["exact sample-linked clinical crosswalk", "sampling dates/age at collection",
                                  "ancestry", "technical run/batch and genotype-batch balance",
                                  "complete medications/dose/timing and interpretation of '/'",
                                  "clinical matching/selection and participant overlap"],
    }, rows


def study_protocols(study: dict) -> dict[str, dict]:
    protocols = {}
    def visit(node):
        if isinstance(node, dict):
            if node.get("type") == "Protocols":
                attr = {x["name"]: x.get("value") for x in node["attributes"]}
                protocols[node["accno"]] = {"type": attr.get("Type"), "description_as_published": attr.get("Description")}
            for value in node.values():
                visit(value)
        elif isinstance(node, list):
            for value in node:
                visit(value)
    visit(study)
    return protocols


def response_evidence_status(entry: dict, data: bytes) -> str:
    if not entry["complete"]:
        return "request_incomplete_or_failed"
    if entry["method"] == "HEAD":
        return "headers_only_not_content_validation"
    if entry["kind"] == "source_supplement":
        if data.startswith(b"%PDF-") and b"%%EOF" in data[-1024:]:
            return "complete_pdf_source"
        return "not_requested_supplement_format_excluded"
    return "complete_source_response"


def verified_source(name: str, ledger: list[dict], raw_dir: Path) -> bytes:
    matches = [x for x in ledger if x["name"] == name and x["complete"] and x["method"] == "GET"]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one complete source: {name}")
    entry = matches[0]
    data = (raw_dir / entry["body_file"]).read_bytes()
    if len(data) != entry["body_bytes"] or sha256(data) != entry["sha256"]:
        raise ValueError("Source cache differs from receipt")
    return data


def verify_head_inventory(ledger: list[dict], inventories: dict[str, dict]) -> dict:
    output = {}
    for accession, inventory in inventories.items():
        files = []
        for name, info in sorted(inventory.items()):
            if not ("_filtered_feature_bc_matrix.tar.gz" in name or name == "count_table.txt"):
                continue
            matches = [x for x in ledger if x["name"] == name + "-head"]
            if len(matches) != 1:
                raise ValueError("Missing or repeated count HEAD record")
            entry = matches[0]
            if entry["method"] != "HEAD" or entry["body_bytes"] != 0:
                raise ValueError("Count body acquisition is not allowed")
            headers = {k.lower(): v for k, v in entry["requests"][-1]["headers"].items()}
            length = int(headers["content-length"]) if "content-length" in headers else None
            files.append({"file": name, "published_bytes": info["size"],
                          "url": entry["url"], "resolved_public_url": entry.get("response_url"),
                          "head_status": entry.get("status"), "head_content_length": length,
                          "head_matches_inventory": entry["complete"] and length == info["size"],
                          "last_modified": headers.get("last-modified"),
                          "etag_not_a_content_hash": headers.get("etag"),
                          "body_downloaded": False, "content_sha256": None})
        output[accession] = {"files": files, "files_count": len(files),
                             "total_published_bytes": sum(f["published_bytes"] for f in files),
                             "all_head_lengths_match_inventory": all(f["head_matches_inventory"] for f in files)}
    return output


def qualify(raw_dir: Path = RAW, work_dir: Path = WORK,
            output_path: Path = ROOT / "data/derived/psc-bach2-input-qualification.json") -> dict:
    plan = json.loads((raw_dir / "plan.json").read_text())
    ledger = read_ledger(raw_dir / "request-ledger.jsonl")
    if sum(x["body_bytes"] for x in ledger) > BUDGET:
        raise ValueError("Response body byte budget exceeded")
    for entry in ledger:
        data = (raw_dir / entry["body_file"]).read_bytes()
        if len(data) != entry["body_bytes"] or sha256(data) != entry["sha256"]:
            raise ValueError("Acquisition cache differs from ledger")
    studies = {}
    for entry in plan["source_reuse"]:
        data = (ROOT / entry["path"]).read_bytes()
        if sha256(data) != entry["sha256"] or len(data) != entry["size_bytes"]:
            raise ValueError("Reused source differs from original receipt")
        if "metadata.body" in entry["path"]:
            study = json.loads(data)
            studies[study["accno"]] = study
    inventories = {acc: file_inventory(study) for acc, study in studies.items()}
    _, sc_records = parse_sdrf(verified_source("E-MTAB-14013.sdrf.txt", ledger, raw_dir).decode())
    liver_header, liver_records = parse_sdrf(verified_source("E-MTAB-14103.sdrf.txt", ledger, raw_dir).decode())
    by_sample = sample_interfaces(sc_records, inventories["E-MTAB-14013"])
    cell_summary, private = summarize_cells(verified_source("scrna_meta.tsv", ledger, raw_dir).decode(), by_sample)
    count_inventory = verify_head_inventory(ledger, inventories)
    text_provenance = json.loads((raw_dir / "supplement-text-provenance.json").read_text())
    pdf = verified_source("publisher-mmc1.pdf", ledger, raw_dir)
    if sha256(pdf) != text_provenance["source_pdf_sha256"] or not pdf.startswith(b"%PDF-"):
        raise ValueError("Supplement PDF provenance failed")
    clinical_text = (ROOT / text_provenance["text_path"]).read_bytes()
    if sha256(clinical_text) != text_provenance["text_sha256"]:
        raise ValueError("PDF text extraction differs from recorded derivative")
    clinical_summary, clinical_rows = parse_clinical_s2(clinical_text.decode())
    json_write(work_dir / "clinical-table-s2-source-values.json", clinical_rows)
    json_write(work_dir / "sample-donor-genotype-and-state-joins.json", private)
    json_write(work_dir / "liver-source-records.json", liver_records)
    source_pins = [{"name": x["name"], "path": x["path"], "url": x["url"],
                    "retrieved_at_utc": x["completed_at_utc"], "bytes": x["size_bytes"],
                    "sha256": x["sha256"], "reused": True}
                   for x in plan["source_reuse"]]
    source_pins += [{"name": x["name"], "path": str((raw_dir / x["body_file"]).relative_to(ROOT)),
                     "url": x["url"], "method": x["method"],
                     "retrieved_at_utc": x["completed_at_utc"], "bytes": x["body_bytes"],
                     "response_body_sha256": x["sha256"], "complete_response": x["complete"],
                     "http_status": x.get("status"), "error": x["error"], "reused": False,
                     "evidence_status": response_evidence_status(x, (raw_dir / x["body_file"]).read_bytes())}
                    for x in ledger]
    result = {
        "schema_version": 1,
        "scope": "Public source-only input qualification; no expression analysis",
        "decision": "E-MTAB-14013 metadata join qualifies; count and feature payload qualification requires separate authorization. E-MTAB-14103 genotype replication is blocked.",
        "ready_for_expression_analysis": False,
        "new_body_bytes": sum(x["body_bytes"] for x in ledger),
        "authorized_new_body_budget_bytes": BUDGET,
        "count_or_read_body_bytes": 0,
        "source_pins": source_pins,
        "count_file_interfaces": count_inventory,
        "E-MTAB-14013": dict(cell_summary, clinical_qualification=clinical_summary),
        "methods_by_accession_as_published": {acc: study_protocols(study) for acc, study in studies.items()},
        "supplement_text_provenance": text_provenance,
        "allele_conventions": {
            "variant_as_published": "rs56258221",
            "source_genotype_labels": ["SNP", "noSNP"],
            "paper_design_labels": ["homozygous carrier", "non-carrier"],
            "genotyping_assay_as_published": "TaqMan C__88670967_10",
            "REF": None, "ALT": None, "risk_nucleotide": None, "genomic_strand": None,
            "genomic_coordinate": None, "variant_build": None,
            "not_inferred_from_RNA_reference": True,
        },
        "E-MTAB-14103": dict(summarize_liver(liver_records), metadata_columns=liver_header),
        "qualification_limits": [
            "Source SNP/noSNP and homozygous/noncarrier labels do not provide nucleotide genotype, REF/ALT, strand or coordinate conventions.",
            "GRCh38 release 93 describes alignment, not a verified rs56258221 coordinate.",
            "All eight CITE-seq donors are male PSC cases. Table S2 reveals genotype-group clinical imbalance, but its patient-number to accession-individual crosswalk is unconfirmed and technical confounding remains unresolved.",
            "Source phenotype discovery and CITE-seq states are reused; this is not an independent cohort or isogenic genotype experiment.",
            "5-prime polyA RNA cannot validate a miRNA translation mechanism. ADT feature contents are not inspected.",
            "Filename conventions and HEAD lengths do not prove original integer UMI quantity, a complete feature universe, barcode membership, or intact archive contents.",
            "Methods 56,209 retained cells and Figure 2 / metadata 55,460 are both preserved; no silent reconciliation.",
            "The public scRNA annotation has 16 major states, not the four later naive activation subclusters.",
            "Metadata paper C15 MAIT maps to source Seurat 16, not 15. Metadata C9 CD4+ TC differs from the Methods C9 NKT label; do not silently correct either.",
            "CITE-seq Methods say TotalSeq A but the key resource table names TotalSeq-C. Actual deposited RNA/ADT feature contents are unverified.",
            "Bulk source Methods mention 47 sequencing-summary samples while the selected liver study has 12; the scopes are unreconciled.",
            "Bulk SDRF has one age 18 against Methods >18 inclusion, and library strand 'not applicable' versus Methods stranded mRNA. Preserve both.",
            "Numeric individual labels are accession-local; no cross-study donor or genotype mapping is inferred.",
        ],
        "next_authorization": {
            "purpose": "Archive integrity, full feature/assay type, original count quantity and retained-barcode membership only, before any separately frozen endpoint effects.",
            "filtered_archive_bytes": count_inventory["E-MTAB-14013"]["total_published_bytes"],
            "filtered_archive_count": count_inventory["E-MTAB-14013"]["files_count"],
            "permission_received": False,
            "not_in_scope": "No full count or raw-read bodies acquired in this qualification.",
        },
    }
    json_write(output_path, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    get = sub.add_parser("acquire")
    get.add_argument("batch", type=Path)
    qualify_parser = sub.add_parser("qualify")
    qualify_parser.add_argument("--output", type=Path, default=ROOT / "data/derived/psc-bach2-input-qualification.json")
    args = parser.parse_args()
    if args.command == "acquire":
        acquire(json.loads(args.batch.read_text()))
    elif args.command == "qualify":
        result = qualify(output_path=args.output)
        print(json.dumps({"decision": result["decision"], "new_body_bytes": result["new_body_bytes"]}))


if __name__ == "__main__":
    main()
