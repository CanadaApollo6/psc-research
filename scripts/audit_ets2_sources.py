"""Audit the public ETS2 benchmark sources and extract source-preserving tables.

This module deliberately stops at summary-level source verification.  It does
not run the released R code, calculate an atlas/program score, or reconstruct
any donor-level count matrix.  Workbook and rank values are retained as the
original strings in the derived tables; normalized Decimal values are used only
for deterministic comparisons and report summaries.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import posixpath
import re
import zipfile
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Iterable, Iterator
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/ets2-benchmark"
DERIVED = ROOT / "data/derived"
REPORTS = ROOT / "reports"

BASE_MANIFEST = ROOT / "config/ets2-benchmark-sources.json"
PROGRAM_MANIFEST = ROOT / "config/ets2-program-sources.json"

ARCHIVE_RELEVANT_MEMBERS = {
    "RNA-seq/RNAseq_CRISPR/TPP_CRISPR_RNAseq_analysis.r",
    "RNA-seq/RNAseq_CRISPR/ko.g1.res_finalSYMBOL_rank.csv",
    "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
    "RNA-seq/RNAseq_MEKi/PD_0.1_MEK.resSYMBOL.rank.csv",
    "RNA-seq/RNAseq_MEKi/PD_0.5_MEK.resSYMBOL.rank.csv",
    "RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv",
    "RNA-seq/RNAseq_biopsies/RNA-seq_colonbiopsies_analysis.R",
    "RNA-seq/RNAseq_overexpression/RNAseq_overexpression_analysis.r",
    "RNA-seq/RNAseq_overexpression/overexpression_250ng_rank.csv",
    "RNA-seq/RNAseq_overexpression/overexpression_500ng_rank.csv",
}

MASTER_MEMBER = "RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv"

MASTER_COLUMNS = [
    "ETS2g1_UP",
    "ETS2g1_DN",
    "CHR21_DN",
    "CHR21_UP",
    "ETS2g2_UP",
    "ETS2g2_DN",
]

S1_CONTRASTS = [
    {
        "id": "s1_g1",
        "label": "ETS2 gRNA 1",
        "left": "ETS2 gRNA 1",
        "right": "expression in unedited cells",
        "start_col": 3,
    },
    {
        "id": "s1_g2",
        "label": "ETS2 gRNA 2",
        "left": "ETS2 gRNA 2",
        "right": "expression in unedited cells",
        "start_col": 9,
    },
    {
        "id": "s1_chr21",
        "label": "chr21q22 deletion",
        "left": "chr21q22 deletion",
        "right": "expression in unedited cells",
        "start_col": 15,
    },
]

S2_CONTRASTS = [
    {
        "id": "s2_oe_250",
        "label": "ETS2 overexpression (250ng)",
        "left": "ETS2 overexpression (250ng)",
        "right": "expression in cells transfected with equivalent amount of mRNA encoding reverse complement of ETS2",
        "start_col": 3,
    },
    {
        "id": "s2_oe_500",
        "label": "ETS2 overexpression (500ng)",
        "left": "ETS2 overexpression (500ng)",
        "right": "expression in cells transfected with equivalent amount of mRNA encoding reverse complement of ETS2",
        "start_col": 9,
    },
    {
        "id": "s2_oe_combined",
        "label": "ETS2 overexpression (250ng and 500ng, dose as covariate)",
        "left": "ETS2 overexpression (250ng and 500ng, dose as covariate)",
        "right": "expression in cells transfected with equivalent amount of mRNA encoding reverse complement of ETS2",
        "start_col": 15,
    },
]

RANK_SPECS = [
    {
        "id": "ko_g1",
        "member": "RNA-seq/RNAseq_CRISPR/ko.g1.res_finalSYMBOL_rank.csv",
        "code_member": "RNA-seq/RNAseq_CRISPR/TPP_CRISPR_RNAseq_analysis.r",
        "source_workbook": "ets2-supplementary-tables.xlsx",
        "source_sheet": "Table S1_loss-of-function",
        "source_contrast": "s1_g1",
        "code_contrast": "g1 - NTC",
        "code_effect_direction": "ETS2 gRNA 1 edited minus non-targeting control",
        "source_comparison_status": "publisher_table_numeric_match_within_tolerance",
        "code_rank_linkage": (
            "R reads ko.g1.res_finalSYMBOL_rank.csv as fGSEA statistics; "
            "the article specifies t-statistic ranking; all unambiguous symbols "
            "match S1 gRNA1 t within 1e-9."
        ),
    },
    {
        "id": "oe_250",
        "member": "RNA-seq/RNAseq_overexpression/overexpression_250ng_rank.csv",
        "code_member": "RNA-seq/RNAseq_overexpression/RNAseq_overexpression_analysis.r",
        "source_workbook": "ets2-supplementary-tables.xlsx",
        "source_sheet": "Table S2_gain-of-function",
        "source_contrast": "s2_oe_250",
        "code_contrast": "ets2_250 - rev_250",
        "code_effect_direction": "ETS2 250ng mRNA minus reverse-complement 250ng control",
        "source_comparison_status": "publisher_table_nonexact_numeric_mismatch",
        "code_rank_linkage": (
            "R reads overexpression_250ng_rank.csv after the split-dose limma "
            "contrast and comments that it is ranked by t statistic; matched "
            "symbols are highly concordant with S2 but not numerically identical."
        ),
    },
    {
        "id": "oe_500",
        "member": "RNA-seq/RNAseq_overexpression/overexpression_500ng_rank.csv",
        "code_member": "RNA-seq/RNAseq_overexpression/RNAseq_overexpression_analysis.r",
        "source_workbook": "ets2-supplementary-tables.xlsx",
        "source_sheet": "Table S2_gain-of-function",
        "source_contrast": "s2_oe_500",
        "code_contrast": "ets2_500 - rev_500",
        "code_effect_direction": "ETS2 500ng mRNA minus reverse-complement 500ng control",
        "source_comparison_status": "publisher_table_nonexact_numeric_mismatch",
        "code_rank_linkage": (
            "R reads overexpression_500ng_rank.csv after the split-dose limma "
            "contrast and comments that it is ranked by t statistic; matched "
            "symbols are highly concordant with S2 but not numerically identical."
        ),
    },
    {
        "id": "meki_100",
        "member": "RNA-seq/RNAseq_MEKi/PD_0.1_MEK.resSYMBOL.rank.csv",
        "code_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "source_workbook": "",
        "source_sheet": "",
        "source_contrast": "",
        "code_contrast": "MEK_100 - ctrl",
        "code_effect_direction": "PD-0325901 100nM minus vehicle control",
        "source_comparison_status": "publisher_differential_table_not_pinned",
        "code_rank_linkage": (
            "R defines MEK_100nM = MEK_100 - ctrl, creates a topTable object, "
            "and then reads the released rank file. The matching write.csv line "
            "is commented/mismatched, so the saved rank's directional linkage "
            "cannot be certified from this release."
        ),
    },
    {
        "id": "meki_500",
        "member": "RNA-seq/RNAseq_MEKi/PD_0.5_MEK.resSYMBOL.rank.csv",
        "code_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "source_workbook": "",
        "source_sheet": "",
        "source_contrast": "",
        "code_contrast": "MEK_500 - ctrl",
        "code_effect_direction": "PD-0325901 500nM minus vehicle control",
        "source_comparison_status": "publisher_differential_table_not_pinned",
        "code_rank_linkage": (
            "R defines MEK_500nM = MEK_500 - ctrl, creates a topTable object, "
            "and then reads the released rank file. The matching write.csv line "
            "is commented/mismatched, so the saved rank's directional linkage "
            "cannot be certified from this release."
        ),
    },
]

DESIGN_SAMPLE_NAMES = {
    "crispr": [
        "1_g1", "1_g4", "2_chr21", "2_g1", "2_g4", "2_NTC",
        "3_chr21", "3_g1", "3_g4", "3_NTC", "4_chr21", "4_NTC",
        "5_chr21", "5_g1", "5_g4", "5_NTC", "6_g1", "6_g4", "6_NTC",
        "7_chr21", "7_g1", "7_g4", "7_NTC", "8_g1", "8_g1_ROX",
        "8_g4", "8_g4_ROX", "8_NTC", "9_g1", "9_g1_ROX", "9_g4",
        "9_g4_ROX", "9_NTC", "10_g1", "10_g1_ROX", "10_g4",
        "10_g4_ROX", "10_NTC",
    ],
    "oe": [
        f"{donor}_{rna}_{dose}"
        for donor in range(1, 9)
        for rna, dose in (
            ("ets2", "d250ng"), ("ets2", "d500ng"),
            ("rev", "d250ng"), ("rev", "d500ng"),
        )
    ],
    "meki": [
        f"{donor}_{dose}" for donor in ("s1", "s2", "s3")
        for dose in ("100", "500", "ctrl")
    ],
}

CRISPR_SELECTED_INDICES = [
    1, 2, 4, 5, 6, 8, 9, 10, 12, 14, 15, 16, 17, 18, 19,
    21, 22, 23, 24, 26, 28, 29, 31, 33, 34, 36, 38,
]
CRISPR_SELECTED_TARGETS = [
    "g1", "g2", "g1", "g2", "NTC", "g1", "g2", "NTC", "NTC",
    "g1", "g2", "NTC", "g1", "g2", "NTC", "g1", "g2", "NTC",
    "g1", "g2", "NTC", "g1", "g2", "NTC", "g1", "g2", "NTC",
]
CRISPR_SELECTED_DONORS = [
    1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6, 6, 6,
    7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10,
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def digest_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def decimal_text(value: Decimal | None) -> str:
    if value is None:
        return ""
    if value == 0:
        return "0"
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


def parse_decimal(value: str) -> Decimal:
    try:
        number = Decimal(value)
    except (InvalidOperation, TypeError) as exc:
        raise ValueError(f"Non-numeric source value: {value!r}") from exc
    if not number.is_finite():
        raise ValueError(f"Non-finite source value: {value!r}")
    return number


def sign_text(value: str) -> str:
    number = parse_decimal(value)
    if number > 0:
        return "positive"
    if number < 0:
        return "negative"
    return "zero"


def column_number(reference: str) -> int:
    letters = "".join(char for char in reference if char.isalpha())
    result = 0
    for char in letters:
        result = result * 26 + ord(char.upper()) - ord("A") + 1
    return result


def column_name(number: int) -> str:
    result = ""
    while number:
        number, remainder = divmod(number - 1, 26)
        result = chr(ord("A") + remainder) + result
    return result


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def verify_manifest(path: Path) -> dict[str, Any]:
    """Verify every local file listed in a source manifest."""

    manifest = read_json(path)
    source_dir = ROOT / manifest["source_directory"]
    records: list[dict[str, Any]] = []
    for source in manifest.get("sources", []):
        source_path = source_dir / source["file"]
        actual_bytes = source_path.stat().st_size
        actual_sha = digest_file(source_path)
        record = dict(source)
        record.update(
            {
                "actual_bytes": actual_bytes,
                "actual_sha256": actual_sha,
                "bytes_verified": actual_bytes == source["bytes"],
                "sha256_verified": actual_sha == source["sha256"],
                "verified": (
                    actual_bytes == source["bytes"]
                    and actual_sha == source["sha256"]
                ),
            }
        )
        if not record["verified"]:
            raise ValueError(f"Pinned source mismatch: {source_path}")
        records.append(record)

    return {
        "path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": digest_file(path),
        "snapshot_date": manifest.get("snapshot_date", ""),
        "source_directory": manifest.get("source_directory", ""),
        "sources": records,
        "published_zip_md5": manifest.get("published_zip_md5", ""),
        "retrieval_failures": manifest.get("retrieval_failures", []),
    }


def archive_contents(archive_path: Path, relative_members: Iterable[str]) -> dict[str, Any]:
    requested = set(relative_members)
    contents: dict[str, bytes] = {}
    records: list[dict[str, Any]] = []
    with zipfile.ZipFile(archive_path) as archive:
        names = archive.namelist()
        root_candidates = [name for name in names if name.endswith("/")]
        if not root_candidates:
            raise ValueError("The public release did not have a root directory")
        root_prefix = root_candidates[0]
        for relative in sorted(requested):
            full_name = root_prefix + relative
            try:
                info = archive.getinfo(full_name)
            except KeyError as exc:
                raise ValueError(f"Missing archive member: {relative}") from exc
            data = archive.read(full_name)
            contents[relative] = data
            records.append(
                {
                    "member": relative,
                    "archive_member": full_name,
                    "bytes": info.file_size,
                    "sha256": sha256_bytes(data),
                    "compress_type": info.compress_type,
                }
            )

        all_names = set(names)
        raw_count_folders: dict[str, bool] = {}
        for relative in sorted(requested):
            folder = relative.rsplit("/", 1)[0]
            raw_count_folders[folder] = root_prefix + folder + "/raw_counts.txt" in all_names

    return {
        "root_prefix": root_prefix,
        "contents": contents,
        "records": records,
        "raw_count_folders": raw_count_folders,
    }


def read_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        data = archive.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(data)
    values = []
    for si in root.iter():
        if local_name(si.tag) == "si":
            values.append("".join(si.itertext()))
    return values


def parse_xlsx_rows(data: bytes, shared_strings: list[str]) -> tuple[list[dict[str, Any]], str]:
    root = ET.fromstring(data)
    dimension = ""
    for child in root:
        if local_name(child.tag) == "dimension":
            dimension = child.attrib.get("ref", "")
            break
    rows: list[dict[str, Any]] = []
    for row_element in root.iter():
        if local_name(row_element.tag) != "row":
            continue
        row_number = int(row_element.attrib["r"])
        values: dict[int, str] = {}
        for cell in row_element:
            if local_name(cell.tag) != "c":
                continue
            reference = cell.attrib.get("r", "")
            if not reference:
                continue
            kind = cell.attrib.get("t", "")
            value_element = next(
                (child for child in cell if local_name(child.tag) == "v"), None
            )
            if kind == "inlineStr":
                inline = next(
                    (child for child in cell if local_name(child.tag) == "is"), None
                )
                value = "" if inline is None else "".join(inline.itertext())
            elif value_element is None:
                value = ""
            else:
                value = value_element.text or ""
                if kind == "s" and value:
                    value = shared_strings[int(value)]
            values[column_number(reference)] = value
        rows.append({"source_row_number": row_number, "cells": values})
    return rows, dimension


def workbook_sheets(path: Path) -> list[dict[str, Any]]:
    with zipfile.ZipFile(path) as archive:
        shared_strings = read_shared_strings(archive)
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map: dict[str, str] = {}
        for relationship in relationships:
            if local_name(relationship.tag) != "Relationship":
                continue
            target = relationship.attrib["Target"]
            rel_map[relationship.attrib["Id"]] = posixpath.normpath(
                posixpath.join("xl", target)
            )
        result = []
        for sheet in workbook.iter():
            if local_name(sheet.tag) != "sheet":
                continue
            rel_id = sheet.attrib.get(
                "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id",
                sheet.attrib.get("r:id", ""),
            )
            xml_member = rel_map[rel_id]
            rows, dimension = parse_xlsx_rows(archive.read(xml_member), shared_strings)
            result.append(
                {
                    "name": sheet.attrib.get("name", ""),
                    "xml_member": xml_member,
                    "dimension": dimension,
                    "rows": rows,
                }
            )
        return result


def extract_workbook_tables(path: Path) -> dict[str, Any]:
    workbook_sha = digest_file(path)
    workbook_bytes = path.stat().st_size
    sheets = workbook_sheets(path)
    output: dict[str, Any] = {
        "path": str(path.relative_to(ROOT)),
        "bytes": workbook_bytes,
        "sha256": workbook_sha,
        "sheets": {},
    }
    for sheet in sheets:
        rows = sheet["rows"]
        valid_rows = [
            row for row in rows if row["cells"].get(1, "").startswith("ENSG")
        ]
        if not valid_rows:
            continue
        workbook_id = path.name
        is_s1 = sheet["name"] == "Table S1_loss-of-function"
        contrasts = S1_CONTRASTS if is_s1 else S2_CONTRASTS
        row_by_ensembl: dict[str, dict[str, Any]] = {}
        for row in valid_rows:
            ensembl = row["cells"].get(1, "")
            if ensembl in row_by_ensembl:
                raise ValueError(f"Duplicate Ensembl ID in {workbook_id}: {ensembl}")
            row_by_ensembl[ensembl] = row
        symbol_rows: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in valid_rows:
            symbol_rows[row["cells"].get(2, "")].append(row)
        last_data_row = max(row["source_row_number"] for row in valid_rows)
        notes = []
        for row in rows:
            if row["source_row_number"] <= last_data_row:
                continue
            for column, value in sorted(row["cells"].items()):
                if value != "":
                    notes.append(
                        {
                            "workbook_id": workbook_id,
                            "worksheet_name": sheet["name"],
                            "source_row_number": row["source_row_number"],
                            "cell_column": column_name(column),
                            "value_original": value,
                            "workbook_sha256": workbook_sha,
                        }
                    )
        output["sheets"][sheet["name"]] = {
            "workbook_id": workbook_id,
            "worksheet_name": sheet["name"],
            "worksheet_xml_member": sheet["xml_member"],
            "worksheet_dimension": sheet["dimension"],
            "workbook_bytes": workbook_bytes,
            "workbook_sha256": workbook_sha,
            "header_group_values": {
                column_name(column): value
                for column, value in sorted(rows[0]["cells"].items())
            },
            "header_field_values": {
                column_name(column): value
                for column, value in sorted(rows[1]["cells"].items())
            },
            "contrasts": contrasts,
            "valid_rows": valid_rows,
            "row_by_ensembl": row_by_ensembl,
            "symbol_rows": dict(symbol_rows),
            "notes": notes,
            "last_data_row": last_data_row,
            "rows_with_xml": len(rows),
        }
    return output


def parse_master(data: bytes, source_record: dict[str, Any]) -> list[dict[str, Any]]:
    text = data.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if reader.fieldnames != MASTER_COLUMNS:
        raise ValueError(f"Unexpected master list columns: {reader.fieldnames}")
    index_by_column = Counter()
    rows: list[dict[str, Any]] = []
    for csv_row_number, row in enumerate(reader, start=2):
        for set_name in MASTER_COLUMNS:
            member = (row.get(set_name) or "").strip()
            if not member:
                continue
            index_by_column[set_name] += 1
            rows.append(
                {
                    "set_name": set_name,
                    "member_index": index_by_column[set_name],
                    "gene_id_original": member,
                    "gene_id_namespace": "Ensembl_gene_id",
                    "gene_id_version_suffix_present": bool_text(
                        bool(re.search(r"\.\d+$", member))
                    ),
                    "master_source_row_number": csv_row_number,
                    "source_member": source_record["member"],
                    "source_member_bytes": source_record["bytes"],
                    "source_member_sha256": source_record["sha256"],
                }
            )
    return rows


def rank_rows(data: bytes) -> list[dict[str, Any]]:
    reader = csv.reader(io.StringIO(data.decode("utf-8-sig"), newline=""))
    rows = []
    for position, row in enumerate(reader, start=1):
        if len(row) != 2:
            raise ValueError(f"Rank row {position} has {len(row)} fields")
        symbol, statistic = row
        if symbol == "":
            raise ValueError(f"Rank row {position} has an empty gene symbol")
        number = parse_decimal(statistic)
        rows.append(
            {
                "rank_position": position,
                "gene_symbol_original": symbol,
                "gene_symbol_namespace": "HGNC_symbol_as_published",
                "statistic_original": statistic,
                "statistic_decimal": number,
                "statistic_sign": sign_text(statistic),
            }
        )
    symbols = [row["gene_symbol_original"] for row in rows]
    if len(symbols) != len(set(symbols)):
        raise ValueError("Rank file has duplicate gene symbols")
    if any(
        rows[index]["statistic_decimal"] < rows[index + 1]["statistic_decimal"]
        for index in range(len(rows) - 1)
    ):
        raise ValueError("Rank file is not sorted in nonincreasing statistic order")
    return rows


def source_stats_by_symbol(
    workbook_tables: dict[str, Any], workbook_name: str, sheet_name: str, contrast_id: str
) -> dict[str, list[dict[str, Any]]]:
    sheet = workbook_tables["sheets"][sheet_name]
    contrast = next(item for item in sheet["contrasts"] if item["id"] == contrast_id)
    start = contrast["start_col"]
    result: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in sheet["valid_rows"]:
        cells = row["cells"]
        result[cells.get(2, "")].append(
            {
                "source_row_number": row["source_row_number"],
                "ensembl_id_original": cells.get(1, ""),
                "gene_symbol_original": cells.get(2, ""),
                "logfc_original": cells.get(start, ""),
                "ave_expr_original": cells.get(start + 1, ""),
                "t_original": cells.get(start + 2, ""),
                "p_value_original": cells.get(start + 3, ""),
                "adj_p_val_original": cells.get(start + 4, ""),
                "b_original": cells.get(start + 5, ""),
            }
        )
    return dict(result)


def differential_rows(
    workbook_tables: dict[str, Any],
    sheet_name: str,
    source_member: str,
    source_record: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    sheet = workbook_tables["sheets"][sheet_name]
    workbook_id = sheet["workbook_id"]
    differential: list[dict[str, Any]] = []
    source_notes: list[dict[str, Any]] = []
    workbook_summary = {
        "workbook_id": workbook_id,
        "worksheet_name": sheet_name,
        "worksheet_xml_member": sheet["worksheet_xml_member"],
        "worksheet_dimension": sheet["worksheet_dimension"],
        "worksheet_body_rows_after_header": (
            int(re.search(r"\d+$", sheet["worksheet_dimension"].split(":")[-1]).group()) - 2
            if re.search(r"\d+$", sheet["worksheet_dimension"].split(":")[-1])
            else None
        ),
        "workbook_bytes": sheet["workbook_bytes"],
        "workbook_sha256": sheet["workbook_sha256"],
        "rows_with_xml": sheet["rows_with_xml"],
        "valid_ensembl_rows": len(sheet["valid_rows"]),
        "first_valid_source_row": min(row["source_row_number"] for row in sheet["valid_rows"]),
        "last_valid_source_row": sheet["last_data_row"],
        "rows_with_any_cell_after_header": sum(
            1
            for row in sheet["valid_rows"]
        ) + len(sheet["notes"]),
        "literal_none_rows": sum(
            row["cells"].get(2, "") == "None" for row in sheet["valid_rows"]
        ),
        "unique_ensembl_ids": len({row["cells"].get(1, "") for row in sheet["valid_rows"]}),
        "unique_non_none_symbols": len(
            {
                row["cells"].get(2, "")
                for row in sheet["valid_rows"]
                if row["cells"].get(2, "") not in ("", "None")
            }
        ),
        "duplicate_non_none_symbols": {
            symbol: len(rows)
            for symbol, rows in sorted(sheet["symbol_rows"].items())
            if symbol not in ("", "None") and len(rows) > 1
        },
        "notes": sheet["notes"],
        "header_group_values": sheet["header_group_values"],
        "header_field_values": sheet["header_field_values"],
    }
    for note in sheet["notes"]:
        source_notes.append(
            {
                "workbook_id": note["workbook_id"],
                "worksheet_name": note["worksheet_name"],
                "source_row_number": note["source_row_number"],
                "cell_column": note["cell_column"],
                "value_original": note["value_original"],
                "workbook_sha256": note["workbook_sha256"],
            }
        )
    for row in sheet["valid_rows"]:
        cells = row["cells"]
        ensembl = cells.get(1, "")
        symbol = cells.get(2, "")
        for contrast in sheet["contrasts"]:
            start = contrast["start_col"]
            logfc = cells.get(start, "")
            t_value = cells.get(start + 2, "")
            differential.append(
                {
                    "workbook_id": workbook_id,
                    "worksheet_name": sheet_name,
                    "worksheet_xml_member": sheet["worksheet_xml_member"],
                    "workbook_sha256": sheet["workbook_sha256"],
                    "worksheet_source_row_number": row["source_row_number"],
                    "source_member": source_member,
                    "source_member_bytes": source_record["bytes"],
                    "source_member_sha256": source_record["sha256"],
                    "contrast_id": contrast["id"],
                    "contrast_label": contrast["label"],
                    "condition_left_as_published": contrast["left"],
                    "condition_right_as_published": contrast["right"],
                    "ensembl_id_original": ensembl,
                    "ensembl_id_namespace": "Ensembl_gene_id",
                    "gene_symbol_original": symbol,
                    "gene_symbol_is_literal_none": bool_text(symbol == "None"),
                    "logfc_original": logfc,
                    "ave_expr_original": cells.get(start + 1, ""),
                    "t_original": t_value,
                    "p_value_original": cells.get(start + 3, ""),
                    "adj_p_val_original": cells.get(start + 4, ""),
                    "b_original": cells.get(start + 5, ""),
                    "logfc_sign": sign_text(logfc),
                    "t_sign": sign_text(t_value),
                    "logfc_t_sign_agree": bool_text(sign_text(logfc) == sign_text(t_value)),
                    "adj_lt_0_05": bool_text(parse_decimal(cells.get(start + 4, "")) < Decimal("0.05")),
                    "adj_lt_0_10": bool_text(parse_decimal(cells.get(start + 4, "")) < Decimal("0.10")),
                    "source_table_footer_present": bool_text(bool(sheet["notes"])),
                }
            )
    return differential, source_notes, workbook_summary


def master_set_metrics(
    master_rows: list[dict[str, Any]],
    s1_sheet: dict[str, Any],
    differential: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    s1_rows = {
        row["cells"][1]: row for row in s1_sheet["valid_rows"]
    }
    s1_by_contrast: dict[str, dict[str, dict[str, str]]] = {}
    for contrast in S1_CONTRASTS:
        values: dict[str, dict[str, str]] = {}
        start = contrast["start_col"]
        for row in s1_sheet["valid_rows"]:
            cells = row["cells"]
            values[cells[1]] = {
                "gene_symbol_original": cells.get(2, ""),
                "logfc_original": cells.get(start, ""),
                "t_original": cells.get(start + 2, ""),
                "adj_p_val_original": cells.get(start + 4, ""),
                "source_row_number": row["source_row_number"],
            }
        s1_by_contrast[contrast["id"]] = values

    set_to_contrast = {
        "ETS2g1_UP": ("s1_g1", "positive"),
        "ETS2g1_DN": ("s1_g1", "negative"),
        "CHR21_DN": ("s1_chr21", "negative"),
        "CHR21_UP": ("s1_chr21", "positive"),
        "ETS2g2_UP": ("s1_g2", "positive"),
        "ETS2g2_DN": ("s1_g2", "negative"),
    }
    rows_by_set: defaultdict[str, list[str]] = defaultdict(list)
    for row in master_rows:
        rows_by_set[row["set_name"]].append(row["gene_id_original"])

    for row in master_rows:
        stats = s1_by_contrast[set_to_contrast[row["set_name"]][0]].get(
            row["gene_id_original"], {}
        )
        symbol = stats.get("gene_symbol_original", "")
        row.update(
            {
                "gene_symbol_s1_original": symbol,
                "gene_symbol_s1_is_literal_none": bool_text(symbol == "None"),
                "s1_source_row_number": stats.get("source_row_number", ""),
                "s1_source_contrast": set_to_contrast[row["set_name"]][0],
                "source_logfc_original": stats.get("logfc_original", ""),
                "source_t_original": stats.get("t_original", ""),
                "source_adj_p_val_original": stats.get("adj_p_val_original", ""),
                "fdr05_membership_reproduced": bool_text(
                    row["gene_id_original"] in s1_by_contrast[set_to_contrast[row["set_name"]][0]]
                    and parse_decimal(stats["adj_p_val_original"]) < Decimal("0.05")
                    and sign_text(stats["logfc_original"]) == set_to_contrast[row["set_name"]][1]
                ),
                "fdr10_membership_reproduced": bool_text(
                    row["gene_id_original"] in s1_by_contrast[set_to_contrast[row["set_name"]][0]]
                    and parse_decimal(stats["adj_p_val_original"]) < Decimal("0.10")
                    and sign_text(stats["logfc_original"]) == set_to_contrast[row["set_name"]][1]
                ),
                "symbol_to_ensembl_count_s1": sum(
                    item["cells"].get(2, "") == symbol
                    for item in s1_sheet["valid_rows"]
                ),
            }
        )

    expected_sets: dict[str, set[str]] = {}
    fdr10_sets: dict[str, set[str]] = {}
    signs_agree: dict[str, int] = {}
    signs_disagree: dict[str, int] = {}
    for set_name, (contrast_id, desired_sign) in set_to_contrast.items():
        expected_sets[set_name] = {
            ensembl
            for ensembl, stats in s1_by_contrast[contrast_id].items()
            if parse_decimal(stats["adj_p_val_original"]) < Decimal("0.05")
            and sign_text(stats["logfc_original"]) == desired_sign
        }
        fdr10_sets[set_name] = {
            ensembl
            for ensembl, stats in s1_by_contrast[contrast_id].items()
            if parse_decimal(stats["adj_p_val_original"]) < Decimal("0.10")
            and sign_text(stats["logfc_original"]) == desired_sign
        }
        agree = disagree = 0
        for stats in s1_by_contrast[contrast_id].values():
            if sign_text(stats["logfc_original"]) == sign_text(stats["t_original"]):
                agree += 1
            else:
                disagree += 1
        signs_agree[contrast_id] = agree
        signs_disagree[contrast_id] = disagree

    actual_sets = {name: set(rows_by_set[name]) for name in MASTER_COLUMNS}
    reproduction = {}
    for name in MASTER_COLUMNS:
        reproduction[name] = {
            "master_count": len(actual_sets[name]),
            "reproduced_fdr05_count": len(expected_sets[name]),
            "exact_member_set": actual_sets[name] == expected_sets[name],
            "fdr10_directional_count": len(fdr10_sets[name]),
        }
    overlaps = {}
    for index, left in enumerate(MASTER_COLUMNS):
        for right in MASTER_COLUMNS[index + 1 :]:
            overlaps[f"{left}__{right}"] = len(actual_sets[left] & actual_sets[right])
    union = set().union(*actual_sets.values())
    occurrence_counts = Counter(member for members in actual_sets.values() for member in members)
    metrics = {
        "columns": {
            name: {
                "rows_in_output": len(rows_by_set[name]),
                "unique_members": len(actual_sets[name]),
                "within_column_duplicate_count": len(rows_by_set[name]) - len(actual_sets[name]),
            }
            for name in MASTER_COLUMNS
        },
        "union_unique_members": len(union),
        "union_members_in_multiple_sets": sum(count > 1 for count in occurrence_counts.values()),
        "union_duplicate_memberships": sum(max(0, count - 1) for count in occurrence_counts.values()),
        "pairwise_overlap_counts": overlaps,
        "fdr_rule": "strict adjusted P value < 0.05; direction determined by source logFC",
        "fdr05_reproduction": reproduction,
        "fdr10_directional_counts": {
            name: len(fdr10_sets[name]) for name in MASTER_COLUMNS
        },
        "source_logfc_t_sign_agreement": {
            contrast_id: {"agree": signs_agree[contrast_id], "disagree": signs_disagree[contrast_id]}
            for contrast_id in signs_agree
        },
        "all_master_members_have_s1_rows": all(
            row["gene_id_original"] in s1_rows for row in master_rows
        ),
    }
    return master_rows, metrics


def make_design_rows(
    code_hashes: dict[str, dict[str, Any]],
    raw_count_folders: dict[str, bool],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def add(
        experiment_id: str,
        code_rel: str,
        sample_index: int,
        sample_name: str,
        target: str,
        donor: str,
        included: bool,
        contrast_id: str,
        formula: str,
        left: str,
        right: str,
        note: str,
    ) -> None:
        code = code_hashes[code_rel]
        rows.append(
            {
                "experiment_id": experiment_id,
                "code_member": code_rel,
                "code_member_bytes": code["bytes"],
                "code_member_sha256": code["sha256"],
                "sample_index": sample_index,
                "sample_name_literal": sample_name,
                "target_as_coded": target,
                "donor_as_coded": donor,
                "included_in_model": bool_text(included),
                "contrast_id": contrast_id,
                "contrast_formula": formula,
                "condition_left_as_coded": left,
                "condition_right_as_coded": right,
                "design_note": note,
            }
        )

    crispr_code = "RNA-seq/RNAseq_CRISPR/TPP_CRISPR_RNAseq_analysis.r"
    selected = dict(
        zip(CRISPR_SELECTED_INDICES, zip(CRISPR_SELECTED_TARGETS, CRISPR_SELECTED_DONORS))
    )
    for index, name in enumerate(DESIGN_SAMPLE_NAMES["crispr"], start=1):
        target_donor = selected.get(index)
        target = target_donor[0] if target_donor else ""
        donor = str(target_donor[1]) if target_donor else ""
        add(
            "crispr_ets2",
            crispr_code,
            index,
            name,
            target,
            donor,
            target_donor is not None,
            "g1_vs_NTC;g2_vs_NTC" if target_donor else "",
            "g1-NTC;g2-NTC" if target_donor else "",
            "g1;g2" if target_donor else "",
            "NTC;NTC" if target_donor else "",
            (
                "R target factor g2 is assigned to literal g4 sample names; "
                "selected model has 10 donor labels, with complete g1/NTC pairs "
                "for donors 2,3,5,6,7,8,9,10."
                if target_donor
                else "Listed in all.names but excluded by the released ko column index vector."
            ),
        )

    oe_code = "RNA-seq/RNAseq_overexpression/RNAseq_overexpression_analysis.r"
    for index, name in enumerate(DESIGN_SAMPLE_NAMES["oe"], start=1):
        donor, rna, dose = name.split("_")
        dose_label = dose.replace("d", "")
        sample = f"{rna}_{dose}"
        add(
            "oe_ets2",
            oe_code,
            index,
            name,
            sample,
            donor,
            True,
            "oe_250_vs_rev_250;oe_500_vs_rev_500;oe_combined_vs_rev",
            "ets2_250-rev_250;ets2_500-rev_500;ets2-rev",
            "ets2_250;ets2_500;ets2",
            "rev_250;rev_500;rev",
            (
                f"Split-dose target is {sample}; {dose_label} dose. Eight complete "
                "paired donors are represented. The combined model includes dose as covariate."
            ),
        )

    meki_code = "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r"
    for index, name in enumerate(DESIGN_SAMPLE_NAMES["meki"], start=1):
        donor, dose = name.split("_")
        target = {"100": "100nM", "500": "500nM", "ctrl": "ctrl"}[dose]
        add(
            "meki",
            meki_code,
            index,
            name,
            target,
            donor,
            True,
            "mek_100_vs_ctrl;mek_500_vs_ctrl",
            "MEK_100-ctrl;MEK_500-ctrl",
            "MEK_100;MEK_500",
            "ctrl;ctrl",
            (
                "Three complete paired donors are represented. The release code "
                "reads the rank files after defining these contrasts but does not "
                "actively export a matching rank CSV."
            ),
        )

    design_summary = {
        "crispr_ets2": {
            "all_names_count": 38,
            "selected_model_sample_count": 27,
            "donor_labels_in_model": 10,
            "complete_g1_ntc_pair_donors": [2, 3, 5, 6, 7, 8, 9, 10],
            "complete_g2_ntc_pair_donors": [2, 3, 5, 6, 7, 8, 9, 10],
            "raw_count_matrix_present_in_same_archive_folder": raw_count_folders["RNA-seq/RNAseq_CRISPR"],
            "contrasts": ["g1-NTC", "g2-NTC"],
        },
        "oe_ets2": {
            "sample_count": 32,
            "donor_count": 8,
            "complete_paired_donors": list(range(1, 9)),
            "raw_count_matrix_present_in_same_archive_folder": raw_count_folders["RNA-seq/RNAseq_overexpression"],
            "contrasts": ["ets2_250-rev_250", "ets2_500-rev_500", "ets2-rev with dose covariate"],
        },
        "meki": {
            "sample_count": 9,
            "donor_count": 3,
            "complete_paired_donors": ["s1", "s2", "s3"],
            "raw_count_matrix_present_in_same_archive_folder": raw_count_folders["RNA-seq/RNAseq_MEKi"],
            "contrasts": ["MEK_100-ctrl", "MEK_500-ctrl"],
        },
    }
    return rows, design_summary


def compare_rank_to_source(
    rank: list[dict[str, Any]],
    source_by_symbol: dict[str, list[dict[str, Any]]],
    tolerance: Decimal = Decimal("1e-9"),
) -> dict[str, Any]:
    rank_symbols = {row["gene_symbol_original"] for row in rank}
    source_non_none = {
        symbol: candidates
        for symbol, candidates in source_by_symbol.items()
        if symbol not in ("", "None")
    }
    source_unique = {
        symbol: candidates[0]
        for symbol, candidates in source_non_none.items()
        if len(candidates) == 1
    }
    shared = rank_symbols & set(source_non_none)
    unambiguous = rank_symbols & set(source_unique)
    ambiguous = shared - unambiguous
    rank_only = rank_symbols - set(source_non_none)
    source_only_unique = set(source_unique) - rank_symbols
    numeric_matches = 0
    sign_mismatches = 0
    differences: list[Decimal] = []
    candidate_matches = 0
    candidate_sign_mismatches = 0
    nearest_sign_mismatches = 0
    nearest_differences: list[Decimal] = []
    for row in rank:
        symbol = row["gene_symbol_original"]
        candidates = source_non_none.get(symbol, [])
        rank_value = row["statistic_decimal"]
        for candidate in candidates:
            source_value = parse_decimal(candidate["t_original"])
            difference = abs(rank_value - source_value)
            if difference <= tolerance:
                candidate_matches += 1
            if sign_text(row["statistic_original"]) != sign_text(candidate["t_original"]):
                candidate_sign_mismatches += 1
        if symbol in source_unique:
            candidate = source_unique[symbol]
            source_value = parse_decimal(candidate["t_original"])
            difference = abs(rank_value - source_value)
            differences.append(difference)
            if difference <= tolerance:
                numeric_matches += 1
            if sign_text(row["statistic_original"]) != sign_text(candidate["t_original"]):
                sign_mismatches += 1
        if candidates:
            nearest = min(
                candidates,
                key=lambda candidate: (
                    abs(rank_value - parse_decimal(candidate["t_original"])),
                    candidate["source_row_number"],
                ),
            )
            nearest_difference = abs(rank_value - parse_decimal(nearest["t_original"]))
            nearest_differences.append(nearest_difference)
            if sign_text(row["statistic_original"]) != sign_text(nearest["t_original"]):
                nearest_sign_mismatches += 1

    return {
        "rank_symbol_count": len(rank_symbols),
        "source_non_none_symbol_count": len(source_non_none),
        "shared_symbol_count": len(shared),
        "unambiguous_one_to_one_symbol_count": len(unambiguous),
        "ambiguous_shared_symbol_count": len(ambiguous),
        "ambiguous_shared_symbols": sorted(ambiguous),
        "rank_only_symbol_count": len(rank_only),
        "source_only_unique_symbol_count": len(source_only_unique),
        "primary_comparison_rule": (
            "Compare only rank symbols with exactly one non-literal-None source row; "
            "ambiguous symbols and rank-only symbols are unavailable for primary identity."
        ),
        "numeric_match_tolerance": decimal_text(tolerance),
        "primary_unambiguous_match_count": len(unambiguous),
        "primary_numeric_match_within_tolerance_count": numeric_matches,
        "primary_numeric_match_within_tolerance_fraction": (
            numeric_matches / len(unambiguous) if unambiguous else None
        ),
        "primary_sign_mismatch_count": sign_mismatches,
        "primary_max_abs_difference": decimal_text(max(differences) if differences else None),
        "primary_mean_abs_difference": decimal_text(
            sum(differences, Decimal(0)) / len(differences) if differences else None
        ),
        "candidate_all_shared_row_count": sum(
            len(source_non_none[symbol]) for symbol in shared
        ),
        "candidate_all_numeric_match_count": candidate_matches,
        "candidate_all_sign_mismatch_count": candidate_sign_mismatches,
        "nearest_candidate_sign_mismatch_count_diagnostic": nearest_sign_mismatches,
        "nearest_candidate_max_abs_difference_diagnostic": decimal_text(
            max(nearest_differences) if nearest_differences else None
        ),
    }


def rank_output_rows(
    spec: dict[str, Any],
    rank_data: list[dict[str, Any]],
    rank_record: dict[str, Any],
    code_record: dict[str, Any],
    source_by_symbol: dict[str, list[dict[str, Any]]] | None,
    source_record: dict[str, Any] | None,
    source_sheet: str,
    source_contrast: str,
) -> list[dict[str, Any]]:
    source_non_none = {
        symbol: candidates
        for symbol, candidates in (source_by_symbol or {}).items()
        if symbol not in ("", "None")
    }
    source_unique = {
        symbol: candidates[0]
        for symbol, candidates in source_non_none.items()
        if len(candidates) == 1
    }
    output = []
    for row in rank_data:
        symbol = row["gene_symbol_original"]
        candidates = source_non_none.get(symbol, [])
        unique = source_unique.get(symbol)
        primary_difference = ""
        primary_within = ""
        primary_sign = ""
        primary_row = ""
        primary_ensembl = ""
        primary_t = ""
        primary_logfc = ""
        primary_adj = ""
        if unique is not None:
            source_value = parse_decimal(unique["t_original"])
            difference = abs(row["statistic_decimal"] - source_value)
            primary_difference = decimal_text(difference)
            primary_within = bool_text(difference <= Decimal("1e-9"))
            primary_sign = bool_text(
                sign_text(row["statistic_original"]) == sign_text(unique["t_original"])
            )
            primary_row = unique["source_row_number"]
            primary_ensembl = unique["ensembl_id_original"]
            primary_t = unique["t_original"]
            primary_logfc = unique["logfc_original"]
            primary_adj = unique["adj_p_val_original"]
        candidate_rows = ";".join(str(item["source_row_number"]) for item in candidates)
        candidate_ensembl = ";".join(item["ensembl_id_original"] for item in candidates)
        candidate_t = ";".join(item["t_original"] for item in candidates)
        candidate_logfc = ";".join(item["logfc_original"] for item in candidates)
        candidate_sign_mismatch = sum(
            sign_text(row["statistic_original"]) != sign_text(item["t_original"])
            for item in candidates
        )
        nearest_row = ""
        nearest_difference = ""
        nearest_sign = ""
        if candidates:
            nearest = min(
                candidates,
                key=lambda item: (
                    abs(row["statistic_decimal"] - parse_decimal(item["t_original"])),
                    item["source_row_number"],
                ),
            )
            nearest_row = nearest["source_row_number"]
            nearest_difference = decimal_text(
                abs(row["statistic_decimal"] - parse_decimal(nearest["t_original"]))
            )
            nearest_sign = bool_text(
                sign_text(row["statistic_original"]) == sign_text(nearest["t_original"])
            )
        output.append(
            {
                "rank_dataset_id": spec["id"],
                "rank_file": spec["member"],
                "rank_member_bytes": rank_record["bytes"],
                "rank_member_sha256": rank_record["sha256"],
                "rank_position": row["rank_position"],
                "gene_symbol_original": symbol,
                "gene_symbol_namespace": row["gene_symbol_namespace"],
                "statistic_original": row["statistic_original"],
                "statistic_decimal": decimal_text(row["statistic_decimal"]),
                "statistic_sign": row["statistic_sign"],
                "rank_order_is_descending": "true",
                "code_member": spec["code_member"],
                "code_member_bytes": code_record["bytes"],
                "code_member_sha256": code_record["sha256"],
                "code_contrast": spec["code_contrast"],
                "code_effect_direction": spec["code_effect_direction"],
                "code_rank_linkage_status": spec["source_comparison_status"],
                "source_workbook_id": spec["source_workbook"],
                "source_worksheet_name": source_sheet,
                "source_contrast": source_contrast,
                "source_workbook_sha256": (source_record or {}).get("sha256", ""),
                "source_match_status": (
                    "one_to_one" if unique is not None
                    else "ambiguous_source_symbol" if candidates
                    else "rank_only"
                ),
                "source_symbol_match_count": len(candidates),
                "source_ensembl_ids_unambiguous": primary_ensembl,
                "source_gene_symbol_original_unambiguous": (
                    unique["gene_symbol_original"] if unique else ""
                ),
                "source_row_number_unambiguous": primary_row,
                "source_t_original_unambiguous": primary_t,
                "source_logfc_original_unambiguous": primary_logfc,
                "source_adj_p_val_original_unambiguous": primary_adj,
                "source_abs_difference_decimal_unambiguous": primary_difference,
                "source_numeric_match_within_1e-9_unambiguous": primary_within,
                "source_sign_match_unambiguous": primary_sign,
                "source_candidate_row_numbers_diagnostic": candidate_rows,
                "source_candidate_ensembl_ids_diagnostic": candidate_ensembl,
                "source_candidate_t_original_diagnostic": candidate_t,
                "source_candidate_logfc_original_diagnostic": candidate_logfc,
                "source_candidate_sign_mismatch_count_diagnostic": candidate_sign_mismatch,
                "source_nearest_row_number_diagnostic": nearest_row,
                "source_nearest_abs_difference_decimal_diagnostic": nearest_difference,
                "source_nearest_sign_match_diagnostic": nearest_sign,
                "source_comparison_status": (
                    "publisher_differential_table_not_pinned"
                    if source_by_symbol is None
                    else spec["source_comparison_status"]
                ),
            }
        )
    return output


def article_evidence(data: bytes) -> dict[str, Any]:
    root = ET.fromstring(data)
    paragraphs = []
    for element in root.iter():
        if local_name(element.tag) == "p":
            paragraphs.append(" ".join("".join(element.itertext()).split()))

    queries = [
        ("ranked_by_t_statistic", "ranked by t-statistic"),
        ("supplementary_differential_tables", "Differential expression results are shown in Supplementary Tables 1 and 2"),
        ("downregulated_g1_program", "downregulated after ETS2 disruption"),
        ("modular_g1_program", "Modular expression of ETS2-regulated genes (downregulated after ETS2 editing, gRNA1)"),
        ("crispr_n", "RNA-seq analysis of differentially expressed genes in ETS2-edited"),
        ("ega_accessions", "EGAD00001011338"),
    ]
    evidence = []
    for evidence_id, phrase in queries:
        match_index = next(
            (index for index, paragraph in enumerate(paragraphs, start=1)
             if phrase.lower() in paragraph.lower()),
            None,
        )
        if match_index is None:
            evidence.append({"id": evidence_id, "phrase": phrase, "paragraph_number": None})
            continue
        excerpt = paragraphs[match_index - 1]
        if len(excerpt) > 900:
            excerpt = excerpt[:897] + "..."
        evidence.append(
            {
                "id": evidence_id,
                "phrase": phrase,
                "paragraph_number": match_index,
                "excerpt": excerpt,
            }
        )
    return {"paragraph_count": len(paragraphs), "evidence": evidence}


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def build_audit(write_outputs: bool = True) -> dict[str, Any]:
    base = verify_manifest(BASE_MANIFEST)
    supplements = verify_manifest(PROGRAM_MANIFEST)
    base_source_records = {
        item["file"]: item for item in base["sources"]
    }
    supplement_records = {
        item["file"]: item for item in supplements["sources"]
    }
    archive_path = RAW / "ets2-public-release.zip"
    archive = archive_contents(archive_path, ARCHIVE_RELEVANT_MEMBERS)
    archive_md5_actual = digest_file(archive_path, "md5")
    if archive_md5_actual != base["published_zip_md5"]:
        raise ValueError("Publisher ZIP MD5 mismatch")
    member_records = {item["member"]: item for item in archive["records"]}

    master = parse_master(archive["contents"][MASTER_MEMBER], member_records[MASTER_MEMBER])
    workbook_tables = extract_workbook_tables(RAW / "ets2-supplementary-tables.xlsx")
    s1_sheet = workbook_tables["sheets"]["Table S1_loss-of-function"]
    s2_sheet = workbook_tables["sheets"]["Table S2_gain-of-function"]
    differential_s1, notes_s1, s1_summary = differential_rows(
        workbook_tables,
        "Table S1_loss-of-function",
        "ets2-supplementary-tables.xlsx",
        supplement_records["ets2-supplementary-tables.xlsx"],
    )
    differential_s2, notes_s2, s2_summary = differential_rows(
        workbook_tables,
        "Table S2_gain-of-function",
        "ets2-supplementary-tables.xlsx",
        supplement_records["ets2-supplementary-tables.xlsx"],
    )
    differential = differential_s1 + differential_s2
    source_notes = notes_s1 + notes_s2
    master, master_metrics = master_set_metrics(master, s1_sheet, differential)

    source_by_symbol_by_contrast = {
        "s1_g1": source_stats_by_symbol(
            workbook_tables, "ets2-supplementary-tables.xlsx", "Table S1_loss-of-function", "s1_g1"
        ),
        "s1_g2": source_stats_by_symbol(
            workbook_tables, "ets2-supplementary-tables.xlsx", "Table S1_loss-of-function", "s1_g2"
        ),
        "s1_chr21": source_stats_by_symbol(
            workbook_tables, "ets2-supplementary-tables.xlsx", "Table S1_loss-of-function", "s1_chr21"
        ),
        "s2_oe_250": source_stats_by_symbol(
            workbook_tables, "ets2-supplementary-tables.xlsx", "Table S2_gain-of-function", "s2_oe_250"
        ),
        "s2_oe_500": source_stats_by_symbol(
            workbook_tables, "ets2-supplementary-tables.xlsx", "Table S2_gain-of-function", "s2_oe_500"
        ),
    }

    rank_rows_all: list[dict[str, Any]] = []
    rank_summaries = []
    code_hashes: dict[str, dict[str, Any]] = {}
    for spec in RANK_SPECS:
        rank_record = member_records[spec["member"]]
        code_record = member_records[spec["code_member"]]
        code_hashes[spec["code_member"]] = code_record
        parsed_rank = rank_rows(archive["contents"][spec["member"]])
        source_by_symbol = source_by_symbol_by_contrast.get(spec["source_contrast"])
        source_record = supplement_records.get("ets2-supplementary-tables.xlsx") if source_by_symbol else None
        source_comparison = (
            compare_rank_to_source(parsed_rank, source_by_symbol)
            if source_by_symbol is not None
            else {
                "primary_comparison_rule": "No pinned publisher differential table for this experiment.",
                "numeric_match_tolerance": "1e-9",
                "primary_unambiguous_match_count": None,
                "primary_numeric_match_within_tolerance_count": None,
                "primary_sign_mismatch_count": None,
            }
        )
        output_rows = rank_output_rows(
            spec,
            parsed_rank,
            rank_record,
            code_record,
            source_by_symbol,
            source_record,
            spec["source_sheet"],
            spec["source_contrast"],
        )
        rank_rows_all.extend(output_rows)
        values = [row["statistic_decimal"] for row in parsed_rank]
        ets2_rows = [row for row in parsed_rank if row["gene_symbol_original"] == "ETS2"]
        rank_summaries.append(
            {
                "rank_dataset_id": spec["id"],
                "rank_file": spec["member"],
                "rank_member_bytes": rank_record["bytes"],
                "rank_member_sha256": rank_record["sha256"],
                "code_member": spec["code_member"],
                "code_member_sha256": code_record["sha256"],
                "code_contrast": spec["code_contrast"],
                "code_effect_direction": spec["code_effect_direction"],
                "code_rank_linkage": spec["code_rank_linkage"],
                "rank_row_count": len(parsed_rank),
                "unique_gene_symbol_count": len({row["gene_symbol_original"] for row in parsed_rank}),
                "finite_statistic_count": len(values),
                "sorted_nonincreasing": all(
                    values[index] >= values[index + 1] for index in range(len(values) - 1)
                ),
                "min_statistic": decimal_text(min(values)),
                "max_statistic": decimal_text(max(values)),
                "ets2_rank_position": ets2_rows[0]["rank_position"] if ets2_rows else None,
                "ets2_statistic_original": ets2_rows[0]["statistic_original"] if ets2_rows else None,
                "source_comparison_status": spec["source_comparison_status"],
                "source_comparison": source_comparison,
            }
        )

    code_hashes = {
        relative: member_records[relative]
        for relative in {
            spec["code_member"] for spec in RANK_SPECS
        }
    }
    design_rows, design_summary = make_design_rows(
        code_hashes,
        archive["raw_count_folders"],
    )

    source_output_paths = {
        "gene_sets": DERIVED / "ets2-source-gene-sets.csv",
        "differential": DERIVED / "ets2-source-differential.csv.gz",
        "ranks": DERIVED / "ets2-source-ranks.csv.gz",
        "design": DERIVED / "ets2-source-design.csv",
        "notes": DERIVED / "ets2-source-notes.csv",
    }
    gene_set_fields = [
        "set_name", "member_index", "gene_id_original", "gene_id_namespace",
        "gene_id_version_suffix_present", "master_source_row_number", "source_member",
        "source_member_bytes", "source_member_sha256", "gene_symbol_s1_original",
        "gene_symbol_s1_is_literal_none", "s1_source_row_number", "s1_source_contrast",
        "source_logfc_original", "source_t_original", "source_adj_p_val_original",
        "fdr05_membership_reproduced", "fdr10_membership_reproduced",
        "symbol_to_ensembl_count_s1",
    ]
    differential_fields = [
        "workbook_id", "worksheet_name", "worksheet_xml_member", "workbook_sha256",
        "worksheet_source_row_number", "source_member", "source_member_bytes",
        "source_member_sha256", "contrast_id", "contrast_label",
        "condition_left_as_published", "condition_right_as_published", "ensembl_id_original",
        "ensembl_id_namespace", "gene_symbol_original", "gene_symbol_is_literal_none",
        "logfc_original", "ave_expr_original", "t_original", "p_value_original",
        "adj_p_val_original", "b_original", "logfc_sign", "t_sign",
        "logfc_t_sign_agree", "adj_lt_0_05", "adj_lt_0_10", "source_table_footer_present",
    ]
    rank_fields = list(rank_rows_all[0]) if rank_rows_all else []
    design_fields = list(design_rows[0]) if design_rows else []
    note_fields = [
        "workbook_id", "worksheet_name", "source_row_number", "cell_column",
        "value_original", "workbook_sha256",
    ]
    if write_outputs:
        write_csv(source_output_paths["gene_sets"], gene_set_fields, master)
        write_gzip_csv(source_output_paths["differential"], differential_fields, differential)
        write_gzip_csv(source_output_paths["ranks"], rank_fields, rank_rows_all)
        write_csv(source_output_paths["design"], design_fields, design_rows)
        write_csv(source_output_paths["notes"], note_fields, source_notes)

    article_path = RAW / "ets2-primary-article.xml"
    article_record = next(
        source for source in base["sources"] if source["file"] == "ets2-primary-article.xml"
    )
    report: dict[str, Any] = {
        "analysis_version": "1.0",
        "snapshot_date": base["snapshot_date"],
        "scope": "Public ETS2 source and semantics audit; no atlas scores or donor count-matrix reconstruction.",
        "source_manifests": [base, supplements],
        "publisher_zip_md5": {
            "expected": base["published_zip_md5"],
            "actual": archive_md5_actual,
            "verified": archive_md5_actual == base["published_zip_md5"],
        },
        "archive": {
            "path": str(archive_path.relative_to(ROOT)),
            "root_prefix": archive["root_prefix"],
            "relevant_members": archive["records"],
            "raw_count_folders": archive["raw_count_folders"],
        },
        "article": {
            "file": article_record["file"],
            "url": article_record["url"],
            "bytes": article_record["actual_bytes"],
            "sha256": article_record["actual_sha256"],
            **article_evidence(article_path.read_bytes()),
        },
        "workbooks": {
            "s1": s1_summary,
            "s2": s2_summary,
            "notes_preserved_in": str(source_output_paths["notes"].relative_to(ROOT)),
        },
        "gene_sets": master_metrics,
        "rank_datasets": rank_summaries,
        "designs": design_summary,
        "output_schemas": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "columns": fields,
            }
            for key, path, fields in (
                ("gene_sets", source_output_paths["gene_sets"], gene_set_fields),
                ("differential", source_output_paths["differential"], differential_fields),
                ("ranks", source_output_paths["ranks"], rank_fields),
                ("design", source_output_paths["design"], design_fields),
                ("notes", source_output_paths["notes"], note_fields),
            )
        },
        "caveats": [
            "The master columns reproduce strict adjusted P value < 0.05 directional memberships from S1; FDR < 0.10 counts are reported separately and are not substituted.",
            "Master members are Ensembl gene IDs without version suffixes; rank files use published gene symbols. Literal None values and duplicate symbol rows remain explicit.",
            "KO gRNA1 ranks numerically match the one-to-one S1 gRNA1 t statistics within Decimal tolerance 1e-9; the raw xlsx and CSV strings differ in some floating-point serialization.",
            "OE rank values are not numerically within 1e-9 of one-to-one S2 t statistics, despite high rank-level concordance. Sign mismatches are reported without choosing a preferred source row by nearest numeric value.",
            "The MEKi release includes code-defined MEK minus control contrasts and rank files but no pinned publisher differential table or active matching export; rank orientation remains unresolved.",
            "Relevant CRISPR, overexpression and MEKi public folders do not contain the raw_counts.txt files referenced by their R scripts. Biopsy raw counts are not a substitute.",
            "R designs describe donor/sample structure only. Saved rank vectors cannot reproduce donor fits, count matrices or held-out-donor validation.",
            "No biological orientation was inferred from expected inflammatory, oxidative-stress or cell-death responses.",
        ],
    }
    if write_outputs:
        report["derived_output_sha256"] = {
            key: digest_file(path)
            for key, path in source_output_paths.items()
        }
        report_path = REPORTS / "ets2-source-audit.json"
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        report["report_path"] = str(report_path.relative_to(ROOT))
    return report


def main() -> None:
    report = build_audit(write_outputs=True)
    print(
        "Verified ETS2 sources; "
        f"master union={report['gene_sets']['union_unique_members']}, "
        f"rank datasets={len(report['rank_datasets'])}."
    )


if __name__ == "__main__":
    main()
