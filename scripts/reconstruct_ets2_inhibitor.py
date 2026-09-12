"""Bounded forensic reconstruction of the released MEK-inhibitor ranks.

The public release does not contain the MEK-inhibitor ``raw_counts.txt`` or an
active export that connects the saved rank CSVs to the code's topTable objects.
This module therefore reports code-supported steps and deterministic checks of
the released vectors.  It never treats the inferred pipeline as recovered
production provenance and never assigns a treatment sign to the saved ranks.
"""

from __future__ import annotations

import csv
import io
import json
import math
import posixpath
import statistics
import zipfile
from pathlib import Path
from typing import Any, Iterable
from xml.etree import ElementTree as ET

from audit_ets2_sources import (
    ARCHIVE_RELEVANT_MEMBERS,
    ROOT,
    RAW,
    RANK_SPECS,
    archive_contents,
    bool_text,
    column_name,
    column_number,
    decimal_text,
    digest_file,
    local_name,
    parse_decimal,
    parse_xlsx_rows,
    rank_rows,
    read_shared_strings,
    verify_manifest,
    write_csv,
)


DERIVED = ROOT / "data/derived"
REPORTS = ROOT / "reports"

ARCHIVE_EXTRA_MEMBERS = {
    "RNA-seq/RNAseq_MEKi/macrophage_genesets_fgsea.csv",
    "RNA-seq/RNAseq_MEKi/fgsea_PD_0.5_mac.pathways.csv",
}

FIG5_SHEET_NAME = "5j"
FIG5_PATHWAY_LABEL_TO_ID = {
    "Macrophage activation": "GOBP_MACROPHAGE_ACTIVATION",
    "Myeloid leukocyte activation": "GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
    "Myeloid leukocyte migration": "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
    "Phagocytosis": "GOBP_PHAGOCYTOSIS",
    "IL-1 production": "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_1_PRODUCTION",
    "IL-6 production": "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_6_PRODUCTION",
    "IL-8 production": "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_8_PRODUCTION",
    "ROS production": "GOBP_POSITIVE_REGULATION_OF_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",
    "TNFSF cytokine production": "GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION",
}

RECONSTRUCTION_STEPS = [
    {
        "step_id": "input_count_matrix",
        "phase": "input",
        "confidence": "high_for_intended_input",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": 'counts=data.matrix(read.table(file="raw_counts.txt"',
        "inference": "The intended input is a nine-column raw_counts.txt file in the MEKi analysis folder.",
        "limitation": "The public MEKi folder has no raw_counts.txt; no count values or donor fits are recovered.",
    },
    {
        "step_id": "sample_names",
        "phase": "sample_annotation",
        "confidence": "high",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": 'names=paste0(rep(c("s1","s2","s3"),each=3),c("_100","_500","_ctrl"))',
        "inference": "The intended nine samples are three donors with 100 nM, 500 nM and control labels.",
        "limitation": "The code labels are preserved as source metadata; biological donor identity is unavailable in the public folder.",
    },
    {
        "step_id": "donor_and_drug_annotation",
        "phase": "sample_annotation",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": 'dge$samples$donor=rep(c("s1","s2","s3"),each=3)',
        "inference": "The intended DGEList metadata assigns donor s1/s2/s3 and drug levels 100nM/500nM/ctrl in three-sample blocks.",
        "limitation": "This identifies the model labels, not independent donor provenance or sample QC in the inaccessible matrix.",
    },
    {
        "step_id": "expression_filter",
        "phase": "filtering",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": "eset.expressed=rowSums(cpm(eset)>0.5)>=9",
        "inference": "Genes with CPM greater than 0.5 in all nine samples were retained before DGEList normalization.",
        "limitation": "The number of genes before/after this filter cannot be audited without raw counts; topTable requests 12,801 rows.",
    },
    {
        "step_id": "library_normalization",
        "phase": "normalization",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": "dge=calcNormFactors(dge)",
        "inference": "edgeR calcNormFactors was applied to the filtered DGEList.",
        "limitation": "Normalization factors and their version-specific implementation cannot be recomputed without counts and the R package environment.",
    },
    {
        "step_id": "donor_aware_design",
        "phase": "design",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": "design=model.matrix(~0+factor(dge$samples$drug)+factor(dge$samples$donor))",
        "inference": "The intended limma design has separate drug levels and donor covariates without an intercept.",
        "limitation": "This is a code design, not evidence that the unavailable public matrix had exactly these columns or QC after filtering.",
    },
    {
        "step_id": "contrasts",
        "phase": "contrast",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": "contrasts=makeContrasts(MEK_100nM=MEK_100-ctrl,MEK_500nM=MEK_500-ctrl",
        "inference": "The code intends 100 nM minus control and 500 nM minus control contrasts.",
        "limitation": "The saved rank files are read after this step, but the release does not expose an active export linking each file to its contrast.",
    },
    {
        "step_id": "voom",
        "phase": "model_fit",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": "voom=voom(dge,design,plot=T)",
        "inference": "The filtered DGEList and donor-aware design are passed to voom.",
        "limitation": "Exact mean-variance trend and weights cannot be recovered without counts and the R package environment.",
    },
    {
        "step_id": "lm_fit",
        "phase": "model_fit",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": "fit=lmFit(voom,design)",
        "inference": "The voom-transformed expression is fitted with limma lmFit.",
        "limitation": "Coefficient and residual values are unavailable without the count matrix.",
    },
    {
        "step_id": "contrasts_fit",
        "phase": "model_fit",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": "fit2=contrasts.fit(fit,contrasts)",
        "inference": "The two drug-control contrasts are applied to the fitted model.",
        "limitation": "The resulting coefficients and standard errors are not in the public release.",
    },
    {
        "step_id": "e_bayes",
        "phase": "model_fit",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": "fit3=eBayes(fit2)",
        "inference": "limma empirical-Bayes moderation is applied before topTable extraction.",
        "limitation": "A moderated t statistic is coefficient divided by a moderated standard error; it cannot recover counts, logFC or variance components by itself.",
    },
    {
        "step_id": "top_table_statistic",
        "phase": "statistic_extraction",
        "confidence": "high_for_intended_code",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": 'MEK_100nM.res=topTable(fit3,number=12801,p.value=1,lfc=0,adjust.method="BH",coef=1)',
        "inference": "The intended output is a BH-adjusted limma topTable with lfc=0, p.value=1 and 12,801 requested rows per coefficient.",
        "limitation": "The public rank vectors contain 12,095 rows; absent counts and code that generated the rank files prevent identifying the 706-row difference.",
    },
    {
        "step_id": "symbol_mapping",
        "phase": "identifier_mapping",
        "confidence": "low_to_medium",
        "evidence_member": "RNA-seq/RNAseq_CRISPR/TPP_CRISPR_RNAseq_analysis.r",
        "evidence_pattern": 'G_list = getBM(filters = "ensembl_gene_id"',
        "inference": "Analogue CRISPR code maps Ensembl row names to HGNC symbols through biomaRt; the MEKi code itself has no mapping call and its rank filename says SYMBOL.",
        "limitation": "Mapping release, one-to-many handling, empty-symbol filtering and any collapse rule for the MEKi rank export are not recoverable.",
    },
    {
        "step_id": "symbol_collapse",
        "phase": "identifier_mapping",
        "confidence": "low",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": 'names(pd100_num)=rownames(pd100)',
        "inference": "The MEKi code assigns the symbols already present in the rank CSV as vector names; it does not show how Ensembl rows were mapped or collapsed before serialization.",
        "limitation": "One-to-many symbol mappings, duplicate handling and whether a symbol was selected by maximum statistic, first row or another rule are non-identifiable from the rank vector.",
    },
    {
        "step_id": "rank_serialization",
        "phase": "rank_consumption",
        "confidence": "medium",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": 'pd100=read.csv(file="PD_0.1_MEK.resSYMBOL.rank.csv"',
        "inference": "The released rank CSV is consumed as a headerless two-column symbol/statistic vector, then named by row names for fGSEA.",
        "limitation": "The generating write step is absent and the active export linkage cannot be certified; sorted order, tie handling, precision and symbol collapse are observed properties, not recovered production rules.",
    },
    {
        "step_id": "fgsea_parameters",
        "phase": "pathway_check",
        "confidence": "medium_to_high",
        "evidence_member": "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
        "evidence_pattern": "fgsea(pathways=lists, stats=pd500_num, eps=0.0, minSize=15, maxSize=500)",
        "inference": "The code specifies fgsea with eps=0, minSize=15 and maxSize=500; the default gseaParam=1 is compatible with the signed ES check.",
        "limitation": "This is an ES compatibility test, not recovery of fgsea's permutation/null calculation or the exact package version.",
    },
]


def _source_fig5_sheet(path: Path, sheet_name: str) -> dict[str, Any]:
    """Read only one sheet from the large source-data workbook."""

    with zipfile.ZipFile(path) as archive:
        shared = read_shared_strings(archive)
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map: dict[str, str] = {}
        for relation in relationships:
            if local_name(relation.tag) != "Relationship":
                continue
            rel_map[relation.attrib["Id"]] = posixpath.normpath(
                posixpath.join("xl", relation.attrib["Target"])
            )
        rel_id = None
        for sheet in workbook.iter():
            if local_name(sheet.tag) == "sheet" and sheet.attrib.get("name") == sheet_name:
                rel_id = sheet.attrib.get(
                    "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id",
                    sheet.attrib.get("r:id", ""),
                )
                break
        if rel_id is None:
            raise ValueError(f"Missing source-data sheet {sheet_name}")
        xml_member = rel_map[rel_id]
        rows, dimension = parse_xlsx_rows(archive.read(xml_member), shared)
        return {
            "sheet_name": sheet_name,
            "xml_member": xml_member,
            "dimension": dimension,
            "rows": rows,
        }


def _rank_spec_by_id(dataset_id: str) -> dict[str, Any]:
    return next(spec for spec in RANK_SPECS if spec["id"] == dataset_id)


def _archive_records_and_contents() -> tuple[dict[str, bytes], dict[str, dict[str, Any]], dict[str, bool], str]:
    selected = set(ARCHIVE_RELEVANT_MEMBERS) | ARCHIVE_EXTRA_MEMBERS
    result = archive_contents(RAW / "ets2-public-release.zip", selected)
    return result["contents"], {item["member"]: item for item in result["records"]}, result["raw_count_folders"], result["root_prefix"]


def _read_gene_sets(data: bytes) -> tuple[dict[str, list[str]], dict[str, Any]]:
    rows = list(csv.reader(io.StringIO(data.decode("utf-8-sig"), newline="")))
    if not rows:
        raise ValueError("The MEKi gene-set file is empty")
    headers = rows[0]
    pathways: dict[str, list[str]] = {}
    duplicates: dict[str, int] = {}
    for index, header in enumerate(headers):
        members = [row[index] for row in rows[1:] if len(row) > index and row[index]]
        pathways[header] = members
        duplicates[header] = len(members) - len(set(members))
    return pathways, {
        "header_count": len(headers),
        "pathway_ids": headers,
        "member_counts": {header: len(pathways[header]) for header in headers},
        "duplicate_member_counts": duplicates,
    }


def _read_fgsea_report(data: bytes) -> dict[str, dict[str, str]]:
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""))
    result = {}
    for row_number, row in enumerate(reader, start=2):
        pathway = row.get("pathway", "")
        if not pathway:
            continue
        row["source_row_number"] = row_number
        result[pathway] = row
    return result


def _fig5_rows(sheet: dict[str, Any], workbook_sha: str) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    pathway_rows: dict[str, dict[str, Any]] = {}
    notes: list[dict[str, Any]] = []
    for row in sheet["rows"]:
        cells = row["cells"]
        pathway = cells.get(1, "")
        if pathway in FIG5_PATHWAY_LABEL_TO_ID:
            pathway_rows[pathway] = {
                "source_data_row_number": row["source_row_number"],
                "pathway_name_reported_original": pathway,
                "pathway_id": FIG5_PATHWAY_LABEL_TO_ID[pathway],
                "source_data_fig5_nes_original": cells.get(2, ""),
                "source_data_fig5_p_minus_log10_original": cells.get(3, ""),
            }
        elif row["source_row_number"] > 1:
            for column, value in sorted(cells.items()):
                if value:
                    notes.append(
                        {
                            "source_data_sheet": FIG5_SHEET_NAME,
                            "source_data_row_number": row["source_row_number"],
                            "cell_column": column_name(column),
                            "value_original": value,
                            "workbook_sha256": workbook_sha,
                        }
                    )
    return pathway_rows, notes


def _signed_weighted_es(rank: list[dict[str, Any]], members: Iterable[str], gsea_param: float = 1.0) -> tuple[float | None, int]:
    """Calculate the signed weighted enrichment score used by fgsea's ES.

    This follows the running-sum definition for a ranked vector and does not
    calculate a permutation null or NES.  A matching ES and NES sign is only a
    compatibility check against the released result.
    """

    member_set = set(members)
    rank_symbols = {row["gene_symbol_original"] for row in rank}
    hits = member_set & rank_symbols
    hit_count = len(hits)
    n = len(rank)
    if not hits or hit_count == n:
        return None, hit_count
    denominator = sum(
        abs(float(row["statistic_decimal"])) ** gsea_param
        for row in rank
        if row["gene_symbol_original"] in hits
    )
    if denominator == 0:
        return None, hit_count
    miss_step = 1.0 / (n - hit_count)
    running = 0.0
    maximum = -math.inf
    minimum = math.inf
    for row in rank:
        if row["gene_symbol_original"] in hits:
            running += abs(float(row["statistic_decimal"])) ** gsea_param / denominator
        else:
            running -= miss_step
        maximum = max(maximum, running)
        minimum = min(minimum, running)
    return (maximum if abs(maximum) >= abs(minimum) else minimum), hit_count


def _pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) < 2:
        return None
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    left_ss = sum((a - left_mean) ** 2 for a in left)
    right_ss = sum((b - right_mean) ** 2 for b in right)
    if left_ss == 0 or right_ss == 0:
        return None
    return numerator / math.sqrt(left_ss * right_ss)


def _rank_summary_rows(
    rank_data: dict[str, list[dict[str, Any]]],
    rank_records: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    rows = []
    for dataset_id in ("meki_100", "meki_500"):
        spec = _rank_spec_by_id(dataset_id)
        parsed = rank_data[dataset_id]
        values = [float(row["statistic_decimal"]) for row in parsed]
        ties = sum(values[index] == values[index + 1] for index in range(len(values) - 1))
        sorted_values = sorted(values)
        q = lambda probability: sorted_values[min(len(sorted_values) - 1, int(probability * (len(sorted_values) - 1)))]
        top = parsed[:5]
        bottom = parsed[-5:]
        rows.append(
            {
                "rank_dataset_id": dataset_id,
                "rank_file": spec["member"],
                "rank_member_bytes": rank_records[spec["member"]]["bytes"],
                "rank_member_sha256": rank_records[spec["member"]]["sha256"],
                "rank_row_count": len(parsed),
                "unique_gene_symbol_count": len({row["gene_symbol_original"] for row in parsed}),
                "finite_statistic_count": sum(math.isfinite(value) for value in values),
                "strict_descending_order": bool_text(all(values[i] > values[i + 1] for i in range(len(values) - 1))),
                "adjacent_tie_count": ties,
                "positive_count": sum(value > 0 for value in values),
                "negative_count": sum(value < 0 for value in values),
                "zero_count": sum(value == 0 for value in values),
                "absolute_statistic_gt_2_count": sum(abs(value) > 2 for value in values),
                "absolute_statistic_gt_5_count": sum(abs(value) > 5 for value in values),
                "absolute_statistic_gt_10_count": sum(abs(value) > 10 for value in values),
                "min_statistic": decimal_text(parse_decimal(str(min(values)))),
                "max_statistic": decimal_text(parse_decimal(str(max(values)))),
                "mean_statistic": decimal_text(parse_decimal(str(statistics.fmean(values)))),
                "median_statistic": decimal_text(parse_decimal(str(statistics.median(values)))),
                "quantile_01_statistic": decimal_text(parse_decimal(str(q(0.01)))),
                "quantile_25_statistic": decimal_text(parse_decimal(str(q(0.25)))),
                "quantile_75_statistic": decimal_text(parse_decimal(str(q(0.75)))),
                "quantile_99_statistic": decimal_text(parse_decimal(str(q(0.99)))),
                "top5_gene_symbols_original": "|".join(row["gene_symbol_original"] for row in top),
                "top5_statistics_original": "|".join(row["statistic_original"] for row in top),
                "bottom5_gene_symbols_original": "|".join(row["gene_symbol_original"] for row in bottom),
                "bottom5_statistics_original": "|".join(row["statistic_original"] for row in bottom),
                "rank_statistic_interpretation": "Candidate t-statistic vector by code/article; no logFC column is present.",
                "t_logfc_plausibility": "Finite two-sided values with both signs and t-like magnitude; direct logFC plausibility is unavailable because no logFC or fitted standard error is retained.",
                "direction_status": "unresolved",
            }
        )
    return rows


def _rank_overlap_rows(rank_data: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    sets = {
        dataset_id: {row["gene_symbol_original"] for row in rows}
        for dataset_id, rows in rank_data.items()
    }
    values = {
        dataset_id: {row["gene_symbol_original"]: float(row["statistic_decimal"]) for row in rows}
        for dataset_id, rows in rank_data.items()
    }
    rows = []
    for left, right in (("meki_100", "meki_500"), ("meki_100", "ko_g1"), ("meki_500", "ko_g1"), ("meki_100", "oe_250"), ("meki_100", "oe_500"), ("meki_500", "oe_500")):
        shared = sorted(sets[left] & sets[right])
        left_values = [values[left][symbol] for symbol in shared]
        right_values = [values[right][symbol] for symbol in shared]
        sign_agree = sum((left_values[i] > 0) == (right_values[i] > 0) for i in range(len(shared)))
        rows.append(
            {
                "left_dataset_id": left,
                "right_dataset_id": right,
                "left_rank_symbol_count": len(sets[left]),
                "right_rank_symbol_count": len(sets[right]),
                "shared_symbol_count": len(shared),
                "left_only_symbol_count": len(sets[left] - sets[right]),
                "right_only_symbol_count": len(sets[right] - sets[left]),
                "shared_fraction_of_left": shared and len(shared) / len(sets[left]),
                "shared_fraction_of_right": shared and len(shared) / len(sets[right]),
                "sign_agree_count": sign_agree,
                "sign_agreement_fraction": sign_agree / len(shared) if shared else None,
                "pearson_statistic_correlation": _pearson(left_values, right_values),
            }
        )
    return rows


def _step_rows(code_contents: dict[str, bytes]) -> list[dict[str, Any]]:
    rows = []
    for step in RECONSTRUCTION_STEPS:
        member = step["evidence_member"]
        text = code_contents[member].decode("utf-8", errors="replace")
        lines = text.splitlines()
        evidence_line_number = None
        evidence_text = ""
        for line_number, line in enumerate(lines, start=1):
            if step["evidence_pattern"] in line:
                evidence_line_number = line_number
                evidence_text = line.strip()
                break
        if evidence_line_number is None:
            raise ValueError(f"Could not find evidence pattern for {step['step_id']}")
        rows.append(
            {
                "step_id": step["step_id"],
                "phase": step["phase"],
                "confidence": step["confidence"],
                "evidence_member": member,
                "evidence_line_number": evidence_line_number,
                "evidence_text_original": evidence_text,
                "inference": step["inference"],
                "limitation": step["limitation"],
            }
        )

    analogue_member = "RNA-seq/RNAseq_CRISPR/TPP_CRISPR_RNAseq_analysis.r"
    analogue_text = code_contents[analogue_member].decode("utf-8", errors="replace")
    for pattern in ('G_list = getBM(filters = "ensembl_gene_id"',):
        if pattern not in analogue_text:
            raise ValueError("Analogue symbol mapping evidence is absent")
    return rows


def write_reconstruction_outputs(report: dict[str, Any], pathway_rows: list[dict[str, Any]], rank_rows_out: list[dict[str, Any]], overlap_rows: list[dict[str, Any]], step_rows: list[dict[str, Any]], note_rows: list[dict[str, Any]]) -> dict[str, str]:
    paths = {
        "ranks": DERIVED / "ets2-inhibitor-reconstruction-ranks.csv",
        "overlaps": DERIVED / "ets2-inhibitor-reconstruction-overlaps.csv",
        "pathways": DERIVED / "ets2-inhibitor-reconstruction-pathways.csv",
        "steps": DERIVED / "ets2-inhibitor-reconstruction-steps.csv",
        "notes": DERIVED / "ets2-inhibitor-reconstruction-notes.csv",
        "report": REPORTS / "ets2-inhibitor-reconstruction.json",
    }
    if rank_rows_out:
        write_csv(paths["ranks"], list(rank_rows_out[0]), rank_rows_out)
    if overlap_rows:
        write_csv(paths["overlaps"], list(overlap_rows[0]), overlap_rows)
    if pathway_rows:
        write_csv(paths["pathways"], list(pathway_rows[0]), pathway_rows)
    if step_rows:
        write_csv(paths["steps"], list(step_rows[0]), step_rows)
    if note_rows:
        write_csv(paths["notes"], list(note_rows[0]), note_rows)
    report["derived_output_sha256"] = {
        key: digest_file(path)
        for key, path in paths.items()
        if key != "report"
    }
    paths["report"].parent.mkdir(parents=True, exist_ok=True)
    paths["report"].write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {key: str(path.relative_to(ROOT)) for key, path in paths.items()}


def build_reconstruction(write_outputs: bool = True) -> dict[str, Any]:
    base_manifest = verify_manifest(ROOT / "config/ets2-benchmark-sources.json")
    program_manifest = verify_manifest(ROOT / "config/ets2-program-sources.json")
    source_data_record = next(
        item for item in program_manifest["sources"] if item["file"] == "ets2-source-data-fig5.xlsx"
    )
    archive_contents_map, archive_records, raw_count_folders, root_prefix = _archive_records_and_contents()
    rank_data: dict[str, list[dict[str, Any]]] = {}
    for spec in RANK_SPECS:
        if spec["id"] in {"meki_100", "meki_500", "ko_g1", "oe_250", "oe_500"}:
            rank_data[spec["id"]] = rank_rows(archive_contents_map[spec["member"]])
    rank_rows_out = _rank_summary_rows(rank_data, archive_records)
    overlap_rows = _rank_overlap_rows(rank_data)

    pathways, pathway_source_metrics = _read_gene_sets(
        archive_contents_map["RNA-seq/RNAseq_MEKi/macrophage_genesets_fgsea.csv"]
    )
    fgsea_report = _read_fgsea_report(
        archive_contents_map["RNA-seq/RNAseq_MEKi/fgsea_PD_0.5_mac.pathways.csv"]
    )
    fig5_sheet = _source_fig5_sheet(RAW / "ets2-source-data-fig5.xlsx", FIG5_SHEET_NAME)
    fig5_pathways, note_rows = _fig5_rows(fig5_sheet, source_data_record["actual_sha256"])
    pathway_rows: list[dict[str, Any]] = []
    pathway_report: list[dict[str, Any]] = []
    meki_500_rank = rank_data["meki_500"]
    for pathway_id in pathway_source_metrics["pathway_ids"]:
        fig5_label = next(
            label for label, mapped_id in FIG5_PATHWAY_LABEL_TO_ID.items() if mapped_id == pathway_id
        )
        fig5 = fig5_pathways[fig5_label]
        published = fgsea_report[pathway_id]
        computed_es, rank_hit_count = _signed_weighted_es(meki_500_rank, pathways[pathway_id])
        sign_flipped_rank = []
        for rank_row in reversed(meki_500_rank):
            flipped = dict(rank_row)
            flipped["statistic_decimal"] = -rank_row["statistic_decimal"]
            flipped["statistic_original"] = decimal_text(-rank_row["statistic_decimal"])
            sign_flipped_rank.append(flipped)
        sign_flipped_es, _ = _signed_weighted_es(sign_flipped_rank, pathways[pathway_id])
        source_es = parse_decimal(published["ES"])
        source_nes = parse_decimal(published["NES"])
        fig5_nes = parse_decimal(fig5["source_data_fig5_nes_original"])
        computed_decimal = parse_decimal(str(computed_es)) if computed_es is not None else None
        es_difference = abs(computed_decimal - source_es) if computed_decimal is not None else None
        fig5_nes_difference = abs(fig5_nes - source_nes)
        row = {
            "pathway_name_reported_original": fig5["pathway_name_reported_original"],
            "pathway_id": pathway_id,
            "source_data_sheet": FIG5_SHEET_NAME,
            "source_data_xml_member": fig5_sheet["xml_member"],
            "source_data_row_number": fig5["source_data_row_number"],
            "source_data_fig5_nes_original": fig5["source_data_fig5_nes_original"],
            "source_data_fig5_p_minus_log10_original": fig5["source_data_fig5_p_minus_log10_original"],
            "source_data_workbook_sha256": source_data_record["actual_sha256"],
            "archive_fgsea_source_row_number": published["source_row_number"],
            "archive_reported_es_original": published["ES"],
            "archive_reported_nes_original": published["NES"],
            "archive_reported_pval_original": published["pval"],
            "archive_reported_size_original": published["size"],
            "archive_fgsea_member_sha256": archive_records["RNA-seq/RNAseq_MEKi/fgsea_PD_0.5_mac.pathways.csv"]["sha256"],
            "gene_set_source_member_count": len(pathways[pathway_id]),
            "gene_set_duplicate_member_count": pathway_source_metrics["duplicate_member_counts"][pathway_id],
            "rank_hit_count": rank_hit_count,
            "computed_es_gsea_param_1_original": decimal_text(computed_decimal),
            "sign_flipped_rank_computed_es_original": decimal_text(
                parse_decimal(str(sign_flipped_es)) if sign_flipped_es is not None else None
            ),
            "computed_es_sign": "positive" if computed_es and computed_es > 0 else "negative" if computed_es and computed_es < 0 else "zero",
            "computed_es_abs_difference_decimal": decimal_text(es_difference),
            "computed_es_within_1e-9_of_archive_es": "true" if es_difference is not None and es_difference <= parse_decimal("1e-9") else "false",
            "reported_nes_sign": "positive" if fig5_nes > 0 else "negative" if fig5_nes < 0 else "zero",
            "computed_es_and_reported_nes_sign_agree": "true" if computed_es is not None and ((computed_es > 0) == (fig5_nes > 0)) else "false",
            "sign_flipped_rank_inverts_es_sign": "true" if computed_es is not None and sign_flipped_es is not None and ((computed_es > 0) != (sign_flipped_es > 0)) else "false",
            "fig5_and_archive_nes_abs_difference_decimal": decimal_text(fig5_nes_difference),
            "fig5_and_archive_nes_within_1e-9": "true" if fig5_nes_difference <= parse_decimal("1e-9") else "false",
            "archive_leading_edge_original": published.get("leadingEdge", ""),
        }
        pathway_rows.append(row)
        pathway_report.append(
            {
                "pathway_id": pathway_id,
                "pathway_name_reported_original": fig5_label,
                "source_gene_set_member_count": len(pathways[pathway_id]),
                "rank_hit_count": rank_hit_count,
                "reported_archive_es": published["ES"],
                "computed_es": row["computed_es_gsea_param_1_original"],
                "reported_fig5_nes": fig5["source_data_fig5_nes_original"],
                "computed_es_within_1e-9": row["computed_es_within_1e-9_of_archive_es"],
                "sign_agree": row["computed_es_and_reported_nes_sign_agree"],
            }
        )

    code_contents = {
        relative: archive_contents_map[relative]
        for relative in {
            "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
            "RNA-seq/RNAseq_CRISPR/TPP_CRISPR_RNAseq_analysis.r",
        }
    }
    step_rows = _step_rows(code_contents)
    rank_source_records = {
        spec["id"]: archive_records[spec["member"]]
        for spec in RANK_SPECS
        if spec["id"] in {"meki_100", "meki_500"}
    }
    report: dict[str, Any] = {
        "analysis_version": "1.0",
        "snapshot_date": base_manifest["snapshot_date"],
        "scope": "Bounded forensic reconstruction of released MEK inhibitor rank vectors; no raw-count recovery or treatment-direction assignment.",
        "source_manifests": [base_manifest, program_manifest],
        "archive": {
            "path": "data/raw/ets2-benchmark/ets2-public-release.zip",
            "root_prefix": root_prefix,
            "raw_count_matrix_present_in_meki_folder": raw_count_folders["RNA-seq/RNAseq_MEKi"],
            "relevant_member_records": [
                archive_records[member]
                for member in sorted(
                    {
                        "RNA-seq/RNAseq_MEKi/MEKi_treated_TPP_analysis.r",
                        "RNA-seq/RNAseq_MEKi/PD_0.1_MEK.resSYMBOL.rank.csv",
                        "RNA-seq/RNAseq_MEKi/PD_0.5_MEK.resSYMBOL.rank.csv",
                        "RNA-seq/RNAseq_MEKi/macrophage_genesets_fgsea.csv",
                        "RNA-seq/RNAseq_MEKi/fgsea_PD_0.5_mac.pathways.csv",
                        "RNA-seq/RNAseq_CRISPR/TPP_CRISPR_RNAseq_analysis.r",
                    }
                )
            ],
        },
        "source_data_fig5": {
            "path": source_data_record["file"],
            "url": source_data_record["url"],
            "bytes": source_data_record["actual_bytes"],
            "sha256": source_data_record["actual_sha256"],
            "sheet": FIG5_SHEET_NAME,
            "worksheet_xml_member": fig5_sheet["xml_member"],
            "worksheet_dimension": fig5_sheet["dimension"],
            "notes": note_rows,
        },
        "rank_vectors": {
            "summaries": rank_rows_out,
            "identity_overlaps": overlap_rows,
            "source_records": rank_source_records,
            "working_context_status": (
                "For 500 nM pathway results, drug-minus-vehicle is a well-supported working "
                "interpretation from the code contrast, source-data sheet note and numerical "
                "vector-to-ES-to-NES chain; exact per-gene historical export orientation remains unverified."
            ),
        },
        "gene_sets": pathway_source_metrics,
        "fgsea_signed_es_compatibility": {
            "rank_dataset_id": "meki_500",
            "algorithm": "weighted running-sum ES with abs(statistic)^1 hit weights and equal miss steps; no NES/null recomputation",
            "pathways": pathway_report,
            "all_nine_es_within_1e-9": all(
                row["computed_es_within_1e-9"] == "true" for row in pathway_report
            ),
            "all_nine_signed_es_agree_with_reported_fig5_nes": all(
                row["sign_agree"] == "true" for row in pathway_report
            ),
            "all_nine_fig5_and_archive_nes_within_1e-9": all(
                row["fig5_and_archive_nes_within_1e-9"] == "true"
                for row in pathway_rows
            ),
            "max_abs_computed_es_minus_archive_es": decimal_text(
                max(
                    parse_decimal(row["computed_es_abs_difference_decimal"])
                    for row in pathway_rows
                )
            ),
            "max_abs_fig5_minus_archive_nes": decimal_text(
                max(
                    parse_decimal(row["fig5_and_archive_nes_abs_difference_decimal"])
                    for row in pathway_rows
                )
            ),
            "all_nine_sign_flipped_rank_es_signs_invert": all(
                row["sign_flipped_rank_inverts_es_sign"] == "true"
                for row in pathway_rows
            ),
            "sign_flip_interpretation": (
                "Reversing the rank order while negating statistics reverses every ES sign, "
                "showing that the reported negative NES signs are compatible with the released "
                "vector orientation. Together with the code contrast and sheet 5j treated-versus-vehicle "
                "note, this supports a drug-minus-vehicle working interpretation for the 500 nM pathway "
                "result, but does not independently certify the historical per-gene export orientation."
            ),
        },
        "t_logfc_plausibility": {
            "observed_rank_properties": [
                "Both vectors are finite, unique-symbol, strictly descending and have no adjacent ties.",
                "Both contain positive and negative values; their observed ranges are compatible with a signed moderated-statistic scale.",
            ],
            "non_identifiability": [
                "The saved vectors contain one statistic per symbol and no logFC, fitted standard error, expression, count or variance column.",
                "A limma moderated t statistic is a coefficient divided by a moderated standard error, so its magnitude cannot recover the underlying logFC or count matrix.",
                "No dose-response, IC50, donor-level effect or clinical efficacy can be inferred from these ranks.",
            ],
        },
        "reconstruction_steps": step_rows,
        "confidence_interpretation": {
            "high_confidence": [
                "The code's intended sample labels, CPM filter, edgeR normalization call, donor-aware design, two drug-control contrasts and limma voom/eBayes sequence.",
                "The released MEKi vectors are finite, unique-symbol, strictly descending two-column vectors with 12,095 rows each.",
                "The nine signed ES values calculated from the 500 nM vector and the nine released GO columns match the released fgsea ES values within 1e-9; their signs agree with source-data Fig. 5 sheet 5j NES signs.",
                "The code-defined drug-versus-control contrast plus the numerical 500 nM vector-to-ES-to-publisher-NES chain makes a drug-minus-vehicle working interpretation well-supported at pathway-result level.",
            ],
            "medium_confidence": [
                "The rank files are intended t-statistic vectors because the article says GSEA lists are ranked by t-statistic and the code consumes them as fGSEA stats.",
                "The fgsea default gseaParam=1 is compatible with the exact signed ES check; package-version and null/NES behavior remain unverified.",
            ],
            "low_confidence_or_unresolved": [
                "The exact count filtering outcome, normalization factors and limma fit cannot be reproduced without controlled-access counts.",
                "The 12,801 topTable request versus 12,095 saved rows cannot be explained from the release.",
                "Symbol mapping, one-to-many collapse, tie handling, rank serialization and the active export that produced each MEKi CSV are absent.",
                "The exact historical per-gene export orientation and 100 nM vector linkage remain unresolved, despite code-defined contrast formulas and 500 nM pathway compatibility.",
            ],
        },
        "limitations": [
            "The 9-sample design from 3 donors is code metadata; saved ranks are not donor-level fits and cannot establish independent replication.",
            "The source-data Fig. 5 NES sheet states treated-versus-vehicle context, but it does not supply the underlying per-gene differential table.",
            "The reconstruction does not download or infer controlled-access EGA count matrices and does not install R packages.",
        ],
    }
    if write_outputs:
        paths = write_reconstruction_outputs(
            report,
            pathway_rows,
            rank_rows_out,
            overlap_rows,
            step_rows,
            note_rows,
        )
        report["output_paths"] = paths
        # The report is rewritten once with its output paths.  The derived table
        # hashes remain unchanged and are already part of its content.
        (ROOT / paths["report"]).write_text(
            json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return report


def main() -> None:
    report = build_reconstruction(write_outputs=True)
    compatibility = report["fgsea_signed_es_compatibility"]
    print(
        "MEKi reconstruction audit: "
        f"vectors={len(report['rank_vectors']['summaries'])}, "
        f"pathway ES matches={compatibility['all_nine_es_within_1e-9']}, "
        f"NES sign matches={compatibility['all_nine_signed_es_agree_with_reported_fig5_nes']}"
    )


if __name__ == "__main__":
    main()
