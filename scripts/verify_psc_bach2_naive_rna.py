#!/usr/bin/env python3
"""Independent, annotation-only BACH2 audit and synthetic mathematical oracle.

Never opens a real matrix. Does not import the production implementation.
Source donor joins and selected keys are written only below ignored work/.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/psc-bach2-naive-rna-independent-review"
TARGET = ("ENSG00000112182", "BACH2", "Gene Expression")
STATES = {"C0: CD4+ TN RTE": ("C0", "0", 9320),
          "C1: CD4+ TN mature": ("C1", "1", 7339)}
FIXED = {
    "work/psc-bach2-followup-review/root-method-decision.json": "f112060bda9cce45aaffe59816b2684e0d0b0bc53753f188a324bda40e161615",
    "reports/psc-bach2-followup-design-review.md": "eef554dbbef39aab08feff1a41b9238eff317c10d5cd62a46b46e77b14e5c368",
    "work/psc-bach2-followup-review/candidate-specification.json": "b54c4fc47a93afc48b6463410902e467e459049cef759dc3ba50d8fb6149d627",
    "work/psc-liver-program-inputs/source-audit/ensembl116-human-gene-crosswalk.reconstructed.tsv": "a60cf0e0ac5782e7c39269a577d1ee24be08d26c100324fe2c38b77a7f81f167",
}
REFERENCE_PATH = next(k for k in FIXED if k.endswith(".reconstructed.tsv"))

RECEIPT_PINS = {'data/raw/psc-bach2-qualification/cohort/sample01_filtered_feature_bc_matrix.tar.gz.receipt.json': '8c7bcab6d748a10e5dda659d1bc749dcc084ce3e873926167645edce2911112c',
 'data/raw/psc-bach2-qualification/cohort/sample02_filtered_feature_bc_matrix.tar.gz.receipt.json': '8b319d197d4335a4f9abbefd318c87b4a69319c855e0e6defefb75331c67d06a',
 'data/raw/psc-bach2-qualification/cohort/sample03_filtered_feature_bc_matrix.tar.gz.receipt.json': 'ee28e8883ece31feb21def1b5179d33daa7d251a4ccd8939f8e80c9b6509256e',
 'data/raw/psc-bach2-qualification/cohort/sample04_filtered_feature_bc_matrix.tar.gz.receipt.json': 'e000a1404e8f0620a7fc20455cda9638b02921db1a8f66bbdb304792cd1e812c',
 'data/raw/psc-bach2-qualification/cohort/sample06_filtered_feature_bc_matrix.tar.gz.receipt.json': 'd8ae0cc044ebd47905730eba4766cd78f993be3a6b2e9cc554ce6aa50d117300',
 'data/raw/psc-bach2-qualification/cohort/sample07_filtered_feature_bc_matrix.tar.gz.receipt.json': 'bc8b40c7668c0ee582d4cb82046ce401d3b9a92c34c4b1268048fa10cb4b9393',
 'data/raw/psc-bach2-qualification/cohort/sample08_filtered_feature_bc_matrix.tar.gz.receipt.json': 'fcc4b6d665ae541aa5264bbc77453dd5c64fd02fcd99ae8c6aebb965df37567f',
 'data/raw/psc-bach2-qualification/counts/acquisition-receipt.json': '9d8d1c9fffee8f830db70d4c7ae868be07b921bdd6cdc7a1da8ce0d99a432511',
 'work/psc-bach2-cohort-qualification/sample01/members/barcodes-integrity.json': '271ec6b37a8bc0ef04dfe1921a1b105b2ebd900e8d4dc84b0dfc7f36327b681e',
 'work/psc-bach2-cohort-qualification/sample01/members/features-integrity.json': 'e40013f124e67968e33ce40b6ead67ff4f14769ce00db43f90e1b3290c8f6d10',
 'work/psc-bach2-cohort-qualification/sample01/members/matrix-integrity.json': '28f430d993b6243f0883dd97ff69066bea645b540dc1583efbfe564e71e6e2ac',
 'work/psc-bach2-cohort-qualification/sample02/members/barcodes-integrity.json': '2dbe409c441e478bb8221be528c0484b784eba4eaf4ae97668fc3d509d1f4195',
 'work/psc-bach2-cohort-qualification/sample02/members/features-integrity.json': 'bc875e2691cc8313ed254a256576c0c2eb32e5c4ad454fb067554cc5f87727c2',
 'work/psc-bach2-cohort-qualification/sample02/members/matrix-integrity.json': 'ee24de7ec06888040b8b10a5fdabd972f2f0dd7fad4b5e0d329bf7693fad7d0a',
 'work/psc-bach2-cohort-qualification/sample03/members/barcodes-integrity.json': '8482fb40cdefc6dea78a6636fa7e05d8b3ba25da5297c4e57f5adeef8835b931',
 'work/psc-bach2-cohort-qualification/sample03/members/features-integrity.json': 'dae645256402c6fc64b68b7e1cafa343f1b4bfdfbd0b3a2732ce42f2969d78ec',
 'work/psc-bach2-cohort-qualification/sample03/members/matrix-integrity.json': '438e50bc2779667c8d81318f22c8bee2fbe139afd5bda9c3a2764a400c4a8724',
 'work/psc-bach2-cohort-qualification/sample04/members/barcodes-integrity.json': 'b3c69eb2b8133f1bd742e5c427790364b8f13c4bac76b74b91cd933d2151a92f',
 'work/psc-bach2-cohort-qualification/sample04/members/features-integrity.json': '6817625035934251b3fdd45463cf4bb65bb4d9f15aa195df75cf8063926c6254',
 'work/psc-bach2-cohort-qualification/sample04/members/matrix-integrity.json': '5199d8e392aec70b0c1bbea43e748d3f689720bb13869ac5ffece56d87fb585e',
 'work/psc-bach2-cohort-qualification/sample05/members/barcodes-integrity.json': 'd9518524a9a540106912709940f2c52d18203592775afb99121efeacda100c68',
 'work/psc-bach2-cohort-qualification/sample05/members/features-integrity.json': '9a591ceb94d46248b6d90757f126ff3391f232dbc3190860f9f866e5bda76d37',
 'work/psc-bach2-cohort-qualification/sample05/members/matrix-integrity.json': 'f1d5a483040aa29c3fb9b23d3caf75f311e89bc5772d4b4724fd1cc52532eb84',
 'work/psc-bach2-cohort-qualification/sample06/members/barcodes-integrity.json': 'b2c58ee65e53d0e5a1792c888a97be61f08eff5e782bafe9421c87cf64e025fd',
 'work/psc-bach2-cohort-qualification/sample06/members/features-integrity.json': '1a3f016dab1457f586a224a09435bfafce9043ff0b5496baaf82d5f25cc31f01',
 'work/psc-bach2-cohort-qualification/sample06/members/matrix-integrity.json': '872e01733b8f481cf055f2a0adebfc0f4da7030b474169392633e98b2e78eced',
 'work/psc-bach2-cohort-qualification/sample07/members/barcodes-integrity.json': '820b7ab66b94ff55f1b0776539d08d9b12e48e60316fdfc765da38c5efc27d5e',
 'work/psc-bach2-cohort-qualification/sample07/members/features-integrity.json': '9891a6b7045175a915b80f97fa4577da7033d9d0a9423b662b6001497e2b2184',
 'work/psc-bach2-cohort-qualification/sample07/members/matrix-integrity.json': '0a736c129bb0ecd854d4e63543ac238974cd5c8c0f46ed8e87e1bc9ec32d80db',
 'work/psc-bach2-cohort-qualification/sample08/members/barcodes-integrity.json': 'ad1b4bc429bfb8dd34551c9f511c512df5bca6b9b90b642eb9cab71428a3d8fe',
 'work/psc-bach2-cohort-qualification/sample08/members/features-integrity.json': '509241330c0615d0da6523b4cd73b0a70d81c3c21d9071ae6b71ef2663a867f8',
 'work/psc-bach2-cohort-qualification/sample08/members/matrix-integrity.json': '0ccdcbd2f85f243bee61cd21f6f31f241b7de3881dedc079a396c402e48f9ad1'}
ARCHIVE_RELEASE_PINS = {'data/raw/psc-bach2-qualification/cohort/sample01_filtered_feature_bc_matrix.tar.gz': {'bytes': 54449959,
                                                                                        'sha256': 'cbede10441fecfc1013077088829ddb9916140d9f8489c371c98efb3492837a6'},
 'data/raw/psc-bach2-qualification/cohort/sample02_filtered_feature_bc_matrix.tar.gz': {'bytes': 44540169,
                                                                                        'sha256': 'b47a0d5322ff6d560e385d5b205b4081b9c88dc56bcf3a7e702863d9a41c64e8'},
 'data/raw/psc-bach2-qualification/cohort/sample03_filtered_feature_bc_matrix.tar.gz': {'bytes': 43425341,
                                                                                        'sha256': 'f0a957fadccf7998a8c366e537483d23bceba1b095b5cbf90e3920402dcf2091'},
 'data/raw/psc-bach2-qualification/cohort/sample04_filtered_feature_bc_matrix.tar.gz': {'bytes': 36719700,
                                                                                        'sha256': '0a7bc817823e121b929a1716abb1f32b132c909b19b0b159a17f2484393f7d52'},
 'data/raw/psc-bach2-qualification/cohort/sample06_filtered_feature_bc_matrix.tar.gz': {'bytes': 36674386,
                                                                                        'sha256': 'e251f2df1148f9ae06c831360295782b681bece649395d2b2bbfe2aaf4b8e826'},
 'data/raw/psc-bach2-qualification/cohort/sample07_filtered_feature_bc_matrix.tar.gz': {'bytes': 29629283,
                                                                                        'sha256': '631ef4435687e966ce3fe743691bd35d601928460c096c1e42ff0701cdf9c9a9'},
 'data/raw/psc-bach2-qualification/cohort/sample08_filtered_feature_bc_matrix.tar.gz': {'bytes': 52982953,
                                                                                        'sha256': '1808db1c3333642526eedd0b880c98cd4df4adcdd5d62866f17dae6d56666e03'},
 'data/raw/psc-bach2-qualification/counts/sample05_filtered_feature_bc_matrix.tar.gz': {'bytes': 26893378,
                                                                                        'sha256': '5c57cc1c8eb48ce0208062b92fcd056b39e0f90e8eb6cb712ccd0d76f95dda29'}}
METADATA_PIN = ("data/raw/psc-bach2-qualification/05-scrna_meta.tsv.body", "0ba598b3edd30dc7f2ba76f0e7779725e60cce50daa01a0e2e04faebcec7c26c", 5914294)
SDRF_PIN = ("data/raw/psc-bach2-qualification/01-E-MTAB-14013.sdrf.txt.body", "58273fa29690bdcaff0060161519ff6bbbee9e2af0df4bef03ef91a94184b10d", 3989)
RNA_AXIS_SHA = "aef32a0841aec54ce17a661c064501697b2b7885bb53dff8a930eef74230ae2d"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_pin(path):
    # Opaque bytes only. No archive or numerical member decoding.
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return {"sha256": digest.hexdigest(), "bytes": size}


def pinned(root, relative, expected, expected_size=None):
    path = root / relative
    require(path.is_file() and not path.is_symlink(), "source missing/link: " + relative)
    result = file_pin(path)
    require(result["sha256"] == expected, "source hash changed: " + relative)
    if expected_size is not None:
        require(result["bytes"] == expected_size, "source size changed: " + relative)
    return {"path": relative, **result}


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, separators=(",", ":"), sort_keys=True,
                                     ensure_ascii=False).encode()).hexdigest()


def source_mapping(path):
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = rows.pop(0)
    scalar = ["Characteristics[individual]", "Characteristics[sex]",
              "Characteristics[disease]", "Characteristics[genotype]", "Factor Value[genotype]"]
    require(all(header.count(key) == 1 for key in scalar), "ambiguous SDRF scalar field")
    derived = [i for i, key in enumerate(header) if key == "Derived Array Data File"]
    require(len(derived) == 3, "changed repeated SDRF file columns")
    result = {}
    for row in rows:
        require(len(row) == len(header), "SDRF width")
        files = [row[i] for i in derived if re.fullmatch(r"sample0[1-8]_filtered_feature_bc_matrix[.]tar[.]gz", row[i])]
        require(len(files) == 1, "missing/ambiguous exact filtered filename")
        sample = files[0].split("_", 1)[0]
        fields = {key: row[header.index(key)] for key in scalar}
        group = fields["Characteristics[genotype]"]
        donor = fields["Characteristics[individual]"]
        require(sample not in result and donor, "duplicate sample or empty donor")
        require(group in ("SNP", "noSNP") and fields["Factor Value[genotype]"] == group,
                "genotype/factor mismatch")
        require(fields["Characteristics[sex]"] == "male" and
                fields["Characteristics[disease]"] == "primary sclerosing cholangitis", "donor scope")
        result[sample] = {"donor": donor, "group": group, "sex": "male", "archive": files[0]}
    require(set(result) == {f"sample0{i}" for i in range(1, 9)}, "eight samples required")
    require(len({r["donor"] for r in result.values()}) == 8, "eight independent source donor labels required")
    require(Counter(r["group"] for r in result.values()) == {"SNP": 4, "noSNP": 4}, "four/four source groups")
    return result


def cell_key(sample, barcode):
    prefix = sample + "_"
    return (sample, barcode[len(prefix):] if barcode.startswith(prefix) else barcode)


def metadata_selection(path, mapping, source_barcodes):
    required = {"sample_name", "sex", "condition", "seurat_clusters", "paper_clusters",
                "paper_clusters_short", "UMAP_1", "UMAP_2", "barcode"}
    keys, selected, totals, states = set(), [], Counter(), Counter()
    positions = {sample: {barcode: i for i, barcode in enumerate(values, 1)}
                 for sample, values in source_barcodes.items()}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(len(reader.fieldnames or []) == len(required) and set(reader.fieldnames) == required,
                "metadata schema changed")
        for row in reader:
            require(None not in row and all(value is not None for value in row.values()), "metadata width")
            sample = row["sample_name"]
            require(sample in mapping, "unknown metadata sample")
            require(row["sex"] == "male" and row["condition"] == mapping[sample]["group"], "metadata source mapping mismatch")
            key = cell_key(sample, row["barcode"])
            require(key not in keys, "duplicate full sample/barcode key")
            require(key[1] in positions[sample], "retained metadata cell absent from its own sample")
            keys.add(key)
            totals[sample] += 1
            # Reject a C0/C1 alias/label discrepancy rather than silently losing cells.
            wanted = row["paper_clusters"] in STATES
            marked = row["paper_clusters_short"] in ("C0", "C1") or row["seurat_clusters"] in ("0", "1")
            require(wanted == marked, "C0/C1 crosswalk scope mismatch")
            if wanted:
                short, seurat, _ = STATES[row["paper_clusters"]]
                require((row["paper_clusters_short"], row["seurat_clusters"]) == (short, seurat), "C0/C1 crosswalk mismatch")
                states[row["paper_clusters"]] += 1
                selected.append({"sample": sample, "donor": mapping[sample]["donor"],
                                 "group": mapping[sample]["group"], "barcode": key[1],
                                 "source_column_1based": positions[sample][key[1]],
                                 "paper_clusters": row["paper_clusters"]})
    return keys, sorted(selected, key=lambda r: (r["sample"], r["source_column_1based"])), totals, states


def reference_identity(path):
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    require(rows[0] == ["Gene stable ID", "NCBI gene (formerly Entrezgene) ID", "HGNC symbol"], "reference header")
    require(rows[-1] == ["[success]"], "reference success footer")
    relations = rows[1:-1]
    require(len(relations) == 91748 and all(len(r) == 3 for r in relations), "reference complete relationships")
    require(len({r[0] for r in relations}) == 86411, "reference ID count")
    forward = {r[2] for r in relations if r[0] == TARGET[0] and r[2]}
    reverse = {r[0] for r in relations if r[2] == TARGET[1]}
    require(forward == {TARGET[1]} and reverse == {TARGET[0]}, "not unique reciprocal exact current gene identity")
    return {"release": "Ensembl 116, June 2026", "exact_ensembl_id": TARGET[0],
            "exact_symbol": TARGET[1], "reciprocal_identity_passed": True,
            "relationship_rows": len(relations), "unique_ensembl_ids": 86411,
            "cache_provenance": "Exact historical bytes reconstructed from pinned original relationship export; not a new retrieval",
            "original_retrieved_at_utc": "2026-09-12T17:40:41Z",
            "source_url": "https://jun2026.archive.ensembl.org/biomart/martservice",
            "genome_build_or_variant_identity_inferred": False}



def receipt_contract(root):
    """Read source integrity/acquisition receipts only; no quantity JSON."""
    pins = [pinned(root,path,sha) for path,sha in RECEIPT_PINS.items()]
    members, archives = {}, {}
    for index in range(1,9):
        sample = f"sample{index:02d}"
        members[sample] = {}
        for kind in ("features","barcodes","matrix"):
            path = f"work/psc-bach2-cohort-qualification/{sample}/members/{kind}-integrity.json"
            row = json.loads((root/path).read_text())
            require(row["decoded_path"] == f"work/psc-bach2-cohort-qualification/{sample}/members/{kind}.txt", "decoded receipt path")
            require(row["nested_gzip_CRC_ISIZE_and_EOF_passed"] is True and row["source_payload_matches_tar"] is True, "closed member integrity")
            members[sample][kind] = row
        archive_path = f"data/raw/psc-bach2-qualification/{'counts' if index == 5 else 'cohort'}/{sample}_filtered_feature_bc_matrix.tar.gz"
        receipt_path = "data/raw/psc-bach2-qualification/counts/acquisition-receipt.json" if index == 5 else archive_path + ".receipt.json"
        acquisition = json.loads((root/receipt_path).read_text())
        fixed = ARCHIVE_RELEASE_PINS[archive_path]
        require(acquisition["complete"] is True and acquisition["sha256"] == fixed["sha256"] and acquisition["body_bytes"] == fixed["bytes"], "closed acquisition pin")
        expected = acquisition["authorized_archive_bytes"] if index == 5 else acquisition["expected_bytes"]
        require(acquisition["body_bytes"] == expected, "closed acquisition bytes")
        archives[sample] = {"path":archive_path,"bytes":fixed["bytes"],"sha256":fixed["sha256"],
                            "source_URL":acquisition["url"],"retrieved_at_UTC":acquisition["completed_at_utc"]}
    return members, archives, pins


def source_audit(root=ROOT):
    pins = [pinned(root, key, value) for key, value in FIXED.items()]
    for path,sha,size in (METADATA_PIN,SDRF_PIN):
        pins.append(pinned(root,path,sha,size))
    members, archives, receipt_pins = receipt_contract(root)
    pins.extend(receipt_pins)
    mapping = source_mapping(root / SDRF_PIN[0])
    barcodes, axes = {}, []
    for sample in sorted(mapping):
        receipt = members[sample]
        archive = archives[sample]
        require(mapping[sample]["archive"] == Path(archive["path"]).name, "SDRF archive identity")
        pins.append(pinned(root,archive["path"],archive["sha256"],archive["bytes"]))
        for kind in ("features", "barcodes"):
            member = receipt[kind]
            pins.append(pinned(root,member["decoded_path"],member["decoded_sha256"],member["decoded_bytes"]))
        feature_path = root / receipt["features"]["decoded_path"]
        with feature_path.open(newline="") as handle:
            axis = [tuple(row) for row in csv.reader(handle, delimiter="\t")]
        require(len(axis) == 36639 and all(len(r) == 3 for r in axis), "full feature axis shape")
        require(len({r[0] for r in axis}) == len(axis), "duplicate feature IDs")
        require(Counter(r[2] for r in axis) == {"Gene Expression": 36601, "Antibody Capture": 38}, "RNA/ADT axes")
        require([i for i, row in enumerate(axis, 1) if row == TARGET] == [12314], "exact unique source triple/one-based position")
        require(sum(r[0] == TARGET[0] for r in axis) == 1 and sum(r[1] == TARGET[1] for r in axis) == 1,
                "target ID/symbol ambiguity")
        rna = [list(r) for r in axis if r[2] == "Gene Expression"]
        require(canonical_digest(rna) == RNA_AXIS_SHA, "RNA axis receipt mismatch")
        axes.append(axis)
        barcode_path = root / receipt["barcodes"]["decoded_path"]
        barcodes[sample] = barcode_path.read_text().splitlines()
        require(all(barcodes[sample]) and len(set(barcodes[sample])) == len(barcodes[sample]), "barcode duplicate/empty")
    require(all(axis == axes[0] for axis in axes), "full ordered source axes differ")
    keys, selected, retained, states = metadata_selection(root / METADATA_PIN[0], mapping, barcodes)
    require(len(keys) == 55460 and len(selected) == 16659, "fixed retained/selected count")
    require(states == {k: v[2] for k, v in STATES.items()}, "fixed state membership counts")
    require(set(r["sample"] for r in selected) == set(mapping), "missing selected donor")
    require(all({r["paper_clusters"] for r in selected if r["sample"] == sample} == set(STATES) for sample in mapping),
            "each donor must contain both states")
    source_cells = sum(len(b) for b in barcodes.values())
    require(source_cells == 59315 and source_cells - len(keys) == 3855, "source-only extras changed")
    group_sizes = Counter(mapping[s]["group"] for s in mapping)
    summary = {"status": "ANNOTATION_ONLY_PASS_NO_EXPRESSION_AUTHORIZATION", "source_pins": pins,
               "source_samples": 8, "source_donors": 8, "group_donors": dict(group_sizes),
               "all_source_donors_male_PSC": True, "retained_cells": len(keys), "selected_cells": len(selected),
               "selected_states": dict(states), "source_only_cells_not_added": source_cells - len(keys),
               "RNA_features": 36601, "ADT_features_kept_separate": 38, "complete_axes_identical": True,
               "target": list(TARGET), "source_position_1based": 12314,
               "canonical_RNA_axis_sha256": canonical_digest([list(r) for r in axes[0] if r[2] == "Gene Expression"]),
               "private_selection_canonical_sha256": canonical_digest(selected),
               "current_reference_identity": reference_identity(root / REFERENCE_PATH),
               "full_sample_plus_barcode_join": True, "clinical_join": False,
               "real_numeric_payloads_parsed": 0, "real_matrices_opened": 0, "network_requests": 0,
               "opaque_archive_bytes_hashed_only": True}
    return summary, {"mapping":mapping,"retained_cell_counts":dict(retained),"source_cell_counts":{s:len(v) for s,v in barcodes.items()},"selected":selected}



def audit_compiled_plan(relative, independent_summary, independent_private, root=ROOT):
    """Compare a compiled author plan to independently decoded annotations.

    Production code is not imported. Matrix hashes are compared to the pinned
    qualification receipt, never obtained by opening a real matrix here.
    """
    require(isinstance(relative,str) and re.fullmatch(r"work/psc-bach2-naive-rna-plan/[A-Za-z0-9][A-Za-z0-9_-]{0,79}",relative)
            and Path(relative).name not in ("runs","independent-review"), "compiled plan must be an authorized preparation-directory role")
    directory = root / relative
    require(not any(p.is_symlink() for p in [directory,*directory.parents]), "compiled plan symlink")
    raw = (directory / "plan.json").read_bytes()
    plan = json.loads(raw)
    selection_raw = (directory / "selection.private.json").read_bytes()
    selection = json.loads(selection_raw)
    schema_raw = (directory / "public-schema.json").read_bytes()
    def encoded(value):
        return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    def encoded_sha(value):
        return hashlib.sha256(encoded(value)).hexdigest()
    require(raw == encoded(plan) and selection_raw == encoded(selection), "noncanonical compiled files")
    require(plan["status"] == "METHOD_PREPARATION_ONLY_NOT_EXECUTION_AUTHORIZATION", "not a preparation plan")
    require(plan["selection_sha256"] == encoded_sha(selection), "compiled selection digest")
    require(plan["source_manifest_sha256"] == encoded_sha(plan["source_manifest"]), "compiled source manifest digest")
    require(schema_raw == encoded(plan["public_schema"]) and plan["public_schema_sha256"] == encoded_sha(plan["public_schema"]), "compiled public schema digest")
    expected_candidate = json.loads((root / "work/psc-bach2-followup-review/candidate-specification.json").read_text())
    require(plan["fixed_candidate"] == expected_candidate, "candidate method changed")
    metadata_path = root / "data/raw/psc-bach2-qualification/05-scrna_meta.tsv.body"
    with metadata_path.open(newline="") as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    retained_keys = [list(cell_key(row["sample_name"],row["barcode"])) for row in metadata]
    require(selection["retained_keys"] == retained_keys, "complete ordered retained keys changed")
    states = {value[0]:key for key,value in STATES.items()}
    normalized = []
    require(len(selection["donors"]) == 8, "compiled donor count")
    for donor in selection["donors"]:
        sample = donor["sample"]
        require(sample in independent_private["mapping"], "compiled sample")
        mapping = independent_private["mapping"][sample]
        require((donor["donor"],donor["group"]) == (mapping["donor"],mapping["group"]), "compiled donor mapping")
        require(donor["retained_count"] == independent_private["retained_cell_counts"][sample], "compiled retained donor count")
        state_counts = Counter(cell["state"] for cell in donor["selected"])
        require(dict(state_counts) == donor["state_counts"], "compiled state counts")
        for cell in donor["selected"]:
            normalized.append({"sample":sample,"donor":donor["donor"],"group":donor["group"],
                               "barcode":cell["barcode"],"source_column_1based":cell["column_1based"],
                               "paper_clusters":states[cell["state"]]})
    normalized.sort(key=lambda r:(r["sample"],r["source_column_1based"]))
    require(normalized == independent_private["selected"], "compiled complete fixed selection disagrees")
    members, archives, _ = receipt_contract(root)
    matrix_pins = plan["source_manifest"]["matrix_expected_pins_NOT_read"]
    archive_pins = plan["source_manifest"]["archive_expected_pins_NOT_read"]
    require(len(matrix_pins) == len(archive_pins) == 8, "compiled source count")
    require({p["sample"] for p in matrix_pins} == {p["sample"] for p in archive_pins} == set(independent_private["mapping"]), "compiled source sample universe")
    for pin in matrix_pins:
        old = members[pin["sample"]]["matrix"]
        require((pin["path"],pin["bytes"],pin["sha256"],pin["rows"],pin["columns"]) ==
                (old["decoded_path"],old["decoded_bytes"],old["decoded_sha256"],36639,independent_private["source_cell_counts"][pin["sample"]]),
                "compiled matrix pin differs from closed integrity receipt and annotation dimensions")
    for pin in archive_pins:
        require({key:pin[key] for key in archives[pin["sample"]]} == archives[pin["sample"]], "compiled archive pin/provenance differs")
    annotation_pins = plan["source_manifest"]["annotation_and_receipt_pins"]
    require(any(pin["path"] == ".gitignore" for pin in annotation_pins), "private ignore policy not pinned")
    for pin in annotation_pins + plan["implementation"]:
        require("/matrix.txt" not in pin["path"] and not pin["path"].endswith(".mtx"), "unexpected real matrix validation read")
        pinned(root,pin["path"],pin["sha256"],pin["bytes"])
    feature_file = root / members["sample01"]["features"]["decoded_path"]
    with feature_file.open(newline="") as handle:
        expected_RNA = [i for i,row in enumerate(csv.reader(handle,delimiter="\t"),1) if row[2] == "Gene Expression"]
    require(plan["RNA_rows_1based"] == expected_RNA, "compiled full RNA-only denominator axis")
    import importlib.metadata
    import sys
    runtime = {"python":sys.version,"cache_tag":sys.implementation.cache_tag,
               "executable":str(Path(sys.executable).absolute()),"executable_sha256":file_pin(Path(sys.executable))["sha256"],
               "packages":{}}
    import base64
    import io
    for package in ("numpy","scipy"):
        distribution = importlib.metadata.distribution(package)
        record = distribution.read_text("RECORD")
        verified = 0
        for file, digest, size in csv.reader(io.StringIO(record)):
            if not digest:
                require(file.endswith(("/RECORD", ".pyc")), "unhashed runtime file")
                continue
            algorithm, expected = digest.split("=", 1)
            require(algorithm == "sha256", "unexpected runtime digest")
            native = Path(distribution.locate_file(file))
            actual = file_pin(native)
            expected_hex = base64.urlsafe_b64decode(expected + "=" * (-len(expected) % 4)).hex()
            require(actual["sha256"] == expected_hex and actual["bytes"] == int(size), "native installed runtime differs")
            verified += 1
        runtime["packages"][package] = {"version":distribution.version,
            "RECORD_sha256":hashlib.sha256(record.encode()).hexdigest(), "verified_installed_files":verified}
    require(runtime == plan["runtime"], "compiled native runtime differs")
    return {"status":"INDEPENDENT_COMPILED_SELECTION_SOURCE_RUNTIME_PASS",
            "plan_directory":relative,"plan_sha256":hashlib.sha256(raw).hexdigest(),
            "selection_sha256":hashlib.sha256(selection_raw).hexdigest(),
            "source_manifest_sha256":encoded_sha(plan["source_manifest"]),
            "runtime_sha256":encoded_sha(runtime),"public_schema_sha256":hashlib.sha256(schema_raw).hexdigest(),
            "implementation":plan["implementation"],"all_retained_keys_compared":55460,
            "selected_keys_and_six_fields_compared":16659,"all_donor_joins_compared":8,
            "matrix_pins_compared_to_immutable_qualification":8,"opaque_archive_pins_compared":8,
            "annotation_and_receipt_hashes_verified":len(annotation_pins),
            "RNA_only_complete_source_positions_compared":36601,
            "real_matrices_opened":0,"expression_authorized":False}


def independent_y(b, r):
    require(isinstance(b, int) and not isinstance(b, bool) and isinstance(r, int) and not isinstance(r, bool), "integer source totals")
    require(0 <= b <= r and r > 0, "invalid full RNA denominator or BACH2 total")
    with localcontext() as ctx:
        ctx.prec = 90
        return float((Decimal(1) + Decimal(1000000) * Decimal(b) / Decimal(r)).ln() / Decimal(2).ln())



def independent_unit_points(units):
    """Exact product identity, independently evaluated by 1100-digit direct ln.

    This does not use the production 320-digit near-one log1p-series. Precision
    here is a computational check, not a statement of biological precision.
    """
    require(set(units) == {"SNP","noSNP"} and all(len(v) == 4 for v in units.values()), "four/four units required")
    factors = {}
    for group, pairs in units.items():
        require(all(type(b) is int and type(r) is int and 0 <= b <= r and r > 0 for b,r in pairs), "invalid source units")
        factors[group] = [Fraction(r + 1000000*b,r) for b,r in pairs]
    products = {group:math.prod(values) for group,values in factors.items()}
    def exact_product_log(q,divisor):
        if q == 1:
            return Decimal(0)
        with localcontext() as context:
            context.prec = 1100
            ratio = Decimal(q.numerator) / Decimal(q.denominator)
            return ratio.ln() / Decimal(2).ln() / Decimal(divisor)
    means = {group:exact_product_log(product,4) for group,product in products.items()}
    delta = exact_product_log(products["SNP"] / products["noSNP"],4)
    leaveouts = []
    for group in ("SNP","noSNP"):
        for factor in factors[group]:
            pa = products["SNP"] / factor if group == "SNP" else products["SNP"]
            pb = products["noSNP"] / factor if group == "noSNP" else products["noSNP"]
            na,nb = (3,4) if group == "SNP" else (4,3)
            leaveouts.append(exact_product_log(pa**nb/pb**na,na*nb))
    return {"means":means,"delta":delta,"leaveouts":leaveouts}


def exact_moments(values):
    require(len(values) >= 2 and all(math.isfinite(v) for v in values), "invalid donor Y")
    q = [Fraction(float(v)) for v in values]
    mean = sum(q, Fraction(0)) / len(q)
    variance = sum(((x - mean) ** 2 for x in q), Fraction(0)) / (len(q) - 1)
    return mean, variance


def fraction_sqrt(q):
    with localcontext() as ctx:
        ctx.prec = 90
        return float((Decimal(q.numerator) / Decimal(q.denominator)).sqrt())


def independent_welch(snp, no_snp):
    """Independent exact-rational moments; incomplete-beta Student-t probability.

    No production helpers and no scipy.stats Welch/t implementation are used.
    Genuine tiny differences remain exact in the rational variance. Conversion
    to floating point occurs only for finite model outputs and special functions.
    """
    from scipy.special import betainc
    from scipy.optimize import brentq
    require(len(snp) == len(no_snp) == 4, "fixed all-eight endpoint")
    ma, sa = exact_moments(snp)
    mb, sb = exact_moments(no_snp)
    va, vb = sa / 4, sb / 4
    delta = float(ma - mb)
    se2 = va + vb
    base = {"delta": delta, "mean_SNP": float(ma), "mean_noSNP": float(mb),
            "se": None, "df": None, "ci95": None, "p_two_sided": None,
            "inference_status": "point_only"}
    if not se2:
        return base
    se = fraction_sqrt(se2)
    df = float(se2 ** 2 / (va ** 2 / 3 + vb ** 2 / 3))
    if not math.isfinite(se) or se <= 0 or not math.isfinite(df) or df <= 0:
        return base
    # P(|T_nu| >= |t|) = I_{nu/(nu+t^2)}(nu/2, 1/2).
    t2 = (ma - mb) ** 2 / se2
    probability_arg = float(Fraction.from_float(df) / (Fraction.from_float(df) + t2))
    p = float(betainc(df / 2, 0.5, probability_arg))
    critical = brentq(lambda t: float(betainc(df / 2, 0.5, df / (df + t*t))) - 0.05,
                      0.0, 100.0, xtol=5e-15)
    return {**base, "se": se, "df": df, "ci95": [delta - critical * se, delta + critical * se],
            "p_two_sided": p, "inference_status": "available"}


def independent_loo(snp, no_snp):
    require(len(snp) == len(no_snp) == 4, "all eight whole-donor omissions required")
    a, b = [Fraction(float(v)) for v in snp], [Fraction(float(v)) for v in no_snp]
    return [float((sum(a) - value) / 3 - sum(b) / 4) for value in a] + [
        float(sum(a) / 4 - (sum(b) - value) / 3) for value in b]


def write_audit(output, summary, private):
    candidate = Path(output)
    require(".." not in candidate.parts, "output traversal")
    candidate = candidate.absolute()
    require(candidate.is_relative_to(WORK) and candidate != WORK, "output must be under owned ignored work")
    require(not any(p.is_symlink() for p in [candidate, *candidate.parents]), "output symlink")
    candidate.mkdir(parents=True, exist_ok=False)
    for name, obj in (("aggregate-audit.json", summary), ("private-selection.json", private)):
        with (candidate / name).open("x") as handle:
            json.dump(obj, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-audit", action="store_true", help="metadata/annotations and opaque archive hashes only")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compiled-plan", help="also compare fixed author plan, without importing production code")
    args = parser.parse_args()
    require(args.source_audit and args.output is not None, "explicit source-audit and fresh ignored output required")
    summary, private = source_audit()
    if args.compiled_plan:
        summary["compiled_plan_audit"] = audit_compiled_plan(args.compiled_plan, summary, private)
    write_audit(args.output, summary, private)
    print(json.dumps({k: summary[k] for k in ("status", "source_donors", "retained_cells", "selected_cells", "real_numeric_payloads_parsed")}))


if __name__ == "__main__":
    main()
