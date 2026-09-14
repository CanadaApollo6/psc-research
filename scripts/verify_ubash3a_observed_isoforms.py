"""Independent source scan using zero-based intervals and a separate parser."""

import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/ubash3a-hypothesis-test"
TOKEN = re.compile(r'(?:^|;)\s*([^;=\s]+)(?:\s+"([^"]*)"|\s+([^;\s]+)|=([^;]*))')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as source:
        for n, line in enumerate(source, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            assert len(fields) == 9
            if fields[2] in {"transcript", "mRNA", "exon"}:
                yield n, fields


def tags(text):
    result = {}
    for key, quoted, simple, gff in TOKEN.findall(text):
        value = quoted if quoted else (simple if simple else unquote(gff))
        if key in {"gene_id", "gene_name", "gene", "gene_symbol", "transcript_id", "ID", "Parent"}:
            assert key not in result or result[key] == value
            result[key] = value
    return result


def ids(feature, values):
    if "transcript_id" in values:
        return [values["transcript_id"]]
    if feature in {"transcript", "mRNA"}:
        return [values["ID"]] if "ID" in values else []
    return values.get("Parent", "").split(",") if "Parent" in values else []


def scan(path, low, high):
    wanted, total = set(), 0
    for _, f in rows(path):
        total += f[2] == "exon"
        overlap = f[0].removeprefix("chr") == "21" and int(f[3]) - 1 < high and int(f[4]) > low
        if not overlap and "UBASH3A" not in f[8] and "ENSG00000160185" not in f[8]:
            continue
        attrs = tags(f[8])
        gene = any(v == "UBASH3A" or v.split(".")[0] == "ENSG00000160185"
                   for k, v in attrs.items() if k in {"gene_id", "gene_name", "gene", "gene_symbol"})
        if overlap or gene:
            wanted.update(ids(f[2], attrs))
    found = {tid: {"exons": [], "chrom": set(), "strand": set(), "lines": []} for tid in wanted}
    for n, f in rows(path):
        if not any(tid in f[8] for tid in wanted):
            continue
        for tid in ids(f[2], tags(f[8])):
            if tid in found:
                row = found[tid]
                row["chrom"].add(f[0].removeprefix("chr"))
                row["strand"].add(f[6])
                row["lines"].append(n)
                if f[2] == "exon":
                    row["exons"].append((int(f[3]) - 1, int(f[4])))
    return found, total


def main():
    inputs = json.loads((ROOT / "config/ubash3a-test-annotation-inputs.json").read_text())
    gene = json.loads((ROOT / "data/raw/ubash3a-consequences/ensembl-gene.json").read_text())
    plan = json.loads((ROOT / "config/ubash3a-hypothesis-test-plan.json").read_text())
    old_plan = json.loads((ROOT / "config/ubash3a-consequence-plan.json").read_text())
    selected = {x["id"] for x in old_plan["expected_transcripts"]}
    refs = []
    for tx in gene["Transcript"]:
        if tx["id"] not in selected:
            continue
        for kind in ["normal", "plus29"]:
            exons = sorted((e["start"] - 1, e["end"] + (29 if kind == "plus29" and e["end"] == 42434954 else 0)) for e in tx["Exon"])
            chain = [(x[1], y[0] + 1) for x, y in zip(exons, exons[1:])]
            pair = (42434983 if kind == "plus29" else 42434954, 42437488)
            refs.append({"id": tx["id"] + "." + str(tx["version"]), "event": kind,
                         "exons": exons, "chain": chain, "suffix": chain[chain.index(pair) + 1:],
                         "coding": tx["Translation"]["start"]})
    observed = list(csv.DictReader((OUT / "observed-transcript-models.csv").open()))
    main_exons = list(csv.DictReader((OUT / "observed-exons.csv").open()))
    main_matches = {(r["source_id"], r["source_transcript_id"], r["event"], r["reference_transcript"]): r for r in csv.DictReader((OUT / "reference-chain-matches.csv").open())}
    summaries = {r["source_id"]: r for r in csv.DictReader((OUT / "annotation-summary.csv").open())}
    seen_matches, checks, source_results = set(), 0, []
    for source in inputs["annotations"]:
        path = ROOT / source["path"]
        assert sha(path) == source["sha256"]
        local, total = scan(path, gene["start"] - 1, gene["end"])
        reported = {r["source_transcript_id"]: r for r in observed if r["source_id"] == source["source_id"]}
        assert set(local) == set(reported)
        assert total == int(summaries[source["source_id"]]["total_source_exon_rows"])
        for tid, item in local.items():
            actual = reported[tid]
            exons = sorted(set(item["exons"]))
            primary_exons = [(int(r["start_1based"]) - 1, int(r["end_1based"])) for r in main_exons if r["source_id"] == source["source_id"] and r["source_transcript_id"] == tid]
            assert exons == primary_exons
            assert str(len(item["exons"]) - len(exons)) == actual["duplicate_exon_rows"]
            assert ";".join(map(str, item["lines"])) == actual["source_line_numbers"]
            eligible = item["chrom"] == {"21"} and item["strand"] == {"+"} and bool(exons)
            chain = [(x[1], y[0] + 1) for x, y in zip(exons, exons[1:])] if eligible else []
            assert chain == [tuple(pair) for pair in json.loads(actual["junction_chain_exon_end_next_start"])]
            event = "ineligible_context" if not eligible else ("plus29" if (42434983, 42437488) in chain else "normal" if (42434954, 42437488) in chain else "other")
            assert event == actual["event"]
            full, downstream = [], []
            for ref in refs:
                if event != ref["event"]:
                    continue
                target = (42434983 if event == "plus29" else 42434954, 42437488)
                suffix = chain[chain.index(target) + 1:]
                f, d = chain == ref["chain"], suffix == ref["suffix"]
                if f:
                    full.append(ref["id"])
                if d:
                    downstream.append(ref["id"])
                if f or d:
                    key = (source["source_id"], tid, event, ref["id"])
                    seen_matches.add(key)
                    primary = main_matches[key]
                    expected = {"full_internal_chain_matches": f, "downstream_chain_matches": d,
                                "reference_coding_start_is_exonic": any(a < ref["coding"] <= b for a, b in exons),
                                "retained_stop_is_fully_exonic": all(any(a < pos <= b for a, b in exons) for pos in range(42442545, 42442548)),
                                "first_exon_start_difference_nt": exons[0][0] - ref["exons"][0][0],
                                "last_exon_end_difference_nt": exons[-1][1] - ref["exons"][-1][1],
                                "exact_all_exon_coordinates": exons == ref["exons"]}
                    for k, v in expected.items():
                        assert primary[k] == str(v), (key, k)
                        checks += 1
            assert set(full) == set(filter(None, actual["full_chain_reference_matches"].split(";")))
            assert set(downstream) == set(filter(None, actual["downstream_chain_reference_matches"].split(";")))
            checks += 7
        source_results.append({"source_id": source["source_id"], "all_exon_rows": total, "locus_models": len(local), "passed": True})
        print(json.dumps(source_results[-1]), flush=True)
    assert seen_matches == set(main_matches)
    result = {"status": "pass", "implementation": "Independent regex attribute reader and zero-based source intervals; no imports from the primary analysis.",
              "script_sha256": sha(Path(__file__)), "sources": source_results, "model_rows": len(observed),
              "exon_rows": len(main_exons), "match_rows": len(main_matches), "scalar_and_structure_checks": checks,
              "input_specification_sha256": sha(ROOT / "config/ubash3a-test-annotation-inputs.json")}
    (ROOT / "reports/ubash3a-hypothesis-test-independent-check.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
