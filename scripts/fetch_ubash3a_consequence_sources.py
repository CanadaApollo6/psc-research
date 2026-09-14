"""Capture or verify the bounded public inputs for the UBASH3A consequence audit.

Ensembl and UniProt requests use the corresponding installed research-skill client.
The saved UTF-8 response text is hashed; compressed transport bytes are not retained.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/ubash3a-consequences"
MANIFEST = ROOT / "config/ubash3a-consequence-sources.json"
PLAN = ROOT / "config/ubash3a-consequence-plan.json"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def skill_client(service, explicit=None):
    if explicit:
        path = Path(explicit)
    else:
        base = Path.home() / ".codex/plugins/cache/openai-curated-remote/life-science-research"
        matches = list(base.glob(f"*/skills/{service}-skill/scripts/rest_request.py"))
        if not matches:
            raise FileNotFoundError(f"Supply --{service}-client for the installed research-skill script")
        path = max(matches, key=lambda p: p.stat().st_mtime_ns)
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def rest(source_id, service, path, params=None):
    return {"source_id": source_id, "service": service, "path": path, "params": params or {}, "file": source_id + ".json"}


def requests_for(stage, plan):
    if stage == "metadata":
        return [
            rest("ensembl-release", "ensembl", "info/data"),
            rest("ensembl-assembly", "ensembl", "info/assembly/homo_sapiens"),
            rest("ensembl-gene", "ensembl", "lookup/id/ENSG00000160185", {"expand": 1}),
            rest("ensembl-variant", "ensembl", "variation/human/rs1893592"),
            rest("uniprot-search", "uniprot", "uniprotkb/search", {"query": "gene:UBASH3A AND organism_id:9606 AND reviewed:true", "size": 10, "format": "json"}),
        ]
    if stage == "sequences":
        gene = json.loads((RAW / "ensembl-gene.json").read_text())
        old = json.loads((ROOT / "data/raw/UBASH3A-transcripts-GRCh38.json").read_text())
        by_id = {t["id"]: t for t in gene["Transcript"]}
        old_by_id = {t["id"]: t for t in old["Transcript"]}
        if gene["id"] != plan["gene_id"] or gene["strand"] != 1 or gene["assembly_name"] != "GRCh38":
            raise ValueError("Gene identity, strand or assembly changed")
        output = [rest("ensembl-gene-sequence", "ensembl", f"sequence/region/human/21:{gene['start']}..{gene['end']}:1")]
        for expected in plan["expected_transcripts"]:
            tid = expected["id"]
            current, prior = by_id[tid], old_by_id[tid]
            for key in ("version", "biotype", "start", "end", "strand", "length", "Translation"):
                if current.get(key) != prior.get(key):
                    raise ValueError(f"Annotation drift for {tid}: {key}")
            exon_key = lambda t: sorted((e["id"], e["start"], e["end"], e["strand"]) for e in t["Exon"])
            if exon_key(current) != exon_key(prior):
                raise ValueError(f"Exon annotation drift for {tid}")
            output.extend([
                rest(tid + "-cdna", "ensembl", "sequence/id/" + tid, {"type": "cdna"}),
                rest(tid + "-cds", "ensembl", "sequence/id/" + tid, {"type": "cds"}),
                rest(tid + "-protein", "ensembl", "sequence/id/" + current["Translation"]["id"], {"type": "protein"}),
            ])
        hits = json.loads((RAW / "uniprot-search.json").read_text())["results"]
        exact = [h for h in hits if h["organism"]["taxonId"] == 9606 and any(g.get("geneName", {}).get("value") == "UBASH3A" for g in h.get("genes", []))]
        if len(exact) != 1:
            raise ValueError(f"Expected one reviewed human UBASH3A record; found {len(exact)}")
        output.append(rest("uniprot-entry", "uniprot", "uniprotkb/" + exact[0]["primaryAccession"], {"format": "json"}))
        return output
    if stage == "literature":
        ids = {"mucaki2020": "PMC7066660", "lindeboom2016": "PMC5045715", "lindeboom2019": "PMC6858879"}
        output = [{"source_id": name, "service": "https", "url": f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML", "file": name + ".xml"} for name, pmc in ids.items()]
        # Europe PMC returned 404 for this article's XML. Reuse the complete,
        # previously pinned primary-paper HTML without overwriting that cache.
        output.append({"source_id": "ge2018", "service": "local_reuse", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6018660/", "file": "ge2018.html", "origin_manifest": "config/mechanism-sources.json", "origin_source_id": "Ge2018-PMC"})
        return output
    raise ValueError(stage)


def fetch_one(spec, existing, clients, offline):
    target = RAW / spec["file"]
    old = existing.get(spec["source_id"])
    if old and target.exists():
        if not target.is_file() or digest(target) != old["sha256"] or target.stat().st_size != old["bytes"]:
            raise ValueError(f"Missing or changed pinned source: {spec['source_id']}")
        return old
    if offline:
        raise FileNotFoundError(f"No pinned source: {spec['source_id']}")
    if target.exists():
        raise ValueError(f"Unpinned file already exists: {target}")
    started = datetime.now(timezone.utc).isoformat()
    record = dict(spec)
    if spec["service"] == "local_reuse":
        parent = json.loads((ROOT / spec["origin_manifest"]).read_text())
        source = next(s for s in parent["sources"] if s["id"] == spec["origin_source_id"])
        origin = ROOT / "data/raw" / source["file"]
        if digest(origin) != source["sha256"]:
            raise ValueError("Previously pinned literature cache changed")
        target.write_bytes(origin.read_bytes())
        record.update(original_retrieved_at_utc=source["retrieved_at_utc"], original_local_cache_path=str(origin.relative_to(ROOT)), representation="Unmodified copy of the earlier pinned primary-paper HTML; current timestamp records local reuse, not a fresh download", source_url=source["url"])
    elif spec["service"] == "https":
        request = urllib.request.Request(spec["url"], headers={"User-Agent": "PSC-public-research/UBASH3A-sequence-audit"})
        with urllib.request.urlopen(request, timeout=45) as response:
            content = response.read(5_000_001)
            if len(content) > 5_000_000:
                raise ValueError("Literature response exceeded the bounded size")
            record.update(status_code=response.status, final_url=response.url, content_type=response.headers.get("content-type"))
        target.write_bytes(content)
        record["representation"] = "HTTP response body as received by urllib; not reformatted"
    else:
        service = spec["service"]
        client = clients[service]
        base_url = "https://rest.ensembl.org" if service == "ensembl" else "https://rest.uniprot.org"
        payload = {"base_url": base_url, "path": spec["path"], "params": spec["params"], "headers": {"Accept": "application/json", "Content-Type": "application/json"}, "response_format": "text", "save_raw": True, "raw_output_path": str(target), "timeout_sec": 45, "max_items": 10}
        run = subprocess.run([sys.executable, str(client)], input=json.dumps(payload), text=True, capture_output=True, timeout=60)
        result = json.loads(run.stdout)
        if not result.get("ok"):
            raise RuntimeError(json.dumps({"source_id": spec["source_id"], "client_result": result}))
        json.loads(target.read_text())
        record.update(url=base_url + "/" + spec["path"] + ("?" + urllib.parse.urlencode(spec["params"]) if spec["params"] else ""), request=payload, status_code=result["status_code"], client_sha256=digest(client), client_result=result)
        record["request"]["raw_output_path"] = str(target.relative_to(ROOT))
        record["client_result"]["raw_output_path"] = str(target.relative_to(ROOT))
        record["representation"] = "Skill client saves requests.response.text as UTF-8 without JSON reformatting; transport decompression/decoding is not a byte-level wire archive"
    if not target.is_file() or target.stat().st_size == 0 or target.stat().st_size > 5_000_000:
        raise ValueError(f"Invalid saved response size: {target}")
    record.update(local_cache_path=str(target.relative_to(ROOT)), retrieval_started_at_utc=started, retrieved_at_utc=datetime.now(timezone.utc).isoformat(), bytes=target.stat().st_size, sha256=digest(target), publisher_checksum=None)
    if old:
        if record["sha256"] != old["sha256"] or record["bytes"] != old["bytes"]:
            raise ValueError(f"Re-downloaded source changed: {spec['source_id']}; do not update the historical pin automatically")
        return old
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["metadata", "sequences", "literature", "all"], default="all")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--ensembl-client")
    parser.add_argument("--uniprot-client")
    args = parser.parse_args()
    plan = json.loads(PLAN.read_text())
    manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {"manifest_id": "ubash3a-consequence-sources-v1", "created_at_utc": datetime.now(timezone.utc).isoformat(), "plan_sha256": digest(PLAN), "sources": [], "failed_attempts": [], "retention": "Raw third-party content remains ignored; hashes and small derived outputs are versioned."}
    if manifest["plan_sha256"] != digest(PLAN):
        raise ValueError("The source manifest's analysis plan has changed")
    clients = {} if args.offline else {s: skill_client(s, getattr(args, s + "_client")) for s in ("ensembl", "uniprot")}
    RAW.mkdir(parents=True, exist_ok=True)
    stages = ["metadata", "sequences", "literature"] if args.stage == "all" else [args.stage]
    failed = False
    for stage in stages:
        existing = {s["source_id"]: s for s in manifest["sources"]}
        with ThreadPoolExecutor(max_workers=3) as pool:
            futures = [(spec, pool.submit(fetch_one, spec, existing, clients, args.offline)) for spec in requests_for(stage, plan)]
            for spec, future in futures:
                try:
                    row = future.result()
                    existing[spec["source_id"]] = row
                    print(json.dumps({"source_id": spec["source_id"], "status": "verified", "bytes": row["bytes"]}), flush=True)
                except Exception as exc:
                    failed = True
                    failure = {"source_id": spec["source_id"], "checked_at_utc": datetime.now(timezone.utc).isoformat(), "request_spec": spec, "error": str(exc)}
                    manifest["failed_attempts"].append(failure)
                    print(json.dumps(failure), flush=True)
        manifest["sources"] = sorted(existing.values(), key=lambda r: r["source_id"])
        if not args.offline:
            MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
        if failed:
            break
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
