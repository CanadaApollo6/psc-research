#!/usr/bin/env python3
"""Read-only, secret-safe PSC model capability check. Never runs predictions.

Use the existing project interpreter. Optional network calls are limited to one
AlphaGenome metadata request and one public AFDB record. No installs, login,
provisioning, auth-response dumps, or sequence/structure requests are performed.
"""
from __future__ import annotations

import argparse
import dataclasses
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import sys
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ("alphagenome", "alphagenome-research", "alphafold", "alphafold3",
            "colabfold", "jax", "jaxlib", "torch")
ENV_NAMES = ("ALPHAGENOME_API_KEY", "ALPHAFOLD_API_KEY", "ALPHAFOLD3_API_KEY",
             "AF3_API_KEY", "AF3_MODEL_DIR", "ALPHAFOLD_DATA_DIR",
             "ALPHAFOLD_MODEL_DIR", "CUDA_VISIBLE_DEVICES")


def now():
    return datetime.now(timezone.utc).isoformat()


def error_code(exc):
    # Never include arbitrary exception text, headers, or auth responses.
    method = getattr(exc, "code", None)
    if callable(method):
        return method().name
    if isinstance(exc, HTTPError):
        return "HTTP_" + str(exc.code)
    return type(exc).__name__


def local_inventory():
    packages = {}
    for name in PACKAGES:
        try:
            dist = importlib.metadata.distribution(name)
            row = {"installed": True, "version": dist.version}
            # direct_url may contain a private URL; retain only commit hash.
            origin = json.loads(dist.read_text("direct_url.json") or "{}")
            commit = origin.get("vcs_info", {}).get("commit_id", "")
            if re.fullmatch(r"[0-9a-f]{40,64}", commit):
                row["commit_id"] = commit
            packages[name] = row
        except importlib.metadata.PackageNotFoundError:
            packages[name] = {"installed": False}
    modules = {}
    for name in ("alphagenome", "alphafold", "alphafold3", "colabfold"):
        modules[name] = importlib.util.find_spec(name) is not None
    relevant = set(ENV_NAMES) | {k for k in os.environ if
        re.match(r"^(ALPHAGENOME|ALPHAFOLD|AFDB|AF[23]_)", k)}
    return {
        "python": sys.version.split()[0],
        "interpreter": ".venv/bin/python" if Path(sys.prefix) == ROOT / ".venv" else "other interpreter",
        "packages": packages,
        "module_discoverable": modules,
        "cli_on_path": {name: shutil.which(name) is not None for name in (
            "alphafold", "alphafold3", "run_alphafold.py", "colabfold_batch", "docker", "nvidia-smi")},
        "environment_key_presence": {key: bool(os.environ.get(key, "")) for key in sorted(relevant)},
        "project_dotenv_present": (ROOT / ".env").is_file(),
        "documented_alphagenome_key_file_present": (ROOT / "data/private/alphagenome_api_key").is_file(),
        "scope": "Current project interpreter, PATH, relevant environment names, and documented key location only. No home/browser/credential-store scan.",
    }


def alphagenome_metadata(enabled):
    from alphagenome_access import key_is_configured
    row = {"credential_configured": key_is_configured(), "request_count": 0,
           "inference_request_count": 0, "sdk_import": "not attempted"}
    if not enabled:
        row["status"] = "disabled"
        return row
    try:
        import alphagenome
        from alphagenome.models import dna_client, variant_scorers
    except Exception as exc:
        row.update(status="sdk_unavailable", sdk_import="failed", import_error_code=error_code(exc))
        return row
    row.update(sdk_import="success", client_version=alphagenome.__version__,
               model_selections=[x.name for x in dna_client.ModelVersion],
               output_types=[x.name for x in dna_client.OutputType],
               recommended_scorers=sorted(variant_scorers.RECOMMENDED_VARIANT_SCORERS))
    if not row["credential_configured"]:
        row["status"] = "credential_missing; no RPC attempted"
        return row
    try:
        from alphagenome_access import load_key
        client = dna_client.create(load_key(), model_version=dna_client.ModelVersion.ALL_FOLDS, timeout=20)
        row["request_count"] = 1
        metadata = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
        row.update(status="metadata_success", selection="ALL_FOLDS", tracks=[])
        for field in dataclasses.fields(metadata):
            frame = getattr(metadata, field.name)
            if frame is not None:
                row["tracks"].append({"output_type": field.name, "track_count": len(frame),
                                      "sha256_csv": hashlib.sha256(frame.to_csv(index=False).encode()).hexdigest()})
        row["boundary"] = "Metadata success is not an inference test or an exact server-build/quota receipt."
    except Exception as exc:
        row.update(status="metadata_failed", error_code=error_code(exc))
    return row


def public_afdb(enabled):
    row = {"url": "https://alphafold.ebi.ac.uk/api/prediction/P57075", "request_count": 0,
           "authentication": "none", "new_structure_inference": False}
    if not enabled:
        row["status"] = "disabled"
        return row
    try:
        row["request_count"] = 1
        request = Request(row["url"], headers={"Accept": "application/json", "User-Agent": "PSC-research-capability-check/1"})
        with urlopen(request, timeout=20) as response:
            payload = response.read(1_048_577)
            row["http_status"] = response.status
        if len(payload) > 1_048_576:
            raise ValueError("Response exceeds cap")
        records = json.loads(payload)
        if not isinstance(records, list):
            raise TypeError("Expected public list")
        row.update(status="success", received_bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(), records=[])
        allowed = ("entryId", "uniprotAccession", "gene", "uniprotStart", "uniprotEnd", "latestVersion", "allVersions",
                   "modelCreatedDate", "sequenceStart", "sequenceEnd", "sequenceVersionDate")
        for record in records:
            item = {k: record[k] for k in allowed if k in record}
            if "uniprotSequence" in record:
                seq = record["uniprotSequence"]
                item.update(sequence_length=len(seq), sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())
            row["records"].append(item)
    except Exception as exc:
        row.update(status="failed", error_code=error_code(exc))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alphagenome-metadata", action="store_true")
    parser.add_argument("--public-afdb", action="store_true")
    parser.add_argument("--output", type=Path, default=ROOT / "work/model-capabilities/capabilities.json")
    args = parser.parse_args()
    output = args.output.resolve()
    work = (ROOT / "work/model-capabilities").resolve()
    if not output.is_relative_to(work):
        parser.error("Output must stay inside ignored work/model-capabilities/")
    if output.exists():
        parser.error("Refusing to overwrite an earlier capability receipt")
    report = {"started_at_utc": now(), "local": local_inventory(),
              "alphagenome": alphagenome_metadata(args.alphagenome_metadata),
              "afdb": public_afdb(args.public_afdb),
              "new_predictions": 0, "finished_at_utc": now()}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
