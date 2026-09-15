"""Use the mandatory service client; add response-header provenance without changing it."""
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
request_path = Path(sys.argv[1])
config = json.loads(request_path.read_text())
service = config.pop("_service")
client_path = Path("/home/riels/.codex/plugins/cache/openai-curated-remote/life-science-research/1.0.3/skills") / (service + "-skill") / "scripts/rest_request.py"
spec = importlib.util.spec_from_file_location("service_client", client_path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
responses = []
original = module.requests.Session.request
def capture(session, *args, **kwargs):
    result = original(session, *args, **kwargs)
    responses.append({
        "url": result.url, "status": result.status_code,
        "headers": {k:v for k,v in result.headers.items() if k.lower() in {
            "content-type", "content-length", "content-encoding", "etag", "last-modified",
            "x-total-results", "x-uniprot-release", "x-uniprot-release-date", "link"
        }},
    })
    return result
module.requests.Session.request = capture
started = datetime.datetime.now(datetime.timezone.utc).isoformat()
try:
    result = module.execute(config)
finally:
    module.requests.Session.request = original
record = {
    "started_at_utc": started,
    "completed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "request": config,
    "client_path": str(client_path),
    "client_sha256": hashlib.sha256(client_path.read_bytes()).hexdigest(),
    "response_headers_captured": responses,
    "result": result,
    "snapshot_format": "FASTA/text decoded as UTF-8 or JSON parsed and reserialized by the unchanged skill client",
}
if result.get("ok") and result.get("raw_output_path"):
    path = Path(result["raw_output_path"])
    record["local_path"] = str(path.relative_to(ROOT)) if path.is_absolute() else str(path)
    record["size_bytes"] = path.stat().st_size
    record["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
record_path = request_path.with_suffix(".receipt.json")
record_path.write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps({"ok": result.get("ok"), "receipt":str(record_path), "size_bytes":record.get("size_bytes"),
                  "response":responses, "error": result.get("error")}, indent=2))
sys.exit(0 if result.get("ok") else 1)
