import datetime, hashlib, importlib.util, json, re, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/"work/ubash3a-proteomics-specificity"
RAW=ROOT/"data/raw/ubash3a-proteomics-specificity"
CLIENT=Path("/home/riels/.codex/plugins/cache/openai-curated-remote/life-science-research/1.0.3/skills/uniprot-skill/scripts/rest_request.py")
spec=importlib.util.spec_from_file_location("uniprot_client",CLIENT)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
original=module.requests.Session.request
observed=[]
def capture(session,*args,**kwargs):
    response=original(session,*args,**kwargs)
    observed.append({"url":response.url,"status":response.status_code,"headers":{
        k.lower():v for k,v in response.headers.items() if k.lower() in {
            "content-type","content-encoding","x-uniprot-release","x-uniprot-release-date","x-total-results","link"}}})
    return response
module.requests.Session.request=capture
first=json.loads((WORK/"uniprot-human-page000.request.receipt.json").read_text())
first_response=first["response_headers_captured"][0]
headers={k.lower():v for k,v in first_response["headers"].items()}
expected=int(headers["x-total-results"]);release=headers["x-uniprot-release"]
assert expected==232837 and release=="2026_03"
pages=[{"page":0,"local_path":first["local_path"],"sha256":first["sha256"],
        "size_bytes":first["size_bytes"],"response":{**first_response,"headers":headers}}]
seen={first_response["url"]}
def next_url(headers):
    value=headers.get("link","")
    if not value:return None
    found=re.findall(r'<([^>]+)>;\s*rel="next"',value)
    assert len(found)==1,value
    return found[0]
url=next_url(headers)
while url:
    assert url not in seen and len(pages)<500
    seen.add(url);index=len(pages)
    path=RAW/f"uniprot-human-page{index:03d}.fasta"
    request={"base_url":"https://rest.uniprot.org","path":url,"response_format":"text",
             "save_raw":True,"raw_output_path":str(path),"timeout_sec":60,"max_items":1}
    observed.clear();started=datetime.datetime.now(datetime.timezone.utc).isoformat()
    result=module.execute(request)
    if not result.get("ok"):
        (WORK/"uniprot-pagination-failure.json").write_text(json.dumps({"page":index,"url":url,"result":result},indent=2)+"\n")
        raise RuntimeError(result)
    response=observed[-1];h=response["headers"]
    assert h["x-uniprot-release"]==release and int(h["x-total-results"])==expected
    data=path.read_bytes();count=sum(line.startswith(b">") for line in data.splitlines())
    assert 1<=count<=500
    pages.append({"page":index,"local_path":str(path.relative_to(ROOT)),"sha256":hashlib.sha256(data).hexdigest(),
                  "size_bytes":len(data),"records":count,"started_at_utc":started,"response":response})
    (WORK/"uniprot-pagination-progress.json").write_text(json.dumps({"pages":len(pages),"last_page":index,"records_so_far":500+sum(x.get("records",0) for x in pages[1:])})+"\n")
    if index%25==0:print(f"UniProt pages received: {len(pages)} / 466",flush=True)
    url=next_url(h)
output=RAW/"uniprot-human-2026_03-paginated-with-isoforms.fasta"
count=0
with output.open("wb") as f:
    for page in pages:
        data=(ROOT/page["local_path"]).read_bytes()
        assert hashlib.sha256(data).hexdigest()==page["sha256"]
        count+=sum(line.startswith(b">") for line in data.splitlines())
        f.write(data)
assert count==expected
receipt={"completed_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
         "client_path":str(CLIENT),"client_sha256":hashlib.sha256(CLIENT.read_bytes()).hexdigest(),
         "release":release,"release_date":headers["x-uniprot-release-date"],"expected_records_including_isoforms":expected,
         "records":count,"page_count":len(pages),"pages":pages,"local_path":str(output.relative_to(ROOT)),
         "size_bytes":output.stat().st_size,"sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
         "snapshot_format":"FASTA text saved by unchanged UniProt skill client, pages concatenated in returned order",
         "complete":True}
(WORK/"uniprot-pagination-completion.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({k:receipt[k] for k in ["complete","records","page_count","size_bytes","sha256"]}),flush=True)
