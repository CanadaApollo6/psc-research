"""Bounded public-source acquisition for the original IBDverse model audit."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/'data/raw/ibdverse-model-audit'
MANIFEST = ROOT/'config/ibdverse-model-audit-sources.json'
PLAN = ROOT/'config/ibdverse-model-audit-plan.json'
BIO_SKILL = Path('/home/riels/.codex/plugins/cache/openai-curated-remote/life-science-research/1.0.3/skills/biostudies-arrayexpress-skill/scripts/rest_request.py')
PMC_SKILL = Path('/home/riels/.codex/plugins/cache/openai-curated-remote/life-sciences-literature/0.1.5/skills/ncbi-pmc-skill/scripts/ncbi_pmc.py')


def now(): return datetime.now(timezone.utc).isoformat()


def sha(path):
    with Path(path).open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def save(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def freeze():
    path = ROOT/'config/ibdverse-model-audit-freeze.json'
    if path.exists():
        record=json.loads(path.read_text())
        if record['plan_sha256'] != sha(PLAN): raise ValueError('Frozen audit plan changed')
        return record
    files=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
    record={'frozen_at_utc':now(),'plan_sha256':sha(PLAN),
            'baseline_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'prior_tracked_sha256':{name:sha(ROOT/name) for name in files},
            'prior_knowledge':'Completed source recovery, df consistency inference and regulatory mechanism comparison; existing main code tree inventory and four small target archive members known.'}
    save(path,record)
    return record


def prior(label):
    record=json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {'sources':[]}
    previous=next((s for s in record['sources'] if s['id']==label),None)
    if previous and previous.get('file'):
        if sha(ROOT/previous['file']) != previous['sha256']: raise ValueError('Cached source changed: '+label)
    return record,previous


def append(record,entry):
    record['sources'].append(entry)
    save(MANIFEST,record)
    print(json.dumps({k:entry.get(k) for k in ['id','status','ok','bytes','error_code']}),flush=True)
    return entry


def fetch(label,url,limit=25000000,headers=None,method='GET'):
    if 'ebi.ac.uk/biostudies/api/' in url: raise ValueError('Use the BioStudies skill wrapper')
    freeze(); record,previous=prior(label)
    if previous:
        if previous['url']!=url or previous.get('request_headers',{})!=(headers or {}): raise ValueError('Source ID reused')
        return previous
    # Archive-member transfers use their separate bounded reader. Include them
    # when checking the stage-wide ceiling for ordinary sources and HDF5 ranges.
    acquired = sum(s.get('bytes',0) or 0 for s in record['sources'])
    for auxiliary in ['ibdverse-model-audit-member-queries.json', 'ibdverse-model-audit-coloc-queries.json']:
        auxiliary_path = ROOT / 'config' / auxiliary
        if auxiliary_path.exists():
            acquired += sum(r['bytes'] for q in json.loads(auxiliary_path.read_text()).get('queries',[]) for r in q['ranges'])
    budget = json.loads(PLAN.read_text())['retrieval']['initial_total_new_download_budget_bytes']
    remaining = budget - acquired
    if remaining <= 0: raise ValueError('Stage-wide acquisition budget exhausted')
    limit = min(limit, remaining)
    entry={'id':label,'url':url,'method':method,'request_headers':headers or {},'maximum_bytes':limit,'requested_at_utc':now()}
    path=RAW/(label+'.data')
    try:
        with requests.request(method,url,headers={'User-Agent':'PSC-public-research/IBDverse-model-audit',**(headers or {})},timeout=45,stream=True) as response:
            entry.update(status=response.status_code,resolved_url=response.url,
                         response_headers={k:response.headers[k] for k in ['Content-Type','Content-Length','Content-Range','ETag','Last-Modified'] if k in response.headers})
            response.raise_for_status()
            if headers and 'Range' in headers and response.status_code!=206: raise ValueError('Requested range was not honored')
            content=bytearray()
            for chunk in response.iter_content(65536):
                content.extend(chunk)
                if len(content)>limit: raise ValueError('Bounded source size exceeded')
            path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(content)
            entry.update(ok=True,file=str(path.relative_to(ROOT)),bytes=len(content),sha256=sha(path))
    except Exception as exc:
        entry.update(ok=False,error_code=type(exc).__name__)
    entry['retrieved_at_utc']=now()
    return append(record,entry)


def skill(label,provider,payload):
    freeze();record,previous=prior(label)
    path=RAW/(label+'.json')
    if provider=='biostudies':
        script=BIO_SKILL
        request={'base_url':'https://www.ebi.ac.uk/biostudies/api/v1','headers':{'Accept':'application/json'},
                 'timeout_sec':45,'save_raw':True,'raw_output_path':str(path),'max_items':2,'max_depth':2,**payload}
    elif provider=='pmc':
        script=PMC_SKILL
        # Retain the compact source/checked-source contract as a local receipt.
        request={'timeout_sec':45,'max_items':2,**payload}
    else: raise ValueError('Unknown source provider')
    request_path=ROOT/'config/ibdverse-model-audit-requests'/(label+'.json')
    if previous:
        if previous['request']!=request: raise ValueError('Skill request changed')
        return previous
    save(request_path,request)
    entry={'id':label,'provider':provider,'request':request,'requested_at_utc':now(),'skill_script_sha256':sha(script)}
    result=subprocess.run([sys.executable,str(script)],input=json.dumps(request),text=True,capture_output=True)
    try: data=json.loads(result.stdout)
    except json.JSONDecodeError: data={'ok':False,'error':{'code':'non_json_skill_output'}}
    receipt=RAW/(label+'-receipt.json');save(receipt,data)
    entry.update(ok=data.get('ok',False),status=data.get('status_code'),error_code=data.get('error',{}).get('code'),
                 receipt=str(receipt.relative_to(ROOT)),receipt_sha256=sha(receipt),retrieved_at_utc=now())
    if provider=='pmc': path=receipt
    if path.exists(): entry.update(file=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path))
    return append(record,entry)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['freeze','initial'])
    args=parser.parse_args();freeze()
    if args.stage=='initial':
        fetch('ibdverse-code-current','https://api.github.com/repos/andersonlab/IBDVerse-sc-eQTL-code/commits/main')
        for accession in ['S-BSST2922','E-MTAB-16986']:
            skill('biostudies-'+accession,'biostudies',{'path':'studies/'+accession})
            skill('biostudies-'+accession+'-info','biostudies',{'path':'studies/'+accession+'/info'})
        skill('pmc-IBDverse','pmc',{'params':{'id':'PMC13441992'}})
