"""Stream complete published BF archives; retain only preselected gene rows."""

import argparse
import gzip
import hashlib
import io
import json
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, gzip_bytes, save_json, sha


class CountedReader(io.RawIOBase):
    def __init__(self, response, cap, label):
        self.response=response;self.cap=cap;self.label=label
        self.digest=hashlib.sha256();self.count=0;self.prefix=bytearray();self.last=time.monotonic()
    def readable(self):return True
    def readinto(self, buffer):
        data=self.response.read(len(buffer))
        self.count+=len(data)
        if self.count>self.cap:raise ValueError('Published BF stream cap exceeded')
        self.digest.update(data)
        if len(self.prefix)<65536:self.prefix.extend(data[:65536-len(self.prefix)])
        buffer[:len(data)]=data
        if time.monotonic()-self.last>15:
            print(json.dumps({'dataset':self.label,'compressed_megabytes_read':round(self.count/1e6,1)}),flush=True)
            self.last=time.monotonic()
        return len(data)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline',action='store_true')
    args=parser.parse_args()
    plan_path=ROOT/'config/coloc-platform-plan.json'
    plan=json.loads(plan_path.read_text())
    for name,digest in plan['frozen_file_sha256'].items():
        if sha((ROOT/name).read_bytes())!=digest:raise ValueError('Platform amendment input changed')
    sources={s['id']:s for s in json.loads((ROOT/'config/coloc-sources.json').read_text())['sources']}
    path=ROOT/'config/coloc-published-bf-queries.json'
    run=json.loads(path.read_text()) if path.exists() else {'plan_sha256':sha(plan_path.read_bytes()),
          'reader_sha256':sha(Path(__file__).read_bytes()),'started_at_utc':datetime.now(timezone.utc).isoformat(),'queries':[]}
    if run['plan_sha256']!=sha(plan_path.read_bytes()) or run['reader_sha256']!=sha(Path(__file__).read_bytes()):raise ValueError('BF reader inputs changed')
    for request in plan['published_qtl_signal_sources']:
        ds=request['dataset_id']
        old=next((q for q in run['queries'] if q['dataset_id']==ds),None)
        if old:
            data=(ROOT/old['file']).read_bytes()
            if sha(data)!=old['sha256'] or sha(gzip.decompress(data))!=old['uncompressed_sha256']:raise ValueError('BF subset changed')
            continue
        if args.offline:raise ValueError('Missing complete published BF query')
        source=sources['coloc-'+ds+'-lbf-header']
        queried_at=datetime.now(timezone.utc).isoformat()
        req=urllib.request.Request(source['url'],headers={'User-Agent':'PSC-public-research/coloc-signal-v01'})
        with urllib.request.urlopen(req,timeout=60) as response:
            size=int(response.headers['Content-Length'])
            if response.status!=200 or size>request['maximum_compressed_bytes']:raise ValueError('Unexpected BF archive size or response')
            reader=CountedReader(response,request['maximum_compressed_bytes'],ds)
            selected=[];scanned=0
            with gzip.GzipFile(fileobj=io.BufferedReader(reader,buffer_size=1024*1024),mode='rb') as stream:
                header=stream.readline()
                if not header.startswith(b'molecular_trait_id\tregion\tvariant\t'):raise ValueError('BF header differs')
                prefix=request['gene_id'].encode()+b'\t'
                for line in stream:
                    scanned+=1
                    if line.startswith(prefix):selected.append(line)
            if reader.count!=size or sha(reader.prefix)!=source['sha256']:raise ValueError('Incomplete stream or changed pinned BF prefix')
        if not selected:raise ValueError('Published signal has no full BF rows')
        plain=header+b''.join(selected)
        data=gzip_bytes(plain)
        filename='data/derived/coloc-query-rows/'+ds+'-PFKFB3-published-bfs.tsv.gz'
        (ROOT/filename).write_bytes(data)
        record={**request,'url':source['url'],'queried_at_utc':queried_at,'finished_at_utc':datetime.now(timezone.utc).isoformat(),
                'compressed_archive_bytes':size,'compressed_archive_sha256':reader.digest.hexdigest(),'archive_eof_verified':True,
                'all_archive_rows_scanned':scanned,'retained_gene_rows':len(selected),'file':filename,
                'sha256':sha(data),'uncompressed_sha256':sha(plain),'bytes':len(data)}
        run['queries'].append(record);save_json(path,run)
        print(json.dumps({'dataset':ds,'status':'complete-archive-scanned','retained_gene_rows':len(selected),'compressed_bytes':size}),flush=True)
        time.sleep(2)
    print(json.dumps({'complete_published_signal_datasets':len(run['queries']),'status':'verified-cache' if args.offline else 'complete'}))


if __name__=='__main__':main()
