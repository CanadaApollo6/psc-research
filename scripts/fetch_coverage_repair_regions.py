"""Bounded public reference and all-gene QTL queries for the seven fixed targets."""

import argparse
import gzip
import hashlib
import json
import time
import zlib
from datetime import datetime,timezone
from pathlib import Path

import pysam

from coloc_common import ROOT,gzip_bytes,save_json,sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['reference','catalogue','alternative'])
    parser.add_argument('--offline',action='store_true')
    args=parser.parse_args()
    plan_path=ROOT/'config/coverage-repair-plan.json';plan=json.loads(plan_path.read_text())
    for name,digest in plan['frozen_file_sha256'].items():
        if sha((ROOT/name).read_bytes())!=digest:raise ValueError('Frozen input changed: '+name)
    sources={s['id']:s for s in json.loads((ROOT/'config/coloc-sources.json').read_text())['sources']}
    extras={s['id']:s for s in json.loads((ROOT/'config/coverage-repair-sources.json').read_text())['sources']}
    path=ROOT/f'config/coverage-repair-{args.stage}-queries.json'
    pins={'plan_sha256':sha(plan_path.read_bytes()),'reader_sha256':sha(Path(__file__).read_bytes()),'pysam_version':pysam.__version__}
    if args.stage=='alternative':
        alt_path=ROOT/'config/coverage-repair-alternative-plan.json';alternative=json.loads(alt_path.read_text())
        pins['alternative_plan_sha256']=sha(alt_path.read_bytes())
    run=json.loads(path.read_text()) if path.exists() else {**pins,'queries':[]}
    if any(run.get(k)!=v for k,v in pins.items()):raise ValueError('Query provenance changed')
    folder=ROOT/'data/derived/coverage-repair-queries';folder.mkdir(exist_ok=True)

    def cached(identifier):
        old=next((q for q in run['queries'] if q['id']==identifier),None)
        if old:
            if sha((ROOT/old['file']).read_bytes())!=old['sha256']:raise ValueError('Changed query cache')
            return True
        if args.offline:raise ValueError('Missing query '+identifier)
        return False

    def save(identifier,plain,file,details):
        data=gzip_bytes(plain);target=ROOT/file
        if target.exists() and target.read_bytes()!=data:raise ValueError('Unpinned output differs')
        target.write_bytes(data)
        run['queries'].append({'id':identifier,'file':file,'sha256':sha(data),'uncompressed_sha256':sha(plain),'bytes':len(data),
                               'finished_at_utc':datetime.now(timezone.utc).isoformat(),**details})
        save_json(path,run);print(json.dumps({'id':identifier,'status':'complete',**{k:v for k,v in details.items() if k in ['rows','scanned_rows']}}),flush=True)

    if args.stage=='reference':
        for gene,chrom in [('PFKFB3','10'),('BCL2L11','2')]:
            identifier='reference-'+gene
            if cached(identifier):continue
            ref=extras['repair-reference-'+gene+'-GRCh37'];start,end=ref['start_1based']-1,ref['end_1based']
            index=sources['match-kgp-chr'+chrom+'-index'];url=index['url'].removesuffix('.tbi')
            if sha((ROOT/'data/raw'/index['file']).read_bytes())!=index['sha256']:raise ValueError('Reference index changed')
            with pysam.TabixFile(url,index=str(ROOT/'data/raw'/index['file'])) as tab:
                header=list(tab.header);rows=list(tab.fetch(chrom,start,end))
            if sum(len(r)+1 for r in rows)>plan['network_bounds']['nominal_uncompressed_bytes_per_query']:raise ValueError('Reference cap exceeded')
            plain=('\n'.join(header+rows)+'\n').encode()
            save(identifier,plain,'data/raw/repair-reference-'+gene+'.vcf.gz',{'source_url':url,'source_index_sha256':index['sha256'],'build':'GRCh37','start_0based':start,'end_exclusive':end,'rows':len(rows),'sample_bearing':True})
            time.sleep(2)
    else:
        datasets=plan['datasets'] if args.stage=='catalogue' else alternative['datasets']
        targets={v['position_grch38'] for v in plan['bcl2l11_variants']}
        start,end=(min(targets)-1001,max(targets)+1000) if args.stage=='catalogue' else (alternative['start_grch38_0based'],alternative['end_grch38_exclusive'])
        for ds in datasets:
            identifier=ds['dataset_id']
            if cached(identifier):continue
            index=(sources[identifier+'-index'] if args.stage=='catalogue' else extras['repair-'+identifier+'-index'])
            head=(sources[identifier+'-header'] if args.stage=='catalogue' else extras['repair-'+identifier+'-header'])
            for source in [index,head]:
                if sha((ROOT/'data/raw'/source['file']).read_bytes())!=source['sha256']:raise ValueError('QTL source changed')
            header=zlib.decompressobj(31).decompress((ROOT/'data/raw'/head['file']).read_bytes()).decode().splitlines()[0]
            columns=header.split('\t');posidx=columns.index('position');geneidx=columns.index('gene_id')
            rows=[];count=0;size=0;digest=hashlib.sha256();seen=set()
            with pysam.TabixFile(ds['url'],index=str(ROOT/'data/raw'/index['file'])) as tab:
                for line in tab.fetch('2',start,end):
                    count+=1;size+=len(line)+1;digest.update((line+'\n').encode())
                    if count>plan['network_bounds']['nominal_rows_per_query'] or size>plan['network_bounds']['nominal_uncompressed_bytes_per_query']:raise ValueError('Nominal transfer cap exceeded')
                    fields=line.split('\t')
                    if len(fields)!=len(columns):raise ValueError('Malformed nominal record')
                    position=int(fields[posidx]);seen.add(position)
                    keep=position in targets if args.stage=='catalogue' else fields[geneidx].split('.')[0]==alternative['gene_id']
                    if keep:rows.append(line)
            plain=(header+'\n'+'\n'.join(rows)+ ('\n' if rows else '')).encode()
            save(identifier,plain,str((folder/(identifier+'-'+args.stage+'.tsv.gz')).relative_to(ROOT)),{'source_url':ds['url'],'index_sha256':index['sha256'],'header_sha256':head['sha256'],'build':'GRCh38','start_0based':start,'end_exclusive':end,'rows':len(rows),'scanned_rows':count,'scanned_bytes':size,'scanned_sha256':digest.hexdigest(),'observed_positions':len(seen),'min_observed_position':min(seen) if seen else None,'max_observed_position':max(seen) if seen else None})
            time.sleep(2)
    expected=2 if args.stage=='reference' else len(datasets)
    if len(run['queries'])!=expected:raise ValueError('Incomplete query family')
    print(json.dumps({'stage':args.stage,'complete_queries':len(run['queries'])}),flush=True)


if __name__=='__main__':main()
