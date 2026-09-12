"""Run/verify one pinned R analysis stage, recording code, inputs and outputs."""

import argparse
import gzip
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, gzip_bytes, save_json, sha, verify_plan


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['baseline','fit','compare'])
    parser.add_argument('--dataset',choices=['GWAS','QTD000031'])
    parser.add_argument('--gene',choices=['PFKFB3','BCL2L11'])
    parser.add_argument('--offline',action='store_true')
    parser.add_argument('--rerun',action='store_true')
    parser.add_argument('--rscript',type=Path,default=ROOT/'work/coloc-runtime/r/bin/Rscript')
    args=parser.parse_args()
    if args.offline and args.rerun: parser.error('Choose offline or rerun')
    if args.stage=='fit' and (not args.dataset or not args.gene): parser.error('Fit requires dataset and gene')
    verify_plan()
    lock_path=ROOT/'config/coloc-analysis-inputs.json'
    lock=json.loads(lock_path.read_text())
    for name,digest in lock['file_sha256'].items():
        if sha((ROOT/name).read_bytes())!=digest: raise ValueError('Frozen analysis input changed: '+name)
    for path in (ROOT/'data/raw/coloc-ld').glob('*.json'):
        record=json.loads(path.read_text())
        if sha(path.with_suffix('.f64').read_bytes())!=record['matrix_sha256']: raise ValueError('LD matrix changed')
    identifier=args.stage+('-'+args.dataset+'-'+args.gene if args.stage=='fit' else '')
    provenance_path=ROOT/'reports/coloc-run-provenance.json'
    provenance=json.loads(provenance_path.read_text()) if provenance_path.exists() else {
        'analysis_input_lock_sha256':sha(lock_path.read_bytes()),'executor_sha256':sha(Path(__file__).read_bytes()),'jobs':[]}
    if provenance['analysis_input_lock_sha256']!=sha(lock_path.read_bytes()) or provenance['executor_sha256']!=sha(Path(__file__).read_bytes()):
        raise ValueError('Analysis lock or executor changed')
    previous=next((r for r in provenance['jobs'] if r['id']==identifier and r['status']=='complete'),None)
    if previous:
        for name,digest in previous['output_sha256'].items():
            if sha((ROOT/name).read_bytes())!=digest: raise ValueError('Saved analysis output changed')
        if not args.rerun:
            print(json.dumps({'id':identifier,'status':'verified-cache'}));return
    elif args.offline:
        raise ValueError('No completed analysis stage '+identifier)
    if args.stage=='baseline':
        names=['reports/coloc-baseline.csv','reports/coloc-baseline-variant-posteriors.csv.gz','reports/coloc-baseline-warnings.json']
    elif args.stage=='compare':
        names=['reports/coloc-reference-comparison.csv','reports/coloc-reference-variant-posteriors.csv.gz']
    else:
        stem=args.dataset+'-'+args.gene
        names=['reports/coloc-LD-diagnostics-'+stem+'.csv','reports/coloc-reference-pips-'+stem+'.csv',
               'reports/coloc-reference-fit-'+stem+'.json','data/raw/coloc-ld/'+stem+'.rds']
    command=[str(args.rscript),'--vanilla','scripts/run_coloc.R',args.stage]
    if args.stage=='fit': command += [args.dataset,args.gene]
    record={'id':identifier,'started_at_utc':datetime.now(timezone.utc).isoformat(),'R':'4.4.3','coloc':'5.2.3','susieR':'0.14.2'}
    log=ROOT/'data/raw'/('coloc-'+identifier+'.log')
    env={**os.environ,'OPENBLAS_NUM_THREADS':'4','OMP_NUM_THREADS':'4'}
    with log.open('w') as stream:
        process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in process.stdout:
            stream.write(line);stream.flush();print(line,end='',flush=True)
        code=process.wait()
    if code:
        raise RuntimeError('R analysis failed; source log retained: '+str(log.relative_to(ROOT)))
    for name in names:
        path=ROOT/name
        if path.exists() and name.endswith('.csv.gz'):
            path.write_bytes(gzip_bytes(gzip.decompress(path.read_bytes())))
    record.update(status='complete',finished_at_utc=datetime.now(timezone.utc).isoformat(),
                  log_file=str(log.relative_to(ROOT)),log_sha256=sha(log.read_bytes()),
                  output_sha256={name:sha((ROOT/name).read_bytes()) for name in names if (ROOT/name).exists()})
    if previous:
        if previous['output_sha256']!=record['output_sha256']: raise ValueError('Replay was not byte-identical')
        print(json.dumps({'id':identifier,'status':'byte-identical-replay'}));return
    provenance['jobs'].append(record);save_json(provenance_path,provenance)
    print(json.dumps({'id':identifier,'status':'complete','outputs':len(record['output_sha256'])}),flush=True)


if __name__=='__main__':
    main()
