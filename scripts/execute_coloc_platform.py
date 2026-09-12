"""Execute and fingerprint the separate, source-supported platform amendment."""

import argparse
import gzip
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, gzip_bytes, save_json, sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['baseline','fit','compare'])
    parser.add_argument('--gene',choices=['PFKFB3','BCL2L11'])
    parser.add_argument('--offline',action='store_true')
    parser.add_argument('--rerun',action='store_true')
    parser.add_argument('--rscript',type=Path,default=ROOT/'work/coloc-runtime/r/bin/Rscript')
    args=parser.parse_args()
    if args.stage=='fit' and not args.gene:parser.error('Fit needs a gene')
    if args.offline and args.rerun:parser.error('Choose offline or rerun')
    lock_path=ROOT/'config/coloc-platform-analysis-inputs.json'
    lock=json.loads(lock_path.read_text())
    for name,digest in lock['file_sha256'].items():
        if sha((ROOT/name).read_bytes())!=digest:raise ValueError('Amendment input changed: '+name)
    for record in json.loads((ROOT/'reports/coloc-platform-input-audit.json').read_text())['datasets']:
        if sha((ROOT/record['matrix_file']).read_bytes())!=record['matrix_sha256']:raise ValueError('Amended LD changed')
    path=ROOT/'reports/coloc-platform-run-provenance.json'
    run=json.loads(path.read_text()) if path.exists() else {'input_lock_sha256':sha(lock_path.read_bytes()),'executor_sha256':sha(Path(__file__).read_bytes()),'jobs':[]}
    if run['input_lock_sha256']!=sha(lock_path.read_bytes()) or run['executor_sha256']!=sha(Path(__file__).read_bytes()):raise ValueError('Amendment executor changed')
    identifier=args.stage+('-'+args.gene if args.stage=='fit' else '')
    old=next((j for j in run['jobs'] if j['id']==identifier),None)
    if old:
        for name,digest in old['output_sha256'].items():
            if sha((ROOT/name).read_bytes())!=digest:raise ValueError('Amendment output changed')
        if not args.rerun:
            print(json.dumps({'id':identifier,'status':'verified-cache'}));return
    elif args.offline:raise ValueError('No completed amendment stage')
    extra={}
    if args.stage=='baseline':
        names=['reports/coloc-platform-baseline.csv','reports/coloc-platform-variant-posteriors.csv.gz','reports/coloc-platform-warnings.json']
    elif args.stage=='compare':
        bf_path=ROOT/'config/coloc-published-bf-queries.json';bf=json.loads(bf_path.read_text())
        if len(bf['queries'])!=2:raise ValueError('Complete published BF queries required')
        for q in bf['queries']:
            if not q['archive_eof_verified'] or sha((ROOT/q['file']).read_bytes())!=q['sha256']:raise ValueError('Published BF input changed or incomplete')
        extra={'published_bf_manifest_sha256':sha(bf_path.read_bytes())}
        names=['reports/coloc-published-signal-comparison.csv','reports/coloc-published-signal-overlap.csv','reports/coloc-published-signal-variant-posteriors.csv.gz']
    else:
        names=['reports/coloc-platform-LD-diagnostics-'+args.gene+'.csv']
        for n in [11386,3504,14890]:
            names += ['reports/coloc-platform-fit-'+args.gene+'-N'+str(n)+'.json',
                      'reports/coloc-platform-pips-'+args.gene+'-N'+str(n)+'.csv',
                      'data/raw/coloc-platform-ld/GWAS-'+args.gene+'-N'+str(n)+'.rds']
    command=[str(args.rscript),'--vanilla','scripts/run_coloc_platform.R',args.stage]+([args.gene] if args.stage=='fit' else [])
    log=ROOT/'data/raw'/('coloc-platform-'+identifier+'.log')
    record={'id':identifier,'started_at_utc':datetime.now(timezone.utc).isoformat(),**extra}
    with log.open('w') as stream:
        proc=subprocess.Popen(command,cwd=ROOT,env={**os.environ,'OPENBLAS_NUM_THREADS':'4','OMP_NUM_THREADS':'4'},stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in proc.stdout:
            stream.write(line);stream.flush();print(line,end='',flush=True)
        status=proc.wait()
    if status:raise RuntimeError('Amendment R stage failed; log retained')
    for name in names:
        target=ROOT/name
        if target.exists() and name.endswith('.csv.gz'):target.write_bytes(gzip_bytes(gzip.decompress(target.read_bytes())))
    record.update(status='complete',finished_at_utc=datetime.now(timezone.utc).isoformat(),log_file=str(log.relative_to(ROOT)),
                  log_sha256=sha(log.read_bytes()),output_sha256={name:sha((ROOT/name).read_bytes()) for name in names if (ROOT/name).exists()})
    if old:
        if old['output_sha256']!=record['output_sha256']:raise ValueError('Amendment replay differed')
        print(json.dumps({'id':identifier,'status':'byte-identical-replay'}));return
    run['jobs'].append(record);save_json(path,run)
    print(json.dumps({'id':identifier,'status':'complete'}),flush=True)


if __name__=='__main__':main()
