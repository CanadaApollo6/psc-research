"""Verify portable coloc artifacts and report separately on ignored local caches."""

import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def verify(root=ROOT, require_cache=False):
    pins={}

    def add(name,digest):
        if name in pins and pins[name]!=digest:
            raise ValueError('Conflicting historical hashes: '+name)
        pins[name]=digest

    for name in ['coloc-plan','coloc-platform-plan','coloc-analysis-inputs','coloc-platform-analysis-inputs']:
        record=json.loads((root/'config'/f'{name}.json').read_text())
        for path,digest in record.get('file_sha256',record.get('frozen_file_sha256',{})).items():add(path,digest)
    for name in ['coloc-gwas-queries','coloc-qtl-queries','coloc-reference-queries','coloc-published-bf-queries']:
        record=json.loads((root/'config'/f'{name}.json').read_text())
        for query in record['queries']:add(query['file'],query['sha256'])
    for name in ['coloc-run-provenance','coloc-platform-run-provenance']:
        run=json.loads((root/'reports'/f'{name}.json').read_text())
        executor='scripts/execute_coloc_platform.py' if name=='coloc-platform-run-provenance' else 'scripts/execute_coloc.py'
        add(executor,run['executor_sha256'])
        for job in run['jobs']:
            if job['status']!='complete':raise ValueError('Unfinished recorded analysis')
            for path,digest in job['output_sha256'].items():add(path,digest)
    summary=json.loads((root/'reports/coloc-summary.json').read_text())
    for key in ['source_sha256','generated_sha256']:
        for path,digest in summary[key].items():add(path,digest)
    for name in ['coloc-sources','coloc-platform-sources']:
        for source in json.loads((root/'config'/f'{name}.json').read_text())['sources']:
            add('data/raw/'+source['file'],source['sha256'])
    checked=0;missing_cache=[]
    for name,digest in sorted(pins.items()):
        path=root/name
        if not path.exists():
            if name.startswith('data/raw/') and not require_cache:
                missing_cache.append(name);continue
            raise ValueError('Missing required pinned file: '+name)
        with path.open('rb') as stream:
            actual=hashlib.file_digest(stream,'sha256').hexdigest()
        if actual!=digest:raise ValueError('Pinned file changed: '+name)
        checked+=1
    return {'status':'all_available_hashes_match','verified_files':checked,'ignored_cache_files_absent':len(missing_cache),
            'missing_cache':missing_cache,'scope':'Portable aggregate data, frozen protocols/code and results. Ignored cache absence is reported, never treated as a reproduced LD fit. Historical logs can differ in incidental warning text on replay and are not used to calculate results.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-cache',action='store_true')
    args=parser.parse_args()
    print(json.dumps(verify(require_cache=args.require_cache),indent=2))
