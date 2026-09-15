"""Verify final provenance, planned failures, prior-file preservation and deliverables."""

import gzip
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

from prkd2_common import ROOT, save_json, sha, verify_plan


def main():
    plan=verify_plan()
    cache={}
    def check(name,digest):
        path=ROOT/name
        if name not in cache:
            cache[name]=sha(path.read_bytes())
        if cache[name]!=digest:
            raise ValueError('Artifact hash mismatch: '+name)
    source=json.loads((ROOT/'config/prkd2-sources.json').read_text())
    for item in source['sources']:
        check('data/raw/'+item['file'],item['sha256'])
    queries_checked=0
    for name in ['gwas','qtl','reference','missing-evidence']:
        record=json.loads((ROOT/f'config/prkd2-{name}-queries.json').read_text())
        for query in record['queries']:
            check(query['file'],query['sha256'])
            if 'uncompressed_sha256' in query:
                assert sha(gzip.decompress((ROOT/query['file']).read_bytes()))==query['uncompressed_sha256']
            queries_checked+=1
        for part in record.get('ranges',[]):
            check('data/raw/'+part['file'],part['sha256'])
    expanded=json.loads((ROOT/'config/prkd2-expanded-sources.json').read_text())
    for item in expanded['files']:
        check(item['file'],item['sha256'])
    archive_bytes=archive_rows=0
    for ds in ['QTD000021','QTD000504']:
        record=json.loads((ROOT/f'config/prkd2-bf-query-{ds}.json').read_text())
        assert record['archive_eof_verified']
        length=int(next(s for s in source['sources'] if s['id']==ds+'-lbf-header')['content_range'].split('/')[-1])
        assert record['compressed_archive_bytes']==length
        check(record['file'],record['sha256'])
        assert sha(gzip.decompress((ROOT/record['file']).read_bytes()))==record['uncompressed_sha256']
        archive_bytes+=length;archive_rows+=record['all_archive_rows_scanned']
    executions=[]
    for prefix in ['prkd2','prkd2-expanded']:
        for stage in ['baseline','fit','compare']:
            name=f'config/{prefix}-{stage}-execution.json'
            record=json.loads((ROOT/name).read_text())
            assert record['status']=='complete' and record['exit_code']==0
            for field in ['input_sha256','output_sha256']:
                for file,digest in record[field].items():check(file,digest)
            check(record['log_file'],record['log_sha256'])
            executions.append(name)
        table=pd.read_csv(ROOT/f'reports/{prefix}-signal-comparison.csv')
        assert len(table)==36 and not table.coverage_gate.any() and table['PP.H4.abf'].isna().all()
        baseline=pd.read_csv(ROOT/f'reports/{prefix}-baseline.csv')
        assert len(baseline)==12 and not baseline.coverage_gate.any()
    verifications=['prkd2-input-verification.json','prkd2-statistical-verification.json','prkd2-normalization-verification.json',
                   'prkd2-replay-verification.json','prkd2-expanded-replay-verification.json']
    for name in verifications:
        assert json.loads((ROOT/'reports'/name).read_text())['all_checks_pass']
    mapping=json.loads((ROOT/'reports/prkd2-missing-allele-verification.json').read_text())
    check(mapping['output'],mapping['output_sha256'])
    figures=json.loads((ROOT/'reports/prkd2-figure-provenance.json').read_text())
    for field in ['inputs','outputs']:
        for name,digest in figures[field].items():check(name,digest)
    log=(ROOT/'work/prkd2/python-tests.log').read_text()
    count=int(re.search(r'Ran (\d+) tests in',log).group(1))
    assert count==250 and log.rstrip().endswith('OK')
    rlog=(ROOT/'work/prkd2/R-qualification.log').read_text()
    assert 'both coverage-threshold sides passed' in rlog

    baseline=plan['baseline_commit']
    allowed={'README.md','docs/research-plan.md','docs/evidence-log.md','docs/data-sources.md','docs/computational-followup-roadmap.md'}
    tree=subprocess.check_output(['git','ls-tree','-r','-z',baseline],cwd=ROOT).split(b'\0')
    changed=[];preserved=0;original_files=0
    for entry in tree:
        if not entry:continue
        metadata,name=entry.split(b'\t',1);mode,kind,oid=metadata.decode().split();name=name.decode()
        assert kind=='blob'
        data=(ROOT/name).read_bytes()
        actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        original_files+=1
        if actual!=oid:
            assert name in allowed, 'Unexpected prior-file change: '+name
            changed.append(name)
        else:preserved+=1
    assert original_files==864 and set(changed)==allowed
    untracked=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines()
    diff=subprocess.check_output(['git','diff','--name-only',baseline],cwd=ROOT,text=True).splitlines()
    deliverables=sorted(set(untracked+diff)-{'reports/prkd2-stage-verification.json'})
    for name in deliverables:
        assert name in allowed or 'prkd2' in Path(name).name.lower() or '/prkd2-' in name, name
        assert not name.startswith(('data/raw/','data/private/','work/','.venv/')),name
    links=0
    for name in deliverables:
        if not name.endswith('.md'):continue
        path=ROOT/name
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',path.read_text()):
            if '://' in target or target.startswith(('#','mailto:')):continue
            destination=(path.parent/unquote(target.split('#')[0])).resolve()
            if destination==ROOT/'reports/prkd2-stage-verification.json':
                links+=1;continue
            assert destination.exists(),f'{name}: {target}'
            links+=1
    subprocess.run(['git','diff','--check'],cwd=ROOT,check=True)
    record={'all_checks_pass':True,'verified_at_utc':datetime.now(timezone.utc).isoformat(),
        'baseline_commit':baseline,'original_tracked_files':original_files,'prior_files_preserved_exactly':preserved,
        'intentional_prior_file_updates':sorted(changed),'source_or_execution_files_hash_checked':len(cache),
        'regional_or_targeted_queries_checked':queries_checked,'complete_published_BF_archive_bytes':archive_bytes,
        'complete_published_BF_archive_rows':archive_rows,'execution_records_verified':executions,
        'independent_verification_records':verifications,'repository_tests_passed':count,'R_qualification_passed':True,
        'catalog_rebuild_preserved_original_files':True,'component_rows':72,'coverage_passing_component_rows':0,
        'baseline_rows':24,'coverage_passing_baseline_rows':0,'relative_document_links_checked':links,
        'genotype_and_raw_cache_files_in_deliverables':False,
        'check_logs':{name:sha((ROOT/name).read_bytes()) for name in ['work/prkd2/python-tests.log','work/prkd2/R-qualification.log','work/prkd2/catalog-rebuild.log']},
        'deliverable_sha256':{name:sha((ROOT/name).read_bytes()) for name in deliverables},
        'verifier_sha256':sha(Path(__file__).read_bytes())}
    save_json(ROOT/'reports/prkd2-stage-verification.json',record)
    print(json.dumps({k:v for k,v in record.items() if k!='deliverable_sha256'},indent=2))


if __name__=='__main__':
    main()
