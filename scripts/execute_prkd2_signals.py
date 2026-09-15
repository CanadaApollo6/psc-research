"""Freeze and execute one PRKD2 statistical stage with input/output hashes."""

import argparse
import gzip
import json
import os
import subprocess
from datetime import datetime, timezone

from prkd2_common import ROOT, gzip_bytes, save_json, sha, verify_plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['compare'])
    args = parser.parse_args()
    verify_plan()
    path = ROOT / ('config/prkd2-' + args.stage + '-execution.json')
    if path.exists():
        previous = json.loads(path.read_text())
        if previous.get('status') == 'complete':
            for name, digest in {**previous['input_sha256'], **previous['output_sha256']}.items():
                assert sha((ROOT / name).read_bytes()) == digest, name
            print(json.dumps({'stage': args.stage, 'status': 'verified-cache'}))
            return
        raise ValueError('An earlier attempt exists; preserve and diagnose it before replay')
    inputs = [ROOT / n for n in ['config/prkd2-plan.json', 'config/prkd2-sources.json',
        'reports/prkd2-input-audit.json', 'scripts/run_prkd2.R', 'scripts/execute_prkd2.py',
        'scripts/run_coloc.R', 'scripts/run_coloc_platform.R', 'scripts/run_prkd2_signals.R', 'scripts/execute_prkd2_signals.py']]
    inputs += sorted((ROOT / 'data/derived/prkd2-inputs').glob('*.csv.gz'))
    inputs += sorted((ROOT / 'data/raw/prkd2-ld').glob('*.f64'))
    if args.stage == 'compare':
        for ds in ['QTD000021', 'QTD000504']:
            record_path = ROOT / ('config/prkd2-bf-query-' + ds + '.json')
            record = json.loads(record_path.read_text())
            assert record['archive_eof_verified']
            data = (ROOT / record['file']).read_bytes()
            assert sha(data) == record['sha256'] and sha(gzip.decompress(data)) == record['uncompressed_sha256']
            inputs += [record_path, ROOT / record['file']]
        inputs += [ROOT / 'data/derived/prkd2-published-credible-sets.csv.gz']
        inputs += sorted((ROOT / 'data/raw/prkd2-ld').glob('*.rds'))
        inputs += sorted((ROOT / 'reports').glob('prkd2-GWAS-fit-*.json'))
    before = {p: sha(p.read_bytes()) for folder in ['reports', 'data/derived', 'data/raw/prkd2-ld']
              for p in (ROOT / folder).glob('prkd2*') if p.is_file()}
    before.update({p: sha(p.read_bytes()) for p in (ROOT / 'data/raw/prkd2-ld').glob('*.rds')})
    record = {'stage': args.stage, 'status': 'started', 'started_at_utc': datetime.now(timezone.utc).isoformat(),
              'input_sha256': {str(p.relative_to(ROOT)): sha(p.read_bytes()) for p in inputs},
              'R': '4.4.3', 'coloc': '5.2.3', 'susieR': '0.14.2',
              'scheduling': 'Each stage pins its complete inputs before execution; published BF transfers may continue during unrelated marginal/GWAS fitting stages.'}
    save_json(path, record)
    log = ROOT / ('work/prkd2/' + args.stage + '.log')
    command = [str(ROOT / 'work/coloc-runtime/r/bin/Rscript'), '--vanilla', 'scripts/run_prkd2_signals.R', args.stage]
    environment = {**os.environ, 'OPENBLAS_NUM_THREADS': '4', 'OMP_NUM_THREADS': '4'}
    with log.open('w') as stream:
        process = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            stream.write(line)
            stream.flush()
            print(line, end='', flush=True)
        code = process.wait()
    changed = []
    for folder in ['reports', 'data/derived']:
        for p in (ROOT / folder).glob('prkd2*'):
            if p.is_file() and (p not in before or sha(p.read_bytes()) != before[p]):
                if p.name.endswith('.csv.gz'):
                    p.write_bytes(gzip_bytes(gzip.decompress(p.read_bytes())))
                changed.append(p)
    changed += [p for p in (ROOT / 'data/raw/prkd2-ld').glob('*.rds') if p not in before or sha(p.read_bytes()) != before[p]]
    record.update(status='complete' if code == 0 else 'failed', exit_code=code,
                  finished_at_utc=datetime.now(timezone.utc).isoformat(), log_file=str(log.relative_to(ROOT)),
                  log_sha256=sha(log.read_bytes()), output_sha256={str(p.relative_to(ROOT)): sha(p.read_bytes()) for p in changed})
    save_json(path, record)
    for name, digest in record['input_sha256'].items():
        if sha((ROOT / name).read_bytes()) != digest:
            raise ValueError('Analysis input changed during execution: ' + name)
    if code:
        raise RuntimeError('Statistical stage failed; full log and partial outputs are preserved')
    print(json.dumps({'stage': args.stage, 'status': 'complete', 'outputs': len(changed)}), flush=True)


if __name__ == '__main__':
    main()
