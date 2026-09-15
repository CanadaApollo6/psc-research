"""Replay every statistical stage in a separate directory without changing originals."""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from prkd2_common import ROOT, gzip_bytes, save_json, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis', choices=['prkd2', 'prkd2-expanded'])
    args = parser.parse_args()
    prefix = args.analysis
    scratch = Path(tempfile.mkdtemp(prefix=prefix+'-replay-', dir=ROOT/'work/prkd2'))
    for folder in ['scripts', 'config', 'reports', 'data/derived', 'data/raw/'+prefix+'-ld']:
        (scratch/folder).mkdir(parents=True, exist_ok=True)
    code = ['run_coloc.R', 'run_coloc_platform.R', 'run_prkd2.R', 'run_prkd2_signals.R', 'run_prkd2_expanded.R']
    for name in code:
        shutil.copy2(ROOT/'scripts'/name, scratch/'scripts'/name)
    shutil.copy2(ROOT/'config/prkd2-plan.json', scratch/'config/prkd2-plan.json')
    for name in [prefix+'-inputs', 'prkd2-query-rows', 'prkd2-expanded-published-bfs']:
        shutil.copytree(ROOT/'data/derived'/name, scratch/'data/derived'/name)
    shutil.copy2(ROOT/'data/derived/prkd2-published-credible-sets.csv.gz', scratch/'data/derived/prkd2-published-credible-sets.csv.gz')
    matrix = 'data/raw/'+prefix+'-ld/GWAS-PRKD2.f64'
    (scratch/matrix).symlink_to(ROOT/matrix)
    environment = {**os.environ, 'OPENBLAS_NUM_THREADS': '4', 'OMP_NUM_THREADS': '4'}
    checked = []
    stage_rows = []
    for stage in ['baseline', 'fit', 'compare']:
        record_name = prefix+'-'+stage+'-execution.json'
        record = json.loads((ROOT/'config'/record_name).read_text())
        for name, digest in record['input_sha256'].items():
            assert sha((ROOT/name).read_bytes()) == digest, name
        code = 'run_prkd2_expanded.R' if prefix == 'prkd2-expanded' else ('run_prkd2_signals.R' if stage == 'compare' else 'run_prkd2.R')
        log = scratch/(stage+'.log')
        with log.open('w') as stream:
            result = subprocess.run([str(ROOT/'work/coloc-runtime/r/bin/Rscript'), '--vanilla', 'scripts/'+code, stage],
                                    cwd=scratch, env=environment, stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError('Replay failed; log preserved at '+str(log))
        for name, digest in record['output_sha256'].items():
            previous, current = ROOT/name, scratch/name
            assert sha(previous.read_bytes()) == digest, name
            data = current.read_bytes()
            if name.endswith('.csv.gz'):
                data = gzip_bytes(gzip.decompress(data))
                current.write_bytes(data)
            actual = sha(data)
            if actual != digest:
                raise ValueError('Replay bytes differ: '+name)
            checked.append({'file': name, 'sha256': digest})
        stage_rows.append({'stage': stage, 'outputs_reproduced': len(record['output_sha256']), 'log_sha256': sha(log.read_bytes())})
        print(json.dumps({'analysis': prefix, 'stage': stage, 'reproduced': True}), flush=True)
    output = {'all_checks_pass': True, 'analysis': prefix, 'isolated_replay_directory': str(scratch.relative_to(ROOT)),
              'stages': stage_rows, 'outputs_reproduced': len(checked), 'verified_files': checked,
              'method': 'Fresh R processes; original code and prepared inputs; all output hashes including fitted RDS objects identical',
              'replay_script_sha256': sha(Path(__file__).read_bytes())}
    save_json(ROOT/'reports'/(prefix+'-replay-verification.json'), output)
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
