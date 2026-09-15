"""Replay new recovery calculations in an isolated work copy and compare bytes."""

import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, save_json, sha


def main():
    reports = ['prkd2-recovery-analysis.json', 'prkd2-recovery-secondary-audit.json',
               'prkd2-recovery-ibdverse-analysis.json', 'prkd2-recovery-ibdverse-format-audit.json']
    files = []
    for name in reports:
        result = json.loads((ROOT / 'reports' / name).read_text())
        files.extend(item['file'] for item in result.get('outputs', [result['output']] if 'output' in result else []))
        files.append('reports/' + name)
    files += ['data/derived/prkd2-recovery-independent-coverage.csv', 'data/derived/prkd2-recovery-independent-format.csv']
    expected = {name: sha((ROOT / name).read_bytes()) for name in sorted(set(files))}
    directory = ROOT / 'work' / ('prkd2-recovery-replay-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    directory.mkdir()
    shutil.copytree(ROOT / 'scripts', directory / 'scripts')
    shutil.copytree(ROOT / 'config', directory / 'config')
    (directory / 'reports').mkdir()
    (directory / 'data/derived').mkdir(parents=True)
    (directory / 'data/raw').symlink_to(ROOT / 'data/raw', target_is_directory=True)
    for name in ('prkd2-expanded-inputs', 'prkd2-query-rows', 'prkd2-recovery-queries'):
        shutil.copytree(ROOT / 'data/derived' / name, directory / 'data/derived' / name)
    shutil.copy2(ROOT / 'data/derived/prkd2-missing-evidence.csv', directory / 'data/derived/prkd2-missing-evidence.csv')
    for source in (ROOT / 'data/derived').glob('prkd2-recovery-*'):
        if source.is_file():
            shutil.copy2(source, directory / 'data/derived' / source.name)
    commands = [[sys.executable, 'scripts/' + name] for name in
                ['analyze_prkd2_recovery.py', 'audit_prkd2_recovery_secondary.py',
                 'analyze_prkd2_recovery_ibdverse.py', 'audit_prkd2_recovery_ibdverse_format.py']]
    commands.append([str(ROOT / 'work/coloc-runtime/r/bin/Rscript'), '--vanilla',
                     'scripts/verify_prkd2_recovery_coverage.R', 'data/derived/prkd2-recovery-independent-coverage.csv'])
    for i, command in enumerate(commands):
        with (directory / f'command-{i}.log').open('w') as log:
            subprocess.run(command, cwd=directory, check=True, stdout=log, stderr=subprocess.STDOUT)
    for name, digest in expected.items():
        if sha((directory / name).read_bytes()) != digest or sha((ROOT / name).read_bytes()) != digest:
            raise ValueError('Recovery replay or original preservation differs: ' + name)
    save_json(ROOT / 'reports/prkd2-recovery-replay-verification.json', {
        'all_checks_pass': True, 'replayed_files': len(expected),
        'replay_work_directory': str(directory.relative_to(ROOT)), 'verifier_sha256': sha(Path(__file__).read_bytes()),
        'completed_at_utc': datetime.now(timezone.utc).isoformat(),
        'method': 'Copied scripts/config/retained input tables; shared read-only-in-use original raw cache; five sequential analysis commands; exact byte comparison.',
        'replayed_file_sha256': expected})
    print(json.dumps({'all_checks_pass': True, 'replayed_files': len(expected), 'work_directory': str(directory.relative_to(ROOT))}))


if __name__ == '__main__':
    main()
