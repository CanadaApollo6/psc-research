"""Summarize completed audit verification, tests, links and preserved outputs."""
import hashlib
import importlib.metadata
import json
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'reports/ibdverse-model-audit-stage-verification.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    verification_path = ROOT / 'reports/ibdverse-model-audit-verification.json'
    replay_path = ROOT / 'reports/ibdverse-model-audit-replay-verification.json'
    verification = json.loads(verification_path.read_text())
    replay = json.loads(replay_path.read_text())
    if not verification['all_checks_pass'] or not replay['all_checks_pass']:
        raise ValueError('Independent verification or replay has not passed')
    if verification['verifier_sha256'] != sha(ROOT / 'scripts/verify_ibdverse_model_audit.py'):
        raise ValueError('Independent verifier changed after its recorded run')
    if replay['replay_script_sha256'] != sha(ROOT / 'scripts/replay_ibdverse_model_audit.py'):
        raise ValueError('Replay script changed after its recorded run')
    for filename, expected in replay['output_sha256'].items():
        if sha(ROOT / filename) != expected:
            raise ValueError('Replayed output changed: ' + filename)
    test_path = ROOT / 'data/raw/ibdverse-model-audit/full-tests.log'
    test_log = test_path.read_text()
    matches = re.findall(r'Ran (\d+) tests in ([0-9.]+)s', test_log)
    if not matches or not test_log.rstrip().endswith('OK'):
        raise ValueError('Full repository test suite did not pass')
    tests, seconds = matches[-1]
    baseline = verification['prior_preservation']['baseline_commit']
    from verify_ibdverse_model_audit import verify_previous_outputs
    preservation = verify_previous_outputs()
    if preservation != verification['prior_preservation']:
        raise ValueError('Prior preservation state changed after verification')
    links = 0
    for filename in ['README.md', 'docs/evidence-log.md', 'docs/ibdverse-model-audit-reproduction.md', 'reports/ibdverse-model-audit.md']:
        path = ROOT / filename
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)', path.read_text()):
            if target.startswith(('http:', 'https:', 'mailto:', 'codex:', '#')):
                continue
            target_path = unquote(target.split('#', 1)[0].strip('<>'))
            resolved = (path.parent / target_path).resolve()
            if resolved != DEST and not resolved.exists():
                raise ValueError('Broken local documentation link: ' + filename + ' -> ' + target)
            links += 1
    modified = set(subprocess.check_output(['git', 'diff', '--name-only', baseline], cwd=ROOT, text=True).splitlines())
    modified.update(subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard'], cwd=ROOT, text=True).splitlines())
    modified.discard(str(DEST.relative_to(ROOT)))
    allowed = {'README.md', 'docs/evidence-log.md'}
    if any('ibdverse_model_' not in name and 'ibdverse-model-audit' not in name and name not in allowed for name in modified):
        raise ValueError('Unrelated file in this audit stage')
    if any(name.startswith(('data/raw/', 'data/private/', 'data/predictions/')) for name in modified):
        raise ValueError('Raw/private source included in versioned change list')
    result = {'all_checks_pass': True, 'stage_verifier_sha256': sha(Path(__file__)),
              'baseline_commit': baseline,
              'unit_tests': {'passed': int(tests), 'seconds': float(seconds), 'log_file': str(test_path.relative_to(ROOT)), 'log_sha256': sha(test_path)},
              'independent_verification_sha256': sha(verification_path), 'replay_verification_sha256': sha(replay_path),
              'original_catalog_outputs_unchanged': True, 'prior_preservation': preservation,
              'local_document_links_checked': links,
              'runtime': {'Python': sys.version, **{name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'h5py', 'requests']}},
              'stage_files_excluding_this_record': len(modified),
              'stage_file_sha256': {name: sha(ROOT / name) for name in sorted(modified)},
              'new_H4_posteriors': 0}
    DEST.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'stage_file_sha256'}, indent=2))


if __name__ == '__main__':
    main()
