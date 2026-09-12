"""Verify the recorded follow-up without requiring ignored source caches."""

import argparse
import hashlib
import json
from pathlib import Path

from verify_coloc_artifacts import verify as verify_previous

ROOT = Path(__file__).resolve().parents[1]


def verify(root=ROOT, require_cache=False):
    summary = json.loads((root / 'reports/coverage-repair-summary.json').read_text())
    pins = {}
    for field in ['versioned_file_sha256', 'ignored_source_sha256']:
        for name, digest in summary[field].items():
            if name in pins and pins[name] != digest:
                raise ValueError('Conflicting source hashes: ' + name)
            pins[name] = digest
    missing, checked = [], 0
    for name, expected in sorted(pins.items()):
        path = root / name
        if not path.exists():
            if name.startswith('data/raw/') and not require_cache:
                missing.append(name)
                continue
            raise ValueError('Missing required file: ' + name)
        with path.open('rb') as stream:
            actual = hashlib.file_digest(stream, 'sha256').hexdigest()
        if actual != expected:
            raise ValueError('Recorded follow-up input or result changed: ' + name)
        checked += 1
    previous = verify_previous(root=root, require_cache=require_cache)
    return {'status': 'all_available_hashes_match', 'followup_verified_files': checked,
            'followup_ignored_cache_files_absent': len(missing), 'missing_cache': missing,
            'previous_colocalisation': previous,
            'scope': 'Versioned public aggregates and recorded local source bytes. Missing ignored caches are reported separately; this is not a fresh source reacquisition or an LD refit.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-cache', action='store_true')
    args = parser.parse_args()
    print(json.dumps(verify(require_cache=args.require_cache), indent=2))
