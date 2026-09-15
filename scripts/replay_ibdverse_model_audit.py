"""Replay the audit from existing, pinned local sources without network access."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    paths = sorted(p for p in (ROOT / 'data/derived/ibdverse-model-audit').glob('*') if p.suffix in {'.csv', '.tsv'})
    paths += [ROOT / 'reports' / ('ibdverse-model-audit' + suffix + '.json')
              for suffix in ['', '-h5ad-schema', '-h5ad-metadata']]
    before = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    ledgers = sorted((ROOT / 'config').glob('ibdverse-model-audit*.json'))
    ledger_hashes = {str(p.relative_to(ROOT)): sha(p) for p in ledgers}
    commands = [
        ['scripts/fetch_ibdverse_model_members.py', '--offline'],
        ['scripts/read_ibdverse_model_metadata.py', '--offline'],
        ['scripts/analyze_ibdverse_model_audit.py'],
        ['scripts/verify_ibdverse_model_audit.py'],
    ]
    for command in commands:
        result = subprocess.run([sys.executable, *command], cwd=ROOT, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError('Offline replay failed: ' + ' '.join(command) + '\n' + result.stderr)
    for name, expected in before.items():
        if sha(ROOT / name) != expected:
            raise ValueError('Audit replay changed ' + name)
    for name, expected in ledger_hashes.items():
        if sha(ROOT / name) != expected:
            raise ValueError('Offline replay changed a frozen acquisition record: ' + name)
    record = {'all_checks_pass': True, 'replay_script_sha256': sha(Path(__file__)),
              'network_requests': 0, 'model_requests': 0, 'new_H4_posteriors': 0,
              'commands': commands, 'generated_files_replayed_exactly': len(before) - 1,
              'additional_coloc_inventory_rebuilt_and_verified_exactly': 1,
              'source_acquisition_records_unchanged': len(ledger_hashes), 'output_sha256': before}
    output = ROOT / 'reports/ibdverse-model-audit-replay-verification.json'
    output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k: v for k, v in record.items() if k != 'output_sha256'}, indent=2))


if __name__ == '__main__':
    main()
