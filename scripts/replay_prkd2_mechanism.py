"""Replay new calculations using cached sources, without model or network requests."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    base = ROOT / 'work' / ('prkd2-mechanism-replay-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    base.mkdir(parents=True, exist_ok=False)
    stages = [
        ('results', 'scripts/summarize_prkd2_mechanism.py', [], ROOT/'data/derived/prkd2-mechanism-results'),
        ('training', 'scripts/audit_prkd2_mechanism_training.py', ['--offline'], ROOT/'data/derived/prkd2-mechanism-training'),
        ('peaks', 'scripts/analyze_prkd2_mechanism_peaks.py', ['--offline'], ROOT/'data/derived/prkd2-mechanism-encode')]
    checked = []
    for name, script, flags, source in stages:
        destination = base/name
        result = subprocess.run([sys.executable, str(ROOT/script), *flags, '--output-directory', str(destination)],
                                cwd=ROOT, text=True, capture_output=True)
        (base/(name+'.log')).write_text(result.stdout+result.stderr)
        if result.returncode: raise RuntimeError('Replay stage failed: '+name)
        for path in sorted(destination.iterdir()):
            original = source/path.name
            if not original.is_file() or sha(original) != sha(path):
                raise AssertionError('Replay differs: '+str(original))
            checked.append({'file':str(original.relative_to(ROOT)), 'sha256':sha(original)})
    result = {'status':'PASS', 'replay_directory':str(base.relative_to(ROOT)), 'network_or_model_requests':0,
              'exact_files':len(checked), 'files':checked, 'replay_script_sha256':sha(Path(__file__))}
    (ROOT/'reports/prkd2-mechanism-replay-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','exact_files':len(checked),'network_or_model_requests':0}))


if __name__ == '__main__': main()
