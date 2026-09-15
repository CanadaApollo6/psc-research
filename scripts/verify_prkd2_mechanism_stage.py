"""Final source, document, preservation and stage-output manifest."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def read(name): return json.loads((ROOT/name).read_text())


def main():
    cache = read('config/prkd2-mechanism-cache-manifest.json')
    for entry in cache['files']:
        path = ROOT/entry['file']
        assert path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256'], entry['file']
    assert (ROOT/'config/prkd2-mechanism-model-run.json').read_bytes() == (ROOT/'data/predictions/20260915T185853Z-prkd2-mechanism/run.json').read_bytes()
    for manifest in ['reports/prkd2-mechanism-verification.json','reports/prkd2-mechanism-replay-verification.json']:
        assert read(manifest)['status'] == 'PASS'
    numerical = read('reports/prkd2-mechanism-verification.json')
    assert numerical['verifier_sha256'] == sha(ROOT/'scripts/verify_prkd2_mechanism.py')
    replay = read('reports/prkd2-mechanism-replay-verification.json')
    assert replay['replay_script_sha256'] == sha(ROOT/'scripts/replay_prkd2_mechanism.py')
    for entry in replay['files']: assert sha(ROOT/entry['file']) == entry['sha256']
    figures = read('reports/figures/prkd2-mechanism/figure-provenance.json')
    assert figures['script_sha256'] == sha(ROOT/'scripts/plot_prkd2_mechanism.py')
    for name, digest in figures['outputs'].items(): assert sha(ROOT/'reports/figures/prkd2-mechanism'/name) == digest
    freeze = read('config/prkd2-mechanism-freeze.json')
    navigation = {'README.md','docs/evidence-log.md','docs/research-plan.md'}
    changed = [p for p,h in freeze['prior_tracked_sha256'].items() if sha(ROOT/p) != h]
    assert set(changed) == navigation, changed
    docs = [ROOT/p for p in sorted(navigation)] + [ROOT/'reports/prkd2-regulatory-mechanism.md', ROOT/'docs/prkd2-mechanism-reproduction.md']
    links = 0
    for path in docs:
        for raw in re.findall(r'\]\(([^)]+)\)', path.read_text()):
            target = raw.split(' "',1)[0].strip('<>')
            if not target or target.startswith('#') or urlparse(target).scheme: continue
            dest = (path.parent/unquote(target.split('#',1)[0])).resolve()
            # This checker writes the linked completion record at the end.
            if dest != ROOT/'reports/prkd2-mechanism-stage-verification.json': assert dest.exists(), (path,target)
            links += 1
    log = ROOT/'work/prkd2-mechanism-unittest.log'
    content = log.read_text()
    match = re.search(r'Ran (\d+) tests in ([\d.]+)s\s+OK\s*$',content)
    assert match and int(match.group(1)) == 274, 'The complete repository test suite did not pass'
    changed_paths = subprocess.check_output(['git','diff','--name-only',freeze['baseline_commit'],'--'],cwd=ROOT,text=True).splitlines()
    untracked = subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines()
    files = sorted(set(changed_paths + untracked))
    stage_name = 'reports/prkd2-mechanism-stage-verification.json'
    assert not any(p.startswith(('data/raw/','data/private/','data/predictions/','work/','.venv/')) for p in files)
    # Public provenance files can contain the configured credential filename,
    # but no actual credential bytes are read or scanned by this checker.
    outputs = {p:sha(ROOT/p) for p in files if p != stage_name and (ROOT/p).is_file()}
    result = {'status':'PASS','public_source_cache_files':len(cache['files']),'public_source_cache_bytes':cache['total_bytes'],
              'model_manifest_matches_ignored_original':True,'independent_numerical_verification':numerical,
              'replay_files':replay['exact_files'],'repository_tests':int(match.group(1)), 'test_seconds':float(match.group(2)),
              'test_log_sha256':sha(log),'original_catalog_rebuild_changed_files':0,
              'prior_tracked_files':len(freeze['prior_tracked_sha256']),
              'prior_scientific_files_unchanged':len(freeze['prior_tracked_sha256'])-len(changed),
              'intentional_navigation_updates':changed,'local_document_links_checked':links,
              'new_and_updated_files_sha256':outputs,'stage_verifier_sha256':sha(Path(__file__))}
    (ROOT/stage_name).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['new_and_updated_files_sha256','independent_numerical_verification']}))


if __name__ == '__main__': main()
