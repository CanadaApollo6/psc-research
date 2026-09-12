"""Replay both ABF baselines in scratch storage, including from a fresh clone."""

import argparse
import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path

from coloc_common import ROOT, gzip_bytes, sha
from verify_coloc_artifacts import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript',type=Path,default=ROOT/'work/coloc-runtime/r/bin/Rscript')
    args=parser.parse_args()
    verify()
    (ROOT/'work').mkdir(exist_ok=True)
    checked=0
    with tempfile.TemporaryDirectory(prefix='coloc-baseline-replay-',dir=ROOT/'work') as directory:
        scratch=Path(directory)
        (scratch/'reports').mkdir();(scratch/'data').mkdir()
        for name in ['config','scripts','data/derived']:
            (scratch/name).symlink_to(ROOT/name,target_is_directory=True)
        for script,manifest in [('run_coloc.R','coloc-run-provenance.json'),('run_coloc_platform.R','coloc-platform-run-provenance.json')]:
            result=subprocess.run([str(args.rscript.resolve()),'--vanilla','scripts/'+script,'baseline'],cwd=scratch,
                                  env={**os.environ,'OPENBLAS_NUM_THREADS':'4','OMP_NUM_THREADS':'4'},capture_output=True,text=True)
            if result.returncode:
                print(result.stderr);raise RuntimeError('R replay failed: '+script)
            record=json.loads((ROOT/'reports'/manifest).read_text())
            job=next(j for j in record['jobs'] if j['id']=='baseline')
            for name,digest in job['output_sha256'].items():
                data=(scratch/name).read_bytes()
                if name.endswith('.csv.gz'):data=gzip_bytes(gzip.decompress(data))
                if sha(data)!=digest:raise ValueError('Replay differed from frozen result: '+name)
                checked+=1
    print(json.dumps({'status':'byte-identical-isolated-replay','baseline_files':checked,'requires_reference_genotypes':False,'historical_files_overwritten':False}))


if __name__=='__main__':main()
