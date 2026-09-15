"""Frozen PRKD2 analysis paths and integrity checks."""

import json
from coloc_common import ROOT, ChainMap, complement, gzip_bytes, reciprocal_base, save_json, sha


def verify_plan():
    plan = json.loads((ROOT / 'config/prkd2-plan.json').read_text())
    for name, digest in plan['frozen_file_sha256'].items():
        if sha((ROOT / name).read_bytes()) != digest:
            raise ValueError('Frozen PRKD2 input changed: ' + name)
    return plan
