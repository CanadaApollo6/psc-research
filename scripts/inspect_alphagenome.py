"""Fetch model metadata without making variant predictions."""

import dataclasses
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import alphagenome
from alphagenome.models import dna_client

from alphagenome_access import load_key

ROOT = Path(__file__).resolve().parents[1]


def main():
    client = dna_client.create(load_key(), model_version=dna_client.ModelVersion.ALL_FOLDS, timeout=30)
    metadata = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
    folder = ROOT / 'data/predictions/metadata'
    folder.mkdir(parents=True, exist_ok=True)
    records = []
    for field in dataclasses.fields(metadata):
        frame = getattr(metadata, field.name)
        if frame is None:
            continue
        path = folder / f'{field.name}.csv'
        frame.to_csv(path, index=False)
        records.append({'output_type': field.name, 'tracks': len(frame), 'columns': list(frame.columns), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    (folder / 'inventory.json').write_text(json.dumps({'retrieved_at_utc': datetime.now(timezone.utc).isoformat(), 'client_version': alphagenome.__version__, 'model_selection': 'ALL_FOLDS', 'variant_predictions_run': False, 'files': records}, indent=2) + '\n')
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        # API exceptions may contain request context; do not print credentials.
        code = exc.code().name if callable(getattr(exc, 'code', None)) else type(exc).__name__
        raise SystemExit(f'AlphaGenome metadata request failed: {code}') from None
