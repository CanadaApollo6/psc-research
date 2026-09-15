"""Read bounded HDF5 observation metadata, never expression or genotypes."""
import argparse
import io
import json
import math

import h5py
import numpy as np

from fetch_ibdverse_model_audit import ROOT, RAW, fetch, freeze, save, sha


class MetadataReader(io.RawIOBase):
    def __init__(self, plan, offline=False):
        self.plan, self.offline, self.position = plan, offline, 0
        self.blocks = {}
        self.size = plan['source_bytes']
        self.manifest_path = ROOT / 'config/ibdverse-model-audit-h5ad-queries.json'
        plan_sha = sha(ROOT / 'config/ibdverse-model-audit-h5ad-plan.json')
        self.record = json.loads(self.manifest_path.read_text()) if self.manifest_path.exists() else {
            'plan_sha256': plan_sha, 'source_url': plan['source_url'], 'blocks': []}
        if self.record['plan_sha256'] != plan_sha:
            raise ValueError('H5AD access plan changed')

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.position

    def seek(self, offset, whence=0):
        self.position = offset if whence == 0 else self.position + offset if whence == 1 else self.size + offset
        if self.position < 0:
            raise ValueError('Negative byte offset')
        return self.position

    def block(self, start):
        if start in self.blocks: return self.blocks[start]
        length = min(self.plan['http_cache_block_bytes'], self.size - start)
        previous = next((r for r in self.record['blocks'] if r['start'] == start), None)
        if previous:
            p = ROOT / previous['file']
            if sha(p) != previous['sha256']: raise ValueError('Cached HDF5 range changed')
            data = p.read_bytes()
        else:
            if self.offline: raise ValueError('Uncached HDF5 byte range')
            if sum(r['bytes'] for r in self.record['blocks']) + length > self.plan['maximum_new_transfer_bytes']:
                raise ValueError('HDF5 transfer cap exceeded')
            if len(self.record['blocks']) >= self.plan['maximum_http_requests']:
                raise ValueError('HDF5 request cap exceeded')
            result = fetch('IBDverse-h5ad-block-' + str(start), self.plan['source_url'], limit=length,
                           headers={'Range': f'bytes={start}-{start + length - 1}'})
            expected = f'bytes {start}-{start + length - 1}/{self.size}'
            if not result.get('ok') or result['response_headers'].get('Content-Range') != expected:
                raise ValueError('HDF5 exact range or remote size verification failed')
            identity = {k: result['response_headers'].get(k) for k in ['ETag', 'Last-Modified']}
            if 'identity' in self.record and self.record['identity'] != identity:
                raise ValueError('HDF5 remote identity changed')
            self.record['identity'] = identity
            self.record['blocks'].append({k: result[k] for k in ['file', 'bytes', 'sha256']} | {'start': start})
            save(self.manifest_path, self.record)
            data = (ROOT / result['file']).read_bytes()
        if len(data) != length: raise ValueError('Cached HDF5 range size differs')
        self.blocks[start] = data
        return data

    def read(self, length=-1):
        if length < 0: raise ValueError('Unbounded HDF5 read prohibited')
        length = min(length, self.size - self.position)
        if length <= 0: return b''
        n = self.plan['http_cache_block_bytes']
        start, end = self.position, self.position + length
        data = b''.join(self.block(a) for a in range(start // n * n, math.ceil(end / n) * n, n))
        self.position = end
        return data[start % n:start % n + length]

    def readinto(self, buffer):
        data = self.read(len(buffer)); buffer[:len(data)] = data
        return len(data)


def read_column(obj):
    if isinstance(obj, h5py.Dataset):
        data = obj[:]
        return np.array([x.decode() if isinstance(x, bytes) else x for x in data], dtype=object)
    if set(obj.keys()) != {'categories', 'codes'}:
        raise ValueError('Unsupported observation encoding')
    cats = read_column(obj['categories'])
    codes = obj['codes'][:]
    if np.any(codes < -1) or np.any(codes >= len(cats)):
        raise ValueError('Invalid category index')
    values = np.empty(len(codes), dtype=object); values[:] = None
    keep = codes >= 0; values[keep] = cats[codes[keep]]
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--schema-only', action='store_true')
    args = parser.parse_args(); freeze()
    plan = json.loads((ROOT / 'config/ibdverse-model-audit-h5ad-plan.json').read_text())
    reader = MetadataReader(plan, args.offline)
    with h5py.File(reader, 'r') as f:
        obs = f['obs']
        schema = [{'name': name, 'encoding': 'categorical' if isinstance(obs[name], h5py.Group) else 'dataset',
                   'shape': list(obs[name]['codes'].shape) if isinstance(obs[name], h5py.Group) and 'codes' in obs[name] else list(obs[name].shape) if isinstance(obs[name], h5py.Dataset) else None}
                  for name in obs.keys()]
        result = {'schema': schema, 'h5py_version': h5py.__version__, 'expression_data_read': False,
                  'actual_executed_model_N_recovered': False, 'actual_selected_PCs_recovered': False}
        save(ROOT / 'reports/ibdverse-model-audit-h5ad-schema.json', result)
        print(json.dumps(result), flush=True)
        if not args.schema_only:
            selected = [c for c in plan['initial_columns'] if c in obs]
            columns = {c: read_column(obs[c]) for c in selected}
            import pandas as pd
            frame = pd.DataFrame(columns)
            result['columns_read'] = selected
            result['observation_rows'] = len(frame)
            result['counts'] = []
            donor = 'Genotyping_ID' if 'Genotyping_ID' in frame else 'individual_id' if 'individual_id' in frame else None
            annotation = 'predicted_labels_tissue' if 'predicted_labels_tissue' in frame else None
            if donor is None and {'sanger_sample_id', 'predicted_labels', 'tissue'} <= set(frame.columns):
                amendment = json.loads((ROOT / 'config/ibdverse-model-audit-h5ad-mapping-amendment.json').read_text())
                if amendment['original_plan_sha256'] != sha(ROOT / 'config/ibdverse-model-audit-h5ad-plan.json'):
                    raise ValueError('Amended HDF5 plan changed')
                if amendment['observed_schema_sha256'] != sha(ROOT / 'reports/ibdverse-model-audit-h5ad-schema.json'):
                    raise ValueError('Observed HDF5 schema changed')
                mapping_path = ROOT / amendment['mapping_file']
                if sha(mapping_path) != amendment['mapping_sha256']:
                    raise ValueError('Published sample mapping changed')
                mapping = pd.read_csv(mapping_path, sep='\t', dtype=str)
                if mapping.sanger_sample_id.duplicated().any():
                    raise ValueError('Sample mapping is not one-to-one')
                frame['individual_id'] = frame.sanger_sample_id.map(mapping.set_index('sanger_sample_id').individual_id)
                frame['predicted_labels_tissue'] = frame.predicted_labels + '_' + frame.tissue
                donor, annotation = 'individual_id', 'predicted_labels_tissue'
                result['mapping_amendment_sha256'] = sha(ROOT / 'config/ibdverse-model-audit-h5ad-mapping-amendment.json')
                result['unmapped_cells_all_annotations'] = int(frame.individual_id.isna().sum())
            if donor and annotation:
                for context in plan['counts_only_for']:
                    part = frame[frame[annotation].eq(context)]
                    counts = part.groupby(donor, dropna=True).size()
                    result['counts'].append({'context': context, 'donor_column': donor,
                        'cells': len(part), 'donor_labels_with_any_cells': len(counts),
                        'donor_labels_with_at_least_five_cells': int((counts >= plan['minimum_cells_per_donor']).sum()),
                        'cells_missing_donor_label': int(part[donor].isna().sum())})
            result['interpretation'] = 'Released-object RNA eligibility counts only. Exact mapping to the executed genotype/covariate model remains unverified.'
            save(ROOT / 'reports/ibdverse-model-audit-h5ad-metadata.json', result)
            print(json.dumps({k: v for k, v in result.items() if k != 'schema'}), flush=True)


if __name__ == '__main__': main()
