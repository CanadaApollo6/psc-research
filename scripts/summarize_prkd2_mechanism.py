"""Reproducible, descriptive summaries of the prospectively fixed PRKD2 run."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from prkd2_mechanism_common import ROOT, sha, write_csv

RUN = ROOT / 'data/predictions/20260915T185853Z-prkd2-mechanism'
OUTPUT = ROOT / 'data/derived/prkd2-mechanism-results'
TARGET = 'ENSG00000105287'
LFC = 'GeneMaskLFCScorer(requested_output=RNA_SEQ)'
ACTIVE = 'GeneMaskActiveScorer(requested_output=RNA_SEQ)'


def one_value(frame, **key):
    selected = frame
    for k, v in key.items():
        selected = selected[selected[k] == v]
    if len(selected) > 1:
        raise ValueError('Ambiguous score identity: ' + str(key))
    return selected.iloc[0].to_dict() if len(selected) else None


def complete_rank(rows, metric, descending=True):
    """A missing fixed-universe value makes a rank unavailable, never smaller-N."""
    if any(r[metric] is None or not np.isfinite(r[metric]) for r in rows):
        return {}
    ordered = sorted(rows, key=lambda r: ((-1 if descending else 1) * r[metric], r['gene_id']))
    return {r['gene_id']: i + 1 for i, r in enumerate(ordered)}


def reference_ld(variants, output):
    ledger = json.loads((ROOT / 'config/prkd2-expanded-sources.json').read_text())
    receipt = next(s for s in ledger['files'] if s['file'] == 'data/raw/prkd2-expanded-EUR-dosages.tsv.gz')
    path = ROOT / receipt['file']
    if sha(path) != receipt['sha256']:
        raise ValueError('Pinned public-reference dosages changed')
    by_pos = {str(v['position_grch37_1based']): v for v in variants}
    arrays, identities = {}, []
    with gzip.open(path, 'rt') as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader)
        if header[:5] != ['chromosome', 'position_grch37', 'rsid', 'ref', 'alt'] or len(header) != 508:
            raise ValueError('Unexpected reference dosage schema')
        for row in reader:
            v = by_pos.get(row[1])
            if v is None:
                continue
            if row[0] != '19' or row[2] not in ['.', v['rsid']] or row[3:5] != [v['reference'], v['alternate']]:
                raise ValueError('Full reference allele identity mismatch')
            if v['rsid'] in arrays:
                raise ValueError('Duplicate reference record')
            arrays[v['rsid']] = np.array([float(x) if x != '.' else np.nan for x in row[5:]])
            identities.append({'selected_rsid': v['rsid'], 'source_rsid_as_published': row[2],
                               'source_chromosome': row[0], 'position_GRCh37_1based': int(row[1]),
                               'reference': row[3], 'alternate': row[4],
                               'match_basis': 'Exact GRCh37 chromosome, position, REF and ALT; source rsID absent if dot'})
    rows = []
    for left in variants:
        for right in variants:
            a, b = arrays.get(left['rsid']), arrays.get(right['rsid'])
            n, r = 0, None
            if a is not None and b is not None:
                keep = np.isfinite(a) & np.isfinite(b)
                n = int(keep.sum())
                if n > 2 and np.std(a[keep]) and np.std(b[keep]):
                    r = float(np.corrcoef(a[keep], b[keep])[0, 1])
            rows.append({'variant_a': left['rsid'], 'variant_b': right['rsid'], 'joint_donors': n,
                         'r_alt_dosage': r, 'r_squared': None if r is None else r*r,
                         'status': 'available' if r is not None else 'unavailable'})
    write_csv(output / 'reference-LD.csv', rows)
    (output / 'reference-LD-provenance.json').write_text(json.dumps({
        'source': receipt, 'calculation': 'Pearson correlation of exact GRCh37 ALT allele dosages, pairwise complete donors; square gives r2.',
        'selected_GRCh37_and_GRCh38_REF_ALT_agree': True, 'allele_identity_crosswalk': identities,
        'interpretation': '503 public 1000 Genomes phase 3 EUR donors. Descriptive reference LD, not study LD, phasing evidence or a colocalisation posterior.',
        'sample_bearing_outputs_written': False}, indent=2) + '\n')


def run(output):
    inputs = json.loads((ROOT / 'config/prkd2-mechanism-inputs.json').read_text())
    for file, digest in inputs['frozen_file_sha256'].items():
        if sha(ROOT / file) != digest:
            raise ValueError('Frozen input changed: ' + file)
    manifest = json.loads((RUN / 'run.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('Incomplete model run')
    variants, genes = inputs['variants'], inputs['fixed_gene_universe']
    tracks = [t for t in inputs['selected_track_metadata'] if t['modality'] == 'rna_seq']
    primary_names = {t['name'] for t in inputs['merged_primary_rna_track_metadata']}
    effect_rows, ranks, regulatory, splicing, outside = [], [], [], [], []
    for variant in variants:
        rsid = variant['rsid']
        entry = next(e for e in manifest['scores'] if e['rsid'] == rsid)
        if sha(RUN / entry['file']) != entry['sha256']:
            raise ValueError('Model tidy scores changed')
        frame = pd.read_csv(RUN / entry['file'], low_memory=False).fillna('')
        if set(frame.variant_id) != {f"chr19:{variant['position_grch38_1based']}:{variant['reference']}>{variant['alternate']}"}:
            raise ValueError('Variant mismatch')
        scored_genes = set(frame[frame.variant_scorer == LFC].gene_id)
        outside.append({'rsid': rsid, 'returned_RNA_genes': len(scored_genes), 'fixed_genes': len(genes),
                        'outside_fixed_universe': ';'.join(sorted(scored_genes - {g['gene_id'] for g in genes}))})
        variant_ranks = []
        for gene in genes:
            primary, primary_active = [], []
            for track in tracks:
                key = {'gene_id': gene['gene_id'], 'track_name': track['name'], 'track_strand': track['strand']}
                row = one_value(frame, variant_scorer=LFC, **key)
                act = one_value(frame, variant_scorer=ACTIVE, **key)
                value = float(row['raw_score']) if row else None
                active = float(act['raw_score']) if act else None
                is_primary = track['name'] in primary_names
                effect_rows.append({'rsid': rsid, 'gene_id': gene['gene_id'], 'ensembl_116_name': gene['gene_name'],
                                    'model_gencode_name': row['gene_name'] if row else '', 'ontology': track['ontology_curie'],
                                    'track_name': track['name'], 'primary': is_primary, 'rna_ln_lfc': value,
                                    'rna_active': active, 'status': 'available' if row else 'missing'})
                if is_primary:
                    primary.append(value)
                    primary_active.append(active)
            complete = len(primary) == len(primary_names) and all(v is not None for v in primary)
            nearest_distance = abs(variant['position_grch38_1based'] - gene['tss_grch38_1based'])
            variant_ranks.append({'rsid': rsid, 'gene_id': gene['gene_id'], 'ensembl_116_name': gene['gene_name'],
                                  'signed_rna_median': float(np.median(primary)) if complete else None,
                                  'absolute_rna_median': float(np.median(np.abs(primary))) if complete else None,
                                  'rna_active_median': float(np.median(primary_active)) if all(v is not None for v in primary_active) else None,
                                  'primary_tracks': len(primary), 'opposed_signs_across_primary_tracks': complete and min(primary) < 0 < max(primary),
                                  'nearest_tss_distance_bp': nearest_distance, 'fixed_gene_count': len(genes), 'status': 'available' if complete else 'missing'})
        rank_map = complete_rank(variant_ranks, 'absolute_rna_median')
        nearest_map = complete_rank(variant_ranks, 'nearest_tss_distance_bp', False)
        for r in variant_ranks:
            r.update(model_rank=rank_map.get(r['gene_id']), nearest_tss_rank=nearest_map[r['gene_id']])
        ranks.extend(variant_ranks)
        for track in inputs['selected_track_metadata']:
            mod = track['modality'].upper()
            if mod not in ['DNASE', 'ATAC', 'CAGE', 'CHIP_HISTONE']:
                continue
            if mod == 'CHIP_HISTONE' and track.get('histone_mark') not in ['H3K27ac', 'H3K4me1', 'H3K4me3']:
                continue
            width = 2001 if mod == 'CHIP_HISTONE' else 501
            for agg in ['DIFF_LOG2_SUM', 'ACTIVE_SUM']:
                scorer = f'CenterMaskScorer(requested_output={mod}, width={width}, aggregation_type={agg})'
                row = one_value(frame, variant_scorer=scorer, track_name=track['name'], track_strand=track['strand'])
                regulatory.append({'rsid': rsid, 'modality': mod, 'ontology': track['ontology_curie'], 'track_name': track['name'],
                                   'strand': track['strand'], 'histone_mark': track.get('histone_mark', ''), 'width_bp': width,
                                   'aggregation': agg, 'score': float(row['raw_score']) if row else None,
                                   'status': 'available' if row else 'missing'})
        splice_tracks = [('SPLICE_SITES', 'donor'), ('SPLICE_SITES', 'acceptor')] + [
            ('SPLICE_SITE_USAGE', t['name']) for t in inputs['selected_track_metadata'] if t['modality'] == 'splice_site_usage']
        for mod, name in splice_tracks:
            row = one_value(frame, gene_id=TARGET, output_type=mod, track_name=name)
            splicing.append({'rsid': rsid, 'gene_id': TARGET, 'modality': mod, 'track_name': name,
                             'max_absolute_score': float(row['raw_score']) if row else None,
                             'status': 'available; unsigned' if row else 'not returned; variant outside target gene mask'})
    output.mkdir(parents=True, exist_ok=True)
    for name, rows in [('all-gene-RNA-effects.csv', effect_rows), ('primary-gene-ranks.csv', ranks),
                       ('regulatory-scores.csv', regulatory), ('target-splicing.csv', splicing), ('gene-universe-accounting.csv', outside)]:
        write_csv(output / name, rows)
    target_rows = []
    for v in variants:
        r = next(r for r in ranks if r['rsid'] == v['rsid'] and r['gene_id'] == TARGET)
        controls = [c for c in variants if c['role'] == 'comparison' and c['matched_anchor'] == v['rsid']]
        values = [next(t['absolute_rna_median'] for t in ranks if t['rsid'] == c['rsid'] and t['gene_id'] == TARGET) for c in controls]
        expected = v['expected_alt_minus_ref_rna_sign']
        observed = None if r['signed_rna_median'] is None else int(np.sign(r['signed_rna_median']))
        target_rows.append({**r, 'role': v['role'], 'substitution': v['reference'] + '>' + v['alternate'],
                            'expected_alt_minus_ref_rna_sign': expected,
                            'source_direction_agrees': observed == expected if expected is not None else None,
                            'selected_comparison_count': len(controls),
                            'absolute_effect_minus_comparison_median': r['absolute_rna_median'] - float(np.median(values)) if values else None})
    write_csv(output / 'target-summary.csv', target_rows)
    reference_ld(variants, output)
    result = {'fixed_genes': len(genes), 'variants': len(variants), 'primary_tracks': len(primary_names),
              'all_fixed_RNA_tracks': len(tracks), 'all_RNA_rows': len(effect_rows),
              'missing_RNA_rows': sum(r['status'] == 'missing' for r in effect_rows),
              'source_run': str(RUN.relative_to(ROOT)), 'source_run_sha256': sha(RUN / 'run.json'),
              'input_manifest_sha256': sha(ROOT / 'config/prkd2-mechanism-inputs.json'),
              'target_summary': target_rows,
              'interpretation': 'Signed model effects and within-window ranks are descriptive. No significance test, new H4 posterior or causal probability.'}
    (output / 'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'fixed_genes': len(genes), 'RNA_rows': len(effect_rows), 'target_ranks': [r['model_rank'] for r in target_rows]}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-directory', type=Path, default=OUTPUT)
    args = parser.parse_args()
    run(args.output_directory)
