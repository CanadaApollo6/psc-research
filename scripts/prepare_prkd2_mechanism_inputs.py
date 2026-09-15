"""Select comparison variants and lock full gene/track universes before prediction."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.ipc as ipc
from alphagenome.models import variant_scoring_utils

from prkd2_mechanism_common import ROOT, RAW, PLAN, freeze, now, sha, save_new, fetch_public, write_csv, clean_records

OUT = ROOT / 'data/derived/prkd2-mechanism'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def mutation_class(reference, alternate):
    return reference + '>' + alternate if reference in 'CT' else reference.translate(COMPLEMENT) + '>' + alternate.translate(COMPLEMENT)


def is_cpg(triplet):
    if len(triplet) != 3:
        raise ValueError('Need exact reference trinucleotide')
    return triplet[1:] == 'CG' or triplet[:2] == 'CG'


def references(plan):
    window = plan['region']
    start0 = window['tss_grch38_1based'] - 1 - window['width'] // 2
    end0 = start0 + window['width']
    prior = json.loads((ROOT / 'config/prkd2-expanded-sources.json').read_text())['files']
    source = next(r for r in prior if r['file'] == 'data/raw/prkd2-GRCh38-region.txt')
    source_path = ROOT / source['file']
    if sha(source_path) != source['sha256']:
        raise ValueError('Previously verified Ensembl reference changed')
    seq = source_path.read_text().strip().upper()
    seq_start0 = source['start_1based'] - 1
    subset = seq[start0 - seq_start0:end0 - seq_start0]
    if len(subset) != window['width']:
        raise ValueError('Incomplete reference interval')
    fai_path = ROOT / 'data/raw/GRCh38.p13.genome.fa.fai'
    old_fai = next(r for r in json.loads((ROOT / 'config/reference-sources.json').read_text())['sources'] if r['file'] == fai_path.name)
    if sha(fai_path) != old_fai['sha256']:
        raise ValueError('Model reference index changed')
    fai_row = next(line.split('\t') for line in fai_path.read_text().splitlines() if line.startswith(window['chromosome'] + '\t'))
    length, offset, bases, line_width = map(int, fai_row[1:5])
    if not 0 <= start0 < end0 <= length:
        raise ValueError('Model reference out of bounds')
    byte_start = offset + (start0 // bases) * line_width + start0 % bases
    last0 = end0 - 1
    byte_end = offset + (last0 // bases) * line_width + last0 % bases
    url = 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/GRCh38.p13.genome.fa'
    # The indexed range must match exactly; the parent FASTA is not downloaded.
    model_path = fetch_public('model-GRCh38.p13-common-window.range', url, max_bytes=1100000,
                              headers={'Range': f'bytes={byte_start}-{byte_end}'})
    record = next(r for r in json.loads((ROOT / 'config/prkd2-mechanism-reference-sources.json').read_text())['sources'] if r['id'] == model_path.name)
    expected_prefix = f'bytes {byte_start}-{byte_end}/'
    if record['status_code'] != 206 or not record['response_headers'].get('Content-Range', '').startswith(expected_prefix):
        raise ValueError('Exact model reference range not honored')
    model_seq = ''.join(model_path.read_text().split()).upper()
    if model_seq != subset:
        raise ValueError('Ensembl and model reference differ in common interval')
    out_seq = RAW / 'references/common-window-sequence.txt'
    if out_seq.exists() and out_seq.read_text() != subset + '\n':
        raise ValueError('Existing decoded model reference changed')
    out_seq.write_text(subset + '\n')
    release_path = fetch_public('ensembl-release.json', 'https://rest.ensembl.org/info/data?content-type=application/json')
    releases = json.loads(release_path.read_text())['releases']
    if 116 not in releases:
        raise ValueError('Protocol requires Ensembl release 116; inspect before changing source version')
    genes_path = fetch_public('ensembl-common-window-genes.json',
                              f'https://rest.ensembl.org/overlap/region/human/19:{start0+1}-{end0}?feature=gene;content-type=application/json')
    gencode = ROOT / 'data/raw/match-model-gencode-v46.feather'
    gencode_receipt = next(r for r in json.loads((ROOT / 'config/psc-comparator-sources.json').read_text())['sources'] if r['file'] == gencode.name)
    if sha(gencode) != gencode_receipt['sha256']:
        raise ValueError('Model GENCODE annotation changed')
    model_gene_ids = set()
    model_transcripts = []
    with pa.memory_map(str(gencode), 'r') as stream:
        reader = ipc.open_file(stream)
        for i in range(reader.num_record_batches):
            batch = reader.get_batch(i).to_pandas()
            part = batch[batch.Chromosome.eq('chr19') & batch.Feature.eq('transcript')]
            for row in part.to_dict('records'):
                # Matches the pinned official model annotation extractor: negative End,
                # positive Start. These are source-model coordinates, not a new TSS definition.
                tss0 = int(row['End'] if row['Strand'] == '-' else row['Start'])
                if start0 <= tss0 < end0:
                    identifier = row['gene_id'].split('.')[0]
                    model_gene_ids.add(identifier)
                    model_transcripts.append({'gene_id': identifier, 'gene_name': row['gene_name'],
                                              'transcript_id': row['transcript_id'], 'strand': row['Strand'],
                                              'start0': int(row['Start']), 'end0': int(row['End']),
                                              'model_tss_query_coordinate0': tss0})
    genes = []
    for gene in json.loads(genes_path.read_text()):
        if gene['seq_region_name'] != '19' or gene['strand'] not in [-1, 1]:
            continue
        tss = gene['start'] if gene['strand'] == 1 else gene['end']
        if start0 <= tss - 1 < end0:
            gid = gene['gene_id'].split('.')[0]
            genes.append({'gene_id': gid, 'gene_name': gene.get('external_name', gid), 'gene_version': gene.get('version'),
                          'biotype': gene['biotype'], 'strand': gene['strand'], 'tss_grch38_1based': tss,
                          'gene_start_grch38_1based': gene['start'], 'gene_end_grch38_1based': gene['end'],
                          'in_model_transcript_tss_universe': gid in model_gene_ids})
    genes.sort(key=lambda r: r['gene_id'])
    if len({r['gene_id'] for r in genes}) != len(genes) or sum(r['gene_id'] == window['target_gene_id'] for r in genes) != 1:
        raise ValueError('Gene annotation identity failure')
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / 'gene-universe.csv', genes)
    write_csv(OUT / 'model-transcripts.csv', sorted(model_transcripts, key=lambda r: (r['gene_id'], r['transcript_id'])))
    save_new(OUT / 'reference-audit.json', {'start0': start0, 'end0': end0, 'width': len(subset),
             'sequence_sha256': hashlib.sha256(subset.encode()).hexdigest(), 'exact_model_ensembl_sequence_match': True,
             'ensembl_reference_source': source, 'model_fai_source': old_fai, 'model_gencode_source': gencode_receipt,
             'frozen_QTL_TSS_1based': window['tss_grch38_1based'],
             'ensembl_target': next(r for r in genes if r['gene_id'] == window['target_gene_id']),
             'model_TSS_coordinate_convention': 'positive Start, negative exclusive End, matching original model extractor'})
    return start0, end0, subset, genes


def select_variants(plan, start0, end0, sequence):
    path = ROOT / 'data/derived/prkd2-expanded-inputs/GWAS-full-evidence.csv.gz'
    baseline = freeze()['prior_tracked_sha256']
    if sha(path) != baseline[str(path.relative_to(ROOT))]:
        raise ValueError('Prior full GWAS changed')
    source_rows = pd.read_csv(path)
    literal_counts = source_rows.SNP.value_counts()
    covariates, exclusions = [], []
    for row in source_rows.to_dict('records'):
        rsid = str(row['SNP'])
        if row['status'] != 'eligible' or not re.fullmatch(r'rs\d+', rsid) or literal_counts[rsid] != 1:
            continue
        pos = int(row['position'])
        if not start0 + 1 < pos < end0:
            continue
        snp_parts = str(row['snp']).split(':')
        if len(snp_parts) != 4 or any(len(a) != 1 or a not in 'ACGT' for a in snp_parts[2:]):
            continue
        reference = sequence[pos - 1 - start0]
        if reference not in snp_parts[2:] or snp_parts[2] == snp_parts[3]:
            exclusions.append({'rsid': rsid, 'reason': 'canonical alleles do not include reference'})
            continue
        alternate = next(a for a in snp_parts[2:] if a != reference)
        triplet = sequence[pos - 2 - start0:pos + 1 - start0]
        maf, pvalue = float(row['control_maf']), float(row['pvalue'])
        if not math.isfinite(maf) or not 0 < maf <= .5 or not math.isfinite(pvalue) or not 0 < pvalue <= 1:
            continue
        covariates.append({'rsid': rsid, 'position_grch38_1based': pos, 'position_grch37_1based': int(row['pos']),
                           'reference': reference, 'alternate': alternate, 'triplet': triplet,
                           'mutation_class': mutation_class(reference, alternate), 'cpg': is_cpg(triplet),
                           'source_control_maf': maf, 'source_psc_pvalue': pvalue,
                           'distance_to_frozen_tss': abs(pos - plan['region']['tss_grch38_1based']),
                           'source_row': int(row['source_row']), 'snp': row['snp']})
    by_id = {r['rsid']: r for r in covariates}
    anchors = []
    for planned in plan['anchors']:
        row = by_id[planned['rsid']]
        if any(row[key] != planned[key] for key in ['position_grch38_1based', 'reference', 'alternate']):
            raise ValueError('Anchor disagrees with verified source/reference')
        anchors.append({**row, 'role': 'anchor', 'matched_anchor': planned['rsid'],
                        'expected_alt_minus_ref_rna_sign': planned['expected_alt_minus_ref_rna_sign']})
    chosen = {r['rsid'] for r in anchors}
    comparisons, audit, statuses = [], [], []
    rules = plan['comparison_selection']
    for anchor in anchors:
        eligible = []
        for row in covariates:
            if row['rsid'] in {a['rsid'] for a in anchors}:
                continue
            maf_diff = abs(row['source_control_maf'] - anchor['source_control_maf'])
            tss_diff = abs(row['distance_to_frozen_tss'] - anchor['distance_to_frozen_tss'])
            separation = abs(row['position_grch38_1based'] - anchor['position_grch38_1based'])
            reasons = []
            if row['source_psc_pvalue'] < .1: reasons.append('PSC P below comparison threshold')
            if row['mutation_class'] != anchor['mutation_class']: reasons.append('mutation class mismatch')
            if row['cpg'] != anchor['cpg']: reasons.append('CpG mismatch')
            if separation > rules['maximum_anchor_distance_bp']: reasons.append('anchor distance')
            if tss_diff > rules['maximum_absolute_tss_distance_difference_bp']: reasons.append('TSS distance difference')
            if maf_diff > rules['maximum_source_control_maf_difference'] + 1e-12: reasons.append('control MAF difference')
            if row['rsid'] in chosen: reasons.append('already selected')
            cost = maf_diff / rules['maximum_source_control_maf_difference'] + tss_diff / rules['maximum_absolute_tss_distance_difference_bp']
            item = {'anchor': anchor['rsid'], 'candidate': row['rsid'], 'eligible': not reasons,
                    'reason': '; '.join(reasons), 'maf_difference': maf_diff, 'tss_distance_difference': tss_diff,
                    'genomic_separation': separation, 'cost': cost}
            audit.append(item)
            if not reasons:
                eligible.append((cost, separation, int(row['rsid'][2:]), row))
        selected = [r[-1] for r in sorted(eligible)[:rules['count_per_anchor']]]
        statuses.append({'anchor': anchor['rsid'], 'eligible_comparisons': len(eligible),
                         'selected_comparisons': len(selected), 'selected_ids': ';'.join(r['rsid'] for r in selected),
                         'status': 'complete pair' if len(selected) == rules['count_per_anchor'] else 'comparison shortfall; anchor retained'})
        for row in selected:
            chosen.add(row['rsid'])
            comparisons.append({**row, 'role': 'comparison', 'matched_anchor': anchor['rsid'],
                                'expected_alt_minus_ref_rna_sign': None})
    inputs = anchors + comparisons
    # Reconcile exact literal rsID and GRCh38 allele identity with a versioned API snapshot.
    variation = fetch_public('selected-ensembl-variants.json', 'https://rest.ensembl.org/variation/human',
                              headers={'Content-Type': 'application/json', 'Accept': 'application/json'},
                              body={'ids': [r['rsid'] for r in inputs]})
    records = json.loads(variation.read_text())
    for row in inputs:
        record = records[row['rsid']]
        mappings = [m for m in record['mappings'] if m['assembly_name'] == 'GRCh38' and m['seq_region_name'] == '19' and m['coord_system'] == 'chromosome']
        if record['name'] != row['rsid'] or len(mappings) != 1:
            raise ValueError('Selected variant has ambiguous or renamed identity')
        m = mappings[0]
        if m['start'] != row['position_grch38_1based'] or m['end'] != m['start'] or m['strand'] != 1:
            raise ValueError('Selected variant coordinate mismatch')
        alleles = m['allele_string'].split('/')
        if alleles[0] != row['reference'] or row['alternate'] not in alleles:
            raise ValueError('Selected variant allele mismatch')
        row.update(chromosome='chr19', interval_start0=start0, interval_end0=end0,
                   ensembl_recorded_alleles=m['allele_string'], reference_verified=True)
    write_csv(OUT / 'comparison-covariates.csv', covariates)
    write_csv(OUT / 'comparison-selection-audit.csv', audit)
    write_csv(OUT / 'comparison-status.csv', statuses)
    write_csv(OUT / 'selected-variants.csv', inputs)
    if exclusions:
        write_csv(OUT / 'comparison-reference-exclusions.csv', exclusions)
    return inputs, statuses


def tracks(plan):
    record = json.loads((ROOT / 'config/prkd2-mechanism-model-metadata.json').read_text())
    all_rows, missing = [], []
    for source in record['files']:
        path = ROOT / source['file']
        if sha(path) != source['sha256']:
            raise ValueError('Fresh model metadata changed')
        frame = pd.read_csv(path)
        if 'ontology_curie' not in frame:
            continue
        if source['modality'] in ['rna_seq', 'splice_site_usage']:
            frame = variant_scoring_utils.merge_stranded_track_metadata(frame)
        frame = frame.copy()
        frame['modality'] = source['modality']
        wanted = [plan['model']['primary_ontology']] + plan['model']['supporting_ontologies']
        selected = frame[frame.ontology_curie.isin(wanted)]
        if source['modality'] == 'chip_histone':
            selected = selected[selected.histone_mark.isin(['H3K27ac', 'H3K4me1', 'H3K4me3'])]
        all_rows.extend(clean_records(selected))
        for ontology in wanted:
            missing.append({'modality': source['modality'], 'ontology': ontology,
                            'available_tracks': int(selected.ontology_curie.eq(ontology).sum())})
    primary = [r for r in all_rows if r['modality'] == 'rna_seq' and r['ontology_curie'] == plan['model']['primary_ontology']]
    if not primary or len({r['name'] for r in primary}) != len(primary):
        raise ValueError('Primary merged RNA metadata missing or duplicated')
    fields = sorted({k for r in all_rows for k in r})
    write_csv(OUT / 'selected-model-tracks.csv', all_rows, fields)
    write_csv(OUT / 'model-context-coverage.csv', missing)
    return all_rows, primary


def main():
    frozen = freeze()
    lock_path = ROOT / 'config/prkd2-mechanism-inputs.json'
    if lock_path.exists():
        locked = json.loads(lock_path.read_text())
        for name, expected in locked['frozen_file_sha256'].items():
            if sha(ROOT / name) != expected:
                raise ValueError('Final input lock changed: ' + name)
        print(json.dumps({'status': 'verified existing input lock', 'variants': len(locked['variants'])}))
        return
    plan = json.loads(PLAN.read_text())
    start0, end0, seq, genes = references(plan)
    variants, statuses = select_variants(plan, start0, end0, seq)
    metadata, primary = tracks(plan)
    if len(variants) > plan['model']['maximum_score_requests']:
        raise ValueError('Prediction request cap exceeded')
    files = [PLAN, ROOT/'config/prkd2-mechanism-freeze.json', ROOT/'config/prkd2-mechanism-model-metadata.json',
             ROOT/'config/prkd2-mechanism-reference-sources.json', Path(__file__), ROOT/'scripts/prkd2_mechanism_common.py',
             ROOT/'data/derived/prkd2-expanded-inputs/GWAS-full-evidence.csv.gz'] + sorted(OUT.glob('*'))
    record = {'locked_at_utc': now(), 'plan_sha256': frozen['plan_sha256'], 'ready_for_inference': True,
              'variants': variants, 'comparison_status': statuses, 'full_gene_annotation': genes,
              'fixed_gene_universe': [r for r in genes if r['in_model_transcript_tss_universe']],
              'selected_track_metadata': metadata, 'merged_primary_rna_track_metadata': primary,
              'track_merging_note': 'Official client removes duplicate unstranded names after strand-specific merge; three raw CD14 RNA rows yield one primary gene score, not three independent assays.',
              'frozen_file_sha256': {str(p.relative_to(ROOT)): sha(p) for p in files}}
    save_new(lock_path, record)
    print(json.dumps({'status': 'locked before new predictions', 'variants': [r['rsid'] for r in variants],
                      'gene_universe': len(record['fixed_gene_universe']), 'primary_rna_tracks': len(primary),
                      'comparisons': statuses, 'locked_at_utc': record['locked_at_utc']}))


if __name__ == '__main__':
    main()
