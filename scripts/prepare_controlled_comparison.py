"""Nominate and reference-check a prospective panel without requesting predictions."""

import csv
import gzip
import hashlib
import io
import json
import re
from collections import defaultdict
from decimal import Decimal
from pathlib import Path

from audit_finemap_archive import csv_bytes
from fetch_sources import validate

ROOT = Path(__file__).resolve().parents[1]


def select_loci(signals, rule):
    regions = defaultdict(list)
    for row in signals:
        regions[row['candidate_gene_as_published']].append(row)
    audit, eligible = [], []
    for name, rows in regions.items():
        first = rows[0]
        top = max(Decimal(r['posterior_probability_as_published']) for r in rows)
        if name in rule['exclude_prior_pilot_regions']:
            reason = 'Prior pilot region; outcomes already inspected'
        elif len(rows) != 1:
            reason = 'Multiple published signals'
        elif not Decimal(str(rule['minimum_published_top_pip'])) <= top < Decimal(str(rule['maximum_published_top_pip_exclusive'])):
            reason = 'Published top probability outside fixed range'
        else:
            reason = 'Eligible'
            eligible.append(first)
        audit.append({'table_label_as_published': name, 'published_signals': len(rows),
                      'highest_published_pip': str(top),
                      'published_set_size_for_single_signal': first['credible_set_size_as_published'] if len(rows) == 1 else '',
                      'eligibility': reason, 'selection_rank': '', 'selected': False})
    ordered = sorted(eligible, key=lambda r: (int(r['credible_set_size_as_published']),
                     -Decimal(r['posterior_probability_as_published']), int(r['chromosome']), int(r['lead_position_b37'])))
    selected = ordered[:rule['number_of_regions']]
    for rank, row in enumerate(ordered, 1):
        entry = next(a for a in audit if a['table_label_as_published'] == row['candidate_gene_as_published'])
        entry['selection_rank'] = rank
        entry['selected'] = rank <= rule['number_of_regions']
    return selected, audit


def mapped_record(record, rsid, build, chromosome):
    if record['name'] != rsid:
        raise ValueError('Renamed rsID requires explicit reconciliation')
    mappings = [m for m in record['mappings'] if m['assembly_name'] == build
                and m['seq_region_name'] == chromosome and m['coord_system'] == 'chromosome']
    if len(mappings) != 1 or mappings[0]['strand'] != 1:
        raise ValueError('Nonunique or nonforward chromosome mapping')
    return mappings[0]


def assess_mapping(mapping37, mapping38):
    alleles37, alleles38 = mapping37['allele_string'].split('/'), mapping38['allele_string'].split('/')
    if alleles37 != alleles38:
        return 'Allele strings differ between genome builds'
    if any(m['start'] != m['end'] for m in (mapping37, mapping38)) or any(len(a) != 1 or a not in 'ACGT' for a in alleles38):
        return 'Indel or other non-SNV record; excluded from this SNV comparison'
    if len(alleles38) != 2 or len(set(alleles38)) != 2:
        return 'Multiallelic or non-biallelic database record; original study alleles unavailable'
    return ''


def noncoding_context(term):
    if term == 'intron_variant':
        return 'intron'
    if term in ('5_prime_UTR_variant', '3_prime_UTR_variant'):
        return term
    if term in ('intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant'):
        return 'intergenic/upstream/downstream'
    return None


def nominate(rows, chromosome, config, read_json, read_base):
    nominated = []
    for source in sorted(rows, key=lambda r: (-Decimal(r['snp_prob_as_published']),
                          int(r['variant_id_as_published'][2:]) if re.fullmatch(r'rs\d+', r['variant_id_as_published']) else float('inf'),
                          r['variant_id_as_published'])):
        rsid = source['variant_id_as_published']
        result = {key: source[key] for key in ['dataset_id', 'variant_id_as_published', 'identifier_type',
                                              'source_row', 'snp_prob_as_published']}
        result.update({'chromosome': chromosome, 'position_grch37_1based': '', 'position_grch38_1based': '',
                       'alleles_grch37_as_database': '', 'alleles_grch38_as_database': '',
                       'reference_bases': '', 'alternate_bases': '', 'most_severe_consequence_grch38': '',
                       'noncoding_context': '', 'reference_eligible': False, 'role': 'excluded', 'reason': '',
                       'source_finemap_alleles_present': False, 'disease_risk_allele': ''})
        if not re.fullmatch(r'rs\d+', rsid):
            result['reason'] = 'No single literal rsID; build/alleles or identifier alias require reconciliation'
        else:
            records = {build: read_json(f'finemap-{build}-{rsid}.json') for build in ('GRCh37', 'GRCh38')}
            mappings = {build: mapped_record(records[build], rsid, build, chromosome) for build in records}
            for build in mappings:
                result[f'position_{build.lower()}_1based'] = mappings[build]['start']
                result[f'alleles_{build.lower()}_as_database'] = mappings[build]['allele_string']
            result['reason'] = assess_mapping(mappings['GRCh37'], mappings['GRCh38'])
            term = records['GRCh38']['most_severe_consequence']
            result['most_severe_consequence_grch38'] = term
            result['noncoding_context'] = noncoding_context(term) or ''
            if not result['reason'] and not result['noncoding_context']:
                result['reason'] = 'Consequence outside fixed noncoding classes'
            if not result['reason']:
                ref, alt = mappings['GRCh38']['allele_string'].split('/')
                bases = [read_json(f'finemap-{build}-{rsid}-sequence.json')['seq'].upper() for build in records]
                bases.append(read_base(f'finemap-model-reference-{rsid}.txt').strip().upper())
                if any(base != ref for base in bases):
                    raise ValueError(f'Reference mismatch for {rsid}')
                result.update(reference_bases=ref, alternate_bases=alt, reference_eligible=True,
                              role='reserve', reason='Verified biallelic SNV; retain as LD-screen reserve if not selected')
        nominated.append(result)
    eligible = [r for r in nominated if r['reference_eligible']]
    for row in eligible[:config['maximum_per_region']]:
        row['role'] = 'anchor_pending_comparators'
        row['reason'] = 'Selected by fixed descending-PIP rule; MAF and comparator eligibility still required'
    return nominated


def main():
    protocol_path = ROOT / 'config/psc-controlled-comparison.json'
    manifest_path = ROOT / 'config/finemap-sources.json'
    lock_path = ROOT / 'config/psc-controlled-comparison-lock.json'
    protocol = json.loads(protocol_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    for source in manifest['sources']:
        validate((ROOT / 'data/raw' / source['file']).read_bytes(), source)
    signal_path = ROOT / 'data/derived/psc-signal-summaries.csv'
    if hashlib.sha256(signal_path.read_bytes()).hexdigest() != protocol['table1_csv_sha256']:
        raise ValueError('Source Table 1 changed after protocol definition')
    archive = next(s for s in manifest['sources'] if s['id'] == 'psc-finemap-archive')
    if archive['sha256'] != protocol['archive_sha256']:
        raise ValueError('Fine-mapping source changed after protocol definition')
    signals = list(csv.DictReader(io.StringIO(signal_path.read_text())))
    selected, selection_audit = select_loci(signals, protocol['locus_selection'])
    if [r['candidate_gene_as_published'] for r in selected] != protocol['locus_selection']['selected_table_labels']:
        raise ValueError('Locus rule does not reproduce the frozen selection')
    table_path = ROOT / 'data/derived/psc-finemap-variants.csv.gz'
    with gzip.open(table_path, 'rt') as stream:
        all_rows = list(csv.DictReader(stream))
    candidates, windows, pools = [], [], []
    read_json = lambda name: json.loads((ROOT / 'data/raw' / name).read_text())
    read_base = lambda name: (ROOT / 'data/raw' / name).read_text()
    for signal, label in zip(selected, protocol['locus_selection']['selected_archive_labels']):
        rows = [r for r in all_rows if r['dataset_id'] == 'gwas/' + label]
        high = [r for r in rows if Decimal(r['snp_prob_as_published']) >= Decimal(str(protocol['anchor_selection']['minimum_archived_marginal_pip']))]
        candidate_rows = nominate(high, signal['chromosome'], protocol['anchor_selection'], read_json, read_base)
        # Coordinate identity checks only where Table 1 directly supplies this rsID.
        for candidate in candidate_rows:
            for identifier_field, position_field in [('lead_gwas_snp', 'lead_position_b37'), ('highest_pp_snp', 'highest_pp_position_b37')]:
                if candidate['position_grch37_1based'] and candidate['variant_id_as_published'] == signal[identifier_field]:
                    if candidate['position_grch37_1based'] != int(signal[position_field]):
                        raise ValueError('Candidate mapping conflicts with the directly tabulated build-37 position')
        candidates.extend(candidate_rows)
        anchors = [r for r in candidate_rows if r['role'] == 'anchor_pending_comparators']
        if anchors:
            center = anchors[0]['position_grch38_1based'] - 1
            start = center - protocol['model']['interval_width'] // 2
            end = start + protocol['model']['interval_width']
            if start < 0 or any(not start <= r['position_grch38_1based'] - 1 < end for r in anchors):
                raise ValueError('Anchor outside common region window')
            windows.append({'archive_region': label, 'chromosome': 'chr' + signal['chromosome'],
                            'center_variant': anchors[0]['variant_id_as_published'], 'interval_start_0based': start,
                            'interval_end_0based_exclusive': end, 'nominated_anchors': len(anchors)})
        low = [r for r in rows if Decimal(r['snp_prob_as_published']) <= Decimal(str(protocol['comparators']['maximum_archived_marginal_pip']))]
        pools.append({'archive_region': label, 'total_source_variants': len(rows), 'high_pip_source_records': len(high),
                      'reference_eligible_high_pip_snvs': sum(r['reference_eligible'] for r in candidate_rows),
                      'nominated_anchors': len(anchors), 'low_pip_source_records': len(low),
                      'low_pip_literal_rsid_records': sum(r['identifier_type'] == 'rsid' for r in low),
                      'fully_matched_comparators': 0})
    outputs = {
        'data/derived/psc-comparison-locus-selection.csv': csv_bytes(selection_audit),
        'data/derived/psc-comparison-candidates.csv': csv_bytes(candidates),
        'data/derived/psc-comparison-windows.csv': csv_bytes(windows),
        'data/derived/psc-comparison-pool-summary.csv': csv_bytes(pools),
    }
    hashes = {'protocol_sha256': hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
              'source_manifest_sha256': hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
              'complete_variant_table_sha256': hashlib.sha256(table_path.read_bytes()).hexdigest(),
              'candidate_table_sha256': hashlib.sha256(outputs['data/derived/psc-comparison-candidates.csv']).hexdigest()}
    if lock_path.exists():
        lock = json.loads(lock_path.read_text())
        for name, value in hashes.items():
            if lock[name] != value:
                raise ValueError(f'Frozen comparison input changed: {name}; review and version explicitly')
    for name, data in outputs.items():
        (ROOT / name).write_bytes(data)
    result = {**hashes, 'protocol_version': protocol['protocol_version'], 'protocol_recorded_at_utc': protocol['recorded_at_utc'],
              'selected_archive_regions': protocol['locus_selection']['selected_archive_labels'],
              'high_pip_source_records': len(candidates),
              'reference_verified_biallelic_snvs': sum(r['reference_eligible'] for r in candidates),
              'nominated_anchors': sum(r['role'] == 'anchor_pending_comparators' for r in candidates),
              'reserves': sum(r['role'] == 'reserve' for r in candidates),
              'excluded_source_records': sum(r['role'] == 'excluded' for r in candidates),
              'fully_matched_comparators': 0, 'model_requests_made': 0, 'ready_for_inference': False,
              'remaining_preparation': ['EUR frequencies and pairwise LD', 'Matched low-PIP comparator identities and covariates',
                                        'Pinned common gene universe and exact current track metadata',
                                        'Final matched input manifest before predictions'],
              'interpretation': 'Candidate inputs and comparison rules are fixed. This script does not perform matching or model inference. Source log discrepancies and incomplete high-PIP LD-screen coverage remain explicit.'}
    (ROOT / 'reports/psc-comparison-preparation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
