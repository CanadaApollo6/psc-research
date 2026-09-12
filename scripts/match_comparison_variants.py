"""Apply the fixed frequency, LD and annotation matching rule before inference."""

import csv
import gzip
import hashlib
import io
import itertools
import json
import math
from pathlib import Path

import numpy as np

from audit_finemap_archive import csv_bytes
from prepare_matching_covariates import static_pair

ROOT = Path(__file__).resolve().parents[1]


def maf(dosages):
    called = np.asarray(dosages, dtype=float)
    called = called[np.isfinite(called)]
    if not len(called):
        return None, 0
    frequency = float(called.sum() / (2 * len(called)))
    return min(frequency, 1 - frequency), len(called)


def dosage_r2(a, b, minimum=400):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError('Different genotype sample ordering/length')
    keep = np.isfinite(a) & np.isfinite(b)
    n = int(keep.sum())
    if n < minimum or np.ptp(a[keep]) == 0 or np.ptp(b[keep]) == 0:
        return None, n
    return float(np.corrcoef(a[keep], b[keep])[0, 1] ** 2), n


def parse_vcf(text, eur_samples, targets):
    header = next((line for line in text.splitlines() if line.startswith('#CHROM\t')), None)
    if header is None:
        raise ValueError('VCF column header missing')
    samples = header.split('\t')[9:]
    if len(set(samples)) != len(samples) or not set(eur_samples) <= set(samples):
        raise ValueError('Missing or duplicate EUR samples')
    indices = [samples.index(sample) for sample in eur_samples]
    keys = {}
    for row in targets:
        key = row['chromosome'], row['position_grch37_1based'], row['reference_bases'], row['alternate_bases']
        if key in keys:
            raise ValueError('Two rsIDs share an exact reference key; explicit alias review required')
        keys[key] = row['rsid']
    output, source_ids = {}, {}
    for line in text.splitlines():
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        if len(fields) != len(samples) + 9:
            raise ValueError('Truncated VCF sample row')
        key = fields[0], int(fields[1]), fields[3], fields[4]
        if key not in keys:
            raise ValueError('Unexpected coordinate/allele row in pinned subset')
        rsid = keys[key]
        if rsid in output:
            raise ValueError('Duplicate exact VCF target')
        if fields[6] not in ('PASS', '.'):
            raise ValueError('Nonpassing population variant filter')
        fmt = fields[8].split(':')
        if 'GT' not in fmt:
            raise ValueError('No GT field')
        gt_index = fmt.index('GT')
        dosages = []
        for index in indices:
            alleles = fields[index + 9].split(':')[gt_index].replace('|', '/').split('/')
            if len(alleles) != 2 or any(a not in ('0', '1', '.') for a in alleles):
                raise ValueError('Expected diploid biallelic GT')
            dosages.append(np.nan if '.' in alleles else sum(map(int, alleles)))
        output[rsid] = np.array(dosages)
        source_ids[rsid] = fields[2]
    return output, source_ids


def choose_pair(ordered, r2, maximum):
    for first, second in itertools.combinations(ordered, 2):
        value, n = r2(first['rsid'], second['rsid'])
        if value is not None and value <= maximum:
            return first, second, value, n
    return None


def load_references(label, manifest):
    result = {}
    for build in ('GRCh37', 'GRCh38', 'model'):
        spec = next(s for s in manifest['sources'] if s['id'] == f'match-{build}-{label}-reference')
        path = ROOT / 'data/raw' / spec['file']
        if hashlib.sha256(path.read_bytes()).hexdigest() != spec['sha256']:
            raise ValueError('Reference changed')
        sequence = (''.join(path.read_text().splitlines()) if build == 'model' else json.loads(path.read_text())['seq']).upper()
        region = spec['region']
        if len(sequence) != region['end0'] - region['start0'] or set(sequence) - set('ACGTN'):
            raise ValueError('Incorrect reference interval length or alphabet')
        result[build] = region['start0'], sequence
    return result


def main():
    protocol = json.loads((ROOT / 'config/psc-controlled-comparison.json').read_text())
    limits = protocol['comparators']
    manifest = json.loads((ROOT / 'config/psc-comparator-sources.json').read_text())
    targets = list(csv.DictReader(io.StringIO((ROOT / 'data/derived/psc-comparator-genotype-targets.csv').read_text())))
    for row in targets:
        for field in ('position_grch37_1based', 'position_grch38_1based', 'nearest_tss_distance_bp'):
            row[field] = int(row[field])
    panel = next(s for s in manifest['sources'] if s['file'] == 'junction-kgp-phase3-panel.data')
    panel_path = ROOT / 'data/raw' / panel['file']
    if hashlib.sha256(panel_path.read_bytes()).hexdigest() != panel['sha256']:
        raise ValueError('Population panel changed')
    samples = list(csv.DictReader(io.StringIO(panel_path.read_text()), delimiter='\t'))
    eur_samples = sorted(s['sample'] for s in samples if s['super_pop'] == 'EUR')
    if len(eur_samples) != len(set(eur_samples)):
        raise ValueError('Duplicate EUR sample in panel')
    query_run = json.loads((ROOT / 'config/psc-comparator-genotype-queries.json').read_text())
    genotypes, source_ids, genotyped_rows = {}, {}, []
    for label in protocol['locus_selection']['selected_archive_labels']:
        rows = [r for r in targets if r['archive_region'] == label]
        query = next(q for q in query_run['queries'] if q['archive_region'] == label)
        path = ROOT / query['file']
        if hashlib.sha256(path.read_bytes()).hexdigest() != query['sha256']:
            raise ValueError('Genotype subset changed')
        text = gzip.decompress(path.read_bytes()).decode()
        if hashlib.sha256(text.encode()).hexdigest() != query['uncompressed_sha256']:
            raise ValueError('Uncompressed genotype subset changed')
        g, ids = parse_vcf(text, eur_samples, rows)
        genotypes.update(g)
        source_ids.update(ids)
        references = load_references(label, manifest)
        for row in rows:
            bases = []
            for build, (start0, sequence) in references.items():
                position = row['position_grch37_1based'] if build == 'GRCh37' else row['position_grch38_1based']
                index = position - 1 - start0
                if not 0 <= index < len(sequence):
                    raise ValueError('Variant outside requested reference')
                bases.append(sequence[index])
            ref_valid = all(base == row['reference_bases'] for base in bases)
            frequency, called = maf(g.get(row['rsid'], []))
            row.update(reference_bases_verified=ref_valid, eur_maf=frequency, eur_called_donors=called,
                       source_vcf_id=ids.get(row['rsid'], ''), genotype_status='exact_alleles_present' if row['rsid'] in g else 'exact_alleles_absent')
            genotyped_rows.append(row)
    ld_cache = {}
    def r2(first, second):
        key = tuple(sorted((first, second)))
        if key not in ld_cache:
            if first not in genotypes or second not in genotypes:
                ld_cache[key] = None, 0
            else:
                ld_cache[key] = dosage_r2(genotypes[first], genotypes[second], limits['minimum_jointly_called_eur_donors_for_ld'])
        return ld_cache[key]
    matches, anchor_results, pair_audit = [], [], []
    for label in protocol['locus_selection']['selected_archive_labels']:
        regional = [r for r in genotyped_rows if r['archive_region'] == label]
        high = [r for r in regional if r['role'] != 'comparator_pool']
        anchors = sorted([r for r in high if r['role'] == 'anchor_pending_comparators'], key=lambda r: (-float(r['pip']), int(r['rsid'][2:])))
        assigned = []
        for anchor in anchors:
            result = {'archive_region': label, 'anchor': anchor['rsid'], 'status': '', 'eur_maf': anchor['eur_maf'],
                      'eur_called_donors': anchor['eur_called_donors'], 'eligible_comparators_at_assignment': 0,
                      'comparator_1': '', 'comparator_2': '', 'comparator_pair_r2': '', 'joint_donors': ''}
            if not anchor['reference_bases_verified'] or anchor['eur_maf'] is None or anchor['eur_maf'] < limits['minimum_eur_maf']:
                result['status'] = 'unmatched: anchor reference, genotype or MAF requirement failed'
                anchor_results.append(result)
                continue
            eligible = []
            for candidate in regional:
                if candidate['role'] != 'comparator_pool':
                    continue
                reason = ''
                if candidate['rsid'] in assigned:
                    reason = 'already_assigned'
                elif not candidate['reference_bases_verified']:
                    reason = 'reference_mismatch'
                elif candidate['eur_maf'] is None:
                    reason = 'missing_exact_genotype'
                elif candidate['eur_maf'] < limits['minimum_eur_maf'] or abs(candidate['eur_maf'] - anchor['eur_maf']) > limits['maximum_absolute_maf_difference']:
                    reason = 'maf_mismatch'
                elif not static_pair(anchor, candidate, limits):
                    reason = 'mutation_context_tss_or_distance_mismatch'
                checks = []
                if not reason:
                    checks = [(other['rsid'], *r2(candidate['rsid'], other['rsid'])) for other in high]
                    checks += [(other, *r2(candidate['rsid'], other)) for other in assigned]
                    if any(value is None for _, value, _ in checks):
                        reason = 'required_ld_unavailable'
                    elif any(value > limits['maximum_ld_r2'] for _, value, _ in checks):
                        reason = 'ld_exceeds_threshold'
                audit = {'archive_region': label, 'anchor': anchor['rsid'], 'comparator': candidate['rsid'],
                         'status': reason or 'eligible', 'matching_cost': '', 'maximum_screen_r2': '', 'minimum_screen_joint_donors': ''}
                if not reason:
                    frequency_difference = abs(candidate['eur_maf'] - anchor['eur_maf'])
                    tss_difference = abs(math.log10(candidate['nearest_tss_distance_bp'] + 1) - math.log10(anchor['nearest_tss_distance_bp'] + 1))
                    physical_distance = abs(candidate['position_grch38_1based'] - anchor['position_grch38_1based'])
                    cost = frequency_difference / .05 + tss_difference / .5 + physical_distance / 250000
                    audit.update(matching_cost=cost, maximum_screen_r2=max(v for _, v, _ in checks), minimum_screen_joint_donors=min(n for _, _, n in checks))
                    eligible.append({**candidate, 'matching_cost': cost, 'maf_difference': frequency_difference,
                                     'log10_tss_distance_difference': tss_difference, 'physical_distance_bp': physical_distance,
                                     'maximum_screen_r2': audit['maximum_screen_r2'], 'minimum_screen_joint_donors': audit['minimum_screen_joint_donors']})
                pair_audit.append(audit)
            eligible.sort(key=lambda r: (r['matching_cost'], int(r['rsid'][2:])))
            result['eligible_comparators_at_assignment'] = len(eligible)
            chosen = choose_pair(eligible, r2, limits['maximum_ld_r2'])
            if not chosen:
                result['status'] = 'unmatched: no mutually eligible comparator pair'
            else:
                first, second, pair_r2, n = chosen
                assigned.extend([first['rsid'], second['rsid']])
                result.update(status='matched', comparator_1=first['rsid'], comparator_2=second['rsid'], comparator_pair_r2=pair_r2, joint_donors=n)
                for order, row in enumerate((first, second), 1):
                    matches.append({'archive_region': label, 'anchor': anchor['rsid'], 'comparator_order': order, **row})
            anchor_results.append(result)
    outputs = {'data/derived/psc-comparator-covariates.csv': genotyped_rows,
               'data/derived/psc-comparator-pair-audit.csv': pair_audit,
               'data/derived/psc-comparator-anchor-status.csv': anchor_results}
    if matches:
        outputs['data/derived/psc-comparator-matches.csv'] = matches
    for name, rows in outputs.items():
        (ROOT / name).write_bytes(csv_bytes(rows))
    ld_rows = [{'variant_1': a, 'variant_2': b, 'eur_dosage_r2': value, 'jointly_called_eur_donors': n} for (a, b), (value, n) in sorted(ld_cache.items())]
    if ld_rows:
        (ROOT / 'data/derived/psc-comparator-ld-audit.csv').write_bytes(csv_bytes(ld_rows))
    summary = {'eur_reference_donors': len(eur_samples), 'genotype_targets': len(targets), 'exact_genotypes_found': len(genotypes),
               'reference_verified_targets': sum(r['reference_bases_verified'] for r in genotyped_rows),
               'anchors_matched': sum(r['status'] == 'matched' for r in anchor_results), 'anchors_total': len(anchor_results),
               'comparators_assigned': len(matches), 'maximum_model_requests_if_complete': len(matches) + sum(r['status'] == 'matched' for r in anchor_results),
               'model_predictions_requested': 0, 'anchor_results': anchor_results,
               'ld_screen_limit': 'Excluded high-PIP anonymous identifiers, indels and multiallelic records remain outside this screen.'}
    (ROOT / 'reports/psc-comparator-matching.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
