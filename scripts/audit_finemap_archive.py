"""Audit every published FINEMAP triplet without inferring missing signal sets."""

import csv
import gzip
import hashlib
import io
import json
import re
import tarfile
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path, PurePosixPath

from fetch_sources import validate

ROOT = Path(__file__).resolve().parents[1]
ZERO = Decimal(0)
PIP_TOLERANCE = Decimal('0.0000501')  # Four-decimal SNP values, plus config rounding.
MASS_TOLERANCE = Decimal('0.00001')
LOG_TOLERANCE = Decimal('0.000005')


def number(value, probability=False):
    result = Decimal(value)
    if not result.is_finite() or (probability and not ZERO <= result <= 1):
        raise ValueError(f'Invalid numeric field: {value}')
    return result


def log_bayes_factor(value):
    # Zero Bayes factors can be serialized as -inf; preserve that source value.
    if value != '-inf':
        number(value)


def identifier_type(identifier):
    if re.fullmatch(r'rs\d+', identifier):
        return 'rsid'
    if re.fullmatch(r'rs\d+(;rs\d+)+', identifier):
        return 'semicolon_rsid_aliases'
    if re.fullmatch(r'chr[0-9XY]+:\d+', identifier):
        return 'coordinate_without_build_or_alleles'
    if re.fullmatch(r'kgp\d+', identifier):
        return 'kgp_identifier'
    if identifier.startswith('SNP_A-'):
        return 'array_probe_identifier'
    return 'other_unresolved_identifier'


def safe_members(archive):
    """Read regular files in memory; never extract paths or follow links."""
    inventory, content, seen = [], {}, set()
    total = 0
    for member in archive.getmembers():
        path = PurePosixPath(member.name)
        if path.is_absolute() or '..' in path.parts or member.name in seen:
            raise ValueError('Unsafe or duplicate archive member')
        seen.add(member.name)
        if not (member.isfile() or member.isdir()):
            raise ValueError('Only regular files and directories are accepted')
        total += member.size
        if member.size > 32 * 1024**2 or total > 128 * 1024**2:
            raise ValueError('Uncompressed archive limit exceeded')
        data = archive.extractfile(member).read() if member.isfile() else b''
        inventory.append({'member': member.name, 'kind': 'file' if member.isfile() else 'directory',
                          'bytes': member.size, 'sha256': hashlib.sha256(data).hexdigest() if member.isfile() else ''})
        if member.isfile():
            content[member.name] = data.decode('utf-8')
    return sorted(inventory, key=lambda x: x['member']), content


def parse_snps(text):
    lines = text.splitlines()
    if not lines or lines[0].split() != ['index', 'snp', 'snp_prob', 'snp_log10bf']:
        raise ValueError('Unexpected SNP table header')
    rows, identifiers, indices = [], set(), set()
    for ordinal, line in enumerate(lines[1:], 1):
        fields = line.split()
        if len(fields) != 4:
            raise ValueError('Malformed SNP row')
        index, identifier, pip, logbf = fields
        if int(index) < 1 or index in indices or identifier in identifiers:
            raise ValueError('Duplicate SNP or source index')
        number(pip, probability=True)
        log_bayes_factor(logbf)
        indices.add(index)
        identifiers.add(identifier)
        rows.append({'source_row': ordinal, 'source_index': index, 'variant_id_as_published': identifier,
                     'identifier_type': identifier_type(identifier), 'snp_prob_as_published': pip,
                     'snp_log10bf_as_published': logbf})
    if not rows:
        raise ValueError('Empty SNP table')
    return rows


def parse_configs(text, identifiers):
    if not text.strip():
        return None
    lines = text.splitlines()
    if lines[0].split() != ['rank', 'config', 'config_prob', 'config_log10bf']:
        raise ValueError('Unexpected configuration header')
    marginals, k_mass, k_count = defaultdict(Decimal), defaultdict(Decimal), Counter()
    seen, singletons = set(), []
    for ordinal, line in enumerate(lines[1:], 1):
        fields = line.split()
        if len(fields) != 4 or int(fields[0]) != ordinal:
            raise ValueError('Malformed or nonsequential configuration row')
        _, label, prob, logbf = fields
        variants = label.split(',')
        key = tuple(sorted(variants))
        if len(set(variants)) != len(variants) or key in seen or not set(variants) <= identifiers:
            raise ValueError('Duplicate configuration or unknown variant')
        seen.add(key)
        value = number(prob, probability=True)
        log_bayes_factor(logbf)
        k = len(variants)
        k_mass[k] += value
        k_count[k] += 1
        for variant in variants:
            marginals[variant] += value
        if k == 1:
            singletons.append((variants[0], value))
    if not seen:
        raise ValueError('Configuration header without rows')
    return {'rows': len(seen), 'mass': sum(k_mass.values(), ZERO), 'marginals': marginals,
            'k_mass': dict(sorted(k_mass.items())), 'k_count': dict(sorted(k_count.items())),
            'singletons': singletons}


def parse_log(text):
    def field(pattern, cast=str):
        match = re.search(pattern, text)
        if not match:
            raise ValueError(f'Missing FINEMAP log field: {pattern}')
        return cast(match[1])

    posterior_text = text.split('- Post-Pr(# of causal SNPs is k)  :', 1)
    if len(posterior_text) != 2:
        raise ValueError('Missing posterior distribution in log')
    posterior = {int(k): number(value, probability=True)
                 for k, value in re.findall(r'^\s*(\d+) -> ([\d.eE+-]+)\s*$', posterior_text[1], re.M)}
    return {'finemap_version': field(r'Welcome to FINEMAP v([\d.]+)'),
            'max_causal_variants': field(r'n-causal-max\s*:\s*(\d+)', int),
            'number_of_snps': field(r'Number of SNPs in region\s*:\s*(\d+)', int),
            'number_of_individuals': field(r'Number of individuals in GWAS\s*:\s*(\d+)', int),
            'posterior_k': posterior}


def singleton_sets(config, snp_count, thresholds=(Decimal('.95'), Decimal('.99'))):
    """Only complete singleton configurations define a one-causal-variant set."""
    if (config is None or set(config['k_count']) != {1} or config['rows'] != snp_count
            or abs(config['mass'] - 1) > MASS_TOLERANCE):
        return []
    ordered = sorted(config['singletons'], key=lambda item: (-item[1], item[0]))
    output = []
    for threshold in thresholds:
        cumulative, included, boundary = ZERO, [], None
        for identifier, probability in ordered:
            if cumulative >= threshold and probability != boundary:
                break
            cumulative += probability
            boundary = probability
            included.append({'variant_id_as_published': identifier, 'config_probability': str(probability),
                             'cumulative_probability': str(cumulative)})
        if cumulative < threshold:
            raise ValueError('Insufficient configuration mass for requested set')
        output.append({'threshold': str(threshold), 'size': len(included), 'mass': str(cumulative),
                       'members': included})
    return output


def audit_dataset(stem, contents):
    rows = parse_snps(contents[stem + '.snp'])
    config = parse_configs(contents[stem + '.config'], {r['variant_id_as_published'] for r in rows})
    log = parse_log(contents[stem + '.log'])
    kind = 'gwas' if '/psc_gwas_fine_mapping/' in stem else 'molecular_qtl'
    label = PurePosixPath(stem).name
    dataset = f'{kind}/{label}'
    issues = []
    if log['number_of_snps'] != len(rows):
        issues.append('log_snp_count_disagrees')
    if log['finemap_version'] != '1.3':
        issues.append('log_version_differs_from_paper_methods_v1.3')
    audit = {'dataset_id': dataset, 'analysis_kind': kind, 'archive_label': label, 'archive_stem': stem,
             'snp_rows': len(rows), 'identifier_types': dict(sorted(Counter(r['identifier_type'] for r in rows).items())),
             'rounded_pip_sum': str(sum((Decimal(r['snp_prob_as_published']) for r in rows), ZERO)),
             'rounded_zero_pip_rows': sum(Decimal(r['snp_prob_as_published']) == 0 for r in rows),
             'negative_infinite_log_bayes_factors': sum(r['snp_log10bf_as_published'] == '-inf' for r in rows),
             'log': {**log, 'posterior_k': {str(k): str(p) for k, p in log['posterior_k'].items()}},
             'config_rows': None, 'config_probability_mass': None, 'config_k_counts': None,
             'config_k_probability_mass': None, 'config_max_causal_variants': None,
             'max_pip_config_difference': None, 'max_posterior_k_log_difference': None,
             'single_causal_reconstructed_sets': [], 'issues': issues}
    if config is None:
        issues.append('empty_configuration_file')
    else:
        if abs(config['mass'] - 1) > MASS_TOLERANCE:
            issues.append('saved_configuration_mass_not_near_one')
        gaps = [abs(Decimal(r['snp_prob_as_published']) - config['marginals'].get(r['variant_id_as_published'], ZERO)) for r in rows]
        posterior_gap = max(abs(config['k_mass'].get(k, ZERO) - log['posterior_k'].get(k, ZERO))
                            for k in set(config['k_mass']) | set(log['posterior_k']))
        audit.update({'config_rows': config['rows'], 'config_probability_mass': str(config['mass']),
                      'config_k_counts': config['k_count'],
                      'config_k_probability_mass': {str(k): str(v) for k, v in config['k_mass'].items()},
                      'config_max_causal_variants': max(config['k_count']),
                      'max_pip_config_difference': str(max(gaps)),
                      'max_posterior_k_log_difference': str(posterior_gap)})
        if max(gaps) > PIP_TOLERANCE:
            issues.append('configuration_marginals_disagree_with_snp_pips')
        if max(config['k_count']) > log['max_causal_variants']:
            issues.append('configuration_exceeds_log_max_causal_variants')
        if posterior_gap > LOG_TOLERANCE:
            issues.append('configuration_posterior_k_disagrees_with_log')
        audit['single_causal_reconstructed_sets'] = singleton_sets(config, len(rows))
    normalized = [{
        'dataset_id': dataset, **r,
        'config_summed_marginal_probability': (str(config['marginals'].get(r['variant_id_as_published'], ZERO))
                                               if config else ''),
    } for r in rows]
    audit['top_variants'] = normalized[:5]
    return audit, normalized


def crosswalk(signals, all_rows):
    """Relate region labels, never silently rename genes or invent rsID aliases."""
    label_map = {'BCL211': 'BCL2L11', 'ETS2': 'PSMG1', 'IL2-IL21': 'IL21'}
    regions = defaultdict(list)
    for row in all_rows:
        regions[row['dataset_id']].append(row)
    result = []
    for signal in signals:
        published_label = signal['candidate_gene_as_published']
        region = label_map.get(published_label, published_label)
        rows = regions['gwas/' + region]
        target = signal['highest_pp_snp']
        coordinate = f"chr{signal['chromosome']}:{signal['highest_pp_position_b37']}"
        matches = [r for r in rows if target in r['variant_id_as_published'].split(';')]
        method = 'exact_rsid' if matches and matches[0]['variant_id_as_published'] == target else 'rsid_in_published_alias_list'
        if not matches:
            matches = [r for r in rows if r['variant_id_as_published'] == coordinate]
            method = 'coordinate_string_candidate_build_and_alleles_unverified' if matches else 'not_resolved'
        if len(matches) > 1:
            raise ValueError('Ambiguous Table 1 crosswalk')
        match = matches[0] if matches else {}
        result.append({
            'source_table_data_row': signal['source_table_data_row'], 'published_region_label': published_label,
            'archive_region_label': region, 'signal_as_published': signal['signal'],
            'highest_pp_snp_as_published': target, 'highest_pp_position_b37_as_published': signal['highest_pp_position_b37'],
            'posterior_probability_as_published': signal['posterior_probability_as_published'],
            'credible_set_size_as_published': signal['credible_set_size_as_published'],
            'archive_variant_identifier': match.get('variant_id_as_published', ''), 'match_method': method,
            'archive_marginal_pip': match.get('snp_prob_as_published', ''),
            'complete_published_signal_membership_available': False,
            'region_link_note': ('Region association only; gene labels are not declared equivalent'
                                 if published_label != region else 'Literal region-label agreement'),
        })
    return result


def csv_bytes(rows):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def main():
    manifest = json.loads((ROOT / 'config/finemap-sources.json').read_text())
    for source in manifest['sources']:
        validate((ROOT / 'data/raw' / source['file']).read_bytes(), source)
    source = next(s for s in manifest['sources'] if s['id'] == 'psc-finemap-archive')
    with tarfile.open(ROOT / 'data/raw' / source['file'], 'r:gz') as archive:
        inventory, contents = safe_members(archive)
    stems = sorted(name[:-4] for name in contents if name.endswith('.snp'))
    for stem in stems:
        if not all(stem + suffix in contents for suffix in ('.config', '.log')):
            raise ValueError('Incomplete FINEMAP triplet')
    unassigned = [name for name in contents if name.endswith(('.config', '.log')) and name.rsplit('.', 1)[0] not in stems]
    if unassigned:
        raise ValueError('Unassigned configuration or log')
    audits, all_rows, sets = [], [], []
    for stem in stems:
        audit, rows = audit_dataset(stem, contents)
        for entry in audit['single_causal_reconstructed_sets']:
            members = entry.pop('members')
            sets.extend({'dataset_id': audit['dataset_id'], 'threshold': entry['threshold'], 'set_size': entry['size'],
                         'set_probability_mass': entry['mass'], **member} for member in members)
        audits.append(audit)
        all_rows.extend(rows)
    signals_path = ROOT / 'data/derived/psc-signal-summaries.csv'
    signals = list(csv.DictReader(io.StringIO(signals_path.read_text())))
    crosswalk_rows = crosswalk(signals, all_rows)
    report = {
        'analysis_version': '1.0', 'source_manifest': 'config/finemap-sources.json',
        'archive_sha256': source['sha256'], 'archive_publisher_md5': source['publisher_checksum'],
        'source_table1_csv_sha256': hashlib.sha256(signals_path.read_bytes()).hexdigest(),
        'archive_members': len(inventory), 'regular_files': len(contents),
        'gwas_regions': sum(a['analysis_kind'] == 'gwas' for a in audits),
        'molecular_qtl_tables': sum(a['analysis_kind'] == 'molecular_qtl' for a in audits),
        'variant_rows_by_analysis_kind': {kind: sum(a['snp_rows'] for a in audits if a['analysis_kind'] == kind)
                                           for kind in ('gwas', 'molecular_qtl')},
        'nonempty_config_tables': sum(a['config_rows'] is not None for a in audits),
        'pip_config_concordant_tables': sum(a['max_pip_config_difference'] is not None and Decimal(a['max_pip_config_difference']) <= PIP_TOLERANCE for a in audits),
        'log_max_causal_conflicts': sum('configuration_exceeds_log_max_causal_variants' in a['issues'] for a in audits),
        'log_posterior_k_conflicts': sum('configuration_posterior_k_disagrees_with_log' in a['issues'] for a in audits),
        'paper_method_finemap_version': '1.3', 'all_log_versions': sorted({a['log']['finemap_version'] for a in audits}),
        'explicit_credible_set_files': [n for n in contents if '.cred' in n],
        'genetic_input_or_ld_files': [n for n in contents if n.endswith(('.z', '.ld', '.bcor'))],
        'extra_gwas_regions_not_in_table1': sorted({a['archive_label'] for a in audits if a['analysis_kind'] == 'gwas'} - {r['archive_region_label'] for r in crosswalk_rows}),
        'model_predictions_run': 0,
        'interpretation': 'Archived marginal probabilities are not complete per-signal credible-set memberships. Config-summed values are retained separately and never replace source PIPs. Missing build/alleles remain missing.',
        'datasets': audits,
    }
    derived, reports = ROOT / 'data/derived', ROOT / 'reports'
    with (derived / 'psc-finemap-variants.csv.gz').open('wb') as stream:
        with gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0) as output:
            output.write(csv_bytes(all_rows))
    (derived / 'psc-finemap-singleton-reconstructions.csv').write_bytes(csv_bytes(sets))
    (derived / 'psc-finemap-table1-crosswalk.csv').write_bytes(csv_bytes(crosswalk_rows))
    (reports / 'psc-finemap-archive-members.csv').write_bytes(csv_bytes(inventory))
    summary_fields = ['dataset_id', 'snp_rows', 'rounded_pip_sum', 'config_rows', 'config_probability_mass',
                      'config_max_causal_variants', 'max_pip_config_difference', 'max_posterior_k_log_difference']
    summary_rows = [{**{key: a[key] for key in summary_fields},
                     'log_max_causal_variants': a['log']['max_causal_variants'],
                     'log_individuals': a['log']['number_of_individuals'], 'issues': ';'.join(a['issues'])} for a in audits]
    (reports / 'psc-finemap-datasets.csv').write_bytes(csv_bytes(summary_rows))
    (reports / 'psc-finemap-audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key != 'datasets'}, indent=2))


if __name__ == '__main__':
    main()
