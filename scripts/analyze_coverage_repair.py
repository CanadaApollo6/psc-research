"""Replay exact deletion matching and missing-variant checks from pinned inputs."""

import csv
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pysam.bcftools

from coloc_common import ROOT, ChainMap, save_json, sha
from coverage_repair_common import (Reference, component_weights, coverage_parts,
                                    exact_qtl_status, normalize, reciprocal_edit,
                                    unordered_edit_candidates)
from prepare_coloc_inputs import canonical_key, save_frame, unique_alias_rows


def read_json(name): return json.loads((ROOT / name).read_text())


def source_rows(name, delimiter='\t'):
    with gzip.open(ROOT / name, 'rt') as stream:
        if delimiter == ' ':
            lines = stream.read().splitlines()
            return [dict(zip(lines[0].lstrip('#').split(), line.split())) for line in lines[1:]]
        return list(csv.DictReader(stream, delimiter=delimiter))


def verify_pins():
    plan = read_json('config/coverage-repair-plan.json')
    for name, digest in plan['frozen_file_sha256'].items():
        if sha((ROOT / name).read_bytes()) != digest: raise ValueError('Frozen input changed: ' + name)
    for stage in ['reference', 'catalogue']:
        run = read_json(f'config/coverage-repair-{stage}-queries.json')
        if run['plan_sha256'] != sha((ROOT / 'config/coverage-repair-plan.json').read_bytes()):
            raise ValueError('Query plan differs')
        for query in run['queries']:
            data = (ROOT / query['file']).read_bytes()
            if sha(data) != query['sha256'] or sha(gzip.decompress(data)) != query['uncompressed_sha256']:
                raise ValueError('Changed query subset')
    sources = {s['id']: s for s in read_json('config/coverage-repair-sources.json')['sources']}
    for spec in sources.values():
        if sha((ROOT / 'data/raw' / spec['file']).read_bytes()) != spec['sha256']:
            raise ValueError('Changed supporting source: ' + spec['id'])
    return plan, sources


def load_references(sources):
    result = {}
    for gene, chromosome in [('PFKFB3', 'chr10'), ('BCL2L11', 'chr2')]:
        for build in ['GRCh37', 'GRCh38']:
            spec = sources[f'repair-reference-{gene}-{build}']
            sequence = (ROOT / 'data/raw' / spec['file']).read_text().strip().upper()
            if set(sequence) - set('ACGT') or len(sequence) != spec['end_1based'] - spec['start_1based'] + 1:
                raise ValueError('Invalid or truncated reference sequence')
            result[gene, build] = Reference(chromosome, spec['start_1based'], sequence)
    return result


def vcf_aggregates(gene, reference):
    panel = pd.read_csv(ROOT / 'data/raw/junction-kgp-phase3-panel.data', sep='\t')
    eur = set(panel.loc[panel.super_pop.eq('EUR'), 'sample'])
    if len(eur) != 503: raise ValueError('Unexpected EUR reference panel')
    records = []
    with gzip.open(ROOT / f'data/raw/repair-reference-{gene}.vcf.gz', 'rt') as stream:
        for line in stream:
            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                samples = [i for i, sample in enumerate(columns) if sample in eur]
                if len(samples) != 503: raise ValueError('Missing public reference samples')
            elif not line.startswith('#'):
                f = line.strip().split('\t')
                if f[3] == '.' or any(set(a) - set('ACGT') for a in [f[3], *f[4].split(',')]): continue
                gt_index = f[8].split(':').index('GT')
                genotypes = [f[i].split(':')[gt_index].replace('|', '/').split('/') for i in samples]
                called = [g for g in genotypes if len(g) == 2 and '.' not in g]
                for j, alt in enumerate(f[4].split(','), 1):
                    try: normalized = normalize(reference, int(f[1]), f[3], alt)
                    except ValueError: continue
                    ac = sum(g.count(str(j)) for g in called)
                    records.append(dict(position_grch37=int(f[1]), source_ref=f[3], source_alt=alt,
                                        source_id=f[2], normalized=normalized, called_donors=len(called),
                                        alt_count=ac, allele_number=2 * len(called), alt_frequency=ac / (2 * len(called))))
    return records


def bcftools_check(reference, edits, label):
    folder = ROOT / 'work/coverage-repair-normalization'
    folder.mkdir(exist_ok=True)
    fasta = folder / (label + '.fa')
    fasta.write_text('>block\n' + reference.sequence + '\n')
    pysam.faidx(str(fasta))
    vcf = folder / (label + '.vcf')
    lines = ['##fileformat=VCFv4.2', f'##contig=<ID=block,length={len(reference.sequence)}>',
             '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
    for identifier, (pos, ref, alt) in edits:
        lines.append(f'block\t{pos - reference.start + 1}\t{identifier}\t{ref}\t{alt}\t.\tPASS\t.')
    vcf.write_text('\n'.join(lines) + '\n')
    output = pysam.bcftools.norm('-f', str(fasta), '-c', 'e', '-Ov', str(vcf), catch_stdout=True)
    checks = []
    expected = {name: normalize(reference, *edit) for name, edit in edits}
    for line in output.splitlines():
        if line.startswith('#'): continue
        f = line.split('\t'); result = (int(f[1]) + reference.start - 1, f[3], f[4])
        if result != expected[f[2]]: raise ValueError('Independent bcftools normalization disagrees')
        checks.append(dict(id=f[2], build=label, bcftools_matches=True, normalized=result))
    if len(checks) != len(edits): raise ValueError('Incomplete independent normalization check')
    return checks


def deletion_analysis(plan, references, forward, reverse):
    ref38, ref37 = references['PFKFB3', 'GRCh38'], references['PFKFB3', 'GRCh37']
    gwas = source_rows('data/derived/coloc-query-rows/GWAS-PFKFB3.txt.gz', ' ')
    public = vcf_aggregates('PFKFB3', ref37)
    matches, audit, edits, dbsnp = [], [], [], []
    for row in gwas:
        if len(row['allele_0']) == len(row['allele_1']) == 1: continue
        position, length = int(row['pos']), max(len(row['allele_0']), len(row['allele_1']))
        if not ref37.start <= position <= ref37.start + len(ref37.sequence) - length: continue
        candidates = unordered_edit_candidates(ref37, position, row['allele_0'], row['allele_1'])
        audit.append({**row, 'possible_reference_edits': candidates})
    for t, source_id in zip(plan['pfkfb3_deletions'], ['rs149556258', 'rs138248884', 'rs140499372']):
        edit38 = (t['position_grch38'], t['ref'], t['alt'])
        mapped, bases, strand = reciprocal_edit(ref38, ref37, forward, reverse, *edit38)
        norm38, norm37 = normalize(ref38, *edit38), normalize(ref37, *mapped)
        edits.append((t['id'], mapped))
        records = [r for r in audit if norm37 in r['possible_reference_edits']]
        refs = [r for r in public if r['normalized'] == norm37]
        if len(records) > 1 or len(refs) > 1: raise ValueError('Ambiguous exact deletion mapping')
        # Explicit phase3 REF/ALT and matching population frequency corroborate
        # the deletion interpretation of otherwise unordered GWAS allele strings.
        freq_difference = None
        if records and refs:
            row, rr = records[0], refs[0]
            if int(row['pos']) != rr['position_grch37'] or {row['allele_0'], row['allele_1']} != {rr['source_ref'], rr['source_alt']}:
                raise ValueError('GWAS and explicit reference record require additional representation audit')
            risk_freq = rr['alt_frequency'] if row['allele_1'] == rr['source_alt'] else 1 - rr['alt_frequency']
            freq_difference = abs(risk_freq - float(row['freq_1_controls']))
            if freq_difference > .15: raise ValueError('Deletion frequency orientation is inconsistent')
        # Exhaustively find all minimal anchored deletion placements within the
        # sequence block. A different repeat placement cannot escape this tract.
        placements = []
        delta = len(t['ref']) - len(t['alt'])
        for pos in range(ref37.start + 1, ref37.start + len(ref37.sequence) - delta):
            ref = ref37.get(pos, delta + 1)
            if ref37.haplotype(pos, ref, ref[0]) == ref37.haplotype(*norm37): placements.append(pos)
        if not placements or min(placements) <= ref37.start + 1 or max(placements) + delta >= ref37.start + len(ref37.sequence) - 1:
            raise ValueError('Repeat-equivalence tract not bounded by reference block')
        # Inspect every original regional GWAS span that could encode this edit,
        # including long padded records starting outside the queried block.
        potential = [r for r in gwas if abs(len(r['allele_0']) - len(r['allele_1'])) == delta
                     and int(r['pos']) <= max(placements) + delta
                     and int(r['pos']) + max(len(r['allele_0']), len(r['allele_1'])) - 1 >= min(placements)]
        if any(r['SNP'] not in {a['SNP'] for a in audit} for r in potential):
            raise ValueError('Potential equivalent GWAS edit lies beyond reference block')
        item = dict(target_id=t['id'], source_variant=f"chr10_{edit38[0]}_{edit38[1]}_{edit38[2]}",
                    deleted_bases=delta, normalized_grch38=norm38, normalized_grch37=norm37,
                    reciprocal_bases=bases, strand=strand, full_flank_sequence_agrees=True,
                    minimal_equivalent_placement_min_grch37=min(placements), minimal_equivalent_placement_max_grch37=max(placements),
                    original_regional_gwas_rows=len(gwas), bounded_gwas_indels_checked=len(audit), potential_equivalent_gwas_rows=len(potential),
                    gwas_matches=[{k: v for k, v in r.items() if k != 'possible_reference_edits'} for r in records],
                    reference_matches=refs, control_vs_reference_risk_frequency_difference=freq_difference,
                    recoverable=bool(records and refs), status='sequence_matched_GWAS_and_reference' if records and refs else 'absent_from_original_regional_GWAS_and_phase3_reference' if not records and not refs else 'incomplete_cross_source_support')
        matches.append(item)
        for build, reference, target in [('GRCh38', ref38, norm38), ('GRCh37', ref37, norm37)]:
            data = read_json(f'data/raw/repair-variant-{source_id}-{build}.json')
            for m in data.get('mappings', []):
                if m['seq_region_name'] != '10' or m['assembly_name'] != build: continue
                alleles = m['allele_string'].split('/')
                exact = []
                for alt in alleles[1:]:
                    try:
                        if normalize(reference, m['start'], alleles[0], alt) == target: exact.append(alt)
                    except ValueError: continue
                dbsnp.append(dict(target_id=t['id'], queried_rsid=source_id, returned_rsid=data['name'], build=build,
                                  mapping=m, exact_target_alt_representations=exact))
    independent = bcftools_check(ref37, edits, 'GRCh37') + bcftools_check(ref38, [(t['id'], (t['position_grch38'], t['ref'], t['alt'])) for t in plan['pfkfb3_deletions']], 'GRCh38')
    return matches, audit, dbsnp, independent


def signal_coverage(matches):
    previous = pd.read_csv(ROOT / 'reports/coloc-published-signal-overlap.csv')
    shared_keys = set(pd.read_csv(ROOT / 'data/derived/coloc-platform-inputs/GWAS-PFKFB3-ld.csv.gz').snp)
    recovered = {m['source_variant'] for m in matches if m['recoverable']}
    rows, variants = [], []
    for ds, components in [('QTD000031', [1, 2]), ('QTD000469', [1])]:
        frame = pd.read_csv(ROOT / f'data/derived/coloc-query-rows/{ds}-PFKFB3-published-bfs.tsv.gz', sep='\t')
        if frame.variant.duplicated().any(): raise ValueError('Duplicate full BF variant')
        keys = []
        for v in frame.variant:
            chrom, pos, ref, alt = v.split('_')
            try: keys.append(canonical_key(chrom, int(pos), ref, alt)[0])
            except ValueError: keys.append(None)
        baseline = np.asarray([k in shared_keys for k in keys])
        added = frame.variant.isin(recovered).to_numpy()
        for component in components:
            weights = component_weights(frame[f'lbf_variable{component}'])
            old, addition, total = coverage_parts(weights, baseline, added)
            prior = previous[previous.dataset_id.eq(ds) & previous.qtl_component.eq(component)]
            if not np.allclose(prior.qtl_signal_mass_retained, old, atol=1e-12, rtol=0): raise ValueError('Baseline signal coverage does not replay')
            rows.append(dict(dataset_id=ds, component=component, full_source_variants=len(frame),
                             baseline_shared_variants=int(baseline.sum()), exact_added_variants=int(added.sum()),
                             original_signal_coverage=old, added_signal_weight=addition, recoverable_signal_coverage=total,
                             remaining_signal_weight=1-total, qtl_90pct_gate=total >= .9,
                             status='support_audit_only_no_new_GWAS_fit_or_posterior'))
            for m in matches:
                one = frame.variant.eq(m['source_variant'])
                if one.sum() != 1: raise ValueError('Frozen deletion absent from full molecular BF vector')
                variants.append(dict(dataset_id=ds, component=component, target_id=m['target_id'], source_variant=m['source_variant'],
                                     component_bf_weight=float(weights[one].sum()), exact_gwas_recoverable=m['recoverable']))
    return rows, variants


def bcl_missingness(plan, references, forward, reverse):
    ref38, ref37 = references['BCL2L11', 'GRCh38'], references['BCL2L11', 'GRCh37']
    public = vcf_aggregates('BCL2L11', ref37)
    gwas = source_rows('data/derived/coloc-query-rows/GWAS-BCL2L11.txt.gz', ' ')
    exact_targets, statuses, effects = [], [], []
    for t in plan['bcl2l11_variants']:
        ref = ref38.get(t['position_grch38']); alt = next(a for a in t['alleles'] if a != ref)
        if ref not in t['alleles']: raise ValueError('BCL target REF mismatch')
        edit, bases, strand = reciprocal_edit(ref38, ref37, forward, reverse, t['position_grch38'], ref, alt)
        source = [r for r in gwas if int(r['pos']) == edit[0] and {r['allele_0'], r['allele_1']} == {edit[1], edit[2]}]
        refs = [r for r in public if r['normalized'] == edit]
        if len(source) != 1 or len(refs) != 1: raise ValueError('BCL identity not unique')
        exact_targets.append(dict(**t, ref=ref, alt=alt, mapped_grch37=edit, reciprocal_bases=bases,
                                  full_flank_sequence_agrees=True, gwas=source[0], reference=refs[0]))
        for ds in plan['datasets']:
            rows = source_rows(f"data/derived/coverage-repair-queries/{ds['dataset_id']}-catalogue.tsv.gz")
            status, at_position, exact, measured = exact_qtl_status(rows, t['position_grch38'], ref, alt, 'ENSG00000153094')
            statuses.append(dict(dataset_id=ds['dataset_id'], sample_group=ds['sample_group'], source_id=t['source_id'],
                                 position_grch38=t['position_grch38'], ref=ref, alt=alt, status=status,
                                 all_gene_position_rows=len(at_position), exact_allele_rows=len(exact), target_gene_rows=len(measured),
                                 observed_alleles=';'.join(sorted({r['ref']+'>'+r['alt'] for r in at_position}))))
            if measured:
                row, aliases, rsids = unique_alias_rows(measured)
                effects.append(dict(dataset_id=ds['dataset_id'], source_id=t['source_id'], ref=ref, alt=alt,
                                    beta_per_alt=float(row['beta']), se=float(row['se']), pvalue=float(row['pvalue']),
                                    maf=float(row['maf']), ma_samples=int(row['ma_samples']), an=int(row['an']), aliases=aliases,
                                    already_present_in_previous_gene_query=True))
                original = source_rows(f"data/derived/coloc-query-rows/{ds['dataset_id']}-BCL2L11.tsv.gz")
                if not any(r == row for r in original): raise ValueError('All-gene result differs from frozen original query')
    return exact_targets, statuses, effects


def main():
    plan, sources = verify_pins()
    references = load_references(sources)
    forward = ChainMap.from_file(ROOT / 'data/raw/coloc-hg38-to-hg19-chain.gz')
    reverse = ChainMap.from_file(ROOT / 'data/raw/coloc-hg19-to-hg38-chain.gz')
    matches, audit, dbsnp, independent = deletion_analysis(plan, references, forward, reverse)
    coverage, weights = signal_coverage(matches)
    targets, statuses, effects = bcl_missingness(plan, references, forward, reverse)
    folder = ROOT / 'data/derived'
    save_json(folder / 'coverage-repair-exact-identities.json', dict(pfkfb3=matches, bcl2l11=targets, ensembl_representation_audit=dbsnp,
              independent_normalization=independent, pysam_version=pysam.__version__, htslib_version=pysam.__samtools_version__))
    save_json(folder / 'coverage-repair-bounded-gwas-indels.json', audit)
    for name, rows in [('pfkfb3-signal-coverage', coverage), ('pfkfb3-deletion-weights', weights),
                       ('bcl2l11-missingness', statuses), ('bcl2l11-blueprint-effects', effects)]:
        save_frame(pd.DataFrame(rows), folder / ('coverage-repair-' + name + '.csv'))
    print(json.dumps({'pfkfb3_exact_recoverable_deletions': sum(m['recoverable'] for m in matches),
                      'pfkfb3_signal_coverage': coverage, 'bcl2l11_statuses': dict(Counter(r['status'] for r in statuses))}))


if __name__ == '__main__': main()
